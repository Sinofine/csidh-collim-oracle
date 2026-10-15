"""Synthesize every analytic phase target and independently verify each word.

Default budget: full dynamic-circuit channel diamond-distance contribution
from phase synthesis at most 2^-64, by telescoping (not postselected distance).
Both ancilla-free and one-reusable-clean-ancilla lowerings are included.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import time
import mpmath as mp
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from exact_angles import Angle
from phase_certificates import positive_angle, verify_rotation, phase_error
from pygridsynth.gridsynth import gridsynth_gates


CLIFFORD_T = {
    0: [],
    1: ["t"],
    2: ["s"],
    3: ["s", "t"],
    4: ["z"],
    5: ["z", "t"],
    6: ["sdg"],
    7: ["tdg"],
}


def synthesize(angle, bits):
    theta = angle.value(180)
    k = int(mp.nint(theta * 4 / mp.pi)) % 8
    entry = dict(
        tokens=CLIFFORD_T[k],
        phase_frame=dict(angle_half=0, pi_quarters=0),
        method="near_exact_clifford_t",
    )
    # Test an exact Clifford+T candidate against the original analytic angle.
    with ctx.workprec(640):
        snapped = phase_error(angle, entry) < arb(2) ** (-bits)
    if not snapped:
        word = gridsynth_gates(
            theta=theta,
            epsilon=mp.mpf(2) ** (-bits - 2),
            dps=180,
            seed=0,
            up_to_phase=False,
        )
        assert set(word) <= set("HSTXW")
        entry = dict(
            tokens=[x.lower() for x in reversed(word) if x != "W"],
            phase_frame=dict(angle_half=1, pi_quarters=word.count("W")),
            method="ross_selinger",
            gridsynth_matrix_order_word=word,
        )
    entry["t_count"] = sum(t in ("t", "tdg") for t in entry["tokens"])
    entry["frobenius_error_upper_display"] = str(
        verify_rotation(str(angle), entry, bits)
    )
    return entry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-bits", type=int, default=64)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "results/phase_synthesis.json"
    )
    args = parser.parse_args()
    if not 1 <= args.budget_bits <= 128:
        parser.error("Supported per-call budget range: 1 to 128 bits.")
    if args.verify_only:
        data = json.loads(args.output.read_text())
        for k, v in data["rotations"].items():
            verify_rotation(k, v, data["rotation_error_bits"])
        print(
            json.dumps(
                dict(
                    verified_rotations=len(data["rotations"]),
                    rotation_error_bits=data["rotation_error_bits"],
                )
            )
        )
        return
    started = time.monotonic()
    revised = json.loads((ROOT / "results/reference/oracle.json").read_text())
    baseline = json.loads((ROOT / "results/bonnetain/oracle.json").read_text())
    ledger = revised["phase_ledger"]
    assert ledger == baseline["phase_ledger"]
    counts = Counter()
    targets = {}
    for row in ledger:
        counts[row["class"]] += row["count"]
        if row["class"] == "GPhase":
            continue
        angle = Angle.from_key(row["angle_key"])
        for a in [angle] + ([angle / 2] if row["class"] == "CPhase" else []):
            a, _ = positive_angle(a)
            targets[str(a)] = a
    assert dict(counts) == revised["resources"]["unresolved_phase_gates"]
    maximum_rotations = counts["Phase"] + 3 * counts["CPhase"]
    bits = args.budget_bits + (maximum_rotations - 1).bit_length() + 1
    mp.mp.dps = 180
    rotations = {}
    for i, (key, angle) in enumerate(sorted(targets.items())):
        rotations[key] = synthesize(angle, bits)
        if (i + 1) % 50 == 0:
            print(json.dumps(dict(completed=i + 1, total=len(targets))), flush=True)

    def cost(a):
        a, _ = positive_angle(a)
        return rotations[str(a)]["t_count"]

    variants = {}
    for policy in ["ancilla_free", "one_clean_ancilla"]:
        rows = []
        total_t = 0
        total_ccx = 0
        nrot = 0
        for source in ledger:
            if source["class"] == "GPhase":
                continue
            a = Angle.from_key(source["angle_key"])
            row = dict(source)
            if source["class"] == "Phase":
                row.update(
                    method="single", t_each=cost(a), extra_ccx_each=0, rotations_each=1
                )
            else:
                free_cost = 3 * cost(a / 2)
                anc_cost = 14 + cost(a)
                if policy == "one_clean_ancilla" and anc_cost < free_cost:
                    row.update(
                        method="and_phase_uncompute",
                        t_each=anc_cost,
                        extra_ccx_each=2,
                        rotations_each=1,
                    )
                else:
                    row.update(
                        method="three_phase_two_cx",
                        t_each=free_cost,
                        extra_ccx_each=0,
                        rotations_each=3,
                    )
            total_t += row["count"] * row["t_each"]
            total_ccx += row["count"] * row["extra_ccx_each"]
            nrot += row["count"] * row["rotations_each"]
            rows.append(row)
        revised_t = revised["resources"]["ccx_7T_conversion_only"] + total_t
        baseline_t = baseline["resources"]["ccx_7T_conversion_only"] + total_t
        variants[policy] = dict(
            phase_t=total_t,
            extra_ccx=total_ccx,
            standalone_phase_t=total_t - 7 * total_ccx,
            lowerings=rows,
            rotation_occurrences=nrot,
            channel_error_bound=f"{2 * nrot} * 2^(-{bits})",
            revised_t=revised_t,
            baseline_t=baseline_t,
            t_reduction_percent=str((1 - mp.mpf(revised_t) / baseline_t) * 100),
        )
        assert 2 * nrot <= 2 ** (bits - args.budget_bits)
    result = dict(
        gate_counts=dict(counts),
        budget_bits=args.budget_bits,
        rotation_error_bits=bits,
        maximum_rotation_occurrences=maximum_rotations,
        seed=0,
        working_decimal_digits=180,
        verification_binary_digits=640,
        rotations=rotations,
        variants=variants,
        elapsed_seconds=time.monotonic() - started,
    )
    args.output.parent.mkdir(exist_ok=True, parents=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: {kk: vv for kk, vv in v.items() if kk != "lowerings"}
                for k, v in variants.items()
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
