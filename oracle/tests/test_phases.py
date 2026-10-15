"""Check controlled-phase wiring, native gadgets, and exact seven-T Toffoli.

Arb checks all controlled-phase isometries at high precision. Native Qarton
checks representative emitted gadgets in double precision, with no dummies.
Neither is a simulation of the cryptographic-size group action.
"""

import argparse
import cmath
import json
import math
from functools import lru_cache
from pathlib import Path
import sys
from flint import arb, acb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from arithmetic import install_before_component_imports

install_before_component_imports()
from exact_angles import Angle
from phase_certificates import (
    angle_ball,
    positive_angle,
    inverse_word,
    matrix_word,
    matmul,
    frobenius_upper,
)
from phase_lowering import lower_phases, LoweredPhase, LoweredCPhase
from qarton.circuit import QuantumSimulator, BitVector


def identity(n):
    return [[acb(int(i == j)) for j in range(1 << n)] for i in range(1 << n)]


def embedded(u, target, n):
    m = [[acb(0) for _ in range(1 << n)] for _ in range(1 << n)]
    for j in range(1 << n):
        for bit in (0, 1):
            i = (j & ~(1 << target)) | (bit << target)
            m[i][j] = u[bit][(j >> target) & 1]
    return m


def controlled_x(controls, target, n):
    m = [[acb(0) for _ in range(1 << n)] for _ in range(1 << n)]
    for j in range(1 << n):
        i = j ^ (1 << target) if all((j >> c) & 1 for c in controls) else j
        m[i][j] = acb(1)
    return m


def certify_controlled(catalogue, policy):
    rows = catalogue["variants"][policy]["lowerings"]

    @lru_cache(maxsize=None)
    def get(key):
        a, sign = positive_angle(Angle.from_key(key))
        e = catalogue["rotations"][str(a)]
        word = e["tokens"] if sign > 0 else inverse_word(e["tokens"])
        beta = angle_ball(a) / 2 if e["phase_frame"]["angle_half"] else arb(0)
        beta += arb.pi() * e["phase_frame"]["pi_quarters"] / 4
        return matrix_word(word), beta * sign

    maximum = arb(0)
    checked = 0
    with ctx.workprec(640):
        for row in rows:
            if row["class"] != "CPhase":
                continue
            a = Angle.from_key(row["angle_key"])
            if row["method"] == "and_phase_uncompute":
                u, beta = get(str(a))
                n = 3
                cx = controlled_x([0, 1], 2, 3)
                actual = matmul(cx, matmul(embedded(u, 2, 3), cx))
            else:
                half = a / 2
                u, beta = get(str(half))
                v, _ = get(str(-half))
                n = 2
                cx = controlled_x([0], 1, 2)
                actual = identity(2)
                for gate in [
                    embedded(u, 0, 2),
                    embedded(u, 1, 2),
                    cx,
                    embedded(v, 1, 2),
                    cx,
                ]:
                    actual = matmul(gate, actual)
            factor = acb(0, beta).exp()
            phase = acb(0, angle_ball(a)).exp()
            # Four columns encode arbitrary two-qubit data and clean scratch.
            diff = [
                [
                    factor * actual[i][j] - (phase if j == 3 else acb(1)) * (i == j)
                    for j in range(4)
                ]
                for i in range(1 << n)
            ]
            bound = frobenius_upper(diff)
            assert bound < 2 * row["rotations_each"] * arb(2) ** (
                -catalogue["rotation_error_bits"]
            )
            if bound > maximum:
                maximum = bound
            checked += 1
    return dict(
        controlled_isometries_checked=checked,
        maximum_frobenius_upper_display=str(maximum),
    )


def native_gadgets(catalogue, policy):
    rows = catalogue["variants"][policy]["lowerings"]
    selected = {}
    for row in rows:
        a = Angle.from_key(row["angle_key"])
        category = (row["class"], row["method"], a.kind, a.sn > 0)
        selected.setdefault(category, row)
    results = []
    with lower_phases(ROOT / "results/phase_synthesis.json", policy):
        for row in selected.values():
            is_cp = row["class"] == "CPhase"
            width = 2 if is_cp else 1
            qc = (LoweredCPhase if is_cp else LoweredPhase)(row["angle_key"])
            assert tuple(list(qc.input_signature)[0].positions) == tuple(range(width))
            amps = [(j + 1) * cmath.exp((j + 1) * 0.271j) for j in range(1 << width)]
            for _ in range(4):
                norm = math.sqrt(sum(abs(a) ** 2 for a in amps))
                amps = [a / norm for a in amps]
            initial = {
                BitVector.from_int(j, length=qc.nbr_bits()): a
                for j, a in enumerate(amps)
            }
            actual = (
                QuantumSimulator(qc, no_dummies=True)
                .simulate_internal(initial)
                .quantum_state()
            )
            theta = float(Angle.from_key(row["angle_key"]))
            wanted = {
                BitVector.from_int(j, length=qc.nbr_bits()): a
                * (cmath.exp(1j * theta) if j == (1 << width) - 1 else 1)
                for j, a in enumerate(amps)
            }
            overlap = sum(a.conjugate() * actual.get(k, 0) for k, a in wanted.items())
            phase = overlap / abs(overlap)
            error = sum(
                abs(actual.get(k, 0) - phase * wanted.get(k, 0)) ** 2
                for k in set(actual) | set(wanted)
            )
            garbage = sum(
                abs(a) ** 2
                for k, a in actual.items()
                if not k.is_zero(qc.output_signature.all_zeroes())
            )
            assert error < 1e-20 and garbage < 1e-20, (row, error, garbage)
            results.append(
                dict(
                    gate_class=row["class"],
                    angle_key=row["angle_key"],
                    method=row["method"],
                    squared_state_error_up_to_global_phase=error,
                    quantum_workspace_nonzero_probability=garbage,
                )
            )
    return results


def exact_toffoli():
    import sympy as s

    root = s.sqrt(2)
    h = s.Matrix([[1, 1], [1, -1]]) / root
    t = s.diag(1, (1 + s.I) / root)
    td = t.conjugate().T

    def one(u, q):
        m = s.zeros(8)
        for j in range(8):
            for bit in (0, 1):
                m[(j & ~(1 << q)) | (bit << q), j] = u[bit, (j >> q) & 1]
        return m

    def cx(c, q):
        m = s.zeros(8)
        for j in range(8):
            m[j ^ (1 << q) if (j >> c) & 1 else j, j] = 1
        return m

    schedule = [
        ("h", 2),
        ("cx", 1, 2),
        ("tdg", 2),
        ("cx", 0, 2),
        ("t", 2),
        ("cx", 1, 2),
        ("tdg", 2),
        ("cx", 0, 2),
        ("t", 1),
        ("t", 2),
        ("h", 2),
        ("cx", 0, 1),
        ("t", 0),
        ("tdg", 1),
        ("cx", 0, 1),
    ]
    u = s.eye(8)
    for g, *args in schedule:
        gate = cx(*args) if g == "cx" else one({"h": h, "t": t, "tdg": td}[g], args[0])
        u = gate * u
    wanted = s.zeros(8)
    for j in range(8):
        wanted[j ^ 4 if j & 3 == 3 else j, j] = 1
    assert all(s.simplify(v) == 0 for v in u - wanted)
    return dict(
        exact_symbolic_matrix_verified=True, T_count=7, chronological_schedule=schedule
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--policy",
        choices=["ancilla_free", "one_clean_ancilla"],
        default="one_clean_ancilla",
    )
    args = parser.parse_args()
    catalogue = json.loads((ROOT / "results/phase_synthesis.json").read_text())
    result = dict(
        policy=args.policy,
        arb=certify_controlled(catalogue, args.policy),
        native_gadgets=native_gadgets(catalogue, args.policy),
        toffoli=exact_toffoli(),
    )
    output = ROOT / "results" / f"phase_validation_{args.policy}.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                output=str(output),
                arb=result["arb"],
                native_gadget_tests=len(result["native_gadgets"]),
                exact_7T_toffoli=True,
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
