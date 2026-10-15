"""Check synthesized phase costs against the constructed Clifford+T circuits."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--policy",
    choices=["ancilla_free", "one_clean_ancilla"],
    default="one_clean_ancilla",
)
args = parser.parse_args()
catalogue = json.loads((ROOT / "results/phase_synthesis.json").read_text())
variant = catalogue["variants"][args.policy]
rows = []
phase_cost = sum(row["count"] * row["t_each"] for row in variant["lowerings"])
assert phase_cost == variant["phase_t"]
assert 2 * variant["rotation_occurrences"] <= 2 ** (
    catalogue["rotation_error_bits"] - catalogue["budget_bits"]
)
for design, folder, model in [
    ("bonnetain", "bonnetain", "baseline_t"),
    ("revised", "reference", "revised_t"),
]:
    original = json.loads((ROOT / f"results/{folder}/oracle.json").read_text())
    result = json.loads(
        (ROOT / f"results/clifford_t_{design}_{args.policy}.json").read_text()
    )
    counts = result["resources"]["primitives"]
    total = 7 * counts["ccx"] + counts.get("t", 0) + counts.get("tdg", 0)
    assert total == result["compiled_T_count"] == variant[model]
    assert total == original["resources"]["ccx_7T_conversion_only"] + phase_cost
    assert counts["ccx"] == original["resources"]["ccx"] + variant["extra_ccx"]
    assert not any(result["resources"]["unresolved_phase_gates"].values())
    assert not result["gate_audit"]["noninteger_append_weights"]
    assert result["space"]["total"] == sum(
        result["space"][k] for k in ["explicit_interface", "clean", "dirty"]
    )
    for layer in result["gate_audit"]["layers"].values():
        assert not layer["unresolved_direct_names"]
        assert layer["known_direct_ccx_sum"] == layer["reported_ccx"]
    rows.append({"design": design, "qubits": result["space"]["total"], "T": total})
print(
    json.dumps(
        {
            "policy": args.policy,
            "rows": rows,
            "phase_T": phase_cost,
            "channel_error_bound": variant["channel_error_bound"],
        },
        indent=2,
    )
)
