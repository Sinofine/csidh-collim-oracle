"""Compare fresh constructions with the same arithmetic and oracle scaffolding."""

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output", type=Path, default=ROOT / "results/comparison/same-design.json"
)
a = parser.parse_args()
paths = {
    "bonnetain": ROOT / "results/bonnetain/oracle.json",
    "revised": ROOT / "results/reference/oracle.json",
}
records = {name: json.loads(path.read_text()) for name, path in paths.items()}
base = records["bonnetain"]
for name, record in records.items():
    assert record["design"] == name
    assert record["parameters"] == base["parameters"]
    assert record["space"]["explicit_interface"] == 1910
    assert record["resources"]["classical_bits"] == 6
    assert (
        record["resources"]["unresolved_phase_gates"]
        == base["resources"]["unresolved_phase_gates"]
    )
    assert (
        record["layer_resources"]["UnitaryBonnetainLadderStep"]
        == base["layer_resources"]["UnitaryBonnetainLadderStep"]
    )
old = base["resources"]["ccx"]
new = records["revised"]["resources"]["ccx"]

rows = [
    dict(
        design=name,
        qubits=r["resources"]["qubits"],
        clean=r["space"]["clean"],
        ccx=r["resources"]["ccx"],
        ccx_7T_only=r["resources"]["ccx_7T_conversion_only"],
        ccx_change_from_bonnetain=r["resources"]["ccx"] - old,
        ccx_change_percent=100 * (r["resources"]["ccx"] / old - 1),
    )
    for name, r in records.items()
]
n = 511
tm = 28 * n * n
ti = 224 * n * n * math.log2(n)
Told = 2600 * (357998 * tm + 2347 * ti)
Tnew = 2600 * (343982 * tm + 2338 * ti)
result = dict(
    rows=rows,
    paper_model=dict(old_T=Told, new_T=Tnew, reduction_percent=100 * (1 - Tnew / Told)),
)
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(rows, indent=2))
print("Paper-model reduction:", result["paper_model"]["reduction_percent"])
