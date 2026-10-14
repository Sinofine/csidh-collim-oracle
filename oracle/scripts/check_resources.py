"""Check recorded circuit resources."""

import argparse
import ast
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--result", type=Path)
parser.add_argument(
    "--compare-to", type=Path, help="Optional independent construction snapshot"
)
args = parser.parse_args()
for file in (ROOT / "code").rglob("*.py"):
    ast.parse(file.read_text(), filename=str(file))
if args.result:
    actual = json.loads(args.result.read_text())
    assert actual["space"]["total"] == actual["resources"]["qubits"]
    assert actual["space"]["total"] == sum(
        actual["space"][key] for key in ["explicit_interface", "clean", "dirty"]
    )
    assert (
        actual["resources"]["ccx_7T_conversion_only"] == 7 * actual["resources"]["ccx"]
    )
    baseline_path = args.compare_to
    if baseline_path is not None:
        baseline = json.loads(baseline_path.read_text())
        assert actual["space"] == baseline["space"], (
            actual["space"],
            baseline["space"],
        )
        assert actual["resources"] == baseline["resources"], (
            "Resource statistics changed"
        )
        for name, saved in baseline["layer_space"].items():
            assert all(
                actual["layer_space"][name][key] == saved[key]
                for key in ["total", "clean", "dirty"]
            ), name

print(
    "Resource accounting checked"
    if args.result
    else "Python sources checked"
)
