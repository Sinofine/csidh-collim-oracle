"""Build the actual phase-lowered Qarton DAG and reconcile its gate counters."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from arithmetic import install_before_component_imports

install_before_component_imports()
from phase_lowering import lower_phases
from construct_oracle import construct

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--design", choices=["revised", "bonnetain"], default="revised")
parser.add_argument(
    "--policy",
    choices=["ancilla_free", "one_clean_ancilla"],
    default="one_clean_ancilla",
)
parser.add_argument(
    "--catalogue", type=Path, default=ROOT / "results/phase_synthesis.json"
)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
output = (
    args.output or ROOT / "results" / f"clifford_t_{args.design}_{args.policy}.json"
)
with lower_phases(args.catalogue, args.policy) as catalogue:
    construct(output, args.design, True)
result = json.loads(output.read_text())
counts = result["resources"]["primitives"]
t = 7 * counts["ccx"] + counts.get("t", 0) + counts.get("tdg", 0)
expected = catalogue["variants"][args.policy][
    "revised_t" if args.design == "revised" else "baseline_t"
]
assert t == expected, (t, expected)
assert not any(result["resources"]["unresolved_phase_gates"].values())
result["compiled_T_count"] = t
result["phase_compilation"] = {
    "policy": args.policy,
    "budget_bits": catalogue["budget_bits"],
    "rotation_error_bits": catalogue["rotation_error_bits"],
}
output.write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(dict(compiled_T_count=t, total_qubits=result["space"]["total"])),
    flush=True,
)
