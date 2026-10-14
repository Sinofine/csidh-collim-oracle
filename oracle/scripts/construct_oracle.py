"""Construct a CSIDH-512 oracle and count its resources."""

import argparse
import json
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
CODE = PROJECT / "code"
sys.path.insert(0, str(CODE))


def construct(output, design="revised", gate_audit=False):
    # Resolve the arithmetic profile before any parent components are imported.
    from arithmetic import install_before_component_imports

    install_before_component_imports()
    use_binary_character = design == "bonnetain"
    use_montgomery_isogeny = design == "bonnetain"
    if use_binary_character or use_montgomery_isogeny:
        import bonnetain

        bonnetain.CHARACTER_CHECKPOINTS = 60
    if use_binary_character:
        import character

        character.append_chain = bonnetain.append_binary_character
        character.UnitaryPositiveCharacterFlag = (
            bonnetain.BonnetainPositiveCharacterFlag
        )
    if use_montgomery_isogeny:
        import conversion

        conversion.UnitaryValueKernelAdapter = bonnetain.BonnetainKernelAdapter
    from group_action import ValueGroupAction
    from parameters import P512, PRIMES
    from points import UNITARY
    from qarton.circuit import Circuit
    from resources import report

    parents = {}
    audit_registry = {}
    audit_noninteger_weights = []
    counter_method = None
    if gate_audit:
        from qarton.circuit.sub_circuits_counter import SubCircuitsCounter

        counter_method = SubCircuitsCounter.add_circuit

        def audited_add(counter, gate, quantity=1):
            audit_registry[gate.name] = gate
            if float(quantity) != int(quantity):
                audit_noninteger_weights.append(
                    {"gate": gate.class_name, "weight": quantity}
                )
            return counter_method(counter, gate, quantity)

        SubCircuitsCounter.add_circuit = audited_add
    original_append = Circuit.append

    def append(parent, gate, *regs, **kwargs):
        result = original_append(parent, gate, *regs, **kwargs)
        if type(parent).__module__ in {
            "group_action",
            "preparation",
            "character",
            "ladder",
            "points",
            "conversion",
            "isogeny",
            "utils",
            "bonnetain",
        }:
            parents[type(parent).__name__] = parent
        return result

    start = time.monotonic()
    Circuit.append = append
    try:
        circuit = ValueGroupAction(P512, PRIMES, 1300, 5, 15, backends=UNITARY)
    finally:
        Circuit.append = original_append
        if gate_audit:
            SubCircuitsCounter.add_circuit = counter_method
    resources = report(circuit)
    resources.pop("native_t", None)
    audit_details = {}
    if gate_audit:
        from qarton.circuit import ResourceReporter

        for name, parent in sorted(parents.items()):
            direct = []
            total = 0
            missing = []
            for (
                gate_name,
                quantity,
            ) in parent.get_sub_circuits_counter()._direct_count_by_name.items():
                gate = audit_registry.get(gate_name)
                if gate is None:
                    missing.append({"name": str(gate_name), "quantity": quantity})
                    continue
                cost = (
                    gate.resources.ccx_count()
                    if isinstance(gate, Circuit)
                    else int(gate.class_name == "ccx")
                )
                total += int(quantity) * cost
                direct.append(
                    {"class": gate.class_name, "quantity": quantity, "ccx_each": cost}
                )
            audit_details[name] = {
                "direct_children": direct,
                "known_direct_ccx_sum": total,
                "reported_ccx": parent.resources.ccx_count(),
                "unresolved_direct_names": missing,
                "class_counts": ResourceReporter(parent).sub_circuits_by_class(),
            }
    result = {
        "design": design,
        "parameters": {"n": 511, "s": 15, "M": 1300, "window": 5},
        "resources": resources,
        "space": {
            "explicit_interface": 1910,
            "clean": circuit.nbr_input_anc(),
            "dirty": circuit.nbr_dirty_anc(),
            "total": circuit.nbr_qubits(),
        },
        "layer_space": {
            name: {
                "total": c.nbr_qubits(),
                "clean": c.nbr_input_anc(),
                "dirty": c.nbr_dirty_anc(),
            }
            for name, c in sorted(parents.items())
        },
        "layer_resources": {name: report(c) for name, c in sorted(parents.items())},
        "elapsed_seconds": time.monotonic() - start,
        "gate_audit": {
            "noninteger_append_weights": audit_noninteger_weights,
            "layers": audit_details,
        }
        if gate_audit
        else None,
    }
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {"output": str(output), "resources": resources, "space": result["space"]},
            indent=2,
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--design", choices=["revised", "bonnetain"], default="revised")
    parser.add_argument(
        "--gate-audit",
        action="store_true",
        help="Record direct child costs and counter weights",
    )
    args = parser.parse_args()
    output = args.output or PROJECT / (
        "results/reference/oracle.json"
        if args.design == "revised"
        else f"results/{args.design}/oracle.json"
    )
    construct(output, args.design, args.gate_audit)


if __name__ == "__main__":
    main()
