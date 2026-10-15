"""Unitary arithmetic backend selection."""

from qarton_compat import install_bitvector_equality_fix

install_bitvector_equality_fix()

from qarton.arithmetic.comparator_adder import CDKMComparatorAdder
from qarton.arithmetic.controlled_adder import ControlledCDKMAdder
from qarton.binary_operations.and_ccx import AndGate, AndGateDecomposed, FakeAndGate
from qarton.circuit.backends import MemoryConstraint, PriorityBackends
from qarton.circuit.basic_gates import MSR

# Count AND operations using their CCX implementation.
AndGate.replace_by_ccx = True
STABLE_UNITARY = MemoryConstraint(
    None,
    priority=PriorityBackends(FakeAndGate, CDKMComparatorAdder, ControlledCDKMAdder),
    exclude=PriorityBackends(MSR, AndGate, AndGateDecomposed),
)


def install_before_component_imports():
    import sys

    forbidden = ("group_action", "preparation", "ladder", "isogeny", "conversion")
    if any(n in sys.modules for n in forbidden):
        raise RuntimeError(
            "Install stable profile in a fresh process before importing components"
        )
    import points as unitary_point_step

    unitary_point_step.UNITARY = STABLE_UNITARY
