"""Small reversible control utilities shared by the reference circuits."""

from qarton.circuit import BoolType, Circuit, memoize


def cswap(circuit, control, a, b):
    for x, y in zip(a, b):
        circuit.cx(x, y)
        circuit.ccx(control, y, x)
        circuit.cx(x, y)


@memoize
class SetFlag(Circuit):
    def __init__(self):
        super().__init__()
        f = self.add_input_output(BoolType())
        self.x(f[0])
