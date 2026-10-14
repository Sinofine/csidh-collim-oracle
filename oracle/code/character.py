"""Fixed-exponent Legendre symbol with a width-five window."""

from points import UNITARY
from qarton.binary_operations import qc_eq
from qarton.circuit import BoolType, Circuit, memoize
from qarton.modular_arithmetic import (
    ModInt,
    ModIntType,
    ModularProductAdd,
    ModularSquareAdd,
)
from schedules import bennett, window


def append_chain(c, x, p, width, checkpoints, root):
    t = ModIntType(p)
    a0, chain, _ = window((p - 1) // 2, width)
    sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
    mul = ModularProductAdd(p, zero_input=True, backends=UNITARY)
    x2 = c.get_anc(t)
    c.append(sq, x, x2)
    tape = [(sq, (x, x2))]
    odd = {1: x}
    for a in range(3, 2**width, 2):
        z = c.get_anc(t)
        c.append(mul, odd[a - 2], x2, z)
        tape.append((mul, (odd[a - 2], x2, z)))
        odd[a] = z

    def apply(i, src, dst, inverse):
        if i == len(chain):
            assert dst is None and not inverse
            root(src)
        else:
            op, a = chain[i]
            g = sq if op == "square" else mul
            args = (src, dst) if op == "square" else (src, odd[a], dst)
            c.append(g.inverse() if inverse else g, *args)

    bennett(len(chain) + 1, checkpoints).bridge(
        0,
        len(chain) + 1,
        checkpoints,
        odd[a0],
        None,
        lambda: c.get_anc(t),
        c.release_anc,
        apply,
    )
    for g, args in reversed(tape):
        c.append(g.inverse(), *args)
        c.release_anc(args[-1])


@memoize
class UnitaryPositiveCharacterFlag(Circuit):
    """x,enable,out -> x,enable,out XOR (enable AND [chi(x)=1])."""

    def __init__(self, p, width=5, checkpoints=44):
        super().__init__()
        x = self.add_input_output(ModIntType(p))
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(BoolType())

        def root(src):
            f = self.get_anc(BoolType())
            qc_eq(ModInt(1, p), src, f)
            self.ccx(enable[0], f[0], out[0])
            qc_eq(ModInt(1, p), src, f)
            self.release_anc(f)

        append_chain(self, x, p, width, checkpoints, root)
