"""Fixed-exponent Legendre symbol with a width-five window."""

from points import UNITARY
from qarton.circuit import BoolType, Circuit, memoize
from qarton.modular_arithmetic import (
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

    # Use the emitted bounded schedule: each recursive prefix may borrow its
    # parent's still-zero destination. The former recursive allocator did not
    # enforce this physical checkpoint bound.
    pool = [c.get_anc(t) for _ in range(checkpoints + 1)]
    dp = bennett(len(chain) + 1, checkpoints)

    def emit(start, length, s, src, dst, free):
        if length == 1:
            return [(start, src, dst, False)]
        k = dp.k[length][s]
        if not k or not free:
            raise ValueError("insufficient checkpoints")
        mid, remaining = free[0], free[1:]
        prefix = emit(start, k, s, src, mid, remaining + [dst])
        suffix = emit(start + k, length - k, s - 1, mid, dst, remaining)
        erase = emit(start, k, s - 1, src, mid, remaining)
        return (
            prefix + suffix + [(i, a, b, not inv) for i, a, b, inv in reversed(erase)]
        )

    events = emit(0, len(chain) + 1, checkpoints, odd[a0], pool[-1], pool[:-1])
    assert len(events) == dp.b[len(chain) + 1][checkpoints]
    for i, src, dst, inverse in events:
        apply(i, src, None if i == len(chain) else dst, inverse)
    c.release_anc(*pool)
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
            # Euler's criterion gives exactly 0, 1, or p-1. Since p is odd,
            # only +1 has low bit 1. No n-bit equality circuit is necessary.
            self.ccx(enable[0], src[0], out[0])

        append_chain(self, x, p, width, checkpoints, root)


@memoize
class UnitaryNegativeCharacterPhase(Circuit):
    """Clean phase (-1)^[chi(x)=-1], including chi(0)=0.

    The high bit distinguishes p-1 from both 0 and 1 on the reachable output
    domain. The terminal marker is Clifford; no TM is charged to this Z.
    """

    def __init__(self, p, width=5, checkpoints=44):
        super().__init__()
        x = self.add_input_output(ModIntType(p))
        append_chain(self, x, p, width, checkpoints, lambda y: self.z(y[-1]))
