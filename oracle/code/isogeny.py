"""Reversible product chains and local isogeny evaluation."""

from points import UNITARY, UnitaryInitialPointStep, build_unitary_point_step
from qarton.arithmetic import qc_lt_uint
from qarton.circuit import (
    BoolType,
    Circuit,
    UInt,
    UIntType,
    dummify,
    memoize,
    qc_bool_reg,
)
from qarton.modular_arithmetic import (
    ModInt,
    ModIntType,
    ModularAdder,
    ModularInverseAdd,
    ModularProductAdd,
    ModularSquareAdd,
)
from qarton.modular_arithmetic.modular_double import ModularDouble
from qarton.modular_arithmetic.modular_inverse import ControlledModularIPMult
from schedules import bennett


def toggle_selected_one(c, bit, x, z):
    c.x(z[0])
    c.cx(bit[0], z[0])
    for a, b in zip(x, z):
        c.ccx(bit[0], a, b)


@memoize
class UnitaryPowerStep(Circuit):
    def __init__(self, p):
        super().__init__()
        t = ModIntType(p)
        x = self.add_input_output(t)
        bit = self.add_input_output(BoolType())
        acc = self.add_input_output(t)
        out = self.add_input_output(t)
        square = self.get_anc(t)
        factor = self.get_anc(t)
        sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
        mul = ModularProductAdd(p, zero_input=True, backends=UNITARY)
        toggle_selected_one(self, bit, x, factor)
        self.append(sq, acc, square)
        self.append(mul, square, factor, out)
        self.append(sq.inverse(), acc, square)
        toggle_selected_one(self, bit, x, factor)
        self.release_anc(square, factor)


@dummify
@memoize
class UnitaryVisitingPower(Circuit):
    def __init__(self, p, bits, checkpoints=15):
        super().__init__()
        assert bits >= 1
        self.p = p
        self.bits = bits
        t = ModIntType(p)
        x = self.add_input_output(t)
        exponent = self.add_input_output(UIntType(bits))
        out = self.add_input_output(t)
        step = UnitaryPowerStep(p)

        def alloc():
            return self.get_anc(t)

        def release(q):
            self.release_anc(q)

        def apply(i, src, dst, inverse):
            initial = None
            if i == 0:
                assert src is None
                initial = alloc()
                toggle_selected_one(self, qc_bool_reg(exponent[bits - 1]), x, initial)
                src = initial
            if i == bits - 1:
                assert dst is None and not inverse
                self.cx_reg(src, out)
            else:
                self.append(
                    step.inverse() if inverse else step,
                    x,
                    qc_bool_reg(exponent[bits - 2 - i]),
                    src,
                    dst,
                )
            if initial is not None:
                toggle_selected_one(self, qc_bool_reg(exponent[bits - 1]), x, initial)
                release(initial)

        bennett(bits, checkpoints).bridge(
            0, bits, checkpoints, None, None, alloc, release, apply
        )

    def dummy_classical_function(self, args):
        x, e, out = args
        v = out.v ^ pow(x.v, int(e), self.p)
        if v >= self.p:
            raise ValueError("XOR target outside canonical domain")
        return x, e, ModInt(v, self.p)

    def dummy_classical_function_inverse(self, args):
        return self.dummy_classical_function(args)


@memoize
class UnitaryValueTail(Circuit):
    def __init__(self, p):
        super().__init__()
        t = ModIntType(p)
        u, N, D, out = [self.add_input_output(t) for _ in range(4)]
        inv = ModularInverseAdd(p, backends=UNITARY)
        mul = ModularProductAdd(p, zero_input=True, backends=UNITARY)
        sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
        tape = []

        def value(g, *args):
            z = self.get_anc(t)
            self.append(g, *args, z)
            tape.append((g, args, z))
            return z

        di = value(inv, D)
        h = value(mul, N, di)
        h2 = value(sq, h)
        h4 = value(sq, h2)
        h8 = value(sq, h4)
        self.append(mul, u, h8, out)
        for g, args, z in reversed(tape):
            self.append(g.inverse(), *args, z)
            self.release_anc(z)


@memoize
class UnitaryWProductUpdate(Circuit):
    def __init__(self, p, ellbits, index):
        super().__init__()
        t = ModIntType(p)
        W = self.add_input_output(t)
        Z = self.add_input_output(t)
        ell = self.add_input_output(UIntType(ellbits))
        N = self.add_input_output(t)
        D = self.add_input_output(t)
        add = ModularAdder(p, backends=UNITARY)
        dbl = ModularDouble(p, backends=UNITARY)
        if index == 0:
            self.cx_reg(W, N)
            self.append(add, Z, N)
            self.cx_reg(Z, D)
            self.append(dbl, D)
            return
        f = self.get_anc(BoolType())
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        mult = ControlledModularIPMult(p, backends=UNITARY)
        self.append(add, Z, W)
        self.append(mult, f, W, N)
        self.append(add.inverse(), Z, W)
        self.append(dbl, Z)
        self.append(mult, f, Z, D)
        self.append(dbl.inverse(), Z)
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        self.release_anc(f)


def append_fused_walk(
    c,
    p,
    primes,
    s,
    coef,
    W,
    Z,
    ell,
    point_class,
    root,
    product_class=UnitaryWProductUpdate,
):
    t = ModIntType(p)
    bits = max(primes).bit_length()
    m = (max(primes) - 1) // 2
    dp = bennett(m, s)
    events = []
    dp.bridge(
        0,
        m,
        s,
        None,
        None,
        lambda: None,
        lambda _: None,
        lambda i, src, dst, inv: events.append((i, inv)),
    )
    first = {}
    last = {}
    for pos, (i, inv) in enumerate(events):
        first.setdefault(i, pos)
        last[i] = pos
    root_pos = first[m - 1]
    assert last[m - 1] == root_pos and last[0] == len(events) - 1
    assert all(first[i] < root_pos < last[i] for i in range(m - 1))
    N = c.get_anc(t)
    D = c.get_anc(t)
    seen = set()
    erased = set()
    cursor = 0

    def alloc():
        return tuple(c.get_anc(t) for _ in range(4))

    def release(rs):
        c.release_anc(*rs)

    def apply(i, src, dst, inverse):
        nonlocal cursor
        pos = cursor
        cursor += 1
        assert events[pos] == (i, inverse)
        initial = None
        if i == 0:
            assert src is None
            initial = alloc()
            c.x(initial[1][0])
            c.cx_reg(W, initial[2])
            c.cx_reg(Z, initial[3])
            src = initial
        product = product_class(p, bits, i, backends=UNITARY)
        if pos == first[i]:
            assert not inverse and i not in seen
            c.append(product, src[2], src[3], ell, N, D)
            seen.add(i)
        if i == m - 1:
            assert dst is None and not inverse and len(seen) == m
            root(N, D)
        if pos == last[i]:
            assert pos >= root_pos and i not in erased
            if i == 0:
                assert len(erased) == m - 1
            c.append(product.inverse(), src[2], src[3], ell, N, D)
            erased.add(i)
        if i != m - 1:
            gate = point_class(p, bits, i)
            c.append(gate.inverse() if inverse else gate, coef, W, Z, ell, *src, *dst)
        if initial is not None:
            c.x(initial[1][0])
            c.cx_reg(W, initial[2])
            c.cx_reg(Z, initial[3])
            release(initial)

    dp.bridge(0, m, s, None, None, alloc, release, apply)
    assert cursor == len(events) and len(erased) == m
    c.release_anc(N, D)


def ledger_point_step(p, bits, index):
    return (
        UnitaryInitialPointStep(p, bits, 0)
        if index == 0
        else build_unitary_point_step(p, bits, index)
    )


@memoize
class UnitaryLocalValueImage(Circuit):
    """d,QW,QZ,ell,0 -> same inputs,dprime; inverse erases matching result.
    Product traversal embeds the image tail at its root. Exponent result is
    explicitly evaluated before the walk and erased after the walk.
    """

    def __init__(self, p, primes, s=15):
        super().__init__()
        t = ModIntType(p)
        bits = max(primes).bit_length()
        d = self.add_input_output(t)
        W = self.add_input_output(t)
        Z = self.add_input_output(t)
        ell = self.add_input_output(UIntType(bits))
        out = self.add_input_output(t)
        u = self.get_anc(t)
        power = UnitaryVisitingPower(p, bits, s)
        self.append(power, d, ell, u)

        def root(N, D):
            self.append(UnitaryValueTail(p), u, N, D, out)

        append_fused_walk(self, p, primes, s, d, W, Z, ell, ledger_point_step, root)
        self.append(power.inverse(), d, ell, u)
        self.release_anc(u)
