"""Reversible Montgomery ladder and cofactor tests."""

from expressions import Expr
from points import UNITARY
from points import UnitaryInitialExpression as ExpressionCircuit
from qarton.binary_operations import qc_eq
from qarton.circuit import BoolType, Circuit, UIntType, dummify, memoize, qc_bool_reg
from qarton.modular_arithmetic import (
    ModInt,
    ModIntType,
    ModularAdder,
    ModularProductAdd,
    ModularSquareAdd,
)
from qarton.modular_arithmetic.modular_neg import ModularNeg
from reference import ladder_reference, ladder_reference_inverse
from schedules import bennett
from utils import cswap


@memoize
class UnitaryBonnetainLadderStep(Circuit):
    def __init__(self, p):
        super().__init__()
        t = ModIntType(p)
        a24 = self.add_input_output(t)
        x = self.add_input_output(t)
        bit = self.add_input_output(BoolType())
        X, Z, U, V = [self.add_input_output(t) for _ in range(4)]
        X4, Z4, X5, Z5 = [self.add_input_output(t) for _ in range(4)]
        a, b = self.get_anc(t), self.get_anc(t)
        add = ModularAdder(p, backends=UNITARY)
        neg = ModularNeg(p, backends=UNITARY)
        mul = ModularProductAdd(p, zero_input=True, backends=UNITARY)
        sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)

        def plus(src, dst):
            self.append(add, src, dst)

        def minus(src, dst):
            self.append(add.inverse(), src, dst)

        def reflect(src, dst):
            self.append(neg, dst)
            plus(src, dst)
            plus(src, dst)

        cswap(self, bit[0], X, U)
        cswap(self, bit[0], Z, V)
        # Y=X-Z, C=U+V; preserve source coordinates reversibly.
        self.append(neg, Z)
        plus(X, Z)
        plus(U, V)
        self.append(mul, Z, V, a)  # CY
        reflect(X, Z)
        reflect(U, V)  # T, D
        self.append(mul, Z, V, b)  # DT
        plus(a, b)
        self.append(sq, b, X5)
        minus(a, b)
        minus(a, b)
        self.append(sq, b, X4)
        self.append(mul, x, X4, Z5)
        self.append(sq.inverse(), b, X4)
        plus(a, b)
        self.append(mul.inverse(), Z, V, b)
        reflect(X, Z)
        reflect(U, V)
        self.append(mul.inverse(), Z, V, a)
        minus(U, V)
        # Doubling: a=Y^2, b=T^2 then E=T^2-Y^2.
        self.append(sq, Z, a)
        reflect(X, Z)
        self.append(sq, Z, b)
        minus(a, b)
        self.append(mul, a24, b, X4)
        plus(a, X4)
        self.append(mul, b, X4, Z4)
        minus(a, X4)
        self.append(mul.inverse(), a24, b, X4)
        plus(a, b)
        self.append(mul, a, b, X4)
        self.append(sq.inverse(), Z, b)
        reflect(X, Z)
        self.append(sq.inverse(), Z, a)
        self.append(neg, Z)
        plus(X, Z)
        cswap(self, bit[0], X, U)
        cswap(self, bit[0], Z, V)
        cswap(self, bit[0], X4, X5)
        cswap(self, bit[0], Z4, Z5)
        self.release_anc(a, b)


def append_ladder(c, p, bits, s, A, x, scalar, root):
    t = ModIntType(p)
    s = max(s, bits.bit_length())
    e = Expr(p, 1)
    a24 = e.scale(e.add(0, e.const(2)), pow(4, -1, p))
    init = ExpressionCircuit(e, [a24])
    ar = c.get_anc(t)
    c.append(init, A, ar)

    def alloc():
        return tuple(c.get_anc(t) for _ in range(4))

    def release(rs):
        c.release_anc(*rs)

    step = UnitaryBonnetainLadderStep(p, backends=UNITARY)

    def apply(i, src, dst, inverse):
        initial = None
        if i == 0:
            assert src is None
            initial = alloc()
            c.x(initial[0][0])
            c.cx_reg(x, initial[2])
            c.x(initial[3][0])
            src = initial
        if i == bits:
            assert dst is None and not inverse
            root(src)
        else:
            c.append(
                step.inverse() if inverse else step,
                ar,
                x,
                qc_bool_reg(scalar[bits - 1 - i]),
                *src,
                *dst,
            )
        if initial is not None:
            c.x(initial[0][0])
            c.cx_reg(x, initial[2])
            c.x(initial[3][0])
            release(initial)

    bennett(bits + 1, s).bridge(0, bits + 1, s, None, None, alloc, release, apply)
    c.append(init.inverse(), A, ar)
    c.release_anc(ar)


@dummify
@memoize
class UnitaryValueLadder(Circuit):
    def __init__(self, p, scalar_bits, checkpoints=15):
        super().__init__()
        self.p = p
        self.scalar_bits = scalar_bits
        t = ModIntType(p)
        A = self.add_input_output(t)
        x = self.add_input_output(t)
        scalar = self.add_input_output(UIntType(scalar_bits))
        X = self.add_input_output(t)
        Z = self.add_input_output(t)

        def root(src):
            self.cx_reg(src[0], X)
            self.cx_reg(src[1], Z)

        append_ladder(self, p, scalar_bits, checkpoints, A, x, scalar, root)

    dummy_classical_function = ladder_reference
    dummy_classical_function_inverse = ladder_reference_inverse


@memoize
class UnitaryCofactorFlag(Circuit):
    def __init__(self, p, checkpoints=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        x = self.add_input_output(t)
        scalar = self.add_input_output(UIntType(p.bit_length()))
        flag = self.add_input_output(BoolType())

        def root(src):
            qc_eq(ModInt(0, p), src[1], flag)
            self.x(flag[0])

        append_ladder(self, p, p.bit_length(), checkpoints, A, x, scalar, root)
