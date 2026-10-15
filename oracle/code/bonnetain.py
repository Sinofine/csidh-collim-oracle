"""Bonnetain--Schrottenloher character and isogeny circuits.

Reference: ePrint 2018/537, Sections 5.1--5.3 and Appendix B.3."""

from expressions import Expr
from isogeny import UnitaryPowerStep, UnitaryVisitingPower
from parameters import P512
from points import UNITARY, MaskedProductIntoZero, UnitaryInitialExpression
from qarton.arithmetic import qc_lt_uint
from qarton.circuit import BoolType, Circuit, UInt, UIntType, memoize
from qarton.modular_arithmetic import (
    ModIntType,
    ModularAdder,
    ModularInverseAdd,
    ModularProductAdd,
    ModularSquareAdd,
)
from qarton.modular_arithmetic.modular_double import ModularDouble
from qarton.modular_arithmetic.modular_inverse import (
    ControlledModularIPMult,
    ModularIPMult,
)
from qarton.modular_arithmetic.modular_neg import ModularNeg
from schedules import bennett

CHARACTER_CHECKPOINTS = 60


def append_binary_character(c, x, p, width, checkpoints, root):
    """Binary square/multiply using the same three-operation bit-step primitive.

    CSIDH-512 uses the published 512-map reference length, including a terminal
    marking callback. Two leading zero bits pad the 510-bit fixed exponent.
    The passed window/checkpoint arguments belong to the shared interface;
    the baseline uses CHARACTER_CHECKPOINTS=4s, selected explicitly by the runner.
    """
    bits = 512 if p == P512 else p.bit_length()
    exponent = (p - 1) // 2
    t = ModIntType(p)
    step = UnitaryPowerStep(p, backends=UNITARY)
    bit = c.get_anc(BoolType())

    def apply(i, src, dst, inverse):
        initial = None
        if i == 0:
            initial = c.get_anc(t)
            if exponent >> (bits - 1) & 1:
                c.cx_reg(x, initial)
            else:
                c.x(initial[0])
            src = initial
        if i == bits - 1:
            root(src)
        else:
            value = exponent >> (bits - 2 - i) & 1
            if value:
                c.x(bit[0])
            c.append(step.inverse() if inverse else step, x, bit, src, dst)
            if value:
                c.x(bit[0])
        if initial is not None:
            if exponent >> (bits - 1) & 1:
                c.cx_reg(x, initial)
            else:
                c.x(initial[0])
            c.release_anc(initial)

    bennett(bits, CHARACTER_CHECKPOINTS).bridge(
        0,
        bits,
        CHARACTER_CHECKPOINTS,
        None,
        None,
        lambda: c.get_anc(t),
        lambda q: c.release_anc(q),
        apply,
    )
    c.release_anc(bit)


@memoize
class BonnetainPositiveCharacterFlag(Circuit):
    """Same x, enable, flag interface as the revised character routine."""

    def __init__(self, p, width=5, checkpoints=44):
        super().__init__()
        x = self.add_input_output(ModIntType(p))
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(BoolType())

        def root(value):
            # Apply the same reachable-domain bit test to both comparators.
            self.ccx(enable[0], value[0], out[0])

        append_binary_character(self, x, p, width, checkpoints, root)


@memoize
class BonnetainPointStep(Circuit):
    """Montgomery iQ -> (i+1)Q in the same four-field checkpoint interface."""

    def __init__(self, p, bits, index):
        super().__init__()
        t = ModIntType(p)
        a24, QX, QZ = [self.add_input_output(t) for _ in range(3)]
        ell = self.add_input_output(UIntType(bits))
        PX, PZ, X, Z = [self.add_input_output(t) for _ in range(4)]
        OPX, OPZ, OX, OZ = [self.add_input_output(t) for _ in range(4)]
        if index == 0:
            # The first transition doubles Q, as in the revised initial step.
            e = Expr(p, 3)
            plus = e.sq(e.add(1, 2))
            minus = e.sq(e.sub(1, 2))
            delta = e.sub(plus, minus)
            nx = e.mul(plus, minus)
            nz = e.mul(delta, e.add(minus, e.mul(0, delta)))
            self.cx_reg(X, OPX)
            self.cx_reg(Z, OPZ)
            self.append(UnitaryInitialExpression(e, [nx, nz]), a24, X, Z, OX, OZ)
            return
        f = self.get_anc(BoolType())
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        a, b = self.get_anc(t), self.get_anc(t)
        temporary = OPX
        add = ModularAdder(p, backends=UNITARY)
        neg = ModularNeg(p, backends=UNITARY)
        mul = ModularProductAdd(p, zero_input=True, backends=UNITARY)
        sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
        cmul = MaskedProductIntoZero(p, borrowed=True)

        def plus(src, dst):
            self.append(add, src, dst)

        def minus(src, dst):
            self.append(add.inverse(), src, dst)

        def reflect(src, dst):
            self.append(neg, dst)
            plus(src, dst)
            plus(src, dst)

        # CY=(QX+QZ)(X-Z), DT=(QX-QZ)(X+Z).
        self.append(neg, Z)
        plus(X, Z)
        plus(QX, QZ)
        self.append(mul, Z, QZ, a)
        reflect(X, Z)
        reflect(QX, QZ)
        self.append(mul, Z, QZ, b)
        plus(a, b)
        self.append(sq, b, temporary)
        self.append(cmul, f, PZ, temporary, OX, OPZ)
        self.append(sq.inverse(), b, temporary)
        minus(a, b)
        minus(a, b)
        self.append(sq, b, temporary)
        self.append(cmul, f, PX, temporary, OZ, OPZ)
        self.append(sq.inverse(), b, temporary)
        plus(a, b)
        self.append(mul.inverse(), Z, QZ, b)
        reflect(X, Z)
        reflect(QX, QZ)
        self.append(mul.inverse(), Z, QZ, a)
        minus(QX, QZ)
        self.append(neg, Z)
        plus(X, Z)
        self.release_anc(a, b)
        self.x(f[0])
        for src, target in [(X, OX), (Z, OZ)]:
            for a, b in zip(src, target):
                self.ccx(f[0], a, b)
        self.x(f[0])
        for old, current, target in [(PX, X, OPX), (PZ, Z, OPZ)]:
            self.cx_reg(old, target)
            for a, b, c in zip(old, current, target):
                self.ccx(f[0], a, c)
                self.ccx(f[0], b, c)
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        self.release_anc(f)


@memoize
class BonnetainProductUpdate(Circuit):
    """Accumulate Y=X-Z and T=X+Z with the shared controlled IP multiplier."""

    def __init__(self, p, bits, index):
        super().__init__()
        t = ModIntType(p)
        X, Z = self.add_input_output(t), self.add_input_output(t)
        ell = self.add_input_output(UIntType(bits))
        N, D = self.add_input_output(t), self.add_input_output(t)
        add = ModularAdder(p, backends=UNITARY)
        if index == 0:
            self.cx_reg(X, N)
            self.append(add.inverse(), Z, N)
            self.cx_reg(X, D)
            self.append(add, Z, D)
            return
        f = self.get_anc(BoolType())
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        mult = ControlledModularIPMult(p, backends=UNITARY)
        self.append(add.inverse(), Z, X)
        self.append(mult, f, X, N)
        self.append(add, Z, X)
        self.append(add, X, Z)
        self.append(mult, f, Z, D)
        self.append(add.inverse(), X, Z)
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        self.release_anc(f)


@memoize
class BonnetainEighthPower(Circuit):
    """Five-square value interface from Appendix B.3, with two work fields."""

    def __init__(self, p):
        super().__init__()
        t = ModIntType(p)
        x, out = self.add_input_output(t), self.add_input_output(t)
        x2, x4 = self.get_anc(t), self.get_anc(t)
        sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
        self.append(sq, x, x2)
        self.append(sq, x2, x4)
        self.append(sq, x4, out)
        self.append(sq.inverse(), x2, x4)
        self.append(sq.inverse(), x, x2)
        self.release_anc(x2, x4)


@memoize
class BonnetainImageTail(Circuit):
    """Appendix B.3's in-place a',d' and inversion/IP-multiply quotient schedule."""

    def __init__(self, p):
        super().__init__()
        t = ModIntType(p)
        ae, de, N, D, out = [self.add_input_output(t) for _ in range(5)]
        y8, t8 = self.get_anc(t), self.get_anc(t)
        power = BonnetainEighthPower(p, backends=UNITARY)
        self.append(power, N, y8)
        self.append(power, D, t8)
        ipmul = ModularIPMult(p, backends=UNITARY)
        add = ModularAdder(p, backends=UNITARY)
        self.append(ipmul, t8, ae)
        self.append(ipmul, y8, de)
        # ae=a'-d'; compute its inverse, then ae=a'+d'.
        self.append(add.inverse(), de, ae)
        inverse = self.get_anc(t)
        inv = ModularInverseAdd(p, backends=UNITARY)
        self.append(inv, ae, inverse)
        self.append(add, de, ae)
        self.append(add, de, ae)
        self.append(ipmul, inverse, ae)
        double = ModularDouble(p, backends=UNITARY)
        self.append(double, ae)
        self.cx_reg(ae, out)
        self.append(double.inverse(), ae)
        self.append(ipmul.inverse(), inverse, ae)
        self.append(add.inverse(), de, ae)
        self.append(add.inverse(), de, ae)
        self.append(inv.inverse(), ae, inverse)
        self.release_anc(inverse)
        self.append(add, de, ae)
        self.append(ipmul.inverse(), y8, de)
        self.append(ipmul.inverse(), t8, ae)
        self.append(power.inverse(), D, t8)
        self.append(power.inverse(), N, y8)
        self.release_anc(y8, t8)


@memoize
class BonnetainKernelAdapter(Circuit):
    """Same A,QX,QZ,ell,enable,out interface as the revised kernel adapter."""

    def __init__(self, p, primes, s=15):
        super().__init__()
        t = ModIntType(p)
        bits = max(primes).bit_length()
        A, QX, QZ = [self.add_input_output(t) for _ in range(3)]
        ell = self.add_input_output(UIntType(bits))
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(UIntType(p.bit_length()))
        e = Expr(p, 1)
        a = e.add(0, e.const(2))
        d = e.sub(0, e.const(2))
        a24 = e.scale(a, pow(4, -1, p))
        initial = UnitaryInitialExpression(e, [a, d, a24])
        ar, dr, a24r = [self.get_anc(t) for _ in range(3)]
        self.append(initial, A, ar, dr, a24r)
        ae, de = self.get_anc(t), self.get_anc(t)
        power = UnitaryVisitingPower(p, bits, 4 * s, backends=UNITARY)
        self.append(power, ar, ell, ae)
        self.append(power, dr, ell, de)
        result = self.get_anc(t)
        walk = BonnetainProductWalk(p, primes, s, backends=UNITARY)
        self.append(walk, a24r, QX, QZ, ell, ae, de, result)
        for a, b, target in zip(A, result, out):
            self.cx(a, target)
            self.ccx(enable[0], a, target)
            self.ccx(enable[0], b, target)
        self.append(walk.inverse(), a24r, QX, QZ, ell, ae, de, result)
        self.release_anc(result)
        self.append(power.inverse(), dr, ell, de)
        self.append(power.inverse(), ar, ell, ae)
        self.release_anc(ae, de)
        self.append(initial.inverse(), A, ar, dr, a24r)
        self.release_anc(ar, dr, a24r)


@memoize
class BonnetainProductWalk(Circuit):
    """Encapsulated shared Bennett walk and tail, with a clean value target."""

    def __init__(self, p, primes, s=15):
        super().__init__()
        t = ModIntType(p)
        a24, QX, QZ = [self.add_input_output(t) for _ in range(3)]
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        ae, de, out = [self.add_input_output(t) for _ in range(3)]
        import isogeny

        isogeny.append_fused_walk(
            self,
            p,
            primes,
            s,
            a24,
            QX,
            QZ,
            ell,
            lambda p, bits, i: BonnetainPointStep(p, bits, i, backends=UNITARY),
            lambda N, D: self.append(BonnetainImageTail(p), ae, de, N, D, out),
            product_class=BonnetainProductUpdate,
        )
