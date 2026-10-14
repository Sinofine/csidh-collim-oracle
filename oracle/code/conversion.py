"""Montgomery/Edwards conversion and kernel-to-curve evaluation."""

from expressions import Expr
from isogeny import UnitaryLocalValueImage
from ladder import UnitaryValueLadder
from points import UNITARY
from points import UnitaryInitialExpression as ExpressionCircuit
from points import UnitaryInitialExpression as UnitaryExpressionCircuit
from qarton.circuit import BoolType, Circuit, UIntType, memoize
from qarton.modular_arithmetic import ModIntType


@memoize
class UnitaryMontgomeryToW(Circuit):
    def __init__(self, p):
        super().__init__()
        e = Expr(p, 3)
        A, X, Z = 0, 1, 2
        xz = e.mul(X, Z)
        W = e.mul(e.mul(e.sub(A, e.const(2)), xz), e.sq(e.sub(X, Z)))
        V = e.mul(e.add(e.add(e.sq(X), e.mul(A, xz)), e.sq(Z)), e.sq(e.add(X, Z)))
        regs = [self.add_input_output(ModIntType(p)) for _ in range(5)]
        self.append(UnitaryExpressionCircuit(e, [W, V]), *regs)


MontgomeryToW = UnitaryMontgomeryToW


@memoize
class UnitaryValueKernelAdapter(Circuit):
    def __init__(self, p, primes, s=15):
        super().__init__()
        t = ModIntType(p)
        A, X, Z = [self.add_input_output(t) for _ in range(3)]
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(UIntType(p.bit_length()))
        W = self.get_anc(t)
        V = self.get_anc(t)
        conv = MontgomeryToW(p)
        self.append(conv, A, X, Z, W, V)
        e = Expr(p, 1)
        d = e.mul(e.sub(0, e.const(2)), e.inv(e.add(0, e.const(2))))
        to_d = ExpressionCircuit(e, [d])
        dr = self.get_anc(t)
        self.append(to_d, A, dr)
        dn = self.get_anc(t)
        local = UnitaryLocalValueImage(p, primes, s)
        self.append(local, dr, W, V, ell, dn)
        e = Expr(p, 1)
        b = e.scale(e.mul(e.add(e.const(1), 0), e.inv(e.sub(e.const(1), 0))), 2)
        to_b = ExpressionCircuit(e, [b])
        B = self.get_anc(t)
        self.append(to_b, dn, B)
        for a, b, o in zip(A, B, out):
            self.cx(a, o)
            self.ccx(enable[0], a, o)
            self.ccx(enable[0], b, o)
        self.append(to_b.inverse(), dn, B)
        self.release_anc(B)
        self.append(local.inverse(), dr, W, V, ell, dn)
        self.release_anc(dn)
        self.append(to_d.inverse(), A, dr)
        self.release_anc(dr)
        self.append(conv.inverse(), A, X, Z, W, V)
        self.release_anc(W, V)


@memoize
class UnitaryValueImageFromPoint(Circuit):
    def __init__(self, p, primes, s=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        x = self.add_input_output(t)
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(UIntType(p.bit_length()))
        X = self.get_anc(t)
        Z = self.get_anc(t)
        ladder = UnitaryValueLadder(p, p.bit_length(), s, backends=UNITARY)
        adapter = UnitaryValueKernelAdapter(p, primes, s)
        self.append(ladder, A, x, cofac, X, Z)
        self.append(adapter, A, X, Z, ell, enable, out)
        self.append(ladder.inverse(), A, x, cofac, X, Z)
        self.release_anc(X, Z)
