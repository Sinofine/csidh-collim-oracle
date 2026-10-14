"""Signed CSIDH actions and the complete group-action circuit."""

from conversion import UnitaryValueImageFromPoint
from points import UNITARY
from preparation import UnitaryArithmeticPrepare, UnitaryArithmeticUnprepare
from qarton.arithmetic import qc_add_uint, qc_lt_uint, qc_sub_uint
from qarton.circuit import BoolType, Circuit, UIntType, memoize, qc_reg_cast
from qarton.modular_arithmetic import ModIntType
from qarton.modular_arithmetic.modular_neg import ModularNeg
from utils import SetFlag, cswap


@memoize
class UnitaryActionSelector(Circuit):
    """e sign/magnitude pairs, round -> same inputs, ell,cofactor,sign,enable.
    Outputs must start zero. Select i with sum_{j<i}|e_j| <= round < sum_{j<=i}|e_j|.
    Promise sum magnitudes <= M. Inactive slots receive a safe valid default degree.
    """

    def __init__(self, p, primes, M):
        super().__init__()
        b = M.bit_length()
        ellbits = max(primes).bit_length()
        n = p.bit_length()
        exps = [
            (self.add_input_output(BoolType()), self.add_input_output(UIntType(b)))
            for _ in primes
        ]
        time = self.add_input_output(UIntType(b))
        ell = self.add_input_output(UIntType(ellbits))
        cofac = self.add_input_output(UIntType(n))
        sign = self.add_input_output(BoolType())
        enable = self.add_input_output(BoolType())

        def compare(a, b, out):
            if len(a) == len(b) == 1:
                self.x(a[0])
                self.ccx(a[0], b[0], out[0])
                self.x(a[0])
            else:
                qc_lt_uint(a, b, out)

        prefixes = []
        for _, mag in exps:
            s = self.get_anc(UIntType(b))
            if prefixes:
                self.cx_reg(prefixes[-1], s)
            qc_add_uint(mag, s)
            prefixes.append(s)
        default = (primes[0], (p + 1) // primes[0])
        for out, value in zip((ell, cofac), default):
            for j in range(len(out)):
                if (value >> j) & 1:
                    self.x(out[j])
        lo = self.get_anc(BoolType())
        hi = self.get_anc(BoolType())
        pick = self.get_anc(BoolType())
        for i, (prime, (sg, mag)) in enumerate(zip(primes, exps)):
            if i:
                compare(time, prefixes[i - 1], lo)
            compare(time, prefixes[i], hi)
            self.cx(lo[0], pick[0])
            self.cx(hi[0], pick[0])
            for out, value, base in zip(
                (ell, cofac), (prime, (p + 1) // prime), default
            ):
                for j in range(len(out)):
                    if ((value ^ base) >> j) & 1:
                        self.cx(pick[0], out[j])
            self.ccx(pick[0], sg[0], sign[0])
            self.cx(pick[0], enable[0])
            self.cx(hi[0], pick[0])
            self.cx(lo[0], pick[0])
            compare(time, prefixes[i], hi)
            if i:
                compare(time, prefixes[i - 1], lo)
        self.release_anc(lo, hi, pick)
        for i in reversed(range(len(prefixes))):
            qc_sub_uint(exps[i][1], prefixes[i])
            if i:
                self.cx_reg(prefixes[i - 1], prefixes[i])
            self.release_anc(prefixes[i])


@memoize
class UnitaryValueOOP(Circuit):
    def __init__(self, p, primes, width=5, s=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        k = self.add_input_output(UIntType(p.bit_length()))
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(UIntType(p.bit_length()))
        x = self.get_anc(UIntType(p.bit_length()))
        pre = self.append(
            UnitaryArithmeticPrepare(p, primes, width, s, backends=UNITARY),
            A,
            ell,
            k,
            x,
        )[-1]
        self.append_controlled(
            UnitaryValueImageFromPoint(p, primes, s),
            A,
            ell,
            k,
            qc_reg_cast(x, t),
            enable,
            out,
            controls=[pre[0]],
        )
        back = self.append_controlled(
            UnitaryArithmeticUnprepare(p, primes, width, s, backends=UNITARY),
            A,
            ell,
            k,
            x,
            controls=[pre[0]],
        )[-1]
        ok = self.get_anc(BoolType())
        self.append_controlled(SetFlag(), ok, controls=[pre[0], back[0]])
        success = self.get_classical(BoolType())
        self.msr(ok[0], success[0])
        self.release_anc(ok)
        tmp = self.get_classical(BoolType())
        for bit in x:
            self.msr(bit, tmp[0])
        self.release_anc(x)
        self.release_classical(tmp, pre, back)
        self.add_output(success)


@memoize
class UnitaryNegIf(Circuit):
    """A,sign -> (-1)^sign A,sign, with clean native workspace."""

    def __init__(self, p):
        super().__init__()
        A = self.add_input_output(ModIntType(p))
        sign = self.add_input_output(BoolType())
        temp = self.get_anc(ModIntType(p))
        self.cx_reg(A, temp)
        neg = ModularNeg(p, backends=UNITARY)
        self.append(neg, temp)
        cswap(self, sign[0], A, temp)
        self.append(neg, temp)
        self.cx_reg(A, temp)
        self.release_anc(temp)


@memoize
class ValueSignedOOP(Circuit):
    """Field target must be zero or F(A); success computes F(A) or erases it."""

    def __init__(self, p, primes, width=5, s=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        k = self.add_input_output(UIntType(p.bit_length()))
        sign = self.add_input_output(BoolType())
        enable = self.add_input_output(BoolType())
        out = self.add_input_output(t)
        neg = UnitaryNegIf(p, backends=UNITARY)
        self.append(neg, A, sign)
        self.append(neg, out, sign)
        ok = self.append(
            UnitaryValueOOP(p, primes, width, s, backends=UNITARY),
            A,
            ell,
            k,
            enable,
            qc_reg_cast(out, UIntType(p.bit_length())),
        )[-1]
        self.append(neg, A, sign)
        self.append(neg, out, sign)
        self.add_output(ok)


@memoize
class ValueEraseWithInverseSignedOOP(Circuit):
    def __init__(self, p, primes, width=5, s=15):
        super().__init__()
        t = ModIntType(p)
        B = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        sign = self.add_input_output(BoolType())
        enable = self.add_input_output(BoolType())
        A = self.add_input_output(t)
        self.x(sign[0])
        ok = self.append(
            ValueSignedOOP(p, primes, width, s), B, ell, cofac, sign, enable, A
        )[-1]
        self.x(sign[0])
        self.add_output(ok)


@memoize
class SwapFields(Circuit):
    def __init__(self, p):
        super().__init__()
        a = self.add_input_output(ModIntType(p))
        b = self.add_input_output(ModIntType(p))
        self.swap_reg(a, b)


@memoize
class ValueSignedInplace(Circuit):
    """Clean on success, reject on failure; two OOP instruments per update."""

    def __init__(self, p, primes, width=5, s=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        sign = self.add_input_output(BoolType())
        enable = self.add_input_output(BoolType())
        B = self.get_anc(t)
        first = self.append(
            ValueSignedOOP(p, primes, width, s), A, ell, cofac, sign, enable, B
        )[-1]
        second = self.append_controlled(
            ValueEraseWithInverseSignedOOP(p, primes, width, s),
            B,
            ell,
            cofac,
            sign,
            enable,
            A,
            controls=[first[0]],
        )[-1]
        self.append_controlled(SwapFields(p), A, B, controls=[first[0], second[0]])
        q = self.get_anc(BoolType())
        self.append_controlled(SetFlag(), q, controls=[first[0], second[0]])
        ok = self.get_classical(BoolType())
        self.msr(q[0], ok[0])
        self.release_anc(q)
        tmp = self.get_classical(BoolType())
        for bit in B:
            self.msr(bit, tmp[0])
        self.release_anc(B)
        self.release_classical(tmp, first, second)
        self.add_output(ok)


ActionSelector = UnitaryActionSelector


@memoize
class ValueGroupAction(Circuit):
    """A,e_1..e_u,zero-output -> same, classical success. Promise sum|e|<=M.
    Per exponent: sign bit then magnitude. Uses exactly M shared-degree slots,
    with classical failure guarding all later isogeny updates. Failure outputs
    carry no action promise; success is an approximate postselected oracle.
    """

    def __init__(self, p, primes, M, width=5, s=15):
        super().__init__()
        t = ModIntType(p)
        b = M.bit_length()
        eb = max(primes).bit_length()
        A = self.add_input_output(t)
        exps = [
            (self.add_input_output(BoolType()), self.add_input_output(UIntType(b)))
            for _ in primes
        ]
        out = self.add_input_output(t)
        self.cx_reg(A, out)
        time = self.get_anc(UIntType(b))
        ell = self.get_anc(UIntType(eb))
        cofac = self.get_anc(UIntType(p.bit_length()))
        sign = self.get_anc(BoolType())
        enable = self.get_anc(BoolType())
        q = self.get_anc(BoolType())
        self.x(q[0])
        alive = self.get_classical(BoolType())
        self.msr(q[0], alive[0])
        self.release_anc(q)
        selector = ActionSelector(p, primes, M, backends=UNITARY)
        step = ValueSignedInplace(p, primes, width, s, backends=UNITARY)
        flat = [v for pair in exps for v in pair]
        last = 0
        for i in range(M):
            for j in range(b):
                if ((i ^ last) >> j) & 1:
                    self.x(time[j])
            last = i
            self.append(selector, *flat, time, ell, cofac, sign, enable)
            success = self.append_controlled(
                step, out, ell, cofac, sign, enable, controls=[alive[0]]
            )[-1]
            q = self.get_anc(BoolType())
            self.append_controlled(SetFlag(), q, controls=[alive[0], success[0]])
            nxt = self.get_classical(BoolType())
            self.msr(q[0], nxt[0])
            self.release_anc(q)
            self.release_classical(alive, success)
            alive = nxt
            self.append(selector.inverse(), *flat, time, ell, cofac, sign, enable)
        for j in range(b):
            if (last >> j) & 1:
                self.x(time[j])
        self.release_anc(time, ell, cofac, sign, enable)
        self.add_output(alive)
