"""Point sampling, character tests, and amplitude amplification."""

import math

from character import UnitaryPositiveCharacterFlag, append_chain
from expressions import Expr
from ladder import UnitaryCofactorFlag, append_ladder
from points import UNITARY, UnitaryInitialExpression
from qarton.arithmetic import qc_eq_uint
from qarton.binary_operations import qc_eq
from qarton.circuit import BoolType, Circuit, UInt, UIntType, memoize, qc_reg_cast
from qarton.modular_arithmetic import ModInt, ModIntType, ModularAdderConstant
from qarton.modular_arithmetic.uniform_sample import UniformSample
from utils import SetFlag


@memoize
class PaperDegreePhase(Circuit):
    def __init__(self, primes):
        super().__init__()
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        f = self.add_input_output(BoolType())
        equal = self.get_anc(BoolType())
        both = self.get_anc(BoolType())
        for l in primes:
            angle = math.acos(-1 / (l - 1))  # g=(1-1/l)/2; cos gamma=1-1/(2g)
            qc_eq_uint(UInt(l), ell, equal)
            self.ccx(equal[0], f[0], both[0])
            self.phase(angle, both[0])
            self.ccx(equal[0], f[0], both[0])
            qc_eq_uint(UInt(l), ell, equal)
        self.release_anc(equal, both)


@memoize
class UnitaryNonzeroUniform(Circuit):
    def __init__(self, p):
        super().__init__()
        raw = self.add_input_output(UIntType(p.bit_length()))
        x = qc_reg_cast(raw, ModIntType(p))
        u = UniformSample(p - 1, backends=UNITARY)
        assert u.nbr_qubits() == p.bit_length() and u.nbr_cbits() == 0
        for operation in u.iterate_operations():
            self.append_operation(operation)
        self.append(ModularAdderConstant(p, ModInt(1, p), backends=UNITARY), x)


def rhs_circuit(p):
    e = Expr(p, 2)
    x2 = e.sq(1)
    rhs = e.add(e.add(e.mul(x2, 1), e.mul(0, x2)), 1)
    return UnitaryInitialExpression(e, [rhs])


@memoize
class UnitaryCharacterFirstGoodPhase(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        x = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        chi = self.get_anc(BoolType())
        poly = rhs_circuit(p)
        character = UnitaryPositiveCharacterFlag(
            p, width, max(checkpoints, 4 * checkpoints - 16), backends=UNITARY
        )

        def character_call(inverse=False):
            r = self.get_anc(t)
            self.append(poly, A, x, r)
            enable = self.get_anc(BoolType())
            self.x(enable[0])
            self.append(character.inverse() if inverse else character, r, enable, chi)
            self.x(enable[0])
            self.release_anc(enable)
            self.append(poly.inverse(), A, x, r)
            self.release_anc(r)

        character_call()

        def terminal(src):
            nz = self.get_anc(BoolType())
            both = self.get_anc(BoolType())
            qc_eq(ModInt(0, p), src[1], nz)
            self.x(nz[0])
            self.ccx(nz[0], chi[0], both[0])
            self.append(PaperDegreePhase(primes, backends=UNITARY), ell, both)
            self.ccx(nz[0], chi[0], both[0])
            self.x(nz[0])
            qc_eq(ModInt(0, p), src[1], nz)
            self.release_anc(nz, both)

        append_ladder(self, p, p.bit_length(), checkpoints, A, x, cofac, terminal)
        character_call(True)
        self.release_anc(chi)


@memoize
class UnitaryDegreeNegativeCharacter(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=44):
        super().__init__()
        x = self.add_input_output(ModIntType(p))
        enable = self.add_input_output(BoolType())
        ell = self.add_input_output(UIntType(max(primes).bit_length()))

        def root(src):
            f = self.get_anc(BoolType())
            both = self.get_anc(BoolType())
            qc_eq(ModInt(p - 1, p), src, f)
            self.ccx(enable[0], f[0], both[0])
            self.append(PaperDegreePhase(primes, backends=UNITARY), ell, both)
            self.ccx(enable[0], f[0], both[0])
            qc_eq(ModInt(p - 1, p), src, f)
            self.release_anc(f, both)

        append_chain(self, x, p, width, checkpoints, root)


@memoize
class UnitaryDegreePositiveCharacter(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=44):
        super().__init__()
        x = self.add_input_output(ModIntType(p))
        enable = self.add_input_output(BoolType())
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        f = self.get_anc(BoolType())
        both = self.get_anc(BoolType())
        qc_eq(ModInt(0, p), x, f)
        self.x(f[0])
        self.ccx(enable[0], f[0], both[0])
        self.append(PaperDegreePhase(primes, backends=UNITARY), ell, both)
        self.ccx(enable[0], f[0], both[0])
        self.x(f[0])
        qc_eq(ModInt(0, p), x, f)
        self.release_anc(f, both)
        self.append(
            UnitaryDegreeNegativeCharacter(
                p, primes, width, checkpoints, backends=UNITARY
            ).inverse(),
            x,
            enable,
            ell,
        )


@memoize
class UnitaryCofactorFirstGoodPhase(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        x = self.add_input_output(t)
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        f = self.get_anc(BoolType())
        k = UnitaryCofactorFlag(
            p, max(checkpoints, p.bit_length().bit_length()), backends=UNITARY
        )
        self.append(k, A, x, cofac, f)
        r = self.get_anc(t)
        poly = rhs_circuit(p)
        self.append(poly, A, x, r)
        self.append(
            UnitaryDegreePositiveCharacter(
                p,
                primes,
                width,
                max(checkpoints, 4 * checkpoints - 16),
                backends=UNITARY,
            ),
            r,
            f,
            ell,
        )
        self.append(poly.inverse(), A, x, r)
        self.release_anc(r)
        self.append(k.inverse(), A, x, cofac, f)
        self.release_anc(f)


@memoize
class UnitaryCofactorFirstPrepare(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=15):
        super().__init__()
        A = self.add_input_output(ModIntType(p))
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        raw = self.add_input_output(UIntType(p.bit_length()))
        x = qc_reg_cast(raw, ModIntType(p))
        u = UnitaryNonzeroUniform(p, backends=UNITARY)
        self.append(u, raw)
        self.append(
            UnitaryCofactorFirstGoodPhase(
                p, primes, width, checkpoints, backends=UNITARY
            ),
            A,
            x,
            ell,
            cofac,
        )
        self.append(u.inverse(), raw)
        f = self.get_anc(BoolType())
        qc_eq_uint(UInt(0), raw, f)
        self.append(PaperDegreePhase(primes, backends=UNITARY), ell, f)
        qc_eq_uint(UInt(0), raw, f)
        self.release_anc(f)
        self.append(u, raw)


@memoize
class UnitaryCoherentPrepare(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=15):
        super().__init__()
        A = self.add_input_output(ModIntType(p))
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        cofac = self.add_input_output(UIntType(p.bit_length()))
        raw = self.add_input_output(UIntType(p.bit_length()))
        x = qc_reg_cast(raw, ModIntType(p))
        u = UnitaryNonzeroUniform(p, backends=UNITARY)
        self.append(u, raw)
        self.append(
            UnitaryCharacterFirstGoodPhase(
                p, primes, width, checkpoints, backends=UNITARY
            ),
            A,
            x,
            ell,
            cofac,
        )
        self.append(u.inverse(), raw)
        f = self.get_anc(BoolType())
        qc_eq_uint(UInt(0), raw, f)
        self.append(PaperDegreePhase(primes, backends=UNITARY), ell, f)
        qc_eq_uint(UInt(0), raw, f)
        self.release_anc(f)
        self.append(u, raw)


@memoize
class UnitaryArithmeticGoodDetector(Circuit):
    def __init__(self, p, width=5, checkpoints=15):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)
        x = self.add_input_output(t)
        cofac = self.add_input_output(UIntType(p.bit_length()))
        f = self.get_anc(BoolType())
        self.append(
            UnitaryCofactorFlag(p, checkpoints, backends=UNITARY), A, x, cofac, f
        )
        nonzero = self.get_classical(BoolType())
        self.msr(f[0], nonzero[0])
        self.release_anc(f)
        r = self.get_anc(t)
        poly = rhs_circuit(p)
        self.append(poly, A, x, r)
        enable = self.get_anc(BoolType())
        self.x(enable[0])
        positive = self.get_anc(BoolType())
        self.append(
            UnitaryPositiveCharacterFlag(
                p, width, max(checkpoints, 4 * checkpoints - 16), backends=UNITARY
            ),
            r,
            enable,
            positive,
        )
        self.x(enable[0])
        self.release_anc(enable)
        char = self.get_classical(BoolType())
        self.msr(positive[0], char[0])
        self.release_anc(positive)
        self.append(poly.inverse(), A, x, r)
        self.release_anc(r)
        ok = self.get_anc(BoolType())
        self.append_controlled(SetFlag(), ok, controls=[nonzero[0], char[0]])
        result = self.get_classical(BoolType())
        self.msr(ok[0], result[0])
        self.release_anc(ok)
        self.release_classical(nonzero, char)
        self.add_output(result)


@memoize
class UnitaryArithmeticPrepare(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=15):
        super().__init__()
        A = self.add_input_output(ModIntType(p))
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        k = self.add_input_output(UIntType(p.bit_length()))
        raw = self.add_input_output(UIntType(p.bit_length()))
        self.append(
            UnitaryCofactorFirstPrepare(
                p, primes, width, checkpoints, backends=UNITARY
            ),
            A,
            ell,
            k,
            raw,
        )
        success = self.append(
            UnitaryArithmeticGoodDetector(p, width, checkpoints, backends=UNITARY),
            A,
            qc_reg_cast(raw, ModIntType(p)),
            k,
        )[-1]
        self.add_output(success)


@memoize
class UnitaryArithmeticUnprepare(Circuit):
    def __init__(self, p, primes, width=5, checkpoints=15):
        super().__init__()
        A = self.add_input_output(ModIntType(p))
        ell = self.add_input_output(UIntType(max(primes).bit_length()))
        k = self.add_input_output(UIntType(p.bit_length()))
        raw = self.add_input_output(UIntType(p.bit_length()))
        self.append(
            UnitaryCoherentPrepare(
                p, primes, width, checkpoints, backends=UNITARY
            ).inverse(),
            A,
            ell,
            k,
            raw,
        )
        f = self.get_anc(BoolType())
        qc_eq_uint(UInt(0), raw, f)
        success = self.get_classical(BoolType())
        self.msr(f[0], success[0])
        self.release_anc(f)
        self.add_output(success)
