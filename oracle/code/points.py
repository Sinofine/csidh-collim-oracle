"""Edwards doubling and differential point steps."""

from expressions import Expr
from qarton.arithmetic import qc_lt_uint
from qarton.arithmetic.comparator_adder import CDKMComparatorAdder
from qarton.binary_operations.and_ccx import AndGate, AndGateDecomposed, FakeAndGate
from qarton.circuit import BoolType, Circuit, UInt, UIntType, memoize
from qarton.circuit.backends import MemoryConstraint, PriorityBackends
from qarton.circuit.basic_gates import MSR
from qarton.modular_arithmetic import (
    ModInt,
    ModIntType,
    ModularAdder,
    ModularConstantProductAdd,
    ModularInverseAdd,
    ModularProductAdd,
    ModularSquareAdd,
)
from qarton.modular_arithmetic.controlled_modular_product import (
    ControlledModularProductAdd,
)

UNITARY = MemoryConstraint(
    None,
    priority=PriorityBackends(FakeAndGate, CDKMComparatorAdder),
    exclude=PriorityBackends(MSR, AndGate, AndGateDecomposed),
)


@memoize
class ExplicitUnitaryWPointStep(Circuit):
    def __init__(self, p, ellbits, index):
        super().__init__()
        t = ModIntType(p)
        A = self.add_input_output(t)  # noqa: F841 - shared point-step interface
        QW = self.add_input_output(t)
        QZ = self.add_input_output(t)
        ell = self.add_input_output(UIntType(ellbits))
        PW, PZ, W, Z = [self.add_input_output(t) for _ in range(4)]
        OPW, OPZ, OW, OZ = [self.add_input_output(t) for _ in range(4)]
        if index == 0:
            raise ValueError(
                "Differential steps require index > 0; use UnitaryInitialPointStep for index 0."
            )
        f = self.get_anc(BoolType())
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        # On disabled branches initialize output to unchanged current pair.
        self.x(f[0])
        for src, out in [(W, OW), (Z, OZ)]:
            for a, b in zip(src, out):
                self.ccx(f[0], a, b)
        self.x(f[0])
        a, b = [self.get_anc(t) for _ in range(2)]
        # OPW remains zero until its final copy: borrow it for the square temporary.
        c = OPW
        mul = ModularProductAdd(p, zero_input=True, backends=UNITARY)
        sq = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
        add = ModularAdder(p, backends=UNITARY)
        # General accumulator required: disabled outputs need not be zero.
        cmul = ControlledModularProductAdd(p, zero_input=False, backends=UNITARY)

        def coordinate(x, y, u, v, scale, out):
            self.append(mul, x, y, a)
            self.append(mul, u, v, b)
            self.append(add.inverse(), b, a)
            self.append(sq, a, c)
            self.append(cmul, f, scale, c, out)
            self.append(sq.inverse(), a, c)
            self.append(add, b, a)
            self.append(mul.inverse(), u, v, b)
            self.append(mul.inverse(), x, y, a)

        coordinate(W, QZ, QW, Z, PZ, OW)
        coordinate(W, QW, Z, QZ, PW, OZ)
        self.release_anc(a, b)
        # Copy previous coordinates selected from old/current pair.
        for old, new, out in [(PW, W, OPW), (PZ, Z, OPZ)]:
            self.cx_reg(old, out)
            for a, b, c in zip(old, new, out):
                self.ccx(f[0], a, c)
                self.ccx(f[0], b, c)
        qc_lt_uint(UInt(2 * index + 1), ell, f)
        self.release_anc(f)


def build_unitary_point_step(p, ellbits, index):
    return ExplicitUnitaryWPointStep(p, ellbits, index, backends=UNITARY)


class UnitaryInitialExpression(Circuit):
    """Compute-copy-uncompute lowering of a field SLP, with output XOR."""

    def __init__(self, expr, outputs):
        super().__init__()
        self.expr = expr
        self.outputs = tuple(outputs)
        p = expr.p
        t = ModIntType(p)
        regs = [self.add_input_output(t) for _ in range(expr.ninputs)]
        outs = [self.add_input_output(t) for _ in outputs]
        tape = []
        for kind, a in expr.nodes:
            z = self.get_anc(t)
            if kind == "const":
                indices = [i for i in range(p.bit_length()) if (a[0] >> i) & 1]
                for i in indices:
                    self.x(z[i])
                tape.append(("const", z, indices))
            else:
                if kind in ("add", "sub"):
                    self.cx_reg(regs[a[0]], z)
                    gate = ModularAdder(p, backends=UNITARY)
                    if kind == "sub":
                        gate = gate.inverse()
                    args = (regs[a[1]], z)
                    self.append(gate, *args)
                    tape.append(("add", z, (gate, args, regs[a[0]])))
                else:
                    if kind == "mul":
                        gate = ModularProductAdd(p, zero_input=True, backends=UNITARY)
                        args = (regs[a[0]], regs[a[1]], z)
                    elif kind == "square":
                        gate = ModularSquareAdd(p, zero_input=True, backends=UNITARY)
                        args = (regs[a[0]], z)
                    elif kind == "scale":
                        gate = ModularConstantProductAdd(
                            p, a[1], zero_input=True, backends=UNITARY
                        )
                        args = (regs[a[0]], z)
                    elif kind == "inv":
                        gate = ModularInverseAdd(p, backends=UNITARY)
                        args = (regs[a[0]], z)
                    else:
                        raise ValueError(kind)
                    self.append(gate, *args)
                    tape.append(("gate", z, (gate, args)))
            regs.append(z)
        for o, i in zip(outs, outputs):
            self.cx_reg(regs[i], o)
        for kind, z, data in reversed(tape):
            if kind == "const":
                for i in data:
                    self.x(z[i])
            elif kind == "add":
                gate, args, src = data
                self.append(gate.inverse(), *args)
                self.cx_reg(src, z)
            else:
                gate, args = data
                self.append(gate.inverse(), *args)
            self.release_anc(z)

    def dummy_classical_function(self, args):
        # Acceleration of an existing native arithmetic DAG, not a dummy-only gate.
        inputs = args[: self.expr.ninputs]
        targets = args[self.expr.ninputs :]
        values = self.expr.evaluate([x.v for x in inputs])
        result = []
        for target, index in zip(targets, self.outputs):
            value = target.v ^ values[index]
            if value >= self.expr.p:
                raise ValueError("XOR target outside canonical field domain")
            result.append(ModInt(value, self.expr.p))
        return tuple(inputs) + tuple(result)

    def dummy_classical_function_inverse(self, args):
        return self.dummy_classical_function(args)


@memoize
class UnitaryInitialPointStep(Circuit):
    def __init__(self, p, bits, index=0):
        if index != 0:
            raise ValueError("Initial step only")
        super().__init__()
        t = ModIntType(p)
        d = self.add_input_output(t)
        QW = self.add_input_output(t)  # noqa: F841 - shared point-step interface
        QZ = self.add_input_output(t)  # noqa: F841 - shared point-step interface
        ell = self.add_input_output(UIntType(bits))  # noqa: F841 - shared point-step interface
        src = [self.add_input_output(t) for _ in range(4)]
        out = [self.add_input_output(t) for _ in range(4)]
        e = Expr(p, 3)
        plus = e.sq(e.add(1, 2))
        minus = e.sq(e.sub(1, 2))
        delta = e.sub(plus, minus)
        dp = e.mul(0, plus)
        nw = e.mul(delta, e.sub(dp, delta))
        nz = e.mul(dp, minus)
        self.cx_reg(src[2], out[0])
        self.cx_reg(src[3], out[1])
        self.append(
            UnitaryInitialExpression(e, [nw, nz]), d, src[2], src[3], out[2], out[3]
        )
