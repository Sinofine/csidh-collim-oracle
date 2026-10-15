# SPDX-License-Identifier: AGPL-3.0-or-later
# Sampling and rotation topology adapted from Qarton, Copyright (C) 2026.
"""Analytic phase targets, retaining integer ratios until synthesis.

The sampler and Y-rotation topology follow the pinned Qarton implementation
(AGPL-3.0-or-later), UniformSample / RotationY / ControlledRotationY /
ControlledHadamard. Only the representation of angles changes. See README.md for the distinction between native double simulation and
high-precision synthesis verification.
"""

from fractions import Fraction
from functools import lru_cache
import json
import mpmath as mp
from qarton.circuit import (
    BoolType,
    Circuit,
    InPlaceCircuit,
    InvolutoryCircuit,
    memoize,
    qc_bool_reg,
    qc_reg_from_bits,
)
from qarton.circuit.encoding import EncodableObject
from qarton.modular_arithmetic import ModIntType


class Angle(EncodableObject):
    """scale * pi, acos(n/d), or acos(sqrt(n/d)); immutable by convention."""

    def __init__(self, kind, n=0, d=1, sn=1, sd=1):
        assert kind in ("pi", "acos", "acos_sqrt")
        v, scale = Fraction(n, d), Fraction(sn, sd)
        self.kind, self.n, self.d = kind, v.numerator, v.denominator
        self.sn, self.sd = scale.numerator, scale.denominator

    def to_data(self):
        return (self.kind, self.n, self.d, self.sn, self.sd)

    @classmethod
    def from_data(cls, data):
        return cls(*data)

    def __str__(self):
        return "A" + json.dumps(self.to_data(), separators=(",", ":"))

    @classmethod
    def from_key(cls, key):
        assert key.startswith("A")
        return cls.from_data(json.loads(key[1:]))

    def __neg__(self):
        return Angle(self.kind, self.n, self.d, -self.sn, self.sd)

    def __truediv__(self, divisor):
        assert isinstance(divisor, int) and divisor
        return Angle(self.kind, self.n, self.d, self.sn, self.sd * divisor)

    def __mul__(self, scalar):
        if isinstance(scalar, int):
            return Angle(self.kind, self.n, self.d, self.sn * scalar, self.sd)
        # Native Qarton simulation uses cmath / double complex, explicitly.
        if isinstance(scalar, complex):
            return float(self) * scalar
        return NotImplemented

    __rmul__ = __mul__

    def value(self, dps=160):
        with mp.workdps(dps):
            if self.kind == "pi":
                base = mp.pi
            elif self.kind == "acos":
                base = mp.acos(mp.mpf(self.n) / self.d)
            else:
                base = mp.acos(mp.sqrt(mp.mpf(self.n) / self.d))
            return base * self.sn / self.sd

    def __float__(self):
        return _angle_float(str(self))


@lru_cache(maxsize=None)
def _angle_float(key):
    return float(Angle.from_key(key).value())


PI = Angle("pi")


@memoize
class ExactRY(InPlaceCircuit):
    def __init__(self, theta):
        super().__init__()
        r = self.add_input_output(BoolType())
        self.theta = theta
        self.phase(-PI / 2, r[0])
        self.h(r[0])
        self.phase(theta, r[0])
        self.h(r[0])
        self.gp(-theta / 2)
        self.phase(PI / 2, r[0])

    def define_inverse(self):
        return ExactRY(-self.theta)


@memoize
class ExactCRY(InPlaceCircuit):
    def __init__(self, theta):
        super().__init__()
        r = self.add_input_output(2)
        self.theta = theta
        self.cp(-PI / 2, r[0], r[1])
        self.h(r[1])
        self.cp(theta, r[0], r[1])
        self.h(r[1])
        self.phase(-theta / 2, r[0])
        self.cp(PI / 2, r[0], r[1])

    def define_inverse(self):
        return ExactCRY(-self.theta)


@memoize
class ExactCH(InvolutoryCircuit):
    def __init__(self):
        super().__init__()
        r = self.add_input_output(2)
        self.cp(-PI / 2, r[0], r[1])
        self.h(r[1])
        self.cp(-PI / 4, r[0], r[1])
        self.cx(r[0], r[1])
        self.cp(PI / 4, r[0], r[1])
        self.h(r[1])
        self.cp(PI / 2, r[0], r[1])


@memoize
class ExactUniformSample(Circuit):
    """Same algorithm/topology as Qarton UniformSample, with analytic ratios."""

    def __init__(self, p):
        super().__init__()
        r = self.add_anc_output(ModIntType(p))
        append_uniform_sample(self, r, p)


def append_uniform_sample(self, r, p):
    """Append onto an existing register, letting Qarton allocate gadget scratch."""
    ones = [i for i in range(p.bit_length()) if (p >> i) & 1]
    mm = 1 << ones[0]
    if len(ones) > 1:
        for i in ones[1:]:
            self.x(r[i])
    for i in range(ones[0]):
        self.h(r[i])
    if len(ones) == 1:
        return
    theta = Angle("acos_sqrt", mm, p, -2)
    self.append(ExactRY(theta), qc_bool_reg(r[ones[1]]))
    self.x(r[ones[1]])
    for i in range(ones[0], ones[1]):
        self.append(ExactCH(), qc_reg_from_bits(r[ones[1]], r[i]))
    self.x(r[ones[1]])
    for m in range(1, len(ones) - 1):
        theta = Angle("acos_sqrt", 1 << ones[m], p - mm, -2)
        self.x(r[ones[m]])
        self.append(ExactCRY(theta), qc_reg_from_bits(r[ones[m]], r[ones[m + 1]]))
        self.x(r[ones[m]])
        self.x(r[ones[m + 1]])
        for i in range(ones[m], ones[m + 1]):
            self.append(ExactCH(), qc_reg_from_bits(r[ones[m + 1]], r[i]))
        self.x(r[ones[m + 1]])
        mm += 1 << ones[m]
