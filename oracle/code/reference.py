"""Classical evaluation callbacks for the Montgomery ladder."""

from qarton.modular_arithmetic import ModInt


def ladder_reference(self, args):
    A, x, scalar, ox, oz = args
    p = self.p
    a24 = (A.v + 2) * pow(4, -1, p) % p
    r0 = (1, 0)
    r1 = (x.v, 1)
    for i in reversed(range(self.scalar_bits)):
        bit = (scalar >> i) & 1
        if bit:
            r0, r1 = r1, r0
        X, Z = r0
        U, V = r1
        a = (X + Z) ** 2 % p
        b = (X - Z) ** 2 % p
        d = (a - b) % p
        doubled = (a * b % p, d * (b + a24 * d) % p)
        u = (X + Z) * (U - V) % p
        v = (X - Z) * (U + V) % p
        added = ((u + v) ** 2 % p, x.v * (u - v) ** 2 % p)
        r0, r1 = doubled, added
        if bit:
            r0, r1 = r1, r0
    values = (ox.v ^ r0[0], oz.v ^ r0[1])
    if max(values) >= p:
        raise ValueError("XOR target outside canonical field domain")
    return A, x, scalar, ModInt(values[0], p), ModInt(values[1], p)


def ladder_reference_inverse(self, args):
    return self.dummy_classical_function(args)
