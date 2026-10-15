"""Independent affine reference plus published x-only / Edwards formulas.

Small parameters only for exhaustive curve enumeration. No cryptographic use.
The Montgomery reference uses affine points, not the project's circuit ladder.
"""

import hashlib
from functools import lru_cache
from math import prod


def inv(x, p):
    if x % p == 0:
        raise ZeroDivisionError("zero field denominator")
    return pow(x, -1, p)


def character(x, p):
    y = pow(x % p, (p - 1) // 2, p)
    return -1 if y == p - 1 else y


def rhs(A, x, p):
    return (x**3 + A * x * x + x) % p


def add(P, Q, A, p, twist=1):
    """Group law on twist*y^2 = x^3 + A*x^2 + x; infinity is None."""
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y + v) % p == 0:
        return None
    if P == Q:
        slope = (3 * x * x + 2 * A * x + 1) * inv(2 * twist * y, p) % p
    else:
        slope = (v - y) * inv(u - x, p) % p
    z = (twist * slope * slope - A - x - u) % p
    return z, (slope * (x - z) - y) % p


def multiply(k, P, A, p, twist=1):
    R = None
    while k:
        if k & 1:
            R = add(R, P, A, p, twist)
        P = add(P, P, A, p, twist)
        k >>= 1
    return R


def prime_factors(n):
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def point_from_x(A, x, p, twist):
    f = rhs(A, x, p) * twist % p  # twist is +1 or -1
    if character(f, p) != 1:
        return None
    y = pow(f, (p + 1) // 4, p)
    assert y * y % p == f
    return x, y


@lru_cache(maxsize=None)
def good_points(A, p, twist):
    """One sign of y for every x of EXACT order p+1; curve and twist."""
    n = p + 1
    result = []
    for x in range(p):
        P = point_from_x(A, x, p, twist)
        if P is not None and multiply(n, P, A, p, twist) is None:
            if all(
                multiply(n // q, P, A, p, twist) is not None for q in prime_factors(n)
            ):
                result.append(P)
    return tuple(result)


def valid_curves(p):
    """Supersingular Montgomery curves with one rational point of order two."""
    result = []
    for A in range(p):
        if character(A * A - 4, p) != -1:
            continue
        count = 1 + sum(1 + character(rhs(A, x, p), p) for x in range(p))
        if count == p + 1:
            result.append(A)
    return result


def xdouble(P, A, p):
    X, Z = P
    a = (X + Z) ** 2 % p
    b = (X - Z) ** 2 % p
    c = (a - b) % p
    return a * b % p, c * (b + (A + 2) * inv(4, p) * c) % p


def xadd(P, Q, D, p):
    X, Z = P
    U, V = Q
    DX, DZ = D
    a = (X + Z) * (U - V) % p
    b = (X - Z) * (U + V) % p
    return DZ * (a + b) ** 2 % p, DX * (a - b) ** 2 % p


def xladder(k, x, A, p):
    R, S, P = (1, 0), (x, 1), (x, 1)
    for bit in bin(k)[2:]:
        if bit == "0":
            R, S = xdouble(R, A, p), xadd(R, S, P, p)
        else:
            R, S = xadd(R, S, P, p), xdouble(S, A, p)
    return R


@lru_cache(maxsize=None)
def select_good_point(A, p, twist):
    """Deterministic functional substitute for good-point state preparation.

    Exhaustive in toy fields; for large p use bounded hash-derived candidates,
    characterize curve/twist and certify full order using x-only ladders.
    This routine is not a quantum state preparation or a resource estimate.
    """
    if p < 4096:
        pts = good_points(A, p, twist)
        if not pts:
            raise ValueError("No full-order point on this curve/twist")
        return pts[0]
    factors = prime_factors(p + 1)
    for attempt in range(1024):
        material = f"csidh-oracle-test:{p}:{A}:{twist}:{attempt}".encode()
        x = (
            int.from_bytes(
                hashlib.shake_256(material).digest((p.bit_length() + 7) // 8), "little"
            )
            % p
        )
        if character(rhs(A, x, p), p) != twist:
            continue
        if any(xladder((p + 1) // q, x, A, p)[1] == 0 for q in factors):
            continue
        X, Z = xladder(p + 1, x, A, p)
        if Z == 0 and X != 0:
            P = point_from_x(A, x, p, twist)
            assert P is not None
            return P
    raise RuntimeError("Bounded deterministic good-point search failed")


def montgomery_image(A, kernel_xs, ell, p):
    """Moriya--Onuki--Takagi (2019/843), equation (4)."""
    ap = pow(A + 2, ell, p) * pow(prod((x + 1) % p for x in kernel_xs), 8, p) % p
    dp = pow(A - 2, ell, p) * pow(prod((x - 1) % p for x in kernel_xs), 8, p) % p
    return 2 * (ap + dp) * inv(ap - dp, p) % p


def w_from_montgomery(A, x, p):
    """w=d*u_E^2*y_E^2, d=(A-2)/(A+2), algebraic birational map.

    u_E can lie in an extension; its square and w are in F_p here.
    This is NOT an assertion of an F_p-isomorphism to the normalized Edwards
    curve for every A. Circuit curve outputs are Montgomery coefficients.
    """
    return ((A - 2) * x * (x - 1) ** 2 * inv((x * x + A * x + 1) * (x + 1) ** 2, p)) % p


def edwards_image(A, kernel_xs, ell, p):
    d = (A - 2) * inv(A + 2, p) % p
    ws = [w_from_montgomery(A, x, p) for x in kernel_xs]
    h = prod((w + 1) * inv(2, p) % p for w in ws) % p
    dnew = pow(d, ell, p) * pow(h, 8, p) % p
    B = 2 * (1 + dnew) * inv(1 - dnew, p) % p
    return B, dnew, ws


def j_montgomery(A, p):
    return 256 * pow(A * A - 3, 3, p) * inv(A * A - 4, p) % p


def velu_j(A, kernel, p, twist=1):
    """Independent short-Weierstrass Velu coefficient sums, all nonzero points.

    Starting from twist*y^2=f(x), use X=twist*(x+A/3), Y=twist^2*y.
    Compare j-invariants; this check alone does not certify a Montgomery sign.
    """
    third = inv(3, p)
    a = twist * twist * (1 - A * A * third) % p
    b = twist**3 * (2 * A**3 * inv(27, p) - A * third) % p
    xs = [twist * (P[0] + A * third) % p for P in kernel]
    t = sum(3 * x * x + a for x in xs) % p
    w = sum(5 * x**3 + 3 * a * x + 2 * b for x in xs) % p
    anew, bnew = (a - 5 * t) % p, (b - 7 * w) % p
    return 1728 * 4 * anew**3 * inv(4 * anew**3 + 27 * bnew * bnew, p) % p


def reference_step(A, ell, direction, p):
    P = select_good_point(A, p, direction)
    Q = multiply((p + 1) // ell, P, A, p, direction)
    half = [multiply(i, Q, A, p, direction)[0] for i in range(1, (ell + 1) // 2)]
    return montgomery_image(A, half, ell, p)


def reference_action(A, exponents, primes, p):
    for e, ell in zip(exponents, primes):
        for _ in range(abs(e)):
            A = reference_step(A, ell, 1 if e > 0 else -1, p)
    return A
