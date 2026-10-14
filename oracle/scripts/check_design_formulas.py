"""Small-field algebra checks for the two isogeny formulas, not gate simulation."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def add(P, Q, A, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y + v) % p == 0:
        return None
    slope = (
        (3 * x * x + 2 * A * x + 1) * pow(2 * y, -1, p)
        if P == Q
        else (v - y) * pow(u - x, -1, p)
    ) % p
    z = (slope * slope - A - x - u) % p
    return z, (slope * (x - z) - y) % p


def mul(k, P, A, p):
    R = None
    while k:
        if k & 1:
            R = add(R, P, A, p)
        P = add(P, P, A, p)
        k >>= 1
    return R


rows = []
for p, degrees in [(11, [3]), (19, [5]), (59, [3, 5])]:
    checked = excluded = 0
    for A in range(p):
        if (A * A - 4) % p == 0:
            continue
        points = [
            (x, y)
            for x in range(p)
            for y in range(p)
            if (y * y - x**3 - A * x * x - x) % p == 0
        ]
        if len(points) + 1 != p + 1:
            continue
        for ell in degrees:
            generators = {mul((p + 1) // ell, P, A, p) for P in points}
            generators.discard(None)
            for Q in generators:
                assert mul(ell, Q, A, p) is None
                Y = T = 1
                h = 1
                try:
                    for i in range(1, (ell + 1) // 2):
                        X = mul(i, Q, A, p)[0]
                        Z = 1
                        Y = Y * (X - Z) % p
                        T = T * (X + Z) % p
                        W = (A - 2) * X * Z * (X - Z) ** 2 % p
                        V = (X * X + A * X * Z + Z * Z) * (X + Z) ** 2 % p
                        h = h * (W + V) * pow(2 * V, -1, p) % p
                    aprime = pow(A + 2, ell, p) * pow(T, 8, p) % p
                    dprime = pow(A - 2, ell, p) * pow(Y, 8, p) % p
                    old = 2 * (aprime + dprime) * pow(aprime - dprime, -1, p) % p
                    d = (A - 2) * pow(A + 2, -1, p) % p
                    dn = pow(d, ell, p) * pow(h, 8, p) % p
                    new = 2 * (1 + dn) * pow(1 - dn, -1, p) % p
                except ValueError:
                    excluded += 1
                    continue
                assert old == new, (p, A, ell, Q, old, new)
                checked += 1
    rows.append(
        dict(p=p, degrees=degrees, checked=checked, excluded_zero_denominator=excluded)
    )
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output", type=Path, default=ROOT / "results/comparison/formula-checks.json"
)
a = parser.parse_args()
a.output.parent.mkdir(parents=True, exist_ok=True)
result = dict(
    rows=rows,
)
a.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
