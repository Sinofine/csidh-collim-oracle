"""Native coherent test of the real 74-prime, M=1300 selector interface.

This checks the cryptographic-parameter selector, not the full 511-bit action.
It does not implement or validate a relation-lattice / Babai front end.
"""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from test_complete import unitary_check
from group_action import UnitaryActionSelector
from parameters import P512, PRIMES
from points import UNITARY
from qarton_compat import collision_regression

M = 1300
configs = []
for index in range(len(PRIMES)):
    for sign in (-1, 1):
        e = [0] * len(PRIMES)
        e[index] = sign * M
        configs.append((tuple(e), M - 1))
mixed = [0] * len(PRIMES)
mixed[0] = 2
mixed[2] = -3
mixed[-1] = 4
configs += [(tuple(mixed), t) for t in (0, 1, 2, 4, 5, 8, 9, M - 1)]
full = [0] * len(PRIMES)
full[0] = -600
full[-1] = 700
configs += [(tuple(full), t) for t in (0, 599, 600, M - 1)]
configs += [(tuple([0] * len(PRIMES)), t) for t in (0, M - 1)]
rows = []
outputs = []
for exponents, time in configs:
    assert sum(map(abs, exponents)) <= M
    flat = tuple(v for x in exponents for v in (x < 0, abs(x)))
    row = (*flat, time, 0, 0, False, False)
    lo = 0
    selected = None
    for ell, e in zip(PRIMES, exponents):
        if lo <= time < lo + abs(e):
            selected = (ell, e < 0, True)
            break
        lo += abs(e)
    ell, sign, enable = selected or (PRIMES[0], False, False)
    rows.append(row)
    outputs.append((*flat, time, ell, (P512 + 1) // ell, sign, enable))
qc = UnitaryActionSelector(P512, PRIMES, M, backends=UNITARY)
result = unitary_check(qc, rows, outputs, "selector_511_bits_74_primes_M1300")
result["bitvector_collision_regression"] = collision_regression()
output = ROOT / "results/selector512_validation.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)
