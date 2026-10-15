"""Independent affine enumeration and ideal-rotation finite-field error audit.

This analysis does not include rotation-synthesis, arithmetic-circuit or hardware
errors. It applies on the promised supersingular cyclic Montgomery domain, with
nonzero x sampled uniformly and cofactor (p+1)/ell.
"""

from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "independent_affine_reference", ROOT / "tests/affine_reference.py"
)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
sys.path.insert(0, str(ROOT / "code"))
from parameters import P512


def success(p, ell):
    g = Fraction((p + 1) * (ell - 1), 2 * ell * (p - 1))
    return 1 - (1 - g) * Fraction(4, (p - 1) ** 2)


rows = []
total_candidates = 0
for p, primes in [(11, (3,)), (59, (3, 5)), (419, (3, 5, 7))]:
    curves = ref.valid_curves(p)
    for ell in primes:
        counts = []
        for A in curves:
            good = 0
            for x in range(1, p):
                P = ref.point_from_x(A, x, p, 1)
                if (
                    P is not None
                    and ref.multiply((p + 1) // ell, P, A, p, 1) is not None
                ):
                    good += 1
                total_candidates += 1
            counts.append(good)
        predicted = (p + 1) * (ell - 1) // (2 * ell)
        assert all(x == predicted for x in counts)
        rows.append(
            dict(
                p=p,
                ell=ell,
                valid_curves=len(curves),
                affine_good_x_counts=counts,
                expected_good_x_count=predicted,
                preparation_success_exact=str(success(p, ell)),
                signed_inplace_success=float(success(p, ell) ** 4),
            )
        )

M = 1300
failure_bound = Fraction(32 * M, 3 * (P512 - 1) ** 2)
result = dict(
    all_affine_candidate_checks=total_candidates,
    rows=rows,
    csidh512=dict(
        n=P512.bit_length(),
        M=M,
        ideal_rotation_failure_bound_numerator=str(failure_bound.numerator),
        ideal_rotation_failure_bound_denominator=str(failure_bound.denominator),
        log2_ideal_rotation_failure_bound=math.log2(failure_bound.numerator)
        - math.log2(failure_bound.denominator),
    ),
)
path = ROOT / "results/preparation_audit.json"
path.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
