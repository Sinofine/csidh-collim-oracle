"""Native small-field validation of the composed, heralded CSIDH oracle.

Qarton's gate evolution is unmodified and every dummy shortcut is disabled.
Only its measurement-outcome sampler is controlled, to cover success and rare
failure paths reproducibly. No zero-probability measurement outcome is forced.
The resulting branch probabilities are recorded; this is not a shot-frequency
estimate and does not certify a cryptographic-size quantum simulation.
"""

import argparse
import cmath
from contextlib import contextmanager
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from arithmetic import install_before_component_imports

install_before_component_imports()
from points import UNITARY
from preparation import (
    UnitaryNonzeroUniform,
    UnitaryCoherentPrepare,
    UnitaryCofactorFirstPrepare,
    UnitaryCharacterFirstGoodPhase,
    UnitaryCofactorFirstGoodPhase,
)
from group_action import UnitaryActionSelector, ValueSignedOOP, ValueGroupAction
from qarton.circuit import QuantumSimulator
from qarton.circuit.bit_vector import BitVector
from qarton.modular_arithmetic import ModInt
from resources import report

spec = importlib.util.spec_from_file_location(
    "independent_affine_reference", ROOT / "tests/affine_reference.py"
)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)


def amplitudes(n):
    a = [(1 + i % 3) * cmath.exp(0.213j * (i + 1)) for i in range(n)]
    for _ in range(4):
        z = math.sqrt(sum(abs(x) ** 2 for x in a))
        a = [x / z for x in a]
        if round(sum(abs(x) ** 2 for x in a), 15) == 1:
            break
    return a


def good_set(A, p, ell):
    return {
        x
        for x in range(1, p)
        if (P := ref.point_from_x(A, x, p, 1)) is not None
        and ref.multiply((p + 1) // ell, P, A, p, 1) is not None
    }


def preparation_success(p, ell):
    """Exact probability, ideal rotations, original degree-only phase angle.

    On the supersingular cyclic Montgomery domain, #good x=(p+1)(ell-1)/(2ell).
    For z=exp(i gamma), cos(gamma)=-1/(ell-1), the bad amplitude factor is
    1+(z-1)*(1+g*(z-1)); its squared magnitude is 4/(p-1)^2.
    """
    g = Fraction((p + 1) * (ell - 1), 2 * ell * (p - 1))
    return 1 - (1 - g) * Fraction(4, (p - 1) ** 2)


class Outcomes:
    def __init__(self, fail_at=None):
        self.fail_at = fail_at
        self.trace = []
        self.ambiguous = 0
        self.failure_selected = False
        self.probability = 1.0

    def choose(self, population, weights, *, k):
        assert k == 1
        pairs = [(int(v), float(w)) for v, w in zip(population, weights)]
        possible = [(v, w) for v, w in pairs if w > 1e-13]
        assert possible
        selected = None
        if len(possible) > 1:
            if self.fail_at == self.ambiguous:
                selected = 0
                self.failure_selected = True
            elif not self.failure_selected:
                selected = 1
            self.ambiguous += 1
        if selected is None:
            selected = max(possible, key=lambda x: x[1])[0]
        weight = dict(possible)[selected]
        probability = weight / sum(weights)
        self.probability *= probability
        self.trace.append(
            dict(outcomes=pairs, selected=selected, conditional_probability=probability)
        )
        return [selected]

    @contextmanager
    def active(self):
        original = random.choices
        random.choices = self.choose
        try:
            yield
        finally:
            random.choices = original


def run_native(qc, rows, amps, policy=None):
    initial = {qc.sig.input_to_bit_list(row): a for row, a in zip(rows, amps)}
    assert len(initial) == len(rows)
    sim = QuantumSimulator(qc, no_dummies=True)
    if policy is None:
        result = sim.simulate_internal(initial).quantum_state()
    else:
        with policy.active():
            result = sim.simulate_internal(initial).quantum_state()
    garbage = sum(
        abs(a) ** 2
        for k, a in result.items()
        if not k.is_zero(qc.output_signature.all_zeroes())
    )
    assert garbage < 1e-20, ("nonzero quantum workspace", garbage)
    assert abs(sum(abs(a) ** 2 for a in result.values()) - 1) < 1e-11
    return result, sim.count_ops, garbage


def compare_state(qc, actual, target_rows, target_amps):
    # Discarded classical scratch records need not be zero. Keep *every*
    # quantum wire, including all workspace, and every declared classical
    # output. Never erase quantum garbage or overlook an incorrect success bit.
    kept = set(qc.sig.all_qubits().to_tuple())
    for reg in qc.output_signature:
        kept.update(reg.positions)
    kept = sorted(kept)
    discarded = set(range(qc.nbr_bits())) - set(kept)
    assert all(len({k[i] for k in actual}) == 1 for i in discarded), (
        "discarded classical records are not classical on this trajectory"
    )

    def key(bits):
        return tuple(bits[i] for i in kept)

    wanted = {
        key(qc.output_signature.to_bit_list(row)): a
        for row, a in zip(target_rows, target_amps)
    }
    assert len(wanted) == len(target_rows)
    got = {key(k): a for k, a in actual.items()}
    assert len(got) == len(actual)
    error = sum(
        abs(got.get(k, 0) - wanted.get(k, 0)) ** 2 for k in set(got) | set(wanted)
    )
    if error >= 1e-20:
        diagnostic = dict(
            error=error,
            actual=[
                (str(qc.sig.bit_list_to_output(k)), str(a)) for k, a in actual.items()
            ],
            expected=[(str(row), str(a)) for row, a in zip(target_rows, target_amps)],
        )
        print("STATE_MISMATCH", json.dumps(diagnostic), flush=True)
    assert error < 1e-20, ("state or relative-phase mismatch", error)
    return error


def unitary_check(qc, rows, outputs, label, transform=None):
    started = time.monotonic()
    a = amplitudes(len(rows))
    expected_rows, expected_amps = (outputs, a) if transform is None else transform(a)
    actual, ops, garbage = run_native(qc, rows, a)
    error = compare_state(qc, actual, expected_rows, expected_amps)
    inverse = qc.inverse()
    norm = math.sqrt(sum(abs(v) ** 2 for v in actual.values()))
    initial = {BitVector(tuple(k)): v / norm for k, v in actual.items()}
    restored = (
        QuantumSimulator(inverse, no_dummies=True)
        .simulate_internal(initial)
        .quantum_state()
    )
    inverse_error = compare_state(inverse, restored, rows, a)
    return dict(
        test=label,
        input_branches=len(rows),
        output_branches=len(actual),
        squared_state_error=error,
        inverse_squared_state_error=inverse_error,
        nonzero_workspace_probability=garbage,
        basic_operations=ops,
        resources=report(qc),
        seconds=time.monotonic() - started,
    )


def uniform_case(p):
    qc = UnitaryNonzeroUniform(p, backends=UNITARY)

    def target(a):
        return list(range(1, p)), [a[0] / math.sqrt(p - 1)] * (p - 1)

    return unitary_check(qc, [0], None, f"uniform_p{p}", target)


def phase_case(p, order):
    primes = (3,) if p == 11 else (3, 5)
    curves = ref.valid_curves(p) if p == 11 else [0, ref.valid_curves(p)[1]]
    rows = [
        (ModInt(A, p), ModInt(x, p), ell, (p + 1) // ell)
        for A in curves
        for ell in primes
        for x in range(1, p)
    ]
    cls = (
        UnitaryCharacterFirstGoodPhase
        if order == "character"
        else UnitaryCofactorFirstGoodPhase
    )
    qc = cls(p, primes, 2 if p == 11 else 3, 3, backends=UNITARY)
    good = {(A, ell): good_set(A, p, ell) for A in curves for ell in primes}

    def target(a):
        return rows, [
            v
            * (
                cmath.exp(1j * math.acos(-1 / (ell - 1)))
                if x.v in good[A.v, ell]
                else 1
            )
            for (A, x, ell, k), v in zip(rows, a)
        ]

    return unitary_check(qc, rows, None, f"{order}_good_phase_p{p}", target)


def prepare_case(p, order):
    primes = (3,) if p == 11 else (3, 5)
    curves = ref.valid_curves(p) if p == 11 else [0]
    rows = [(ModInt(A, p), ell, (p + 1) // ell, 0) for A in curves for ell in primes]
    cls = (
        UnitaryCoherentPrepare if order == "character" else UnitaryCofactorFirstPrepare
    )
    qc = cls(p, primes, 2 if p == 11 else 3, 3, backends=UNITARY)

    def target(a):
        outputs, weights = [], []
        for row, amp in zip(rows, a):
            A, ell, k, _ = row
            good = good_set(A.v, p, ell)
            g = len(good) / (p - 1)
            z = cmath.exp(1j * math.acos(-1 / (ell - 1)))
            for x in range(1, p):
                factor = (z if x in good else 1) + (z - 1) * (1 + g * (z - 1))
                outputs.append((A, ell, k, x))
                weights.append(amp * factor / math.sqrt(p - 1))
        return outputs, weights

    result = unitary_check(qc, rows, None, f"{order}_preparation_p{p}", target)
    result["exact_preparation_success"] = {
        str(ell): str(preparation_success(p, ell)) for ell in primes
    }
    return result


def selector_case():
    p, primes, M = 59, (3, 5), 2
    rows, outputs = [], []
    for a in range(M + 1):
        for b in range(M + 1 - a):
            for sa in (False, True):
                for sb in (False, True):
                    for t in range(M):
                        row = (sa, a, sb, b, t, 0, 0, False, False)
                        if t < a:
                            ell, sign, en = 3, sa, True
                        elif t < a + b:
                            ell, sign, en = 5, sb, True
                        else:
                            ell, sign, en = 3, False, False
                        rows.append(row)
                        outputs.append(row[:5] + (ell, 60 // ell, sign, en))
    return unitary_check(
        UnitaryActionSelector(p, primes, M, backends=UNITARY),
        rows,
        outputs,
        "selector_p59_M2",
    )


def oop_case():
    p, primes = 11, (3,)
    rows, outputs = [], []
    for A in ref.valid_curves(p):
        for sign in (False, True):
            for enable in (False, True):
                B = ref.reference_step(A, 3, -1 if sign else 1, p) if enable else A
                for old in (0, B):
                    row = (ModInt(A, p), 3, 4, sign, enable, ModInt(old, p))
                    if row in rows:
                        continue
                    rows.append(row)
                    outputs.append(row[:-1] + (ModInt(0 if old == B else B, p), True))
    qc = ValueSignedOOP(p, primes, 2, 3, backends=UNITARY)
    policy = Outcomes()
    a = amplitudes(len(rows))
    started = time.monotonic()
    actual, ops, garbage = run_native(qc, rows, a, policy)
    error = compare_state(qc, actual, outputs, a)
    success = float(preparation_success(p, 3) ** 2)
    assert abs(policy.probability - success) < 1e-11
    return dict(
        test="signed_oop_compute_and_erase_p11",
        input_branches=len(rows),
        squared_state_error=error,
        nonzero_workspace_probability=garbage,
        success_probability=policy.probability,
        exact_success_probability=str(preparation_success(p, 3) ** 2),
        measurement_trace=policy.trace,
        basic_operations=ops,
        resources=report(qc),
        seconds=time.monotonic() - started,
    )


def group_case(p, M, fail_at=None):
    primes = (3,) if p == 11 else (3, 5)
    qc = ValueGroupAction(p, primes, M, 2 if p == 11 else 3, 3, backends=UNITARY)
    rows, exps = [], []
    if fail_at is not None:
        configs = [(0, (M,))]
    elif p == 11:
        configs = [
            (A, (e,))
            for A in (ref.valid_curves(p) if M == 1 else [0, 5])
            for e in range(-M, M + 1)
        ]
    else:
        # The p=11 suite covers zero and both signs exhaustively. Here use a
        # complex superposition across different degrees and opposite signs,
        # which also exposes degree-dependent heralding probabilities.
        configs = [(0, e) for e in ((1, 0), (0, -1))]
    for A, e in configs:
        flat = tuple(v for x in e for v in (x < 0, abs(x)))
        rows.append((ModInt(A, p), *flat, ModInt(0, p)))
        exps.append(e)
    a = amplitudes(len(rows))
    policy = Outcomes(fail_at)
    started = time.monotonic()
    actual, ops, garbage = run_native(qc, rows, a, policy)
    result = dict(
        test=f"group_p{p}_M{M}" + ("" if fail_at is None else f"_failure{fail_at}"),
        input_branches=len(rows),
        measurement_trace=policy.trace,
        path_probability=policy.probability,
        nonzero_workspace_probability=garbage,
        basic_operations=ops,
        resources=report(qc),
        seconds=time.monotonic() - started,
    )
    if fail_at is None:
        outputs = []
        factors = []
        for row, e in zip(rows, exps):
            B = ref.reference_action(row[0].v, e, primes, p)
            outputs.append(row[:-1] + (ModInt(B, p), True))
            degrees = [l for l, x in zip(primes, e) for _ in range(abs(x))]
            degrees += [primes[0]] * (M - len(degrees))
            factors.append(
                math.prod(float(preparation_success(p, l) ** 2) for l in degrees)
            )
        norm = math.sqrt(sum(abs(v * f) ** 2 for v, f in zip(a, factors)))
        expected = [v * f / norm for v, f in zip(a, factors)]
        result["squared_state_error"] = compare_state(qc, actual, outputs, expected)
        assert abs(policy.probability - norm * norm) < 1e-11
        overlap = sum(v.conjugate() * w for v, w in zip(a, expected))
        result["ideal_action_conditional_trace_distance"] = math.sqrt(
            max(0, 1 - abs(overlap) ** 2)
        )
        result["predicted_success_probability"] = norm * norm
        assert len(policy.trace) == 1 + M * (12 + 3 * p.bit_length()), len(policy.trace)
    else:
        assert policy.failure_selected
        for k in actual:
            decoded = qc.sig.bit_list_to_output(k)
            assert not decoded[-1], ("failure misreported as success", decoded)
            assert decoded[:-2] == rows[0][:-1], (
                "input exponent or curve corrupted",
                decoded,
            )
            if fail_at % 6 < 3:
                completed = fail_at // 6
                expected_A = ref.reference_action(rows[0][0].v, (completed,), primes, p)
                assert decoded[-2] == ModInt(expected_A, p), (
                    "failed OOP changed the last completed action output",
                    decoded,
                )
        assert len(policy.trace) < 1 + M * (12 + 3 * p.bit_length()), (
            "later action slot was not skipped"
        )
        result["later_action_slots_skipped"] = True
        result["failure_selected_at_ambiguous_measurement"] = fail_at
    return result


CASES = {
    "uniform11": lambda: uniform_case(11),
    "uniform59": lambda: uniform_case(59),
    "phase11_character": lambda: phase_case(11, "character"),
    "phase11_cofactor": lambda: phase_case(11, "cofactor"),
    "prepare11_character": lambda: prepare_case(11, "character"),
    "prepare11_cofactor": lambda: prepare_case(11, "cofactor"),
    "prepare59_character": lambda: prepare_case(59, "character"),
    "prepare59_cofactor": lambda: prepare_case(59, "cofactor"),
    "selector": selector_case,
    "signed_oop": oop_case,
    "group11_M1": lambda: group_case(11, 1),
    "group11_M2": lambda: group_case(11, 2),
    "group59_M1": lambda: group_case(59, 1),
    **{f"failure{i}": (lambda i=i: group_case(11, 2, i)) for i in range(6)},
    "failure_mid_chain": lambda: group_case(11, 3, 6),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=list(CASES), required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = args.output or ROOT / "results/complete_validation" / f"{args.case}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    print(json.dumps({"starting": args.case}), flush=True)
    result = CASES[args.case]()
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps({k: result[k] for k in ("test", "input_branches", "seconds")}),
        flush=True,
    )


if __name__ == "__main__":
    main()
