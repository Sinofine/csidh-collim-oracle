"""Native Qarton coherent component tests; every classical dummy is disabled."""

import cmath
import importlib.util
import json
import math
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from arithmetic import install_before_component_imports

install_before_component_imports()
from bonnetain import BonnetainPointStep
from character import UnitaryNegativeCharacterPhase, UnitaryPositiveCharacterFlag
from conversion import UnitaryValueKernelAdapter
from points import UNITARY, ExplicitUnitaryWPointStep, MaskedProductIntoZero
from qarton.circuit import QuantumSimulator
from qarton.circuit.bit_vector import BitVector
from qarton.modular_arithmetic import ModInt
from resources import report


def check(qc, inputs, outputs, label, phases=None):
    resources = report(qc)
    assert (
        resources["classical_bits"] == 0 and resources["primitives"].get("msr", 0) == 0
    ), resources
    raw = [(1 + i % 3) * cmath.exp(0.213j * (i + 1)) for i in range(len(inputs))]
    norm = math.sqrt(sum(abs(a) ** 2 for a in raw))
    amps = [a / norm for a in raw]
    initial = {qc.sig.input_to_bit_list(row): a for row, a in zip(inputs, amps)}
    target = {
        tuple(qc.output_signature.to_bit_list(row)): a
        * (1 if phases is None else phases[i])
        for i, (row, a) in enumerate(zip(outputs, amps))
    }
    assert len(initial) == len(inputs), "duplicate test inputs"
    started = time.monotonic()
    sim = QuantumSimulator(qc, no_dummies=True)
    actual = sim.simulate_internal(initial).quantum_state()
    # Compare immutable full bit tuples, including every ancilla. BitVector is
    # mutable, so never use a simulator-owned key as the comparison identity.
    actual_bits = {tuple(k): a for k, a in actual.items()}
    error = sum(
        abs(actual_bits.get(k, 0) - target.get(k, 0)) ** 2
        for k in set(actual_bits) | set(target)
    )
    garbage = sum(
        abs(a) ** 2
        for k, a in actual.items()
        if not k.is_zero(qc.output_signature.all_zeroes())
    )
    if error >= 1e-20:
        print(
            "Mismatch",
            label,
            "actual",
            [(qc.sig.bit_list_to_output(k), a) for k, a in list(actual.items())[:3]],
            "expected",
            list(zip(outputs, amps))[:3],
            flush=True,
        )
    assert error < 1e-20 and garbage < 1e-20, (label, error, garbage)
    # Apply the actual inverse to the actual complete output, including work.
    measured_norm = sum(abs(a) ** 2 for a in actual.values())
    assert abs(measured_norm - 1) < 1e-12
    # quantum_state() rounds amplitudes to 15 decimals; normalize that returned
    # representation before the simulator's strict input normalization check.
    inverse = qc.inverse()
    assert inverse.input_signature == qc.output_signature
    renormalized = {
        BitVector(tuple(k)): a / math.sqrt(measured_norm) for k, a in actual.items()
    }
    restored = (
        QuantumSimulator(inverse, no_dummies=True)
        .simulate_internal(renormalized)
        .quantum_state()
    )
    restored = {tuple(k): a for k, a in restored.items()}
    original_in_inverse_layout = {
        tuple(inverse.output_signature.to_bit_list(row)): a
        for row, a in zip(inputs, amps)
    }
    inv_error = sum(
        abs(restored.get(k, 0) - original_in_inverse_layout.get(k, 0)) ** 2
        for k in set(restored) | set(original_in_inverse_layout)
    )
    if inv_error >= 1e-20:
        for k, a in list(restored.items())[:2]:
            wanted = next(iter(original_in_inverse_layout))
            print(
                "inverse mismatch",
                label,
                "actual",
                inverse.sig.bit_list_to_output(BitVector(k)),
                a,
                "expected input first",
                inputs[0],
                amps[0],
                "diff first",
                [i for i, (u, v) in enumerate(zip(k, wanted)) if u != v],
                flush=True,
            )
    assert inv_error < 1e-20, (label, inv_error)
    row = dict(
        test=label,
        branches=len(inputs),
        squared_state_error=error,
        nonzero_workspace_probability=garbage,
        inverse_squared_error=inv_error,
        basic_operations=sim.count_ops,
        output_norm_squared=measured_norm,
        resources=report(qc),
        seconds=time.monotonic() - started,
    )
    print(json.dumps({k: row[k] for k in ["test", "branches", "seconds"]}), flush=True)
    return row


def main():
    rng = random.Random(20261015)
    results = []
    for p, w, s in [(11, 2, 2), (59, 3, 3)]:
        m = lambda x: ModInt(x % p, p)
        xs = list(range(p)) if p == 11 else [0, 1, 2, 3, 5, 17, p - 1]
        rows = [(m(x), bool(e), bool(t)) for x in xs for e in (0, 1) for t in (0, 1)]
        expected = [
            (x, e, t ^ (e and pow(x.v, (p - 1) // 2, p) == 1)) for x, e, t in rows
        ]
        results.append(
            check(
                UnitaryPositiveCharacterFlag(p, w, s, backends=UNITARY),
                rows,
                expected,
                f"flag_p{p}",
            )
        )
        rows = [m(x) for x in xs]
        phases = [-1 if pow(x, (p - 1) // 2, p) == p - 1 else 1 for x in xs]
        results.append(
            check(
                UnitaryNegativeCharacterPhase(p, w, s, backends=UNITARY),
                rows,
                rows,
                f"phase_p{p}",
                phases,
            )
        )
        pairs = (
            [(a, b) for a in range(p) for b in range(p)]
            if p == 11
            else [(0, 0), (0, 1), (1, 0)]
            + [(rng.randrange(p), rng.randrange(p)) for _ in range(12)]
        )
        rows = [(bool(f), m(a), m(b), m(0)) for a, b in pairs for f in (0, 1)]
        expected = [(f, a, b, m(a.v * b.v if f else 0)) for f, a, b, _ in rows]
        results.append(
            check(
                MaskedProductIntoZero(p, backends=UNITARY),
                rows,
                expected,
                f"masked_product_p{p}",
            )
        )
        for design, cls in [
            ("w", ExplicitUnitaryWPointStep),
            ("montgomery", BonnetainPointStep),
        ]:
            rows, expected = [], []
            for j in range(12):
                coef, qw, qz, pw, pz, x, z = [rng.randrange(p) for _ in range(7)]
                ell = 3 if j % 2 else 5
                active = ell > 3
                if design == "w":
                    nx = 4 * pz * (x * qz - qw * z) ** 2 % p
                    nz = 4 * pw * (x * qw - z * qz) ** 2 % p
                else:
                    a = (qw + qz) * (x - z)
                    b = (qw - qz) * (x + z)
                    nx = pz * (a + b) ** 2 % p
                    nz = pw * (a - b) ** 2 % p
                row = (
                    m(coef),
                    m(qw),
                    m(qz),
                    ell,
                    m(pw),
                    m(pz),
                    m(x),
                    m(z),
                    m(0),
                    m(0),
                    m(0),
                    m(0),
                )
                result = row[:8] + (
                    (m(x), m(z), m(nx), m(nz)) if active else (m(pw), m(pz), m(x), m(z))
                )
                rows.append(row)
                expected.append(result)
            results.append(
                check(
                    cls(p, 3, 1, backends=UNITARY),
                    rows,
                    expected,
                    f"{design}_point_p{p}",
                )
            )
    spec = importlib.util.spec_from_file_location(
        "independent_affine_reference", ROOT / "tests/affine_reference.py"
    )
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    for p, primes in [(11, (3,)), (59, (3, 5)), (419, (7,))]:
        rows, expected = [], []
        for ell in primes:
            for direction in (-1, 1):
                P = ref.select_good_point(0, p, direction)
                Q = ref.multiply((p + 1) // ell, P, 0, p, direction)
                B = ref.reference_step(0, ell, direction, p)
                for enable in (True,) if p == 419 else (False, True):
                    for target in (3,) if p == 419 else (0, 3):
                        row = (
                            ModInt(0, p),
                            ModInt(Q[0], p),
                            ModInt(1, p),
                            ell,
                            enable,
                            target,
                        )
                        rows.append(row)
                        expected.append(row[:-1] + (target ^ (B if enable else 0),))
        results.append(
            check(
                UnitaryValueKernelAdapter(p, primes, 4, backends=UNITARY),
                rows,
                expected,
                f"kernel_adapter_p{p}",
            )
        )
    output = ROOT / "results/component-tests.json"
    output.write_text(json.dumps({"tests": results}, indent=2) + "\n")


if __name__ == "__main__":
    main()
