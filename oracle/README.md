# CSIDH circuits in Qarton

Reversible circuits and resource estimates for the CSIDH-512 group-action oracle. Two constructions are provided: the Bonnetain–Schrottenloher algorithm and the revised algorithm with fixed-window character computation and Edwards-coordinate isogeny evaluation.

Both use the same arithmetic backend and outer action circuit, including the same masked zero-target products and terminal character-bit tests. The implementation does not yet reproduce the papers' resource bounds; the differences are listed below.

## Requirements

Python 3.14 or later and the Qarton revision pinned in `requirements.txt`:

```bash
python3.14 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The scripts run directly; the project does not need to be installed.

## Resource estimates

Run from the project directory:

```bash
python scripts/construct_oracle.py --design bonnetain
python scripts/construct_oracle.py --design revised
python scripts/compare_designs.py
```

Construction counts circuit resources without simulating a quantum state. Results are written to `results/bonnetain/`, `results/reference/`, and `results/comparison/`. The parameters are `n=511`, `s=15`, `M=1300`, and window width 5.

| Construction | Total qubits | CCX gates |
|---|---:|---:|
| Bonnetain | 44,351 | 3,159,759,945,292,600 |
| Revised | 42,814 | 3,055,649,556,172,200 |

Qubit totals include the 1,910-bit interface. Six classical bits are reported separately. The JSON files also report `7 * CCX`; this conversion excludes synthesis of the remaining phase operations.

## Code

- `code/character.py`, `points.py`, `isogeny.py`: revised character and local isogeny circuits.
- `code/bonnetain.py`: Bonnetain character and local isogeny circuits.
- `code/ladder.py`, `preparation.py`, `conversion.py`, `group_action.py`: cofactor multiplication, point preparation, coordinate conversion, and action composition.
- `code/arithmetic.py`, `expressions.py`, `schedules.py`: arithmetic selection, field expressions, and reversible schedules.
- `code/parameters.py`, `resources.py`, `utils.py`, `reference.py`: parameters, counting, controls, and classical ladder callbacks.
- `scripts/`: construction, comparison, and verification commands.
- `results/`: recorded resource counts and formula checks.
- `tests/`: native coherent component tests and an independent affine reference.

The Bonnetain construction follows ePrint 2018/537, Sections 5.1–5.3 and Appendix B. The revised character and local isogeny circuits follow Sections 3.1 and 3.2 of the supplied revised manuscript, respectively.

## Checks

```bash
python scripts/check_resources.py --result results/reference/oracle.json
python scripts/check_resources.py --result results/bonnetain/oracle.json
python scripts/check_design_formulas.py
python tests/test_components.py
```

The first two commands check resource accounting. The formula check compares the two isogeny expressions on small classical curves; it does not verify the complete quantum circuit.

The component tests disable Qarton's classical callbacks and check complex amplitudes, workspace cleanup, and inverse circuits. They cover character flags/phases, masked products, differential point steps, and the kernel-to-curve adapter at small primes.

For direct child-circuit counts:

```bash
python scripts/construct_oracle.py --design bonnetain --gate-audit --output results/gate-audit/bonnetain.json
python scripts/construct_oracle.py --design revised --gate-audit --output results/gate-audit/revised.json
python scripts/audit_gate_counts.py
```

The resulting JSON files contain the direct child counts for each recorded circuit component.

## Differences from the paper estimates

The revised construction uses 1,537 fewer qubits and 3.29488% fewer CCX gates than the Bonnetain construction with the common arithmetic changes. The paper's arithmetic model predicts a 2.78288% gate reduction. These are different accounting models, and the current implementation has unresolved schedule differences:

- The character circuit saves 60.157% CCX. It makes 5772 versus 2300 field-primitive calls, with a single terminal CCX counted separately. The paper charges 5775 versus 2301, including terminal marking.
- Both differential point steps make 10 multiplication/square calls including cleanup and masked outputs. Neither currently matches the paper's `7 TM` step estimate. The Edwards step shares two products between its output coordinates; their common factor of four does not change the projective point. Per out-of-place action, character computation saves 39,934,240,416 CCX and local mapping saves 108,216,938 CCX. These savings are charged over 2600 calls.
- The ladder implements the listed 16 multiplication/square operations; the lemma states 15. The cofactor scalar uses 511 bits while the reference chain has 512 steps.
- The shared local traversal does not yet reproduce the original `B(m+1,s)` versus revised `B(m,s)` accounting. Register lifetimes and the cited arithmetic workspace bounds also remain to be reconciled.
- Preparation returns a classical success flag. Failure handling, successful-branch amplitudes, exceptional points, zero denominators, and complete cleanup require full-circuit verification. Character phase and controlled Boolean-flag costs are not interchangeable.

The current results establish construction and resource counts. Small-field formula checks cover 72 valid odd-prime-order kernels over fields of size 11, 19, and 59. The native component suite covers 13 tests and 436 coherent branches. Full-oracle quantum functionality and synthesis of the remaining phase gates have not been verified.
