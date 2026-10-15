# CSIDH quantum algorithms

Code for CSIDH group-action circuits and collimation-sieve experiments.

| Directory | Purpose |
|---|---|
| [oracle/](oracle/) | Qarton constructions and resource counts for the Bonnetain and revised CSIDH-512 group-action oracles |
| [sieve/haskell/](sieve/haskell/) | Peikert's collimation-sieve simulator with binary and four-list collimation |
| [sieve/rust/](sieve/rust/) | Rust implementation with fixed random seeds and binary/four-list collimation |

The three programs run independently. The sieve programs simulate classical phase-vector processing and count modeled oracle queries. They do not call the Qarton circuit or measure its quantum gate costs.

## Group-action oracle

Requires Python 3.14 or later. From `oracle/`:

```sh
python3.14 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python scripts/construct_oracle.py --design bonnetain
python scripts/construct_oracle.py --design revised
python scripts/compare_designs.py
```

Recorded results are in `oracle/results/`. The current constructions use 44,351 and 42,814 qubits respectively. They do not yet reproduce the papers' resource estimates; see the [oracle README](oracle/README.md) for the accounting differences.

## Haskell sieve

Requires Stack. From `sieve/haskell/`:

```sh
stack build
stack test
stack run -- 64 10 10 0.25 4 --build-map
```

The arguments are `log N`, `log L`, `log S`, discard threshold, and arity. Arity is 2 or 4. Use `log N = 0` for the exact CSIDH-512 group order. `--build-map` constructs the root reindexing tables. See the [Haskell README](sieve/haskell/README.md).

## Rust sieve

Requires Rust and Cargo. From `sieve/rust/`:

```sh
cargo test --release
cargo build --release
./target/release/collimation-sieve-rs 80 10 10 0.25 4 --build-map --seed 1
```

The positional arguments match the Haskell program. The Rust program also accepts a random seed. See the [Rust README](sieve/rust/README.md) for experiment options.

## Sources and licenses

The oracle code follows Bonnetain–Schrottenloher, ePrint 2018/537, and the supplied revised CSIDH manuscript. The Haskell sieve extends Chris Peikert's reference implementation; the Rust sieve is a separate implementation. Both sieve directories retain their BSD-3-Clause licenses and copyright notices. Those licenses do not assign a license to the oracle code. Qarton is an external dependency pinned in `oracle/requirements.txt`.
