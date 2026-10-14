# Supplementary Rust implementation

This is a dependency-free Rust reimplementation of Peikert's generalized
collimation sieve. The Haskell implementation in this artifact is the primary
simulator used by the manuscript; this Rust implementation provides fixed-seed
reproducibility and an independent cross-check of the binary and four-ary
paths. 

Implemented features:

- binary and 4-ary recursive collimation;
- adaptive child lengths and discard/retry from Section 3.2.2;
- high-bit quotient measurement from Algorithm 2;
- streaming 4-list matching using two pair-sum heaps;
- explicit `inv[J]` and reverse `pi` at the root;
- exact CSIDH-512 order (`log N = 0`) and power-of-two test groups;
- actual/model oracle-query counts, maximum vector length, discard rate,
  full regularization, punctured states, and numerical inverse-QFT success
  estimates for the 1, 2, 4, and 8 closest outcomes;
- a fixed seed for repeatable experiments.

The `--build-map` and `--print-inv` modes materialize `inv` and `pi` at the
root. Internal recursive nodes use the phase-vector path. 

## Build and test

```sh
cargo test --release
cargo build --release
```

## Examples

```sh
# Refined 4-ary experiment, N=2^80 and L=S=2^10
./target/release/collimation-sieve-rs 80 10 10 0.25 4 --seed 1

# Construct inv and pi, but do not charge terminal output to construction time
./target/release/collimation-sieve-rs 80 10 10 0.25 4 --build-map --seed 1

# Construct and print every inv[J]; printing time is reported separately
./target/release/collimation-sieve-rs 80 10 10 0.25 4 --print-inv --seed 1

# Print a progress line after every 5000 entered sieve nodes
./target/release/collimation-sieve-rs 0 10 10 0.25 4 --progress 5000 --seed 1
```
