#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root_dir"
cargo build --release

seeds="${SEEDS:-$(seq 1 50)}"
configs=(
  "full 10 4"
  "full 11 5"
  "full 12 6"
  "punctured 6 6"
  "punctured 8 8"
  "punctured 10 10"
)

for spec in "${configs[@]}"; do
  read -r experiment log_l log_s <<< "$spec"
  for seed in $seeds; do
    echo "===== $experiment N=2^80 logL=$log_l logS=$log_s r=4 seed=$seed ====="
    target/release/collimation-sieve-rs \
      80 "$log_l" "$log_s" 0.25 4 \
      --build-map \
      --leaf-rounding nearest \
      --keep-short-root \
      --seed "$seed"
  done
done
