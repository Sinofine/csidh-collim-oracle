# Collimation Sieve

This Haskell implementation is the primary simulator accompanying the
artifact. It simulates and benchmarks a generalization of
Kuperberg's quantum *collimation sieve* algorithm for arbitrary cyclic
groups. We have used it to estimate the quantum (in)security of CSIDH
("Commutative SIDH"), a proposed "post quantum" drop-in replacement
for noninteractive Diffie-Hellman-style key agreement and
encryption. See
[our paper](https://web.eecs.umich.edu/~cpeikert/pubs/csidh-sieve.pdf)
for further details.

## Building

1. Make sure you have a recent version of
[Haskell stack installed](https://docs.haskellstack.org/en/stable/install_and_upgrade/). (Warning:
pre-1.0 versions of stack will crash with a parsing error.)

2. Run `stack build` from this directory. The pinned Stack snapshot may need
   to download its compiler and dependencies on the first build.

3. Run the test suite with `stack test`.

## Execution

* Run the sieve with

        stack run <log N> <log L> <log S> [threshold] [arity] [mapping-mode]
  where `N` is the group order (use `log N=0` for the exact CSIDH-512
  group order), `L` is the desired phase vector length, `S` is the
  desired range size, and `threshold < 1` (optional, defaults to 0.25)
  is the factor that determines whether a phase vector is too short
  (i.e., any vector shorter than `threshold * L` is discarded).
  `arity` is optional and may be `2` (the default binary sieve) or `4`
  (the four-list streaming interval-merge sieve).
  For arity 4, `mapping-mode` may be `--build-map` or `--print-inv`.

* Small parameters like `log N=150`, `log L=16`, `log S=16` will cause
  the sieve to finish relatively quickly (and use relatively little
  memory), and will give you an idea of how it works and how to
  interpret the output.

### Four-list collimation

The arity-4 implementation generates the sums of the first and last
pairs of phase vectors as ascending and descending streams. It uses
heap-based multi-way merging and a retained sliding interval, so the
two lists of all pairwise sums are never materialized. For four
length-`L` inputs, generating the pair streams takes `O(L^2 log L)` time.
The implementations sort the retained output, giving total time
`O(L^2 log L + |J| log |J|)`. The two heaps use `O(L)` memory; the retained
matching window has expected `O(L)` size under the uniform-range sieve
schedule but can contain `O(L^2)` entries in the worst case. The current
implementation also groups equal left-pair sums, whose maximum group size is
an additional worst-case working-memory term.

The library exposes two variants. `collimate4At` and `collimate4` retain only
the output phase multipliers for lower-memory simulation. The corresponding
`collimate4AtWithMap` and `collimate4WithMap` functions additionally return
`inverseTable[J] = (j1,j2,j3,j4)` and the reverse `permutationTable`, thereby
materializing the classical `inv` and `pi` tables required to reindex the
surviving tuples. These tables use `O(|J|)` additional memory.

To materialize the complete `inv` and `pi` tables at the root without printing
their entries, use `--build-map`, for example:

        stack run -- 64 10 10 0.25 4 --build-map

This is the appropriate mode for timing the root mapping construction. The
mapping is fully evaluated inside the reported CPU time, rather than being
left lazy until the reporting phase. Use external process timing when a
wall-clock measurement is required.

To additionally print the complete `inv` table, use `--print-inv`, for example:

        stack run -- 20 3 6 0.25 4 --print-inv

The output contains one line `inv[J] = (j1,j2,j3,j4)` for every surviving
root tuple. Printing occurs after the reported sieve CPU time, so terminal I/O is
not counted as mapping construction; its duration is reported separately as
`inv printing time`. Internal nodes continue to use the phase-only path in
both modes.

The arity-4 interval schedule uses the calibrated constant
`kappa_4 = 151/315` for the sum of four uniform random variables. This
is the four-list analogue of the binary implementation's `2/3`
constant. Thus `S_(i+1)` is approximately
`(151/315) * S_i * L^3`, and the product of the four requested child
lengths is approximately `(315/151) * L * S_(i+1) / S_i`. This keeps
the expected output length close to the requested length. For example:

        stack run 150 16 16 0.25 4

* `stack run` will use as many CPU cores as it deems appropriate,
  which may not give the best performance. To specify the number of
  cores and other options, use `stack exec` with RTS
  options. (Caution: `stack run` silently ignores RTS options.) For
  example,

        stack exec -- collimation-sieve-exe <args> +RTS -N4 -s -RTS
  runs with four runtime capabilities and outputs memory and
  garbage-collection statistics. On workloads with no useful parallel sparks,
  `+RTS -N1 -s -RTS` can be faster and gives a cleaner reproducibility
  baseline.

## Measurement scope

This program is a classical phase-vector simulator. Its query counters model
calls to the hidden-shift oracle, and its `inv` and `pi` objects are classical
tables. They are not measurements of a coherent QROM/QRACM circuit, T gates,
or quantum workspace. In `--build-map` and `--print-inv` modes the complete
mapping is materialized at the root; internal recursive nodes retain the
lower-memory phase-vector path.

* **WARNING:** a length limit in Haskell's `vector` package
  effectively imposes a length limit of 2^30 on the sieve's vectors,
  so the sieve will not run reliably for `log L >= 26` unless
  `threshold` is increased somewhat, e.g., to `0.4` (and even this is
  not a guarantee). We are working on a tweak to address this, but in
  any case note that `log L >= 26` can use 100 GB of RAM or more.
