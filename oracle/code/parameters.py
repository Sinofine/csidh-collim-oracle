"""Public CSIDH-512 field and degree parameters."""

from math import prod

PRIMES = tuple(
    i for i in range(3, 374, 2) if all(i % d for d in range(2, int(i**0.5) + 1))
) + (587,)
P512 = 4 * prod(PRIMES) - 1
