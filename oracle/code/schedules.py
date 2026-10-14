"""Fixed-window schedule and exact Bennett dynamic programming."""

from functools import lru_cache
from math import inf

from parameters import P512


def window(exponent, width=5):
    if exponent <= 0 or width < 1:
        raise ValueError("positive exponent and width required")
    bits = bin(exponent)[2:]
    i = 0
    windows = []
    gap = 0
    while i < len(bits):
        if bits[i] == "0":
            gap += 1
            i += 1
            continue
        j = min(i + width, len(bits))
        while bits[j - 1] == "0":
            j -= 1
        a = int(bits[i:j], 2)
        gap += j - i
        windows.append((gap, a))
        gap = 0
        i = j
    if gap:
        windows.append((gap, 0))
    a0 = windows[0][1]
    chain = []
    value = a0
    for g, a in windows[1:]:
        chain.extend([("square", 0)] * g)
        if a:
            chain.append(("multiply", a))
        value = (value << g) + a
    assert value == exponent
    return a0, tuple(chain), tuple(windows)


class Bennett:
    def __init__(self, t, s):
        self.b = [[inf] * (s + 1) for _ in range(t + 1)]
        self.k = [[0] * (s + 1) for _ in range(t + 1)]
        for j in range(s + 1):
            self.b[1][j] = 1
        for j in range(1, s + 1):
            for n in range(2, t + 1):
                cost, k = min(
                    (self.b[k][j] + self.b[k][j - 1] + self.b[n - k][j - 1], k)
                    for k in range(1, n)
                )
                self.b[n][j] = cost
                self.k[n][j] = k

    def bridge(
        self, start, length, s, source, target, allocate, release, apply, inverse=False
    ):
        """Clean evaluation with target retained, including a different inverse schedule.
        A call's output is external to its s temporary checkpoints.
        """
        if length == 1:
            apply(start, source, target, inverse)
            return
        if s < 1:
            raise ValueError("insufficient checkpoint budget")
        k = self.k[length][s]
        tmp = allocate()
        if not inverse:
            self.bridge(start, k, s, source, tmp, allocate, release, apply)
            self.bridge(
                start + k, length - k, s - 1, tmp, target, allocate, release, apply
            )
            self.bridge(start, k, s - 1, source, tmp, allocate, release, apply, True)
        else:
            self.bridge(start, k, s - 1, source, tmp, allocate, release, apply)
            self.bridge(
                start + k,
                length - k,
                s - 1,
                tmp,
                target,
                allocate,
                release,
                apply,
                True,
            )
            self.bridge(start, k, s, source, tmp, allocate, release, apply, True)
        release(tmp)


@lru_cache(None)
def bennett(t, s):
    return Bennett(t, s)


def certificate(p=P512, width=5):
    a0, chain, windows = window((p - 1) // 2, width)
    sq = 1 + sum(op == "square" for op, a in chain)
    mul = (2 ** (width - 1) - 1) + sum(op == "multiply" for op, a in chain)
    return dict(
        p=p,
        field_bits=p.bit_length(),
        exponent=(p - 1) // 2,
        width=width,
        a0=a0,
        windows=windows,
        chain_length=len(chain),
        forward_squarings=sq,
        forward_multiplications=mul,
        full_storage_field_registers=2 + (2 ** (width - 1)) + 1 + len(chain),
    )
