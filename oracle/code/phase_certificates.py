"""Independent Arb ball-arithmetic verifier of chronological Clifford+T words."""

from flint import arb, acb, ctx
from exact_angles import Angle


INVERSE = {"h": "h", "x": "x", "z": "z", "s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t"}


def inverse_word(tokens):
    return [INVERSE[t] for t in reversed(tokens)]


def positive_angle(angle):
    return (-angle, -1) if angle.sn < 0 else (angle, 1)


def angle_ball(angle):
    if angle.kind == "pi":
        value = arb.pi()
    elif angle.kind == "acos":
        value = (arb(angle.n) / angle.d).acos()
    else:
        value = (arb(angle.n) / angle.d).sqrt().acos()
    return value * angle.sn / angle.sd


def matmul(a, b):
    return [
        [
            sum((a[i][k] * b[k][j] for k in range(len(b))), acb(0))
            for j in range(len(b[0]))
        ]
        for i in range(len(a))
    ]


def matrix_word(tokens):
    rt = arb(2).sqrt()
    w = acb(0, arb.pi() / 4).exp()
    I = [[acb(1), acb(0)], [acb(0), acb(1)]]
    gates = {
        "h": [[acb(1) / rt, acb(1) / rt], [acb(1) / rt, -acb(1) / rt]],
        "x": [[acb(0), acb(1)], [acb(1), acb(0)]],
        "z": [[acb(1), acb(0)], [acb(0), acb(-1)]],
        "s": [[acb(1), acb(0)], [acb(0), acb(0, 1)]],
        "sdg": [[acb(1), acb(0)], [acb(0), acb(0, -1)]],
        "t": [[acb(1), acb(0)], [acb(0), w]],
        "tdg": [[acb(1), acb(0)], [acb(0), w.conjugate()]],
    }
    result = I
    for token in tokens:
        result = matmul(gates[token], result)
    return result


def frobenius_upper(difference):
    # Sum upper magnitudes to avoid dependence on a numerical eigensolver.
    total = arb(0)
    for row in difference:
        for value in row:
            total += value.abs_upper() ** 2
    return total.sqrt().upper()


def phase_error(angle, entry, precision=640):
    """Certify ||exp(i beta) U_word - P(theta)||_2 <= Frobenius bound."""
    with ctx.workprec(precision):
        theta = angle_ball(angle)
        beta = arb(0)
        if entry["phase_frame"]["angle_half"]:
            beta += theta / 2
        beta += arb.pi() * entry["phase_frame"]["pi_quarters"] / 4
        factor = acb(0, beta).exp()
        u = matrix_word(entry["tokens"])
        wanted = [[acb(1), acb(0)], [acb(0), acb(0, theta).exp()]]
        return frobenius_upper(
            [[factor * u[i][j] - wanted[i][j] for j in range(2)] for i in range(2)]
        )


def verify_rotation(key, entry, bits):
    with ctx.workprec(640):
        bound = phase_error(Angle.from_key(key), entry)
        assert bound < arb(2) ** (-bits), (key, bound, bits)
        assert entry["t_count"] == sum(t in ("t", "tdg") for t in entry["tokens"])
        return bound
