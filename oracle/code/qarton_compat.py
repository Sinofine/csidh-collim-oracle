"""Minimal, explicit compatibility repair for the pinned Qarton revision.

Qarton 8c305f6 compares BitVectors by hash equality. Python integer hashes
collide, for example at 1 and 2**sys.hash_info.modulus.bit_length() on this
runtime. Distinct large-register basis states must compare by bits and length.
No gate evolution, arithmetic circuit or resource counter is replaced.
"""

from qarton.circuit.bit_vector import BitVector


def exact_bitvector_equal(self, other):
    return (
        isinstance(other, BitVector)
        and len(self) == len(other)
        and self.to_int() == other.to_int()
    )


def install_bitvector_equality_fix():
    BitVector.__eq__ = exact_bitvector_equal


def collision_regression():
    import sys

    shift = sys.hash_info.modulus.bit_length()
    a = BitVector.from_int(1, shift + 1)
    b = BitVector.from_int(1 << shift, shift + 1)
    assert hash(a) == hash(b), "Expected deterministic Python integer hash collision"
    assert a.to_int() != b.to_int()
    assert a != b and len({a, b}) == 2
    assert a == a.copy() and a != BitVector.from_int(1, shift + 2)
    return dict(
        shift=shift,
        distinct_states_preserved=True,
        collision_pair=[str(a.to_int()), str(b.to_int())],
    )
