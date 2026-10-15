"""Install actual Qarton Clifford+T gadgets from the verified phase catalogue.

Run each catalogue/policy in a fresh process: Qarton memoizes constructions.
Global phases are discarded only for the present, classically controlled
instrument. This context must not be used to manufacture a quantum-controlled
version of an already synthesized circuit.
"""

from contextlib import contextmanager
import json
from pathlib import Path
from qarton.circuit import (
    BoolType,
    Circuit,
    InPlaceCircuit,
    qc_bool_reg,
    qc_reg_from_bits,
)
from qarton.circuit.basic_gates import Phase, CPhase, GPhase
from exact_angles import Angle
from phase_certificates import inverse_word, positive_angle

_data = None
_policy = None
_rows = None


def tokens_for(key):
    angle, sign = positive_angle(Angle.from_key(key))
    word = _data["rotations"][str(angle)]["tokens"]
    return word if sign > 0 else inverse_word(word)


def emit(parent, bit, tokens):
    for token in tokens:
        assert token in ("h", "x", "z", "t", "tdg", "s", "sdg")
        getattr(parent, token)(bit)


class LoweredPhase(InPlaceCircuit):
    def __init__(self, key):
        super().__init__()
        reg = self.add_input_output(BoolType())
        emit(self, reg[0], tokens_for(key))


class LoweredCPhase(InPlaceCircuit):
    def __init__(self, key):
        super().__init__()
        reg = self.add_input_output(2)
        row = _rows[("CPhase", key)]
        if row["method"] == "and_phase_uncompute":
            anc = self.get_anc(BoolType())
            self.ccx(reg[0], reg[1], anc[0])
            self.append(LoweredPhase(key), anc)
            self.ccx(reg[0], reg[1], anc[0])
            self.release_anc(anc)
        else:
            angle = Angle.from_key(key) / 2
            # P_c(a) P_t(a) CX P_t(-a) CX = CP(2a), up to the
            # recorded, state-independent phases of the synthesized P gates.
            self.append(LoweredPhase(str(angle)), qc_bool_reg(reg[0]))
            self.append(LoweredPhase(str(angle)), qc_bool_reg(reg[1]))
            self.cx(reg[0], reg[1])
            self.append(LoweredPhase(str(-angle)), qc_bool_reg(reg[1]))
            self.cx(reg[0], reg[1])


@contextmanager
def lower_phases(catalogue, policy="ancilla_free"):
    global _data, _policy, _rows
    assert _data is None, "Use a fresh process for each compiled construction."
    _data = json.loads(Path(catalogue).read_text())
    _policy = policy
    _rows = {
        (r["class"], r["angle_key"]): r for r in _data["variants"][policy]["lowerings"]
    }
    original = Circuit.append_basic_gate

    def append(parent, gate, *bits):
        if isinstance(gate, GPhase):
            return
        if isinstance(gate, CPhase):
            assert isinstance(gate.phi, Angle)
            return parent.append(LoweredCPhase(str(gate.phi)), qc_reg_from_bits(*bits))
        if isinstance(gate, Phase):
            assert isinstance(gate.phi, Angle)
            return parent.append(LoweredPhase(str(gate.phi)), qc_reg_from_bits(*bits))
        return original(parent, gate, *bits)

    Circuit.append_basic_gate = append
    try:
        yield _data
    finally:
        Circuit.append_basic_gate = original
