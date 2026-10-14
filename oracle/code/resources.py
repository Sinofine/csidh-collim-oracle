"""Qubit, gate, and unresolved phase-operation counts."""

from qarton.circuit import ResourceReporter


def report(circuit):
    r = ResourceReporter(circuit)
    counts = r.sub_circuits_by_class()
    return dict(
        qubits=r.nbr_qubits(),
        classical_bits=circuit.nbr_cbits(),
        ccx=r.ccx_count(),
        native_t=r.t_count(),
        cx=counts.get("cx", 0),
        unresolved_phase_gates={
            k: counts.get(k, 0) for k in ["Phase", "CPhase", "GPhase"]
        },
        ccx_7T_conversion_only=7 * r.ccx_count(),
        primitives={
            k: v
            for k, v in counts.items()
            if k in ["x", "cx", "ccx", "h", "t", "tdg", "msr"]
        },
    )
