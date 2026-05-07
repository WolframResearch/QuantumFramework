"""
Generate Stim fixtures for cross-package testing.

For each scenario in `cases`, this script runs Stim's TableauSimulator
and records the resulting stabilizer-generator strings (with signs).
The output JSON is consumed by Tests/CrossPackage_Stim.wlt's TIER G.

Run from repo root:
    python3 Tests/fixtures/generate_stim_fixtures.py

This file is local-only; the JSON it produces is committed so the WL test
can run without requiring stim at test time.
"""

import json
import sys
import os
import stim


def stabilizers_after(gates):
    """Apply gates (list of (name, *targets)) to TableauSimulator
    initialized in the all-|0> state, then return the list of stabilizer
    generator strings.

    Returns a list of strings like "+XX", "-Z" -- positive stabilizers
    use a leading "+", negative use "-".
    """
    sim = stim.TableauSimulator()
    n_qubits = 0
    for g in gates:
        name, *targets = g
        for t in targets:
            n_qubits = max(n_qubits, t + 1)
    # Pre-allocate the qubits Stim needs by issuing a no-op identity on the
    # highest-indexed qubit; this ensures the tableau has space.
    if n_qubits > 0:
        sim.do(stim.CircuitInstruction("I", [n_qubits - 1]))

    for g in gates:
        name, *targets = g
        sim.do(stim.CircuitInstruction(name, targets))

    # current_inverse_tableau().inverse() gives the forward tableau; its
    # row-stabilizers are the post-state stabilizers.
    tableau = sim.current_inverse_tableau().inverse()
    n = len(tableau)
    stabs = []
    for i in range(n):
        # A Tableau encodes Z_i -> P_i conjugation rules. The stabilizers
        # of the post-state |psi> = U|0...0> are U Z_i U^dag, which is
        # exactly tableau.z_output(i).
        stab = tableau.z_output(i)
        stabs.append(str(stab))
    return stabs


cases = {
    # Single-qubit identity start: |0> stabilizer is +Z.
    "single_qubit_zero": {
        "n_qubits": 1,
        "gates": [],
        "expected_stabilizers": None,  # filled in by Stim
    },
    # H|0> = |+> stabilizer is +X.
    "single_qubit_plus_via_H": {
        "n_qubits": 1,
        "gates": [("H", 0)],
    },
    # Bell state via H + CX.
    "bell_phi_plus": {
        "n_qubits": 2,
        "gates": [("H", 0), ("CX", 0, 1)],
    },
    # Bell PhiMinus via H + CX + Z on q1.
    "bell_phi_minus": {
        "n_qubits": 2,
        "gates": [("H", 0), ("CX", 0, 1), ("Z", 0)],
    },
    # GHZ-3 state.
    "ghz_3": {
        "n_qubits": 3,
        "gates": [("H", 0), ("CX", 0, 1), ("CX", 1, 2)],
    },
    # Linear 4-qubit cluster state: H on each + CZ between adjacent.
    "cluster_4": {
        "n_qubits": 4,
        "gates": [("H", 0), ("H", 1), ("H", 2), ("H", 3),
                  ("CZ", 0, 1), ("CZ", 1, 2), ("CZ", 2, 3)],
    },
    # X gate then Z measurement-prep: |1> stabilizer is -Z.
    "single_qubit_one_via_X": {
        "n_qubits": 1,
        "gates": [("X", 0)],
    },
    # |+i> = SH|0> stabilizer is +Y.
    "single_qubit_plus_i_via_SH": {
        "n_qubits": 1,
        "gates": [("H", 0), ("S", 0)],
    },
    # |-i> = S^dag H|0> stabilizer is -Y.
    "single_qubit_minus_i_via_SdagH": {
        "n_qubits": 1,
        "gates": [("H", 0), ("S_DAG", 0)],
    },
    # CNOT^2 = identity check (back to |00> with stabilizers ZI, IZ).
    "cnot_squared_identity": {
        "n_qubits": 2,
        "gates": [("CX", 0, 1), ("CX", 0, 1)],
    },
    # SWAP exchanges |+0>: H on q0, then SWAP. Expected stabilizers IX, ZI.
    "swap_after_H": {
        "n_qubits": 2,
        "gates": [("H", 0), ("SWAP", 0, 1)],
    },
}


def main():
    output = {"_metadata": {"stim_version": stim.__version__,
                            "generator_script": "Tests/fixtures/generate_stim_fixtures.py"}}
    for case_name, case in cases.items():
        gates = case["gates"]
        try:
            stabs = stabilizers_after(gates)
            output[case_name] = {
                "n_qubits": case["n_qubits"],
                "gates": gates,
                "stabilizers": stabs,
            }
        except Exception as e:
            output[case_name] = {"error": str(e), "n_qubits": case["n_qubits"], "gates": gates}

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stim_fixtures.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, sort_keys=True)
    print(f"Wrote {len(cases)} cases to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
