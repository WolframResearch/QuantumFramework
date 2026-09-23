BeginTestSection["QuantumPartialTrace - circuit bent pair"]

(* A trace pair {a, b} joins output wire a to input wire b.  On a circuit this is
   drawn as a cup feeding the input wire and a cap closing the output wire.  The
   two were placed on the wrong wires, so for a != b the loop never closed: both
   wires stayed open and the value was a full operator instead of the trace. *)

VerificationTest[
    QuantumPartialTrace[QuantumCircuitOperator[{QuantumOperator["X", {{1}, {2}}]}], {{1, 2}}]["Order"],
    {{}, {}},
    TestID -> "Circuit-bent-pair-closes-order"
]

VerificationTest[
    Normal @ QuantumPartialTrace[QuantumCircuitOperator[{QuantumOperator["X", {{1}, {2}}]}], {{1, 2}}]["QuantumOperator"]["Computational"]["Matrix"],
    {{0}},
    TestID -> "Circuit-bent-pair-value-is-trace-of-X"
]

(* Nontrivial, non-scalar bent trace: contract output 1 against input 2 of CNOT.
   The circuit route must equal the operator route, which is the ground truth,
   and both equal the hand-computed contraction Sum_k <k|_o1 CNOT |k>_i2. *)
VerificationTest[
    Normal @ QuantumPartialTrace[QuantumCircuitOperator[{QuantumOperator["CNOT"]}], {{1, 2}}]["QuantumOperator"]["Computational"]["Matrix"],
    Normal @ QuantumPartialTrace[QuantumOperator["CNOT"], {{1, 2}}]["Computational"]["Matrix"],
    TestID -> "Circuit-bent-pair-matches-operator"
]

VerificationTest[
    Normal @ QuantumPartialTrace[QuantumCircuitOperator[{QuantumOperator["CNOT"]}], {{1, 2}}]["QuantumOperator"]["Computational"]["Matrix"],
    {{1, 1}, {0, 0}},
    TestID -> "Circuit-bent-pair-matches-hand-contraction"
]

(* A same-wire pair already closed correctly and must keep doing so. *)
VerificationTest[
    QuantumPartialTrace[QuantumCircuitOperator[{QuantumOperator["X"]}], {{1, 1}}]["Order"],
    {{}, {}},
    TestID -> "Circuit-same-wire-still-closes"
]

EndTestSection[]


BeginTestSection["QuantumPartialTrace - repeated wire"]

(* A wire can be traced at most once.  A repeated pair used to emit
   TensorContract::lvreps and, on a multi-wire operator, return an operator that
   reported ValidQ True with the contraction left unevaluated, while a one-wire
   operator returned a Failure.  Both must now be rejected the same way. *)

VerificationTest[
    Head @ QuantumPartialTrace[QuantumOperator["CNOT"], {{1, 1}, {1, 1}}],
    Failure,
    TestID -> "Repeated-wire-multiwire-operator-rejected"
]

VerificationTest[
    Head @ QuantumPartialTrace[QuantumOperator["X"], {{1, 1}, {1, 1}}],
    Failure,
    TestID -> "Repeated-wire-onewire-operator-rejected"
]

VerificationTest[
    Head @ QuantumPartialTrace[QuantumState["PhiPlus"], {{1, 1}, {1, 1}}],
    Failure,
    TestID -> "Repeated-wire-state-rejected"
]

(* The guard is not a false positive: distinct pairs on a multi-wire operator
   still trace both wires (the full trace of the two-qubit identity is 4). *)
VerificationTest[
    Normal @ QuantumPartialTrace[QuantumOperator[IdentityMatrix[4], {1, 2}], {{1, 1}, {2, 2}}]["Matrix"],
    {{4}},
    TestID -> "Distinct-pairs-still-trace"
]

EndTestSection[]
