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


BeginTestSection["QuantumPartialTrace - operator basis mismatch on traced wire"]

(* On an operator whose traced wire carries different output and input bases, the trace is that of
   the represented map A = F_out . C . F_in^-1, not of the stored coefficients C.  A stored identity
   with a Pauli-X output basis and a computational input basis is the Hadamard, whose trace is 0;
   the plain index contraction used to return Tr[stored I] = 2. *)
VerificationTest[
    QuantumPartialTrace[
        QuantumOperator[IdentityMatrix[2], QuantumBasis[QuditBasis["PauliX"], QuditBasis[2]]],
        {{1, 1}}
    ]["Number"],
    0,
    TestID -> "Mismatched-basis-traces-the-map-not-the-coefficients"
]

(* The trace of a map does not depend on the basis: the same Hadamard, in three encodings, traces
   to 0.  The first stores it in the computational basis, the second re-expresses it in the X basis. *)
VerificationTest[
    QuantumPartialTrace[QuantumOperator[{{1, 1}, {1, -1}}/Sqrt[2]], {{1, 1}}]["Number"] == 0 &&
    QuantumPartialTrace[QuantumOperator[QuantumOperator[{{1, 1}, {1, -1}}/Sqrt[2]], "PauliX"], {{1, 1}}]["Number"] == 0,
    True,
    TestID -> "Trace-is-basis-independent-across-encodings"
]

(* A non-orthonormal traced basis reconciles through the exact inverse F_in^-1, not the adjoint.
   Stored {{2,3},{5,7}} with output basis columns {{1,0},{1,1}} represents {{7,10},{5,7}}, trace 14. *)
VerificationTest[
    QuantumPartialTrace[
        QuantumOperator[{{2, 3}, {5, 7}},
            QuantumBasis[QuditBasis[Transpose[{{1, 1}, {0, 1}}]], QuditBasis[2]]],
        {{1, 1}}
    ]["Number"],
    14,
    TestID -> "Non-orthonormal-traced-wire-uses-inverse"
]

(* Partial trace keeping a wire: reconciling only the traced wire equals converting the whole
   operator to the computational basis first (the ground-truth route). *)
VerificationTest[
    Normal @ QuantumPartialTrace[
        QuantumOperator[Partition[Range[16], 4], QuantumBasis[QuditBasis[{"PauliX", "PauliX"}], QuditBasis[{2, 2}]]],
        {{1, 1}}
    ]["Computational"]["Matrix"],
    Normal @ QuantumPartialTrace[
        QuantumOperator[Partition[Range[16], 4], QuantumBasis[QuditBasis[{"PauliX", "PauliX"}], QuditBasis[{2, 2}]]]["Computational"],
        {{1, 1}}
    ]["Computational"]["Matrix"],
    TestID -> "Local-reconciliation-matches-computational-route"
]

(* A bent pair a != b with a mismatched traced wire matches the whole-operator route too. *)
VerificationTest[
    Normal @ QuantumPartialTrace[
        QuantumOperator[Partition[Range[16], 4], QuantumBasis[QuditBasis[{"PauliX", 2}], QuditBasis[{2, 2}]]],
        {{1, 2}}
    ]["Computational"]["Matrix"],
    Normal @ QuantumPartialTrace[
        QuantumOperator[Partition[Range[16], 4], QuantumBasis[QuditBasis[{"PauliX", 2}], QuditBasis[{2, 2}]]]["Computational"],
        {{1, 2}}
    ]["Computational"]["Matrix"],
    TestID -> "Bent-pair-mismatched-basis-matches-computational-route"
]

EndTestSection[]
