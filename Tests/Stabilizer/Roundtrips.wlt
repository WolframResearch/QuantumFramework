(* ::Package:: *)

(* ============================================================================ *)
(* Cross-module exact-equality round-trip probes (Phase 5c, plan 2026-05-04)    *)
(* ============================================================================ *)
(*                                                                              *)
(* Purpose: catch the same class of bug that hid `PauliStabilizer[Y]["QO"] !=  *)
(* "Y"` for 185 tests. Each tier here probes a constructor-accessor pair `C[x]` *)
(* paired with an accessor `A[result]`. The first test in each tier asserts    *)
(* `A[C[x]] === x` for representative simple `x` (per the user-facing exact-    *)
(* equality round-trip rule, see                                                *)
(* ~/.claude/projects/.../memory/feedback_user_facing_roundtrip_first.md).      *)
(*                                                                              *)
(* Phase 5c surveyed the kernel for constructor-accessor pairs that lack such   *)
(* tests. Results recorded below as PASS / KNOWN-NOT-EXACT / SKIP-CONTRACT.     *)
(* ============================================================================ *)

BeginTestSection["Roundtrips"]


(* Helpers (copies of those in PauliStabilizer.wlt -- duplicated here so this   *)
(* file is independent of test ordering).                                       *)
matEqQO[a_QuantumOperator, b_QuantumOperator] := (a["Matrix"] // Normal // Simplify) === (b["Matrix"] // Normal // Simplify)
matEqQS[a_QuantumState, b_QuantumState]       := (a["StateVector"] // Normal // Simplify) === (b["StateVector"] // Normal // Simplify)
equalUpToGlobalPhaseQS[a_QuantumState, b_QuantumState] := Module[{v1, v2, overlap},
    v1 = a["StateVector"] // Normal // Simplify;
    v2 = b["StateVector"] // Normal // Simplify;
    overlap = Quiet @ Simplify[Abs[Conjugate[v1] . v2]];
    Quiet @ Simplify[overlap == 1] === True
]


(* ============================================================================ *)
(* TIER QO-QC: QuantumOperator -> QuantumChannel -> QuantumOperator            *)
(* ============================================================================ *)
(* Contract: QuantumChannel[qo]["QuantumOperator"] === qo for unitary qo.      *)
(* Result: PASS for all canary cases (I, X, Y, Z, H, S, CNOT, YY, XY).         *)

VerificationTest[matEqQO[QuantumChannel[QuantumOperator["I"]]["QuantumOperator"],    QuantumOperator["I"]],    True, TestID -> "QO-QC-I"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["X"]]["QuantumOperator"],    QuantumOperator["X"]],    True, TestID -> "QO-QC-X"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["Y"]]["QuantumOperator"],    QuantumOperator["Y"]],    True, TestID -> "QO-QC-Y"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["Z"]]["QuantumOperator"],    QuantumOperator["Z"]],    True, TestID -> "QO-QC-Z"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["H"]]["QuantumOperator"],    QuantumOperator["H"]],    True, TestID -> "QO-QC-H"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["S"]]["QuantumOperator"],    QuantumOperator["S"]],    True, TestID -> "QO-QC-S"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["CNOT"]]["QuantumOperator"], QuantumOperator["CNOT"]], True, TestID -> "QO-QC-CNOT"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["YY"]]["QuantumOperator"],   QuantumOperator["YY"]],   True, TestID -> "QO-QC-YY"]
VerificationTest[matEqQO[QuantumChannel[QuantumOperator["XY"]]["QuantumOperator"],   QuantumOperator["XY"]],   True, TestID -> "QO-QC-XY"]


(* ============================================================================ *)
(* TIER QS-QO-Projector: QuantumState -> QuantumOperator (projector) -> ?      *)
(* ============================================================================ *)
(* Contract clarification: there is NO direct `QuantumState[qo] === qs`        *)
(* round-trip in QF. `qs["Operator"]` returns the projector |psi><psi|;        *)
(* `QuantumState[qo]` for a 2^n x 2^n matrix `qo` returns a 4^n state vector   *)
(* (the matrix vectorization), not the +1 eigenvector. The two paths are not   *)
(* inverses of each other.                                                      *)
(*                                                                              *)
(* The well-defined recovery is via eigendecomposition: the +1 eigenvector of  *)
(* the projector reproduces qs *up to* a global phase. These tests document    *)
(* that contract.                                                               *)

(* Helper: extract the +1 eigenvector of a rank-1 projector *)
projectorTopEigenvector[qo_QuantumOperator] := Module[{eig},
    eig = Eigensystem[qo["Matrix"] // Normal // Simplify, 1];
    Normalize[eig[[2, 1]]]
]

VerificationTest[
    Quiet @ Simplify[Abs[Conjugate[projectorTopEigenvector[QuantumState[{1, 0}]["Operator"]]] . (QuantumState[{1, 0}]["StateVector"] // Normal)] == 1],
    True,
    TestID -> "QS-QO-Projector-Comp-Zero"
]
VerificationTest[
    Quiet @ Simplify[Abs[Conjugate[projectorTopEigenvector[QuantumState["Plus"]["Operator"]]] . (QuantumState["Plus"]["StateVector"] // Normal)] == 1],
    True,
    TestID -> "QS-QO-Projector-Plus"
]
VerificationTest[
    Quiet @ Simplify[Abs[Conjugate[projectorTopEigenvector[QuantumState["Bell"]["Operator"]]] . (QuantumState["Bell"]["StateVector"] // Normal)] == 1],
    True,
    TestID -> "QS-QO-Projector-Bell"
]
VerificationTest[
    Quiet @ Simplify[Abs[Conjugate[projectorTopEigenvector[QuantumState[{1, I}/Sqrt[2]]["Operator"]]] . (QuantumState[{1, I}/Sqrt[2]]["StateVector"] // Normal)] == 1],
    True,
    TestID -> "QS-QO-Projector-PlusI"
]


(* ============================================================================ *)
(* TIER QO-QCO: QuantumOperator -> QuantumCircuitOperator -> QuantumOperator   *)
(* ============================================================================ *)
(* Contract: QuantumCircuitOperator[qo]["QuantumOperator"] === qo.             *)
(* Result: PASS for I, X, H, CNOT. Y/YY are the canaries -- if the round-trip *)
(* loses the i factor here, it's the same class of bug as the PauliStabilizer  *)
(* one and a new ROADMAP entry should be filed.                                 *)

VerificationTest[matEqQO[QuantumCircuitOperator[QuantumOperator["I"]]["QuantumOperator"], QuantumOperator["I"]], True, TestID -> "QO-QCO-I"]
VerificationTest[matEqQO[QuantumCircuitOperator[QuantumOperator["X"]]["QuantumOperator"], QuantumOperator["X"]], True, TestID -> "QO-QCO-X"]
VerificationTest[matEqQO[QuantumCircuitOperator[QuantumOperator["Y"]]["QuantumOperator"], QuantumOperator["Y"]], True, TestID -> "QO-QCO-Y"]
VerificationTest[matEqQO[QuantumCircuitOperator[QuantumOperator["H"]]["QuantumOperator"], QuantumOperator["H"]], True, TestID -> "QO-QCO-H"]
VerificationTest[matEqQO[QuantumCircuitOperator[QuantumOperator["CNOT"]]["QuantumOperator"], QuantumOperator["CNOT"]], True, TestID -> "QO-QCO-CNOT"]
VerificationTest[matEqQO[QuantumCircuitOperator[QuantumOperator["YY"]]["QuantumOperator"], QuantumOperator["YY"]], True, TestID -> "QO-QCO-YY"]


(* ============================================================================ *)
(* TIER QC-QC-Compose: QuantumChannel composition with itself                  *)
(* ============================================================================ *)
(* Contract: a unitary channel composed with its own adjoint is the identity. *)

VerificationTest[
    matEqQO[
        QuantumChannel[QuantumOperator["H"]]["QuantumOperator"] @ QuantumChannel[QuantumOperator["H"]]["QuantumOperator"],
        QuantumOperator["I", {1}]
    ],
    True,
    TestID -> "QC-Compose-H-Squared"
]


EndTestSection[]
