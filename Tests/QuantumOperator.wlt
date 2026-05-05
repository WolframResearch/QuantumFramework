(* QuantumOperator tests. *)

BeginTestSection["QuantumOperator - constructors"]

VerificationTest[QuantumOperator[]["Dimensions"], {2, 2}, TestID -> "Empty"];

VerificationTest[QuantumOperator["X"]["Dimensions"], {2, 2}, TestID -> "X-bare"];

VerificationTest[QuantumOperator["X"[3]]["Dimensions"], {3, 3}, TestID -> "X-3"];

VerificationTest[QuantumOperator["Hadamard"]["Dimensions"], {2, 2}, TestID -> "Hadamard"];

VerificationTest[QuantumOperator["CNOT"]["Dimensions"], {2, 2, 2, 2}, TestID -> "CNOT"];

(* CNOT[3] gives a 2-dim control with a 3-dim target (the qutrit "shift" gate),
   matching pre-refactor semantics. *)
VerificationTest[QuantumOperator["CNOT"[3]]["Dimensions"], {2, 3, 2, 3}, TestID -> "CNOT-3"];

VerificationTest[QuantumOperator["Toffoli"]["Dimensions"], {2, 2, 2, 2, 2, 2}, TestID -> "Toffoli"];

VerificationTest[QuantumOperator["Fourier"[3]]["Dimensions"], {3, 3}, TestID -> "Fourier-3"];

VerificationTest[QuantumOperator["XRotation"[Pi/3]]["Dimensions"], {2, 2}, TestID -> "XRotation"];

VerificationTest[QuantumOperator["Phase"[Pi/4]]["Dimensions"], {2, 2}, TestID -> "Phase"];

VerificationTest[QuantumOperator["Permutation"[{2, 2}, Cycles[{{1, 2}}]]]["Dimensions"], {2, 2, 2, 2}, TestID -> "Permutation"];

VerificationTest[QuantumOperator["Spider"]["Dimensions"], {2, 2}, TestID -> "Spider-bare"];

VerificationTest[QuantumOperator["Curry"]["Dimensions"], {4, 2, 2}, TestID -> "Curry-default"];

VerificationTest[QuantumOperator["XSpider"[Pi/2], {{1, 2}, {3}}]["Dimensions"], {2, 2, 2}, TestID -> "XSpider-with-order"];

VerificationTest[
  QuantumOperator["Liouvillian"[QuantumOperator["X"]]]["Dimensions"],
  {2, 2},
  TestID -> "Liouvillian"
];

(* explicit matrix *)
VerificationTest[
  QuantumOperator[{{0, 1}, {1, 0}}]["Dimensions"],
  {2, 2},
  TestID -> "Explicit-matrix"
];

VerificationTest[QuantumOperator["XYZ"]["Dimensions"], {2, 2, 2, 2, 2, 2}, TestID -> "PauliString-XYZ"];

EndTestSection[]


BeginTestSection["QuantumOperator - properties"]

VerificationTest[QuantumOperator["X"]["UnitaryQ"], True, TestID -> "X-Unitary"];

VerificationTest[QuantumOperator["Hadamard"]["UnitaryQ"], True, TestID -> "Hadamard-Unitary"];

VerificationTest[QuantumOperator["CNOT"]["Arity"], 2, TestID -> "CNOT-Arity"];

VerificationTest[
  QuantumOperator["X"]["MatrixRepresentation"] - {{0, 1}, {1, 0}} // Total // Total,
  0,
  TestID -> "X-Matrix"
];

EndTestSection[]


BeginTestSection["QuantumOperator - composition"]

VerificationTest[
  (QuantumOperator["X"] @ QuantumOperator["X"])["MatrixRepresentation"] - IdentityMatrix[2] // Total // Total,
  0,
  TestID -> "XX=I"
];

VerificationTest[
  Chop[Flatten[QuantumOperator["X"][QuantumState["0"]]["StateVector"]] - {0, 1}] // Total,
  0,
  TestID -> "X-on-Zero"
];

EndTestSection[]


BeginTestSection["QuantumOperator - shortcut roundtrip"]

(* For each common named gate, QuantumShortcut should emit the new call-form
   shorthand, and feeding it through FromOperatorShorthand should reconstruct
   an operator with the same matrix. *)

roundtripQ[op_] := Norm[Flatten[
  op["Matrix"] - Head[op][First[QuantumShortcut[op]]]["Matrix"]
]] // Chop // # === 0 &;

VerificationTest[roundtripQ[QuantumOperator["X"]], True, TestID -> "Roundtrip-X"];
VerificationTest[roundtripQ[QuantumOperator["X"[3]]], True, TestID -> "Roundtrip-X-3"];
VerificationTest[roundtripQ[QuantumOperator["Y"]], True, TestID -> "Roundtrip-Y"];
VerificationTest[roundtripQ[QuantumOperator["Z"]], True, TestID -> "Roundtrip-Z"];
VerificationTest[roundtripQ[QuantumOperator["Hadamard"]], True, TestID -> "Roundtrip-Hadamard"];
VerificationTest[roundtripQ[QuantumOperator["NOT"]], True, TestID -> "Roundtrip-NOT"];
VerificationTest[roundtripQ[QuantumOperator["S"]], True, TestID -> "Roundtrip-S"];
VerificationTest[roundtripQ[QuantumOperator["T"]], True, TestID -> "Roundtrip-T"];
VerificationTest[roundtripQ[QuantumOperator["Phase"[Pi/4]]], True, TestID -> "Roundtrip-Phase"];
VerificationTest[roundtripQ[QuantumOperator["PhaseShift"[2]]], True, TestID -> "Roundtrip-PhaseShift"];
VerificationTest[roundtripQ[QuantumOperator["U"[Pi/3, Pi/4, Pi/5]]], True, TestID -> "Roundtrip-U"];
VerificationTest[roundtripQ[QuantumOperator["U2"[Pi/4, Pi/3]]], True, TestID -> "Roundtrip-U2"];

(* QuantumShortcut emits the new call form, not the legacy list form, for any
   parameterized gate. *)
VerificationTest[
  QuantumShortcut[QuantumOperator["U"[Pi/3, Pi/4, Pi/5]]],
  {"U"[Pi/3, Pi/4, Pi/5] -> {1}},
  TestID -> "Shortcut-U-callform"
];

VerificationTest[
  QuantumShortcut[QuantumOperator["X"[3]]],
  {"X"[3] -> {1}},
  TestID -> "Shortcut-X3-callform"
];

EndTestSection[]


BeginTestSection["QuantumOperator - failure"]

VerificationTest[
  Quiet @ QuantumOperator["NotAnActualGate"[3]],
  Failure["InvalidName", _],
  SameTest -> MatchQ,
  TestID -> "InvalidName-call-form"
];

EndTestSection[]
