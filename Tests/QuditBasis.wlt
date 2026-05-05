(* QuditBasis tests.

   Covers:
     - all public construction schemes after the v2.0 list -> call-form refactor
     - basic property access (Names, Dimension, Elements, Qudits)
     - tensor-product composition
     - the Failure path for unknown names
*)

BeginTestSection["QuditBasis - constructors"]

VerificationTest[QuditBasis["Bell"]["Dimension"], 4, TestID -> "Bell-bare"];

VerificationTest[QuditBasis["Fourier"]["Dimension"], 2, TestID -> "Fourier-default"];

VerificationTest[QuditBasis["Fourier"[3]]["Dimension"], 3, TestID -> "Fourier-call-3"];

VerificationTest[QuditBasis["PauliX"[3]]["Dimension"], 3, TestID -> "PauliX-3"];

VerificationTest[QuditBasis["GellMann"[3]]["Dimension"], 9, TestID -> "GellMann-3"];

VerificationTest[QuditBasis["HesseSIC"]["Dimension"], 9, TestID -> "HesseSIC-bare"];

VerificationTest[QuditBasis["HoggarSIC"]["Dimension"], 64, TestID -> "HoggarSIC-bare"];

(* tensor-product LIST stays as a list (semantic composition) *)
VerificationTest[QuditBasis["Computational"[{2, 3, 4}]]["Dimension"], 24, TestID -> "Computational-tensor-list"];

(* alternative dispatch: Pauli string of single-letter names *)
VerificationTest[QuditBasis["XYZ"]["Dimension"], 8, TestID -> "PauliString-XYZ"];

VerificationTest[QuditBasis[3]["Dimension"], 3, TestID -> "Integer-3"];

VerificationTest[QuditBasis[{"a", "b"}, {{1, 0}, {0, 1}}]["Dimension"], 2, TestID -> "Explicit-elements"];

EndTestSection[]


BeginTestSection["QuditBasis - properties"]

VerificationTest[QuditBasis["Pauli"]["Qudits"], 1, TestID -> "Pauli-Qudits"];

VerificationTest[Length[QuditBasis["Pauli"]["Names"]], 4, TestID -> "Pauli-NameCount"];

VerificationTest[QuditBasis["Computational"[2]]["ElementDimensions"], {2}, TestID -> "Comp2-ElementDims"];

VerificationTest[QuditBasis["GellMann"[2]]["Dimension"], 4, TestID -> "GellMann-2-dim"];

EndTestSection[]


BeginTestSection["QuditBasis - tensor product"]

VerificationTest[
  QuantumTensorProduct[QuditBasis["Pauli"], QuditBasis["Pauli"]]["Dimension"],
  16,
  TestID -> "Pauli-tensor-Pauli"
];

VerificationTest[
  QuditBasis["Computational"[{2, 3}]]["Dimensions"],
  {2, 3},
  TestID -> "Comp-tensor-dimensions"
];

EndTestSection[]


BeginTestSection["QuditBasis - failure"]

(* The cache wrapper turns InvalidName into a ConfirmationFailed wrapping
   the InvalidName, so check for either. *)
VerificationTest[
  Quiet @ QuditBasis["NotAnActualName"["bar"]],
  Failure["InvalidName" | "ConfirmationFailed", _],
  SameTest -> MatchQ,
  TestID -> "InvalidName-call-form"
];

EndTestSection[]
