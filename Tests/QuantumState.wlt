(* QuantumState tests. *)

BeginTestSection["QuantumState - constructors"]

VerificationTest[QuantumState[]["Dimension"], 2, TestID -> "Empty-default-0"];

VerificationTest[QuantumState["Bell"]["Dimension"], 4, TestID -> "Bell-bare"];

VerificationTest[QuantumState["GHZ"]["Dimension"], 8, TestID -> "GHZ-default"];

VerificationTest[QuantumState["GHZ"[4]]["Dimension"], 16, TestID -> "GHZ-4"];

VerificationTest[QuantumState["W"[3]]["Dimension"], 8, TestID -> "W-3"];

VerificationTest[QuantumState["Werner"]["Dimension"], 4, TestID -> "Werner-default"];

VerificationTest[QuantumState["Werner"[1/3, 2]]["Dimension"], 4, TestID -> "Werner-call"];

VerificationTest[QuantumState["Dicke"[4, 2]]["Dimension"], 16, TestID -> "Dicke-4-2"];

VerificationTest[QuantumState["BlochVector"[{0, 1, 0}]]["Dimension"], 2, TestID -> "BlochVector"];

VerificationTest[QuantumState["UniformSuperposition"[3]]["Dimension"], 8, TestID -> "UniformSuperposition-3"];

VerificationTest[QuantumState["UniformMixture"[2]]["Dimension"], 4, TestID -> "UniformMixture-2"];

VerificationTest[QuantumState["Plus"]["Dimension"], 2, TestID -> "Plus-bare"];

VerificationTest[QuantumState["+"]["Dimension"], 2, TestID -> "+-shorthand"];

VerificationTest[QuantumState["101"]["Dimension"], 8, TestID -> "Digit-string-101"];

VerificationTest[QuantumState[{1/Sqrt[2], 1/Sqrt[2]}]["Dimension"], 2, TestID -> "Explicit-vector"];

EndTestSection[]


BeginTestSection["QuantumState - properties"]

VerificationTest[QuantumState["Bell"]["Qudits"], 2, TestID -> "Bell-Qudits"];

VerificationTest[QuantumState["Bell"]["Purity"] == 1, True, TestID -> "Bell-PurePurity"];

VerificationTest[Chop[QuantumState["UniformMixture"[2]]["Purity"] - 1/4], 0, TestID -> "UnifMix-Purity"];

VerificationTest[QuantumState["Bell"]["MatrixRepresentation"] // Dimensions, {4, 4}, TestID -> "Bell-MatrixDim"];

EndTestSection[]


BeginTestSection["QuantumState - bipartite operations"]

VerificationTest[
  QuantumPartialTrace[QuantumState["Bell"], {1}]["Dimension"],
  2,
  TestID -> "Bell-PartialTrace-1"
];

VerificationTest[
  Chop[QuantityMagnitude @ UnitConvert[
    QuantumPartialTrace[QuantumState["Bell"], {1}]["VonNeumannEntropy"], "Bits"
  ] - 1],
  0,
  TestID -> "Bell-EntropyHalf"
];

EndTestSection[]


BeginTestSection["QuantumState - failure"]

VerificationTest[
  Quiet @ QuantumState["NotAnActualName"["bar"]],
  Failure["InvalidName", _],
  SameTest -> MatchQ,
  TestID -> "InvalidName-call-form"
];

EndTestSection[]
