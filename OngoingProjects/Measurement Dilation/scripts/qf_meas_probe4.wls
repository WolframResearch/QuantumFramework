#!/usr/bin/env wolframscript
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

SetAttributes[timeMs, HoldFirst];
timeMs[expr_] := (ClearSystemCache[]; Round[1000. First[AbsoluteTiming[expr]], 0.1]);

circuit[n_] := QuantumCircuitOperator[Join[
   Table["H" -> i, {i, n}],
   Table["CNOT" -> {i, i + 1}, {i, n - 1}],
   Table["RY"[1.0/(i + 1)] -> i, {i, n}]
]];

obs[n_] := QuantumOperator[StringJoin[Table["Z", {n}]]];

(* X-basis probabilities two ways: (a) dense state rebase  (b) append H gates, read computational *)
Do[
  Module[{psi, pGate, tGate, pRebase, tRebase},
    psi = circuit[n][];
    (* (b) rotate to X with H gates inside the circuit, then computational probs *)
    tGate = timeMs[
      pGate = (circuit[n] /* QuantumCircuitOperator[Table["H" -> i, {i, n}]])[]["ProbabilitiesList"]];
    Print["xProb_viaGates n=", n, " : ", <|"len" -> Length[pGate], "tMs" -> tGate|> // InputForm];
    If[n <= 8,
      tRebase = timeMs[pRebase = QuantumState[psi, QuantumBasis["PauliX", n]]["ProbabilitiesList"]];
      Print["xProb_viaRebase n=", n, " : ",
        <|"tMs" -> tRebase, "maxDiff" -> Max[Abs[Sort[pGate] - Sort[pRebase]]]|> // InputForm]];
  ],
  {n, {4, 6, 8, 10, 12}}
];

(* bra-op-ket expectation of Z^n stays bounded where QMO route died *)
Do[
  Module[{psi = circuit[n][], o = obs[n], eB, tB},
    tB = timeMs[eB = (psi["Dagger"][o[psi]])["Scalar"]];
    Print["braExp n=", n, " : ", <|"val" -> Re[eB], "tMs" -> tB|> // InputForm];
  ],
  {n, {8, 10, 12}}
];
