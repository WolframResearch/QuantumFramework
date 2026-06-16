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

(* expectation of Z^{xn}: three routes, printed per-n *)
Do[
  Module[{psi = circuit[n][], o = obs[n], eT, eB, eQ, tT, tB, tQ},
    tT = timeMs[eT = Tr[o["MatrixRepresentation"] . psi["DensityMatrix"]]];
    tB = timeMs[eB = (psi["Dagger"][o[psi]])["Scalar"]];
    tQ = timeMs[eQ = QuantumMeasurementOperator[o][psi]["Mean"]];
    Print["expVal n=", n, " : ",
      <|"Tr" -> Re[eT], "bra" -> Re[eB], "qmo" -> Re[eQ],
        "tMs_Tr" -> tT, "tMs_bra" -> tB, "tMs_qmo" -> tQ,
        "maxDiff" -> Max[Abs[{eT - eB, eT - eQ}]]|> // InputForm];
  ],
  {n, {4, 6, 8, 10}}
];

(* single-outcome probability of |0..0>: full-dist vs projector inner product *)
Do[
  Module[{psi = circuit[n][], bit = StringJoin[Table["0", {n}]], pF, pP, tF, tP, pl},
    tF = timeMs[pl = psi["ProbabilitiesList"]]; pF = First[pl];
    tP = timeMs[pP = Abs[(QuantumState[bit]["Dagger"][psi])["Scalar"]]^2];
    Print["singleProb n=", n, " : ",
      <|"pFull" -> pF, "pProj" -> pP, "diff" -> Abs[pF - pP],
        "tMs_full" -> tF, "tMs_proj" -> tP|> // InputForm];
  ],
  {n, {4, 6, 8, 10, 12}}
];

(* explicit dilated superoperator-applied state dimension (measure all n) *)
Do[
  Module[{psi = circuit[n][], qmo = QuantumMeasurementOperator[Range[n]], s},
    s = qmo["SuperOperator"][psi];
    Print["dilatedSuperState n=", n, " : ",
      <|"stateDim" -> s["Dimension"], "2^(2n)" -> 2^(2 n)|> // InputForm];
  ],
  {n, {4, 6, 8}}
];
