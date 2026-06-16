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

(* probabilities in a NON-computational (X) basis directly from the state, no measurement op *)
Do[
  Module[{psi = circuit[n][], pX, tX},
    tX = timeMs[pX = QuantumState[psi, QuantumBasis["PauliX", n]]["ProbabilitiesList"]];
    Print["xBasisProbs n=", n, " : ", <|"len" -> Length[pX], "tMs" -> tX, "tot" -> Total[pX]|> // InputForm];
  ],
  {n, {4, 6, 8, 10}}
];
