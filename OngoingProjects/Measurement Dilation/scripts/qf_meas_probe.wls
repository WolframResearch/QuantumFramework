#!/usr/bin/env wolframscript
(* Probe: cost of computing probabilities / expectations via
   (1) measurement dilation (QuantumMeasurementOperator),
   (2) state -> probabilities,
   (3) inner-product bra-op-ket (QCO with Dagger).
   Working-tree paclet. No Quiet, no Print; final expression is the result Association. *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

SetAttributes[timeMs, HoldFirst];
timeMs[expr_] := (ClearSystemCache[]; Round[1000. First[AbsoluteTiming[expr]], 0.1]);

(* structured circuit: H wall, CNOT ladder, RY rotations -> nontrivial Born probs *)
circuit[n_] := QuantumCircuitOperator[Join[
   Table["H" -> i, {i, n}],
   Table["CNOT" -> {i, i + 1}, {i, n - 1}],
   Table["RY"[1.0/(i + 1)] -> i, {i, n}]
]];

(* ---- A. dimension blow-up of dilation: measure all n qubits as n SEPARATE single-qubit measurements ---- *)
(* sequential single-qubit measurements each append one pointer qudit -> 2^(n+k) *)
seqMeasCirc[n_, k_] := circuit[n] /* QuantumCircuitOperator[List /@ Range[k]];

(* superoperator output dimension for measuring all n qubits (one joint register) *)
supDim[n_] := Module[{qmo = QuantumMeasurementOperator[Range[n]]},
  <|"sysDim" -> 2^n, "supOutDim" -> qmo["SuperOperator"]["OutputDimension"], "2^(2n)" -> 2^(2 n)|>];

(* ---- run for a range of n ---- *)
ns = {4, 6, 8, 10, 12};

psiCache = Association[];
getPsi[n_] := Lookup[psiCache, n, psiCache[n] = circuit[n][]];

(* Route 2: state -> full probability distribution *)
route2[n_] := Module[{psi, t, pl},
  psi = getPsi[n];
  t = timeMs[pl = psi["ProbabilitiesList"]];
  <|"n" -> n, "len" -> Length[pl], "tMs" -> t, "p0" -> First[pl]|>
];

(* Route 1: full measurement applied (dilation), all n qubits, as ONE joint register *)
route1joint[n_] := Module[{qm, t, pl},
  t = timeMs[qm = (circuit[n] /* QuantumCircuitOperator[{Range[n]}])[]];
  pl = qm["ProbabilitiesList"];
  <|"n" -> n, "head" -> Head[qm], "tMs" -> t, "p0" -> First[pl]|>
];

(* Route 1b: n SEPARATE single-qubit measurements (pointer accumulation) *)
route1seq[n_] := Module[{qm, t},
  t = timeMs[qm = seqMeasCirc[n, n][]];
  <|"n" -> n, "head" -> Head[qm], "tMs" -> t|>
];

(* agreement check: do route1 and route2 give the same distribution? *)
agree[n_] := Module[{psi = getPsi[n], pl2, pl1},
  pl2 = psi["ProbabilitiesList"];
  pl1 = (circuit[n] /* QuantumCircuitOperator[{Range[n]}])[]["ProbabilitiesList"];
  Max[Abs[pl1 - pl2]]
];

Print["supDimTable: ", Association[(# -> supDim[#]) & /@ {4, 6, 8}] // InputForm];
Print["route2_stateProbs: ", (route2 /@ ns) // InputForm];
Print["route1_jointMeas: ", (route1joint /@ {4, 6, 8, 10}) // InputForm];
Print["route1b_seqMeas: ", (route1seq /@ {4, 6, 8, 10}) // InputForm];
Print["agreement_maxAbsDiff: ", Association[(# -> agree[#]) & /@ {4, 6, 8}] // InputForm];
