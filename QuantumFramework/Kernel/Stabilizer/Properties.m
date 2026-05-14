Package["Wolfram`QuantumFramework`"]



(* ============================================================================ *)
(* Direct property handlers                                                     *)
(* ============================================================================ *)

ps_PauliStabilizer["Qudits" | "Qubits"] := Dimensions[ps["Tableau"]][[2]]

ps_PauliStabilizer["GeneratorCount"] := Dimensions[ps["Tableau"]][[3]] / 2

ps_PauliStabilizer["StabilizerSigns"] := Drop[ps["Signs"], ps["GeneratorCount"]]
ps_PauliStabilizer["StabilizerTableau" | "Stabilizer"] := Map[Drop[#, ps["GeneratorCount"]] &, ps["Tableau"], {2}]
ps_PauliStabilizer["StabilizerX"] := ps["Stabilizer"][[1]]
ps_PauliStabilizer["StabilizerZ"] := ps["Stabilizer"][[2]]

ps_PauliStabilizer["DestabilizerSigns"] := Take[ps["Signs"], ps["GeneratorCount"]]
ps_PauliStabilizer["DestabilizerTableau" | "Destabilizer"] := Map[Take[#, ps["GeneratorCount"]] &, ps["Tableau"], {2}]
ps_PauliStabilizer["DestabilizerX"] := ps["Destabilizer"][[1]]
ps_PauliStabilizer["DestabilizerZ"] := ps["Destabilizer"][[2]]

ps_PauliStabilizer["Phase"] := (1 - ps["Signs"]) / 2
ps_PauliStabilizer["X"] := ps["Tableau"][[1]]
ps_PauliStabilizer["Z"] := ps["Tableau"][[2]]


(* ============================================================================ *)
(* Flat 2n*2n binary representations                                            *)
(* ============================================================================ *)

ps_PauliStabilizer["Matrix"] := ArrayReshape[Transpose[ps["Tableau"], {2, 3, 1}], 2 {ps["GeneratorCount"], ps["Qubits"]}]

ps_PauliStabilizer["p"] := With[{n = ps["Qudits"]},
    Diagonal[ps["Matrix"] . PadLeft[identityMatrix[n], {-2 n, 2 n}] . Transpose[ps["Matrix"]]]
]

ps_PauliStabilizer["TableauPhase"] := Join[ps["Matrix"], ArrayReshape[ps["Phase"], {2 ps["GeneratorCount"], 1}], 2]


(* ============================================================================ *)
(* Display-form property handlers                                               *)
(* (Definitions of PauliForm, TableauForm live in Formatting.m.)                *)
(* ============================================================================ *)

ps_PauliStabilizer["PauliForm" | "Generators" | "Stabilizers", n_ : Infinity] := PauliForm[ps, n]
ps_PauliStabilizer["Destabilizers", n_ : Infinity] := PauliForm[ps, n, False]
ps_PauliStabilizer["PauliStrings", n_ : Infinity] := Join[ps["Stabilizers", n], ps["Destabilizers", n]]
ps_PauliStabilizer["PauliSymbols", n_ : Infinity] := If[StringStartsQ[#, "-"], - StringDrop[#, 1], #] & /@ Join[ps["Stabilizers", n], ps["Destabilizers", n]]

ps_PauliStabilizer["StabilizerTableauForm", n_ : Infinity] := TableauForm[ps, n]
ps_PauliStabilizer["TableauForm", n_ : Infinity] := TableauForm[Take[ps["Signs"], UpTo[n]], Map[Take[#, UpTo[n]] &, ps["Tableau"], {2}], False]


(* ============================================================================ *)
(* Method-grade symbolic-measurement operations.                                *)
(* Implementations live in Stabilizer/SymbolicMeasure.m as PackageScope helpers.*)
(* ============================================================================ *)

ps_PauliStabilizer["SymbolicMeasure", q_Integer] := symbolicMeasure[ps, q]
ps_PauliStabilizer["SymbolicMeasure", qudits : {___Integer}] := symbolicMeasure[ps, qudits]

ps_PauliStabilizer["SubstituteOutcomes", rules_] := substituteOutcomes[ps, rules]

ps_PauliStabilizer["SampleOutcomes"] := sampleOutcomes[ps]
ps_PauliStabilizer["SampleOutcomes", n_Integer ? Positive] := sampleOutcomes[ps, n]


(* ============================================================================ *)
(* Method-grade InnerProduct + Expectation.                                     *)
(* Implementations live in Stabilizer/InnerProduct.m as PackageScope helpers.   *)
(* ============================================================================ *)

ps_PauliStabilizer["InnerProduct", other_, opts : OptionsPattern[{Method -> "Direct"}]] := stabilizerInnerProduct[ps, other, opts]

ps_PauliStabilizer["Expectation", pauli_String] := stabilizerExpectation[ps, pauli]
