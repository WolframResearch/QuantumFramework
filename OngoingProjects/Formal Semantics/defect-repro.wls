(* The QuantumPositiveTransform defect of Appendix D (item 1), then the two problems fixed      *)
(* since 56118a44: a circuit applied to a tableau, and a one-qubit measurement on a tableau.    *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
emit[label_, val_] := WriteString[$Output, label, " -> ", ToString[val, InputForm], "\n"];

(* Defect 1 (Appendix D). QuantumPositiveTransform: the state and operator wrappers fail, the array form works. *)
qs = QuantumState["0", 3];
emit["1a  Through[{Head, ArrayQ}[QuantumPositiveTransform[qs]]] (Schrodinger picture)", Through[{Head, ArrayQ}[QuantumPositiveTransform[qs]]]];
emit["1b  Through[{Head, ArrayQ}[QuantumPositiveTransform[QuantumOperator[\"X\"[3]]]]]", Through[{Head, ArrayQ}[QuantumPositiveTransform[QuantumOperator["X"[3]]]]]];
emit["1b' Head @ QuantumPositiveTransform[QuantumWignerTransform[qs]] (phase-space picture)", Head @ QuantumPositiveTransform[QuantumWignerTransform[qs]]];
emit["1c  Head @ QuantumPositiveTransform[QuditBasis[\"Wigner\"[3]]] (should be QuditBasis)", Head @ QuantumPositiveTransform[QuditBasis["Wigner"[3]]]];
w = QuantumPhaseSpaceTransform[qs]["StateVector"];
emit["1d  Dimensions @ QuantumPositiveTransform[w] (array form)", Dimensions @ QuantumPositiveTransform[w]];
emit["1e  sectors recombine to w", With[{p = Normal @ QuantumPositiveTransform[w]}, p[[1]] - p[[2]] == Normal[w]]];

(* Fixed since 56118a44: a circuit applied to a tableau now returns the resulting PauliStabilizer. *)
qco = QuantumCircuitOperator["Bell"];
ps = PauliStabilizer[2];
emit["2a  MatchQ[qco[ps], _[_PauliStabilizer]]", MatchQ[qco[ps], _[_PauliStabilizer]]];
emit["2b  MatchQ[qco[ps, Method -> \"Stabilizer\"], _[_PauliStabilizer, ___]]", MatchQ[qco[ps, Method -> "Stabilizer"], _[_PauliStabilizer, ___]]];
emit["2c  MatchQ[ps[qco], _[_QuantumCircuitOperator]]", MatchQ[ps[qco], _[_QuantumCircuitOperator]]];
emit["2d  works: QuantumOperator[\"CNOT\"][ps[{\"H\" -> 1}]][\"Stabilizers\"]", QuantumOperator["CNOT"][ps[{"H" -> 1}]]["Stabilizers"]];
emit["2e  works: ps[{\"H\" -> 1, \"CNOT\" -> {1, 2}}][\"Stabilizers\"]", ps[{"H" -> 1, "CNOT" -> {1, 2}}]["Stabilizers"]];
emit["2f  works: qco[Method -> \"Stabilizer\"][\"Stabilizers\"]", qco[Method -> "Stabilizer"]["Stabilizers"]];

(* Fixed since 56118a44: a one-qubit measurement operator on a many-qubit tableau now returns the outcome Association. *)
ps3 = PauliStabilizer[3][{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}];
emit["3a  Head @ QuantumMeasurementOperator[\"Z\", {2}][ps3]", Head @ QuantumMeasurementOperator["Z", {2}][ps3]];
emit["3b  works: Keys @ QuantumMeasurementOperator[\"IZI\"][ps3]", Keys @ QuantumMeasurementOperator["IZI"][ps3]];
emit["3c  works: Keys @ ps3[\"M\", \"IZI\"]", Keys @ ps3["M", "IZI"]];
WriteString[$Output, "DONE\n"];
