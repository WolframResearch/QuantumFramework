(* section3-audit-check.wls: the combining rules of section 4 and the circuit-and-operator
   defects of Appendix D (items 2 to 7), run against the working tree. Each line is the
   call followed by what the kernel returns; messages are left in the output. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
out[label_, expr_] := WriteString[$Output, label, ": ", ToString[expr, InputForm], "\n"];

(* Application: padding by qubits to the full input order, spectators, the local-dimension check *)
out["A1 X on wire 3 applied to |1>", Normal @ QuantumOperator["X", {3}][QuantumState["1"]]["StateVector"]];
out["A1b full input order of that one-element circuit", QuantumCircuitOperator[{QuantumOperator["X", {3}]}]["FullInputOrder"]];
out["A1c X[3] on wire 2 applied to |1>, a missing qutrit", QuantumOperator["X"[3], {2}][QuantumState["1"]]];
out["A1d X on wire 1 applied to a qubit-qutrit state: dimensions of the result", QuantumOperator["X", {1}][QuantumState[{1, 0, 0, 0, 0, 0}, {2, 3}]]["Dimensions"]];

(* Composition: operands sorted first, the boundary formula, disjoint wires, mixed pictures, parameter intervals *)
c31 = QuantumOperator["CNOT", {3, 1}];
out["A2 CNOT{3,1}[CNOT{3,1}] order", c31[c31]["Order"]];
out["A2b CNOT{2,4}[CNOT{1,2}] order", QuantumOperator["CNOT", {2, 4}][QuantumOperator["CNOT", {1, 2}]]["Order"]];
r = QuantumOperator["X", {1}][QuantumOperator["Z", {2}]];
out["A2c X1[Z2] order, matrix == X(x)Z", {r["Order"], Normal @ r["Matrix"] === KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]}];
m = QuantumOperator["X", QuantumBasis["Picture" -> "Heisenberg"]][QuantumOperator["Z"]];
out["A2d Heisenberg X composed with Schrodinger Z: head, picture (defect 2)", {Head[m], m["Picture"]}];
out["A2e tensor product across pictures: head", Head @ QuantumTensorProduct[QuantumOperator["X", QuantumBasis["Picture" -> "Heisenberg"]], QuantumOperator["Z"]]];
p1 = QuantumOperator["RX"[\[Theta]], "Parameters" -> {{\[Theta], 0, 1}}];
p2 = QuantumOperator["RZ"[\[Theta]], "Parameters" -> {{\[Theta], 2, 3}}];
out["A2f parameter intervals on the direct path", p1[p2]["ParameterSpec"]];
p3 = QuantumOperator[QuantumOperator["RZ"[\[Theta]]]["Matrix"], QuantumBasis["PauliX"], "Parameters" -> {{\[Theta], 2, 3}}];
out["A2g parameter intervals on the fallback forced by a basis mismatch", p1[p3]["ParameterSpec"]];

(* Tensor product: the empty product and a mixed pair *)
out["A3 QuantumTensorProduct[]", QuantumTensorProduct[]];
out["A3b operator with channel: head", Head @ QuantumTensorProduct[QuantumOperator["X"], QuantumChannel["BitFlip"[1/10]]]];

(* Dagger: coordinates versus the adjoint of the map; the channel adjoint keeps its order *)
f = QuantumOperator[{{1, I}, {2, 3}}, QuantumBasis["Output" -> QuditBasis[{{1, 0}, {0, 2}}], "Input" -> QuditBasis[2]]];
out["A4 rep(f), rep(f dagger), rep(f)^dagger", {Normal @ f["MatrixRepresentation"], Normal @ f["Dagger"]["MatrixRepresentation"], ConjugateTranspose @ Normal @ f["MatrixRepresentation"]}];
out["A4b after Computational, rep(f dagger) == rep(f)^dagger", Normal @ f["Computational"]["Dagger"]["MatrixRepresentation"] === ConjugateTranspose @ Normal @ f["MatrixRepresentation"]];
bf = QuantumChannel["BitFlip"[1/10]];
out["A4c BitFlip order, Adjoint order, SuperDagger order", {bf["Order"], bf["Adjoint"]["Order"], SuperDagger[bf]["Order"]}];

(* Partial trace: state positions; an operator pair without a change of basis (defect 4); circuit labels (defect 3); a repeated pair (defect 7) *)
pt = QuantumPartialTrace[QuantumState["PhiPlus"], {1}];
out["A5 trace of qudit 1 of PhiPlus: matrix, InputQudits, StateType", {Normal @ pt["DensityMatrix"], pt["InputQudits"], pt["StateType"]}];
h = QuantumOperator[IdentityMatrix[2], QuantumBasis["Output" -> QuditBasis["PauliX"], "Input" -> QuditBasis[2]]];
out["A5b stored I2 with PauliX output and computational input: MatrixRepresentation", Normal @ h["MatrixRepresentation"]];
out["A5c its pair trace, then the same after Computational (defect 4)", {Normal @ QuantumPartialTrace[h, {{1, 1}}]["Matrix"], Normal @ QuantumPartialTrace[h["Computational"], {{1, 1}}]["Matrix"]}];
qcx = QuantumCircuitOperator[{QuantumOperator["X", {{1}, {2}}]}];
out["A5d circuit X from wire 2 to wire 1 traced {{1,2}}: order, then label and order of each element (defect 3)", {QuantumPartialTrace[qcx, {{1, 2}}]["Order"], Map[{#["Label"], #["Order"]} &, QuantumPartialTrace[qcx, {{1, 2}}]["Elements"]]}];
out["A5e circuit X traced {{1,1}}: order", QuantumPartialTrace[QuantumCircuitOperator[{QuantumOperator["X"]}], {{1, 1}}]["Order"]];
rep = QuantumPartialTrace[QuantumOperator["CNOT"], {{1, 1}, {1, 1}}];
out["A5f repeated pair: head, ValidQ, contraction left unevaluated (defect 7)", {Head[rep], rep["ValidQ"], ! FreeQ[rep, TensorContract]}];

(* Sum and scalar: union of wires; the coordinate diagonal; the density level (defect 5) *)
s = QuantumOperator["X", {1}] + QuantumOperator["Z", {2}];
out["A6 X1 + Z2: order, matrix == X(x)I + I(x)Z", {s["Order"], Normal @ s["Matrix"] === KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]] + KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]]}];
out["A6b rep(1 + h) - rep(h) for PauliX output and computational input", Normal[(1 + h)["MatrixRepresentation"]] - Normal[h["MatrixRepresentation"]]];
out["A6c 2 + X[ToMatrix]: head (defect 5)", Head[2 + QuantumOperator["X"]["ToMatrix"]]];

(* Measurement: PostMeasurementState is the target qudits; measuring circuits; the value re-enters a circuit; direct application shifts *)
qm = QuantumMeasurementOperator["Z", {1}][QuantumState["PhiPlus"]];
out["A7 PostMeasurementState qudits and matrix; qudits of each branch state", {qm["PostMeasurementState"]["OutputQudits"], Normal @ qm["PostMeasurementState"]["DensityMatrix"], Map[#["OutputQudits"] &, qm["States"]]}];
out["A7b measuring circuit: qco[] and qco[qs] heads", {Head @ QuantumCircuitOperator[{"H", {1}}][], Head @ QuantumCircuitOperator[{"H", {1}}][QuantumState["0"]]}];
out["A7c circuit applied to a QuantumMeasurement: head", Head @ QuantumCircuitOperator[{"X" -> 1}][QuantumMeasurementOperator["Z", {1}][QuantumState["+"]]]];
out["A7d Z on wire 3 applied to |1> directly, then inside a circuit: probabilities", {Normal @ QuantumMeasurementOperator["Z", {3}][QuantumState["1"]]["ProbabilitiesList"], Normal @ QuantumCircuitOperator[{QuantumMeasurementOperator["Z", {3}]}][QuantumState["1"]]["ProbabilitiesList"]}];

(* Wrappers and circuits: a barrier lost on append; a channel applied to a circuit compiles (defect 6) *)
qcb = QuantumCircuitOperator[{"H", "Barrier", "X"}];
out["A8 qc[Z] elements, Z[qc] elements", {Length @ qcb[QuantumOperator["Z"]]["Elements"], Length @ QuantumOperator["Z"][qcb]["Elements"]}];
out["A8b chan[qc]: head, with its message (defect 6)", Head @ bf[QuantumCircuitOperator[{"H"}]]];

(* Evolution: None in place of the state gives the propagator *)
out["A9 QuantumEvolve[Z, None, t]: head", Head @ QuantumEvolve[QuantumOperator["Z"], None, t]];
