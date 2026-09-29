PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

emit[label_, val_] := WriteString[$Output, label, " -> ", ToString[val, InputForm], "\n"];
zeroQ[m_] := Max[Abs[Flatten[N[Normal[m]]]]] < 10^-9;

(* ---------- the two storage shapes ---------- *)
ps0 = PauliStabilizer[3];
emit["register: head, keys, stabilizers", {Head[ps0], Keys[First[ps0]], ps0["Stabilizers"]}];
ps1 = ps0["H", 1];
emit["one gate: head, keys (packed form)", {Head[ps1], Keys[First[ps1]]}];
ps2 = ps0[{"H", "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}];
emit["gate list: head, keys, stabilizers", {Head[ps2], Keys[First[ps2]], ps2["Stabilizers"]}];
emit["packed and canonical agree on the tableau", ps1["Tableau"] == ps0["H" -> 1]["Tableau"]];
emit["StabilizerStateQ on ps, on a dense Bell state", {StabilizerStateQ[ps2], StabilizerStateQ[QuantumState["Bell"]]}];

(* ---------- the two conversions between the representations ---------- *)
qs2 = QuantumState[ps2];
emit["QuantumState[ps]: head, dimension, equals dense GHZ", {Head[qs2], qs2["Dimension"], qs2 == QuantumState["GHZ"[3]]}];
back = PauliStabilizer[qs2];
emit["PauliStabilizer[QuantumState[ps]]: head, keys, same stabilizers", {Head[back], Keys[First[back]], back["Stabilizers"] == ps2["Stabilizers"]}];
emit["round trip keeps the global phase: PauliStabilizer[qs]['State'] == qs", PauliStabilizer[I qs2]["State"] == I qs2];
emit["PauliStabilizer of a non-stabilizer state", Head[PauliStabilizer[QuantumState[{1, Exp[I Pi/4]}/Sqrt[2]]]]];
emit["QuantumCircuitOperator[ps]: head, gates; QuantumOperator[ps]: head", {Head[QuantumCircuitOperator[ps2]], QuantumCircuitOperator[ps2]["Gates"], Head[QuantumOperator[ps2]]}];
emit["synthesized circuit reproduces the state", QuantumCircuitOperator[ps2][] == qs2];

(* ---------- application through the tensor-representation heads (upvalues) ---------- *)
emit["qo[ps] Clifford: head", Head[QuantumOperator["H"][ps0]]];
emit["qo[ps] non-Clifford T: head", Head[QuantumOperator["T"][ps0]]];
emit["qo[ps] generic unitary: result", QuantumOperator["RandomUnitary"[2]][ps0]];
emit["qco[ps] Clifford circuit: head", Head[QuantumCircuitOperator["Bell"][ps0]]];
emit["qco[] with Method Stabilizer: head", Head[QuantumCircuitOperator["Bell"][Method -> "Stabilizer"]]];
emit["qco[qs] with Method Stabilizer: head", Head[QuantumCircuitOperator["Bell"][QuantumState["00"], Method -> "Stabilizer"]]];
emit["Method Stabilizer on H,T circuit: head", Head[QuantumCircuitOperator[{"H" -> 1, "T" -> 1}][Method -> "Stabilizer"]]];
emit["Method Stabilizer on a generic unitary circuit: result", QuantumCircuitOperator[{QuantumOperator["RandomUnitary"[2]]}][Method -> "Stabilizer"]];

(* ---------- measurement in the tableau ---------- *)
m1 = ps2[{1}];
emit["ps[{1}] measurement: head, keys, value heads", {Head[m1], Keys[m1], Head /@ Values[m1]}];
m2 = ps0[{1}];
emit["deterministic measurement on |000>: keys", Keys[m2]];
emit["qmo[ps] Z on qubit 1: head, keys", With[{r = QuantumMeasurementOperator["Z", {1}][ps2]}, {Head[r], If[AssociationQ[r], Keys[r], r]}]];
emit["ps expectation ZZI, XXX, ZII", {ps2["Expectation", "ZZI"], ps2["Expectation", "XXX"], ps2["Expectation", "ZII"]}];
emit["inner product <000|GHZ>", ps0["InnerProduct", ps2]];

(* ---------- the non-Clifford boundary: frames ---------- *)
sf = ps2["T", 1];
emit["frame: head, keys, length, qubits", {Head[sf], Keys[First[sf]], sf["Length"], sf["Qubits"]}];
emit["frame state == dense T on GHZ", zeroQ[Normal[sf["StateVector"]] - Normal[QuantumOperator["T", {1}][qs2]["StateVector"]]]];
sf2 = sf["T", 2];
emit["second T doubles the components", sf2["Length"]];
emit["frame + frame, c frame: heads and lengths", {Head[sf + sf2], (sf + sf2)["Length"], Head[2 sf], (2 sf)["Coefficients"]}];
emit["StabilizerStateQ on a frame of length 1 and 2", {StabilizerStateQ[StabilizerFrame[ps2]], StabilizerStateQ[sf]}];
emit["frame expectation ZZI", sf["Expectation", "ZZI"]];
emit["qmo[frame] result head", Head[QuantumMeasurementOperator["Z", {1}][sf]]];

(* ---------- channels on the tableau ---------- *)
bf = QuantumChannel["BitFlip"[1/10]][ps2];
emit["BitFlip[ps]: head, structure", {Head[bf], {#1, Head[#2]} & @@@ bf}];
emit["AmplitudeDamping[ps]: head", Head[QuantumChannel["AmplitudeDamping"[1/10]][ps2]]];
cc = CliffordChannel[ps2];
emit["CliffordChannel[ps]: head, keys", {Head[cc], Keys[First[cc]]}];

(* ---------- graph states ---------- *)
gs = GraphState[CycleGraph[3]];
emit["GraphState: head, stabilizers, PauliStabilizer head", {Head[gs], gs["Stabilizers"], Head[gs["PauliStabilizer"]]}];
emit["GraphState[ps] on GHZ (not graph form)", GraphState[ps2]];
emit["LocalComplement head", Head[LocalComplement[gs, 1]]];

WriteString[$Output, "DONE\n"];
