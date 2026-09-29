(* section4-audit-check.wls: the doubling, phase-space, maps-on-density-matrices, and
   measurement-versus-channel facts (now sections 3, 5, and 7), run against the working tree,
   and kernel defect 8 of Appendix D. Each line is the call followed by what the kernel returns. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
out[label_, expr_] := WriteString[$Output, label, ": ", ToString[expr, InputForm], "\n"];

(* D1 "Bend" is the doubling A (x) conj(A), not a map-to-state bending *)
xb = QuantumOperator["X"]["Bend"];
out["D1 X Bend: order, matrix dimensions, equals X (x) conj(X)", {xb["Order"], Dimensions @ Normal @ xb["Matrix"], Normal @ xb["Matrix"] === KroneckerProduct[PauliMatrix[1], Conjugate @ PauliMatrix[1]]}];
out["D1b (-X) Bend equals X Bend", Normal @ QuantumOperator[-PauliMatrix[1]]["Bend"]["Matrix"] === Normal @ xb["Matrix"]];
out["D1c |+> Bend is the density vector", Normal @ QuantumState["+"]["Bend"]["StateVector"] === Flatten @ Normal @ QuantumState["+"]["DensityMatrix"]];
out["D1d the stored state of X: output and input qudits", {QuantumOperator["X"]["State"]["OutputQudits"], QuantumOperator["X"]["State"]["InputQudits"]}];

(* D2 the classical-structure operators and their shapes on a qubit *)
out["D2 matrix dimensions of Copy, Decohere, Discard, Measure, Encode", Map[Dimensions @ Normal @ QuantumOperator[#]["Matrix"] &, {"Copy", "Decohere", "Discard", "Measure", "Encode"}]];
out["D2b Copy, Decohere, Discard matrices", Map[Normal @ QuantumOperator[#]["Matrix"] &, {"Copy", "Decohere", "Discard"}]];
w = QuantumOperator["WSpider"];
out["D2c WSpider: order, m.m^dagger", {w["Order"], Normal[w["Matrix"] . ConjugateTranspose[w["Matrix"]]]}];
out["D2d Multiplexer[] dimensions", Dimensions @ Normal @ QuantumOperator["Multiplexer"[]]["Matrix"]];

(* D3 cup and cap: the unnormalized Bell state and the snake identity *)
out["D3 Cup state vector", Normal @ QuantumOperator["Cup"]["State"]["StateVector"]];
out["D3b snake: Cap after Cup on a spectator wire equals the identity", Normal @ QuantumCircuitOperator[{"Cup"[2] -> {2, 3}, "Cap"[2] -> {1, 2}}]["QuantumOperator"]["Computational"]["Matrix"] === IdentityMatrix[2]];

(* D4 the phase-space bases in even and odd dimension *)
els[d_] := Normal /@ QuantumBasis["Wigner"[d]]["Output"]["Elements"];
gram[d_] := Simplify @ Outer[Tr[ConjugateTranspose[#1] . #2] &, els[d], els[d], 1];
fzx[d_, q_, p_] := MatrixPower[Normal @ QuantumOperator["Fourier"[d]]["Matrix"], 2] . MatrixPower[Normal @ QuantumOperator["Z"[d]]["Matrix"], q] . MatrixPower[Normal @ QuantumOperator["X"[d]]["Matrix"], p];
out["D4 d=4, q=p=1: 2 e^{i pi qp/d} F^2 Z X is an element; 2 e^{2 pi i qp/d} F^2 Z X is not", {MemberQ[Simplify @ els[4], Simplify[2 Exp[I Pi/4] fzx[4, 1, 1]]], MemberQ[Simplify @ els[4], Simplify[2 Exp[2 Pi I/4] fzx[4, 1, 1]]]}];
out["D4b Gram matrices: d=2 is 8 I, d=4 is 16 I, d=3 is 3 I", {gram[2] === 8 IdentityMatrix[4], gram[4] === 16 IdentityMatrix[16], gram[3] === 3 IdentityMatrix[9]}];
out["D4c all elements Hermitian for d=2,3,4", AllTrue[Join[els[2], els[3], els[4]], HermitianMatrixQ]];
out["D4d traces of the elements, d=2 and d=4", {Simplify[Tr /@ els[2]], Simplify[Tr /@ els[4]]}];
out["D4e d=2 elements", Simplify @ els[2]];

(* D5 the inverse takes square roots and returns density matrices *)
qt = QuantumWignerTransform[QuantumState[{1, 0, 0}]];
out["D5 qutrit: dimension after the transform and after the inverse", {qt["Dimension"], QuantumWeylTransform[qt]["Dimension"]}];
inv = QuantumWeylTransform[QuantumWignerTransform[QuantumState[Exp[I a] {1, 1}/Sqrt[2]]]];
out["D5b inverse of the transform of e^{i a}|+>: state type, equals |+><+| for real a", {inv["StateType"], Assuming[Element[a, Reals], FullSimplify[Normal @ inv["DensityMatrix"] - Normal @ QuantumState[{1, 1}/Sqrt[2]]["DensityMatrix"]]] === {{0, 0}, {0, 0}}}];
psi4 = QuantumState[{1, I, -1/2, (1 + I)/3}];
out["D5c d=4 complex state: inverse of the transform equals its density matrix", Simplify[Normal @ QuantumWeylTransform[QuantumWignerTransform[psi4]]["DensityMatrix"] - Normal @ psi4["DensityMatrix"]] === ConstantArray[0, {4, 4}]];

(* D6 the rate matrix is the matrix of the Liouvillian; column sums *)
rate[m_] := Normal @ With[{r = HamiltonianTransitionRate[QuantumOperator[m]]}, If[MatrixQ[r], r, r["Matrix"]]];
wOf[m_] := Normal @ QuantumWignerTransform[QuantumState[m]]["StateVector"];
rho0 = (IdentityMatrix[2] + PauliMatrix[2])/2;
lio = -I h (PauliMatrix[1] . rho0 - rho0 . PauliMatrix[1]);
out["D6 qubit H = h X: R, R.W(rho) == W(-i[H,rho]), == W(-i[H,rho])/i, column sums of R", {rate[h PauliMatrix[1]], Simplify[rate[h PauliMatrix[1]] . wOf[rho0] - wOf[lio]] === {0, 0, 0, 0}, Simplify[rate[h PauliMatrix[1]] . wOf[rho0] - wOf[lio]/I] === {0, 0, 0, 0}, Simplify @ Total @ rate[h PauliMatrix[1]]}];
h3m = Normal[QuantumOperator["X"[3]]["Matrix"]]; h3 = h3m + ConjugateTranspose[h3m];
rho3 = DiagonalMatrix[{1/2, 1/3, 1/6}];
out["D6b qutrit Hermitian H: R.W(rho) == W(-i[H,rho]), column sums", {Simplify[rate[h3] . wOf[rho3] - wOf[-I (h3 . rho3 - rho3 . h3)]] === ConstantArray[0, 9], Simplify @ Total @ rate[h3]}];

(* D7 phase-space evolution: symbolic time fails, an interval works (defect 8) *)
psPlus = QuantumPhaseSpaceTransform[QuantumState[{1, 1}/Sqrt[2]]];
out["D7 QuantumEvolve[Z, phase-space |+>, t]: head", Head @ QuantumEvolve[QuantumOperator["Z"], psPlus, t]];
out["D7b QuantumEvolve[Z, phase-space |+>, {t, 0, 1}]: head, picture", With[{e = QuantumEvolve[QuantumOperator["Z"], psPlus, {t, 0, 1}]}, {Head[e], e["Picture"]}]];
out["D7c the same with a dense state and symbolic t: head", Head @ QuantumEvolve[QuantumOperator["Z"], QuantumState[{1, 1}/Sqrt[2]], t]];

(* D8 the two transforms on a channel, and the wrong-size target basis *)
out["D8 QuantumWignerTransform[chan], QuantumPhaseSpaceTransform[chan]: heads", {Head @ QuantumWignerTransform[QuantumChannel["BitFlip"[1/10]]], Head @ QuantumPhaseSpaceTransform[QuantumChannel["BitFlip"[1/10]]]}];
out["D8b QuditBasis form with a wrong-size target: head and dimension; QuantumBasis form: head", {Head @ QuantumPhaseSpaceTransform[QuditBasis[2], QuantumBasis["Wigner"[3]]], QuantumPhaseSpaceTransform[QuditBasis[2], QuantumBasis["Wigner"[3]]]["Dimension"], Head @ QuantumPhaseSpaceTransform[QuantumBasis[2], QuantumBasis["Wigner"[3]]]}];

(* D9 QuantumOperator represents maps on density matrices on the doubled space (section 3) *)
xm = PauliMatrix[1];
zero00 = {{1, 0}, {0, 0}};
out["D9 Left[X] == X (x) I; Right[X] == I (x) X^T", {Normal @ QuantumOperator["Left"["X"]]["MatrixRepresentation"] === KroneckerProduct[xm, IdentityMatrix[2]], Normal @ QuantumOperator["Right"["X"]]["MatrixRepresentation"] === KroneckerProduct[IdentityMatrix[2], Transpose[xm]]}];
out["D9b Liouvillian[X] == -i(X (x) I - I (x) X^T)", Normal @ QuantumOperator["Liouvillian"["X"]]["MatrixRepresentation"] === -I (KroneckerProduct[xm, IdentityMatrix[2]] - KroneckerProduct[IdentityMatrix[2], Transpose[xm]])];
out["D9c Liouvillian[X] applied to |0><0| == -i[X, |0><0|]", Normal @ QuantumOperator["Liouvillian"["X"]][QuantumState[zero00]]["DensityMatrix"] === -I (xm . zero00 - zero00 . xm)];
out["D9d bit-flip channel Sum M rho M^dag on |0><0| == diag(3/4, 1/4)", Normal @ QuantumChannel[{Sqrt[3/4] IdentityMatrix[2], Sqrt[1/4] xm}][QuantumState[zero00]]["DensityMatrix"] === {{3/4, 0}, {0, 1/4}}];

(* D10 measurement list = effects (square-rooted to measurement operators); channel list = Kraus (section 5) *)
e0 = {{3/4, 0}, {0, 1/4}}; e1 = {{1/4, 0}, {0, 3/4}};
qmoEff = QuantumMeasurementOperator[{e0, e1}];
out["D10 measurement POVMElements are the effects {e0, e1}", (Normal /@ qmoEff["POVMElements"]) === {e0, e1}];
out["D10b measurement operators equal Sqrt[e0], Sqrt[e1]", Simplify[(Normal @ #["MatrixRepresentation"] & /@ Values @ qmoEff["Operators"]) - {MatrixPower[e0, 1/2], MatrixPower[e1, 1/2]}] === {ConstantArray[0, {2, 2}], ConstantArray[0, {2, 2}]}];
out["D10c channel from the same list keeps e0, e1 as Kraus: E(|0><0|) == e0.|0><0|.e0 + e1.|0><0|.e1", Normal @ QuantumChannel[{e0, e1}][QuantumState[zero00]]["DensityMatrix"] === e0 . zero00 . ConjugateTranspose[e0] + e1 . zero00 . ConjugateTranspose[e1]];
