PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

emit[label_, val_] := WriteString[$Output, label, " -> ", ToString[val, InputForm], "\n"];
zeroQ[m_] := Max[Abs[Flatten[N[Normal[m]]]]] < 10^-9;
permutationMatrixQ[m_] := With[{n = N[Normal[m]]},
    zeroQ[n (1 - n)] && zeroQ[Total[n] - 1] && zeroQ[Total[n, {2}] - 1]
];
vec[qs_] := Chop[N[Normal[qs["StateVector"]]]];
mat[qo_] := Chop[N[Normal[qo["Matrix"]]]];
stateQ = Wolfram`QuantumFramework`PackageScope`QuantumStateQ;

(* ---------- odd d: the Wigner basis in d = 3 ---------- *)
wb = QuditBasis["Wigner"[3]];
els = Normal /@ wb["Elements"];
emit["Wigner[3] element count, dimension", {Length[els], wb["Dimension"]}];
emit["elements are 3x3 Hermitian with trace 1", {AllTrue[els, Dimensions[#] == {3, 3} &], AllTrue[els, HermitianMatrixQ[Chop[N[#]]] &], Chop[N[Tr /@ els]]}];
emit["Gram matrix Tr[A_a A_b] == d Identity", zeroQ[Outer[Tr[#1 . #2] &, els, els, 1] - 3 IdentityMatrix[9]]];
emit["sum of elements == d Identity", zeroQ[Total[els] - 3 IdentityMatrix[3]]];

(* ---------- states ---------- *)
q0 = QuantumState["0", 3];
w0 = QuantumPhaseSpaceTransform[q0];
emit["W(|0>) picture, dimensions, qudits", {w0["Picture"], w0["Dimensions"], w0["Qudits"]}];
emit["W(|0>) vector", vec[w0]];
emit["W(|0>) == Tr[A_a rho]/d", zeroQ[vec[w0] - Chop[N[Tr[# . Normal[q0["DensityMatrix"]]] & /@ els]/3]]];
emit["W(|0>) grid", Chop[N[Normal[w0["PhaseSpace"]]]]];
strange = QuantumState[{0, 1, -1}/Sqrt[2], 3];
ws = QuantumPhaseSpaceTransform[strange];
emit["W(strange) min entry, total, sum|W|", {Min[vec[ws]], Chop[Total[vec[ws]]], Total[Abs[vec[ws]]]}];

(* ---------- operators: functoriality, equivariance, Cliffords, inverse ---------- *)
F3 = QuantumOperator["Fourier"[3]]; Z3 = QuantumOperator["Z"[3]]; X3 = QuantumOperator["X"[3]];
S3 = QuantumOperator[DiagonalMatrix[{1, 1, Exp[2 Pi I/3]}], 3, "Label" -> "S3"];
R3 = QuantumOperator["RandomUnitary"[3]];
{WF, WZ, WX, WS, WR} = QuantumPhaseSpaceTransform /@ {F3, Z3, X3, S3, R3};
emit["W(Fourier) dims, real", {Dimensions[mat[WF]], zeroQ[Im[mat[WF]]]}];
emit["Cliffords F, X, Z, S3 are permutations", permutationMatrixQ /@ {mat[WF], mat[WX], mat[WZ], mat[WS]}];
emit["random unitary: permutation, real, column sums 1", {permutationMatrixQ[mat[WR]], zeroQ[Im[mat[WR]]], zeroQ[Total[mat[WR]] - 1]}];
emit["functoriality W(F.S) == W(F).W(S)", zeroQ[mat[QuantumPhaseSpaceTransform[F3[S3]]] - mat[WF] . mat[WS]]];
emit["functoriality W(R.F) == W(R).W(F)", zeroQ[mat[QuantumPhaseSpaceTransform[R3[F3]]] - mat[WR] . mat[WF]]];
emit["equivariance W(R rho) == W(R).W(rho)", zeroQ[vec[QuantumPhaseSpaceTransform[R3[strange]]] - mat[WR] . vec[ws]]];
emit["equivariance through the heads W(R)[W(rho)]", zeroQ[vec[WR[ws]] - vec[QuantumPhaseSpaceTransform[R3[strange]]]]];
emit["inverse on states: Weyl[W(rho)] == rho", QuantumWeylTransform[ws] == strange];
emit["inverse on operators: Weyl[W(F)] dims", Dimensions[mat[QuantumWeylTransform[WF]]]];
emit["inverse on operators lands on the doubled operator F (x) conj F", zeroQ[mat[QuantumWeylTransform[WF]] - mat[F3["Double"]]]];
emit["W kills the global phase: W(I F) == W(F)", zeroQ[mat[QuantumPhaseSpaceTransform[I F3]] - mat[WF]]];

(* ---------- measurement ---------- *)
qmo = QuantumMeasurementOperator[{1}, 3];
born = Chop[N[Values[qmo[strange]["Probabilities"]]]];
Wqmo = QuantumPhaseSpaceTransform[qmo];
qmW = Wqmo[ws];
emit["Born probabilities of the strange state", born];
emit["W(qmo): head, POVM dims, targets", {Head[Wqmo], Dimensions[mat[Wqmo["POVM"]]], Wqmo["Targets"]}];
emit["W(qmo)[W(rho)]: head, outcome count", {Head[qmW], Length[qmW["Outcomes"]]}];
emit["W(qmo)[W(rho)] outcomes", qmW["Outcomes"]];
probs9 = Chop[N[Values[qmW["Probabilities"]]]];
emit["W(qmo)[W(rho)] probabilities", probs9];
pdiag = QuantumState[DiagonalMatrix[{0, 1/2, 1/2}], 3];
emit["W(diag Born) vector", vec[QuantumPhaseSpaceTransform[pdiag]]];
emit["probabilities == W(diag Born)", zeroQ[probs9 - vec[QuantumPhaseSpaceTransform[pdiag]]]];
emit["probabilities == Tr[A_a diag Born] (no 1/d)", zeroQ[probs9 - Chop[N[Tr[# . DiagonalMatrix[{0, 1/2, 1/2}]] & /@ els]]]];
emit["probabilities == coordinates of the doubled register marginal", zeroQ[probs9 - vec[QuantumPhaseSpaceTransform[QuantumPartialTrace[qmo[strange]["State"], {2}]]]]];

(* ---------- SIC dual frame: coordinates are POVM probabilities ---------- *)
SeedRandom[1]; rho2 = QuantumState["RandomMixed"];
sic = QuantumPhaseSpaceTransform[rho2, "QBismSIC"];
emit["SIC coordinates, nonnegative, sum", {vec[sic], Min[vec[sic]] >= 0, Chop[Total[vec[sic]]]}];
emit["SIC coordinates == QBismSICPOVM probabilities", zeroQ[vec[sic] - Values[QuantumMeasurementOperator["QBismSICPOVM"][rho2]["Probabilities"]]]];
emit["SIC transform of H: real, min entry", {zeroQ[Im[mat[QuantumPhaseSpaceTransform[QuantumOperator["H"], "QBismSIC"]]]], Min[mat[QuantumPhaseSpaceTransform[QuantumOperator["H"], "QBismSIC"]]]}];

(* ---------- positive transform: the sign-sector construction on tensors ---------- *)
pv = Normal[QuantumPositiveTransform[vec[ws]]];
emit["Positive(W): dims (sector, point), min >= 0", {Dimensions[pv], Min[pv] >= 0}];
emit["Positive(W): plus - minus == W, total == sum|W|", {zeroQ[pv[[1]] - pv[[2]] - vec[ws]], Chop[Total[pv, 2] - Total[Abs[vec[ws]]]] == 0}];
pm = Normal[QuantumPositiveTransform[mat[WR]]];
emit["Positive(W(R)): dims (sector, point, sector, point), min >= 0", {Dimensions[pm], Min[pm] >= 0}];
emit["Positive(W(R)) blocks are T+, T-, T-, T+", {zeroQ[pm[[1, All, 1, All]] - Ramp[mat[WR]]], zeroQ[pm[[1, All, 2, All]] - Ramp[-mat[WR]]], zeroQ[pm[[2, All, 1, All]] - Ramp[-mat[WR]]], zeroQ[pm[[2, All, 2, All]] - Ramp[mat[WR]]]}];
out = ArrayReshape[pm, {18, 18}] . Flatten[pv];
emit["sectors of Positive(W(R)).Positive(W) recover W(R).W", zeroQ[out[[;; 9]] - out[[10 ;;]] - mat[WR] . vec[ws]]];
(* measurement at the level of the map: transform, apply, invert *)
jointBack = QuantumWeylTransform[qmW["State"]];
jointDense = qmo[strange]["State"];
emit["measurement map: dims after inverse vs dense post-measurement joint", {jointBack["Dimensions"], jointDense["Dimensions"]}];
emit["measurement map is equivariant (inverse of transformed joint == dense joint)", zeroQ[Normal[jointBack["DensityMatrix"]] - Normal[jointDense["DensityMatrix"]]]];
emit["QuditBasis has a 'Canonical' property", MemberQ[QuditBasis[2]["Properties"], "Canonical"]];
posState = QuantumPositiveTransform[ws];
emit["QuantumPositiveTransform on a state: head, valid state", {Head[posState], stateQ[posState]}];

(* ---------- generator: the Liouvillian as a rate matrix ---------- *)
Hm = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}};
Hop = QuantumOperator[Hm, 3, "Label" -> "H"];
Rt = HamiltonianTransitionRate[Hop];
emit["rate matrix dims, real, column sums 0", {Dimensions[Rt], zeroQ[Im[N[Rt]]], zeroQ[Total[N[Rt]]]}];
t0 = 7/10;
rhot = QuantumState[MatrixExp[-I Hm t0] . Normal[strange["DensityMatrix"]] . MatrixExp[I Hm t0], 3];
emit["exp(R t).W(0) == W(rho(t))", zeroQ[MatrixExp[N[Rt] t0] . vec[ws] - vec[QuantumPhaseSpaceTransform[rhot]]]];

(* ---------- channel on a qubit (even d) ---------- *)
chan = QuantumChannel["BitFlip"[1/10]];
Wchan = QuantumPhaseSpaceTransform[chan];
emit["W(channel): head, dims, real", {Head[Wchan], Dimensions[mat[Wchan]], zeroQ[Im[mat[Wchan]]]}];
emit["W(channel).W(rho) == W(channel[rho])", zeroQ[mat[Wchan] . vec[QuantumPhaseSpaceTransform[rho2]] - vec[QuantumPhaseSpaceTransform[chan[rho2]]]]];

(* ---------- even d: the qubit convention ---------- *)
wb2 = QuditBasis["Wigner"[2]]; els2 = Normal /@ wb2["Elements"];
emit["Wigner[2] elements", Chop[N[els2]]];
emit["Wigner[2] traces, Gram matrix, sum", {Chop[N[Tr /@ els2]], Chop[N[Outer[Tr[#1 . #2] &, els2, els2, 1]]], Chop[N[Total[els2]]]}];
w0q = QuantumPhaseSpaceTransform[QuantumState["0"]];
emit["qubit W(|0>) vector, total", {vec[w0q], Total[vec[w0q]]}];
emit["qubit W(|0>) == Tr[A_a rho]/d ?", zeroQ[vec[w0q] - Chop[N[Tr[# . {{1, 0}, {0, 0}}] & /@ els2]/2]]];
grid = Chop[N[Normal[w0q["PhaseSpace"]]]];
emit["qubit grid 2d x 2d, total, Total[grid][[;;;;2]]", {grid, Total[Flatten[grid]], Total[grid][[;; ;; 2]]}];
wH = QuantumPhaseSpaceTransform[QuantumOperator["H"]];
emit["qubit W(H): matrix, column sums", {mat[wH], Total[mat[wH]]}];
emit["qubit equivariance W(H).W(|0>) == W(H|0>)", zeroQ[mat[wH] . vec[w0q] - vec[QuantumPhaseSpaceTransform[QuantumOperator["H"][QuantumState["0"]]]]]];
emit["qubit inverse Weyl[W(|0>)] == |0>", QuantumWeylTransform[w0q] == QuantumState["0"]];
emit["W(channel) column sums (even d)", Total[mat[Wchan]]];

(* ---------- tensor power of a one-qudit basis; global-qudit QuasiProbability ---------- *)
pair = QuantumPhaseSpaceTransform[QuantumTensorProduct[rho2, rho2], "QBismSIC"];
emit["2-qubit SIC: dims, product of single-qubit coordinates", {pair["Dimensions"], zeroQ[vec[pair] - Flatten[KroneckerProduct[vec[sic], vec[sic]]]]}];
bell = QuantumState["PhiPlus"];
emit["QuasiProbability dims (one d=4 qudit) vs head-wise W(bell) grid dims", {Dimensions[Normal[bell["QuasiProbability"]]], Dimensions[Normal[QuantumPhaseSpaceTransform[bell]["PhaseSpace"]]]}];

(* ---------- wrong-dimension bases ---------- *)
emit["qutrit with the qubit SIC basis", Head[QuantumPhaseSpaceTransform[q0, "QBismSIC"]]];
emit["qutrit with QuditBasis['QBismSIC']", Head[QuantumPhaseSpaceTransform[q0, QuditBasis["QBismSIC"]]]];

WriteString[$Output, "DONE\n"];
