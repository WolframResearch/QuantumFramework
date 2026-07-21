(* Physics/math tests for the parametric named operators:
   XRotation/RX, YRotation/RY, ZRotation/RZ, R, U/U3, U2, U1/Phase/P, PhaseShift,
   GlobalPhase, Diagonal, FlipSign, WignerD.

   Rotations are pinned by their generator and angle convention at a SYMBOLIC
   angle (so the whole one-parameter family is checked, not one sampled value),
   plus unitarity. The U family is pinned against its OpenQASM/IBM closed form.
   The qudit d>2 rotations are locked with a characterization test: they are
   currently non-unitary because the generalized Pauli generators are unitary but
   not Hermitian. *)

cm[x_] := Normal[QuantumOperator[x]["Computational"]["Matrix"]]

matEq[a_, b_] := Dimensions[a] === Dimensions[b] && Chop[N[a - b], 10^-10] === ConstantArray[0, Dimensions[a]]
symEq[a_, b_] := Dimensions[a] === Dimensions[b] && FullSimplify[a - b] === ConstantArray[0, Dimensions[a]]
unitaryQ[m_] := matEq[m . ConjugateTranspose[m], IdentityMatrix[Length[m]]]

ClearAll[t, th, ph, lam];
PX = {{0, 1}, {1, 0}};
PY = {{0, -I}, {I, 0}};
PZ = {{1, 0}, {0, -1}};
jZ[j_] := Normal[QuantumOperator["JZ"[j]]["MatrixRepresentation"]];
jY[j_] := Normal[QuantumOperator["JY"[j]]["MatrixRepresentation"]];


BeginTestSection["Rotations - generator and angle convention (symbolic)"]

(* R_k(t) = exp(-i t P_k / 2), the standard half-angle qubit rotation. *)
VerificationTest[cm["RX"[t]], MatrixExp[-I t / 2 PX], SameTest -> symEq, TestID -> "RX-generator"]
VerificationTest[cm["RY"[t]], MatrixExp[-I t / 2 PY], SameTest -> symEq, TestID -> "RY-generator"]
VerificationTest[cm["RZ"[t]], MatrixExp[-I t / 2 PZ], SameTest -> symEq, TestID -> "RZ-generator"]
VerificationTest[cm["XRotation"[t]], cm["RX"[t]], SameTest -> symEq, TestID -> "XRotation=RX"]
VerificationTest[cm["YRotation"[t]], cm["RY"[t]], SameTest -> symEq, TestID -> "YRotation=RY"]
VerificationTest[cm["ZRotation"[t]], cm["RZ"[t]], SameTest -> symEq, TestID -> "ZRotation=RZ"]

(* Default angle is pi/2. *)
VerificationTest[cm["RX"[]], cm["RX"[Pi / 2]], SameTest -> matEq, TestID -> "RX-default-pi/2"]
VerificationTest[cm["RY"[]], cm["RY"[Pi / 2]], SameTest -> matEq, TestID -> "RY-default-pi/2"]
VerificationTest[cm["RZ"[]], cm["RZ"[Pi / 2]], SameTest -> matEq, TestID -> "RZ-default-pi/2"]

EndTestSection[]


BeginTestSection["Rotations - unitarity at d=2 and d>2 characterization"]

Do[
    VerificationTest[unitaryQ[cm[g[1.3]]], True, TestID -> "Unitary-" <> g <> "-d2"],
    {g, {"RX", "RY", "RZ"}}
]

(* CHARACTERIZATION (present behavior, not an endorsement): the qudit rotations
   R_k(t, d) are built as exp(-i t X_d / 2) from the generalized Pauli X_d, Y_d,
   Z_d, which for d > 2 are unitary but not Hermitian. The exponential of a
   non-Hermitian generator is not unitary, so these "rotations" are non-unitary
   at d > 2 and no message is raised. A conventional qudit rotation would use a
   Hermitian generator (e.g. the spin-j operators "JX"/"JY"/"JZ"). If the qudit
   generator is ever made Hermitian, these three tests flip to True. *)
Do[
    VerificationTest[unitaryQ[cm[g[1.3, 3]]], False, TestID -> "NonUnitary-" <> g <> "-d3"],
    {g, {"RX", "RY", "RZ"}}
]

EndTestSection[]


BeginTestSection["R - multi-Pauli exponential"]

VerificationTest[cm["R"[t, "X"]], MatrixExp[-I t / 2 PX], SameTest -> symEq, TestID -> "R-X=RX"]
VerificationTest[cm["R"[t, "XX"]], MatrixExp[-I t / 2 KroneckerProduct[PX, PX]], SameTest -> symEq, TestID -> "R-XX"]
VerificationTest[cm["R"[t, "ZZ"]], MatrixExp[-I t / 2 KroneckerProduct[PZ, PZ]], SameTest -> symEq, TestID -> "R-ZZ"]
VerificationTest[cm["R"[t, "XY"]], MatrixExp[-I t / 2 KroneckerProduct[PX, PY]], SameTest -> symEq, TestID -> "R-XY"]

EndTestSection[]


BeginTestSection["U family - OpenQASM/IBM closed forms"]

u3[a_, b_, c_] := {{Cos[a / 2], -Exp[I c] Sin[a / 2]}, {Exp[I b] Sin[a / 2], Exp[I (b + c)] Cos[a / 2]}};
VerificationTest[cm["U3"[th, ph, lam]], u3[th, ph, lam], SameTest -> symEq, TestID -> "U3-closed-form"]
VerificationTest[cm["U"[th, ph, lam]], u3[th, ph, lam], SameTest -> symEq, TestID -> "U=U3"]
VerificationTest[cm["U2"[ph, lam]], {{1, -Exp[I lam]}, {Exp[I ph], Exp[I (ph + lam)]}} / Sqrt[2], SameTest -> symEq, TestID -> "U2-closed-form"]
VerificationTest[cm["U1"[t]], {{1, 0}, {0, Exp[I t]}}, SameTest -> symEq, TestID -> "U1=phase"]
VerificationTest[cm["Phase"[t]], {{1, 0}, {0, Exp[I t]}}, SameTest -> symEq, TestID -> "Phase-closed-form"]
VerificationTest[cm["P"[t]], {{1, 0}, {0, Exp[I t]}}, SameTest -> symEq, TestID -> "P-closed-form"]
(* U3[] defaults to U(0, pi, pi) = Identity. *)
VerificationTest[cm["U3"[]], IdentityMatrix[2], SameTest -> matEq, TestID -> "U3-default=I"]
(* Phase[]/P[] default to pi, giving Z. *)
VerificationTest[cm["Phase"[]], PZ, SameTest -> matEq, TestID -> "Phase-default=Z"]

(* Qudit phase gate: a single phase on the top basis state. *)
VerificationTest[cm["Phase"[t, 3]], DiagonalMatrix[{1, 1, Exp[I t]}], SameTest -> symEq, TestID -> "Phase-d3"]

EndTestSection[]


BeginTestSection["PhaseShift and GlobalPhase"]

(* PhaseShift[k] = diag(1, exp(i sign(k) 2 pi / 2^|k|)), the QFT R_k rotation;
   negative k gives the inverse. *)
VerificationTest[cm["PhaseShift"[1]], PZ, SameTest -> matEq, TestID -> "PhaseShift1=Z"]
VerificationTest[cm["PhaseShift"[2]], {{1, 0}, {0, I}}, SameTest -> matEq, TestID -> "PhaseShift2=S"]
VerificationTest[cm["PhaseShift"[3]], {{1, 0}, {0, Exp[I Pi / 4]}}, SameTest -> matEq, TestID -> "PhaseShift3=T"]
VerificationTest[cm["PhaseShift"[-2]], {{1, 0}, {0, -I}}, SameTest -> matEq, TestID -> "PhaseShift-2=Sdag"]
Do[
    VerificationTest[cm["PhaseShift"[k]], DiagonalMatrix[{1, Exp[I 2 Pi / 2^k]}], SameTest -> matEq, TestID -> "PhaseShift-Rk-" <> ToString[k]],
    {k, {1, 2, 3, 4}}
]

VerificationTest[cm["GlobalPhase"[t, 3]], Exp[I t] IdentityMatrix[3], SameTest -> symEq, TestID -> "GlobalPhase-d3"]
VerificationTest[cm["GlobalPhase"[t]], {{Exp[I t]}}, SameTest -> symEq, TestID -> "GlobalPhase-default-d1"]

EndTestSection[]


BeginTestSection["Diagonal and FlipSign"]

VerificationTest[cm["Diagonal"[{2, 3}]], DiagonalMatrix[{2, 3}], SameTest -> matEq, TestID -> "Diagonal-list"]
VerificationTest[cm["Diagonal"[t]], t IdentityMatrix[2], SameTest -> symEq, TestID -> "Diagonal-scalar=cI"]
VerificationTest[cm["Diagonal"[{1, -1, -1, 1}]], DiagonalMatrix[{1, -1, -1, 1}], SameTest -> matEq, TestID -> "Diagonal-2qubit"]

(* FlipSign[digits, d] places -1 at the basis state with the given digits. *)
VerificationTest[cm["FlipSign"[{1}]], DiagonalMatrix[{1, -1}], SameTest -> matEq, TestID -> "FlipSign-1"]
VerificationTest[cm["FlipSign"[{0}]], DiagonalMatrix[{-1, 1}], SameTest -> matEq, TestID -> "FlipSign-0"]
VerificationTest[cm["FlipSign"[{1, 1}]], DiagonalMatrix[{1, 1, 1, -1}], SameTest -> matEq, TestID -> "FlipSign-11=CZ"]
VerificationTest[cm["FlipSign"[{1, 0}]], DiagonalMatrix[{1, 1, -1, 1}], SameTest -> matEq, TestID -> "FlipSign-10"]
VerificationTest[cm["FlipSign"["11"]], cm["FlipSign"[{1, 1}]], SameTest -> matEq, TestID -> "FlipSign-string"]
(* Qutrit target. *)
VerificationTest[cm["FlipSign"[{2}, 3]], DiagonalMatrix[{1, 1, -1}], SameTest -> matEq, TestID -> "FlipSign-qutrit"]

EndTestSection[]


BeginTestSection["WignerD - Euler-angle spin-j rotation"]

(* WignerD[j, {a, b, c}] = exp(-i a Jz) exp(-i b Jy) exp(-i c Jz), the active
   Euler-angle rotation on the spin-j representation (descending m, |0> = m=+j).
   Reference built from QF's own "JZ"/"JY" spin matrices. *)
wd[j_, a_, b_, c_] := MatrixExp[-I a jZ[j]] . MatrixExp[-I b jY[j]] . MatrixExp[-I c jZ[j]];
{aa, bb, cc} = {2/5, 7/10, 11/10};
Do[
    VerificationTest[
        Normal[QuantumOperator["WignerD"[j, {aa, bb, cc}]]["Computational"]["Matrix"]],
        wd[j, aa, bb, cc],
        SameTest -> matEq, TestID -> "WignerD-Rz-Ry-Rz-" <> ToString[j]
    ];
    VerificationTest[
        Normal[QuantumOperator["WignerD"[j, bb]]["Computational"]["Matrix"]],
        MatrixExp[-I bb jY[j]],
        SameTest -> matEq, TestID -> "WignerD-single=Ry-" <> ToString[j]
    ];
    VerificationTest[
        unitaryQ[Normal[QuantumOperator["WignerD"[j, {aa, bb, cc}]]["Computational"]["Matrix"]]],
        True, TestID -> "WignerD-unitary-" <> ToString[j]
    ],
    {j, {1/2, 1, 3/2}}
]

(* Default j = 1/2. *)
VerificationTest[
    Normal[QuantumOperator["WignerD"[{aa, bb, cc}]]["Computational"]["Matrix"]],
    wd[1/2, aa, bb, cc],
    SameTest -> matEq, TestID -> "WignerD-default-j"
]

EndTestSection[]
