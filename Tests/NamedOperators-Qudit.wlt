(* Physics/math tests for the qudit (Weyl-Heisenberg) named operators at d=3,5. *)

cm[x_] := Normal[QuantumOperator[x]["Computational"]["Matrix"]]
cmO[x_, o_] := Normal[QuantumOperator[x, o]["Computational"]["Matrix"]]
matEq[a_, b_] := Dimensions[a] === Dimensions[b] && Chop[N[a - b], 10^-10] === ConstantArray[0, Dimensions[a]]
unitaryQ[m_] := matEq[m . ConjugateTranspose[m], IdentityMatrix[Length[m]]]
shiftRef[d_] := Table[If[Mod[i - 1, d] == Mod[j, d], 1, 0], {i, d}, {j, d}];
clockRef[d_] := DiagonalMatrix[Table[Exp[-2 Pi I n / d], {n, 0, d - 1}]];
stdClock[d_] := DiagonalMatrix[Table[Exp[2 Pi I n / d], {n, 0, d - 1}]];
CNOTm = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};

BeginTestSection["Qudit Pauli - shift, inverse clock, and Y"]
(* X_d cyclic shift; Z_d QF inverse clock diag(w^-n); Y_d = -i Z_d X_d. *)
Do[
    VerificationTest[cm["Shift"[d]], shiftRef[d], SameTest -> matEq, TestID -> "Shift-" <> ToString[d]];
    VerificationTest[cm["X"[d]], shiftRef[d], SameTest -> matEq, TestID -> "X-shift-" <> ToString[d]];
    VerificationTest[cm["PauliX"[d]], shiftRef[d], SameTest -> matEq, TestID -> "PauliX-shift-" <> ToString[d]];
    VerificationTest[cm["ShiftPhase"[d]], clockRef[d], SameTest -> matEq, TestID -> "ShiftPhase-" <> ToString[d]];
    VerificationTest[cm["Z"[d]], clockRef[d], SameTest -> matEq, TestID -> "Z-invclock-" <> ToString[d]];
    VerificationTest[cm["PauliZ"[d]], clockRef[d], SameTest -> matEq, TestID -> "PauliZ-invclock-" <> ToString[d]];
    VerificationTest[cm["Y"[d]], -I cm["Z"[d]] . cm["X"[d]], SameTest -> matEq, TestID -> "Y=-iZX-" <> ToString[d]];
    VerificationTest[cm["PauliY"[d]], cm["Y"[d]], SameTest -> matEq, TestID -> "PauliY=Y-" <> ToString[d]];
    VerificationTest[cm["NOT"[d]], cm["X"[d]], SameTest -> matEq, TestID -> "NOT=X-" <> ToString[d]];
    VerificationTest[matEq[MatrixPower[cm["RootNOT"[d]], 2], cm["X"[d]]], True, TestID -> "RootNOT^2=X-" <> ToString[d]];
    VerificationTest[matEq[MatrixPower[cm["X"[d]], d], IdentityMatrix[d]], True, TestID -> "X^d=I-" <> ToString[d]];
    VerificationTest[matEq[MatrixPower[cm["Z"[d]], d], IdentityMatrix[d]], True, TestID -> "Z^d=I-" <> ToString[d]];
    VerificationTest[unitaryQ[cm["X"[d]]], True, TestID -> "X-unitary-" <> ToString[d]];
    VerificationTest[unitaryQ[cm["Y"[d]]], True, TestID -> "Y-unitary-" <> ToString[d]];
    VerificationTest[unitaryQ[cm["Z"[d]]], True, TestID -> "Z-unitary-" <> ToString[d]],
    {d, {3, 5}}
]
EndTestSection[]

BeginTestSection["Qudit Weyl relation X Z = w Z X"]
Do[
    VerificationTest[matEq[cm["X"[d]] . cm["Z"[d]], Exp[2 Pi I / d] cm["Z"[d]] . cm["X"[d]]], True, TestID -> "Weyl-XZ=wZX-" <> ToString[d]],
    {d, {2, 3, 5}}
]
EndTestSection[]

BeginTestSection["HeisenbergWeyl - displacement operators"]
(* HW[d,i,a] = X_d^i Z_std^a; HW[d,0,1] = standard clock, conjugate of QF "Z"[d]. *)
Do[
    VerificationTest[cm["HeisenbergWeyl"[d, 1, 0]], cm["X"[d]], SameTest -> matEq, TestID -> "HW-shift-" <> ToString[d]];
    VerificationTest[cm["HeisenbergWeyl"[d, 0, 1]], stdClock[d], SameTest -> matEq, TestID -> "HW-stdclock-" <> ToString[d]];
    VerificationTest[matEq[cm["HeisenbergWeyl"[d, 0, 1]], cm["Z"[d]]], False, TestID -> "HW-clock-differs-from-Z-" <> ToString[d]];
    VerificationTest[cm["HeisenbergWeyl"[d, 1, 1]], cm["X"[d]] . stdClock[d], SameTest -> matEq, TestID -> "HW-XZ-" <> ToString[d]];
    VerificationTest[unitaryQ[cm["HeisenbergWeyl"[d, 1, 1]]], True, TestID -> "HW-unitary-" <> ToString[d]],
    {d, {3, 5}}
]
EndTestSection[]

BeginTestSection["SUM - qudit generalized CNOT"]
sumRef[d_] := Table[With[{i = IntegerDigits[in - 1, d, 2], o = IntegerDigits[out - 1, d, 2]}, If[i[[1]] == o[[1]] && o[[2]] == Mod[Total[i], d], 1, 0]], {out, d^2}, {in, d^2}];
VerificationTest[cm["SUM"[]], CNOTm, SameTest -> matEq, TestID -> "SUM[2]=CNOT"]
Do[
    VerificationTest[cm["SUM"[d]], sumRef[d], SameTest -> matEq, TestID -> "SUM-" <> ToString[d]];
    VerificationTest[unitaryQ[cm["SUM"[d]]], True, TestID -> "SUM-unitary-" <> ToString[d]],
    {d, {3, 5}}
]
EndTestSection[]

BeginTestSection["Fourier / InverseFourier / Hadamard"]
Do[
    VerificationTest[cm["Fourier"[d]], N[FourierMatrix[d]], SameTest -> matEq, TestID -> "Fourier=DFT-" <> ToString[d]];
    VerificationTest[cm["InverseFourier"[d]], ConjugateTranspose[N[FourierMatrix[d]]], SameTest -> matEq, TestID -> "InvFourier=Fdag-" <> ToString[d]];
    VerificationTest[unitaryQ[cm["Fourier"[d]]], True, TestID -> "Fourier-unitary-" <> ToString[d]];
    VerificationTest[matEq[ConjugateTranspose[N[FourierMatrix[d]]] . cm["X"[d]] . N[FourierMatrix[d]], cm["Z"[d]]], True, TestID -> "FdagXF=Z-" <> ToString[d]];
    VerificationTest[cm["Hadamard"[d]], cm["Fourier"[d]], SameTest -> matEq, TestID -> "Hadamard=Fourier-" <> ToString[d]],
    {d, {2, 3, 5}}
]
VerificationTest[cmO["Fourier"[2], {1, 2}], KroneckerProduct[N[FourierMatrix[2]], N[FourierMatrix[2]]], SameTest -> matEq, TestID -> "Fourier-multi-tensor"]
VerificationTest[cmO["Hadamard"[3], {1, 2}], KroneckerProduct[N[FourierMatrix[3]], N[FourierMatrix[3]]], SameTest -> matEq, TestID -> "Hadamard-multi-tensor"]
EndTestSection[]
