(* Physics/math tests for the single-qubit and Pauli-string named operators:
   X, Y, Z, PauliX/Y/Z, Shift, ShiftPhase, H, Hadamard, NOT, RootNOT, S, T, V,
   SX, "0", "1", and the case-insensitive aliases.

   These assert matrix CONTENT and the defining gate identities, not just
   dimensions. A shape-only check (or a serialize/reconstruct roundtrip) is blind
   to a wrong matrix of the right shape: that is the class of gap that let the JX
   diagonal-Jz bug survive. Every gate here is pinned by its closed form and by
   an algebraic relation that fails for any other 2x2 of the same shape. *)

cm[x_] := Normal[QuantumOperator[x]["Computational"]["Matrix"]]

(* numeric matrix equality, tolerant of floating dust and radicals/exponentials *)
matEq[a_, b_] := Dimensions[a] === Dimensions[b] && Chop[N[a - b], 10^-10] === ConstantArray[0, Dimensions[a]]
unitaryQ[m_] := matEq[m . ConjugateTranspose[m], IdentityMatrix[Length[m]]]

PX = {{0, 1}, {1, 0}};
PY = {{0, -I}, {I, 0}};
PZ = {{1, 0}, {0, -1}};
Hm = {{1, 1}, {1, -1}} / Sqrt[2];
Sm = {{1, 0}, {0, I}};
Tm = {{1, 0}, {0, Exp[I Pi / 4]}};
Vm = {{1 + I, 1 - I}, {1 - I, 1 + I}} / 2;
I2 = IdentityMatrix[2];


BeginTestSection["Pauli - qubit closed forms"]

VerificationTest[cm["X"], PX, SameTest -> matEq, TestID -> "X-matrix"]
VerificationTest[cm["Y"], PY, SameTest -> matEq, TestID -> "Y-matrix"]
VerificationTest[cm["Z"], PZ, SameTest -> matEq, TestID -> "Z-matrix"]
VerificationTest[cm["PauliX"], PX, SameTest -> matEq, TestID -> "PauliX-matrix"]
VerificationTest[cm["PauliY"], PY, SameTest -> matEq, TestID -> "PauliY-matrix"]
VerificationTest[cm["PauliZ"], PZ, SameTest -> matEq, TestID -> "PauliZ-matrix"]
VerificationTest[cm["Shift"], PX, SameTest -> matEq, TestID -> "Shift=X"]
VerificationTest[cm["ShiftPhase"], PZ, SameTest -> matEq, TestID -> "ShiftPhase=Z"]
VerificationTest[cm["H"], Hm, SameTest -> matEq, TestID -> "H-matrix"]
VerificationTest[cm["Hadamard"], Hm, SameTest -> matEq, TestID -> "Hadamard-matrix"]
VerificationTest[cm["NOT"], PX, SameTest -> matEq, TestID -> "NOT=X"]
VerificationTest[cm["S"[]], Sm, SameTest -> matEq, TestID -> "S-matrix"]
VerificationTest[cm["T"[]], Tm, SameTest -> matEq, TestID -> "T-matrix"]
VerificationTest[cm["V"[]], Vm, SameTest -> matEq, TestID -> "V-matrix"]

(* "0" and "1" are the sign-flip (phase-oracle) operators on the named basis
   state, not the state-injection kets that the same strings denote in circuit
   shorthand. As named operators they are -Z and Z. *)
VerificationTest[cm["0"[]], -PZ, SameTest -> matEq, TestID -> "0=-Z"]
VerificationTest[cm["1"[]], PZ, SameTest -> matEq, TestID -> "1=Z"]

EndTestSection[]


BeginTestSection["Pauli - Hermiticity and unitarity"]

Do[
    VerificationTest[HermitianMatrixQ[cm[g]], True, TestID -> "Herm-" <> g],
    {g, {"X", "Y", "Z"}}
]

Do[
    VerificationTest[unitaryQ[cm[g]], True, TestID -> "Unitary-" <> StringReplace[g, {"[" -> "", "]" -> ""}]],
    {g, {"X", "Y", "Z", "H", "NOT"}}
]
Do[
    VerificationTest[unitaryQ[cm[g[]]], True, TestID -> "Unitary-br-" <> g],
    {g, {"S", "T", "V", "RootNOT"}}
]

EndTestSection[]


BeginTestSection["Pauli - defining algebraic relations"]

VerificationTest[matEq[cm["X"] . cm["X"], I2], True, TestID -> "X^2=I"]
VerificationTest[matEq[cm["Y"] . cm["Y"], I2], True, TestID -> "Y^2=I"]
VerificationTest[matEq[cm["Z"] . cm["Z"], I2], True, TestID -> "Z^2=I"]

(* Pauli multiplication table: X Y = i Z and cyclic. *)
VerificationTest[matEq[cm["X"] . cm["Y"], I cm["Z"]], True, TestID -> "XY=iZ"]
VerificationTest[matEq[cm["Y"] . cm["Z"], I cm["X"]], True, TestID -> "YZ=iX"]
VerificationTest[matEq[cm["Z"] . cm["X"], I cm["Y"]], True, TestID -> "ZX=iY"]

(* Hadamard conjugation swaps X and Z. *)
VerificationTest[matEq[cm["H"] . cm["X"] . cm["H"], cm["Z"]], True, TestID -> "HXH=Z"]
VerificationTest[matEq[cm["H"] . cm["Z"] . cm["H"], cm["X"]], True, TestID -> "HZH=X"]
VerificationTest[matEq[cm["H"] . cm["H"], I2], True, TestID -> "H^2=I"]

(* Phase-gate tower: S^2 = Z, T^2 = S, T^4 = Z. *)
VerificationTest[matEq[cm["S"[]] . cm["S"[]], PZ], True, TestID -> "S^2=Z"]
VerificationTest[matEq[MatrixPower[cm["T"[]], 2], cm["S"[]]], True, TestID -> "T^2=S"]
VerificationTest[matEq[MatrixPower[cm["T"[]], 4], PZ], True, TestID -> "T^4=Z"]

(* V and RootNOT are the principal square root of X (the IBM SX gate). *)
VerificationTest[matEq[cm["V"[]] . cm["V"[]], PX], True, TestID -> "V^2=X"]
VerificationTest[matEq[MatrixPower[cm["RootNOT"[]], 2], PX], True, TestID -> "RootNOT^2=X"]
VerificationTest[matEq[cm["RootNOT"[]], Vm], True, TestID -> "RootNOT=V"]
VerificationTest[matEq[cm["SX"[]], Vm], True, TestID -> "SXbracket=V"]

EndTestSection[]


BeginTestSection["Pauli - case-insensitive aliases"]

(* The registry keeps a case-folding layer so a lowercase spelling dispatches to
   the registered gate. These guard that layer: if the case-insensitive lookup
   regresses to a dead branch, each lowercase name fails to construct and the
   matrix comparison errors out. *)
VerificationTest[cm["x"], PX, SameTest -> matEq, TestID -> "alias-x"]
VerificationTest[cm["y"], PY, SameTest -> matEq, TestID -> "alias-y"]
VerificationTest[cm["z"], PZ, SameTest -> matEq, TestID -> "alias-z"]
VerificationTest[cm["h"], Hm, SameTest -> matEq, TestID -> "alias-h"]
VerificationTest[cm["not"], PX, SameTest -> matEq, TestID -> "alias-not"]
VerificationTest[cm["hadamard"], Hm, SameTest -> matEq, TestID -> "alias-hadamard"]

EndTestSection[]


BeginTestSection["Pauli - SX bare name is two-qubit S tensor X (characterization)"]

(* "SX" is the only registered name whose characters all lie in the Pauli-string
   alphabet {I,X,Y,Z,H,S,T,V,P}, so the bare string is parsed by the chain rule
   as the two-qubit product S tensor X, not the single-qubit sqrt(X). The
   single-qubit gate is the bracketed "SX"[] or "V". *)
VerificationTest[QuantumOperator["SX"]["Arity"], 2, TestID -> "SX-arity-2"]
VerificationTest[cm["SX"], KroneckerProduct[Sm, PX], SameTest -> matEq, TestID -> "SX=S(x)X"]
VerificationTest[QuantumOperator["SX"[]]["Arity"], 1, TestID -> "SXbracket-arity-1"]

EndTestSection[]


BeginTestSection["Pauli strings - chain rule"]

(* A multi-character string over the Pauli-string alphabet is the tensor product
   of the per-character single-qubit operators. *)
VerificationTest[cm["XZ"], KroneckerProduct[PX, PZ], SameTest -> matEq, TestID -> "XZ=X(x)Z"]
VerificationTest[cm["IZ"], KroneckerProduct[I2, PZ], SameTest -> matEq, TestID -> "IZ=I(x)Z"]
VerificationTest[cm["XYZ"], KroneckerProduct[PX, PY, PZ], SameTest -> matEq, TestID -> "XYZ=X(x)Y(x)Z"]

EndTestSection[]
