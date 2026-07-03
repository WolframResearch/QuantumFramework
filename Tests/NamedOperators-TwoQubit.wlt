(* Physics/math tests for the multi-qubit named operators:
   the controlled family (C/Controlled, C0/Controlled0, CX, CY, CZ, CH, CS, CT,
   CPHASE, CNOT, C0NOT), SWAP, RootSWAP, CSWAP/Fredkin, C0SWAP, Toffoli, Deutsch,
   Multiplexer/BlockDiagonal, and Braid.

   Each gate is pinned by its explicit computational matrix or truth table and,
   where the semantics matter, by the control-subspace structure. The general
   Controlled/Controlled0 machinery is checked against the standard CNOT/C0NOT so
   that a wiring or block-placement regression is caught, not only a dimension
   change. *)

cm[x_] := Normal[QuantumOperator[x]["Computational"]["Matrix"]]
cmO[x_, o_] := Normal[QuantumOperator[x, o]["Computational"]["Matrix"]]

matEq[a_, b_] := Dimensions[a] === Dimensions[b] && Chop[N[a - b], 10^-10] === ConstantArray[0, Dimensions[a]]
symEq[a_, b_] := Dimensions[a] === Dimensions[b] && FullSimplify[a - b] === ConstantArray[0, Dimensions[a]]
unitaryQ[m_] := matEq[m . ConjugateTranspose[m], IdentityMatrix[Length[m]]]

ClearAll[t];
I2 = IdentityMatrix[2];
Hm = {{1, 1}, {1, -1}} / Sqrt[2];
Sm = {{1, 0}, {0, I}};
PX = {{0, 1}, {1, 0}};
CNOTm = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
C0NOTm = {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
SWAPm = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};
(* 3-qubit permutation-of-basis references, qubit 1 the high bit *)
toffM = ReplacePart[IdentityMatrix[8], {{7, 7} -> 0, {8, 8} -> 0, {7, 8} -> 1, {8, 7} -> 1}];
fredM = ReplacePart[IdentityMatrix[8], {{6, 6} -> 0, {7, 7} -> 0, {6, 7} -> 1, {7, 6} -> 1}];
c0swapM = ReplacePart[IdentityMatrix[8], {{2, 2} -> 0, {3, 3} -> 0, {2, 3} -> 1, {3, 2} -> 1}];


BeginTestSection["Controlled - named gates, explicit matrices"]

VerificationTest[cm["CNOT"], CNOTm, SameTest -> matEq, TestID -> "CNOT-matrix"]
VerificationTest[cm["CX"], CNOTm, SameTest -> matEq, TestID -> "CX=CNOT"]
VerificationTest[cm["CZ"], DiagonalMatrix[{1, 1, 1, -1}], SameTest -> matEq, TestID -> "CZ-matrix"]
VerificationTest[cm["CY"], {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -I}, {0, 0, I, 0}}, SameTest -> matEq, TestID -> "CY-matrix"]
VerificationTest[cm["CH"], ArrayFlatten[{{I2, 0}, {0, Hm}}], SameTest -> matEq, TestID -> "CH-matrix"]
VerificationTest[cm["CS"], DiagonalMatrix[{1, 1, 1, I}], SameTest -> matEq, TestID -> "CS-matrix"]
VerificationTest[cm["CT"], DiagonalMatrix[{1, 1, 1, Exp[I Pi / 4]}], SameTest -> matEq, TestID -> "CT-matrix"]
VerificationTest[cm["CPHASE"[t]], DiagonalMatrix[{1, 1, 1, Exp[I t]}], SameTest -> symEq, TestID -> "CPHASE-matrix"]
VerificationTest[cm["C0NOT"], C0NOTm, SameTest -> matEq, TestID -> "C0NOT-matrix"]

VerificationTest[matEq[cm["CNOT"] . cm["CNOT"], IdentityMatrix[4]], True, TestID -> "CNOT^2=I"]
Do[
    VerificationTest[unitaryQ[cm[g]], True, TestID -> "Unitary-" <> g],
    {g, {"CNOT", "CX", "CY", "CZ", "CH", "CS", "CT", "C0NOT"}}
]

EndTestSection[]


BeginTestSection["Controlled - general machinery reproduces the standard gates"]

(* Controlled[op -> target] builds the control on the complementary wire with the
   op acting in the |1> control subspace; Controlled0 uses the |0> subspace. *)
VerificationTest[cm["Controlled"["PauliX" -> {2}]], CNOTm, SameTest -> matEq, TestID -> "Controlled-PauliX=CNOT"]
VerificationTest[cm["Controlled0"["PauliX" -> {2}]], C0NOTm, SameTest -> matEq, TestID -> "Controlled0-PauliX=C0NOT"]
VerificationTest[cm["Controlled"["H" -> {2}]], cm["CH"], SameTest -> matEq, TestID -> "Controlled-H=CH"]
VerificationTest[cm["Controlled"["S"[] -> {2}]], cm["CS"], SameTest -> matEq, TestID -> "Controlled-S=CS"]
VerificationTest[cm["Controlled"["Phase"[t] -> {2}]], cm["CPHASE"[t]], SameTest -> symEq, TestID -> "Controlled-Phase=CPHASE"]
(* Two controls: op fires only when both control qubits are |1>. *)
VerificationTest[cmO["Controlled"["NOT", {1, 2}], {3}], toffM, SameTest -> matEq, TestID -> "Controlled-2ctrl=Toffoli"]

EndTestSection[]


BeginTestSection["SWAP family"]

VerificationTest[cm["SWAP"], SWAPm, SameTest -> matEq, TestID -> "SWAP-matrix"]
VerificationTest[unitaryQ[cm["SWAP"]], True, TestID -> "SWAP-unitary"]
VerificationTest[matEq[cm["SWAP"] . cm["SWAP"], IdentityMatrix[4]], True, TestID -> "SWAP^2=I"]
VerificationTest[matEq[MatrixPower[cm["RootSWAP"[]], 2], SWAPm], True, TestID -> "RootSWAP^2=SWAP"]
VerificationTest[unitaryQ[cm["RootSWAP"[]]], True, TestID -> "RootSWAP-unitary"]

EndTestSection[]


BeginTestSection["Toffoli / Fredkin / C0SWAP truth tables"]

VerificationTest[cmO["Toffoli", {1, 2, 3}], toffM, SameTest -> matEq, TestID -> "Toffoli-truth-table"]
VerificationTest[cm["Fredkin"[]], fredM, SameTest -> matEq, TestID -> "Fredkin-truth-table"]
VerificationTest[cm["CSWAP"[]], fredM, SameTest -> matEq, TestID -> "CSWAP=Fredkin"]
VerificationTest[cm["C0SWAP"[]], c0swapM, SameTest -> matEq, TestID -> "C0SWAP-truth-table"]
Do[
    VerificationTest[unitaryQ[cmO["Toffoli", {1, 2, 3}]], True, TestID -> "Toffoli-unitary"];
    VerificationTest[unitaryQ[cm["Fredkin"[]]], True, TestID -> "Fredkin-unitary"],
    {1}
]

EndTestSection[]


BeginTestSection["Deutsch - half-angle convention (characterization)"]

(* QF's Deutsch[theta] doubly-controlled block is i R_x(theta) = i exp(-i theta X
   / 2), so the Toffoli is reached at theta = pi (not the literature's pi/2). A
   user transcribing theta from a reference that reaches Toffoli at pi/2 gets a
   gate off by a factor of two. These lock QF's actual convention. *)
VerificationTest[cmO["Deutsch"[Pi], {1, 2, 3}], toffM, SameTest -> matEq, TestID -> "Deutsch-Pi=Toffoli"]
VerificationTest[
    cmO["Deutsch"[t], {1, 2, 3}][[7 ;; 8, 7 ;; 8]],
    I MatrixExp[-I t / 2 PX],
    SameTest -> symEq, TestID -> "Deutsch-block=i-Rx"
]
VerificationTest[matEq[cmO["Deutsch"[Pi / 2], {1, 2, 3}], toffM], False, TestID -> "Deutsch-Pi/2-not-Toffoli"]

EndTestSection[]


BeginTestSection["Multiplexer / BlockDiagonal"]

(* Multiplexer[A, B, ...] is the direct sum (uniformly-controlled / select). *)
VerificationTest[cm["Multiplexer"[QuantumOperator["X"], QuantumOperator["Z"]]], ArrayFlatten[{{PX, 0}, {0, {{1, 0}, {0, -1}}}}], SameTest -> matEq, TestID -> "Multiplexer-XZ"]
VerificationTest[cm["BlockDiagonal"[QuantumOperator["X"], QuantumOperator["Z"]]], cm["Multiplexer"[QuantumOperator["X"], QuantumOperator["Z"]]], SameTest -> matEq, TestID -> "BlockDiagonal=Multiplexer"]
VerificationTest[cm["Multiplexer"[QuantumOperator["X"], QuantumOperator["Z"], QuantumOperator["H"]]], ArrayFlatten[{{PX, 0, 0}, {0, {{1, 0}, {0, -1}}, 0}, {0, 0, Hm}}], SameTest -> matEq, TestID -> "Multiplexer-3ops"]

EndTestSection[]


BeginTestSection["Braid - unitary and Yang-Baxter"]

(* Both braid variants are unitary solutions of the Yang-Baxter equation
   (R (x) I)(I (x) R)(R (x) I) = (I (x) R)(R (x) I)(I (x) R). *)
ybe[m_] := Module[{r12 = KroneckerProduct[m, I2], r23 = KroneckerProduct[I2, m]}, matEq[r12 . r23 . r12, r23 . r12 . r23]];
Do[
    VerificationTest[unitaryQ[cm["Braid"[v]]], True, TestID -> "Braid-unitary-" <> v];
    VerificationTest[ybe[cm["Braid"[v]]], True, TestID -> "Braid-YBE-" <> v],
    {v, {"Bell", "SwapPhase"}}
]
VerificationTest[cm["Braid"[]], cm["Braid"["Bell"]], SameTest -> matEq, TestID -> "Braid-default=Bell"]

EndTestSection[]
