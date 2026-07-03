(* Physics/math tests for the open-system, random, and ladder named operators:
   the oscillator ladder ("+", "-", Up, Down, I+, I-), RandomUnitary,
   RandomHermitian, Switch, Channel, Measurement, Liouvillian, Hamiltonian.

   Each is pinned by the property that defines it: ladder operators build the
   number operator, random operators keep their ensemble invariant, the quantum
   switch is the two-order superposition, channels act as the stated CPTP map,
   and the Liouvillian is trace-preserving with Hamiltonian = i * Liouvillian. *)

cm[x_] := Normal[QuantumOperator[x]["Computational"]["Matrix"]]
matEq[a_, b_] := Dimensions[a] === Dimensions[b] && Chop[N[a - b], 10^-10] === ConstantArray[0, Dimensions[a]]
unitaryQ[m_] := matEq[m . ConjugateTranspose[m], IdentityMatrix[Length[m]]]
PX = {{0, 1}, {1, 0}};
PZ = {{1, 0}, {0, -1}};
Sm = {{1, 0}, {0, I}};
I2 = IdentityMatrix[2];

BeginTestSection["Oscillator ladder (+, -, Up, Down)"]
(* "+"[d] is the truncated annihilation operator a (superdiagonal sqrt(n));
   "-"[d] is the creation a^dagger. a^dagger a is the number operator. *)
VerificationTest[cm["+"[2]], {{0, 1}, {0, 0}}, SameTest -> matEq, TestID -> "plus-2=sigma+"]
VerificationTest[cm["+"[3]], {{0, 1, 0}, {0, 0, Sqrt[2]}, {0, 0, 0}}, SameTest -> matEq, TestID -> "plus-3=annihilation"]
VerificationTest[cm["-"[3]], ConjugateTranspose[{{0, 1, 0}, {0, 0, Sqrt[2]}, {0, 0, 0}}], SameTest -> matEq, TestID -> "minus-3=creation"]
VerificationTest[ConjugateTranspose[cm["+"[3]]], cm["-"[3]], SameTest -> matEq, TestID -> "minus=dagger-plus"]
VerificationTest[cm["-"[3]] . cm["+"[3]], DiagonalMatrix[{0, 1, 2}], SameTest -> matEq, TestID -> "adag.a=number"]
VerificationTest[cm["+"[3]] . cm["-"[3]], DiagonalMatrix[{1, 2, 0}], SameTest -> matEq, TestID -> "a.adag-truncated"]
VerificationTest[cm["Up"[3]], cm["+"[3]], SameTest -> matEq, TestID -> "Up=plus"]
VerificationTest[cm["Down"[3]], cm["-"[3]], SameTest -> matEq, TestID -> "Down=minus"]
VerificationTest[cm["I+"[3]], cm["+"[3]], SameTest -> matEq, TestID -> "Iplus=plus"]
VerificationTest[cm["I-"[3]], cm["-"[3]], SameTest -> matEq, TestID -> "Iminus=minus"]
EndTestSection[]

BeginTestSection["RandomUnitary / RandomHermitian - ensemble invariants"]
SeedRandom[1234];
VerificationTest[unitaryQ[cm["RandomUnitary"]], True, TestID -> "RandomUnitary-d2"]
VerificationTest[unitaryQ[Normal[QuantumOperator["RandomUnitary"[3]]["Computational"]["Matrix"]]], True, TestID -> "RandomUnitary-d3"]
VerificationTest[unitaryQ[Normal[QuantumOperator["RandomUnitary"[{2, 2}]]["Computational"]["Matrix"]]], True, TestID -> "RandomUnitary-2x2"]
VerificationTest[HermitianMatrixQ[cm["RandomHermitian"]], True, TestID -> "RandomHermitian-d2"]
VerificationTest[HermitianMatrixQ[Normal[QuantumOperator["RandomHermitian"[QuantumBasis[3]]]["Computational"]["Matrix"]]], True, TestID -> "RandomHermitian-d3"]
EndTestSection[]

BeginTestSection["Switch - quantum switch"]
(* Switch[A, B] applies A B in the |0> control branch and B A in the |1> branch
   (control on qubit 1), i.e. block-diag(A B, B A). *)
VerificationTest[
    Normal[QuantumOperator["Switch"[QuantumOperator["X"], QuantumOperator["S"[]]], {1, 2}]["Computational"]["Matrix"]],
    ArrayFlatten[{{PX . Sm, 0}, {0, Sm . PX}}], SameTest -> matEq, TestID -> "Switch=blockdiag(AB,BA)"
]
EndTestSection[]

BeginTestSection["Channel - bit-flip CPTP action and delegation"]
(* Channel[BitFlip[p]] : rho -> (1-p) rho + p X rho X, a doubled-space
   superoperator. It delegates to QuantumChannel[...]["DiscardExtraQudits"]. *)
mbf = cm["Channel"["BitFlip"[3/10]]];
VerificationTest[ArrayReshape[mbf . Flatten[{{1, 0}, {0, 0}}], {2, 2}], {{7/10, 0}, {0, 3/10}}, SameTest -> matEq, TestID -> "BitFlip-on-0"]
VerificationTest[ArrayReshape[mbf . Flatten[{{0, 0}, {0, 1}}], {2, 2}], {{3/10, 0}, {0, 7/10}}, SameTest -> matEq, TestID -> "BitFlip-on-1"]
VerificationTest[mbf, Normal[QuantumChannel["BitFlip"[3/10]]["DiscardExtraQudits"]["Computational"]["Matrix"]], SameTest -> matEq, TestID -> "Channel-delegation"]
EndTestSection[]

BeginTestSection["Measurement - delegation"]
(* Measurement[args] delegates to QuantumMeasurementOperator[args] with the extra
   pointer qudits discarded. *)
VerificationTest[cm["Measurement"["Z"]], Normal[QuantumMeasurementOperator["Z"]["DiscardExtraQudits"]["Computational"]["Matrix"]], SameTest -> matEq, TestID -> "Measurement-delegation"]
EndTestSection[]

BeginTestSection["Liouvillian / Hamiltonian - Lindblad generator"]
(* Liouvillian[H, {L}, {g}] is the vectorized Lindblad generator; it is
   trace-preserving (vec(I)^T L = 0). With no dissipation it is the commutator
   generator -i(H (x) I - I (x) H^T), and Hamiltonian[...] = i Liouvillian[...]. *)
liou = Normal[QuantumOperator["Liouvillian"[QuantumOperator[PZ], {QuantumOperator["X"]}, {1/2}]]["Computational"]["Matrix"]];
VerificationTest[{1, 0, 0, 1} . liou, {0, 0, 0, 0}, SameTest -> matEq, TestID -> "Liouvillian-trace-preserving"]
liouH = Normal[QuantumOperator["Liouvillian"[QuantumOperator[PZ]]]["Computational"]["Matrix"]];
VerificationTest[liouH, -I (KroneckerProduct[PZ, I2] - KroneckerProduct[I2, Transpose[PZ]]), SameTest -> matEq, TestID -> "Liouvillian-Hamiltonian-only"]
VerificationTest[Normal[QuantumOperator["Hamiltonian"[QuantumOperator[PZ]]]["Computational"]["Matrix"]], I liouH, SameTest -> matEq, TestID -> "Hamiltonian=i-Liouvillian"]
(* The rule form [H, {L} -> {g}] agrees with the positional form [H, {L}, {g}]. *)
VerificationTest[
    Normal[QuantumOperator["Liouvillian"[QuantumOperator[PZ], {QuantumOperator["X"]} -> {1/2}]]["Computational"]["Matrix"]],
    liou, SameTest -> matEq, TestID -> "Liouvillian-rule=positional"
]
EndTestSection[]
