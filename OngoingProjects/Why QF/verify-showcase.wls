(* Consolidated, boolean-asserted re-verification of QF-Showcase.md against the
   live working-tree kernel (QuantumFramework 2.0.2, NOT the installed paclet).
   Each check either compares QF to an INDEPENDENT oracle (textbook matrices, a
   frame-correct trace, the parity operator, a closed form) or reproduces a
   specific showcase cell verbatim (a QF-internal identity the document itself
   asserts), decided by a boolean; no eyeballing.  Failure modes named in each key. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];

(* Regimes covered (each showcase cell re-verified at its own fixed dimension):
     general symbolic  : Rabi P(t), Lindblad rho(t), QAOA landscape, purity loss, Werner(p), ZFS(d,e)
     exactly solvable  : <Jx>=1, textbook spin-1 H_truth oracle, CHSH (2+Sqrt2)/4, [X,Y]=2iZ
     asymptotic/limit  : ZFS e->0 degenerate {0,2d,2d}; S(t)->0 product points (Limit); large-detuning Rabi -> 0
     numeric reference : Wigner/Husimi NIntegrate; S(pi/8)=1; Grover integer arithmetic
     failure/edge      : |+_x> mislabel gotcha; removable singularity of S via Limit; Werner boundary p=1/2 *)

results = <||>;
check[name_, got_, exp_] := (results[name] = {TrueQ[got === exp], got, exp});
checkTol[name_, got_, exp_, tol_ : 10^-2] := (results[name] = {NumericQ[got] && Abs[got - exp] < tol, got, exp});

(* ---------- Claim 2: version + faithful spin-1 angular momentum ---------- *)
check["version 2.0.2", PacletObject["Wolfram/QuantumFramework"]["Version"], "2.0.2"];

mx  = QuantumState[{0, 0, 1}, "JX"[1]];
mxZ = QuantumState[mx, "Computational"[3]];
jxOp = QuantumOperator["JX"[1]];
check["<Jx>=1 in JX basis (L417)", mx["Dagger"][jxOp[mx]]["Scalar"], 1];
check["<Jx>=1 in computational (L417)", mxZ["Dagger"][jxOp[mxZ]]["Scalar"], 1];
(* JX basis is real: the imaginary part of the whole basis matrix vanishes symbolically *)
check["JX basis real (symbolic)",
   FullSimplify[Im[Normal[QuditBasis["JX"[1]]["Matrix"]]]], ConstantArray[0, {3, 3}]];

(* ZFS spectrum vs an INDEPENDENT textbook spin-1 build, not the doc literal *)
{jx, jy, jz} = QuantumOperator /@ {"JX"[1], "JY"[1], "JZ"[1]};
hzfs = 2 \[FormalD] jz^2 + \[FormalE] (jx^2 - jy^2);
zfsSpec = QuantumOperator[hzfs, QuantumBasis[hzfs["Eigenvectors"]]]["Matrix"] // Normal // FullSimplify;
Jzt = DiagonalMatrix[{1, 0, -1}];
Jxt = (1/Sqrt[2]) {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}};
Jyt = (1/Sqrt[2]) {{0, -I, 0}, {I, 0, -I}, {0, I, 0}};
Htruth = 2 \[FormalD] Jzt . Jzt + \[FormalE] (Jxt . Jxt - Jyt . Jyt);
check["ZFS spectrum == textbook oracle (L438)",
   Sort[Eigenvalues[zfsSpec]] // FullSimplify, Sort[Eigenvalues[Htruth]] // FullSimplify];
(* limiting case: axial e -> 0 gives the degenerate doublet {0, 2d, 2d} *)
check["ZFS e->0 degenerate {0,2d,2d}",
   Sort[Eigenvalues[zfsSpec /. \[FormalE] -> 0]] // FullSimplify, Sort[{0, 2 \[FormalD], 2 \[FormalD]}]];
(* invariant: hzfs is Hermitian in the computational frame, symbolic over real d,e *)
hzfsComp = Normal[QuantumOperator[hzfs, "Computational"[3]]["Matrix"]];
check["ZFS Hermitian in comp frame (symbolic real d,e)",
   FullSimplify[ComplexExpand[hzfsComp - ConjugateTranspose[hzfsComp]], {\[FormalD] \[Element] Reals, \[FormalE] \[Element] Reals}], ConstantArray[0, {3, 3}]];

(* ---------- Two formerly-open kernel defects are FIXED at HEAD ---------- *)
(* cross-basis composition *)
xB = QuantumOperator["X", "PauliX"]; zB = QuantumOperator["Z", "PauliZ"];
check["cross-basis X@Z == sigma_x.sigma_z",
   Normal[QuantumOperator[xB[zB], "Computational"[2]]["Matrix"]] // FullSimplify,
   PauliMatrix[1] . PauliMatrix[3]];
(* complex-basis bra route: FRAME-CORRECT oracle (computational DM) vs bra-op-ket *)
Ycomp = Normal[QuantumOperator["Y"]["Matrix"]];
oracle[st_] := Tr[Normal[st["Computational"]["DensityMatrix"]] . Ycomp] // FullSimplify;
brk[st_] := st["Dagger"][QuantumOperator["Y"][st]]["Scalar"];
yInY = QuantumState[QuantumState[{1, I}/Sqrt[2]], "Y"];   (* true |y+> rebased into complex Y basis *)
check["complex-bra route == trace oracle (|y+>)", brk[yInY], oracle[yInY]];
check["complex-bra route == +1 (|y+>)", brk[yInY], 1];
(* the construction gotcha, correctly labelled: {1,I}/Sqrt2 tagged Y is |+_x>, not |y+> *)
plusX = QuantumState[{1, I}/Sqrt[2], "Y"];
check["{1,I}/Sqrt2 in Y is |+_x> (comp vector)",
   Normal[plusX["Computational"]["StateVector"]] // FullSimplify, Exp[I Pi/4] {1, 1}/Sqrt[2] // FullSimplify];
check["|+_x> has <Y>=0 (brk==oracle)", brk[plusX], oracle[plusX]];
check["|+_x> has <Y>=0", brk[plusX], 0];

(* ---------- Flagship symbolic reproductions (Claims 1,2,3,4,5,7) ---------- *)
ClearAll[h];
h[t_] := \[Omega]0/2 QuantumOperator["Z"] + \[CapitalOmega]/2 (Cos[\[Omega] t] QuantumOperator["X"] + Sin[\[Omega] t] QuantumOperator["Y"]);
urot = QuantumEvolve[h[t], None, t];
Prabi = Assuming[{\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, t \[Element] Reals},
   FullSimplify[ComplexExpand[Abs[Normal[urot["Matrix"]][[2, 1]]]^2]]];
docRabi = \[CapitalOmega]^2/(\[CapitalOmega]^2 + (\[Omega] - \[Omega]0)^2) Sin[t/2 Sqrt[\[CapitalOmega]^2 + (\[Omega] - \[Omega]0)^2]]^2;
check["Rabi P(t) closed form (L83)",
   FullSimplify[Prabi - docRabi, {\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, t \[Element] Reals}], 0];
(* invariant: the propagator is unitary (symbolic, on the physical domain) *)
um = Normal[urot["Matrix"]];
check["Rabi propagator U(t) unitary (L60)",
   FullSimplify[ComplexExpand[ConjugateTranspose[um] . um], {\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, t \[Element] Reals}], IdentityMatrix[2]];
(* asymptotic limit: the transition-probability envelope collapses as the detuning runs away *)
check["Rabi peak suppressed at large detuning (L83)",
   Limit[\[CapitalOmega]^2/(\[CapitalOmega]^2 + (\[Omega] - \[Omega]0)^2), \[Omega] -> \[Infinity]], 0];

Block[{$Assumptions = \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0 && \[CapitalOmega]2 \[Element] Reals && t >= 0},
  \[Rho]t = QuantumEvolve[\[CapitalOmega]2/2 QuantumOperator["Z"],
     {Sqrt[\[Gamma]1] QuantumOperator[{{0, 1}, {0, 0}}], Sqrt[\[Gamma]\[Phi]] QuantumOperator["Z"]},
     QuantumState[Array[\[Rho], {2, 2}]], t];
  rtMat = FullSimplify /@ Normal[\[Rho]t["Computational"]["DensityMatrix"]];
  docRt = {{\[Rho][1, 1] + (1 - E^(-\[Gamma]1 t)) \[Rho][2, 2], \[Rho][1, 2] E^(-((\[Gamma]1/2 + 2 \[Gamma]\[Phi] + I \[CapitalOmega]2) t))},
     {\[Rho][2, 1] E^(-((\[Gamma]1/2 + 2 \[Gamma]\[Phi] - I \[CapitalOmega]2) t)), \[Rho][2, 2] E^(-\[Gamma]1 t)}};
  check["Lindblad rho(t) matrix (L272)", FullSimplify[rtMat - docRt], {{0, 0}, {0, 0}}];
  (* invariant: the Lindblad map is trace-preserving *)
  check["Lindblad map trace-preserving (L272)", FullSimplify[Tr[rtMat] - (\[Rho][1, 1] + \[Rho][2, 2])], 0];
  (* invariant: it is Hermiticity-preserving on a Hermitian input rho *)
  rtH = rtMat /. {\[Rho][1, 1] -> r1, \[Rho][2, 2] -> r2, \[Rho][1, 2] -> ra + I rb, \[Rho][2, 1] -> ra - I rb};
  check["Lindblad map Hermiticity-preserving (L272)",
     FullSimplify[rtH - ConjugateTranspose[rtH], {r1 \[Element] Reals, r2 \[Element] Reals, ra \[Element] Reals, rb \[Element] Reals}], {{0, 0}, {0, 0}}];
];

(* Heisenberg S(t): J=1 concrete for a clean curve; Limit at removable singularities, Bell point numeric *)
hheis = QuantumOperator["XX"] + QuantumOperator["YY"] + QuantumOperator["ZZ"];
sMag = QuantityMagnitude[QuantumPartialTrace[QuantumEvolve[hheis, QuantumState["01"], t], {2}]["VonNeumannEntropy"]];
check["S(t)->0 at product point t=0 (Limit)", Simplify[Limit[sMag, t -> 0]], 0];
check["S(t)->0 at product point t=pi/4 (Limit)", Simplify[Limit[sMag, t -> Pi/4]], 0];
check["S(pi/8)=1 at Bell point (exact)", Simplify[sMag /. t -> Pi/8], 1];

graph = CompleteGraph[4]; nq = VertexCount[graph];
hcost = QuantumOperator[Total[With[{id = QuantumOperator[StringJoin[ConstantArray["I", nq]]]/2},
     (id - QuantumOperator["ZZ"/2, #] &) /@ (List @@@ EdgeList[graph])]], "Label" -> "Hamiltonian"];
oneLayer = QuantumCircuitOperator[{
     QuantumCircuitOperator[{StringJoin[ConstantArray["0", nq]], "H" -> Range[nq]}],
     QuantumCircuitOperator[("R"[\[Gamma], "ZZ"] -> # &) /@ (List @@@ EdgeList[graph]), "Parameters" -> {\[Gamma]}],
     QuantumCircuitOperator[{"RX"[\[Theta]] -> VertexList[graph]}, "Parameters" -> {\[Theta]}]}, "Parameters" -> {\[Theta], \[Gamma]}];
landscape = FullSimplify[(oneLayer /* hcost /* oneLayer["Dagger"])[]["Scalar"], Element[{\[Gamma], \[Theta]}, Reals]];
docLandscape = 3/8 (7 + Cos[2 \[Theta]] + 2 Cos[4 \[Gamma]] Sin[\[Theta]]^2 - 8 Cos[\[Gamma]]^2 Sin[\[Gamma]] Sin[2 \[Theta]]);
check["QAOA landscape closed form (L183)", FullSimplify[landscape - docLandscape], 0];

psiB = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}];
(* invariant: the constructed Bloch state is normalized <psi|psi> = 1 *)
check["Bloch state normalized <psi|psi>=1", FullSimplify[psiB["Norm"], {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}], 1];
ploss = FullSimplify[psiB["Purity"] - QuantumChannel["AmplitudeDamping"[\[Gamma]]][psiB]["Purity"],
   {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals, 0 <= \[Gamma] <= 1}];
check["amplitude-damping purity loss (L246)", FullSimplify[ploss - 2 \[Gamma] (1 - \[Gamma]) Sin[\[Theta]/2]^4], 0];

sicProbs = FullSimplify[ComplexExpand[Tr[psiB["DensityMatrix"] . Normal[#]] & /@
      QuantumMeasurementOperator["TetrahedronSICPOVM"]["POVMElements"]]];
check["SIC probabilities sum to 1 (L294)", FullSimplify[Total[sicProbs]], 1];
(* the DEFINING SIC property: rank-one, equiangular Gram Tr[Ei.Ej] = (delta_ij + 1/3)/4 *)
sicE = Normal /@ QuantumMeasurementOperator["TetrahedronSICPOVM"]["POVMElements"];
check["SIC equiangular Gram Tr[Ei.Ej] (diag 1/4, off 1/12)",
   FullSimplify[Table[Tr[sicE[[i]] . sicE[[j]]], {i, 4}, {j, 4}]],
   (ConstantArray[1, {4, 4}] + 2 IdentityMatrix[4])/12];

aOp = AnnihilationOperator[];
uosc = QuantumEvolve[\[Omega] (aOp["Dagger"] @ aOp + 1/2), None, t];
(* QF's state == is phase-blind (it verifies the label alpha -> alpha e^{-i w t} but NOT
   the vacuum phase e^{-i w t/2}); test the phase directly on the amplitude list. *)
cohLHS = uosc @ CoherentState[];
cohRHS = E^(-I \[Omega] t/2) CoherentState[][\[FormalAlpha] E^(-I \[Omega] t)];
check["coherent-state evolution, phase-resolved amplitudes (L342)",
   FullSimplify[ComplexExpand[cohLHS["AmplitudesList"] - cohRHS["AmplitudesList"]],
     {\[FormalAlpha] \[Element] Reals, \[Omega] \[Element] Reals, t \[Element] Reals}], ConstantArray[0, 16]];

cluster3 = QuantumCircuitOperator[{"H" -> Range[3], "CZ", "CZ" -> {2, 3}}];
check["stabilizer === statevector (L748)",
   QuantumState[cluster3[Method -> "Stabilizer"]]["StateVector"] === cluster3[]["StateVector"], True];
stab5 = QuantumCircuitOperator[Join[Thread["H" -> Range[5]],
     {"CZ" -> {1, 2}, "CZ" -> {2, 3}, "CZ" -> {3, 4}, "CZ" -> {4, 5}}]][Method -> "Stabilizer"];
check["cluster expectations {XZIII,IZXZI,XXIII} (L704)",
   {stab5["Expectation", "XZIII"], stab5["Expectation", "IZXZI"], stab5["Expectation", "XXIII"]}, {1, 1, 0}];

check["[X,Y]=2iZ vs Pauli oracle (L794)",
   Normal[Commutator[QuantumOperator["X"], QuantumOperator["Y"]]["Matrix"]],
   PauliMatrix[1] . PauliMatrix[2] - PauliMatrix[2] . PauliMatrix[1]];

psiRY = QuantumOperator["RY"[\[Theta]]][QuantumState["0"]];
(* reproduces showcase L831 verbatim: the manual inner product form, not the bra-op-ket route *)
energy = Simplify[Conjugate[psiRY["StateVector"]] . ((QuantumOperator["X"] + QuantumOperator["Z"])[psiRY])["StateVector"], \[Theta] \[Element] Reals];
check["Maximize<X+Z>=Sqrt2 (L835)",
   First[Assuming[\[Theta] \[Element] Reals, Maximize[{energy, 0 <= \[Theta] <= 2 Pi}, \[Theta]] // FullSimplify]], Sqrt[2]];

wernerReduce = Reduce[Simplify[QuantumEntanglementMonotone[QuantumState["Werner"[p, 2]], {1}, "Realignment"], 0 < p < 1] > 0 && 0 < p < 1, p];
check["Werner threshold 0<p<1/2 (L845)", Resolve[ForAll[p, Equivalent[wernerReduce, 0 < p < 1/2]], Reals], True];

chsh = QuantumCircuitOperator["OPENQASM 3.0;
qubit[4] q;
bit[4] c;
reset q[0];
reset q[3];
U(Pi/2, 0, Pi) q[3];
U(Pi/2, 0, Pi) q[0];
swap q[0] q[3];
ctrl(0) @ negctrl(1) @ U(Pi/2, Pi, Pi) q[0] q[3];
ctrl(1) @ negctrl(0) @ U(Pi/2, 0, 0) q[0] q[3];
reset q[1];
U(Pi/2, 0, Pi) q[1];
reset q[2];
U(Pi/2, 0, Pi) q[2];
U(0, 0, 0) q[0];
ctrl(1) @ negctrl(0) @ U(Pi/2, 0, Pi) q[1] q[0];
U(Pi/4, 0, 0) q[3];
ctrl(0) @ negctrl(1) @ U(Pi/2, 0, Pi) q[2] q[3];
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];"];
chshWin = Probability[BitAnd[x, y] == BitXor[a, b], {a, x, y, b} \[Distributed] chsh[]["MultivariateDistribution"]] // FullSimplify;
check["CHSH Tsirelson (2+Sqrt2)/4 (L954)", FullSimplify[chshWin - (2 + Sqrt[2])/4], 0];

(* ---------- Doc-fix computations (Reduce forms, Wigner/Husimi measure, Grover count) ---------- *)
check["Reduce equality -> gphi==0 (L281 fix)",
   Reduce[1/(\[FormalG]1/2 + 2 \[FormalG]p) == 2/\[FormalG]1, \[FormalG]p], \[FormalG]p == 0 && \[FormalG]1 != 0];
check["ForAll-both equality is False (L281)",
   Reduce[ForAll[{\[FormalG]1, \[FormalG]p}, \[FormalG]1 > 0 && \[FormalG]p >= 0, 1/(\[FormalG]1/2 + 2 \[FormalG]p) == 2/\[FormalG]1]], False];
check["ForAll-both <= is True (L279)",
   Reduce[ForAll[{\[FormalG]1, \[FormalG]p}, \[FormalG]1 > 0 && \[FormalG]p >= 0, 1/(\[FormalG]1/2 + 2 \[FormalG]p) <= 2/\[FormalG]1]], True];

wig = WignerRepresentation[CatState[][2, 0], {-6, 6}, {-5, 5}, "GridSize" -> 160];
hus = HusimiQRepresentation[CatState[][2, 0], {-6, 6}, {-5, 5}, "GridSize" -> 160];
(* PrecisionGoal -> 2 matches the integration to this row's 0.02 tolerance, so the
   oscillatory cat-fringe integrand does not trip NIntegrate::slwcon chasing machine precision *)
checkTol["Wigner integrates to 1 over dx dp (L387)", NIntegrate[wig[x, p], {x, -6, 6}, {p, -5, 5}, PrecisionGoal -> 2], 1, 0.02];
checkTol["Husimi integrates to 2 over dx dp (L387 fix)", NIntegrate[hus[x, p], {x, -6, 6}, {p, -5, 5}], 2, 0.05];
(* W(0,0)=1/pi is EXACT via the parity eigenvalue, not "approaching": <Parity>=+1 for the
   even cat, alpha-independent.  Assert <parity> = +1 for symbolic real alpha. *)
parity = Normal[PhaseShiftOperator[Pi]["Matrix"]];  (* e^{i pi a\[Dagger]a}, tracks $FockSize *)
check["even-cat parity == +1 (symbolic alpha, W(0,0)=1/pi exact)",
   FullSimplify[Tr[Normal[CatState[][\[FormalAlpha], 0]["DensityMatrix"]] . parity], \[FormalAlpha] \[Element] Reals], 1];

check["Grover Floor[(Pi/4)Sqrt256]==12 (L622 fix)", Floor[(Pi/4) Sqrt[256]], 12];
check["Grover exact-angle round == 12 (L638/L640 code)", Round[Pi/(4 ArcSin[Sqrt[1/2^8]]) - 1/2], 12];
check["Grover naive Round == 13 (the trap, not used)", Round[(Pi/4) Sqrt[256]], 13];

(* ---------- tally ---------- *)
fails = Select[results, ! First[#] &];
Print["=============================================================="];
KeyValueMap[
  Print["  [", If[First[#2], "PASS", "FAIL"], "]  ", StringPadRight[#1, 46],
        If[First[#2], "", "  got=" <> ToString[#2[[2]]] <> "  exp=" <> ToString[#2[[3]]]]] &, results];
Print["=============================================================="];
Print[ToString[Count[results, {True, __}]], " / ", ToString[Length[results]], " PASS"];
Print[If[Length[fails] == 0, "ALL PASS", "FAILURES: " <> ToString[Keys[fails]]]];
