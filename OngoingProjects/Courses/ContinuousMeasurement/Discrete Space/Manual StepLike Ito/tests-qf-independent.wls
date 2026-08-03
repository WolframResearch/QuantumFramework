#!/usr/bin/env wolframscript
(* QF-free test harness for manual-implementation-ito-qf-independent.wl .
   Sections:
     0  validate the independent Lindblad reference on analytic amplitude damping
     1  exact reproduction: new T vs a matrix-mode transcription of the notebook
        algorithm, on every notebook example (d=2, eta=None, V=None)
     2  positivity + trace preserved along every trajectory
     3  ensemble average -> Lindblad master equation (physics)
     4  edge/audit discriminators: d=4 two-qubit (paper), eta<1 qutrit, V!=None qutrit
     5  structure preservation at large dt *)

out[s_] := WriteString[$Output, s, "\n"];
$ProgressReporting = False;
fails = 0;
check[label_, cond_, extra_ : ""] := (
   out["  [" <> If[TrueQ[cond], "PASS", fails++; "FAIL"] <> "] " <> label <> If[extra === "", "", "  (" <> extra <> ")"]];
   TrueQ[cond]);

Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/ContinuousMeasurement/Manual StepLike Ito/manual-implementation-ito-qf-independent.wl"];

(* ---- standard operators/states, built by hand (no QF) ---- *)
X = PauliMatrix[1]; Y = PauliMatrix[2]; Z = PauliMatrix[3];
Jm = (X - I Y)/2;                       (* qubit lowering = {{0,0},{1,0}} *)
ket = <|"0" -> {{1, 0}, {0, 0}}, "+" -> {{1, 1}, {1, 1}}/2.|>;
id2 = IdentityMatrix[2];

(* ---- notebook algorithm, matrix mode (idFn selects IdentityMatrix dim) ---- *)
origR[rho_, L_, dt_, dw_, eta_] :=
  MapThread[Sqrt[#3] Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &, {L, dw, eta}];
origAlg[rho0_, Hmat_, Lmat_, etaSpec_, Vmat_, dt_, tf_, iddim_] :=
 Module[{dw, Heff, etam},
  etam = If[etaSpec === None, ConstantArray[1, Length[Lmat]], etaSpec];
  Heff = IdentityMatrix[iddim] - dt (I Hmat
       + .5 Sum[l . l + ConjugateTranspose[l] . l, {l, Lmat}]
       + If[Vmat === None, 0, .5 Sum[l . l + ConjugateTranspose[l] . l, {l, Vmat}]]);
  dw = RandomVariate[NormalDistribution[0, Sqrt[dt]], {Length[Lmat], Floor[tf/dt]}];
  FoldList[
   Module[{M, st, Ups},
     Ups = Total[Sqrt[etam] Lmat Sow@origR[#1, Lmat, dt, #2, etam]];
     M = Heff + .5 Ups . Ups + Ups;
     st = M . #1 . ConjugateTranspose[M]
        + If[Vmat === None, 0, dt Sum[l . #1 . ConjugateTranspose[l], {l, Vmat}]]
        + dt Sum[(1 - etam[[l]]) Lmat[[l]] . #1 . ConjugateTranspose[Lmat[[l]]], {l, Range[Length@Lmat]}];
     st/Tr[st]] &,
   rho0, Transpose[dw]]];

(* ---- independent Lindblad reference: rho(t) = unvec(Exp[Lsup t] vec(rho0)) ---- *)
lindbladExpm[rho0_, H_, collapse_, tf_] := Module[{d = Length[H], idm, Lsup},
  idm = IdentityMatrix[d];
  Lsup = -I (KroneckerProduct[H, idm] - KroneckerProduct[idm, Transpose[H]])
    + Total[Map[KroneckerProduct[#, Conjugate[#]]
       - (1/2) KroneckerProduct[ConjugateTranspose[#] . #, idm]
       - (1/2) KroneckerProduct[idm, Transpose[ConjugateTranspose[#] . #]] &, collapse]];
  ArrayReshape[MatrixExp[N[Lsup] tf] . Flatten[N[rho0]], {d, d}]];

tdist[a_, b_] := Max[Abs[Flatten[a - b]]];
minEig[rho_] := Min[Re@Eigenvalues[rho]];
hermErr[rho_] := Max[Abs[Flatten[rho - ConjugateTranspose[rho]]]];   (* 0 iff Hermitian *)

out["=================================================================="];
out["SECTION 0  independent Lindblad reference vs analytic amplitude damping"];
out["=================================================================="];
Module[{kappa = 0.7, tt = 1.3, rho, an00, an01},
  rho = lindbladExpm[ket["+"], 0. id2, {Sqrt[kappa] Jm}, tt];
  an00 = 0.5 Exp[-kappa tt]; an01 = 0.5 Exp[-kappa tt/2];
  check["rho00(t) = (1/2)e^{-k t}", Abs[rho[[1, 1]] - an00] < 1*^-9, "|d|=" <> ToString[Abs[rho[[1, 1]] - an00]]];
  check["rho01(t) = (1/2)e^{-k t/2}", Abs[rho[[1, 2]] - an01] < 1*^-9, "|d|=" <> ToString[Abs[rho[[1, 2]] - an01]]];
];

out[""];
out["=================================================================="];
out["SECTION 1  exact reproduction: new T vs notebook algorithm (d=2)"];
out["=================================================================="];
examples = {
  <|"name" -> "gamma>>Omega", "H" -> 1.0 X, "L" -> {Sqrt[2.0] Z}, "rho0" -> ket["0"], "dt" -> 0.005, "tf" -> 20.|>,
  <|"name" -> "gamma<<Omega", "H" -> 1.0 X, "L" -> {Sqrt[0.1] Z}, "rho0" -> ket["0"], "dt" -> 0.005, "tf" -> 20.|>,
  <|"name" -> "diffusion(J-)", "H" -> 0.5 (2.0 X + 1.0 Z), "L" -> {Sqrt[0.2] Jm}, "rho0" -> ket["+"], "dt" -> 0.01, "tf" -> 10.|>,
  <|"name" -> "WM Fig4.6", "H" -> 0.5 (3.0 Y), "L" -> {Sqrt[1.0] Jm}, "rho0" -> ket["+"], "dt" -> 0.01, "tf" -> 10.|>,
  <|"name" -> "two-point(trunc)", "H" -> 0.5 (3.5 Y), "L" -> {Sqrt[0.7] Jm}, "rho0" -> ket["0"], "dt" -> 0.01, "tf" -> 8.|>};
Do[
  Module[{ex = examples[[i]], seed = 100 + i, oT, oR, nT, nR, dTraj, dRec},
   SeedRandom[seed]; {oT, oR} = Reap[origAlg[ex["rho0"], ex["H"], ex["L"], None, None, ex["dt"], ex["tf"], 2]];
   SeedRandom[seed]; {nT, nR} = Reap[\[ScriptCapitalT][ex["rho0"], ex["H"], ex["L"], ex["dt"], ex["tf"]]];
   dTraj = Max[Abs[Flatten[oT - nT]]];
   dRec = Max[Abs[Flatten[oR[[1]] - nR[[1]]]]];
   check[ex["name"] <> ": trajectory matches notebook", dTraj < 1*^-10, "max|d|=" <> ToString[dTraj]];
   check[ex["name"] <> ": readouts (Sow) match notebook", dRec < 1*^-10, "max|d|=" <> ToString[dRec]];
   ], {i, Length[examples]}];

out[""];
out["=================================================================="];
out["SECTION 2  positivity + trace along the new trajectory"];
out["=================================================================="];
Do[
  Module[{ex = examples[[i]], traj},
   SeedRandom[200 + i]; traj = \[ScriptCapitalT][ex["rho0"], ex["H"], ex["L"], ex["dt"], ex["tf"]];
   check[ex["name"] <> ": all states PSD", Min[minEig /@ traj] > -1*^-9, "min eig=" <> ToString[Min[minEig /@ traj]]];
   check[ex["name"] <> ": trace = 1", Max[Abs[Tr[#] - 1 & /@ traj]] < 1*^-9, "max|Tr-1|=" <> ToString[Max[Abs[Tr[#] - 1 & /@ traj]]]];
   check[ex["name"] <> ": Hermitian", Max[hermErr /@ traj] < 1*^-12, "max|rho-rho^dag|=" <> ToString[Max[hermErr /@ traj]]];
   ], {i, Length[examples]}];

out[""];
out["=================================================================="];
out["SECTION 1b  the public readout \[ScriptCapitalR] reproduces the integrator's sown dy"];
out["=================================================================="];
Module[{ex = examples[[3]], seed = 400, traj, rec, dwReco, dyR},
  SeedRandom[seed]; {traj, rec} = Reap[\[ScriptCapitalT][ex["rho0"], ex["H"], ex["L"], ex["dt"], ex["tf"]]];
  (* invert dy_r = Sqrt[1] Tr[(L_r+L_r^dag) rho0] dt + dw_r for dw, then re-run \[ScriptCapitalR] *)
  dwReco = MapThread[#1 - Tr[(#2 + ConjugateTranspose[#2]) . ex["rho0"]] ex["dt"] &, {rec[[1, 1]], ex["L"]}];
  dyR = \[ScriptCapitalR][ex["rho0"], ex["L"], ex["dt"], dwReco, ConstantArray[1., Length[ex["L"]]]];
  check["\[ScriptCapitalR] round-trips the first sown readout", Max[Abs[dyR - rec[[1, 1]]]] < 1*^-12, "max|d|=" <> ToString[Max[Abs[dyR - rec[[1, 1]]]]]];
];

out[""];
out["=================================================================="];
out["SECTION 3  ensemble average -> Lindblad master equation (d=2, J-)"];
out["=================================================================="];
Module[{kappa = 0.6, dt = 0.02, tf = 6., K = 800, H, L, rho0, avg, ref, td},
  H = 0.5 X; L = {Sqrt[kappa] Jm}; rho0 = ket["+"];
  SeedRandom[321];
  avg = Mean[Table[Last@\[ScriptCapitalT][rho0, H, L, dt, tf], K]];
  ref = lindbladExpm[rho0, H, L, tf];
  td = tdist[avg, ref];
  check["ensemble mean (K=" <> ToString[K] <> ") ~ Lindblad", td < 0.03, "max|d|=" <> ToString[td] <> " (MC ~1/Sqrt[K]=" <> ToString[N[1/Sqrt[K]]] <> ")"];
];

out[""];
out["=================================================================="];
out["SECTION 4  edge / audit discriminators"];
out["=================================================================="];

out["-- (4a) d=4 two-qubit (paper's own example): dimension fix --"];
Module[{om = 1., kap = 0.01, k1 = 0.005, k2 = 0.005, H, L, eta, rho0, crashed, nT, avg, ref, td},
  H = (om/2) (KroneckerProduct[X, id2] + KroneckerProduct[id2, X]) + kap KroneckerProduct[Z, Z];
  L = {Sqrt[2 k1] KroneckerProduct[Z, id2], Sqrt[2 k2] KroneckerProduct[id2, Z]};
  eta = {0.85, 0.85}; rho0 = IdentityMatrix[4]/4.;
  (* the notebook's Heff line, verbatim, at d=4: IdentityMatrix[2] - dt(...) cannot form.
     Checked structurally (no FoldList) to avoid the nested-held-expression blow-up. *)
  crashed = ! MatrixQ[IdentityMatrix[2] - 0.02 (I H + 0.5 Total[(# . # + ConjugateTranspose[#] . #) & /@ L])];
  check["notebook code (IdentityMatrix[2]) cannot form Heff at d=4", crashed];
  SeedRandom[1]; nT = \[ScriptCapitalT][rho0, H, L, eta, 0.02, 1.];
  check["new T runs at d=4, PSD + trace", Min[minEig /@ nT] > -1*^-9 && Max[Abs[Tr[#] - 1 & /@ nT]] < 1*^-9];
  (* physics: Z(x)I and I(x)Z are involutory so eta-bug is normalized away here; average -> Lindblad *)
  SeedRandom[9]; avg = Mean[Table[Last@\[ScriptCapitalT][rho0, H, L, eta, 0.02, 3.], 300]];
  ref = lindbladExpm[rho0, H, L, 3.];
  td = tdist[avg, ref];
  check["d=4 ensemble mean ~ Lindblad", td < 0.03, "max|d|=" <> ToString[td]];
];

out["-- (4b) eta<1 qutrit (spin-1 J-, L^2 != 0): efficiency-factor fix --"];
Module[{Jz, Jm3, gam = 0.6, eta = 0.5, H, L, rho0, dt = 0.02, tf = 5., K = 400, refL, aNew, aOrig, dNew, dOrig},
  Jz = DiagonalMatrix[{1., 0., -1.}];
  Jm3 = Sqrt[2.] {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};            (* spin-1 lowering, (Jm3)^2 != 0 *)
  H = 1.0 Jz; L = {Sqrt[gam] Jm3}; rho0 = DiagonalMatrix[{1., 0., 0.}];
  refL = lindbladExpm[rho0, H, L, tf];                          (* eta-independent unconditional ME *)
  SeedRandom[55]; aNew = Mean[Table[Last@\[ScriptCapitalT][rho0, H, L, {eta}, dt, tf], K]];
  SeedRandom[55]; aOrig = Mean[Table[Last@origAlg[rho0, H, L, {eta}, None, dt, tf, 3], K]];
  dNew = tdist[aNew, refL]; dOrig = tdist[aOrig, refL];
  out["    new-mean vs Lindblad = " <> ToString[dNew] <> " ;  notebook(eta-unweighted)-mean vs Lindblad = " <> ToString[dOrig]];
  check["new T matches eta-independent Lindblad", dNew < 0.04, "max|d|=" <> ToString[dNew]];
  check["notebook eta-handling farther from Lindblad (decisive magnitude in tests_disc.wls)", dOrig > dNew, "dOrig=" <> ToString[dOrig] <> " dNew=" <> ToString[dNew]];
];

out["-- (4c) V != None qutrit (spin-1 J- as unmonitored channel): spurious-V^2 fix --"];
Module[{Jz, Jm3, gl = 0.3, gv = 0.4, H, L, V, rho0, dt = 0.02, tf = 5., K = 400, refL, aNew, aOrig, dNew, dOrig},
  Jz = DiagonalMatrix[{1., 0., -1.}];
  Jm3 = Sqrt[2.] {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  H = 1.0 Jz; L = {Sqrt[gl] Jm3}; V = {Sqrt[gv] Jm3}; rho0 = DiagonalMatrix[{1., 0., 0.}];
  refL = lindbladExpm[rho0, H, Join[L, V], tf];                 (* both channels dissipate *)
  SeedRandom[77]; aNew = Mean[Table[Last@\[ScriptCapitalT][rho0, H, L, None, V, dt, tf], K]];
  SeedRandom[77]; aOrig = Mean[Table[Last@origAlg[rho0, H, L, None, V, dt, tf, 3], K]];
  dNew = tdist[aNew, refL]; dOrig = tdist[aOrig, refL];
  out["    new-mean vs Lindblad = " <> ToString[dNew] <> " ;  notebook(spurious V^2)-mean vs Lindblad = " <> ToString[dOrig]];
  check["new T matches Lindblad with V channel", dNew < 0.04, "max|d|=" <> ToString[dNew]];
  check["notebook V-handling farther from Lindblad (decisive magnitude in tests_disc.wls)", dOrig > dNew, "dOrig=" <> ToString[dOrig] <> " dNew=" <> ToString[dNew]];
];

out[""];
out["=================================================================="];
out["SECTION 5  structure preservation at large dt (positivity for any step)"];
out["=================================================================="];
Module[{H, L, rho0, traj},
  H = 0.5 (3.0 Y); L = {Sqrt[1.0] Jm}; rho0 = ket["+"];
  SeedRandom[5]; traj = \[ScriptCapitalT][rho0, H, L, 0.5, 20.];  (* huge dt=0.5 *)
  check["dt=0.5: all states PSD", Min[minEig /@ traj] > -1*^-9, "min eig=" <> ToString[Min[minEig /@ traj]]];
  check["dt=0.5: trace = 1", Max[Abs[Tr[#] - 1 & /@ traj]] < 1*^-9];
];

out[""];
out["=================================================================="];
out["SECTION 6  two NON-COMMUTING monitored operators (d=2): L1 L2 != L2 L1"];
out["=================================================================="];
Module[{H = 0.5 Z, L = {Sqrt[0.5] Jm, Sqrt[0.3] X}, r0 = ket["+"], dt = 0.01, tf = 5., seed = 909, oT, nT, dTraj, refL, avg, td},
  check["operators do not commute", (L[[1]] . L[[2]]) =!= (L[[2]] . L[[1]])];
  SeedRandom[seed]; oT = origAlg[r0, H, L, None, None, dt, tf, 2];
  SeedRandom[seed]; nT = \[ScriptCapitalT][r0, H, L, dt, tf];
  dTraj = Max[Abs[Flatten[oT - nT]]];
  check["non-commuting: reproduces notebook algorithm", dTraj < 1*^-9, "max|d|=" <> ToString[dTraj]];
  refL = lindbladExpm[r0, H, L, tf];
  SeedRandom[7]; avg = Mean[Table[Last@\[ScriptCapitalT][r0, H, L, dt, tf], 800]];
  td = tdist[avg, refL];
  check["non-commuting: ensemble mean ~ Lindblad", td < 0.03, "max|d|=" <> ToString[td]];
];

out[""];
out["=================================================================="];
out["SECTION 7  convergence: ensemble mean tracks Lindblad under dt refinement"];
out["=================================================================="];
Module[{H = 0.5 X, L = {Sqrt[0.6] Jm}, r0 = ket["+"], tf = 4., K = 800, ref, d1, d2},
  ref = lindbladExpm[r0, H, L, tf];
  SeedRandom[11]; d1 = tdist[Mean[Table[Last@\[ScriptCapitalT][r0, H, L, 0.04, tf], K]], ref];
  SeedRandom[11]; d2 = tdist[Mean[Table[Last@\[ScriptCapitalT][r0, H, L, 0.01, tf], K]], ref];
  out["    mean vs Lindblad at dt=0.04: " <> ToString[d1] <> " ; at dt=0.01: " <> ToString[d2]];
  check["ensemble mean ~ Lindblad at coarse dt=0.04", d1 < 0.05, "d1=" <> ToString[d1]];
  check["ensemble mean ~ Lindblad at fine dt=0.01 (refinement does not degrade)", d2 < 0.04, "d2=" <> ToString[d2]];
];

out[""];
out["=================================================================="];
out["SECTION 8  accuracy breakdown at large dt: state stays physical but the mean is biased"];
out["=================================================================="];
Module[{H = 0.5 X, L = {Sqrt[0.6] Jm}, r0 = ket["+"], tf = 4., K = 1200, ref, biases, bigTraj},
  ref = lindbladExpm[r0, H, L, tf];
  biases = Table[SeedRandom[13]; tdist[Mean[Table[Last@\[ScriptCapitalT][r0, H, L, dtv, tf], K]], ref], {dtv, {0.02, 0.2, 0.5, 1.0}}];
  out["    |mean - Lindblad| at dt={0.02,0.2,0.5,1.0} = " <> ToString[biases]];
  bigTraj = (SeedRandom[1]; \[ScriptCapitalT][r0, H, L, 1.0, tf]);
  check["at dt=1.0 the state is still physical (structure preserved)", Min[minEig /@ bigTraj] > -1*^-9 && Max[hermErr /@ bigTraj] < 1*^-12];
  check["mean is accurate at small dt=0.02", biases[[1]] < 0.02, "bias=" <> ToString[biases[[1]]]];
  check["first-order bias grows with dt and breaks accuracy at dt=1.0", biases[[-1]] > 3 biases[[1]] && biases[[-1]] > 0.04, "bias(dt=1)=" <> ToString[biases[[-1]]]];
];

out[""];
out["=================================================================="];
out["SECTION 9  guarded signature: 5/6/7-arg bind, malformed call fails fast"];
out["=================================================================="];
Module[{H = 0.5 X, L = {Sqrt[0.6] Jm}, r0 = ket["+"], t5, t6, t7, badEta, badVs, badDt},
  t5 = \[ScriptCapitalT][r0, H, L, 0.1, 1.];
  t6 = \[ScriptCapitalT][r0, H, L, {0.7}, 0.1, 1.];
  t7 = \[ScriptCapitalT][r0, H, L, {0.7}, {Sqrt[0.2] Jm}, 0.1, 1.];
  badEta = \[ScriptCapitalT][r0, H, L, 0.7, 0.1, 1.];             (* scalar in eta slot *)
  badVs = \[ScriptCapitalT][r0, H, L, {0.7}, 0.5, 0.05, 0.2];     (* scalar in Vs slot *)
  badDt = \[ScriptCapitalT][r0, H, L, "notNumeric", 1.];         (* non-numeric dt *)
  check["5-arg returns Floor[tf/dt]+1 states", MatrixQ[Last[t5]] && Length[t5] === Floor[1./0.1] + 1];
  check["6-arg (eta) returns a trajectory", MatrixQ[Last[t6]]];
  check["7-arg (eta, V) returns a trajectory", MatrixQ[Last[t7]]];
  check["scalar in eta slot fails fast (unevaluated)", Head[badEta] === \[ScriptCapitalT] && ! MatrixQ[badEta]];
  check["scalar in Vs slot fails fast (unevaluated)", Head[badVs] === \[ScriptCapitalT] && ! MatrixQ[badVs]];
  check["non-numeric dt fails fast (unevaluated)", Head[badDt] === \[ScriptCapitalT] && ! MatrixQ[badDt]];
];

out[""];
out["=================================================================="];
out["TOTAL FAILURES: " <> ToString[fails]];
out["=================================================================="];
If[fails > 0, Exit[1]];
