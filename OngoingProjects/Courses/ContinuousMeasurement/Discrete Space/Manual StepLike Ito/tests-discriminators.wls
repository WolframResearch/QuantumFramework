#!/usr/bin/env wolframscript
(* Decisive discriminators for findings 2 (eta<1) and 3 (V!=None).
   Key idea: run the corrected and the notebook integrators on the SAME random
   stream (common seed), so the Monte-Carlo noise cancels in the difference and
   dsys = |mean_notebook - mean_new| is a clean, MC-free measure of the bug. *)
out[s_] := WriteString[$Output, s, "\n"];
$ProgressReporting = False;
fails = 0;
check[label_, cond_, extra_ : ""] := (out["  [" <> If[TrueQ[cond], "PASS", fails++; "FAIL"] <> "] " <> label <> If[extra === "", "", "  (" <> extra <> ")"]];);

Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/ContinuousMeasurement/Manual StepLike Ito/manual-implementation-ito-qf-independent.wl"];

origR[rho_, L_, dt_, dw_, eta_] :=
  MapThread[Sqrt[#3] Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &, {L, dw, eta}];
origAlg[rho0_, Hmat_, Lmat_, etaSpec_, Vmat_, dt_, tf_, iddim_] :=
 Module[{dw, Heff, etam},
  etam = If[etaSpec === None, ConstantArray[1, Length[Lmat]], etaSpec];
  Heff = IdentityMatrix[iddim] - dt (I Hmat
       + .5 Sum[l . l + ConjugateTranspose[l] . l, {l, Lmat}]
       + If[Vmat === None, 0, .5 Sum[l . l + ConjugateTranspose[l] . l, {l, Vmat}]]);
  dw = RandomVariate[NormalDistribution[0, Sqrt[dt]], {Length[Lmat], Floor[tf/dt]}];
  FoldList[Module[{M, st, Ups},
     Ups = Total[Sqrt[etam] Lmat Sow@origR[#1, Lmat, dt, #2, etam]];
     M = Heff + .5 Ups . Ups + Ups;
     st = M . #1 . ConjugateTranspose[M]
        + If[Vmat === None, 0, dt Sum[l . #1 . ConjugateTranspose[l], {l, Vmat}]]
        + dt Sum[(1 - etam[[l]]) Lmat[[l]] . #1 . ConjugateTranspose[Lmat[[l]]], {l, Range[Length@Lmat]}];
     st/Tr[st]] &, rho0, Transpose[dw]]];

lindbladExpm[rho0_, H_, collapse_, tf_] := Module[{d = Length[H], idm, Lsup},
  idm = IdentityMatrix[d];
  Lsup = -I (KroneckerProduct[H, idm] - KroneckerProduct[idm, Transpose[H]])
    + Total[Map[KroneckerProduct[#, Conjugate[#]]
       - (1/2) KroneckerProduct[ConjugateTranspose[#] . #, idm]
       - (1/2) KroneckerProduct[idm, Transpose[ConjugateTranspose[#] . #]] &, collapse]];
  ArrayReshape[MatrixExp[N[Lsup] tf] . Flatten[N[rho0]], {d, d}]];
tdist[a_, b_] := Max[Abs[Flatten[a - b]]];

Jz = DiagonalMatrix[{1., 0., -1.}];
Jm3 = Sqrt[2.] {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};   (* spin-1 lowering, (Jm3)^2 = 2|-1><+1| != 0 *)
rho0 = DiagonalMatrix[{1., 0., 0.}];
dt = 0.01; tf = 2.5; K = 1500;

out["=== (4b) eta<1 qutrit: efficiency-factor fix (spin-1 J-, L^2 != 0) ==="];
Module[{gam = 0.8, eta = 0.2, H = 1.0 Jz, L, refL, aNew, aOrig, dNew, dOrig, dsys},
  L = {Sqrt[gam] Jm3};
  refL = lindbladExpm[rho0, H, L, tf];                       (* eta-independent unconditional ME *)
  SeedRandom[55]; aNew = Mean[Table[Last@\[ScriptCapitalT][rho0, H, L, {eta}, dt, tf], K]];
  SeedRandom[55]; aOrig = Mean[Table[Last@origAlg[rho0, H, L, {eta}, None, dt, tf, 3], K]];
  dNew = tdist[aNew, refL]; dOrig = tdist[aOrig, refL]; dsys = tdist[aOrig, aNew];
  out["    new vs Lindblad          = " <> ToString[dNew]];
  out["    notebook vs Lindblad     = " <> ToString[dOrig]];
  out["    systematic |orig - new|  = " <> ToString[dsys] <> "   <- MC-free measure of the bug"];
  check["new T matches eta-independent Lindblad (MC-limited)", dNew < 0.03, "dNew=" <> ToString[dNew]];
  check["bug is real and systematic (dsys >> new residual)", dsys > 5 dNew, "dsys=" <> ToString[dsys] <> " vs 5 dNew=" <> ToString[5 dNew]];
  check["corrected is closer to Lindblad than notebook", dOrig > dNew, "dOrig=" <> ToString[dOrig] <> " > dNew=" <> ToString[dNew]];
];

out[""];
out["=== (4c) V!=None qutrit: spurious-V^2 fix (spin-1 J- unmonitored) ==="];
Module[{gl = 0.3, gv = 0.8, H = 1.0 Jz, L, V, refL, aNew, aOrig, dNew, dOrig, dsys},
  L = {Sqrt[gl] Jm3}; V = {Sqrt[gv] Jm3};
  refL = lindbladExpm[rho0, H, Join[L, V], tf];
  SeedRandom[77]; aNew = Mean[Table[Last@\[ScriptCapitalT][rho0, H, L, None, V, dt, tf], K]];
  SeedRandom[77]; aOrig = Mean[Table[Last@origAlg[rho0, H, L, None, V, dt, tf, 3], K]];
  dNew = tdist[aNew, refL]; dOrig = tdist[aOrig, refL]; dsys = tdist[aOrig, aNew];
  out["    new vs Lindblad          = " <> ToString[dNew]];
  out["    notebook vs Lindblad     = " <> ToString[dOrig]];
  out["    systematic |orig - new|  = " <> ToString[dsys] <> "   <- MC-free measure of the bug"];
  check["new T matches Lindblad with V channel (MC-limited)", dNew < 0.03, "dNew=" <> ToString[dNew]];
  check["bug is real and systematic (dsys >> new residual)", dsys > 5 dNew, "dsys=" <> ToString[dsys] <> " vs 5 dNew=" <> ToString[5 dNew]];
  check["corrected is closer to Lindblad than notebook", dOrig > dNew, "dOrig=" <> ToString[dOrig] <> " > dNew=" <> ToString[dNew]];
];

out[""];
out["TOTAL FAILURES: " <> ToString[fails]];
