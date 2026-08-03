#!/usr/bin/env wolframscript
(* Benchmark: notebook algorithm (matrix mode) vs the new QF-independent T.
   Both are pure WL (no FFI), so this is a straight kernel comparison.
   min-over-repeats with ClearSystemCache before each (warm-cache lies without it). *)

out[s_] := WriteString[$Output, s, "\n"];
$ProgressReporting = False;
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/ContinuousMeasurement/Manual StepLike Ito/manual-implementation-ito-qf-independent.wl"];

X = PauliMatrix[1]; Y = PauliMatrix[2]; Z = PauliMatrix[3]; Jm = (X - I Y)/2;
id2 = IdentityMatrix[2];

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

SetAttributes[bestOf, HoldFirst];
bestOf[expr_, n_] := Min[Table[ClearSystemCache[]; First@AbsoluteTiming[expr], {n}]];

SetAttributes[run, HoldRest];
run[label_, oExpr_, nExpr_, reps_] := Module[{to, tn},
  to = bestOf[oExpr, reps]; tn = bestOf[nExpr, reps];
  out["  " <> label];
  out["     notebook algorithm : " <> ToString[to] <> " s"];
  out["     new T              : " <> ToString[tn] <> " s"];
  out["     speedup            : " <> ToString[to/tn] <> "x"];
];

out["=== single qubit, 1 operator ==="];
run["d=2, 4000 steps (dt=0.005, tf=20)",
  origAlg[{{1., 0.}, {0., 0.}}, 1.0 X, {Sqrt[2.0] Z}, None, None, 0.005, 20., 2],
  \[ScriptCapitalT][{{1., 0.}, {0., 0.}}, 1.0 X, {Sqrt[2.0] Z}, 0.005, 20.], 5];

run["d=2, 40000 steps (dt=0.01, tf=400) [J-]",
  origAlg[{{0.5, 0.5}, {0.5, 0.5}}, 0.5 (3.5 Y), {Sqrt[0.7] Jm}, None, None, 0.01, 400., 2],
  \[ScriptCapitalT][{{0.5, 0.5}, {0.5, 0.5}}, 0.5 (3.5 Y), {Sqrt[0.7] Jm}, 0.01, 400.], 3];

out["=== two qubits, 2 operators (paper example, d=4) ==="];
Module[{H, L, r0},
  H = (1./2) (KroneckerProduct[X, id2] + KroneckerProduct[id2, X]) + 0.01 KroneckerProduct[Z, Z];
  L = {Sqrt[0.01] KroneckerProduct[Z, id2], Sqrt[0.01] KroneckerProduct[id2, Z]};
  r0 = IdentityMatrix[4]/4.;
  run["d=4, 5000 steps (dt=0.01, tf=50), eta=0.85",
   origAlg[r0, H, L, {0.85, 0.85}, None, 0.01, 50., 4],
   \[ScriptCapitalT][r0, H, L, {0.85, 0.85}, 0.01, 50.], 4];
];

out["=== the notebook's actual two-point run: new T alone at 400,000 steps ==="];
Module[{t},
  t = bestOf[\[ScriptCapitalT][{{1., 0.}, {0., 0.}}, 0.5 (3.5 Y), {Sqrt[0.7] Jm}, 0.01, 4000.], 3];
  out["  new T, d=2, 400000 steps (dt=0.01, tf=4000): " <> ToString[t] <> " s"];
];
