#!/usr/bin/env wolframscript
(* Verify + benchmark the prototypes in sde-schemes-prototype.wl against the
   density-matrix Rouchon-Ralph baseline (manual-implementation-ito-qf-independent.wl).
   QF-free: plain matrices/vectors, runs from a bare wolframscript.

   Sections:
     A  algebraic identity + numeric equivalence: \[ScriptCapitalP] reproduces \[ScriptCapitalT] at eta=1
     B  norm / positivity preservation of \[ScriptCapitalP] (any dt)
     C  ensemble mean of \[ScriptCapitalP] -> Lindblad master equation
     D  positivity trade: Euler-Maruyama \[ScriptCapitalE] goes non-positive where \[ScriptCapitalT]/\[ScriptCapitalP] stay PSD
     E  cost benchmark: \[ScriptCapitalP] (O(d^2)/step) vs \[ScriptCapitalT] (O(d^3)/step) across dimension
     F  Talay-Tubaro extrapolation: first-order Rouchon mean -> weak order 2, common noise *)

out[s_] := WriteString[$Output, s, "\n"];
sf[x_] := ToString[x, InputForm];              (* single-line exact form for tiny reals *)
nf[x_] := ToString[NumberForm[N[x], 4]];       (* short 4-sig-fig form for table magnitudes *)
$ProgressReporting = False;
fails = 0;
check[label_, cond_, extra_ : ""] := (
   out["  [" <> If[TrueQ[cond], "PASS", fails++; "FAIL"] <> "] " <> label <>
     If[extra === "", "", "  (" <> extra <> ")"]];
   TrueQ[cond]);

base = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/ContinuousMeasurement/Manual StepLike Ito/";
Get[base <> "manual-implementation-ito-qf-independent.wl"];
Get[base <> "sde-schemes-prototype.wl"];

(* ---- independent Lindblad reference (vectorized Liouvillian), from the baseline harness ---- *)
lindbladExpm[rho0_, H_, collapse_, tf_] := Module[{d = Length[H], idm, Lsup},
  idm = IdentityMatrix[d];
  Lsup = -I (KroneckerProduct[H, idm] - KroneckerProduct[idm, Transpose[H]])
    + Total[Map[KroneckerProduct[#, Conjugate[#]]
       - (1/2) KroneckerProduct[ConjugateTranspose[#] . #, idm]
       - (1/2) KroneckerProduct[idm, Transpose[ConjugateTranspose[#] . #]] &, collapse]];
  ArrayReshape[MatrixExp[N[Lsup] tf] . Flatten[N[rho0]], {d, d}]];

tdist[a_, b_] := Max[Abs[Flatten[a - b]]];
minEig[rho_] := Min[Re@Eigenvalues[rho]];
proj[psi_] := KroneckerProduct[psi, Conjugate[psi]];   (* |psi><psi| *)

X = PauliMatrix[1]; Y = PauliMatrix[2]; Z = PauliMatrix[3]; Jm = (X - I Y)/2;

(* spin-j operators in the basis m = j, j-1, ..., -j (index i <-> m = j-(i-1)) *)
spinOps[j_] := Module[{d = 2 j + 1, mm = Range[j, -j, -1], jz, jp, jmv},
  jz = DiagonalMatrix[N@mm];
  jp = SparseArray[{i_, k_} /; i == k - 1 :> N@Sqrt[j (j + 1) - mm[[k]] (mm[[k]] + 1)], {d, d}];
  jmv = Transpose[jp];                                  (* J- = (J+)^dag, real *)
  {jz, Normal[jmv]}];

out["=================================================================="];
out["SECTION A  \[ScriptCapitalP] reproduces \[ScriptCapitalT] at eta=1 (algebra + numerics)"];
out["=================================================================="];

(* A0: the exact algebraic identity normalize(M |psi><psi| M^dag) === |phi><phi|,
   phi = M psi/||M psi||. Proven 2x2 symbolic; then a random d=5 spot check to show it
   is dimension-general (the identity is M rho M^dag = (M psi)(M psi)^dag on rank one). *)
Module[{Msym, psisym, resid, Mr, pr},
  Msym = {{m11, m12}, {m21, m22}}; psisym = {p1, p2};
  resid = FullSimplify[
    With[{num = Msym . proj[psisym] . ConjugateTranspose[Msym]}, num/Tr[num]]
      - With[{ph = Msym . psisym}, proj[ph]/(Conjugate[ph] . ph)],
    Assumptions -> Element[{m11, m12, m21, m22, p1, p2}, Complexes]];
  check["algebra: N[M rho M^dag] === proj(M psi/||M psi||), rho=|psi><psi| (2x2 symbolic)",
     resid === {{0, 0}, {0, 0}}];
  SeedRandom[17]; Mr = RandomComplex[{-1 - I, 1 + I}, {5, 5}]; pr = RandomComplex[{-1 - I, 1 + I}, 5];
  check["same identity at d=5 (dimension-general, numeric)",
     tdist[With[{n = Mr . proj[pr] . ConjugateTranspose[Mr]}, n/Tr[n]], proj[Normalize[Mr . pr]]] < 1*^-13];
];

(* Weak order 1 of the ensemble mean (E_dW[N[M rho M^dag]] = rho + dt Lindblad(rho) +
   O(dt^2)) is the defining property of the Rouchon-Ralph CP map for any monitored L
   (their theorem). It is demonstrated NUMERICALLY, free of Monte-Carlo noise, by the
   multilevel coupled differences of Section F, and confirmed against the exact Lindblad
   reference in Section C. A closed-form symbolic average of the *normalized* nonlinear
   one-step map is deferred: the naive dt->eps^2, dW->eps g moment expansion mis-handles
   the 1/Tr normalization for some operators (it is exact for a lowering operator but
   carries a spurious O(dt) term for a raising one), so it is not a trustworthy proof and
   is not asserted here. *)

(* A1: numeric equivalence, shared seed => same dW. Includes multi-operator (nL=2)
   non-commuting cases, which alone exercise the double sum 1/2 Sum_{r,s} dy_r dy_s L_r L_s
   with L_r L_s != L_s L_r (dormant at nL=1). *)
Module[{cases, noncommQ},
  noncommQ[Ls_] := Length[Ls] >= 2 && Max[Abs[Flatten[Ls[[1]] . Ls[[2]] - Ls[[2]] . Ls[[1]]]]] > 1*^-6;
  cases = {
    <|"name" -> "qubit d=2, nL=1 (H=X, L=Sqrt2 Z)", "psi0" -> {1., 0.},
      "H" -> 1.0 X, "L" -> {Sqrt[2.0] Z}, "dt" -> 0.005, "tf" -> 10.|>,
    <|"name" -> "qubit d=2, nL=1 (J-)", "psi0" -> {1., 1.}/Sqrt[2.],
      "H" -> 0.5 (2.0 X + 1.0 Z), "L" -> {Sqrt[0.6] Jm}, "dt" -> 0.01, "tf" -> 8.|>,
    <|"name" -> "spin-1 d=3, nL=1 (Jz, J-)", "psi0" -> {1., 0., 0.},
      "H" -> 1.0 spinOps[1][[1]], "L" -> {Sqrt[0.5] spinOps[1][[2]]}, "dt" -> 0.01, "tf" -> 6.|>,
    <|"name" -> "qubit d=2, nL=2 non-commuting (J-, Z)", "psi0" -> {1., 0.},
      "H" -> 0.5 Z, "L" -> {Sqrt[0.5] Jm, Sqrt[0.3] Z}, "dt" -> 0.01, "tf" -> 6.|>,
    <|"name" -> "spin-1 d=3, nL=2 non-commuting (J-, Jz)", "psi0" -> {0., 1., 0.},
      "H" -> 0.4 spinOps[1][[1]], "L" -> {Sqrt[0.5] spinOps[1][[2]], Sqrt[0.2] spinOps[1][[1]]},
      "dt" -> 0.01, "tf" -> 5.|>};
  check["the nL=2 cases genuinely use non-commuting operators (double-sum ordering live)",
     AllTrue[Select[cases, Length[#["L"]] >= 2 &], noncommQ[#["L"]] &]];
  MapIndexed[
   Function[{ex, idx},
    Module[{seed = 30 + First[idx], dm, dmRec, ps, fm, dTraj, dRec, dFm},
     SeedRandom[seed]; {dm, dmRec} = Reap[\[ScriptCapitalT][proj[ex["psi0"]], ex["H"], ex["L"], ex["dt"], ex["tf"]]];
     SeedRandom[seed]; ps = Reap[\[ScriptCapitalP][ex["psi0"], ex["H"], ex["L"], ex["dt"], ex["tf"]]];
     SeedRandom[seed]; fm = \[ScriptCapitalF][ex["psi0"], ex["H"], ex["L"], ex["dt"], ex["tf"]];
     dTraj = Max[MapThread[tdist[#1, proj[#2]] &, {dm, ps[[1]]}]];
     dRec = tdist[dmRec[[1]], ps[[2, 1]]];
     dFm = Max[MapThread[tdist[#1, proj[#2]] &, {dm, fm}]];
     check[ex["name"] <> ": matvec \[ScriptCapitalP] matches density-matrix trajectory", dTraj < 1*^-11,
        "max|d|=" <> sf[dTraj]];
     check[ex["name"] <> ": form-M \[ScriptCapitalF] matches density-matrix trajectory", dFm < 1*^-11,
        "max|d|=" <> sf[dFm]];
     check[ex["name"] <> ": Sow-ed readout matches baseline", dRec < 1*^-11,
        "max|dy diff|=" <> sf[dRec]]]],
   cases];
];

out[""];
out["=================================================================="];
out["SECTION B  \[ScriptCapitalP] preserves the norm (=> PSD, trace 1) at any dt"];
out["=================================================================="];
With[{H = 0.5 (3.0 Y), L = {Sqrt[1.0] Jm}, psi0 = {1., 1.}/Sqrt[2.], dts = {0.01, 0.5, 5.}},
  Scan[
   Function[dtv,
    Module[{traj, worst},
     SeedRandom[7]; traj = \[ScriptCapitalP][psi0, H, L, dtv, 40.];
     worst = Max[Abs[(Conjugate[#] . #) - 1. & /@ traj]];
     check["dt=" <> sf[dtv] <> ": ||psi||=1 along trajectory", worst < 1*^-12,
        "max||psi||^2-1|=" <> sf[worst]]]],
   dts]];

out[""];
out["=================================================================="];
out["SECTION C  ensemble mean of \[ScriptCapitalP] -> Lindblad master equation"];
out["=================================================================="];
Module[{kappa = 0.6, dt = 0.02, tf = 6., K = 800, H, L, psi0, avg, ref, td},
  H = 0.5 X; L = {Sqrt[kappa] Jm}; psi0 = {1., 1.}/Sqrt[2.];
  SeedRandom[321];
  avg = Mean[Table[proj[Last@\[ScriptCapitalP][psi0, H, L, dt, tf]], K]];
  ref = lindbladExpm[proj[psi0], H, L, tf];
  td = tdist[avg, ref];
  check["mean of |psi><psi| (K=" <> sf[K] <> ") ~ Lindblad", td < 0.03,
     "max|d|=" <> sf[td] <> " (MC 1/Sqrt[K]=" <> sf[N[1/Sqrt[K]]] <> ")"];
];
(* exactly-solvable closed form, independent of the vectorized-Liouvillian reference:
   |0> under L = Sqrt[gamma] |1><0|, H=0 (amplitude damping) has rho00(t) = e^{-gamma t}. *)
Module[{gamma = 0.7, tf = 1.3, K = 1500, psi0 = {1., 0.}, avg, analytic},
  SeedRandom[99];
  avg = Mean[Table[proj[Last@\[ScriptCapitalP][psi0, 0. Z, {Sqrt[gamma] Jm}, 0.01, tf]], K]];
  analytic = Exp[-gamma tf];
  check["ensemble mean rho00(t) matches analytic e^{-gamma t} (amplitude damping, closed form)",
     Abs[Re[avg[[1, 1]]] - analytic] < 0.03,
     "|d|=" <> sf[Abs[Re[avg[[1, 1]]] - analytic]] <> " analytic=" <> sf[analytic]];
];

out[""];
out["=================================================================="];
out["SECTION D  positivity trade: Euler-Maruyama \[ScriptCapitalE] vs Rouchon \[ScriptCapitalP]/\[ScriptCapitalT]"];
out["=================================================================="];
With[{H = 0.5 X, L = {Sqrt[1.0] Jm}, psi0 = {1., 0.}, seed = 4, dts = {0.02, 0.2, 0.5}},
 With[{r0 = proj[psi0]},
  out["    minimum eigenvalue over the trajectory (negative = left the physical set):"];
  Scan[
   Function[dtv,
    Module[{eE, eT},
     SeedRandom[seed]; eE = Min[minEig /@ \[ScriptCapitalE][r0, H, L, dtv, 20.]];
     SeedRandom[seed]; eT = Min[minEig /@ (proj /@ \[ScriptCapitalP][psi0, H, L, dtv, 20.])];
     out["      dt=" <> sf[dtv] <> ":  Euler-Maruyama min eig = " <> sf[eE]
        <> " ;  Rouchon min eig = " <> sf[eT]]]],
   dts];
  Module[{eE, eT},
   SeedRandom[seed]; eE = Min[minEig /@ \[ScriptCapitalE][r0, H, L, 0.5, 20.]];
   SeedRandom[seed]; eT = Min[minEig /@ (proj /@ \[ScriptCapitalP][psi0, H, L, 0.5, 20.])];
   check["Rouchon stays PSD at dt=0.5", eT > -1*^-9, "min eig=" <> sf[eT]];
   check["Euler-Maruyama goes non-positive at dt=0.5 (structure lost)", eE < -1*^-3,
      "min eig=" <> sf[eE]]]]];
(* the trade is structure, not accuracy: the Euler-Maruyama ensemble MEAN still tracks the
   Lindblad master equation at small dt, even as its individual states leave the physical set. *)
Module[{H = 0.5 X, L = {Sqrt[0.6] Jm}, r0 = proj[{1., 1.}/Sqrt[2.]], dt = 0.01, tf = 4., K = 400, ref, emAvg, td},
  ref = lindbladExpm[r0, H, L, tf];
  SeedRandom[71]; emAvg = Mean[Table[Last@\[ScriptCapitalE][r0, H, L, dt, tf], K]];
  td = tdist[emAvg, ref];
  check["Euler-Maruyama mean still tracks Lindblad at small dt (it loses structure, not accuracy)",
     td < 0.06, "max|d|=" <> sf[td] <> " (MC~" <> sf[N[1/Sqrt[K]]] <> ")"];
];

out[""];
out["=================================================================="];
out["SECTION E  cost: \[ScriptCapitalP] (state vector) vs \[ScriptCapitalT] (density matrix) across dimension"];
out["=================================================================="];
SetAttributes[bestOf, HoldFirst];
bestOf[expr_, n_] := Min[Table[ClearSystemCache[]; First@AbsoluteTiming[expr], {n}]];
(* time the two public pure-state maps of the deliverable, matvec \[ScriptCapitalP]
   (O(d^2)) and form-M \[ScriptCapitalF] (one matmul), against the density matrix. *)
Module[{js = {1/2, 3/2, 7/2, 15/2, 31/2}, gamma = 0.3, dt = 0.01, tf = 10., reps = 5},
  out["    spin-j: H=Jz, L=Sqrt[0.3] J-, 1000 steps, min over " <> sf[reps] <> " (ClearSystemCache)."];
  out["    d capped at 32: the dense full trajectory is O(d^2) memory (a d=64 trajectory is"];
  out["    ~65 MB vs the ket's ~1 MB) and its timing turns memory-bound and noisy past that;"];
  out["    the pure-state advantage only grows beyond, as O(d) arithmetic and O(d) memory."];
  out["      " <> StringPadRight["d", 6] <> StringPadRight["dense(s)", 13] <>
     StringPadRight["pure mv(s)", 13] <> StringPadRight["pure fM(s)", 13] <>
     StringPadRight["mv x", 11] <> "fM x"];
  Scan[
   Function[j,
    Module[{d = 2 j + 1, jz, jmv, H, L, psi0, r0, tD, tMv, tFm},
     {jz, jmv} = spinOps[j];
     H = jz; L = {Sqrt[gamma] jmv}; psi0 = UnitVector[d, 1]; r0 = proj[psi0];
     tD  = bestOf[\[ScriptCapitalT][r0, H, L, dt, tf], reps];
     tMv = bestOf[\[ScriptCapitalP][psi0, H, L, dt, tf], reps];
     tFm = bestOf[\[ScriptCapitalF][psi0, H, L, dt, tf], reps];
     out["      " <> StringPadRight[sf[d], 6] <> StringPadRight[nf[tD], 13] <>
        StringPadRight[nf[tMv], 13] <> StringPadRight[nf[tFm], 13] <>
        StringPadRight[nf[tD/tMv] <> "x", 11] <> nf[tD/tFm] <> "x"]]],
   js];
];

out[""];
out["SECTION F  weak order: base Rouchon is order 1; Talay-Tubaro extrapolation removes the leading O(dt) error"];
out["=================================================================="];
(* Rouchon (eta=1, no V) folded over externally supplied dW (rows = steps, cols = ops);
   density-matrix form, reusing the deliverable's rouchonHeffLcorr for the constants. *)
rouchonNoise[rho0_, H_, Ls_, dt_, dwRows_] := Module[
  {heff, lcorr, lsum, step},
  {heff, lcorr} = rouchonHeffLcorr[H, Ls, dt];
  lsum = # + ConjugateTranspose[#] & /@ Ls;
  step = Function[{rho, dwv},
    With[{dy = MapThread[Re[Tr[#1 . rho]] dt + #2 &, {lsum, dwv}]},
     With[{s = dy . Ls},
      With[{m = heff + s + (1/2) s . s - lcorr},
       With[{num = m . rho . ConjugateTranspose[m]}, num/Tr[num]]]]]];
  Fold[step, rho0, dwRows]];

(* pin rouchonNoise to the baseline before it carries the whole Section F conclusion:
   fed the dW that \[ScriptCapitalT] draws under a shared seed, it must reproduce it. *)
Module[{H = 0.5 Z, L = {Sqrt[0.4] Jm, Sqrt[0.25] X}, rho0 = proj[{1., 1.}/Sqrt[2.]], dt = 0.02, ns = 200, dwm, rn, rt},
  SeedRandom[808]; dwm = RandomVariate[NormalDistribution[0, Sqrt[dt]], {Length[L], ns}];
  rn = rouchonNoise[rho0, H, L, dt, Transpose[dwm]];
  SeedRandom[808]; rt = Last@\[ScriptCapitalT][rho0, H, L, dt, ns dt];
  check["rouchonNoise (the Section F integrator) matches baseline T on a shared dW",
     tdist[rn, rt] < 1*^-11, "max|d|=" <> sf[tdist[rn, rt]]];
];

(* Weak error e(dt) = E[f(rho_dt)] - f_exact has an expansion c1 dt + c2 dt^2 + ...
   Comparing a mean to f_exact is Monte-Carlo-floor limited. Instead drive one
   common Brownian path at dt/8 and coarsen it to dt/4, dt/2, dt; every reported
   quantity is a per-path combination whose coefficients sum to zero, so the O(1)
   trajectory value cancels and the estimator is bias- (not MC-) limited:
     base coupled difference   D_k  = E[f_dt - f_{dt/2}]     ~ c1 dt/2  (halves: ratio 2)
     Richardson  R(dt) = 2 f_{dt/2} - f_dt ;  DR_k = E[R(dt) - R(dt/2)] ~ dt^2 (quarters: ratio 4) *)
Module[{H, L, rho0, tf = 3.2, K = 4000, obs, exact, nc = 20, dt, data,
        mf, D1, D2, D3, DR1, DR2},
  H = 0.5 X; L = {Sqrt[0.6] Jm}; rho0 = proj[{1., 1.}/Sqrt[2.]];
  obs[rho_] := Re[rho[[1, 1]]];                              (* rho00, a smooth observable *)
  exact = obs[lindbladExpm[rho0, H, L, tf]];
  dt = tf/nc;                                                (* coarsest dt = 0.16 *)
  out["    observable rho00 at t=" <> sf[tf] <> ", exact = " <> sf[exact]
     <> " ; K=" <> sf[K] <> " common-path trajectories, dt in {0.16,0.08,0.04,0.02}"];
  SeedRandom[2024];
  data = Table[
    Module[{dw8, dw4, dw2, dw1},
     dw8 = RandomVariate[NormalDistribution[0, Sqrt[dt/8]], {8 nc, 1}];
     dw4 = Total /@ Partition[dw8, 2];
     dw2 = Total /@ Partition[dw4, 2];
     dw1 = Total /@ Partition[dw2, 2];
     obs[rouchonNoise[rho0, H, L, #1, #2]] & @@@
       {{dt, dw1}, {dt/2, dw2}, {dt/4, dw4}, {dt/8, dw8}}],  (* {f1,f2,f4,f8}, common path *)
    K];
  mf = Mean /@ Transpose[data];                              (* raw means {f1,f2,f4,f8} *)
  D1 = Mean[(#[[1]] - #[[2]]) &  /@ data];                   (* E[f_dt   - f_{dt/2}] *)
  D2 = Mean[(#[[2]] - #[[3]]) &  /@ data];                   (* E[f_{dt/2}- f_{dt/4}] *)
  D3 = Mean[(#[[3]] - #[[4]]) &  /@ data];                   (* E[f_{dt/4}- f_{dt/8}] *)
  DR1 = Mean[((2 #[[2]] - #[[1]]) - (2 #[[3]] - #[[2]])) & /@ data];  (* E[R(dt)  - R(dt/2)] *)
  DR2 = Mean[((2 #[[3]] - #[[2]]) - (2 #[[4]] - #[[3]])) & /@ data];  (* E[R(dt/2)- R(dt/4)] *)
  out["    raw weak bias E[f]-exact at the four dt (MC-floor limited, ~" <>
     sf[N[Sqrt[exact (1 - exact)/K]]] <> "): " <> sf[mf - exact]];
  out["    base coupled diffs   D = {" <> sf[D1] <> ", " <> sf[D2] <> ", " <> sf[D3] <> "}"];
  out["      ratios D1/D2, D2/D3 (weak order 1 => ~2): " <> sf[{D1/D2, D2/D3}]];
  out["    Richardson R(dt)=2 f_{dt/2}-f_dt ; coupled diffs DR = {" <> sf[DR1] <> ", " <> sf[DR2] <> "}"];
  out["      |DR1|/|D1| = " <> sf[Abs[DR1/D1]] <> " , |DR2|/|D2| = " <> sf[Abs[DR2/D2]]
     <> "  (leading O(dt) term removed; O(dt^2) residual sits below the MC floor here)"];
  check["base Rouchon mean is weak order 1 (coupled diff halves: D1/D2 ~ 2)", 1.6 < D1/D2 < 2.6,
     "D1/D2=" <> sf[D1/D2]];
  check["base order 1 confirmed at the next level too (D2/D3 ~ 2)", 1.6 < D2/D3 < 2.6,
     "D2/D3=" <> sf[D2/D3]];
  check["Talay-Tubaro removes the leading O(dt) error at the coarse level (|DR1| << |D1|)",
     Abs[DR1] < 0.1 Abs[D1], "|DR1/D1|=" <> sf[Abs[DR1/D1]]];
  check["extrapolant much more accurate than base at the next level too (|DR2| < |D2|)",
     Abs[DR2] < 0.5 Abs[D2], "|DR2/D2|=" <> sf[Abs[DR2/D2]]];
];

out[""];
out["=================================================================="];
out["TOTAL FAILURES: " <> sf[fails]];
out["=================================================================="];
If[fails > 0, Exit[1]];
