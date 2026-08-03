(* ============================================================================
   Consolidated verification: NDSolve/WhenEvent Euler-Maruyama reproduction of
   TransmonLab's ItoProcess Cartesian-Bloch quantum-trajectory SDE.
   Layers: (S) symbolic invariants incl. exactly-solvable and Zeno limits,
   (A) real-rate Ito<->Lindblad, (B) dimensionless EM==Ito==Lindblad with r<=1,
   (E) transverse-sector and eta<1 mixed-state cases, (G) Zeno suppression.
   No Quiet, no Print.  Run: wolframscript -file verify_all.wls
   ============================================================================ *)
SetDirectory["/private/tmp/claude-10838/-Users-mohammadb-Documents-GitHub-QuantumMachines-TransmonLab--claude-worktrees-compassionate-banach-e7090b/ae3164c5-9a0b-4ab5-bd0b-7adb86cf987e/scratchpad"];
Get["common.wl"];

report = <||>;
add[k_, v_] := (report[k] = v; Echo[v, k <> ": "]);
ok[vt_] := TrueQ[vt["Outcome"] == "Success"];
within[val_, target_, ntraj_] := Abs[val - target] < 4/Sqrt[ntraj];   (* 4-sigma MC bound *)

(* ---------- (S) SYMBOLIC LAYER (arbitrary symbolic rates) ---------- *)
Clear[dep, g1, wx, x0, y0, z0, X, Y, Z, x, y, z, gci, gba, gd, gphi, eta, th, ph, xx, yy, zz];

(* S1: affine mean ODE closes; DSolve residual 0 *)
meanSol = DSolve[{X'[t] == -dep X[t], Y'[t] == -dep Y[t] - wx Z[t],
    Z'[t] == g1 (1 - Z[t]) + wx Y[t], X[0] == x0, Y[0] == y0, Z[0] == z0},
   {X, Y, Z}, t] // First;
S1 = VerificationTest[
   FullSimplify[{D[X[t] /. meanSol, t] - (-dep (X[t] /. meanSol)),
      D[Y[t] /. meanSol, t] - (-dep (Y[t] /. meanSol) - wx (Z[t] /. meanSol)),
      D[Z[t] /. meanSol, t] - (g1 (1 - (Z[t] /. meanSol)) + wx (Y[t] /. meanSol))}],
   {0, 0, 0}, TestID -> "S1 affine-mean-ODE closed form"];

(* S2: dW-coefficient of d(r^2) = 2 sqrt(gci) z (1-r^2) -> 0 at r=1 *)
bV = diffV[{x, y, z}, {gci, gba, gd, gphi, g1, wx}];
S2 = VerificationTest[
   FullSimplify[2 ({x, y, z} . bV) - 2 Sqrt[gci] z (1 - (x^2 + y^2 + z^2))],
   0, TestID -> "S2 r^2 noise vanishes at r=1"];

(* S3: Ito generator L[r^2] <= 0 on the sphere for efficiency eta<=1 *)
xyz = {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]};
Lr2sphere = FullSimplify[(2 (xyz . driftV[xyz, {gci, gba, gd, gphi, g1, wx}])
    + diffV[xyz, {gci, gba, gd, gphi, g1, wx}] . diffV[xyz, {gci, gba, gd, gphi, g1, wx}]) /.
   gba -> 2 eta gd - gci, Assumptions -> {0 <= th <= Pi, 0 <= ph < 2 Pi}];
S3 = VerificationTest[
   Reduce[Lr2sphere > 0 && 0 <= th <= Pi && g1 >= 0 && gphi >= 0 && gd >= 0
      && gci >= 0 && 0 <= eta <= 1, {th}],
   False, TestID -> "S3 r<=1 invariant: no sphere point has L[r^2]>0"];

(* S4: factored 2x2-Lindblad Bloch rates == affine drift (independent machinery) *)
S4 = VerificationTest[
   Simplify[lindblad2x2Bloch[{xx, yy, zz}, {0, 0, gd, gphi, g1, wx}]
      - driftV[{xx, yy, zz}, {0, 0, gd, gphi, g1, wx}],
    {xx, yy, zz, wx, g1, gphi, gd} \[Element] Reals],
   {0, 0, 0}, TestID -> "S4 2x2-Lindblad == affine drift"];

(* S5: 2x2 Lindblad generator invariants -- trace preservation & Hermiticity,
   and rho >= 0  <=>  r<=1 (positivity IS the Bloch-ball invariant). *)
drhoSym = lindblad2x2Gen[(IdentityMatrix[2] + xx PauliMatrix[1] + yy PauliMatrix[2]
     + zz PauliMatrix[3])/2, {0, 0, gd, gphi, g1, wx}];
S5a = VerificationTest[Simplify[Tr[drhoSym], Reals], 0, TestID -> "S5a Tr[drho]=0 (trace preserved)"];
S5b = VerificationTest[
   FullSimplify[ComplexExpand[drhoSym - ConjugateTranspose[drhoSym]]],
   {{0, 0}, {0, 0}}, TestID -> "S5b drho Hermitian"];
(* rho positivity IS the Bloch-ball invariant: Tr=1, Det=(1-r^2)/4 >= 0 iff r<=1 *)
rhoP = (IdentityMatrix[2] + xx PauliMatrix[1] + yy PauliMatrix[2] + zz PauliMatrix[3])/2;
S5c = VerificationTest[
   FullSimplify[{Tr[rhoP] - 1, Det[rhoP] - (1 - (xx^2 + yy^2 + zz^2))/4}],
   {0, 0}, TestID -> "S5c rho: Tr=1, Det=(1-r^2)/4 (>=0 iff r<=1)"];

(* S6: exactly-solvable special case -- pure measurement (wx=g1=gphi=0):
   <z> constant (QND), <x>,<y> decay at gd. *)
pmBloch = #[t] & /@ lindbladBloch[{gci, gba, gd, 0, 0, 0}, {x0, y0, z0}];
S6 = VerificationTest[FullSimplify[pmBloch - {x0 Exp[-gd t], y0 Exp[-gd t], z0}],
   {0, 0, 0}, TestID -> "S6 pure-measurement QND special case"];

(* S7: Zeno suppression -- take the (y,z) block of the DRIFT Jacobian FROM driftV
   (so a drift typo would be caught), then its slow eigenvalue in the strong-
   measurement limit dep->Infinity is -(g1 + wx^2/dep): the Rabi rate is suppressed
   by the measurement (Quantum Zeno effect). *)
Clear[dep, uu, vv];
jacYZ = D[driftV[{0, uu, vv}, {gci, gba, gd, gphi, g1, wx}][[{2, 3}]], {{uu, vv}}] /.
   gd -> dep - g1/2 - gphi;                        (* = {{-dep,-wx},{wx,-g1}} *)
lamSlow = SelectFirst[Eigenvalues[jacYZ],
   FreeQ[Limit[#, dep -> Infinity], DirectedInfinity | Infinity] &];  (* the slow root -> -g1 *)
zenoSeries = Normal@Series[lamSlow, {dep, Infinity, 1}];
S7 = VerificationTest[Simplify[zenoSeries - (-g1 - wx^2/dep)], 0,
   TestID -> "S7 Zeno: (y,z) drift-Jacobian slow eigenvalue -> -(g1 + wx^2/dep)"];

add["S symbolic all pass?", AllTrue[{S1, S2, S3, S4, S5a, S5b, S5c, S6, S7}, ok]];

(* ---------- (A) REAL TransmonLab equivalence-case rates ---------- *)
ratesA = {11737.880001524305, 0., 5868.9400007621525,
          33333.33333333333, 33333.333333333336, 2 Pi*10^6};
initA = N[{Sin[0.001], 0., Cos[0.001]}]; tfA = 1.*^-5;
{gciA, gbaA, gdA, gphiA, g1A, wxA} = ratesA;
add["A eta=(gci+gba)/(2gd)", (gciA + gbaA)/(2 gdA)];
add["A dep-(1/T2+gd)", (g1A/2 + gphiA + gdA) - (1/20.*^-6 + gdA)];
lindZA = lindbladZfun[ratesA, initA][tfA];         (* DSolve, exact *)
lind2xA = lindblad2x2Zfun[ratesA, initA][tfA];     (* MatrixExp propagator, exact *)
add["A Lindblad DSolve vs Liouvillian-exp |diff|", Abs[lindZA - lind2xA]];
add["A Lindblad z(tf)", lindZA];
itoA = itoEnsemble[ratesA, initA, tfA, tfA/500, 1000, 42];
add["A Ito meanZ(tf) [nSteps=500,N=1000]", Last[itoA["meanZ"]]];
add["A |Ito-Lindblad| < 4/sqrt(N)?", within[Last[itoA["meanZ"]], lindZA, 1000]];

(* D2: deterministic explicit-EM converges only O(dt) on the fast-Rabi real case.
   Report the relative error honestly (still converging, not converged), and that
   it shrinks with dt.  Below the stability threshold noProj blows past r=1. *)
ratesDet = {0., 0., 0., gphiA, g1A, wxA};
lindDet = lindbladZfun[ratesDet, initA][tfA];
add["A det stability threshold nSteps*", wxA^2 tfA/((g1A/2 + gphiA) + g1A)];
add["A det EM noProj z(tf) @1000 (below thr, blows up)",
   Last@First@emPath[ratesDet, initA, tfA, tfA/1000, False]];
detErr = Table[
   Abs[Last@First@emPath[ratesDet, initA, tfA, tfA/ns, True] - lindDet]/Abs[lindDet],
   {ns, {8000, 16000, 32000}}];
add["A det EM proj REL ERROR @{8000,16000,32000} (O(dt), shrinking, NOT converged)", detErr];
add["A det EM rel error monotonically shrinks?", detErr[[1]] > detErr[[2]] > detErr[[3]]];

(* ---------- equivalence-case scaffold (DRY): Lindblad + EM-proj + native Ito ---------- *)
equivCase[rates_, init_, tf_, dt_, n_, seed_] := <|
  "lind" -> lindbladZfun[rates, init][tf],
  "em"   -> emEnsembleVec[rates, init, tf, dt, n, seed, True],
  "ito"  -> itoEnsemble[rates, init, tf, dt, n, seed]|>;
emIto[c_, n_] := within[Last[c["em"]["meanZ"]], c["lind"], n] &&
   within[Last[c["ito"]["meanZ"]], c["lind"], n];
stdFin[a_] := StandardDeviation[a["zTrajs"][[All, -1]]];

(* ---------- (B) DIMENSIONLESS: EM==Ito==Lindblad, r<=1 ---------- *)
ratesB = {3., 0., 1.5, 0., 1., 3.}; initB = {0., 0., 1.}; tfB = 3.; nB = 1000;
cB = equivCase[ratesB, initB, tfB, tfB/600, nB, 2024];
emBnp = emEnsembleVec[ratesB, initB, tfB, tfB/150, nB, 2024, False];
add["B {Lindblad, EMproj, Ito} meanZ(tf)", {cB["lind"], Last[cB["em"]["meanZ"]], Last[cB["ito"]["meanZ"]]}];
add["B EM & Ito reproduce Lindblad < 4/sqrt(N)?", emIto[cB, nB]];
add["B EMproj maxr <= 1 (r<=1 every path)?", cB["em"]["maxr"] <= 1 + 1.*^-9];
add["B EMnoProj overshoots r>1?", emBnp["maxr"] > 1.05];

(* ---------- (D) SECOND MOMENT (teeth): the affine-drift MEAN is diffusion-blind
   -- a wrong diffusion still matches the Lindblad mean -- so also check Std[z(tf)].
   EM reproduces the native-Ito Std; a half-strength diffusion keeps the mean but
   breaks the Std -- the teeth the diffusion-blind mean lacks. ---------- *)
halfProcB = ItoProcess[{driftV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, ratesB],
    List /@ (diffV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, ratesB]/2),
    {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}}, {{\[FormalX], \[FormalY], \[FormalZ]}, initB}, t];
halfStdB = StandardDeviation[(BlockRandom[SeedRandom[2024];
   RandomFunction[halfProcB, {0., tfB, tfB/600}, nB,
     Method -> "StochasticRungeKuttaScalarNoise"]]["ValueList"])[[All, -1, 3]]];
add["D Std[z(tf)] {Ito, EMproj, half-diffusion}", {stdFin[cB["ito"]], stdFin[cB["em"]], halfStdB}];
add["D EM reproduces Ito 2nd moment AND half-diffusion caught (teeth)?",
   Abs[stdFin[cB["em"]] - stdFin[cB["ito"]]] < 0.06 && Abs[halfStdB - stdFin[cB["ito"]]] > 0.1];

(* ---------- (E) transverse (equatorial) init and eta<1 (mixed state) ---------- *)
(* E1: equatorial init {1,0,0} exercises the <x>,<y> sector (never near the pole). *)
cE1 = equivCase[ratesB, {1., 0., 0.}, tfB, tfB/600, nB, 77];
add["E1 equatorial: EM & Ito reproduce Lindblad, EMproj r<=1?",
   emIto[cE1, nB] && cE1["em"]["maxr"] <= 1 + 1.*^-9];

(* E2: efficiency eta<1 (gci+gba < 2 gd) => MIXED trajectories, Bloch radius
   strictly inside the ball -- the invariant is tested off its boundary.
   EM-proj still holds r<=1 exactly; Ito-SRK overshoots slightly (not norm-preserving). *)
ratesEta = {1., 0., 1.5, 0., 1., 3.};   (* eta = 1/3 *)
cE2 = equivCase[ratesEta, initB, tfB, tfB/600, nB, 88];
add["E2 {eta, Ito mean final Bloch radius (mixed<1)}",
   {(ratesEta[[1]] + ratesEta[[2]])/(2 ratesEta[[3]]), Mean[cE2["ito"]["finalR"]]}];
add["E2 eta<1: EM & Ito reproduce Lindblad; mixed r<1; EMproj r<=1?",
   emIto[cE2, nB] && Mean[cE2["ito"]["finalR"]] < 0.99 && cE2["em"]["maxr"] <= 1 + 1.*^-9];

(* ---------- (G) ZENO: strong measurement suppresses coherent Rabi ---------- *)
(* dep >> wx: (y,z) subspace overdamped (discriminant>0, no oscillation), slow
   longitudinal rate ~ g1 + wx^2/dep.  Ensemble reproduces the (Zeno) Lindblad. *)
ratesZ = {20., 0., 10., 0., 0.5, 2.}; initZ = {0., 0., 1.}; tfZ = 4.; nZ = 800;
depZ = depOf[ratesZ[[5]], ratesZ[[4]], ratesZ[[3]]];    (* = 10.25 *)
cZ = equivCase[ratesZ, initZ, tfZ, tfZ/1600, nZ, 5];
add["G discriminant>0 (overdamped/Zeno, no oscillation)?",
   (depZ - ratesZ[[5]])^2 - 4 ratesZ[[6]]^2 > 0];
jacZ = D[driftV[{0, uu, vv}, ratesZ][[{2, 3}]], {{uu, vv}}];   (* (y,z) Jacobian from driftV *)
add["G exact slow eigenvalue vs leading Zeno -(g1+wx^2/dep)",
   {Max[Eigenvalues[jacZ]], -(ratesZ[[5]] + ratesZ[[6]]^2/depZ)}];
add["G strong-meas ensemble reproduces Lindblad, EMproj r<=1?",
   emIto[cZ, nZ] && cZ["em"]["maxr"] <= 1 + 1.*^-9];

(* ---------- SUMMARY ---------- *)
add["ALL SYMBOLIC PASS", AllTrue[{S1, S2, S3, S4, S5a, S5b, S5c, S6, S7}, ok]];
add["ALL NUMERIC PASS", And[
   within[Last[itoA["meanZ"]], lindZA, 1000],
   detErr[[1]] > detErr[[2]] > detErr[[3]],
   emIto[cB, nB], cB["em"]["maxr"] <= 1 + 1.*^-9, emBnp["maxr"] > 1.05,
   Abs[stdFin[cB["em"]] - stdFin[cB["ito"]]] < 0.06, Abs[halfStdB - stdFin[cB["ito"]]] > 0.1,
   emIto[cE1, nB], cE1["em"]["maxr"] <= 1 + 1.*^-9,
   emIto[cE2, nB], Mean[cE2["ito"]["finalR"]] < 0.99, cE2["em"]["maxr"] <= 1 + 1.*^-9,
   emIto[cZ, nZ], cZ["em"]["maxr"] <= 1 + 1.*^-9]];
Export["verify_all_result.m", report];
"verify_all-done"
