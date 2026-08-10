(* Validation battery for engine-2d3d.wl against the kernel-verified 1D essay
   (../measurement-in-a-potential.md) and closed forms.

   Run:  wolframscript -file validation.wls
   Every line prints value :: PASS/FAIL against the tolerance stated in its label;
   the last line counts failures. Statistical checks carry loose tolerances and fixed
   seeds. Section D3 asserts that a deliberately under-resolved grid IS visibly wrong. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "engine-2d3d.wl"}]];

w[s_] := WriteString[$Output, s, "\n"];
fails = 0;
chk[label_, val_, tol_] := w[label <> " :: " <> ToString[val, InputForm] <> " :: " <>
    If[TrueQ[Abs[val] <= tol], "PASS", fails++; "FAIL"]];
chkRatio[label_, val_, tol_] := chk[label, val - 1, tol];
note[label_, val_] := w[label <> " :: " <> ToString[val, InputForm]];

(* The 1D steady moments, one symbolic definition specialized everywhere below. Section D0
   DERIVES this object from the Riccati and selects the physical branch, so these lines
   are the named form of a derived result, not an independent transcription: the free and
   harmonic steady variances are its om -> 0 and general specializations. *)
sxSSsym[lam_, om_, ms_, hs_] := Sqrt[(-ms om^2 + Sqrt[ms^2 om^4 + 4 lam^2 hs^2])/(8 lam^2 ms)];
modeSSsym[lam_, om_, ms_, hs_] := With[{sx = sxSSsym[lam, om, ms, hs]},
   With[{cc = 2 lam ms sx^2}, {sx, cc, ms^2 om^2 sx + 4 lam ms sx cc}]];
freeSSvar[lam_] := sxSSsym[lam, 0, mass, hbar];
harmSSvar[lam_, om_] := sxSSsym[lam, om, mass, hbar];

w["==== A: transform conventions (manufactured plane wave, non-square grids) ===="]
nx = 16; lx = 13.; ny = 8; ly = 7.; nz = 4; lz = 5.;
gA = grid2[{nx, ny}, {lx, ly}];
{ax, ay} = {3, -2};
psiPW = Exp[I ((2. Pi/lx) ax gA["X"] + (2. Pi/ly) ay gA["Y"])]/Sqrt[nx ny];
cPW = ft[psiPW];
ix = First@FirstPosition[Round[gA["px"] lx/(2. Pi)], ax];
iy = First@FirstPosition[Round[gA["py"] ly/(2. Pi)], ay];
chk["A1 2D plane wave lands on its own bin, weight - 1 (tol 1e-12)", Abs[cPW[[ix, iy]]]^2 - 1, 1.*^-12]
chk["A2 off-bin weight (tol 1e-12)", Total[Abs[cPW]^2, 2] - Abs[cPW[[ix, iy]]]^2, 1.*^-12]
chk["A3 Parseval, flat 2-norm ratio - 1 (tol 1e-12)", stateNorm[cPW]/stateNorm[psiPW] - 1, 1.*^-12]
chk["A4 round trip (tol 1e-12)", Max@Abs[ift[cPW] - psiPW], 1.*^-12]
gA3 = grid3[{nx, ny, nz}, {lx, ly, lz}];
az = 1;
psiPW3 = Exp[I ((2. Pi/lx) ax gA3["X"] + (2. Pi/ly) ay gA3["Y"] + (2. Pi/lz) az gA3["Z"])]/Sqrt[nx ny nz];
cPW3 = ft[psiPW3];
iz = First@FirstPosition[Round[gA3["pz"] lz/(2. Pi)], az];
chk["A5 3D plane wave bin weight - 1 (tol 1e-12)", Abs[cPW3[[ix, iy, iz]]]^2 - 1, 1.*^-12]
chk["A6 3D round trip (tol 1e-12)", Max@Abs[ift[cPW3] - psiPW3], 1.*^-12]
(* A7: leakage guards fire on the right axis, so a swapped marginal orientation cannot
   pass silently *)
gL = grid2[{64, 64}, {16., 16.}];
chk["A7 contained packet: worst edge/tail leakage (tol 1e-6)",
  Max[{edgeMass2[gL, #], tailMass2[gL, #]} &[gaussian2[gL, {0., 0.}, {0., 0.}, {1., 1.}]]], 1.*^-6]
edgeWall = edgeMass2[gL, gaussian2[gL, {7.2, 0.}, {0., 0.}, {0.7, 0.7}]];
chk["A7 packet parked at the x wall: loud on x, quiet on y (PASS = x > 0.1 and y < 1e-4)",
  Boole[! (edgeWall[[1]] > 0.1 && edgeWall[[2]] < 1.*^-4)], 0]
tailFast = tailMass2[gL, gaussian2[gL, {0., 0.}, {0., 11.5}, {0.7, 0.7}]];
chk["A7 packet near the py ladder top: loud on py, quiet on px (PASS = py > 0.1 and px < 1e-4)",
  Boole[! (tailFast[[2]] > 0.1 && tailFast[[1]] < 1.*^-4)], 0]

w["==== B: engine vs closed forms ===="]
(* B1: free spreading, full 4x4 covariance *)
gB = grid2[{128, 96}, {40., 30.}];
{x0, y0, px0, py0, sx, sy} = {-2., 1., 0.7, -0.4, 1.2, 0.8};
psiB0 = gaussian2[gB, {x0, y0}, {px0, py0}, {sx, sy}];
dtB = 0.005; ntB = 200; tfB = ntB dtB;
psiBT = Last[trajectory2[gB, 0 &, {}, dtB, ntB][psiB0]];
mB = moments2[gB, psiBT];
(* the exact Gaussian flow is MatrixExp of the filter's own drift matrix, so this line
   tests driftMatrix2 and the grid dynamics against each other *)
expAB = MatrixExp[driftMatrix2[0. IdentityMatrix[2]] tfB];
sig0B = DiagonalMatrix[{sx^2, sy^2, hbar^2/(4 sx^2), hbar^2/(4 sy^2)}];
chk["B1 free-spreading mean vs MatrixExp of the drift (tol 1e-10)",
  Max@Abs[mB["mean"] - expAB . {x0, y0, px0, py0}], 1.*^-10]
chk["B1 free-spreading 4x4 covariance vs the MatrixExp flow (tol 1e-10)",
  Max@Abs[mB["cov"] - expAB . sig0B . Transpose[expAB]], 1.*^-10]
chk["B1 norm error (tol 1e-12)", normError2[psiBT], 1.*^-12]

(* B2: harmonic breathing at lambda = 0, Strang error O(dt^2) *)
omB = 1.3; gB2 = grid2[{128, 128}, {20., 20.}];
psiB20 = gaussian2[gB2, {0., 0.}, {0., 0.}, {0.9, 0.5}];
dtB2 = 0.002; ntB2 = 600; tB2 = ntB2 dtB2;
psiB2T = Last[trajectory2[gB2, 0.5 mass omB^2 (#1^2 + #2^2) &, {}, dtB2, ntB2][psiB20]];
mB2 = moments2[gB2, psiB2T];
expAH = MatrixExp[driftMatrix2[mass omB^2 IdentityMatrix[2]] tB2];
sig0H = DiagonalMatrix[{0.9^2, 0.5^2, hbar^2/(4 0.9^2), hbar^2/(4 0.5^2)}];
chk["B2 breathing covariance vs the MatrixExp flow (tol 1e-4, Strang O(dt^2))",
  Max@Abs[mB2["cov"] - expAH . sig0H . Transpose[expAH]], 1.*^-4]
chk["B2 means stay at rest (tol 1e-10)", Max@Abs[mB2["mean"]], 1.*^-10]

(* B3: one measurement update vs the exact Gaussian posterior *)
gM = grid2[{128, 128}, {20., 20.}];
{mx0, my0, sxM, syM} = {0.3, -0.2, 0.8, 1.1};
psiM = gaussian2[gM, {mx0, my0}, {0., 0.}, {sxM, syM}];
dtM = 0.01;
chA = axisChannels2[gM, {2., 0.5}];
outsA = {1.0, -0.7};
mA = moments2[gM, measUpdate2[gM, chA, dtM, outsA][psiM]];
postVar[v_, lam_, dtq_] := 1/(1/v + 4 lam dtq);
postMean[m_, v_, lam_, xb_, dtq_] := postVar[v, lam, dtq] (m/v + 4 lam dtq xb);
chk["B3a two-channel posterior variance error (tol 1e-12)",
  Max@Abs[{mA["cov"][[1, 1]] - postVar[sxM^2, 2., dtM], mA["cov"][[2, 2]] - postVar[syM^2, 0.5, dtM]}], 1.*^-12]
chk["B3a two-channel posterior mean error (tol 1e-12)",
  Max@Abs[{mA["mean"][[1]] - postMean[mx0, sxM^2, 2., 1.0, dtM], mA["mean"][[2]] - postMean[my0, syM^2, 0.5, -0.7, dtM]}], 1.*^-12]
th = 0.6; uR = {Cos[th], Sin[th]}; lamR = 1.7; xbR = 0.9;
mR = moments2[gM, measUpdate2[gM, setupChannels2[gM, {channel[uR, lamR]}], dtM, {xbR}][psiM]];
sig0 = {{sxM^2, 0}, {0, syM^2}}; m0 = {mx0, my0};
kden = 1 + 4 lamR dtM uR . sig0 . uR;
chk["B3b rotated-channel posterior covariance error (tol 1e-12)",
  Max@Abs[mR["cov"][[1 ;; 2, 1 ;; 2]] - (sig0 - 4 lamR dtM Outer[Times, sig0 . uR, sig0 . uR]/kden)], 1.*^-12]
chk["B3b rotated-channel posterior mean error (tol 1e-12)",
  Max@Abs[mR["mean"][[1 ;; 2]] - (m0 + (4 lamR dtM (xbR - uR . m0)/kden) sig0 . uR)], 1.*^-12]
postD = Abs[measUpdate2[gM, chA, dtM, outsA][psiM]]^2;
bayesRaw = Exp[-2 dtM Total[MapThread[#1["lambda"] (#1["c"] - #2)^2 &, {chA, outsA}]]] Abs[psiM]^2;
chk["B3c Bayes verbatim: posterior density vs likelihood times prior (tol 1e-14)",
  Max@Abs[postD - bayesRaw/Total[bayesRaw, 2]], 1.*^-14]

(* B3d: simultaneous NON-orthogonal channels, the case orthogonal-axis tests cannot see:
   the posterior precision adds 4 dt lambda_i u_i u_i^T per channel *)
uA = {1., 0.}; uB = {Cos[0.6], Sin[0.6]}; {lA, lB} = {2., 1.3};
chNO = setupChannels2[gM, {channel[uA, lA], channel[uB, lB]}];
outsNO = {0.7, -0.4};
mNO = moments2[gM, measUpdate2[gM, chNO, dtM, outsNO][psiM]];
sig0pos = DiagonalMatrix[{sxM^2, syM^2}];
sigPostNO = Inverse[Inverse[sig0pos] + 4 dtM (lA Outer[Times, uA, uA] + lB Outer[Times, uB, uB])];
mPostNO = sigPostNO . (Inverse[sig0pos] . {mx0, my0} + 4 dtM (lA outsNO[[1]] uA + lB outsNO[[2]] uB));
chk["B3d non-orthogonal pair posterior covariance (tol 1e-12)",
  Max@Abs[mNO["cov"][[1 ;; 2, 1 ;; 2]] - sigPostNO], 1.*^-12]
chk["B3d non-orthogonal pair posterior mean (tol 1e-12)",
  Max@Abs[mNO["mean"][[1 ;; 2]] - mPostNO], 1.*^-12]

(* B3e: two channels along the SAME direction are one channel at the summed strength *)
chk["B3e same-direction strengths add in the Kraus exponent (tol 1e-14)",
  Max@Abs[measUpdate2[gM, setupChannels2[gM, {channel[uB, 0.9], channel[uB, 1.4]}], dtM, {0.5, 0.5}][psiM] -
     measUpdate2[gM, setupChannels2[gM, {channel[uB, 2.3]}], dtM, {0.5}][psiM]], 1.*^-14]

(* B3f: the strength guard returns a Failure for a symbolic strength *)
ClearAll[symLam];
chk["B3f setup rejects a symbolic strength with a Failure (0 = True)",
  Boole[! FailureQ[setupChannels2[gM, {channel[{1, 0}, symLam]}]]], 0]

(* B4: sampler statistics on a correlated state (batched draws, statistical tolerances) *)
phi = 0.5; c1 = (gM["X"] - 0.5) Cos[phi] + (gM["Y"] + 0.3) Sin[phi];
c2 = -(gM["X"] - 0.5) Sin[phi] + (gM["Y"] + 0.3) Cos[phi];
psiC = normalizeState[Exp[-c1^2/(4 1.5^2) - c2^2/(4 0.6^2)]];
mC = moments2[gM, psiC];
SeedRandom[42];
samples = RandomChoice[Flatten[Abs[psiC]^2] -> gM["rflat"], 60000];
mS = Mean[samples];
ds = Map[# - mS &, samples];
chk["B4 sample mean vs state mean (tol 0.03, 60000 draws)",
  Max@Abs[mS - mC["mean"][[1 ;; 2]]], 0.03]
chk["B4 sample covariance vs state covariance, correlated state (tol 0.05)",
  Max@Abs[Transpose[ds] . ds/Length[samples] - mC["cov"][[1 ;; 2, 1 ;; 2]]], 0.05]

(* B4b: simultaneous outcomes are correlated through the SHARED position sample:
   Cov(out_1, out_2) = u_1 . Sigma_pos . u_2 at ANY step size, since the noises are
   independent. The draws below batch the engine's own sampling law (the same weighted
   choice off rflat that jointSample2 makes) with the per-channel noises added exactly
   as drawOutcomes2 adds them, at a step long enough that the signal is not buried
   under the outcome noise (the estimator's statistical floor sits well inside the
   tolerance). Per-marginal sampling would give zero cross-covariance here. *)
SeedRandom[8];
dtB4 = 0.25;
rs = RandomChoice[Flatten[Abs[psiC]^2] -> gM["rflat"], 40000];
o1 = rs . uA + RandomVariate[NormalDistribution[0, 1/Sqrt[4 lA dtB4]], 40000];
o2 = rs . uB + RandomVariate[NormalDistribution[0, 1/Sqrt[4 lB dtB4]], 40000];
chk["B4b outcome cross-covariance minus u1.Sigma.u2 (tol 0.08, 40000 draws, ~4.5 sigma)",
  (Mean[o1 o2] - Mean[o1] Mean[o2]) - uA . mC["cov"][[1 ;; 2, 1 ;; 2]] . uB, 0.08]

(* B5: 3D core: free spreading against the 6x6 MatrixExp flow, posterior, conditional steps *)
g3 = grid3[{64, 48, 32}, {24., 20., 16.}];
s3 = {1.1, 0.9, 0.7};
psi30 = gaussian3[g3, {0.5, -0.4, 0.2}, {0.3, -0.2, 0.5}, s3];
dt3 = 0.01; nt3 = 50; t3 = nt3 dt3;
m3T = moments3[g3, Nest[evolveStep3[g3, 0 &, dt3], psi30, nt3]];
expA3 = MatrixExp[driftMatrix2[0. IdentityMatrix[3]] t3];
sig03 = DiagonalMatrix[Join[s3^2, hbar^2/(4 s3^2)]];
chk["B5 3D free-spreading mean vs MatrixExp of the 6x6 drift (tol 1e-10)",
  Max@Abs[m3T["mean"] - expA3 . {0.5, -0.4, 0.2, 0.3, -0.2, 0.5}], 1.*^-10]
chk["B5 3D free-spreading 6x6 covariance vs the MatrixExp flow (tol 1e-10)",
  Max@Abs[m3T["cov"] - expA3 . sig03 . Transpose[expA3]], 1.*^-10]
ch3 = axisChannels3[g3, {2., 0.5, 1.}];
m3u = moments3[g3, measUpdate3[g3, ch3, dtM, {0.8, -0.6, 0.4}][psi30]];
chk["B5 3D three-channel posterior variance error (tol 1e-12)",
  Max@Abs[Diagonal[m3u["cov"]][[1 ;; 3]] -
     MapThread[postVar[#1^2, #2, dtM] &, {s3, {2., 0.5, 1.}}]], 1.*^-12]
SeedRandom[7];
chk["B5 3D one conditional step norm error (tol 1e-12)",
  normError2[step3[g3, 0 &, ch3, dt3, drawOutcomes3[g3, ch3, dt3][psi30]][psi30]], 1.*^-12]
SeedRandom[11];
chk["B5 3D sampler mean vs packet center (tol 0.03, 20000 draws)",
  Max@Abs[Mean[RandomChoice[Flatten[Abs[psi30]^2] -> g3["rflat"], 20000]] - {0.5, -0.4, 0.2}], 0.03]
(* a real 3D conditional run: one channel with all three components nonzero, on a state
   correlated across the axes, in an isotropic trap; the conditional state must stay pure
   in the 6x6 symplectic sense and localize the measured combination *)
g3r = grid3[{64, 64, 64}, {16., 16., 16.}];
u3 = Normalize[{1., 1., 1.}];
ch3r = setupChannels3[g3r, {channel[u3, 1.5]}];
psi3r0 = normalizeState[Exp[-(g3r["X"]^2 + g3r["Y"]^2 + g3r["Z"]^2)/(4 0.9^2) + 0.3 g3r["X"] g3r["Y"]/(2 0.9^2)]];
psi3rT = BlockRandom[SeedRandom[13];
   Nest[Function[p, step3[g3r, 0.5 mass (#1^2 + #2^2 + #3^2) &, ch3r, 0.002,
       drawOutcomes3[g3r, ch3r, 0.002][p]][p]], psi3r0, 800]];
m3r = moments3[g3r, psi3rT];
omega6 = ArrayFlatten[{{0, IdentityMatrix[3]}, {-IdentityMatrix[3], 0}}];
chk["B5 3D monitored run stays pure: max |symplectic eigenvalue - hbar/2| (tol 1e-8)",
  Max@Abs[Abs[Eigenvalues[I omega6 . m3r["cov"]]] - hbar/2], 1.*^-8]
chk["B5 3D measured combination localizes below its initial variance (0 = True)",
  Boole[! (posEmbed2[u3] . m3r["cov"] . posEmbed2[u3] <
      posEmbed2[u3] . moments3[g3r, psi3r0]["cov"] . posEmbed2[u3])], 0]
chk["B5 3D run norm error (tol 1e-12)", normError2[psi3rT], 1.*^-12]
chk["B5 3D run stays contained: worst per-axis edge/tail guard (tol 1e-8)",
  Max[{edgeMass3[g3r, psi3rT], tailMass3[g3r, psi3rT]}], 1.*^-8]

w["==== C: separable 2D runs against the essay's verified 1D limits ===="]
(* C1: free particle, both axes monitored at different strengths *)
gF = grid2[{128, 128}, {50., 50.}];
chF = axisChannels2[gF, {1., 2.}];
psiF = BlockRandom[SeedRandom[11];
   Last[trajectory2[gF, 0 &, chF, 0.005, 1200][gaussian2[gF, {0., 0.}, {0., 0.}, {2., 2.}]]]];
mF = moments2[gF, psiF];
chkRatio["C1 Sxx / free steady, lambda_x = 1 (tol 0.02, O(dt) bias)", mF["cov"][[1, 1]]/freeSSvar[1.], 0.02]
chkRatio["C1 Syy / free steady, lambda_y = 2 (tol 0.02)", mF["cov"][[2, 2]]/freeSSvar[2.], 0.02]
chk["C1 cross-block covariance: product state stays product (tol 1e-10)",
  Max@Abs[{mF["cov"][[1, 2]], mF["cov"][[1, 4]], mF["cov"][[2, 3]], mF["cov"][[3, 4]]}], 1.*^-10]
chk["C1 per-axis purity SxSp - C^2 - hbar^2/4 (tol 1e-3)",
  Max@Abs[{mF["cov"][[1, 1]] mF["cov"][[3, 3]] - mF["cov"][[1, 3]]^2,
     mF["cov"][[2, 2]] mF["cov"][[4, 4]] - mF["cov"][[2, 4]]^2} - hbar^2/4], 1.*^-3]
chk["C1 final-state leakage after 1200 measured steps (tol 1e-4)",
  Max[{edgeMass2[gF, psiF], tailMass2[gF, psiF]}], 1.*^-4]

(* C2: harmonic trap, different frequencies and strengths per axis *)
gH = grid2[{128, 128}, {16., 16.}];
{omx, omy, lhx, lhy} = {1., 1.5, 2., 0.7};
psiH = BlockRandom[SeedRandom[23];
   Last[trajectory2[gH, 0.5 mass (omx^2 #1^2 + omy^2 #2^2) &, axisChannels2[gH, {lhx, lhy}], 0.002, 1500][
     gaussian2[gH, {0., 0.}, {0., 0.}, {Sqrt[hbar/(2 mass omx)], Sqrt[hbar/(2 mass omy)]}]]]];
mHv = moments2[gH, psiH];
chkRatio["C2 Sxx / harmonic Riccati steady (tol 0.01)", mHv["cov"][[1, 1]]/harmSSvar[lhx, omx], 0.01]
chkRatio["C2 Syy / harmonic Riccati steady (tol 0.01)", mHv["cov"][[2, 2]]/harmSSvar[lhy, omy], 0.01]
chk["C2 squeezed below both zero points (0 = both True)",
  Boole[! (mHv["cov"][[1, 1]] < hbar/(2 mass omx) && mHv["cov"][[2, 2]] < hbar/(2 mass omy))], 0]
chk["C2 final-state leakage after 1500 measured steps (tol 1e-8)",
  Max[{edgeMass2[gH, psiH], tailMass2[gH, psiH]}], 1.*^-8]

(* C3a: exact 2D Lindblad on 32^2: additive heating to machine precision *)
g32 = grid2[{32, 32}, {16., 16.}];
{l32x, l32y} = {1., 0.6}; dt32 = 0.002;
fm = FourierMatrix[32, FourierParameters -> {0, -1}];
F2 = KroneckerProduct[fm, fm];
psi32 = gaussian2[g32, {0., 0.}, {0., 0.}, {1., 1.}];
chk["C3a flattened transform matrix matches ft (tol 1e-12)",
  Max@Abs[F2 . Flatten[psi32] - Flatten[ft[psi32]]], 1.*^-12]
p2op = ConjugateTranspose[F2] . (Flatten[g32["p2"]] F2);
uniM = Transpose[Map[Flatten[evolveStep2[g32, 0 &, dt32][Partition[#, 32]]] &, N[IdentityMatrix[1024]]]];
chk["C3a materialized Strang propagator is unitary (tol 1e-12)",
  Max@Abs[uniM . ConjugateTranspose[uniM] - IdentityMatrix[1024]], 1.*^-12]
Xf = Flatten[g32["X"]]; Yf = Flatten[g32["Y"]];
decHalf = Exp[-(dt32/2) ((l32x/2) Outer[Subtract, Xf, Xf]^2 + (l32y/2) Outer[Subtract, Yf, Yf]^2)];
lstep[rho_] := decHalf (uniM . (decHalf rho) . ConjugateTranspose[uniM]);
rhos = NestList[lstep, With[{v = Flatten[psi32]}, KroneckerProduct[v, Conjugate[v]]], 40];
incs = Differences[Map[Re[Total[# Transpose[p2op], 2]] &, rhos]];
chk["C3a per-step <p^2> gain minus (lx + ly) hbar^2 dt, max dev (tol 1e-10)",
  Max@Abs[incs - (l32x + l32y) hbar^2 dt32], 1.*^-10]
chk["C3a trace drift (tol 1e-12)", Max@Abs[Map[Re[Tr[#]] - 1 &, rhos]], 1.*^-12]
(* positivity is the invariant with teeth (the essay's phrase): the sandwich structure
   preserves the trace even for a wrong-signed kernel, the spectrum would not survive *)
chk["C3a rho stays a state: most negative eigenvalue over 5 sampled steps (tol 1e-12)",
  Max[0., -Min[Map[Min[Re[Eigenvalues[#]]] &, rhos[[1 ;; 41 ;; 10]]]]], 1.*^-12]
chk["C3a rho stays Hermitian (tol 1e-12)",
  Max[Map[Max[Abs[# - ConjugateTranspose[#]]] &, rhos[[1 ;; 41 ;; 10]]]], 1.*^-12]
chk["C3a p2 operator Hermitian (tol 1e-12)", Max@Abs[p2op - ConjugateTranspose[p2op]], 1.*^-12]

(* C3b: trajectory-averaged SSE heating, free particle, 60 seeds *)
gE = grid2[{64, 64}, {24., 24.}];
chE = axisChannels2[gE, {1., 0.5}];
psiE0 = gaussian2[gE, {0., 0.}, {0., 0.}, {1., 1.}];
p2of[mom_] := mom["cov"][[3, 3]] + mom["mean"][[3]]^2 + mom["cov"][[4, 4]] + mom["mean"][[4]]^2;
eTraj[seed_] := BlockRandom[SeedRandom[seed];
   Map[p2of[moments2[gE, #]]/(2 mass) &,
    trajectory2[gE, 0 &, chE, 0.002, 400][psiE0][[1 ;; -1 ;; 50]]]];
eMean = Mean[Map[eTraj, Range[60]]];
chkRatio["C3b ensemble dE/dt over (lx + ly) hbar^2/2m (tol 0.1, 60 seeds)",
  ((Last[eMean] - First[eMean])/0.8)/((1. + 0.5) hbar^2/(2 mass)), 0.1]

(* C4: freeze the y axis, replay the 1D essay record: exact product at every step *)
grid1[n_, l_] := With[{dx = l/n},
   <|"n" -> n, "dx" -> dx, "x" -> dx Range[-n/2, n/2 - 1],
     "p" -> (2 Pi hbar/l) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
gaussian1[g_, xc_, pc_, s_] := Normalize[Exp[-(g["x"] - xc)^2/(4 s^2) + I pc g["x"]/hbar]];
drawOut1[g_, lam_, dt_][psi_] := RandomChoice[Abs[psi]^2 -> g["x"]] +
   RandomVariate[NormalDistribution[0, 1/Sqrt[4 lam dt]]];
measUpdate1[g_, lam_, dt_, xb_][psi_] := Normalize[Exp[-lam dt (g["x"] - xb)^2] psi];
potentialStep1[g_, Vfun_, dt_][psi_] := Exp[-I Vfun[g["x"]] dt/(2 hbar)] psi;
unitaryStep1[g_, dt_][psi_] := ift[Exp[-I g["p"]^2 dt/(2 mass hbar)] ft[psi]];
evolveStep1[g_, Vfun_, dt_] := potentialStep1[g, Vfun, dt] @* unitaryStep1[g, dt] @* potentialStep1[g, Vfun, dt];
stepV1[g_, Vfun_, lam_, dt_, xb_][psi_] := evolveStep1[g, Vfun, dt][measUpdate1[g, lam, dt, xb][psi]];
trajectory1[g_, Vfun_, lam_, dt_, nt_][psi0_] :=
   NestList[stepV1[g, Vfun, lam, dt, Sow@drawOut1[g, lam, dt][#]][#] &, psi0, nt];
g1 = grid1[128, 16.];
{states1, {outs1}} = BlockRandom[SeedRandom[21];
   Reap[trajectory1[g1, 0.5 mass #^2 &, 2., 0.002, 500][gaussian1[g1, 0., 0., Sqrt[hbar/(2 mass)]]]]];
g2f = grid2[{128, 32}, {16., 8.}];
psiy0 = Normalize[N[ConstantArray[1., 32]]];
states2 = replay2[g2f, 0.5 mass #1^2 &, setupChannels2[g2f, {channel[{1, 0}, 2.]}], 0.002,
    Map[List, outs1]][Outer[Times, gaussian1[g1, 0., 0., Sqrt[hbar/(2 mass)]], psiy0]];
chk["C4 frozen-axis replay: max over 500 steps of |psi2D - psi1D x flat| (tol 1e-12)",
  Max[MapThread[Max[Abs[#2 - Outer[Times, #1, psiy0]]] &, {states1, states2}]], 1.*^-12]

w["==== D: multimode Gaussian filter ===="]
(* D0: derive the 1D steady state from the Riccati (the essay's route), prove there is
   exactly one physical branch and that it attracts, take its limits, and only then use
   it as the target of the multimode fixed-point identities. *)
ClearAll[lamXs, lamYs, omXs, omYs, ms, hs, thS, omS, lamS, om1S, om2S, phS, SxV, CcV, SpV];
$Assumptions = lamXs > 0 && lamYs > 0 && omXs > 0 && omYs > 0 && ms > 0 && hs > 0 &&
   omS > 0 && lamS > 0 && 0 < thS < Pi/2 && om1S > 0 && om2S > 0 && 0 < phS < Pi/2;
ric1D = {2 CcV/ms - 4 lamS SxV^2, SpV/ms - ms omS^2 SxV - 4 lamS SxV CcV,
   -2 ms omS^2 CcV + lamS hs^2 - 4 lamS CcV^2};
roots1D = SolveValues[Thread[ric1D == 0], {SxV, CcV, SpV}];
physQ[sol_] := TrueQ[Resolve[ForAll[{ms, lamS, omS, hs},
     ms > 0 && lamS > 0 && omS >= 0 && hs > 0,
     Element[sol[[1]] | sol[[3]], Reals] && sol[[1]] > 0 && sol[[3]] > 0], Reals]];
chk["D0 exactly one physical root among the algebraic roots (0 = unique)",
  Boole[Count[Map[physQ, roots1D], True] != 1], 0]
rootPhys = SelectFirst[roots1D, physQ];
chk["D0 derived physical root equals the named closed form (exact 0)",
  Max@Abs[Simplify[rootPhys - modeSSsym[lamS, omS, ms, hs]]], 0]
chk["D0 purity at the steady state: Sx Sp - C^2 - hs^2/4 (exact 0)",
  Simplify[rootPhys[[1]] rootPhys[[3]] - rootPhys[[2]]^2 - hs^2/4], 0]
(* stability, symbolically: Routh-Hurwitz on the linearized flow at the root, for ALL
   positive parameters, the same ForAll shape as the branch selection *)
jacSym = Simplify[D[ric1D, {{SxV, CcV, SpV}}] /. Thread[{SxV, CcV, SpV} -> modeSSsym[lamS, omS, ms, hs]]];
ClearAll[zRH];
cpoly = CoefficientList[CharacteristicPolynomial[jacSym, zRH], zRH];
{a0RH, a1RH, a2RH} = Simplify[Most[cpoly]/Last[cpoly]];
chk["D0 stability: Routh-Hurwitz for the physical root, all positive parameters (0 = proved)",
  Boole[! TrueQ[Resolve[ForAll[{ms, lamS, omS, hs}, ms > 0 && lamS > 0 && omS > 0 && hs > 0,
       a2RH > 0 && a0RH > 0 && a2RH a1RH - a0RH > 0], Reals]]], 0]
chk["D0 limit om -> 0 recovers the free-particle steady variance (exact 0)",
  Simplify[Limit[rootPhys[[1]], omS -> 0, Direction -> "FromAbove",
     Assumptions -> lamS > 0 && ms > 0 && hs > 0] - Sqrt[hs/(4 lamS ms)]], 0]
chk["D0 classical limit hs -> 0: the conditional variance collapses (exact 0)",
  Simplify[Limit[rootPhys[[1]], hs -> 0, Direction -> "FromAbove",
    Assumptions -> lamS > 0 && ms > 0 && omS > 0]], 0]
chk["D0 strong-measurement limit lam -> Infinity (exact 0)",
  Simplify[Limit[rootPhys[[1]], lamS -> Infinity,
    Assumptions -> omS > 0 && ms > 0 && hs > 0]], 0]
(* detection efficiency: the essay's rescaling lam -> eta lam, hbar -> hbar/Sqrt[eta] on
   the information terms gives the eta-Riccati fixed point and the mixed-state purity *)
ClearAll[etaS];
etaSSm = modeSSsym[etaS lamS, omS, ms, hs/Sqrt[etaS]];
rhsEta = Block[{hbar = hs, mass = ms}, riccatiRHS2[{{ms omS^2}},
     {<|"u" -> {1}, "lambda" -> lamS, "eta" -> etaS|>}][
    {{etaSSm[[1]], etaSSm[[2]]}, {etaSSm[[2]], etaSSm[[3]]}}]];
chk["D0 eta-filter fixed point at the essay's rescaled steady state (exact 0)",
  Max@Abs[Simplify[rhsEta, {lamS > 0, omS > 0, ms > 0, hs > 0, 0 < etaS <= 1}]], 0]
chk["D0 eta purity: Sx Sp - C^2 - hs^2/(4 eta) at the rescaled steady state (exact 0)",
  Simplify[etaSSm[[1]] etaSSm[[3]] - etaSSm[[2]]^2 - hs^2/(4 etaS),
   {lamS > 0, omS > 0, ms > 0, hs > 0, 0 < etaS <= 1}], 0]
{sxA, ccA, spA} = modeSSsym[lamXs, omXs, ms, hs];
{syA, ccB, spB} = modeSSsym[lamYs, omYs, ms, hs];
rhsAxes = Block[{hbar = hs, mass = ms},
   riccatiRHS2[ms DiagonalMatrix[{omXs^2, omYs^2}],
     {<|"u" -> {1, 0}, "lambda" -> lamXs|>, <|"u" -> {0, 1}, "lambda" -> lamYs|>}][
    {{sxA, 0, ccA, 0}, {0, syA, 0, ccB}, {ccA, 0, spA, 0}, {0, ccB, 0, spB}}]];
chk["D0a axis-aligned symbolic fixed point (exact 0)", Max@Abs[Simplify[rhsAxes]], 0]
{sxR, ccR, spR} = modeSSsym[lamS, omS, ms, hs];
r2m = {{Cos[thS], -Sin[thS]}, {Sin[thS], Cos[thS]}};
r4m = ArrayFlatten[{{r2m, 0}, {0, r2m}}];
sigRot = r4m . {{sxR, 0, ccR, 0}, {0, hs/(2 ms omS), 0, 0}, {ccR, 0, spR, 0},
     {0, 0, 0, hs ms omS/2}} . Transpose[r4m];
rhsRot = Block[{hbar = hs, mass = ms},
   riccatiRHS2[ms omS^2 IdentityMatrix[2],
     {<|"u" -> {Cos[thS], Sin[thS]}, "lambda" -> lamS|>}][sigRot]];
chk["D0b rotated-channel symbolic fixed point, theta symbolic (exact 0)",
  Max@Abs[Simplify[rhsRot]], 0]
ClearAll[k11, k12, k22, c1, s1, c2, s2, l1, l2,
  mm11, mm12, mm13, mm14, mm22, mm23, mm24, mm33, mm34, mm44];
mSym = {{mm11, mm12, mm13, mm14}, {mm12, mm22, mm23, mm24},
   {mm13, mm23, mm33, mm34}, {mm14, mm24, mm34, mm44}};
kSym = {{k11, k12}, {k12, k22}};
aSym = Block[{mass = ms}, driftMatrix2[kSym]];
hQuad = ArrayFlatten[{{kSym/2, 0}, {0, IdentityMatrix[2]/(2 ms)}}];
chk["D0c quadratic drift conserves energy for any moments (exact 0)",
  Simplify[Tr[(aSym . mSym + mSym . Transpose[aSym]) . hQuad]], 0]
dSym = hs^2 l1 Outer[Times, {0, 0, c1, s1}, {0, 0, c1, s1}] +
   hs^2 l2 Outer[Times, {0, 0, c2, s2}, {0, 0, c2, s2}];
chk["D0c diffusion feeds energy at hbar^2 (l1 + l2)/2m (exact 0)",
  Simplify[Tr[dSym . hQuad] - hs^2 (l1 + l2)/(2 ms), {c1^2 + s1^2 == 1, c2^2 + s2^2 == 1}], 0]

(* theta -> 0 collapses the rotated channel onto the x channel at the level of the flow
   itself, for ARBITRARY symmetric covariance, not only at the fixed point *)
rhsThetaArb = Block[{hbar = hs, mass = ms}, riccatiRHS2[ms omS^2 IdentityMatrix[2],
     {<|"u" -> {Cos[thS], Sin[thS]}, "lambda" -> lamS|>}][mSym]];
rhsAxisArb = Block[{hbar = hs, mass = ms}, riccatiRHS2[ms omS^2 IdentityMatrix[2],
     {<|"u" -> {1, 0}, "lambda" -> lamS|>}][mSym]];
chk["D0 theta -> 0 collapse of the rotated flow, arbitrary Sigma (exact 0)",
  Max@Abs[Simplify[(rhsThetaArb /. thS -> 0) - rhsAxisArb]], 0]
(* the D2 conduit: a rotated trap couples the axes only while the mode frequencies differ *)
kRotSym = ms {{Cos[phS], -Sin[phS]}, {Sin[phS], Cos[phS]}} . DiagonalMatrix[{om1S^2, om2S^2}] .
    Transpose[{{Cos[phS], -Sin[phS]}, {Sin[phS], Cos[phS]}}];
chk["D0 degenerate-trap limit om2 -> om1 closes the cross-stiffness conduit (exact 0)",
  Simplify[kRotSym[[1, 2]] /. om2S -> om1S], 0]

(* the filter's mode-count generality, certified at 3 modes with everything symbolic:
   three axis channels at three strengths in a diagonal 3D trap, the per-axis 1D steady
   states embedded in the 6x6, and the SAME riccatiRHS2 code *)
ClearAll[l1s, l2s, l3s, o1s, o2s, o3s];
$Assumptions = $Assumptions && l1s > 0 && l2s > 0 && l3s > 0 && o1s > 0 && o2s > 0 && o3s > 0;
mode3 = Map[modeSSsym[#[[1]], #[[2]], ms, hs] &, {{l1s, o1s}, {l2s, o2s}, {l3s, o3s}}];
sig6 = ArrayFlatten[{{DiagonalMatrix[mode3[[All, 1]]], DiagonalMatrix[mode3[[All, 2]]]},
    {DiagonalMatrix[mode3[[All, 2]]], DiagonalMatrix[mode3[[All, 3]]]}}];
rhs6 = Block[{hbar = hs, mass = ms}, riccatiRHS2[ms DiagonalMatrix[{o1s^2, o2s^2, o3s^2}],
     {<|"u" -> {1, 0, 0}, "lambda" -> l1s|>, <|"u" -> {0, 1, 0}, "lambda" -> l2s|>,
      <|"u" -> {0, 0, 1}, "lambda" -> l3s|>}][sig6]];
chk["D0 three-mode symbolic fixed point of the same filter code (exact 0)",
  Max@Abs[Simplify[rhs6]], 0]

(* the filter's 1D restriction IS the essay's five-number filter: run both on one
   synthetic outcome stream and compare every component every step *)
ess1D[lam_, om_, dtq_][{mx_, mp_, sx_, cc_, sp_}, xb_] :=
   With[{dWf = 2 Sqrt[lam] dtq (xb - mx)},
    Join[{mx + mp/mass dtq + 2 Sqrt[lam] sx dWf, mp - mass om^2 mx dtq + 2 Sqrt[lam] cc dWf},
     {sx, cc, sp} + dtq {2 cc/mass - 4 lam sx^2, sp/mass - mass om^2 sx - 4 lam sx cc,
        -2 mass om^2 cc + lam hbar^2 - 4 lam cc^2}]];
SeedRandom[19];
outs1f = RandomVariate[NormalDistribution[0, 3], 400];
essRun = FoldList[ess1D[2., 1., 0.002], {0., 0., 0.5, 0., 0.5}, outs1f];
oneRun = filterRun2[{{mass 1.^2}}, {<|"u" -> {1}, "lambda" -> 2.|>}, 0.002, Map[List, outs1f]][
   {{0., 0.}, {{0.5, 0.}, {0., 0.5}}}];
chk["D0 one-mode filterRun2 equals the essay's five-number filter over 400 steps (tol 1e-12)",
  Max@Abs[Map[Flatten[{#[[1]], #[[2, 1, 1]], #[[2, 1, 2]], #[[2, 2, 2]]}] &, oneRun] - essRun], 1.*^-12]

(* D1: isotropic trap watched through one rotated quadrature: grid vs filter *)
omD = 1.; lamD = 2.; thD = Pi/5.; dtD = 0.002; ntD = 1500;
uD = {Cos[thD], Sin[thD]}; uPerp = {-Sin[thD], Cos[thD]};
gD = grid2[{128, 128}, {16., 16.}];
chD = setupChannels2[gD, {channel[uD, lamD]}];
psiD0 = gaussian2[gD, {0., 0.}, {0., 0.}, {Sqrt[hbar/(2 mass omD)], Sqrt[hbar/(2 mass omD)]}];
{psiDf, {trkD}} = BlockRandom[SeedRandom[31];
   Reap[momentsRun2[gD, 0.5 mass omD^2 (#1^2 + #2^2) &, chD, dtD, ntD][psiD0]]];
outsD = trkD[[All, 1]];
momD = trkD[[All, 2]];
posCovs = Map[#["cov"][[1 ;; 2, 1 ;; 2]] &, momD];
sigThGrid = Map[uD . # . uD &, posCovs];
sigPerpGrid = Map[uPerp . # . uPerp &, posCovs];
filtD = filterRun2[mass omD^2 IdentityMatrix[2], chD, dtD, outsD][
   {{0., 0., 0., 0.}, 0.5 IdentityMatrix[4]}];
sigThFilt = Map[uD . #[[2, 1 ;; 2, 1 ;; 2]] . uD &, Rest[filtD]];
chk["D1 grid vs filter, measured-quadrature variance, max over transient (tol 0.005)",
  Max@Abs[sigThGrid - sigThFilt], 0.005]
chkRatio["D1 final measured variance / 1D Riccati steady (tol 0.01)",
  Last[sigThGrid]/sxSSsym[lamD, omD, mass, hbar], 0.01]
chk["D1 dark mode: grid Sigma_perp minus ground, max over run (tol 1e-6, the Strang ripple is (hbar/2m om)(om dt)^2/4 = 5e-7)",
  Max@Abs[sigPerpGrid - hbar/(2 mass omD)], 1.*^-6]
meanGridD = Map[#["mean"][[1 ;; 2]] &, momD];
meanFiltD = Map[#[[1, 1 ;; 2]] &, Rest[filtD]];
chk["D1 mean tracking, max gap over position range (tol 0.005)",
  Max@Abs[meanGridD - meanFiltD]/(Max[meanGridD] - Min[meanGridD]), 0.005]
chk["D1 dark-mode mean pinned at zero (tol 1e-10)", Max@Abs[Map[uPerp . # &, meanGridD]], 1.*^-10]
(* the one invariant that exists in 2D but not 1D: the multimode uncertainty relation.
   At unit efficiency the conditional state is PURE, so every symplectic eigenvalue must
   EQUAL hbar/2, not merely exceed it; the two-sided check would catch a state that
   quietly went mixed *)
omegaSym = {{0, 0, 1, 0}, {0, 0, 0, 1}, {-1, 0, 0, 0}, {0, -1, 0, 0}};
chk["D1 pure-state symplectic spectrum on the grid: max |nu - hbar/2| over the run (tol 1e-9)",
  Max[Map[Max[Abs[Abs[Eigenvalues[I omegaSym . #["cov"]]] - hbar/2]] &, momD]], 1.*^-9]
nuMinFilt = Min[Map[Min[Abs[Eigenvalues[I omegaSym . #[[2]]]]] &, Rest[filtD]]];
chk["D1 multimode uncertainty in the filter (tol 0.005, Euler bias)",
  Max[hbar/2 - nuMinFilt, 0.], 0.005]
chk["D1 filter covariance stays a covariance (0 = final Sigma positive semidefinite)",
  Boole[! PositiveSemidefiniteMatrixQ[Last[filtD][[2]]]], 0]
chk["D1 final-state leakage after 1500 measured steps (tol 1e-8)",
  Max[{edgeMass2[gD, psiDf], tailMass2[gD, psiDf]}], 1.*^-8]
(* the dark mode as a theorem: the pair (drift, measured direction) has observability
   rank 2 of 4 in the isotropic trap, so the orthogonal mode is invisible to the record
   no matter how the filter is run. The rank is taken at SYMBOLIC frequency, mass, and
   channel angle, and the row count follows the phase-space dimension. *)
obsRank[K_, u_] := MatrixRank[NestList[# . driftMatrix2[K] &, posEmbed2[u], 2 Length[K] - 1]];
chk["D1 observability rank of the isotropic trap through one tilted channel is 2, symbolically (0 = True)",
  Boole[Block[{mass = ms}, obsRank[ms omS^2 IdentityMatrix[2], {Cos[thS], Sin[thS]}]] != 2], 0]
(* the dark-mode ripple pinned to its closed form: the Strang map's invariant ellipse
   sits (hbar/2m om)(om dt)^2/8 above the ground variance, and a run started on the
   ground circle orbits it, so the peak excursion is twice that, (hbar/2m om)(om dt)^2/4.
   Measured with the measurement OFF (pure unitary ground-state run), so the number is
   deterministic; it must match the closed form and vanish quadratically with the step *)
rippleD1[dtq_] := Max[Map[Abs[moments2[gD, #]["cov"][[1, 1]] - hbar/(2 mass omD)] &,
    NestList[evolveStep2[gD, 0.5 mass omD^2 (#1^2 + #2^2) &, dtq], psiD0, Round[1.6/dtq]]]];
chkRatio["D1 measurement-off ripple over (hbar/2m om)(om dt)^2/4 (tol 1e-4)",
  rippleD1[dtD]/(hbar/(2 mass omD) (omD dtD)^2/4), 1.*^-4]
chk["D1 ripple vanishes quadratically: doubled-step ratio minus 4 (tol 0.05)",
  rippleD1[2 dtD]/rippleD1[dtD] - 4, 0.05]
note["D1 measured-mode mean wanders (single realization, heating lives here)",
  Max@Abs[Map[uD . # &, meanGridD]]]
note["D1 single-run energy gain beside the ensemble ramp lam hbar^2 t/2m",
  With[{eOf = Function[mom,
       (mom["cov"][[3, 3]] + mom["mean"][[3]]^2 + mom["cov"][[4, 4]] + mom["mean"][[4]]^2)/(2 mass)
        + 0.5 mass omD^2 (mom["cov"][[1, 1]] + mom["mean"][[1]]^2 + mom["cov"][[2, 2]] + mom["mean"][[2]]^2)]},
   {eOf[Last[momD]] - eOf[First[momD]], lamD hbar^2 (ntD - 1) dtD/(2 mass)}]]

(* D2: rotated anisotropic trap, x-only record localizes y through cross-covariance;
   same step dtD and length ntD as D1 *)
{om1, om2, ph2} = {1., 2., Pi/6.};
rPh = {{Cos[ph2], -Sin[ph2]}, {Sin[ph2], Cos[ph2]}};
kD2 = mass rPh . DiagonalMatrix[{om1^2, om2^2}] . Transpose[rPh];
gD2 = grid2[{128, 128}, {16., 16.}];
chD2 = setupChannels2[gD2, {channel[{1, 0}, 2.]}];
sBroad = 1.2;
{psiD2f, {trkD2}} = BlockRandom[SeedRandom[57];
   Reap[momentsRun2[gD2, 0.5 (kD2[[1, 1]] #1^2 + 2 kD2[[1, 2]] #1 #2 + kD2[[2, 2]] #2^2) &,
      chD2, dtD, ntD][gaussian2[gD2, {0., 0.}, {0., 0.}, {sBroad, sBroad}]]]];
momD2 = trkD2[[All, 2]];
filtD2 = filterRun2[kD2, chD2, dtD, trkD2[[All, 1]]][
   {{0., 0., 0., 0.}, DiagonalMatrix[{sBroad^2, sBroad^2, hbar^2/(4 sBroad^2), hbar^2/(4 sBroad^2)}]}];
syyGrid = Map[#["cov"][[2, 2]] &, momD2];
syyFilt = Map[#[[2, 2, 2]] &, Rest[filtD2]];
sxyGrid = Map[#["cov"][[1, 2]] &, momD2];
sxyFilt = Map[#[[2, 1, 2]] &, Rest[filtD2]];
chk["D2 grid vs filter, Sigma_yy gap over the transient swing (tol 0.05, relative)",
  Max@Abs[syyGrid - syyFilt]/(First[syyGrid] - Last[syyGrid]), 0.05]
chk["D2 grid vs filter, Sigma_xy gap over its range (tol 0.05, relative)",
  Max@Abs[sxyGrid - sxyFilt]/(Max[sxyGrid] - Min[sxyGrid]), 0.05]
chk["D2 y localized without ever being measured: final/initial Sigma_yy (tol: below 0.5)",
  Boole[! (Last[syyGrid]/First[syyGrid] < 0.5)], 0]
note["D2 Sigma_yy initial -> final", {First[syyGrid], Last[syyGrid]}]
note["D2 Sigma_xy range (the information conduit)", {Min[sxyGrid], Max[sxyGrid]}]
myGrid = Map[#["mean"][[2]] &, momD2];
myFilt = Map[#[[1, 2]] &, Rest[filtD2]];
chk["D2 <y> tracked from the x record alone, max gap over range (tol 0.06)",
  Max@Abs[myGrid - myFilt]/(Max[myGrid] - Min[myGrid]), 0.06]
chk["D2 final-state leakage after 1500 measured steps (tol 1e-8)",
  Max[{edgeMass2[gD2, psiD2f], tailMass2[gD2, psiD2f]}], 1.*^-8]
chk["D2 pure-state symplectic spectrum on the grid (tol 1e-8)",
  Max[Map[Max[Abs[Abs[Eigenvalues[I omegaSym . #["cov"]]] - hbar/2]] &, momD2]], 1.*^-8]
(* the filter's Euler step, started from a broad pure state, transiently dips below the
   symplectic floor: an integrator artifact of the coarse step, bounded here and absent
   from the grid, whose covariance the check above pins to the pure floor throughout *)
nuMinFiltD2 = Min[Map[Min[Abs[Eigenvalues[I omegaSym . #[[2]]]]] &, Rest[filtD2]]];
chk["D2 filter symplectic dip below hbar/2 stays bounded (tol 0.08, Euler artifact)",
  Max[hbar/2 - nuMinFiltD2, 0.], 0.08]
note["D2 filter min symplectic eigenvalue over the transient", nuMinFiltD2]
(* one channel sees everything here: the rotated anisotropic pair (drift, x direction)
   is fully observable at symbolic frequencies and tilt, the theorem-level converse of
   D1's dark mode; and the physical degeneracy statement joining the two demos is that
   this rank COLLAPSES back to 2 when the mode frequencies degenerate or the tilt
   vanishes, which is when the record's reach shrinks to one mode *)
chk["D2 observability rank of the rotated anisotropic trap through x is 4, symbolically (0 = True)",
  Boole[Block[{mass = ms}, obsRank[kRotSym, {1, 0}]] != 4], 0]
chk["D2 rank collapses to 2 at degenerate frequencies om2 -> om1 (0 = True)",
  Boole[Block[{mass = ms}, obsRank[kRotSym /. om2S -> om1S, {1, 0}]] != 2], 0]
chk["D2 rank collapses to 2 at zero tilt ph -> 0 (0 = True)",
  Boole[Block[{mass = ms}, obsRank[kRotSym /. phS -> 0, {1, 0}]] != 2], 0]
(* the algebraic steady state anchors the demo: RiccatiSolve's stabilizing solution
   solves the engine's own Riccati, is pure, and the engine's Riccati flow reaches it;
   the t = 3 demo above is deliberately MID-TRANSIENT (its final moments differ from
   the steady values), which is why the grid-vs-filter comparison runs through the
   transient rather than at the fixed point *)
sSteadyD2 = RiccatiSolve[{Transpose[driftMatrix2[kD2]], Transpose[{posEmbed2[{1., 0.}]}]},
   {hbar^2 2. Outer[Times, momEmbed2[{1., 0.}], momEmbed2[{1., 0.}]], {{1/(4 2.)}}}];
chk["D2 RiccatiSolve steady state solves the engine's Riccati (tol 1e-12)",
  Max@Abs[riccatiRHS2[kD2, chD2][sSteadyD2]], 1.*^-12]
chk["D2 algebraic steady state is pure (tol 1e-9)",
  Max@Abs[Abs[Eigenvalues[I omegaSym . sSteadyD2]] - hbar/2], 1.*^-9]
sig0D2m = DiagonalMatrix[{sBroad^2, sBroad^2, hbar^2/(4 sBroad^2), hbar^2/(4 sBroad^2)}];
chk["D2 the engine's Riccati flow integrated to t = 30 lands on the algebraic steady (tol 5e-3)",
  Max@Abs[Nest[# + dtD riccatiRHS2[kD2, chD2][#] &, sig0D2m, 15000] - sSteadyD2], 5.*^-3]
note["D2 mid-transient vs steady {Syy, Sxy}",
  {{Last[syyGrid], sSteadyD2[[2, 2]]}, {Last[sxyGrid], sSteadyD2[[1, 2]]}}]

(* D3: failure regime: the steady width must stay resolved, several points per sigma *)
lamStrong = 20.;
freeSSD3 = freeSSvar[lamStrong];
runD3[nn_] := With[{gg = grid2[{nn, nn}, {16., 16.}]},
   BlockRandom[SeedRandom[5];
    moments2[gg, Last[trajectory2[gg, 0 &, axisChannels2[gg, {lamStrong, lamStrong}], 0.001, 1000][
        gaussian2[gg, {0., 0.}, {0., 0.}, {1., 1.}]]]]["cov"][[1, 1]]]];
chkRatio["D3 Sxx/steady at 2.7 points per steady sigma (tol 0.02)", runD3[128]/freeSSD3, 0.02]
chk["D3 under-resolved grid is visibly wrong, 0.67 points per sigma (PASS = deviation > 0.05)",
  Boole[! (Abs[runD3[32]/freeSSD3 - 1] > 0.05)], 0]

(* D4: where the Gaussian filter is the WRONG object: a double well splits the monitored
   conditional state into lobes, which no 4x4 Gaussian can carry. The channel is TILTED
   (Pi/6), so the bimodal x physics writes itself into the x-y correlation through the
   shared record and the state is genuinely two-dimensional while it lasts. The grid
   stays perfectly healthy while the filter, handed the surrogate quadratic trap with
   the true well-bottom curvature (V'' at x = 3 is 4 a (3 x^2 - b) = 5.76) and fed the
   same record, mistracks. Here the grid is the instrument and the filter is the
   approximation that breaks. Same step dtD and length ntD as D1. *)
gW = grid2[{128, 128}, {20., 20.}];
chW = setupChannels2[gW, {channel[{Cos[Pi/6.], Sin[Pi/6.]}, 0.5]}];
{psiWf, {trkW}} = BlockRandom[SeedRandom[71];
   Reap[momentsRun2[gW, 0.08 (#1^2 - 9.)^2 + 0.5 mass #2^2 &, chW, dtD, ntD][
     gaussian2[gW, {0., 0.}, {0., 0.}, {1.6, Sqrt[0.5]}]]]];
momW = trkW[[All, 2]];
filtW = filterRun2[DiagonalMatrix[{5.76, 1.}], chW, dtD, trkW[[All, 1]]][
   {{0., 0., 0., 0.}, DiagonalMatrix[{1.6^2, 0.5, hbar^2/(4 1.6^2), 0.5}]}];
mxG = Map[#["mean"][[1]] &, momW];
mxF = Map[#[[1, 1]] &, Rest[filtW]];
chk["D4 grid stays healthy in the double well (leakage tol 1e-8)",
  Max[{edgeMass2[gW, psiWf], tailMass2[gW, psiWf]}], 1.*^-8]
chk["D4 the lobes write themselves into Sigma_xy through the tilted record (PASS = max > 0.2)",
  Boole[! (Max@Abs[Map[#["cov"][[1, 2]] &, momW]] > 0.2)], 0]
chk["D4 the Gaussian filter mistracks the mean it was built to report (PASS = gap > 0.2 of range)",
  Boole[! (Max@Abs[mxG - mxF]/(Max[mxG] - Min[mxG]) > 0.2)], 0]
note["D4 bimodal transient: max grid Sigma_xx beside the surrogate-well ground variance",
  {Max[Map[#["cov"][[1, 1]] &, momW]], hbar/(2 mass Sqrt[5.76])}]
note["D4 final Sigma_xx, grid vs filter", {Last[Map[#["cov"][[1, 1]] &, momW]],
   Last[Map[#[[2, 1, 1]] &, Rest[filtW]]]}]

(* D5: a run driven past the small grid's ceilings: the 2D form of the essay's closing
   caveat. The continuum heats along the ramp (lx + ly) hbar^2 t forever; the small ring
   cannot follow, and the betrayal is episodic rather than monotonic: whenever the packet
   reaches a wall or the momentum ladder top, wrap-around aliases fast components into
   slow ones and the seam-crossing mask re-localizes the packet, so the conditional
   <p^2> collapses instead of climbing. The guards go loud along the way and the
   endpoint sits far below the ramp: from there on the grid is lying about the
   continuum. *)
gSat = grid2[{32, 32}, {8., 8.}];
chSat = axisChannels2[gSat, {2., 2.}];
satStates = BlockRandom[SeedRandom[3];
   FoldList[Function[{p, k}, Nest[Function[q, step2[gSat, 0 &, chSat, 0.002,
        drawOutcomes2[gSat, chSat, 0.002][q]][q]], p, 5000]],
    gaussian2[gSat, {0., 0.}, {0., 0.}, {0.8, 0.8}], Range[10]]];
satP2 = Map[With[{m = moments2[gSat, #]},
    m["cov"][[3, 3]] + m["mean"][[3]]^2 + m["cov"][[4, 4]] + m["mean"][[4]]^2] &, Rest[satStates]];
satGuards = Map[Max[{edgeMass2[gSat, #], tailMass2[gSat, #]}] &, Rest[satStates]];
chk["D5 the guards go loud during the run (PASS = worst checkpoint reading > 0.1)",
  Boole[! (Max[satGuards] > 0.1)], 0]
chk["D5 the moment collapses episodically instead of climbing (PASS = biggest checkpoint drop < -20)",
  Boole[! (Min[Differences[satP2]] < -20)], 0]
note["D5 <p^2> at t = 10, 20, ..., 100 beside the ramp 4t (the endpoint is wherever the wandering left it)",
  {satP2, 4. Range[10, 100, 10]}]
note["D5 worst guard reading per checkpoint", satGuards]

w["==== summary ===="]
w["OPEN FAILURES: " <> ToString[fails]]
If[fails > 0, Exit[1]];
