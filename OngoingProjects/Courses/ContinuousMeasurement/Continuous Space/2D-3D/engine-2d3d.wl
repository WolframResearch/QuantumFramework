(* Continuous position measurement on a 2D and 3D grid.

   Companion prototype to ../measurement-in-a-potential.md, whose conventions are binding:
   hbar = mass = 1; ft/ift with FourierParameters -> {0, -1}; the state is a normalized
   rank-d array of amplitudes with Sqrt[dx dy (dz)] absorbed, so Abs[psi]^2 is probability
   per grid point and every expectation is a dot product; one measurement channel applies
   the Gaussian Kraus reweighting Exp[-lambda dt (c - out)^2] in its measured coordinate c
   and renormalizes; its outcome is the measured coordinate of ONE joint position sample
   drawn from Abs[psi]^2, plus Gaussian noise of standard deviation 1/Sqrt[4 lambda dt].

   A channel measures the coordinate combination u . r at strength lambda, with u a unit
   vector: axis-aligned u recovers independent per-axis monitoring, a tilted u measures a
   single rotated quadrature. All simultaneous channels share the SAME joint position
   sample per step (their outcomes are correlated through the position), with independent
   per-channel noises; the per-channel Kraus masks commute, so their exponents add.

   Grid-shaped symbols carry a 2 (rank-2 state) or 3 (rank-3 state) suffix so this file can
   be loaded next to the 1D essay's definitions without pattern collisions; hbar, mass, ft,
   and ift are shared names redefined identically to the essay's. *)

hbar = 1; mass = 1;

(* Rank-agnostic transform pair and state norm. Fourier on a rank-d array is the
   d-dimensional transform, one 1D transform along each axis, so the essay's ft/ift carry
   over verbatim. Normalize does not accept rank-2 arrays (and Norm on a matrix is the
   spectral norm, not the amplitude 2-norm), so the renormalization is spelled out on the
   flattened amplitude list. *)
ft[psi_] := Fourier[psi, FourierParameters -> {0, -1}];
ift[c_] := InverseFourier[c, FourierParameters -> {0, -1}];
stateNorm[psi_] := Norm[Flatten[psi]];
normalizeState[psi_] := psi/stateNorm[psi];
normError2[psi_] := Abs[stateNorm[psi] - 1];

(* Probability that has reached the outer tenth of an axis, fed a coordinate list and the
   matching marginal; identical to the essay's helper. *)
outerTenth2[coord_, probs_] := Total@Pick[probs, UnitStep[Abs[coord] - 0.9 Max@Abs[coord]], 1];

(* ------------------------------ 2D grid ------------------------------ *)

(* Per-axis coordinates and momenta exactly as the 1D grid builds them, momenta stored in
   Fourier's rotated order axis by axis. The rank-2 arrays X, Y, PX, PY put each per-axis
   list into full-grid orientation: dimension 1 is the x index, dimension 2 the y index,
   and p2[[i, j]] = px[[i]]^2 + py[[j]]^2 is the kinetic quadratic form in that same
   orientation. *)
grid2[{nx_, ny_}, {lxIn_, lyIn_}] := With[{lx = N[lxIn], ly = N[lyIn]}, With[{
    xs = (lx/nx) Range[-nx/2, nx/2 - 1], ys = (ly/ny) Range[-ny/2, ny/2 - 1],
    px = (2 Pi hbar/lx) RotateLeft[Range[-nx/2, nx/2 - 1], nx/2],
    py = (2 Pi hbar/ly) RotateLeft[Range[-ny/2, ny/2 - 1], ny/2]},
   With[{xa = Outer[Plus, xs, 0. ys], ya = Outer[Plus, 0. xs, ys]},
    <|"n" -> {nx, ny}, "dx" -> {lx/nx, ly/ny}, "x" -> xs, "y" -> ys, "px" -> px, "py" -> py,
      "X" -> xa, "Y" -> ya, "rflat" -> Transpose[{Flatten[xa], Flatten[ya]}],
      "PX" -> Outer[Plus, px, 0. py], "PY" -> Outer[Plus, 0. px, py],
      "p2" -> Outer[Plus, px^2, py^2]|>]]];

(* Product Gaussian packet with per-axis centers, mean momenta, and position widths. *)
gaussian2[g_, {x0_, y0_}, {px0_, py0_}, {sx_, sy_}] := normalizeState[
   Exp[-(g["X"] - x0)^2/(4 sx^2) - (g["Y"] - y0)^2/(4 sy^2) + I (px0 g["X"] + py0 g["Y"])/hbar]];

(* ------------------------------ channels ------------------------------ *)

(* channel[u, lambda] measures u . r at strength lambda. setupChannels2 normalizes u,
   precomputes the measured coordinate as a full-grid array, and drops lambda = 0 entries:
   an unmonitored combination contributes no mask and no outcome (its outcome noise
   1/Sqrt[4 lambda dt] would be infinite). Strengths must be numeric and provably
   nonnegative real (NumericQ, Element[.., Reals], >= 0; exact constants like Pi
   qualify); an undetermined symbol, a negative value, or a complex value returns a
   Failure instead of being dropped like a zero (undetermined strengths belong in raw
   channel associations fed to the filter's riccatiRHS2, which never goes through this
   setup). axisChannels2 is the common per-axis case. *)
channel[u_, lambda_] := <|"u" -> Normalize[N[u]], "lambda" -> lambda|>;
setupChannels2[g_, chans_] := Enclose[
   ConfirmAssert[VectorQ[Lookup[chans, "lambda"], NumericQ[#] && Element[#, Reals] && # >= 0 &],
    "channel strengths must be nonnegative real numbers"];
   Map[Append[#, "c" -> #["u"][[1]] g["X"] + #["u"][[2]] g["Y"]] &,
    Select[chans, #["lambda"] > 0 &]]];
axisChannels2[g_, {lx_, ly_}] := setupChannels2[g, {channel[{1, 0}, lx], channel[{0, 1}, ly]}];

(* ------------------------------ outcomes ------------------------------ *)

(* One joint position sample from Abs[psi]^2: the flattened density weights the cached
   list of coordinate tuples directly, so no index decode exists to get wrong. *)
jointSample2[g_][psi_] := RandomChoice[Flatten[Abs[psi]^2] -> g["rflat"]];

drawOutcomes2[g_, chans_List, dt_][psi_] := With[{r = jointSample2[g][psi]},
   Map[#["u"] . r + RandomVariate[NormalDistribution[0, 1/Sqrt[4 #["lambda"] dt]]] &, chans]];

(* ------------------------------ measurement update ------------------------------ *)

(* All channels at once: the per-channel Gaussian exponents add, then one Exp, one
   pointwise product, one renormalization. The chans_List pattern keeps a Failure from
   setupChannels2 from flowing onward as if it were an empty channel list. *)
measUpdate2[g_, chans_List, dt_, outs_][psi_] := normalizeState[
   Exp[-dt Total[MapThread[#1["lambda"] (#1["c"] - #2)^2 &, {chans, outs}]]] psi];

(* ------------------------------ unitary step (Strang) ------------------------------ *)

(* The potential phase is pointwise on the array; Vfun must accept the full X and Y arrays
   (any arithmetic/listable formula does). The kinetic phase is pointwise in momentum
   space. Half a potential kick around a full kinetic drift, as in the essay. *)
potentialStep2[g_, Vfun_, dt_][psi_] := Exp[-I Vfun[g["X"], g["Y"]] dt/(2 hbar)] psi;
kineticStep2[g_, dt_][psi_] := ift[Exp[-I g["p2"] dt/(2 mass hbar)] ft[psi]];
evolveStep2[g_, Vfun_, dt_] := potentialStep2[g, Vfun, dt] @* kineticStep2[g, dt] @* potentialStep2[g, Vfun, dt];

(* ------------------------------ conditional step and runners ------------------------------ *)

step2[g_, Vfun_, chans_, dt_, outs_][psi_] := evolveStep2[g, Vfun, dt][measUpdate2[g, chans, dt, outs][psi]];

(* trajectory2 draws fresh outcomes each step and Sows the per-channel outcome vector;
   called inside Reap it returns {states, {outcomes}}. replay2 consumes a given list of
   per-step outcome vectors instead: the conditional state is a deterministic functional
   of the record, which is what makes record-driven cross-checks exact. *)
trajectory2[g_, Vfun_, chans_, dt_, nt_][psi0_] :=
   NestList[step2[g, Vfun, chans, dt, Sow[drawOutcomes2[g, chans, dt][#]]][#] &, psi0, nt];
replay2[g_, Vfun_, chans_, dt_, outsList_][psi0_] :=
   FoldList[step2[g, Vfun, chans, dt, #2][#1] &, psi0, outsList];

(* Long runs that only need moments: Nest the state, Sow {outcomes, moments2} per step. *)
momentsRun2[g_, Vfun_, chans_, dt_, nt_][psi0_] := Nest[
   Function[p, With[{outs = drawOutcomes2[g, chans, dt][p]},
     With[{next = step2[g, Vfun, chans, dt, outs][p]}, Sow[{outs, moments2[g, next]}]; next]]],
   psi0, nt];

(* ------------------------------ conditional moments ------------------------------ *)

(* Means and the full covariance of a quadrature vector R: Sigma_ab = Re<R_a R_b> -
   <R_a><R_b>, the symmetrized covariance. Each operator is applied to psi once, the
   position ones pointwise and the momentum ones in the transform; flattening turns the
   whole moment table into one Gram matrix of applied states, at any rank and for any
   number of quadratures. *)
momentsOf[psi_, ops_] := With[{fp = Flatten[psi], opsF = Map[Flatten, ops]},
   With[{mus = Re[Conjugate[fp] . Transpose[opsF]]},
    <|"mean" -> mus,
      "cov" -> Re[Conjugate[opsF] . Transpose[opsF]] - Outer[Times, mus, mus]|>]];
moments2[g_, psi_] := With[{c = ft[psi]},
   momentsOf[psi, {g["X"] psi, g["Y"] psi, ift[g["PX"] c], ift[g["PY"] c]}]];
moments3[g_, psi_] := With[{c = ft[psi]},
   momentsOf[psi, {g["X"] psi, g["Y"] psi, g["Z"] psi,
     ift[g["PX"] c], ift[g["PY"] c], ift[g["PZ"] c]}]];

(* Per-axis edge and momentum-tail leakage, the essay's two grid health checks read off
   the per-axis marginals. *)
edgeMass2[g_, psi_] := With[{P = Abs[psi]^2},
   {outerTenth2[g["x"], Total[P, {2}]], outerTenth2[g["y"], Total[P]]}];
tailMass2[g_, psi_] := With[{Pc = Abs[ft[psi]]^2},
   {outerTenth2[g["px"], Total[Pc, {2}]], outerTenth2[g["py"], Total[Pc]]}];

(* ------------------------------ multimode Gaussian filter ------------------------------ *)

(* Phase-space order z = (x, y, px, py); the potential is the quadratic form
   V = (1/2) {x, y} . K . {x, y}, so the drift matrix is A = {{0, 1/m}, {-K, 0}} in
   2x2 blocks. A channel along u enters twice: its position embedding ut = (u, 0) carries
   the information gain, its momentum embedding vt = (0, u) carries the backaction
   diffusion. The covariance obeys the deterministic matrix Riccati flow

      dSigma/dt = A.Sigma + Sigma.At + Sum_i lambda_i hbar^2 vt_i vt_i^T
                  - Sum_i 4 lambda_i (Sigma.ut_i)(Sigma.ut_i)^T,

   record-blind, while the outcomes drive only the means,

      dz = A.z dt + Sum_i 4 lambda_i (Sigma.ut_i) (out_i - ut_i . z) dt.

   Restricted to one axis these are verbatim the essay's Riccati and its filterStep.
   The constructors read the mode count off their arguments, so the same filter code
   serves any number of modes (a 3x3 stiffness matrix gives the 6x6 drift). A channel
   may carry a detection efficiency "eta" (default 1): the information terms, the
   -Sigma M Sigma contraction and the innovation gain, scale as eta lambda, while the
   backaction diffusion keeps the bare lambda, the essay's Part IV rescaling; at
   eta < 1 the filter covariance is a mixed state's, above the pure floor. One record
   caveat: drawOutcomes2/3 generate unit-efficiency records (outcome noise
   1/Sqrt[4 lambda dt]); an eta < 1 filter is meant for a record whose outcome noise is
   1/Sqrt[4 eta lambda dt], i.e. a lossy detector's stream, which the pure-state grid
   engine does not generate. *)
driftMatrix2[K_] := ArrayFlatten[{{0, IdentityMatrix[Length[K]]/mass}, {-K, 0}}];
posEmbed2[u_] := Join[u, 0 u];
momEmbed2[u_] := Join[0 u, u];
etaOf[ch_] := Lookup[ch, "eta", 1];
riccatiRHS2[K_, chans_][sig_] := With[{a = driftMatrix2[K]},
   a . sig + sig . Transpose[a]
    + Total[Map[hbar^2 #["lambda"] Outer[Times, momEmbed2[#["u"]], momEmbed2[#["u"]]] &, chans]]
    - Total[Map[4 etaOf[#] #["lambda"] Outer[Times, sig . posEmbed2[#["u"]], sig . posEmbed2[#["u"]]] &, chans]]];
filterStep2[K_, chans_, dt_][{z_, sig_}, outs_] := {
   z + driftMatrix2[K] . z dt
     + Total[MapThread[4 etaOf[#1] #1["lambda"] (sig . posEmbed2[#1["u"]]) (#2 - posEmbed2[#1["u"]] . z) dt &,
        {chans, outs}]],
   sig + riccatiRHS2[K, chans][sig] dt};
filterRun2[K_, chans_, dt_, outsList_][{z0_, sig0_}] :=
   FoldList[filterStep2[K, chans, dt], {z0, sig0}, outsList];

(* ------------------------------ 3D core step ------------------------------ *)

(* The same construction one axis longer. Every rank-2 statement above generalizes by
   adding the z axis to the Outer products, the channel combination, the sample decode,
   and the kinetic quadratic form; ft/ift and the masks are rank-agnostic already. *)
grid3[{nx_, ny_, nz_}, {lxIn_, lyIn_, lzIn_}] := With[{lx = N[lxIn], ly = N[lyIn], lz = N[lzIn]}, With[{
    xs = (lx/nx) Range[-nx/2, nx/2 - 1], ys = (ly/ny) Range[-ny/2, ny/2 - 1],
    zs = (lz/nz) Range[-nz/2, nz/2 - 1],
    px = (2 Pi hbar/lx) RotateLeft[Range[-nx/2, nx/2 - 1], nx/2],
    py = (2 Pi hbar/ly) RotateLeft[Range[-ny/2, ny/2 - 1], ny/2],
    pz = (2 Pi hbar/lz) RotateLeft[Range[-nz/2, nz/2 - 1], nz/2]},
   With[{xa = Outer[Plus, xs, 0. ys, 0. zs], ya = Outer[Plus, 0. xs, ys, 0. zs],
     za = Outer[Plus, 0. xs, 0. ys, zs]},
    <|"n" -> {nx, ny, nz}, "dx" -> {lx/nx, ly/ny, lz/nz},
      "x" -> xs, "y" -> ys, "z" -> zs, "px" -> px, "py" -> py, "pz" -> pz,
      "X" -> xa, "Y" -> ya, "Z" -> za,
      "rflat" -> Transpose[{Flatten[xa], Flatten[ya], Flatten[za]}],
      "PX" -> Outer[Plus, px, 0. py, 0. pz], "PY" -> Outer[Plus, 0. px, py, 0. pz],
      "PZ" -> Outer[Plus, 0. px, 0. py, pz],
      "p2" -> Outer[Plus, px^2, py^2, pz^2]|>]]];

gaussian3[g_, r0_, p0_, s_] := normalizeState[Exp[
    -(g["X"] - r0[[1]])^2/(4 s[[1]]^2) - (g["Y"] - r0[[2]])^2/(4 s[[2]]^2)
     - (g["Z"] - r0[[3]])^2/(4 s[[3]]^2)
     + I (p0[[1]] g["X"] + p0[[2]] g["Y"] + p0[[3]] g["Z"])/hbar]];

setupChannels3[g_, chans_] := Enclose[
   ConfirmAssert[VectorQ[Lookup[chans, "lambda"], NumericQ[#] && Element[#, Reals] && # >= 0 &],
    "channel strengths must be nonnegative real numbers"];
   Map[Append[#, "c" -> #["u"][[1]] g["X"] + #["u"][[2]] g["Y"] + #["u"][[3]] g["Z"]] &,
    Select[chans, #["lambda"] > 0 &]]];
axisChannels3[g_, {lx_, ly_, lz_}] := setupChannels3[g,
   {channel[{1, 0, 0}, lx], channel[{0, 1, 0}, ly], channel[{0, 0, 1}, lz]}];

jointSample3[g_][psi_] := RandomChoice[Flatten[Abs[psi]^2] -> g["rflat"]];
drawOutcomes3[g_, chans_List, dt_][psi_] := With[{r = jointSample3[g][psi]},
   Map[#["u"] . r + RandomVariate[NormalDistribution[0, 1/Sqrt[4 #["lambda"] dt]]] &, chans]];

measUpdate3[g_, chans_List, dt_, outs_][psi_] := normalizeState[
   Exp[-dt Total[MapThread[#1["lambda"] (#1["c"] - #2)^2 &, {chans, outs}]]] psi];

potentialStep3[g_, Vfun_, dt_][psi_] := Exp[-I Vfun[g["X"], g["Y"], g["Z"]] dt/(2 hbar)] psi;
kineticStep3[g_, dt_][psi_] := ift[Exp[-I g["p2"] dt/(2 mass hbar)] ft[psi]];
evolveStep3[g_, Vfun_, dt_] := potentialStep3[g, Vfun, dt] @* kineticStep3[g, dt] @* potentialStep3[g, Vfun, dt];

step3[g_, Vfun_, chans_, dt_, outs_][psi_] := evolveStep3[g, Vfun, dt][measUpdate3[g, chans, dt, outs][psi]];

(* Per-axis edge and momentum-tail guards for the rank-3 state, the essay's honesty
   conditions one axis longer. *)
edgeMass3[g_, psi_] := With[{P = Abs[psi]^2},
   {outerTenth2[g["x"], Total[P, {2, 3}]], outerTenth2[g["y"], Total[Total[P], {2}]],
    outerTenth2[g["z"], Total[P, {1, 2}]]}];
tailMass3[g_, psi_] := With[{Pc = Abs[ft[psi]]^2},
   {outerTenth2[g["px"], Total[Pc, {2, 3}]], outerTenth2[g["py"], Total[Total[Pc], {2}]],
    outerTenth2[g["pz"], Total[Pc, {1, 2}]]}];
