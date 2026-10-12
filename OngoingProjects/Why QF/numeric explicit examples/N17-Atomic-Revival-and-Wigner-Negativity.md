---
Template: Default
Title: "How much cavity Wigner negativity can the atomic population record guarantee?"
---

# How much cavity Wigner negativity can the atomic population record guarantee?

A two-level atom exchanges excitations with a single cavity mode. If the field begins in a coherent state, its photon-number components drive the atom at different Rabi frequencies. The atomic oscillations first lose contrast as those components slip out of phase, then recover when some of their phases realign. This is collapse and revival.

In a cavity experiment the cheapest thing to record is the atom's population: prepare the atom and field, let them interact for a time $t$, detect whether the atom is excited, and repeat. The result is the inversion $z(t)=P_e-P_g$ as a function of time. The cavity field itself is harder to see. Its Wigner function can take negative values, a signature that no classical distribution over phase space reproduces, and measuring it requires phase-space tomography of the field.

So the question of this essay is: **given only the record $z(t)$, how much Wigner negativity of the cavity field can we guarantee?** We measure negativity by the negative volume of the Wigner function and ask for the largest value that every field compatible with the record must reach. Throughout, "the field" means the **unconditional cavity state**: we discard the atom without selecting a measurement outcome.

The answer comes in three steps, each computed below. First, a symmetry of the coupling makes $z(t)$ blind to everything about the field except its photon-number distribution. Second, a field with the same photon-number distribution but a random phase is a real preparation with the same record, and its negative volume is a lower bound for every field with that distribution. Third, the record determines the photon-number distribution, so that lower bound is exactly what the record guarantees.

## The atom and field: fixing the physical model

Let $\lvert e\rangle$ and $\lvert g\rangle$ span the atomic space $\mathcal H_a=\mathbb C^2$. The cavity has Fock states $\lvert n\rangle$, with $n$ photons. Numerically we retain $n=0,\ldots,d-1$, so the joint space is $\mathcal H_a\otimes\mathcal H_f^{(d)}=\mathbb C^2\otimes\mathbb C^d$. The atom is always the first subsystem and the field the second.

The operators $\hat\sigma_-=\lvert g\rangle\langle e\rvert$ and $\hat\sigma_+=\hat\sigma_-^\dagger$ lower and raise the atom. The field operator obeys $\hat a\lvert n\rangle=\sqrt n\lvert n-1\rangle$. In the rotating-wave approximation, the laboratory Hamiltonian is

$$
\hat H=\hbar\omega_c\hat I_a\otimes\hat a^\dagger\hat a
+\frac{\hbar\omega_a}{2}\hat\sigma_z\otimes\hat I_f
+\hbar g\left(\hat\sigma_+\otimes\hat a+\hat\sigma_-\otimes\hat a^\dagger\right),
\qquad
\hat\sigma_z=\lvert e\rangle\langle e\rvert-\lvert g\rangle\langle g\rvert.
$$

The first two terms are the uncoupled field and atomic energies. The interaction converts $\lvert e,n\rangle$ into $\lvert g,n+1\rangle$ and back, and the population in that pair oscillates at angular frequency $2g\sqrt{n+1}$. We assume resonance, $\omega_a=\omega_c$, and $g\sqrt{\bar n+1}\ll\omega_c$, as required for the rotating-wave model. In the interaction picture the free energies disappear, leaving $\hat H_I=\hbar g\hat h$ with

$$
\hat h=\hat\sigma_-\otimes\hat a^\dagger+\hat\sigma_+\otimes\hat a.
$$

This is the one-atom specialization of the Hamiltonian in [Jarvis et al., Eq. (1)](https://arxiv.org/html/0812.3766v2). We measure time by $\tau=gt$.

Load the checked-out QF paclet, then its quantum and bosonic constructors:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"]
```

Load the main context:

```wl
Needs["Wolfram`QuantumFramework`"]
```

Expose the bosonic context:

```wl
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

A coherent state has amplitudes $c_n=e^{-|\alpha|^2/2}\alpha^n/\sqrt{n!}$ and mean photon number $\bar n=|\alpha|^2$. Set $\alpha=2$, so $\bar n=4$, and retain 24 Fock levels; a larger cutoff will check this choice at the end:

```wl
{alpha, d} = {2., 24}
```

Build QF's normalized truncated coherent state:

```wl
coherent = CoherentState[d][alpha]
```

Embed photon annihilation on the full product space. This operator will also describe a photon leaving the cavity:

```wl
a = QuantumTensorProduct[QuantumOperator["I"], AnnihilationOperator[d]]
```

Use the computational labels $\lvert0\rangle_a=\lvert e\rangle$ and $\lvert1\rangle_a=\lvert g\rangle$, and construct atomic lowering on subsystem 1 in this excited-first convention:

```wl
sigmaMinus = QuantumOperator[{{0, 0}, {1, 0}}, {1}, 2]
```

For operators, QF's `f[g]` means the composition $fg$; QF aligns the addressed subsystems and carries the untouched subsystem through. Form the dimensionless interaction $\hat h$; `SuperDagger` takes the adjoint:

```wl
h = sigmaMinus[SuperDagger[a]] + SuperDagger[sigmaMinus][a]
```

Confirm that $\hat h$ is Hermitian:

```wl
h["HermitianQ"]
```

## What the population record can see: a conserved excitation number

The interaction only trades an atomic excitation for a photon. So the total number of excitations,

$$
\hat N=\hat\sigma_+\hat\sigma_-\otimes\hat I_f+\hat I_a\otimes\hat a^\dagger\hat a,
$$

is conserved. In other words, $\hat h$ never connects $\lvert e,n\rangle$ to anything except $\lvert g,n+1\rangle$, and both carry $N=n+1$ excitations.

Build $\hat N$ on the joint space:

```wl
excitations = SuperDagger[sigmaMinus][sigmaMinus] + SuperDagger[a][a]
```

Verify that $[\hat N,\hat h]=0$; the largest matrix element of the commutator is:

```wl
Max[Abs[(excitations[h] - h[excitations])["Matrix"]]]
```

Photon loss is not conserving: it removes one excitation. We couple the cavity to a memoryless (Markov) vacuum bath and neglect atomic decay and dephasing. Write $\hat A=\hat I_a\otimes\hat a$ and

$$
\frac{d\hat\rho}{d\tau}=-i[\hat h,\hat\rho]+r\,\mathcal D[\hat A]\hat\rho,
\qquad
\mathcal D[\hat A]\hat\rho=\hat A\hat\rho\hat A^\dagger
-\frac12\left(\hat A^\dagger\hat A\hat\rho+\hat\rho\hat A^\dagger\hat A\right),
\qquad r=\kappa/g.
$$

This is the zero-temperature loss equation in [Meunier et al., Eq. (26)](https://arxiv.org/pdf/quant-ph/0603191); their loss rate $\gamma$ is our $\kappa$, while their coupling is twice our $g$. Here $\kappa$ is the photon-number decay rate: without the atom, $\langle\hat n\rangle$ decays as $e^{-\kappa t}$.

The jump $\hat A$ lowers $N$ by exactly one, $[\hat N,\hat A]=-\hat A$. Verify it; the largest matrix element of $[\hat N,\hat A]+\hat A$ is:

```wl
Max[Abs[(excitations[a] - a[excitations] + a)["Matrix"]]]
```

Now split the joint density matrix into blocks $\hat\rho_{NM}$, the part that maps excitation number $M$ to excitation number $N$. The commutator with $\hat h$ and the terms with $\hat A^\dagger\hat A$, which both conserve $N$, keep every block where it is, and $\hat A\hat\rho\hat A^\dagger$ moves block $(N,M)$ to $(N-1,M-1)$. So the diagonal blocks $\hat\rho_{NN}$ evolve among themselves, never fed by the off-diagonal ones. The observable $\hat\sigma_z\otimes\hat I_f$ is diagonal in $N$, so $z(\tau)$ reads only those diagonal blocks.

In other words: **$z(\tau)$ depends only on the diagonal blocks of the initial joint state.** With the atom starting in $\lvert e\rangle$, the diagonal block with $N$ excitations is $p_{N-1}\lvert e,N-1\rangle\langle e,N-1\rvert$, where $p_n=|c_n|^2$ is the photon-number distribution. Everything else about the field, its phase and every coherence between different photon numbers, sits in the off-diagonal blocks, which the record never sees. The truncated $\hat a$ still lowers $n$ by exactly one, so this holds exactly in the truncated model, not only in the limit $d\to\infty$.

The assumption that the atom starts in $\lvert e\rangle$ matters; we return to it in the last section.

## Two field preparations, one atomic record

The symmetry argument predicts that any two fields with the same $p_n$ produce the same record. The simplest partner of a coherent state is a **phase-randomized coherent state**: a coherent state whose phase is unknown and uniformly distributed, as for a laser with no phase reference. Averaging $\lvert\alpha e^{i\phi}\rangle\langle\alpha e^{i\phi}\rvert$ over $\phi$ removes every coherence $\rho_{mn}$ with $m\neq n$ and keeps the Poisson weights:

$$
\hat\rho_{\text{rand}}=\sum_n p_n\lvert n\rangle\langle n\rvert,\qquad p_n=e^{-\bar n}\frac{\bar n^n}{n!}.
$$

Read the photon-number distribution of the truncated coherent state:

```wl
photonProbabilities = coherent["ProbabilitiesList"]
```

Build the phase-randomized field as a density matrix diagonal in the Fock basis:

```wl
dephased = QuantumState[DiagonalMatrix[photonProbabilities], d]
```

Prepare the atom in $\lvert e\rangle$ next to each field:

```wl
inputs = QuantumTensorProduct[QuantumState["0"], #] & /@ {coherent, dephased}
```

Choose five loss rates, from a lossless cavity to a photon lifetime $g/\kappa=5$:

```wl
losses = {0., .02, .05, .1, .2}
```

Evolve both preparations at every loss rate. QF takes the generator, the jump operator with its rate, and the initial state directly; the rate multiplies the dissipator, so it is equivalent to a jump $\sqrt r\,\hat A$:

```wl
trajectories = Table[
    QuantumEvolve[h, {a} -> {r}, input, {tau, 0., 16.}, AccuracyGoal -> 9, PrecisionGoal -> 9],
    {input, inputs}, {r, losses}
];
```

The large-$\bar n$ estimate of the first revival is $\tau_r\simeq2\pi\sqrt{\bar n}=4\pi$, so the window $0\leq\tau\leq16$ contains the collapse and the first revival. Sample it every $0.1$:

```wl
times = N[Range[0, 160]/10];
```

Evaluate every trajectory on that time grid:

```wl
jointStates = Map[Function[evolution, evolution /@ times], trajectories, {2}];
```

The atomic state is $\hat\rho_a=\operatorname{Tr}_f\hat\rho$; in QF the listed subsystem is the one discarded, so trace out subsystem 2. With the excited-first labels, $\hat\sigma_z$ is `PauliMatrix[3]`. Compute $z=\operatorname{Tr}(\hat\rho_a\hat\sigma_z)$ for every preparation, loss rate and time:

```wl
inversion = Map[Re @ Tr[PauliMatrix[3] . QuantumPartialTrace[#, {2}]["DensityMatrix"]] &, jointStates, {3}];
```

Compare the two records: for each loss rate, the largest difference between the coherent and phase-randomized $z(\tau)$ is:

```wl
AssociationThread[losses, Max /@ Abs[inversion[[1]] - inversion[[2]]]]
```

As one can see, the two records agree to the accuracy of the time integration at every loss rate, including the strong-loss ones where the field loses most of its energy. Plot them for a lossless and a lossy cavity; solid lines are the coherent field, dashed lines the phase-randomized one:

```wl
ListLinePlot[
    Flatten[Table[Transpose[{times, inversion[[input, k]]}], {k, {1, 3}}, {input, 2}], 1],
    PlotStyle -> {
        Directive[StandardBlue, AbsoluteThickness[4], Opacity[.45]], Directive[StandardBlue, Dashed, AbsoluteThickness[1.5]],
        Directive[StandardOrange, AbsoluteThickness[4], Opacity[.45]], Directive[StandardOrange, Dashed, AbsoluteThickness[1.5]]
    },
    PlotLegends -> Placed[{"coherent, \[Kappa]/g = 0", "phase-randomized, \[Kappa]/g = 0", "coherent, \[Kappa]/g = 0.05", "phase-randomized, \[Kappa]/g = 0.05"}, Right],
    Frame -> True, FrameLabel -> {"g t", "z = \!\(\*SubscriptBox[\(P\), \(e\)]\) - \!\(\*SubscriptBox[\(P\), \(g\)]\)"}, PlotRange -> {-1, 1},
    PlotLabel -> "Atomic population record for two field preparations", ImageSize -> 720, AspectRatio -> 1/3
]
```

![The atomic inversion for the coherent and the phase-randomized field, lossless and lossy.](img/n17/same-record.png)

The collapse, the revival, and the way loss erodes the revival are all identical. An experimenter holding this record cannot tell which field was in the cavity.

The fields themselves differ. For the unconditional cavity state $\hat\rho_f=\operatorname{Tr}_a\hat\rho$, discard subsystem 1 and keep the density matrices:

```wl
fieldMatrices = Map[Normal @ QuantumPartialTrace[#, {1}]["DensityMatrix"] &, jointStates, {3}];
```

Recall that a field element $\langle n\rvert\hat\rho_f\lvert m\rangle=\sum_s\langle s,n\rvert\hat\rho\lvert s,m\rangle$ with $n\neq m$ comes from blocks with different $N$. So the phase-randomized field at every later time should be exactly the Fock-diagonal part of the coherent field. Verify it; the largest difference over all loss rates and times is:

```wl
Max[Abs[fieldMatrices[[2]] - Map[DiagonalMatrix @* Diagonal, fieldMatrices[[1]], {2}]]]
```

In short, the record carries the photon-number distribution $p_n(\tau)$ of the field and nothing of its photon-number coherences.

## The Wigner function in angular harmonics: what the record misses

Use dimensionless quadratures $\hat x=(\hat a+\hat a^\dagger)/\sqrt2$ and $\hat p=(\hat a-\hat a^\dagger)/(i\sqrt2)$, normalize the Wigner function by $\int W\,dx\,dp=1$, so that vacuum has $W_0(x,p)=e^{-(x^2+p^2)}/\pi$, and switch to polar coordinates $x=r\cos\theta$, $p=r\sin\theta$. Writing $\hat\rho_f=\sum_{m,n}\rho_{mn}\lvert m\rangle\langle n\rvert$ and using the Wigner function of each $\lvert m\rangle\langle n\rvert$, the Wigner function separates into angular harmonics:

$$
W(r,\theta)=\sum_{k\geq0}c_k\,\mathrm{Re}\!\left[e^{-ik\theta}\,W_k(r)\right],
\qquad
W_k(r)=\sum_n\rho_{n+k,n}\,R_{nk}(r),
\qquad c_0=1,\ c_{k\geq1}=2,
$$

$$
R_{nk}(r)=\frac{(-1)^n}{\pi}\sqrt{\frac{n!}{(n+k)!}}\,(\sqrt2\,r)^k\,L_n^{(k)}(2r^2)\,e^{-r^2},
$$

with $L_n^{(k)}$ the generalized Laguerre polynomial and $(\sqrt2\,r)^0=1$, including at the origin. In words: the $k$-th harmonic, the part of $W$ that winds $k$ times around the origin, is built only from the coherences $\rho_{n+k,n}$ between photon numbers that differ by $k$. The populations build $W_0(r)$ alone, which is rotationally symmetric. Averaging $W$ over the angle keeps exactly $W_0$, and averaging the field over its phase does the same thing: the Wigner function of the phase-randomized field is the $k=0$ harmonic of the original one.

Define the radial functions $R_{nk}$ at a list of radii, grouped by $k$:

```wl
radialFunctions[r_] := Table[
    N[(-1)^n Sqrt[n!/(n + k)!] If[k == 0, 1, (Sqrt[2] r)^k] LaguerreL[n, k, 2 r^2] Exp[-r^2]/Pi],
    {k, 0, d - 1}, {n, 0, d - 1 - k}
]
```

Define the angular factors $c_ke^{-ik\theta}$ at a list of angles:

```wl
harmonicWeights[theta_] := Table[If[k == 0, 1, 2] Exp[-I k theta], {k, 0, d - 1}]
```

Collect the harmonics $W_k(r)$ of a field density matrix from its diagonals $\rho_{n+k,n}$:

```wl
angularModes[rho_, radial_] := MapThread[Dot, {Table[Diagonal[rho, -k], {k, 0, Length[radial] - 1}], radial}]
```

Define the Wigner function at a list of phase-space points $\{x,p\}$:

```wl
wignerAt[rho_, points_] := With[{r = Norm /@ points, theta = Arg[#1 + I #2] & @@@ points},
    Re[Total[angularModes[rho, radialFunctions[r]] harmonicWeights[theta]]]
]
```

Before trusting this, compare it with QF's own Wigner grid. Take the lossless coherent field at $\tau=6.3$ and build QF's representation on a $40\times40$ grid over $[-5,5]^2$:

```wl
qfWigner = WignerRepresentation[QuantumState[fieldMatrices[[1, 1, 64]], d], {-5., 5.}, {-5., 5.}, "GridSize" -> 40]
```

Evaluate the harmonic sum at QF's grid points; the largest difference from QF's values is:

```wl
Max[Abs[Flatten[qfWigner["ValuesOnGrid"]] - wignerAt[fieldMatrices[[1, 1, 64]], Flatten[qfWigner["Grid"], 1]]]]
```

As expected, the two agree to machine precision: the harmonic sum is the same Wigner function, organized by angular winding.

The negative volume

$$
\mathcal N_W=\int\max(0,-W)\,dx\,dp=\int_0^\infty 2\pi r\,dr\,\big\langle\max(0,-W)\big\rangle_\theta
$$

is half the negativity indicator $\delta$ of [Kenfack and Życzkowski, Eq. (2.2)](https://arxiv.org/html/quant-ph/0406015); here $\langle\cdot\rangle_\theta$ is the average over the angle. It vanishes for a nonnegative Wigner function. That does not imply a classical state: squeezed Gaussian states, for example, also have nonnegative Wigner functions.

Integrating on a polar grid suits the harmonic form. The radial functions are computed once, and each density matrix then costs one matrix product. Use the midpoint rule in $r$ and equally spaced angles; the angular sum is exact for every harmonic the truncated field can carry. Build the grid with its radial weights $2\pi r\,\Delta r$:

```wl
wignerGrid[radius_, radialPoints_, angularPoints_] := With[
    {
        radii = radius/radialPoints (Range[radialPoints] - 1/2),
        angles = N[2 Pi Range[0, angularPoints - 1]/angularPoints]
    },
    <|
        "RadialWeights" -> 2 Pi radius/radialPoints radii,
        "Radial" -> radialFunctions[radii],
        "Harmonics" -> Developer`ToPackedArray[harmonicWeights[angles], Complex]
    |>
]
```

Take radii up to $7$ with 350 points and 256 angles; the checks at the end vary all three:

```wl
grid = wignerGrid[7., 350, 256];
```

Evaluate the Wigner function on the whole grid, radii by angles, as one matrix product:

```wl
wigner[grid_][rho_] := Re[Transpose[angularModes[rho, grid["Radial"]]] . grid["Harmonics"]]
```

Check the normalization $\int W\,dx\,dp=1$ for the same field:

```wl
grid["RadialWeights"] . (Mean /@ wigner[grid][fieldMatrices[[1, 1, 64]]])
```

Define the negative volume: clip the negative part with `Ramp`, then average over the angle:

```wl
negativeVolume[grid_][rho_] := grid["RadialWeights"] . (Mean /@ Ramp[-wigner[grid][rho]])
```

Define the negative volume of the phase-averaged field. It differs from the line above only in the order of the two operations: average over the angle first, then clip:

```wl
phaseAveragedNegativeVolume[grid_][rho_] := grid["RadialWeights"] . Ramp[-(Mean /@ wigner[grid][rho])]
```

That difference in order is the whole bound. The function $u\mapsto\max(0,-u)$ is convex, so the clipped average can never exceed the average of the clipped values:

$$
\max\!\big(0,-\langle W\rangle_\theta\big)\leq\big\langle\max(0,-W)\big\rangle_\theta
\quad\Longrightarrow\quad
\mathcal N_W[W_0]\leq\mathcal N_W[W].
$$

In other words, averaging over the phase can cancel negative regions against positive ones but can never create negativity. So **every field with the photon-number distribution $p_n$ has at least the negative volume of the phase-averaged field**, and the phase-randomized preparation reaches that value exactly.

## The guaranteed negativity against the true negativity

Compute the negative volume of both fields at every loss rate and time:

```wl
negativity = Map[negativeVolume[grid], fieldMatrices, {3}];
```

Confirm that the phase-randomized field's negative volume is the phase-averaged negative volume of the coherent field; the largest difference is:

```wl
Max[Abs[negativity[[2]] - Map[phaseAveragedNegativeVolume[grid], fieldMatrices[[1]], {2}]]]
```

Confirm the bound everywhere; the smallest value of $\mathcal N_W(\text{coherent})-\mathcal N_W(\text{phase-randomized})$ over all loss rates and times is:

```wl
Min[negativity[[1]] - negativity[[2]]]
```

As expected, it is zero up to rounding: the bound never fails. Plot both negative volumes for each loss rate, solid for the coherent field and dashed for the phase-randomized one; the shaded band is the negativity that the record cannot guarantee:

```wl
GraphicsRow[
    MapThread[
        Function[{k, color}, ListLinePlot[
            {Transpose[{times, negativity[[1, k]]}], Transpose[{times, negativity[[2, k]]}]},
            PlotStyle -> {Directive[color, AbsoluteThickness[2]], Directive[color, Dashed, AbsoluteThickness[2]]},
            Filling -> {1 -> {2}}, FillingStyle -> Directive[color, Opacity[.15]],
            Frame -> True, FrameLabel -> {"g t", If[k == 1, "negative volume", None]},
            PlotRange -> {{0, 16}, {0, .37}}, PlotLabel -> Row[{"\[Kappa]/g = ", losses[[k]]}], ImageSize -> 260, AspectRatio -> 1
        ]],
        {Range[5], {StandardBlue, StandardOrange, StandardGreen, StandardRed, StandardPurple}}
    ],
    Spacings -> 0, ImageSize -> 1300
]
```

![Negative volume of the coherent field and of the phase-randomized field for five loss rates.](img/n17/negativity-bound.png)

Have you noticed where the band is widest? In the lossless cavity it is broadest between the collapse and the revival, a little after the large-$\bar n$ estimate $\tau_r/2=2\pi$; at $\bar n=4$ these asymptotic times run early, and the revival in the record also peaks later than $4\pi$. The middle of that interval is where the large-$\bar n$ picture puts the field's two branches farthest apart. For large $\bar n$ the joint state stays close to a sum of two products, each an atomic state times a coherent state whose phase advances by $\pm\tau/(2\sqrt{\bar n})$, so the field splits into two components that rotate in opposite directions around the origin. Near $\tau_r/2$ each has turned by a quarter circle, they sit on opposite sides of the origin, and the field is close to a superposition of two coherent states. Its interference fringes are the negativity, and they are made of coherences between photon numbers: exactly what the record cannot see. Near $\tau_r$ the two components meet again on the far side of phase space, the atom's oscillations return, and the fringe pattern collapses into a smaller structure.

Look at the two fields near the large-$\bar n$ estimates of those two times. Pick the samples nearest $2\pi$ and $4\pi$:

```wl
snapshotIndices = 1 + Round[10 {2 Pi, 4 Pi}]
```

Build QF's Wigner representation of the lossless coherent and phase-randomized fields at those times:

```wl
snapshots = Table[
    WignerRepresentation[QuantumState[fieldMatrices[[input, 1, j]], d], {-5., 5.}, {-5., 5.}, "GridSize" -> 121],
    {input, 2}, {j, snapshotIndices}
];
```

Plot all four on one color scale. A Wigner function never exceeds $1/\pi$ in magnitude, so the scale is $\pm1/\pi$; blue is negative and red positive:

```wl
Grid[
    MapThread[
        Function[{w, label}, ListDensityPlot[
            Flatten[MapThread[Append, {w["Grid"], w["ValuesOnGrid"]}, 2], 1],
            ColorFunction -> (Blend[{StandardBlue, White, StandardRed}, Rescale[#, {-1/Pi, 1/Pi}]] &), ColorFunctionScaling -> False,
            PlotRange -> {-1/Pi, 1/Pi}, PlotLegends -> BarLegend[{(Blend[{StandardBlue, White, StandardRed}, Rescale[#, {-1/Pi, 1/Pi}]] &), {-1/Pi, 1/Pi}}, LegendLabel -> "W"],
            FrameLabel -> {"x", "p"}, PlotLabel -> label, ImageSize -> 330
        ]],
        {snapshots, Outer[Row[{#1, ", g t = ", #2}] &, {"coherent", "phase-randomized"}, times[[snapshotIndices]]]},
        2
    ],
    Spacings -> {1, 1}
]
```

![Wigner functions of the coherent and phase-randomized fields near half the revival time and near the revival.](img/n17/snapshots.png)

The top row shows the branch picture directly. Near $\tau_r/2$ the coherent field has two lobes above and below the origin with fringes between them. Near $\tau_r$ the lobes have merged on the side opposite to where the field started. The bottom row shows the same fields averaged over the angle: rings, whose negative ones come from the photon-number distribution alone. Those rings are all the record can vouch for.

Loss shrinks both curves, but not equally. Compute, for each loss rate, the fraction of the time-integrated negative volume that the record guarantees:

```wl
AssociationThread[losses, (Total /@ negativity[[2]])/(Total /@ negativity[[1]])]
```

As one can see, the guaranteed share falls as loss grows: the negativity that survives photon loss is mostly carried by coherences, and the Fock-diagonal part disappears sooner. Why the populations lose their negativity faster is a question this essay leaves open.

## Does the record determine the photon-number distribution?

The bound above needs $p_n$, and the experimenter is given $z(\tau)$. For zero loss each excitation pair evolves on its own. Starting from the excited atom,

$$
\lvert\Psi(\tau)\rangle=\sum_{n}c_n\left[
\cos(\sqrt{n+1}\,\tau)\lvert e,n\rangle
-i\sin(\sqrt{n+1}\,\tau)\lvert g,n+1\rangle
\right],
\qquad
z_0(\tau)=\sum_n p_n\cos\!\left(2\sqrt{n+1}\,\tau\right).
$$

So the record is a linear function of $p_n$, a sum of cosines whose frequencies $2\sqrt{n+1}$ are all different. Build the matrix of those cosines on the sampled times:

```wl
design = Table[Cos[2 Sqrt[n + 1] t], {t, times}, {n, 0, d - 1}];
```

Check the closed form against the lossless QF record; the largest difference is:

```wl
Max[Abs[inversion[[1, 1]] - design . photonProbabilities]]
```

This also checks the atomic sign, the coupling convention and the lossless trajectory over the whole window. Now invert: recover $p_n$ from the record alone by linear least squares. The largest error in the recovered distribution is:

```wl
Max[Abs[LeastSquares[design, inversion[[1, 1]]] - photonProbabilities]]
```

How long must the record be? Repeat the inversion with records that stop at $\tau=4$, at half the large-$\bar n$ revival time $2\pi$, and at $\tau=8,12,16$; the largest error in each recovered distribution is:

```wl
AssociationThread[#, Function[tmax, With[{rows = Flatten[Position[times, t_ /; t <= tmax]]}, Max[Abs[LeastSquares[design[[rows]], inversion[[1, 1, rows]]] - photonProbabilities]]]] /@ #] &[{4., N[2 Pi], 8., 12., 16.}]
```

Have you noticed the pattern? A record that ends during the collapse cannot determine the distribution at all: during the collapse the cosines are too nearly dependent. At half the revival time the recovery is still rough; by $\tau=8$ the whole distribution is recovered accurately, and the recovery keeps improving toward the revival. These records are free of noise, so the only limit is the conditioning of the cosine matrix against the accuracy of the time integration. A measured record carries shot noise, and then the natural scale is the frequency spacing: near $n=\bar n$ neighbouring Rabi frequencies differ by about $1/\sqrt{\bar n}$, and resolving them by their phase difference takes a record of length about $2\pi\sqrt{\bar n}$, the revival time.

With loss, the diagonal blocks still evolve linearly and on their own, so $z(\tau)$ is still a linear function of the initial $p_n$, and with $\kappa$ known the same inversion applies. I demonstrate the inversion only for the lossless record; for the lossy records the statement rests on this linearity and on knowing $\kappa$.

## Where the blindness ends: an atom with its own phase

The symmetry argument used one more assumption: the atom starts in $\lvert e\rangle$. Suppose instead it starts in a superposition $(\lvert e\rangle+e^{i\varphi}\lvert g\rangle)/\sqrt2$. The diagonal block with $N$ excitations then contains $\lvert e,N-1\rangle$ and $\lvert g,N\rangle$ together, and the coherence between them, $\propto e^{-i\varphi}c_{N-1}c_N^*$, involves the field's coherence between neighbouring photon numbers. Within that block $\hat h$ acts as $\sqrt N$ times a Pauli-$X$ matrix, which rotates the populations into the imaginary part of the coherence. So the record picks up the field phase only through $\mathrm{Im}\big(e^{-i\varphi}c_{N-1}c_N^*\big)$. That requires a relative phase between the atomic dipole and the field.

Define a function that evolves an atomic state next to the coherent and the phase-randomized fields in the lossless cavity and returns the largest difference between their records:

```wl
atomRecordGap[atom_] := Max[Abs[Subtract @@ Map[
    Function[field, With[
        {evolution = QuantumEvolve[h, {a} -> {0.}, QuantumTensorProduct[QuantumState[atom], field], {tau, 0., 16.}, AccuracyGoal -> 9, PrecisionGoal -> 9]},
        Re @ Tr[PauliMatrix[3] . QuantumPartialTrace[evolution[#], {2}]["DensityMatrix"]] & /@ times
    ]],
    {coherent, dephased}
]]]
```

Compare an atom in $(\lvert e\rangle+\lvert g\rangle)/\sqrt2$, QF's `"+"`, with one in $(\lvert e\rangle+i\lvert g\rangle)/\sqrt2$, QF's `"R"`:

```wl
AssociationThread[{"+", "R"}, atomRecordGap /@ {"+", "R"}]
```

As one can see, a superposition alone is not enough. With a real $\alpha$, the atom in `"+"` has no relative phase to the field and its record is still blind; the atom in `"R"` is a quarter turn out of phase with the field, and its record tells the two fields apart. The population record becomes sensitive to the field's coherences only when the atom brings a phase reference of its own. A classical drive would also break the blindness, because it does not conserve $\hat N$. The counter-rotating terms beyond the rotating-wave approximation change $N$ by two, so they conserve only its parity; they would let the record see coherences between photon numbers of equal parity.

## Checking the numerical resolution

Check trace conservation over every evolved joint state; the largest deviation of the trace from one is:

```wl
Max[Abs[Map[#["Norm"] &, jointStates, {3}] - 1]]
```

Refine the phase-space grid: radii up to $8$, twice the radial points and twice the angles. The largest change in the lossless negative volumes of both fields is:

```wl
Max[Abs[Map[negativeVolume[wignerGrid[8., 700, 512]], fieldMatrices[[All, 1]], {2}] - negativity[[All, 1]]]]
```

Raise the Fock cutoff from 24 to 28 and recompute the coherent field at $\kappa/g=0.05$. Every $d$-dependent object is rebuilt inside the `Block`:

```wl
cutoffNegativity = Block[{d = 28},
    With[{field = CoherentState[d][alpha], loss = QuantumTensorProduct[QuantumOperator["I"], AnnihilationOperator[d]]},
        With[
            {
                evolution = QuantumEvolve[sigmaMinus[SuperDagger[loss]] + SuperDagger[sigmaMinus][loss], {loss} -> {.05}, QuantumTensorProduct[QuantumState["0"], field], {tau, 0., 16.}, AccuracyGoal -> 9, PrecisionGoal -> 9],
                cutoffGrid = wignerGrid[7., 350, 256]
            },
            negativeVolume[cutoffGrid][Normal @ QuantumPartialTrace[evolution[#], {1}]["DensityMatrix"]] & /@ times
        ]
    ]
];
```

The largest change in the negative volume is:

```wl
Max[Abs[cutoffNegativity - negativity[[1, 3]]]]
```

These are sensitivity checks at the stated settings, not rigorous bounds on the Wigner function over the whole plane. They are small compared with the differences between the two fields that the argument rests on.

## The answer

The population record of an atom that starts in $\lvert e\rangle$ is blind to everything about the cavity field except its photon-number distribution, because the coupling conserves the excitation number and photon loss only lowers it. A record that extends past the collapse determines that distribution. The phase-averaged field is then a real preparation with the same record, and convexity makes its negative volume a lower bound for every field with the same distribution. So the record guarantees exactly the negative volume of the phase-averaged field, the part of the negativity that the populations produce on their own. The rest, carried by coherences between photon numbers, is invisible to it. That rest is most of the negativity near half the revival time, where the field is close to a superposition of two coherent states, and an increasing share of it as photon loss grows.

<!-- Evaluated cell by cell in Wolfram Language 15.0.1 with QF 2.1.3 (QuantumFramework HEAD f46ec628), 2026-10-05. -->

## Undermind research references

Collection-level research context: N17 was selected from the twenty-question collection associated with the search below. This is not a dedicated search or validation of the result in this essay.

- [Compact computational questions across quantum science](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Compact%20computational%20questions%20across%20quantum%20science): saved candidate-discovery search for the parent collection.

Workspace: [Roadmap for Wolfram quantum](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3); workspace ID: `bcedc042-7a0d-4f3d-b4ca-4e1c138571a3`.
