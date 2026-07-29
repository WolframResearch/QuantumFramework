---
Template: Default
---

# Continuous Measurement in a Potential

**Generalizing the stochastic Schrödinger equation of continuous position measurement from the free particle to an arbitrary one-dimensional Hamiltonian $\hat H = \hat p^2/2m + V(\hat x)$. The only measurement-independent piece of the split-operator scheme is the unitary evolution, so a symmetric (Strang) split-operator step carries any potential, and the free particle is the $V=0$ special case, recovered exactly. Two potentials mark the range. In the harmonic trap the conditional state stays Gaussian: the measurement squeezes its position below the zero-point value, and the five-number moment filter still closes and agrees with the grid. In the double well the moment hierarchy does not close, only the grid split-operator applies, and continuous position monitoring localizes the particle in one well and suppresses coherent tunneling by the quantum Zeno effect. Throughout, the measurement backaction heats the momentum at the state-independent rate $\lambda\hbar^2$, which drives both the trapped oscillator's unbounded heating and, in the double well, the anti-Zeno escape over the barrier that competes with the Zeno suppression.**

Mads Bahrami (last updated: July 28, 2026)

### Setting the Stage

This continues *Trajectories, the Record, and the Master Equation*, which simulates the stochastic Schrödinger equation (SSE) of a free particle under continuous weak position measurement. That document ends by flagging one extension: add a potential to $\hat H$ and study the same record-conditioned dynamics in a trap, the setting of cavity and levitated optomechanics. This is that extension.

The observation that makes it a small change is that the measurement is independent of the Hamiltonian. Over one time step the conditional state is acted on by the Kraus operator $\hat M(\bar x)\propto e^{-\lambda\,dt\,(\hat x-\bar x)^2}$, drawn and renormalized, and then evolved unitarily under $\hat H$ for $dt$. Only the second factor knows about $V$. So generalizing the simulation means generalizing one function, the unitary step, and leaving the measurement machinery untouched. Everything else, the outcome distribution, the renormalization, the conditional-unconditional split, the backaction, transfers verbatim.

Every number below is produced by a cell you can rerun and change; evaluate top to bottom and vary the potential, the measurement strength $\lambda$, and the initial state. As before $\hbar$, the mass $m$, and $\lambda$ stay symbolic wherever the algebra is symbolic, and take numerical values (`hbar`, `mass`) only in simulation.

## Part I: One Function Carries the Potential

We reuse the free-particle machinery. Fix units $\hbar=m=1$, build the position grid with its wrapped momenta, a Gaussian packet, the free (spectral) unitary step, the two conditional moments, the exact outcome draw, and the Kraus update:

```wl
hbar = 1.; mass = 1.;
grid[n_, ll_] := With[{dx = ll/n},
   <|"n" -> n, "dx" -> dx, "x" -> dx Range[-n/2, n/2 - 1],
     "p" -> (2 Pi hbar/ll) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
gaussian[g_, x0_, p0_, s_] := Normalize[Exp[-(g["x"] - x0)^2/(4 s^2) + I p0 g["x"]/hbar]];
unitaryStep[g_, dt_][ψ_] := InverseFourier[
   Exp[-I g["p"]^2 dt/(2 mass hbar)] Fourier[ψ, FourierParameters -> {0, -1}],
   FourierParameters -> {0, -1}];
xMean[g_, ψ_] := g["x"] . Abs[ψ]^2;
xVar[g_, ψ_] := g["x"]^2 . Abs[ψ]^2 - xMean[g, ψ]^2;
drawExact[g_, λ_, dt_][ψ_] := RandomChoice[Abs[ψ]^2 -> g["x"]] +
   RandomVariate[NormalDistribution[0, 1/Sqrt[4 λ dt]]];
measUpdate[g_, λ_, dt_, xb_][ψ_] := Normalize[Exp[-λ dt (g["x"] - xb)^2] ψ];
normError[ψ_] := Abs[ψ . Conjugate[ψ] - 1];
projector[ψ_] := Outer[Times, ψ, Conjugate@ψ];
```

For $\hat H = \hat p^2/2m + V(\hat x)$ the potential is diagonal in position, so its propagator $e^{-iV(\hat x)\,dt/\hbar}$ is a pointwise phase. Splitting the exponential symmetrically, $e^{-i\hat H\,dt/\hbar}\approx e^{-iV\,dt/2\hbar}\,e^{-i\hat T\,dt/\hbar}\,e^{-iV\,dt/2\hbar}$, gives a second-order (Strang) step: a half phase kick, a full free (kinetic) step, another half phase kick. The half kick is one line, and the full evolution is the composition, with the free `unitaryStep` reused unchanged as the kinetic factor:

```wl
potentialStep[g_, Vfun_, dt_][ψ_] := Exp[-I Vfun[g["x"]] dt/(2 hbar)] ψ;
evolveStep[g_, Vfun_, dt_] :=
   potentialStep[g, Vfun, dt] @* unitaryStep[g, dt] @* potentialStep[g, Vfun, dt];
```

The conditional step is measurement then evolution, exactly as before, but with `evolveStep` in place of `unitaryStep`, and one run reaps the record while nesting the trajectory:

```wl
stepV[g_, Vfun_, λ_, dt_, xb_][ψ_] := evolveStep[g, Vfun, dt][measUpdate[g, λ, dt, xb][ψ]];
trajectory[g_, Vfun_, λ_, dt_, nt_][ψ0_] :=
   NestList[stepV[g, Vfun, λ, dt, drawExact[g, λ, dt][#]][#] &, ψ0, nt];
exptV[g_, Vfun_, λ_, dt_, nt_][ψ0_] := Reap[
   NestList[stepV[g, Vfun, λ, dt, Sow@drawExact[g, λ, dt][#]][#] &, ψ0, nt]];
```

The potential step is a pure phase, $|e^{-iV\,dt/2\hbar}|=1$, so it preserves the norm exactly regardless of how large or stiff $V$ is; the whole `evolveStep` is unitary to machine precision on any potential. The consistency check is the zero potential: with $V\equiv 0$ each phase kick is the identity and `evolveStep` must reduce to the free `unitaryStep` bit for bit.

```wl
gChk = grid[128, 40.]; ψChk = gaussian[gChk, 1.5, 0.7, 1.3];
{"max |evolveStep(V=0) - unitaryStep|" -> Max@Abs[evolveStep[gChk, 0 &, 0.01][ψChk] - unitaryStep[gChk, 0.01][ψChk]],
 "identical expression" -> (evolveStep[gChk, 0 &, 0.01][ψChk] === unitaryStep[gChk, 0.01][ψChk]),
 "‖evolveStep(V=3x²)ψ‖ error" -> normError[evolveStep[gChk, 3 #^2 &, 0.02][ψChk]]}
```

The difference is zero and the two are identical expressions: the free particle is the $V=0$ special case, recovered without approximation. The third entry checks the general claim, that `evolveStep` is unitary for any potential: the potential factor is a pure phase, so the norm is preserved to machine precision even on a stiff $V$. So every result the free-particle document established survives verbatim when we pass `0 &` for the potential, and the generalization is safe to build on.

## Part II: The Free Particle Survives the Generalization

Rerun the headline free-particle results through the generalized code, passing the zero potential. A broad packet under measurement localizes to the steady conditional width $(\hbar/4\lambda m)^{1/4}$, and the split-operator holds its norm at one to machine precision:

```wl
gF = grid[256, 50.]; λR = 1.; dt = 0.005; nt = 600;
statesF = BlockRandom[SeedRandom[11]; First@exptV[gF, 0 &, λR, dt, nt][gaussian[gF, 0., 0., 2.]]];
widthsF = Sqrt[xVar[gF, #] & /@ statesF];
AssociationThread[{"\!\(\*SqrtBox[\"Vx\"]\): start -> end", "predicted steady", "final ‖ψ‖ error"},
  {{First@widthsF, Last@widthsF}, (hbar/(4 λR mass))^(1/4), normError[Last@statesF]}]
```

The conditional width relaxes from $2$ to $\approx 0.71$, the fourth root $(\hbar/4\lambda m)^{1/4}$. At unit detection efficiency the conditional state is pure, which is what lets us carry it as a single state vector; the split-operator's substantive guarantee, which a general integrator does not give, is that the norm stays one, here to $\sim 10^{-16}$. The state mixes only in the ensemble average of Part IV. Backaction heating is also unchanged: on the position grid the unconditional (Lindblad) state integrates by the same operator split, the double-commutator dissipator damping each density-matrix element $\rho_{ij}$ by $e^{-\frac{\lambda}{2}(x_i-x_j)^2 dt}$ and the Hamiltonian acting through `evolveStep` assembled as a matrix. Build the pieces once, as reusable helpers, since the trap and the well will need them too:

```wl
decohereMat[g_, λ_, dt_] := Exp[-(λ/2) Outer[(#1 - #2)^2 &, g["x"], g["x"]] dt];
uniMat[g_, Vfun_, dt_] := Transpose[evolveStep[g, Vfun, dt] /@ IdentityMatrix[g["n"]]];
lindblad[g_, Vfun_, λ_, dt_] := With[{u = uniMat[g, Vfun, dt], dec = decohereMat[g, λ, dt]},
   Function[ρ, u . (ρ dec) . ConjugateTranspose[u]]];
p2Op[g_] := Transpose[InverseFourier[g["p"]^2 Fourier[#, FourierParameters -> {0, -1}],
      FourierParameters -> {0, -1}] & /@ IdentityMatrix[g["n"]]];
```

Propagate the free master equation from a pure packet and read the momentum variance, which the moments predict to grow linearly as $\langle\hat p^2\rangle(t) = \langle\hat p^2\rangle_0 + \lambda\hbar^2 t$:

```wl
gS = grid[32, 16.]; dtS = 0.0005; tf = 0.4; ntS = Round[tf/dtS]; ψ0 = gaussian[gS, 0., 0., 1.];
ρfree = Nest[lindblad[gS, 0 &, λR, dtS], projector[ψ0], ntS];
{"⟨p²⟩(0)" -> Re[Conjugate[ψ0] . p2Op[gS] . ψ0], "⟨p²⟩(tf) master" -> Re@Tr[ρfree . p2Op[gS]],
 "⟨p²⟩(0) + λℏ²tf" -> Re[Conjugate[ψ0] . p2Op[gS] . ψ0] + λR hbar^2 tf}
```

The momentum variance rises from $0.25$ to $0.65 = 0.25 + \lambda\hbar^2 t_f$: continuous position measurement heats the momentum at the constant rate $\lambda\hbar^2$, independent of the state and the Hamiltonian. This rate is the one piece of physics common to every case below.

The generalized code also accepts any initial wavefunction, not just a Gaussian. Take a two-packet superposition, a cat state delocalized over $x=\pm 6$, and watch continuous position measurement collapse it:

```wl
gC = grid[256, 50.]; ψCat = Normalize[gaussian[gC, -6., 0., 1.] + gaussian[gC, 6., 0., 1.]];
statesC = BlockRandom[SeedRandom[5]; First@exptV[gC, 0 &, 1., 0.005, 600][ψCat]];
{"initial Vx" -> xVar[gC, ψCat], "final Vx" -> xVar[gC, Last@statesC],
 "final ‖ψ‖ error" -> normError[Last@statesC]}
```

The initial variance is $\approx 37$, dominated by the $\pm 6$ separation of the two lobes; the final variance is $\approx 0.5$, the single-packet steady value, with the norm preserved. Continuous position measurement resolves which lobe the particle is in and localizes it there. Displaying the conditional density $|\psi(x,t)|^2$ against time makes the collapse visible: the two ridges of the initial cat compete, and one wins.

```wl
band = Flatten@Position[UnitStep[10 - Abs[gC["x"]]], 1];
ArrayPlot[Reverse[(Abs[statesC]^2)[[All, band]]],
 DataRange -> {gC["x"][[{First@band, Last@band}]], {0, 600 dt}}, ColorFunction -> "SunsetColors",
 AspectRatio -> 1.3, ImageSize -> 300, Frame -> True, FrameLabel -> {"x", "t"},
 PlotLabel -> "A cat state collapses to one packet"]
```

## Part III: A Trap That Closes, the Harmonic Oscillator

The harmonic oscillator $V(x) = \tfrac12 m\omega^2 x^2$ is the exact setting of cavity and levitated optomechanics, where the position of a mechanical mode is monitored by scattered or cavity light. Because $\hat H$ is quadratic, a Gaussian conditional state stays Gaussian under both the measurement and the evolution, so its five moments close among themselves. The conditional covariances obey a deterministic matrix Riccati flow, the free-particle equations of the companion document plus the trap terms $-m\omega^2 V_x$ and $-2m\omega^2 C$:

```wl
ClearAll[Vx, Cc, Vp, λ, ω, m, ℏ];
riccati = {2 Cc/m - 4 λ Vx^2, Vp/m - m ω^2 Vx - 4 λ Vx Cc, -2 m ω^2 Cc + λ ℏ^2 - 4 λ Cc^2};
steady = Solve[Thread[riccati == 0], {Vx, Cc, Vp}];
posRoot[rule_] := With[{v = {Vx, Cc, Vp} /. rule /. {λ -> 1., ω -> 1., m -> 1., ℏ -> 1.}},
   Max@Abs@Im[v] < 10^-10 && Re[v[[1]]] > 0 && Re[v[[3]]] > 0];
{VxSS, CcSS, VpSS} = Simplify[{Vx, Cc, Vp} /. First@Select[steady, posRoot], {λ > 0, ω > 0, m > 0, ℏ > 0}];
VxSS
```

The steady conditional position variance is $V_x^{ss} = \dfrac{1}{2\sqrt2}\sqrt{\dfrac{-m\omega^2+\sqrt{4\hbar^2\lambda^2+m^2\omega^4}}{\lambda^2 m}}$. Two limits check it. The steady conditional state is pure, so its covariances saturate the uncertainty product, and the free particle is $\omega\to 0$:

```wl
{"purity Vx Vp - C² (expect ℏ²/4)" -> Simplify[VxSS VpSS - CcSS^2, {λ > 0, ω > 0, m > 0, ℏ > 0}],
 "ω→0 (expect free Sqrt[ℏ/4λm])" -> Simplify[Limit[VxSS, ω -> 0], {λ > 0, m > 0, ℏ > 0}]}
```

The conditional purity invariant $V_x V_p - C^2$ equals $\hbar^2/4$ exactly, so the steady conditional state is a pure Gaussian, and as $\omega\to 0$ the steady variance returns to the free-particle $\sqrt{\hbar/4\lambda m}$. Now compare it to the ground-state zero-point variance $V_x^{zp} = \hbar/2m\omega$. In the single dimensionless measurement strength $\mu = 2\hbar\lambda/m\omega^2$ the ratio collapses to a clean form, below $1$ for every $\mu>0$:

```wl
μ = 2 ℏ λ/(m ω^2);
{"Vx_ss / Vx_zp" -> Simplify[VxSS/(ℏ/(2 m ω)), {λ > 0, ω > 0, m > 0, ℏ > 0}],
 "closed form Sqrt[2(Sqrt[1+μ²]-1)]/μ" -> Simplify[VxSS/(ℏ/(2 m ω)) == Sqrt[2 (Sqrt[1 + μ^2] - 1)]/μ,
    {λ > 0, ω > 0, m > 0, ℏ > 0}]}
```

The ratio matches the closed form $\sqrt{2(\sqrt{1+\mu^2}-1)}/\mu$, and its two regimes follow by expanding it:

```wl
Block[{q}, {"weak, μ→0: ratio →" -> Simplify[Normal@Series[Sqrt[2 (Sqrt[1 + q^2] - 1)]/q, {q, 0, 2}], q > 0],
   "strong, μ→∞: ratio·√μ →" -> Limit[Sqrt[q] Sqrt[2 (Sqrt[1 + q^2] - 1)]/q, q -> Infinity]}]
```

So the conditional variance dips just below the zero-point value for weak measurement, $V_x^{ss}/V_x^{zp}\approx 1-\mu^2/8$, and falls off as $\sqrt{2/\mu}$ for strong. Continuous position measurement squeezes the conditional position variance below the ground-state zero-point value at *any* nonzero strength, deepening monotonically with $\lambda$: the measurement extracts position information continuously, and the conditional state is a pure squeezed state, narrower in position than the vacuum. Fix numbers in a strong regime and read the squeezing off the closed form:

```wl
ωR = 1.; λR = 2.;
VxZP = hbar/(2 mass ωR); VpZP = hbar mass ωR/2;
harmNum = {ℏ -> hbar, m -> mass, ω -> ωR, λ -> λR};
{"Vx_zp" -> VxZP, "Vx_ss" -> (VxSS /. harmNum), "Vx_ss/Vx_zp" -> (VxSS /. harmNum)/VxZP,
 "Vp_ss/Vp_zp (anti-squeezed)" -> (VpSS /. harmNum)/VpZP, "μ" -> (μ /. harmNum)}
```

At $\mu=4$ the conditional position is squeezed to $0.625$ of the zero-point variance, and the momentum is anti-squeezed to $2.58$ times its zero-point value, the backaction paid for the position information. The grid split-operator scheme reproduces this without assuming Gaussianity. Evolve one trajectory from the ground state under measurement in the trap and track the conditional variance:

```wl
Vharm = 0.5 mass ωR^2 #^2 &;
gH = grid[256, 24.]; dtH = 0.002; tfH = 8.; ntH = Round[tfH/dtH];
runH = BlockRandom[SeedRandom[11]; trajectory[gH, Vharm, λR, dtH, ntH][gaussian[gH, 0., 0., Sqrt[VxZP]]]];
vxH = xVar[gH, #] & /@ runH;
{"grid Vx (late-time mean)" -> Mean[vxH[[Round[0.6 ntH] ;;]]], "Riccati Vx_ss" -> (VxSS /. harmNum),
 "Vx_zp" -> VxZP, "final ‖ψ‖ error" -> normError[Last@runH]}
```

The grid conditional variance relaxes to $\approx 0.313$, matching the Riccati steady value to grid resolution and sitting well below the zero-point $0.5$, with the norm held to machine precision. The relaxation is the same in every run, because the conditional variance obeys a deterministic (noiseless) equation while only the mean is driven by the record:

```wl
ListLinePlot[{Transpose[{dtH Range[0, ntH], vxH}]},
 Frame -> True, ImageSize -> 500, AspectRatio -> 1/2, PlotRange -> {0, 0.55},
 GridLines -> {None, {{VxZP, Directive[Gray, Dashed]}, {VxSS /. harmNum, Directive[Red, Dashed]}}},
 FrameLabel -> {"time t", "conditional variance Vx"},
 PlotLabel -> "Measurement squeezes Vx below the zero-point line (gray) to the Riccati value (red)"]
```

Because the moments close, the same conditional dynamics is captured by five real numbers instead of a grid: the means $\langle\hat x\rangle,\langle\hat p\rangle$, driven by the record, and the covariances $V_x,C,V_p$, following the deterministic Riccati flow. This is the continuous quantum Kalman-Bucy filter, the identical algorithm the optomechanics experiments run on a field-programmable gate array to reconstruct the mechanical state in real time. Hand it to the built-in stochastic integrator with the harmonic drifts in place:

```wl
ClearAll[mx, mp, vx, vc, vp];
momentProc = ItoProcess[
   {{mp[t]/mass, -mass ωR^2 mx[t], 2 vc[t]/mass - 4 λR vx[t]^2,
     vp[t]/mass - mass ωR^2 vx[t] - 4 λR vx[t] vc[t], -2 mass ωR^2 vc[t] + λR hbar^2 - 4 λR vc[t]^2},
    {{2 Sqrt[λR] vx[t]}, {2 Sqrt[λR] vc[t]}, {0}, {0}, {0}},
    {mx[t], mp[t], vx[t], vc[t], vp[t]}},
   {{mx, mp, vx, vc, vp}, {0., 0., VxZP, 0., VpZP}}, {t, 0}];
tdM = BlockRandom[SeedRandom[7]; RandomFunction[momentProc, {0., tfH, 0.01}, 300,
    Method -> "StochasticRungeKuttaScalarNoise"]];
finM = tdM["ValueList"][[All, -1]];
{"moment-filter Vx(tf) mean" -> Mean[#[[3]] & /@ finM], "Vx spread across runs" -> StandardDeviation[#[[3]] & /@ finM],
 "⟨x⟩ spread across runs" -> StandardDeviation[#[[1]] & /@ finM],
 "grid Vx (late)" -> Mean[vxH[[Round[0.6 ntH] ;;]]], "Riccati Vx_ss" -> (VxSS /. harmNum)}
```

The moment filter returns the same steady conditional variance as the grid, and its spread across realizations is exactly zero: the covariance is deterministic. The conditional mean, by contrast, spreads by several units across runs, so the filter's stochasticity lives entirely in the means, which the record drives through the gains $2\sqrt{\lambda}V_x$ and $2\sqrt{\lambda}C$. The five-number filter and the 256-point grid describe the same conditional state, because for a quadratic Hamiltonian the Gaussian description is exact. This is the experimental correspondence made quantitative: the conditional position uncertainty a levitated nanoparticle experiment reconstructs by Kalman filtering is the same $\sqrt{V_x}$ computed here, its scale set through $\mu$ by the ratio of the measurement rate to the mechanical frequency. At the unit efficiency assumed here the steady value sits below the zero-point; the above-zero-point figure a real experiment reports (Magrini et al., $1.3$ times zero-point) reflects its finite detection efficiency $\eta<1$, which mixes the conditional state and lifts its uncertainty, and is the extension deferred to the close.

## Part IV: Backaction in the Trap Heats Without Bound

The conditional state reaches a steady squeezed pure state. The unconditional state, averaged over the record, does not reach any steady state at all. The trap has no way to remove the energy the measurement injects, so the closed oscillator heats indefinitely. The unconditional second moments close on their own, and their energy has a constant source:

```wl
ClearAll[XX, SS, PP];
usol = First@DSolve[{XX'[t] == 2 SS[t]/m, SS'[t] == PP[t]/m - m ω^2 XX[t],
    PP'[t] == -2 m ω^2 SS[t] + λ ℏ^2, XX[0] == x0, SS[0] == 0, PP[0] == p0}, {XX, SS, PP}, t];
energy = Simplify[PP[t]/(2 m) + (1/2) m ω^2 XX[t] /. usol];
{"E(t)" -> energy, "dE/dt" -> Simplify[D[energy, t]]}
```

The energy is exactly linear in time, $E(t) = E_0 + \lambda\hbar^2 t/2m$, with no oscillatory part: the breathing of $\langle\hat x^2\rangle$ and $\langle\hat p^2\rangle$ against each other cancels in the sum, leaving only the steady heating $dE/dt = \lambda\hbar^2/2m$. In oscillator quanta this is a mean occupation rising at $dn/dt = \lambda\hbar/2m\omega$, a thermal-like mixed Gaussian of ever-increasing temperature. The grid master equation confirms it and shows the individual moments growing:

```wl
gHL = grid[96, 24.]; dtHL = 0.002;
stepHL = lindblad[gHL, Vharm, λR, dtHL];
ρ0HL = projector[gaussian[gHL, 0., 0., Sqrt[VxZP]]];
x2Mean[g_, ρ_] := Re[Diagonal[ρ] . g["x"]^2];
Emeas[ρ_] := Re@Tr[ρ . p2Op[gHL]]/(2 mass) + 0.5 mass ωR^2 x2Mean[gHL, ρ];
tsHL = {0., 0.5, 1., 1.5, 2.};
ρsHL = FoldList[Nest[stepHL, #1, Round[(#2[[2]] - #2[[1]])/dtHL]] &, ρ0HL, Partition[tsHL, 2, 1]];
TableForm[{Emeas /@ ρsHL, x2Mean[gHL, #] & /@ ρsHL, Re@Tr[#.p2Op[gHL]] & /@ ρsHL,
   Re@Tr[#.#] & /@ ρsHL},
 TableHeadings -> {{"E", "⟨x²⟩", "⟨p²⟩", "purity"}, tsHL}]
```

The energy climbs in equal steps, $\{0.5, 1, 1.5, 2, 2.5\}$ over $t\in[0,2]$, a slope of exactly $\lambda\hbar^2/2m = 1$; $\langle\hat x^2\rangle$ and $\langle\hat p^2\rangle$ both grow, and the purity falls from $1$ as the unconditional state mixes. There is no steady state to plot. The contrast with the conditional state is starker in the trap than for the free particle: one realization sits in a fixed pure squeezed packet whose variance never moves, while the average over realizations heats without bound.

```wl
ListLinePlot[Transpose[{tsHL, Emeas /@ ρsHL}],
 Frame -> True, ImageSize -> 460, AspectRatio -> 1/2, PlotRange -> {0, All}, Mesh -> All,
 FrameLabel -> {"time t", "unconditional energy E(t)"},
 PlotLabel -> "The trapped unconditional state heats linearly: dE/dt = λℏ²/2m"]
```

A genuine thermal *steady* state requires an added mechanical damping channel to balance the heating; that is exactly why the optomechanics experiments cool with feedback rather than watching the bare measurement, whose only effect on the unconditional state is to heat it. Restoring damping is a one-line addition to the master equation and is the natural next extension.

## Part V: A Trap That Does Not Close, the Double Well

Introduce a quartic potential and the Gaussian description fails. For a double well $V(x) = b(x^2-a^2)^2$ the force is cubic, so the equation for the mean momentum needs the third moment, which for a non-Gaussian state is not fixed by the mean and the variance:

```wl
ClearAll[x, a, b];
{"V'(x)" -> D[b (x^2 - a^2)^2, x], "so d⟨p⟩/dt = -⟨V'⟩ = -4b(⟨x³⟩ - a²⟨x⟩) needs ⟨x³⟩" -> True,
 "harmonic contrast: V'=" -> D[(1/2) m ω^2 x^2, x]}
```

Where the harmonic force $m\omega^2 x$ closes the hierarchy at the first moment, the cubic force $4b(x^3-a^2 x)$ pulls in $\langle\hat x^3\rangle$, that pulls in $\langle\hat x^4\rangle$, and no finite set of moments closes. Only the grid split-operator scheme, which never assumes a shape for the conditional state, applies. Build the well $V=b(x^2-a^2)^2$, with two minima at $x=\pm a$ separated by a barrier, and a helper that diagonalizes it on any grid to return the symmetric and antisymmetric doublet and the one-well combination:

```wl
aa = 1.2; bb = 1.0; Vdw = bb (#^2 - aa^2)^2 &;
leftMask[g_] := Boole[Negative[g["x"]]];
wellDoublet[g_] := Block[{ham, ev, vec, ord, sym, anti},
   ham = p2Op[g]/(2 mass) + DiagonalMatrix[Vdw[g["x"]]];
   {ev, vec} = Eigensystem[ham]; ord = Ordering[Re@ev];
   sym = Normalize[vec[[ord[[1]]]]]; anti = Normalize[vec[[ord[[2]]]]];
   <|"H" -> ham, "E" -> Re@ev[[ord[[1 ;; 2]]]], "sym" -> sym, "anti" -> anti,
     "left" -> With[{plus = Normalize[sym + anti], minus = Normalize[sym - anti]},
        If[leftMask[g] . Abs[plus]^2 >= 0.5, plus, minus]]|>];
gD = grid[256, 14.]; dw = wellDoublet[gD];
ψSym = dw["sym"]; ψAnti = dw["anti"]; Ttun = 2 Pi hbar/(dw["E"][[2]] - dw["E"][[1]]);
leftP[ψ_] := leftMask[gD] . Abs[ψ]^2;
{"E0, E1" -> dw["E"], "barrier" -> bb aa^4, "tunneling period Ttun" -> Ttun,
 "Hermiticity |H - H†|" -> Max@Abs[dw["H"] - ConjugateTranspose[dw["H"]]]}
```

The two lowest states form a doublet split by $\Delta = E_1-E_0 \approx 0.40$, below the barrier $\approx 2.07$; their symmetric and antisymmetric combinations are the states localized in the left and right wells, and a particle prepared in one well tunnels to the other with period $T_{\rm tun} = 2\pi\hbar/\Delta \approx 15.8$. The Hamiltonian is Hermitian to machine precision, so the spectrum is real and the `Re` above only strips rounding. First confirm that the conditional state under measurement is genuinely non-Gaussian, which is what forecloses the moment filter. Track the skewness and excess kurtosis of $|\psi(x)|^2$ along a measured trajectory started in one well:

```wl
central[ψ_, k_] := With[{d = Abs[ψ]^2, mu = xMean[gD, ψ]}, (gD["x"] - mu)^k . d];
skew[ψ_] := central[ψ, 3]/xVar[gD, ψ]^(3/2); exkurt[ψ_] := central[ψ, 4]/xVar[gD, ψ]^2 - 3;
trajDW = BlockRandom[SeedRandom[3]; trajectory[gD, Vdw, 2., 0.002, 1500][gaussian[gD, -aa, 0., 0.45]]];
{"max |skewness|" -> Max@Abs[skew /@ trajDW[[1 ;; ;; 100]]],
 "max |excess kurtosis|" -> Max@Abs[exkurt /@ trajDW[[1 ;; ;; 100]]]}
```

Both run to order one, far from the Gaussian value of zero: a five-number filter forcing skewness and kurtosis to zero cannot track this state, and only the grid survives. Now the two physical effects of watching a particle in a double well. First, measurement-induced localization. Prepare the symmetric eigenstate, delocalized equally over both wells, and monitor it: each run collapses into one well, chosen at random, and stays pure.

```wl
locRuns = BlockRandom[SeedRandom[10]; Table[Last@trajectory[gD, Vdw, 8., 0.002, 1500][ψSym], {8}]];
{"start P_left (symmetric)" -> leftP[ψSym], "P_left of 8 runs" -> (leftP /@ locRuns),
 "max ‖ψ‖ error over runs" -> Max[normError /@ locRuns]}
```

The symmetric state starts balanced at $P_{\rm left}=0.5$; after measurement each run has collapsed toward $0$ or $1$, localized in the right or left well, and every conditional state stays a normalized pure state (only the ensemble average, as in Part IV, is mixed). Different records send the particle to different wells, and the ensemble of outcomes is the position measurement resolving a delocalized state into a definite, if random, location. Two runs side by side show the density collapsing to opposite wells:

```wl
bandD = Flatten@Position[UnitStep[3.5 - Abs[gD["x"]]], 1];
showRun[seed_] := ArrayPlot[Reverse[(Abs[#]^2 & /@ BlockRandom[SeedRandom[seed];
       trajectory[gD, Vdw, 8., 0.002, 1500][ψSym]])[[All, bandD]]],
   DataRange -> {gD["x"][[{First@bandD, Last@bandD}]], {0, 3}}, ColorFunction -> "SunsetColors",
   AspectRatio -> 1.3, ImageSize -> 230, Frame -> True, FrameLabel -> {"x", "t"}];
GraphicsRow[{showRun[4], showRun[1]}, ImageSize -> 500]
```

Second, the quantum Zeno effect. Watching the particle suppresses the coherent tunneling between wells. The cleanest way to see it is the unconditional (Lindblad) population, which is deterministic: the measurement dephases the left-right coherence at rate $2\lambda a^2$, converting coherent oscillation into slow incoherent equilibration. On a coarser grid the dense master equation can afford, take the same one-well state from the doublet helper and propagate it at a few measurement strengths, reading the left-well population over one tunneling period:

```wl
gDL = grid[128, 12.]; dtDL = 0.005;
ψLwell = wellDoublet[gDL]["left"];
leftPL[ρ_] := Re[Diagonal[ρ] . leftMask[gDL]];
EmeasL[ρ_] := Re@Tr[ρ . p2Op[gDL]]/(2 mass) + Re[Diagonal[ρ] . Vdw[gDL["x"]]];
smp = Round[Ttun/dtDL/8];
zenoStates[λ_] := NestList[Nest[lindblad[gDL, Vdw, λ, dtDL], #, smp] &, projector[ψLwell], 8];
zTimes = Range[0, 8] smp dtDL; zStates = zenoStates /@ {0., 0.5, 4.};
TableForm[Transpose[Prepend[Map[leftPL, zStates, {2}], zTimes]],
 TableHeadings -> {None, {"time t", "P_left, λ=0", "λ=0.5", "λ=4"}}]
```

Without measurement the left-well population oscillates coherently, $0.96\to 0.03$ at the half period, where the particle has fully tunneled to the other well, and back. At $\lambda=0.5$ that coherent oscillation is gone: the measurement dephases the left-right superposition, so the population no longer swings across and returns but relaxes monotonically toward equipartition at $0.5$. The coherent tunneling transfer is suppressed. This is the quantum Zeno signature of a continuous position measurement: watching the particle destroys the phase coherence the tunneling needs.

The backaction is the other half of the story. The same $\lambda\hbar^2$ heating that has no steady state in the harmonic trap is present here too, and over one tunneling period it carries the particle over the barrier at every strength shown:

```wl
{"barrier height" -> bb aa^4, "energy at one tunneling period, λ = 0, 0.5, 4" -> (EmeasL[Last@#] & /@ zStates)}
```

By one period the energy has already passed the barrier ($\approx 2.07$) at $\lambda=0.5$ ($\approx 5.4$), and reaches $\approx 33$ at $\lambda=4$; the stronger the measurement, the faster the population equilibrates over the top rather than through the barrier. So watching does two things at once, both growing with $\lambda$: it dephases the wells, killing the coherent oscillation, and it heats the particle, eventually over the barrier.

```wl
ListLinePlot[Transpose[{zTimes, leftPL /@ #}] & /@ zStates,
 PlotLegends -> {"λ=0 (unwatched)", "λ=0.5 (dephased)", "λ=4 (heated over barrier)"}, Frame -> True, ImageSize -> 500,
 AspectRatio -> 1/2, PlotRange -> {0, 1}, PlotMarkers -> Automatic,
 FrameLabel -> {"time t", "left-well population"},
 PlotLabel -> "Watching kills the coherent oscillation (λ=0.5) and, harder still, heats over the barrier (λ=4)"]
```

Holding the particle in one well would take the dephasing without the heating, and unit-efficiency position measurement cannot separate them: the one coupling delivers both. This is the honest reach of the Zeno effect under continuous position measurement, and the reason a real experiment steers with feedback rather than relying on the measurement alone.

## Part VI: The Division of Labor

The generalization draws a clean line between the two ways to simulate a monitored quantum system. When the conditional state lives in a small closed parametrization, the built-in stochastic integrator suffices and the geometry keeps the state physical on its own: a qubit is three Bloch numbers confined to the ball, and a linear (Gaussian) system, the free particle or the harmonic trap, is five moments whose covariance flow has an attracting pure-state fixed point. Both close because the Hamiltonian is at most quadratic, and for both the Kalman-Bucy moment filter is exact and is literally the algorithm the optomechanics experiments run.

Otherwise, hand-write the step whose structure you can exploit. A potential that is anything beyond quadratic breaks the moment closure, the conditional state is no longer Gaussian, and only the grid split-operator scheme applies: the measurement as an exact Kraus operator with its renormalization, independent of $\hat H$, and the Hamiltonian as a Strang split whose kinetic factor is spectral and whose potential factor is a pointwise phase. The double well is the smallest example that forces this choice, and it rewards it with physics the moment filter cannot represent, the non-Gaussian conditional state, the collapse into one well, and the Zeno suppression of coherent tunneling set against the backaction heating that drives the particle over the barrier.

### Where This Leaves Us

The measurement is the same in every case; only the unitary changes, and one Strang step carries any one-dimensional potential. The free particle is recovered exactly at $V=0$. In the harmonic trap the conditional state is a pure squeezed Gaussian, its position driven below the zero-point value by the measurement, while the unconditional state heats without bound at $dE/dt = \lambda\hbar^2/2m$, the same backaction rate throughout. In the double well the conditional state goes non-Gaussian, localizes into one well, and has its coherent tunneling suppressed by observation, while the same backaction heats it over the barrier. The natural extensions are to reduce the detection efficiency below unity, so a single record no longer purifies the conditional state, or to add mechanical damping, so the trapped unconditional state reaches the thermal steady state the bare measurement never does.
