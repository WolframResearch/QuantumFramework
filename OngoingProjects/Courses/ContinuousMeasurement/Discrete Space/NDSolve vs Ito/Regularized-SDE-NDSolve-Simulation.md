---
Template: Default
---

# Handing NDSolve a regularized Stratonovich equation: the state and readout of a monitored qubit

**This document shows how to simulate a continuously monitored qubit with `NDSolve`. We have implemented this approach for a model that follows four quantities: three numbers that describe the qubit state (i.e., the components of the Bloch vector) and the detector's accumulated readout. The dynamical equations for these four quantities are usually written as a nonlinear stochastic differential equation in Itô form. The obstacle is that `NDSolve` cannot use noise as an ordinary function of time—or, more generally, any variable that must be sampled from a distribution. Following the essence of the Wong–Zakai theorem, we therefore convert the original Itô equation to its equivalent Stratonovich form, which is suited to smooth approximations of the noise, and replace the rough Wiener path with a continuous piecewise-cubic path. `NDSolve` can then solve the resulting random ordinary differential equation. As the noise grid is refined, this regularized equation approaches the original stochastic model under the conditions of the Wong–Zakai theorem. Symbolic checks recover the standard master-equation evolution exactly, while finite-sample simulations compare the mean and spread with independent references. These numerical comparisons support the method but do not by themselves prove convergence.**

Mads Bahrami (last updated: August 12, 2026)

## The trajectory equations: the state and the readout in one SDE

A continuously monitored quantum state can be understood schematically as $d\rho_c=\mathcal L(\rho_c)\,dt+\mathcal M(\rho_c)\,dW$: the drift $\mathcal L(\rho_c)\,dt$ describes deterministic Hamiltonian evolution, relaxation, dephasing, and the unconditional disturbance caused by measurement, while $\mathcal M(\rho_c)\,dW$ is the stochastic, record-dependent measurement update. The detector record obeys $dQ=s(\rho_c)\,dt+dW$, where $s(\rho_c)$ is the expected signal, so $dW=dQ-s(\rho_c)\,dt$ is the innovation—the unpredictable part of the observed result—and the same Wiener increment updates both the record and the state. Because the stochastic average $\mathbb E[dW]=0$ and $dW^2=dt$, $dW$ is of order $\sqrt{dt}$; thus $d\rho_c$ is not an ordinary time derivative but a deterministic prediction plus a random correction conditioned on what the detector reports. Averaging over all possible records removes the stochastic term and recovers the master equation.

Continuous homodyne monitoring of a driven qubit produces two coupled equations driven by the same noise: one for the observer's running estimate of the Bloch vector $\{x,y,z\}=\{\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle\}$, in the convention $\hat\sigma_z=|g\rangle\langle g|-|e\rangle\langle e|$ that puts the ground state at $z=+1$ (some literature uses the opposite sign), and one for the readout $Q$ the detector accumulates while measuring $\hat\sigma_z$. Stacking the readout as a fourth component, the pair is a single scalar-noise Itô stochastic differential equation, a drift $\vec a$ times $dt$ plus a diffusion $\vec b$ times one Wiener increment $dW$:

$$
d\begin{pmatrix} x \\ y \\ z \\ Q \end{pmatrix}
= \underbrace{\begin{pmatrix} -\Gamma_2^{\mathrm{eff}}\,x \\ -\Gamma_2^{\mathrm{eff}}\,y-\Omega_x\,z \\ \gamma_1(1-z)+\Omega_x\,y \\ \sqrt{\Gamma_{CI}}\,z \end{pmatrix}}_{\displaystyle \vec a}\,dt
\;+\; \underbrace{\begin{pmatrix} -\sqrt{\Gamma_{CI}}\,x z-\sqrt{\Gamma_{BA}}\,y \\ -\sqrt{\Gamma_{CI}}\,y z+\sqrt{\Gamma_{BA}}\,x \\ \sqrt{\Gamma_{CI}}\,(1-z^2) \\ 1 \end{pmatrix}}_{\displaystyle \vec b}\,dW .
$$

This is the target: the process whose mean we will check against the Lindblad master equation and whose spread we will check against a native integrator. The top three rows, the Bloch state, close on their own, nothing in them depends on $Q$; the fourth row is the record, an integral whose drift $\sqrt{\Gamma_{CI}}\,z$ carries the $\hat\sigma_z$ signal and whose noise is the same $dW$ that kicks the state.

The six parameters have units of inverse time: $\Gamma_{CI}$ is the informative (coherent-information) rate, which sharpens the conditional state toward a $\hat\sigma_z$ eigenstate; $\Gamma_{BA}$ is the backaction rate of the orthogonal homodyne quadrature, which rotates the transverse phase; $\Gamma_d$ is the measurement-induced ensemble dephasing; $\gamma_\phi$ is intrinsic pure dephasing; $\gamma_1=1/T_1$ is relaxation; and $\Omega_x$ is the Rabi frequency. One combination recurs, the effective transverse decay $\Gamma_2^{\mathrm{eff}}=\gamma_1/2+\gamma_\phi+\Gamma_d$. When $\Gamma_d>0$, the detector efficiency is $\eta=(\Gamma_{CI}+\Gamma_{BA})/(2\Gamma_d)$, with a quantum-limited detector satisfying $0\le\eta\le1$. The benchmark ensembles below use $\eta=1$; one deliberately unphysical example violates the admissibility condition. Numerically, time and all rates are nondimensionalized by a chosen reference rate (the representative example takes that reference to be $\gamma_1$). Throughout, the code takes a rate vector $r=\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$.

For a single scalar noise, the Itô equation above is identical to a Stratonovich equation with the same diffusion and a drift lowered by one exact term (the Itô-Stratonovich conversion).

$$
d\begin{pmatrix} x \\ y \\ z \\ Q \end{pmatrix}
=
\begin{pmatrix}
\big(\tfrac{\Gamma_{BA}+\Gamma_{CI}}{2}-\Gamma_2^{\mathrm{eff}}-\Gamma_{CI}\,z^2\big)\,x-\sqrt{\Gamma_{BA}\Gamma_{CI}}\,y z \\
\big(\tfrac{\Gamma_{BA}+\Gamma_{CI}}{2}-\Gamma_2^{\mathrm{eff}}-\Gamma_{CI}\,z^2\big)\,y+\sqrt{\Gamma_{BA}\Gamma_{CI}}\,x z-\Omega_x\,z \\
\gamma_1(1-z)+\Gamma_{CI}\,z(1-z^2)+\Omega_x\,y \\
\sqrt{\Gamma_{CI}}\,z
\end{pmatrix} dt
\;+\;
\begin{pmatrix}
-\sqrt{\Gamma_{CI}}\,x z-\sqrt{\Gamma_{BA}}\,y \\
-\sqrt{\Gamma_{CI}}\,y z+\sqrt{\Gamma_{BA}}\,x \\
\sqrt{\Gamma_{CI}}\,(1-z^2) \\
1
\end{pmatrix} \circ dW .
$$

We feed the Stratonovich form to the random ODE because continuous, piecewise-smooth approximations of a one-dimensional Wiener path converge under the usual [Wong-Zakai hypotheses](https://arxiv.org/abs/1403.7281) to the Stratonovich solution. The fact that the driving noise is scalar matters: with several noncommuting noise channels, the limiting equation can depend on how the missing iterated integrals are approximated. At any nonzero grid spacing the ODE is only a regularized model; the equality with the target SDE is a limit statement, not an identity at finite resolution. The readout row is the same in both calculi because its diffusion coefficient is the constant one, so its conversion correction vanishes. After nondimensionalizing time, $Q$ is dimensionless; before nondimensionalization it carries the same units as $W$, namely $\sqrt{\text{time}}$.

## Encoding the Stratonovich vector field

Encode the Stratonovich drift, all four rows and name the entries of the rate vector $r=\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$:

```wl
ClearAll[driftStrat];
driftStrat[{x_, y_, z_, q_}, r_] :=
  With[{\[CapitalGamma]CI = r[[1]], \[CapitalGamma]BA = r[[2]], \[CapitalGamma]d = r[[3]],
    \[Gamma]\[Phi] = r[[4]], \[Gamma]1 = r[[5]], \[CapitalOmega]x = r[[6]]},
   With[{\[CapitalGamma]2 = \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d},
    {((\[CapitalGamma]BA + \[CapitalGamma]CI)/2 - \[CapitalGamma]2 - \[CapitalGamma]CI z^2) x - Sqrt[\[CapitalGamma]BA \[CapitalGamma]CI] y z,
     ((\[CapitalGamma]BA + \[CapitalGamma]CI)/2 - \[CapitalGamma]2 - \[CapitalGamma]CI z^2) y + Sqrt[\[CapitalGamma]BA \[CapitalGamma]CI] x z - \[CapitalOmega]x z,
     \[Gamma]1 (1 - z) + \[CapitalGamma]CI z (1 - z^2) + \[CapitalOmega]x y,
     Sqrt[\[CapitalGamma]CI] z}]];
```

Encode the Stratonovich diffusion, the same $\vec b$ as the Itô form:

```wl
ClearAll[diffStrat];
diffStrat[{x_, y_, z_, q_}, r_] :=
  With[{\[CapitalGamma]CI = r[[1]], \[CapitalGamma]BA = r[[2]]},
   {-Sqrt[\[CapitalGamma]CI] x z - Sqrt[\[CapitalGamma]BA] y,
    -Sqrt[\[CapitalGamma]CI] y z + Sqrt[\[CapitalGamma]BA] x,
    Sqrt[\[CapitalGamma]CI] (1 - z^2),
    1}];
```

To read these back as formulas, write the six rates as symbols, using string subscripts so the labels stay $\Gamma_{CI},\Gamma_{BA},\dots$ and never collide with the Bloch $x$ or fuse into a $C\cdot I$ product. The readbacks and the symbolic conversion below use $x,y,z,q,Q,t$, a formal noise $w$, the starting values, and the rate symbols themselves as unassigned variables, so clear them first:

```wl
ClearAll[x, y, z, q, Q, t, w, x0, y0, z0, q0, \[CapitalGamma], \[Gamma], \[CapitalOmega]];
ratesSym = {Subscript[\[CapitalGamma], "CI"], Subscript[\[CapitalGamma], "BA"],
   Subscript[\[CapitalGamma], "d"], Subscript[\[Gamma], "\[Phi]"],
   Subscript[\[Gamma], "1"], Subscript[\[CapitalOmega], "x"]};
```

Verify that the encoded drift is the Stratonovich drift written above, by reading it back on the symbolic rates:

```wl
driftStrat[{x, y, z, Q}, ratesSym] // MatrixForm
```

As one can see, the four rows are term for term the drift in the display equation, with Wolfram ordering each sum its own way: the state rows carry the $-\Gamma_{CI}z^2$ sharpening and the $\sqrt{\Gamma_{BA}\Gamma_{CI}}$ cross term, and the readout row is the bare $\sqrt{\Gamma_{CI}}\,z$ signal. Confirm the diffusion the same way:

```wl
diffStrat[{x, y, z, Q}, ratesSym] // MatrixForm
```

As expected, the three state rows are $\vec b$ and the readout row is the constant one, the additive noise that makes $Q$ the same in both stochastic calculi. The code and the equations are one object, so from here on we compute with the code.

Fix a representative normalized benchmark, rescaled by $\gamma_1$; its readout is quantum-limited ($\eta=1$), $\Gamma_{CI}+\Gamma_{BA}=2\Gamma_d$ exactly:

```wl
ratesB = {0.252, 0.1, 0.176, 1., 1., 3.};
```

Start from the ground state with the readout zeroed, a four-vector: the three Bloch numbers followed by the record $Q$:

```wl
initB = {0., 0., 1., 0.};
```

Fix a short time window and the trajectory count:

```wl
tfB = 3.; nB = 600;
```

For the trajectory to stay a density matrix, the measurement must fit inside the total dephasing, $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$, an admissibility condition derived below. For $\Gamma_d>0$ it reads $\eta\le1+\gamma_\phi/\Gamma_d$; a physical detector with $\eta\le1$ and nonnegative intrinsic dephasing satisfies it. Confirm this set passes:

```wl
ratesB[[1]] + ratesB[[2]] <= 2 (ratesB[[3]] + ratesB[[4]])
```

The quantum-limit claim is one division away; define the efficiency and read this set's:

```wl
ClearAll[efficiency];
efficiency[r_] := (r[[1]] + r[[2]])/(2 r[[3]]);
efficiency[ratesB]
```

The efficiency is one at the input precision, by construction.

The reference for the mean is the Lindblad master equation. Because the Itô drift $\vec a$ is affine in $\{x,y,z\}$, the ensemble mean obeys that drift as a closed linear ODE. Encode $\vec a$'s three state rows once; the symbolic conversion check below must land on this same object, so one encoding serves both:

```wl
ClearAll[driftIto];
driftIto[{x_, y_, z_}, r_] :=
  With[{\[CapitalGamma]d = r[[3]], \[Gamma]\[Phi] = r[[4]], \[Gamma]1 = r[[5]], \[CapitalOmega]x = r[[6]]},
   With[{\[CapitalGamma]2 = \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d},
    {-\[CapitalGamma]2 x, -\[CapitalGamma]2 y - \[CapitalOmega]x z, \[Gamma]1 (1 - z) + \[CapitalOmega]x y}]];
```

Calling that encoding "the Lindblad master equation" requires a derivation. Build the generator itself, the drive $H=\tfrac{\Omega_x}{2}\sigma_x$ with relaxation $\gamma_1\,\mathcal{D}[\sigma_-]$ toward the $z=+1$ ground state and total dephasing $\tfrac{\gamma_\phi+\Gamma_d}{2}\,\mathcal{D}[\sigma_z]$, where $\mathcal{D}[L]\rho=L\rho L^\dagger-\tfrac12(L^\dagger L\rho+\rho L^\dagger L)$; project it onto the Pauli components and subtract `driftIto`:

```wl
With[{sx = PauliMatrix[1], sy = PauliMatrix[2], sz = PauliMatrix[3],
  sm = {{0, 1}, {0, 0}}, id = IdentityMatrix[2]},
 Module[{\[Rho], diss, gen},
  \[Rho] = (id + x sx + y sy + z sz)/2;
  diss[L_] := L . \[Rho] . ConjugateTranspose[L] -
    (ConjugateTranspose[L] . L . \[Rho] + \[Rho] . ConjugateTranspose[L] . L)/2;
  gen = -I (ratesSym[[6]]/2) (sx . \[Rho] - \[Rho] . sx) +
    ratesSym[[5]] diss[sm] + ((ratesSym[[4]] + ratesSym[[3]])/2) diss[sz];
  Simplify[(Tr[gen . #] & /@ {sx, sy, sz}) - driftIto[{x, y, z}, ratesSym]]]]
```

Every component cancels, so the affine drift is the Lindblad master equation in Bloch coordinates. Relaxation carries $\gamma_1$, while the $\sigma_z$ dissipator carries intrinsic and measurement-induced dephasing. Solve that affine system with `DSolve`, building the equations from `driftIto` itself so the reference and the checked drift stay one object, and return $\langle\vec v\rangle(t)$:

```wl
ClearAll[lindbladBloch];
lindbladBloch[r_, {x0_, y0_, z0_}] :=
  Module[{sol, xf, yf, zf, t},
   sol = DSolve[Join[MapThread[#1'[t] == #2 &, {{xf, yf, zf}, driftIto[{xf[t], yf[t], zf[t]}, r]}],
       {xf[0] == x0, yf[0] == y0, zf[0] == z0}], {xf, yf, zf}, t] // First;
   Function[tt, Re[{xf[tt], yf[tt], zf[tt]} /. sol]]];
```

Two numerical scales enter the simulation: the spacing used to sample the Wiener path and the integration step chosen by `NDSolve`. A convenient dimensional estimate includes the drive, longitudinal relaxation, transverse decay, and measurement strength. It is only a starting mesh heuristic; multiplicative noise supplies no universal accuracy guarantee from rates alone.

```wl
ClearAll[maxRate];
maxRate[r_] :=
  With[{\[CapitalGamma]CI = r[[1]], \[CapitalGamma]BA = r[[2]], \[CapitalGamma]d = r[[3]],
    \[Gamma]\[Phi] = r[[4]], \[Gamma]1 = r[[5]], \[CapitalOmega]x = r[[6]]},
   Max[Abs[\[CapitalOmega]x], \[Gamma]1, \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d,
    \[CapitalGamma]CI + \[CapitalGamma]BA]];
```

Use forty target intervals across the shortest estimated timescale as the working mesh. This is a reproducible default, not a certified tolerance:

```wl
ClearAll[dtGenOf];
dtGenOf[r_] := 1/(40 maxRate[r]);
```

## Regularizing the noise: a Wiener path NDSolve can integrate

A white-noise increment $dW$ is not an ordinary function of time, so it cannot simply appear as a classical forcing term in `NDSolve`. Sample a Wiener path on a grid and interpolate those samples by a continuous piecewise-cubic function $W_h(t)$. Away from the knots, $W_h'(t)$ is an ordinary polynomial forcing, and

$$
d\vec v=\vec a_{\mathrm{Strat}}\,dt+\vec b\circ dW
\quad\longrightarrow\quad
\dot{\vec v}_h=\vec a_{\mathrm{Strat}}(\vec v_h)+\vec b(\vec v_h)W_h'(t).
$$

For Wolfram Language's [default order-three interpolation of value data](https://reference.wolfram.com/language/ref/Interpolation.html), $W_h$ is continuous but its first derivative need not be continuous at a knot. Thus the replacement is a piecewise-smooth random ODE, not a globally smooth forcing. The forcing changes on the mesh scale $h$; calling $h$ a literal correlation time would be too strong because the interpolated derivative is neither stationary nor independent from interval to interval.

One grid carries both the noise and the read-off. The input `dtMax` is an upper bound: `Subdivide` lands exactly on $t_f$, so the realized spacing can be slightly smaller. At least three intervals give the cubic interpolant four points:

```wl
ClearAll[noiseGrid, noiseStep];
noiseGrid[tf_, dtMax_] := Subdivide[0., tf, Max[3, Ceiling[tf/dtMax]]];
noiseStep[tf_, dtMax_] := tf/Max[3, Ceiling[tf/dtMax]];
```

On that grid, draw independent Gaussian increments with variance equal to the realized spacing, accumulate them, and interpolate the sampled path:

```wl
ClearAll[wienerFun];
wienerFun[tf_, dtMax_, seed_] := BlockRandom[SeedRandom[seed];
   With[{grid = noiseGrid[tf, dtMax]},
    With[{incr = RandomVariate[NormalDistribution[0, Sqrt[noiseStep[tf, dtMax]]], Length[grid] - 1]},
     Interpolation[Transpose[{grid, Prepend[Accumulate[incr], 0.]}], InterpolationOrder -> 3]]]];
```

Visualize one realization of the path and its piecewise-defined derivative, the classical forcing that stands in for white noise:

```wl
With[{Wf = wienerFun[tfB, dtGenOf[ratesB], 11],
  g = noiseGrid[tfB, dtGenOf[ratesB]]},
 With[{s = Subdivide[0., tfB, 2 (Length[g] - 1)]},
  Row[{
    ListLinePlot[Transpose[{s, Wf[s]}], Frame -> True, GridLines -> Automatic,
     PlotLabel -> "regularized Wiener path W(t)", FrameLabel -> {"t", "W"},
     ImageSize -> Medium], "   ",
    ListLinePlot[Transpose[{s, Wf'[s]}], Frame -> True, GridLines -> Automatic,
     PlotLabel -> "its derivative W'(t): piecewise-smooth noise", FrameLabel -> {"t", "W'"},
     ImageSize -> Medium]}]]]
```

The realized noise spacing $h$ selects the regularized random ODE and controls its bias relative to the white-noise limit. `MaxStepSize` controls only the numerical integration of that fixed ODE. Capping it below $h$ makes `NDSolve` resolve each interpolation interval, but neither the cap nor the rate heuristic establishes convergence by itself.

Assemble one trajectory. Build the piecewise-smooth forcing, form the random ODE $\dot{\vec v}=\vec a_{\mathrm{Strat}}(\vec v)+\vec b_{\mathrm{Strat}}(\vec v)W_h'(t)$ for the full four-vector, and hand it to `NDSolve`. If no cap is supplied, use one quarter of the realized noise-grid spacing:

```wl
ClearAll[oneTraj];
oneTraj[r_, init_, tf_, dtGen_, seed_, maxstep_ : Automatic, goal_ : Automatic] :=
  Module[{Wf, xx, yy, zz, QQ, t, rhs, eqs, stepCap},
   Wf = wienerFun[tf, dtGen, seed];
   stepCap = Replace[maxstep, Automatic -> noiseStep[tf, dtGen]/4];
   rhs = driftStrat[{xx[t], yy[t], zz[t], QQ[t]}, r] +
     diffStrat[{xx[t], yy[t], zz[t], QQ[t]}, r] Wf'[t];
   eqs = Join[MapThread[#1'[t] == #2 &, {{xx, yy, zz, QQ}, rhs}],
     {xx[0] == init[[1]], yy[0] == init[[2]], zz[0] == init[[3]], QQ[0] == init[[4]]}];
   NDSolveValue[eqs, {xx, yy, zz, QQ}, {t, 0, tf}, MaxStepSize -> stepCap,
    AccuracyGoal -> goal, PrecisionGoal -> goal]];
```

Show one trajectory two ways: the monitored component $z(t)$ against the exact Lindblad mean, and the readout $Q(t)$ it accumulates, both built from the single four-vector solution:

```wl
With[{traj = oneTraj[ratesB, initB, tfB, dtGenOf[ratesB], 7, dtGenOf[ratesB]/4],
   lb = lindbladBloch[ratesB, initB[[1 ;; 3]]]},
 Row[{
   Plot[{traj[[3]][t], lb[t][[3]]}, {t, 0, tfB},
    PlotLegends -> Placed[{"one trajectory z(t)", "Lindblad \[LeftAngleBracket]z\[RightAngleBracket](t)"}, Below],
    Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "z"},
    PlotLabel -> "state: the monitored component", ImageSize -> Medium], "   ",
   Plot[traj[[4]][t], {t, 0, tfB},
    Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "Q"},
    PlotLabel -> "readout: the accumulated record Q(t)", ImageSize -> Medium]}]]
```

As one can see, the state path fluctuates around the Lindblad mean, as a single record must, while the readout $Q$ climbs steeply where $z$ is near one and bends as $z$ falls, its $\sqrt{\Gamma_{CI}}\,z$ signal buried under noise. The two panels share one $dW$, and that is the physical content of a measurement record: the noise written into $Q$ is the same noise that kicked $z$.

## The ensemble: matching the mean and the spread

A single path tells us little about the law. We therefore retain full trajectories from two finite-step simulations: the built-in `StratonovichProcess` integrator and the regularized `NDSolve` route. Their means can be compared with the closed Lindblad solution, while their spreads probe information that the mean cannot see. Agreement is evidence at the stated time steps and sample size; it is not by itself a convergence theorem.

Before any simulation, let's verify the Itô-Stratonovich conversion. If the explicit Stratonovich drift truly carries the correction, its Itô form must be the affine drift $\vec a$ we encoded as `driftIto` for the Lindblad reference, with the readout row untouched. Convert the full four-component field and test it against that target, under the physical assumption that the rates are nonnegative (needed only to fold $\sqrt{\Gamma_{CI}}\,\sqrt{\Gamma_{BA}}$ into a single root):

```wl
Assuming[{Subscript[\[CapitalGamma], "CI"] >= 0, Subscript[\[CapitalGamma], "BA"] >= 0},
 FullSimplify[
  (ItoProcess[StratonovichProcess[{driftStrat[{x, y, z, q}, ratesSym],
        List /@ diffStrat[{x, y, z, q}, ratesSym], {x, y, z, q}},
       {{x, y, z, q}, {x0, y0, z0, q0}}, t]]["Drift"] /. {x[t] -> x, y[t] -> y, z[t] -> z, q[t] -> q})
   == Append[driftIto[{x, y, z}, ratesSym], Sqrt[Subscript[\[CapitalGamma], "CI"]] z]]]
```

The equality holds row by row: the Itô form of the Stratonovich SDE we hand `NDSolve` is exactly the affine drift on all three Bloch rows, with the readout row unchanged, so the correction is carried in full and $Q$ is the same in both calculi.

The trajectories are independent integrations, so the regularized generator below spreads them across local subkernels with `ParallelTable`, which launches the configured kernels on its first use and ships the `Global` context definitions to them automatically. The wall-clock gain is set by the machine's core mix and scheduling, typically a factor of a few rather than the kernel count; with no subkernels available the same cells simply evaluate sequentially.

Generate `ntraj` full trajectories from the built-in `StratonovichProcess` integrator. Feed it the same four-row Stratonovich field and request the [documented order-$3/2$ scalar-noise method](https://reference.wolfram.com/language/ref/StratonovichProcess.html) rather than the default order-$1/2$ Euler-Maruyama method. The method order describes asymptotic discretization behavior; it does not remove finite-step bias:

```wl
ClearAll[stratRun];
stratRun[r_, init_, tf_, dt_, ntraj_, seed_] :=
  Module[{t, v, proc, td},
   v = {\[FormalX][t], \[FormalY][t], \[FormalZ][t], \[FormalQ][t]};
   proc = StratonovichProcess[{driftStrat[v, r], List /@ diffStrat[v, r], v},
     {{\[FormalX], \[FormalY], \[FormalZ], \[FormalQ]}, init}, t];
   td = BlockRandom[SeedRandom[seed];
     RandomFunction[proc, {0., tf, dt}, ntraj, Method -> "StochasticRungeKuttaScalarNoise"]];
   td["ValueList"]];
```

Generate the same count of regularized trajectories, each `NDSolve` solution read off on a shared time grid, keeping all four components so the record stays beside the state that produced it. Three arguments carry defaults: the integration-step cap is one quarter of the realized noise spacing, the read-off grid is the noise grid itself, and the solver goals pass to `NDSolve` unchanged:

```wl
ClearAll[regRun];
regRun[r_, init_, tf_, dtGen_, ntraj_, baseSeed_, maxstep_ : Automatic, grid_ : Automatic, goal_ : Automatic] :=
  With[{ms = maxstep /. Automatic -> noiseStep[tf, dtGen]/4,
    g = grid /. Automatic -> noiseGrid[tf, dtGen]},
   ParallelTable[
    With[{tr = oneTraj[r, init, tf, dtGen, baseSeed + i, ms, goal]},
     Transpose[#[g] & /@ tr]],
    {i, ntraj}]];
```

For the harder regimes and the step-size study below we will want only the final values, so reduce each generator to them. The native final Bloch vectors:

```wl
ClearAll[stratRef];
stratRef[r_, init_, tf_, dt_, ntraj_, seed_] := stratRun[r, init, tf, dt, ntraj, seed][[All, -1, 1 ;; 3]];
```

And the regularized final Bloch vectors, read off at $t_f$ alone, the cap and the goals defaulting the same way:

```wl
ClearAll[ensemble];
ensemble[r_, init_, tf_, dtGen_, ntraj_, baseSeed_, maxstep_ : Automatic, goal_ : Automatic] :=
  regRun[r, init, tf, dtGen, ntraj, baseSeed, maxstep, {tf}, goal][[All, -1, 1 ;; 3]];
```

Run both on the representative case, on the shared grid set by `dtGenOf[ratesB]`. First the native ensemble:

```wl
stratRunB = stratRun[ratesB, initB, tfB, dtGenOf[ratesB], nB, 2024];
```

Then the regularized ensemble on the same grid. Its displayed statistics resolve differences only at Monte-Carlo scale, so the production ensembles pass `AccuracyGoal` and `PrecisionGoal` of five to `NDSolve`; the pathwise identity checks below, where solver error itself is on display, keep the defaults:

```wl
regRunB = regRun[ratesB, initB, tfB, dtGenOf[ratesB], nB, 9100, Automatic, Automatic, 5];
```

Average each ensemble over its trajectories and compare the two mean paths of $\langle z\rangle$ against the single exact Lindblad trajectory:

```wl
With[{g = noiseGrid[tfB, dtGenOf[ratesB]], lb = lindbladBloch[ratesB, initB[[1 ;; 3]]]},
 ListLinePlot[{
    Transpose[{g, Mean[stratRunB][[All, 3]]}],
    Transpose[{g, Mean[regRunB][[All, 3]]}],
    Transpose[{g, lb[#][[3]] & /@ g}]},
   PlotStyle -> {Automatic, Dashed, Thick},
   PlotLegends -> {"StratonovichProcess mean", "regularized mean", "Lindblad \[LeftAngleBracket]z\[RightAngleBracket](t)"},
   Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "\[LeftAngleBracket]z\[RightAngleBracket]"},
   PlotLabel -> "averaging the trajectories rebuilds the master equation"]]
```

The plot compares both sampled ensemble means with the exact Lindblad curve over the whole stored grid. In the white-noise theory the ensemble average obeys the master equation exactly; finite ensembles fluctuate around it.

Quantify the largest gap between the regularized mean of $z$ and the Lindblad value anywhere on the grid, in Monte-Carlo standard errors, skipping $t=0$ (where all paths coincide and the standard error vanishes):

```wl
With[{g = Rest[noiseGrid[tfB, dtGenOf[ratesB]]], lb = lindbladBloch[ratesB, initB[[1 ;; 3]]],
  zc = regRunB[[All, 2 ;;, 3]]},
 Max[Abs[Mean[zc] - (lb[#][[3]] & /@ g)]/(StandardDeviation[zc]/Sqrt[nB])]]
```

This is the largest pointwise standardized gap on the stored grid. Because it is a maximum over many correlated times, ordinary one-, two-, or three-standard-error thresholds do not calibrate it as a simultaneous test. Use it as a scale diagnostic; a formal whole-curve test would require the sampling distribution of the maximum, for example from an independent bootstrap.

The mean sees only the drift, so it is a weak test. The measurement backaction lives in the standard deviation of the final Bloch vector, its spread across the ensemble, compared across the two routes:

```wl
TableForm[
 Transpose[{StandardDeviation[stratRunB[[All, -1, 1 ;; 3]]],
   StandardDeviation[regRunB[[All, -1, 1 ;; 3]]]}],
 TableHeadings -> {{"Spread in x", "Spread in y", "Spread in z"}, {"Stratonovich", "Regularized"}}]
```

The displayed spreads provide a diffusion-sensitive comparison between the two finite-step routes. All three components fluctuate: $x$ starts at zero and its mean stays there, yet the backaction quadrature feeds its fluctuations from $y$.

The sample spread is itself uncertain. A large-sample delta-method estimate is $\mathrm{SE}(s)\approx\sqrt{(m_4-s^4)/(4s^2n)}$, where $m_4$ is the fourth central moment. The plug-in code below estimates the moments from the same sample. It is an asymptotic scale estimate, not an exact finite-sample error bar.

Encode $\mathrm{SE}$, the error bar of one spread, both moments read off the sample:

```wl
ClearAll[sdErr];
sdErr[x_] := With[{v = CentralMoment[x, 2]}, Sqrt[(CentralMoment[x, 4] - v^2)/(4 v Length[x])]];
```

Compare the Stratonovich and regularized spreads by dividing their absolute difference by their combined standard error,

$$
G=
\frac{\left|s_{\mathrm{Strat}}-s_{\mathrm{reg}}\right|}
{\sqrt{\mathrm{SE}(s_{\mathrm{Strat}})^2+
       \mathrm{SE}(s_{\mathrm{reg}})^2}},
$$

so $G=1$ means that the two spreads differ by one combined standard error. Encode this statistic:

```wl
ClearAll[spreadGap];
spreadGap[x1_, x2_] := Abs[StandardDeviation[x1] - StandardDeviation[x2]]/
   Sqrt[sdErr[x1]^2 + sdErr[x2]^2];
```

Read each component's gap in those units:

```wl
AssociationThread[{"Spread gap in x", "Spread gap in y",
  "Spread gap in z"},
 spreadGap[stratRunB[[All, -1, 1 ;; 3]], regRunB[[All, -1, 1 ;; 3]]]]
```

The returned numbers express each spread difference in approximate combined-standard-error units. They are useful diagnostics, but they are not a distribution-free hypothesis test and do not justify the phrase “statistically indistinguishable” without a calibrated sampling analysis.

The conditional state must stay in the unit Bloch ball. On its boundary the radial diffusion vanishes,
$(x,y,z)\cdot\vec b=\sqrt{\Gamma_{CI}}z(1-x^2-y^2-z^2)=0$.
The remaining stochastic-invariance condition is that the Itô generator applied to $r^2=x^2+y^2+z^2$ point inward or tangent there. Its boundary value is $2(x,y,z)\cdot\vec a+\vec b\cdot\vec b$:

$$\left.\mathcal L(r^2)\right|_{r=1}=(1-z^2)\,(\Gamma_{CI}+\Gamma_{BA}-2(\Gamma_d+\gamma_\phi))-\gamma_1(1-z)^2 ,$$

where $\mathcal L$ is the Itô generator.

Verify that closed form, substituting a point on the sphere and subtracting the claim:

```wl
With[{b3 = diffStrat[{x, y, z, q}, ratesSym][[1 ;; 3]]},
 Simplify[((2 {x, y, z} . driftIto[{x, y, z}, ratesSym] + b3 . b3) /.
     {x -> Sqrt[1 - z^2] Cos[\[Phi]], y -> Sqrt[1 - z^2] Sin[\[Phi]]}) -
   ((1 - z^2) (ratesSym[[1]] + ratesSym[[2]] - 2 (ratesSym[[3]] + ratesSym[[4]])) -
     ratesSym[[5]] (1 - z)^2)]]
```

The subtraction cancels identically. The relaxation term is nonpositive and vanishes quadratically at the ground pole, whereas the measurement term is linear in $1-z$ near that pole. Consequently the generator is inward or tangent at every boundary point if and only if $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$. Under nonnegative rates and the usual existence assumptions, this is the Bloch-ball invariance condition for this SDE. Every ensemble rate set below satisfies it; one single-path example deliberately violates it. The following maximum is only a numerical spot-check of the regularized solver:

```wl
Max[Map[Norm, regRunB[[All, All, 1 ;; 3]], {2}]]
```

The radius sits at one to within a hair, with no projection step. Read this as a spot-check of positivity, not a proof of it: it samples one ensemble on a fixed grid from a state that began on the surface. Violate the admissibility condition, which for a detector means $\eta>1+\gamma_\phi/\Gamma_d$, more information extracted than the dephasing paid for, and the equation stops being a valid unraveling: the trajectory can leave the ball, and its coordinates no longer describe a density matrix. That failure is one cell away. Fix rates that fail the check, twice the information the dephasing supplies, start on the sphere away from the pole where the closed form above makes the flux point outward, and read the largest radius one path reaches:

```wl
With[{ratesBad = {2., 0., 0.5, 0., 1., 0.}},
 With[{tr = oneTraj[ratesBad, {0.8, 0., 0.6, 0.}, 1., dtGenOf[ratesBad], 5, dtGenOf[ratesBad]/4]},
  Max[Norm /@ Transpose[#[Subdivide[0., 1., 200]] & /@ tr[[1 ;; 3]]]]]]
```

The radius runs past one. At the chosen initial boundary point the radial diffusion is zero while the radial generator points outward, so the model immediately violates the invariance condition. The equation still integrates; what fails is its interpretation as a density-matrix trajectory.

The boundary is just as sharp from the inside. Saturate the condition with no relaxation and no intrinsic dephasing, $\gamma_1=\gamma_\phi=0$ with both quadratures on and the drive on too, $\Gamma_{CI}+\Gamma_{BA}=2\Gamma_d$ exactly, and the closed form above vanishes identically: the sphere itself is invariant, so a pure state stays pure along every single trajectory, not merely on average. Start a path on the sphere and read how far it ever leaves:

```wl
With[{rPure = {0.352, 0.4, 0.376, 0., 0., 3.}},
 With[{tr = oneTraj[rPure, {0.6, 0., 0.8, 0.}, 3., dtGenOf[rPure], 5, dtGenOf[rPure]/4],
   s = Subdivide[0., 3., 400]},
  Max[Abs[Norm /@ Transpose[#[s] & /@ tr[[1 ;; 3]]] - 1]]]]
```

The radius sticks to the sphere at solver error, the sharpest test of the diffusion's scale in the document: a noise term wrong by any factor would walk the path off the sphere at once. On the boundary, positivity is not an inequality respected but an equality enforced.

## The readout: the average signal and the shared noise

We carried the readout $Q$ beside the state through every run but have only looked at the state. Turn to the record itself. Taking the mean of $dQ=\sqrt{\Gamma_{CI}}\,z\,dt+dW$ and using $\langle dW\rangle=0$, the average record is the integrated signal,

$$\langle Q(t)\rangle=\sqrt{\Gamma_{CI}}\int_0^t\langle z(s)\rangle\,ds,$$

a closed prediction we already hold, since $\langle z\rangle$ is the Lindblad mean. And the records are already in hand: every stored trajectory carries $Q$ as its fourth component, so slice the record ensemble out of the same run whose state we just averaged:

```wl
recQB = regRunB[[All, All, 4]];
```

The prediction is itself a closed form, not a quadrature: append $u'=\sqrt{\Gamma_{CI}}\,\langle z\rangle$ to the affine system and let `DSolve` integrate it exactly:

```wl
ClearAll[qpredB];
qpredB = Module[{sol, xf, yf, zf, uf, t},
   sol = DSolve[Join[MapThread[#1'[t] == #2 &, {{xf, yf, zf}, driftIto[{xf[t], yf[t], zf[t]}, ratesB]}],
       {uf'[t] == Sqrt[ratesB[[1]]] zf[t], xf[0] == initB[[1]], yf[0] == initB[[2]],
        zf[0] == initB[[3]], uf[0] == 0}], {xf, yf, zf, uf}, t] // First;
   Function[tt, Re[uf[tt] /. sol]]];
```

Overlay the ensemble-averaged record on that prediction:

```wl
With[{g = noiseGrid[tfB, dtGenOf[ratesB]]},
 ListLinePlot[{Transpose[{g, Mean[recQB]}], Transpose[{g, qpredB /@ g}]},
  PlotStyle -> {Automatic, Thick},
  PlotLegends -> {"ensemble mean \[LeftAngleBracket]Q\[RightAngleBracket](t)", "\!\(\*SqrtBox[\(\[CapitalGamma]CI\)]\) \[Integral]\[LeftAngleBracket]z\[RightAngleBracket] dt"},
  Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "\[LeftAngleBracket]Q\[RightAngleBracket]"},
  PlotLabel -> "the average record is the integrated signal"]]
```

As expected, the mean record follows the integrated mean signal. Each trajectory satisfies

$$
Q_i(t)=\sqrt{\Gamma_{CI}}\int_0^t z_i(s)\,ds+W_i(t),
$$

so averaging $n_{\mathrm{traj}}$ trajectories gives

$$
\bar Q(t)=\sqrt{\Gamma_{CI}}\int_0^t\bar z(s)\,ds+\bar W(t).
$$

In the infinite-ensemble limit, $\bar z\to\langle z\rangle$ and $\bar W\to0$, recovering the exact prediction $\langle Q(t)\rangle=\sqrt{\Gamma_{CI}}\int_0^t\langle z(s)\rangle\,ds$. At finite $n_{\mathrm{traj}}$, the averaged Wiener path has standard deviation $\sqrt{t/n_{\mathrm{traj}}}$, so it leaves a Brownian-like residual wander rather than cancelling exactly. For example, at $t=3$ with $n_{\mathrm{traj}}=600$, this direct noise contribution has a standard deviation of about $0.071$. The complete discrepancy also includes sampling fluctuations in $\bar z$, so its appropriate error scale is the standard error of the mean record, $s_Q(t)/\sqrt{n_{\mathrm{traj}}}$. The signal survives averaging, while the statistical fluctuations shrink as $1/\sqrt{n_{\mathrm{traj}}}$.

Quantify the agreement at each nonzero grid time $t_j$ with the standardized residual

$$
Z_Q(t_j)=
\frac{\left|\bar Q(t_j)-Q_{\mathrm{pred}}(t_j)\right|}
{s_Q(t_j)/\sqrt{n_{\mathrm{traj}}}},
$$

where $Q_{\mathrm{pred}}(t_j)=\sqrt{\Gamma_{CI}}\int_0^{t_j}\langle z(s)\rangle\,ds$ and $s_Q(t_j)$ is the sample standard deviation of the records at that time. The denominator therefore includes both the direct Wiener noise and the trajectory-to-trajectory variation of the integrated signal. Omit $t=0$, where every record is fixed at zero and the standard error vanishes, and report the largest standardized residual on the sampled grid:

```wl
With[{g = Rest[noiseGrid[tfB, dtGenOf[ratesB]]], qc = recQB[[All, 2 ;;]]},
 Max[Abs[Mean[qc] - (qpredB /@ g)]/(StandardDeviation[qc]/Sqrt[nB])]]
```

The result measures the largest pointwise discrepancy in standard-error units. Because the maximum is taken over correlated times, pointwise Gaussian thresholds do not give its false-alarm probability. It is a diagnostic of scale, not a calibrated simultaneous test and not a claim about unsampled times.

The record's spread deserves the same two-route treatment as the state's, and the native ensemble has been carrying the fourth row all along; read its final-$Q$ spread against the regularized ensemble's, in the same units as before:

```wl
spreadGap[stratRunB[[All, -1, 4]], recQB[[All, -1]]]
```

The final-record spread supplies the corresponding diffusion-sensitive comparison for the fourth row. Here the two finite samples differ on the scale estimated above.

The average removed the noise, but a single record still carries it, and carries a specific one. Since $dQ-\sqrt{\Gamma_{CI}}\,z\,dt=dW$ exactly, subtracting the integrated signal from one record must return the very Wiener path that drove the state. Build one trajectory, strip its signal, and overlay the remainder on its driving path $W(t)$:

```wl
With[{r = ratesB, dt = dtGenOf[ratesB], seed = 7},
 With[{Wf = wienerFun[tfB, dt, seed], traj = oneTraj[r, initB, tfB, dt, seed, dt/4]},
  Module[{u, t, zint},
   zint = NDSolveValue[{u'[t] == traj[[3]][t], u[0] == 0.}, u, {t, 0, tfB}];
   Plot[{traj[[4]][t] - Sqrt[r[[1]]] zint[t], Wf[t]}, {t, 0, tfB},
    PlotStyle -> {Thick, Dashed},
    PlotLegends -> Placed[{"record minus signal: Q(t) - \!\(\*SqrtBox[\(\[CapitalGamma]CI\)]\) \[Integral]z dt", "driving path W(t)"}, Below],
    Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "W"},
    PlotLabel -> "stripping the signal from the record returns the driving noise"]]]]
```

Notice that the stripped record lands on the driving path to the width of the line. This is the physical content of a measurement made literal: the fluctuations in the record are not merely correlated with the noise that kicked the state, they are that noise. The observer reads $\sqrt{\Gamma_{CI}}\int z$ as signal; the remainder is the single innovation $dW$ the detector wrote into the record and the backaction at once.

Rather than trust the line width, put a number on it, the largest gap between the stripped record and the driving path across the window; define the residual once, its solver goals adjustable, and read it at the default tolerance:

```wl
ClearAll[recordResidual];
recordResidual[goal_] := With[{r = ratesB, dt = dtGenOf[ratesB], seed = 7},
   With[{Wf = wienerFun[tfB, dt, seed], traj = oneTraj[r, initB, tfB, dt, seed, dt/4, goal]},
    Module[{u, t, zint},
     zint = NDSolveValue[{u'[t] == traj[[3]][t], u[0] == 0.}, u, {t, 0, tfB},
       AccuracyGoal -> goal, PrecisionGoal -> goal];
     With[{s = Subdivide[0., tfB, 1000]},
      Max[Abs[traj[[4]][s] - Sqrt[r[[1]]] zint[s] - Wf[s]]]]]]];
recordResidual[Automatic]
```

The same ODE defines $Q'=\sqrt{\Gamma_{CI}}z+W'$, so this is a self-consistency check, not an independent validation of the route. Its residual combines the numerical errors of the four-component solve and the auxiliary integral. Rerun the same trajectory with tighter goals:

```wl
recordResidual[10]
```

The residual decreases when the requested goals are tightened, as expected for integration error. The exact regularized ODE satisfies the identity; the reported nonzero residual is numerical.

One regime makes the tie between the monitored population and its record explicit. Switch off relaxation and drive, $\gamma_1=\Omega_x=0$. Phase backaction and transverse dephasing still affect $x$ and $y$, but not $z$. Dividing the Stratonovich $z$ row by $1-z^2$ and using the ordinary chain rule gives $d\,\mathrm{arctanh}\,z=\Gamma_{CI}z\,dt+\sqrt{\Gamma_{CI}}\,dW=\sqrt{\Gamma_{CI}}\,dQ$. Therefore

$$z(t)=\tanh\!\left(\mathrm{arctanh}\,z_0+\sqrt{\Gamma_{CI}}\,[Q(t)-Q(0)]\right).$$

The formula determines the monitored component $z$, not by itself the entire Bloch vector. First check the differential identity on the symbolic flow, with the formal $w$ standing in for the regularized noise:

```wl
With[{rFree = ReplacePart[ratesSym, {5 -> 0, 6 -> 0}]},
 Simplify[(driftStrat[{x, y, z, q}, rFree][[3]] + diffStrat[{x, y, z, q}, rFree][[3]] w)/(1 - z^2) -
   Sqrt[rFree[[1]]] (driftStrat[{x, y, z, q}, rFree][[4]] + diffStrat[{x, y, z, q}, rFree][[4]] w)]]
```

The identity holds with $\Gamma_{BA},\Gamma_d,\gamma_\phi$ symbolic: only relaxation or a drive breaks this closed relation between $z$ and $Q$. The numerical example uses $Q(0)=0$ and a bare-measurement rate set on the admissibility boundary, $\Gamma_{CI}=2\Gamma_d$ with $\gamma_\phi=0$:

```wl
ClearAll[rG, tanhResidual];
rG = {0.352, 0., 0.176, 0., 0., 0.};
tanhResidual[goal_] := With[{z0 = 0.6},
   With[{tr = oneTraj[rG, {0., 0., z0, 0.}, tfB, dtGenOf[rG], 11, dtGenOf[rG]/4, goal],
     s = Subdivide[0., tfB, 1000]},
    Max[Abs[tr[[3]][s] - Tanh[ArcTanh[z0] + Sqrt[rG[[1]]] tr[[4]][s]]]]]];
tanhResidual[Automatic]
```

Tighten the goals and the violation must follow the solver, not the noise grid:

```wl
tanhResidual[10]
```

It falls with the tolerance. For each regularized path, the ordinary chain rule makes the identity exact independently of the white-noise limit. In this regime the monitored coordinate $z$ is the record read through a $\tanh$.

## The two time steps: which one sets the answer

The simulation uses two distinct time increments: the noise-grid spacing $\Delta t_{\mathrm{gen}}$ and `NDSolve`'s maximum integration step. Which one controls the numerical result? Test the integration step first by holding the noise grid, the noise realizations, and the solver goals fixed, then comparing the final-$z$ mean and spread for two step caps: $\Delta t_{\mathrm{gen}}/4$ and $10\Delta t_{\mathrm{gen}}$. The first row uses the stored ensemble; the second recomputes the same trajectories with the relaxed step cap:

```wl
TableForm[
 Map[{Mean[#][[3]], StandardDeviation[#][[3]]} &,
  {regRunB[[All, -1, 1 ;; 3]],
   ensemble[ratesB, initB, tfB, dtGenOf[ratesB], nB, 9100, 10 dtGenOf[ratesB], 5]}],
 TableHeadings -> {{"step \!\(\*SubscriptBox[\(\[CapitalDelta]t\), \(gen\)]\)/4",
    "step 10 \!\(\*SubscriptBox[\(\[CapitalDelta]t\), \(gen\)]\)"}, {"mean z", "spread z"}}]
```

The same seeds generate the same interpolated paths in both rows. For this parameter set and these observables, relaxing the cap moves the mean and the spread by far less than their Monte-Carlo uncertainties; the residual difference between the rows is solver error at the production goals, not a noise-discretization effect. This is a case-specific sensitivity check, not evidence that `MaxStepSize` is generally irrelevant.

The regularization dial is $\Delta t_{\mathrm{gen}}$, and the native `StratonovichProcess` ensemble supplies a finite-step reference. Read its final-$z$ mean:

```wl
Mean[stratRunB[[All, -1, 3]]]
```

and its final-$z$ spread, the number the regularized route must reproduce as the grid tightens:

```wl
StandardDeviation[stratRunB[[All, -1, 3]]]
```

Check the native integrator at three step sizes. The runs use the same integer seed, but changing the step changes how random numbers are consumed, so they should be treated as separate Monte Carlo samples rather than pathwise-coupled trajectories:

```wl
natZB = Prepend[Table[stratRef[ratesB, initB, tfB, dtGenOf[ratesB]/f, nB, 2024][[All, 3]],
    {f, {2, 4}}], stratRunB[[All, -1, 3]]];
TableForm[Transpose[{dtGenOf[ratesB]/{1, 2, 4}, StandardDeviation /@ natZB}],
 TableHeadings -> {None, {"dt", "spread z"}}]
```

Read the two refined rows as distances from the working-step row, in units of the estimates' sampling errors in quadrature:

```wl
TableForm[{spreadGap[#, First[natZB]]} & /@ Rest[natZB],
 TableHeadings -> {{"dt/2", "dt/4"}, {"Spread gap in z"}}]
```

The refined estimates show no resolved trend relative to their Monte Carlo uncertainty. That supports using the working step as a practical reference at this sample size, but it does not prove convergence.

Now sweep $\Delta t_{\mathrm{gen}}$ from coarse to fine, the integrator step held below it, and tabulate the regularized mean and spread of the final $z$ against the grid step:

```wl
nSweep = 300;
sweepZB = Table[ensemble[ratesB, initB, tfB, tfB/k, nSweep, 6000, Automatic, 5][[All, 3]], {k, {6, 12, 25, 100, 400}}];
TableForm[Transpose[{N[tfB/{6, 12, 25, 100, 400}], Mean /@ sweepZB, StandardDeviation /@ sweepZB}],
 TableHeadings -> {None, {"\!\(\*SubscriptBox[\(\[CapitalDelta]t\), \(gen\)]\)", "mean z", "spread z"}}]
```

The table shows how the estimated mean and spread vary across the sampled noise meshes. Put the mean column in pointwise Monte Carlo standard-error units relative to the exact Lindblad value:

```wl
With[{lbz = lindbladBloch[ratesB, initB[[1 ;; 3]]][tfB][[3]]},
 TableForm[
  Transpose[{N[tfB/{6, 12, 25, 100, 400}],
    Abs[Mean[#] - lbz]/(StandardDeviation[#]/Sqrt[nSweep]) & /@ sweepZB}],
  TableHeadings -> {None, {"noise step", "gap of the mean from Lindblad, in standard errors"}}]]
```

The table resolves no systematic mean bias at this sample size. The spread is more sensitive to the regularization mesh, as expected because the affine mean is diffusion-blind. A finite sweep can bound effects at its sampled meshes; it cannot establish the $h\to0$ limit.

## A strong readout: the same route where the drift is large

The gap between the Itô and Stratonovich drifts scales with the measurement rate, so a strong informative readout is the stringent test of feeding `NDSolve` the Stratonovich form. Fix a strong readout from a tilted start, again a four-vector with the readout zeroed; the rates are again quantum-limited ($\eta=1$), with the intrinsic dephasing supplying the slack in the admissibility check:

```wl
ratesS = {7., 1., 4., 1., 1., 3.}; initS = {0.6, 0., 0.6, 0.}; tfS = 1.5; nS = 600;
ratesS[[1]] + ratesS[[2]] <= 2 (ratesS[[3]] + ratesS[[4]])
```

Read its efficiency:

```wl
efficiency[ratesS]
```

Run the native reference for the strong case:

```wl
stratS = stratRef[ratesS, initS, tfS, dtGenOf[ratesS], nS, 3030];
```

Run the regularized ensemble for the strong case:

```wl
regS = ensemble[ratesS, initS, tfS, dtGenOf[ratesS], nS, 5300, Automatic, 5];
```

Compare the mean of the final Bloch vector across the three routes:

```wl
TableForm[
  Transpose[{lindbladBloch[ratesS, initS[[1 ;; 3]]][tfS], Mean[stratS],
   Mean[regS]}],
  TableHeadings -> {{"Mean of x", "Mean of y",
    "Mean of z"}, {"Lindblad", "Stratonovich", "regularized"}}]
```

As one can see, the regularized mean still tracks the exact value in every row where the correction is large; the $x$ row has decayed from its initial tilt toward zero under the transverse decay, and all three routes agree on the small residual. Compare the spread across the two routes:

```wl
TableForm[
  Transpose[{StandardDeviation[stratS], StandardDeviation[regS]}],
  TableHeadings -> {{"SD of x", "SD of y", "SD of z"}, {"Stratonovich", "regularized"}}]
```

As expected, the spreads agree closely where the Stratonovich correction is largest. Put the sampling scale on all three components again:

```wl
TableForm[{spreadGap[stratS, regS]},
 TableHeadings -> {{"Spread gap"}, {"x", "y", "z"}}]
```

All three gaps sit at the sampling scale even at strong coupling. See the final-$z$ marginal agree, the regularized histogram against the native Stratonovich reference:

```wl
Histogram[{regS[[All, 3]], stratS[[All, 3]]}, 24, "PDF",
  Frame -> True, FrameLabel -> {"final z", "density"},
  PlotLabel -> "Strong readout (final value of z)", 
 PlotLegends -> {"Regularized route", "Stratonovich ensemble"}]
```

This confirms that the regularized histogram sits on top of the native Stratonovich reference: the route reproduces the final-$z$ marginal, not merely its mean, even where the drift correction is strong.

Repeat the finite mesh-sensitivity check in the strong-readout case, where the conversion correction is larger:

```wl
sweepZS = Table[ensemble[ratesS, initS, tfS, tfS/k, nSweep, 6600, Automatic, 5][[All, 3]], {k, {120, 240, 480}}];
TableForm[Transpose[{N[tfS/{120, 240, 480}], StandardDeviation /@ sweepZS}],
 TableHeadings -> {None, {"\!\(\*SubscriptBox[\(\[CapitalDelta]t\), \(gen\)]\)", "SD in z"}}]
```

Now read each row as a distance from the native strong-case spread, in units of the two estimates' sampling errors in quadrature:

```wl
TableForm[
 Transpose[{N[tfS/{120, 240, 480}], spreadGap[#, stratS[[All, 3]]] & /@ sweepZS}],
 TableHeadings -> {None, {"noise step", "gap of the spread from the native run, in sampling errors"}}]
```

These three estimates show no resolved trend beyond their Monte Carlo scale. That supports the working mesh for the displayed observable and sample size; it does not certify the forty-interval heuristic for other rates, initial states, observables, or accuracy targets.

## What is now true

The monitored qubit and its accumulated readout form one four-component SDE driven by one scalar noise. To make an equation `NDSolve` can accept, first convert the target Itô drift to the equivalent Stratonovich drift, then replace the Wiener path by a continuous piecewise-cubic approximation $W_h$. `NDSolve` receives the random ODE $\dot v_h=a_{\mathrm{Strat}}(v_h)+b(v_h)W_h'$. Under the scalar-noise Wong-Zakai conditions, its solution approaches the target Itô process through the equivalent Stratonovich representation as $h\to0$.

The exact algebra is stronger than the numerical evidence: the Itô conversion, Lindblad generator, readout identity, and Bloch-ball flux all reduce symbolically to the claimed formulas. The simulations are finite-resolution checks. Their means and spreads agree with the stated references at the displayed scales, but regularization bias, stochastic-integrator bias, ODE error, and Monte Carlo uncertainty remain distinct error sources.

Integrating the readout beside the state costs one extra equation and preserves the shared innovation explicitly. For $Q(0)=0$, its ensemble mean is $\langle Q(t)\rangle=\sqrt{\Gamma_{CI}}\int_0^t\langle z(s)\rangle\,ds$; subtracting the pathwise signal returns the same regularized Wiener path that drives the state; and with drive and relaxation absent, the monitored coordinate is a closed $\tanh$ function of the record. The Bloch-vector interpretation requires $\Gamma_{CI}+\Gamma_{BA}\le2(\Gamma_d+\gamma_\phi)$; outside that condition the equation may still integrate, but its trajectories need not represent density matrices.
