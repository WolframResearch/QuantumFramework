---
Template: Default
---

# Handing NDSolve the explicit Stratonovich SDE: the state and the readout of a monitored qubit

**A companion to *Watching a Qubit: Continuous Measurement and Quantum Trajectories* and a direct sequel to *Regularizing the Noise so NDSolve Can Integrate It*. There we smoothed the white noise so `NDSolve` could integrate the trajectory, then spent the document deriving the single drift term that separates the Itô and Stratonovich readings, so we could subtract it before handing the equation over. Here we drop that bookkeeping: the trajectory has an explicit Stratonovich form whose drift already carries the correction, so we write the drift and the diffusion out as they stand and feed `NDSolve` the equation directly. And we integrate all four components, the three Bloch numbers together with the readout $Q$ the detector accumulates, so a sample record rides alongside the state that produced it. What it buys is checked, not asserted: the built-in Stratonovich-to-Itô conversion certifies the explicit drift symbolically, the ensemble mean lands on the exact Lindblad master equation, and the ensemble spread lands on an independent native `StratonovichProcess` run, for every component of the Bloch vector, not just the monitored one.**

Mads Bahrami (last updated: August 9, 2026)

I strongly believe in a computation-first narrative: if I cannot compute it, I cannot claim to understand it, so every claim below is turned into a short cell you can run and change. The environment is a live Wolfram notebook; evaluate the cells from top to bottom, since later ones use objects the earlier ones define. Change a rate, reseed the noise, lengthen the window, and watch what moves: my suggestion is to read each output and its meaning first, before unpacking the input that produced it.

Let's start!

## The trajectory equations: the state and the readout in one SDE

Continuous homodyne monitoring of a driven qubit produces two coupled equations driven by the same noise: one for the observer's running estimate of the Bloch vector $\{x,y,z\}=\{\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle\}$, and one for the readout $Q$ the detector accumulates while measuring $\hat\sigma_z$. Stacking the readout as a fourth component, the pair is a single scalar-noise Itô stochastic differential equation, a drift $\vec a$ times $dt$ plus a diffusion $\vec b$ times one Wiener increment $dW$:

$$
d\begin{pmatrix} x \\ y \\ z \\ Q \end{pmatrix}
= \underbrace{\begin{pmatrix} -\Gamma_2^{\mathrm{eff}}\,x \\ -\Gamma_2^{\mathrm{eff}}\,y-\Omega_x\,z \\ \gamma_1(1-z)+\Omega_x\,y \\ \sqrt{\Gamma_{CI}}\,z \end{pmatrix}}_{\displaystyle \vec a}\,dt
\;+\; \underbrace{\begin{pmatrix} -\sqrt{\Gamma_{CI}}\,x z-\sqrt{\Gamma_{BA}}\,y \\ -\sqrt{\Gamma_{CI}}\,y z+\sqrt{\Gamma_{BA}}\,x \\ \sqrt{\Gamma_{CI}}\,(1-z^2) \\ 1 \end{pmatrix}}_{\displaystyle \vec b}\,dW .
$$

This is the target: the process whose mean we will check against the Lindblad master equation and whose spread we will check against a native integrator. The top three rows, the Bloch state, close on their own, nothing in them depends on $Q$; the fourth row is the record, an integral whose drift $\sqrt{\Gamma_{CI}}\,z$ carries the $\hat\sigma_z$ signal and whose noise is the same $dW$ that kicks the state.

Every parameter is a rate, in units where $\gamma_1=1$: $\Gamma_{CI}$ is the informative (coherent-information) rate, which sharpens the state toward a $\hat\sigma_z$ eigenstate; $\Gamma_{BA}$ is the backaction rate of the orthogonal quadrature, which only kicks the phase; $\Gamma_d$ is the measurement dephasing $\Gamma_m/2$; $\gamma_\phi$ is the intrinsic pure dephasing; $\gamma_1=1/T_1$ is the relaxation rate; and $\Omega_x$ is the Rabi frequency. One combination recurs, the effective transverse decay $\Gamma_2^{\mathrm{eff}}=\gamma_1/2+\gamma_\phi+\Gamma_d$. Throughout, the code takes a rate vector $r=\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$ and names its entries, so each vector field below reads term for term like these formulas.

Now the same process written the other way. For a single scalar noise, the Itô equation above is identical to a Stratonovich equation with the same diffusion and a drift lowered by one exact term (the Itô-Stratonovich conversion, unambiguous for scalar noise, with no Lévy-area subtlety). We will not compute that term. Its effect is already written into the explicit Stratonovich drift below, and that is the drift we hand `NDSolve`:

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

The reason this is the form to feed a smooth-noise solver is one sentence: a white noise replaced by a differentiable curve of vanishing correlation time integrates to the Stratonovich reading, not the Itô one, so `NDSolve`, which sees only a smooth forcing, must be given the Stratonovich drift. That is the entire design of this route, and it is why the correction never has to be computed: the smoothing puts it back. The readout row is the same in both forms, because its noise is additive, $dQ=\sqrt{\Gamma_{CI}}\,z\,dt+dW$, with nothing to reinterpret.

## Encoding the Stratonovich vector field

Encode the Stratonovich drift, all four rows. Name the entries of the rate vector $r=\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$ and fold the transverse rate $\Gamma_2^{\mathrm{eff}}$ into a local name, so each row reads term for term like the display equation:

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

Encode the Stratonovich diffusion, the same $\vec b$ as the Itô form, with the additive readout row equal to one:

```wl
ClearAll[diffStrat];
diffStrat[{x_, y_, z_, q_}, r_] :=
  With[{\[CapitalGamma]CI = r[[1]], \[CapitalGamma]BA = r[[2]]},
   {-Sqrt[\[CapitalGamma]CI] x z - Sqrt[\[CapitalGamma]BA] y,
    -Sqrt[\[CapitalGamma]CI] y z + Sqrt[\[CapitalGamma]BA] x,
    Sqrt[\[CapitalGamma]CI] (1 - z^2),
    1}];
```

To read these back as formulas, write the six rates as symbols, using string subscripts so the labels stay $\Gamma_{CI},\Gamma_{BA},\dots$ and never collide with the Bloch $x$ or fuse into a $C\cdot I$ product:

```wl
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

Fix a representative dispersive transmon, rescaled by $\gamma_1$, starting from the ground state with the readout zeroed, over a short window. The initial condition is a four-vector, the three Bloch numbers followed by the record $Q$:

```wl
ratesB = {0.352, 0., 0.176, 1., 1., 3.}; initB = {0., 0., 1., 0.}; tfB = 3.; nB = 600;
```

The reference for the mean is the Lindblad master equation. Because the Itô drift $\vec a$ is affine in $\{x,y,z\}$, the ensemble mean obeys that drift as a closed linear ODE, which `DSolve` solves analytically, with no time-stepping error; return $\langle\vec v\rangle(t)$:

```wl
ClearAll[lindbladBloch];
lindbladBloch[r_, {x0_, y0_, z0_}] :=
  With[{\[CapitalGamma]d = r[[3]], \[Gamma]\[Phi] = r[[4]], \[Gamma]1 = r[[5]], \[CapitalOmega]x = r[[6]]},
   Module[{\[CapitalGamma]2 = \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d, sol, xf, yf, zf},
    sol = DSolve[{xf'[t] == -\[CapitalGamma]2 xf[t], yf'[t] == -\[CapitalGamma]2 yf[t] - \[CapitalOmega]x zf[t],
        zf'[t] == \[Gamma]1 (1 - zf[t]) + \[CapitalOmega]x yf[t], xf[0] == x0, yf[0] == y0, zf[0] == z0},
       {xf, yf, zf}, t] // First;
    Function[tt, Re[{xf[tt], yf[tt], zf[tt]} /. sol]]]];
```

Two time steps will drive every simulation below, and they should follow the physics, not the trajectory count $n$: a step size sets accuracy, $n$ sets statistics, and they are independent axes. Fix them from the fastest rate in each case. Let $\Lambda$ be the largest of the Rabi frequency, the transverse decay $\Gamma_2^{\mathrm{eff}}$, and the total measurement rate $\Gamma_{CI}+\Gamma_{BA}$:

```wl
ClearAll[maxRate];
maxRate[r_] :=
  With[{\[CapitalGamma]CI = r[[1]], \[CapitalGamma]BA = r[[2]], \[CapitalGamma]d = r[[3]],
    \[Gamma]\[Phi] = r[[4]], \[Gamma]1 = r[[5]], \[CapitalOmega]x = r[[6]]},
   Max[\[CapitalOmega]x, \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d, \[CapitalGamma]CI + \[CapitalGamma]BA]];
```

Set the noise grid $\Delta t_{\mathrm{gen}}$ to resolve that fastest timescale with forty knots, fine enough to sit in the white-noise limit. Everything else follows: the `NDSolve` step is a quarter of it, and the native reference below runs on the same grid:

```wl
ClearAll[dtGenOf];
dtGenOf[r_] := 1/(40 maxRate[r]);
```

## Regularizing the noise: a Wiener path NDSolve can integrate

A white-noise increment $dW$ is not a function of time, it is a distribution, and `NDSolve` cannot integrate a distribution. But if we draw the noise on a grid and connect the dots into a smooth curve $W(t)$, then $W'(t)$ is an ordinary function, and $d\vec v=\vec a\,dt+\vec b\,dW$ becomes the ordinary differential equation $\dot{\vec v}=\vec a(\vec v)+\vec b(\vec v)\,W'(t)$. In other words, we replace white noise by a smooth approximation with a short but finite correlation time.

Concretely: draw i.i.d. Gaussian increments of variance $\Delta t_{\mathrm{gen}}$ on a grid of spacing $\Delta t_{\mathrm{gen}}$, accumulate them into a Wiener path, and interpolate with a cubic so $W(t)$ is a continuous curve and its derivative $W'(t)$ an ordinary function of time, no longer the distribution that white noise is:

```wl
ClearAll[wienerFun];
wienerFun[tf_, dtGen_, seed_] := BlockRandom[SeedRandom[seed];
   Module[{n = Max[1, Round[tf/dtGen]], grid, h, incr},
    grid = Subdivide[0., tf, n]; h = tf/n;
    incr = RandomVariate[NormalDistribution[0, Sqrt[h]], n];
    Interpolation[Transpose[{grid, Prepend[Accumulate[incr], 0.]}], InterpolationOrder -> 3]]];
```

Visualize one realization of the path and its derivative, the ordinary function that stands in for white noise:

```wl
With[{Wf = wienerFun[tfB, dtGenOf[ratesB], 11]},
 GraphicsRow[{
   Plot[Wf[t], {t, 0, tfB}, Frame -> True, GridLines -> Automatic,
    PlotLabel -> "regularized Wiener path W(t)", FrameLabel -> {"t", "W"}],
   Plot[Wf'[t], {t, 0, tfB}, Frame -> True, GridLines -> Automatic,
    PlotLabel -> "its derivative W'(t): the smooth noise", FrameLabel -> {"t", "W'"}]},
  ImageSize -> 640]]
```

As one can see, $W(t)$ is a continuous random walk and $W'(t)$ is a jagged but genuine function, varying on the scale $\Delta t_{\mathrm{gen}}$. There are two time increments here, and they play different roles. The first, $\Delta t_{\mathrm{gen}}$, is the grid the Gaussian increments are drawn on: it sets the correlation time of the smoothed noise, and it is the dial that fixes which stochastic equation we actually solve, converging to the ideal white-noise SDE only as $\Delta t_{\mathrm{gen}}\to 0$. The second, the `NDSolve` step, carries no physics; it only has to resolve the smooth forcing $W'(t)$ between noise knots, and we hold it below $\Delta t_{\mathrm{gen}}$ so the solver is guaranteed to see every wiggle. Which of the two actually moves the answer is a question we put to the machine below.

Assemble one trajectory. Build a smooth noise, form the random ODE $\dot{\vec v}=\vec a_{\mathrm{Strat}}(\vec v)+\vec b_{\mathrm{Strat}}(\vec v)\,W'(t)$ for the full four-vector, and hand it to `NDSolve` with its step held below $\Delta t_{\mathrm{gen}}$:

```wl
ClearAll[oneTraj];
oneTraj[r_, init_, tf_, dtGen_, seed_, maxstep_] :=
  Module[{Wf, xx, yy, zz, QQ, t, rhs, eqs},
   Wf = wienerFun[tf, dtGen, seed];
   rhs = driftStrat[{xx[t], yy[t], zz[t], QQ[t]}, r] +
     diffStrat[{xx[t], yy[t], zz[t], QQ[t]}, r] Wf'[t];
   eqs = Join[MapThread[#1'[t] == #2 &, {{xx, yy, zz, QQ}, rhs}],
     {xx[0] == init[[1]], yy[0] == init[[2]], zz[0] == init[[3]], QQ[0] == init[[4]]}];
   NDSolveValue[eqs, {xx, yy, zz, QQ}, {t, 0, tf}, MaxStepSize -> maxstep]];
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
   Plot[traj[[4]][t], {t, 0, tfB}, PlotStyle -> ColorData[97, 2],
    Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "Q"},
    PlotLabel -> "readout: the accumulated record Q(t)", ImageSize -> Medium]}]]
```

As one can see, the state path fluctuates around the Lindblad mean, as a single record must, while the readout $Q$ climbs steeply where $z$ is near one and bends as $z$ falls, its $\sqrt{\Gamma_{CI}}\,z$ signal buried under noise. The two panels share one $dW$, and that is the physical content of a measurement record: the noise written into $Q$ is the same noise that kicked $z$.

## The ensemble: matching the mean and the second moment

A single path tells us little; the mean lives in the ensemble and the second moment in how widely the paths spread. Rather than compare final values, we keep the whole trajectory from each route, average many of them, and watch the two ensemble means track the exact Lindblad curve across the entire interval. The two routes are the built-in `StratonovichProcess` integrator and the regularized `NDSolve` route, both solving the same Stratonovich SDE, and the yardstick is a single Lindblad solution from `DSolve`.

The exact mean obeys the affine Itô drift $\vec a$, which is the Lindblad master equation, and $\vec a$ is also the target the built-in conversion should return from our Stratonovich field. Encode it on the three state rows:

```wl
ClearAll[driftIto];
driftIto[{x_, y_, z_}, r_] :=
  With[{\[CapitalGamma]d = r[[3]], \[Gamma]\[Phi] = r[[4]], \[Gamma]1 = r[[5]], \[CapitalOmega]x = r[[6]]},
   With[{\[CapitalGamma]2 = \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d},
    {-\[CapitalGamma]2 x, -\[CapitalGamma]2 y - \[CapitalOmega]x z, \[Gamma]1 (1 - z) + \[CapitalOmega]x y}]];
```

Before any simulation, the central claim can be settled exactly, by the same built-in we are keeping. We never hand-compute the Itô-Stratonovich correction, but the conversion `ItoProcess[StratonovichProcess[...]]` applies it for us, rewriting our Stratonovich vector field in Itô form. If the explicit Stratonovich drift truly carries the correction, its Itô form must be the affine drift $\vec a$ on the three state rows, with the readout row untouched. Convert the full four-component field and subtract that target, under the physical assumption that the rates are nonnegative (needed only to fold $\sqrt{\Gamma_{CI}}\,\sqrt{\Gamma_{BA}}$ into a single root):

```wl
Assuming[{Subscript[\[CapitalGamma], "CI"] >= 0, Subscript[\[CapitalGamma], "BA"] >= 0},
 FullSimplify[
  (ItoProcess[StratonovichProcess[{driftStrat[{x, y, z, q}, ratesSym],
        List /@ diffStrat[{x, y, z, q}, ratesSym], {x, y, z, q}},
       {{x, y, z, q}, {x0, y0, z0, q0}}, t]]["Drift"] /. {x[t] -> x, y[t] -> y, z[t] -> z, q[t] -> q})
   - Append[driftIto[{x, y, z}, ratesSym], Sqrt[Subscript[\[CapitalGamma], "CI"]] z]]]
```

A bare vector of zeros, one per row $x,y,z,Q$: the Itô form of the Stratonovich SDE we hand `NDSolve` is exactly the affine drift on all three Bloch rows, and the readout row is unchanged, so the correction is carried in full and $Q$ is the same in both calculi. This settles the central claim at no numerical cost and more decisively than any Monte-Carlo match; the ensemble below only confirms what the algebra already proves.

Start the worker kernels:

```wl
LaunchKernels[]
```

Share the definitions the regularized trajectory needs, so `ParallelTable` can spread the paths across cores:

```wl
DistributeDefinitions[driftStrat, diffStrat, wienerFun, oneTraj]
```

Generate `ntraj` full trajectories from the built-in `StratonovichProcess` integrator, keeping every step of each path. Feed it the same Stratonovich drift and diffusion the `NDSolve` route uses, its three state rows, and ask for the order-3/2 scalar-noise method: the default is Euler-Maruyama, order 1/2, whose average drifts off the exact curve at any finite step. The process variables are Wolfram's formal symbols, unbound by construction, so they are safe placeholders (the fourth slot, the readout, is a passive `0` since the state rows never use it):

```wl
ClearAll[stratRun];
stratRun[r_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_] :=
  Module[{v = {\[FormalX][t], \[FormalY][t], \[FormalZ][t], 0}, proc, td},
   proc = StratonovichProcess[{driftStrat[v, r][[1 ;; 3]], List /@ diffStrat[v, r][[1 ;; 3]],
      {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}},
     {{\[FormalX], \[FormalY], \[FormalZ]}, {x0, y0, z0}}, t];
   td = BlockRandom[SeedRandom[seed];
     RandomFunction[proc, {0., tf, dt}, ntraj, Method -> "StochasticRungeKuttaScalarNoise"]];
   td["ValueList"]];
```

Generate the same count of regularized trajectories, each `NDSolve` solution read off on a shared time grid. Two arguments carry defaults so the one real relationship lives in one place: the integrator step is $\Delta t_{\mathrm{gen}}/4$, held a factor of four below the noise grid, and the read-off grid is the noise grid itself:

```wl
ClearAll[regRun];
regRun[r_, init_, tf_, dtGen_, ntraj_, baseSeed_, maxstep_ : Automatic, grid_ : Automatic] :=
  With[{ms = maxstep /. Automatic -> dtGen/4, g = grid /. Automatic -> Range[0., tf, dtGen]},
   ParallelTable[
    With[{tr = oneTraj[r, init, tf, dtGen, baseSeed + i, ms]},
     Table[{tr[[1]][s], tr[[2]][s], tr[[3]][s]}, {s, g}]],
    {i, ntraj}]];
```

For the harder regimes and the step-size study below we will want only the final values, so reduce each generator to them. The native final Bloch vectors:

```wl
ClearAll[stratRef];
stratRef[r_, init_, tf_, dt_, ntraj_, seed_] := stratRun[r, init, tf, dt, ntraj, seed][[All, -1]];
```

And the regularized finals paired with the largest Bloch radius each path reached, its integrator step defaulting the same way:

```wl
ClearAll[ensemble];
ensemble[r_, init_, tf_, dtGen_, ntraj_, baseSeed_, maxstep_ : Automatic] :=
  With[{p = regRun[r, init, tf, dtGen, ntraj, baseSeed, maxstep, Subdivide[0., tf, 100]]},
   {p[[All, -1]], Max[Norm /@ #] & /@ p}];
```

Run both on the representative case, on the shared grid set by `dtGenOf[ratesB]`. First the native ensemble:

```wl
stratRunB = stratRun[ratesB, initB[[1 ;; 3]], tfB, dtGenOf[ratesB], nB, 2024];
```

Then the regularized ensemble on the same grid:

```wl
regRunB = regRun[ratesB, initB, tfB, dtGenOf[ratesB], nB, 9100];
```

Average each ensemble over its trajectories and compare the two mean paths of $\langle z\rangle$ against the single exact Lindblad trajectory:

```wl
With[{g = Range[0., tfB, dtGenOf[ratesB]], lb = lindbladBloch[ratesB, initB[[1 ;; 3]]]},
 ListLinePlot[{
    Transpose[{g, Mean[stratRunB][[All, 3]]}],
    Transpose[{g, Mean[regRunB][[All, 3]]}],
    Transpose[{g, lb[#][[3]] & /@ g}]},
   PlotStyle -> {ColorData[97, 1], Directive[Dashed, ColorData[97, 2]], Directive[Thick, Black]},
   PlotLegends -> {"StratonovichProcess mean", "regularized mean", "Lindblad \[LeftAngleBracket]z\[RightAngleBracket](t)"},
   Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "\[LeftAngleBracket]z\[RightAngleBracket]"},
   PlotLabel -> "averaging the trajectories rebuilds the master equation", ImageSize -> 520]]
```

As one can see, both ensemble means, the built-in `StratonovichProcess` integrator and the regularized route, lie on the exact Lindblad curve across the whole interval, not only at the final time. Averaging over the trajectories rebuilds the deterministic master equation, which is the sense in which the smooth Lindblad flow was an average over jagged records all along.

The mean sees only the drift, so it is a weak test. The measurement backaction lives in the second moment, the spread of the final Bloch vector, rows $x,y,z$ and each row native Stratonovich then regularized:

```wl
Transpose[{StandardDeviation[stratRunB[[All, -1]]], StandardDeviation[regRunB[[All, -1]]]}]
```

As expected, the spreads agree component by component, so the regularized route reproduces the noise, not just the drift. The $x$ row is flat zero: with the backaction quadrature off ($\Gamma_{BA}=0$) and no initial $x$-coherence, nothing sources $x$, and the substance sits in $y$ and $z$.

Read the largest Bloch radius any regularized path reached. The state stays inside the unit ball because the trajectory is a positivity-preserving unraveling of a Lindblad equation: for a physically admissible measurement, one whose informative and backaction rates stay within the decoherence the channels supply so the map is completely positive, the conditional state is a genuine density matrix at every step. The noise being tangent to the sphere, $(x,y,z)\cdot\vec b=\sqrt{\Gamma_{CI}}\,z\,(1-x^2-y^2-z^2)$, is part of that but not the whole of it, since the Itô term couples the drift to the noise. All three rate sets here are admissible, so the maximum should sit at one to within the small overshoot of an unprojected scheme:

```wl
Max[Map[Norm, regRunB, {2}]]
```

The radius sits at one to within a hair, with no projection step. Read this as a spot-check of positivity, not a proof of it: it samples one ensemble on a fixed grid from a state that began on the surface. Push the measurement rates far past the decoherence budget, beyond the efficiency any real detector can have, and the equation stops being a valid unraveling: the trajectory runs off the ball and the integration diverges.

## The readout: the average signal and the shared noise

We carried the readout $Q$ beside the state through every run but have only looked at the state. Turn to the record itself. Taking the mean of $dQ=\sqrt{\Gamma_{CI}}\,z\,dt+dW$ and using $\langle dW\rangle=0$, the average record is the integrated signal,

$$\langle Q(t)\rangle=\sqrt{\Gamma_{CI}}\int_0^t\langle z(s)\rangle\,ds,$$

a closed prediction we already hold, since $\langle z\rangle$ is the Lindblad mean. Define a generator that keeps the readout $Q$ instead of dropping it, reading each full trajectory's record off the shared grid:

```wl
ClearAll[recRun];
recRun[r_, init_, tf_, dtGen_, ntraj_, baseSeed_] :=
  With[{g = Range[0., tf, dtGen], ms = dtGen/4},
   ParallelTable[
    With[{tr = oneTraj[r, init, tf, dtGen, baseSeed + i, ms]},
     Table[tr[[4]][s], {s, g}]], {i, ntraj}]];
```

Run it on the representative case, on the same grid as before:

```wl
recQB = recRun[ratesB, initB, tfB, dtGenOf[ratesB], nB, 8200];
```

Overlay the ensemble-averaged record on the closed prediction $\sqrt{\Gamma_{CI}}\int_0^t\langle z\rangle\,ds$, integrating the Lindblad $\langle z\rangle$:

```wl
With[{g = Range[0., tfB, dtGenOf[ratesB]], lb = lindbladBloch[ratesB, initB[[1 ;; 3]]]},
 Module[{u, t, qpred},
  qpred = NDSolveValue[{u'[t] == Sqrt[ratesB[[1]]] lb[t][[3]], u[0] == 0.}, u, {t, 0, tfB}];
  ListLinePlot[{Transpose[{g, Mean[recQB]}], Transpose[{g, qpred /@ g}]},
   PlotStyle -> {ColorData[97, 2], Directive[Thick, Black]},
   PlotLegends -> {"ensemble mean \[LeftAngleBracket]Q\[RightAngleBracket](t)", "\[Sqrt]\[CapitalGamma]CI \[Integral]\[LeftAngleBracket]z\[RightAngleBracket] dt"},
   Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "\[LeftAngleBracket]Q\[RightAngleBracket]"},
   PlotLabel -> "the average record is the integrated signal", ImageSize -> 520]]]
```

As one can see, the mean record climbs along the integrated $\langle z\rangle$, the two curves tracking with a residual wander that is the ensemble-averaged noise, of order $\sqrt{t/n_{\mathrm{traj}}}$ and not yet fully cancelled at this sample size. The signal is what survives averaging; the noise is what does not.

The average removed the noise, but a single record still carries it, and carries a specific one. Since $dQ-\sqrt{\Gamma_{CI}}\,z\,dt=dW$ exactly, subtracting the integrated signal from one record must return the very Wiener path that drove the state. Build one trajectory, strip its signal, and overlay the remainder on its driving path $W(t)$:

```wl
With[{r = ratesB, dt = dtGenOf[ratesB], seed = 7},
 With[{Wf = wienerFun[tfB, dt, seed], traj = oneTraj[r, initB, tfB, dt, seed, dt/4]},
  Module[{u, t, zint},
   zint = NDSolveValue[{u'[t] == traj[[3]][t], u[0] == 0.}, u, {t, 0, tfB}];
   Plot[{traj[[4]][t] - Sqrt[r[[1]]] zint[t], Wf[t]}, {t, 0, tfB},
    PlotStyle -> {Directive[Thick, ColorData[97, 2]], Directive[Dashed, Black]},
    PlotLegends -> Placed[{"record minus signal: Q(t) - \[Sqrt]\[CapitalGamma]CI \[Integral]z dt", "driving path W(t)"}, Below],
    Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "W"},
    PlotLabel -> "stripping the signal from the record returns the driving noise", ImageSize -> 520]]]]
```

As one can see, the stripped record lands on the driving path to the width of the line. This is the physical content of a measurement made literal: the fluctuations in the record are not merely correlated with the noise that kicked the state, they are that noise. The observer reads $\sqrt{\Gamma_{CI}}\int z$ as signal; the remainder is the single innovation $dW$ the detector wrote into the record and the backaction at once.

## The two time steps: which one sets the answer

The route carries two time increments, and the section title is a real question: which of them moves the answer? Put the integrator step to the test first. Hold the noise grid fixed and read the final-$z$ mean and spread with the `NDSolve` step held far below $\Delta t_{\mathrm{gen}}$, then again with it pushed well above, each row $\{\text{mean } z,\text{ spread } z\}$:

```wl
Map[{Mean[First[#]][[3]], StandardDeviation[First[#]][[3]]} &,
 {ensemble[ratesB, initB, tfB, dtGenOf[ratesB], nB, 9100, dtGenOf[ratesB]/4],
  ensemble[ratesB, initB, tfB, dtGenOf[ratesB], nB, 9100, 10 dtGenOf[ratesB]]}]
```

As one can see, with the same noise realizations under both, the two rows coincide to solver precision: for this parameter set and `NDSolve`'s default tolerances, its adaptive error control refines the step wherever $W'(t)$ varies, so relaxing `MaxStepSize` from a quarter of $\Delta t_{\mathrm{gen}}$ to ten times it leaves the observables unmoved. We still cap it below $\Delta t_{\mathrm{gen}}$ as a cheap guarantee, not because this run needs it, and it is not the binding dial.

The binding dial is $\Delta t_{\mathrm{gen}}$. Read the native reference the shrinking grid should approach, the mean and spread of the final $z$, as $\{\text{mean } z,\text{ spread } z\}$:

```wl
{Mean[stratRunB[[All, -1, 3]]], StandardDeviation[stratRunB[[All, -1, 3]]]}
```

Now sweep $\Delta t_{\mathrm{gen}}$ from coarse to fine, the integrator step held below it, and tabulate the regularized mean and spread of the final $z$ with the standard error of the mean, each row $\{\Delta t_{\mathrm{gen}},\text{ mean } z,\text{ s.e. of mean},\text{ spread } z\}$:

```wl
Table[With[{fz = First[ensemble[ratesB, initB, tfB, tfB/k, 300, 6000]][[All, 3]]},
   {N[tfB/k], Mean[fz], StandardDeviation[fz]/Sqrt[Length[fz]], StandardDeviation[fz]}],
  {k, {6, 12, 25, 100, 400}}]
```

As one can see, the mean stays put within its standard error across the sweep while the spread descends from an inflated coarse value and flattens toward the native reference as $\Delta t_{\mathrm{gen}}$ shrinks. In the white-noise limit the affine drift fixes the mean exactly; at finite $\Delta t_{\mathrm{gen}}$ any regularization bias in the mean sits below the Monte-Carlo scatter here, while the bias in the second moment is large enough to watch converge out. So of the two increments it is $\Delta t_{\mathrm{gen}}$ that sets the answer, and at this sampling accuracy the spread is where we can see it move.

## Phase backaction: the x channel comes alive

So far the backaction quadrature has been off ($\Gamma_{BA}=0$), which decoupled $x$ and left it trivial. Turn it on and the $x$ channel comes alive, but watch where: the affine Itô drift still decouples $x$'s mean, so $\langle x\rangle$ decays on its own, while the diffusion's $\sqrt{\Gamma_{BA}}\,y$ term now feeds $x$'s fluctuations from $y$. So the mean is the weak test of the coupling and the spread is the real one. Fix a case with both quadratures on, from a tilted start:

```wl
ratesX = {0.4, 0.4, 0.4, 0., 1., 2.}; initX = {0.3, 0., 0.6, 0.}; tfX = 2.; nX = 600;
```

Run the native reference:

```wl
stratX = stratRef[ratesX, initX[[1 ;; 3]], tfX, dtGenOf[ratesX], nX, 4040];
```

Run the regularized ensemble:

```wl
regX = ensemble[ratesX, initX, tfX, dtGenOf[ratesX], nX, 7100];
```

Compare the mean of the final Bloch vector, rows $x,y,z$ and each row exact Lindblad, native Stratonovich, regularized:

```wl
Transpose[{lindbladBloch[ratesX, initX[[1 ;; 3]]][tfX], Mean[stratX], Mean[regX[[1]]]}]
```

As one can see, the $x$ row is now genuinely nonzero and the three routes agree on it, its mean decaying under the drift that decouples it. Now the spread, where the backaction actually lives, rows $x,y,z$ and each row native Stratonovich then regularized:

```wl
Transpose[{StandardDeviation[stratX], StandardDeviation[regX[[1]]]}]
```

As expected, the $x$ spread is now substantial, fed by the backaction noise the $y$ channel injects, and the regularized route matches it. This is the full-Bloch-vector check the representative case could not give: with the coupling on, every component carries a nonzero second moment, and the regularized Stratonovich route reproduces all three.

## A strong readout: the same route where the drift is large

The gap between the Itô and Stratonovich drifts scales with the measurement rate, so a strong informative readout is the stringent test of feeding `NDSolve` the Stratonovich form. Fix a strong readout from a tilted start, again a four-vector with the readout zeroed:

```wl
ratesS = {8., 0., 4., 1., 1., 3.}; initS = {0.6, 0., 0.6, 0.}; tfS = 1.5; nS = 600;
```

Run the native reference for the strong case:

```wl
stratS = stratRef[ratesS, initS[[1 ;; 3]], tfS, dtGenOf[ratesS], nS, 3030];
```

Run the regularized ensemble for the strong case:

```wl
regS = ensemble[ratesS, initS, tfS, dtGenOf[ratesS], nS, 5300];
```

Compare the mean of the final Bloch vector, rows $x,y,z$ and each row exact Lindblad, native Stratonovich, regularized:

```wl
Transpose[{lindbladBloch[ratesS, initS[[1 ;; 3]]][tfS], Mean[stratS], Mean[regS[[1]]]}]
```

As one can see, the regularized mean still tracks the exact value in every row where the correction is large; the $x$ row has decayed from its initial tilt toward zero, since with $\Gamma_{BA}=0$ nothing rotates it back, and all three routes agree on the small residual. Compare the spread, rows $x,y,z$ and each row native Stratonovich then regularized:

```wl
Transpose[{StandardDeviation[stratS], StandardDeviation[regS[[1]]]}]
```

As expected, the spreads still agree row by row where the Stratonovich correction is largest. See the whole distribution agree, the regularized final-$z$ histogram against the native Stratonovich reference:

```wl
Histogram[{regS[[1]][[All, 3]], stratS[[All, 3]]}, 24, "PDF",
 ChartLegends -> {"regularized", "Stratonovich ref"},
 Frame -> True, FrameLabel -> {"final z", "density"},
 PlotLabel -> "strong readout: the regularized route matches the Stratonovich ensemble", ImageSize -> 520]
```

As one can see, the regularized distribution sits on top of the native Stratonovich reference, so the route reproduces the ensemble, not merely its average, even where the drift correction is strong.

## What is now true

The trajectory of a continuously monitored qubit is one four-component stochastic differential equation: the three Bloch numbers and the readout the detector accumulates, all driven by a single noise. Written in Stratonovich form, its drift already carries the term that separates the two stochastic calculi, which the built-in Stratonovich-to-Itô conversion confirms exactly, so a smooth-noise integrator needs nothing added: replace the white noise by a cubic-interpolated Wiener path, hand `NDSolve` the Stratonovich drift and diffusion as they stand, and the solution is the target Itô process. That target is what the two references pin down, because the Itô drift is affine the ensemble mean is the Lindblad master equation for every Bloch component, and because the diffusion is exactly $\vec b$ the ensemble spread is the native `StratonovichProcess` spread. Integrating the readout beside the state costs one more equation and buys a sample record: averaged over the ensemble it is the integrated signal $\sqrt{\Gamma_{CI}}\int\langle z\rangle$, and stripped of that signal a single record returns the very Wiener path that kicked the state, which is the physical content of a measurement, the record and the backaction are one draw. Once the backaction quadrature is on, the coupled $x$ channel it feeds carries its own second moment, and the route reproduces that too, so the full Bloch vector is matched, not just the monitored component. Of the two time increments, only the noise-generation step sets the answer, the dial that fixes which SDE is solved and is sent to zero, while the integrator step, for the tolerances used here, is absorbed by the solver's own error control.
