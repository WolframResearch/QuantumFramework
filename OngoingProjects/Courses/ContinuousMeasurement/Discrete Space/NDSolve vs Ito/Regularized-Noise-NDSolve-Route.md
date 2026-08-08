---
Template: Default
---

# Regularizing the noise so NDSolve can integrate it: the Wong-Zakai route for a continuously measured qubit

**A second companion to *Trajectories, the Record, and the Master Equation*, and a sequel to *NDSolve/WhenEvent vs ItoProcess*. There we admitted that `NDSolve` has no native Wiener input and hand-rolled Euler-Maruyama inside a `WhenEvent`. Here we take the other road the reader keeps suggesting: instead of teaching `NDSolve` to take discrete noise steps, we make the noise *smooth*, so the stochastic differential equation becomes an ordinary differential equation with a wiggly time-dependent forcing that `NDSolve` already knows how to integrate. The plan is simple and the payoff is real, but there is a subtlety that decides everything about the statistics: smoothing the noise silently changes the stochastic calculus from Ito to Stratonovich. This essay builds the method from primitives, shows exactly where that change enters, repairs it with one exact drift term, and then measures what the method buys us: it stays in the Bloch ball, it reproduces the ensemble mean and second moment, and it beats the hand-rolled scheme wherever the drive is stiff. Every number below is produced by a cell you can rerun, and you are encouraged to change the rates and rerun.**

Mads Bahrami (last updated: August 7, 2026)

I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. So this document is a live Wolfram notebook. Evaluate the cells top to bottom, read each output before you read the input that produced it, and treat the whole thing as one continuous experiment: some cells depend on earlier ones. The rate names in the code are ordinary letters, `{GCI, GBA, Gd, gph, g1, Ox}`, standing for the six physical rates $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$ (measurement coherent-information, backaction dephasing, measurement dephasing $\Gamma_m/2$, pure dephasing, $1/T_1$, and Rabi), so the code reads like the formulas.

## Setting up: the same qubit, the same SDE

Recall the object from the companion. A continuously monitored qubit is fully described by its Bloch vector $\{x,y,z\}=\{\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle\}$, and under homodyne detection of $\hat\sigma_z$ (with $T_1$ decay, pure dephasing, and a Rabi drive) it follows a single-noise Ito stochastic differential equation

$$
d\vec v=\vec a(\vec v)\,dt+\vec b(\vec v)\,dW,\qquad \vec v=\{x,y,z\},
$$

with drift $\vec a$ and diffusion $\vec b$. In words, the $dt$ part is the deterministic Lindblad flow and the single $dW$ part is the measurement backaction along the record's noise. Define the drift, with the transverse rate $\Gamma_2^{\mathrm{eff}}=\gamma_1/2+\gamma_\phi+\Gamma_d$ folded in:

```wl
ClearAll[driftV, diffV];
driftV[{x_, y_, z_}, {GCI_, GBA_, Gd_, gph_, g1_, Ox_}] :=
  With[{G2 = g1/2 + gph + Gd}, {-G2 x, -G2 y - Ox z, g1 (1 - z) + Ox y}];
```

And the diffusion, the single $dW$ column that all three components share:

```wl
diffV[{x_, y_, z_}, {GCI_, GBA_, Gd_, gph_, g1_, Ox_}] :=
  {-Sqrt[GCI] x z - Sqrt[GBA] y, -Sqrt[GCI] y z + Sqrt[GBA] x, Sqrt[GCI] (1 - z^2)};
```

Throughout, our clean worked case is a representative dispersive transmon with a gentle drive, rescaled by $\gamma_1$ so the physics is dimensionless. Fix it once:

```wl
ratesB = {0.352, 0., 0.176, 1., 1., 3.}; initB = {0., 0., 1.}; tfB = 3.;
```

The exact reference for the mean is the Lindblad master equation, which for this affine drift closes in a `DSolve` form with no integration error. Wrap it so its third component is $\langle z\rangle(t)$:

```wl
ClearAll[lindbladBloch];
lindbladBloch[rates_, {x0_, y0_, z0_}] :=
  Module[{GCI, GBA, Gd, gph, g1, Ox, G2, sol, xf, yf, zf},
   {GCI, GBA, Gd, gph, g1, Ox} = rates;
   G2 = g1/2 + gph + Gd;
   sol = DSolve[{xf'[t] == -G2 xf[t], yf'[t] == -G2 yf[t] - Ox zf[t],
       zf'[t] == g1 (1 - zf[t]) + Ox yf[t], xf[0] == x0, yf[0] == y0, zf[0] == z0},
      {xf, yf, zf}, t] // First;
   Function[tt, Re[{xf[tt], yf[tt], zf[tt]} /. sol]]];
```

Let's start.

## The random ODE: turning a Wiener path into a forcing NDSolve can eat

Here is the whole idea in one sentence. A white-noise increment $dW$ is not a function of time, it is a distribution, and `NDSolve` cannot integrate a distribution; but if we draw the noise on a grid and connect the dots into a smooth curve $W(t)$, then $W'(t)$ *is* an ordinary function, and $d\vec v=\vec a\,dt+\vec b\,dW$ becomes the ordinary differential equation $\dot{\vec v}=\vec a(\vec v)+\vec b(\vec v)\,W'(t)$. In other words, we replace the true white noise by a smooth approximation with a short but finite correlation time, and hand the resulting random ODE to the solver we already trust.

The regularization is nothing more than an interpolation. Draw independent Gaussian increments $\Delta W_n\sim\mathcal N(0,\sqrt{\Delta t_{\mathrm{gen}}})$ on a grid of spacing $\Delta t_{\mathrm{gen}}$, accumulate them into a Wiener path, and interpolate with a cubic so the result is twice differentiable. Define a function that returns one such smooth path:

```wl
ClearAll[wienerFun];
wienerFun[tf_, dtGen_, seed_] := BlockRandom[SeedRandom[seed];
   Module[{grid = Range[0., tf, dtGen], incr},
    incr = RandomVariate[NormalDistribution[0, Sqrt[dtGen]], Length[grid] - 1];
    Interpolation[Transpose[{grid, Prepend[Accumulate[incr], 0.]}],
     InterpolationOrder -> 3]]];
```

Visualize one realization of the path and its derivative, the smooth stand-in for white noise:

```wl
With[{Wf = wienerFun[tfB, tfB/300, 11]},
 GraphicsRow[{
   Plot[Wf[t], {t, 0, tfB}, Frame -> True, GridLines -> Automatic,
    PlotLabel -> "regularized Wiener path W(t)", FrameLabel -> {"t", "W"}],
   Plot[Wf'[t], {t, 0, tfB}, Frame -> True, GridLines -> Automatic,
    PlotLabel -> "its derivative W'(t): the smooth noise", FrameLabel -> {"t", "W'"}]},
  ImageSize -> 620]]
```

As one can see, $W(t)$ is a continuous random walk and $W'(t)$ is a jagged but genuine function, varying on the scale $\Delta t_{\mathrm{gen}}$. That scale is the one dial that matters, and we will return to it. For now, notice the promise we have just made: as $\Delta t_{\mathrm{gen}}\to 0$ the curve $W'(t)$ becomes ever more violent and approaches white noise. Whether that limit reproduces our Ito SDE is the whole question of the next section.

With the noise in hand, one trajectory is a single `NDSolve` call. Feed the random ODE, with a slot for a correction function we will motivate in a moment (pass `None` for now, which subtracts nothing):

```wl
ClearAll[oneTraj];
oneTraj[rates_, init_, tf_, dtGen_, seed_, corr_, maxstep_] :=
  Module[{Wf, cF, xx, yy, zz, t, rhs, eqs},
   Wf = wienerFun[tf, dtGen, seed];
   cF = If[corr =!= None, corr, (0 {##} &)];
   rhs = driftV[{xx[t], yy[t], zz[t]}, rates] - cF[xx[t], yy[t], zz[t]] +
     diffV[{xx[t], yy[t], zz[t]}, rates] Wf'[t];
   eqs = Join[MapThread[#1'[t] == #2 &, {{xx, yy, zz}, rhs}],
     {xx[0] == init[[1]], yy[0] == init[[2]], zz[0] == init[[3]]}];
   NDSolveValue[eqs, {xx, yy, zz}, {t, 0, tf}, MaxStepSize -> maxstep]];
```

Notice the `MaxStepSize`: it forces `NDSolve`'s internal step below $\Delta t_{\mathrm{gen}}$ so the solver actually resolves the wiggles of $W'(t)$ rather than stepping over them. That is the reader's own suggestion, and it is correct; we will make it precise in the section on the two time steps. Show one trajectory's $z(t)$ against the exact Lindblad mean:

```wl
With[{traj = oneTraj[ratesB, initB, tfB, tfB/300, 7, None, tfB/1200]},
 Plot[{traj[[3]][t], lindbladBloch[ratesB, initB][t][[3]]}, {t, 0, tfB},
  PlotStyle -> {ColorData[97, 1], Directive[Thick, Black]},
  PlotLegends -> {"one smoothed trajectory z(t)", "Lindblad <z>(t)"},
  Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "z"},
  PlotLabel -> "a single regularized-noise trajectory", ImageSize -> 480]]
```

As expected, a single path does not track the mean; it fluctuates around it, because each measurement record pushes this particular qubit one way or the other. The mean lives in the ensemble, and the second moment lives in how widely these paths spread. The next question is the sharp one: if we average many such paths, do we recover the Ito ensemble, or something else?

## The hidden calculus: smoothing lands on Stratonovich, not Ito

Here is the fact that decides everything. If you approximate white noise by a smooth process and take the correlation time to zero, the solution of the resulting ODE does not converge to the Ito solution of $d\vec v=\vec a\,dt+\vec b\,dW$. It converges to the *Stratonovich* solution. This is the Wong-Zakai theorem, and it is not a numerical artifact: it survives in the exact limit. Equivalently, smoothing the noise quietly reinterprets the multiplication $\vec b\,dW$ from the Ito rule to the Stratonovich rule, and the two rules differ by a definite drift.

For a single scalar noise, that difference is clean, with none of the Levy-area ambiguity that plagues several independent noises. The Stratonovich SDE $d\vec v=\vec a\,dt+\vec b\circ dW$ is identical to the Ito SDE whose drift is $\vec a+\tfrac12(\vec b\cdot\nabla)\vec b$. In short, the smoothed route with drift $\vec a$ secretly solves the Ito equation with an *extra* drift $\vec c=\tfrac12(\vec b\cdot\nabla)\vec b$. Compute that extra drift symbolically, as half the Jacobian of $\vec b$ contracted with $\vec b$, for general rates and a general Bloch point:

```wl
ClearAll[x, y, z, GCI, GBA, Gd, gph, g1, Ox];
corrGeneral = With[{v = {x, y, z}, b = diffV[{x, y, z}, {GCI, GBA, Gd, gph, g1, Ox}]},
   FullSimplify[(1/2) Grad[b, v] . b]]
```

As one can see, the correction is nonzero and, reading the third component, $c_z=\Gamma_{CI}\,z(z^2-1)$: a pure cubic in $z$ scaled by the measurement rate. So the smoothed route does not solve our SDE unless we cancel this term. The repair is exact algebra, not a smaller step: subtract $\vec c$ from the drift so the smoothed route's hidden Ito drift becomes $\vec a$ again. Turn the correction into a fast pure function of a Bloch point:

```wl
ClearAll[corrFn];
corrFn[rates_] := Module[{v = {\[FormalA], \[FormalB], \[FormalC]}, c},
   c = (1/2) Grad[diffV[v, rates], v] . diffV[v, rates];
   Function[{X, Y, Z}, Evaluate[c /. Thread[v -> {X, Y, Z}]]]];
```

The correction has a physical reading worth pausing on. It depends only on the measurement rates, never on the drive or the dissipation. Confirm that $\vec c$ is free of $\Omega_x,\gamma_1,\gamma_\phi,\Gamma_d$:

```wl
FreeQ[corrGeneral, Ox | g1 | gph | Gd]
```

This confirms that the Ito-Stratonovich gap here is purely a measurement-backaction effect: it is the drift the state acquires from a detector of finite bandwidth, and it vanishes when there is nothing to measure. Which raises the natural limit: does turning the measurement off remove it?

## The invariants survive: the ball, the mean, and the honest limit

A method is only trustworthy if the things that must not move do not move. Three invariants govern this SDE, and we check each by computation rather than assertion.

First, the naive limit, stated carefully so we do not fool ourselves. One might expect $\vec c\to 0$ as the informative rate $\Gamma_{CI}\to 0$. Take that limit and read the result:

```wl
Simplify[corrGeneral /. GCI -> 0]
```

Note that it is *not* zero: with the informative rate off, the phase backaction still drives a transverse Stratonovich term set by $\Gamma_{BA}$ alone, so the honest no-measurement limit is $\vec b\equiv 0$, both quadratures off, not $\Gamma_{CI}\to 0$. Confirm that killing the whole diffusion kills the correction:

```wl
Simplify[corrGeneral /. {GCI -> 0, GBA -> 0}]
```

As expected, with no diffusion there is no calculus to reinterpret, so Ito and Stratonovich coincide and the correction is exactly zero. This is the honest statement of the limit, and it is a good example of why we compute the thing instead of reciting it.

Second, the ball $r\le 1$. Positivity of the qubit is exactly $r=|\{x,y,z\}|\le 1$, and the reason the continuous flow respects it is that the noise is tangent to the Bloch sphere. Show that the noise's radial component, the projection $\vec v\cdot\vec b=\sqrt{\Gamma_{CI}}\,z(1-r^2)$, vanishes on $r=1$:

```wl
FullSimplify[{x, y, z} . diffV[{x, y, z}, {GCI, GBA, Gd, gph, g1, Ox}] -
   Sqrt[GCI] z (1 - (x^2 + y^2 + z^2))]
```

A bare zero: the identity holds, so a noise step slides the state along the sphere and never off it. We will confirm below that the discrete integrator inherits this to machine precision.

Third, the mean. The drift is affine in $\{x,y,z\}$, and an affine drift makes the ensemble mean obey a closed, deterministic ODE, exactly the Lindblad equation. The cleanest fingerprint of "affine" is that the second derivative of the drift vanishes identically:

```wl
Simplify[D[driftV[{x, y, z}, {GCI, GBA, Gd, gph, g1, Ox}], {{x, y, z}, 2}]]
```

An all-zero Hessian, so $d\langle\vec v\rangle/dt=\vec a(\langle\vec v\rangle)$ and the mean is the `DSolve` line with no approximation. This is why matching the mean is only a weak test, and why we will insist on the second moment.

## Does it reproduce the ensemble? The mean, the spread, and the ball

Now we run the method and compare it, corrected and uncorrected, against two references: the exact Lindblad mean, and the native `ItoProcess` route integrated by the SDE-native scalar-noise stochastic Runge-Kutta. The native route is the independent yardstick for the second moment; the Lindblad line is the exact yardstick for the mean.

The ensemble is many trajectories on independent noise realizations. Using a per-index seed makes the corrected and uncorrected ensembles share the same noise, so their difference is the correction alone, not Monte-Carlo scatter. (`LaunchKernels[]` starts worker kernels and `DistributeDefinitions` copies our functions to them, so `ParallelTable` spreads the trajectories across cores; drop the `Parallel` prefix to run serially.) Build the ensemble, returning each path's final $z$ and its largest Bloch radius so we can test $r\le 1$ at the same time:

```wl
LaunchKernels[];
DistributeDefinitions[driftV, diffV, corrFn, wienerFun, oneTraj];
ClearAll[ensemble];
ensemble[rates_, init_, tf_, dtGen_, ntraj_, baseSeed_, corrected_, maxstep_] :=
  With[{cF = If[corrected, corrFn[rates], None]},
   Transpose@ParallelTable[
     Module[{tr = oneTraj[rates, init, tf, dtGen, baseSeed + i, cF, maxstep]},
      {tr[[3]][tf], Max[Table[Norm[Through[tr[s]]], {s, Subdivide[0., tf, 100]}]]}],
     {i, ntraj}]];
```

The native Ito reference returns the final $z$ of each path directly:

```wl
ClearAll[itoZfinal];
itoZfinal[rates_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_] :=
  Module[{proc, td},
   proc = ItoProcess[{driftV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, rates],
      List /@ diffV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, rates],
      {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}},
     {{\[FormalX], \[FormalY], \[FormalZ]}, {x0, y0, z0}}, t];
   td = BlockRandom[SeedRandom[seed];
     RandomFunction[proc, {0., tf, dt}, ntraj,
      Method -> "StochasticRungeKuttaScalarNoise"]];
   (td["ValueList"])[[All, -1, 3]]];
```

Run all three on the weak-readout case: native Ito, smoothed uncorrected, smoothed corrected:

```wl
nB = 600;
itoB = itoZfinal[ratesB, initB, tfB, tfB/600, nB, 2024];
uncB = ensemble[ratesB, initB, tfB, tfB/300, nB, 9100, False, tfB/1200];
corB = ensemble[ratesB, initB, tfB, tfB/300, nB, 9100, True, tfB/1200];
```

Compare the exact mean against the three ensemble means, in order Lindblad, native Ito, smoothed uncorrected, smoothed corrected:

```wl
{lindbladBloch[ratesB, initB][tfB][[3]], Mean[itoB], Mean[uncB[[1]]], Mean[corB[[1]]]}
```

As one can see, in this weak-readout regime all the means sit close together: the correction $\vec c\propto\Gamma_{CI}$ is small, so the uncorrected route looks almost right. Because the corrected and uncorrected ensembles ran on the same noise, their mean difference isolates the correction alone, with no Monte-Carlo scatter. Read it:

```wl
Mean[corB[[1]] - uncB[[1]]]
```

The shift is faint, smaller than the Monte-Carlo scatter on either separate mean, so only the shared-noise pairing makes it visible at all. This is exactly the trap. A method that is quietly solving the wrong equation can hide when the wrong term is small, and the companion's own regime is precisely where it hides. To give the test teeth we turn the measurement up.

Take a strong informative readout and a tilted start, where the cubic correction is large. Run the native reference and the two smoothed ensembles:

```wl
ratesS = {8., 0., 4., 1., 1., 3.}; initS = {0.6, 0., 0.6}; tfS = 1.5; nS = 600;
itoS = itoZfinal[ratesS, initS, tfS, tfS/1200, nS, 3030];
uncS = ensemble[ratesS, initS, tfS, tfS/400, nS, 5300, False, tfS/1600];
corS = ensemble[ratesS, initS, tfS, tfS/400, nS, 5300, True, tfS/1600];
```

Compare the exact mean against the two smoothed means, in order Lindblad, smoothed uncorrected, smoothed corrected:

```wl
{lindbladBloch[ratesS, initS][tfS][[3]], Mean[uncS[[1]]], Mean[corS[[1]]]}
```

As one can see, the corrected mean lands on the exact Lindblad value while the uncorrected mean sits well below it. Measure that gap in units of the Monte-Carlo standard error, so it is not a matter of taste:

```wl
(lindbladBloch[ratesS, initS][tfS][[3]] - Mean[uncS[[1]]])/(StandardDeviation[uncS[[1]]]/Sqrt[nS])
```

Now the failure is loud: the uncorrected mean misses the exact Lindblad value by many standard errors. The Stratonovich drift points the state toward the equator, partially cancelling the measurement's collapse toward the poles, so it under-measures. That same under-measurement shows up in the spread, which is the quantity the companion used as its sharp test. Read the standard deviation of the final $z$ for the three routes, in order native Ito, smoothed uncorrected, smoothed corrected:

```wl
{StandardDeviation[itoS], StandardDeviation[uncS[[1]]], StandardDeviation[corS[[1]]]}
```

As expected, the uncorrected spread is visibly narrower, while the corrected spread matches the native Ito reference. The correction restores both the first and the second moment. See the same thing as a picture, the final-$z$ distributions side by side:

```wl
Histogram[{uncS[[1]], corS[[1]], itoS}, 24, "PDF",
 ChartLegends -> {"smoothed uncorrected", "smoothed corrected", "Ito ref"},
 Frame -> True, FrameLabel -> {"final z", "density"},
 PlotLabel -> "strong readout: the correction restores the spread",
 ImageSize -> 520]
```

As one can see, the uncorrected histogram is the narrow one, while the corrected distribution sits squarely on top of the native Ito reference.

Finally the ball. Read the largest Bloch radius reached by any path in the strong-readout ensembles, uncorrected then corrected:

```wl
{Max[uncS[[2]]], Max[corS[[2]]]}
```

Both sit at $r=1$ to within a hair, the same tiny non-norm-preserving overshoot the native stochastic Runge-Kutta shows. The continuous flow keeps the state in the ball because the noise is tangent, and the smooth integrator inherits that without any projection step, unlike the explicit Euler-Maruyama of the companion, which needed one.

## The two time steps: which one sets the answer

The reader who proposed this method asked the right question: there are two time increments, and the `NDSolve` step should be smaller than the one used to generate the noise. That instinct is correct, and it is worth stating precisely, because the two increments do not play symmetric roles.

The first, $\Delta t_{\mathrm{gen}}$, is the grid on which we draw the Wiener increments. It sets the correlation time of the regularized noise, and it is the knob that controls *which SDE we are solving and how accurately*: as $\Delta t_{\mathrm{gen}}\to 0$ the smoothed process approaches white noise and the trajectory converges to the Stratonovich SDE. The second, the `NDSolve` step $h$, only has to resolve the smooth forcing between noise knots. Because $W'(t)$ varies on the scale $\Delta t_{\mathrm{gen}}$, we need $h<\Delta t_{\mathrm{gen}}$, or `NDSolve` steps over the wiggles and integrates a noise we never generated. In short: $h$ is the ODE-accuracy dial, $\Delta t_{\mathrm{gen}}$ is the SDE-limit dial, and they must be sent to zero in that order, $h$ first.

Confirm that $\Delta t_{\mathrm{gen}}$ is the dial that moves the statistics by halving it and watching the corrected mean and spread settle toward the native Ito reference:

```wl
Grid[Prepend[
  Table[With[{z = First@ensemble[ratesS, initS, tfS, tfS/k, 400, 6000, True, tfS/(4 k)]},
     {N[tfS/k], Mean[z], StandardDeviation[z]}], {k, {100, 200, 400}}],
  {"dtGen", "corrected mean", "corrected std"}], Frame -> All]
```

As one can see, the corrected mean and spread are stable as $\Delta t_{\mathrm{gen}}$ shrinks with $h$ held safely below it, which tells us the answer has converged in the SDE-limit dial and is not an artifact of the step. The one thing neither dial can do is remove the Ito-Stratonovich gap: shrinking $h$ and $\Delta t_{\mathrm{gen}}$ makes the smoothed route converge to Stratonovich, and only subtracting $\vec c$ moves it to Ito. That is why the correction is algebra, not refinement.

## Performance: where regularizing pays for itself

We paid a price and we should say what we bought. Per trajectory, this route is more expensive than the native `ItoProcess`, which integrates the whole ensemble in one vectorized call, while here each path is its own `NDSolve`. So on a gentle, well-conditioned problem the native route wins on cost, and this one earns its keep as an *independent* cross-check: it shares no discretization with the stochastic Runge-Kutta, so agreement between the two is meaningful.

The place it wins outright is the stiff regime, which is exactly where the companion's hand-rolled scheme broke. Explicit Euler-Maruyama treats the fast Rabi rotation explicitly, so it is stable only when the step is small enough, below a threshold $n^*=\Omega_x^2 t_f/(\Gamma_2^{\mathrm{eff}}+\gamma_1)$ in step count. Turn the measurement off so the single path is deterministic and the exact Lindblad answer is available, and push the drive to the true transmon strength:

```wl
ratesF = {0., 0., 0., 1., 1., 60 Pi}; initF = N[{Sin[0.001], 0., Cos[0.001]}]; tfF = 1./3;
```

Read the exact final $z$ this deterministic path must reach:

```wl
lindbladBloch[ratesF, initF][tfF][[3]]
```

Compute the explicit-Euler stability threshold, the step count below which the fast rotation is amplified rather than resolved:

```wl
nStar = ratesF[[6]]^2 tfF/((ratesF[[5]]/2 + ratesF[[4]] + ratesF[[3]]) + ratesF[[5]])
```

The threshold is several thousand steps. Watch explicit Euler-Maruyama cross it: below $n^*$ it amplifies the rotation and leaves the ball entirely. Define the bare forward-Euler update and sweep the step count:

```wl
ClearAll[emDet];
emDet[rates_, init_, tf_, nsteps_] := Module[{dt = tf/nsteps, path},
   path = NestList[(# + driftV[#, rates] dt) &, init, nsteps];
   {"z(tf)" -> Last[path][[3]], "max r" -> Max[Norm /@ path]}];
Grid[Prepend[{#, Sequence @@ Values[emDet[ratesF, initF, tfF, #]]} & /@ {1000, 4000, 8000},
   {"nSteps", "EM z(tf)", "EM max r"}], Frame -> All]
```

As one can see, at a thousand steps Euler-Maruyama has blown up past the ball entirely, and only well above $n^*$ does it become bounded, and even then it is still visibly off the exact value, because it is only first order. Now the smoothed route, which with the measurement off is simply `NDSolve` on the drift, integrated with an adaptive step capped at $t_f/1000$. Read its final $z$, its largest radius, and the number of steps it actually took:

```wl
ClearAll[ndDet];
ndDet[rates_, init_, tf_, maxstep_] := Module[{xx, yy, zz, t, eqs, sol, gr},
   eqs = Join[MapThread[#1'[t] == #2 &,
      {{xx, yy, zz}, driftV[{xx[t], yy[t], zz[t]}, rates]}],
     {xx[0] == init[[1]], yy[0] == init[[2]], zz[0] == init[[3]]}];
   sol = NDSolve[eqs, {xx, yy, zz}, {t, 0, tf}, MaxStepSize -> maxstep][[1]];
   gr = Subdivide[0., tf, 3000];
   {"z(tf)" -> (zz[tf] /. sol),
    "max r" -> Max[Table[Norm[{xx[s], yy[s], zz[s]} /. sol], {s, gr}]],
    "steps" -> Length[(xx /. sol)["Grid"]]}];
ndDet[ratesF, initF, tfF, tfF/1000]
```

The contrast is the whole point: `NDSolve` reaches the exact answer to machine-level accuracy in about a thousand steps and stays in the ball, while Euler-Maruyama needs several times more just to avoid diverging and is still inaccurate when it gets there. The reason is that decoupling the noise grid from the integration step frees us to use a stable, adaptive integrator on the drift, instead of being pinned to the explicit scheme's stability window. So the regularized route is not only an independent statistical cross-check; it is the better tool wherever the deterministic part is stiff.

## What is now true

The regularized-noise route is a legitimate way to integrate a continuously monitored qubit's trajectory with `NDSolve`, but it changes the stochastic calculus under the hood: smoothing a single scalar noise and taking its correlation time to zero converges to the Stratonovich reading, so the route silently solves the Ito SDE with an extra drift $\tfrac12(\vec b\cdot\nabla)\vec b$, a pure measurement-backaction term that in the third Bloch component is $\Gamma_{CI}z(z^2-1)$ and that no step size can remove. Subtracting it is exact algebra and restores the Ito ensemble: the mean returns to the Lindblad line, the second moment returns to the native stochastic Runge-Kutta value, and the state stays in the Bloch ball to machine precision without any projection, because the noise was tangent all along. Of the two time increments, the noise-generation step is the dial that fixes the SDE limit and must be sent to zero, while the integrator step must merely sit below it to resolve the noise, and keeping them separate is exactly what lets a stable adaptive integrator handle a stiff drive that the explicit Euler-Maruyama of the companion cannot. The correction is invisible when the measurement is weak and worth many standard errors when it is strong, which is the sense in which a quiet failure and a loud one are the same bug seen at two readout strengths.
