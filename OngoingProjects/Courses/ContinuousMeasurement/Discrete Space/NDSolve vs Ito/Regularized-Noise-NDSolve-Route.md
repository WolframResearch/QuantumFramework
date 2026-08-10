---
Template: Default
---

# Regularizing the noise so NDSolve can integrate it: a continuously measured qubit's trajectory

**A companion to *Trajectories, the Record, and the Master Equation* and a sequel to *NDSolve/WhenEvent vs ItoProcess*. There we admitted that `NDSolve` has no native Wiener input and hand-rolled Euler-Maruyama inside a `WhenEvent`. Here we take the other road: make the noise smooth, so the stochastic differential equation becomes an ordinary differential equation with a wiggly time-dependent forcing that `NDSolve` already knows how to integrate. The equation has two readings, Itô and Stratonovich, that differ by one exact drift term; a smooth noise integrates to the Stratonovich reading, so the correct route feeds `NDSolve` the Itô drift minus that term. This document writes both forms out in full, builds the corrected route from primitives, and measures what it buys: it keeps the state in the Bloch ball, reproduces the ensemble mean and second moment, and beats the hand-rolled scheme wherever the drive is stiff.**

Mads Bahrami (last updated: August 8, 2026)

## The trajectory equations: state and readout, in Itô and Stratonovich form

Continuous measurement produces two coupled stochastic equations driven by the same noise: one for the observer's running estimate of the state, one for the readout the detector writes down. In the most generic diffusive (homodyne) form, the conditional state $\rho$ obeys a stochastic master equation, and each monitored channel emits a current:

$$
d\rho = \Big(-\tfrac{i}{\hbar}[H,\rho] + \sum_k \mathcal{D}[L_k]\rho\Big)dt \;+\; \sum_k \mathcal{H}[c_k]\rho\;dW_k ,
\qquad
dQ_k = \sqrt{\eta_k}\,\big\langle c_k + c_k^\dagger\big\rangle\,dt + dW_k ,
$$

with the dissipator $\mathcal{D}[L]\rho = L\rho L^\dagger - \tfrac12\{L^\dagger L,\rho\}$ and the measurement (innovation) superoperator $\mathcal{H}[c]\rho = c\rho + \rho c^\dagger - \operatorname{Tr}[(c+c^\dagger)\rho]\,\rho$. The one link that matters: the *same* increment $dW_k$ drives both equations, so the noise in the record is the noise that kicks the state. Here $H$ is the Hamiltonian, the $L_k$ are the Lindblad operators (dissipation), the $c_k$ are the measured operators, the $\eta_k$ are the detector efficiencies with $0\le\eta_k\le1$, and the $dW_k$ are independent Wiener increments with zero mean and $dW_j\,dW_k=\delta_{jk}\,dt$.

Specialize to our qubit: a Rabi drive $H=\tfrac{\Omega_x}{2}\hat\sigma_x$, relaxation $\gamma_1\mathcal{D}[\hat\sigma_-]$, dephasing $\tfrac{\gamma_\phi+\Gamma_d}{2}\mathcal{D}[\hat\sigma_z]$, and one homodyne channel measuring $\hat\sigma_z$ through $c=\tfrac12(\sqrt{\Gamma_{CI}}-i\sqrt{\Gamma_{BA}})\hat\sigma_z$. Because $c+c^\dagger=\sqrt{\Gamma_{CI}}\,\hat\sigma_z$, the record reads out $\langle\hat\sigma_z\rangle=z$. Writing $\rho=\tfrac12(I+x\hat\sigma_x+y\hat\sigma_y+z\hat\sigma_z)$ collapses the operator equation onto the Bloch vector $\{x,y,z\}=\{\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle\}$ and the current onto one scalar. Stacking the record as a fourth component $Q$, the integrated homodyne signal, the state and its readout form one Itô SDE, a drift $\vec a$ times $dt$ plus a diffusion $\vec b$ times one scalar increment $dW$:

$$
d\begin{pmatrix} x \\ y \\ z \\ Q \end{pmatrix}
= \underbrace{\begin{pmatrix} -\Gamma_2^{\mathrm{eff}}\,x \\[2pt] -\Gamma_2^{\mathrm{eff}}\,y-\Omega_x\,z \\[2pt] \gamma_1(1-z)+\Omega_x\,y \\[2pt] \sqrt{\Gamma_{CI}}\,z \end{pmatrix}}_{\displaystyle \vec a}\,dt
\;+\; \underbrace{\begin{pmatrix} -\sqrt{\Gamma_{CI}}\,x z-\sqrt{\Gamma_{BA}}\,y \\[2pt] -\sqrt{\Gamma_{CI}}\,y z+\sqrt{\Gamma_{BA}}\,x \\[2pt] \sqrt{\Gamma_{CI}}\,(1-z^2) \\[2pt] 1 \end{pmatrix}}_{\displaystyle \vec b}\,dW .
$$

The fourth row is the record $Q$: its drift $\sqrt{\Gamma_{CI}}\,z$ carries the $\hat\sigma_z$ signal and its noise is the same $dW$. The top three rows, the Bloch state, close on their own (nothing depends on $Q$), so the code below encodes that state block and $Q$ is a passive integral we can reconstruct from a trajectory and its noise.

Every parameter is a rate, in the units where $\gamma_1=1$: $\Gamma_{CI}$ is the informative (coherent-information) rate, which sharpens the state toward a $\hat\sigma_z$ eigenstate; $\Gamma_{BA}$ is the backaction rate of the orthogonal quadrature, which only kicks the phase; $\Gamma_d$ is the measurement-induced dephasing $\Gamma_m/2$; $\gamma_\phi$ is the intrinsic pure dephasing; $\gamma_1=1/T_1$ is the relaxation rate; and $\Omega_x$ is the Rabi frequency. Two combinations recur: the transverse decay $\Gamma_2^{\mathrm{eff}}=\gamma_1/2+\gamma_\phi+\Gamma_d=1/T_2^{\mathrm{eff}}$, and the detector efficiency $\eta=(\Gamma_{CI}+\Gamma_{BA})/(2\Gamma_d)$, with $0\le\eta\le1$, which we set to one (optimal homodyne phase) throughout. The $dt$ column of $\vec a$ is the Lindblad flow and the $dW$ column of $\vec b$ is the measurement backaction.

For a single scalar noise, the Itô equation above is the *same process* as a Stratonovich equation with the same diffusion and a drift lowered by one term (this is the Itô-Stratonovich conversion, unambiguous for scalar noise, with no Levy-area subtlety):

$$
d\begin{pmatrix} x \\ y \\ z \\ Q \end{pmatrix}
= \Big(\vec a-\tfrac12(\vec b\cdot\nabla)\vec b\Big)\,dt \;+\; \vec b\circ dW .
$$

which becomes:

$$
d\begin{pmatrix} x \\ y \\ z \\ Q \end{pmatrix}
=
\begin{pmatrix}
\big(\tfrac{\Gamma_{BA}+\Gamma_{CI}}{2}-\Gamma_2^{\mathrm{eff}}-\Gamma_{CI}\,z^2\big)\,x-\sqrt{\Gamma_{BA}\Gamma_{CI}}\,y z \\[4pt]
\big(\tfrac{\Gamma_{BA}+\Gamma_{CI}}{2}-\Gamma_2^{\mathrm{eff}}-\Gamma_{CI}\,z^2\big)\,y+\sqrt{\Gamma_{BA}\Gamma_{CI}}\,x z-\Omega_x\,z \\[4pt]
\gamma_1(1-z)+\Gamma_{CI}\,z(1-z^2)+\Omega_x\,y \\[4pt]
\sqrt{\Gamma_{CI}}\,z
\end{pmatrix} dt
\;+\;
\begin{pmatrix}
-\sqrt{\Gamma_{CI}}\,x z-\sqrt{\Gamma_{BA}}\,y \\[2pt]
-\sqrt{\Gamma_{CI}}\,y z+\sqrt{\Gamma_{BA}}\,x \\[2pt]
\sqrt{\Gamma_{CI}}\,(1-z^2) \\[2pt]
1
\end{pmatrix} \circ dW .
$$

In short, writing the correction as $\vec c\equiv\tfrac12(\vec b\cdot\nabla)\vec b$, the Stratonovich drift is $\vec a-\vec c$ component by component, while the diffusion $\vec b$ is untouched, now read as a Stratonovich product $\circ\,dW$. The readout is the same in both forms, because its noise is additive, so only the state equation carries the correction. This is the whole reason the document exists: a smooth noise integrates to the Stratonovich equation, so to solve the Itô equation we will hand `NDSolve` the drift $\vec a-\vec c$, and the smoothing puts $\vec c$ back. We will compute $\vec c$ symbolically below. Let's start.

## Encoding the drift and the diffusion

Encode the state drift, the top three rows of $\vec a$, reading the rate vector positionally as $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$:

```wl
ClearAll[driftV];
driftV[{x_, y_, z_}, r_] :=
  With[{\[CapitalGamma]2 = r[[5]]/2 + r[[4]] + r[[3]]},
   {-\[CapitalGamma]2 x, -\[CapitalGamma]2 y - r[[6]] z, r[[5]] (1 - z) + r[[6]] y}];
```

Encode the state diffusion, the top three rows of $\vec b$:

```wl
ClearAll[diffV];
diffV[{x_, y_, z_}, r_] :=
  {-Sqrt[r[[1]]] x z - Sqrt[r[[2]]] y, -Sqrt[r[[1]]] y z + Sqrt[r[[2]]] x, Sqrt[r[[1]]] (1 - z^2)};
```

Fix a representative dispersive transmon, rescaled by $\gamma_1$, from a ground-state start over a short window:

```wl
ratesB = {0.352, 0., 0.176, 1., 1., 3.}; initB = {0., 0., 1.}; tfB = 3.; nB = 600;
```

The exact reference for the mean is the Lindblad master equation, which for this affine drift closes in a `DSolve` form with no integration error. Define it as a function returning $\langle\vec v\rangle(t)$, its third component the mean $\langle z\rangle$:

```wl
ClearAll[lindbladBloch];
lindbladBloch[r_, {x0_, y0_, z0_}] :=
  Module[{\[CapitalGamma]2 = r[[5]]/2 + r[[4]] + r[[3]], sol, xf, yf, zf},
   sol = DSolve[{xf'[t] == -\[CapitalGamma]2 xf[t], yf'[t] == -\[CapitalGamma]2 yf[t] - r[[6]] zf[t],
       zf'[t] == r[[5]] (1 - zf[t]) + r[[6]] yf[t], xf[0] == x0, yf[0] == y0, zf[0] == z0},
      {xf, yf, zf}, t] // First;
   Function[tt, Re[{xf[tt], yf[tt], zf[tt]} /. sol]]];
```

## The random ODE: turning a Wiener path into a forcing NDSolve can eat

A white-noise increment $dW$ is not a function of time, it is a distribution, and `NDSolve` cannot integrate a distribution. But if we draw the noise on a grid and connect the dots into a smooth curve $W(t)$, then $W'(t)$ is an ordinary function, and $d\vec v=\vec a\,dt+\vec b\,dW$ becomes the ordinary differential equation $\dot{\vec v}=\vec a(\vec v)+\vec b(\vec v)\,W'(t)$. In other words, we replace white noise by a smooth approximation with a short but finite correlation time.

Define a function that draws Gaussian increments on a grid of spacing $\Delta t_{\mathrm{gen}}$, accumulates them into a Wiener path, and interpolates with a cubic so the result is twice differentiable:

```wl
ClearAll[wienerFun];
wienerFun[tf_, dtGen_, seed_] := BlockRandom[SeedRandom[seed];
   Module[{grid = Range[0., tf, dtGen], incr},
    incr = RandomVariate[NormalDistribution[0, Sqrt[dtGen]], Length[grid] - 1];
    Interpolation[Transpose[{grid, Prepend[Accumulate[incr], 0.]}], InterpolationOrder -> 3]]];
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

As one can see, $W(t)$ is a continuous random walk and $W'(t)$ is a jagged but genuine function, varying on the scale $\Delta t_{\mathrm{gen}}$. Recall from the two forms above that this smoothing lands on the Stratonovich reading, so the drift we must feed `NDSolve` is $\vec a-\vec c$. To get $\vec c$ in closed form, write the six rates as symbols:

```wl
ratesSym = {Subscript[\[CapitalGamma], "CI"], Subscript[\[CapitalGamma], "BA"],
   Subscript[\[CapitalGamma], "d"], Subscript[\[Gamma], "\[Phi]"],
   Subscript[\[Gamma], "1"], Subscript[\[CapitalOmega], "x"]};
```

Compute the correction $\vec c=\tfrac12(\vec b\cdot\nabla)\vec b$, as half the Jacobian of $\vec b$ contracted with $\vec b$:

```wl
corrGeneral = With[{v = {x, y, z}}, FullSimplify[(1/2) Grad[diffV[v, ratesSym], v] . diffV[v, ratesSym]]]
```

As one can see, this is the correction $\vec c$ on the three state rows (its record row is zero); reading its third component, $c_z=\Gamma_{CI}\,z(z^2-1)$, it is a cubic in $z$ set by the informative rate. Turn it into a fast pure function of a Bloch point, so we can subtract it inside the solver:

```wl
ClearAll[corrFn];
corrFn[r_] := Module[{v = {\[FormalA], \[FormalB], \[FormalC]}, c},
   c = (1/2) Grad[diffV[v, r], v] . diffV[v, r];
   Function[{X, Y, Z}, Evaluate[c /. Thread[v -> {X, Y, Z}]]]];
```

Assemble one trajectory: build a smooth noise, subtract the correction from the drift, and hand the random ODE to `NDSolve` with its step held below $\Delta t_{\mathrm{gen}}$ so it resolves the wiggles of $W'(t)$:

```wl
ClearAll[oneTraj];
oneTraj[r_, init_, tf_, dtGen_, seed_, maxstep_] :=
  Module[{Wf, cF, xx, yy, zz, t, rhs, eqs},
   Wf = wienerFun[tf, dtGen, seed];
   cF = corrFn[r];
   rhs = driftV[{xx[t], yy[t], zz[t]}, r] - cF[xx[t], yy[t], zz[t]] +
     diffV[{xx[t], yy[t], zz[t]}, r] Wf'[t];
   eqs = Join[MapThread[#1'[t] == #2 &, {{xx, yy, zz}, rhs}],
     {xx[0] == init[[1]], yy[0] == init[[2]], zz[0] == init[[3]]}];
   NDSolveValue[eqs, {xx, yy, zz}, {t, 0, tf}, MaxStepSize -> maxstep]];
```

Show one trajectory's $z(t)$ against the exact Lindblad mean, building the reference curve once in the `With` so the plot variable stays out of its solver:

```wl
With[{traj = oneTraj[ratesB, initB, tfB, tfB/300, 7, tfB/1200],
   lb = lindbladBloch[ratesB, initB]},
 Plot[{traj[[3]][t], lb[t][[3]]}, {t, 0, tfB},
  PlotStyle -> {ColorData[97, 1], Directive[Thick, Black]},
  PlotLegends -> {"one smoothed trajectory z(t)", "Lindblad <z>(t)"},
  Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "z"},
  PlotLabel -> "a single regularized-noise trajectory", ImageSize -> 480]]
```

As expected, a single path does not track the mean; it fluctuates around it, because each measurement record pushes this particular qubit one way or the other. The mean lives in the ensemble, and the second moment in how widely the paths spread.

## The invariants: the ball, the mean, and what the correction touches

The correction should be a pure measurement effect, with no trace of the drive or the dissipation. Confirm that $\vec c$ is free of $\Omega_x,\gamma_1,\gamma_\phi,\Gamma_d$:

```wl
FreeQ[corrGeneral, Subscript[\[CapitalOmega], "x"] | Subscript[\[Gamma], "1"] |
   Subscript[\[Gamma], "\[Phi]"] | Subscript[\[CapitalGamma], "d"]]
```

This confirms the Itô-Stratonovich gap is set entirely by the measurement rates, so it vanishes when there is nothing to measure. Take that no-measurement limit, both quadratures off:

```wl
Simplify[corrGeneral /. {Subscript[\[CapitalGamma], "CI"] -> 0, Subscript[\[CapitalGamma], "BA"] -> 0}]
```

As expected, with no diffusion there is no calculus to reinterpret, so the correction is exactly zero. Note that turning off only the informative rate is not enough, because the phase backaction carries its own term:

```wl
Simplify[corrGeneral /. Subscript[\[CapitalGamma], "CI"] -> 0]
```

Now the ball $r\le 1$. Positivity of the qubit is exactly $r=|\{x,y,z\}|\le 1$, and the continuous flow respects it because the noise is tangent to the Bloch sphere. Show that the noise's radial component, the projection $\vec v\cdot\vec b=\sqrt{\Gamma_{CI}}\,z(1-r^2)$, vanishes on $r=1$:

```wl
FullSimplify[{x, y, z} . diffV[{x, y, z}, ratesSym] - Sqrt[Subscript[\[CapitalGamma], "CI"]] z (1 - (x^2 + y^2 + z^2))]
```

A bare zero: a noise step slides the state along the sphere and never off it. Finally the mean. The drift is affine in $\{x,y,z\}$, so the ensemble mean obeys the Lindblad ODE exactly; the fingerprint of "affine" is a vanishing second derivative:

```wl
Simplify[D[driftV[{x, y, z}, ratesSym], {{x, y, z}, 2}]]
```

An all-zero Hessian, so $d\langle\vec v\rangle/dt=\vec a(\langle\vec v\rangle)$ and matching the mean is only a weak test. The measurement backaction lives in the second moment, which we check next.

## The ensemble: matching the mean and the second moment

We compare the corrected route against two references: the exact Lindblad mean, and the native `ItoProcess` route integrated by the SDE-native scalar-noise stochastic Runge-Kutta. The Lindblad line is the exact yardstick for the mean; the native route is an independent yardstick for the spread, built on different machinery, so agreement between them is meaningful.

Define the native Itô reference, returning the final $z$ of each path:

```wl
ClearAll[itoZfinal];
itoZfinal[r_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_] :=
  Module[{proc, td},
   proc = ItoProcess[{driftV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, r],
      List /@ diffV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, r],
      {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}},
     {{\[FormalX], \[FormalY], \[FormalZ]}, {x0, y0, z0}}, t];
   td = BlockRandom[SeedRandom[seed];
     RandomFunction[proc, {0., tf, dt}, ntraj, Method -> "StochasticRungeKuttaScalarNoise"]];
   (td["ValueList"])[[All, -1, 3]]];
```

Start the worker kernels:

```wl
LaunchKernels[]
```

Share the definitions with them, so `ParallelTable` can spread the trajectories across cores:

```wl
DistributeDefinitions[driftV, diffV, corrFn, wienerFun, oneTraj]
```

Define the ensemble as many corrected trajectories on independent noise (a per-index seed), returning each path's final $z$ and its largest Bloch radius:

```wl
ClearAll[ensemble];
ensemble[r_, init_, tf_, dtGen_, ntraj_, baseSeed_, maxstep_] :=
  Transpose@ParallelTable[
    Module[{tr = oneTraj[r, init, tf, dtGen, baseSeed + i, maxstep]},
     {tr[[3]][tf], Max[Table[Norm[Through[tr[s]]], {s, Subdivide[0., tf, 100]}]]}],
    {i, ntraj}];
```

Run the native reference on the weak-readout case:

```wl
itoB = itoZfinal[ratesB, initB, tfB, tfB/600, nB, 2024];
```

Run the corrected ensemble on the same case:

```wl
corB = ensemble[ratesB, initB, tfB, tfB/300, nB, 9100, tfB/1200];
```

Compare the mean of the final $z$, in order exact Lindblad, native Itô, corrected smoothed:

```wl
{lindbladBloch[ratesB, initB][tfB][[3]], Mean[itoB], Mean[corB[[1]]]}
```

As one can see, the corrected mean lands on the exact Lindblad value, within the Monte-Carlo scatter $1/\sqrt{n}$ shared by both ensembles. That certifies the drift; now the noise. Compare the standard deviation of the final $z$, native Itô against corrected smoothed:

```wl
{StandardDeviation[itoB], StandardDeviation[corB[[1]]]}
```

As expected, the two spreads agree, so the noise term is right, not just the drift. Read the largest Bloch radius any corrected path reached:

```wl
Max[corB[[2]]]
```

The radius sits at one to within a hair, the same tiny overshoot the native stochastic Runge-Kutta shows, and no projection step was needed to hold it.

The correction scales with the measurement rate, so a strong readout is the stringent test. Fix a strong informative readout from a tilted start:

```wl
ratesS = {8., 0., 4., 1., 1., 3.}; initS = {0.6, 0., 0.6}; tfS = 1.5; nS = 600;
```

Run the native reference for the strong case:

```wl
itoS = itoZfinal[ratesS, initS, tfS, tfS/1200, nS, 3030];
```

Run the corrected ensemble for the strong case:

```wl
corS = ensemble[ratesS, initS, tfS, tfS/400, nS, 5300, tfS/1600];
```

Compare the mean of the final $z$, in order exact Lindblad, native Itô, corrected smoothed:

```wl
{lindbladBloch[ratesS, initS][tfS][[3]], Mean[itoS], Mean[corS[[1]]]}
```

As one can see, the corrected mean still tracks the exact value where the correction is large. See the whole distribution agree, the corrected final-$z$ histogram against the native Itô reference:

```wl
Histogram[{corS[[1]], itoS}, 24, "PDF",
 ChartLegends -> {"corrected smoothed", "Ito ref"},
 Frame -> True, FrameLabel -> {"final z", "density"},
 PlotLabel -> "strong readout: corrected route matches the Ito ensemble", ImageSize -> 520]
```

As one can see, the corrected distribution sits on top of the native Itô reference, so the route reproduces the ensemble, not merely its average.

## The two time steps: which one sets the answer

There are two time increments, and they do not play symmetric roles. The first, $\Delta t_{\mathrm{gen}}$, is the grid the Wiener increments are drawn on; it sets the correlation time of the regularized noise and is the dial that fixes which SDE we solve, converging to it as $\Delta t_{\mathrm{gen}}\to 0$. The second, the `NDSolve` step $h$, only has to resolve the smooth forcing between noise knots, so it must sit below $\Delta t_{\mathrm{gen}}$, or the solver steps over the wiggles and integrates a noise we never generated. In short, $h$ is the accuracy dial and $\Delta t_{\mathrm{gen}}$ is the SDE-limit dial, sent to zero in that order.

Read the native reference the sweep should approach, mean then spread:

```wl
{Mean[itoS], StandardDeviation[itoS]}
```

Halve $\Delta t_{\mathrm{gen}}$ with $h$ held below it, and tabulate the corrected mean and spread against it, as step, mean, spread:

```wl
Grid[Prepend[
  Table[With[{z = First@ensemble[ratesS, initS, tfS, tfS/k, 500, 6000, tfS/(4 k)]},
     {N[tfS/k], Mean[z], StandardDeviation[z]}], {k, {100, 200, 400}}],
  {"dtGen", "corrected mean", "corrected std"}], Frame -> All]
```

As one can see, the corrected mean and spread are stable as $\Delta t_{\mathrm{gen}}$ shrinks and sit at the native values, so the answer has converged in the SDE-limit dial and is not an artifact of the step.

## Performance: where regularizing pays for itself

Per trajectory this route costs more than the native `ItoProcess`, which integrates the whole ensemble in one vectorized call, so on a gentle problem the native route wins and the corrected smoothed route earns its keep as an independent cross-check. The place it wins outright is the stiff regime, where the companion's explicit scheme broke: Euler-Maruyama treats the fast Rabi rotation explicitly, so it is stable only above a step count $n^*=\Omega_x^2 t_f/(\Gamma_2^{\mathrm{eff}}+\gamma_1)$.

Turn the measurement off, so the single path is deterministic and the exact answer is available, and push the drive to the true transmon strength:

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

Define the bare forward-Euler update, returning the final $z$ and the largest radius:

```wl
ClearAll[emDet];
emDet[r_, init_, tf_, nsteps_] := Module[{dt = tf/nsteps, path},
   path = NestList[(# + driftV[#, r] dt) &, init, nsteps];
   {"z(tf)" -> Last[path][[3]], "max r" -> Max[Norm /@ path]}];
```

Sweep the step count across $n^*$, tabulated as steps, final $z$, largest radius:

```wl
Grid[Prepend[{#, Sequence @@ Values[emDet[ratesF, initF, tfF, #]]} & /@ {1000, 4000, 8000},
   {"nSteps", "EM z(tf)", "EM max r"}], Frame -> All]
```

As one can see, at a thousand steps Euler-Maruyama has blown past the ball entirely, and only well above $n^*$ does it become bounded, and even then it is still off the exact value, because it is only first order. Define the smoothed route with the measurement off, which is simply `NDSolve` on the drift, returning the final $z$, the largest radius, and the number of steps it took:

```wl
ClearAll[ndDet];
ndDet[r_, init_, tf_, maxstep_] := Module[{xx, yy, zz, t, eqs, sol, gr},
   eqs = Join[MapThread[#1'[t] == #2 &, {{xx, yy, zz}, driftV[{xx[t], yy[t], zz[t]}, r]}],
     {xx[0] == init[[1]], yy[0] == init[[2]], zz[0] == init[[3]]}];
   sol = NDSolve[eqs, {xx, yy, zz}, {t, 0, tf}, MaxStepSize -> maxstep][[1]];
   gr = Subdivide[0., tf, 3000];
   {"z(tf)" -> (zz[tf] /. sol),
    "max r" -> Max[Table[Norm[{xx[s], yy[s], zz[s]} /. sol], {s, gr}]],
    "NDSolve steps" -> Length[(xx /. sol)["Grid"]]}];
```

Run it at a step capped at $t_f/1000$:

```wl
ndDet[ratesF, initF, tfF, tfF/1000]
```

The contrast is the whole point: `NDSolve` reaches the exact answer to machine-level accuracy in about a thousand steps and stays in the ball, while Euler-Maruyama needs several times more just to avoid diverging. Decoupling the noise grid from the integration step frees us to use a stable, adaptive integrator on the drift, instead of being pinned to the explicit scheme's stability window.

## What is now true

The trajectory equation has two readings that differ by the single drift term $\vec c=\tfrac12(\vec b\cdot\nabla)\vec b$, a pure measurement-backaction vector whose third component is $\Gamma_{CI}\,z(z^2-1)$. A smooth noise integrates to the Stratonovich reading, so the correct way to put this SDE through `NDSolve` is to feed it the Itô drift minus $\vec c$: the smoothing adds $\vec c$ back, and the two cancel. With the term carried, the route reproduces the exact Lindblad mean, matches the native Itô route's second moment and full distribution in both weak and strong readout, and keeps the state in the Bloch ball to machine precision without any projection, because the noise was tangent all along. Of the two time increments, the noise-generation step is the dial that fixes the SDE limit and is sent to zero, while the integrator step merely sits below it to resolve the noise, and keeping them separate is what lets a stable adaptive integrator carry a stiff drive that the explicit Euler-Maruyama of the companion cannot.
