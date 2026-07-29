---
Template: Default
---

# NDSolve/WhenEvent vs ItoProcess: a continuously measured qubit's trajectory

**A companion to *Trajectories, the Record, and the Master Equation*. There the free particle is watched on a grid; here the object is a continuously monitored qubit, whose complete state is the three Bloch numbers $\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle$ living in the ball $r=|\{x,y,z\}|\le 1$. The trajectory is a Cartesian-Bloch Ito SDE, integrated natively by the built-in `RandomFunction[ItoProcess[...]]`. `NDSolve` has no native Wiener input, so this document hand-rolls Euler-Maruyama inside `NDSolve` via `WhenEvent`, and asks a sharp question: does the hand-rolled route reproduce the built-in ensemble and the Lindblad master equation, while keeping $r\le 1$, and at what cost? Every number below is produced by a cell you can rerun.**

Mads Bahrami (last updated: July 29, 2026)

This document is nearly self-contained: evaluate the cells top to bottom, and the only external ingredient is one Wolfram Function Repository function, `ItoProcessChangeVariables`, which performs the Ito change of variables for the $r\le 1$ proof and is fetched automatically on first use. The six rates $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$ (measurement coherent-information, backaction dephasing, measurement dephasing $\Gamma_m/2$, pure dephasing, $1/T_1$, and Rabi) stay symbolic where the algebra is symbolic and become numbers only when we simulate. The code uses these same Greek letters for the rates, so it reads like the formulas.

## From the density matrix to the Bloch SDE

The trajectory follows from continuous homodyne detection of $\hat\sigma_z$ (Gambetta07), with $T_1$ decay, pure dephasing, and a Rabi drive. *Homodyne detection* beats the cavity's output field against a strong local-oscillator reference at the same frequency, so the recorded photocurrent reads out one chosen field quadrature, here the one carrying the qubit's $\hat\sigma_z$ (which-state) information. Write the effective transverse rate as $\Gamma_2^{\mathrm{eff}}=\gamma_1/2+\gamma_\phi+\Gamma_d=1/T_2^{\mathrm{eff}}$. The object actually being tracked is the conditional density matrix $\rho$, which obeys a stochastic master equation, a deterministic Lindblad drift plus a state-dependent measurement backaction along the record's noise:

$$
d\rho = \Big(\underbrace{-i[H,\rho] + \gamma_1\,\mathcal{D}[\hat\sigma_-]\rho + \tfrac{\gamma_\phi+\Gamma_d}{2}\,\mathcal{D}[\hat\sigma_z]\rho}_{\text{Lindblad drift}}\Big)\,dt \;+\; \underbrace{\mathcal{H}[c]\rho}_{\text{homodyne backaction}}\,dW,
$$

with the Rabi Hamiltonian $H=\tfrac{\Omega_x}{2}\hat\sigma_x$, the GKSL dissipator $\mathcal{D}[L]\rho = L\rho L^\dagger-\tfrac12\{L^\dagger L,\rho\}$, and the homodyne measurement (innovation) superoperator $\mathcal{H}[c]\rho = c\rho+\rho c^\dagger-\operatorname{Tr}[(c+c^\dagger)\rho]\,\rho$ (Wiseman-Milburn; the calligraphic $\mathcal{H}$ is the measurement superoperator, not the Hamiltonian $H$), with measurement operator $c=\tfrac12\big(\sqrt{\Gamma_{CI}}-i\sqrt{\Gamma_{BA}}\big)\hat\sigma_z$. Here $\Gamma_{CI}$ weights the informative quadrature that sharpens the state toward a $\hat\sigma_z$ eigenstate, and $\Gamma_{BA}$ the orthogonal quadrature that only kicks the phase.

Projected onto the Bloch vector $\{x,y,z\}=\{\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle\}$, the Lindblad drift becomes the $dt$ terms and the backaction the $dW$ terms of a single scalar-noise Ito SDE:

$$
d\begin{pmatrix} x \\ y \\ z \end{pmatrix}
= \underbrace{\begin{pmatrix} -\Gamma_2^{\mathrm{eff}}\,x \\[2pt] -\Gamma_2^{\mathrm{eff}}\,y-\Omega_x z \\[2pt] \gamma_1(1-z)+\Omega_x y \end{pmatrix}}_{\displaystyle \vec a}\,dt
\;+\; \underbrace{\begin{pmatrix} -\sqrt{\Gamma_{CI}}\,xz-\sqrt{\Gamma_{BA}}\,y \\[2pt] -\sqrt{\Gamma_{CI}}\,yz+\sqrt{\Gamma_{BA}}\,x \\[2pt] \sqrt{\Gamma_{CI}}\,(1-z^2) \end{pmatrix}}_{\displaystyle \vec b}\,dW.
$$

This is the standard split $d\vec v=\vec a(\vec v)\,dt+\vec b(\vec v)\,dW$ with $\vec v=\{x,y,z\}$: the two labeled columns are the drift vector $\vec a$ (the $dt$ part) and the diffusion vector $\vec b$ (the $dW$ part, a single column since one scalar noise $dW$ drives all three components). In code, the drift, with the transverse rate $\Gamma_2^{\mathrm{eff}}$ folded in:

```wl
driftV[{x_, y_, z_}, {\[CapitalGamma]CI_, \[CapitalGamma]BA_, \[CapitalGamma]d_, \[Gamma]\[Phi]_, \[Gamma]1_, \[CapitalOmega]x_}] :=
  With[{\[CapitalGamma]2eff = \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d},
   {-\[CapitalGamma]2eff x, -\[CapitalGamma]2eff y - \[CapitalOmega]x z, \[Gamma]1 (1 - z) + \[CapitalOmega]x y}];
```

And the diffusion $\vec b$, the second column above:

```wl
diffV[{x_, y_, z_}, {\[CapitalGamma]CI_, \[CapitalGamma]BA_, \[CapitalGamma]d_, \[Gamma]\[Phi]_, \[Gamma]1_, \[CapitalOmega]x_}] :=
  {-Sqrt[\[CapitalGamma]CI] x z - Sqrt[\[CapitalGamma]BA] y,
   -Sqrt[\[CapitalGamma]CI] y z + Sqrt[\[CapitalGamma]BA] x,
   Sqrt[\[CapitalGamma]CI] (1 - z^2)};
```

The *native* route hands $\{\vec a,\vec b\}$ to `ItoProcess` and integrates with `RandomFunction` using the SDE-native scalar-noise stochastic Runge-Kutta. The *hand-rolled* route freezes the ODE (`x'[t]==0`) inside `NDSolve` and injects one explicit Euler-Maruyama step at each $dt$ through a `WhenEvent`, adding $\vec a\,dt+\vec b\,dW$ to the current state. Two mistakes to avoid: do not integrate the drift in the ODE and re-add it in the event (double count), and the event must ADD the increment to the state, not replace the state with it.

## The mean is the Lindblad equation, exactly

The drift is *affine* in $\{x,y,z\}$: $\vec a(\vec v)=A\vec v+\vec c$. Taking the Ito expectation kills the $dW$ term ($\mathbb{E}[\vec b\,dW]=0$), so the mean $\langle \vec v\rangle=\mathbb{E}[\vec v]$ obeys a closed, deterministic linear ODE,

$$
\frac{d\langle x\rangle}{dt}=-\Gamma_2^{\mathrm{eff}}\langle x\rangle,\quad
\frac{d\langle y\rangle}{dt}=-\Gamma_2^{\mathrm{eff}}\langle y\rangle-\Omega_x\langle z\rangle,\quad
\frac{d\langle z\rangle}{dt}=\gamma_1(1-\langle z\rangle)+\Omega_x\langle y\rangle,
$$

which is exactly the Bloch form of the Lindblad drift above: transverse decay at $\Gamma_2^{\mathrm{eff}}$, relaxation to the ground state at $\gamma_1$, and the Rabi rotation $\Omega_x$. Solve it symbolically with `DSolve`:

```wl
Clear[\[CapitalGamma]2eff, \[Gamma]1, \[CapitalOmega]x, x0, y0, z0];
meanSol = DSolve[{X'[t] == -\[CapitalGamma]2eff X[t], Y'[t] == -\[CapitalGamma]2eff Y[t] - \[CapitalOmega]x Z[t],
    Z'[t] == \[Gamma]1 (1 - Z[t]) + \[CapitalOmega]x Y[t], X[0] == x0, Y[0] == y0, Z[0] == z0},
   {X, Y, Z}, t] // First;
```

The transverse component decouples and decays as a pure exponential, $\langle x\rangle(t)=x_0\,e^{-\Gamma_2^{\mathrm{eff}}t}$, while $\langle y\rangle,\langle z\rangle$ form a damped Rabi pair:

```wl
X[t] /. meanSol
```

So the ensemble mean of the trajectory SDE is the Lindblad master equation, identically. The mean sees only the drift; it is blind to any diffusion that leaves the drift alone, which makes matching it a weak test. We sharpen it with the second moment below.

## The ball $r\le 1$ is an invariant, proven

Positivity of the qubit density matrix is exactly $r\le 1$, because $\rho=(I+x\hat\sigma_x+y\hat\sigma_y+z\hat\sigma_z)/2$ has $\det\rho=(1-r^2)/4$, which is non-negative iff $r\le 1$. The cleanest way to see that the SDE respects this is to change coordinates from Cartesian $\{x,y,z\}$ to spherical $\{r,\theta,\phi\}$ and read off the radius's own one-dimensional Ito equation $dr=a_r\,dt+b_r\,dW$: its two coefficients settle the invariant directly. Build the Cartesian process symbolically first:

```wl
Clear[x, y, z, x0, y0, z0, \[CapitalGamma]CI, \[CapitalGamma]BA, \[CapitalGamma]d, \[Gamma]\[Phi], \[Gamma]1, \[CapitalOmega]x];
ratesSym = {\[CapitalGamma]CI, \[CapitalGamma]BA, \[CapitalGamma]d, \[Gamma]\[Phi], \[Gamma]1, \[CapitalOmega]x};
procSym = ItoProcess[{driftV[{x, y, z}, ratesSym], List /@ diffV[{x, y, z}, ratesSym], {x, y, z}},
   {{x, y, z}, {x0, y0, z0}}, t];
```

The Wolfram Function Repository's `ItoProcessChangeVariables` applies Ito's rule, including the second-order correction, to rewrite the SDE in spherical coordinates. This is the one external ingredient in the document, fetched automatically on first use:

```wl
sphProc = ResourceFunction["ItoProcessChangeVariables"][procSym, {r, \[Theta], \[Phi]}, "Cartesian" -> "Spherical"];
```

First the noise. The radial diffusion coefficient $b_r$ is the top entry of the transformed diffusion:

```wl
FullSimplify[sphProc["Diffusion"][[1, 1]], Assumptions -> {r[t] > 0}]
```

It is proportional to $(1-r^2)$, namely $b_r=\sqrt{\Gamma_{CI}}\,\cos\theta\,(1-r^2)$, so it vanishes on the surface $r=1$ for every choice of rates. The measurement noise is tangent to the Bloch sphere: a noise step slides the state along the sphere and never off it.

Second the drift. Evaluate the radial drift $a_r$ on the surface $r=1$:

```wl
Clear[\[Eta]];
radialDrift = FullSimplify[sphProc["Drift"][[1]] /. r[t] -> 1, Assumptions -> {0 <= \[Theta][t] <= Pi}];
```

Write the measurement efficiency as $\eta=(\Gamma_{CI}+\Gamma_{BA})/(2\Gamma_d)$ and ask `Reduce` whether any point on the sphere, for any physical efficiency $\eta\in[0,1]$ and non-negative rates, can make the radial drift point outward, $a_r>0$:

```wl
Reduce[(radialDrift /. \[CapitalGamma]BA -> 2 \[Eta] \[CapitalGamma]d - \[CapitalGamma]CI) > 0 &&
   0 <= \[Theta][t] <= Pi && \[Gamma]1 >= 0 && \[Gamma]\[Phi] >= 0 && \[CapitalGamma]d >= 0 &&
   \[CapitalGamma]CI >= 0 && 2 \[Eta] \[CapitalGamma]d - \[CapitalGamma]CI >= 0 && 0 <= \[Eta] <= 1, {\[Theta][t]}]
```

`False`: no such point exists, so $a_r\le 0$ everywhere on $r=1$. The noise cannot push the state off the sphere and the drift can only pull it inward, so the exact continuous flow never leaves the ball. Discrete integrators are another matter, which is the whole point of what follows.

## The three integrators

The native route builds an `ItoProcess` from the drift, the diffusion (as a one-column matrix, since the noise is scalar), and the state variables, then lets `RandomFunction` integrate `ntraj` paths on the grid $\{0,t_f,dt\}$ with the SDE-native scalar-noise stochastic Runge-Kutta. The process variables are Wolfram's formal symbols, unbound by construction, so they are safe placeholders:

```wl
itoEnsemble[rates_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_] :=
  Module[{proc, td, vals},
   proc = ItoProcess[{driftV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, rates],
      List /@ diffV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, rates],
      {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}},
     {{\[FormalX], \[FormalY], \[FormalZ]}, {x0, y0, z0}}, t];
   td = BlockRandom[SeedRandom[seed];
     RandomFunction[proc, {0., tf, dt}, ntraj, Method -> "StochasticRungeKuttaScalarNoise"]];
   vals = td["ValueList"];   (* full data: ntraj x npts x 3, the {x,y,z} of every path *)
   <|"trajectories" -> vals, "times" -> td["Times"],
     "meanZ" -> Mean[vals[[All, All, 3]]], "zTrajs" -> vals[[All, All, 3]],
     "maxr" -> Max[Map[Norm, vals, {2}]]|>];
```

Both ensemble builders keep the full simulated data: the returned association carries a `trajectories` key, an `ntraj`$\times$`npts`$\times 3$ array of the Bloch vector $\{x,y,z\}$ at every path and step, and `meanZ`, `zTrajs`, and `maxr` are convenience reductions of it.

The hand-rolled route is explicit Euler-Maruyama. For the SDE $d\vec v=\vec a\,dt+\vec b\,dW$ the update is

$$
\vec v_{n+1}=\vec v_n+\vec a(\vec v_n)\,\Delta t+\vec b(\vec v_n)\,\Delta W_n,\qquad \Delta W_n\sim\mathcal{N}(0,\sqrt{\Delta t}),
$$

with an optional projection $\vec v\mapsto \vec v/\max(1,|\vec v|)$ that snaps any overshoot back onto the ball. One `NDSolve` carries three length-`ntraj` vectors; a single `WhenEvent` at each grid point draws `ntraj` Wiener increments $\Delta W$ (one per trajectory), reuses `driftV`/`diffV`, and applies the update. The midpoint event phase `Mod[t+dt/2,dt]==0` fires exactly $N=t_f/\Delta t$ times and samples cleanly on the grid:

```wl
emEnsemble[rates_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_, project_] :=
  Module[{xs, ys, zs, sol, grid, xyzg},
   sol = BlockRandom[SeedRandom[seed];
     NDSolve[{xs'[t] == ConstantArray[0., ntraj], ys'[t] == ConstantArray[0., ntraj],
       zs'[t] == ConstantArray[0., ntraj],
       WhenEvent[Mod[t + dt/2, dt] == 0,
        Block[{\[CapitalDelta]W = RandomVariate[NormalDistribution[0, Sqrt[dt]], ntraj], xyz, new, capr},
         xyz = {xs[t], ys[t], zs[t]};
         new = xyz + MapThread[#1 dt + #2 \[CapitalDelta]W &, {driftV[xyz, rates], diffV[xyz, rates]}];
         capr = If[project, Clip[Sqrt[new[[1]]^2 + new[[2]]^2 + new[[3]]^2], {1., Infinity}],
                   ConstantArray[1., ntraj]];
         {xs[t], ys[t], zs[t]} -> (#/capr & /@ new)]],
       xs[0] == ConstantArray[N@x0, ntraj], ys[0] == ConstantArray[N@y0, ntraj],
       zs[0] == ConstantArray[N@z0, ntraj]},
      {xs, ys, zs}, {t, 0, tf}, MaxStepSize -> dt][[1]]];
   grid = Range[0., tf, dt];
   xyzg = ({xs[#], ys[#], zs[#]} /. sol) & /@ grid;   (* npts x 3 x ntraj *)
   <|"trajectories" -> Transpose[xyzg, {2, 3, 1}],   (* full data: ntraj x npts x 3 *)
     "times" -> grid, "meanZ" -> Mean /@ xyzg[[All, 3]],
     "zTrajs" -> Transpose[xyzg[[All, 3]]],
     "maxr" -> Max[Sqrt[xyzg[[All, 1]]^2 + xyzg[[All, 2]]^2 + xyzg[[All, 3]]^2]]|>];
```

The exact reference is the mean ODE from the first section, which closes in a `DSolve` closed form (no integration error). Wrapped as a function of the rates, initial state, and time, it returns the full mean Bloch vector $\langle\vec v\rangle(t)=\{\langle x\rangle,\langle y\rangle,\langle z\rangle\}(t)$; its third component is $\langle z\rangle$:

```wl
lindbladBloch[rates_, {x0_, y0_, z0_}] :=
  Module[{\[CapitalGamma]CI, \[CapitalGamma]BA, \[CapitalGamma]d, \[Gamma]\[Phi], \[Gamma]1, \[CapitalOmega]x, \[CapitalGamma]2eff, sol},
   {\[CapitalGamma]CI, \[CapitalGamma]BA, \[CapitalGamma]d, \[Gamma]\[Phi], \[Gamma]1, \[CapitalOmega]x} = rates;
   \[CapitalGamma]2eff = \[Gamma]1/2 + \[Gamma]\[Phi] + \[CapitalGamma]d;
   sol = DSolve[{\[FormalX]'[t] == -\[CapitalGamma]2eff \[FormalX][t],
      \[FormalY]'[t] == -\[CapitalGamma]2eff \[FormalY][t] - \[CapitalOmega]x \[FormalZ][t],
      \[FormalZ]'[t] == \[Gamma]1 (1 - \[FormalZ][t]) + \[CapitalOmega]x \[FormalY][t],
      \[FormalX][0] == x0, \[FormalY][0] == y0, \[FormalZ][0] == z0},
     {\[FormalX], \[FormalY], \[FormalZ]}, t] // First;
   Function[tt, Re[{\[FormalX][tt], \[FormalY][tt], \[FormalZ][tt]} /. sol]]];   (* {<x>,<y>,<z>}(tt) *)
```

## Reproduction on a clean dimensionless case

Rescaling every rate by $\gamma_1$ makes the physics dimensionless without changing it: the trajectory in $\tau=\gamma_1 t$ is identical, and the dimensionless ratios that survive are what matter. Take the decoherence structure of a representative dispersive transmon ($30/20\,\mu s$ for $T_1/T_2$, a weak readout), rescaled by $\gamma_1$, which gives $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1\}/\gamma_1=\{0.352,0,0.176,1,1\}$ with efficiency $\eta=(\Gamma_{CI}+\Gamma_{BA})/(2\Gamma_d)=1$ (optimal homodyne phase). Pair it with a deliberately gentle drive $\Omega_x/\gamma_1=3$, from a ground-state start $\{0,0,1\}$. The drive is gentled on purpose so the explicit scheme sits well inside its stability window; the true transmon drive is about sixty times faster, and that stiff case is the honest catch below. The rate vector is positional, $\{\Gamma_{CI},\Gamma_{BA},\Gamma_d,\gamma_\phi,\gamma_1,\Omega_x\}$:

```wl
ratesB = {0.352, 0., 0.176, 1., 1., 3.}; initB = {0., 0., 1.}; tfB = 3.; nB = 1000;
```

First the exact reference: the affine mean ODE closed in `DSolve` form, evaluated as a function of time (its third component is $\langle z\rangle$):

```wl
lindB = lindbladBloch[ratesB, initB];
```

Next the native route, integrating $n$ trajectories of the `ItoProcess` with the SDE-native scalar-noise stochastic Runge-Kutta:

```wl
itoB = itoEnsemble[ratesB, initB, tfB, tfB/600, nB, 2024];
```

Then the hand-rolled route, the same ensemble through the projected Euler-Maruyama inside `NDSolve`/`WhenEvent`:

```wl
emB = emEnsemble[ratesB, initB, tfB, tfB/600, nB, 2024, True];
```

Compare the exact final $\langle z\rangle$ against the two ensemble means:

```wl
{lindB[tfB][[3]], Last[itoB["meanZ"]], Last[emB["meanZ"]]}
```

Both ensemble means land within the Monte-Carlo scatter $1/\sqrt{n}$ of the exact Lindblad, which is the sense in which they agree. As curves, the two ensemble means track the exact black line:

```wl
ListLinePlot[{Transpose[{itoB["times"], itoB["meanZ"]}],
   Transpose[{emB["times"], emB["meanZ"]}],
   Transpose[{itoB["times"], lindB[#][[3]] & /@ itoB["times"]}]},
  PlotStyle -> {ColorData[97, 1], Directive[Dashed, ColorData[97, 2]], Directive[Thick, Black]},
  PlotLegends -> {"ItoProcess (native)", "NDSolve/WhenEvent EM", "Lindblad (DSolve)"},
  Frame -> True, FrameLabel -> {"t", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
  ImageSize -> 460, PlotLabel -> "EM == Ito == Lindblad"]
```

Now the invariant. Run one more ensemble with the projection turned off, at a coarser step so any overshoot is visible:

```wl
emBnp = emEnsemble[ratesB, initB, tfB, tfB/50, nB, 2024, False];
```

Compare the largest Bloch radius reached along the way. The projected Euler-Maruyama holds the ball to machine precision, the native SRK overshoots it a little (it is not norm-preserving), and the unprojected Euler-Maruyama overshoots more, since off the sphere the noise is no longer tangent and nothing pulls the state back:

```wl
{"EM proj maxr" -> emB["maxr"], "Ito maxr" -> itoB["maxr"], "EM noProj maxr" -> emBnp["maxr"]}
```

The projected value is $1$ to machine precision (the projection holds the ball exactly); the native SRK sits a hair above $1$, its small non-norm-preserving overshoot; and the unprojected value sits clearly above $1$, since without the projection the off-sphere overshoot is never corrected. Here the case is gentle, so the excursion stays bounded; the outright runaway is the stiff regime below, where a coarse step falls below the stability threshold.

## Giving the check teeth: the second moment

Because the mean is diffusion-blind, matching it certifies only the drift. The measurement backaction lives in the *second* moment, the spread of the final $z$ across paths. A small helper reads that spread off an ensemble:

```wl
stdFin[a_] := StandardDeviation[a["zTrajs"][[All, -1]]];
```

To show the spread has teeth, build a deliberately wrong SDE: keep the drift but halve the diffusion, $\vec b\mapsto \vec b/2$. The mean is unchanged (it only sees the drift), so this cannot be caught by $\langle z\rangle$:

```wl
halfProc = ItoProcess[{driftV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, ratesB],
    List /@ (diffV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, ratesB]/2),
    {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}}, {{\[FormalX], \[FormalY], \[FormalZ]}, initB}, t];
```

Integrate this wrong process and keep the final $z$ of each path:

```wl
halfFin = (BlockRandom[SeedRandom[2024];
    RandomFunction[halfProc, {0., tfB, tfB/600}, nB,
      Method -> "StochasticRungeKuttaScalarNoise"]]["ValueList"])[[All, -1, 3]];
```

Compare the standard deviations. The Euler-Maruyama spread matches the native Ito spread, while the half-diffusion spread is clearly smaller, even though its mean still sits on the Lindblad value:

```wl
{"Std Ito" -> stdFin[itoB], "Std EM" -> stdFin[emB],
 "Std half-diffusion" -> StandardDeviation[halfFin],
 "half-diffusion mean (still matches)" -> Mean[halfFin]}
```

The two correct spreads agree, while the half-diffusion spread is markedly smaller (halving the noise amplitude narrows the distribution), even though its mean is unmoved. The standard deviation is what distinguishes a correct measurement model from a wrong one.

## The honest catch: Euler-Maruyama is only strong order 0.5

The catch is that the hand-rolled route is explicit Euler-Maruyama, and it treats the fast Rabi rotation explicitly. Linearizing the drift about a point, the $(y,z)$ block is a damped rotation with eigenvalues $-\tfrac{\Gamma_2^{\mathrm{eff}}+\gamma_1}{2}\pm i\,\Omega_x$ (for $\Omega_x$ large). Explicit Euler is stable only when $|1+\lambda\,\Delta t|\le 1$, i.e. $\Delta t\lesssim 2\,\mathrm{Re}(-\lambda)/|\lambda|^2\approx(\Gamma_2^{\mathrm{eff}}+\gamma_1)/\Omega_x^2$, so the step count must exceed

$$
n^{*}=\frac{\Omega_x^2\,t_f}{\Gamma_2^{\mathrm{eff}}+\gamma_1}.
$$

Turn off the measurement ($\Gamma_{CI}=\Gamma_{BA}=0$, so the diffusion vanishes and the single path is deterministic) and keep only dephasing, relaxation, and the drive, now at the true transmon strength. In the same dimensionless units the real drive-to-relaxation ratio is $\Omega_x/\gamma_1=2\pi\,(1\,\mathrm{MHz})(30\,\mu s)=60\pi\approx 188$, and the final time $\tau_f=\gamma_1 t_f=1/3$ is the dimensionless form of $10\,\mu s$:

```wl
ratesFast = {0., 0., 0., 1., 1., 60 Pi}; initFast = N[{Sin[0.001], 0., Cos[0.001]}]; tfFast = 1./3;
```

The exact answer and the threshold $n^{*}=\Omega_x^2 t_f/(\Gamma_2^{\mathrm{eff}}+\gamma_1)$ are (positions 5, 4, 3 of the rate vector are $\gamma_1,\gamma_\phi,\Gamma_d$):

```wl
{"exact Lindblad z(tf)" -> lindbladBloch[ratesFast, initFast][tfFast][[3]],
 "threshold n*" -> ratesFast[[6]]^2 tfFast/((ratesFast[[5]]/2 + ratesFast[[4]] + ratesFast[[3]]) + ratesFast[[5]])}
```

The exact $\langle z\rangle$ lands partway up the damped-Rabi envelope, and $n^{*}$ is a few thousand steps. Now sweep the step count across the threshold. One row of the sweep runs a single deterministic path at a given step count, with and without projection, and reads off the final $\langle z\rangle$:

```wl
detRow[ns_] := {ns, Last[emEnsemble[ratesFast, initFast, tfFast, tfFast/ns, 1, 1, False]["meanZ"]],
   Last[emEnsemble[ratesFast, initFast, tfFast, tfFast/ns, 1, 1, True]["meanZ"]]};
```

Tabulate the sweep and watch the transition across $n^{*}$:

```wl
Grid[Prepend[detRow /@ {1000, 4000, 8000, 16000, 32000},
   {"nSteps", "EM noProj z(tf)", "EM proj z(tf)"}], Frame -> All]
```

Below $n^{*}$ the explicit scheme amplifies: the unprojected value blows past one, and projection prevents blowup but pins the state at the pole, which is wrong physics. Above $n^{*}$ it converges to the Lindblad, but only at first order in $\Delta t$, so it is still visibly above the exact value even at the finest step counts in the sweep. The native SRK reaches the same answer at a thousand steps, because its stability region is far larger. This is the division of labor: for a stiff, fast-Rabi qubit, the native `ItoProcess` route is the right tool, and the projected `NDSolve/WhenEvent` scheme is a demonstration of what one can build when a native Wiener input is absent, not the tool of choice.

## What is now true

The Cartesian-Bloch trajectory SDE carries its own positivity: the noise is tangent to the sphere and the drift points inward at $r=1$, so the exact flow never leaves the ball, and because that drift is affine the ensemble mean is the Lindblad master equation identically. `RandomFunction` on `ItoProcess` integrates this in the coordinates that keep it stable; the hand-rolled Euler-Maruyama in `NDSolve/WhenEvent` reaches the same ensemble, and the same second moment, but only as a strong-order-0.5 explicit scheme whose stability against the fast Rabi is bought with either a projection or a small step, so its cost tracks $\Omega_x^2$ while the SRK's does not. The projection is what makes the discrete Cartesian scheme stay in the ball that the continuous SDE never left on its own.
