---
Template: Default
---

# Watching a Qubit: Continuous Measurement and Quantum Trajectories

<!-- #| style: Subtitle -->
A computation-first tour of the stochastic master equation, quantum trajectories, and the detector record, built from plain matrices in the Wolfram Language.

<!-- #| style: Author -->
Mads Bahrami (last updated: August 7, 2026)

<!-- #| style: Affiliation -->
Wolfram Research, Inc.

### Setting the Stage: How This Notebook Flows

This notebook is a computation-first tour of what happens when you watch a quantum system instead of leaving it alone: how the smooth Lindblad master equation is really an average over detector records nobody looked at, how conditioning on one record turns that smooth equation into a jagged stochastic one, how a single step of that stochastic equation can be written so the state stays a state, how strong monitoring freezes a qubit while weak monitoring lets it drift, how monitored emission sends a pure state wandering over the Bloch sphere, and how the detector current, all by itself, carries a spectral peak at the qubit's oscillation frequency. Every one of these statements is turned into a short piece of code you can run and change.

I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. Real understanding shows up when you can do the thing, simulate it, predict it, estimate it, rather than repeat words about it. So the rhythm here is always the same: state an idea, say it again in plainer words, compute it, then read the output aloud and let it raise the next question.

The environment you are reading is a live Wolfram notebook. Evaluate the cells from top to bottom; some cells define objects that later cells use, so order matters. The story is meant to be read like a movie, one continuous take with a few headings for breath, and my suggestion is to look at each output and its meaning first, before unpacking every line of the input that produced it. You are not locked into the code as given: change a rate, reseed the noise, monitor a different operator, and watch what moves. That is the whole point of a notebook.

One habit is worth fixing before we start, because it decides what every plot in this notebook is allowed to mean. We will constantly produce three different kinds of object, and they do not have the same standing:

| Object | Where it comes from | What it is |
|---|---|---|
| A detector current | the instrument | observed data |
| A state path computed from a record | record, model, initial state | inferred, not directly observed |
| A current and state path from `simulate` | model plus random numbers | synthetic, useful for insight, not evidence |
| The average of many `simulate` runs matched to the master equation | two calculations of one model | a check on the code, not a test of nature |

In other words, computation plays three roles here: it infers a state from data, it illustrates a model by simulation, and it checks code against an independent calculation. Keeping those apart is what keeps us honest. We set $\hbar=1$ throughout, so every Hamiltonian is measured in inverse time.

### Prerequisites and How to Read This

You should be comfortable with density matrices (a state is a positive, unit-trace, Hermitian matrix $\rho$), with the matrix exponential, and with the idea of a Wiener increment $dW$: a zero-mean Gaussian random kick with variance $dt$, the mathematical version of one tick of white noise. Nothing beyond that is assumed. The pace steepens only in Part VIII and the appendix, where we ask how small the time step must be; you can read the physics parts first and come back to those.

Let's start!

## Part I: Two Equations, One Monitored Qubit

Leave a quantum system connected to a large environment (the modes of the electromagnetic field, a transmission line, a cavity) but look at none of what leaks out. Averaged over everything the environment carries away, the state $\rho$ obeys the Lindblad master equation

$$
\frac{d\rho}{dt}=-i[H,\rho]+\sum_k\mathcal D[c_k]\rho,
\qquad
\mathcal D[c]\rho=c\rho c^\dagger-\tfrac12\{c^\dagger c,\rho\}.
$$

In other words, this is the smooth, deterministic equation you get when you throw the detector record away. It has no random term, because averaging over all possible records has washed the randomness out. Each $\mathcal D[c]$ is a dissipator: it is how one decay or dephasing channel $c$ nudges the state.

A physical state must keep unit trace, so the first thing to check is that the dissipator moves probability around without creating or destroying any. Build $\mathcal D[c]$ from plain matrices and confirm it is trace-free for a completely general channel and state:

```wl
dissipator[channel_, state_] :=
  channel . state . ConjugateTranspose[channel] -
   (ConjugateTranspose[channel] . channel . state +
      state . ConjugateTranspose[channel] . channel)/2;

FullSimplify[Tr[dissipator[{{c11, c12}, {c21, c22}}, {{r11, r12}, {r21, r22}}]]]
```

As one can see, the trace of the dissipator vanishes identically, with no assumption on $c$ or $\rho$. That is exactly why the master equation can be trusted to keep $\operatorname{Tr}\rho=1$: dissipation only reshuffles probability.

Now put a detector on some of those channels. Monitor operators $L_r$ with efficiencies $0\le\eta_r\le1$ ($\eta_r=1$ means every photon is caught, $\eta_r=0$ means none is), and leave the rest, call them $V_j$, unwatched. The record you read is random, so the state conditioned on that record now obeys a stochastic master equation:

$$
d\rho_c=-i[H,\rho_c]\,dt
+\sum_j\mathcal D[V_j]\rho_c\,dt
+\sum_r\mathcal D[L_r]\rho_c\,dt
+\sum_r\sqrt{\eta_r}\,\mathcal H[L_r]\rho_c\,dW_r,
\qquad
\mathcal H[c]\rho=c\rho+\rho c^\dagger-\operatorname{Tr}[(c+c^\dagger)\rho]\,\rho.
$$

Equivalently, the smooth master equation has grown one extra term per monitored channel, each carrying a Wiener increment $dW_r$. That term, $\mathcal H[L_r]$, is the measurement backaction: it is how learning something from the detector kicks the state. What the detector actually hands you each step is a real readout increment

$$
dR_r=\sqrt{\eta_r}\,\operatorname{Tr}[(L_r+L_r^\dagger)\rho_c]\,dt+dW_r .
$$

In words, the readout is a weak signal (proportional to the expectation of $L_r+L_r^\dagger$) buried in a fresh burst of white noise $dW_r$. The current an experimenter plots is $I_r=dR_r/dt$.

The backaction term also has to leave the state normalized: a random kick that changed the trace would not preserve probability. That subtraction $-\operatorname{Tr}[(c+c^\dagger)\rho]\,\rho$ is precisely what keeps it trace-free. Confirm it for a general channel and a unit-trace state:

```wl
measurementBackaction[channel_, state_] :=
  channel . state + state . ConjugateTranspose[channel] -
   Tr[(channel + ConjugateTranspose[channel]) . state] state;

FullSimplify[
 Tr[measurementBackaction[{{c11, c12}, {c21, c22}},
    {{r11, r12}, {Conjugate[r12], 1 - r11}}]], r11 \[Element] Reals]
```

As expected, the trace of the backaction also vanishes. So both new pieces, the deterministic dissipator and the random backaction, conserve probability, and $\rho$ stays a valid state whether or not we keep the record. That is the whole logical content of Part I; the rest of the notebook is about computing with these two equations. Two limits are worth remembering as checks, and we will meet both later: with every $\eta_r=0$ the record carries no information and the conditional state falls back onto the plain master equation, and with every loss channel monitored at $\eta_r=1$ a pure state stays pure.

## Part II: A Step That Stays a State

The conditional equation is a stochastic differential equation, so the obvious way to march it forward is one Euler step: add the drift times $dt$ and the backaction times the noise. Let us try exactly that and watch it break.

First, the qubit toolkit we will use for the rest of the notebook. The lowering operator $J_-=(X-iY)/2$ sends $|0\rangle$ to $|1\rangle$, so in this notebook $|0\rangle$ is the excited state that decays:

```wl
X = PauliMatrix[1]; Y = PauliMatrix[2]; Z = PauliMatrix[3];
identity2 = IdentityMatrix[2];
jMinus = (X - I Y)/2;
densityMatrix[stateVector_] := Outer[Times, stateVector, Conjugate[stateVector]];
blochVector[state_] := Re[Tr[state . #] & /@ {X, Y, Z}];
purity[state_] := Re@Tr[state . state];
ket0 = {1, 0}; ket1 = {0, 1}; ketPlus = {1, 1}/Sqrt[2];
```

`blochVector[state]` returns the three expectation values $\langle X\rangle,\langle Y\rangle,\langle Z\rangle$, and `purity[state]` returns $\operatorname{Tr}(\rho^2)$, which equals one exactly for a pure state and drops below one for a mixed one. When the state is inferred from a record, both of these are inferred too; they are not detector readings.

Now the naive Euler step of the conditional equation, for a single monitored operator, written straight from the two superoperators above:

```wl
naiveStep[hamiltonian_, channel_, dt_][state_, noise_] :=
  state + (-I (hamiltonian . state - state . hamiltonian) +
      dissipator[channel, state]) dt +
   measurementBackaction[channel, state] noise;
```

Take the equal superposition $|+\rangle$, monitor $Z$ with no Hamiltonian, and feed in a single noise kick a little larger than one standard deviation. At $|+\rangle$ the signal $\operatorname{Tr}[(Z+Z^\dagger)\rho]$ vanishes, so this kick is pure noise. Look at the smallest eigenvalue of the result:

```wl
naiveResult = naiveStep[0 identity2, Z, 0.1][densityMatrix[ketPlus], 0.5];
Min@Eigenvalues[naiveResult]
```

The smallest eigenvalue is negative. A negative eigenvalue means a negative probability, so the Euler step has walked the state clean off the set of valid density matrices. The reason is structural: the backaction enters linearly in the noise, and a large enough kick overshoots the Bloch ball. Shrinking $dt$ makes it rarer, but a rare unphysical state is still a broken simulation. This forces a better step.

The fix, due to Rouchon and Ralph ([arXiv:1410.5345](https://arxiv.org/abs/1410.5345)), is to write the update as a sum of terms of the form $A\rho A^\dagger$ and then renormalize:

$$
\rho_{n+1}=\frac{\widetilde\rho_{n+1}}{\operatorname{Tr}\widetilde\rho_{n+1}},
\qquad
\widetilde\rho_{n+1}=M\rho_nM^\dagger
+dt\sum_jV_j\rho_nV_j^\dagger
+dt\sum_r(1-\eta_r)L_r\rho_nL_r^\dagger,
$$

$$
M=1-\Big(iH+\tfrac12\sum_rL_r^\dagger L_r
+\tfrac12\sum_jV_j^\dagger V_j\Big)dt
+S+\tfrac12\Big(S^2-dt\sum_r\eta_rL_r^2\Big),
\qquad
S=\sum_r\sqrt{\eta_r}\,L_r\,dR_r .
$$

The idea is worth saying plainly: if $\rho_n$ is positive, then every $A\rho_nA^\dagger$ is positive, a sum of positive matrices is positive, and dividing by a positive trace keeps it positive. So a valid state stays valid by construction, for any $dt\ge0$ and any $0\le\eta_r\le1$. (The compact $S^2$ form of the last term is an exact rewrite of the Milstein double sum $\tfrac12\sum_{r,s}\sqrt{\eta_r\eta_s}L_rL_s(dR_rdR_s-\delta_{rs}dt)$; we test its consequences below rather than take it on faith.)

Build that step from primitives. It reads the fixed model once, then returns a function that maps a state and a readout increment to the next state:

```wl
rouchonStep[hamiltonian_, monitored_, efficiencies_, unmonitored_, dt_] :=
 With[{identity = IdentityMatrix[Length[hamiltonian]],
   drift = I hamiltonian +
     Total[ConjugateTranspose[#] . # & /@ Join[monitored, unmonitored]]/2,
   itoCorrection =
     Total[MapThread[#1 (#2 . #2) &, {efficiencies, monitored}]]/2},
  Function[{state, readout},
   Module[{signal, m, updated},
    signal = Total@MapThread[Sqrt[#2] #1 #3 &, {readout, efficiencies, monitored}];
    m = identity - drift dt + signal + signal . signal/2 - itoCorrection dt;
    updated = m . state . ConjugateTranspose[m] +
      dt Total[# . state . ConjugateTranspose[#] & /@ unmonitored] +
      dt Total[MapThread[(1 - #2) #1 . state . ConjugateTranspose[#1] &,
         {monitored, efficiencies}]];
    updated/Re[Tr[updated]]]]];
```

Now revisit the exact kick that broke the Euler step, the same $|+\rangle$, the same $dt$, the same noise (which at $|+\rangle$ is the whole readout increment), and read the smallest eigenvalue again:

```wl
rouchonResult =
  rouchonStep[0 identity2, {Z}, {1.}, {}, 0.1][densityMatrix[ketPlus], {0.5}];
{Min@Eigenvalues[rouchonResult], Tr[rouchonResult]}
```

The smallest eigenvalue is now zero to roundoff and the trace is one: the step that overshot the ball now lands on its surface. Positivity is no longer a matter of luck with the step size.

With one physical step in hand, two small functions finish the toolkit. The first generates a synthetic record: at each step it draws a Wiener increment, turns it into a readout $dR$ through the signal-plus-noise rule, and advances the state. It returns the state path, the record, and the noise it drew:

```wl
readoutIncrement[state_, monitored_, efficiencies_, dt_, noise_] :=
  MapThread[
   Sqrt[#3] Re@Tr[(#1 + ConjugateTranspose[#1]) . state] dt + #2 &,
   {monitored, noise, efficiencies}];

simulate[initialState_, hamiltonian_, monitored_, efficiencies_, unmonitored_,
   dt_, finalTime_] := Module[
  {stepCount = Round[finalTime/dt], step, noise, states, record},
  step = rouchonStep[hamiltonian, monitored, efficiencies, unmonitored, dt];
  noise = RandomVariate[NormalDistribution[0, Sqrt[dt]],
    {stepCount, Length[monitored]}];
  states = FoldList[
    Function[{state, dW},
     step[state, readoutIncrement[state, monitored, efficiencies, dt, dW]]],
    initialState, noise];
  record = MapThread[readoutIncrement[#1, monitored, efficiencies, dt, #2] &,
    {Most[states], noise}];
  <|"Times" -> dt Range[0, stepCount], "States" -> states,
   "Record" -> record, "Noise" -> noise|>];

simulate[initialState_, hamiltonian_, monitored_, dt_, finalTime_] :=
  simulate[initialState, hamiltonian, monitored,
   ConstantArray[1., Length[monitored]], {}, dt, finalTime];
```

The second uses a record you supply and never draws a random number. This is the function to point at calibrated experimental data:

```wl
infer[initialState_, hamiltonian_, monitored_, efficiencies_, unmonitored_,
   dt_, record_] :=
  <|"Times" -> dt Range[0, Length[record]],
   "States" -> FoldList[
     rouchonStep[hamiltonian, monitored, efficiencies, unmonitored, dt],
     initialState, record], "Record" -> record|>;

infer[initialState_, hamiltonian_, monitored_, dt_, record_] :=
  infer[initialState, hamiltonian, monitored,
   ConstantArray[1., Length[monitored]], {}, dt, record];
```

The short call form of each function assumes unit efficiency and no unwatched channel; we turn those knobs in Part VII. Notice the deliberate split: `simulate` contains the only random-number generator in the notebook and its outputs are synthetic, while `infer` is deterministic given its record. That split is the honest boundary between a model and a measurement, and Part III makes it concrete. For a hardened, dimension-general version of exactly this integrator with full input checking, see the companion package `manual-implementation-ito-qf-independent.wl` in this folder; here we keep the step bare so the machinery stays visible.

## Part III: One Trajectory, and the Average Is the Master Equation

This is the heart of the notebook, so we build it around one model and look at it from every side. Drive a qubit with $H=\Omega X$, start it in the excited state $|0\rangle$, and monitor its $Z$ component through $\sqrt\gamma\,Z$:

```wl
omega = 1.0; gamma = 0.4;
hamiltonian = omega X;
measuredOperators = {Sqrt[gamma] Z};
initialState = densityMatrix[ket0];
dt = 0.005; finalTime = 10.0;
```

To have an independent standard to compare against, solve the smooth master equation directly with a differential-equation solver. This is a completely separate calculation from the trajectory stepper: no Rouchon step, no random numbers, just `NDSolveValue` on the Lindblad equation:

```wl
masterEquationSolution[initialState_, hamiltonian_, channels_, {t_, t0_, t1_}] :=
  NDSolveValue[{
    state'[t] == -I (hamiltonian . state[t] - state[t] . hamiltonian) +
      Total[dissipator[#, state[t]] & /@ channels],
    state[t0] == initialState}, state, {t, t0, t1}];

masterSolution =
  masterEquationSolution[initialState, hamiltonian, measuredOperators,
   {t, 0, finalTime}];
```

Generate one record and its state path, and plot the monitored component $\langle Z\rangle$ against the smooth master-equation curve:

```wl
oneRun = BlockRandom[SeedRandom[1];
  simulate[initialState, hamiltonian, measuredOperators, dt, finalTime]];
timeGrid = oneRun["Times"];
masterStates = masterSolution /@ timeGrid;

ListLinePlot[{
   Transpose[{timeGrid, blochVector[#][[3]] & /@ oneRun["States"]}],
   Transpose[{timeGrid, blochVector[#][[3]] & /@ masterStates}]},
  PlotStyle -> {Directive[Opacity[0.6], ColorData[97, 1]], Directive[Thick, Red]},
  PlotLegends -> {"one record", "master equation"},
  Frame -> True, FrameLabel -> {"time", "\[LeftAngleBracket]Z\[RightAngleBracket]"},
  GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 540,
  AspectRatio -> 1/2, PlotLabel -> "One jagged record against the smooth average"]
```

As one can see, the single record jitters constantly and never sits on the smooth curve. That is the correct picture, not a numerical defect: the smooth curve is the average over all possible records, and no individual record follows an average, just as no single coin flip equals one half.

Before trusting that jagged path, confirm it is genuinely inferable and not an artifact of the random draw. Feed the record `simulate` produced back into `infer` and compare the two state paths:

```wl
inferred = infer[initialState, hamiltonian, measuredOperators, dt, oneRun["Record"]];
Max[Norm[#1 - #2, "Frobenius"] & @@@
  Transpose[{oneRun["States"], inferred["States"]}]]
```

The two paths coincide exactly, because they fold the same update rule over the same record. This is the sim-to-experiment bridge in miniature: with real data you would drop the record `simulate` invented and hand `infer` the calibrated detector increments instead, and the state path would come out the same way.

Now the central claim. If the smooth curve really is the average over records, then averaging many independent trajectories should reproduce it. Average one hundred and fifty runs and measure the largest gap to the master equation over the whole time interval:

```wl
manyStates = BlockRandom[SeedRandom[42];
  Table[
   simulate[initialState, hamiltonian, measuredOperators, dt, finalTime]["States"],
   {150}]];
meanStates = Mean[manyStates];
Max[Norm[#1 - #2, "Frobenius"] & @@@ Transpose[{meanStates, masterStates}]]
```

The gap is small, and it is Monte Carlo scatter: it shrinks like one over the square root of the number of runs, not because the two calculations are secretly the same code, but because the random backaction has exactly the right mean. Plot the averaged $\langle Z\rangle$ on top of the master-equation curve to see them lie together:

```wl
ListLinePlot[{
   Transpose[{timeGrid, blochVector[#][[3]] & /@ meanStates}],
   Transpose[{timeGrid, blochVector[#][[3]] & /@ masterStates}]},
  PlotStyle -> {Directive[Thick, ColorData[97, 1]], Directive[Dashed, Red]},
  PlotLegends -> {"mean of 150 records", "master equation"},
  Frame -> True, FrameLabel -> {"time", "\[LeftAngleBracket]Z\[RightAngleBracket]"},
  GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 540,
  AspectRatio -> 1/2, PlotLabel -> "Averaging the records rebuilds the smooth equation"]
```

Averaging over records is the whole meaning of the word unravelling: the deterministic master equation is what you get by forgetting which record occurred.

That forgetting has a sharp signature in the purity. With the channel monitored at unit efficiency, every single record keeps the state pure, while the average over records is mixed. Compare the purity of one record with the purity of the averaged state, along the whole path:

```wl
ListLinePlot[{
   Transpose[{timeGrid, purity /@ oneRun["States"]}],
   Transpose[{timeGrid, purity /@ meanStates}]},
  PlotStyle -> {Directive[Thick, ColorData[97, 1]], Directive[Dashed, Red]},
  PlotLegends -> {"one record", "average over records"},
  Frame -> True, FrameLabel -> {"time", "purity  Tr[\[Rho]^2]"},
  GridLines -> Automatic, PlotRange -> {0.45, 1.03}, ImageSize -> 540,
  AspectRatio -> 1/2, PlotLabel -> "A kept record stays pure; the average decoheres"]
```

Here is the payoff of the whole part: decoherence is not something that happens to a monitored qubit. Each record leaves it perfectly pure. Decoherence is what you see when you average over the records you chose not to keep. Ordinary open-system decay is monitored dynamics with the detector switched off.

Before moving on, let us summarize the most important points we have computed so far:

- The dissipator $\mathcal D[c]$ and the backaction $\mathcal H[c]$ are both trace-free, so $\rho$ stays a unit-trace state with or without the record.
- A naive Euler step of the conditional equation can produce a negative eigenvalue; the Rouchon step writes the update as $\sum A\rho A^\dagger$ and stays positive by construction.
- One record gives a jagged $\langle Z\rangle$ that never sits on the smooth master-equation curve.
- Feeding a simulated record into `infer` reproduces the simulated state path exactly.
- The average of many records reproduces the master equation, with only Monte Carlo scatter left over.
- Each unit-efficiency record stays pure; the average is mixed. Decoherence is the average over discarded records.

## Part IV: Strong and Weak Monitoring

What does monitoring strength do to the qubit? For our model, the record-averaged $z=\langle Z\rangle$ obeys a damped-oscillator equation that we can read off and check by hand:

$$
\ddot z+2\gamma\dot z+4\Omega^2z=0,
\qquad z(0)=1,\quad\dot z(0)=0,
\qquad\text{exponents } -\gamma\pm\sqrt{\gamma^2-4\Omega^2}.
$$

In other words, the drive $\Omega$ wants to rotate $Z$ (the $4\Omega^2$ restoring term) and the monitoring $\gamma$ wants to damp it. Before choosing any numbers, verify the equation, both initial conditions, and the strong-monitoring rate symbolically:

```wl
symbolicZ = Exp[-symbolicGamma symbolicTime] (
   Cosh[Sqrt[symbolicGamma^2 - 4 symbolicOmega^2] symbolicTime] +
    symbolicGamma Sinh[Sqrt[symbolicGamma^2 - 4 symbolicOmega^2] symbolicTime]/
     Sqrt[symbolicGamma^2 - 4 symbolicOmega^2]);

FullSimplify[{
   D[symbolicZ, {symbolicTime, 2}] +
    2 symbolicGamma D[symbolicZ, symbolicTime] + 4 symbolicOmega^2 symbolicZ,
   symbolicZ /. symbolicTime -> 0,
   D[symbolicZ, symbolicTime] /. symbolicTime -> 0,
   Limit[symbolicGamma (symbolicGamma -
       Sqrt[symbolicGamma^2 - 4 symbolicOmega^2]), symbolicGamma -> Infinity]},
  Assumptions -> symbolicGamma > 2 symbolicOmega > 0]
```

The four entries come back as the equation residual, the value at zero, the derivative at zero, and the strong-monitoring slow-rate coefficient. The first three confirm $z$ solves the problem; the last says that when $\gamma\gg\Omega$ the slow exponent is $\gamma-\sqrt{\gamma^2-4\Omega^2}\simeq 2\Omega^2/\gamma$. Read that rate plainly: stronger monitoring makes $Z$ change more slowly. Watching the qubit freezes it, an instance of the quantum Zeno effect.

Package the exact solution so we can plot it in either regime:

```wl
zDephasing[t_, omega_, gamma_] := With[{q = gamma^2 - 4 omega^2},
  Which[
   q > 0, Exp[-gamma t] (Cosh[Sqrt[q] t] + gamma Sinh[Sqrt[q] t]/Sqrt[q]),
   q < 0, Exp[-gamma t] (Cos[Sqrt[-q] t] + gamma Sin[Sqrt[-q] t]/Sqrt[-q]),
   True, Exp[-gamma t] (1 + gamma t)]];
```

Choose one strong rate and one weak rate, well separated, and plot the two exact curves:

```wl
omegaD = 1.0; gammaStrong = 20.0; gammaWeak = 0.1;
timeGridD = Range[0., 8.0, 0.0025];

ListLinePlot[{
   Transpose[{timeGridD, zDephasing[timeGridD, omegaD, gammaStrong]}],
   Transpose[{timeGridD, zDephasing[timeGridD, omegaD, gammaWeak]}]},
  PlotStyle -> {Directive[Thick, ColorData[97, 1]], Directive[Thick, ColorData[97, 2]]},
  PlotLegends -> {"strong monitoring: barely moves", "weak monitoring: oscillates"},
  Frame -> True, FrameLabel -> {"time", "\[LeftAngleBracket]Z\[RightAngleBracket]"},
  GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 540,
  AspectRatio -> 1/2, PlotLabel -> "Strong monitoring freezes Z; weak monitoring lets it swing"]
```

As expected, the strongly monitored qubit clings to its initial $\langle Z\rangle=1$ while the weakly monitored one swings freely under the drive. Confirm that the independent master-equation solver agrees with the exact curve in both regimes:

```wl
strongMaster = masterEquationSolution[densityMatrix[ket0], omegaD X,
   {Sqrt[gammaStrong] Z}, {t, 0, 8.0}];
weakMaster = masterEquationSolution[densityMatrix[ket0], omegaD X,
   {Sqrt[gammaWeak] Z}, {t, 0, 8.0}];

{Max@Abs[(blochVector[strongMaster[#]][[3]] & /@ timeGridD) -
     zDephasing[timeGridD, omegaD, gammaStrong]],
 Max@Abs[(blochVector[weakMaster[#]][[3]] & /@ timeGridD) -
     zDephasing[timeGridD, omegaD, gammaWeak]]}
```

Both gaps are at solver tolerance. Now the individual records, which the smooth curves hide. Run one trajectory in each regime and plot the inferred $\langle Z\rangle$:

```wl
strongRun = BlockRandom[SeedRandom[7];
  simulate[densityMatrix[ket0], omegaD X, {Sqrt[gammaStrong] Z}, 0.0025, 8.0]];
weakRun = BlockRandom[SeedRandom[7];
  simulate[densityMatrix[ket0], omegaD X, {Sqrt[gammaWeak] Z}, 0.0025, 8.0]];

ListLinePlot[{
   Transpose[{strongRun["Times"], blochVector[#][[3]] & /@ strongRun["States"]}],
   Transpose[{weakRun["Times"], blochVector[#][[3]] & /@ weakRun["States"]}]},
  PlotStyle -> {ColorData[97, 1], ColorData[97, 2]},
  PlotLegends -> {"strong: one record", "weak: one record"},
  Frame -> True, FrameLabel -> {"time", "inferred \[LeftAngleBracket]Z\[RightAngleBracket]"},
  GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 540,
  AspectRatio -> 1/2, PlotLabel -> "Individual records: pinned near an eigenstate versus free to drift"]
```

The strongly monitored record stays pinned near a $Z$ eigenstate for long stretches, with occasional abrupt hops: the continuous limit of quantum jumps. The weakly monitored record wanders with the drive. These paths illustrate the model; only a detector record from a real device could support a claim about a real device.

## Part V: Monitored Spontaneous Emission

So far we monitored a Hermitian operator, $Z$. Monitoring emission is different, because the jump operator $J_-$ is not Hermitian. In homodyne detection the instrument beats the emitted light against a reference and records one real quadrature of the field. With a complete, unit-efficiency record and a pure initial state, the inferred state stays pure, so it lives on the surface of the Bloch sphere and we can watch it move there.

Set up a driven, detuned qubit whose emission $\sqrt{\gamma_S}\,J_-$ is monitored by homodyne detection, starting from $|+\rangle$:

```wl
omegaS = 2.0; detuningS = 1.0;
hamiltonianS = (omegaS X + detuningS Z)/2;
gammaS = 0.2; measuredS = {Sqrt[gammaS] jMinus};
initialStateS = densityMatrix[ketPlus];
dtS = 0.01; finalTimeS = 10.0;
```

Generate five independent records and turn each state path into a curve of Bloch vectors:

```wl
emissionRuns = BlockRandom[SeedRandom[3];
  Table[simulate[initialStateS, hamiltonianS, measuredS, dtS, finalTimeS], {5}]];
emissionBlochPaths = (blochVector /@ #["States"]) & /@ emissionRuns;
```

Draw a bare Bloch sphere and overlay the five inferred paths:

```wl
blochSphere = Graphics3D[{
    {Opacity[0.12], LightBlue, Sphere[]},
    {Gray, Line[{{-1.1, 0, 0}, {1.1, 0, 0}}],
     Line[{{0, -1.1, 0}, {0, 1.1, 0}}], Line[{{0, 0, -1.1}, {0, 0, 1.1}}]}},
   Boxed -> False, PlotRange -> 1.15, ImageSize -> 420];

Show[blochSphere,
 Graphics3D[{Thick,
   MapIndexed[{ColorData[97, First@#2], Line@#1} &, emissionBlochPaths]}]]
```

Each curve is one possible inferred state path for one possible record, wandering on the surface as the detector delivers its noisy quadrature. Now watch what averaging does: run many records and compare the mean state with the master equation:

```wl
emissionMaster = masterEquationSolution[initialStateS, hamiltonianS, measuredS,
   {t, 0, finalTimeS}];
emissionTimeGrid = emissionRuns[[1]]["Times"];
emissionMeanStates = Mean@BlockRandom[SeedRandom[11];
    Table[simulate[initialStateS, hamiltonianS, measuredS, dtS, finalTimeS]["States"],
     {150}]];

Max[Norm[#1 - #2, "Frobenius"] & @@@
  Transpose[{emissionMeanStates, emissionMaster /@ emissionTimeGrid}]]
```

The mean tracks the master equation to Monte Carlo scatter. So the surface-wandering pure records, averaged together, sink into the ball and decay toward the ground state: ordinary spontaneous decay is the averaged, unwatched version of these jittering pure trajectories, exactly the lesson of Part III now carried over from a Hermitian channel to an emission channel.

## Part VI: What the Detector Record Shows

The state paths we have plotted are inferred objects. The detector current is not: it is the raw thing the instrument produces, and its statistics are observable without any model. From one long current we can compute its mean, how readings separated in time move together, and how its power is spread over frequency. The last two are the autocovariance and the power spectrum.

Take a driven, emitting qubit and record for a long time so the statistics settle:

```wl
omegaC = 3.5; gammaC = 0.7;
hamiltonianC = omegaC Y/2;
measuredC = {Sqrt[gammaC] jMinus};
dtC = 0.01; finalTimeC = 1200.0;

longRun = BlockRandom[SeedRandom[1];
  simulate[densityMatrix[ket0], hamiltonianC, measuredC, dtC, finalTimeC]];
simulatedCurrent = longRun["Record"][[All, 1]]/dtC;
```

A single current sample is almost all noise. Confirm that its standard deviation is set by the white-noise level $1/\sqrt{dt}$, not by the signal:

```wl
{Mean[simulatedCurrent], StandardDeviation[simulatedCurrent], 1/Sqrt[dtC]}
```

The mean is small and the spread matches $1/\sqrt{dt}$: one reading tells you next to nothing. The coherent motion is hiding in the correlations between readings. The autocovariance measures how readings separated by a delay move together; compute it with a fast Fourier transform, dividing each delay by its own number of data pairs:

```wl
autocovariance[data_, maximumLag_Integer, dt_, stride_Integer : 1] := Module[
  {centered = N[data - Mean[data]], count = Length[data], paddedLength, sums, lags},
  lags = Range[0, Min[maximumLag, count - 1], stride];
  paddedLength = 2^Ceiling[Log[2, 2 count - 1]];
  sums = Re@InverseFourier[
     Abs[Fourier[PadRight[centered, paddedLength],
        FourierParameters -> {1, -1}]]^2, FourierParameters -> {1, -1}];
  Transpose[{dt lags, sums[[lags + 1]]/(count - lags)}]];

currentCovariance = Rest@autocovariance[simulatedCurrent, Floor[8/dtC], dtC, 4];

ListLinePlot[currentCovariance, Frame -> True, ImageSize -> 540,
 AspectRatio -> 1/2, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"lag", "current autocovariance"},
 PlotLabel -> "A driven oscillation survives inside the noisy current"]
```

We drop the zero-lag point because it holds the large white-noise variance; away from it, the autocovariance oscillates at the driven frequency. So the coherent motion that no single reading revealed is plainly there once we correlate readings across time. The same information, sorted by frequency instead of by delay, is the power spectrum. Estimate it by splitting the record into segments, removing each segment's mean, and averaging the segment spectra:

```wl
averagedPeriodogram[data_, dt_, segmentLength_Integer] := Module[
  {segments = Partition[N[data], segmentLength, segmentLength], twoSided,
   meanSpectrum, indices, frequencies},
  twoSided = (dt/segmentLength) Abs[
       Fourier[# - Mean[#], FourierParameters -> {1, -1}]]^2 & /@ segments;
  meanSpectrum = Mean[twoSided];
  indices = Range[2, Ceiling[segmentLength/2]];
  frequencies = 2 Pi Range[1, Length[indices]]/(segmentLength dt);
  Transpose[{frequencies, 2 meanSpectrum[[indices]]}]];

spectrum = Select[averagedPeriodogram[simulatedCurrent, dtC, 4096], First[#] <= 8 &];
oscillationFrequency = Sqrt[omegaC^2 - gammaC^2/16];

ListLinePlot[spectrum, Frame -> True, ImageSize -> 540, AspectRatio -> 1/2,
 PlotRange -> All, GridLines -> {{oscillationFrequency}, {2}},
 FrameLabel -> {"angular frequency", "one-sided power spectrum"},
 PlotLabel -> "A peak rises above the flat white-noise floor"]
```

Under this one-sided convention ideal white noise sits at a flat level of two, and above it a peak rises in the driven band. Where should it sit? The record-averaged coherence oscillates at $\omega_{\mathrm{osc}}=\sqrt{\Omega_C^2-\gamma_C^2/16}$ and decays at rate $3\gamma_C/4$, so the line is centered at $\omega_{\mathrm{osc}}$ with a half width of $3\gamma_C/4$. A single noisy periodogram cannot pin a broad line to one bin, so instead of reading the argmax, compare the mean power in that predicted band with a control band away from it, against the white-noise floor of two:

```wl
halfWidth = 3 gammaC/4;
bandPower[spectrum_, center_, width_] :=
  Mean@Select[spectrum, center - width <= First[#] <= center + width &][[All, 2]];

{bandPower[spectrum, oscillationFrequency, halfWidth], bandPower[spectrum, 1.0, halfWidth]}
```

The predicted band carries power well above two while the control band sits at the floor: the coherent line lives exactly where the model puts it. That is the honest reach of this section. The mean, the autocovariance, and the excess power are genuine observed features of the current, computed from the record alone, while turning that line into a precise frequency or a decay rate needs a longer record or a model fit, not a width read off by eye.

## Part VII: Efficiency and Dimension

Two knobs on the toolkit have been sitting at their defaults. Detection efficiency $\eta$ changes the record-based state estimate, though not the master equation you get after averaging. Turn it below one and the conditional state can no longer stay pure: you are catching only a fraction of the emitted information, and the fraction you miss mixes your estimate. Compare one record at full efficiency with one at half efficiency, driven along the same noise seed:

```wl
omegaE = 2.0; gammaE = 1.0;
hamiltonianE = omegaE X/2;
measuredE = {Sqrt[gammaE] jMinus};
initialStateE = densityMatrix[ketPlus];
dtE = 0.01; finalTimeE = 6.0;

fullEfficiency = BlockRandom[SeedRandom[5];
  simulate[initialStateE, hamiltonianE, measuredE, {1.0}, {}, dtE, finalTimeE]];
halfEfficiency = BlockRandom[SeedRandom[5];
  simulate[initialStateE, hamiltonianE, measuredE, {0.5}, {}, dtE, finalTimeE]];

{Max[Abs[purity[#] - 1] & /@ fullEfficiency["States"]],
 Mean[purity /@ halfEfficiency["States"]]}
```

The full-efficiency record stays pure to roundoff; the half-efficiency record is mixed, its mean purity well below one. In short, efficiency is how much of the environment's information you actually hold, and holding less of it leaves you with a more uncertain, more mixed estimate.

Now confirm that averaging over records erases $\eta$ entirely. Take three efficiencies, and for each average the final state over many records and compare with the single master-equation solution:

```wl
efficiencyMaster = masterEquationSolution[initialStateE, hamiltonianE, measuredE,
   {t, 0, finalTimeE}];

AssociationThread[{1.0, 0.5, 0.0},
 MapThread[
  Function[{efficiency, seed},
   Norm[Mean@BlockRandom[SeedRandom[seed];
        Table[simulate[initialStateE, hamiltonianE, measuredE, {efficiency}, {},
           dtE, finalTimeE]["States"][[-1]], {100}]] -
      efficiencyMaster[finalTimeE], "Frobenius"]],
  {{1.0, 0.5, 0.0}, {101, 102, 103}}]]
```

All three efficiencies land on the same master-equation state within Monte Carlo scatter. So $\eta$ reshapes every individual record but leaves their average untouched: the unwatched master equation cannot tell how efficiently you were watching.

The second knob is dimension. Nothing in `rouchonStep` assumed a qubit; it read the dimension from the Hamiltonian. Monitor a spin-one system through its lowering operator $J_-$ and drive it with $J_x=(J_-+J_+)/2$:

```wl
jMinus3 = {{0, 0, 0}, {Sqrt[2], 0, 0}, {0, Sqrt[2], 0}};
jX3 = (jMinus3 + ConjugateTranspose[jMinus3])/2;
hamiltonian3 = 1.5 jX3;
measured3 = {Sqrt[0.6] jMinus3};
initialState3 = densityMatrix[{1, 0, 0}];
dt3 = 0.01; finalTime3 = 8.0;
```

Run one qutrit record at efficiency $\eta=0.4$ and check it stays a physical state by reading its smallest eigenvalue along the path:

```wl
qutritRun = BlockRandom[SeedRandom[2];
  simulate[initialState3, hamiltonian3, measured3, {0.4}, {}, dt3, finalTime3]];
Min[Flatten[Eigenvalues[(# + ConjugateTranspose[#])/2] & /@ qutritRun["States"]]]
```

The smallest eigenvalue never goes negative: the positivity argument never mentioned the dimension, so it holds here unchanged. Finally, average the full qutrit density matrix over many records and compare with its master equation:

```wl
qutritMaster = masterEquationSolution[initialState3, hamiltonian3, measured3,
   {t, 0, finalTime3}];
qutritMean = Mean@BlockRandom[SeedRandom[9];
    Table[simulate[initialState3, hamiltonian3, measured3, {0.4}, {}, dt3,
       finalTime3]["States"], {120}]];

Max[Norm[#1 - #2, "Frobenius"] & @@@
  Transpose[{qutritMean, qutritMaster /@ qutritRun["Times"]}]]
```

The full three-by-three mean matches the master equation, populations and coherences alike. Monitoring, unravelling, and the master-equation-as-average picture were never facts about qubits; they are facts about open quantum systems of any size.

## Part VIII: Can We Trust This Step Size?

Every plot above used a fixed $dt$, and the Rouchon step earns its keep by guaranteeing a valid state at any $dt$. But a valid state is not the same as an accurate one. In other words, positivity is a physicality check, not an error estimate: a positive density matrix can still be the wrong positive density matrix if the step is too coarse. So the step size needs a separate test, and this section runs the two cheapest ones. (The stricter pathwise and detector-statistics tests live in the appendix.)

The first test isolates the deterministic error. Set the efficiency to zero: the channel still dephases the qubit, but no information enters, so every record is identical and the trajectory reduces to the master equation with no random scatter. Define the trace distance between two states, halve the step repeatedly, and watch the distance from the master-equation answer:

```wl
traceDistance[firstState_, secondState_] :=
  Total[SingularValueList[N[firstState - secondState]]]/2;

observedOrders[errors_Association] := AssociationThread[Most@Keys[errors],
   Log[2, Most[Values[errors]]/Rest[Values[errors]]]];

deterministicSteps = {0.08, 0.04, 0.02, 0.01};
deterministicErrors = AssociationThread[deterministicSteps,
   traceDistance[
      Last@simulate[initialState, hamiltonian, measuredOperators, {0.}, {},
         #, finalTime]["States"],
      masterSolution[finalTime]] & /@ deterministicSteps];

<|"errors" -> deterministicErrors, "orders" -> observedOrders[deterministicErrors]|>
```

As one can see, the error falls in a regular way and the estimated order sits near one: halving the step roughly halves the error, the mark of a first-order scheme. This test matters most when the Hamiltonian is fast, but it says nothing about the random part.

The second test is the one we already passed in Part III, stated as a criterion: the mean over many records must approach the master equation, and the gap must be consistent with Monte Carlo scatter rather than a fixed bias. We saw the gap shrink toward the solver curve; increasing the number of records shrinks it further. Positivity, deterministic convergence, and the correct ensemble mean answer three different questions, and a step you trust has to pass all three. There is no universal good value of $dt$: a step is adequate only relative to a stated model, time interval, detector bandwidth, and accuracy goal.

## Appendix: Certifying a Step Size

The two checks in Part VIII are enough for the physics above. Certifying a step size for research needs two stricter tools, gathered here so the main story stays uncluttered.

The first is pathwise convergence. Two independently sampled trajectories represent different detector histories and must not be compared point by point; to refine honestly, sample the finest Wiener increments once and build every coarser record by summing adjacent fine increments. This shares a single underlying noise path across resolutions:

```wl
coarsenIncrements[increments_, blockSize_] :=
  Total /@ Partition[N[increments], blockSize];

coupledFamily[initialState_, hamiltonian_, monitored_, efficiencies_, finestStep_,
   finalTime_, blockSizes_, seed_] := Module[
  {stepCount = Round[finalTime/finestStep], fineSteps, fineNoise},
  fineSteps = ConstantArray[finestStep, stepCount];
  fineNoise = BlockRandom[SeedRandom[seed];
    RandomVariate[NormalDistribution[0, Sqrt[finestStep]],
     {stepCount, Length[monitored]}]];
  AssociationThread[N[finestStep blockSizes],
   Function[block,
     With[{steps = coarsenIncrements[fineSteps, block],
       noise = coarsenIncrements[fineNoise, block]},
      Module[{step = rouchonStep[hamiltonian, monitored, efficiencies, {}, steps[[1]]],
        states},
       states = FoldList[
         Function[{state, dW},
          step[state, readoutIncrement[state, monitored, efficiencies, steps[[1]], dW]]],
         initialState, noise];
       <|"States" -> states, "Noise" -> noise|>]]] /@ blockSizes]];

adjacentErrors[family_Association] := Module[
  {stepValues = Keys[family], trajectories = Values[family], ratios},
  ratios = Round[Most[stepValues]/Rest[stepValues]];
  AssociationThread[Most[stepValues],
   MapThread[
    Function[{pair, ratio},
     Max@MapThread[traceDistance,
       {pair[[1]]["States"], pair[[2]]["States"][[1 ;; ;; ratio]]}]],
    {Partition[trajectories, 2, 1], ratios}]]];

rootMeanSquare[errorList_] := With[{keys = Keys@First[errorList]},
   AssociationThread[keys, Sqrt[Mean[(Values /@ errorList)^2]]]];
```

Run four resolutions across sixty-four independent noise paths, sharing noise within each family, and estimate the order from the root-mean-square error over paths:

```wl
families = coupledFamily[initialState, hamiltonian, measuredOperators, {1.0},
     0.0025, 2.0, {8, 4, 2, 1}, #] & /@ Range[101, 164];
maximumErrors = rootMeanSquare[adjacentErrors /@ families];

<|"errors" -> maximumErrors, "orders" -> observedOrders[maximumErrors]|>
```

The error decreases through the refinements and the order is close to one, so for this single monitored channel the scheme is first order pathwise, not merely positive.

The second tool checks the detector record against the state that generated it. For a simulated record, subtract the predicted signal and divide by $\sqrt{dt}$; this residual of the record against its own prediction (its innovation, in the language of filtering) should reproduce the standard-normal samples the simulator drew,

$$
z_{r,n}=\frac{dR_{r,n}-\sqrt{\eta_r}\operatorname{Tr}[(L_r+L_r^\dagger)\rho_n]\,dt}{\sqrt{dt}} .
$$

Form these residuals for a simulated record and read their mean, variance, and lag-one correlation:

```wl
normalizedInnovations[run_Association, monitored_, efficiencies_] := MapThread[
   Function[{state, record, dt},
    (record - readoutIncrement[state, monitored, efficiencies, dt,
        ConstantArray[0., Length[monitored]]])/Sqrt[dt]],
   {Most[run["States"]], run["Record"], Differences[run["Times"]]}];

innovationRun = BlockRandom[SeedRandom[55];
  simulate[initialState, hamiltonian, measuredOperators, 0.0025, 5.0]];
innovations = Flatten@normalizedInnovations[innovationRun, measuredOperators, {1.0}];

<|"mean" -> Mean[innovations], "variance" -> Variance[innovations],
 "lag-one correlation" -> Correlation[Most[innovations], Rest[innovations]]|>
```

The residuals have mean near zero, variance near one, and negligible lag-one correlation: the record generator draws fresh, independent noise, exactly as the model demands. This check validates the record and its pairing with the pre-step state; it does not replace the pathwise test.

A last caution, because it is a real trap. The compact Rouchon step was built from a first-order argument stated for commuting monitored operators, and $X$ and $Z$ do not commute. The step stays positive by construction here, as Part II guaranteed for any channels, but that positivity does not certify the convergence order. Monitor $X$ and $Z$ at once and measure the order instead of assuming it:

```wl
ncMeasured = {Sqrt[0.2] X, Sqrt[0.3] Z};
ncFamilies = coupledFamily[densityMatrix[ketPlus], 0.3 Y, ncMeasured, {0.8, 0.5},
     0.0025, 1.0, {8, 4, 2, 1}, #] & /@ Range[301, 332];
ncErrors = rootMeanSquare[adjacentErrors /@ ncFamilies];

<|"errors" -> ncErrors, "orders" -> observedOrders[ncErrors]|>
```

The trace distances still fall, but the measured order is visibly below one. The lesson is not to copy a convergence order from one model to another: when the monitored operators fail to commute, the step still gives a physical state, but you must measure its accuracy for yourself.

## Where This Leaves Us

We built, from plain matrices, a positivity-preserving step for the conditional master equation, a `simulate` that manufactures synthetic records and their state paths and an `infer` that reconstructs a state path from a record it did not generate, an independent master-equation solver to check against, and a small kit for reading a detector current. With those, we watched a single record jitter while its average rebuilt the smooth master equation, saw strong monitoring freeze a qubit and weak monitoring free it, sent monitored emission wandering over the Bloch sphere, pulled a spectral peak out of a noisy current, dialed efficiency down to mix the conditional state, and ran the whole story unchanged in three dimensions.

The natural next steps live alongside this notebook. To do something with the record rather than only observe it, the sibling notebook on record analysis and feedback estimates the state in real time and feeds it back. To leave the discrete Hilbert space for a particle on a line, the continuous-space notebooks replace the qubit with a monitored position. And to close the loop with a real experiment, replace the synthetic readout here with calibrated detector increments in `infer`, and reserve part of the data for an independent test: continuous records have been used to estimate mechanical quantum states in real time and checked against later or separate measurements ([Magrini et al., Nature 595, 373 (2021)](https://www.nature.com/articles/s41586-021-03602-3); [Rossi et al., Phys. Rev. Lett. 123, 163601 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.163601)). A synthetic trajectory in this notebook does not inherit that status; it earns only the right to be checked, which is what we did throughout.

### Acknowledgment

The state update follows the positivity-preserving scheme of Rouchon and Ralph. The examples use standard models of monitored $Z$ noise and homodyne-detected emission. If you have suggestions or questions, reach out at [quantum@wolfram.com](mailto:quantum@wolfram.com).
