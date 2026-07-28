---
Template: TechNote
Name: TimeEvolution
Title: Time Evolution
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/TimeEvolution
Keywords: [QuantumEvolve, Schrodinger equation, Lindblad, Liouvillian, Kossakowski, open system, propagator, Rabi oscillation, decoherence, Bloch vector, Heisenberg picture, superradiance]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [GettingStarted, QuantumMachineLearning]
---

[QuantumEvolve]() solves the equation of motion of a quantum system. Which equation that is depends on what you hand it: the Schrodinger equation for a closed system, the Lindblad master equation once jump operators are present, the Heisenberg equation when the thing being evolved is an operator, and the Kossakowski equation when the jump operators are correlated.

The call is the same shape in every case - a Hamiltonian, optionally some dissipation, an initial condition, and a time - and the result is always an object that depends on time. Supply a time to it and you get an ordinary state or operator back, which answers every property the framework knows about.

Two choices control everything else. A bare symbol for the time solves the equation *exactly*, through [DSolve](); a `{t, t0, t1}` interval solves it *numerically*, through [NDSolve](). And an initial condition of [None]() evolves the propagator instead of any particular state, which can then be applied to whatever state you like.

```wl
<< Wolfram`QuantumFramework`
```

## A Driven Qubit

A [PauliX]() Hamiltonian drives a qubit between $|0\rangle$ and $|1\rangle$. Evolving $|0\rangle$ over one full period:

```wl
psi = QuantumEvolve[QuantumOperator["PauliX"], QuantumState["0"], {t, 0, 2 Pi}]
```

<!-- => the QuantumState summary box: a pure state of dimension 2, carrying the parameter t -->

---

The result is a state as a function of time, and `t` is recorded on it as a parameter:

```wl
psi["Parameters"]
```

<!-- => {t} -->

---

Supplying a time gives the state at that time - an ordinary [QuantumState]() with nothing left symbolic:

```wl
psi[1.]
```

<!-- => the QuantumState summary box: a pure state of dimension 2, no parameters -->

---

From there every state property is available as usual:

```wl
psi[1.]["ProbabilitiesList"]
```

<!-- => {0.2919266077143294, 0.7080733922856707} -->

---

Sampling the interval shows the population sloshing between the two levels and back:

```wl
Table[Round[psi[t]["ProbabilitiesList"], 0.0001], {t, 0, 3}]
```

<!-- => {{1., 0.}, {0.2919, 0.7081}, {0.1732, 0.8268}, {0.9801, 0.0199}} -->

---

Plotting it gives the Rabi oscillation:

```wl
Plot[psi[t]["ProbabilitiesList"][[2]], {t, 0, 2 Pi},
    AxesLabel -> {"t", "population of |1\[RightAngleBracket]"}]
```

<!-- => the plot: a sin^2 curve rising from 0 to 1 and back over the interval -->

## What the Solution Actually Is

A numerically evolved state is not a table of samples and it is not a symbolic expression. Its amplitudes *are* the solver's own [InterpolatingFunction](), held as the state's array and applied to the time parameter:

```wl
psi["State"]
```

<!-- => the InterpolatingFunction summary box, Domain {{0., 6.283185307179586}} and
   Output dimensions {2}, applied to t -->

---

That single object carries the whole two-component solution, so the state still knows its own shape without interpolating anything:

```wl
psi["Dimensions"]
```

<!-- => {2} -->

---

Nothing is computed until a time is supplied; `psi[1.]` is what triggers an interpolation. This matters for large systems, where the alternative - one scalar interpolation per amplitude, each carrying its own copy of the solver's time grid - costs a copy of that grid per dimension. Asking for it explicitly with `"Expand" -> True` shows what is being avoided:

```wl
expanded = QuantumEvolve[QuantumOperator["PauliX"], QuantumState["0"], {t, 0, 2 Pi}, "Expand" -> True];
Normal[expanded["State"]]
```

<!-- => a list of two separate InterpolatingFunction summary boxes, each applied to t -->

---

The two forms give the same numbers:

```wl
expanded[1.]["ProbabilitiesList"]
```

<!-- => {0.29192656493030605, 0.7080734350696939} -->

## Exact Solutions

Passing a bare symbol instead of an interval solves the equation symbolically rather than numerically:

```wl
exact = QuantumEvolve[QuantumOperator["PauliX"], QuantumState["0"], t];
Simplify[Normal[exact["StateVector"]]]
```

<!-- => {Cos[t], -I Sin[t]} -->

---

That is the textbook solution, and its populations are what the plot above traced:

```wl
Simplify[exact["ProbabilitiesList"], Assumptions -> t \[Element] Reals]
```

<!-- => {Cos[t]^2, Sin[t]^2} -->

---

The initial state may be symbolic too, in which case so is the answer:

```wl
Simplify[Normal[QuantumEvolve[QuantumOperator["X"], QuantumState[{\[Alpha], \[Beta]}], t]["StateVector"]]]
```

<!-- => {α Cos[t] - I β Sin[t], β Cos[t] - I α Sin[t]} -->

---

Omit the initial state and the time entirely and the register state $|0\rangle$ is used, with the formal symbol $t$ as the time parameter. This is the shortest thing that says anything:

```wl
QuantumEvolve[QuantumOperator["X"]]
```

<!-- => the QuantumState summary box: Cos[t] |0> - I Sin[t] |1>, in the formal symbol t -->

---

The numerical solution agrees with the closed form to the solver's tolerance:

```wl
Max[Abs[psi[1.]["ProbabilitiesList"] - {Cos[1.] ^ 2, Sin[1.] ^ 2}]]
```

<!-- => 2.5987900498236627*^-8 -->

## Time-Dependent Hamiltonians

Nothing requires the Hamiltonian to be constant. A symbol in it that matches the evolution parameter simply makes the generator time dependent. This is a spin with splitting $\omega_0$ driven transversally at the same frequency - a resonant drive:

```wl
\[Omega]0 = 2 Pi;
\[Omega]p = Pi/10;
tf = 2 Pi/\[Omega]p;
driven = \[Omega]0 QuantumOperator["JZ"] + 2 \[Omega]p Cos[\[Omega]0 t] QuantumOperator["JX"];
```

---

With no initial state given, the register state is used, so the Hamiltonian and the time interval are the whole call:

```wl
rabi = QuantumEvolve[driven, {t, 0, tf}]
```

<!-- => the QuantumState summary box: a pure state of dimension 2, parameter t -->

---

The Bloch vector traces the slow Rabi nutation underneath the fast Larmor precession. Over the interval $2\pi/\omega_p$ it makes one full circuit from the north pole and back:

```wl
Plot[Evaluate[Re[rabi[t]["BlochVector"]]], {t, 0, tf}, PlotRange -> All,
    PlotLegends -> {Subscript[r, x], Subscript[r, y], Subscript[r, z]}]
```

<!-- => the plot: three components, z sweeping from +1 down to -1 and back, x and y small -->

---

Halfway through, the spin has been driven to the south pole - a $\pi$ pulse:

```wl
Round[Chop[Re[rabi[tf/2]["BlochVector"]]], 0.0001]
```

<!-- => {-0.025, -0.0002, -0.9997} -->

---

The evolution is unitary throughout, so the state stays pure and the Bloch vector stays on the sphere:

```wl
rabi[tf/2]["Purity"]
```

<!-- => 1 -->

## Observables and Trajectories

An evolved state answers every ordinary state property once a time is given, so an observable is measured as it would be on any state:

```wl
QuantumMeasurement[QuantumMeasurementOperator[QuantumOperator["PauliZ"]] @ psi[1.]]["Mean"]
```

<!-- => -0.4161467845713413 -->

---

For a qubit the Bloch vector gives the motion directly. Under [PauliX]() the state precesses about the $x$ axis, so the $x$ component stays zero:

```wl
Chop[psi[1.]["BlochCartesianCoordinates"]]
```

<!-- => {0, -0.9092974506127978, -0.4161467845713415} -->

---

Drawing that as a curve on the Bloch sphere shows the precession as the great circle it is:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"],
    ParametricPlot3D[Evaluate[Re[psi[t]["BlochVector"]]], {t, 0, 2 Pi}]]
```

<!-- => the Bloch sphere with a great circle in the y-z plane traced on it -->

## Open Systems: the Lindblad Equation

Jump operators turn the Schrodinger equation into the Lindblad master equation,

$$\partial_t \rho = -i[H, \rho] + \sum_i \gamma_i \left( L_i \rho L_i^\dagger - \tfrac{1}{2} \{ L_i^\dagger L_i, \rho \} \right),$$

with $\rho$ the density matrix, $L_i$ the jump operators and $\gamma_i$ their rates. Passing them as `Ls -> rates` is enough to switch equations. Here a [PauliZ]() Hamiltonian with [PauliX]() dephasing at rate 0.6:

```wl
noisy = QuantumEvolve[QuantumOperator["PauliZ"], {QuantumOperator["PauliX"]} -> {0.6},
    QuantumState["0"], {t, 0, 3}]
```

<!-- => the QuantumState summary box: a mixed state of dimension 2, parameter t -->

---

The state is a density matrix now rather than a vector, because dissipation does not preserve purity. Purity falls from 1 toward the maximally mixed value of 1/2:

```wl
Table[noisy[t]["Purity"], {t, 0, 3}]
```

<!-- => {1., 0.5453589757953925, 0.5041148713352459, 0.5003732925729159} -->

---

The populations equalize as it goes:

```wl
noisy[3.]["ProbabilitiesList"]
```

<!-- => {0.5136618551616506, 0.48633814483834936} -->

## Relaxing to a Thermal State

Dephasing drives a qubit to the maximally mixed state because it has no preferred direction. Amplitude damping does have one, and at finite temperature the two directions compete. With a mean occupation $n$ of the bath, decay happens at rate $\gamma(n+1)$ and absorption at rate $\gamma n$:

```wl
\[CapitalOmega] = 50; \[Gamma] = 1; n = 3;
H = \[CapitalOmega]/2 QuantumOperator["Z"];
\[Sigma]m = QuantumOperator[{{0, 1}, {0, 0}}];
Ls = {\[Sigma]m, \[Sigma]m["Dagger"]};
\[Gamma]s = \[Gamma] {n + 1, n};
```

---

Starting from a pure state tilted off the axis:

```wl
\[Rho]0 = QuantumState[{Cos[Pi/8], Exp[I Pi/4] Sin[Pi/8]}]
```

<!-- => the QuantumState summary box: a pure state of dimension 2 -->

---

```wl
\[Rho]t = QuantumEvolve[H, Ls -> \[Gamma]s, \[Rho]0, {t, 0, 1}]
```

<!-- => the QuantumState summary box: a mixed state of dimension 2, parameter t -->

---

By $t = 1$ the populations have all but reached the thermal ratio $\{(n+1)/(2n+1),\, n/(2n+1)\}$, which for $n = 3$ is $\{4/7, 3/7\}$:

```wl
\[Rho]t[1.]["ProbabilitiesList"]
```

<!-- => {0.5716858359636007, 0.4283141640363992} -->

---

```wl
N[{n + 1, n}/(2 n + 1)]
```

<!-- => {0.5714285714285714, 0.42857142857142855} -->

---

The trajectory spirals in from the surface of the Bloch sphere toward that point, the fast [PauliZ]() precession winding around the slow decay:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"],
    ParametricPlot3D[Evaluate[Re[\[Rho]t[t]["BlochVector"]]], {t, 0, 1}]]
```

<!-- => the Bloch sphere with an inward spiral traced on it -->

---

The rates can equally be folded into the operators as $\sqrt{\gamma_i} L_i$, in which case the jump operators are given as a plain list and the result is identical:

```wl
\[Rho]t2 = QuantumEvolve[H, Sqrt[\[Gamma]s] Ls, \[Rho]0, {t, 0, 1}];
Max[Abs[Normal[\[Rho]t2[1.]["DensityMatrix"]] - Normal[\[Rho]t[1.]["DensityMatrix"]]]]
```

<!-- => 0. -->

## The Liouvillian and the Propagator

An open-system generator can be assembled on its own as a Liouvillian superoperator $\mathcal{L}$, the object with $\mathrm{d}\rho/\mathrm{d}t = \mathcal{L}\rho$:

```wl
\[ScriptCapitalL] = QuantumOperator["Liouvillian"[QuantumOperator["Z"], {\[Sigma]m}, {0.3}]]
```

<!-- => the QuantumOperator summary box: a superoperator on one qubit -->

---

It acts on a density matrix flattened into a vector, so for a qubit its matrix is 4x4 rather than 2x2:

```wl
Dimensions[\[ScriptCapitalL]["Matrix"]]
```

<!-- => {4, 4} -->

---

Passing [None]() in place of an initial state evolves the *propagator* rather than any particular state. Given a Liouvillian this produces $e^{\mathcal{L}t}$, the map that takes any initial density matrix to its state at time $t$:

```wl
U = QuantumEvolve[I \[ScriptCapitalL], None, {t, 0, 5}]
```

<!-- => the QuantumOperator summary box: a superoperator on one qubit, parameter t -->

---

Because it is the map and not a trajectory, it can be applied to any state without solving again:

```wl
\[Rho]start = QuantumState[{Cos[Pi/8], Sin[Pi/8]}];
U[2.][\[Rho]start]["ProbabilitiesList"]
```

<!-- => {0.9196283974343614, 0.08037160256563854} -->

---

Solving the same problem directly as a Lindblad equation gives the same answer, to the solver's tolerance:

```wl
QuantumEvolve[QuantumOperator["Z"], {\[Sigma]m} -> {0.3}, \[Rho]start, {t, 0, 5}][2.]["ProbabilitiesList"]
```

<!-- => {0.919628396623777, 0.08037160337622301} -->

---

The propagator also exists for closed systems, where it is the familiar unitary. Solved exactly it is a closed form:

```wl
Simplify[Normal[QuantumEvolve[QuantumOperator["PauliX"], None, t]["Matrix"]]]
```

<!-- => {{Cos[t], -I Sin[t]}, {-I Sin[t], Cos[t]}} -->

## The Heisenberg Picture

The same superoperator run backwards evolves *observables* instead of states. Taking its adjoint and applying it to an operator - passed through [QuantumState]() as a matrix state, then read back as an operator - gives $A(t)$:

```wl
A0 = QuantumOperator["Z"];
At = U["Dagger"][A0["MatrixQuantumState"]]["Operator"];
```

---

The two pictures are two ways of writing the same number, $\operatorname{Tr}[A_0 \rho(t)] = \operatorname{Tr}[A(t) \rho_0]$. In the Schrodinger picture the state moves:

```wl
Re[Tr[Normal[A0["Matrix"]] . Normal[U[2.][\[Rho]start]["DensityMatrix"]]]]
```

<!-- => 0.8392567948687228 -->

---

In the Heisenberg picture the observable moves and the state does not:

```wl
Re[Tr[Normal[At[2.]["Matrix"]] . Normal[\[Rho]start["DensityMatrix"]]]]
```

<!-- => 0.839256794868723 -->

---

They agree to machine precision across the interval, since both read the same solved propagator:

```wl
Max[Table[Abs[Tr[Normal[A0["Matrix"]] . Normal[U[t][\[Rho]start]["DensityMatrix"]]] -
    Tr[Normal[At[t]["Matrix"]] . Normal[\[Rho]start["DensityMatrix"]]]], {t, 0, 5, 0.5}]]
```

<!-- => 1.1102230246251565*^-16 -->

## The Kossakowski Equation

The Lindblad equation assumes the jump operators are uncorrelated. Dropping that assumption gives the Kossakowski equation,

$$\partial_t \rho = -i[H, \rho] + \sum_{ij} \beta_{ij} \left( L_i \rho L_j^\dagger - \tfrac{1}{2} \{ L_j^\dagger L_i, \rho \} \right),$$

where the rates become a positive-semidefinite matrix $\beta_{ij}$ whose off-diagonal entries encode correlations between channels. A diagonal $\beta$ recovers the Lindblad equation. Building the Liouvillian with a matrix in place of a list of rates is all it takes.

Two atoms decaying into a shared field is the standard example. The Hamiltonian is just their bare energies:

```wl
\[Omega]1 = 0.9; \[Omega]2 = 0.7;
Hk = QuantumOperator[\[Omega]1/2 "ZI" + \[Omega]2/2 "IZ"];
```

---

Each atom decays through its own lowering operator:

```wl
jumpOps = QuantumTensorProduct[QuantumOperator[#], QuantumOperator[#2]] & @@@ {{"-", "I"}, {"I", "-"}}
```

<!-- => {QuantumOperator summary box, QuantumOperator summary box}: two 2-qubit operators -->

---

The correlation between the two channels is the off-diagonal $\gamma_{12}$; the matrix is positive semidefinite as long as $\gamma \geq \gamma_{12}$:

```wl
\[Gamma]d = 0.5; \[Gamma]12 = 0.45;
\[Beta] = {{\[Gamma]d, \[Gamma]12}, {\[Gamma]12, \[Gamma]d}};
\[ScriptCapitalL]k = QuantumOperator["Liouvillian"[Hk, jumpOps, \[Beta]]];
Dimensions[\[ScriptCapitalL]k["Matrix"]]
```

<!-- => {16, 16} -->

---

Evolve three initial states with one excitation between them - a product state, the singlet, and the triplet. When the Hamiltonian is already a superoperator, the initial state comes directly after it:

```wl
states = QuantumState /@ {"01", "PsiMinus", "PsiPlus"};
evolved = QuantumEvolve[I \[ScriptCapitalL]k, #, {t, 0, 5}] & /@ states;
```

---

Count the excitations left at each time with the total-number operator:

```wl
totalExcitation = Total @ With[{mat = {{1, 0}, {0, 0}}}, Array[QuantumOperator[mat, {#}] &, 2]];
excitations = Table[Re[Tr[Normal[#[t]["DensityMatrix"]] . Normal[totalExcitation["Matrix"]]]],
    {t, 0, 5, 0.05}] & /@ evolved;
```

---

All three start with one excitation and all three decay, but at three different rates. The correlation splits the symmetric and antisymmetric combinations apart: the singlet is subradiant, decaying at $\gamma - \gamma_{12}$, and the triplet is superradiant, decaying at $\gamma + \gamma_{12}$ - here a factor of nineteen between them:

```wl
ListLinePlot[excitations, DataRange -> {0, 5}, PlotRange -> {0, 1.05},
    PlotLegends -> {"|01\[RightAngleBracket]", "singlet", "triplet"},
    AxesLabel -> {"t", "excitations"}]
```

<!-- => the plot: three curves from 1, the triplet falling fastest and the singlet slowest -->

---

At the end of the interval they are far apart:

```wl
Round[excitations[[All, -1]], 0.0001]
```

<!-- => {0.3705, 0.7082, 0.0327} -->

## Watching for Events

`"AdditionalEquations"` injects extra equations or specifications into the solve, which is how hybrid and piecewise systems are expressed - and how events are detected. Here a strong pulsed drive with a constant detuning:

```wl
rabiDrive = 10. Sin[Pi t] QuantumOperator["X"] - QuantumState["1"]["Operator"];
pulsed = QuantumEvolve[rabiDrive, {t, 0, 2}];
Plot[pulsed[t]["Probability"][[1]], {t, 0, 2}, PlotRange -> {0, 1}]
```

<!-- => the plot: the ground-state population oscillating rapidly within a pulsed envelope -->

---

To find the times at which the population crosses one half, solve again with a [WhenEvent](). The dependent variable of the internal differential equation is the formal symbol $s$, so the condition is written on `Indexed[\[FormalS][t], 1]`, and [Sow]()/[Reap]() collect the crossings:

```wl
{sol, points} = Reap @ QuantumEvolve[rabiDrive, {t, 0, 5},
    "AdditionalEquations" -> WhenEvent[Abs[Indexed[\[FormalS][t], 1]] ^ 2 == 0.5, Sow[t]]];
Length[First[points]]
```

<!-- => 20 -->

---

```wl
Round[Take[First[points], 4], 0.0001]
```

<!-- => {0.2287, 0.4156, 0.5747, 0.7573} -->

---

The crossings really are crossings - evaluating the returned solution at them gives one half back:

```wl
Max[Abs[(sol[#]["Probability"][[1]] & /@ First[points]) - 0.5]]
```

<!-- => 1.783993921589122*^-6 -->

---

Shading the plot between them marks out the intervals where the qubit is mostly excited:

```wl
Plot[sol[t]["Probability"][[1]], {t, 0, 2}, Mesh -> First[points],
    MeshShading -> {Green, Orange}, PlotRange -> {0, 1}]
```

<!-- => the plot: the same curve, alternately shaded green and orange between crossings -->

## The Equations Themselves

`"ReturnEquations" -> True` stops short of solving and hands back the system instead, in exactly the form the solver would have been called with - equations, dependent variable, and time specification:

```wl
QuantumEvolve[QuantumOperator["XX"], {t, 0, 1}, "ReturnEquations" -> True]
```

<!-- => {{Derivative[1][\[FormalS]][t] == SparseArray[<4>] . \[FormalS][t],
        \[FormalS][0] == SparseArray[<1>]}, \[FormalS], {t, 0, 1}} -->

---

The generator is right there in the equation, and it is exactly $-iH$ for the Hamiltonian that went in - here $-i\,X \otimes X$:

```wl
eqs = QuantumEvolve[QuantumOperator["XX"], {t, 0, 1}, "ReturnEquations" -> True];
MatrixForm[Normal[eqs[[1, 1, 2, 1]]]]
```

<!-- => the 4x4 matrix with -I on the antidiagonal and 0 elsewhere -->

---

Reaching for this is worthwhile when the system needs to be handed to a solver directly - to couple it to equations of a different kind, to pick a [Method]() by hand, or simply to see what is being integrated.

The triple is in [NDSolveValue]()'s argument order, but two things in it are for the framework's benefit rather than the solver's: the matrices are sparse, which [NDSolveValue]() will not accept as an initial condition, and the dependent variable is a formal symbol. Make both ordinary and it runs:

```wl
solved = NDSolveValue @@ (eqs /. {sa_SparseArray :> Normal[sa], \[FormalS] -> s})
```

<!-- => the InterpolatingFunction summary box, Domain {{0., 1.}} and Output dimensions {4} -->

---

What comes back is the same solution [QuantumEvolve]() would have handed over wrapped in a state - $|00\rangle$ rotating into $|11\rangle$ as $\cos t$ and $-i \sin t$:

```wl
Chop[solved[1.]]
```

<!-- => {0.5403023235460563, 0, 0, 0. - 0.8414709538478929 I} -->

## Details

- The time argument is either a symbol, which solves the equation exactly with [DSolve](), or `{t, t0, t1}`, which solves it numerically with [NDSolve](). Options not recognized by [QuantumEvolve]() are passed through to whichever solver runs, so `AccuracyGoal`, `Method`, `MaxSteps` and the rest all work.
- The initial condition selects the equation. A state gives the Schrodinger or Lindblad equation; an operator gives the Heisenberg equation; [None]() gives the propagator. Omitting it uses the register state $|0\ldots0\rangle$, and omitting the time along with it uses the formal symbol $t$.
- Jump operators are given either as `{L1, L2, ...} -> {\[Gamma]1, \[Gamma]2, ...}` or as a plain list with the rates folded in as $\sqrt{\gamma_i} L_i$. Correlated channels need the Kossakowski form, which is reached by building `` QuantumOperator["Liouvillian"[H, Ls, \[Beta]]] `` with a rate matrix and evolving `I \[ScriptCapitalL]`.
- With a Hamiltonian that is already a superoperator, the state argument comes directly after it: `` QuantumEvolve[\[ScriptCapitalL], state, {t, 0, 2}] ``.
- The evolution parameter is recorded on the returned object, so `state["Parameters"]` names it and `state[t]` supplies it.
- A numerically evolved state holds the solver's [InterpolatingFunction]() as its amplitudes, so it answers questions about its shape without evaluating anything, and interpolation happens once a time is supplied. `"Expand" -> True` instead stores one scalar interpolation per amplitude, each carrying its own copy of the solver's time grid.
- `"MergeInterpolatingFunctions"` controls what happens when the Hamiltonian itself contains interpolating functions of time. The default merges them into a single sparse interpolation, which is substantially faster than letting each enter the dynamical equation separately.

## Possible Issues

Ask for amplitudes before supplying a time and the answer is a function of the parameter rather than numbers, which is rarely what is wanted from a numerical solve. Supply the time first:

```wl
psi[0.5]["ProbabilitiesList"]
```

<!-- => {0.7701511674960476, 0.2298488325039524} -->

---

The numerical and exact routes agree to the solver's tolerance, not to machine precision, and the difference above is around the eighth decimal. Tighten it with the usual [NDSolve]() options.

---

Defining a time-dependent Hamiltonian as a function of its time argument does not work the way it looks like it should. The operator's matrix is a [SparseArray](), and a pattern variable is not substituted into it, so the symbol in the stored expression is never replaced:

```wl
badH[x_] = Sin[x] QuantumOperator["X"];
Normal[badH[0.3]["Matrix"]]
```

<!-- => {{0, Sin[x]}, {Sin[x], 0}}: the argument was ignored -->

---

Write the Hamiltonian as a plain expression in the evolution parameter instead, and let [QuantumEvolve]() handle the substitution:

```wl
Normal[(Sin[t] QuantumOperator["X"])["Matrix"]]
```

<!-- => {{0, Sin[t]}, {Sin[t], 0}} -->

---

If the exact solver cannot find a closed form, the output is the unevaluated solver call it gave up on, along with a message. Here a Rabi problem with a general detuning, which has no elementary solution:

```wl
QuantumEvolve[1/2 Subscript[\[Omega], 0] QuantumOperator["Z"] +
    \[Gamma]r Cos[Subscript[\[Omega], p] t] QuantumOperator["X"], t]
```

<!-- => QuantumEvolve::error, then the unevaluated DSolveValue[...] call -->

---

Numbers in the Hamiltonian make the same problem solvable numerically, so the fallback when a closed form does not exist is to give the parameters values and pass an interval instead.
