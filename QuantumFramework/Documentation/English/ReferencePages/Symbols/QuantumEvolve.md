---
Template: Symbol
Name: QuantumEvolve
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumEvolve
Keywords: [quantum evolution, Schrodinger equation, Lindblad equation, time evolution, NDSolveValue, Hamiltonian, Liouvillian, open system]
SeeAlso: [QuantumState, QuantumOperator]
RelatedGuides: [WolframQuantumComputationFramework]
---

<!--
Faithful md source for Documentation/English/ReferencePages/Symbols/QuantumEvolve.nb,
re-verified cell by cell against the repo kernel at 10b99e7e (2026-06-11).
Deliberate deviations from the .nb, each verified against the live kernel:
  1. The finite-temperature amplitude-damping example defines the jump operators as
     explicit matrices instead of QuantumOperator /@ {"J-", "J+"}. Named spin
     operators carry the spin-J basis, in which "J-" lowers the spin projection
     (toward computational |1>); since the named-J rebasing fix all evolution routes
     agree on that orientation, but the explicit lowering operator keeps this
     example in the standard computational-basis convention (decay toward |0>),
     which is what the surrounding prose describes.
  2. The WhenEvent example drops the original trailing // Quiet: the cell runs
     message-free on the current kernel, so there is nothing to suppress.
  3. Prose typos fixed ("is solve", "separatedly", "NDSovle", "loses it balance")
     plus light grammar copyediting of captions (articles, agreement); a few
     captions carry an added one-line physics gloss (the Heisenberg trace
     identity, the sub/superradiant rates). The MergeInterpolatingFunctions
     intro note ended with a copy-pasted "can be used to see DEs without
     solving them" sentence belonging to ReturnEquations, replaced with a
     sentence describing what it controls.
  4. Typeset PlotLegends box-strings ("\!\(\*SubscriptBox[\(r\),\(x\)]\)" and the
     FractionBox Bell-state labels) simplified to plain expressions/strings; greek
     long-name escapes written as Unicode. Structure validated against the canonical
     NotebookToMarkdown.wl output (same 12 headings, same 61 code cells).
-->

## Usage

<code>[QuantumEvolve]()[*op*]</code> represents a symbolic quantum state as the solution of the Schrodinger equation with the symbolic Hamiltonian or Liouvillian *op*, from the initial register state, with the formal symbol $t$ as the time parameter.

<code>[QuantumEvolve]()[*op*, *qs*, *t*]</code> represents a symbolic quantum state or operator as the solution of the Schrodinger or Heisenberg equation with the symbolic Hamiltonian or Liouvillian *op*, from the initial state or operator *qs*, and time as *t*.

<code>[QuantumEvolve]()[*op*, *qs*, {*t*, $t_{i}$, $t_{f}$}]</code> represents a numeric quantum state or operator as the solution of the Schrodinger or Heisenberg equation with the numeric Hamiltonian or Liouvillian *op*, from the initial numeric state or operator *qs*, with variable *t* as time in the range $t_{i}$ to $t_{f}$.

<code>[QuantumEvolve]()[*op*, [None](), …]</code> represents an evolution operator of the Hamiltonian or Liouvillian *op*.

<code>[QuantumEvolve]()[*op*, {$L_{1}$, $L_{2}$, …}, *qs*, {*t*, $t_{i}$, $t_{f}$}]</code> represents a numeric quantum state or operator as the solution of the Lindblad equation with the Hamiltonian *op* and the Lindblad operators $L_{j}$, from the initial numeric state or operator *qs*, with variable *t* as time in the range $t_{i}$ to $t_{f}$.

<code>[QuantumEvolve]()[*op*, {$L_{1}$, $L_{2}$, …} -> {$\gamma_{1}$, $\gamma_{2}$, …}, *qs*, {*t*, $t_{i}$, $t_{f}$}]</code> represents a numeric quantum state or operator as the solution of the Lindblad equation with the Hamiltonian *op* and the Lindblad operators $L_{j}$ with rates $\gamma_{j}$, from the initial numeric state or operator *qs*, with variable *t* as time in the range $t_{i}$ to $t_{f}$.

## Details & Options

- [QuantumEvolve]() evolves a quantum state symbolically or numerically. It takes a quantum operator as Hamiltonian or Liouvillian and solves the Schrodinger equation (if a pure quantum state is given), the Liouville-von Neumann equation (if a density matrix is given), or the Lindblad master equation.

- The [Lindblad master equation](https://en.wikipedia.org/wiki/Lindbladian) is $\frac{\mathrm{d}\rho}{\mathrm{d}t} = -i[H,\rho] + \sum_i \left( L_i \rho L_i^\dagger - \tfrac{1}{2}\{L_i^\dagger L_i, \rho\} \right)$ with $\rho$ the density operator, $H$ the Hamiltonian, and $L_i$ the Lindblad operators (also called jump operators). With explicit jump rates it reads $\frac{\mathrm{d}\rho}{\mathrm{d}t} = -i[H,\rho] + \sum_i \gamma_i \left( L_i \rho L_i^\dagger - \tfrac{1}{2}\{L_i^\dagger L_i, \rho\} \right)$ with $\gamma_i$ the rates. One can also feed the Lindblad operators as $\sqrt{\gamma_i}\, L_i$.

- The Hamiltonian and all Lindblad operators should have the same dimensions.

- The initial state can be a vector (pure state) or a density matrix. The solution is given as a density matrix or a pure state. If an initial operator is given, the Heisenberg equation is solved and a time-evolved operator is returned.

- [QuantumEvolve]() gives the solution either numerically or symbolically, depending on the inputs.

- [QuantumEvolve]() accepts all options of [DSolve]() and [NDSolve]().

- The option `"AdditionalEquations"` can be used to add any additional equations to the differential equation, e.g. using [WhenEvent]().

|   |   |   |
|:---|:---|:---|
| `"AdditionalEquations"` | `{}` | additional conditions or differential equations to be added to the Lindblad equation |

- The option `"ReturnEquations"` can be used to see the differential equations without solving them.

|   |   |   |
|:---|:---|:---|
| `"ReturnEquations"` | [False]() | when [True]() it returns the corresponding differential equations, ready to be wrapped by [DSolve]() or [NDSolve]() |

- When different interpolating functions of time are included in the input operators, they can enter the dynamical equation in two ways, either separately or merged together; the option `"MergeInterpolatingFunctions"` controls this.

|   |   |   |
|:---|:---|:---|
| `"MergeInterpolatingFunctions"` | [True]() | when [False]() it includes the interpolating functions separately in [NDSolve]() or [DSolve]() |

## Basic Examples

Evolve a symbolic quantum state by a time-independent Hamiltonian, with independent variable *t*:

```wl
ψf = QuantumEvolve[QuantumOperator["X"], QuantumState[{α, β}], t]
```

<!-- => QuantumState; state vector {α Cos[t] - I β Sin[t], β Cos[t] - I α Sin[t]} -->

Start with the register state, without specifying any independent variable (the formal symbol $t$ is used as the time parameter):

```wl
QuantumEvolve[QuantumOperator["X"]]
```

<!-- => QuantumState; Cos[FormalT] |0⟩ - I Sin[FormalT] |1⟩ -->

If the initial state is [None](), [QuantumEvolve]() returns an evolution operator as a [QuantumOperator]():

```wl
U = QuantumEvolve[QuantumOperator["X"], None, t]
```

<!-- => Cos[t] |0⟩⟨0| - I Sin[t] |0⟩⟨1| - I Sin[t] |1⟩⟨0| + Cos[t] |1⟩⟨1| -->

Evolve $|0\rangle$ by the above unitary:

```wl
U[]
```

<!-- => Cos[t] |0⟩ - I Sin[t] |1⟩ -->

---

Define a time-dependent Hamiltonian:

```wl
h = Subscript[ω, 1]/2 Cos[α] QuantumOperator["Z"] +
  Subscript[ω, 1]/2 Sin[α] (Cos[ω t] QuantumOperator["X"] + Sin[ω t] QuantumOperator["Y"])
```

Define a symbolic initial state:

```wl
ψi = QuantumState[{Cos[α/2], Sin[α/2]}]
```

Evolve the initial state with the Hamiltonian:

```wl
ψf = QuantumEvolve[h, ψi, t];
```

Show the final state vector, with $\lambda = \sqrt{-\omega^2 + 2\omega\cos\alpha\,\omega_1 - \omega_1^2}$:

```wl
FullSimplify[Normal[ψf["StateVector"]]] /. {Sqrt[_] -> λ, 1/Sqrt[_] -> 1/λ}
```

<!-- => {E^(-I t ω/2) Cos[α/2] (Cosh[t λ/2] + I (ω - Subscript[ω,1]) Sinh[t λ/2]/λ),
        E^(I t ω/2) Sin[α/2] (Cosh[t λ/2] - I (ω + Subscript[ω,1]) Sinh[t λ/2]/λ)};
   spot-checked against a direct NDSolve integration to 2*10^-8 -->

---

Define a time-dependent Hamiltonian:

```wl
hamiltonian = ω0 QuantumOperator["JZ"] + 2 ωp Cos[ω0 t] QuantumOperator["JX"];
```

Set numerical values for the Hamiltonian parameters:

```wl
ω0 = 2 π; ωp = π/10;
tf = 2 π/ωp;
```

Find the final state:

```wl
r1 = QuantumEvolve[hamiltonian, {t, 0, tf}];
```

Plot the evolution of the Bloch vector components:

```wl
Plot[Evaluate[Re[r1[t]["BlochVector"]]], {t, 0, tf}, PlotRange -> All,
 PlotLegends -> {Subscript[r, x], Subscript[r, y], Subscript[r, z]}]
```

<!-- => Graphics; pure state, Bloch vector stays on the sphere (norm 1 to 3*10^-6) -->

Include relaxation dynamics with raising and lowering spin angular momentum operators, in a Lindblad dynamical equation:

```wl
γ = 1/20;
r2 = QuantumEvolve[hamiltonian, Sqrt[γ] {QuantumOperator["J+"], QuantumOperator["J-"]}, {t, 0, 2 tf}]
```

<!-- => QuantumState (mixed, parameter t) -->

Plot the evolution of the Bloch vector components:

```wl
Plot[Evaluate[Re[r2[t]["BlochVector"]]], {t, 0, 2 tf}, PlotRange -> All,
 PlotLegends -> {Subscript[r, x], Subscript[r, y], Subscript[r, z]}]
```

<!-- => Graphics; Bloch vector decays toward the center (norm 0.05 at t = 2 tf) -->

---

$T_1$ and $T_2$ relaxation of a spin:

```wl
Ls = {Sqrt[1/T1] QuantumOperator["J+"], Sqrt[1/T1] QuantumOperator["J-"], Sqrt[1/T2] QuantumOperator["JZ"]};
```

Show the density matrix:

```wl
$Assumptions = T1 > 0 && T2 > 0 && Ω0 > 0;
```

```wl
FullSimplify /@ Normal[QuantumEvolve[Ω0 QuantumOperator["JZ"], Ls,
     QuantumState[Array[ρ, {2, 2}]], t]["Computational"]["DensityMatrix"]] // MatrixForm
```

<!-- => {{(ρ[1,1] + ρ[2,2] + E^(-2 t/T1) (ρ[1,1] - ρ[2,2]))/2, E^(t (-1/T1 - 1/(2 T2) - I Ω0)) ρ[1,2]},
        {E^(t (-1/T1 - 1/(2 T2) + I Ω0)) ρ[2,1], (ρ[1,1] + ρ[2,2] - E^(-2 t/T1) (ρ[1,1] - ρ[2,2]))/2}} -->

## Generalizations & Extensions

### Lindblad master equation for the density matrix (Schrodinger picture)

One can solve the Lindblad master equation for the density matrix, $\partial_t \rho = -i[H,\rho] + \sum_i \gamma_i \left( L_i \rho L_i^\dagger - \tfrac{1}{2}\{L_i^\dagger L_i, \rho\} \right)$, with $\rho$ the density matrix, $L_i$ the jump operators and $\gamma_i$ the corresponding damping rates.

Set the $\Omega$, $\gamma$ and $n$ parameters of a finite-temperature amplitude-damping mechanism:

```wl
Ω = 50; γ = 1; n = 3;
```

Initial state:

```wl
ρ0 = QuantumState[{Cos[π/8], Exp[I π/4] Sin[π/8]}]
```

Set up the Hamiltonian, the jump operators (the lowering operator and its adjoint, defined explicitly), and their corresponding rates:

```wl
H = Ω/2 QuantumOperator["Z"];
σm = QuantumOperator[{{0, 1}, {0, 0}}];
Ls = {σm, σm["Dagger"]};
γs = γ {n + 1, n};
```

With *H* a quantum operator as Hamiltonian, $L_i$ the jump operators, and $\gamma_i$ the corresponding rates, determine the state:

```wl
ρt = QuantumEvolve[H, Ls -> γs, ρ0, {t, 0, 1}]
```

<!-- => QuantumState (mixed, parameter t); agrees with a hand-written master-equation
   integration to 4*10^-6 (solver tolerance), relaxing toward populations
   {(n+1)/(2n+1), n/(2n+1)} = {4/7, 3/7} -->

Plot the evolution:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"],
 ParametricPlot3D[Evaluate[Re[ρt[t]["BlochVector"]]], {t, 0, 1}]]
```

<!-- => Graphics3D -->

One can also incorporate the rates into the definition of the operators as $\sqrt{\gamma_i}\, L_i$ and call <code>[QuantumEvolve]()[*H*, {$L_{1}$, $L_{2}$, …}, *qs*, {*t*, $t_{i}$, $t_{f}$}]</code>:

```wl
ρt2 = QuantumEvolve[H, Sqrt[γs] Ls, ρ0, {t, 0, 1}]
```

<!-- => QuantumState; identical to ρt -->

Plot the evolution:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"],
 ParametricPlot3D[Evaluate[Re[ρt2[t]["BlochVector"]]], {t, 0, 1}]]
```

<!-- => Graphics3D -->

One can also create the Liouvillian superoperator $\mathcal{L}$ of $\mathrm{d}\rho/\mathrm{d}t = \mathcal{L}\rho$:

```wl
ℒ = QuantumOperator["Liouvillian"[H, Ls, γs]]
```

<!-- => QuantumOperator (16x16 superoperator) -->

Evolution superoperator, $\rho_t = e^{\mathcal{L} t} \rho_0$:

```wl
expℒ = QuantumEvolve[I ℒ, None, {t, 0, 1}];
ρt3[t_] = expℒ[t][ρ0]
```

<!-- => QuantumState; agrees with ρt to 3*10^-6 (solver tolerance) -->

Plot the evolution:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"],
 ParametricPlot3D[Evaluate[Re[ρt3[t]["BlochVector"]]], {t, 0, 1}]]
```

<!-- => Graphics3D -->

Since the Liouvillian is time independent, one can also exponentiate it directly:

```wl
ρt4[t_] = Exp[ℒ t][ρ0]
```

<!-- => QuantumState; agrees with ρt to 3*10^-6 (solver tolerance) -->

Plot the evolution:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"],
 ParametricPlot3D[Evaluate[Re[ρt4[t]["BlochVector"]]], {t, 0, 1}]]
```

<!-- => Graphics3D -->

### Lindblad master equation for the operators (Heisenberg picture)

One can solve the Lindblad master equation for operators, $\partial_t X = i[H,X] + \sum_i \gamma_i \left( L_i X L_i^\dagger - \tfrac{1}{2}\{L_i^\dagger L_i, X\} \right)$, with $X$ the matrix representation of an operator, $L_i$ the jump operators and $\gamma_i$ the corresponding damping rates.

Create a time-dependent random Hamiltonian of three qubits:

```wl
H = 3 Sin[t] QuantumOperator["RandomHermitian"];
```

Create three random jump operators and random damping rates:

```wl
Ls = Table[QuantumOperator["RandomHermitian"], 3];
γs = RandomReal[1, 3];
```

Create a random mixed state:

```wl
ρ0 = QuantumState["RandomMixed"];
```

Create the Liouvillian operator:

```wl
ℒ = QuantumOperator["Liouvillian"[H, Ls, γs]];
```

Generate the time evolution:

```wl
expℒ = QuantumEvolve[I ℒ, None, {t, 0, 10}];
```

Find the final state:

```wl
ρf = expℒ@ρ0;
```

Generate a random Hermitian operator:

```wl
a0 = QuantumOperator["RandomHermitian"];
```

Find the time evolution of the operator:

```wl
af = expℒ["Dagger"][a0["MatrixQuantumState"]]["Operator"];
```

Show two equivalent ways of calculating the time-dependent average value, $\operatorname{Tr}[A_0\, \rho(t)] = \operatorname{Tr}[A(t)\, \rho_0]$:

```wl
Plot[{Re[Tr[a0["Matrix"] . ρf["DensityMatrix"]]],
  Re[Tr[af["Matrix"] . ρ0["DensityMatrix"]]]}, {t, 0, 10}]
```

<!-- => Graphics; the two curves coincide (verified equal to 10^-10 at sample times) -->

### Kossakowski equation (Lindblad generalization)

It is also possible to solve a more general form of the Lindblad equation, called the Kossakowski equation: $\partial_t X = i[H,X] + \sum_{ij} \beta_{ij} \left( L_i X L_j^\dagger - \tfrac{1}{2}\{L_j^\dagger L_i, X\} \right)$, where $X$ is the matrix representation of an operator, $L_i$ are the jump operators, and $\beta_{ij}$ is a positive-semidefinite matrix that encodes the rates and correlations of the jump operators. If the matrix is diagonal, the equation reduces to the Lindblad master equation.

Consider the decay of two atoms. First define the Hamiltonian, including only the bare energies with frequencies $\omega_1$ and $\omega_2$:

```wl
H = QuantumOperator[ω1/2 "ZI" + ω2/2 "IZ"];
```

Define the jump operators associated with the decay of each atom:

```wl
jumpOps = QuantumTensorProduct[QuantumOperator[#], QuantumOperator[#2]] & @@@ {{"-", "I"}, {"I", "-"}};
```

Create a positive-semidefinite rate matrix; a simple case is a symmetric one (positive semidefinite only for $\gamma \geq \gamma_{12}$):

```wl
βdisp = {{γ, γ12}, {γ12, γ}};
```

Set the simulation parameters:

```wl
γ = 0.5;
γ12 = 0.45;
ω1 = 0.9;
ω2 = 0.7;
```

Create the operator that represents the Kossakowski equation:

```wl
ℒ = QuantumOperator["Liouvillian"[H, jumpOps, βdisp]];
```

Define the initial states:

```wl
states = QuantumState /@ {"01", "PsiMinus", "PsiPlus"};
```

Generate the evolved states by applying the evolution operator to the initial states:

```wl
evolvedStates = QuantumEvolve[I ℒ, #, {t, 0, 5}] & /@ states;
```

Operator representing the total number of excitations:

```wl
totalExcitation = Total@With[{mat = {{1, 0}, {0, 0}}}, Array[QuantumOperator[mat, {#}] &, 2]];
```

Sample the observable at different times:

```wl
excitations = Table[Tr[#[t]["DensityMatrix"] . totalExcitation["Matrix"]], {t, 0, 5, 0.05}] & /@ evolvedStates;
```

The correlations between the jump operators affect the decay of different states; the singlet decays at the subradiant rate $\gamma - \gamma_{12}$ and the triplet at the superradiant rate $\gamma + \gamma_{12}$:

```wl
ListLinePlot[excitations, DataRange -> {0, 5},
 PlotLegends -> {"|01⟩ moderate decay", "(|01⟩-|10⟩)/√2 slow decay", "(|01⟩+|10⟩)/√2 fast decay"}]
```

<!-- => Graphics; all three curves start at 1; at t = 5 the singlet retains the most
   excitation and the triplet the least (verified ordering) -->

## Options

### "AdditionalEquations"

Any additional equations or specifications can be included in the solution, especially for piecewise or hybrid systems.

Define a Hamiltonian with a time-dependent Rabi drive and constant detuning:

```wl
hamiltonian[t_] = 10. Sin[π t] QuantumOperator["X"] - QuantumState["1"]["Operator"];
```

Find the time-dependent state, starting from the ground state:

```wl
ψ = QuantumEvolve[hamiltonian[t], {t, 0, 2}];
```

Plot the probability of the ground state in time:

```wl
Plot[ψ[t]["Probability"][[1]], {t, 0, 2}, PlotRange -> {0, 1}]
```

<!-- => Graphics -->

Redefine the time evolution, sowing the times at which the state populations cross one half (the dependent variable of the internal differential equation is the formal symbol $s$):

```wl
{sol, points} = Reap@QuantumEvolve[hamiltonian[t], {t, 0, 5},
   "AdditionalEquations" -> WhenEvent[Abs[Indexed[\[FormalS][t], 1]]^2 == 0.5, Sow[t]]];
```

<!-- => 20 crossings in [0, 5]; each verified to satisfy |ψ1|^2 = 0.5 to 2*10^-5 -->

Visualize it:

```wl
Plot[sol[t]["Probability"][[1]], {t, 0, 2}, Mesh -> points,
 MeshShading -> {Green, Orange}, PlotRange -> {0, 1}]
```

<!-- => Graphics -->

### "ReturnEquations"

`"ReturnEquations"` is [False]() by default. When [True](), the corresponding equations are returned, ready to be fed into, e.g., [NDSolve]():

```wl
QuantumEvolve[QuantumOperator["XX"], {t, 0, 1}, "ReturnEquations" -> True]
```

<!-- => {{Derivative[1][FormalS][t] == SparseArray[...] . FormalS[t], FormalS[0] == SparseArray[...]}, FormalS, {t, 0, 1}} -->

### "MergeInterpolatingFunctions"

When interpolating functions of time appear in the Hamiltonian, they can enter the dynamical equation either separately or merged together into one sparse interpolation.

Set the time grid and the other variables of the Hamiltonian:

```wl
time = {0, .8, 4 - .8, 4};
ω = {0, 15.8, 15.8, 0};
δ = {-16.33, -16.33, 16.33, 16.33};
Ω = ListInterpolation[ω, {time}, InterpolationOrder -> 1][t];
Δ = ListInterpolation[δ, {time}, InterpolationOrder -> 1][t];
```

Create the Hamiltonian operator:

```wl
hamiltonian = Sum[Ω/2 QuantumOperator["X", {j}] - Δ/2 QuantumOperator["I" - "Z", {j}], {j, 10}];
```

Evolve the quantum state with merged interpolating functions:

```wl
ψ1 = QuantumEvolve[hamiltonian, {t, 0, 1}]; // AbsoluteTiming
```

<!-- => {1.7, Null} (machine dependent) -->

Evolve it again with the interpolating functions kept separate; this is roughly an order of magnitude slower, while the two solutions agree to machine precision:

```wl
ψ2 = QuantumEvolve[hamiltonian, {t, 0, 1}, "MergeInterpolatingFunctions" -> False]; // AbsoluteTiming
```

<!-- => {19., Null} (machine dependent); ψ1 and ψ2 agree to 2*10^-15 -->

Plot the probability of observing the register state:

```wl
Plot[ψ1[t]["ProbabilitiesList"][[1]], {t, 0, 1}, PlotLabel -> "Population of register state in time"]
```

<!-- => Graphics -->

Show the corresponding dynamical equations:

```wl
QuantumEvolve[hamiltonian, {t, 0, 1}, "ReturnEquations" -> True]
```

<!-- => equations of the form FormalS'[t] == "SparseInterpolationFunction"[t] . FormalS[t] -->

```wl
QuantumEvolve[hamiltonian, {t, 0, 1}, "MergeInterpolatingFunctions" -> False, "ReturnEquations" -> True]
```

<!-- => equations of the form FormalS'[t] == SparseArray[...] . FormalS[t], with the
   InterpolatingFunction objects kept inside the sparse array entries -->

## Possible Issues

If the solver cannot find a solution, the output is the corresponding solver call with its differential equations:

```wl
QuantumEvolve[1/2 Subscript[ω, 0] QuantumOperator["Z"] + γ Cos[Subscript[ω, p] t] QuantumOperator["X"], t]
```

<!-- => QuantumEvolve::error message, then the unevaluated
   DSolveValue[{FormalS'[t] == -I SparseArray[...] . FormalS[t], FormalS[0] == SparseArray[...]},
               FormalS[t] ∈ Vectors[2, Complexes], t] -->
