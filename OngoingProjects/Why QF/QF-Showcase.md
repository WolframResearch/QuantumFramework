# What Makes Wolfram QuantumFramework Unique: A Verified Showcase

Wolfram QuantumFramework is a quantum laboratory built into the Wolfram Language. The objects of quantum theory, states, operators, channels, measurements, circuits, are exact symbolic expressions: they carry free parameters, compose with one another in closed form, and answer questions about themselves through one uniform interface. A computation here typically ends in a formula, the propagator $U(t)$ of a time-dependent Hamiltonian, an optimization landscape you can differentiate, an entanglement entropy $S(t)$ you can solve, rather than a table of numbers. The same objects run in any finite dimension, qubits, qutrits, spin-$j$ multiplets, $N$-level systems, and they reach outward as well, to OpenQASM, to device transpilation, to real quantum hardware. Its natural habitat is the computation whose deliverable is insight: a closed form, an exact threshold, a verified physical statement. These are not estimates that tighten with more shots or a finer grid; they are *exact*, so where a numeric run or a finite measurement returns a value, QF returns the closed form, the exact integer, or the algebraic root that the value only approximates.

This document makes that concrete as seven claims, each developed through worked examples:

1. **Symbolic-first, exact algebra by default.** Time-dependent Hamiltonians solved to closed-form propagators, an exact QAOA landscape, a many-body entanglement entropy as a formula in time.
2. **One object model spans all of quantum theory.** Channels, POVMs, bosonic Fock space, and phase-space quasi-probabilities are the same kind of object as a gate; a master equation solves symbolically for every state and every rate at once; a basis rides on every object; the model runs in any finite dimension, from a spin-1 qutrit to a five-site ring; and it closes with a full state-tomography loop.
3. **Measurement records are quantum wires.** A measurement appends a pointer qudit instead of collapsing the state, so feedforward is a controlled gate on a record: teleportation runs without a classical channel, a coherent error is digitized and corrected at codeword fidelity exactly $1$, and a measurement can even be run backwards.
4. **One object model, three computational engines.** The same circuit runs as exact state-vector algebra, as a steered tensor-network contraction, or as a compiled stabilizer tableau that carries a thousand qubits with ease, and the representations convert into one another exactly.
5. **A quantum object is an ordinary Wolfram Language expression.** Native functions operate *on* the object: `Exp` exponentiates a generator into a propagator, `D` differentiates a gate into another operator, `Plus` and `Commutator` build a Hamiltonian, and the closed forms they yield flow on into solvers.
6. **Every object draws itself.** Bloch spheres, amplitude charts that carry the phases, measurement histograms, circuit diagrams, and the tensor network behind a circuit, each a one-call property of the object it depicts.
7. **An interoperability hub, down to real hardware.** OpenQASM in and out natively; transpilation to a device target; submission to a real QPU, with IBM as the worked example.

The claims build on one another, symbolic objects at the base and real hardware at the summit, so we begin with the algebra everything else inherits.

---

## Claim 1: Symbolic-first, exact algebra by default

**Key feature:** symbolic algebra pervades the *whole* object model, not one function. The building blocks you reach for, a **state** (`QuantumState`), a **Hamiltonian** (`QuantumOperator`), its **time evolution** (`QuantumEvolve`), and a **circuit** (`QuantumCircuitOperator`), are each fully symbolic, and they compose into exact closed-form results. We take them in that order, starting with the simplest: a state whose amplitudes are symbols.

### A symbolic state
A `QuantumState` can carry free symbols. The object to build is the general qubit $|\psi\rangle = \cos\alpha\,|0\rangle + e^{i\beta}\sin\alpha\,|1\rangle$: not one state but the whole two-parameter family at once. Build it:
```wolfram
psi = QuantumState[{Cos[\[Alpha]], Sin[\[Alpha]] Exp[I \[Beta]]}]
```
The object carries both amplitudes as symbols. Read them back as a plain vector:
```wolfram
psi["StateVector"] // Normal
```
The amplitudes return exactly as $\{\cos\alpha,\ e^{i\beta}\sin\alpha\}$, symbols intact.

The same object reports its exact Bloch geometry, the expectation triple $\vec r = \big(\langle X\rangle, \langle Y\rangle, \langle Z\rangle\big)$ that places the qubit on the unit sphere (`ComplexExpand` treats $\alpha$ and $\beta$ as real):
```wolfram
psi["BlochVector"] // ComplexExpand // FullSimplify
```
The Bloch vector returns exactly as $\bigl(\sin 2\alpha\,\cos\beta,\ \sin 2\alpha\,\sin\beta,\ \cos 2\alpha\bigr)$: the full geometry as functions of $\alpha$ and $\beta$, without committing to numbers.

### A time-dependent symbolic Hamiltonian
Hamiltonians need not be static. The following is a spin-$\tfrac12$ in a magnetic field that **rotates** in the $xy$-plane at frequency $\omega$, on top of a static $z$-field, the standard model of magnetic resonance and of a driven qubit:
$$
H(t) = \frac{\omega_0}{2}\,Z + \frac{\Omega}{2}\big(\cos\omega t\,X + \sin\omega t\,Y\big)
$$
It is **explicitly time-dependent**, built as a function `h[t]` whose operator carries the time `t`:
```wolfram
ClearAll[h]
h[t_] := \[Omega]0/2 QuantumOperator["Z"] +
   \[CapitalOmega]/2 (Cos[\[Omega] t] QuantumOperator["X"] + Sin[\[Omega] t] QuantumOperator["Y"])
```
Because the field direction turns, $H$ does **not commute with itself at different times**, so the propagator is not simply $e^{-i\int H\,dt}$; it needs genuine time-ordering. Check the non-commutativity directly, letting the kernel's `Commutator` act on the two operator objects $h(t_1)$ and $h(t_2)$ themselves:
```wolfram
Commutator[h[t1], h[t2]]["Matrix"] // Normal // FullSimplify
```
The commutator is nonzero, its diagonal carrying $\mp\tfrac{i}{2}\Omega^2\sin[(t_1-t_2)\omega]$: this is the genuinely hard, time-ordered case, not a disguised static problem.

### Its propagator, in closed form
`QuantumEvolve` is QF's one entry point for time evolution: it evolves a state forward or returns the propagator itself, under closed-system (Schrodinger) or open-system (Lindblad) dynamics, numerically or, when the inputs carry symbols, in closed form. Hand it the Hamiltonian `h[t]`, `None` in place of an initial state, and the time variable `t`, and it returns the *propagator* $U(t)$ solving $i\,\partial_t U = H(t)\,U$. (Verified here: $U(0)=I$ and $i\,U'(t)-H(t)\,U(t)=0$ exactly.)
```wolfram
urot = QuantumEvolve[h[t], None, t]
```
Simplify its matrix under real parameters and $\Omega > 0$:
```wolfram
Assuming[{\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, t \[Element] Reals},
  FullSimplify[ComplexExpand[Normal[urot["Matrix"]]]]] // MatrixForm
```
The matrix returns exactly as the lab-frame Rabi propagator, written here with the detuning $\delta \equiv \omega - \omega_0$ and the generalized Rabi frequency $\Omega_R \equiv \sqrt{\Omega^2 + \delta^2}$:
$$
U(t) = \begin{pmatrix}
e^{-i\omega t/2}\!\left(\cos\dfrac{\Omega_R t}{2} + i\dfrac{\delta}{\Omega_R}\sin\dfrac{\Omega_R t}{2}\right) & -\,i\,e^{-i\omega t/2}\,\dfrac{\Omega}{\Omega_R}\sin\dfrac{\Omega_R t}{2} \\[3mm]
-\,i\,e^{i\omega t/2}\,\dfrac{\Omega}{\Omega_R}\sin\dfrac{\Omega_R t}{2} & e^{i\omega t/2}\!\left(\cos\dfrac{\Omega_R t}{2} - i\dfrac{\delta}{\Omega_R}\sin\dfrac{\Omega_R t}{2}\right)
\end{pmatrix}
$$
The $e^{\mp i\omega t/2}$ factors are the rotating-frame phases; the $\cos$ and $\sin$ of $\Omega_R t/2$ are the dressed-state oscillations. The physically meaningful readout is the transition probability $|\langle 1|U(t)|0\rangle|^2$. Pull the lower-left matrix entry, the $0\to1$ amplitude:
```wolfram
amp = Normal[urot["Matrix"]][[2, 1]];
```
Then square its magnitude under the same real-parameter assumptions:
```wolfram
Assuming[{\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, t \[Element] Reals},
  FullSimplify[ComplexExpand[Abs[amp]^2]]]
```
The result is the textbook Rabi formula, $P_{0\to1}(t) = \frac{\Omega^2}{\Omega^2 + \delta^2}\,\sin^2\!\big(\tfrac{t}{2}\sqrt{\Omega^2 + \delta^2}\big)$. On resonance ($\omega=\omega_0$, so $\delta=0$) this is $\sin^2(\Omega t/2)$: complete Rabi flopping between $|0\rangle$ and $|1\rangle$. The oscillation rate $\sqrt{\Omega^2+\delta^2}$ is the generalized Rabi frequency, the dressed-state gap in the rotating frame, here recovered as a closed-form function of the drive detuning and strength.

### When the closed form needs special functions: Landau-Zener
The rotating-field propagator was elementary, but symbolic evolution is not confined to elementary functions. Sweep the energy bias *linearly* through an avoided crossing at rate $v$ with a constant gap $\Delta$, the Landau-Zener model $H(t) = \tfrac{v t}{2}\,Z + \tfrac{\Delta}{2}\,X$, the universal question of whether a system dragged through a level crossing stays in its energy branch or jumps the gap, again time-dependent and non-commuting:
```wolfram
hlz = v t/2 QuantumOperator["Z"] + \[CapitalDelta]/2 QuantumOperator["X"]
```
Ask `QuantumEvolve` for the propagator:
```wolfram
ulz = QuantumEvolve[hlz, None, t]
```
Then simplify its full matrix:
```wolfram
lzmat = FullSimplify[Normal[ulz["Matrix"]]]
```
The matrix returns as a genuinely special-function object. Its first entry alone reads $e^{-i v t^2/4}\,{}_1F_1\!\big(\tfrac{i\Delta^2}{8v};\ \tfrac12;\ \tfrac{i}{2}\,v t^2\big)$, and across the four entries Kummer ${}_1F_1$, fractional-order Hermite, and parabolic cylinder functions $D_\nu$ appear side by side, glued by $\Gamma$-function coefficients, the single dimensionless ratio $\Delta^2/v$ controlling every index and phase. They are one family written several ways: every piece solves the same Weber equation, here in the raw form the symbolic solver found it. The matrix checks out as a genuine propagator (numerically $U(0)=I$, solves the Schrodinger equation, unitary), and in the long-time limit its asymptotics collapse to the famous Landau-Zener transition probability $P = e^{-\pi\Delta^2/2v}$: symbolic evolution reaching past elementary functions to the special-function solution the physics demands.

The payoff of having the *whole* propagator, not just its asymptote, is sharpest as a picture. First freeze the special-function matrix into a numeric parametric operator: fix $v = 1$ and $\Delta = \sqrt{g}$ so a single dimensionless knob $g = \Delta^2/v$ tunes the adiabaticity, and declare $g$ and the time $t$ as parameters:
```wolfram
op = N@QuantumOperator[lzmat /. {v -> 1., \[CapitalDelta] -> Sqrt[g]},
   "Parameters" -> {g, t}]
```
Now `op[g, t]` is the Landau-Zener propagator at adiabaticity $g$ and time $t$, evaluated numerically on demand.

The survival is the probability of staying in the diabatic ground state $|0\rangle$ across the sweep. Propagate $|0\rangle$ from well before the crossing at $t = -10$ to time $t$, with $U(t,-10) = \mathrm{op}[g,t]\,\mathrm{op}[g,-10]^\dagger$, and read off $|\langle 0|\psi(t)\rangle|^2$:
```wolfram
ClearAll[survival];
survival[g_, t_] = With[{c = op[g, -10]["Dagger"][]},
   Abs[(QuantumState["0"]["Dagger"]@op[g, t][c])["Scalar"]]^2];
```
`survival[g, t]` is now a closed-form function of adiabaticity and time, assembled entirely by composing QF objects, with no differential solver in the loop.

Take three regimes: a fast (diabatic) sweep, an intermediate one, and a slow (adiabatic) sweep, set by the single ratio $\Delta^2/v$:
```wolfram
ratios = {0.1, 0.6, 1.6};
```
Plot the survival across the crossing for each, with the asymptotic Landau-Zener formula $e^{-\pi\Delta^2/2v}$ drawn as dashed lines:
```wolfram
Plot[Evaluate@Table[survival[g, t], {g, ratios}], {t, -10, 10},
 Frame -> True, PlotRange -> {0, 1.05},
    FrameLabel -> {"time  t  (avoided crossing at t = 0)",
        "Survival  |\[LeftAngleBracket]0|\[Psi](t)\[RightAngleBracket]|^2"},
    PlotLegends -> Placed[{"\[CapitalDelta]\.b2/v = 0.1  (diabatic)", "0.6  (intermediate)",
    "1.6  (adiabatic)"}, {0.78, 0.6}],
    Epilog -> {Gray, Dashed, Sequence @@ (Line[{{-10, Exp[-Pi #/2]}, {10, Exp[-Pi #/2]}}] & /@ ratios)},
    GridLines -> {{0}, None}]
```
The dashed lines are the textbook scalar; the solid curves are the exact survival from the special-function propagator. Each plunges through the crossing at $t = 0$ and then *rings* around its asymptote with a slowly decaying envelope, the Landau-Zener-Stückelberg oscillations that the parabolic cylinder functions carry and the bare formula cannot show. The asymptotic result is nothing but the long-time average of the exact answer, and the single ratio $\Delta^2/v$ orders the whole family: a fast (diabatic) sweep barely dips and stays in its branch, a slow (adiabatic) sweep transfers almost completely. The closed form is the entire time-resolved crossing, ringing and all, where the textbook gives only the final number.

### A parametric circuit: QAOA, built entirely from a graph
QAOA is the canonical parametric circuit for optimization: one layer alternates a cost unitary $e^{-i\gamma H_C}$ with a mixer $e^{-i\theta B}$, and we tune the two angles to maximize a classical objective. The objective here is MaxCut: color the vertices of a graph with two colors $s_i = \pm1$ so that as many edges as possible join different colors, $\max_{s\,\in\,\{\pm1\}^n} \sum_{(i,j)\in E} \tfrac{1 - s_i s_j}{2}$. Everything below is *derived from the problem graph*, so nothing is hardcoded. Take the complete graph on four vertices, $K_4$, a frustrated instance where you cannot cut all six edges:
```wolfram
graph = CompleteGraph[4]
```
Read off the classical optimum, its best two-coloring of the vertices:
```wolfram
FindMaximumCut[graph]
```
The classical answer is a cut of $4$ of the $6$ edges, splitting the four vertices into two pairs. Record the vertex count for reuse:
```wolfram
n = VertexCount[graph]
```

The MaxCut cost Hamiltonian $H_C = \sum_{(i,j)}\tfrac{1 - Z_i Z_j}{2}$ is assembled straight from the edge list, one $(1-Z_iZ_j)/2$ term per edge, and carries a label that will name its block in the diagram below:
```wolfram
hcost = QuantumOperator[Total[With[{id = QuantumOperator[StringJoin[ConstantArray["I", n]]]/2},
   (id - QuantumOperator["ZZ"/2, #] &) /@ (List @@@ EdgeList[graph])]], "Label" -> "Hamiltonian"]
```
Change `graph` and this operator adapts to any graph.

The ansatz is built one labeled piece at a time. State prep puts every qubit in $|+\rangle$, an equal superposition over all cuts:
```wolfram
stateprep = QuantumCircuitOperator[{StringJoin[ConstantArray["0", n]], "H" -> Range[n]}, "Label" -> "Prep"]
```
The cost layer applies one $ZZ(\gamma)$ rotation per edge, carrying the angle $\gamma$:
```wolfram
qcost = QuantumCircuitOperator[("R"[\[Gamma], "ZZ"] -> # &) /@ (List @@@ EdgeList[graph]), "Parameters" -> {\[Gamma]}, "Label" -> "Cost"]
```
The mixer applies an $R_X(\theta)$ to every vertex:
```wolfram
mixer = QuantumCircuitOperator[{"RX"[\[Theta]] -> VertexList[graph]}, "Parameters" -> {\[Theta]}, "Label" -> "Mixer"]
```
Compose the three into one parametric layer:
```wolfram
oneLayer = QuantumCircuitOperator[{stateprep, qcost, mixer}, "Parameters" -> {\[Theta], \[Gamma]}]
```
Not a single qubit index was written by hand.

The cost landscape is the expectation $\langle\psi|H_C|\psi\rangle$, and here QF leans on a feature worth naming: **a `QuantumCircuitOperator` can hold states as elements**. Composing the prep circuit (the ket), the operator $H_C$, and the prep circuit's `["Dagger"]` (the bra) builds $\langle\psi|H_C|\psi\rangle$ as a *single object*:
```wolfram
meanvalue = (oneLayer /* hcost /* oneLayer["Dagger"])
```
That object draws itself:
```wolfram
meanvalue["Diagram"]
```
The diagram shows the bra-ket sandwich block by block, the labeled `Prep`/`Cost`/`Mixer` layer, the `Hamiltonian` in the middle, and the daggered layer mirrored on the right: the expectation value as a picture. Calling the same object with `[]["Scalar"]` reads off the number, no manual inner product:
```wolfram
landscape = FullSimplify[meanvalue[]["Scalar"], Element[{\[Gamma], \[Theta]}, Reals]]
```
The exact QAOA energy landscape returns as $\langle H_C\rangle(\gamma,\theta) = \frac{3}{8}\big(7 + \cos 2\theta + 2\cos 4\gamma\,\sin^2\theta - 8\cos^2\gamma\,\sin\gamma\,\sin 2\theta\big)$, a closed form in the two angles.

Because the landscape is a formula, the whole optimization problem can be looked at before it is solved. Plot the surface over a full period of both angles:
```wolfram
Plot3D[landscape, {\[Gamma], 0, 2 \[Pi]}, {\[Theta], 0, 2 \[Pi]},
  PlotRange -> All,
  AxesLabel -> {"\[Gamma]", "\[Theta]", "cost"},
  ColorFunction -> "TemperatureMap",
  PlotPoints -> 50]
```
The surface is the variational problem in full view: smooth, with period $\pi$ in $\theta$ but $2\pi$ in $\gamma$ (the cost angle enters through $\sin\gamma$, not only $4\gamma$, so $\gamma$ needs a full $2\pi$ to repeat), and flat at the value $3$ along the line $\gamma = 0$, where the cost layer is off and the equal superposition cuts $3$ of the $6$ edges on average. What an experiment estimates point by point, with shot noise on every sample, is here one exact surface, and its peaks are not read off a grid but certified algebraically. Take the maximum exactly:
```wolfram
Maximize[landscape, {\[Theta], \[Gamma]}]
```
The maximum returns as a closed-form algebraic number, a root of a cubic equal to $\langle H_C\rangle_{\max}$, with the optimal angles likewise exact `Root` expressions. The lesson is the honest one: a single QAOA layer does *not* reach the true cut of $4$, and the exact landscape certifies that ceiling as a cubic root rather than an estimate. Falling short is the generic behavior of one layer; an instance where $p=1$ happens to saturate the cut is the lucky special case, and climbing the rest of the way to $4$ is what deeper QAOA, more layers, buys. What survives is the showcase point: the variational optimum is an exact symbolic quantity, a root of a cubic with closed-form angles, not a number sampled near the answer.

### The payoff: a full many-body computation in closed form
The same machinery scales to a genuine many-body computation. Build the 2-qubit Heisenberg Hamiltonian $H = J\,(XX + YY + ZZ)$ symbolically (Hamiltonian simulation, the workhorse application of a quantum computer):
```wolfram
hheis = J (QuantumOperator["XX"] + QuantumOperator["YY"] + QuantumOperator["ZZ"])
```
Quench the product state $|01\rangle$ under it, that is, compute $|\psi(t)\rangle = e^{-iHt}\,|01\rangle$:
```wolfram
psit = QuantumEvolve[hheis, QuantumState["01"], t]
```
Display the evolved state at time $t$:
```wolfram
FullSimplify[psit] // TraditionalForm
```
The evolved ket displays in closed form: the excitation oscillates coherently between $|01\rangle$ and $|10\rangle$ at angular frequency $2J$, inside an overall phase $e^{iJt}$. Now read off the **entanglement entropy of one qubit as a function of time**, $S(t) = -\operatorname{Tr}\big[\rho_1\log_2\rho_1\big]$ with $\rho_1 = \operatorname{Tr}_2\,|\psi(t)\rangle\langle\psi(t)|$, the standard measure of how entangled the two spins have become, and a derived, nonlinear quantity:
```wolfram
s = FullSimplify[QuantumPartialTrace[psit, {2}]["VonNeumannEntropy"], t J \[Element] Reals]
```
QF returns the **entire curve as a formula in $t$ and $J$**: the binary entropy of the two oscillating populations, zero whenever the state is a product, a full bit at the Bell points, periodic with the exchange oscillation. Because it is a formula, it can be differentiated, expanded, or solved in closed form, and the same approach extends to larger symbolic spin chains.

Plot it (at $J = 1$):
```wolfram
Plot[s /. {J -> 1}, {t, 0, Pi},
  Frame -> True, FrameLabel -> {"time  J t", "entanglement entropy  S(t)  [bits]"}, PlotRange -> {0, 1.05}]
```

---

## Claim 2: One object model spans all of quantum theory

**Key feature:** if Claim 1 went *deep* on one object's symbolic dynamics, this claim goes *wide*. QF gives you whole *categories* of object beyond unitary gates: open systems, generalized measurements, continuous-variable optics, phase space, and state estimation. Two generalities run underneath the whole model: every object carries its basis, and every object lives in any finite dimension, qutrits and $N$-level systems as readily as qubits. The point is *scope*, the range of quantum theory QF represents as objects you compute with directly.

### An open-system object: a quantum channel
A noise process is its own object, a `QuantumChannel`, not a gate. Take a generic pure qubit on the Bloch sphere, polar angle $\theta$ and azimuth $\phi$, and confirm it is pure through the purity $\operatorname{Tr}[\rho^2]$, which equals $1$ exactly for a pure state and drops below $1$ for a mixed one. Build the state:
```wolfram
psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
```
Read its purity:
```wolfram
psi["Purity"] // FullSimplify
```
The purity is $1$, as it must be for any pure state.

Send it through amplitude damping, the decay channel $\rho \mapsto K_0\rho K_0^\dagger + K_1\rho K_1^\dagger$ with $K_1 = \sqrt{\gamma}\,|0\rangle\langle1|$ and $K_0 = |0\rangle\langle0| + \sqrt{1-\gamma}\,|1\rangle\langle1|$: the excitation is lost with probability $\gamma$, the model of spontaneous emission. The question: how much purity does that cost, $\operatorname{Tr}[\rho^2] - \operatorname{Tr}\big[\mathcal{E}_\gamma(\rho)^2\big]$, as an exact function of the damping rate $\gamma$ and the state:
```wolfram
FullSimplify[psi["Purity"] - QuantumChannel["AmplitudeDamping"[\[Gamma]]][psi]["Purity"],
  {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals, 0 <= \[Gamma] <= 1}]
```
The purity loss returns exactly as $2\,\gamma(1-\gamma)\,\sin^4\tfrac{\theta}{2}$. The pure input has become mixed: irreversible, non-unitary dynamics captured directly as a channel object. The loss depends only on $\gamma$ and the excited-state weight $\sin^2(\theta/2)$, never on the azimuth $\phi$, and it vanishes at $\gamma=0$ (no damping), at $\gamma=1$ (fully relaxed to the pure ground state), and at $\theta=0$ ($|0\rangle$, which damping leaves untouched).

### The master equation, solved once for every state and every rate
A channel is one snapshot of decoherence; the full dynamics is the Lindblad master equation, $\dot\rho = -i[H,\rho] + \sum_i \gamma_i \big(L_i\rho L_i^\dagger - \tfrac12\{L_i^\dagger L_i,\rho\}\big)$, and `QuantumEvolve` solves it *symbolically*. Hand it a density matrix whose four entries are symbols and rates that stay symbols, and one solve answers for every qubit experiment at once. Set the physical region: positive rates, a real level splitting, and forward time:
```wolfram
$Assumptions = \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0 && \[CapitalOmega] \[Element] Reals && t >= 0;
```
Write the decay operator $\sigma_- = |0\rangle\langle 1|$ as its matrix, so it empties the excited state:
```wolfram
\[Sigma]m = QuantumOperator[{{0, 1}, {0, 0}}]
```
Take the all-states-at-once initial density matrix $\rho_0$, every entry a symbol:
```wolfram
\[Rho]0 = QuantumState[Array[\[Rho], {2, 2}]]
```
Solve the master equation for a qubit with level splitting $\Omega$, energy relaxation at $\gamma_1$ (jump operator $\sqrt{\gamma_1}\,\sigma_-$), and pure dephasing at $\gamma_\phi$ (jump operator $\sqrt{\gamma_\phi}\,Z$):
```wolfram
\[Rho]t = QuantumEvolve[\[CapitalOmega]/2 QuantumOperator["Z"],
   {Sqrt[\[Gamma]1] \[Sigma]m, Sqrt[\[Gamma]\[Phi]] QuantumOperator["Z"]}, \[Rho]0, t]
```
Read its density matrix in the computational basis:
```wolfram
FullSimplify /@ Normal[\[Rho]t["Computational"]["DensityMatrix"]] // MatrixForm
```
The matrix returns exactly as the complete relaxation law for an arbitrary qubit:
$$
\rho(t) = \begin{pmatrix}
\rho_{11} + \left(1 - e^{-\gamma_1 t}\right)\rho_{22} & \rho_{12}\, e^{-\left(\frac{\gamma_1}{2} + 2\gamma_\phi + i\Omega\right) t} \\[1mm]
\rho_{21}\, e^{-\left(\frac{\gamma_1}{2} + 2\gamma_\phi - i\Omega\right) t} & \rho_{22}\, e^{-\gamma_1 t}
\end{pmatrix}
$$
The excited population empties at $\gamma_1$, the lost weight lands on the ground state, and the coherence precesses at $\Omega$ while decaying at $\Gamma_2 = \tfrac{\gamma_1}{2} + 2\gamma_\phi$: every experiment on every qubit, in one matrix. With $T_1 = 1/\gamma_1$ and $T_2 = 1/\Gamma_2$ read off the exponents, the textbook bound between the two times is now a theorem about this solution:
```wolfram
Reduce[ForAll[{\[Gamma]1, \[Gamma]\[Phi]}, \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0, 1/(\[Gamma]1/2 + 2 \[Gamma]\[Phi]) <= 2/\[Gamma]1]]
```
The reduction returns `True`: $T_2 \le 2T_1$ holds for every admissible rate pair, handed over by the computation rather than quoted from a textbook, and the same `Reduce` asked for equality returns $\gamma_\phi = 0$ (verified): the famous factor of two is saturated exactly when dephasing vanishes.

### A generalized measurement: a POVM
Real detectors are not always projective. A POVM is a set of positive effects $\{E_i\}$ summing to the identity, $\sum_i E_i = I$, with outcome probabilities given by the Born rule $p_i = \operatorname{Tr}[\rho\,E_i]$; build one from your own effects with `QuantumMeasurementOperator[{E1, E2, ...}]`. A *symmetric informationally complete* POVM is the maximally symmetric case: $d^2$ rank-one effects $E_i = \tfrac1d|\psi_i\rangle\langle\psi_i|$ with equal pairwise overlaps $|\langle\psi_i|\psi_j\rangle|^2 = \tfrac{1}{d+1}$, four outcomes on a two-level system, pointing along the vertices of a tetrahedron in the Bloch ball. Measure a generic Bloch qubit $\psi(\theta,\phi)$, the Born rule applied symbolically to all four effects. Build the state:
```wolfram
psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
```
Build the tetrahedron SIC POVM:
```wolfram
sic = QuantumMeasurementOperator["TetrahedronSICPOVM"]
```
Apply the Born rule $p_i = \operatorname{Tr}[\rho\,E_i]$ to all four effects at once:
```wolfram
probs = FullSimplify[ComplexExpand[Tr[psi["DensityMatrix"] . Normal[#]] & /@ sic["POVMElements"]]]
```
The four probabilities return exactly as $p_i = \tfrac14\big(1 + \hat r \cdot \hat n_i\big)$, the first reading $\tfrac{1+\cos\theta}{4}$ and the others $\tfrac{1}{12}\big(3 - \cos\theta + \sqrt2\,\sin\theta\,(\ldots)\big)$ along the remaining tetrahedron directions $\hat n_i$, with $\hat r$ the state's Bloch vector; they sum to $1$ identically (verified). Pin the state to the north pole and the general formulas collapse to numbers:
```wolfram
probs /. \[Theta] -> 0
```
The distribution collapses to simple rationals: $|0\rangle$ puts a half on the tetrahedron vertex nearest the pole and a sixth on each of the other three. And because the SIC is *informationally complete*, these four probabilities are not merely statistics of the state, they **are** the state: $\theta$ and $\phi$ can be solved back from them, which is what makes SIC measurements the natural alphabet of tomography.

The SIC is not only a symbolic Born-rule calculation; it is a measurement you drop straight into a circuit, alongside the gates, on whichever wire you choose. Place the tetrahedron SIC as a measurement on the first qubit:
```wolfram
sicQubit1 = QuantumMeasurementOperator["TetrahedronSICPOVM", {1}]
```
Build a small entangling circuit that ends in that POVM:
```wolfram
circ = QuantumCircuitOperator[{"H" , "RY"[1.05] -> 2,
       "CNOT" , "RZ"[2.04], sicQubit1}]
```
Draw it:
```wolfram
circ["Diagram"]
```
The diagram shows the four-outcome POVM sitting at the end of the first wire like any other measurement. Applying the circuit returns a `QuantumMeasurement`, and that object draws its own outcome statistics:
```wolfram
circ[]["ProbabilitiesPlot"]
```
Four bars summing to one, one tetrahedron outcome suppressed and another enhanced: the SIC reading qubit $1$. The twist is that qubit $1$ is *entangled* with qubit $2$, so its reduced state is mixed, its Bloch vector pulled inside the ball (length less than $1$), and the four probabilities are correspondingly drawn in from the pure-state extremes toward the uniform $\tfrac14$. The same informationally complete alphabet reports a mixed, entangled qubit as faithfully as a pure one, which is exactly what tomography of a real, noisy qubit demands.

### Continuous-variable quantum optics
Bosonic Fock space, not qubits. The state of an ideal laser mode is the coherent state $|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\tfrac{\alpha^n}{\sqrt{n!}}\,|n\rangle$, defined as the eigenstate of the photon annihilation operator, $\hat a\,|\alpha\rangle = \alpha\,|\alpha\rangle$, with a *complex* eigenvalue $\alpha$ because $\hat a$ is not Hermitian. QF builds these objects directly in a truncated Fock space (sixteen levels by default). Check the eigenvalue equation with $\alpha$ fully symbolic. Build the annihilation operator:
```wolfram
aOp = AnnihilationOperator[]
```
Build the coherent state, left unnormalized so the symbolic algebra stays clean:
```wolfram
alphaKet = CoherentState["Normalized" -> False]
```
Subtract $\alpha\,|\alpha\rangle$ from $\hat a\,|\alpha\rangle$ and read the amplitudes (the truncation edge excluded):
```wolfram
Simplify[Most[(aOp[alphaKet] - \[FormalAlpha] alphaKet)["AmplitudesList"]]]
```
Every amplitude of $\hat a|\alpha\rangle - \alpha|\alpha\rangle$ returns identically zero in the symbolic $\alpha$ (the truncation-edge amplitude excluded): the defining equation, verified for the whole family of coherent states at once.

The property that earns these states the title "most classical states of light" is dynamical: under the oscillator Hamiltonian $H = \omega\,(\hat a^\dagger\hat a + \tfrac12)$ a coherent state *stays* a coherent state, $e^{-iHt}\,|\alpha\rangle = e^{-i\omega t/2}\,|\alpha\,e^{-i\omega t}\rangle$, its label circling the phase-space origin exactly like a classical oscillator amplitude, with only the vacuum phase out front. Ask `QuantumEvolve` for the oscillator propagator, symbolically:
```wolfram
uosc = QuantumEvolve[\[Omega] (aOp["Dagger"] @ aOp + 1/2), None, t]
```
Then check the coherent-state identity it implies, with $\alpha$, $\omega$, and $t$ all symbols:
```wolfram
uosc @ CoherentState[] == E^(-I \[Omega] t/2) CoherentState[][\[FormalAlpha] E^(-I \[Omega] t)]
```
`True`, with $\alpha$, $\omega$, and $t$ all symbols: one evolution settles the dynamics of every coherent state at every frequency, an exact statement about the family rather than a simulation of one member.

What still separates this most-classical light from its neighbors is photon statistics. The second-order coherence at zero delay, $g^{(2)}(0) = \dfrac{\langle a^\dagger a^\dagger a\, a\rangle}{\langle a^\dagger a\rangle^2} = \dfrac{\langle n(n-1)\rangle}{\langle n\rangle^2}$, is the normalized chance of detecting two photons at once; it sorts light by photon statistics, antibunched ($<1$, nonclassical), coherent ($=1$), and bunched thermal ($>1$).
```wolfram
g2s = <|"Fock" -> G2Coherence[FockState[1]], "Coherent" -> G2Coherence[CoherentState[20][1.5]],
  "Thermal" -> G2Coherence[ThermalState[1.0, 20]]|>
```
Each class lands where the formula predicts: zero for the single photon, one for the laser, two for thermal light. The three values chart directly:
```wolfram
BarChart[g2s, ChartLabels -> {"Fock |1\[RightAngleBracket]", "Coherent", "Thermal"},
  PlotLabel -> "Second-order coherence g2(0)", Frame -> True, PlotRange -> {0, 2.2}]
```
The ordering on the chart, antibunched below $1$, coherent at $1$, bunched above, is the photon-statistics fingerprint that tells these light sources apart where a single intensity measurement cannot.

### Phase-space quasi-probabilities: Wigner and Husimi
The Wigner function is the phase-space quasi-probability $W(x,p) = \dfrac{1}{\pi}\displaystyle\int \langle x - y\,|\,\rho\,|\,x + y\rangle\, e^{2ipy}\,dy$. It integrates to $1$ like a probability density but, unlike one, can go **negative**, and that negativity is a signature of nonclassicality with no classical analogue. The sharpest specimen is a Schrodinger cat state, the even superposition $(|\alpha\rangle + |{-\alpha}\rangle)/\mathcal{N}$ of two coherent states of opposite phase, here with $\alpha = 2$. Build its Wigner function on a phase-space window:
```wolfram
wigner = WignerRepresentation[CatState[][2, 0], {-6, 6}, {-2, 2}, "GridSize" -> 90]
```
Read it at two telling points, the origin and half a fringe up in momentum:
```wolfram
{wigner[0, 0], wigner[0, 1/2]}
```
The first value is positive, the second negative. The first is the identification that matters: the midpoint between the two coherent components, where a classical mixture of them would put essentially nothing, instead carries the strongest feature of the entire function, an interference peak approaching the ceiling $1/\pi$ that no Wigner function can exceed; and the neighboring fringe is firmly negative.

The Wigner function is one of a family of phase-space quasi-probabilities, and its smoothed sibling makes that negativity vanish. The Husimi-$Q$ function $Q(x,p) = \tfrac1\pi\langle\alpha|\rho|\alpha\rangle$, with $\alpha = (x+ip)/\sqrt2$, is the Wigner function convolved with a coherent-state Gaussian; that smoothing forces it non-negative everywhere, at the cost of blurring away the very nonclassicality the Wigner reveals. Build it over the same window, a touch taller in $p$ since the Husimi spreads wider:
```wolfram
husimi = HusimiQRepresentation[CatState[][2, 0], {-6, 6}, {-3, 3}, "GridSize" -> 90]
```
Read it at the origin where the Wigner peaked:
```wolfram
husimi[0, 0]
```
The Husimi returns next to nothing there, against the Wigner's sharp peak: the Gaussian average has cancelled the oscillating fringes. Plot the two surfaces side by side:
```wolfram
Row[{
   Plot3D[wigner[x, p], {x, -6, 6}, {p, -2, 2}, PlotRange -> All, Mesh -> None,
     PlotPoints -> 50, ColorFunction -> "SunsetColors",
     AxesLabel -> {"x", "p", "W"}, PlotLabel -> "Wigner W(x,p)", ImageSize -> Small],
   Plot3D[husimi[x, p], {x, -6, 6}, {p, -3, 3}, PlotRange -> All, Mesh -> None,
     PlotPoints -> 50, ColorFunction -> "SunsetColors",
     AxesLabel -> {"x", "p", "Q"}, PlotLabel -> "Husimi Q(x,p)", ImageSize -> Small]}]
```
The same cat state, two quasi-probabilities. The Wigner on the left carries the two Gaussian lobes at $x = \pm 2\sqrt2$, the coherent components, *plus* the interference ridge oscillating above and below zero between them; that negativity is the signature of the superposition. The Husimi on the right keeps only the two smooth, strictly non-negative lobes, each near $1/2\pi$, the fringes smoothed flat. Both integrate to $1$, but only the Wigner can go negative, which is exactly why Wigner negativity, and not the always-positive Husimi, is the witness of nonclassicality. Decoherence erases the Wigner fringes and leaves the lobes, collapsing the left picture toward the right.

### Every object carries its basis
The scope extends inward too: every state, operator, and measurement above is stored *together with the basis it is written in*. A basis is a stored matrix $B$ whose rows are the basis vectors in computational coordinates, a change of basis is the similarity $A_{\mathcal B} = B^{-1} A B$, and objects written in different bases reconcile automatically through the computational frame. About three dozen named bases ship (Bell, Fourier, Schwinger, Dirac, Gell-Mann, Wigner, MUB, SIC), and the same engine lowers non-computational measurements for OpenQASM export and supplies the Bell-derived basis that two-qubit KAK gate compilation runs in. Watch the representation move while the physics stays put: take the qubit $|+\rangle$, the spread-out superposition $(|0\rangle+|1\rangle)/\sqrt2$ in the computational basis, and write it in the $X$ basis instead:
```wolfram
QuantumState[QuantumState["+"], "X"] // TraditionalForm
```
The state renders as the single ket $|x_+\rangle$: in its own eigenbasis $|+\rangle$ is one basis vector, not a superposition. Same state, different coordinates.

The measurable content is untouched by any such rewriting, and nothing here is qubit-bound. Take a spin-1 particle in its $m_x = +1$ eigenstate, defined directly in the $J_x$ basis:
```wolfram
mx = QuantumState[{0, 0, 1}, "JX"[1]]
```
Display it:
```wolfram
mx // TraditionalForm
```
The state renders as the single ket $|J_1^x\rangle$. Re-express it in the computational ($J_z$) basis:
```wolfram
mxZ = QuantumState[mx, "Computational"[3]]
```
Display it again:
```wolfram
mxZ // TraditionalForm
```
The same state spreads into a three-term superposition over the $J_z$ levels: one ket in its own basis, a sum in the other. The named operator $J_x$ is itself stored in its *own* eigenbasis, so the expectation $\langle J_x\rangle = \langle\psi|J_x|\psi\rangle$ mixes three declared bases at once; compute it from both representations of the state:
```wolfram
<|"JX basis" -> mx["Dagger"][QuantumOperator["JX"[1]][mx]]["Scalar"],
  "Computational basis" -> mxZ["Dagger"][QuantumOperator["JX"[1]][mxZ]]["Scalar"]|>
```
Both entries return as $\langle J_x\rangle = 1$: the state in either coordinate system, the operator in a third, and QF reconciled them all; the invariant survives every change of coordinates.

### Any finite dimension: a spin-1 qutrit and a five-site ring
The same spin-1 system carries honest three-level *dynamics*, not just a change of basis, and the representation extends to any dimension at all. Its zero-field-splitting Hamiltonian $H = 2d\,J_z^2 + e\,(J_x^2 - J_y^2)$, the model of an NV-center defect, is built from QF's native spin-1 operators with both couplings left symbolic. First the three spin-1 operators:
```wolfram
{jx, jy, jz} = QuantumOperator /@ {"JX"[1], "JY"[1], "JZ"[1]}
```
Assemble the zero-field-splitting Hamiltonian from them:
```wolfram
hzfs = 2 \[FormalD] jz^2 + \[FormalE] (jx^2 - jy^2)
```
Display its matrix:
```wolfram
hzfs["Matrix"] // MatrixForm
```
The operator returns fully symbolic, and its matrix displays in the basis the algebra carried along, the $J_x$ eigenbasis: the couplings stay symbols, and the coordinates stay attached to the object.

A Hamiltonian is simplest in its energy eigenbasis, and the operator hands over its own eigenvectors. Build a basis from them and read $H$ there:
```wolfram
QuantumOperator[hzfs, QuantumBasis[hzfs["Eigenvectors"]]]["Matrix"] // Normal // FullSimplify
```
The matrix returns exactly as $\operatorname{diag}\big(0,\ 2d - e,\ 2d + e\big)$: the zero-field-split spectrum in closed form, symbols intact. Building a basis from an object's own eigenvectors works in any dimension, and the spin-$j$ operators ship natively for every $j$.

Dimension need not be a power of two either. A particle hopping on an $N$-site ring lives in an $N$-dimensional Hilbert space, here $N = 5$, with the tight-binding Hamiltonian $H = -\sum_{j=1}^{N}\big(|j\rangle\langle j{+}1| + |j{+}1\rangle\langle j|\big)$, site $N{+}1$ identified with site $1$. Fix the ring size:
```wolfram
n = 5;
```
Build the hopping matrix, $-1$ for every pair of neighboring sites:
```wolfram
ring = -Table[
    Boole[Mod[i - j, n] == 1 || Mod[j - i, n] == 1], {i, n}, {j, n}];
```
Look at it in the site basis:
```wolfram
ring // MatrixForm
```
The matrix displays as a circulant: $-1$ on the two neighbor diagonals and in the far corners, the bond that closes the chain into a ring, and zero everywhere else. In the site basis the physics is local hopping; by Bloch's theorem the same operator should be *diagonal* in the momentum (Fourier) basis. Rewrite it there and read the diagonal:
```wolfram
Normal @ Chop @
  N @ Diagonal[
    QuantumOperator[QuantumOperator[ring], QuantumBasis["Fourier"[n]]][
     "Matrix"]]
```
The diagonal returns exactly the tight-binding dispersion $E(k) = -2\cos(2\pi k / N)$ at the five Bloch momenta, and the off-diagonal terms vanish: a five-level system, diagonalized by a named basis the same way a qubit would be.

### State estimation: tomography from measured counts
Tomography is the inverse of measurement: recover an unknown state from the statistics of measurements in several bases, and it is part of the same object model. The three qubit Pauli bases $X, Y, Z$ are mutually unbiased and together informationally complete: their outcome frequencies estimate the three expectations that pin down the density matrix through $\rho = \tfrac12\big(I + \langle X\rangle X + \langle Y\rangle Y + \langle Z\rangle Z\big)$. QF runs the whole loop. First, pick an unknown state to reconstruct:
```wolfram
true = QuantumState[N @ {Cos[Pi/6], Sin[Pi/6] Exp[I Pi/4]}]
```
Simulate $2000$ measurement shots in each of the three Pauli bases (seeded for reproducibility):
```wolfram
SeedRandom[42];
data = QuantumMeasurementSimulation[true, QuantumMeasurementOperator /@ {"X", "Y", "Z"}, 2000];
```
Look at the raw counts:
```wolfram
Values[data]
```
Three pairs of counts return, one per basis: multinomially sampled tallies scattered around the Born-rule proportions, different on every fresh sample, and all a real experiment ever sees.

Now hand the counts to the estimator, which does linear inversion followed by maximum-likelihood:
```wolfram
est = QuantumStateEstimate[data]
```
Read the reconstructed density matrix:
```wolfram
est["MaximumLikelihoodState"]["DensityMatrix"] // Normal // Chop
```
The recovered matrix sits on top of the true one, whose exact entries are $\rho_{11} = \tfrac34$ and $\rho_{12} = \tfrac{\sqrt{3}}{4}\,e^{-i\pi/4}$, with entry-wise residuals at the scale shot noise sets, $\sim 1/\sqrt{2000}$; which entries deviate, and by how much, changes with every fresh sample. To score the reconstruction in one number, compute its fidelity $F(\rho,\sigma) = \operatorname{Tr}\sqrt{\sqrt{\rho}\,\sigma\,\sqrt{\rho}}$ to the state we started from, equal to $1$ exactly when the two states coincide:
```wolfram
QuantumSimilarity[true, est["MaximumLikelihoodState"]]
```
The fidelity returns very close to unity: the unknown state recovered from finite, noisy synthetic data, using nothing but counts in three bases. Past these unitary measurement bases the same basis engine continues into the over-complete operator frames, where $B^{-1}$ becomes a dual frame and a state turns into a probability vector, the mechanism behind the SIC measurement above.

---

## Claim 3: Measurement records are quantum wires, and feedforward is a controlled gate

**Key feature:** Claim 2 used measurement to *estimate* a state; here measurement itself comes apart. QF implements a projective measurement as von Neumann's measurement *interaction*, not as a collapse: it appends a pointer qudit and entangles it with the system, $V\,|\psi\rangle = \sum_k \big(P_k|\psi\rangle\big)\otimes|k\rangle$. The record of the measurement is therefore itself a quantum wire, numbered $0, -1, -2, \dots$ alongside the system wires, every branch stays in the state, and *feedforward*, acting on the system conditioned on an earlier outcome, is nothing but a gate controlled on a record wire. Protocols that mix measurement with conditioned gates, teleportation, error correction, measurement-based computing, stay inside one coherent symbolic object.

### Teleportation without a classical channel
The cleanest place to see records-as-wires earn their keep is teleportation. In the textbook protocol Alice measures her two qubits in the Bell basis and phones the two bits to Bob, who applies $X^{m_2} Z^{m_1}$ to his half of an entangled pair. The deferred-measurement principle says the phone call is physically optional: condition the corrections on the records themselves. Build exactly that, qubit $1$ carrying the unknown state, qubits $2$ and $3$ in a Bell pair:
```wolfram
tele = QuantumCircuitOperator[{
    "H" -> 2, "CNOT" -> {2, 3},
    "CNOT" -> {1, 2}, "H" -> 1,
    {1}, {2},
    "CNOT" -> {-1, 3},
    "CZ" -> {0, 3}}]
```
Draw the protocol:
```wolfram
tele["Diagram"]
```
The diagram is the whole protocol. The two measurements, written as the bare qubit lists `{1}` and `{2}`, drop their outcomes onto the record wires the intro described, numbered $0$ and $-1$ below the system wires. The corrections then reach back from those record wires: `"CNOT" -> {-1, 3}` is an $X$ on Bob's qubit controlled by the second record, and `"CZ" -> {0, 3}` a $Z$ controlled by the first, both ordinary controlled gates with a record wire in the control slot. No classical bit ever leaves the object; the phone call has become two wires.

Teleport a generic Bloch qubit $\psi(\theta,\phi) = \cos\tfrac\theta2\,|0\rangle + e^{i\phi}\sin\tfrac\theta2\,|1\rangle$. The protocol succeeds if Bob ends up holding exactly the input, $\rho_{\text{Bob}} = |\psi\rangle\langle\psi|$, which for a qubit means the two Bloch vectors coincide. First build the input state:
```wolfram
psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
```
Read its Bloch vector:
```wolfram
FullSimplify[psi["BlochVector"], {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}]
```
It returns the textbook $\big(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta\big)$. Now run the protocol:
```wolfram
out = tele[QuantumTensorProduct[psi, QuantumState["00"]]]
```
Read Bob's Bloch vector, the partial trace discarding Alice's two qubits and the two records, the four wires before Bob's:
```wolfram
FullSimplify[QuantumPartialTrace[out, {1, 2, 3, 4}]["BlochVector"], {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}]
```
Bob's Bloch vector returns as the *same* $\big(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta\big)$, identically in $\theta$ and $\phi$: every pure state teleports exactly, the recovered state *is* the input, with no classical channel anywhere in the computation.

### A coherent error, digitized and corrected
The three-qubit repetition code encodes $\alpha|0\rangle + \beta|1\rangle$ as $\alpha|000\rangle + \beta|111\rangle$ and protects it against bit flips. The classic subtlety is that real errors are not discrete flips; here the error will be a *coherent* rotation $R_X(2\theta)$ on the middle qubit, partway between nothing ($\theta = 0$) and a full flip ($\theta = \pi/2$). The founding theorem of error correction says that measuring the syndrome *digitizes* this continuous error into those two discrete cases. Watch it happen symbolically.

First the syndrome measurement. The parity $Z_1 Z_2$ must be measured without learning anything else about the state: a *degenerate* projective measurement with two rank-two projectors $P_\pm = (I \pm Z \otimes Z)/2$, which is precisely what `QuantumMeasurementOperator` accepts as a list of projectors. First the single-qubit identity:
```wolfram
id2 = IdentityMatrix[2];
```
The two-qubit parity operator $Z \otimes Z$:
```wolfram
zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
```
The syndrome measurement on a pair of qubits, the two rank-two projectors $P_\pm = (I \pm Z \otimes Z)/2$:
```wolfram
syn[q1_, q2_] := QuantumMeasurementOperator[
   {QuantumOperator[(KroneckerProduct[id2, id2] + zz)/2, {q1, q2}],
    QuantumOperator[(KroneckerProduct[id2, id2] - zz)/2, {q1, q2}]}, {q1, q2}];
```
Then the full code: encode, inflict the coherent error, measure both parities, and decode each of the four syndrome patterns into its correction. The decoder is three feedforward gates with mixed polarities, `"C"["X", ones, zeros]` conditions an $X$ on some records reading $1$ and others reading $0$. First the encoder:
```wolfram
enc = QuantumCircuitOperator[{"CNOT" , "CNOT" -> {1, 3}},
   "Label" -> "Encode"]
```
Then the full code, encoder through feedforward correction:
```wolfram
qec = QuantumCircuitOperator[{enc,
        "RX"[2 \[Theta]] -> 2,
        syn[1, 2], syn[2, 3],
        "C"["X", {0}, {-1}],
        "C"["X" -> 2, {0, -1}, {}],
        "C"["X" -> 3, {-1}, {0}]}]
```
Draw it:
```wolfram
qec["Diagram"]
```
The diagram lays out the whole code on five wires: the labeled `Encode` block spreading the data qubit across the three code qubits, the $R_X(2\theta)$ error on qubit $2$, the two parity measurements dropping their syndromes onto record wires $0$ and $-1$, and the three corrections reaching back from those records, each an $X$ controlled on one record and anti-controlled on the other. Run it on a fully symbolic data qubit:
```wolfram
coded = qec[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]]
```
Read the joint distribution of the two syndrome records, the diagonal of their reduced density matrix (the records sit on the first two of the five wires, the data on the last three):
```wolfram
Assuming[{\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1},
  FullSimplify[ComplexExpand[Diagonal[Normal[QuantumPartialTrace[coded, {3, 4, 5}]["DensityMatrix"]]]]]]
```
The four syndrome weights return as $\{\cos^2\theta,\ 0,\ 0,\ \sin^2\theta\}$: the continuous rotation has been *digitized*. The code never sees a partial error; it sees "no error" with probability $\cos^2\theta$ or "qubit $2$ flipped" with probability $\sin^2\theta$, and the weights are independent of $\alpha$ and $\beta$, which is exactly what "measuring the syndrome without learning the state" means.

The feedforward gates then act inside each branch. The theorem to check: the corrected state matches the ideal codeword $\alpha|000\rangle + \beta|111\rangle$ at fidelity $F = 1$, for *every* input and *every* error angle. Compute $F$ symbolically:
```wolfram
Assuming[{\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1},
  FullSimplify[ComplexExpand[QuantumSimilarity[QuantumPartialTrace[coded, {1, 2}],
    enc[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]], "Fidelity"]]]]
```
The fidelity returns as exactly $1$, for every input $(\alpha, \beta)$ and every error angle $\theta$: not a high number for one sampled case, but the digitization-of-errors theorem itself, handed over by the computation.

### Run the measurement backwards
Because no branch was discarded, the measurement is a *reversible* interaction, and reversible in the laboratory, not only on paper. That a measurement can be undone while its record stays coherent, before it irreversibly decoheres, is an experimental fact: the reversal, or 'uncollapsing', of a partial (weak) measurement has been demonstrated on superconducting qubits, returning the qubit to its pre-measurement state. QF exposes its underlying dilation as `m["SuperOperator"]`: the isometry $V$ that carries out the measurement by appending the pointer qudit, $V\,|\psi\rangle = \sum_k \big(P_k|\psi\rangle\big)\otimes|k\rangle$, mapping the qubit's two-dimensional space into the $2\times2$ system-and-pointer space. Because $V$ is an isometry it satisfies $V^\dagger V = I$, so it has a left inverse $V^\dagger$ that un-appends the pointer and undoes the measurement. Wrap that inverse as a labeled operator and run the two back to back, measure then un-measure. Build the input state:
```wolfram
psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
```
The projective measurement on qubit $1$:
```wolfram
m = QuantumMeasurementOperator[{1}]
```
Its dilation isometry, daggered, is the un-measurement:
```wolfram
unmeasure = QuantumOperator[m["SuperOperator"]["Dagger"], "Label" -> "Un-measure"]
```
Compose the two into a measure-then-unmeasure circuit:
```wolfram
backward = QuantumCircuitOperator[{m, unmeasure}]
```
Draw it:
```wolfram
backward["Diagram"]
```
The diagram is two boxes on one wire: the measurement dropping its outcome onto the pointer record, then the labeled `Un-measure` block reaching back to absorb it. Apply the whole circuit to the input and read the state that comes out:
```wolfram
FullSimplify[ComplexExpand[Normal[backward[psi]["StateVector"]]], {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}]
```
The state returns as exactly $\{\cos\tfrac\theta2,\ e^{i\phi}\sin\tfrac\theta2\}$, the input recovered phase and all: the measurement is undone exactly. A projective measurement is irreversible only once its record is discarded; while the record stays a wire, it is reversible physics, $V^\dagger V = I$ exactly as von Neumann wrote it down.

A record, then, is a wire you can condition gates on, trace out, read jointly with the system, or undo. The cost is equally explicit: each measurement appends a pointer qudit, so a circuit with $k$ measurements on $n$ qubits evolves a state of dimension $2^{n+k}$. That growth is not overhead; it is the price of keeping every branch alive, and it is precisely what the teleportation identity, the digitization theorem, and the un-measurement above are made of. Keeping every branch alive is, of course, one computational strategy among several, which raises the next claim: QF has three.

---

## Claim 4: One object model, three computational engines

**Key feature:** the objects of Claims 1 through 3 never specified *how* they are computed. Underneath a circuit application sit three interchangeable engines: exact state-vector algebra (`Method -> "Schrodinger"`), tensor-network contraction (the default), and a compiled stabilizer tableau (`Method -> "Stabilizer"`), each suited to a different regime of symbols, structure, and scale. The exact state-vector path is the dense baseline, a literal $2^n$ matrix-vector product available on demand; the two engines that genuinely change the *representation*, and so repay a closer look, are the tensor-network default and the stabilizer tableau, which we take in turn. Because all three share one object model, their representations then convert into one another exactly.

### The default engine: a steered tensor-network contraction
Applying a circuit is, by default, a tensor-network contraction, and the *order* of the pairwise contractions decides the cost: every step produces an intermediate tensor, the expense is set by the largest intermediate the order ever creates (in matrix-product language, the bond dimension it forces you to carry), and a poor order can be exponentially more expensive than a good one on the very same network. QF exposes that choice as a `Method` sub-option, `"Path" -> "Greedy"` or `"Optimal"`, or an explicit list you compute and inspect yourself with `GreedyContractionPath` / `OptimalContractionPath`. Run a recognizable algorithm through it: Grover's search for the single input a Boolean oracle accepts among $N = 2^8 = 256$ candidates. The named `"GroverPhase"` circuit is *one* Grover iteration $G = D\,O$, the diffusion $D$ after the phase oracle $O$; stack the peak number of them, $k = \big\lfloor\tfrac{\pi}{4}\sqrt N\big\rceil = 12$, on the uniform superposition and contract the whole search through a greedily chosen path. The search register size:
```wolfram
n = 8;
```
The marked input, the single accepted candidate among $2^8 = 256$:
```wolfram
marked = 200;
```
Build the Boolean oracle that accepts only that input:
```wolfram
oracle = BooleanFunction[2^marked, n];
```
One Grover iteration $G = D\,O$, the diffusion after the phase oracle:
```wolfram
step = QuantumCircuitOperator["GroverPhase"[oracle]]
```
The peak iteration count $k = \big\lfloor\tfrac{\pi}{4}\sqrt N\big\rceil$:
```wolfram
kopt = Round[Pi/(4 ArcSin[Sqrt[1/2^n]]) - 1/2];
```
The full search: a Hadamard layer, then $k$ stacked iterations:
```wolfram
grover = QuantumCircuitOperator[Join[Thread["H" -> Range[n]], ConstantArray[step, kopt]]]
```
Draw the flattened circuit, gate labels off since the diffusion block repeats $k$ times:
```wolfram
grover["Flatten"]["Diagram", "ShowGateLabels" -> False]
```
Apply it to $|0\dots0\rangle$ through a greedily chosen contraction path:
```wolfram
found = grover[Method -> {"TensorNetwork", "Path" -> "Greedy"}]
```
Read the largest outcome probability:
```wolfram
Max[found["Probabilities"]]
```
Because the contraction stayed exact, the probability comes back as an exact rational, a ratio of large integers. Read it as a decimal:
```wolfram
% // N
```
The success probability is all but certain, sitting on the analytic Grover peak $\sin^2\!\big((2k+1)\arcsin\sqrt{1/N}\big)$ that the iteration count $k = 12$ is chosen to hit; the most probable bitstring is the marked input itself, $200 = 11001000$ (verified). A search of a few hundred gates, contracted in one steered pass.

The contraction plan is itself an inspectable object. Pull the circuit's tensor network:
```wolfram
net = grover["TensorNetwork"]
```
Compute the greedy path and draw it as a tree, its branching the order in which the network is folded together:
```wolfram
ContractionTree[
 TensorNetworkContraction[net, GreedyContractionPath[net]], 
 TreeElementLabel -> (All -> None), 
 TreeElementShapeFunction -> (All -> None)]
```
The tree is the engine's execution plan: the leaves are the circuit's gate and state tensors, every internal node is one pairwise contraction, and the branching structure is the schedule the greedy path chose, which contractions nest inside which. The depth and balance of that tree are what the path tunes to keep each intermediate tensor small.

### The stabilizer engine: the state as its symmetry group
For Clifford circuits there is a representation that never touches amplitudes at all: store the state as the group of Pauli operators that fix it, $\{P : P\,|\psi\rangle = +|\psi\rangle\}$, whose $n$ commuting generators determine the $n$-qubit state uniquely and fit in a binary symplectic tableau. This is the Gottesman-Knill representation, in which a Clifford gate is a constant-size rewrite of the tableau rather than a $2^n \times 2^n$ matrix product. The engine is one option away. Build a five-qubit linear cluster state, the entangled resource behind measurement-based computing, as a layer of Hadamards followed by nearest-neighbor controlled-$Z$ gates:
```wolfram
cluster5 = QuantumCircuitOperator[
   Join[Thread["H" -> Range[5]],
    {"CZ" -> {1, 2}, "CZ" -> {2, 3}, "CZ" -> {3, 4}, "CZ" -> {4, 5}}]]
```
Draw the circuit:
```wolfram
cluster5["Diagram"]
```
The diagram is an ordinary gate-model circuit: a Hadamard on every wire, then controlled-$Z$ gates down the chain $1$-$2$-$3$-$4$-$5$. Nothing about it is stabilizer-specific yet; the engine is a choice made at application time.

Apply it with the stabilizer engine instead of the default:
```wolfram
stab5 = cluster5[Method -> "Stabilizer"]
```
A purely Clifford circuit that would otherwise carry a $2^5$ amplitude vector returns instead a `PauliStabilizer`: the state stored as its symmetry group, not as its amplitudes.

Read the generators of that group:
```wolfram
stab5["PauliForm"]
```
Five Pauli strings that determine the state completely. They are exactly the cluster-state generators $K_v = X_v \prod_{w \sim v} Z_w$ along the chain, the $X$ on each site flanked by $Z$ on its neighbors.

Pauli expectations are closed-form algebra on the generators, with no state vector anywhere:
```wolfram
<|"XZIII" -> stab5["Expectation", "XZIII"], "IZXZI" -> stab5["Expectation", "IZXZI"], "XXIII" -> stab5["Expectation", "XXIII"]|>
```
The two generators $K_1$ and $K_3$ return deterministically $+1$, while $XX$ on the first two sites, which anticommutes with $K_1 = X_1 Z_2$, is exactly unbiased, all read off by $\mathbb{F}_2$ phase tracking.

The payoff of the representation is scale. Stream twenty thousand random Clifford gates onto a thousand qubits, a state whose amplitude vector would need $2^{1000}$ entries. The qubit count:
```wolfram
n = 1000;
```
Build the random gate stream, twenty thousand Hadamard, CNOT, and phase gates (seeded for reproducibility):
```wolfram
SeedRandom[1];
gates = Flatten @ Table[{"H" -> RandomInteger[{1, n}], "CNOT" -> RandomSample[Range[n], 2], "S" -> RandomInteger[{1, n}]}, {6666}];
```
Apply them to the $n$-qubit stabilizer state and time the evolution:
```wolfram
First @ AbsoluteTiming[big = PauliStabilizer[n]["ApplyCircuit", gates];]
```
The whole gate stream folds almost instantly: `ApplyCircuit` runs the tableau updates through one compiled kernel, so the per-gate cost is small and independent of the $2^n$ that never gets built. The bipartite entanglement entropy across the half-chain cut is then a closed form, $S(A) = \operatorname{rank}_{\mathbb{F}_2}(\text{generators restricted to } A) - |A|$, exact at any qubit count:
```wolfram
big["Entropy", Range[n/2]]
```
The entropy returns just shy of its half-chain maximum: the random circuit has driven the chain to near-maximal, Page-like entanglement, and whatever gate sequence the sampler draws, the value is an exact integer from an $\mathbb{F}_2$ rank, not an estimate. The formalism's boundary is first-class too: a non-Clifford gate ($T$, or a general phase) steps outside the stabilizer group, and QF answers with a `StabilizerFrame`, a superposition of stabilizer states, rather than a wrong answer.

### The representations convert exactly
The engines are not silos; one object model means the representations interconvert. Start from an ordinary three-qubit linear cluster circuit, a Hadamard on every wire followed by two nearest-neighbor $CZ$ bonds:
```wolfram
cluster3 = QuantumCircuitOperator[{"H" -> Range[3], "CZ", "CZ" -> {2, 3}}]
```
Draw it:
```wolfram
cluster3["Diagram"]
```
The bare `"CZ"` rides on the default first pair and `"CZ" -> {2, 3}` on the next, so the diagram shows the bond chain $1\!-\!2\!-\!3$: nothing here is stabilizer-specific, just gates.

Apply that same circuit through the stabilizer engine, which stores the state as its symmetry group rather than its amplitudes:
```wolfram
stab3 = cluster3[Method -> "Stabilizer"]
```
A `PauliStabilizer` returns, the state held as the three cluster generators that fix it. Now expand that tableau back into an exact state, that is, solve $S_i|\psi\rangle = |\psi\rangle$ for the unique joint $+1$ eigenstate of the three generators:
```wolfram
QuantumState[stab3] // TraditionalForm
```
The ket returns as $\tfrac{1}{\sqrt8}\sum_{x}(-1)^{x_1 x_2 + x_2 x_3}\,|x_1 x_2 x_3\rangle$: a uniform superposition over all eight bitstrings whose sign is fixed by the two cluster bonds, negative exactly on $011$ and $110$, the tableau expanded into exact amplitudes. The agreement with the default engine is not approximate; check it with `SameQ`, the language's structural identity:
```wolfram
QuantumState[stab3]["StateVector"] === cluster3[]["StateVector"]
```
`True`: the two engines produced the *same exact expression*, global phase included.

The conversion also runs back to circuits: the tableau synthesizes a Clifford circuit that prepares it.
```wolfram
stab3["Circuit"]["Diagram"]
```
A five-gate Clifford circuit returns, interleaved Hadamards and controlled-NOTs, and applying it to $|000\rangle$ reproduces the state vector exactly (verified). Circuit to tableau, tableau to state, tableau back to circuit: three representations of one object, each available where it is cheapest, with symbols, structure, or scale deciding which engine runs.

---

## Claim 5: A quantum object is an ordinary Wolfram Language expression

**Key feature:** this is the substrate under everything above, the reason Claim 1's algebra came out symbolic and Claim 4's three engines share one model: a QF object is a symbolic expression like any other in the language, so the kernel's own functions operate *on* the object and hand back another object. `Exp` exponentiates a generator into a propagator, `D` differentiates a parametrized gate into another operator, `Plus` and `Commutator` build a Hamiltonian from named Paulis, and the closed-form scalars these objects yield flow on into solvers like `Maximize` and `Reduce`. Nothing is converted, extracted, or exported first: the same `Exp`, `D`, `Plus` you would call on a number-matrix act unchanged on a `QuantumOperator` and return a `QuantumOperator`. The overload set is real but finite, `Exp` (and `MatrixExp`), `D`, `Plus`, `Times`, `Power`, `Tr`, `Commutator` and the analytic functions act on operators; for `Eigenvalues` or `Det` you still reach through `obj["Matrix"]`.

### Superposition by native arithmetic
Write a state the way a physicist writes it, a sum of kets scaled by amplitudes, and check it against the named one:
```wolfram
(QuantumState["0"] + QuantumState["1"])/Sqrt[2] == QuantumState["+"]
```
The comparison returns `True`: `Plus` and the scalar division act on the state objects, assembling $(|0\rangle+|1\rangle)/\sqrt2$ as a genuine `QuantumState`, and QF's overloaded `==` confirms it is $|+\rangle$. No amplitudes were pulled out by hand, and the same arithmetic scales up to the heavier overloads below.

### The exponential map: a generator becomes a propagator
The propagator of a Hamiltonian $H$ is the matrix exponential $e^{-iHt}$. Write it exactly the way a physicist does, with the kernel's `Exp` acting directly on the generator object:
```wolfram
uexp = Exp[(-I t) QuantumOperator["X"]]
```
Read its label and its matrix (trig-simplified):
```wolfram
{uexp["Label"], 
 MatrixForm@FullSimplify @ ExpToTrig @ Normal @ uexp["Matrix"]}
```
The result is a new `QuantumOperator` whose label reads $e^{-iXt}$ and whose matrix is the $X$-rotation $\cos t\,I - i\sin t\,X = \bigl(\begin{smallmatrix}\cos t & -i\sin t\\ -i\sin t & \cos t\end{smallmatrix}\bigr)$. `Exp` performed a genuine matrix exponential on the operator, mixing the two levels through the off-diagonal $-i\sin t$ rather than exponentiating entries one by one, and the generator was never lowered to a number-matrix first. (`MatrixExp[...]` returns the identical operator.)

### Operator algebra in one line
Ordinary `Plus` and scalar multiplication assemble a Hamiltonian on the object, and the Pauli decomposition reads the coefficients straight back out:
```wolfram
hab = a QuantumOperator["X"] + b QuantumOperator["Z"]
```
Read its label and its Pauli decomposition:
```wolfram
{hab["Label"], hab["PauliDecompose"]}
```
The operator carries the human-readable label $a\,X + b\,Z$, and `hab["PauliDecompose"]` returns `<|X -> a, Z -> b|>`: a round-trip through the algebra, coefficients in by arithmetic and back out by decomposition, on a symbolic operator rather than a numeric lookup table. The algebra runs deeper than arithmetic. The kernel's `Commutator` delivers the $su(2)$ relation as an operator identity:
```wolfram
Commutator[QuantumOperator["X"], QuantumOperator["Y"]] == 2 I QuantumOperator["Z"]
```
`True`: $[X, Y] = 2iZ$, computed on the operators and machine-checked against $2iZ$.

### Operator calculus: recover $H(t)$ from a target $U(t)$
Claim 1 solved Hamiltonians forward into their propagators; the inverse is pure operator calculus. Choose a unitary trajectory $U(t)$, and the Hamiltonian that drives it is $H(t) = i\,\dot U(t)\,U^\dagger(t)$, the engine behind counterdiabatic driving and shortcuts to adiabaticity. Compose the target $U(t) = R_Z(t^2)\,R_X(t)$ from named gates:
```wolfram
utgt = QuantumOperator["RZ"[t^2]] @ QuantumOperator["RX"[t]]
```
Differentiate it *as an operator* with the kernel's `D` and close it with the operator's own `Dagger`, forming $H = i\,\dot U\,U^\dagger$:
```wolfram
hdrive = I D[utgt, t] @ utgt["Dagger"]
```
Read the control field off each axis:
```wolfram
FullSimplify[ComplexExpand[#]] & /@ hdrive["PauliDecompose"]
```
The controls return axis by axis as the drive $H(t) = t\,Z + \tfrac12\big(\cos t^2\, X + \sin t^2\, Y\big)$, with a vanishing identity component, traceless as a control should be. `D` operated on the `QuantumOperator` and returned a `QuantumOperator` (not a bare matrix), and composition `@`, `Dagger`, and the Pauli readout chained on the same object; the defining identity $i\,\dot U = H\,U$ holds identically. The textbook inverse-design formula runs on the operator itself, with no differential solver involved.

### Calculus on states too
The derivative overload is not operator-only. Take a state whose amplitudes are functions of time. Build it:
```wolfram
psi = QuantumState[{Cos[t], Sin[t]}]
```
Differentiate it twice in place and add it back:
```wolfram
Simplify[D[psi, {t, 2}]["StateVector"] + psi["StateVector"]] // Normal
```
The sum collapses to the zero vector: the object satisfies the harmonic identity $|\psi''\rangle = -|\psi\rangle$ through WL's own calculus, the second derivative taken on the `QuantumState` directly.

### The scalars flow into solvers
When a QF object yields a symbolic scalar, the language's solvers carry on from there. The expectation of $X + Z$ in the trial state $R_Y(\theta)|0\rangle$ is a closed form. Build the trial state:
```wolfram
psi = QuantumOperator["RY"[\[Theta]]][QuantumState["0"]]
```
Form the energy $\langle\psi|(X+Z)|\psi\rangle$:
```wolfram
energy = Simplify[Conjugate[psi["StateVector"]] . ((QuantumOperator["X"] + QuantumOperator["Z"])[psi])["StateVector"], \[Theta] \[Element] Reals];
```
Hand it to `Maximize`:
```wolfram
Assuming[\[Theta] \[Element] Reals, Maximize[{energy, 0 <= \[Theta] <= 2 Pi}, \[Theta]] // FullSimplify]
```
The energy is $\cos\theta + \sin\theta$, and `Maximize` returns $\sqrt2$ at $\theta = \pi/4$: the exact extremal value, recognizable as the operator norm of $X + Z$, a complete variational eigenvalue calculation done by composing objects with native arithmetic and optimization. The same flow certifies a physical threshold. The Werner state's realignment criterion, a quantity positive exactly when it witnesses entanglement, is a symbolic function of the mixing $p$:
```wolfram
crit = Simplify[QuantumEntanglementMonotone[QuantumState["Werner"[p, 2]], "Realignment"], 0 < p < 1];
```
Turn the criterion into the exact entanglement threshold:
```wolfram
Reduce[crit > 0 && 0 < p < 1, p]
```
The criterion is $\tfrac{|3-4p|-1}{2}$, and `Reduce` turns it into the exact threshold: the Werner state is entangled precisely for $0 < p < \tfrac{1}{2}$, a physical fact certified by handing a QF-derived expression to a native solver.

---

## Claim 6: Every object draws itself

**Key feature:** Claim 5 let the kernel's functions read an object; here the object draws *itself*. Visualization is property dispatch, not a separate plotting library: each object exposes plot properties that render through the Wolfram Language graphics stack and forward its options, so a publication figure is one call on the object it depicts. The catalog is wide. States draw as Bloch vectors, amplitude and probability charts, and qudit pie and sector wheels; measurements as outcome histograms; stabilizer states as their binary symplectic tableau; circuits as gate diagrams, tensor-network graphs, and multiway and causal views. Adiabatic evolutions add spectrum, gap, and coupling plots (in the `QuantumOptimization` sub-package), the underlying `TensorNetwork` renders as a hypergraph through the TensorNetworks paclet, and the phase-space objects of Claim 2 plot as Wigner surfaces. Below, one example per object kind.

### A state on the Bloch sphere
A single qubit draws its own Bloch vector:
```wolfram
QuantumState[{Cos[Pi/5], Sin[Pi/5] Exp[I Pi/3]}]["BlochPlot"]
```
The output is a `Graphics3D` of the sphere with the state's arrow at polar angle $2\pi/5$ and azimuth $\pi/3$: the geometry of the state, drawn by the state.

### Amplitudes that keep their phases
A bar chart of probabilities discards the phases; `AmplitudesChart` draws the complex amplitudes $\psi_x = \langle x|\psi\rangle$ themselves, heights for the magnitudes $|\psi_x|$ and colors for the phases $\arg\psi_x$. Feed a quantum Fourier transform the two-term superposition $\tfrac{1}{\sqrt2}(|000\rangle + |101\rangle)$, which the transform spreads into a phase-rich pattern across all eight basis states:
```wolfram
QuantumCircuitOperator["Fourier"[3]][
  QuantumState[Normalize[{1, 0, 0, 0, 0, 1, 0, 0}]]]["AmplitudesChart", 
 ChartLegends -> Automatic, AspectRatio -> 1/3, 
 PlotRange -> {Automatic, 1}]
```
The eight bars return with magnitudes tracing an interference comb, a clean null at $|100\rangle$, and the surviving bars tinted by a wheel of distinct phase colors keyed by the legend: the QFT's phase ramp made visible alongside the magnitudes, the full complex state in one figure.

### A measurement draws its statistics
A `QuantumMeasurement` object charts its own outcome distribution, the Born probabilities $p(x) = |\langle x|\psi\rangle|^2$. Run the Bernstein-Vazirani circuit for a hidden string $s = 1011$, a single oracle query that encodes $s$ in the interference pattern, and let it draw the statistics of its own final measurement:
```wolfram
QuantumCircuitOperator[
   "BernsteinVazirani"[{1, 0, 1, 1}]][]["ProbabilitiesPlot", 
 "LabelsAngle" -> 90 Degree, AspectRatio -> 1/3]
```
The histogram collapses to a single bar of height $1$ at $1011$, every other outcome exactly zero: one query has pulled the entire hidden string out of the Born statistics, the deterministic signature of the algorithm. The hardware result in Claim 7 is the same kind of object, which is why simulated and measured statistics compare so directly there.

### A stabilizer state shows its tableau
The stabilizer states of Claim 4 have a native picture of their own: the binary symplectic tableau, the actual data structure the formalism computes on. Feed the five-cycle graph $C_5$ straight into `PauliStabilizer` to build its ring cluster state, and draw the tableau:
```wolfram
PauliStabilizer[QuantumCircuitOperator["Graph"[CycleGraph[5]]]]["TableauForm"]
```
The display is the sign column beside a ten-by-ten binary grid, dividers separating the $x$ bits from the $z$ bits and the destabilizer rows from the stabilizer rows; the stabilizer rows are the five ring generators $K_v = X_v \prod_{w \sim v} Z_w$, so the $z$ block is literally the adjacency matrix of the cycle while the $x$ block is the identity. This view *is* the algebra: the entanglement geometry of the ring is drawn in the binary pattern, and a Clifford gate acts by rewriting a few columns of precisely this grid.

### The circuit, and the network behind it
A circuit draws itself as a gate diagram. Take quantum phase estimation, the subroutine inside Shor's factoring and the HHL linear solver, on three counting qubits:
```wolfram
qpe = QuantumCircuitOperator["PhaseEstimation"[3]]
```
Draw it as a gate diagram:
```wolfram
qpe["Diagram"]
```
The wire-and-gate picture lays out the algorithm: Hadamards on the counting register, the controlled powers of the unitary, and the inverse Fourier transform that reads the phase back. The same object also draws its computational anatomy, the tensor network the engine contracts:
```wolfram
qpe["TensorNetworkGraph"]
```
A `Graph` whose vertices are the gate and state tensors and whose edges are the contracted indices: the diagram shows the physics, the network shows the computation, and Claim 4 steered the order in which such a network is contracted.

The same network renders in the dual picture too, as a hypergraph in which the roles swap, every index a vertex and every tensor one hyperedge spanning its indices:
```wolfram
qpe["TensorNetwork"]["Hypergraph"]
```
The result is a `Hypergraph` with one hyperedge per tensor of the phase-estimation circuit, drawn over the network's indices: the tensor-contraction picture of the algorithm, in which contracting two tensors is the merging of two hyperedges.

---

## Claim 7: An interoperability hub, down to real hardware

**Key feature:** the objects so far stayed inside QF; this claim sends them out. A QF object can *leave* QF for the rest of the quantum ecosystem, and foreign circuits can *enter* it; `QuantumQASM` is the single OpenQASM door, and it runs in pure Wolfram Language with no external account. The door opens wider at each step below: text export, foreign import, device-aware transpilation, and finally a physical quantum processor.

### Export to OpenQASM 3.0 (native, no Python)
`QuantumQASM` turns a circuit into the standard text format read by Qiskit, Cirq, and most hardware stacks.
```wolfram
QuantumQASM[QuantumCircuitOperator["Toffoli"]]
```
The output is a complete OpenQASM 3.0 program for the three-qubit Toffoli gate: rather than a black-box `ccx`, the emitter writes the standard Clifford+$T$ decomposition, six CNOTs (each a `ctrl @ U`) interleaved with the non-Clifford $T$ and $T^\dagger$ phases and a pair of Hadamards bracketing the target, every gate a concrete `U` rotation. The native emitter writes OpenQASM 3.0, and the importer below reads OpenQASM 2.0 or 3.0 alike. A 2.0 *dump* is one option away through the qiskit path, `QuantumQASM[circuit, "Version" -> 2]` (that route needs Python; the 3.0 export and all imports do not). Export needs concrete angles: a circuit with a symbolic parameter returns a `Failure` rather than a QASM string, so resolve parameters to numbers at this boundary while keeping them symbolic everywhere else in QF (Claims 1 and 5).

### Import foreign OpenQASM and compute with it (native, no Python)
The same door runs the other way: hand `QuantumCircuitOperator` an OpenQASM string authored anywhere and get back a QF circuit, parsed in pure WL, `reset`, controlled gates, mid-circuit `measure` and all. The importer reads each gate's angle as a Wolfram Language expression, so the angles can stay *symbolic*: write `Pi/2` and `Pi`, the exact constant, in place of their decimal expansions, and the imported circuit is exact rather than a floating-point copy. Import a four-qubit circuit that runs the optimal quantum strategy for the CHSH game:
```wolfram
chsh = QuantumCircuitOperator["OPENQASM 3.0;
qubit[4] q;
bit[4] c;
reset q[0];
reset q[3];
U(Pi/2, 0, Pi) q[3];
U(Pi/2, 0, Pi) q[0];
swap q[0] q[3];
ctrl(0) @ negctrl(1) @ U(Pi/2, Pi, Pi) q[0] q[3];
ctrl(1) @ negctrl(0) @ U(Pi/2, 0, 0) q[0] q[3];
reset q[1];
U(Pi/2, 0, Pi) q[1];
reset q[2];
U(Pi/2, 0, Pi) q[2];
U(0, 0, 0) q[0];
ctrl(1) @ negctrl(0) @ U(Pi/2, 0, Pi) q[1] q[0];
U(Pi/4, 0, 0) q[3];
ctrl(0) @ negctrl(1) @ U(Pi/2, 0, Pi) q[2] q[3];
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];"]
```
Note that QF reads `U(Pi/2, 0, Pi)`: the angle is the Wolfram constant `Pi`. A spec-strict OpenQASM 3 parser like qiskit knows only the lowercase `pi` and rejects `Pi` as an undefined symbol.

Once imported it is an ordinary QF object, exact end to end because every angle stayed a symbol.

In the CHSH game the players win when their answers satisfy $x \wedge y = a \oplus b$, and no classical strategy wins more than $\tfrac34$ of the time; entanglement beats that bound. Read the circuit's joint measurement distribution, label the four outcomes $(a, x, y, b)$, and ask Wolfram Language for the winning probability $\Pr[\,x \wedge y = a \oplus b\,]$, in closed form:
```wolfram
Probability[BitAnd[x, y] == BitXor[a, b], {a, x, y, b} \[Distributed] chsh[]["MultivariateDistribution"]] // FullSimplify
```
Because the angles never left exact arithmetic, the winning probability returns as a closed form, $\tfrac{2+\sqrt2}{4} = \cos^2(\pi/8)$, exactly the **Tsirelson bound**, the quantum maximum for CHSH and the algebraic number a sampled estimate only converges toward. A circuit authored in another tool became an exact probability we read off with native WL statistics.

### A hardware-aware transpilation target
`QiskitTarget` carries a device spec (qubit count, native gate set, connectivity), validated at construction against qiskit's own `Target` schema through the Python bridge. Build a 5-qubit linear-chain device:
```wolfram
target = QiskitTarget[<|"NumQubits" -> 5, "BasisGates" -> {"cx", "rz", "sx", "x"},
   "CouplingMap" -> {{0,1},{1,2},{2,3},{3,4}}|>]
```
Read its coupling map back:
```wolfram
target["CouplingMap"]
```
Hand the target to `QuantumQASM` as an option, and the circuit is transpiled to that backend's native gates and connectivity:
```wolfram
QuantumQASM[QuantumCircuitOperator["Fourier"[3]], "Target" -> target]
```
The QFT comes back rewritten for the device: its controlled-phase rotations decomposed into the native `rz`/`sx`/`x` set, the entangling steps as `cx` on physical qubits, every two-qubit gate routed along an edge the coupling map allows. Transpiling a Fourier transform, with its all-to-all controlled phases, is the realistic version of this task, where decomposition and routing both have real work to do.

### Submit to a real IBM QPU, and compare to the exact result
An IBM account connects through `ServiceConnect["IBMQuantumPlatform"]`, after which `IBMJobSubmit` sends a circuit to the hardware asynchronously and returns an `IBMJob` handle. Create an API key and copy your instance's cloud resource name (CRN) from your IBM Quantum Platform account, then connect, the placeholder `"xxx"` values standing in for your own:
```wolfram
api = "xxx";
LocalSymbol["ibm_crn"] = "xxx";
ServiceConnect["IBMQuantumPlatform", "New", Authentication -> {"apikey" -> api}]
```
The call returns a `ServiceObject` connection handle. A QPU job needs a circuit that *ends in a measurement*; in QF a list of qubit indices placed in the circuit is exactly that, so `{1, 2, 3}` measures all three qubits. Build the measured GHZ circuit:
```wolfram
ghz = QuantumCircuitOperator[{"GHZ"[3], {1, 2, 3}}]
```
Submit it to the backend, which returns a handle immediately:
```wolfram
job = IBMJobSubmit[ghz, "ibm_fez", "Shots" -> 4096]
```
Check the status:
```wolfram
job["Status"]
```
The handle returns immediately with status `"Queued"` while the job waits in line on the device. The same circuit runs exactly in the Wolfram Language by just applying it:
```wolfram
wl = ghz[]
```
The result is a `QuantumMeasurement` with probabilities $\tfrac12$ at $000$ and $\tfrac12$ at $111$. Once the job completes, the hardware outcome comes back as the **same kind of `QuantumMeasurement` object**, decoded into the same qubit order:
```wolfram
qpu = job["Refresh"][]
```
This run returned $4096$ shots from `ibm_fez`. Because the two results are the same kind of object, hardware noise is one chart away, the exact and measured probabilities paired per basis state:
```wolfram
BarChart[Transpose[Values @* KeySort /@ {wl["Probabilities"], qpu["Probabilities"]}],
  AspectRatio -> 1/2, Frame -> True, ChartLegends -> {"Exact WL", "QPU"},
  PlotLabels -> Keys[wl["Probabilities"]]]
```
The exact bars are $\tfrac12$ at $000$ and $111$ and zero elsewhere; the device's two tall bars land visibly short of $\tfrac12$, with the missing weight spread thinly over the six bitstrings the exact state forbids. That gap, read bar by bar off one figure, *is* the device noise.

---

## Where This Leaves Us

Seven claims, each cashed out in executed code. Symbolic algebra carried dynamics to closed form: the Rabi propagator of a rotating drive, a Landau-Zener sweep expressed in parabolic cylinder functions, a QAOA landscape whose exact maximum is a closed-form cubic root, the entanglement entropy of a quench as a function of time. One object model held the rest of the theory: channels, POVMs, Fock space, and phase space; a master equation solved once for every initial state and every rate, with $T_2 \le 2T_1$ arriving as a theorem rather than a remembered rule; a basis and a dimension on every object, a spin-1 Hamiltonian and a five-site ring each diagonalized by its own named basis; a state reconstructed from measured counts at a fidelity a hair below unity. Measurement stayed inside the physics: records became quantum wires, teleportation ran without a classical channel, a coherent error was digitized and corrected at codeword fidelity exactly $1$, and a measurement ran backwards. Underneath, three engines computed the same objects: exact state-vector algebra, a steered tensor-network contraction, and a stabilizer tableau that carried a thousand qubits to near-maximal half-chain entanglement computed exactly, the representations converting into one another exactly. Because these objects are ordinary Wolfram Language expressions, the kernel's own functions transformed them directly, exponentiating a generator into its propagator and differentiating a gate to recover the Hamiltonian that drives it, and each of them drew itself. And the same objects left the system: out to OpenQASM and back, through a transpiler to a device's native gates, and onto a physical quantum processor, whose $4096$-shot histogram stands in this document next to the exact distribution it approximates.

Behind all seven claims is one mechanism: quantum objects that are exact symbolic expressions compose, with each other and with the entire language. One consequence ran through every example: the deliverable was *exact*, the algebraic QAOA optimum, the closed-form relaxation law, the integer stabilizer entropy, the exactly-recovered teleported state, never a value sampled near the answer but the ground truth such a sample converges toward. That is the quiet numerical claim under the symbolic one: where a numeric simulation or a finite measurement hands you an estimate, QF hands you the exact object it estimates. The distance from writing a Hamiltonian to reading its closed form, or to running it on hardware, is a handful of composable calls. A result that is a formula can be differentiated, reduced, or pushed to a limit; a result that is a circuit is already in the format the rest of the quantum ecosystem speaks.

The next steps are the reader's: change the graph under the QAOA section, raise the spin past $1$, swap damping for dephasing, tilt the repetition code's coherent error off the $X$ axis, point the last section at a different backend. The document re-derives itself.

---

**Setup.** Evaluate once before the cells; this installs the development builds of the QuantumFramework and TensorNetworks paclets and loads them:

```wolfram
PacletInstall["https://www.wolfr.am/DevWQCF", 
 ForceVersionInstall -> True];
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
PacletInstall["https://www.wolfr.am/DevWTN", 
 ForceVersionInstall -> True];
Needs["Wolfram`TensorNetworks`"];
```

Every cell was executed as shown (Wolfram Language 15.0, `Wolfram/QuantumFramework` 2.0.0) and every stated output is the captured result, including the hardware run on `ibm_fez`; reproducing the transpilation cells needs a Python installation with qiskit, and the hardware cells an IBM Quantum account.

Taken together, these features, and above all the way symbolic algebra, numerical simulation, and visualization live inside one object, make QF a particularly effective tool for teaching and learning quantum theory at every level, from a first qubit to graduate research: a student can build a state, evolve it, draw it, and read its closed form in a single expression, the computation-first habit this document follows throughout. Wolfram develops that direction directly, in the book [*Quantum Computing by Model and Simulation: A Computational Approach*](https://www.wolfram-media.com/products/quantum-computing-by-model-and-simulation/) (Wolfram Media), which teaches the subject through modeling and simulation in the framework, and in [Wolfram U](https://www.wolfram.com/wolfram-u/), whose free interactive courses and study groups include quantum computing and related topics.
