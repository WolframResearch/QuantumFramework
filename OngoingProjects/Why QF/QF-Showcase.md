# What Makes Wolfram QuantumFramework Unique: A Verified Showcase

Wolfram QuantumFramework is a quantum laboratory built into the Wolfram Language. The objects of quantum theory, states, operators, channels, measurements, circuits, are exact symbolic expressions: they carry free parameters, compose with one another in closed form, and answer questions about themselves through one uniform interface. A computation here typically ends in a formula, the propagator $U(t)$ of a time-dependent Hamiltonian, an optimization landscape you can differentiate, an entanglement entropy $S(t)$ you can solve, rather than a table of numbers. The same objects run in any finite dimension, qubits, qutrits, spin-$j$ multiplets, $N$-level systems, and they reach outward as well, to OpenQASM, to device transpilation, to real quantum hardware. Its natural habitat is the computation whose deliverable is insight: a closed form, an exact threshold, a verified physical statement.

This document makes that concrete as seven claims, each developed through worked examples:

1. **Symbolic-first, exact algebra by default.** Time-dependent Hamiltonians solved to closed-form propagators, a control Hamiltonian recovered from a target evolution, an exact QAOA landscape, a many-body entanglement entropy as a formula in time.
2. **One object model spans all of quantum theory.** Channels, POVMs, bosonic Fock space, and phase-space quasi-probabilities are the same kind of object as a gate, and a master equation solves symbolically for every state and every rate at once; the model stores the basis on every object, runs in any finite dimension, and closes with a full state-tomography loop.
3. **Measurement records are quantum wires.** A measurement appends a pointer qudit instead of collapsing the state, so feedforward is a controlled gate on a record: teleportation runs without a classical channel, a coherent error is digitized and corrected at codeword fidelity exactly $1$, and a measurement can even be run backwards.
4. **One object model, three computational engines.** The same circuit runs as exact state-vector algebra, as a steered tensor-network contraction, or as a compiled stabilizer tableau that carries a thousand qubits with ease, and the representations convert into one another exactly.
5. **A quantum object is an ordinary Wolfram Language expression.** Built-ins like `Plot` and `Reduce` apply to it directly, precisely because nothing about it is special.
6. **Every object draws itself.** Bloch spheres, amplitude charts that carry the phases, measurement histograms, circuit diagrams, and the tensor network behind a circuit, each a one-call property of the object it depicts.
7. **An interoperability hub, down to real hardware.** OpenQASM in and out natively; transpilation to a device target; submission to a real QPU, with IBM as the worked example.

The claims build on one another, symbolic objects at the base and real hardware at the summit, so we begin with the algebra everything else inherits.

---

## Claim 1: Symbolic-first, exact algebra by default

**Key feature:** symbolic algebra pervades the *whole* object model, not one function. The building blocks you reach for, a **state** (`QuantumState`), a **Hamiltonian** (`QuantumOperator`), its **time evolution** (`QuantumEvolve`), and a **circuit** (`QuantumCircuitOperator`), are each fully symbolic, and they compose into exact closed-form results. We take them in that order, starting with the simplest: a state whose amplitudes are symbols.

### A symbolic state
A `QuantumState` can carry free symbols. The object to build is the general qubit $|\psi\rangle = \cos\alpha\,|0\rangle + e^{i\beta}\sin\alpha\,|1\rangle$: not one state but the whole two-parameter family at once. Build it and read the amplitudes back:
```wolfram
psi = QuantumState[{Cos[\[Alpha]], Sin[\[Alpha]] Exp[I \[Beta]]}];
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
It is **explicitly time-dependent**: `\[FormalT]` is QF's time variable, and it appears inside the operator.
```wolfram
H = \[Omega]0/2 QuantumOperator["Z"] +
    \[CapitalOmega]/2 (Cos[\[Omega] \[FormalT]] QuantumOperator["X"] + Sin[\[Omega] \[FormalT]] QuantumOperator["Y"]);
```
Because the field direction turns, $H$ does **not commute with itself at different times**, so the propagator is not simply $e^{-i\int H\,dt}$; it needs genuine time-ordering. Check the non-commutativity directly, computing $[H(t_1), H(t_2)]$:
```wolfram
Hm[t_] := Normal[H["Matrix"]] /. \[FormalT] -> t;
Commutator[Hm[t1], Hm[t2], Dot] // FullSimplify
```
The commutator is nonzero, its diagonal carrying $\mp\tfrac{i}{2}\Omega^2\sin[(t_1-t_2)\omega]$: this is the genuinely hard, time-ordered case, not a disguised static problem.

### Its propagator, in closed form
`QuantumEvolve` integrates the time-dependent Schrodinger equation $i\,\partial_t U = H(t)\,U$ and returns the propagator $U(t)$ symbolically. (Verified here: $U(0)=I$ and $i\,U'(t)-H(t)\,U(t)=0$ exactly.)
```wolfram
U = QuantumEvolve[H, None];
```
Simplify its matrix under real parameters and $\Omega > 0$:
```wolfram
Assuming[{\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, \[FormalT] \[Element] Reals},
  FullSimplify[ComplexExpand[Normal[U["Matrix"]]]]] // MatrixForm
```
The matrix returns exactly as the lab-frame Rabi propagator, written here with the detuning $\delta \equiv \omega - \omega_0$ and the generalized Rabi frequency $\Omega_R \equiv \sqrt{\Omega^2 + \delta^2}$:
$$
U(t) = \begin{pmatrix}
e^{-i\omega t/2}\!\left(\cos\dfrac{\Omega_R t}{2} + i\dfrac{\delta}{\Omega_R}\sin\dfrac{\Omega_R t}{2}\right) & -\,i\,e^{-i\omega t/2}\,\dfrac{\Omega}{\Omega_R}\sin\dfrac{\Omega_R t}{2} \\[3mm]
-\,i\,e^{i\omega t/2}\,\dfrac{\Omega}{\Omega_R}\sin\dfrac{\Omega_R t}{2} & e^{i\omega t/2}\!\left(\cos\dfrac{\Omega_R t}{2} - i\dfrac{\delta}{\Omega_R}\sin\dfrac{\Omega_R t}{2}\right)
\end{pmatrix}
$$
The $e^{\mp i\omega t/2}$ factors are the rotating-frame phases; the $\cos$ and $\sin$ of $\Omega_R t/2$ are the dressed-state oscillations. The physically meaningful readout is the transition probability $|\langle 1|U(t)|0\rangle|^2$:
```wolfram
amp = Normal[U["Matrix"]][[2, 1]];
Assuming[{\[Omega] \[Element] Reals, \[Omega]0 \[Element] Reals, \[CapitalOmega] > 0, \[FormalT] \[Element] Reals},
  FullSimplify[ComplexExpand[Abs[amp]^2]]]
```
The result is the textbook Rabi formula, $P_{0\to1}(t) = \frac{\Omega^2}{\Omega^2 + \delta^2}\,\sin^2\!\big(\tfrac{t}{2}\sqrt{\Omega^2 + \delta^2}\big)$. On resonance ($\omega=\omega_0$, so $\delta=0$) this is $\sin^2(\Omega t/2)$: complete Rabi flopping between $|0\rangle$ and $|1\rangle$. The oscillation rate $\sqrt{\Omega^2+\delta^2}$ is the generalized Rabi frequency, the dressed-state gap in the rotating frame, here recovered as a closed-form function of the drive detuning and strength.

### When the closed form needs special functions: Landau-Zener
The rotating-field propagator was elementary, but symbolic evolution is not confined to elementary functions. Sweep the energy bias *linearly* through an avoided crossing at rate $v$ with a constant gap $\Delta$, the Landau-Zener model $H(t) = \tfrac{v t}{2}\,Z + \tfrac{\Delta}{2}\,X$, the universal question of whether a system dragged through a level crossing stays in its energy branch or jumps the gap, again time-dependent and non-commuting:
```wolfram
HLZ = v \[FormalT]/2 QuantumOperator["Z"] + \[CapitalDelta]/2 QuantumOperator["X"];
```
Ask `QuantumEvolve` for the propagator and simplify its full matrix:
```wolfram
ULZ = QuantumEvolve[HLZ, None];
FullSimplify[Normal[ULZ["Matrix"]]]
```
The matrix returns as a genuinely special-function object. Its first entry alone reads $e^{-i v t^2/4}\,{}_1F_1\!\big(\tfrac{i\Delta^2}{8v};\ \tfrac12;\ \tfrac{i}{2}\,v t^2\big)$, and across the four entries Kummer functions ${}_1F_1$, Hermite functions of fractional order, and parabolic cylinder functions $D_\nu$ appear side by side, glued by $\Gamma$-function coefficients. The cylinder functions carry purely imaginary indices $\nu = \pm\,i\,\Delta^2/4v$ and complex-rotated time arguments $e^{i\pi/4}\sqrt{v}\,t$ and $e^{3i\pi/4}\sqrt{v}\,t$: the hallmark of the Landau-Zener problem, with the single dimensionless combination $\Delta^2/v$ controlling every index, coefficient, and phase. And it is one function family written three ways: every piece solves the same Weber equation and interconverts by standard identities, here in the raw form the symbolic solver found it. The matrix checks out as a genuine propagator (numerically $U(0)=I$, solves the Schrodinger equation, and is unitary), and in the long-time limit its asymptotics collapse to the famous Landau-Zener transition probability $P = e^{-\pi\Delta^2/2v}$. This is symbolic evolution reaching past elementary functions to the special-function solution the physics demands.

### Designing the dynamics: recover $H(t)$ from a target $U(t)$
The previous examples *solved* given Hamiltonians. The inverse is just as symbolic: choose a unitary trajectory $U(t)$, and the Hamiltonian that drives it is $H(t) = i\,\dot U(t)\,U^\dagger(t)$, the engine behind counterdiabatic driving and shortcuts to adiabaticity. Compose the target trajectory $U(t) = R_Z(t^2)\,R_X(t)$ from QF's named parametric gates, a quadratically chirped $Z$ rotation following an $X$ rotation:
```wolfram
U = QuantumOperator["RZ"[\[FormalT]^2]] @ QuantumOperator["RX"[\[FormalT]]];
```
Differentiate the trajectory *as an operator*, `D` acts on the parametric `QuantumOperator` directly, close it with the operator's own `Dagger`, and read the control field off each axis with the Pauli decomposition:
```wolfram
H = I D[U, \[FormalT]] @ U["Dagger"];
FullSimplify[ComplexExpand[#]] & /@ H["PauliDecompose"]
```
The controls return axis by axis: a linearly growing field $t$ on $Z$, the quadrature pair $\tfrac12\cos t^2$ and $\tfrac12\sin t^2$ on $X$ and $Y$, and a vanishing identity component, traceless as a control should be. Together they are the drive $H(t) = t\,Z + \tfrac12\big(\cos t^2\, X + \sin t^2\, Y\big)$ that realizes the chosen evolution, and $i\,\dot U = H\,U$ holds identically (verified). You specify the evolution gate by gate; the algebra returns the fields to apply on each axis, a closed-form inverse with no differential solver involved.

### A parametric circuit: QAOA, built entirely from a graph
QAOA is the canonical parametric circuit for optimization: one layer alternates a cost unitary $e^{-i\gamma H_C}$ with a mixer $e^{-i\theta B}$, and we tune the two angles to maximize a classical objective. The objective here is MaxCut: color the vertices of a graph with two colors $s_i = \pm1$ so that as many edges as possible join different colors, $\max_{s\,\in\,\{\pm1\}^n} \sum_{(i,j)\in E} \tfrac{1 - s_i s_j}{2}$. Everything below is *derived from the problem graph*, so nothing is hardcoded. Take the triangle, a frustrated instance where you cannot cut all three edges, and read off the classical optimum:
```wolfram
g = CompleteGraph[3];
FindMaximumCut[g]
```
The classical answer is a cut of $2$ of the $3$ edges, one vertex against the other two. Record the vertex count for reuse:
```wolfram
n = VertexCount[g]
```

The MaxCut cost Hamiltonian $H_C = \sum_{(i,j)}\tfrac{1 - Z_i Z_j}{2}$ is assembled straight from the edge list, one $(1-Z_iZ_j)/2$ term per edge, and carries a label that will name its block in the diagram below:
```wolfram
Hc = QuantumOperator[Total[With[{id = QuantumOperator[StringJoin[ConstantArray["I", n]]]/2},
   (id - QuantumOperator["ZZ"/2, #] &) /@ (List @@@ EdgeList[g])]], "Label" -> "Hamiltonian"];
```
Change `g` and this operator adapts to any graph.

The ansatz is built one labeled piece at a time. State prep puts every qubit in $|+\rangle$, an equal superposition over all cuts:
```wolfram
stateprep = QuantumCircuitOperator[{StringJoin[ConstantArray["0", n]], "H" -> Range[n]}, "Label" -> "Prep"];
```
The cost layer applies one $ZZ(\gamma)$ rotation per edge, carrying the angle $\gamma$:
```wolfram
qcost = QuantumCircuitOperator[("R"[\[Gamma], "ZZ"] -> # &) /@ (List @@@ EdgeList[g]), "Parameters" -> {\[Gamma]}, "Label" -> "Cost"];
```
The mixer applies an $R_X(\theta)$ to every vertex:
```wolfram
mixer = QuantumCircuitOperator[{"RX"[\[Theta]] -> VertexList[g]}, "Parameters" -> {\[Theta]}, "Label" -> "Mixer"];
```
Compose the three into one parametric layer:
```wolfram
oneLayer = QuantumCircuitOperator[{stateprep, qcost, mixer}, "Parameters" -> {\[Theta], \[Gamma]}];
```
Not a single qubit index was written by hand.

The cost landscape is the expectation $\langle\psi|H_C|\psi\rangle$, and here QF leans on a feature worth naming: **a `QuantumCircuitOperator` can hold states as elements**. Composing the prep circuit (the ket), the operator $H_C$, and the prep circuit's `["Dagger"]` (the bra) builds $\langle\psi|H_C|\psi\rangle$ as a *single object*, and that object draws itself:
```wolfram
meanvalue = (oneLayer /* Hc /* oneLayer["Dagger"]);
meanvalue["Diagram"]
```
The diagram shows the bra-ket sandwich block by block, the labeled `Prep`/`Cost`/`Mixer` layer, the `Hamiltonian` in the middle, and the daggered layer mirrored on the right: the expectation value as a picture. Calling the same object with `[]["Scalar"]` reads off the number, no manual inner product:
```wolfram
landscape = FullSimplify[meanvalue[]["Scalar"], Element[{\[Gamma], \[Theta]}, Reals]]
```
The exact QAOA energy landscape returns as $\langle H_C\rangle(\gamma,\theta) = \frac{3}{8}\big(3 + \cos 2\gamma + 2\cos 2\theta\,\sin^2\gamma - 2\sin 2\gamma\,\sin 2\theta\big)$, a closed form in the two angles.

Because the landscape is a formula, the whole optimization problem can be looked at before it is solved. Plot the surface over a full period of both angles:
```wolfram
Plot3D[landscape, {\[Gamma], 0, 2 \[Pi]}, {\[Theta], 0, 2 \[Pi]},
  PlotRange -> All,
  AxesLabel -> {"\[Gamma]", "\[Theta]", "cost"},
  ColorFunction -> "TemperatureMap",
  PlotPoints -> 50]
```
The surface is the variational problem in full view: smooth and doubly periodic with period $\pi$ in each angle (every term enters through $2\gamma$ and $2\theta$, so the $2\pi$ window tiles the unit cell four times over), and flat at the value $\tfrac32$ along the line $\gamma = 0$, where the cost layer is off and the equal superposition cuts $\tfrac32$ edges on average. What an experiment estimates point by point, with shot noise on every sample, is here one exact surface, and its peaks are not read off a grid but certified algebraically. Take the maximum exactly:
```wolfram
Maximize[landscape, {\[Theta], \[Gamma]}]
```
The maximum returns as exactly $2$, with the optimal angles as closed-form algebraic numbers, and it equals `FindMaximumCut[g]`: because the landscape is exact, the variational optimum *provably* equals the true cut value, the same way a symbolic VQE expectation lands exactly on a Hamiltonian's eigenvalue rather than merely near it.

### The payoff: a full many-body computation in closed form
The same machinery scales to a genuine many-body computation. Build the 2-qubit Heisenberg Hamiltonian $H = J\,(XX + YY + ZZ)$ symbolically (Hamiltonian simulation, the workhorse application of a quantum computer):
```wolfram
H = J (QuantumOperator["XX"] + QuantumOperator["YY"] + QuantumOperator["ZZ"]);
```
Quench the product state $|01\rangle$ under it, that is, compute $|\psi(t)\rangle = e^{-iHt}\,|01\rangle$, and display the state at time $t$:
```wolfram
st = QuantumEvolve[H, QuantumState["01"]];
FullSimplify[st] // TraditionalForm
```
The evolved ket displays in closed form: the excitation oscillates coherently between $|01\rangle$ and $|10\rangle$ at angular frequency $2J$, inside an overall phase $e^{iJt}$. Now read off the **entanglement entropy of one qubit as a function of time**, $S(t) = -\operatorname{Tr}\big[\rho_1\log_2\rho_1\big]$ with $\rho_1 = \operatorname{Tr}_2\,|\psi(t)\rangle\langle\psi(t)|$, the standard measure of how entangled the two spins have become, and a derived, nonlinear quantity:
```wolfram
s = FullSimplify[QuantumPartialTrace[st, {2}]["VonNeumannEntropy"], \[FormalT] J \[Element] Reals]
```
QF returns the **entire curve as a formula in $t$ and $J$**: the binary entropy of the two oscillating populations, zero whenever the state is a product, a full bit at the Bell points, periodic with the exchange oscillation. Because it is a formula, it can be differentiated, expanded, or solved in closed form, and the same approach extends to larger symbolic spin chains.

Plot it (at $J = 1$):
```wolfram
Plot[s /. {J -> 1, \[FormalT] -> t}, {t, 0, Pi},
  Frame -> True, FrameLabel -> {"time  J t", "entanglement entropy  S(t)  [bits]"}, PlotRange -> {0, 1.05}]
```

---

## Claim 2: One object model spans all of quantum theory

**Key feature:** QF gives you whole *categories* of object beyond unitary gates: open systems, generalized measurements, continuous-variable optics, phase space, and state estimation. Two quiet generalities run underneath the whole model: every object carries its basis, and every object lives in any finite dimension. The point is *scope*, the range of quantum theory QF represents as objects you compute with directly.

### An open-system object: a quantum channel
A noise process is its own object, a `QuantumChannel`, not a gate. Take a generic pure qubit on the Bloch sphere, polar angle $\theta$ and azimuth $\phi$, and confirm it is pure through the purity $\operatorname{Tr}[\rho^2]$, which equals $1$ exactly for a pure state and drops below $1$ for a mixed one:
```wolfram
psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}];
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
A channel is one snapshot of decoherence; the full dynamics is the Lindblad master equation, $\dot\rho = -i[H,\rho] + \sum_i \gamma_i \big(L_i\rho L_i^\dagger - \tfrac12\{L_i^\dagger L_i,\rho\}\big)$, and `QuantumEvolve` solves it *symbolically*. Hand it a density matrix whose four entries are symbols and rates that stay symbols, and one solve answers for every qubit experiment at once. Set the physical region, the decay operator $\sigma_- = |0\rangle\langle 1|$ (written as its matrix so it empties the excited state), and the all-states-at-once $\rho_0$:
```wolfram
$Assumptions = \[Gamma]1 > 0 && \[Gamma]\[Phi] >= 0 && \[CapitalOmega] \[Element] Reals && t >= 0;
\[Sigma]m = QuantumOperator[{{0, 1}, {0, 0}}];
\[Rho]0 = QuantumState[Array[\[Rho], {2, 2}]];
```
Solve the master equation for a qubit with level splitting $\Omega$, energy relaxation at $\gamma_1$ (jump operator $\sqrt{\gamma_1}\,\sigma_-$), and pure dephasing at $\gamma_\phi$ (jump operator $\sqrt{\gamma_\phi}\,Z$):
```wolfram
\[Rho]t = QuantumEvolve[\[CapitalOmega]/2 QuantumOperator["Z"],
   {Sqrt[\[Gamma]1] \[Sigma]m, Sqrt[\[Gamma]\[Phi]] QuantumOperator["Z"]}, \[Rho]0, t];
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
Real detectors are not always projective. A POVM is a set of positive effects $\{E_i\}$ summing to the identity, $\sum_i E_i = I$, with outcome probabilities given by the Born rule $p_i = \operatorname{Tr}[\rho\,E_i]$; build one from your own effects with `QuantumMeasurementOperator[{E1, E2, ...}]`. A *symmetric informationally complete* POVM is the maximally symmetric case: $d^2$ rank-one effects $E_i = \tfrac1d|\psi_i\rangle\langle\psi_i|$ with equal pairwise overlaps $|\langle\psi_i|\psi_j\rangle|^2 = \tfrac{1}{d+1}$, four outcomes on a two-level system, pointing along the vertices of a tetrahedron in the Bloch ball. Measure the generic Bloch qubit $\psi$ from the channel example, the Born rule applied symbolically to all four effects:
```wolfram
sic = QuantumMeasurementOperator["TetrahedronSICPOVM"];
probs = FullSimplify[ComplexExpand[Tr[psi["DensityMatrix"] . Normal[#]] & /@ sic["POVMElements"]]]
```
The four probabilities return exactly as $p_i = \tfrac14\big(1 + \hat r \cdot \hat n_i\big)$, the first reading $\tfrac{1+\cos\theta}{4}$ and the others $\tfrac{1}{12}\big(3 - \cos\theta + \sqrt2\,\sin\theta\,(\ldots)\big)$ along the remaining tetrahedron directions $\hat n_i$, with $\hat r$ the state's Bloch vector; they sum to $1$ identically (verified). Pin the state to the north pole and the general formulas collapse to numbers:
```wolfram
probs /. \[Theta] -> 0
```
The distribution collapses to simple rationals: $|0\rangle$ puts a half on the tetrahedron vertex nearest the pole and a sixth on each of the other three. And because the SIC is *informationally complete*, these four probabilities are not merely statistics of the state, they **are** the state: $\theta$ and $\phi$ can be solved back from them, which is what makes SIC measurements the natural alphabet of tomography.

### Continuous-variable quantum optics
Bosonic Fock space, not qubits. The state of an ideal laser mode is the coherent state $|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\tfrac{\alpha^n}{\sqrt{n!}}\,|n\rangle$, defined as the eigenstate of the photon annihilation operator, $\hat a\,|\alpha\rangle = \alpha\,|\alpha\rangle$, with a *complex* eigenvalue $\alpha$ because $\hat a$ is not Hermitian. QF builds these objects directly in a truncated Fock space (sixteen levels by default). Check the eigenvalue equation with $\alpha$ fully symbolic:
```wolfram
aOp = AnnihilationOperator[];
alphaKet = CoherentState["Normalized" -> False];
Simplify[Most[(aOp[alphaKet] - \[FormalAlpha] alphaKet)["AmplitudesList"]]]
```
Every amplitude of $\hat a|\alpha\rangle - \alpha|\alpha\rangle$ returns identically zero in the symbolic $\alpha$ (the truncation-edge amplitude excluded): the defining equation, verified for the whole family of coherent states at once.

The property that earns these states the title "most classical states of light" is dynamical: under the oscillator Hamiltonian $H = \omega\,(\hat a^\dagger\hat a + \tfrac12)$ a coherent state *stays* a coherent state, $e^{-iHt}\,|\alpha\rangle = e^{-i\omega t/2}\,|\alpha\,e^{-i\omega t}\rangle$, its label circling the phase-space origin exactly like a classical oscillator amplitude, with only the vacuum phase out front. Ask `QuantumEvolve` for the whole identity, symbolically:
```wolfram
U = QuantumEvolve[\[Omega] (aOp["Dagger"] @ aOp + 1/2), None, t];
U @ CoherentState[] == E^(-I \[Omega] t/2) CoherentState[][\[FormalAlpha] E^(-I \[Omega] t)]
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

### A phase-space quasi-probability
The Wigner function is the phase-space quasi-probability $W(x,p) = \dfrac{1}{\pi}\displaystyle\int \langle x - y\,|\,\rho\,|\,x + y\rangle\, e^{2ipy}\,dy$. It integrates to $1$ like a probability density but, unlike one, can go **negative**, and that negativity is a signature of nonclassicality with no classical analogue. The sharpest specimen is a Schrodinger cat state, the even superposition $(|\alpha\rangle + |{-\alpha}\rangle)/\mathcal{N}$ of two coherent states of opposite phase, here with $\alpha = 2$. Build its Wigner function on a phase-space window and read it at two telling points, the origin and half a fringe up in momentum:
```wolfram
w = WignerRepresentation[CatState[][2, 0], {-6, 6}, {-2, 2}, "GridSize" -> 90];
{w[0, 0], w[0, 1/2]}
```
The two values return as $+0.318$ and $-0.236$. The first is the identification that matters: the midpoint between the two coherent components, where a classical mixture of them would put essentially nothing, instead carries the strongest feature of the entire function, an interference peak at $\approx 1/\pi$, the ceiling no Wigner function can exceed; and the neighboring fringe is firmly negative. The whole surface plots directly:
```wolfram
Plot3D[w[x, p], {x, -6, 6}, {p, -2, 2}, PlotRange -> All,
  AxesLabel -> {"x", "p", "W"}, ColorFunction -> "TemperatureMap",
  PlotPoints -> 50]
```
The two smooth Gaussian lobes at $x = \pm 2\sqrt{2}$ are the coherent components, each topping out near $1/2\pi$ as a half-weight classical blob should; the ridge between them, oscillating above and below zero along $p$, is the coherence of the superposition itself. Decoherence erases exactly those fringes and leaves the lobes, which is why this surface is the standard portrait of a quantum superposition of macroscopically distinct states.

### Every object carries its basis
The scope extends inward too: every state, operator, and measurement above is stored *together with the basis it is written in*. A basis is a stored matrix $B$ whose rows are the basis vectors in computational coordinates, a change of basis is the similarity $A_{\mathcal B} = B^{-1} A B$, and objects written in different bases reconcile automatically through the computational frame. About three dozen named bases ship (Bell, Fourier, Schwinger, Dirac, Gell-Mann, Wigner, MUB, SIC), and the same engine lowers non-computational measurements for OpenQASM export and supplies the magic (Bell) basis for two-qubit KAK gate compilation. Watch the representation move while the physics stays put: take the qubit $|+\rangle$, the spread-out superposition $(|0\rangle+|1\rangle)/\sqrt2$ in the computational basis, and write it in the $X$ basis instead:
```wolfram
QuantumState[QuantumState["+"], "X"] // TraditionalForm
```
The state renders as the single ket $|x_+\rangle$: in its own eigenbasis $|+\rangle$ is one basis vector, not a superposition. Same state, different coordinates.

The measurable content is untouched by any such rewriting, and nothing here is qubit-bound. Take a spin-1 particle in its $m_x = +1$ eigenstate, defined directly in the $J_x$ basis:
```wolfram
mx = QuantumState[{0, 0, 1}, "JX"[1]];
mx // TraditionalForm
```
The state renders as the single ket $|J_1^x\rangle$. Re-express it in the computational ($J_z$) basis:
```wolfram
mxZ = QuantumState[mx, "Computational"[3]];
mxZ // TraditionalForm
```
The same state spreads into a three-term superposition over the $J_z$ levels: one ket in its own basis, a sum in the other. The named operator $J_x$ is itself stored in its *own* eigenbasis, so the expectation $\langle J_x\rangle = \langle\psi|J_x|\psi\rangle$ mixes three declared bases at once; compute it from both representations of the state:
```wolfram
<|"JX basis" -> mx["Dagger"][QuantumOperator["JX"[1]][mx]]["Scalar"],
  "Computational basis" -> mxZ["Dagger"][QuantumOperator["JX"[1]][mxZ]]["Scalar"]|>
```
Both entries return as $\langle J_x\rangle = 1$: the state in either coordinate system, the operator in a third, and QF reconciled them all; the invariant survives every change of coordinates.

### Any finite dimension: a spin-1 qutrit and a five-site ring
A spin-1 particle, an NV-center defect for instance, is a genuine three-level system. Build its zero-field-splitting Hamiltonian $H = 2d\,J_z^2 + e\,(J_x^2 - J_y^2)$ from QF's native spin-1 operators, with both couplings left symbolic:
```wolfram
{Jx, Jy, Jz} = QuantumOperator /@ {"JX"[1], "JY"[1], "JZ"[1]};
H = 2 \[FormalD] Jz^2 + \[FormalE] (Jx^2 - Jy^2);
H["Matrix"] // MatrixForm
```
The operator returns fully symbolic, and its matrix displays in the basis the algebra carried along, the $J_x$ eigenbasis: the couplings stay symbols, and the coordinates stay attached to the object.

A Hamiltonian is simplest in its energy eigenbasis, and the operator hands over its own eigenvectors. Build a basis from them and read $H$ there:
```wolfram
QuantumOperator[H, QuantumBasis[H["Eigenvectors"]]]["Matrix"] // Normal // FullSimplify
```
The matrix returns exactly as $\operatorname{diag}\big(0,\ 2d - e,\ 2d + e\big)$: the zero-field-split spectrum in closed form, symbols intact. Building a basis from an object's own eigenvectors works in any dimension, and the spin-$j$ operators ship natively for every $j$.

Dimension need not be a power of two either. A particle hopping on an $N$-site ring lives in an $N$-dimensional Hilbert space, here $N = 5$, with the tight-binding Hamiltonian $H = -\sum_{j=1}^{N}\big(|j\rangle\langle j{+}1| + |j{+}1\rangle\langle j|\big)$, site $N{+}1$ identified with site $1$. Build it, $-1$ for every pair of neighboring sites, and look at it in the site basis:
```wolfram
n = 5;
ring = -Table[
    Boole[Mod[i - j, n] == 1 || Mod[j - i, n] == 1], {i, n}, {j, n}];
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
Tomography is the inverse of measurement: recover an unknown state from the statistics of measurements in several bases, and it is part of the same object model. The three qubit Pauli bases $X, Y, Z$ are mutually unbiased and together informationally complete: their outcome frequencies estimate the three expectations that pin down the density matrix through $\rho = \tfrac12\big(I + \langle X\rangle X + \langle Y\rangle Y + \langle Z\rangle Z\big)$ ([Adamson and Steinberg 2010](https://arxiv.org/abs/0808.0944)). QF runs the whole loop. First, simulate $2000$ measurement shots of an unknown state in each of the three bases:
```wolfram
SeedRandom[42];
true = QuantumState[N @ {Cos[Pi/6], Sin[Pi/6] Exp[I Pi/4]}];
data = QuantumMeasurementSimulation[true, QuantumMeasurementOperator /@ {"X", "Y", "Z"}, 2000];
Values[data]
```
Three pairs of counts return, one per basis: multinomially sampled tallies scattered around the Born-rule proportions, different on every fresh sample, and all a real experiment ever sees.

Now hand the counts to the estimator, which does linear inversion followed by maximum-likelihood, and read the reconstructed density matrix:
```wolfram
est = QuantumStateEstimate[data];
est["MaximumLikelihoodState"]["DensityMatrix"] // Normal // Chop
```
The recovered matrix sits on top of the true one, whose exact entries are $\rho_{11} = \tfrac34$ and $\rho_{12} = \tfrac{\sqrt{3}}{4}\,e^{-i\pi/4}$, with entry-wise residuals at the scale shot noise sets, $\sim 1/\sqrt{2000}$; which entries deviate, and by how much, changes with every fresh sample. To score the reconstruction in one number, compute its fidelity $F(\rho,\sigma) = \operatorname{Tr}\sqrt{\sqrt{\rho}\,\sigma\,\sqrt{\rho}}$ to the state we started from, equal to $1$ exactly when the two states coincide:
```wolfram
QuantumSimilarity[true, est["MaximumLikelihoodState"]]
```
The fidelity returns within a part in a thousand of unity, and typically much closer: the unknown state recovered from finite, noisy synthetic data, using nothing but counts in three bases. Past these unitary measurement bases the same basis engine continues into the over-complete operator frames, where $B^{-1}$ becomes a dual frame and a state turns into a probability vector, the mechanism behind the SIC measurement above.

---

## Claim 3: Measurement records are quantum wires, and feedforward is a controlled gate

**Key feature:** QF implements a projective measurement as von Neumann's measurement *interaction*, not as a collapse: it appends a pointer qudit and entangles it with the system, $V\,|\psi\rangle = \sum_k \big(P_k|\psi\rangle\big)\otimes|k\rangle$. The record of the measurement is therefore itself a quantum wire, numbered $0, -1, -2, \dots$ alongside the system wires, every branch stays in the state, and *feedforward*, acting on the system conditioned on an earlier outcome, is nothing but a gate controlled on a record wire. Protocols that mix measurement with conditioned gates, teleportation, error correction, measurement-based computing, stay inside one coherent symbolic object.

### The record is a wire
Measure all three qubits of a Fourier-transformed GHZ state, then act on qubit $2$ with a CNOT controlled by the *second measurement's record*; the circuit draws what that means:
```wolfram
ff = QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1}, {2}, {3}, "CNOT" -> {-1, 2}}];
ff["Diagram"]
```
The diagram shows the three measurements dropping their records onto new wires labeled $0$, $-1$, $-2$, and the final CNOT reaching from record wire $-1$ to qubit $2$: a gate conditioned on a measurement outcome, drawn and composed like any other gate. Nothing classical ever leaves the object; the record is part of the quantum state, and the circuit keeps computing with it.

### Teleportation without a classical channel
In the textbook protocol Alice measures her two qubits in the Bell basis and phones the two bits to Bob, who applies $X^{m_2} Z^{m_1}$ to his half of an entangled pair. The deferred-measurement principle says the phone call is physically optional: condition the corrections on the records themselves. Build exactly that, qubit $1$ carrying the unknown state, qubits $2$ and $3$ in a Bell pair:
```wolfram
tele = QuantumCircuitOperator[{
    "H" -> 2, "CNOT" -> {2, 3},
    "CNOT" -> {1, 2}, "H" -> 1,
    {1}, {2},
    "CNOT" -> {-1, 3},
    "CZ" -> {0, 3}}];
tele["Diagram"]
```
The diagram is the whole protocol: the Bell pair, the Bell-basis rotation, the two measurements, and then the corrections, an $X$ on Bob's qubit controlled by the second record and a $Z$ controlled by the first, both ordinary controlled gates.

Teleport the generic Bloch qubit $\psi(\theta,\phi)$ from Claim 2. The protocol succeeds if Bob ends up holding exactly the input, $\rho_{\text{Bob}} = |\psi\rangle\langle\psi|$, which for a qubit means the two Bloch vectors coincide; compare them symbolically (the partial trace discards Alice's two qubits and the two records, the four wires before Bob's):
```wolfram
out = tele[QuantumTensorProduct[psi, QuantumState["00"]]];
FullSimplify[ComplexExpand[
   QuantumPartialTrace[out, {1, 2, 3, 4}]["BlochVector"] - psi["BlochVector"]],
  {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}]
```
The difference returns as $\{0, 0, 0\}$ identically in $\theta$ and $\phi$: every pure state teleports exactly, proved in one expression, with no classical channel anywhere in the computation.

### A coherent error, digitized and corrected
The three-qubit repetition code encodes $\alpha|0\rangle + \beta|1\rangle$ as $\alpha|000\rangle + \beta|111\rangle$ and protects it against bit flips. The classic subtlety is that real errors are not discrete flips; here the error will be a *coherent* rotation $R_X(2\theta)$ on the middle qubit, partway between nothing ($\theta = 0$) and a full flip ($\theta = \pi/2$). The founding theorem of error correction says that measuring the syndrome *digitizes* this continuous error into those two discrete cases. Watch it happen symbolically.

First the syndrome measurement. The parity $Z_1 Z_2$ must be measured without learning anything else about the state: a *degenerate* projective measurement with two rank-two projectors $P_\pm = (I \pm Z \otimes Z)/2$, which is precisely what `QuantumMeasurementOperator` accepts as a list of projectors:
```wolfram
id2 = IdentityMatrix[2];
zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
syn[q1_, q2_] := QuantumMeasurementOperator[
   {QuantumOperator[(KroneckerProduct[id2, id2] + zz)/2, {q1, q2}],
    QuantumOperator[(KroneckerProduct[id2, id2] - zz)/2, {q1, q2}]}, {q1, q2}];
```
Then the full code: encode, inflict the coherent error, measure both parities, and decode each of the four syndrome patterns into its correction. The decoder is three feedforward gates with mixed polarities, `"C"[op, ones, zeros]` conditions on some records reading $1$ and others reading $0$:
```wolfram
enc = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}}];
qec = QuantumCircuitOperator[{enc,
    QuantumOperator["RX"[2 \[Theta]], {2}],
    syn[1, 2], syn[2, 3],
    QuantumOperator["C"[QuantumOperator["X", {1}], {0}, {-1}]],
    QuantumOperator["C"[QuantumOperator["X", {2}], {0, -1}, {}]],
    QuantumOperator["C"[QuantumOperator["X", {3}], {-1}, {0}]]}];
```
Run it on a fully symbolic data qubit and read the joint distribution of the two syndrome records, the diagonal of their reduced density matrix (the records sit on the first two of the five wires, the data on the last three):
```wolfram
st = qec[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]];
Assuming[{\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1},
  FullSimplify[ComplexExpand[Diagonal[Normal[QuantumPartialTrace[st, {3, 4, 5}]["DensityMatrix"]]]]]]
```
The four syndrome weights return as $\{\cos^2\theta,\ 0,\ 0,\ \sin^2\theta\}$: the continuous rotation has been *digitized*. The code never sees a partial error; it sees "no error" with probability $\cos^2\theta$ or "qubit $2$ flipped" with probability $\sin^2\theta$, and the weights are independent of $\alpha$ and $\beta$, which is exactly what "measuring the syndrome without learning the state" means.

The feedforward gates then act inside each branch. The theorem to check: the corrected state matches the ideal codeword $\alpha|000\rangle + \beta|111\rangle$ at fidelity $F = 1$, for *every* input and *every* error angle. Compute $F$ symbolically:
```wolfram
Assuming[{\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals, \[Alpha]^2 + \[Beta]^2 == 1},
  FullSimplify[ComplexExpand[QuantumSimilarity[QuantumPartialTrace[st, {1, 2}],
    enc[QuantumTensorProduct[QuantumState[{\[Alpha], \[Beta]}], QuantumState["00"]]], "Fidelity"]]]]
```
The fidelity returns as exactly $1$, for every input $(\alpha, \beta)$ and every error angle $\theta$: not a high number for one sampled case, but the digitization-of-errors theorem itself, handed over by the computation.

### Run the measurement backwards
Because no branch was discarded, the measurement interaction $V$ is an isometry, and an isometry has a left inverse: $V^\dagger V = I$. Measure the Bloch qubit $\psi$, then append the dagger of the very same interaction:
```wolfram
m = QuantumMeasurementOperator[{1}];
FullSimplify[ComplexExpand[
   Normal[QuantumCircuitOperator[{m, m["SuperOperator"]["Dagger"]}][psi]["StateVector"]] -
   Normal[psi["StateVector"]]], {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}]
```
The difference is $\{0, 0\}$: the measurement is undone exactly, phase and all. A projective measurement is irreversible only once its record is discarded; while the record stays a wire, it is reversible physics, exactly as von Neumann wrote it down.

A record, then, is a wire you can condition gates on, trace out, read jointly with the system, or undo. The cost is equally explicit: each measurement appends a pointer qudit, so a circuit with $k$ measurements on $n$ qubits evolves a state of dimension $2^{n+k}$. That growth is not overhead; it is the price of keeping every branch alive, and it is precisely what the teleportation identity, the digitization theorem, and the un-measurement above are made of. Keeping every branch alive is, of course, one computational strategy among several, which raises the next claim: QF has three.

---

## Claim 4: One object model, three computational engines

**Key feature:** the objects of Claims 1 through 3 never specified *how* they are computed. Underneath a circuit application sit three interchangeable engines: exact state-vector algebra (`Method -> "Schrodinger"`), tensor-network contraction (the default), and a compiled stabilizer tableau (`Method -> "Stabilizer"`). The engines cover different regimes, symbols, structure, and scale, and because they share one object model, their representations convert into one another exactly.

### The default engine: a steered tensor-network contraction
Applying a circuit is, by default, a tensor-network contraction, and the *order* of the pairwise contractions decides the cost: every step produces an intermediate tensor, the expense is set by the largest intermediate the order ever creates (in matrix-product language, the bond dimension it forces you to carry), and a poor order can be exponentially more expensive than a good one on the very same network. QF exposes that choice as a `Method` sub-option: `"Path" -> "Greedy"` or `"Optimal"` (with the cost objectives `"flops"`, `"size"`, `"write"`, `"combo"`), or an explicit list of pairwise contractions. Apply a 10-qubit GHZ circuit through a greedily chosen path and read the nonzero amplitudes:
```wolfram
ghz = QuantumCircuitOperator["GHZ"[10]];
ghz10 = ghz[QuantumState["Register"[10]], Method -> {"TensorNetwork", "Path" -> "Greedy"}];
Most @ ArrayRules @ ghz10["StateVector"]
```
The state returns with exactly two nonzero amplitudes, $\{1\} \to \tfrac{1}{\sqrt2}$ and $\{1024\} \to \tfrac{1}{\sqrt2}$, that is $|0\rangle^{\otimes 10}$ and $|1\rangle^{\otimes 10}$: the GHZ state, computed through an explicitly steered contraction.

The plan itself is an inspectable object. Pull the circuit's tensor network, compute the greedy path, and draw the contraction as a tree, each node an intermediate tensor labeled by its dimensions:
```wolfram
net = ghz["TensorNetwork"];
ContractionTree[TensorNetworkContraction[net, GreedyContractionPath[net]], "Labels" -> "Dimensions"]
```
The tree is the engine's execution plan: the leaves are the circuit's gate and state tensors, every internal node is one pairwise contraction, and the dimension labels trace the intermediate sizes the greedy choice keeps small, exactly the quantity a good path minimizes and a bad one lets blow up.

### The stabilizer engine: the state as its symmetry group
For Clifford circuits there is a representation that never touches amplitudes at all: store the state as the group of Pauli operators that fix it, $\{P : P\,|\psi\rangle = +|\psi\rangle\}$, whose $n$ commuting generators determine the $n$-qubit state uniquely and fit in a binary symplectic tableau. This is the Gottesman-Knill representation, in which a Clifford gate is a constant-size rewrite of the tableau rather than a $2^n \times 2^n$ matrix product. The engine is one option away; apply the familiar GHZ circuit with `Method -> "Stabilizer"` and read the generators of the result:
```wolfram
ps = QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"];
ps["PauliForm"]
```
The same circuit that produced an eight-amplitude ket elsewhere in this document now returns a `PauliStabilizer`: three Pauli strings that determine the state completely.

Pauli expectations are closed-form algebra on the generators, with no state vector anywhere:
```wolfram
<|"XXX" -> ps["Expectation", "XXX"], "ZZI" -> ps["Expectation", "ZZI"], "ZII" -> ps["Expectation", "ZII"]|>
```
The two stabilizers return deterministically $+1$, while the unstabilized single-qubit $Z$ is exactly unbiased, all read off by $\mathbb{F}_2$ phase tracking.

The payoff of the representation is scale. Stream twenty thousand random Clifford gates onto a thousand qubits, a state whose amplitude vector would need $2^{1000}$ entries, and time the evolution:
```wolfram
SeedRandom[1]; n = 1000;
gates = Flatten @ Table[{"H" -> RandomInteger[{1, n}], "CNOT" -> RandomSample[Range[n], 2], "S" -> RandomInteger[{1, n}]}, {6666}];
First @ AbsoluteTiming[big = PauliStabilizer[n]["ApplyCircuit", gates];]
```
The whole gate stream folds in a few hundredths of a second ($0.06$ s in this run): `ApplyCircuit` runs the tableau updates through one compiled kernel, so the cost per gate is about three microseconds and independent of the $2^n$ that never gets built. The bipartite entanglement entropy across the half-chain cut is then a closed form, $S(A) = \operatorname{rank}_{\mathbb{F}_2}(\text{generators restricted to } A) - |A|$ ([Fattal et al. 2004](https://arxiv.org/abs/quant-ph/0406168)), exact at any qubit count:
```wolfram
big["Entropy", Range[n/2]]
```
The entropy returns a few ebits shy of the $500$-ebit ceiling: the random circuit has driven the chain to near-maximal, Page-like entanglement, and whatever gate sequence the sampler draws, the value is an exact integer from an $\mathbb{F}_2$ rank, not an estimate. The formalism's boundary is first-class too: a non-Clifford gate ($T$, or a general phase) steps outside the stabilizer group, and QF answers with a `StabilizerFrame`, a superposition of stabilizer states, rather than a wrong answer.

### The representations convert exactly
The engines are not silos; one object model means the representations interconvert. Expand the three-qubit tableau back into an exact state, that is, solve $S_i|\psi\rangle = |\psi\rangle$ for the unique joint $+1$ eigenstate of the three generators:
```wolfram
QuantumState[ps] // TraditionalForm
```
The ket returns as $\tfrac{1}{\sqrt2}\big(|000\rangle + |111\rangle\big)$: the tableau expanded into exact amplitudes. The agreement with the default engine is not approximate; check it with `===`, the language's structural identity:
```wolfram
QuantumState[ps]["StateVector"] === QuantumCircuitOperator["GHZ"[3]][]["StateVector"]
```
`True`: the two engines produced the *same exact expression*, global phase included.

The conversion also runs back to circuits: the tableau synthesizes a Clifford circuit that prepares it.
```wolfram
ps["Circuit"]["Diagram"]
```
A four-gate Clifford circuit returns, and applying it to $|000\rangle$ reproduces the state vector exactly (verified). Circuit to tableau, tableau to state, tableau back to circuit: three representations of one object, each available where it is cheapest, with symbols, structure, or scale deciding which engine runs.

---

## Claim 5: A quantum object is an ordinary Wolfram Language expression

**Key feature:** a QF object is a symbolic expression like any other in the language, with no separate runtime or data format behind it. That is the whole trick: the same `obj["..."]` dispatch and the same built-in functions (`Plot`, `Reduce`, `NProbability`, ...) apply to it directly, with nothing to convert or export first. Each example below is a different WL operation flowing through a QF object.

### Property dispatch: one object, many questions
Ask one object many questions, get answers, no extra packages: here the purity $\operatorname{Tr}[\rho^2]$, the von Neumann entropy $S = -\operatorname{Tr}[\rho\log_2\rho]$, whether the state is entangled, and the concurrence $C$ (a $0$-to-$1$ entanglement measure, maximal for a Bell state).
```wolfram
bell = QuantumState["Bell"];
<|"Purity" -> bell["Purity"], "Entropy" -> bell["VonNeumannEntropy"],
  "EntangledQ" -> QuantumEntangledQ[bell], "Concurrence" -> QuantumEntanglementMonotone[bell, "Concurrence"]|>
```
The answers return self-labeled: the state is pure, carries zero global entropy, is entangled, and has maximal concurrence, four facts off one object.

The same dispatch returns **closed forms** when the state carries a symbol. A Werner state is the singlet diluted with noise, $\rho_W(p) = (1-p)\,|\Psi^-\rangle\langle\Psi^-| + \tfrac{p}{3}\big(I - |\Psi^-\rangle\langle\Psi^-|\big)$, the standard model of an entangled pair degraded by a noisy channel. Ask for its purity and entropy as exact functions of the mixing $p$:
```wolfram
w = QuantumState["Werner"[p, 2]];
Simplify /@ <|"Purity" -> w["Purity"], "Entropy" -> w["VonNeumannEntropy"]|>
```
Both come back as closed forms in the mixing parameter: the purity is $\left|\,1 - 2p + \tfrac{4p^2}{3}\,\right|$ and the entropy is $S(p) = \big[(p-1)\log(1-p) - p\log\tfrac{p}{3}\big]/\log 2$ bits.

### Dispatch returns WL data: a Hamiltonian as an `Association` of Pauli strings
A property can hand back a structured Wolfram Language container. A variational quantum eigensolver writes a Hamiltonian as a sum of Pauli strings, $H = \sum_\alpha h_\alpha P_\alpha$ with $P_\alpha \in \{I, X, Y, Z\}^{\otimes n}$, and measures each string in its own eigenbasis; commuting strings share a measurement setting ([Gokhale et al. 2020](https://arxiv.org/abs/2012.04001)). Take the hydrogen-molecule Hamiltonian in its two-qubit form as a plain numeric matrix and ask for that decomposition:
```wolfram
H2 = QuantumOperator[{{0., 0, 0, 0}, {0, -0.2738, 0.182, 0},
                      {0, 0.182, -1.8302, 0}, {0, 0, 0, 0.1824}}];
H2["PauliDecompose"]
```
The dense matrix comes back as the `Association` $-0.4804\,II + 0.3435\,ZI - 0.4347\,IZ + 0.5716\,ZZ + 0.091\,XX + 0.091\,YY$: the six observables a VQE estimates, in an ordinary key-value container, so the measurement grouping itself (the $Z$ strings share one setting, $XX/YY$ another) is a `Keys`/`Select` computation away.

### Straight into `Plot`
A QF call goes directly inside a built-in. The concurrence of $\cos\theta\,|00\rangle+\sin\theta\,|11\rangle$ is computed by QF and drawn by the language's own `Plot`.
```wolfram
conc[t_?NumericQ] := QuantumEntanglementMonotone[
  QuantumState[{Cos[t], 0, 0, Sin[t]}, 2], {1}, "Concurrence"];
Plot[conc[t], {t, 0, Pi/2}, Frame -> True,
 FrameLabel -> {"\[Theta]", "Concurrence"}, PlotTheme -> "Detailed"]
```
QF also returns the closed form: concurrence $= \sin(2\theta)$.

### Straight into `Reduce`
The physical question: for how much noise does the Werner state's entanglement survive? Get the realignment entanglement criterion, a quantity that is positive exactly when it certifies entanglement, as a symbolic function of $p$:
```wolfram
crit = Simplify[QuantumEntanglementMonotone[QuantumState["Werner"[p, 2]], "Realignment"], 0 < p < 1]
```
Feed that QF expression straight into `Reduce` to solve for the entanglement threshold:
```wolfram
Reduce[crit > 0 && 0 < p < 1, p]
```
The criterion comes back as $\tfrac{|3-4p|-1}{2}$, and `Reduce` turns it into the exact threshold: the Werner state is entangled precisely for $0 < p < \tfrac{1}{2}$.

---

## Claim 6: Every object draws itself

**Key feature:** visualization is property dispatch, not a separate plotting library. Each object exposes plot properties that render through the Wolfram Language graphics stack and forward its options, so a publication figure is one call on the object it depicts. The catalog is wide: states render as Bloch vectors, amplitude and probability charts, and qudit pie and sector wheels (`"PieChart"`, `"SectorChart"`); measurements as outcome histograms; stabilizer states as their binary symplectic tableau (`"TableauForm"`); circuits as diagrams, tensor-network graphs, and multiway and causal-graph views; adiabatic evolutions as spectrum, gap, and coupling plots (`"EnergySpectrumPlot"`, `"SpectralGapPlot"`, `"AdiabaticCouplingPlot"`, `"AdiabaticPathPlot"` in the `QuantumOptimization` sub-package); the underlying `TensorNetwork` also renders as a hypergraph through the TensorNetworks paclet; and the phase-space objects of Claim 2 plot as Wigner surfaces. Below, one example per object kind.

### A state on the Bloch sphere
A single qubit draws its own Bloch vector:
```wolfram
QuantumState[{Cos[Pi/5], Sin[Pi/5] Exp[I Pi/3]}]["BlochPlot"]
```
The output is a `Graphics3D` of the sphere with the state's arrow at polar angle $2\pi/5$ and azimuth $\pi/3$: the geometry of the state, drawn by the state.

### Amplitudes that keep their phases
A bar chart of probabilities discards the phases; `AmplitudesChart` draws the complex amplitudes $\psi_x = \langle x|\psi\rangle$ themselves, heights for the magnitudes $|\psi_x|$ and colors for the phases $\arg\psi_x$. Chart a GHZ state pushed through a quantum Fourier transform, which scrambles two real amplitudes into a phase-rich pattern:
```wolfram
QuantumCircuitOperator[{"GHZ"[3], "Fourier"[3]}][]["AmplitudesChart", ChartLegends -> Automatic]
```
The eight bars return with magnitudes tracing a cosine envelope, the bar at $|100\rangle$ exactly zero, and the surviving bars tinted by a wheel of distinct phase colors keyed by the legend: magnitude *and* phase in one figure, the full complex state on display.

### A measurement draws its statistics
A `QuantumMeasurement` object charts its own outcome distribution, the Born probabilities $p(x) = |\langle x|\psi\rangle|^2$:
```wolfram
QuantumCircuitOperator[{"GHZ"[3], {1, 2, 3}}][]["ProbabilitiesPlot"]
```
Two bars of height $\tfrac12$ at $000$ and $111$: the GHZ correlations as a histogram. The hardware result in Claim 7 is the same kind of object, which is why simulated and measured statistics compare so directly there.

### A stabilizer state shows its tableau
The stabilizer states of Claim 4 have a native picture of their own: the binary symplectic tableau, the actual data structure the formalism computes on. Draw it for the GHZ state:
```wolfram
PauliStabilizer["GHZ"]["TableauForm"]
```
The display is the sign column beside a binary grid, dividers separating the $x$ bits from the $z$ bits and the destabilizer rows from the stabilizer rows; the bottom three rows are exactly the generators $XXX$, $ZZI$, $IZZ$ in symplectic coordinates. This view *is* the algebra: a Clifford gate acts on the state by rewriting a few columns of precisely this grid.

### The circuit, and the network behind it
A circuit draws itself as a gate diagram:
```wolfram
qcv = QuantumCircuitOperator[{"GHZ"[3], "Fourier"[3]}];
qcv["Diagram"]
```
The familiar wire-and-gate picture, Hadamard and CNOT chain into the QFT block. The same object also draws its computational anatomy, the tensor network the engine contracts:
```wolfram
qcv["TensorNetworkGraph"]
```
A `Graph` whose vertices are the gate and state tensors and whose edges are the contracted indices: the diagram shows the physics, the network shows the computation, and Claim 4 steered the order in which such a network is contracted.

The same network renders in the dual picture too, as a hypergraph in which the roles swap, every index a vertex and every tensor one hyperedge spanning its indices:
```wolfram
qcv["TensorNetwork"]["Hypergraph"]
```
The result is a `Hypergraph` with thirteen hyperedges, one per tensor, the three $|0\rangle$ states and the ten gates, drawn over the network's indices: the tensor-contraction picture of the same circuit, in which contracting two tensors is the merging of two hyperedges.

---

## Claim 7: An interoperability hub, down to real hardware

**Key feature:** a QF object can *leave* QF for the rest of the quantum ecosystem, and foreign circuits can *enter* it; `QuantumQASM` is the single OpenQASM door, and it runs in pure Wolfram Language with no external account. The door opens wider at each step below: text export, foreign import, device-aware transpilation, and finally a physical quantum processor.

### Export to OpenQASM 3.0 (native, no Python)
`QuantumQASM` turns a circuit into the standard text format read by Qiskit, Cirq, and most hardware stacks.
```wolfram
QuantumQASM[QuantumCircuitOperator["GHZ"[3]]]
```
The output is a complete OpenQASM 3.0 program: a `qubit[3]` register, a Hadamard on `q[0]` written as a `U` rotation, and the two CNOTs as `ctrl @ U` gates. The native emitter writes OpenQASM 3.0, and the importer below reads OpenQASM 2.0 or 3.0 alike. A 2.0 *dump* is one option away through the qiskit path, `QuantumQASM[qc, "Version" -> 2]` (that route needs Python; the 3.0 export and all imports do not). Export needs concrete angles: a circuit with a symbolic parameter returns a `Failure` rather than a QASM string, so resolve parameters to numbers at this boundary while keeping them symbolic everywhere else in QF (Claims 1 and 5).

### Import foreign OpenQASM and compute with it (native, no Python)
The same door runs the other way: hand `QuantumCircuitOperator` an OpenQASM string authored anywhere and get back a QF circuit, parsed in pure WL, `reset`, controlled gates, mid-circuit `measure` and all. Import a four-qubit circuit that runs the optimal quantum strategy for the CHSH game:
```wolfram
qc = QuantumCircuitOperator["OPENQASM 3.0;
qubit[4] q;
bit[4] c;
reset q[0];
reset q[3];
U(1.5707963267948966, 0., 3.141592653589793) q[3];
U(1.5707963267948966, 0., 3.141592653589793) q[0];
swap q[0] q[3];
ctrl(0) @ negctrl(1) @ U(1.5707963267948966, 3.141592653589793, 3.141592653589793) q[0] q[3];
ctrl(1) @ negctrl(0) @ U(1.5707963267948966, 0., 0.) q[0] q[3];
reset q[1];
U(1.5707963267948966, 0., 3.141592653589793) q[1];
reset q[2];
U(1.5707963267948966, 0., 3.141592653589793) q[2];
U(0., 0., 0.) q[0];
ctrl(1) @ negctrl(0) @ U(1.5707963267948966, 0., 3.141592653589793) q[1] q[0];
U(0.7853981633974483, 0., 0.) q[3];
ctrl(0) @ negctrl(1) @ U(1.5707963267948966, 0., 3.141592653589793) q[2] q[3];
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];"];
```
Once imported it is an ordinary QF object. In the CHSH game the players win when their answers satisfy $x \wedge y = a \oplus b$, and no classical strategy wins more than $\tfrac34$ of the time; entanglement beats that bound. Read the circuit's joint measurement distribution, label the four outcomes $(a, x, y, b)$, and ask Wolfram Language for the winning probability $\Pr[\,x \wedge y = a \oplus b\,]$ in one line:
```wolfram
NProbability[BitAnd[x, y] == BitXor[a, b], {a, x, y, b} \[Distributed] qc[]["MultivariateDistribution"]]
```
The answer is $0.8536 = \cos^2(\pi/8) = \tfrac{2+\sqrt2}{4}$, exactly the **Tsirelson bound**, the quantum maximum for CHSH. A circuit authored in another tool became a probability distribution we query with native WL statistics.

### A hardware-aware transpilation target
`QiskitTarget` carries a device spec (qubit count, native gate set, connectivity), validated at construction against qiskit's own `Target` schema through the Python bridge. Build a 5-qubit linear-chain device and read its coupling map back:
```wolfram
target = QiskitTarget[<|"NumQubits" -> 5, "BasisGates" -> {"cx", "rz", "sx", "x"},
   "CouplingMap" -> {{0,1},{1,2},{2,3},{3,4}}|>];
target["CouplingMap"]
```
Hand the target to `QuantumQASM` as an option, and the circuit is transpiled to that backend's native gates and connectivity:
```wolfram
QuantumQASM[QuantumCircuitOperator["GHZ"[3]], "Target" -> target]
```
The GHZ circuit comes back rewritten for the device: the Hadamard compiled into `rz`/`sx` rotations, the CNOTs as `cx` on physical qubits `$0, $1, $2`, every move allowed by the coupling map.

### Submit to a real IBM QPU, and compare to the exact result
With an IBM account connected through `ServiceConnect["IBMQuantumPlatform"]`, `IBMJobSubmit` sends a circuit to the hardware asynchronously and returns an `IBMJob` handle. A QPU job needs a circuit that *ends in a measurement*; in QF a list of qubit indices placed in the circuit is exactly that, so `{1, 2, 3}` measures all three qubits. Build the measured GHZ circuit, submit it, and check the status:
```wolfram
qc = QuantumCircuitOperator[{"GHZ"[3], {1, 2, 3}}];
job = IBMJobSubmit[qc, "ibm_fez", "Shots" -> 4096];
job["Status"]
```
The handle returns immediately with status `"Queued"` while the job waits in line on the device. The same circuit runs exactly in the Wolfram Language by just applying it:
```wolfram
wl = qc[]
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
The exact bars are $\tfrac12$ at $000$ and $111$ and zero elsewhere; the device's two tall bars land visibly short of $\tfrac12$, a few percent each in this run, with the missing weight spread thinly over the six bitstrings the exact state forbids. That gap, read bar by bar off one figure, *is* the device noise.

---

## Where This Leaves Us

Seven claims, each cashed out in executed code. Symbolic algebra carried dynamics to closed form: the Rabi propagator of a rotating drive, a Landau-Zener sweep expressed in parabolic cylinder functions, a control Hamiltonian recovered exactly from its target evolution, a QAOA landscape whose maximum provably equals the true cut, the entanglement entropy of a quench as a function of time. One object model held the rest of the theory: channels, POVMs, Fock space, and phase space; a master equation solved once for every initial state and every rate, with $T_2 \le 2T_1$ arriving as a theorem rather than a remembered rule; a basis and a dimension on every object; a state reconstructed from measured counts at a fidelity a hair below unity. Measurement stayed inside the physics: records became quantum wires, teleportation ran without a classical channel, a coherent error was digitized and corrected at codeword fidelity exactly $1$, and a measurement ran backwards. Underneath, three engines computed the same objects: exact state-vector algebra, a steered tensor-network contraction, and a stabilizer tableau that carried a thousand qubits to near-maximal half-chain entanglement computed exactly, the representations converting into one another exactly. Because these objects are ordinary Wolfram Language expressions, the language computed with them directly, and each of them drew itself. And the same objects left the system: out to OpenQASM and back, through a transpiler to a device's native gates, and onto a physical quantum processor, whose $4096$-shot histogram stands in this document next to the exact distribution it approximates.

Behind all seven claims is one mechanism: quantum objects that are exact symbolic expressions compose, with each other and with the entire language. The distance from writing a Hamiltonian to reading its closed form, or to running it on hardware, is a handful of composable calls. A result that is a formula can be differentiated, reduced, or pushed to a limit; a result that is a circuit is already in the format the rest of the quantum ecosystem speaks.

The next steps are the reader's: change the graph under the QAOA section, raise the spin past $1$, swap damping for dephasing, tilt the repetition code's coherent error off the $X$ axis, point the last section at a different backend. The document re-derives itself.

---

**Setup.** Evaluate once before the cells:

```wolfram
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
Needs["Wolfram`TensorNetworks`"];
```

Every cell was executed as shown (Wolfram Language 15.0, `Wolfram/QuantumFramework` 2.0.0) and every stated output is the captured result, including the hardware run on `ibm_fez`; reproducing the transpilation cells needs a Python installation with qiskit, and the hardware cells an IBM Quantum account.
