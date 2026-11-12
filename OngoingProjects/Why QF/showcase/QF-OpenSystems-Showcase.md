# Open Quantum Systems in QuantumFramework: A Verified Showcase

An open quantum system exchanges energy and information with an environment it cannot track, so its state $\rho$ obeys not the Schrödinger equation but a master equation, in Lindblad form

$$
\frac{d\rho}{dt} = -i[H,\rho] + \sum_i \gamma_i \left( L_i \rho L_i^\dagger - \tfrac{1}{2}\{L_i^\dagger L_i, \rho\} \right),
$$

with jump operators $L_i$ and rates $\gamma_i$. Every mainstream tool integrates this equation numerically: you fix numbers for the rates and the initial state, and you get a table of density matrices. QuantumFramework solves the same equation **symbolically** when the inputs are symbolic, and its solution is an ordinary `QuantumState` that the rest of the object algebra (channels, operators, circuits, entanglement monotones, phase-space transforms) and the rest of the Wolfram Language (`Reduce`, `Eigensystem`, `WhenEvent`) consume directly. This document is a series of six worked examples built on exactly that combination, laddered from the recognizable to the consequential:

1. **Characterize a qubit**: one solve answers for *all* initial states and *all* rates, and the famous bound $T_2 \le 2T_1$ falls out as a theorem.
2. **Locate a transition**: the entanglement-sudden-death time of a two-qubit state, exactly, with its own validity condition.
3. **Exploit structure**: correlated decay of two atoms, with superradiant and subradiant rates $\gamma \pm \gamma_{12}$ in closed form and a decoherence-free subspace as an exact statement.
4. **Close the algebra**: the propagator of a Liouvillian is an operator object; it equals the named amplitude-damping channel; it composes with gates into a Ramsey sequence with a closed-form fringe.
5. **Couple to classical control**: event-triggered feedback (`WhenEvent`) *inside* the master equation.
6. **Change representation**: the same master equation as a real flow on discrete-Wigner quasi-probabilities, and the exact time at which a magic state loses its negativity.

Each example ends with a one-line *Elsewhere* note locating what the same question costs in a numeric-first tool, with QuTiP (the de-facto standard for open systems) as the reference point; that comparison is the explicit purpose of this document. Every cell was executed against the repo working tree (paclet 2.0.0, Wolfram Language 15.0) and is shown with its verified result. The companion design log is [QF-OpenSystems-Pipeline.md](QF-OpenSystems-Pipeline.md).

Setup:

```wolfram
Needs["Wolfram`QuantumFramework`"]
```

One convention used throughout: the decay operator is written explicitly as

```wolfram
σm = QuantumOperator[{{0, 1}, {0, 0}}];
```

so that $\sigma_- = |0\rangle\langle 1|$ empties the excited state $|1\rangle$, the textbook direction. (QF's named `"J-"` lowers *spin*, and its spin-up state is $|0\rangle$, so it transfers population the other way; defining the matrix once avoids the convention trap.)

---

## Example 1: Qubit relaxation, solved once for every initial state and every rate

**The question.** A qubit with level splitting $\Omega$ decays at rate $\gamma_1$ and dephases at rate $\gamma_\phi$. What does *any* initial state do?

**The QF move.** Feed the master equation a fully symbolic density matrix `Array[ρ, {2, 2}]` and fully symbolic rates. One `QuantumEvolve` call returns the general solution; every concrete experiment is a substitution into it.

Declare the physical region once, so simplifications know the rates are rates:

```wolfram
$Assumptions = γ1 > 0 && γφ >= 0 && Ω ∈ Reals && t >= 0;
```

A density matrix whose four entries are symbols stands for every qubit state at once:

```wolfram
ρ0 = QuantumState[Array[ρ, {2, 2}]]
```

The two decoherence mechanisms are jump operators with symbolic rates, $\sqrt{\gamma_1}\,\sigma_-$ for energy relaxation and $\sqrt{\gamma_\phi}\,Z$ for pure dephasing:

```wolfram
Ls = {Sqrt[γ1] σm, Sqrt[γφ] QuantumOperator["Z"]};
```

Solve the Lindblad equation once, with the Hamiltonian $\tfrac{\Omega}{2}Z$ and time left symbolic:

```wolfram
ρt = QuantumEvolve[Ω/2 QuantumOperator["Z"], Ls, ρ0, t];
```

Its density matrix is the general solution:

```wolfram
FullSimplify /@ Normal[ρt["Computational"]["DensityMatrix"]] // MatrixForm
```

The matrix returns as the complete relaxation law for an arbitrary qubit:

$$
\rho(t) = \begin{pmatrix}
\rho_{11} + \left(1 - e^{-\gamma_1 t}\right)\rho_{22} & \rho_{12}\, e^{-\left(\frac{\gamma_1}{2} + 2\gamma_\phi + i\Omega\right) t} \\[1mm]
\rho_{21}\, e^{-\left(\frac{\gamma_1}{2} + 2\gamma_\phi - i\Omega\right) t} & \rho_{22}\, e^{-\gamma_1 t}
\end{pmatrix}
$$

The excited population empties at $\gamma_1$, the lost weight lands on the ground state, and the coherence precesses at $\Omega$ while decaying at $\frac{\gamma_1}{2} + 2\gamma_\phi$, every experiment on every qubit, in one matrix.

Two physical time scales sit in the exponents: the population relaxation rate $\Gamma_1 = \gamma_1$ and the coherence decay rate $\Gamma_2 = \tfrac{\gamma_1}{2} + 2\gamma_\phi$, i.e. $T_1 = 1/\Gamma_1$ and $T_2 = 1/\Gamma_2$. The textbook bound between them is now a one-line theorem about this solution:

```wolfram
Reduce[ForAll[{γ1, γφ}, γ1 > 0 && γφ >= 0, 1/(γ1/2 + 2 γφ) <= 2/γ1]]
```

The reduction returns `True`: $T_2 \le 2T_1$ holds for every admissible rate pair, a theorem about the solution rather than a remembered rule. Equality holds exactly on the pure-relaxation line:

```wolfram
Reduce[1/(γ1/2 + 2 γφ) == 2/γ1 && γ1 > 0 && γφ >= 0, γφ]
```

This returns exactly as $\gamma_1 > 0 \,\wedge\, \gamma_\phi = 0$: the famous factor of two is saturated precisely when dephasing vanishes.

**Elsewhere.** A numeric solver answers this question one $(\rho_0, \gamma_1, \gamma_\phi)$ triple at a time; in QuTiP the bound $T_2 \le 2T_1$ is something you quote from a textbook, not something the computation hands you.

---

## Example 2: Entanglement sudden death, exactly

**The question.** Entanglement under local noise can vanish at a *finite* time, abruptly, while every local observable still decays smoothly (Yu and Eberly's "entanglement sudden death"). At what time, exactly, and for which states?

**The QF move.** Solve the two-qubit master equation with the mixing angle $\theta$ and the decay rate $\gamma$ symbolic; the concurrence of the solution is a closed form, and `Reduce` extracts the death time *together with its own validity condition*.

The initial state interpolates between product and maximally entangled:

```wolfram
$Assumptions = γ > 0 && t >= 0 && 0 < θ < π/2;
ψ0 = QuantumState[{Cos[θ], 0, 0, Sin[θ]}]
```

Each qubit decays locally; the jump operators are $\sigma_-$ on either factor:

```wolfram
{L1, L2} = {QuantumTensorProduct[σm, QuantumOperator["I"]], QuantumTensorProduct[QuantumOperator["I"], σm]};
```

Solve the master equation in the interaction frame (the free Hamiltonian only rotates phases that drop out of the concurrence):

```wolfram
ρt = QuantumEvolve[0 QuantumOperator["ZZ"], {L1, L2} -> {γ, γ}, ψ0, t];
```

The solution is an X-state, the four corners and the diagonal:

```wolfram
dm = FullSimplify /@ Normal[ρt["Computational"]["DensityMatrix"]]; MatrixForm[dm]
```

The matrix returns with the entangling corner and the single-excitation populations it competes against, writing $E \equiv e^{-\gamma t}$:

$$
\rho(t) = \begin{pmatrix}
\cos^2\theta + (1-E)^2\sin^2\theta & 0 & 0 & E\,\cos\theta\sin\theta \\
0 & E(1-E)\sin^2\theta & 0 & 0 \\
0 & 0 & E(1-E)\sin^2\theta & 0 \\
E\,\cos\theta\sin\theta & 0 & 0 & E^2\sin^2\theta
\end{pmatrix}
$$

For an X-state the concurrence is $C = 2\max\left(0,\ |\rho_{14}| - \sqrt{\rho_{22}\rho_{33}}\right)$; the contest between the entangling corner and the populations it must beat comes out in closed form:

```wolfram
conc = FullSimplify[2 (Abs[dm[[1, 4]]] - Sqrt[dm[[2, 2]] dm[[3, 3]]])]
```

The concurrence returns exactly as $C(t) = 2\sin\theta\, e^{-2\gamma t}\left(e^{\gamma t}(\cos\theta - \sin\theta) + \sin\theta\right)$: the prefactor $2\sin\theta\,e^{-2\gamma t}$ is strictly positive, so entanglement lives or dies with the bracket. `Reduce` on the bracket returns the death time together with the condition for sudden death to exist at all:

```wolfram
Reduce[E^(γ t) (Cos[θ] - Sin[θ]) + Sin[θ] == 0 && t > 0, t, Reals]
```

The reduction returns exactly as $t\gamma = \ln\frac{\sin\theta}{\sin\theta - \cos\theta} \,\wedge\, \sin\theta > \cos\theta$: the death time $t^* = \frac{1}{\gamma}\ln\frac{1}{1 - \cot\theta}$, valid precisely when $\theta > \pi/4$, both the formula and its domain delivered by the same computation.

For $\theta \le \pi/4$ (more weight on $|00\rangle$ than on the doubly excited $|11\rangle$) entanglement only decays asymptotically; past $\pi/4$ it dies at the finite $t^*$ above. The closed form agrees with QF's own numeric concurrence:

```wolfram
{QuantumEntanglementMonotone[QuantumState[dm /. {θ -> Pi/3, γ -> 1, t -> 0.3}], {{1}, {2}}, "Concurrence"], Max[0, conc /. {θ -> Pi/3, γ -> 1, t -> 0.3}]}
```

Both entries return as $0.353558$, agreeing to $10^{-15}$: the object model's numeric entanglement monotone and the hand-read X-state formula are the same number.

**Elsewhere.** In QuTiP you would scan a $(\theta, t)$ grid of `mesolve` runs, compute `concurrence` at each grid point, and read the death time off a contour plot; the boundary $\theta > \pi/4$ would be an empirical observation, not a derived condition.

---

## Example 3: Correlated decay, superradiance, and a decoherence-free subspace

**The question.** Two atoms radiating into a *shared* environment do not decay independently. The general master equation for that situation is Kossakowski's form,

$$
\frac{d\rho}{dt} = -i[H,\rho] + \sum_{ij} \beta_{ij}\left( L_i \rho L_j^\dagger - \tfrac12\{L_j^\dagger L_i, \rho\}\right),
$$

where the positive-semidefinite matrix $\beta$ carries both the individual rates (diagonal) and the environment-induced correlations (off-diagonal). How fast do collective states decay, and can a state hide from the environment entirely?

**The QF move.** The matrix-rate form is a constructor argument: `QuantumOperator["Liouvillian"[H, {L1, L2}, β]]` accepts the full $\beta$ matrix (shown numerically below; for symbolic rates QF works in the eigenbasis of $\beta$, whose eigenvectors are the physical collective modes). Because operators are algebra, the collective modes are one line, and the sub- and superradiant decay rates come out in closed form.

The local decay operators and the correlated rate matrix, with $r = \gamma_{12}/\gamma \in [0, 1)$ and time measured in units of $1/\gamma$:

```wolfram
{L1, L2} = {QuantumTensorProduct[σm, QuantumOperator["I"]], QuantumTensorProduct[QuantumOperator["I"], σm]};
```

The environment couples not to the atoms but to the eigenmodes of $\beta$, and that diagonalization is itself a symbolic one-liner:

```wolfram
Eigensystem[{{1, r}, {r, 1}}]
```

The eigensystem returns exactly as $\{\{1-r,\ 1+r\},\ \{(-1,1),\ (1,1)\}\}$: the antisymmetric combination of the two atoms decays at the *subradiant* rate $1-r$, the symmetric one at the *superradiant* rate $1+r$. Because QF operators form an algebra, the collective jump operators are the same one line:

```wolfram
{Lp, Lm} = {(L1 + L2)/Sqrt[2], (L1 - L2)/Sqrt[2]};
```

Now run the master equation at a strongly correlated environment, $r = 9/10$, kept *exact*. Evolve the antisymmetric Bell state $|S\rangle = (|01\rangle - |10\rangle)/\sqrt2$ and read off its excitation number $\langle N\rangle = \operatorname{Tr}[\rho(t) N]$:

```wolfram
ρS = QuantumEvolve[0 QuantumOperator["ZZ"], {Lp, Lm} -> {1 + 9/10, 1 - 9/10}, QuantumState[{0, 1, -1, 0}/Sqrt[2]], t];
Simplify[Tr[Normal[ρS["Computational"]["DensityMatrix"]] . DiagonalMatrix[{0, 1, 1, 2}]]]
```

The excitation returns exactly as $e^{-t/10}$: the singlet radiates at precisely the subradiant rate $1 - r = 1/10$, an exact exponential with an exact rational rate, not a fitted decay constant. The symmetric state $|T\rangle = (|01\rangle + |10\rangle)/\sqrt2$ under the identical environment:

```wolfram
ρT = QuantumEvolve[0 QuantumOperator["ZZ"], {Lp, Lm} -> {1 + 9/10, 1 - 9/10}, QuantumState[{0, 1, 1, 0}/Sqrt[2]], t];
Simplify[Tr[Normal[ρT["Computational"]["DensityMatrix"]] . DiagonalMatrix[{0, 1, 1, 2}]]]
```

This one returns exactly as $e^{-19t/10}$, the superradiant rate $1 + r = 19/10$: same two atoms, same environment, decay nineteen times faster, on nothing but a relative sign in the state.

At perfect correlation $r = 1$ the subradiant rate is not small, it is *zero*, and QF can sit exactly on that singular point. Evolving the singlet with rates $\{2, 0\}$:

```wolfram
ρDFS = QuantumEvolve[0 QuantumOperator["ZZ"], {Lp, Lm} -> {2, 0}, QuantumState[{0, 1, -1, 0}/Sqrt[2]], t];
Simplify[Normal[ρDFS["Computational"]["DensityMatrix"]] - Normal[QuantumState[{0, 1, -1, 0}/Sqrt[2]]["DensityMatrix"]]]
```

The difference returns as the exact zero matrix, for symbolic $t$: the singlet is a decoherence-free subspace, frozen identically for all time, not to within solver tolerance.

Finally, the Kossakowski matrix itself is a first-class input: the named `"Liouvillian"` constructor takes $\beta$ whole, and the resulting numeric evolution agrees with the collective-mode solution:

```wolfram
LK = QuantumOperator["Liouvillian"[0 QuantumOperator["ZZ"], {L1, L2}, {{1., 0.9}, {0.9, 1.}}]];
ρnum = QuantumEvolve[I LK, QuantumState[{0, 1, -1, 0}/Sqrt[2]], {t, 0, 3}];
{Re[Tr[Normal[ρnum[3.]["Computational"]["DensityMatrix"]] . DiagonalMatrix[{0, 1, 1, 2}]]], N[E^(-3/10)]}
```

Both give $0.740818$, the numeric matrix-rate path and the exact $e^{-3/10}$ agreeing to solver precision. (One honest boundary, logged in the design doc: the matrix-rate constructor currently requires numeric entries, and the fully-symbolic-in-$r$ master-equation solve exceeds `DSolve`'s patience; the symbolic content above comes from the exact $\beta$ eigendecomposition plus exact-rational runs, which is also how the physics is actually organized.)

**Elsewhere.** QuTiP's `mesolve` takes a flat list of collapse operators; a Kossakowski matrix must be hand-diagonalized into effective jump operators before the solver ever sees it, and the decoherence-free subspace at $\gamma_{12} = \gamma$ can only be approached as a sequence of numeric runs.

---

## Example 4: The master equation, the channel, and the circuit are one algebra

**The question.** The Lindblad generator, its propagator, the CPTP channel it integrates to, and the gates of a circuit are, mathematically, citizens of one algebra. Does the software know that?

**The QF move.** Each of those four things is a QF object, and they compose. The generator is a named operator; `QuantumEvolve[..., None, t]` turns it into the propagator *as an operator object*; the propagator provably equals the named amplitude-damping channel; and the evolved state drops into a circuit between two Hadamard gates to give a closed-form Ramsey fringe.

The Lindblad generator of pure decay, packed as one operator object:

```wolfram
$Assumptions = γ > 0 && γφ >= 0 && δ ∈ Reals && t >= 0 && t1 >= 0 && t2 >= 0;
gen = QuantumOperator["Hamiltonian"[None, {σm}, {γ}]];
```

With `None` in the state slot, `QuantumEvolve` returns the propagator itself, a parametric `QuantumOperator`:

```wolfram
U = QuantumEvolve[gen, None, t];
Normal[U["Matrix"]] // MatrixForm
```

The matrix returns as the amplitude-damping propagator on the vectorized density matrix:

$$
S(t) = \begin{pmatrix}
1 & 0 & 0 & 1 - e^{-\gamma t} \\
0 & e^{-\gamma t/2} & 0 & 0 \\
0 & 0 & e^{-\gamma t/2} & 0 \\
0 & 0 & 0 & e^{-\gamma t}
\end{pmatrix}
$$

with the populations relaxing at $\gamma$ and the coherences at $\gamma/2$, written once for all times.

A propagator of a time-independent generator must form a semigroup, $S(t_1)S(t_2) = S(t_1{+}t_2)$, and this one does, symbolically:

```wolfram
FullSimplify[(Normal[U["Matrix"]] /. t -> t1) . (Normal[U["Matrix"]] /. t -> t2) == (Normal[U["Matrix"]] /. t -> t1 + t2)]
```

The check returns `True` for symbolic $t_1$, $t_2$: the CPTP semigroup law, verified as an identity rather than sampled.

Integrating the master equation for time $t$ should be *the same map* as the textbook amplitude-damping channel at damping parameter $\gamma(t) = 1 - e^{-\gamma t}$. QF has that channel as a named object, so the claim is checkable on a fully symbolic state:

```wolfram
viaME = FullSimplify /@ Normal[QuantumEvolve[0 QuantumOperator["Z"], {Sqrt[γ] σm}, ρ0, t]["Computational"]["DensityMatrix"]];
```

```wolfram
viaChannel = FullSimplify /@ Normal[QuantumChannel["AmplitudeDamping"[1 - E^(-γ t)]][ρ0]["Computational"]["DensityMatrix"]];
```

```wolfram
FullSimplify[viaME == viaChannel]
```

The equality returns `True` on the fully symbolic $\rho_0$ from Example 1: integrating the master equation for time $t$ *is* the Kraus channel at damping $1 - e^{-\gamma t}$, for every state and every rate, as an identity of the algebra.

The continuous-time and discrete-channel pictures of decoherence are literally the same object, with the time dependence of the damping parameter derived rather than assumed.

Finally, the circuit picture. A Ramsey experiment is a Hadamard, a free evolution of duration $t$ under detuning $\delta$ with decay and dephasing, and a second Hadamard; the evolved open-system state is an ordinary `QuantumState`, so the gates apply to it directly:

```wolfram
ψplus = QuantumOperator["H"][QuantumState["0"]];
```

```wolfram
ρwait = QuantumEvolve[δ/2 QuantumOperator["Z"], {Sqrt[γ] σm, Sqrt[γφ] QuantumOperator["Z"]}, ψplus, t];
```

```wolfram
P0 = FullSimplify[ComplexExpand[Normal[QuantumOperator["H"][ρwait]["Computational"]["DensityMatrix"]][[1, 1]]]]
```

The probability returns exactly as $P_0(t) = \tfrac12\left(1 + e^{-(\gamma + 4\gamma_\phi)t/2}\cos\delta t\right)$: the Ramsey fringe at the detuning frequency inside the decoherence envelope $e^{-\Gamma_2 t}$, $\Gamma_2 = \tfrac{\gamma}{2} + 2\gamma_\phi$, the workhorse formula of qubit characterization, derived end-to-end by the object algebra.

**Elsewhere.** QuTiP has numeric propagators (`propagator`) and numeric channel conversions (`to_kraus`, `to_choi`), but no named symbolic channel to equate against and no symbolic state for gates to act on; each arrow of this example exists only pointwise in the rates and time.

---

## Example 5: Event-triggered control inside the master equation

**The question.** Real experiments are hybrid: a classical controller watches the quantum system and acts on it, when a population crosses a threshold, a drive switches off, a pulse fires, a setpoint changes. Can the *solver* express that, rather than the user stitching together piecewise runs?

**The QF move.** `QuantumEvolve` accepts `"AdditionalEquations"`, and its numeric path is `NDSolve`, so the entire `WhenEvent`/`DiscreteVariables` machinery of the language applies *inside* the master equation. The drive amplitude becomes a discrete variable of the ODE system; an event on the quantum state's matrix elements flips it.

The system: a qubit under a resonant Rabi drive of amplitude $2\pi$, decaying at $\gamma = 0.5$. First, *calibration*: an event that fires when the excited population crosses $0.6$ and `Sow`s the crossing time, in a single integration:

```wolfram
{sol, pts} = Reap[QuantumEvolve[π QuantumOperator["X"], {Sqrt[0.5] σm}, QuantumState["0"], {\[FormalT], 0, 2}, "AdditionalEquations" -> WhenEvent[Re[Indexed[\[FormalS][\[FormalT]], {2, 2}]] > 0.6, Sow[\[FormalT]]]]]; pts
```

The crossing times return as $\{\{0.296047,\ 1.31377\}\}$: the damped Rabi oscillation pierces the threshold on its first rise and once more on its second, each instant located by the integrator's event machinery.

The event condition watches a *matrix element of the evolving density matrix* (`\[FormalS]` is the dependent variable of the master equation; `Indexed[..., {2, 2}]` is the excited population), and the solver reports the crossing to event-location precision, with no scanning and no restart.

Second, *control*: make the drive amplitude itself a discrete variable of the ODE system, switched off at exactly the crossing time just found, so the drive fills the qubit to the target population and then releases it. `DiscreteVariables` is an `NDSolve` option, and `QuantumEvolve` forwards it:

```wolfram
H = QuantumOperator[drive[\[FormalT]]/2 {{0, 1}, {1, 0}}, "ParameterSpec" -> \[FormalT]];
```

```wolfram
f = QuantumEvolve[H, {Sqrt[0.5] σm}, QuantumState["0"], {\[FormalT], 0, 3}, "AdditionalEquations" -> {WhenEvent[\[FormalT] > pts[[1, 1]], drive[\[FormalT]] -> 0], drive[0] == 2 π}, DiscreteVariables -> {drive}, "ReturnSolution" -> True];
```

```wolfram
Table[Re[f[τ][[2, 2]]], {τ, {0.2, 0.5, 1.0, 2.0, 3.0}}]
```

The samples return as $\{0.329,\ 0.542,\ 0.422,\ 0.256,\ 0.155\}$: rising toward the threshold before the switch at $t = 0.296$, then relaxing at exactly the bare decay rate ($0.6\,e^{-0.5\,(3 - 0.296)} \approx 0.155$, matching the last sample). Any time slice wraps back into the object algebra as `QuantumState[f[τ]]`. (`"ReturnSolution" -> True` hands back the raw `InterpolatingFunction`: solutions with event discontinuities currently trip QF's symbolic re-packaging step, a boundary logged in the design doc.)

Plot the controlled trajectory (code only, as in the main showcase):

```wolfram
Plot[Re[f[τ][[2, 2]]], {τ, 0, 3}, PlotRange -> {0, 1}, AxesLabel -> {"t", "P₁"}]
```

**Elsewhere.** QuTiP's solvers have no event mechanism; the standard workaround is to integrate to a guessed crossing, detect it in post-processing, then restart `mesolve` with new arguments, manually, once per event.

---

## Example 6: The master equation in discrete phase space

**The question.** A qutrit state can be written as a discrete Wigner quasi-probability vector $W$, nine real numbers on a $3\times3$ phase-space grid, with $\sum_k W_k = 1$. In odd dimensions, negative entries of $W$ are *the* resource for quantum computation: a state with non-negative $W$ is efficiently classically simulable, and Wigner negativity is a magic monotone (Howard et al. 2014). How does decoherence look in this representation, and when exactly does a magic state's negativity die?

**The QF move.** The Wigner basis is an ordinary `QuantumBasis`; a state expressed in it carries `Picture -> "PhaseSpace"`, and `QuantumEvolve` recognizes that picture and writes the master equation as a *real* flow on the quasi-probabilities, a classical-looking rate equation. The negativity death time then comes out exactly.

The "strange state", the most negative qutrit state:

```wolfram
$Assumptions = γ > 0 && t >= 0;
ψS = QuantumState[{0, 1, -1}/Sqrt[2]];
```

Its discrete Wigner representation is one basis change:

```wolfram
wS = QuantumState[ψS["Double"], QuantumBasis["Wigner"[3]]];
Simplify[Normal[wS["StateVector"]]]
```

The quasi-probabilities return exactly as $\{-\tfrac13, \tfrac16, \tfrac16, \tfrac16, \tfrac16, \tfrac16, \tfrac16, \tfrac16, \tfrac16\}$: a true probability-like vector summing to $1$, except for the $-\tfrac13$ sitting at the origin of phase space, the largest negativity any qutrit state can have.

The noise: full depolarizing, built from all eight non-identity Heisenberg-Weyl unitaries $X^a Z^b$ as jump operators at rate $\gamma/9$:

```wolfram
{X3, Z3} = Normal[QuantumOperator[#]["Matrix"]] & /@ {"X"[3], "Z"[3]};
Ls = Sqrt[γ/9] QuantumOperator[#] & /@ Rest[Flatten[Table[MatrixPower[X3, a] . MatrixPower[Z3, b], {a, 0, 2}, {b, 0, 2}], 1]];
```

Ask `QuantumEvolve` for the equations it would solve in the phase-space picture:

```wolfram
eqs = QuantumEvolve[0 QuantumOperator["Z"[3]], Ls, wS, t, "ReturnEquations" -> True];
```

The first element of `eqs` is literally $\dot W = G\,W$ with a *real* $9\times9$ generator $G$ acting on the nine quasi-probabilities; no imaginary unit appears anywhere. Sum each column of $G$:

```wolfram
Total[eqs[[1, 1, 2, 1]]]
```

The column sums return as $\{0,0,0,0,0,0,0,0,0\}$.

Every column of the generator sums to zero: in this representation the Lindblad equation *is* a continuous-time Markov master equation on nine phase-space points, probability conserved, classical in form, quantum only through the negative entry it transports. Solve the master equation and transform the solution back into the Wigner basis:

```wolfram
ρt = QuantumEvolve[0 QuantumOperator["Z"[3]], Ls, ψS, t];
wt = Simplify[Normal[QuantumState[ρt["Double"], QuantumBasis["Wigner"[3]]]["StateVector"]]]
```

The quasi-probabilities return exactly as $W(t) = \left\{\tfrac19 - \tfrac49 e^{-\gamma t},\ \tfrac{2 + e^{-\gamma t}}{18} \times 8\right\}$: the negativity sits in the first component, relaxing toward the uniform distribution $\tfrac19$ of the maximally mixed state. It crosses zero at a finite, exact time:

```wolfram
Reduce[Min[wt] == 0 && t > 0 && γ > 0, t, Reals]
```

The reduction returns exactly as $\gamma > 0 \,\wedge\, t = \frac{2\ln 2}{\gamma}$.

At $t^* = \frac{2\ln 2}{\gamma}$, exactly, the strange state stops being a magic resource: before this instant no classical sampler can reproduce it, after it the Wigner function is a genuine probability distribution. A phase transition in computational power, time-stamped in closed form by the master equation.

**Elsewhere.** QuTiP computes continuous-variable Wigner functions for *display* of a numerically evolved state; a discrete-Wigner representation in which the master equation itself becomes a classical Markov flow, with an exact classicalization time, has no counterpart.

---

## Where This Leaves Us

Six questions about open quantum systems, six answers that are formulas rather than tables. One symbolic solve characterized a qubit for every initial state and every rate pair at once, and handed back $T_2 \le 2T_1$ as a theorem. The death of entanglement under local noise arrived with its own exact time stamp, $t^* = \frac{1}{\gamma}\ln\frac{1}{1-\cot\theta}$, and its own validity condition $\theta > \pi/4$. Two atoms in a shared environment decayed at exactly $1 \pm r$ along the eigenmodes of the Kossakowski matrix, and at $r = 1$ the singlet froze identically, not approximately. The propagator of a Liouvillian proved itself a semigroup, proved itself equal to the named amplitude-damping channel at $\gamma(t) = 1 - e^{-\gamma t}$, and then sat between two Hadamard gates to derive the Ramsey fringe. A classical control law fired inside the master equation. And the same master equation, re-expressed in a discrete phase-space basis, became a real flow of quasi-probabilities whose negativity, the resource that powers quantum computation, expired at exactly $t^* = \frac{2\ln 2}{\gamma}$.

The common mechanism is the one the main showcase documents across all of quantum theory: QF objects are exact symbolic expressions, the master-equation solver is one method on that algebra rather than the centerpiece of the tool, and whatever comes out, a state, a propagator, a channel, is immediately a citizen of both the object model and the Wolfram Language. Where a numeric-first package answers an open-system question with a trajectory, this one answers with the function the trajectory samples. Change the jump operators, the correlation $r$, the control threshold, or the phase-space dimension, and the document re-derives itself.
