# Talk outline (50 minutes): what the Wolfram Language and QuantumFramework let a physicist do

Eight topics, about 25 slides. Each topic opens with a slide that asks one physics question and names the computational idea that answers it. Its result slides each show one short visible input, a picture or a formula, and a sentence saying what was found. The long setup code lives in `TalkSetup.wl`, loaded by one visible `Get` on the title slide, so nothing on a slide is hidden.

Visuals marked **(existing)** are plots already in the verified examples. Visuals marked **(new)** are a few lines of code built from objects the examples already compute: each would be written and checked in a kernel before the deck is built. The diagram, plot and Wigner properties named below (`"Diagram"`, `"ProbabilitiesPlot"`, `"BlochSpherePlot"`, `WignerFunction`, `GraphState`) were confirmed to exist in the current QuantumFramework.

## Opening (2 slides)

Title slide, then an agenda slide listing the eight questions below, each with its computational idea in a few words.

## 1. What does a driven, decaying qubit do, for every rate at once?

Core idea: solve the master equation once with every rate and every initial-state entry left as a symbol.

What is found: the complete relaxation law in closed form, and $T_2\le2T_1$ derived from it by `Reduce`, with equality exactly when pure dephasing vanishes.

Slides and visuals:
- The symbolic density matrix (existing formula).
- Bloch-vector trajectories from the symbolic solution for a few rate ratios, spiralling into the steady state (new: `ParametricPlot3D` of the solution).
- $T_2\le2T_1$ from `Reduce` (existing).

## 2. Does watching a qubit speed it up or freeze it?

Core idea: the Liouvillian of a driven, measured qubit has a closed-form spectrum, so the crossover can be found exactly.

What is found: below the exceptional point $k=\Omega$ measurement speeds relaxation; above it the qubit freezes at the Zeno rate $\Omega^2/2k$. The same Liouvillian gives the Mollow triplet and fluorescence squeezed below the vacuum.

Slides and visuals:
- Eigenvalues moving in the complex plane as $k$ grows, colliding at $k=\Omega$ and splitting onto the real axis (new: `ListPlot` of the symbolic eigenvalues).
- Rabi frequency and slowest rate against $k$ (existing plot).
- Mollow triplet, weak against strong drive (existing plot).

## 3. Can a continuous error be fixed by a discrete correction?

Core idea: measurement records are wires in the circuit, so feedforward correction is part of one symbolic object.

What is found: parity measurements turn an arbitrary rotation into "no error" with probability $\cos^2\theta$ or "one bit flipped" with $\sin^2\theta$, and the correction restores every input at every angle with fidelity exactly 1.

Slides and visuals:
- The full circuit with its record wires and classically controlled corrections (new: `qec["Diagram"]`, already checked to render).
- The two syndrome probabilities against the error angle (new: `Plot` of the symbolic result).
- Fidelity 1 for symbolic data and angle (existing).

## 4. Can watching a particle squeeze it below its zero-point motion?

Core idea: the conditional covariances obey a Riccati equation that `SolveValues` and `Resolve` handle with $\lambda,\omega,m,\hbar$ symbolic; a grid simulation checks it.

What is found: yes, at every measurement strength; a detector of efficiency $\eta$ leaves purity exactly $\sqrt\eta$; and five numbers (a quantum Kalman filter) track what 256 amplitudes compute.

Slides and visuals:
- Squeezing ratio against measurement strength with its weak and strong limits overlaid (new: `Plot` of the closed form).
- The conditional probability density $|\psi(x,t)|^2$ of the watched particle as a space-time density plot, jittering as readings arrive (new: from the stored grid states).
- Kalman filter on top of the full wavefunction's mean (existing plot).

## 5. What does a beam splitter do to $n$ photons, for every $n$?

Core idea: `QuantumInterferometer` keeps photon numbers, angles and squeezing as symbols, and draws itself.

What is found: the amplitude is a Wigner $D$-function (a beam splitter is a spin rotation); odd splits of $|n,n\rangle$ vanish for every $n$; four-photon cancellation needs an unbalanced splitter at $R=(3\pm\sqrt3)/6$; a three-mode network flips the sign of $|2\rangle$ alone (the KLM gate).

Slides and visuals:
- Interferometer diagram and the symbolic amplitude (new diagram, existing formula).
- Output distributions of $|n,n\rangle$ for $n=1,\ldots,4$ as bar charts, the odd outputs empty (new: `BarChart` of the existing formula).
- Return probability of $|2,2\rangle$ against reflectivity, touching zero twice (new: `Plot`).
- KLM gate diagram and its heralded amplitudes (new diagram, existing result).
- Wigner function of a squeezed or cat state (new: `WignerFunction` with `DensityPlot`).

## 6. How far can exact simulation go without a state vector?

Core idea: stabilizer tableaux turn Clifford circuits into bit rewrites; T gates add terms to a stabilizer frame.

What is found: the five-qubit code's fifteen errors have fifteen distinct syndromes; a 1000-qubit GHZ state has one bit of entanglement across any cut; the cost of magic grows with the number of T gates, not qubits.

Slides and visuals:
- Syndrome table as a 15-by-4 array picture (new: `ArrayPlot` of the existing table).
- A graph state drawn as its graph (new: `GraphState`, to check).
- Frame size against T count beside $2^n$ on a log scale (new: plot of the existing numbers).

## 7. How short can a gate on a real superconducting qubit be?

Core idea: a five-level transmon with a DRAG term, evolved and calibrated in one parametrized `QuantumOperator`.

What is found: short pulses leak out of the qubit; DRAG suppresses leakage near its weak-drive value, but completing the flip needs amplitude, detuning and DRAG tuned together; decoherence, not leakage, sets the speed limit.

Slides and visuals:
- Level populations during an 8 ns pulse with and without DRAG (new: `Plot` from `QuantumEvolve`).
- Leakage and transfer error against the DRAG coefficient (existing plot).
- Optional, if TransmonLab is loaded: the readout spirals and single-shot IQ clouds from the Qolab demo (existing `IQPlot`).

## 8. Why is sampling photons hard for a classical computer?

Core idea: amplitudes are permanents; their cost depends on how photons occupy the modes.

What is found: bosons bunch, classical particles split, fermions antibunch, set only by the exchange sign; generic amplitudes cost exponentially more with photon number while bunched ones stay cheap.

Slides and visuals:
- Bosons, classical particles and fermions on one splitter as three bar charts (new: from the existing table).
- Amplitude time against photon number for the three occupation types on a log scale (new: the existing timing call over more values of $n$).

## Closing (1 slide)

The eight questions in one table against what answered them: closed-form solution, exact proof by `Reduce`/`Resolve`, measurement as circuit wires, symbolic photon numbers, stabilizer algebra, a parametrized pulse model, permanents.

Dropped from the 50-minute version, available as backup: continuous QEC filtering, Wong-Zakai with NDSolve, the gravitational soliton, Dirac matrices, Snyder spacetime, Shor's discrete log, the formal semantics.
