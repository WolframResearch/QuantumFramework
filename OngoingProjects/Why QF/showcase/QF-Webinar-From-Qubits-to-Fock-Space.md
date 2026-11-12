---
Template: ComputationalEssay
ResourceType: ComputationalEssay
Name: From Qubits to Fock Space with Wolfram Quantum Framework
Author: Mads Bahrami
Date: 2026
Description: A compact live tour of symbolic quantum computation across circuits, open systems, basis changes, tensor networks, and continuous-variable optics
Abstract: Wolfram Quantum Framework treats states, operators, channels, measurements, and circuits as symbolic objects that compose. This webinar follows one idea across several physical settings: parameters should survive the computation. We begin with a symbolic qubit, derive an exact QAOA landscape, and watch a Heisenberg interaction generate entanglement. The same object model then describes amplitude damping, a generalized measurement, and a spin-1 basis transformation. We finish in bosonic Fock space with coherent and cat states, then expose the tensor network behind a circuit. The goal is not a catalog of functions, but a compact demonstration of one algebra spanning distinct quantum problems.
Keywords: [quantum computing, symbolic computation, Wolfram Language, open quantum systems, Fock space, tensor networks, quantum optics]
Sources: ["[Wolfram Quantum Framework](https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/)", "[QF verified showcase](QF-Showcase.md)"]
Links: ["[What Is a Computational Essay?](https://writings.stephenwolfram.com/2017/11/what-is-a-computational-essay/)"]
---

Quantum software often separates qubit circuits, open systems, spin models, and quantum optics into different computational worlds. Wolfram Quantum Framework takes the opposite view: these are different objects in one algebra.

The decisive feature is not merely that symbols are accepted as input. Symbols survive composition. A parameter introduced in a state, gate, channel, or Hamiltonian can propagate through an entire workflow and return as an exact physical statement.

The examples target Wolfram Quantum Framework 2.1.0 and Wolfram Language 14.3 or later. Participants should install the [paclet](https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/) before the session.

Load the framework and its continuous-variable extension:

```wl
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

## A state can represent a whole family

Begin with the general pure qubit
$|\psi\rangle=\cos\alpha|0\rangle+e^{i\beta}\sin\alpha|1\rangle$.
The state is stored once, with both parameters intact.

Construct the symbolic state:

```wl
psi = QuantumState[{Cos[\[Alpha]], Exp[I \[Beta]] Sin[\[Alpha]]}]
```

Its Bloch vector is an exact map from amplitudes to geometry:

```wl
FullSimplify[ComplexExpand[psi["BlochVector"]]]
```

The result is
$\big(\sin(2\alpha)\cos\beta,\sin(2\alpha)\sin\beta,\cos(2\alpha)\big)$.
One object represents every pure qubit and exposes its physical geometry without numerical sampling.

## A circuit can become an exact landscape

The same symbolic propagation reaches through a variational circuit. For MaxCut on the complete graph $K_4$, one QAOA layer has two angles, $\gamma$ for the cost layer and $\theta$ for the mixer.

Define the graph and derive its dimensions from the graph itself:

```wl
graph = CompleteGraph[4];
n = VertexCount[graph];
edges = List @@@ EdgeList[graph];
```

Build $H_C=\sum_{(i,j)\in E}(1-Z_iZ_j)/2$ directly from those edges:

```wl
hCost = QuantumOperator[Total[With[{id = QuantumOperator[StringRepeat["I", n]]/2},
    (id - QuantumOperator["ZZ"/2, #] &) /@ edges]], "Label" -> "Hamiltonian"]
```

Prepare $|+\rangle^{\otimes n}$:

```wl
prep = QuantumCircuitOperator[{StringRepeat["0", n], "H" -> Range[n]},
    "Label" -> "Prep"]
```

Make one cost rotation for every graph edge:

```wl
cost = QuantumCircuitOperator[("R"[\[Gamma], "ZZ"] -> # &) /@ edges,
    "Parameters" -> {\[Gamma]}, "Label" -> "Cost"]
```

Make the transverse-field mixer:

```wl
mixer = QuantumCircuitOperator[{"RX"[\[Theta]] -> Range[n]},
    "Parameters" -> {\[Theta]}, "Label" -> "Mixer"]
```

Compose the ansatz without reducing any block to a matrix:

```wl
oneLayer = QuantumCircuitOperator[{prep, cost, mixer},
    "Parameters" -> {\[Theta], \[Gamma]}]
```

The expectation value is itself a composable bra-operator-ket circuit:

```wl
meanValue = oneLayer /* hCost /* oneLayer["Dagger"];
meanValue["Diagram"]
```

Evaluating that object returns the entire variational objective as one formula:

```wl
landscape = FullSimplify[meanValue[]["Scalar"],
    Element[{\[Gamma], \[Theta]}, Reals]]
```

The result is
$$
\langle H_C\rangle=\frac{3}{8}\left(7+\cos2\theta+2\cos4\gamma\sin^2\theta-8\cos^2\gamma\sin\gamma\sin2\theta\right).
$$
An experiment estimates points on this surface. Here the surface is the result.

Plot the exact objective before optimizing it:

```wl
Plot3D[landscape, {\[Gamma], 0, 2 Pi}, {\[Theta], 0, 2 Pi},
    PlotRange -> All, AxesLabel -> {"\[Gamma]", "\[Theta]", "cost"},
    ColorFunction -> "TemperatureMap", PlotPoints -> 50]
```

Because the objective remains symbolic, its optimum can be certified algebraically rather than read from a grid:

```wl
Maximize[landscape, {\[Theta], \[Gamma]}]
```

The maximum and the optimizing angles are exact algebraic numbers. This is the useful meaning of symbolic variational computation: not symbolic gates in isolation, but an analytic objective at the end of the circuit.

## Evolution can produce a formula for entanglement

Now replace the circuit by a many-body Hamiltonian. Let two spins interact through
$H=J(XX+YY+ZZ)$, beginning in the product state $|01\rangle$.

Construct the Hamiltonian and evolve the state exactly:

```wl
hHeisenberg = J Total[QuantumOperator /@ {"XX", "YY", "ZZ"}];
psiT = QuantumEvolve[hHeisenberg, QuantumState["01"], t]
```

Trace out the second spin and ask for the von Neumann entropy:

```wl
entropy = FullSimplify[
    QuantumPartialTrace[psiT, {2}]["VonNeumannEntropy"], t J \[Element] Reals]
```

The answer is the binary entropy of the two exchange populations. It is zero at product states and reaches one bit at the Bell points.

Plot the complete entanglement history for $J=1$:

```wl
Plot[entropy /. J -> 1, {t, 0, Pi}, Frame -> True,
    FrameLabel -> {"time", "entanglement entropy [bits]"},
    PlotRange -> {0, 1.05}]
```

The state, its reduced density matrix, and a nonlinear entanglement measure have all remained exact functions of time.

## Noise and measurement use the same object model

Closed evolution is only one part of quantum mechanics. A channel is a first-class object, so a symbolic state can pass directly through a nonunitary process.

Construct a generic Bloch-sphere state:

```wl
qubit = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
```

Compute the purity lost under amplitude damping with probability $\gamma$:

```wl
purityLoss = FullSimplify[
    qubit["Purity"] - QuantumChannel["AmplitudeDamping"[\[Gamma]]][qubit]["Purity"],
    {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals, 0 <= \[Gamma] <= 1}]
```

The exact result is
$2\gamma(1-\gamma)\sin^4(\theta/2)$.
It immediately reveals what numerical cases obscure: the azimuth is irrelevant, the ground state is fixed, and full relaxation is pure again.

Generalized measurements fit the same algebra. A tetrahedral SIC POVM has four outcomes and is informationally complete for a qubit.

Construct the measurement and apply the Born rule symbolically to all four effects:

```wl
sic = QuantumMeasurementOperator["TetrahedronSICPOVM"];
probabilities = FullSimplify[ComplexExpand[
    Tr[qubit["DensityMatrix"] . Normal[#]] & /@ sic["POVMElements"]]]
```

Verify normalization for every $\theta$ and $\phi$ at once:

```wl
FullSimplify[Total[probabilities],
    Element[{\[Theta], \[Phi]}, Reals]]
```

The result is `1`. States, channels, and measurements do not require separate workflows. They compose and expose the same property-based interface.

## Change the basis, keep the physics

A useful basis can turn a difficult description into a transparent one. The framework stores the basis with each object and reconciles different bases during composition.

Define a spin-1 particle in the $m_x=+1$ eigenstate:

```wl
mx = QuantumState[{0, 0, 1}, "JX"[1]]
```

Rewrite the same state in the three-dimensional computational basis:

```wl
mxZ = QuantumState[mx, "Computational"[3]];
{mx, mxZ} // TraditionalForm
```

One representation is a single $J_x$ eigenket. The other is a three-term superposition. Test the invariant in both representations:

```wl
jx = QuantumOperator["JX"[1]];
{mx["Dagger"][jx[mx]]["Scalar"], mxZ["Dagger"][jx[mxZ]]["Scalar"]}
```

Both expectations are `1`. The coordinates changed; the state and its measurable content did not.

## The same algebra reaches Fock space

Qubits are not the boundary. The continuous-variable extension supplies bosonic operators and states in a truncated Fock basis.

Construct the annihilation operator and a symbolic coherent state:

```wl
a = AnnihilationOperator[];
alphaKet = CoherentState["Normalized" -> False];
```

Check the defining equation $a|\alpha\rangle=\alpha|\alpha\rangle$ away from the truncation edge:

```wl
Simplify[Most[(a[alphaKet] - \[FormalAlpha] alphaKet)["AmplitudesList"]]]
```

Every returned amplitude is zero. The code is finite-dimensional, but the symbolic identity survives throughout the represented sector.

Photon correlations then distinguish three physically different states of light:

```wl
g2 = <|"Fock" -> G2Coherence[FockState[1]],
    "Coherent" -> G2Coherence[CoherentState[20][1.5]],
    "Thermal" -> G2Coherence[ThermalState[1., 20]]|>
```

Display antibunched, coherent, and bunched statistics together:

```wl
BarChart[g2, ChartLabels -> Automatic, Frame -> True,
    FrameLabel -> {None, "g2(0)"}, PlotRange -> {0, 2.2}]
```

The values are $0$, numerically $1$, and approximately $2$. The tiny deviations from the ideal coherent and thermal values are the visible cost of Fock-space truncation, not statistical noise.

A cat state makes the phase-space distinction more vivid. Its Wigner function can be negative, while its Gaussian-smoothed Husimi function cannot.

Compute both representations of the same cat state before the live session:

```wl
cat = CatState[][2, 0];
wigner = WignerRepresentation[cat, {-6, 6}, {-2, 2}, "GridSize" -> 90];
husimi = HusimiQRepresentation[cat, {-6, 6}, {-3, 3}, "GridSize" -> 90];
```

Show the two phase-space surfaces side by side:

```wl
Row[{Plot3D[wigner[x, p], {x, -6, 6}, {p, -2, 2}, Mesh -> None,
    PlotRange -> All, ColorFunction -> "SunsetColors", PlotLabel -> "Wigner"],
  Plot3D[husimi[x, p], {x, -6, 6}, {p, -3, 3}, Mesh -> None,
    PlotRange -> All, ColorFunction -> "SunsetColors", PlotLabel -> "Husimi"]}]
```

The Wigner interference fringes dip below zero. The Husimi surface retains the two coherent lobes but smooths those fringes away. Both pictures came from the same state object.

## A symbolic circuit can change its representation

A circuit remains symbolic underneath its diagram. This makes exact circuit identities available to the ordinary algebra of Wolfram Language.

Compose two opposite rotations without assigning a value to $\phi$:

```wl
cancellation = QuantumCircuitOperator[{"RZ"[\[Phi]], "RZ"[-\[Phi]]},
    "Parameters" -> {\[Phi]}]
```

Convert the circuit to its operator and simplify under the physical assumption that $\phi$ is real:

```wl
FullSimplify[Normal[QuantumOperator[cancellation]["Matrix"]],
    \[Phi] \[Element] Reals]
```

The result is the $2\times2$ identity matrix. This is semantic simplification: the equality is proved for the whole parameter family, not inferred from selected angles or from the appearance of a diagram.

The abstraction also separates a circuit from the engine that evaluates it. The default circuit path is a tensor-network contraction, and the network is inspectable.

Expose the tensor network behind the QAOA circuit:

```wl
oneLayer["TensorNetworkGraph"]
```

The gate diagram describes the algorithm. This graph describes the contraction problem that computes it.

For Clifford circuits, the same object can instead compile to an exact stabilizer tableau. A three-qubit cluster state is a compact check.

Build the cluster circuit:

```wl
cluster = QuantumCircuitOperator[
    {"H" -> Range[3], "CZ", "CZ" -> {2, 3}}]
```

Evaluate it through two representations and compare the resulting state vectors:

```wl
stabilizer = cluster[Method -> "Stabilizer"];
QuantumState[stabilizer]["StateVector"] === cluster[]["StateVector"]
```

The result is `True`. The physical object stays fixed while the computational representation changes.

## What the framework changes

The unifying idea is now concrete. A symbolic parameter can begin in a state or circuit and end in a Bloch vector, optimization landscape, entropy, purity loss, or measurement distribution. The same objects survive changes of basis, changes of physical model, and changes of computational engine.

This matters because the final expression is often the scientific result. A formula can be differentiated, optimized, simplified, checked in limiting cases, and reused as input to the rest of Wolfram Language. The workflow remains compact because the objects retain the structure of the physics.

## What this talk does not show

This tour omits several substantial parts of the framework: symbolic Lindblad master equations, measurement records and feedforward, state tomography, explicit contraction-path control, large stabilizer calculations, OpenQASM interchange, device transpilation, and hardware execution. Those belong in focused follow-up sessions.

The live emphasis should remain on the seven transitions shown here. They make one claim repeatedly, in different physics: the framework is broad because its algebra is coherent, not because its function list is long.

## References

[1] Wolfram Research, [*Wolfram Quantum Framework*](https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/), Wolfram Paclet Repository.

[2] E. Farhi, J. Goldstone, and S. Gutmann, [*A Quantum Approximate Optimization Algorithm*](https://arxiv.org/abs/1411.4028), 2014.

[3] S. Wolfram, [*What Is a Computational Essay?*](https://writings.stephenwolfram.com/2017/11/what-is-a-computational-essay/), Wolfram Writings, 2017.
