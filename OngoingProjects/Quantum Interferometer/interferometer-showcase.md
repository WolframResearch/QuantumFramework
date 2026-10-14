---
Template: Default
Name: QuantumInterferometer Showcase
Keywords: [QuantumInterferometer, linear optics, boson sampling, Hong-Ou-Mandel, Gaussian boson sampling, beam splitter, Mach-Zehnder]
---

# QuantumInterferometer: five worked examples

`QuantumInterferometer` describes a linear-optical network: light enters several modes,
passes through beam splitters and phase shifters (and, optionally, squeezers), and photons
are counted at the output. The object stores the mode-to-mode transfer matrix and returns
photon-counting amplitudes exactly, as permanents of its submatrices for a passive network
and as hafnians once squeezers make the network active, with no cutoff imposed on photon
number. Each device below also draws itself with `"Diagram"`: beam splitters are crossings
shaded by reflectivity, phase shifters and squeezers are labeled boxes. Every example starts
from a standard effect in quantum optics.

Load the framework and its second-quantization layer.

```wl
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

## 1. Two-photon interference on a single beam splitter

Send one photon into each port of a beam splitter. Classically each photon leaves through
either port independently, so half the time the two would come out separately, one per port.
Quantum mechanics interferes the two paths that end in a coincidence, and on a balanced
splitter that amplitude cancels exactly: the photons always leave together. This is the
Hong-Ou-Mandel effect.

Build a single beam splitter, keeping the reflectivity angle `θ` (so `θ = π/4` is balanced) and
the relative phase `φ` symbolic.

```wl
hom = QuantumInterferometer[{"BS"[θ, φ] -> {1, 2}}, 2];
```

The diagram is one crossing joining the two modes, labeled with the angle and phase.

```wl
hom["Diagram", "ShowParameters" -> True]
```

Send one photon into each mode; the output is a superposition over the three two-photon
patterns.

```wl
hom[{1, 1}]
```

Read the probability that the two photons come out in separate modes, the coincidence, in
closed form.

```wl
coincidence = FullSimplify[hom["Probability", {1, 1} -> {1, 1}], (θ | φ) ∈ Reals]
```

This is $\cos^2(2\theta)$: independent of the phase, and zero when the splitter is balanced.
Set the splitter to balanced and read the output directly.

```wl
balanced = QuantumInterferometer[{"BS"[π/4, 0] -> {1, 2}}, 2];
```

The picture is the same crossing, now fixed at `θ = π/4`.

```wl
balanced["Diagram", "ShowParameters" -> True]
```

With one photon in each mode the coincidence term is gone: the state is an equal superposition
of both photons in mode 1 and both in mode 2.

```wl
balanced[{1, 1}]
```

*Reference:* C. K. Hong, Z. Y. Ou, L. Mandel, Phys. Rev. Lett. **59**, 2044 (1987).

## 2. Realizing an arbitrary device: three meshes for one unitary

Any passive linear-optical transformation on `m` modes is an `m×m` unitary, and it can be
built physically from two-mode beam splitters and phase shifters laid out in a mesh. Reck's
layout is triangular; Clements' is rectangular and half as deep, which matters on a real chip
because a photon then passes through fewer lossy elements. `QuantumInterferometer` also offers
a cosine-sine recursion.

Draw a Haar-random six-mode unitary to serve as the target transformation.

```wl
SeedRandom[6];
w = RandomVariate[CircularUnitaryMatrixDistribution[6]];
```

Decompose it three ways into meshes of beam splitters and phase shifters.

```wl
meshes = AssociationMap[QuantumInterferometer[w, Method -> #] &, {"Clements", "Reck", "CosineSine"}];
```

The three layouts place the same elements in different patterns.

```wl
Column[KeyValueMap[Labeled[#2["Diagram"], #1, Top] &, meshes]]
```

Each uses the same number of beam splitters and output phases. Compare the depth and confirm
every mesh rebuilds the target to numerical zero.

```wl
Dataset @ Map[
    <|
        "BeamSplitters" -> #["BeamSplitterCount"],
        "Depth" -> #["Depth"],
        "RebuildError" -> Chop @ Max @ Abs[QuantumInterferometer[#["Elements"], 6]["Unitary"] - w]
    |> &,
    meshes
]
```

The rectangular mesh reaches the target in the fewest layers.

*References:* M. Reck et al., Phys. Rev. Lett. **73**, 58 (1994); W. R. Clements et al., Optica
**3**, 1460 (2016); I. Dhand, S. K. Goyal, Phys. Rev. A **92**, 043813 (2015).

## 3. Many-photon interference and forbidden outputs

The two-photon cancellation generalizes. Send one photon into every mode of a balanced
multiport, the discrete Fourier transform (the three-mode case is the "tritter"), and the
symmetry of the Fourier matrix extinguishes whole families of output patterns by exact
destructive interference. Because each probability is a permanent of a submatrix of an exact
matrix, the surviving ones come out as exact rationals rather than floating-point noise.

Build the balanced three-mode splitter from the Fourier matrix.

```wl
tritter = QuantumInterferometer[FourierMatrix[3]];
```

Its diagram is three crossings followed by a column of output phases.

```wl
tritter["Diagram", "ShowParameters" -> True]
```

Send one photon into each mode and read the exact output distribution.

```wl
tritter["Probabilities", {1, 1, 1}]
```

Only the fully bunched outputs (all three photons in one mode) and the fully spread output
survive; every pattern that splits as two-plus-one has probability exactly zero. The effect
sharpens with size. Build the four-mode Fourier transform.

```wl
fourmode = QuantumInterferometer[FourierMatrix[4]];
```

Its mesh is deeper, but built from the same two parts.

```wl
fourmode["Diagram", "ShowParameters" -> True]
```

Form the full output distribution for one photon per mode.

```wl
p = fourmode["Probabilities", {1, 1, 1, 1}];
```

Count how many of the output patterns are forbidden.

```wl
{Count[Values[p], 0], Length[p]}
```

Most of the patterns are extinguished. These exact zeros are a stringent signature of genuine
multi-photon interference, and computing them is exactly the hard permanent problem behind
boson sampling.

*Reference:* M. C. Tichy, K. Mayer, A. Buchleitner, K. Mølmer, Phys. Rev. Lett. **113**,
020502 (2014).

## 4. Squeezed light and Gaussian boson sampling

Adding squeezers makes the network active: it no longer conserves photon number, and the
input is the vacuum. A two-mode squeezer emits photons in pairs, so the counts in the two
modes are perfectly correlated.

Squeeze two modes oppositely and mix them on a balanced splitter.

```wl
tmsv = QuantumInterferometer[{"S"[0.5] -> 1, "S"[-0.5] -> 2, "BS"[π/4, 0] -> {1, 2}}, 2];
```

The diagram shows the two squeezer boxes feeding a beam splitter.

```wl
tmsv["Diagram", "ShowParameters" -> True]
```

Starting from vacuum, read the output distribution and keep the patterns with real weight.

```wl
Select[tmsv["Probabilities", {0, 0}, "Cutoff" -> 6], # > 10^-12 &]
```

Only the balanced patterns `{n, n}` appear: the photons come in pairs. Their weights follow
the thermal law $\tanh^{2n}(r)\,/\cosh^2(r)$ at squeezing `r = 1/2`.

```wl
Table[Tanh[0.5]^(2 n)/Cosh[0.5]^2, {n, 0, 3}]
```

Now place one squeezer on each of four modes.

```wl
SeedRandom[1];
squeezers = QuantumInterferometer[Thread["S"[0.6] -> Range[4]]];
```

The diagram is four squeezer boxes, nothing else yet.

```wl
squeezers["Diagram", "ShowParameters" -> True]
```

Feed them into a random interferometer to make a Gaussian network.

```wl
gbs = QuantumInterferometer["Random"[4]][squeezers];
```

The diagram shows the squeezers followed by the mesh.

```wl
gbs["Diagram"]
```

For such a network the amplitude to detect one photon in every mode is a hafnian of a matrix
read off the network, which is the basis of Gaussian boson sampling.

```wl
hafnian = gbs["Amplitude", {0, 0, 0, 0} -> {1, 1, 1, 1}];
```

Compute the same amplitude a second way: build the Fock-space circuit at six levels per mode
and overlap it with the target state.

```wl
circuit = (SuperDagger[FockState[{1, 1, 1, 1}, 6]] @ gbs["CircuitOperator", 6][])["Scalar"];
```

The two independent routes agree to machine precision.

```wl
{hafnian, circuit, Abs[hafnian - circuit]}
```

*References:* C. S. Hamilton et al., Phys. Rev. Lett. **119**, 170501 (2017); hafnian
algorithm, A. Björklund, B. Gupt, N. Quesada, arXiv:1805.12498.

## 5. A programmable Mach-Zehnder, checked against the full simulation

The building block of programmable photonics is the Mach-Zehnder: two balanced beam splitters
with a phase `φ` between them. The phase steers a single photon from one output toward the
other, a tunable switch.

Build the device, keeping the internal phase `φ` symbolic.

```wl
mz = QuantumInterferometer[{"BS"[π/4, 0] -> {1, 2}, "PS"[φ] -> 1, "BS"[π/4, 0] -> {1, 2}}, 2];
```

The diagram is two crossings with a phase box on the upper arm.

```wl
mz["Diagram", "ShowParameters" -> True]
```

Send a photon into the first mode; the output amplitudes depend on the phase.

```wl
mz[{1, 0}]
```

The routing probabilities are $\sin^2(\varphi/2)$ to stay and $\cos^2(\varphi/2)$ to cross, so
the phase sweeps the photon continuously from one port to the other.

```wl
{
    FullSimplify[mz["Probability", {1, 0} -> {1, 0}], φ ∈ Reals],
    FullSimplify[mz["Probability", {1, 0} -> {0, 1}], φ ∈ Reals]
}
```

Fix the phase to a number to check the many-photon statistics against the full simulation.

```wl
mz1 = QuantumInterferometer[{"BS"[π/4, 0] -> {1, 2}, "PS"[1.0] -> 1, "BS"[π/4, 0] -> {1, 2}}, 2];
```

The diagram now carries the value instead of the symbol.

```wl
mz1["Diagram", "ShowParameters" -> True]
```

Send one photon into each mode and read the output distribution from the permanent.

```wl
permanent = mz1["Probabilities", {1, 1}]
```

Compute the same distribution from the truncated Fock-space circuit.

```wl
circuit = mz1["CircuitOperator", 3][FockState[{1, 1}, 3]]["Probability"]
```

The two distributions agree to machine precision.

```wl
Max @ Abs @ Values @ Merge[{KeyMap[#["Name"] &, permanent], -KeyMap[#["Name"] &, circuit]}, Total]
```

Finally, the dagger of the device inverts it on every mode.

```wl
FullSimplify[mz["Dagger"][mz]["Unitary"], φ ∈ Reals] == IdentityMatrix[2]
```
