---
Template: Default
Name: QuantumInterferometer Showcase
Title: "Quantum Optics with QuantumInterferometer"
Author: Mads Bahrami
Keywords: [QuantumInterferometer, linear optics, coherent states, boson sampling, Hong-Ou-Mandel, squeezed light, beam splitter, Mach-Zehnder, heralding, quantum walks]
---

# Quantum Optics with QuantumInterferometer

Beam splitters, coherent light, multiphoton interference, squeezed states, attenuation,
waveguide arrays, and conditional optical gates: a computation-first tour with the Wolfram Language.

Mads Bahrami (last updated: October 24, 2026)

## Introduction: from optical inputs to measured outputs

These fourteen experiments connect optical inputs to photon-counting probabilities,
coherent amplitudes, and conditional states. Each calculation tests a physical prediction.
Start a fresh Wolfram Language kernel and evaluate the cells in order; later examples
reuse earlier definitions. Change the inputs and parameters to explore the results.

### Prerequisites and setup

Basic complex numbers, matrices, and quantum states are assumed. A **mode** is an independent
optical field component, such as a spatial path. A **Fock state** has definite photon numbers.
Its **occupation pattern**, such as `{2, 0, 1}`, lists the photon number in each mode.
A **coherent state** has a complex coherent amplitude $\alpha$ and mean photon number $|\alpha|^2$.
We assume ideal photon-number-resolving detectors throughout.

`QuantumInterferometer` (QI) represents an optical network of beam splitters, phase shifters,
and optionally squeezers. Passive networks conserve photon number and have a unitary
transfer matrix. Their Fock-state transition amplitudes are normalized permanents;
the squeezed-vacuum examples use hafnians. Transition probabilities are squared moduli
of amplitudes.

Individual transition queries need no Fock-space cutoff. Listing coherent-state or
squeezed-state probabilities requires a total photon-number cutoff. A Fock-space circuit
instead requires a number of levels per mode. Passive `"Amplitude"` and `"Probability"`
queries accept symbolic occupations; `"Probabilities"` requires explicit occupations.
Squeezed networks require explicit occupations but allow symbolic squeezing parameters.
Symbolic occupation formulas retain their nonnegative-integer domain; explicit assumptions
let `Refine` display the formula on that domain.

Load the quantum-state and operator tools:

```wl
Needs["Wolfram`QuantumFramework`"]
```

Load the second-quantization tools:

```wl
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

### Conventions: transfer matrices, beam splitter phases, and mode ordering

For a passive network, $W_{jk}=\langle1_j|U|1_k\rangle$: column $k$ gives the
single-photon amplitudes for input mode $k$. We use **transmissivity** $T=\cos^2\theta$
and **reflectivity** $R=\sin^2\theta$. These are probabilities; transmission and reflection
amplitudes also carry phases. Beam splitter angles and phase shifts are real;
coherent amplitudes may be complex. A phase shifter contributes $e^{i\phi}$ per photon.

Build a beam splitter with zero reflection phase:

```wl
conventionBS = QuantumInterferometer[{"BS"[theta, 0] -> {1, 2}}, 2];
```

Draw the beam splitter with its parameters:

```wl
conventionBS["Diagram", "ShowParameters" -> True]
```

Read its transfer matrix:

```wl
conventionBS["Unitary"]
```

The reflected amplitudes have opposite real signs. A reflection phase of `Pi/2`
gives the symmetric convention used in many textbooks. Build that beam splitter:

```wl
symmetricConventionBS = QuantumInterferometer[{"BS"[theta, Pi/2] -> {1, 2}}, 2];
```

Draw the beam splitter in the symmetric convention:

```wl
symmetricConventionBS["Diagram", "ShowParameters" -> True]
```

Read its transfer matrix:

```wl
symmetricConventionBS["Unitary"]
```

The symmetric convention has imaginary reflection amplitudes. We state any convention
change explicitly; some texts also reverse the phase-shift sign. All optical parameters
here are dimensionless. Waveguide propagation lengths are measured in inverse coupling units.

## Part I: Optical components and photon statistics

### 1. Two-photon interference: indistinguishability and coincidences

Send one photon into each input mode. Assume the photons have identical polarization,
spectrum, and temporal shape. At a balanced beam splitter, the two amplitudes for one photon
in each output cancel. Both photons leave in the same output mode: the Hong-Ou-Mandel effect.
Distinguishable photons instead have coincidence probability $1/2$.

Build a beam splitter with symbolic angle `theta` and reflection phase `phi`:

```wl
hom = QuantumInterferometer[{"BS"[theta, phi] -> {1, 2}}, 2];
```

Draw the beam splitter:

```wl
hom["Diagram", "ShowParameters" -> True]
```

Compute the output Fock-state amplitudes for one photon in each input:

```wl
hom[{1, 1}]
```

The coincidence outcome is `{1, 1}`. Compute its probability with
[`FullSimplify`](https://reference.wolfram.com/language/ref/FullSimplify.html), assuming real parameters:

```wl
coincidence = FullSimplify[hom["Probability", {1, 1} -> {1, 1}], (theta | phi) ∈ Reals]
```

The probability is $\cos^2(2\theta)$, independent of the reflection phase.
Build the balanced beam splitter at $\theta=\pi/4$:

```wl
balanced = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2];
```

Draw the balanced beam splitter:

```wl
balanced["Diagram", "ShowParameters" -> True]
```

The coincidence amplitude vanishes. Compute the two surviving bunching amplitudes:

```wl
balanced[{1, 1}]
```

Indistinguishable alternatives interfere as amplitudes. Distinguishable alternatives
contribute separate probabilities. Verify the indistinguishable coincidence probability:

```wl
balanced["Probability", {1, 1} -> {1, 1}] == 0
```

For the distinguishable reference, read the transfer matrix:

```wl
balancedMatrix = balanced["Unitary"]
```

Add the probabilities for both photons being transmitted or both being reflected:

```wl
Abs[balancedMatrix[[1, 1]] balancedMatrix[[2, 2]]]^2 +
 Abs[balancedMatrix[[1, 2]] balancedMatrix[[2, 1]]]^2 == 1/2
```

Distinguishability restores coincidences. A delay scan requires a model of the temporal
overlap, which is absent from these fixed-mode inputs.

*Reference:* C. K. Hong, Z. Y. Ou, L. Mandel, Phys. Rev. Lett. **59**, 2044 (1987).

### 2. Optical meshes: three realizations of one unitary

A passive network on $m$ modes has an $m\times m$ unitary transfer matrix. Beam splitters
and phase shifters can realize it in several layouts: triangular Reck, rectangular Clements,
or recursive cosine-sine meshes. Compare their component counts, depths, and reconstruction errors.

Set a reproducible random seed:

```wl
SeedRandom[6];
```

Generate a six-mode Haar-random target unitary:

```wl
targetUnitary = RandomVariate[CircularUnitaryMatrixDistribution[6]];
```

Decompose the target into the three meshes:

```wl
meshes = AssociationMap[
 QuantumInterferometer[targetUnitary, Method -> #] &,
 {"Clements", "Reck", "CosineSine"}];
```

Draw the three meshes:

```wl
Column[KeyValueMap[Labeled[#2["Diagram"], #1, Top] &, meshes]]
```

Count the beam splitters in each mesh:

```wl
Map[#["BeamSplitterCount"] &, meshes]
```

Depth counts beam splitter layers. Compare the mesh depths:

```wl
Map[#["Depth"] &, meshes]
```

Check the decomposition by rebuilding each mesh from its optical elements.
This tests the elements rather than the stored target matrix. Build the reconstructed networks:

```wl
rebuiltMeshes = Map[QuantumInterferometer[#["Elements"], 6] &, meshes];
```

Draw the reconstructed meshes:

```wl
Column[KeyValueMap[Labeled[#2["Diagram"], #1, Top] &, rebuiltMeshes]]
```

[`Chop`](https://reference.wolfram.com/language/ref/Chop.html) replaces small numerical
residuals with zero. Compute the largest matrix-entry error for each reconstruction:

```wl
Map[Chop[Max[Abs[#["Unitary"] - targetUnitary]]] &, rebuiltMeshes]
```

All three meshes reproduce the target. The Clements mesh uses the fewest layers in this example.

*References:* M. Reck et al., Phys. Rev. Lett. **73**, 58 (1994); W. R. Clements et al., Optica
**3**, 1460 (2016); I. Dhand, S. K. Goyal, Phys. Rev. A **92**, 043813 (2015).

### 3. Fourier multiports: symmetry and suppressed outputs

A Fourier multiport realizes the discrete Fourier transform. With one photon per input,
Fourier symmetry makes certain transition amplitudes vanish. Exact matrix entries let us
identify these suppressed occupation patterns without numerical tolerances.

Build the three-mode Fourier multiport, called a tritter, with
[`FourierMatrix`](https://reference.wolfram.com/language/ref/FourierMatrix.html):

```wl
tritter = QuantumInterferometer[FourierMatrix[3]];
```

Draw the tritter:

```wl
tritter["Diagram", "ShowParameters" -> True]
```

Compute the exact output distribution for one photon per input:

```wl
tritter["Probabilities", {1, 1, 1}]
```

Only `{1, 1, 1}` and the three fully bunched outcomes have nonzero probability.
All permutations of `{2, 1, 0}` are suppressed. Build the four-mode Fourier multiport:

```wl
fourier4 = QuantumInterferometer[FourierMatrix[4]];
```

Draw the four-mode network:

```wl
fourier4["Diagram", "ShowParameters" -> True]
```

Compute its output distribution for one photon per input:

```wl
fourierProbabilities = fourier4["Probabilities", {1, 1, 1, 1}];
```

Count the suppressed occupation patterns:

```wl
Count[Values[fourierProbabilities], 0]
```

Count all occupation patterns with four photons in four modes:

```wl
Length[fourierProbabilities]
```

Fourier symmetry suppresses 25 of the 35 patterns. These exact zeros expose interference
that is harder to analyze in a generic network.

*Reference:* M. C. Tichy, K. Mayer, A. Buchleitner, K. Mølmer, Phys. Rev. Lett. **113**,
020502 (2014).

### 4. Squeezed light: pair creation and Gaussian photon statistics

Squeezing makes a network active: energy exchange with a pump allows photon-pair creation
and annihilation. Opposite single-mode squeezers followed by a balanced beam splitter
prepare a two-mode squeezed vacuum. The two output modes have equal photon numbers.

Build this preparation with a symbolic squeezing parameter `r`:

```wl
symbolicPair = QuantumInterferometer[
 {"S"[r] -> 1, "S"[-r] -> 2, "BS"[Pi/4, 0] -> {1, 2}}, 2];
```

Draw the squeezed-vacuum preparation:

```wl
symbolicPair["Diagram", "ShowParameters" -> True]
```

Compute the probabilities of zero, one, and two pairs for positive `r`:

```wl
pairLawValues = FullSimplify[Table[
 symbolicPair["Probability", {0, 0} -> {n, n}], {n, 0, 2}], r > 0]
```

The pair-number distribution is geometric: $P(n,n)=\tanh^{2n}r/\cosh^2r$.
Verify the first three probabilities:

```wl
FullSimplify[
 pairLawValues == Table[Tanh[r]^(2 n)/Cosh[r]^2, {n, 0, 2}], r > 0]
```

Verify that the infinite pair-number distribution is normalized:

```wl
FullSimplify[Sum[Tanh[r]^(2 n)/Cosh[r]^2, {n, 0, Infinity}], r > 0] == 1
```

Check that zero squeezing gives vacuum:

```wl
Limit[pairLawValues, r -> 0] == {1, 0, 0}
```

Set the squeezing parameter to `0.5` for a numerical comparison:

```wl
tmsv = QuantumInterferometer[
 {"S"[0.5] -> 1, "S"[-0.5] -> 2, "BS"[Pi/4, 0] -> {1, 2}}, 2];
```

Draw this squeezed-vacuum preparation:

```wl
tmsv["Diagram", "ShowParameters" -> True]
```

List probabilities through six total photons, omitting numerical residuals below `10^-12`:

```wl
Select[tmsv["Probabilities", {0, 0}, "Cutoff" -> 6], # > 10^-12 &]
```

Only patterns `{n, n}` occur. Either mode alone has a thermal photon-number distribution.
Compute the geometric probabilities for the retained pair numbers:

```wl
Table[Tanh[0.5]^(2 n)/Cosh[0.5]^2, {n, 0, 3}]
```

For a four-mode example, prepare one squeezed vacuum in each mode:

```wl
squeezers = QuantumInterferometer[Thread["S"[0.6] -> Range[4]]];
```

Draw the four squeezers:

```wl
squeezers["Diagram", "ShowParameters" -> True]
```

Set the random seed for the passive network:

```wl
SeedRandom[1];
```

Build a random four-mode passive network:

```wl
randomMixer = QuantumInterferometer["Random"[4]];
```

Draw the passive network:

```wl
randomMixer["Diagram"]
```

Apply the passive network after the squeezers:

```wl
gbs = randomMixer[squeezers];
```

Draw the combined network:

```wl
gbs["Diagram"]
```

The vacuum-to-four-photon amplitude is a normalized hafnian. Such amplitudes determine
the probabilities in Gaussian boson sampling. Compute the amplitude for `{1, 1, 1, 1}`:

```wl
gaussianAmplitude = gbs["Amplitude", {0, 0, 0, 0} -> {1, 1, 1, 1}]
```

Build the Fock-space circuit with six levels per mode:

```wl
gaussianFockCircuit = gbs["CircuitOperator", 6];
```

Apply the circuit to vacuum:

```wl
gaussianFockOutput = gaussianFockCircuit[];
```

Compute the overlap with the four-photon target:

```wl
fockAmplitude = (SuperDagger[FockState[{1, 1, 1, 1}, 6]] @ gaussianFockOutput)["Scalar"]
```

Compare the two amplitudes:

```wl
Abs[gaussianAmplitude - fockAmplitude]
```

The difference is at machine precision. The hafnian calculation avoids constructing
the full Fock-space circuit.

*References:* C. S. Hamilton et al., Phys. Rev. Lett. **119**, 170501 (2017); hafnian
algorithm, A. Björklund, B. Gupt, N. Quesada, arXiv:1805.12498.

### 5. Mach-Zehnder interferometry: controlling the output with a phase

Two balanced beam splitters and an internal phase shift form a Mach-Zehnder interferometer.
The phase controls the single-photon output probabilities. Build it with symbolic phase `phi`:

```wl
mz = QuantumInterferometer[
 {"BS"[Pi/4, 0] -> {1, 2}, "PS"[phi] -> 1, "BS"[Pi/4, 0] -> {1, 2}}, 2];
```

Draw the Mach-Zehnder interferometer:

```wl
mz["Diagram", "ShowParameters" -> True]
```

Compute the output amplitudes for a photon in the first input:

```wl
mz[{1, 0}]
```

The first-output probability is $\sin^2(\phi/2)$. Compute it for a real phase:

```wl
FullSimplify[mz["Probability", {1, 0} -> {1, 0}], phi ∈ Reals]
```

The second-output probability is $\cos^2(\phi/2)$. Compute it:

```wl
FullSimplify[mz["Probability", {1, 0} -> {0, 1}], phi ∈ Reals]
```

Set the phase to `1.0` for a Fock-space comparison:

```wl
mz1 = QuantumInterferometer[
 {"BS"[Pi/4, 0] -> {1, 2}, "PS"[1.0] -> 1, "BS"[Pi/4, 0] -> {1, 2}}, 2];
```

Draw the interferometer with this phase:

```wl
mz1["Diagram", "ShowParameters" -> True]
```

Compute the two-photon output distribution using permanents:

```wl
mzProbabilities = mz1["Probabilities", {1, 1}]
```

Build the Fock-space circuit with three levels per mode:

```wl
mzFockCircuit = mz1["CircuitOperator", 3];
```

Propagate one photon in each input mode:

```wl
mzFockOutput = mzFockCircuit[FockState[{1, 1}, 3]];
```

Read the output probabilities:

```wl
fockProbabilities = mzFockOutput["Probability"]
```

Subtract the distributions by matching their occupation patterns:

```wl
mzProbabilityDifferences = Merge[{mzProbabilities, -fockProbabilities}, Total];
```

Compute the largest absolute probability difference:

```wl
Max[Abs[Values[mzProbabilityDifferences]]]
```

The distributions agree to machine precision. Construct the inverse interferometer:

```wl
inverseMZ = mz["Dagger"];
```

Draw the inverse, which reverses the element order and inverts each element:

```wl
inverseMZ["Diagram", "ShowParameters" -> True]
```

Apply the inverse after the original interferometer:

```wl
roundTripMZ = inverseMZ[mz];
```

Draw the complete sequence:

```wl
roundTripMZ["Diagram", "ShowParameters" -> True]
```

Verify that the combined transfer matrix is the identity:

```wl
FullSimplify[roundTripMZ["Unitary"], phi ∈ Reals] == IdentityMatrix[2]
```

## Part II: Input states and multiphoton interference

Changing the input state changes the interference and photon statistics, even in the same network.

### 6. Coherent inputs: photon statistics and dark ports

A coherent state with $\alpha=1$ and a one-photon Fock state both have mean photon number one.
The coherent state's photon number fluctuates; the Fock state's does not.

`"Coherent"[{a1, a2}]` specifies coherent amplitudes. Their squared moduli sum to the mean
total photon number, so this vector need not have unit norm. A plain list specifies a
Fock-state occupation pattern instead.

Compute the output coherent amplitudes for input `"Coherent"[{1, 0}]`:

```wl
balanced["Coherent"[{1, 0}]]
```

Compare the output distribution for the one-photon Fock input `{1, 0}`:

```wl
balanced["Probabilities", {1, 0}]
```

The coherent input can produce vacuum. Verify its vacuum probability:

```wl
balanced["Probability", "Coherent"[{1, 0}] -> {0, 0}] == Exp[-1]
```

Compute the three output probabilities with exactly two photons:

```wl
coherentTwoPhotonProbabilities =
 balanced["Probability", "Coherent"[{1, 0}] -> #] & /@ {{2, 0}, {1, 1}, {0, 2}};
```

Add them to obtain the total two-photon probability:

```wl
twoPhotonProbability = Total[coherentTwoPhotonProbabilities]
```

Verify the two-photon probability for a Poisson distribution with mean one:

```wl
twoPhotonProbability == Exp[-1]/2
```

Photon-number conservation forbids vacuum from the one-photon Fock input. Check it:

```wl
balanced["Probability", {1, 0} -> {0, 0}] == 0
```

It also forbids a two-photon output. Check one such outcome:

```wl
balanced["Probability", {1, 0} -> {2, 0}] == 0
```

A coherent-state distribution has an infinite tail. A finite cutoff omits probability;
it does not justify renormalization. List the probabilities through six total photons:

```wl
coherentProbabilities = balanced["Probabilities", "Coherent"[{1, 0}], "Cutoff" -> 6];
```

Verify that the retained probability equals the Poisson cumulative probability:

```wl
Total[coherentProbabilities] == CDF[PoissonDistribution[1], 6]
```

Equal coherent amplitudes interfere destructively in one output, making it a **dark port**.
Compute the output coherent amplitudes:

```wl
balanced["Coherent"[{1, 1}]]
```

Verify that the dark port cannot contain a photon:

```wl
balanced["Probability", "Coherent"[{1, 1}] -> {1, 0}] == 0
```

The Fock input `{1, 1}` instead bunches into either output. Compute the probability
that both photons leave through the port that was dark for the coherent input:

```wl
balanced["Probability", {1, 1} -> {2, 0}]
```

Coherent-state interference depends on the relative input phase. Give both inputs coherent
amplitude `alphaLaser` and add phase `laserPhase` to the second. Compute the output amplitudes:

```wl
coherentOutput = balanced["Coherent"[{alphaLaser, alphaLaser Exp[I laserPhase]}]]
```

`First` extracts the coherent amplitude vector from the `"Coherent"` expression. Its squared moduli
are the mean output photon numbers, $|\alpha|^2\{1-\cos\chi,1+\cos\chi\}$, for relative
phase $\chi$. Verify this law:

```wl
FullSimplify[Abs[First[coherentOutput]]^2 ==
 Abs[alphaLaser]^2 {1 - Cos[laserPhase], 1 + Cos[laserPhase]}, laserPhase ∈ Reals]
```

The total mean photon number stays fixed. Compare the output amplitudes at relative phases zero and `Pi`:

```wl
AssociationMap[FullSimplify[coherentOutput /. laserPhase -> #] &, {0, Pi}]
```

Changing the relative phase exchanges the dark port while preserving the input mean photon number.

*Source:* Agarwal, *Quantum Optics*, Sections 5.3-5.5, pp. 105-107; Gerry and Knight,
*Introductory Quantum Optics*, Section 6.2.

### 7. Twin Fock inputs: why odd photon numbers are suppressed

A twin Fock input has equal photon numbers in the two modes. Start with three photons
per input of the balanced beam splitter. Compute the output distribution:

```wl
twinProbabilities = balanced["Probabilities", {3, 3}]
```

Every odd output photon number is suppressed. To examine arbitrary inputs, use
`{nTwin, nTwin}` and the photon-number-conserving output `{mTwin, 2 nTwin - mTwin}`.
The physical domain is integer `nTwin >= 0` and `0 <= mTwin <= 2 nTwin`.

Store the symbolic transition amplitude with its occupation domain:

```wl
symbolicTwinAmplitude = balanced[
 "Amplitude", {nTwin, nTwin} -> {mTwin, 2 nTwin - mTwin}];
```

[`Refine`](https://reference.wolfram.com/language/ref/Refine.html) removes the domain
conditions when the occupations are assumed physical. Display the amplitude on this domain:

```wl
Refine[symbolicTwinAmplitude,
 Element[{nTwin, mTwin}, Integers] && nTwin >= 0 && 0 <= mTwin <= 2 nTwin]
```

This is the [`WignerD`](https://reference.wolfram.com/language/ref/WignerD.html) rotation
matrix element $d^n_{0,m-n}(\pi/2)$ in Wolfram's convention. The stored expression still
returns zero outside its nonnegative-integer domain. Check negative and fractional substitutions:

```wl
AssociationMap[symbolicTwinAmplitude /. {nTwin -> #, mTwin -> 0} &, {-1, 1/2}]
```

Substitute three photons per input and order the amplitudes by the first output's photon number:

```wl
twinSymbolicSix = Table[
 symbolicTwinAmplitude /. {nTwin -> 3, mTwin -> q}, {q, 0, 6}]
```

The amplitudes retain relative signs. Verify that their squared moduli give the direct QI probabilities:

```wl
Abs[twinSymbolicSix]^2 == Table[
 balanced["Probability", {3, 3} -> {q, 6 - q}], {q, 0, 6}]
```

The general odd-output `WignerD` expression does not simplify automatically in this kernel.
A creation-operator argument gives the parity rule directly.

Let $z_1,z_2$ represent the output creation operators. One input pair contributes
$(z_1+z_2)(z_2-z_1)/2=(z_2^2-z_1^2)/2$. Apart from input Fock normalization,
$n$ pairs therefore give $(z_2^2-z_1^2)^n/2^n$. Only even powers occur,
so neither output can have an odd photon number.

Expand the six-photon polynomial:

```wl
Expand[(z2^2 - z1^2)^3/2^3]
```

Verify the odd-output suppression in the QI distribution:

```wl
AllTrue[KeySelect[twinProbabilities, OddQ[First[#["Name"]]] &], PossibleZeroQ]
```

Verify normalization of the complete distribution:

```wl
Total[twinProbabilities] == 1
```

For input `{n, n}`, the coefficient of $z_1^{2k}z_2^{2n-2k}$ is
$(-1)^k\binom nk/2^n$. Multiplying by $\sqrt{(2k)!(2n-2k)!}/n!$ gives the normalized
Fock-state amplitude. Its squared modulus is $\binom{2k}{k}\binom{2(n-k)}{n-k}/4^n$.

Form the normalized amplitude for the even-output occupation pattern:

```wl
twinEvenAmplitude = (-1)^kTwin Binomial[nTwin, kTwin] Sqrt[
 (2 kTwin)! (2 nTwin - 2 kTwin)!]/(2^nTwin nTwin!);
```

Compute its probability on the physical integer domain:

```wl
twinEvenProbability = FullSimplify[Abs[twinEvenAmplitude]^2,
 Element[{nTwin, kTwin}, Integers] && 0 <= kTwin <= nTwin];
```

Verify the binomial form of this probability:

```wl
FullSimplify[
 twinEvenProbability == Binomial[2 kTwin, kTwin] Binomial[2 (nTwin - kTwin), nTwin - kTwin]/4^nTwin,
 Element[{nTwin, kTwin}, Integers] && 0 <= kTwin <= nTwin]
```

Verify normalization for every nonnegative integer input using
[`Sum`](https://reference.wolfram.com/language/ref/Sum.html):

```wl
Sum[Binomial[2 kTwin, kTwin] Binomial[2 (nTwin - kTwin), nTwin - kTwin]/4^nTwin,
 {kTwin, 0, nTwin}, Assumptions -> Element[nTwin, Integers] && nTwin >= 0] == 1
```

Define the even-output probability on its physical domain:

```wl
twinLaw[n_Integer?NonNegative, k_Integer?NonNegative] /; k <= n :=
 Binomial[2 k, k] Binomial[2 (n - k), n - k]/4^n
```

Check this law against the six-photon QI probabilities:

```wl
Table[balanced["Probability", {3, 3} -> {2 k, 6 - 2 k}] == twinLaw[3, k],
 {k, 0, 3}]
```

Check the symbolic QI amplitude against the law for zero through five photons per input:

```wl
Table[
 Abs[symbolicTwinAmplitude /. {nTwin -> n, mTwin -> 2 k}]^2 == twinLaw[n, k],
 {n, 0, 5}, {k, 0, n}]
```

Check normalization independently using the complete QI distributions over this range:

```wl
Table[Total[balanced["Probabilities", {n, n}]] == 1, {n, 0, 5}]
```

The range includes vacuum and the original two-photon case.

*Source:* Agarwal, Section 5.7, p. 109. The accompanying errata corrects the paths in Figure 5.5.

### 8. Four-photon interference: tuning suppression with reflectivity

At a balanced beam splitter, `{1, 3}` cannot produce `{2, 2}`. Verify this suppression:

```wl
balanced["Probability", {1, 3} -> {2, 2}] == 0
```

For input `{2, 2}`, suppressing the same output requires a different reflectivity.
Build a beam splitter with symbolic angle:

```wl
fourBS = QuantumInterferometer[{"BS"[theta, 0] -> {1, 2}}, 2];
```

Draw the beam splitter:

```wl
fourBS["Diagram", "ShowParameters" -> True]
```

Compute the return amplitude for `{2, 2}`:

```wl
fourAmplitude = FullSimplify[
 fourBS["Amplitude", {2, 2} -> {2, 2}], theta ∈ Reals]
```

Verify its expression in terms of transmission and reflection amplitudes:

```wl
FullSimplify[fourAmplitude == 1 - 6 Sin[theta]^2 Cos[theta]^2, theta ∈ Reals]
```

Writing $R=\sin^2\theta$ gives return amplitude $1-6R(1-R)$.
Find the physical reflectivities at which it vanishes:

```wl
suppressionSettings = reflectivity /. Solve[
 1 - 6 reflectivity (1 - reflectivity) == 0 && 0 < reflectivity < 1,
 reflectivity, Reals]
```

Build a beam splitter at each suppression setting:

```wl
suppressionDevices = QuantumInterferometer[
 {"BS"[ArcSin[Sqrt[#]], 0] -> {1, 2}}, 2] & /@ suppressionSettings;
```

Draw both beam splitters:

```wl
Column[#["Diagram", "ShowParameters" -> True] & /@ suppressionDevices]
```

Verify zero return probability at both settings:

```wl
RootReduce[#["Probability", {2, 2} -> {2, 2}]] & /@ suppressionDevices
```

The settings exchange transmissivity and reflectivity. At either endpoint, the return
probability is one. Plot the return probability across the full angular range:

```wl
Plot[Evaluate[fourAmplitude^2], {theta, 0, Pi/2},
 Frame -> True, FrameLabel -> {"Beam splitter angle", "P({2,2} to {2,2})"},
 PlotLabel -> "Four-photon suppression at two unbalanced settings",
 GridLines -> Automatic, PlotRange -> {0, 1}, AspectRatio -> 1/2, ImageSize -> 480]
```

The zeros depend on the input photon number. Compute the return amplitude for symbolic
twin Fock input `{nReturn, nReturn}`:

```wl
returnAmplitudeFamily = FullSimplify[
 fourBS["Amplitude", {nReturn, nReturn} -> {nReturn, nReturn}],
 Element[nReturn, Integers] && nReturn >= 0 && theta ∈ Reals]
```

Evaluate the amplitude for one, two, and three photons per input:

```wl
AssociationMap[
 FullSimplify[returnAmplitudeFamily /. nReturn -> #, theta ∈ Reals] &, {1, 2, 3}]
```

The different angular polynomials have different zeros. Verify the four-photon member:

```wl
FullSimplify[(returnAmplitudeFamily /. nReturn -> 2) == fourAmplitude, theta ∈ Reals]
```

*Source:* Agarwal, Exercises 5.6 and 5.9, pp. 134-135.

### 9. A photon and coherent light: preparing a state by heralding

Coherent inputs remain coherent under passive optics. Adding a Fock-state input allows
nonclassical state preparation by **heralding**: a specified detection outcome in one mode
selects a conditional state in the other.

Build a beam splitter for one photon in the first input and coherent amplitude `beta` in the second:

```wl
mixedBS = QuantumInterferometer[{"BS"[theta, 0] -> {1, 2}}, 2];
```

Draw the beam splitter; the first output will be measured for the heralding outcome:

```wl
mixedBS["Diagram", "ShowParameters" -> True]
```

Keep `beta` complex. For one photon detected in the first output, compute the vacuum
and one-photon amplitudes in the second:

```wl
heraldLeading = FullSimplify[
 Table[mixedBS["Amplitude", {1, "Coherent"[beta]} -> {1, k}], {k, 0, 1}], theta ∈ Reals]
```

Verify the two unnormalized conditional amplitudes:

```wl
FullSimplify[heraldLeading == Exp[-Abs[beta]^2/2] {Cos[theta], beta Cos[2 theta]},
 theta ∈ Reals]
```

The coherent amplitude controls the relative magnitude and phase of these coefficients.
Verify that a phase shift of `beta` rotates only the one-photon coefficient:

```wl
FullSimplify[(heraldLeading /. beta -> beta Exp[I heraldDrivePhase]) ==
 heraldLeading {1, Exp[I heraldDrivePhase]}, Element[{theta, heraldDrivePhase}, Reals]]
```

Small transmission and coherent amplitudes suppress higher photon numbers. Set the
transmission amplitude to `1/5`:

```wl
heraldBS = QuantumInterferometer[{"BS"[ArcCos[1/5], 0] -> {1, 2}}, 2];
```

Draw the preparation network:

```wl
heraldBS["Diagram", "ShowParameters" -> True]
```

Set `beta = 1/5` and list joint probabilities through six total photons:

```wl
mixedProbabilities = heraldBS["Probabilities", {1, "Coherent"[1/5]}, "Cutoff" -> 6];
```

Select the joint outcomes with one photon in the first output:

```wl
N[KeySelect[mixedProbabilities, First[#["Name"]] == 1 &]]
```

These are joint probabilities before conditioning. To construct the conditional state,
retain the corresponding amplitudes for zero through eight photons in the second output:

```wl
heraldAmplitudes = Table[
 heraldBS["Amplitude", {1, "Coherent"[1/5]} -> {1, k}], {k, 0, 8}];
```

Their squared norm is the retained success probability. Compute it:

```wl
heraldSuccess = Total[Abs[heraldAmplitudes]^2];
```

Evaluate the retained success probability numerically:

```wl
N[heraldSuccess]
```

Normalize the selected amplitudes to construct the conditional state:

```wl
heraldState = QuantumState[heraldAmplitudes/Sqrt[heraldSuccess], 9];
```

Verify normalization of the conditional state:

```wl
FullSimplify[heraldState["Norm"]] == 1
```

The leading coefficients approximate $(|0\rangle-|1\rangle)/\sqrt2$.
Compute the fidelity (squared overlap) with this target state:

```wl
N[Abs[({1, -1} . Take[Normal[heraldState["StateVector"]], 2])/Sqrt[2]]^2]
```

The fidelity applies to the retained state. Bound the omitted joint probability
using the coherent input's Poisson tail:

```wl
N[1 - CDF[PoissonDistribution[1/25], 8], 20]
```

Photon-number conservation gives this bound: with one herald photon, an omitted output
requires more than eight photons from the coherent input. Divide by the retained success
probability for an upper bound on the conditional tail.

At zero coherent amplitude, the heralding outcome leaves vacuum in the second output.
Check its joint probability:

```wl
heraldBS["Probability", {1, "Coherent"[0]} -> {1, 0}] == 1/25
```

Verify that no photon can remain in the second output:

```wl
heraldBS["Probability", {1, "Coherent"[0]} -> {1, 1}] == 0
```

*Source:* Agarwal, Section 5.9 and Figure 5.6, pp. 110-111. Relative signs follow our
real beam splitter convention.

### 10. Multiphoton interferometry: fringes and factorial moments

A fringe is the phase dependence of a specified outcome for a specified input.
Compute the single-photon return probability in the Mach-Zehnder interferometer:

```wl
fringeOne = FullSimplify[mz["Probability", {1, 0} -> {1, 0}], phi ∈ Reals]
```

Verify the single-photon fringe:

```wl
FullSimplify[fringeOne == Sin[phi/2]^2, phi ∈ Reals]
```

Compute the coincidence fringe for one photon per input:

```wl
fringeTwo = FullSimplify[mz["Probability", {1, 1} -> {1, 1}], phi ∈ Reals]
```

Verify the two-photon fringe:

```wl
FullSimplify[fringeTwo == Cos[phi]^2, phi ∈ Reals]
```

Compute the probability of output `{3, 1}` for input `{2, 2}`:

```wl
fringeFour = FullSimplify[mz["Probability", {2, 2} -> {3, 1}], phi ∈ Reals]
```

Verify the four-photon fringe:

```wl
FullSimplify[fringeFour == 3 Sin[2 phi]^2/8, phi ∈ Reals]
```

The four-photon fringe has a shorter period and a lower detection probability.
A shorter period alone does not establish better phase sensitivity; fluctuations
and photon resources also matter. Plot the probabilities without rescaling:

```wl
Plot[Evaluate[{fringeOne, fringeTwo, fringeFour}], {phi, 0, 2 Pi},
 Frame -> True, FrameLabel -> {"Internal phase", "Outcome probability"},
 PlotLegends -> {"{1,0} to {1,0}", "{1,1} to {1,1}", "{2,2} to {3,1}"},
 PlotLabel -> "Multiphoton outcome probabilities", GridLines -> Automatic,
 PlotRange -> {0, 1}, AspectRatio -> 1/2, ImageSize -> 540]
```

For a three-photon example, use reflectivity $2/3$ at both beam splitters,
the symmetric reflection-phase convention, and phase factor $e^{-i\varphi}$ in the second arm.
Build the interferometer:

```wl
threeMZ = QuantumInterferometer[{
 "BS"[ArcSin[Sqrt[2/3]], Pi/2] -> {1, 2},
 "PS"[-phase3] -> 2,
 "BS"[ArcSin[Sqrt[2/3]], Pi/2] -> {1, 2}}, 2];
```

Draw the three-photon interferometer:

```wl
threeMZ["Diagram", "ShowParameters" -> True]
```

Compute the probability of `{2, 1}` becoming `{1, 2}`:

```wl
threeFringe = FullSimplify[
 threeMZ["Probability", {2, 1} -> {1, 2}], phase3 ∈ Reals]
```

Verify the fringe with three times the internal phase:

```wl
FullSimplify[threeFringe == 16 (1 + Cos[3 phase3])/81, phase3 ∈ Reals]
```

The factorial moment $\langle n_1n_2(n_2-1)\rangle$ counts ordered detection triples.
In this three-photon sector, only `{1, 2}` contributes, with factor two.
Compute the complete output distribution:

```wl
threeProbabilities = threeMZ["Probabilities", {2, 1}];
```

Verify normalization before using the probabilities:

```wl
FullSimplify[ComplexExpand[Total[threeProbabilities]]] == 1
```

Compute the factorial moment from the occupation patterns:

```wl
factorialMoment = Total[KeyValueMap[
 With[{n = #1["Name"]}, #2 n[[1]] n[[2]] (n[[2]] - 1)] &, threeProbabilities]];
```

Verify that the moment is twice the selected outcome probability:

```wl
FullSimplify[ComplexExpand[factorialMoment - 2 threeFringe]] == 0
```

*Source:* Agarwal, Sections 5.12, 5.15-5.16, pp. 113-123, and Exercise 5.10.

## Part III: State preparation, loss, and photonic computation

We now reverse a preparation, discard an output mode, and compose larger optical networks.

### 11. Squeezed-state interference: disentangling two modes

`symbolicPair` prepares a two-mode squeezed vacuum by mixing opposite single-mode squeezed
vacua. Reversing its beam splitter should recover the independent squeezed vacua.
Keep `r > 0` symbolic and construct the inverse beam splitter:

```wl
inverseBalanced = balanced["Dagger"];
```

Draw the inverse beam splitter:

```wl
inverseBalanced["Diagram", "ShowParameters" -> True]
```

Apply it after the two-mode squeezed-vacuum preparation:

```wl
separatedSqueezers = inverseBalanced[symbolicPair];
```

Draw the full sequence:

```wl
separatedSqueezers["Diagram", "ShowParameters" -> True]
```

Build the direct preparation with the same squeezing parameter:

```wl
directSqueezers = QuantumInterferometer[{"S"[r] -> 1, "S"[-r] -> 2}, 2];
```

Draw the direct preparation:

```wl
directSqueezers["Diagram", "ShowParameters" -> True]
```

Choose output occupations zero through six in each mode:

```wl
squeezePatterns = Tuples[Range[0, 6], 2];
```

Subtract the preparation amplitudes for each occupation pattern:

```wl
squeezeAmplitudeDifferences =
 (separatedSqueezers["Amplitude", {0, 0} -> #] -
  directSqueezers["Amplitude", {0, 0} -> #]) & /@ squeezePatterns;
```

Find the largest absolute difference for positive squeezing:

```wl
separationError = FullSimplify[Max[Abs[squeezeAmplitudeDifferences]], r > 0]
```

Verify exact agreement for positive squeezing:

```wl
separationError == 0
```

The independent squeezed vacua have factorized probabilities. Each has
$P(2k)=\binom{2k}{k}\tanh^{2k}r/(4^k\cosh r)$ and $P(2k+1)=0$.
Define this single-mode probability, including the vacuum case:

```wl
singleSqueezedProbability[n_Integer?NonNegative, r_] :=
 If[EvenQ[n], Binomial[n, n/2] If[n == 0, 1, Tanh[r]^n]/(2^n Cosh[r]), 0]
```

Verify the zero-squeezing vacuum probability:

```wl
singleSqueezedProbability[0, 0] == 1
```

Compute the single-mode probabilities for zero through four photons:

```wl
singleSqueezedValues = singleSqueezedProbability[#, r] & /@ Range[0, 4];
```

Independent modes have joint probabilities given by an outer product. Form that reference matrix:

```wl
independentSqueezedProbabilities = Outer[Times, singleSqueezedValues, singleSqueezedValues];
```

Compute the joint probabilities after separation, with rows for mode 1 and columns for mode 2:

```wl
separatedSqueezedProbabilities = Table[
 separatedSqueezers["Probability", {0, 0} -> {j, k}], {j, 0, 4}, {k, 0, 4}];
```

Verify factorization for arbitrary positive `r`:

```wl
FullSimplify[separatedSqueezedProbabilities == independentSqueezedProbabilities, r > 0]
```

The total photon-number cutoff differs from the number of Fock levels per mode.
Return to `tmsv` at `r = 0.5`: cutoff $C$ retains pairs through $\lfloor C/2\rfloor$.
Compute the retained probabilities at three cutoffs:

```wl
pairRetainedProbabilities = AssociationMap[
 Total[tmsv["Probabilities", {0, 0}, "Cutoff" -> #]] &, {4, 8, 12}]
```

Compute the omitted probabilities:

```wl
Map[1 - # &, pairRetainedProbabilities]
```

Verify the retained probabilities against the geometric tail formula:

```wl
KeyValueMap[
 Abs[#2 - (1 - Tanh[0.5]^(2 (Floor[#1/2] + 1)))] < 10^-12 &,
 pairRetainedProbabilities]
```

The retained probability approaches one without renormalizing the listed outcomes.
QI computes squeezed-network probabilities directly, but `"Samples"` currently supports
only passive networks. Gaussian boson sampling here refers to the probability model.

*Source:* Agarwal, Section 5.8, p. 109, and Sections 3.1-3.2; Gerry and Knight, Section 7.7.

### 12. Optical attenuation: tracing out an unobserved output

A beam splitter with vacuum in its second input models attenuation when the second output
is unobserved. The joint evolution conserves photon number; tracing out the loss mode gives
the attenuated state. Build a beam splitter with symbolic angle:

```wl
attenuator = QuantumInterferometer[{"BS"[lossAngle, 0] -> {1, 2}}, 2];
```

Draw the beam splitter with its observed and loss outputs:

```wl
attenuator["Diagram", "ShowParameters" -> True]
```

Send three photons into the first input and vacuum into the second:

```wl
jointLossState = attenuator[FockState[{3, 0}, 4]];
```

Trace out the second output to obtain the attenuated state:

```wl
attenuatedState = QuantumPartialTrace[jointLossState, {2}];
```

Each photon is transmitted with probability $T=\cos^2\theta$.
Compute the binomial probabilities for zero through three transmitted photons:

```wl
lossProbabilities = PDF[BinomialDistribution[3, Cos[lossAngle]^2], Range[0, 3]]
```

Verify that the reduced density matrix is diagonal with these probabilities:

```wl
FullSimplify[Normal[attenuatedState["DensityMatrix"]] ==
 DiagonalMatrix[lossProbabilities], lossAngle ∈ Reals]
```

The joint state is pure, but the reduced state is generally mixed.
Verify normalization of its photon-number distribution:

```wl
FullSimplify[Total[lossProbabilities] == 1, lossAngle ∈ Reals]
```

Verify that the mean transmitted photon number is `3 T`:

```wl
FullSimplify[Range[0, 3] . lossProbabilities == 3 Cos[lossAngle]^2, lossAngle ∈ Reals]
```

Check perfect transmission:

```wl
(lossProbabilities /. lossAngle -> 0) == {0, 0, 0, 1}
```

Check perfect reflection:

```wl
(lossProbabilities /. lossAngle -> Pi/2) == {1, 0, 0, 0}
```

Check the reduced purity at equal transmissivity and reflectivity:

```wl
FullSimplify[attenuatedState["Purity"] /. lossAngle -> Pi/4] == 5/16
```

For arbitrary input photon number `nLoss`, let `kLoss` photons be transmitted.
The loss output contains `nLoss - kLoss` photons. Compute the amplitude on the physical domain:

```wl
lossAmplitudeFamily = Refine[
 attenuator["Amplitude", {nLoss, 0} -> {kLoss, nLoss - kLoss}],
 Element[{nLoss, kLoss}, Integers] && 0 <= kLoss <= nLoss]
```

The coefficient of $z_1^kz_2^{n-k}$ in $(\cos\theta\,z_1+\sin\theta\,z_2)^n$,
with Fock normalization, is $\sqrt{\binom nk}\cos^k\theta\sin^{n-k}\theta$.
Its squared modulus is binomial. Compute the probability on the physical integer domain:

```wl
lossCountLaw = FullSimplify[PDF[BinomialDistribution[nLoss, tLoss], kLoss],
 Element[{nLoss, kLoss}, Integers] && 0 <= kLoss <= nLoss && 0 < tLoss < 1]
```

Substitute explicit occupations into the QI amplitude while keeping the angle symbolic.
Verify the binomial law for zero through four input photons and $0<\theta<\pi/2$:

```wl
Table[FullSimplify[
 Abs[lossAmplitudeFamily /. {nLoss -> n, kLoss -> k}]^2 ==
  (lossCountLaw /. {nLoss -> n, kLoss -> k, tLoss -> Cos[lossAngle]^2}),
 0 < lossAngle < Pi/2], {n, 0, 4}, {k, 0, n}]
```

Verify normalization for every nonnegative integer input photon number:

```wl
Sum[lossCountLaw, {kLoss, 0, nLoss},
 Assumptions -> Element[nLoss, Integers] && nLoss >= 0 && 0 < tLoss < 1] == 1
```

Compute the mean transmitted photon number from the binomial law:

```wl
transmittedMean = Sum[kLoss lossCountLaw, {kLoss, 0, nLoss},
 Assumptions -> Element[nLoss, Integers] && nLoss >= 0 && 0 < tLoss < 1];
```

Verify that the mean is the input photon number times the transmissivity:

```wl
FullSimplify[transmittedMean == nLoss tLoss,
 Element[nLoss, Integers] && nLoss >= 0 && 0 < tLoss < 1]
```

Two attenuators with transmissivity $1/2$ have net transmissivity $1/4$.
Use the finite `QuantumState` route to propagate a mixed state through the second stage.
First, propagate three photons and a vacuum loss input through the balanced beam splitter:

```wl
firstLossJoint = balanced[FockState[{3, 0}, 4]];
```

Trace out the first loss output:

```wl
halfLossState = QuantumPartialTrace[firstLossJoint, {2}];
```

Add a fresh vacuum loss mode for the second stage:

```wl
secondLossInput = QuantumTensorProduct[halfLossState, FockState[0, 4]];
```

Propagate the mixed input through the second beam splitter:

```wl
secondLossJoint = balanced[secondLossInput];
```

Trace out the second loss output:

```wl
twiceAttenuated = QuantumPartialTrace[secondLossJoint, {2}];
```

Compute the binomial reference with transmissivity `1/4`:

```wl
twiceProbabilities = PDF[BinomialDistribution[3, 1/4], Range[0, 3]]
```

Verify the density matrix after both stages:

```wl
FullSimplify[Normal[twiceAttenuated["DensityMatrix"]] == DiagonalMatrix[twiceProbabilities]]
```

Four levels per mode suffice because the total photon number is at most three.
The attenuation comes from discarding a mode after unitary evolution.

*Source:* Agarwal, Section 5.10, pp. 111-112.

### 13. Waveguide arrays: transport, correlations, and quantum walks

Coupled waveguides exchange light. Propagation through an array is a passive mode
transformation, so QI can represent it as an equivalent optical network.

#### Coupled waveguides: continuous photon transport

Set the nearest-neighbor coupling to one and measure propagation length in inverse coupling
units. The propagator is $W(z)=e^{-izH}$, where $H$ is the single-particle coupling Hamiltonian.
Define $H$ for two waveguides:

```wl
twoGuideH = {{0, 1}, {1, 0}};
```

Compute the propagator with
[`MatrixExp`](https://reference.wolfram.com/language/ref/MatrixExp.html):

```wl
twoGuidePropagator = MatrixExp[-I couplingLength twoGuideH]
```

Build the equivalent beam splitter with reflection phase `-Pi/2`:

```wl
twoGuideSplitter = QuantumInterferometer[
 {"BS"[couplingLength, -Pi/2] -> {1, 2}}, 2];
```

Draw the equivalent beam splitter:

```wl
twoGuideSplitter["Diagram", "ShowParameters" -> True]
```

Verify agreement between the Hamiltonian propagator and beam splitter transfer matrix:

```wl
FullSimplify[twoGuidePropagator == twoGuideSplitter["Unitary"], couplingLength ∈ Reals]
```

The negative reflection phase matches the sign in $e^{-izH}$.
For a larger array, compute the propagator numerically before constructing QI's optical mesh.
Symbolic parameters remain available in element-based constructions.

Construct the nearest-neighbor Hamiltonian for five waveguides:

```wl
arrayHamiltonian = SparseArray[
 {Band[{1, 2}] -> 1, Band[{2, 1}] -> 1}, {5, 5}];
```

Verify that the array Hamiltonian is Hermitian:

```wl
HermitianMatrixQ[arrayHamiltonian]
```

Define the equivalent optical network at numerical propagation length `length`:

```wl
arrayDevice[length_?NumericQ] :=
 QuantumInterferometer[N[MatrixExp[-I length arrayHamiltonian]]]
```

Construct the network at unit propagation length:

```wl
arrayAtOne = arrayDevice[1];
```

Draw its equivalent beam splitter mesh:

```wl
arrayAtOne["Diagram"]
```

Read the transfer matrix:

```wl
arrayMatrix = arrayAtOne["Unitary"];
```

Verify that the array transfer matrix is unitary:

```wl
UnitaryMatrixQ[arrayMatrix]
```

Construct the zero-length network:

```wl
arrayAtZero = arrayDevice[0];
```

Draw its uncoupled modes:

```wl
arrayAtZero["Diagram"]
```

Verify the identity transfer matrix at zero length:

```wl
Max[Abs[arrayAtZero["Unitary"] - IdentityMatrix[5]]] < 10^-12
```

Compute the single-photon output probabilities for the central input:

```wl
arraySingleProbabilities = Table[
 arrayAtOne["Probability", {0, 0, 1, 0, 0} -> UnitVector[5, j]], {j, 5}]
```

Verify agreement with the squared moduli of the third transfer-matrix column:

```wl
Max[Abs[arraySingleProbabilities - Abs[arrayMatrix[[All, 3]]]^2]] < 10^-12
```

Reflection symmetry requires a symmetric distribution. Verify it:

```wl
Max[Abs[arraySingleProbabilities - Reverse[arraySingleProbabilities]]] < 10^-12
```

Verify normalization of the single-photon output distribution:

```wl
Abs[Total[arraySingleProbabilities] - 1] < 10^-12
```

A coherent input reproduces the same spatial profile in its mean output photon numbers.
Send symbolic coherent amplitude `arrayAlpha` into the central waveguide:

```wl
arrayClassicalField = arrayAtOne["Coherent"[arrayAlpha UnitVector[5, 3]]]
```

The mean output photon numbers equal the single-photon probabilities times $|\alpha|^2$.
Verify this relation, using `Chop` for numerical residuals:

```wl
Chop[FullSimplify[Abs[First[arrayClassicalField]]^2 -
 Abs[arrayAlpha]^2 arraySingleProbabilities]] == ConstantArray[0, 5]
```

Mean photon numbers do not describe two-photon correlations. The normally ordered
correlation $G_{jk}=\langle n_j(n_k-\delta_{jk})\rangle$ counts ordered detection pairs.
Compute the output distribution for photons entering waveguides two and four:

```wl
arrayPairProbabilities = arrayAtOne["Probabilities", {0, 1, 0, 1, 0}];
```

Compute the correlation matrix by averaging over occupation patterns:

```wl
arrayCorrelations = Total[KeyValueMap[
 With[{n = #1["Name"]}, #2 (Outer[Times, n, n] - DiagonalMatrix[n])] &,
 arrayPairProbabilities]];
```

Two photons give two ordered pairs. Verify that the matrix entries sum to two:

```wl
Abs[Total[arrayCorrelations, 2] - 2] < 10^-12
```

Plot the two-photon correlation matrix:

```wl
ArrayPlot[arrayCorrelations, Frame -> True, FrameLabel -> {"Output mode", "Output mode"},
 PlotLabel -> "Two photons propagating through five coupled waveguides",
 PlotLegends -> Automatic, ImageSize -> 380]
```

QI can sample without first enumerating the full distribution. Here we compare sample
frequencies with the distribution already calculated. Set a reproducible seed:

```wl
SeedRandom[23];
```

Draw three thousand samples from the two-photon input:

```wl
arraySampleCounts = arrayAtOne["SampleCounts", {0, 1, 0, 1, 0}, 3000];
```

Verify the number of samples:

```wl
Total[arraySampleCounts] == 3000
```

Convert the sample counts to relative frequencies:

```wl
arrayFrequencies = arraySampleCounts/3000;
```

Subtract the frequencies from the predicted probabilities by occupation pattern:

```wl
arrayFrequencyDifferences = Merge[{arrayPairProbabilities, -arrayFrequencies}, Total];
```

Compute the largest absolute frequency difference:

```wl
Max[Abs[Values[arrayFrequencyDifferences]]]
```

This error fluctuates between runs. Increase the sample size to explore its dependence
on the number of samples.

#### Discrete-time quantum walks: coupling position and direction

A discrete-time quantum walk alternates a **coin operation**, which mixes two internal
states, with a shift conditioned on the coin state. Here two optical modes at each position
encode the coin states: one moves left and the other right.

Use nine positions with periodic boundaries. From the center, four steps cannot wrap around,
so the result also describes four steps on an infinite line. Define the positions:

```wl
walkPositions = Range[-4, 4];
```

Define the Hadamard coin operation:

```wl
walkCoin = {{1, 1}, {1, -1}}/Sqrt[2];
```

Construct the conditional shift with periodic boundaries:

```wl
walkShift = KroneckerProduct[RotateLeft[IdentityMatrix[9]], {{1, 0}, {0, 0}}] +
 KroneckerProduct[RotateRight[IdentityMatrix[9]], {{0, 0}, {0, 1}}];
```

Apply the coin operation before the shift to form one step:

```wl
walkStep = walkShift . KroneckerProduct[IdentityMatrix[9], walkCoin];
```

Verify that one step is unitary:

```wl
UnitaryMatrixQ[walkStep]
```

Modes 9 and 10 encode the two coin states at position zero. Prepare a relative phase of $\pi/2$
to obtain a symmetric position distribution. Build the preparation beam splitter:

```wl
walkPreparation = QuantumInterferometer[
 {"BS"[Pi/4, Pi/2] -> {9, 10}}, 18];
```

Draw the preparation network:

```wl
walkPreparation["Diagram", "ShowParameters" -> True]
```

Compute the transfer matrix for preparation followed by four walk steps:

```wl
walkTransfer = MatrixPower[walkStep, 4] . walkPreparation["Unitary"];
```

Construct the equivalent optical network from this matrix:

```wl
walkDevice = QuantumInterferometer[N[walkTransfer]];
```

Draw the eighteen-mode network:

```wl
walkDevice["Diagram", ImageSize -> 1000]
```

Compute the output distribution for a photon entering mode 9:

```wl
walkProbabilities = walkDevice["Probabilities", UnitVector[18, 9]];
```

Sum over the two coin states at each position to obtain the position distribution:

```wl
walkPositionProbabilities = Total[KeyValueMap[
 #2 (Total /@ Partition[#1["Name"], 2]) &, walkProbabilities]]
```

Verify normalization after summing over the coin states:

```wl
Abs[Total[walkPositionProbabilities] - 1] < 10^-12
```

Verify reflection symmetry of the walk's position distribution:

```wl
Max[Abs[walkPositionProbabilities - Reverse[walkPositionProbabilities]]] < 10^-12
```

Compute the position distribution for a four-step classical unbiased random walk:

```wl
classicalPositionProbabilities =
 PDF[BinomialDistribution[4, 1/2], (4 + walkPositions)/2]
```

The classical walk adds path probabilities; the quantum walk adds path amplitudes.
Pair each position with its probability for both walks:

```wl
walkPlotData = Transpose[{walkPositions, #}] & /@
 {walkPositionProbabilities, classicalPositionProbabilities};
```

Compare the two position distributions:

```wl
ListLinePlot[walkPlotData, Frame -> True, FrameLabel -> {"Position", "Probability"},
 PlotLegends -> {"Quantum walk", "Classical walk"},
 PlotLabel -> "Four-step position distributions", PlotMarkers -> Automatic,
 GridLines -> Automatic, PlotRange -> All, AspectRatio -> 1/2, ImageSize -> 520]
```

Summing over the coin states preserves the interference accumulated during propagation.
Increase the position range when adding steps to avoid the periodic boundary.

*Source:* Agarwal, Sections 19.1-19.3 and 19.5, pp. 432-442 and 447-452.

### 14. Conditional optical gates: implementing a nonlinear sign change

A phase shifter multiplies an $n$-photon Fock state by $e^{in\phi}$.
Build a one-mode phase shifter with a symbolic phase:

```wl
ordinaryPhase = QuantumInterferometer[{"PS"[phaseNS] -> 1}, 1];
```

Draw the phase shifter:

```wl
ordinaryPhase["Diagram", "ShowParameters" -> True]
```

Compute its amplitude, assuming a nonnegative integer photon number:

```wl
phaseAmplitude = Refine[
 ordinaryPhase["Amplitude", {phasePhotons} -> {phasePhotons}],
 Element[phasePhotons, Integers] && phasePhotons >= 0]
```

Verify the photon-number phase law:

```wl
FullSimplify[phaseAmplitude == Exp[I phasePhotons phaseNS],
 Element[phasePhotons, Integers] && phasePhotons >= 0 && phaseNS ∈ Reals]
```

At phase `Pi/2`, read the vacuum, one-photon, and two-photon amplitudes:

```wl
Table[phaseAmplitude /. {phaseNS -> Pi/2, phasePhotons -> n}, {n, 0, 2}]
```

The two-photon amplitude changes sign, but the one-photon amplitude also changes.
A **nonlinear sign gate** instead maps
$c_0|0\rangle+c_1|1\rangle+c_2|2\rangle$ to
$c_0|0\rangle+c_1|1\rangle-c_2|2\rangle$.

Linear optics can implement this gate conditionally with auxiliary photons and heralding.
Mode 1 carries the signal; modes 2 and 3 contain one photon and vacuum. Success requires
the same auxiliary output occupation pattern `{1, 0}`. Construct the three-mode transfer matrix:

```wl
nsMatrix = {
 {1 - Sqrt[2], 2^(-1/4), Sqrt[3/Sqrt[2] - 2]},
 {2^(-1/4), 1/2, (1 - Sqrt[2])/2},
 {Sqrt[3/Sqrt[2] - 2], (1 - Sqrt[2])/2, Sqrt[2] - 1/2}};
```

The matrix is real. Verify unitarity using its transpose:

```wl
RootReduce[nsMatrix . Transpose[nsMatrix]] == IdentityMatrix[3]
```

The selected transition amplitude is $\langle n,1,0|U|n,1,0\rangle$.
Photon-number conservation preserves the signal photon number on success.
Build the network from the transfer matrix:

```wl
nsGate = QuantumInterferometer[nsMatrix];
```

Draw the network before heralding:

```wl
nsGate["Diagram"]
```

Compute the selected transition amplitudes for signal photon numbers zero, one, and two:

```wl
nsAmplitudes = RootReduce[Table[
 nsGate["Amplitude", {n, 1, 0} -> {n, 1, 0}], {n, 0, 2}]]
```

Verify the amplitude vector for the nonlinear sign gate:

```wl
nsAmplitudes == {1, 1, -1}/2
```

The common amplitude $1/2$ gives success probability $1/4$ for every normalized input
in the gate subspace. Verify this for arbitrary complex coefficients:

```wl
FullSimplify[Total[Abs[nsAmplitudes {c0, c1, c2}]^2] == 1/4,
 Total[Abs[{c0, c1, c2}]^2] == 1]
```

Check the full state evolution and heralding projection. Four levels per mode retain all
sectors because the total input photon number is at most three. Prepare the signal with
$|c_0|^2+|c_1|^2+|c_2|^2=1$, declaring its coefficients for later substitution:

```wl
nsSignal = QuantumState[{c0, c1, c2, 0}, 4, "Parameters" -> {c0, c1, c2}];
```

Add the one-photon and vacuum auxiliary inputs:

```wl
nsInput = QuantumTensorProduct[nsSignal, FockState[1, 4], FockState[0, 4]];
```

Propagate the full input state:

```wl
nsOutput = nsGate[nsInput];
```

Verify that propagation preserves the parameter declarations:

```wl
nsOutput["ParameterSpec"] === nsInput["ParameterSpec"]
```

Reshape the output amplitudes so the three indices label the three modes:

```wl
nsOutputTensor = ArrayReshape[Normal[nsOutput["StateVector"]], {4, 4, 4}];
```

Auxiliary photon numbers one and zero select indices two and one.
Extract the unnormalized conditional signal vector:

```wl
nsSelected = nsOutputTensor[[All, 2, 1]];
```

Verify the conditional vector before normalization:

```wl
FullSimplify[nsSelected == {c0, c1, -c2, 0}/2]
```

Compute the success probability:

```wl
nsSuccess = Total[Abs[nsSelected]^2];
```

Verify success probability `1/4` for every normalized input in the gate subspace:

```wl
FullSimplify[nsSuccess == 1/4, Total[Abs[{c0, c1, c2}]^2] == 1]
```

Normalize the selected vector and carry its parameter declarations into the conditional state:

```wl
nsConditionalState = QuantumState[nsSelected/Sqrt[nsSuccess], 4,
 "ParameterSpec" -> nsOutput["ParameterSpec"]];
```

Read the conditional state vector under the input normalization constraint:

```wl
FullSimplify[Normal[nsConditionalState["StateVector"]], Total[Abs[{c0, c1, c2}]^2] == 1]
```

Verify normalization of the conditional state:

```wl
FullSimplify[nsConditionalState["Norm"] == 1, Total[Abs[{c0, c1, c2}]^2] == 1]
```

Bind a normalized complex input through the conditional state's parameter interface:

```wl
nsExampleState = nsConditionalState[
 <|c0 -> 1/Sqrt[3], c1 -> I/Sqrt[3], c2 -> -1/Sqrt[3]|>];
```

Verify that the two-photon coefficient changes sign:

```wl
FullSimplify[Normal[nsExampleState["StateVector"]]] == {1, I, 1, 0}/Sqrt[3]
```

The output coefficients are proportional to `{1, I, 1}`, as required.
The gate is defined only on signal photon numbers zero through two.
Compute the selected transition amplitude for a three-photon signal:

```wl
nsBeyondAmplitude = RootReduce[nsGate["Amplitude", {3, 1, 0} -> {3, 1, 0}]]
```

Verify that it differs from both gate amplitudes:

```wl
nsBeyondAmplitude != 1/2 && nsBeyondAmplitude != -1/2
```

*Source:* Kok and Lovett, *Introduction to Optical Quantum Information Processing*,
Section 6.2.2, Exercise 6.4 and Figure 6.5, pp. 183-186. QI computes the unitary evolution;
projection onto the heralding outcome gives the conditional gate.

### Consistency checks: normalization, limits, and independent calculations

When changing an example, revisit the relevant checks:

| Regime | Examples |
|---|---|
| General symbolic | Optical phases, coherent amplitudes, twin Fock inputs, squeezed-state factorization, binomial attenuation, gate superpositions |
| Exact special case | Hong-Ou-Mandel suppression, geometric pair-number distribution, nonlinear sign gate |
| Limit | Vacuum, zero squeezing, zero coherent amplitude, perfect transmission or reflection, zero propagation |
| Numerical comparison | Hafnian and permanent amplitudes versus Fock-space circuits; sample frequencies versus probabilities |
| Boundary | Omitted probability, distinguishable photons, factorial moments, periodic walk boundaries, gate inputs above two photons |

### Scope and further experiments

QI supports symbolic coherent amplitudes in passive networks, symbolic passive Fock-state
transition queries, and finite `QuantumState` superpositions and mixtures. Passive propagation
preserves declared state parameters for later substitution. Squeezed networks
allow symbolic optical parameters with explicit occupations. A Gaussian state specified only
by its mean vector and covariance matrix is not a native input; its Fock-space representation
can use the finite-state route.

Homodyne detection, temporal wave-packet overlap, detector imperfections, and adaptive
feed-forward require additional modeling. Change one input, optical parameter, or heralding
outcome at a time, then rerun the normalization and limiting-case checks.
