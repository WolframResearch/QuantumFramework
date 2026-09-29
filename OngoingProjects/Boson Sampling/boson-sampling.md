---
Template: Default
Title: Boson Sampling Explained
Author: Mads Bahrami
---

# Boson Sampling Explained

<!-- #| style: Subtitle -->
Hong-Ou-Mandel Interference, Multiport Interferometers, Permanents, and the Classical Cost of Simulating Identical Photons: A Guided, Computation-First Narrative with the Wolfram Language and QuantumFramework

<!-- #| style: Author -->
Mads Bahrami (last updated: Sep 27, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

### Setting the Stage: How This Notebook Flows

This notebook is a computation-first tour of boson sampling: how two identical photons meeting at a beam splitter leave together (the Hong-Ou-Mandel effect), how the photons' distinguishability switches that effect off, how any passive optical network is assembled from beam splitters and phase shifters, how the amplitude for a pattern of photon detections becomes the permanent of a matrix, why some patterns are exactly forbidden, and why the cost of simulating identical photons on a classical computer grows so quickly. Along the way we simulate the photons directly in their Fock space with the second-quantization tools of [QuantumFramework](https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/) (Fock states, beam splitters, phase shifters, quantum circuits, and a symbolic algebra of creation and annihilation operators), and we cross-check every result against the permanent formula and against two computations written in plain Wolfram Language.

In other words, I've tried to build a catalogue of computational experiments on identical photons, tools and examples that let you run the experiments yourself and learn directly from what you observe. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it.

Before we start, pay attention to a few things. The environment you see is a [Wolfram (Mathematica) notebook](https://www.wolfram.com/notebooks/). You should evaluate the cells from top to bottom. Although I did my best to make the input cells independent, several functions defined early (for example the one that builds an interferometer) are used until the very end, so some cells depend on previous ones. A few cells in Part IV time computations on purpose; they take up to about a minute each, and the absolute times depend on your computer, so only their trends matter.

Additionally, the story is presented as a continuous sequence, like a movie. I have added a few headings to help with transitions from one topic to another, but I have avoided breaking the narrative into rigid sections. Sometimes a feature is introduced and used before we explain why it applies. This is intentional, because the ability to apply an idea often matters more than following an abstract proof first. Overall, the rhythm is: concept, then computation, then interpretation.

Inevitably, some of the code will look complicated at first. My suggestion is to focus on the output and its meaning before worrying about every detail of the input code. Once you understand what a piece of code is doing conceptually, go back and unpack how it is written.

Remember that you are not locked into the code as given. You can (and should) modify it, try your own variations (other angles, other interferometers, more photons), and run your own numerical experiments.

If you have any suggestions or questions, please reach out to us at [quantum@wolfram.com](mailto:quantum@wolfram.com)

### Prerequisites and How to Read This

To follow the computations comfortably you should be at ease with Dirac notation, tensor products, and the creation and annihilation operators of a harmonic oscillator, $\hat a^\dagger|n\rangle=\sqrt{n+1}\,|n+1\rangle$ and $\hat a|n\rangle=\sqrt{n}\,|n-1\rangle$. No background in complexity theory is assumed; the few facts we need are stated where they appear. The notebook uses QuantumFramework 2.1 or later, whose second-quantization sub-package provides the Fock states, the optical elements, and a symbolic algebra of ladder operators that needs version 15 of the Wolfram Language.

Parts I through III build the physics, from two photons to a random interferometer. Part IV asks how the cost of simulating the photons grows, and it is the part where cells take longer to evaluate.

Let's start!

## Part I: Two Photons and One Beam Splitter

### Setting Up: QuantumFramework and Its Second-Quantization Tools

QuantumFramework keeps its bosonic tools (Fock states, annihilation operators, displacement, squeezing, beam splitters, phase shifters) in a separate sub-package. Loading the main package alone leaves those names undefined, so we load both.

Load QuantumFramework together with its second-quantization sub-package:

```wl
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

### Fock States: Photons as Occupation Numbers

For light, a mode is one way a photon can travel, for example one arm of an interferometer. A state of several photons in several modes is specified by how many photons occupy each mode, $|n_1,n_2,\ldots\rangle$, which is called a Fock (or occupation-number) state. A computer cannot store the infinitely many levels $n=0,1,2,\ldots$ of a mode, so each mode is truncated to a finite number of levels $d$, holding between zero and $d-1$ photons. In other words, a Fock state lists how many photons sit in each mode, and the truncation fixes the largest number of photons a single mode can hold.

Create the two-mode Fock state $|1,1\rangle$ (one photon in each mode), keeping three levels (zero, one, or two photons) in each mode:

```wl
in = FockState[{1, 1}, 3];
in["Formula"]
```

Show the dimensions of this truncated two-mode Fock space:

```wl
in["Dimensions"]
```

So the space has $3\times 3=9$ basis states $|n_1,n_2\rangle$, but only three of them, $|2,0\rangle$, $|1,1\rangle$ and $|0,2\rangle$, contain exactly two photons. We will see in a moment why the other six can never be populated in this experiment.

QuantumFramework orders these basis states by reading the occupations as the digits of a base-$d$ number, mode 1 being the most significant digit: the position of $|n_1,n_2,\ldots\rangle$ in the state vector is that number plus one.

Define a function that returns the position of an occupation pattern in the state vector:

```wl
ClearAll[fockIndex]
fockIndex[occ_List, d_Integer] := FromDigits[occ, d] + 1
```

Verify that the state vector of $|1,1\rangle$ has its only nonzero entry at that position:

```wl
Normal[in["StateVector"]] == UnitVector[9, fockIndex[{1, 1}, 3]]
```

### The Beam Splitter: A Two-Mode Unitary That Conserves Photons

A beam splitter is a partially reflecting mirror that mixes two modes. QuantumFramework's `BeamSplitterOperator[{θ, φ}, d, {1, 2}]` represents the unitary $\hat B(\theta,\phi)=\exp\left(\theta\,(e^{i\phi}\hat a_1\hat a_2^\dagger-e^{-i\phi}\hat a_1^\dagger\hat a_2)\right)$ acting on modes 1 and 2, each truncated to $d$ levels. The mixing angle θ controls how much light crosses from one mode to the other: for $\phi=0$, a photon that enters mode 1 stays with amplitude $\cos\theta$ and crosses into mode 2 with amplitude $\sin\theta$, while a photon that enters mode 2 crosses into mode 1 with amplitude $-\sin\theta$. So a single photon stays in its mode with probability $\cos^2\theta$, and $\theta=\pi/4$ is the balanced (50:50) beam splitter. The angle φ is a relative phase between the two paths. There are, of course, different conventions: many texts parametrize a beam splitter by its transmissivity $T=\cos^2\theta$, and the balanced splitter is often written with a factor $i$ on the reflected path, which corresponds to $\phi=\pi/2$ here. Nothing below depends on that choice. The option `Method -> "Recurrence"` builds each matrix element of the untruncated beam splitter from a recurrence and keeps the ones that fit below the cutoff; for symbolic angles and small cutoffs like ours, it is faster than exponentiating the truncated generator.

Define a beam splitter with a symbolic mixing angle θ (and $\phi=0$) acting on modes 1 and 2:

```wl
bs = BeamSplitterOperator[{\[Theta], 0}, 3, {1, 2}, Method -> "Recurrence"];
```

Look at the generator of $\hat B$: the term $\hat a_1\hat a_2^\dagger$ removes a photon from mode 1 and puts it into mode 2, and its partner does the reverse. So the generator moves photons between the modes without creating or destroying any, which means $\hat B$ commutes with the total photon number $\hat n_1+\hat n_2=\hat a_1^\dagger\hat a_1+\hat a_2^\dagger\hat a_2$. In other words, a beam splitter redirects photons; it never adds or removes them.

We can check this without any truncation. QuantumFramework's symbolic bosonic algebra works with the ladder operators themselves, as noncommuting formal symbols obeying $[\hat a,\hat a^\dagger]=1$. `FieldVariables` creates them, and `BosonicNormalOrder` rewrites a product of them with every creation operator to the left of every annihilation operator, treating the symbols listed in its `"Scalars"` option as ordinary numbers. The built-in [Commutator](https://reference.wolfram.com/language/ref/Commutator.html)`[x, y]` gives the commutator `x ** y - y ** x` of two such products. QuantumFramework displays the annihilation operators of our two modes as the formal symbols a and b.

Create the ladder operators of the two modes, the generator $G=\theta\,(e^{i\phi}\hat a_1\hat a_2^\dagger-e^{-i\phi}\hat a_1^\dagger\hat a_2)$ of the beam splitter, and the total photon number:

```wl
fieldVars = Join[FieldVariables[\[FormalA]], FieldVariables[\[FormalB]]];
{a1, a1Dag, a2, a2Dag} = fieldVars;
generator = \[Theta] (Exp[I \[Phi]] a1 ** a2Dag - Exp[-I \[Phi]] a1Dag ** a2);
totalNumber = a1Dag ** a1 + a2Dag ** a2;
```

Compute the commutator of the generator with the total photon number:

```wl
BosonicNormalOrder[Commutator[generator, totalNumber], fieldVars, "Scalars" -> {\[Theta], \[Phi]}]
```

The truncated operator we built inherits this property.

Verify that the truncated beam splitter commutes with the truncated total photon number:

```wl
With[{a1 = AnnihilationOperator[3, {1}], a2 = AnnihilationOperator[3, {2}]},
 With[{ntot = a1["Dagger"] @ a1 + a2["Dagger"] @ a2},
  FullSimplify[Normal[(bs @ ntot - ntot @ bs)["Matrix"]]] == ConstantArray[0, {9, 9}]]]
```

This is why only the three two-photon states can appear after the beam splitter: the input $|1,1\rangle$ carries two photons, and the photon number is conserved.

Photon-number conservation also means that the beam splitter acts separately on the states with zero, one, two, and more photons in total. With three levels per mode, the sectors with at most two photons are complete, while a state like $|1,2\rangle$ belongs to a three-photon sector that the truncation cuts short: it has no room for $|3,0\rangle$. The `"Recurrence"` construction is unitary on the complete sectors, those with at most $d-1$ photons, and those are the only ones our photons will visit.

Verify that the beam splitter is unitary on the six states with at most two photons:

```wl
With[{mat = Normal[bs["Matrix"]], sel = fockIndex[#, 3] & /@ Select[Tuples[Range[0, 2], 2], Total[#] <= 2 &]},
 FullSimplify[mat[[sel, sel]] . ConjugateTranspose[mat[[sel, sel]]] == IdentityMatrix[Length[sel]], \[Theta] \[Element] Reals]]
```

### Hong-Ou-Mandel Interference: Why Two Photons Leave Together

Now let's send the two photons through the beam splitter, one in each input port.

Send $|1,1\rangle$ through the beam splitter and show the output state:

```wl
psi = FullSimplify[bs[in], \[Theta] \[Element] Reals];
psi["Formula"]
```

As one can see, the photons leave in a superposition of three patterns. The amplitude of $|1,1\rangle$, one photon in each output port, is $\cos 2\theta=\cos^2\theta-\sin^2\theta$. It is the sum of two histories that end in the same detection event: both photons are transmitted (amplitude $\cos\theta\cdot\cos\theta$), or both photons are reflected (amplitude $-\sin\theta\cdot\sin\theta$). Because the photons are identical, nothing in the final state records which history happened, so the two amplitudes are added before squaring.

Compute the norm of the output state:

```wl
FullSimplify[psi["Norm"], \[Theta] \[Element] Reals]
```

QuantumFramework's `"Probability"` divides the squared amplitudes by their total, so a lost norm never shows up in it; that is why we check the norm itself.

Compute the probability of each detection pattern:

```wl
FullSimplify[#, \[Theta] \[Element] Reals] & /@ psi["Probability"]
```

Verify that these probabilities do not depend on the phase φ of the beam splitter:

```wl
FullSimplify[Values[BeamSplitterOperator[{\[Theta], \[Phi]}, 3, {1, 2}, Method -> "Recurrence"][in]["Probability"]] == Values[psi["Probability"]], (\[Theta] | \[Phi]) \[Element] Reals]
```

The pattern $|1,1\rangle$, one photon in each output, is what a pair of detectors records as a coincidence, and its probability depends on θ. Let's find where it vanishes.

Find the mixing angles between 0 and $\pi/2$ at which the coincidence amplitude $\cos 2\theta$ is zero:

```wl
Solve[Cos[2 \[Theta]] == 0 && 0 < \[Theta] < \[Pi]/2, \[Theta]]
```

Therefore, at a balanced beam splitter the two histories cancel exactly, and two identical photons that enter through different ports always leave through the same port, together. This is the Hong-Ou-Mandel effect, named after [Hong, Ou, and Mandel](https://doi.org/10.1103/PhysRevLett.59.2044), whose 1987 experiment used this fourth-order (two-photon) interference to measure the time interval between two photons. In short, identical photons bunch.

### Distinguishable Photons: Switching the Interference Off

The cancellation above needs the two histories to be indistinguishable. What happens if the photons carry a label that tells them apart? A natural label is polarization. Give every path two polarization modes, horizontal (H) and vertical (V), and order the four modes as $1_H, 2_H, 1_V, 2_V$. The photon in path 1 is horizontally polarized, and the photon in path 2 is polarized at an angle χ, $\cos\chi\,|H\rangle+\sin\chi\,|V\rangle$. At $\chi=0$ the photons are identical; at $\chi=\pi/2$ their polarizations are orthogonal, and a measurement of polarization would reveal which photon went where. Both photons may end up in the same horizontal mode, so each mode again needs three levels.

Prepare one horizontal photon in path 1 and one photon polarized at angle χ in path 2 (modes ordered $1_H, 2_H, 1_V, 2_V$):

```wl
inPol = Cos[\[Chi]] FockState[{1, 1, 0, 0}, 3] + Sin[\[Chi]] FockState[{1, 0, 0, 1}, 3];
```

A beam splitter that does not care about polarization, called non-polarizing, mixes the two horizontal modes and, identically, the two vertical modes.

Build the non-polarizing beam splitter as a circuit, one copy acting on the horizontal modes and one on the vertical modes:

```wl
bsPol = QuantumCircuitOperator[{
    BeamSplitterOperator[{\[Theta], 0}, 3, {1, 2}, Method -> "Recurrence"],
    BeamSplitterOperator[{\[Theta], 0}, 3, {3, 4}, Method -> "Recurrence"]}];
```

Send the two photons through the non-polarizing beam splitter:

```wl
outPol = bsPol[inPol];
```

Compute the norm of the output state:

```wl
FullSimplify[outPol["Norm"], (\[Theta] | \[Chi]) \[Element] Reals]
```

The detectors count photons in each path but ignore polarization. A coincidence is therefore any pattern with exactly one photon in path 1 ($n_{1_H}+n_{1_V}=1$) and exactly one photon in path 2 ($n_{2_H}+n_{2_V}=1$).

Compute the coincidence probability, adding up the probabilities of all patterns with one photon in each path ([ComplexExpand](https://reference.wolfram.com/language/ref/ComplexExpand.html) treats the angles θ and χ as real):

```wl
coincidence = FullSimplify @ ComplexExpand @ Total @ KeySelect[KeyMap[First, outPol["Probability"]],
    #[[1]] + #[[3]] == 1 && #[[2]] + #[[4]] == 1 &]
```

Verify that the coincidence probability equals $\cos^4\theta+\sin^4\theta-\tfrac{1}{2}\sin^2 2\theta\cos^2\chi$:

```wl
FullSimplify[coincidence == Cos[\[Theta]]^4 + Sin[\[Theta]]^4 - Sin[2 \[Theta]]^2 Cos[\[Chi]]^2/2, (\[Theta] | \[Chi]) \[Element] Reals]
```

Look at the two limits. For identical photons ($\chi=0$) the formula collapses to $\cos^2 2\theta$, the Hong-Ou-Mandel result.

Verify that at $\chi=0$ the coincidence probability is $\cos^2 2\theta$:

```wl
FullSimplify[(coincidence /. \[Chi] -> 0) == Cos[2 \[Theta]]^2, \[Theta] \[Element] Reals]
```

For orthogonal polarizations ($\chi=\pi/2$) the formula becomes $\cos^4\theta+\sin^4\theta$, which is exactly what two independent classical particles give: both transmitted or both reflected, with probabilities multiplied rather than amplitudes added. In other words, how far the coincidence probability at the balanced splitter falls below the value for distinguishable photons measures how indistinguishable the photons are, through the squared overlap $\cos^2\chi$ of their polarization states.

Plot the coincidence probability against the mixing angle θ, for photons ranging from identical ($\chi=0$) to fully distinguishable ($\chi=\pi/2$):

```wl
Plot[Evaluate[coincidence /. \[Chi] -> {0, \[Pi]/6, \[Pi]/3, \[Pi]/2}], {\[Theta], 0, \[Pi]/2},
 PlotLegends -> {"\[Chi] = 0", "\[Chi] = \[Pi]/6", "\[Chi] = \[Pi]/3", "\[Chi] = \[Pi]/2"},
 Frame -> True, GridLines -> Automatic, AspectRatio -> 1/2, ImageSize -> Large,
 FrameLabel -> {"mixing angle \[Theta]", "coincidence probability"},
 PlotLabel -> "Coincidence probability for photons with polarization overlap Cos[\[Chi]]"]
```

As you can see, the balanced beam splitter ($\theta=\pi/4$) is where identical and distinguishable photons differ most: the classical curve never drops below one half, while the identical photons never arrive in coincidence.

Compute the coincidence probability at the balanced splitter as a function of χ:

```wl
FullSimplify[coincidence /. \[Theta] -> \[Pi]/4, \[Chi] \[Element] Reals]
```

As one can see, it falls from one half for orthogonal polarizations to zero for identical photons. This drop is the Hong-Ou-Mandel dip. In experiments the dip is usually recorded by delaying one photon, which makes the two photons distinguishable by their arrival times rather than by their polarizations.

### Too Few Levels: What Truncation Does to Bunching

Recall that we kept three levels in each mode, because the two photons can bunch into one mode. Let's see what happens if we are stingy and keep only two levels, enough for zero or one photon per mode.

Apply the same beam splitter with only two levels per mode, and compute the norm of the output state:

```wl
FullSimplify[BeamSplitterOperator[{\[Theta], 0}, 2, {1, 2}, Method -> "Recurrence"][FockState[{1, 1}, 2]]["Norm"], \[Theta] \[Element] Reals]
```

With two levels per mode the two-photon sector is no longer complete: the bunched components $|2,0\rangle$ and $|0,2\rangle$ have nowhere to go, so this construction simply drops them and the state loses norm; at the balanced splitter nothing at all is left. The default `"MatrixExp"` construction fails differently. It exponentiates the truncated generator, which keeps the operator unitary, but in the truncated space $\hat a_1\hat a_2^\dagger$ can no longer put a second photon into mode 2.

Repeat with the default `"MatrixExp"` construction, and show the output state:

```wl
FullSimplify[BeamSplitterOperator[{\[Theta], 0}, 2, {1, 2}][FockState[{1, 1}, 2]], \[Theta] \[Element] Reals]["Formula"]
```

Now the norm is preserved, but the physics is wrong: $|1,1\rangle$ is frozen, and the photons pass through untouched for every θ. In other words, a simulation of $n$ photons needs at least $n+1$ levels in every mode, because all $n$ photons may bunch into one of them. A truncation that is one level too small either leaks probability or silently erases the interference. This requirement is what makes the truncated Fock space grow so quickly as photons are added.

### The Permanent Appears: The Coincidence Amplitude of a 2×2 Matrix

Everything a beam splitter does to a single photon fits in a $2\times 2$ matrix $W$, the beam splitter's transfer matrix: its entry $W_{jk}=\langle 1_j|\hat B|1_k\rangle$ is the amplitude for a photon that enters mode $k$ to leave from mode $j$ (here $|1_k\rangle$ is one photon in mode $k$ and none elsewhere).

Define a function that reads the transfer matrix $W_{jk}=\langle 1_j|\hat U|1_k\rangle$ of an $m$-mode operator by sending one photon into each mode:

```wl
ClearAll[transferMatrix]
transferMatrix[op_, m_Integer, d_Integer] := Transpose @ Map[
   Normal[op[FockState[#, d]]["StateVector"][[fockIndex[#, d] & /@ IdentityMatrix[m]]]] &,
   IdentityMatrix[m]]
```

Compute the transfer matrix of our beam splitter:

```wl
w = FullSimplify[transferMatrix[bs, 2, 3], \[Theta] \[Element] Reals];
w // MatrixForm
```

As one can see, a beam splitter with $\phi=0$ acts on a single photon as an ordinary rotation matrix, whose columns hold the amplitudes $\cos\theta$, $\sin\theta$ and $-\sin\theta$ we started from. For two photons, the coincidence amplitude we found above adds the two histories $W_{11}W_{22}$ (both transmitted) and $W_{12}W_{21}$ (both reflected) with the same sign. This combination is the permanent of $W$. The permanent of an $n\times n$ matrix is defined like the determinant, as a sum over all permutations σ of products $\prod_i A_{i\,\sigma(i)}$, but with every term added with a plus sign: $\mathrm{Perm}(A)=\sum_{\sigma}\prod_{i=1}^{n}A_{i\,\sigma(i)}$. The Wolfram Language computes it with [Permanent](https://reference.wolfram.com/language/ref/Permanent.html).

Compute the permanent of $W$:

```wl
FullSimplify[Permanent[w], \[Theta] \[Element] Reals]
```

Verify that the permanent of $W$ is exactly the amplitude of $|1,1\rangle$ in the output state:

```wl
FullSimplify[Permanent[w] == psi["StateVector"][[fockIndex[{1, 1}, 3]]], \[Theta] \[Element] Reals]
```

Two identical fermions in the same situation pick up a minus sign when their histories are exchanged, so their coincidence amplitude is the determinant $W_{11}W_{22}-W_{12}W_{21}$ instead of the permanent.

Compute the determinant of $W$, the corresponding amplitude for two identical fermions:

```wl
FullSimplify[Det[w], \[Theta] \[Element] Reals]
```

This means that two identical fermions always leave through different ports, the exact opposite of the bunching bosons: the exchange sign turns the destructive interference into constructive interference. We will return to this contrast at the very end, because the difference between a permanent and a determinant is also the difference between a hard and an easy computation.

For distinguishable photons, the probabilities of the two histories add instead of the amplitudes, which amounts to the permanent of the matrix of single-photon probabilities $|W_{jk}|^2$.

Compute the permanent of $|W_{jk}|^2$, the coincidence probability for distinguishable photons:

```wl
FullSimplify[Permanent[Abs[w]^2], \[Theta] \[Element] Reals]
```

Verify that it equals the coincidence probability of the polarization model at $\chi=\pi/2$:

```wl
FullSimplify[Permanent[Abs[w]^2] == (coincidence /. \[Chi] -> \[Pi]/2), \[Theta] \[Element] Reals]
```

Before we move on, let us summarize the most important points we have learned so far:

- A Fock state $|n_1,n_2,\ldots\rangle$ lists the photons in each mode, and `FockState[{n1, n2, ...}, d]` keeps $d$ levels per mode.

- A beam splitter commutes with the total photon number $\hat n_1+\hat n_2$, in the untruncated algebra and in the truncated space, so it only redistributes photons, and it is unitary on every photon-number sector that the truncation holds completely.

- For $|1,1\rangle$ the coincidence amplitude is $\cos 2\theta$, which vanishes at the balanced splitter: identical photons leave together (the Hong-Ou-Mandel effect).

- With a polarization overlap $\cos\chi$, the coincidence probability is $\cos^4\theta+\sin^4\theta-\tfrac{1}{2}\sin^2 2\theta\cos^2\chi$, which interpolates between identical and distinguishable photons; at the balanced splitter it is $\tfrac{1}{2}\sin^2\chi$, the Hong-Ou-Mandel dip.

- A simulation of $n$ photons needs $n+1$ levels per mode; fewer levels either lose probability or erase the interference.

- The coincidence amplitude is $\mathrm{Perm}(W)$ of the transfer matrix; for fermions it is $\det W$, and for distinguishable photons the coincidence probability is $\mathrm{Perm}(|W|^2)$.

## Part II: Three Photons in a Three-Mode Interferometer

### Any Interferometer Is a Mesh of Beam Splitters: The Reck Decomposition

A passive linear-optical network with $m$ modes (mirrors, beam splitters, phase shifters, no photons created and no light lost) acts on a single photon as an $m\times m$ unitary matrix $W$. The remarkable converse, given by [Reck, Zeilinger, Bernstein, and Bertani in 1994](https://doi.org/10.1103/PhysRevLett.73.58) as a recursive algorithm, is that every $m\times m$ unitary factorizes into a sequence of two-mode beam-splitter transformations, so any interferometer can be built in the laboratory. We will build such a factorization ourselves, with $m$ single-mode phase shifters followed by $m(m-1)/2$ beam splitters, in the order the light meets them. Let's first count degrees of freedom to see why that many is the right number. A unitary in $U(m)$ carries $m^2$ real parameters. Each beam splitter carries two ($\theta$ and $\phi$), which gives $m(m-1)$, and the $m$ phase shifters add the remaining $m$. So the count matches exactly. In other words, beam splitters and phase shifters are all we ever need.

First we need the transfer matrix of a QuantumFramework beam splitter placed on modes $p$ and $q$ of an $m$-mode network. With the phase φ restored, it is the identity except for the $2\times 2$ block $\begin{pmatrix}\cos\theta & -e^{-i\phi}\sin\theta \\ e^{i\phi}\sin\theta & \cos\theta\end{pmatrix}$ on rows and columns $p$ and $q$.

Define the $2\times 2$ block of a beam splitter with angles θ and φ:

```wl
ClearAll[bsBlock]
bsBlock[{t_, f_}] := {{Cos[t], -Exp[-I f] Sin[t]}, {Exp[I f] Sin[t], Cos[t]}}
```

Verify that this block is the transfer matrix that QuantumFramework's `BeamSplitterOperator` produces, for symbolic θ and φ:

```wl
FullSimplify[bsBlock[{\[Theta], \[Phi]}] == transferMatrix[BeamSplitterOperator[{\[Theta], \[Phi]}, 2, {1, 2}, Method -> "Recurrence"], 2, 2], (\[Theta] | \[Phi]) \[Element] Reals]
```

Now the decomposition itself. Take the unitary $W$ and look at two neighboring rows $i-1$ and $i$ of column $j$, with entries $a=W_{i-1,j}$ and $b=W_{ij}$. A beam splitter on modes $i-1$ and $i$ mixes only these two rows. Multiplying them from the left by the inverse of the block with $\tan\theta=|b|/|a|$ and $\phi=\arg(b\,a^*)$, the phase of $b$ relative to $a$, turns the entry $b$ into zero. Sweeping column by column, from the bottom row upward, zeros every entry below the diagonal after $m(m-1)/2$ steps. What remains is an upper-triangular unitary matrix, and an upper-triangular unitary matrix is diagonal: a list of phases $D$. Writing $B_k$ for the transfer matrix of the $k$-th beam splitter of the sweep and $K=m(m-1)/2$ for their number, this means $B_K^{-1}\ldots B_1^{-1}\,W=D$, so $W=B_1\ldots B_K\,D$: in the laboratory, the light first meets the phase shifters $D$ and then the beam splitters in reverse order. This is essentially a QR decomposition by Givens rotations, carried out with physical beam splitters.

Define one nulling step, which records the beam-splitter angles it used and updates the two rows it mixes:

```wl
ClearAll[nullingStep]
nullingStep[w_, {i_, j_}] := With[{a = w[[i - 1, j]], b = w[[i, j]]},
  With[{angles = If[PossibleZeroQ[b], {0, 0}, Simplify[{ArcTan[Abs[a], Abs[b]], Arg[b Conjugate[a]]}]]},
   Sow[{angles, {i - 1, i}}];
   ReplacePart[w, Thread[{i - 1, i} -> ConjugateTranspose[bsBlock[angles]] . w[[{i - 1, i}]]]]]]
```

Define the decomposition, which folds the nulling steps over the entries below the diagonal, column by column and from the bottom up, and returns the remaining diagonal together with the recorded beam splitters:

```wl
ClearAll[reckDecomposition]
reckDecomposition[u_?SquareMatrixQ] := With[{m = Length[u]},
  Reap[Fold[nullingStep, u, Catenate @ Table[{i, j}, {j, m - 1}, {i, m, j + 1, -1}]]]]
```

Define a function that tabulates recorded beam splitters by their angles and the pair of modes each one acts on:

```wl
ClearAll[splitterTable]
splitterTable[splitters_] := TableForm[Append @@@ splitters, TableHeadings -> {None, {"\[Theta]", "\[Phi]", "modes"}}]
```

Let's test it on the $3\times 3$ discrete Fourier matrix $F_{rs}=\omega^{(r-1)(s-1)}/\sqrt{3}$ with $\omega=e^{2\pi i/3}$, which the Wolfram Language provides as [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html).

Decompose the $3\times 3$ Fourier matrix and show the diagonal of phases that remains:

```wl
{fourierDiagonal, {fourierSplitters}} = reckDecomposition[FourierMatrix[3]];
FullSimplify[fourierDiagonal] // MatrixForm
```

Show the angles θ and φ and the mode pairs of the three beam splitters:

```wl
splitterTable[fourierSplitters]
```

As one can see, three beam splitters are enough, with exact angles. The decomposition uses nothing but algebra, so it also accepts symbolic entries; let's try it on a single beam splitter with symbolic angles. The nulling step returns mixing angles between 0 and $\pi/2$ and phases between $-\pi$ and $\pi$, so we give the symbolic angles those ranges.

Decompose the transfer matrix of a beam splitter with symbolic angles α and β, and show the beam splitter it records:

```wl
{symbolicDiagonal, {symbolicSplitters}} = FullSimplify[reckDecomposition[bsBlock[{\[Alpha], \[Beta]}]], 0 < \[Alpha] < \[Pi]/2 && -\[Pi] < \[Beta] < \[Pi]];
splitterTable[symbolicSplitters]
```

Verify that no phases are left over:

```wl
symbolicDiagonal == IdentityMatrix[2]
```

As one can see, the decomposition returns exactly the beam splitter we started from. A phase shifter on mode $p$ is QuantumFramework's `PhaseShiftOperator[α, d, {p}]`, the operator $e^{i\alpha\,\hat n_p}$, which multiplies a photon in mode $p$ by $e^{i\alpha}$. So we can now turn the decomposition into a QuantumFramework circuit.

Define a function that turns any unitary into a QuantumFramework circuit of phase shifters followed by beam splitters, with $d$ levels in every mode:

```wl
ClearAll[interferometer]
interferometer[u_?UnitaryMatrixQ, d_Integer, opts : OptionsPattern[BeamSplitterOperator]] :=
 With[{r = reckDecomposition[u]},
  QuantumCircuitOperator @ Join[
    MapIndexed[PhaseShiftOperator[#1, d, #2] &, Arg @ Diagonal @ First[r]],
    BeamSplitterOperator[#1, d, #2, opts] & @@@ Reverse[Catenate[Last[r]]]]]
```

### The Fourier Tritter: Six Dark Outputs

The $3\times 3$ Fourier matrix describes a symmetric three-port beam splitter, often called a tritter: a single photon entering any port leaves through each of the three ports with probability one third. We will send three photons into it, one per port, so each mode needs four levels.

Build the Fourier tritter with four levels per mode and draw its circuit:

```wl
tritter = interferometer[FourierMatrix[3], 4, Method -> "Recurrence"];
tritter["Diagram"]
```

Verify that the circuit acts on a single photon exactly as the Fourier matrix (two levels per mode are enough to read the transfer matrix):

```wl
FullSimplify[transferMatrix[interferometer[FourierMatrix[3], 2, Method -> "Recurrence"], 3, 2] == FourierMatrix[3]]
```

Three photons can be arranged in three modes in $\binom{3+3-1}{3}=10$ ways, from $|3,0,0\rangle$ to $|1,1,1\rangle$. For distinguishable photons each photon would pick its exit port independently, and every one of these patterns would occur. Let's see what identical photons do.

Send one photon into each mode, $|1,1,1\rangle$, and compute the probabilities of the detection patterns that occur:

```wl
tritterOut = tritter[FockState[{1, 1, 1}, 4]];
tritterProbabilities = KeyMap[First, FullSimplify /@ tritterOut["Probability"]]
```

Verify that the output state is normalized:

```wl
FullSimplify[tritterOut["Norm"]] == 1
```

Have you noticed something interesting? Only four of the ten patterns ever occur: either all three photons leave together through one port, or they leave one per port. The six patterns with exactly two photons in one port are completely dark, for any number of repetitions of the experiment. This is a many-photon generalization of the Hong-Ou-Mandel effect, the zero-transmission (or suppression) law of [Tichy and collaborators (2010)](https://arxiv.org/abs/1002.5038): when $n$ photons enter the $n$ ports of the $n$-mode Fourier interferometer, one per port, every output pattern for which the sum of the output mode labels of all the photons, counted from 0, is not a multiple of $n$ is exactly dark. The law gives a sufficient condition; as its authors note, the converse does not hold in general, so a pattern that passes the test can still be dark. For $n=2$ the coincidence pattern has labels $0+1=1$, which is odd, and we recover the Hong-Ou-Mandel effect.

Verify that the coincidence amplitude of two photons in the $2\times 2$ Fourier matrix, which is a balanced beam splitter up to phases, is zero:

```wl
Permanent[FourierMatrix[2]] == 0
```

Now let's test the law on every pattern of our three photons.

The occupation lists of $n$ photons in $m$ modes are the solutions of $n_1+\ldots+n_m=n$ in nonnegative integers, which [FrobeniusSolve](https://reference.wolfram.com/language/ref/FrobeniusSolve.html) finds.

Define a function that lists all ways of placing $n$ photons in $m$ modes, as occupation lists:

```wl
ClearAll[photonPatterns]
photonPatterns[n_, m_] := ReverseSort @ FrobeniusSolve[ConstantArray[1, m], n]
```

Verify that there are $\binom{n+m-1}{n}$ patterns for three photons in three modes:

```wl
Length[photonPatterns[3, 3]] == Binomial[3 + 3 - 1, 3]
```

Verify that, for three photons, the law accounts for every dark pattern: a pattern is dark exactly when the sum of its photons' mode labels (counted from 0) is not a multiple of 3:

```wl
AllTrue[photonPatterns[3, 3], PossibleZeroQ[Lookup[tritterProbabilities, Key[#], 0]] === ! Divisible[Range[0, 2] . #, 3] &]
```

Recall that $n$ photons need $n+1$ levels in every mode. With three photons we can now see what one level too few does to a suppression law.

Send $|1,1,1\rangle$ through the tritter built with three levels per mode and the default construction, and compute the probabilities of the patterns that occur:

```wl
KeyMap[First, FullSimplify /@ interferometer[FourierMatrix[3], 3][FockState[{1, 1, 1}, 3]]["Probability"]]
```

As one can see, the patterns with two photons in one port, which the law forbids, now carry most of the probability, and the probability of $|1,1,1\rangle$ drops below its value of one third.

Compute the squared norm of the output state when the tritter is built with the recurrence construction instead:

```wl
FullSimplify[interferometer[FourierMatrix[3], 3, Method -> "Recurrence"][FockState[{1, 1, 1}, 3]]["Norm"]^2]
```

So with one level too few, the default construction invents detection events that the law forbids, and the recurrence construction loses more than half of the probability.

### Every Probability Is a Permanent

Recall that the coincidence amplitude of two photons was the permanent of their transfer matrix. The same holds for any number of photons. For $n$ photons entering with occupations $S=(s_1,\ldots,s_m)$ and leaving with occupations $T=(t_1,\ldots,t_m)$, build the $n\times n$ matrix $W_{T,S}$ by taking row $j$ of $W$ as many times as there are photons leaving mode $j$, and column $k$ as many times as there are photons entering mode $k$. Then the amplitude is $\langle T|\hat U|S\rangle=\mathrm{Perm}(W_{T,S})/\sqrt{\prod_k s_k!\,\prod_j t_j!}$, and for distinguishable photons the probability is $\mathrm{Perm}(|W_{T,S}|^2)/\prod_j t_j!$. In words: every way of routing the $n$ photons to the detected modes contributes one product of single-photon amplitudes, and identical photons add all these contributions before squaring.

Define the row (or column) list of an occupation pattern, and the two probabilities:

```wl
ClearAll[modeList, bosonProbability, distinguishableProbability]
modeList[occ_List] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occ]
bosonProbability[w_, s_, t_] := Abs[Permanent[w[[modeList[t], modeList[s]]]]]^2/(Times @@ (s!) Times @@ (t!))
distinguishableProbability[w_, s_, t_] := Permanent[Abs[w[[modeList[t], modeList[s]]]]^2]/Times @@ (t!)
```

For distinguishable photons, each photon picks its exit on its own: the history in which photon $k$ leaves through port $j_k$, for every $k$, has probability $\prod_k|W_{j_k k}|^2$, and the histories that end in the same pattern add up.

Verify, for a general complex $3\times 3$ matrix, that adding up these histories pattern by pattern gives the formula for distinguishable photons:

```wl
With[{wGeneral = Array[\[FormalW], {3, 3}]},
 With[{byPattern = GroupBy[Tuples[Range[3], 3], Lookup[Counts[#], Range[3], 0] &,
     Total[Times @@ MapIndexed[Abs[wGeneral[[#1, First[#2]]]]^2 &, #] & /@ #] &]},
  AllTrue[photonPatterns[3, 3], Expand[byPattern[#] - distinguishableProbability[wGeneral, {1, 1, 1}, #]] === 0 &]]]
```

Verify that every tritter probability computed in Fock space by QuantumFramework equals the permanent formula:

```wl
AllTrue[photonPatterns[3, 3], FullSimplify[Lookup[tritterProbabilities, Key[#], 0] == bosonProbability[FourierMatrix[3], {1, 1, 1}, #]] &]
```

The permanent formula holds for every interferometer, not only for the Fourier one. Take a mesh of three beam splitters with symbolic angles on three modes. The transfer matrix of each splitter is the identity with its $2\times 2$ block placed on its two modes, and the light meets the first splitter first, so the mesh has the transfer matrix $W=B_3B_2B_1$.

Verify that QuantumFramework's Fock-space amplitudes for $|1,1,1\rangle$ through this symbolic mesh are the permanents of $W$:

```wl
With[{splitters = {{{\[Theta]1, \[Phi]1}, {1, 2}}, {{\[Theta]2, \[Phi]2}, {2, 3}}, {{\[Theta]3, \[Phi]3}, {1, 2}}}},
 With[{out = QuantumCircuitOperator[BeamSplitterOperator[#1, 4, #2, Method -> "Recurrence"] & @@@ splitters][FockState[{1, 1, 1}, 4]]["StateVector"],
   w = Dot @@ Reverse[ReplacePart[IdentityMatrix[3], Thread[Tuples[#2, 2] -> Flatten[bsBlock[#1]]]] & @@@ splitters]},
  AllTrue[photonPatterns[3, 3], Simplify[out[[fockIndex[#, 4]]] - Permanent[w[[modeList[#], {1, 2, 3}]]]/Sqrt[Times @@ (#!)]] === 0 &]]]
```

The formula also covers inputs with more than one photon in a mode, through the factor $\prod_k s_k!$.

Send $|2,1,0\rangle$, two photons into the first port and one into the second, through the tritter and compute the probabilities of the patterns that occur:

```wl
probabilities210 = KeyMap[First, FullSimplify /@ tritter[FockState[{2, 1, 0}, 4]]["Probability"]]
```

Verify that every one of these probabilities equals the permanent formula with $S=(2,1,0)$:

```wl
AllTrue[photonPatterns[3, 3], FullSimplify[Lookup[probabilities210, Key[#], 0] == bosonProbability[FourierMatrix[3], {2, 1, 0}, #]] &]
```

As one can see, this input spreads evenly over nine patterns and leaves $|1,1,1\rangle$ dark, a dark pattern about which the law for one photon per port says nothing.

Compare identical photons with distinguishable photons, pattern by pattern:

```wl
TableForm[
 FullSimplify[{bosonProbability[FourierMatrix[3], {1, 1, 1}, #], distinguishableProbability[FourierMatrix[3], {1, 1, 1}, #]} & /@ photonPatterns[3, 3]],
 TableHeadings -> {QuditName /@ photonPatterns[3, 3], {"bosons", "distinguishable"}}]
```

As you may have noticed, distinguishable photons populate all ten patterns, and together the six patterns with two photons in one port carry most of their probability. Identical photons empty exactly those patterns and move their weight to the fully bunched patterns and to $|1,1,1\rangle$.

With the permanent formula we can also test the suppression law for more photons, without building any Fock space. Recall that the law only guarantees darkness; it does not promise that every pattern it allows is bright.

Find the patterns of six photons in the six-mode Fourier interferometer that are dark although the law allows them:

```wl
Select[photonPatterns[6, 6], Divisible[Range[0, 5] . #, 6] && PossibleZeroQ[bosonProbability[FourierMatrix[6], ConstantArray[1, 6], #]] &]
```

As one can see, at six photons interference darkens more patterns than the law predicts. The law is a sufficient condition for darkness, not a necessary one.

Before we move on, let us summarize the most important points we have learned so far:

- Any $m$-mode interferometer is $m(m-1)/2$ beam splitters and $m$ phase shifters, which matches the $m^2$ real parameters of $U(m)$.

- `interferometer[u, d]` builds that mesh from QuantumFramework's `BeamSplitterOperator` and `PhaseShiftOperator`, and its transfer matrix reproduces `u` exactly.

- Three photons in the Fourier tritter populate only the patterns whose mode labels add up to a multiple of 3; the Hong-Ou-Mandel effect is the two-photon case of the same law. For six photons some patterns the law allows are dark as well.

- With one level too few per mode, the tritter either shows patterns the law forbids or loses more than half of its probability, depending on how the beam splitters are built.

- Every Fock-space probability equals $|\mathrm{Perm}(W_{T,S})|^2/\prod s!\prod t!$, also for inputs with several photons in a mode, and distinguishable photons follow $\mathrm{Perm}(|W_{T,S}|^2)/\prod t!$ instead.

## Part III: Boson Sampling with a Random Interferometer

### From Galton Boards to Photons: What Boson Sampling Asks

Picture a Galton board: balls dropped from the top bounce left or right at every peg, and pile up in bins at the bottom. Simulating it is easy: follow one ball at a time, flip a coin at every peg, and record the bin it lands in; the pile of many balls just adds up these independent single-ball outcomes. Distinguishable photons in an interferometer behave exactly like those balls: each photon picks its exit port on its own, with the probabilities $|W_{jk}|^2$, and the pattern of detections adds up the photons' independent choices. Recall that we checked exactly this for three photons: the formula for distinguishable photons is the sum over the histories in which each photon picks its exit on its own. Identical photons cannot be followed one at a time. The amplitude for a detection pattern adds up every way of routing the photons to the detected modes, $n!$ of them, and these ways interfere, as we saw in the tritter. That sum is the permanent.

Boson sampling, proposed by [Aaronson and Arkhipov](https://arxiv.org/abs/1011.3245), is the task of producing detection patterns with the probabilities of $n$ identical photons sent through an $m$-mode interferometer. In other words, a boson sampler is a Galton board for identical photons, and the question is whether a classical computer can imitate it. Their hardness argument, the evidence that a classical computer cannot do this efficiently, takes the interferometer to be chosen at random. A random interferometer means a unitary drawn uniformly from $U(m)$ (the Haar measure), which the Wolfram Language samples with [CircularUnitaryMatrixDistribution](https://reference.wolfram.com/language/ref/CircularUnitaryMatrixDistribution.html). Unlike the Fourier matrix, it has no symmetry that could enforce a suppression law.

Generate a Haar-random $6\times 6$ unitary with a fixed seed, and verify that it is unitary:

```wl
SeedRandom[2026];
u6 = RandomVariate[CircularUnitaryMatrixDistribution[6]];
UnitaryMatrixQ[u6]
```

Build the corresponding QuantumFramework interferometer for three photons (four levels per mode), and count its gates:

```wl
sampler = interferometer[u6, 4];
sampler["GateCount"]
```

These are the $6\cdot 5/2$ beam splitters and the six phase shifters of the mesh.

Verify that the interferometer's transfer matrix reproduces the random unitary to machine precision:

```wl
Max[Abs[transferMatrix[interferometer[u6, 2], 6, 2] - u6]] < 10^-12
```

Send three photons into the first three modes, and collect the probability of every detection pattern:

```wl
inputPattern = {1, 1, 1, 0, 0, 0};
samplerOut = sampler[FockState[inputPattern, 4]];
probabilities = KeyMap[First, samplerOut["Probability"]];
```

Verify that every one of the $\binom{3+6-1}{3}$ possible patterns occurs:

```wl
Length[probabilities] == Binomial[3 + 6 - 1, 3]
```

Verify that the output state is normalized:

```wl
Chop[samplerOut["Norm"] - 1] == 0
```

Here [Chop](https://reference.wolfram.com/language/ref/Chop.html) sets the tiny machine-precision residue of the norm to zero, so the comparison tests the physics rather than the rounding.

Verify that every probability agrees with the permanent formula:

```wl
Max[Abs[Values[probabilities] - (bosonProbability[u6, inputPattern, #] & /@ Keys[probabilities])]] < 10^-12
```

So the QuantumFramework simulation in Fock space and the permanent formula are two representations of the same physics, and they agree to machine precision.

Compute the probability of every pattern for distinguishable photons in the same interferometer:

```wl
distinguishable = AssociationMap[distinguishableProbability[u6, inputPattern, #] &, Keys[probabilities]];
```

Plot the probabilities of the patterns for identical and for distinguishable photons, sorted by the probability for identical photons:

```wl
With[{order = Ordering[Values[probabilities], All, Greater]},
 ListPlot[{Values[probabilities][[order]], Values[distinguishable][[order]]}, Joined -> True, PlotMarkers -> Automatic,
  PlotLegends -> {"identical photons", "distinguishable photons"}, Frame -> True, GridLines -> Automatic,
  AspectRatio -> 1/2, ImageSize -> Large, FrameLabel -> {"pattern (sorted by the identical-photon probability)", "probability"},
  PlotLabel -> "Three photons in a Haar-random six-mode interferometer"]]
```

As one can see, identical photons do not follow the distinguishable ones: pattern by pattern, interference raises some probabilities above the classical value and pushes others below it, and without a symmetry there is no simple rule telling which is which. Every probability is its own permanent. One rule does survive in every interferometer. When all $n$ photons leave through the same mode $j$, the matrix $W_{T,S}$ has $n$ identical rows, so its permanent is $n!\,\prod_k W_{jk}$, while its determinant would vanish. Squaring and dividing by $n!$ gives $n!\,\prod_k|W_{jk}|^2$ for identical photons, exactly $n!$ times the classical value $\prod_k|W_{jk}|^2$.

Verify, for a general complex $6\times 6$ matrix and for two to five photons entering the first modes, that every fully bunched pattern is exactly $n!$ times more likely for identical photons than for distinguishable ones:

```wl
With[{wGeneral = Array[\[FormalW], {6, 6}]},
 AllTrue[Range[2, 5], Function[n,
   AllTrue[n IdentityMatrix[6], FullSimplify[bosonProbability[wGeneral, PadRight[ConstantArray[1, n], 6], #] == n! distinguishableProbability[wGeneral, PadRight[ConstantArray[1, n], 6], #]] &]]]]
```

Define a function that returns the total probability of the bunched patterns, those with two or more photons in some mode, for identical and for distinguishable photons entering the interferometer $u$ with occupations $s$:

```wl
ClearAll[bunchedProbabilities]
bunchedProbabilities[u_, s_] := With[{bunched = Select[photonPatterns[Total[s], Length[s]], Max[#] >= 2 &]},
  Map[f |-> Total[f[u, s, #] & /@ bunched], <|"identical" -> bosonProbability, "distinguishable" -> distinguishableProbability|>]]
```

Compute them for our random interferometer:

```wl
bunchedProbabilities[u6, inputPattern]
```

In this interferometer, as in the Hong-Ou-Mandel effect, identical photons prefer to share modes. Is that a property of this particular interferometer? Averaged over Haar-random interferometers, the bunching of identical photons has an exact answer. The average of $\hat U|S\rangle\langle S|\hat U^\dagger$ over the Haar measure commutes with every interferometer, because the measure does not change when all its unitaries are multiplied by a fixed one. The $n$-photon states span a space on which the interferometers act with no smaller invariant subspace, so by Schur's lemma that average is the identity divided by the dimension $\binom{n+m-1}{n}$: every pattern is equally likely on average. Only $\binom{m}{n}$ of the patterns are collision-free, with no mode receiving more than one photon, so the average probability of the bunched patterns is $1-\binom{m}{n}/\binom{n+m-1}{n}$.

Compute this average for three photons as a function of the number of modes $m$:

```wl
bunchedAverage = FullSimplify[1 - Binomial[m, 3]/Binomial[m + 2, 3], m > 2]
```

Evaluate it for six modes:

```wl
bunchedAverage /. m -> 6
```

Draw a thousand Haar-random six-mode interferometers and compute the bunched probabilities in each (this cell takes a few seconds):

```wl
SeedRandom[2027];
haarSamples = Table[bunchedProbabilities[RandomVariate[CircularUnitaryMatrixDistribution[6]], inputPattern], 1000];
```

Compute the averages over the thousand draws:

```wl
Mean[haarSamples]
```

Compute the standard errors of these averages:

```wl
StandardDeviation[haarSamples]/Sqrt[1000]
```

As one can see, the estimate for identical photons lies within its standard error of the exact average, and distinguishable photons land in bunched patterns far less often. So bunching is not an accident of one interferometer: on average, identical photons share modes more often than distinguishable ones.

### Sampling: What a Boson Sampler Actually Outputs

An experiment never shows the probabilities themselves. Each run ends with detectors clicking, which yields one pattern; repeating the run produces a list of patterns drawn from the distribution above. That list is the output of a boson sampler.

Draw ten detection events from the distribution, the record that ten runs of the experiment would produce:

```wl
RandomChoice[Values[probabilities] -> (QuditName /@ Keys[probabilities]), 10]
```

Every fresh evaluation of this cell gives a different list, exactly like a new batch of runs in the laboratory. The first boson-sampling experiments, reported at the end of 2012, worked at the same scale as this example: [Broome and collaborators](https://doi.org/10.1126/science.1231440) sent two and three photons through a six-mode network, [Spring and collaborators](https://doi.org/10.1126/science.1231692) three and four photons through a six-mode chip, and [Tillmann and collaborators](https://doi.org/10.1038/nphoton.2013.102) three photons through a five-mode circuit.

## Part IV: The Classical Cost of Simulating Identical Photons

### Counting States: The Truncated Fock Space and the Photon Sector

Let's recall how QuantumFramework represented the three photons. Every mode needed $n+1=4$ levels, so the state lives in a truncated Fock space of dimension $(n+1)^m=4^6$. Because the interferometer conserves photon number, all the amplitude stays inside the $n$-photon sector, the $\binom{n+m-1}{n}$ patterns we have been listing, and QuantumFramework stores the state as a sparse array holding only those entries.

Verify that the output state stores exactly one amplitude per three-photon pattern:

```wl
Length[samplerOut["StateVector"]["NonzeroValues"]] == Binomial[3 + 6 - 1, 3]
```

How many modes should we take? Exact sampling means drawing detection patterns from the ideal distribution itself; approximate sampling means drawing them from some distribution close to it, which is the most a real experiment, with its noise and losses, can do. Aaronson and Arkhipov carried out their hardness argument for approximate sampling with many modes, of order $n^5\log^2 n$, and they suspect that an improved analysis could bring this down to $m$ of order $n^2$. Their argument needs at least that many, because with fewer modes almost every detection event puts two photons into the same mode. Recall that, averaged over Haar-random interferometers, every pattern of identical photons is equally likely, so the average probability of the collision-free patterns is $\binom{m}{n}/\binom{n+m-1}{n}$.

Compute this probability for as many modes as photons, $m=n$, as the number of photons grows:

```wl
Limit[Binomial[n, n]/Binomial[2 n - 1, n], n -> Infinity]
```

Compute it for $m=c\,n^2$ modes, with $c$ a positive constant:

```wl
Limit[Binomial[c n^2, n]/Binomial[c n^2 + n - 1, n], n -> Infinity, Assumptions -> c > 0]
```

So with as many modes as photons almost every detection event has a collision, while with a number of modes of order $n^2$ a fixed fraction of the events stays collision-free, and that fraction approaches one as $c$ grows. This is why $n^2$ is the natural scale for the number of modes.

Tabulate both counts for $n$ photons in $m=n^2$ modes:

```wl
TableForm[Table[{n, n^2, (n + 1)^(n^2), Binomial[n + n^2 - 1, n]}, {n, 2, 8}],
 TableHeadings -> {None, {"n", "m", "(n+1)^m", "Binomial[n+m-1,n]"}}]
```

As one can see, the truncated Fock space grows far faster than the photon sector, because it reserves room for every mode to hold every photon at once. The sparse storage keeps QuantumFramework away from the first count, but it still labels each basis state by one index over the whole truncated space.

Check whether the truncated Fock space for five photons in twenty-five modes still fits in the range of machine integers used for that index:

```wl
(5 + 1)^25 <= Developer`$MaxMachineInteger
```

This means that at five photons in twenty-five modes the Fock-space state cannot even be created, although only $\binom{29}{5}$ of its amplitudes would ever be nonzero. The photon sector itself is still modest there. We will come back to it by computing the amplitudes in the photon sector directly.

### Timing QuantumFramework as the Photon Number Grows

To see how the cost grows, we take $n$ photons in $n$ modes, where the truncated space stays indexable up to large $n$, and time the two stages of the simulation separately: building the mesh of beam splitters for a Haar-random unitary, and applying it to $|1,\ldots,1\rangle$. [MaxMemoryUsed](https://reference.wolfram.com/language/ref/MaxMemoryUsed.html)`[expr]` reports the largest amount of memory used while `expr` is evaluated. QuantumFramework also caches the operators it builds, so evaluating these timing cells a second time in the same session measures less work; to repeat a measurement, start a fresh kernel.

Define a function that returns, for $n$ photons in $m$ modes, the time to build the interferometer, the time to apply it, and the peak memory (in megabytes) used while applying it:

```wl
ClearAll[qfCost]
qfCost[n_, m_] := Block[{apply},
  SeedRandom[n m];
  With[{u = RandomVariate[CircularUnitaryMatrixDistribution[m]], input = FockState[PadRight[ConstantArray[1, n], m], n + 1]},
   ClearSystemCache[];
   With[{timedBuild = AbsoluteTiming[interferometer[u, n + 1]]},
    ClearSystemCache[];
    With[{bytes = MaxMemoryUsed[apply = First @ AbsoluteTiming[Last[timedBuild][input]]]},
     {First[timedBuild], apply, N[bytes/2^20]}]]]]
```

Measure the cost for $n$ photons in $n$ modes, $n=2,\ldots,11$ (this cell takes about half a minute):

```wl
costs = Table[qfCost[n, n], {n, 2, 11}];
```

Tabulate the build and apply times:

```wl
TableForm[MapThread[Prepend, {costs[[All, ;; 2]], Range[2, 11]}], TableHeadings -> {None, {"n", "build (s)", "apply (s)"}}]
```

As you can see, building takes most of the time from three photons on. It makes $m(m-1)/2$ beam splitters, each a matrix on the $(n+1)^2$ levels of its two modes, and its cost grows polynomially. The apply time stays small while the photon sector is small and starts to climb only at the end of the range. Applying the gates carries the whole photon sector through the mesh, and the sector grows by nearly a factor of four with every added photon.

Compute the limit of the ratio of successive sector sizes, $\binom{2n+1}{n+1}/\binom{2n-1}{n}$:

```wl
Limit[Binomial[2 n + 1, n + 1]/Binomial[2 n - 1, n], n -> Infinity]
```

Compute the peak memory per amplitude of the photon sector, in bytes, dividing by the sector size $\binom{2n-1}{n}$:

```wl
TableForm[Transpose[{Range[2, 11], 2^20 costs[[All, 3]]/Binomial[2 Range[2, 11] - 1, Range[2, 11]]}], TableHeadings -> {None, {"n", "bytes per amplitude"}}]
```

As one can see, the memory per amplitude falls while the sector is small, where a fixed amount of memory dominates, and once the sector passes a few thousand amplitudes it stops falling and starts to rise: from there on the memory grows at least as fast as the sector does. Extend the range to twelve or thirteen photons in a fresh kernel and watch whether the apply time overtakes the build time, keeping an eye on the memory.

### The Same Physics Without the Truncated Fock Space

The truncated Fock space is bookkeeping: the physics lives in the photon sector. So let's describe the output state directly there, in plain Wolfram Language, in two ways. Both start from one observation: the interferometer transforms each creation operator linearly with its transfer matrix, $\hat U\hat a_k^\dagger\hat U^\dagger=\sum_j W_{jk}\hat a_j^\dagger$.

Let's check this for a single beam splitter, in the untruncated algebra. Since $\hat B=e^{G}$ with the generator $G$ we built earlier, $\hat B\hat a^\dagger\hat B^\dagger=e^{G}\hat a^\dagger e^{-G}=\hat a^\dagger+[G,\hat a^\dagger]+\tfrac{1}{2}[G,[G,\hat a^\dagger]]+\ldots$. If the commutator with $G$ maps the two creation operators linearly into each other, the whole series is the exponential of the $2\times 2$ matrix of that linear map.

Compute the commutators of the generator with the two creation operators:

```wl
commutators = BosonicNormalOrder[Commutator[generator, #], fieldVars, "Scalars" -> {\[Theta], \[Phi]}] & /@ {a1Dag, a2Dag}
```

As one can see, the commutator turns each creation operator into a multiple of the other one.

Verify that the exponential of the matrix of this linear map is the transfer matrix of the beam splitter:

```wl
FullSimplify[MatrixExp[Transpose[Coefficient[#, {a1Dag, a2Dag}] & /@ commutators]] == bsBlock[{\[Theta], \[Phi]}], (\[Theta] | \[Phi]) \[Element] Reals]
```

So a beam splitter rotates the creation operators with its transfer matrix, and a mesh of beam splitters composes these rotations into the transfer matrix $W$ of the whole interferometer. The interferometer also leaves the vacuum alone, so the output state is $\prod_k\big(\sum_j W_{jk}\hat a_j^\dagger\big)|0\rangle$, with one factor for each input photon. Creation operators commute with each other, so this product is an ordinary polynomial in the variables $x_j$ standing for $\hat a_j^\dagger$, and $x^T=\prod_j x_j^{t_j}$ applied to the vacuum is $\sqrt{\prod_j t_j!}\,|T\rangle$. In other words, the output amplitude of a pattern is the coefficient of its monomial, times $\sqrt{\prod t!}$.

Define a function that computes all output amplitudes as the coefficients of $\prod_{k=1}^{n}\big(\sum_j W_{jk}\,x_j\big)$, for $n$ photons entering the first $n$ modes:

```wl
ClearAll[polynomialAmplitudes]
polynomialAmplitudes[w_, n_] := With[{x = Array[\[FormalX], Length[w]]},
  Association[(#1 -> #2 Sqrt[Times @@ (#1!)]) & @@@ CoefficientRules[Expand[Times @@ (x . w[[All, Range[n]]])], x]]]
```

Verify, for a general symbolic $4\times 4$ matrix, that the coefficient of every pattern of four photons in four modes gives the permanent formula:

```wl
With[{wGeneral = Array[\[FormalW], {4, 4}]},
 With[{amps = polynomialAmplitudes[wGeneral, 4]},
  AllTrue[photonPatterns[4, 4], Expand[amps[#] - Permanent[wGeneral[[modeList[#], Range[4]]]]/Sqrt[Times @@ (#!)]] === 0 &]]]
```

Verify that it reproduces every QuantumFramework amplitude of the six-mode example:

```wl
With[{amps = polynomialAmplitudes[u6, 3], v = samplerOut["StateVector"]},
 Max[Abs[(amps[#] - v[[fockIndex[#, 4]]]) & /@ Keys[amps]]] < 10^-12]
```

The second way is the first-quantized picture: $n$ identical photons form a totally symmetric wavefunction of rank $n$ over the $m$ modes. Its component at mode labels $(j_1,\ldots,j_n)$ is $\mathrm{Perm}(W_{T,S})/\sqrt{n!}$ for the pattern $T$ those labels describe (with $S$ one photon in each of the first $n$ modes), and the factor $1/\sqrt{n!}$ normalizes it over all ordered lists of labels. A symmetric tensor is fixed by its components with $j_1\le j_2\le\ldots\le j_n$, one for each pattern, and the Wolfram Language represents exactly this with the structured array [SymmetrizedArray](https://reference.wolfram.com/language/ref/SymmetrizedArray.html): given a rule for the components and the symmetry `Symmetric[All]`, it evaluates the rule once per independent component and stores nothing else.

Define the symmetric wavefunction of $n$ photons entering the first $n$ modes as a `SymmetrizedArray`, whose independent components are permanents:

```wl
ClearAll[symmetricAmplitudes]
symmetricAmplitudes[w_, n_] := SymmetrizedArray[{j__} :> Permanent[w[[{j}, Range[n]]]]/Sqrt[n!], ConstantArray[Length[w], n], Symmetric[All]]
```

Verify that a symmetric tensor of rank 3 over six modes has exactly one independent component per three-photon pattern:

```wl
Length[SymmetrizedIndependentComponents[{6, 6, 6}, Symmetric[All]]] == Binomial[3 + 6 - 1, 3]
```

Verify that the wavefunction of the six-mode example is normalized over all ordered lists of labels:

```wl
Chop[Total[Abs[Flatten[Normal[symmetricAmplitudes[u6, 3]]]]^2] - 1] == 0
```

A pattern $T$ is described by $n!/\prod t!$ ordered lists of labels, the multinomial coefficient that [Multinomial](https://reference.wolfram.com/language/ref/Multinomial.html) computes, so its Fock amplitude is the component times $\sqrt{n!/\prod t!}$.

Verify that this reproduces every QuantumFramework amplitude of the six-mode example:

```wl
With[{sa = symmetricAmplitudes[u6, 3], v = samplerOut["StateVector"]},
 Max[Abs[(sa[[Sequence @@ modeList[#]]] Sqrt[Multinomial @@ #] - v[[fockIndex[#, 4]]]) & /@ photonPatterns[3, 6]]] < 10^-12]
```

Now let's time the three methods on the same interferometers, with few photons in many modes and with as many photons as modes. The QuantumFramework time includes building the gates, as before.

Define a function that times each method on a Haar-random interferometer for $n$ photons in $m$ modes, skipping QuantumFramework where its Fock-space index would overflow:

```wl
ClearAll[methodTimes]
methodTimes[n_, m_] := Block[{u},
  SeedRandom[n + 100 m];
  u = RandomVariate[CircularUnitaryMatrixDistribution[m]];
  ClearSystemCache[];
  {If[(n + 1)^m > Developer`$MaxMachineInteger, Missing["TooLarge"],
     First @ AbsoluteTiming[interferometer[u, n + 1][FockState[PadRight[ConstantArray[1, n], m], n + 1]]]],
   First @ AbsoluteTiming[polynomialAmplitudes[u, n]],
   First @ AbsoluteTiming[symmetricAmplitudes[u, n]]}]
```

Time the three methods (in seconds) for four photons in sixteen modes, five photons in twenty-five modes, eight photons in eight modes, and ten photons in ten modes (this cell takes about a minute):

```wl
TableForm[Join[#, methodTimes @@ #] & /@ {{4, 16}, {5, 25}, {8, 8}, {10, 10}},
 TableHeadings -> {None, {"n", "m", "QuantumFramework", "polynomial", "SymmetrizedArray"}}]
```

As one can see, with few photons in many modes, the regime of the hardness argument, both the polynomial expansion and the `SymmetrizedArray` beat QuantumFramework by orders of magnitude, and they handle five photons in twenty-five modes, where the Fock-space index overflows. There the photon sector is small, while QuantumFramework still builds all $m(m-1)/2$ beam splitters, each on the $(n+1)^2$ levels of its two modes, however few photons there are. As the number of photons approaches the number of modes, the gap narrows: the photon sector becomes large, QuantumFramework's sparse state carries it through the mesh efficiently once the gates are built, and at ten photons in ten modes the `SymmetrizedArray` falls behind both, because it computes one full permanent for every pattern while the polynomial expansion shares partial products between patterns. So the representation matters by orders of magnitude, but no representation removes the growth of the photon sector: every method that produces the whole distribution handles all $\binom{n+m-1}{n}$ amplitudes.

### One Amplitude Is One Permanent: Why the Cost Doubles with Every Photon

Suppose we give up on the whole distribution and ask for a single probability. Even then we need a permanent. Computing the permanent is #P-complete already for matrices of zeros and ones, a 1979 theorem of Valiant on which [Aaronson and Arkhipov](https://arxiv.org/abs/1011.3245) build, and #P-complete problems are believed to be far beyond efficient computation. `Permanent` offers several exact methods, among them Ryser's formula, which adds up one term for each subset of the columns, and Glynn's formula, which adds up one term for each pattern of signs. The cost of both grows like $2^n$ times a polynomial in $n$; Aaronson and Arkhipov quote about $2^{n+1}n^2$ floating-point operations for Ryser's algorithm. Let's measure Glynn's formula.

Define a function that times Glynn's formula for the permanent of the $n\times n$ block of a Haar-random $n^2$-mode unitary, the matrix behind one output amplitude of $n$ photons, using [RepeatedTiming](https://reference.wolfram.com/language/ref/RepeatedTiming.html), which evaluates an expression several times and returns a trimmed mean of the times:

```wl
ClearAll[permanentTime]
permanentTime[n_] := Block[{a},
  SeedRandom[n];
  a = RandomVariate[CircularUnitaryMatrixDistribution[n^2]][[;; n, ;; n]];
  ClearSystemCache[];
  First @ RepeatedTiming[Permanent[a, Method -> "Glynn"]]]
```

Time it for $n=10,\ldots,24$ photons (this cell takes about a minute):

```wl
permanentTimes = Table[{n, permanentTime[n]}, {n, 10, 24}];
```

Plot the timings on a logarithmic scale:

```wl
ListLogPlot[permanentTimes, Joined -> True, PlotMarkers -> Automatic, Frame -> True, GridLines -> Automatic,
 AspectRatio -> 1/2, ImageSize -> Large, FrameLabel -> {"photons n", "seconds"},
 PlotLabel -> "Time to compute one n\[Times]n permanent"]
```

As one can see, the points fall on a straight line, and on a logarithmic scale a straight line means that the time grows exponentially with $n$.

Compute the average factor by which the time grows with each added photon, as the geometric mean over all fourteen steps:

```wl
(permanentTimes[[-1, 2]]/permanentTimes[[1, 2]])^(1/14)
```

Therefore, every added photon roughly doubles the cost of a single amplitude, as the $2^n$ growth of Ryser's and Glynn's formulas predicts. At that rate, fifty photons would put one amplitude far beyond any practical waiting time on a laptop, and a hundred photons beyond any computer. This is also why the exact sampling algorithm of [Clifford and Clifford (2018)](https://arxiv.org/abs/1706.01260), which draws each sample through a chain of conditional probabilities without ever computing the whole distribution, still pays about the price of two $n\times n$ permanents for every sample it draws, when the number of modes grows polynomially with $n$.

### Permanents Versus Determinants: Why Fermions Are Easy

Recall that two identical fermions have the determinant, not the permanent, as their amplitude. The two definitions differ only in the signs of the terms, yet the signs change everything: the alternating signs let the determinant be computed by Gaussian elimination, in a number of operations of order $n^3$, while no such shortcut is known for the permanent.

Time the determinant of a random complex $1000\times 1000$ matrix:

```wl
Block[{a = RandomComplex[{-1 - I, 1 + I}, {1000, 1000}]}, ClearSystemCache[]; First @ AbsoluteTiming[Det[a]]]
```

As one can see, the amplitude of a thousand identical fermions costs less than the permanent of two dozen photons we timed above. As [Aaronson and Arkhipov](https://arxiv.org/abs/1011.3245) point out, this is why noninteracting fermions are easy to simulate classically: the determinant can be computed efficiently, while the permanent is #P-complete. For identical photons they showed that an efficient classical algorithm sampling exactly from the output distribution would collapse the polynomial hierarchy to its third level, a consequence complexity theorists consider very unlikely, and that the same holds even for approximate sampling from random interferometers with many modes, if two conjectures about permanents of Gaussian random matrices hold. The only difference between the two particles is one sign under exchange.

### Where This Leaves Us (and What Comes Next)

You now have a complete, computation-first toolkit for boson sampling: Fock states and beam splitters in QuantumFramework, together with its symbolic algebra of ladder operators, the Hong-Ou-Mandel effect and its dip as the photons become distinguishable, the truncation rule of $n+1$ levels per mode, the Reck decomposition that builds any interferometer from beam splitters and phase shifters, the suppression law of the Fourier interferometer and the dark patterns it misses, the permanent formula for every detection probability, sampling from a random interferometer and the exact Haar average of bunching, two descriptions of the same output state that live in the photon sector alone (a polynomial in creation operators and a `SymmetrizedArray` of permanents), and the three growing costs: the truncated Fock space, the photon sector, and the permanent behind every single amplitude. With these ideas unified, the notebook's message becomes practical: the physics of identical photons is simple to state and simple to simulate for a few photons, and the same interference that makes the photons bunch is what makes them expensive to imitate. From here, the most natural continuations are Gaussian boson sampling, where squeezed light replaces single photons (QuantumFramework's `SqueezeOperator` is the starting point) and hafnians replace permanents; realistic imperfections such as photon loss and partial distinguishability, which we met briefly through the polarization angle χ; and the validation of boson samplers, where suppression laws like the one we found for the tritter help certify that an experiment behaves quantum mechanically.
