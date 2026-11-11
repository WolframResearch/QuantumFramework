---
Template: Default
Title: Why Boson Sampling Is Hard
Author: Mads Bahrami
---

# Why Boson Sampling Is Hard

<!-- #| style: Subtitle -->
Two-photon interference, normalized permanent amplitudes, classical computation, and the assumptions behind sampling hardness: a computation-first narrative with the Wolfram Language and QuantumFramework

<!-- #| style: Author -->
Mads Bahrami (last updated: Nov 2, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

## Modes, Occupation Lists, and Assumptions

An **occupation list** $\{n_1,n_2,\ldots,n_m\}$ records the photon count in each of the $m$ modes, with one entry per mode. Here, each input or output mode corresponds to a port of the optical network. An **input occupation list** specifies how many photons enter each input port; an **output occupation list** specifies how many photons are detected at each output port. The occupation list $\{n_1,n_2,\ldots,n_m\}$ labels the Fock state $|n_1,n_2,\ldots,n_m\rangle$.

A **collision-free** input or output occupation list has at most one photon per mode. For example, $\{1,0,1\}$ describes two photons, one in mode 1 and one in mode 3.

We will follow two photons through a beam splitter, derive the normalized permanent amplitude, compare ways of computing it, and examine what this says about classical sampling. Computing one amplitude, collecting every probability, and drawing one sample are different tasks; keeping them separate is essential to the hardness argument.

The examples assume a passive lossless network and number-resolving detectors. The photons may enter different ports, but their internal wave packets must match in arrival time, spectrum, polarization, and unresolved spatial structure to interfere fully. Partial distinguishability reduces that interference. Loss changes the detected photon-number sectors, while threshold detectors merge different occupation lists into the same click pattern; neither alone makes the photons distinguishable.

Load QuantumFramework and its second-quantization tools, then evaluate the cells in order:

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

Let's start!

## Two-Photon Interference at a Beam Splitter

Start with the smallest interesting case: two identical photons, one into each port of a single balanced beam splitter. A balanced beam splitter sends each incoming photon out either side with equal probability, so a classical intuition says the two photons land on the same side half the time and on opposite sides half the time. The quantum answer is stranger: the two photons never come out on opposite sides. They always leave together. In other words, the output occupation list $\{1,1\}$, one photon in each output mode, has probability zero. This is the Hong-Ou-Mandel effect, and it is the whole subject in miniature.

Build a two-mode interferometer that is a single balanced beam splitter (the angle $\theta=\pi/4$ makes the reflection and transmission amplitudes equal in size):

```wl
qbs = QuantumInterferometer[{"BS"[\[Pi]/4, 0] -> {1, 2}}, 2]
```

Draw it, so we can see the one element the photons meet:

```wl
qbs["Diagram"]
```

Everything a passive network of $m$ modes does to a single photon is held in one $m\times m$ matrix $W$, its transfer matrix: the entry $W_{jk}$ is the amplitude for a photon entering port $k$ to leave port $j$. For the beam splitter $m=2$, and we store its matrix as `Wbs`.

Read the transfer matrix:

```wl
Wbs = qbs["TransferMatrix"];
Wbs // MatrixForm
```

As one can see, the four entries are equal in magnitude but not in sign: one of them carries a minus. That sign is not decoration. It is the seed of the interference, and we are about to watch it cancel an entire outcome.

Ask `QuantumInterferometer` for the full distribution over output occupation lists when one photon enters each port, giving the input occupation list $\{1,1\}$:

```wl
qbs["Probabilities", {1, 1}]
```

As one can see, the output occupation list $\{1,1\}$ has probability zero, while the two bunched output occupation lists $\{2,0\}$ and $\{0,2\}$ split the probability between them. The photons leave together, every time. The next question is where that zero comes from, and the cleanest way to see it is to compute the same three numbers for three different kinds of particle and read them side by side.

### Occupation Lists and the Permanent Formula

To compare the three particle types, we need the **permanent** of a square matrix. It is defined like the determinant, as a sum over all ways of pairing rows with columns, except that every term is added with a plus sign where the determinant alternates the signs.

The probability helpers below assume valid input and output occupation lists with the same positive total photon count.

Define a helper, `modeList`, that converts an input or output occupation list into a **mode assignment list**: a list containing the mode label once for each photon in that mode. Replace each count $n_j$ with $n_j$ copies of the mode label $j$, then join the groups into one list. For example,

$$
\{2,0,1\}\;\longrightarrow\;\{\{1,1\},\{\},\{3\}\}\;\longrightarrow\;\{1,1,3\}.
$$

The occupation list has one entry per mode; the mode assignment list has one entry per photon. These entries are bookkeeping for the occupied modes and do not make the photons distinguishable.

```wl
ClearAll[modeList];
modeList[occ_List] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occ]
```

Name the three output occupation lists we are comparing:

```wl
outs = {{2, 0}, {1, 1}, {0, 2}};
```

For indistinguishable bosons, let $S=\{s_1,\ldots,s_m\}$ and $T=\{t_1,\ldots,t_m\}$ be the input and output occupation lists, with $\sum_k s_k=\sum_j t_j=n$, since a passive lossless network neither creates nor destroys photons. Define `bosonProb` to compute the probability of observing $T$ given $S$ and the transfer matrix:

```wl
ClearAll[bosonProb];
bosonProb[w_, s_, t_] := Abs[Permanent[w[[modeList[t], modeList[s]]]]]^2 / Times @@ Factorial[Join[s, t]]
```

To unpack this calculation, start with the device's transfer matrix $W$, whose entries are single-photon amplitudes. The input mode assignment list selects columns, and the output mode assignment list selects rows: repeat column $k$ $s_k$ times and row $j$ $t_j$ times. Call the resulting $n\times n$ matrix $M=W_{T,S}$.

Each column of $M$ represents an incoming photon, and each row an outgoing detection. A complete connection between them uses every row and column once. Because the photons do not interact, its contribution is the product of the corresponding single-photon amplitudes. For example, one contribution is $M_{11}M_{22}\cdots M_{nn}$.

Identical bosons require us to add these contributions with no extra exchange signs. Summing over all $n!$ permutations $\sigma$ of the columns gives precisely the permanent:

$$
\operatorname{Per}(M)=\sum_{\sigma}\prod_{a=1}^{n}M_{a,\sigma(a)}.
$$

Normalizing the input and output Fock states then gives the amplitude:

$$
A(T\leftarrow S)=\frac{\operatorname{Per}(M)}{\sqrt{\prod_k s_k!\,\prod_j t_j!}}.
$$

The factorials account for the normalization of states with multiple photons in the same mode. To obtain the probability, take the squared magnitude of the amplitude:

$$
P(T\mid S)=|A(T\leftarrow S)|^2=\frac{|\operatorname{Per}(M)|^2}{\prod_k s_k!\,\prod_j t_j!}.
$$

This is what `bosonProb` evaluates: `Abs[...]^2` takes the squared magnitude of the permanent, and `Times @@ Factorial[Join[s, t]]` multiplies the factorials of all input and output occupation counts.

### Comparing Particle Statistics

For classical distinguishable particles the amplitudes never interfere. Each photon crosses the network on its own, leaving input mode $k$ for output mode $j$ with probability $|W_{jk}|^2$, so the probability of an output occupation list is the permanent of the block of entrywise squared magnitudes $|W|^2_{T,S}$, selected by the same repeated rows and columns, divided only by the output factorials:

$$ P_{\text{dist}}(T\mid S) = \frac{\operatorname{Per}\bigl(|W|^2_{T,S}\bigr)}{\prod_j t_j!} $$

The permanent counts each assignment of photons to output modes $\prod_j t_j!$ times, once for every way of matching the photons in mode $j$ to the $t_j$ copies of row $j$, and the output factorials remove that overcount. There is no $s_k!$, because photons that start in the same mode are still different particles. Define that:

```wl
ClearAll[distProb];
distProb[w_, s_, t_] := Permanent[(Abs[w]^2)[[modeList[t], modeList[s]]]] / Times @@ Factorial[t]
```

For indistinguishable noninteracting fermions, odd permutations contribute a minus sign, giving a **determinant**. A valid input occupation list has only zero-or-one counts; output collisions have probability zero, and the definition returns that zero directly instead of relying on a numerically computed determinant of a matrix with a repeated row. Define the probability for such inputs:

```wl
ClearAll[fermionProb];
fermionProb[w_, s_, t_] := If[Max[t] > 1, 0, Abs[Det[w[[modeList[t], modeList[s]]]]]^2]
```

Put the three particle types in one table, with one row per output occupation list:

```wl
Grid[
 Join[
  {Style[#, Bold] & /@ {"output", "indistinguishable bosons", "distinguishable bosons", "indistinguishable fermions"}},
  Table[Prepend[Through[{bosonProb, distProb, fermionProb}[Wbs, {1, 1}, t]], Row[{"|", Row[t, ","], "\[RightAngleBracket]"}]], {t, outs}]],
 Frame -> All, Spacings -> {2, 1}]
```

The distinguishable probabilities are a quarter, a half, a quarter. Bosonic coincidence vanishes because its two contributions are $+\tfrac12$ and $-\tfrac12$: positive exchange coefficients do not prevent destructive interference between signed optical amplitudes. Fermions instead leave one particle in each output mode with certainty. A bunched output repeats a matrix row, making its determinant zero, as required by Pauli exclusion.

Check the coincidence amplitude at a variable beam-splitter angle:

```wl
ClearAll[\[Theta]];
With[{splitter = QuantumInterferometer[{"BS"[\[Theta], 0] -> {1, 2}}, 2]},
 FullSimplify[splitter["Amplitude", {1, 1} -> {1, 1}] == Cos[2 \[Theta]],
  \[Theta] \[Element] Reals]
]
```

At $\theta=0$ this amplitude is one, restoring the input; at $\theta=\pi/4$ it vanishes. We can now derive the amplitude for an arbitrary passive network.

## Many-Photon Amplitudes as Normalized Permanents

The same amplitude formula follows directly from how a passive network acts on the creation operators. Each network sends the operator that creates a photon in input port $k$ to a superposition over the outputs,

$$ a_k^\dagger \;\longmapsto\; \sum_j W_{jk}\, a_j^\dagger . $$

Take an input occupation list $S=\{s_1,\ldots,s_m\}$, with $s_k$ photons entering port $k$ and $n=\sum_k s_k$ photons in total. Its Fock state is $|S\rangle=\prod_k \frac{(a_k^\dagger)^{s_k}}{\sqrt{s_k!}}\,|0\rangle$. Substituting the transformed operator for every $a_k^\dagger$ and expanding the product gives

$$
\prod_k\left(\sum_j W_{jk}\,a_j^\dagger\right)^{s_k}
      =\sum_T c_T\prod_j (a_j^\dagger)^{t_j},
$$

where the sum runs over all output occupation lists $T=\{t_1,\ldots,t_m\}$ with $\sum_j t_j=n$. The creation operators commute with one another, so this expansion is ordinary polynomial algebra, and the coefficient $c_T$ collects all contributions to $T$. Acting on the vacuum, $\prod_j (a_j^\dagger)^{t_j}|0\rangle=\sqrt{\prod_j t_j!}\,|T\rangle$, so the input normalization and this factor together give

$$
A(T\leftarrow S)=c_T\sqrt{\frac{\prod_j t_j!}{\prod_k s_k!}},
\qquad
\operatorname{Per}(W_{T,S})=c_T\prod_j t_j!.
$$

Thus coefficient extraction and the permanent formula compute the same complex amplitude. Because the creation operators commute, the same expansion can be carried out with ordinary variables $x_j$ in place of $a_j^\dagger$, in the polynomial $F_S(x_1,\ldots,x_m)=\prod_k(\sum_j W_{jk}x_j)^{s_k}$. Its coefficients give every amplitude at once, which is how `QuantumInterferometer` computes the full output distribution below.

To watch this on a network with no special structure, draw one at random. A random unitary is the natural "generic" network, and `QuantumInterferometer["Random"[m]]` builds an $m$-mode one and factors it into beam splitters and phase shifters.

Build a random five-mode interferometer with a fixed seed, so the run is reproducible:

```wl
SeedRandom[2026];
qi = QuantumInterferometer["Random"[5]]
```

Draw the mesh the photons will cross:

```wl
qi["Diagram"]
```

Read its transfer matrix:

```wl
W = qi["TransferMatrix"];
W // MatrixForm
```

Now use input occupation list $\{1,1,1,0,0\}$, with one photon entering each of ports 1, 2, and 3, and ask for the amplitude of output occupation list $\{1,0,1,1,0\}$, with one photon leaving each of ports 1, 3, and 4. The corresponding mode assignment lists select columns $\{1,2,3\}$ and rows $\{1,3,4\}$ of $W$. Both occupation lists are collision-free, so all factorials in the normalization are one and the amplitude is the permanent of the selected $3\times3$ matrix.

Compute that amplitude directly from `QuantumInterferometer`:

```wl
qi["Amplitude", {1, 1, 1, 0, 0} -> {1, 0, 1, 1, 0}]
```

Confirm it equals the permanent of the matrix selected by the same input and output occupation lists:

```wl
With[{in = {1, 1, 1, 0, 0}, out = {1, 0, 1, 1, 0}},
 qi["Amplitude", in -> out] ===
  Permanent[W[[modeList[out], modeList[in]]]]
]
```

For this collision-free example, `"Amplitude"` uses the direct permanent route, so this checks selection rather than independently verifying the formula.

Now put two photons in port 1 and one in port 3: $S=\{2,0,1,0,0\}$. Choose $T=\{1,2,0,0,0\}$. The selected columns are $\{1,1,3\}$ and rows $\{1,2,2\}$, so both input and output collisions enter the normalization: $\sqrt{2!\,2!}=2$.

Check the amplitude against the normalized permanent:

```wl
With[{in = {2, 0, 1, 0, 0}, out = {1, 2, 0, 0, 0}},
 Abs[qi["Amplitude", in -> out] -
   Permanent[W[[modeList[out], modeList[in]]]]/2] < 10^-12
]
```

The formula does not prescribe an algorithm. For approximate transfer matrices, `QuantumInterferometer` may group repeated occupations, using only occupied rows and columns and choosing the orientation with fewer summands. The grouped sum has

$$
K=\min\left\{\prod_k(s_k+1),\;\prod_j(t_j+1)\right\}
$$

terms: six in this example, versus $2^3=8$ for collision-free input and output. Exact transfer matrices retain the direct `Permanent` route. These are evaluations of the same quantity; `bosonProb` remains our simple reference.

## Computing and Checking the Full Output Distribution

A physical boson sampler produces output occupation lists. To inspect the probabilities underlying those samples, collect the whole distribution for input occupation list $\{1,1,1,0,0\}$:

```wl
probs = qi["Probabilities", {1, 1, 1, 0, 0}]
```

Show the five most likely output occupation lists:

```wl
TakeLargest[probs, 5]
```

There is one probability for each output occupation list. For $n$ photons in $m$ modes there are $\binom{n+m-1}{n}$ such lists, counting every way to place the photons in the modes.

Confirm the distribution has exactly that many entries:

```wl
Length[probs] == Binomial[3 + 5 - 1, 3]
```

Confirm it is a genuine probability distribution, summing to one (`Chop` sets the tiny machine-precision residue to zero so the comparison tests the physics, not the rounding; see [Chop](https://reference.wolfram.com/language/ref/Chop.html)):

```wl
Chop[Total[probs] - 1] == 0
```

For passive Fock inputs, `"Probabilities"` obtains all amplitudes by expanding $F_S$ and extracting its coefficients, then squares their magnitudes. `FrobeniusSolve` only enumerates the possible output occupation lists; it does not calculate amplitudes. This is a different algorithm from separately evaluating selected permanents.

First compare the complex amplitudes, so relative-phase errors cannot disappear upon squaring. The call `qi[in]` returns the coefficient-based amplitude association:

```wl
With[{in = {1, 1, 1, 0, 0}}, {amps = qi[in]},
   Max[Abs[Lookup[amps, QuditName[#], 0] -
              Permanent[W[[modeList[#], modeList[in]]]]/
                Sqrt[Times @@ Factorial[Join[in, #]]]] & /@
         FrobeniusSolve[ConstantArray[1, 5], Total[in]]] < 10^-12
  ]
```

Verify the probabilities independently as well:

```wl
Max[Abs[Lookup[probs, QuditName[#], 0] - bosonProb[W, {1, 1, 1, 0, 0}, #]] & /@
   FrobeniusSolve[{1, 1, 1, 1, 1}, 3]] < 10^-12
```

Both routes agree within numerical precision. Polynomial expansion shares work across outputs, but it is not automatically a polynomial-time algorithm: a generic expansion has $\binom{n+m-1}{n}$ coefficients, which can grow exponentially when the photon and mode counts grow together. A separate permanent evaluation for every output is one approach, not a requirement. What does one generic permanent cost?

## The Cost of One Amplitude

How does the cost of `qi["Amplitude", S -> T]` grow with photon count? Compare collision-free, paired, and fully bunched occupation lists in random networks with $m=2n$ modes. The input occupies the first half of the ports; the output repeats its occupation pattern in the second half.

For six photons in 12 modes, the three input occupation lists are:

- **Collision-free:** each photon occupies its own mode. For six photons, $\{1,1,1,1,1,1,0,0,0,0,0,0\}$.
- **Paired:** photons occupy modes in groups of two, with one singly occupied mode if the total is odd. For six photons, $\{2,2,2,0,0,0,0,0,0,0,0,0\}$.
- **Fully bunched:** all photons occupy one mode, $\{6,0,0,0,0,0,0,0,0,0,0,0\}$.

Define the three input occupation lists. For pairs, `Quotient` gives the number of doubly occupied ports and `Mod` supplies any remaining photon; `PadRight` fills unused ports with zeros:

```wl
ClearAll[occupations];
occupations[n_] := PadRight[#, 2 n] & /@ <|
 "collision-free" -> ConstantArray[1, n],
 "pairs" -> Append[ConstantArray[2, Quotient[n, 2]], Mod[n, 2]],
 "bunched" -> {n}
|>
```

Display the input occupation lists for six photons:

```wl
occupations[6]
```

Build a reproducible network for each photon count. `BlockRandom` keeps the seed local, so this does not change the random state of the surrounding calculation:

```wl
ClearAll[randomInterferometer];
randomInterferometer[n_] := BlockRandom[
 SeedRandom[1107 + 1000 n];
 QuantumInterferometer["Random"[2 n]]
]
```

Each measurement uses one `AbsoluteTiming` evaluation after `ClearSystemCache[]`. Construct the network and output occupation list before starting the clock. The result is `{seconds, amplitude}`:

```wl
ClearAll[timeAmplitude];
timeAmplitude[q_, s_List] := With[{t = RotateRight[s, Total[s]]},
 ClearSystemCache[];
 AbsoluteTiming[q["Amplitude", s -> t]]
]
```

Use the same network for all three occupation patterns and keep their names as association keys:

```wl
ClearAll[amplitudeBenchmark];
amplitudeBenchmark[n_] := With[{q = randomInterferometer[n]},
 timeAmplitude[q, #] & /@ occupations[n]
]
```

Use photon counts up to 24 for the plot; shorten the list for a quicker run:

```wl
photonCounts = {2, 6, 10, 14, 18, 22, 24};
```

Run the benchmark, indexing the results first by photon count and then by occupation pattern:

```wl
amplitudeTimings = AssociationMap[amplitudeBenchmark, photonCounts];
```

Group the `{photon count, seconds}` pairs by occupation pattern. `KeyValueMap` retains each photon count, and `Merge` collects the three curves:

```wl
plotData = Merge[
 KeyValueMap[Function[{n, results}, {n, First[#]} & /@ results], amplitudeTimings],
 Identity];
```

Plot the cost of one QI amplitude:

```wl
ListLogPlot[Values[plotData],
 Joined -> True, PlotMarkers -> Automatic, Frame -> True, GridLines -> Automatic,
 AspectRatio -> 1/2, ImageSize -> Large,
 FrameLabel -> {"photons n (m = 2 n modes)", "seconds per QI amplitude"},
 PlotLegends -> Keys[plotData],
 PlotLabel -> "Occupation structure changes the cost of one amplitude"]
```

The collision-free and paired timings rise rapidly, while QI exploits repeated occupations and occupied modes to keep the fully bunched case cheap. Here every selected entry is $W_{n+1,1}$, so the normalized amplitude reduces to $W_{n+1,1}^{\,n}$. This structured case does not represent a generic permanent. Submillisecond single-call timings are noisy; focus on the large separation in cost. For the collision-free and paired cases, QI calls the built-in `Permanent` directly, so those two curves measure the cost of a general permanent and track each other.

For a general dense complex $n\times n$ matrix, algorithms such as Ryser's and Glynn's take $O(n2^n)$ arithmetic operations to evaluate its permanent, while Gaussian elimination evaluates its determinant in $O(n^3)$. These timings characterize particular implementations and occupation patterns. Computing a general permanent exactly is #P-hard (Valiant, 1979), a worst-case result that permits easy structured instances.

## From Amplitude Cost to Sampling Hardness

A hard probability calculation does not by itself establish hard sampling. A sampler can avoid constructing the full distribution, and even a distribution written using permanents may admit easy samples.

Distinguishable photons supply a counterexample. Their output probabilities involve permanents of entrywise $|W|^2$, yet sampling needs only independent draws from its columns: draw an output for each incoming photon, then count the results.

Confirm that each column of $|W|^2$ is a genuine probability distribution over the output modes, by checking that every column sums to one:

```wl
Chop[Total[Abs[W]^2] - 1]
```

This verifies the normalized single-photon distributions used by that sampler, without evaluating any permanent.

Define that sampler. Each photon in the input picks its output mode from its own column of $|W|^2$, and counting the picks gives an output occupation list:

```wl
ClearAll[distinguishableSample];
distinguishableSample[w_, s_] := With[{m = Length[w]},
 Lookup[Counts[RandomChoice[Abs[w[[All, #]]]^2 -> Range[m]] & /@ modeList[s]], Range[m], 0]
]
```

Indistinguishable bosons have no such shortcut. The simplest exact sampler enumerates every output occupation list, weights each by `bosonProb`, and draws from that list:

```wl
ClearAll[bosonSamples];
bosonSamples[w_, s_, k_] := With[{outputs = FrobeniusSolve[ConstantArray[1, Length[w]], Total[s]]},
 RandomChoice[(bosonProb[w, s, #] & /@ outputs) -> outputs, k]
]
```

Draw ten thousand samples of each kind for input occupation list $\{1,1,1,0,0\}$ and record how often at least two photons leave through the same port:

```wl
SeedRandom[2027];
With[{in = {1, 1, 1, 0, 0}, k = 10^4},
 N @ <|
   "indistinguishable bosons" -> Mean[Boole[Max[#] > 1] & /@ bosonSamples[W, in, k]],
   "distinguishable bosons" -> Mean[Boole[Max[#] > 1] & /@ Table[distinguishableSample[W, in], k]]
 |>
]
```

Indistinguishable bosons leave through a shared port much more often than distinguishable photons sent through the same network: the bunching seen at the beam splitter survives in a generic network and shows up in the samples themselves, with no amplitude in sight. The boson sampler paid for this by computing every probability first, enumerating all $\binom{n+m-1}{n}$ output occupation lists; a useful classical sampler would have to avoid exactly that cost.

The missing bridge is a complexity result, and its mechanism fits in a few steps. A classical sampler is a program that turns random bits into an output occupation list, so the probability of any chosen output is the fraction of random-bit strings that produce it. Stockmeyer's approximate counting estimates such a fraction to within a constant factor in randomized polynomial time with an NP oracle, a class that lies in the third level of the polynomial hierarchy. If a polynomial-time classical algorithm sampled exactly from the boson distribution, this would estimate $|\operatorname{Per}(M)|^2$ to within a constant factor, a task that is itself #P-hard. By Toda's theorem every problem in the polynomial hierarchy reduces to a #P problem, so the whole hierarchy would collapse to its third level, contrary to a standard complexity assumption ([Aaronson and Arkhipov, 2010](https://arxiv.org/abs/1011.3245)).

Real devices are imperfect, so the claim that matters concerns approximate sampling: producing samples from a distribution $Q$ whose total-variation distance from the ideal distribution $P$ satisfies

$$
d_{\mathrm{TV}}(P,Q)=\tfrac12\sum_T|P(T)-Q(T)|\le\varepsilon,
$$

in time polynomial in $n$ and $1/\varepsilon$. Such a sampler may put all of its error on the outputs we ask about, so the argument hides the matrix of interest as a random $n\times n$ block of a larger network, and the sampler cannot tell which block matters. That guarantees only a small additive error on most outputs. Two further assumptions finish the argument, each with its own job. Permanent anticoncentration says typical values of $|\operatorname{Per}|^2$ are not much smaller than their mean, which turns the additive error into a constant-factor one. Average-case hardness says that estimating the permanents of Gaussian random matrices to that accuracy remains #P-hard. Under these assumptions, polynomial-time approximate sampling would again collapse the polynomial hierarchy; they do not prove that every classical algorithm takes exponential time.

The hiding step is where the original argument needs a dilute regime, with many more modes than photons: only then does a small block of a random network look like a matrix of independent Gaussian entries. That is not a universal restriction: later work extends hardness evidence to a linear number of modes, where collisions are common ([Bouland et al., 2026](https://journals.aps.org/prx/abstract/10.1103/xc7b-sjm5)).

Classical algorithms also improve on full-distribution enumeration. Clifford and Clifford give exact sampling in $O(n2^n+\operatorname{poly}(m,n))$ time ([Clifford and Clifford, 2018](https://arxiv.org/abs/1706.01260)), with further average-case improvements when the number of modes grows only in proportion to the photon count ([Clifford and Clifford, 2020](https://arxiv.org/abs/2005.04214)). Algorithmic improvements and conditional hardness coexist: the former exploit structure or share work, while the latter concern suitable growing families of instances.

## Summary and Next Questions

We have connected two-photon interference to normalized permanent amplitudes, checked those amplitudes and probabilities against an independent coefficient calculation, measured how photon count and occupation structure affect QI amplitude calculations, and drawn samples in which bosonic bunching is visible directly. The computations explain the optical probabilities; the cited complexity results supply the separate, conditional argument about classical sampling.

The next questions concern partial distinguishability, loss, and detector resolution, each of which changes the experiment differently. Replacing single-photon inputs by squeezed vacuum leads to Gaussian boson sampling, where hafnians take the role played here by permanents.
