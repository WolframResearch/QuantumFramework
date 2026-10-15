---
Template: Default
Title: Why Boson Sampling Is Hard
Author: Mads Bahrami
---

# Why Boson Sampling Is Hard

<!-- #| style: Subtitle -->
Two-photon interference, permanent amplitudes, the cost of one permanent against its sign-flipped twin the determinant, and an honestly cited bridge to sampling hardness: a computation-first narrative with the Wolfram Language and QuantumFramework

<!-- #| style: Author -->
Mads Bahrami (last updated: Oct 14, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

### Setting the Stage: How This Notebook Flows

This notebook follows one thread from a single beam splitter to the edge of computational complexity: two identical photons meeting at a beam splitter refuse to come out apart, the amplitude for any detection pattern turns out to be a permanent of the network's transfer matrix, one permanent is empirically exponential to compute while its sign-flipped twin the determinant is cheap, and the distance between those two is the entire distance between fermions and bosons. I want to be careful about one thing throughout, because it is the place where this subject is most often oversold: computing one probability, computing the whole distribution, and drawing a sample from that distribution are three different tasks, and only the first two are what the permanent directly controls.

In other words, I've tried to build a catalogue of computational experiments that earn each step by running it. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. So I will compute the interference, compute the permanents, time them, and only then state plainly which part of "boson sampling is hard" the computation has actually established and which part rests on a conjecture I can cite but not run.

The environment you see is a live [Wolfram notebook](https://www.wolfram.com/notebooks/). Evaluate the cells from top to bottom; a few later cells use a network and a distribution defined earlier, so be mindful of those dependencies. Two cells in the second half time a computation on purpose and take up to about a minute; only their trend matters, not the absolute seconds on your machine.

The story runs as a continuous sequence, like a movie, with a few headings only as transitions. My suggestion is to focus on each output and its meaning first, and worry about how the input cell is written afterward. And remember you are not locked into the code: change the number of photons, the number of modes, the seed, and watch the argument bend or hold.

If you have any suggestions or questions, please reach out to us at [quantum@wolfram.com](mailto:quantum@wolfram.com)

### Prerequisites and How to Read This

You only need to be comfortable with the idea of photons occupying modes: a state $|n_1,n_2,\ldots\rangle$ says how many photons sit in each mode, and a detector reports one such pattern. Everything optical is handled by `QuantumInterferometer`, a QuantumFramework tool that builds a passive optical network and hands back its amplitudes and outcome probabilities directly. The two combinatorial objects at the heart of the story, the permanent and the determinant, are built into the Wolfram Language.

Two physical assumptions sit underneath everything below, and it is worth naming them now because relaxing them is exactly where real experiments live. "Identical" photons means photons indistinguishable in every internal degree of freedom that a detector does not resolve: arrival time, spectrum, polarization, spatial mode. And the occupation distributions assume number-resolving detectors and a lossless network. Drop either assumption and the interference degrades toward the classical case, which is one reason loss and distinguishability are studied as more than nuisances.

Let's start!

## Two Photons and One Beam Splitter: Where the Coincidences Go

Start with the smallest interesting case: two identical photons, one into each port of a single balanced beam splitter. A balanced beam splitter sends each incoming photon out either side with equal probability, so a classical intuition says the two photons land on the same side half the time and on opposite sides half the time. The quantum answer is stranger: the two photons never come out on opposite sides. They always leave together. In other words, the coincidence outcome, one photon in each output, has probability zero. This is the Hong-Ou-Mandel effect, and it is the whole subject in miniature.

Build a two-mode interferometer that is a single balanced beam splitter (the angle $\theta=\pi/4$ makes the reflection and transmission amplitudes equal in size):

```wl
qbs = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]
```

Draw it, so we can see the one element the photons meet:

```wl
qbs["Diagram"]
```

Everything the beam splitter does to a single photon is held in one $2\times 2$ matrix $W$, its transfer matrix: the entry $W_{jk}$ is the amplitude for a photon entering port $k$ to leave port $j$.

Read the transfer matrix:

```wl
Wbs = qbs["TransferMatrix"];
Wbs // MatrixForm
```

As one can see, the four entries are equal in magnitude but not in sign: one of them carries a minus. That sign is not decoration. It is the seed of the interference, and we are about to watch it cancel an entire outcome.

Ask `QuantumInterferometer` for the full distribution of detection patterns when one photon enters each port, the input $|1,1\rangle$:

```wl
qbs["Probabilities", {1, 1}]
```

As one can see, the coincidence pattern $|1,1\rangle$ has probability zero, while the two "bunched" patterns $|2,0\rangle$ and $|0,2\rangle$ split the probability between them. The photons leave together, every time. The next question is where that zero comes from, and the cleanest way to see it is to compute the same three numbers for three different kinds of particle and read them side by side.

To do that by hand we need the combinatorial object that assembles a many-photon amplitude from the single-photon matrix. For several identical photons the amplitude of a detection pattern is not one entry of $W$ but a sum over every way of routing the photons to the detected modes. That sum is the **permanent**: the permanent of a square matrix is defined exactly like the determinant, as a sum over all ways of pairing rows with columns, except that every term is added with a plus sign where the determinant alternates the signs.

Define a helper that turns an occupation list into the rows (or columns) of the matrix a pattern selects, repeating mode $j$ once for every photon in it:

```wl
ClearAll[modeList];
modeList[occ_List] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occ]
```

Name the three output patterns we are comparing:

```wl
outs = {{2, 0}, {1, 1}, {0, 2}};
```

For indistinguishable bosons the probability of output $T$ from input $S$ is the squared permanent of the block of $W$ selected by $T$ and $S$, divided by a normalization that counts photons sharing a mode. Define that:

```wl
ClearAll[bosonProb];
bosonProb[t_] := Abs[Permanent[Wbs[[modeList[t], modeList[{1, 1}]]]]]^2 / (Times @@ (t!) Times @@ ({1, 1}!))
```

For classical distinguishable particles the amplitudes never interfere; the probability is instead the permanent of the matrix of squared magnitudes $|W|^2$, which is a plain bookkeeping of independent paths. Define that:

```wl
ClearAll[distProb];
distProb[t_] := Permanent[(Abs[Wbs]^2)[[modeList[t], modeList[{1, 1}]]]] / (Times @@ (t!))
```

For identical fermions the exchange sign is negative, so the permanent becomes a **determinant**. Define that, with the squared magnitude giving the probability:

```wl
ClearAll[fermionProb];
fermionProb[t_] := Abs[Det[Wbs[[modeList[t], modeList[{1, 1}]]]]]^2
```

Put the three particle types in one table, output pattern by output pattern:

```wl
Grid[
 Join[
  {Style[#, Bold] & /@ {"output", "bosons", "distinguishable", "fermions"}},
  Table[{Row[{"|", Row[t, ","], "\[RightAngleBracket]"}], bosonProb[t], distProb[t], fermionProb[t]}, {t, outs}]],
 Frame -> All, Spacings -> {2, 1}]
```

This small table holds the whole physics. The distinguishable column is the flat classical expectation: a quarter, a half, a quarter. The boson column has lost its coincidence entirely, and that loss is destructive interference, even though the permanent adds all its terms with a plus sign. The point is subtle and worth saying plainly: the exchange coefficients are positive, but the optical amplitudes they multiply carry signs, and here the two routes to the coincidence outcome arrive as $+\tfrac12$ and $-\tfrac12$ and cancel. The fermion column is the sharpest: it is exactly the complementary pattern, one particle on each side with certainty. And notice why the bunched outcomes vanish for fermions without any extra rule. When two fermions would occupy one mode, the determinant selects two identical rows of $W$, and a determinant with a repeated row is zero. The Pauli exclusion principle is not imposed here; it falls out of the minus sign.

So a single beam splitter already shows that identical particles produce correlated counts that independent particles do not, and that one sign under exchange decides the shape of those correlations. The next question is whether this survives scaling up, and whether the permanent is really the object doing the work when the network is large and generic.

## From One Photon to Many: The Amplitude Is a Permanent

The rule we used by hand for the beam splitter is not special to two modes. It follows from how a passive network acts on the creation operators. Each network sends the operator that creates a photon in input port $k$ to a superposition over the outputs,

$$ a_k^\dagger \;\longmapsto\; \sum_j W_{jk}\, a_j^\dagger . $$

In other words, $W$ is the single-photon network, and a many-photon input state is built by applying several of these creation operators to the vacuum. Expanding the input Fock state $|S\rangle$, applying the substitution to every creation operator, and collecting the terms that land on a given output occupation $T$ produces

$$ A(T \leftarrow S) \;=\; \frac{\operatorname{Per}(W_{T,S})}{\sqrt{\prod_k s_k!\,\prod_j t_j!}}, \qquad P(T\mid S) = |A(T \leftarrow S)|^2, $$

where $W_{T,S}$ repeats row $j$ once per photon leaving mode $j$ and column $k$ once per photon entering mode $k$. The permanent and the factorials are consequences of the physics, not rules supplied by the software: the permanent is the sum over photon routings, and the factorials count the indistinguishable orderings when photons share a mode.

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

Read its transfer matrix, the five-mode analogue of the $2\times 2$ block above:

```wl
W = qi["TransferMatrix"];
W // MatrixForm
```

Now send three photons in, one into each of the first three ports, and ask for the amplitude of one specific pattern, all three coming out again one per port. By the formula this single number should be the permanent of the top-left $3\times 3$ block of $W$, since that input and output are collision-free and the factorials are all one.

Compute that amplitude directly from `QuantumInterferometer`:

```wl
qi["Amplitude", {1, 1, 1, 0, 0} -> {1, 1, 1, 0, 0}]
```

Confirm it equals the permanent of the selected block:

```wl
qi["Amplitude", {1, 1, 1, 0, 0} -> {1, 1, 1, 0, 0}] === Permanent[W[[{1, 2, 3}, {1, 2, 3}]]]
```

This confirms the formula on one pattern, but it is a weak check on its own, because `"Amplitude"` computes the permanent internally, so both sides of the equality run essentially the same calculation. The honest test is to compare two genuinely different routes to the numbers across the whole distribution, which is the next section.

## Reading the Whole Distribution: Two Routes to the Same Permanents

A boson sampler does not hand you one amplitude; it produces whole detection patterns, drawn with the probabilities the physics assigns. Collect that entire distribution for our three-photon input.

Collect the outcome distribution for the input $|1,1,1,0,0\rangle$:

```wl
probs = qi["Probabilities", {1, 1, 1, 0, 0}];
```

Show the five most likely detection patterns:

```wl
TakeLargest[probs, 5]
```

There is one probability for every way to place three photons in five modes, which is $\binom{n+m-1}{n}$ patterns.

Confirm the distribution has exactly that many entries:

```wl
Length[probs] == Binomial[3 + 5 - 1, 3]
```

Confirm it is a genuine probability distribution, summing to one (`Chop` sets the tiny machine-precision residue to zero so the comparison tests the physics, not the rounding; see [Chop](https://reference.wolfram.com/language/ref/Chop.html)):

```wl
Chop[Total[probs] - 1] == 0
```

Here is the real check. `QuantumInterferometer` reaches these probabilities by one route: it expands the product of output operators as a polynomial and reads off coefficients. Our formula reaches them by a different route: the squared permanent of a selected block, divided by the occupation factorials. If both routes agree on every one of the patterns, the permanent formula is confirmed by something independent, not by itself.

Verify that every probability in the collected distribution equals the squared permanent of $W$ over its normalization, taking the largest disagreement across all patterns:

```wl
Max[Abs[(Lookup[probs, QuditName[#], 0] -
      Abs[Permanent[W[[modeList[#], {1, 2, 3}]]]]^2/(Times @@ (#!))) & /@
    FrobeniusSolve[{1, 1, 1, 1, 1}, 3]]] < 10^-12
```

As one can see, the two routes agree to machine precision across the whole distribution. The entire table of detection probabilities is a table of permanents, confirmed by a calculation that never formed a permanent on the first route and never expanded a polynomial on the second. So to imitate the sampler by computing its probabilities, a classical computer must evaluate a permanent for every pattern it wants. The next question is the one that decides the cost: how expensive is a single permanent?

## The Cost of One Amplitude: A Permanent Against Its Sign-Flipped Twin

The permanent and the determinant are built from the same sum over pairings; they differ only in the sign carried by the odd permutations. That one difference is the difference between the boson amplitude and the fermion amplitude, and it is also, as far as anyone knows, the difference between hard and easy. The determinant's alternating signs are exactly what let Gaussian elimination collapse the sum, so a determinant costs on the order of $n^3$. No such cancellation is known for the permanent, and the best known exact methods sum on the order of $2^n$ terms.

The clean way to see this is not to race two different-sized matrices, but to hand the same matrices to both functions and watch their costs diverge. Define a function that builds a random $n$-mode network and times its transfer matrix's determinant and permanent on the identical matrix (`RepeatedTiming` evaluates several times and returns a trimmed mean; `ClearSystemCache[]` keeps a warm cache from flattering the second call):

```wl
ClearAll[costs];
costs[n_] := Block[{w},
   SeedRandom[n];
   w = QuantumInterferometer["Random"[n]]["Unitary"];
   ClearSystemCache[];
   {n, First @ RepeatedTiming[Det[w]], First @ RepeatedTiming[Permanent[w, Method -> "Glynn"]]}]
```

Time both from 2 through 20 modes (this cell takes about a minute):

```wl
timings = Table[costs[n], {n, 2, 20}];
```

Plot the two costs on the same logarithmic scale:

```wl
ListLogPlot[{timings[[All, {1, 2}]], timings[[All, {1, 3}]]},
 Joined -> True, PlotMarkers -> Automatic, Frame -> True, GridLines -> Automatic,
 AspectRatio -> 1/2, ImageSize -> Large,
 FrameLabel -> {"matrix size n", "seconds"},
 PlotLegends -> {"determinant (fermions)", "permanent (bosons)"},
 PlotLabel -> "Cost of one amplitude on the same n\[Times]n unitaries"]
```

As one can see, the permanent's points climb a straight line on a logarithmic scale, which means the time grows exponentially with the number of photons, while the determinant's stay almost flat across the whole range. Same matrices, two combinatorial functions, and only the exchange sign changed between them.

Measure the average factor by which each cost grows per added mode over the upper range, as the geometric mean of the consecutive ratios:

```wl
AssociationThread[{"permanent", "determinant"},
 GeometricMean[Ratios[Select[timings, 12 <= #[[1]] <= 20 &][[All, #]]]] & /@ {3, 2}]
```

The permanent's factor sits near two: each added photon roughly doubles the cost of one bosonic amplitude, exactly the $2^n$ growth the theory predicts. The determinant's factor sits near one: adding a mode barely moves the cost of the fermionic amplitude. So the plus sign is not just a label on a formula; it is the line between a cost that grows with every photon and a cost that essentially does not.

One honest caveat belongs right here, because it is easy to overclaim from a plot. A timing experiment shows that this particular algorithm is exponential; it cannot show that every possible algorithm must be. The determinant is itself a sum over the same enormous set of pairings, yet it is cheap, so "many terms" is not what makes a computation hard. What is true is sharper and more modest: no polynomial-time exact algorithm for the general permanent is known, and computing the permanent of a general matrix is #P-hard (Valiant, 1979). That is a statement about structure, not about the size of a sum.

## From Hard Probabilities to Hard Sampling: What the Computation Does and Does Not Show

Here is the trap the title invites, and the place to be most careful. We have shown that each probability is a permanent and that one permanent is expensive. It is tempting to conclude immediately that no classical computer can keep up with the sampler. But producing a sample is not the same task as computing a probability, and the gap between them is not a technicality.

The cleanest way to feel the gap is to find a case where the probabilities are still permanents and yet sampling is trivially easy. Distinguishable photons are exactly that case. Their output probabilities are permanents of the entrywise $|W|^2$, which we used in the beam-splitter table. Yet nobody needs a permanent to sample them, because each photon moves independently: a photon entering mode $k$ leaves with the probability distribution given by column $k$ of $|W|^2$, and you simply draw each photon's output on its own and count.

Confirm that each column of $|W|^2$ is a genuine probability distribution over the output modes, by checking that every column sums to one:

```wl
Chop[Total[Abs[W]^2] - 1]
```

As one can see, every column sums to one, so the distinguishable sampler is just independent draws, one per photon, no permanent required. This is the counterexample that forbids the easy inference: "the probabilities are permanents" cannot by itself make sampling hard, because here it does not.

So what does make boson sampling hard? The missing step is a complexity-theoretic one, and it is the part of this story that I can cite but cannot earn by running a cell. Aaronson and Arkhipov showed that an efficient classical algorithm sampling the exact boson-sampling distribution would let one estimate permanents well enough to collapse the polynomial hierarchy, which is believed impossible; the stronger and more realistic case of approximate sampling rests further on two conjectures, that estimating the permanent of a Gaussian random matrix is hard on average and that the output probabilities anticoncentrate ([Aaronson and Arkhipov, 2010](https://arxiv.org/abs/1011.3245); [a recent account of the remaining assumptions](https://arxiv.org/abs/2411.04566)). The hardness also lives in a specific regime, many more modes than photons, so that photons rarely collide and the relevant amplitudes are permanents of Haar-random submatrices, which is why the collision-free sampler, not the $m=n$ full matrix, is the object of the theorems. And even there, "hard" means a conjectured exponential scaling, not an impossibility: the best known classical samplers are exponential but markedly cheaper than computing one amplitude per pattern ([Clifford and Clifford, 2020](https://arxiv.org/abs/2005.04214)).

That is the honest shape of the argument, and it is worth stating as one sentence because it gives the computation-first part a precise job: indistinguishable-photon interference makes each detection probability a normalized permanent, so computing a single probability is already hard in the exact sense; a separate, conjecture-dependent complexity result turns that into strong evidence that reproducing the device's samples is classically hard too.

## Where This Leaves Us (and What Comes Next)

You now have a computation-first path through boson sampling that keeps its promises: a balanced beam splitter whose coincidence outcome vanishes by interference, a table in which one sign under exchange turns the correlated boson counts into the complementary fermion counts and the determinant enforces Pauli exclusion for free, the general rule that a detection amplitude is a normalized permanent of the transfer matrix, verified on a random five-mode network by two independent routes that agree to machine precision, a timing of the same matrices through the permanent and the determinant that shows the exchange sign separating an exponential cost from a flat one, and a careful account of why none of this by itself proves sampling is hard, with the complexity-theoretic bridge cited rather than claimed.

From here the natural next steps are the ones that separate a theorem from a laboratory. What happens when the photons are not perfectly identical, or when some are lost, since a real experiment is never ideal, and how much interference, and therefore how much of the hardness, survives. And what changes when single photons are replaced by squeezed light, the setting called Gaussian boson sampling, where the permanent gives way to a close relative, the hafnian. The permanent is why the exact probabilities are hard; what remains is the more delicate question of when the samples are hard too.
