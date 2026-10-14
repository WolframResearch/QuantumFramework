---
Template: Default
Title: Why Boson Sampling Is Hard
Author: Mads Bahrami
---

# Why Boson Sampling Is Hard

<!-- #| style: Subtitle -->
Identical Photons, Permanents, and the One Sign That Separates Easy from Hard: A Short, Computation-First Argument with the Wolfram Language and QuantumFramework

<!-- #| style: Author -->
Mads Bahrami (last updated: Oct 12, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

### Setting the Stage: How This Notebook Flows

This notebook makes one argument, and only one: when identical photons pass through an optical network, the probability of each detection pattern is a permanent of a matrix, and a permanent is expensive to compute, so a classical computer cannot keep up with the photons. We will send a few photons into a random interferometer, read the distribution of outcomes straight from `QuantumInterferometer`, show that each probability is a permanent, time that permanent as photons are added, and finish with the one small difference (a sign) that makes the fermionic cousin of the same problem easy.

In other words, I want to earn the claim "boson sampling is hard" by computing it, not by repeating it. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. Everything below is run, not asserted.

Before we start, a few things. The environment you see is a [Wolfram (Mathematica) notebook](https://www.wolfram.com/notebooks/). Evaluate the cells from top to bottom; a few later cells use a network and a distribution defined earlier, so order matters. Two cells in the second half time a computation on purpose and take up to about a minute; only their trend matters, not the absolute seconds on your machine.

The story runs as a continuous sequence, like a movie, with a few headings only as transitions. My suggestion is to read each output and its meaning first, and worry about how the input cell is written afterward. And remember you are not locked into the code: change the number of photons, the number of modes, the seed, and watch the argument bend or hold.

If you have any suggestions or questions, please reach out to us at [quantum@wolfram.com](mailto:quantum@wolfram.com)

### Prerequisites and How to Read This

You only need to be comfortable with the idea of photons occupying modes: a state $|n_1,n_2,\ldots\rangle$ says how many photons sit in each mode, and a detector reports one such pattern. Everything optical is handled by `QuantumInterferometer`, a QuantumFramework tool that builds a passive optical network and hands back its outcome probabilities directly. The combinatorial objects at the end, the permanent and the determinant, are built into the Wolfram Language.

Let's start!

## A Boson Sampler: Photons Into a Random Network

A passive optical network is a box of mirrors, beam splitters and phase shifters: it moves photons between modes without creating, destroying or interacting them. A boson sampler feeds $n$ identical photons into such a box, one per input port, and lets detectors at the output report which modes the photons left in. In other words, the machine's whole output is a list of detection patterns, drawn with whatever probabilities the physics assigns.

Load QuantumFramework together with its second-quantization sub-package, where the optical tools live:

```wl
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

The hardness argument wants a network with no special structure, so we take one at random, a unitary drawn uniformly over all optical networks of its size. `QuantumInterferometer["Random"[m]]` builds such an $m$-mode network and factors it into beam splitters and phase shifters for us.

Build a random five-mode interferometer with a fixed seed:

```wl
SeedRandom[2026];
qi = QuantumInterferometer["Random"[5]]
```

Draw the network, so we can see the mesh of beam splitters and phase shifters the photons will meet:

```wl
qi["Diagram"]
```

Now send three photons in, one into each of the first three ports, and ask for the probability of every detection pattern at the output. `QuantumInterferometer` returns the whole distribution at once, keyed by the output pattern in Dirac notation.

Collect the outcome distribution for the input $|1,1,1,0,0\rangle$:

```wl
probs = qi["Probabilities", {1, 1, 1, 0, 0}];
```

Show the five most likely detection patterns:

```wl
TakeLargest[probs, 5]
```

The sampler's job is to produce patterns with exactly these probabilities. There is one number for every way to place three photons in five modes, which is $\binom{n+m-1}{n}$ patterns.

Confirm that the distribution has one probability per pattern:

```wl
Length[probs] == Binomial[3 + 5 - 1, 3]
```

Confirm that it is a genuine probability distribution, summing to one (`Chop` sets the tiny machine-precision residue to zero so the comparison tests the physics, not the rounding):

```wl
Chop[Total[probs] - 1] == 0
```

So far this is just a machine producing a distribution. The real question is where these numbers come from, and that is where the permanent appears.

## Every Probability Is a Permanent

Everything the network does to a single photon is held in one $m \times m$ matrix $W$, its transfer matrix: the entry $W_{jk}$ is the amplitude for a photon that enters port $k$ to leave port $j$. In other words, $W$ is the network seen by one photon at a time.

Read the transfer matrix of our network:

```wl
W = qi["TransferMatrix"];
W // MatrixForm
```

With several identical photons, the amplitude for a detection pattern is not one entry of $W$ but a combination of many, one for every way of routing the photons to the detected modes. That combination is the permanent. The permanent of a square matrix is defined like the determinant, as a sum over all ways of pairing rows with columns, but with every term added with a plus sign (the determinant alternates the signs). Build the relevant submatrix by taking the rows and columns of $W$ for the occupied output and input modes.

Compute, directly from `QuantumInterferometer`, the amplitude for all three photons to come out one per port again, in modes 1, 2 and 3:

```wl
qi["Amplitude", {1, 1, 1, 0, 0} -> {1, 1, 1, 0, 0}]
```

That single number is the permanent of the three-by-three block of $W$ on the occupied modes.

Confirm that the amplitude equals that permanent exactly:

```wl
qi["Amplitude", {1, 1, 1, 0, 0} -> {1, 1, 1, 0, 0}] === Permanent[W[[{1, 2, 3}, {1, 2, 3}]]]
```

This is not special to that one pattern. For a general input $S$ and output $T$, the rule is the same: repeat row $j$ of $W$ once for every photon leaving mode $j$, repeat column $k$ once for every photon entering mode $k$, take the permanent of that matrix, and divide by a normalization that counts photons sharing a mode.

Define the list that repeats mode $j$ once for every photon in it, so a pattern becomes the rows (or columns) of its matrix:

```wl
ClearAll[modeList]
modeList[occ_List] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occ]
```

Verify that every probability in the distribution we collected is the squared permanent of $W$, divided by that normalization:

```wl
Max[Abs[(Lookup[probs, Ket[#], 0] -
      Abs[Permanent[W[[modeList[#], {1, 2, 3}]]]]^2/(Times @@ (#!))) & /@
    FrobeniusSolve[{1, 1, 1, 1, 1}, 3]]] < 10^-12
```

As one can see, the whole distribution is a table of permanents. `QuantumInterferometer` computes it exactly this way, straight from the transfer matrix, without ever building and evolving a giant list of every way the photons could be arranged. So to imitate the sampler, a classical computer must evaluate a permanent for every pattern it wants a probability for. The next question is the only one that matters: how expensive is one permanent?

## The Cost: One Permanent Doubles With Every Photon

Take the hardest honest case, $n$ photons entering $n$ ports of an $n$-mode network, so the matrix behind one amplitude is the full $n \times n$ transfer matrix. The permanent has no shortcut like the one the determinant enjoys: the best known exact methods add up a number of terms that grows like $2^n$. In other words, every photon added roughly doubles the work for a single amplitude.

Define a function that builds a random $n$-mode network and times the permanent of its transfer matrix, the matrix behind the amplitude of $n$ photons entering one per port (`RepeatedTiming` evaluates the permanent several times and returns a trimmed mean):

```wl
ClearAll[permanentTime]
permanentTime[n_] := Block[{w},
   SeedRandom[n];
   w = QuantumInterferometer["Random"[n]]["Unitary"];
   ClearSystemCache[];
   First @ RepeatedTiming[Permanent[w, Method -> "Glynn"]]]
```

Time it for 10 through 20 photons (this cell takes about a minute):

```wl
permanentTimes = Table[{n, permanentTime[n]}, {n, 10, 20}];
```

Plot the timings on a logarithmic scale:

```wl
ListLogPlot[permanentTimes, Joined -> True, PlotMarkers -> Automatic,
 Frame -> True, GridLines -> Automatic, AspectRatio -> 1/2, ImageSize -> Large,
 FrameLabel -> {"photons n", "seconds"},
 PlotLabel -> "Time to compute one n\[Times]n permanent"]
```

As one can see, the points fall on a straight line, and a straight line on a logarithmic scale means the time grows exponentially with the number of photons.

Compute the average factor by which the time grows with each added photon, as the geometric mean over the range:

```wl
(permanentTimes[[-1, 2]]/permanentTimes[[1, 2]])^(1/(Length[permanentTimes] - 1))
```

The factor sits near two: each added photon roughly doubles the cost of a single amplitude, exactly as the $2^n$ growth predicts. At that rate a few dozen photons put one amplitude beyond a laptop, and a hundred beyond any computer, even though the physics is nothing more than photons sliding between modes.

## Fermions Are Easy: The Sign Under Exchange

Here is the sharpest way to see that the permanent is the whole difficulty. Identical fermions, electrons rather than photons, obey the same linear optics, but their histories pick up a minus sign when two are exchanged. That single sign turns the permanent into a determinant. And the determinant, because its alternating signs allow elimination, is cheap: its cost grows like $n^3$, not $2^n$.

Time the determinant of a thousand-by-thousand random matrix, the fermionic amplitude for a thousand particles:

```wl
Block[{a = RandomComplex[{-1 - I, 1 + I}, {1000, 1000}]},
 ClearSystemCache[]; First @ AbsoluteTiming[Det[a]]]
```

As one can see, the amplitude of a thousand identical fermions costs less than the permanent of a couple of dozen photons we timed just above. The physics of the two particles is identical up to one sign under exchange, and that sign is the entire distance between a computation that finishes before you lift your hands from the keyboard and one that outlasts the universe.

## Where This Leaves Us

The argument closes on itself: a boson sampler produces detection patterns whose probabilities are permanents of its transfer matrix, `QuantumInterferometer` hands us those probabilities directly, each one is a permanent we verified by hand, the cost of a permanent doubles with every photon, and the fermionic version with its one extra sign collapses to an easy determinant. Boson sampling is hard for exactly one reason, and now you have computed it: the permanent has no shortcut, and the identical photons are computing one every time they land.

From here the natural next step is to ask what happens when the photons are not perfectly identical or when some are lost, since a real experiment is never ideal, and whether squeezed light in place of single photons (Gaussian boson sampling, where hafnians replace permanents) changes the story. But the core is already in your hands, and it is a single sentence: the permanent is why it is hard.
