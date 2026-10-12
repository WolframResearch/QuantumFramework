---
Template: Default
Title: "A nonclassical steady state from balanced drive and loss: the driven-dissipative Kerr oscillator"
---

# A nonclassical steady state from balanced drive and loss: the driven-dissipative Kerr oscillator

A single cavity mode holds light. Left alone it is a harmonic oscillator, and every photon in it costs the same energy. Add a Kerr medium and that stops being true: the energy to add one more photon now depends on how many are already there, so the cavity's resonance shifts with its own intensity. Drive the mode with a laser and photons flow in; let the mode leak through an imperfect mirror and photons flow out. Drive and loss do not cancel into nothing. They settle into a steady state, a density matrix that no longer changes in time, and the nonlinearity makes that steady state nonclassical.

Two pieces of this model are standard and appear separately in the quantum-optics literature. The Kerr term $\chi(\hat a^\dagger\hat a)^2$ generates unitary evolution of a coherent state: its photon-number components dephase and rephase, with no drive and no loss. A coherently driven, leaking cavity with a linear Hamiltonian $\hat H=\hbar G(\hat a+\hat a^\dagger)$ relaxes to a coherent state $\lvert\beta\rangle$ with $\beta=-2iG/\kappa$, a steady state that is still as classical as the drive that made it. This essay puts the nonlinearity, the drive, and the loss together, and asks what the combination relaxes to.

The model has every ingredient of those two, so the question is what changes when they act at once. The answer, computed below: the steady state is unique, it is approached slowly on a timescale set by the slowest decay mode of the dynamics, and its photon-number statistics are sub-Poissonian, a feature no coherent state has. Only numerics settle the first two, and the third is the signature of the Kerr term surviving into the steady state.

## The model: a Kerr cavity that is driven and leaks

Work in the frame rotating at the drive frequency, so the free evolution is slow and the drive is time independent. The cavity mode has annihilation operator $\hat a$ with $\hat a\lvert n\rangle=\sqrt n\lvert n-1\rangle$, and we retain Fock states $\lvert0\rangle,\ldots,\lvert d-1\rangle$. The Hamiltonian is

$$
\hat H=-\Delta\,\hat a^\dagger\hat a+\frac{K}{2}\,\hat a^\dagger\hat a^\dagger\hat a\hat a+F\left(\hat a+\hat a^\dagger\right).
$$

The first term is the detuning $\Delta=\omega_{\text{drive}}-\omega_{\text{cavity}}$ between the drive and the bare cavity. The second is the Kerr term in normal-ordered form; $\hat a^\dagger\hat a^\dagger\hat a\hat a=\hat n(\hat n-1)$ costs $K$ for each pair of photons, so the energy to add the $n$-th photon grows with $n$. This is the same nonlinearity as the unitary $\chi(\hat a^\dagger\hat a)^2$: the two differ by $\chi\,\hat a^\dagger\hat a$, a term proportional to $\hat n$ that only renames the detuning. The third term is the coherent drive of strength $F$. Photon loss is not in $\hat H$: it is a jump operator $\hat a$ with rate $\kappa$, coupling the mode to a zero-temperature bath, so the open evolution obeys the Lindblad equation

$$
\frac{d\hat\rho}{dt}=-i[\hat H,\hat\rho]+\kappa\left(\hat a\hat\rho\hat a^\dagger-\tfrac12\hat a^\dagger\hat a\hat\rho-\tfrac12\hat\rho\hat a^\dagger\hat a\right).
$$

Load the checked-out QF paclet:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"]
```

Load the main context:

```wl
Needs["Wolfram`QuantumFramework`"]
```

Expose the bosonic constructors:

```wl
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

Fix the parameters: a detuned, moderately nonlinear cavity driven harder than it leaks, truncated to 24 Fock levels. The convergence sweep at the end checks the cutoff:

```wl
{detuning, kerr, drive, kappa, d} = {0.5, 0.5, 1.0, 0.2, 24};
```

Build photon annihilation on the truncated mode, and its adjoint:

```wl
a = AnnihilationOperator[d]; ad = a["Dagger"];
```

Assemble the Hamiltonian from those operators; QF's `@` composes operators, so `ad @ a` is $\hat a^\dagger\hat a$:

```wl
H = -detuning (ad @ a) + kerr/2 (ad @ ad @ a @ a) + drive (a + ad)
```

Confirm that $\hat H$ is Hermitian:

```wl
H["HermitianQ"]
```

## Letting it settle: evolving from the vacuum

Start the cavity empty and watch it fill. `QuantumEvolve` takes the Hamiltonian, the jump operator with its rate as a rule, the initial state, and the time span, and integrates the Lindblad equation directly. Evolve the vacuum out to a time long compared with the photon lifetime $1/\kappa=5$:

```wl
solution = QuantumEvolve[H, {a} -> {kappa}, FockState[0, d], {t, 0, 30}];
```

The mean photon number is $\langle\hat n\rangle=\operatorname{Tr}(\hat\rho\,\hat a^\dagger\hat a)$. Define it on a density-matrix state, reading the trace against the number operator's matrix:

```wl
number = ad @ a; numberMatrix = Normal[number["Matrix"]];
meanPhoton[state_] := Re @ Tr[Normal[state["DensityMatrix"]] . numberMatrix]
```

The evolution result is callable in time; evaluate the mean photon number on a grid and plot its approach:

```wl
times = N[Range[0, 120]/4];
ListLinePlot[
    Transpose[{times, meanPhoton[solution[#]] & /@ times}],
    Frame -> True, FrameLabel -> {"t", "\[LeftAngleBracket]n\[RightAngleBracket]"},
    PlotLabel -> "Mean photon number relaxing toward steady state",
    GridLines -> Automatic, PlotRange -> All, ImageSize -> 640, AspectRatio -> 1/3
]
```

As one can see, the photon number does not climb smoothly to a plateau: it overshoots, dips, and rings down toward a limit, and even at $t=30$, six photon lifetimes, it is still visibly drifting. Reading the asymptote off a trajectory is slow and never exact. The next section gets the limit without evolving at all.

## The steady state without waiting: the null space of the Liouvillian

At steady state $d\hat\rho/dt=0$. Written with the Lindblad super-operator $\mathcal L$ as $\dot{\hat\rho}=\mathcal L[\hat\rho]$, that is the eigenvalue equation $\mathcal L[\hat\rho_{\text{ss}}]=0$: the steady state is whatever $\mathcal L$ annihilates, its null space. In other words, no time integration is needed. We vectorize the density matrix, build $\mathcal L$ as an ordinary matrix, and read off its kernel.

Build the Lindblad super-operator for this model; QF's `"Liouvillian"` takes the Hamiltonian, the jump operators, and their rates:

```wl
liouvillian = QuantumOperator["Liouvillian"[H, {a}, {kappa}]];
```

Confirm it is a super-operator, an object that acts on density matrices rather than on states:

```wl
liouvillian["Type"]
```

Find its null space from the matrix representation. The number of basis vectors is the number of independent steady states:

```wl
nullVectors = NullSpace[N @ Normal @ liouvillian["Matrix"]]; Length[nullVectors]
```

As one can see, the null space is one dimensional: at these parameters the steady state is unique. A Kerr cavity can, at other drive and detuning, carry two long-lived branches at once, the bistability that underlies many dispersive readout and switching schemes; that is a statement about the dimension of this null space, and here it is one. No feature of the model announces this in advance; the computation settles it.

The null vector is the steady density matrix stacked into a column, up to an arbitrary complex scale. Reshape it to a $d\times d$ matrix, make it exactly Hermitian to fix that phase, and set its trace to one:

```wl
raw = Partition[nullVectors[[1]], d];
steady = QuantumState[(raw + ConjugateTranspose[raw])/Tr[raw + ConjugateTranspose[raw]], d];
```

Read its mean photon number:

```wl
meanPhoton[steady]
```

This is the horizontal limit the evolving trajectory was crawling toward, obtained in one linear-algebra step instead of a long integration. Confirm the two agree by overlaying the null-space value on the trajectory; the dashed line is the steady state:

```wl
ListLinePlot[
    Transpose[{times, meanPhoton[solution[#]] & /@ times}],
    Epilog -> {Dashed, Red, Line[{{0, meanPhoton[steady]}, {30, meanPhoton[steady]}}]},
    Frame -> True, FrameLabel -> {"t", "\[LeftAngleBracket]n\[RightAngleBracket]"},
    PlotLabel -> "Evolution (blue) approaching the null-space steady state (red)",
    GridLines -> Automatic, PlotRange -> All, ImageSize -> 640, AspectRatio -> 1/3
]
```

As expected, the ringing trajectory is settling onto the dashed line but has not reached it. How slowly it does so is set by $\mathcal L$ itself.

## How long is long: the Liouvillian gap

The eigenvalues of $\mathcal L$ have non-positive real parts; the zero eigenvalue is the steady state, and every other mode decays as $e^{-|\mathrm{Re}\,\lambda|t}$. The slowest of these, the smallest nonzero decay rate, is the Liouvillian gap, and its reciprocal is the time the slowest transient takes to die. Collect the decay rates and take the smallest nonzero one:

```wl
decayRates = Sort[Select[-Re[Eigenvalues[N @ Normal @ liouvillian["Matrix"]]], # > 10^-8 &]];
gap = decayRates[[1]]
```

The reciprocal of the gap is the time the slowest transient takes to die:

```wl
1/gap
```

In short, this relaxation time is a sizeable fraction of the window the evolution ran out to, $t=30$, so a transient that needs several such times to vanish is still alive at the end of it. This is why reading the steady state off the trajectory was so slow, and why the null space is the right tool: it does not care how small the gap is.

## Is the cutoff big enough: a convergence sweep

Every number above lives in a 24-level truncation of an infinite Fock space. Rather than assert that 24 is enough, recompute the steady-state photon number at a range of cutoffs and watch it stop moving. Package the null-space construction as a function of the cutoff:

```wl
steadyMeanPhoton[cutoff_] := Module[{b = AnnihilationOperator[cutoff], bd, ham, vec, rho},
    bd = b["Dagger"];
    ham = -detuning (bd @ b) + kerr/2 (bd @ bd @ b @ b) + drive (b + bd);
    vec = NullSpace[N @ Normal @ QuantumOperator["Liouvillian"[ham, {b}, {kappa}]]["Matrix"]][[1]];
    rho = Partition[vec, cutoff]; rho = (rho + ConjugateTranspose[rho])/Tr[rho + ConjugateTranspose[rho]];
    Re @ Tr[rho . Normal[(bd @ b)["Matrix"]]]
]
```

Sweep the cutoff from a box too small to hold the state up to the working value:

```wl
AssociationThread[#, steadyMeanPhoton /@ #] &[{8, 12, 16, 20, 24}]
```

As one can see, the photon number is still moving at a cutoff of 8, where the box clips the tail of the distribution, but it has stopped changing by 12 and agrees digit for digit from there to 24. The truncation is not a free parameter of the physics; it is a numerical control, and the sweep shows where it stops mattering.

## A nonclassical balance: sub-Poissonian light

The steady state looks, in the plots above, like an ordinary blob of a few photons. What makes it nonclassical is hidden in the shape of its photon-number distribution. Read the populations $p_n=\langle n\rvert\hat\rho_{\text{ss}}\lvert n\rangle$:

```wl
populations = steady["ProbabilitiesList"];
```

A coherent state of the same mean has a Poisson distribution, whose variance equals its mean. The Mandel $Q$ parameter measures the departure from that, $Q=(\langle\Delta n^2\rangle-\langle n\rangle)/\langle n\rangle$, and is negative exactly when the distribution is narrower than Poisson, a property no classical field has. Compute it from the populations:

```wl
With[{n = Range[0, d - 1]}, meanN = n . populations; (n^2 . populations - meanN^2 - meanN)/meanN]
```

As one can see, $Q$ is negative: the steady state is sub-Poissonian. The Kerr term is responsible. Once a photon is in the cavity it shifts the resonance away from the drive, so the drive is less able to add the next one, and the photon number is squeezed below the Poissonian spread the drive alone would produce.

The same physics shows up in the second-order coherence $g^{(2)}(0)=\langle\hat a^\dagger\hat a^\dagger\hat a\hat a\rangle/\langle\hat a^\dagger\hat a\rangle^2=\langle n(n-1)\rangle/\langle n\rangle^2$, the chance of finding two photons at once relative to a coherent field. Compute it, again from the populations alone:

```wl
With[{n = Range[0, d - 1]}, ((n (n - 1)) . populations)/(n . populations)^2]
```

As expected, $g^{(2)}(0)<1$: the light is antibunched, photons arrive spread apart rather than clumped, which is the time-domain face of the same sub-Poissonian statistics. Plot the steady-state populations against the Poisson distribution of equal mean to see the narrowing directly:

```wl
With[{meanN = Range[0, d - 1] . populations},
    ListPlot[
        {Callout[Transpose[{Range[0, d - 1], populations}], "Kerr steady state"],
         Callout[Transpose[{Range[0, d - 1], PDF[PoissonDistribution[meanN], Range[0, d - 1]]}], "Poisson, same mean"]},
        Filling -> Axis, PlotRange -> {{0, 12}, All}, Frame -> True,
        FrameLabel -> {"n", "\!\(\*SubscriptBox[\(p\), \(n\)]\)"},
        PlotLabel -> "Steady-state photon statistics versus Poisson",
        ImageSize -> 560
    ]
]
```

The Kerr distribution sits higher at its center and falls off faster than the Poisson curve: less spread for the same mean. This nonclassicality is in the photon-number statistics; it need not show in phase space. Check the Wigner function of the steady state, whose most negative value is a witness of phase-space nonclassicality:

```wl
Min[Flatten[WignerRepresentation[steady, {-5., 5.}, {-5., 5.}, "GridSize" -> 60]["ValuesOnGrid"]]]
```

As one can see, the minimum is zero to machine precision: at these parameters the Wigner function stays non-negative, so the nonclassicality lives entirely in the number statistics, not in phase-space negativity. The two are different witnesses, and a state can carry one without the other.

To see that the Kerr term is what does this, remove it. With $K=0$ the Hamiltonian is linear and the model is the driven, leaking cavity whose steady state is a coherent state. Build that steady state the same way, with a loss rate chosen to keep its amplitude inside the same 24-level box:

```wl
linear = With[{vec = NullSpace[N @ Normal @ QuantumOperator["Liouvillian"[drive (a + ad), {a}, {1.0}]]["Matrix"]][[1]]},
    Module[{rho = Partition[vec, d]}, QuantumState[(rho + ConjugateTranspose[rho])/Tr[rho + ConjugateTranspose[rho]], d]]
];
```

Its Mandel $Q$ is zero to rounding, Poissonian:

```wl
With[{n = Range[0, d - 1], p = linear["ProbabilitiesList"]}, meanN = n . p; (n^2 . p - meanN^2 - meanN)/meanN]
```

And it is exactly the coherent state $\lvert\beta\rangle$ with $\beta=-2F/\kappa\,i$, the amplitude predicted for the linear cavity; the distance between them is:

```wl
QuantumDistance[linear, CoherentState[d][-2. drive/1.0 I]]
```

As expected, the distance vanishes: with the nonlinearity off, drive and loss balance into a coherent state, as classical as the laser that drives it. The Kerr term is the whole difference between that and the sub-Poissonian state found above.

## Where this leaves us

The driven-dissipative Kerr cavity relaxes to a single steady state, and at the parameters here its photon-number distribution is narrower than the Poisson distribution of the same mean, so the light is sub-Poissonian and antibunched. That nonclassicality is the Kerr nonlinearity surviving the balance of drive and loss: the intensity-dependent resonance resists adding each next photon, and removing the Kerr term returns the coherent, classical steady state of the linear cavity. The steady state is reached slowly, on the timescale set by the smallest nonzero decay rate of the Liouvillian, which makes evolving in time an unreliable way to find it and the null space of the super-operator the direct one. Both the uniqueness of the steady state and the slowness of its approach are properties of the Liouvillian spectrum that only the computation reveals; neither is visible in the Hamiltonian alone.

<!-- Smoke-tested cell by cell in Wolfram Language with QF loaded from the checked-out paclet, 2026-10-12. The .md source is the deliverable; no .nb was built. -->
