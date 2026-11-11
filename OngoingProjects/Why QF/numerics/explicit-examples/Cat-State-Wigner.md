---
Template: Default
Title: "A Schrodinger-cat state and the negativity of its Wigner function"
---

# A Schrodinger-cat state and the negativity of its Wigner function

A coherent state $\lvert\alpha\rangle$ is the most classical state a single mode of light can be in: its Wigner function is a Gaussian bump sitting at the point $(x,p)$ that the classical amplitude $\alpha$ marks out, and it is non-negative everywhere, so it behaves like an ordinary probability density over phase space. A Schrodinger-cat state is the superposition of two such classical-like fields pulled far apart,

$$
\lvert\psi\rangle=\mathcal N\left(\lvert\alpha\rangle+e^{i\phi}\lvert-\alpha\rangle\right),
$$

with $\mathcal N$ fixing the norm and $\phi$ the relative phase. When $\lvert\alpha\rvert$ is large the two pieces are macroscopically distinguishable, which is the "dead and alive" picture Schrodinger used to needle the Copenhagen interpretation.

The question this example answers is what happens to the Wigner function in between the two lobes. The two coherent pieces are classical on their own, but their superposition is not, and phase space records the difference in one sharp way: the Wigner function develops oscillating fringes between the lobes, and those fringes dip below zero. A negative Wigner function has no classical probability interpretation, so it is the cleanest signature that the state is nonclassical. We build the state in the Wolfram QuantumFramework (QF), read its photon-number structure, compute the Wigner function, and watch the negativity appear exactly where the two branches interfere.

Throughout, $\lvert0\rangle=\{1,0,0,\dots\}$ is the vacuum and the Fock states $\lvert n\rangle$ carry $n$ photons; we truncate the infinite ladder at a finite cutoff, which the last section checks.

## Building the cat: a superposition of two coherent states

Setting $\phi=0$ makes an **even cat**, $\phi=\pi$ an **odd cat**, and $\phi=\pi/2$ the Yurke-Stoler state. Load the checked-out QF paclet:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"]
```

Load the main context:

```wl
Needs["Wolfram`QuantumFramework`"]
```

Expose the bosonic (second-quantization) context, which holds the Fock-space states and the phase-space functions:

```wl
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

The bosonic states live on a truncated Fock space of size `$FockSize`; set a cutoff of 40 levels, comfortably above the $\bar n=\lvert\alpha\rvert^2\approx 8$ photons we will put in:

```wl
SetFockSpaceSize[40]
```

`CatState[size][α, φ]` builds the normalized superposition on a space of the given size; with `CatState[][α, φ]` the size defaults to `$FockSize`. Build an even cat with a real amplitude $\alpha=2.8$:

```wl
even = CatState[][2.8, 0]
```

And the odd cat, with $\phi=\pi$:

```wl
odd = CatState[][2.8, Pi]
```

Both carry about eight photons and sit as two lobes on the $x$-axis at $x=\pm\sqrt2\,\alpha$, since a real $\alpha$ displaces along position only.

## Even and odd photon numbers: the interference in the number basis

The superposition is already visible before we reach phase space. Expanding $\lvert\pm\alpha\rangle$ in Fock states, the $\lvert-\alpha\rangle$ amplitudes carry a factor $(-1)^n$, so adding the two coherent states ($\phi=0$) cancels every odd $n$ and doubles every even $n$, while subtracting them ($\phi=\pi$) does the reverse. Read the first ten photon-number probabilities of the even cat; `Chop` sets the tiny machine-precision residue to zero so the comb is readable:

```wl
Chop[even["ProbabilitiesList"][[1 ;; 10]]]
```

Only even photon numbers survive. Do the same for the odd cat:

```wl
Chop[odd["ProbabilitiesList"][[1 ;; 10]]]
```

As one can see, the two combs are interleaved: the even cat lives on $n=0,2,4,\dots$ and the odd cat on $n=1,3,5,\dots$. That strict alternation is the interference between the two coherent branches, written in the photon-number basis.
The same even-or-odd structure makes a cat an eigenstate of $\hat a^2$: removing two photons maps even to even and odd to odd, and returns the same state. Build the annihilation operator:

```wl
a = AnnihilationOperator[]
```

Confirm that $\hat a^2$ sends the even cat to a multiple of itself, by checking that the normalized result is the original state:

```wl
QuantumSimilarity[(a[a][even])["Normalize"], even]
```

The similarity is one, so $\hat a^2\lvert\psi\rangle\propto\lvert\psi\rangle$. Read off the eigenvalue as the expectation $\langle\psi\rvert\hat a^2\lvert\psi\rangle$:

```wl
(SuperDagger[even] @ a @ a @ even)["Scalar"]
```

The eigenvalue is $\alpha^2$: both coherent pieces $\lvert\pm\alpha\rangle$ are eigenstates of $\hat a$ with eigenvalues $\pm\alpha$, so each is an eigenstate of $\hat a^2$ with the common eigenvalue $\alpha^2$, and so is any superposition of them.

## The Wigner function: two lobes, fringes, and a dip below zero

`WignerRepresentation[state, {xmin,xmax}, {ymin,ymax}]` returns the Wigner function as an `InterpolatingFunction` of $(x,p)$ over the requested window. Compute it for the even cat:

```wl
wig = WignerRepresentation[even, {-5, 5}, {-5, 5}]
```

Plot it over the phase plane:

```wl
Plot3D[wig[x, p], {x, -5, 5}, {p, -5, 5},
  PlotRange -> All, ColorFunction -> "TemperatureMap",
  PlotPoints -> 80, Mesh -> 15, MeshStyle -> Opacity[.2],
  AxesLabel -> {"x", "p", "W(x,p)"}, LabelStyle -> 13,
  PlotLabel -> "Wigner function of the even cat state (\[Alpha]=2.8)"]
```
Two Gaussian lobes sit at $x=\pm\sqrt2\,\alpha$, the phase-space images of $\lvert\alpha\rangle$ and $\lvert-\alpha\rangle$. Between them runs a ridge of fringes that oscillate along $p$, with no classical counterpart: they are the interference of the two branches. The smallest value on the grid is the deepest fringe:

```wl
Min[Flatten[wig["ValuesOnGrid"]]]
```

It is negative. Two separated Gaussians alone would make a non-negative function; the negativity is created entirely by the cross term between $\lvert\alpha\rangle$ and $\lvert-\alpha\rangle$, and it lives in the region between them, not on the lobes.

## The origin value and the photon-number parity

The deepest and the highest points of the fringe pattern both sit at the origin, halfway between the lobes, and the value there is fixed by a symmetry. The Wigner function at the origin is proportional to the expectation of the photon-number parity $\hat\Pi=(-1)^{\hat a^\dagger\hat a}$, which assigns $+1$ to even Fock states and $-1$ to odd ones. The even cat has only even photons and the odd cat only odd ones, so each is an eigenstate of parity:

```wl
Chop[Tr[DiagonalMatrix[(-1)^Range[0, $FockSize - 1]] . #["DensityMatrix"]] & /@ {even, odd}]
```

The even cat is an even-parity eigenstate and the odd cat an odd-parity one, opposite eigenvalues of $\hat\Pi$, so their Wigner functions sit at opposite extremes of the origin value. Evaluate the even cat's representation at the origin:

```wl
wig[0, 0]
```

This is $+1/\pi$, the largest value a Wigner function can take in this convention: `WignerRepresentation` normalizes the vacuum to $W(0,0)=1/\pi$ and bounds $\lvert W\rvert\le1/\pi$, so a parity of $+1$ pins the even cat to the top of that range and a parity of $-1$ pins the odd cat to $-1/\pi$.

A note on conventions, because QF offers a second route that scales differently. The closed-form point-value function `WignerFunction[ρ, α]` works in the complex amplitude $\alpha=(x+ip)/\sqrt2$ and normalizes the vacuum to $2/\pi$. Its origin value for the even cat is:

```wl
WignerFunction[even, 0]
```

This is $+2/\pi$, and the odd cat gives $-2/\pi$. The factor of two against `wig[0, 0]` is just the Jacobian between $d^2\alpha$ and $dx\,dp$; the physics, parity times the vacuum peak, is the same in both.

The contrast with a classical state is immediate. Build the coherent state with the same amplitude:

```wl
coh = CoherentState[][2.8]
```

Compute its Wigner function:

```wl
wigCoh = WignerRepresentation[coh, {-5, 5}, {-5, 5}]
```

Its smallest value on the grid is:

```wl
Min[Flatten[wigCoh["ValuesOnGrid"]]]
```

Zero to machine precision: the coherent Wigner function is a single non-negative Gaussian, with none of the fringes. Placing the two side by side shows what the superposition adds:

```wl
With[{cf = (Blend[{Blue, White, Red}, Rescale[#, {-1/Pi, 1/Pi}]] &)},
  GraphicsRow[
    MapThread[
      DensityPlot[#1[x, p], {x, -5, 5}, {p, -5, 5}, PlotPoints -> 120,
        ColorFunction -> cf, ColorFunctionScaling -> False, PlotRange -> {-1/Pi, 1/Pi},
        FrameLabel -> {"x", "p"}, PlotLabel -> #2,
        PlotLegends -> BarLegend[{cf, {-1/Pi, 1/Pi}}, LegendLabel -> "W"]] &,
      {{wig, wigCoh}, {"even cat: fringes dip below zero", "coherent state: non-negative Gaussian"}}],
    Spacings -> 20]]
```
## Checking the truncation

Everything above lives on a Fock space cut at 40 levels. A cat with $\bar n\approx8$ photons has negligible population that high, but confirm it rather than assume it. The total weight above level 30 is:

```wl
Total[even["ProbabilitiesList"][[31 ;;]]]
```

It is vanishingly small, so the truncation does not distort the state, the fringes, or the origin value; raising `$FockSize` leaves them unchanged.

## What the Wigner function shows

A Schrodinger-cat state is a superposition of two coherent fields, and that superposition leaves two marks. In the photon-number basis it is a strict even-or-odd comb, which makes the cat an eigenstate of $\hat a^2$ with eigenvalue $\alpha^2$. In phase space it is a ridge of fringes between the two coherent lobes, and those fringes go negative. The negativity is not an artifact of the lobes, which are classical Gaussians on their own; it is the interference between them, and its value at the origin is locked to the photon-number parity, $+1/\pi$ for the even cat and $-1/\pi$ for the odd one. A coherent state of the same amplitude shows none of this: its Wigner function is a single non-negative Gaussian. The dip below zero is the phase-space fingerprint of the superposition, and it is the reason a cat state is a resource no classical field can supply.
