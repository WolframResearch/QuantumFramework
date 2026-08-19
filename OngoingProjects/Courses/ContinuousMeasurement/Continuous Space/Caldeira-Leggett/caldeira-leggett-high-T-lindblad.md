---
Template: Default
Title: The High-Temperature Lindblad Limit of Caldeira-Leggett
Author: Mads Bahrami
---

# The High-Temperature Lindblad Limit of Caldeira-Leggett

<!-- #| style: Subtitle -->
Keeping one term past white noise turns the non-positive Caldeira-Leggett master equation into a guaranteed Lindblad generator; carrying the fluctuation-path expansion to all orders reveals its Bernoulli-number structure and why exactly two terms are special. Worked out and verified with the Wolfram Language.

<!-- #| style: Author -->
Mads Bahrami (last updated: August 17, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

<!-- #| style: Abstract -->
The textbook high-temperature limit of quantum Brownian motion collapses the Feynman-Vernon decoherence kernel to white noise, and the resulting Caldeira-Leggett master equation is famously not completely positive: its Kossakowski matrix carries a negative eigenvalue. Following Pleasance, Aurell and Petruccione (arXiv:2508.14262), this essay treats the kernel as a tempered distribution and keeps exactly one term past white noise, a $\delta''(\tau)$ piece whose weight is fixed by the single integral $\int_0^\infty x^2/\sinh^2 x\, dx = \pi^2/6$. That one term fills the momentum slot the standard treatment leaves empty, lifts the Kossakowski determinant from negative to positive, and does so with a generator that is completely independent of the ultraviolet cutoff $\Omega$. We build the kernel algebra from primitives, assemble the Kossakowski matrix and read its spectrum, then watch the correction rescue positivity dynamically: for a quadratic potential the master equation closes on the $2\times 2$ covariance, and a position-squeezed state that pierces the Robertson-Schrödinger uncertainty floor under Caldeira-Leggett stays above it once the new term is present. Carrying the fluctuation-path expansion to all orders, we find its coefficients are the Bernoulli numbers and that the whole tower resums to the exact Ohmic kernel; the first higher term, a $\delta''''(\tau)$, carries a potential-dependent operator that leaves the two-operator Lindblad family and is anti-diffusive because the next Bernoulli number is negative, yet for a harmonic potential it sharpens the agreement with the quantum Gibbs state from order $\hbar^2$ to order $\hbar^4$, and the expansion converges down to the first Matsubara frequency, where the bath memory turns genuinely non-Markovian.

### Setting the Stage: How This Notebook Flows

This notebook is a computation-first tour of one idea in open quantum systems: the correct Markovian limit of the Caldeira-Leggett model of quantum Brownian motion. We start from the Feynman-Vernon influence functional and its decoherence kernel, watch the standard high-temperature limit throw that kernel down to white noise, and then recover the one term it discards by treating the kernel as a tempered distribution rather than an ordinary function. From that single correction we assemble the Markovian generator, read its Kossakowski matrix to decide whether it is a legitimate Lindblad generator, and, for a harmonic potential, reduce the whole master equation to a drift-plus-diffusion flow of the covariance matrix and watch positivity be lost or saved in a single plot. The last stretch goes past the paper: it carries the fluctuation-path expansion to all orders, resums it into the exact kernel, and asks what the higher terms do to the operator content, to positivity, and to how far down in temperature the master equation can reach.

In other words, I have tried to build a catalogue of computational experiments around one paper's central claim, so that every abstract statement about positivity and cutoff-independence is something you can evaluate, perturb, and see for yourself. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. Here that conviction has teeth, because the entire dispute between the standard master equation and the corrected one is settled by the sign of a $2\times 2$ determinant and by whether one curve dips below zero.

The environment you see is a live Wolfram notebook. Evaluate the cells from top to bottom. I kept the symbolic parameters $\hbar$, $M$, $\omega_0$, $\gamma$, $k_BT$ genuinely symbolic throughout the algebra, and only substituted numbers where a plot demands them, so most cells return closed forms you can inspect. Some later cells reuse definitions from earlier ones, so be mindful of order.

The story is presented as a continuous sequence, like a movie. I have added headings to help with transitions, but the narrative runs straight through: concept, then computation, then interpretation. Some code will look dense at first. My suggestion is to read the output and its meaning before worrying about every detail of the input, then go back and unpack the code once you know what it produced. And remember that you are not locked into the code as given: change the temperature, change the squeezing, add the cross-diffusion term, and watch the verdict move.

### Prerequisites and How to Read This

This tour is self-contained but moves at the pace of a research paper. It assumes comfort with the density-operator description of open systems, the Lindblad (GKSL) form of a Markovian master equation, canonical commutators, and a little distribution theory (the Dirac delta and its derivatives as tempered distributions). The final part uses the covariance-matrix description of Gaussian states and the Lyapunov equation for its steady state; both are built up here from the moment equations, so no prior exposure is strictly required, though the pace there is steeper.

Let's start!

## Part I: The Caldeira-Leggett Model and the Feynman-Vernon Kernel

Fix the object of study. A quantum particle of mass $M$, position $\hat X$, momentum $\hat P$, moving in a potential $V(\hat X)$, is coupled bilinearly to a bath of harmonic oscillators held at temperature $T$ (write $\beta = 1/k_BT$). The total Hamiltonian is $\hat H = \hat H_S + \hat H_B + \hat H_I + \hat H_C$, with system Hamiltonian $\hat H_S = \hat P^2/2M + V(\hat X)$, bath $\hat H_B = \sum_k (\hat p_k^2/2m_k + \tfrac12 m_k\omega_k^2\hat x_k^2)$, bilinear coupling $\hat H_I = -\hat X\sum_k c_k\hat x_k$, and a counter-term $\hat H_C = \hat X^2\sum_k c_k^2/(2m_k\omega_k^2)$ that keeps the Hamiltonian translation-invariant. This is the Caldeira-Leggett model, the workhorse of quantum Brownian motion.

Integrating out the bath under a path-integral treatment leaves the reduced dynamics of the particle governed by the Feynman-Vernon influence functional $\mathcal{F}[X,X'] = \exp(-\tfrac{i}{\hbar}S_i - \tfrac{1}{\hbar}S_r)$, a double integral over a forward path $X$ and a backward path $X'$. The imaginary part $S_i$ is dissipation (friction); the real part $S_r$ is decoherence, and it is the one with memory. Writing the paths in terms of the fluctuation coordinate $\Delta X = X - X'$, the decoherence action is

$$
S_r[\Delta X] = \int_{t_i}^{t_f}\!dt\int_{t_i}^{t}\!ds\; k_r(t-s)\,\Delta X(t)\,\Delta X(s),
$$

so the entire memory of decoherence lives in one object, the real kernel $k_r(\tau)$. Equivalently, $k_r$ is the force-force correlation function of a classical Gaussian noise with which the bath is, for the purpose of decoherence, indistinguishable. Everything in this essay is a statement about this one function.

For an Ohmic bath with spectral density $J(\omega) = 2\eta\,\omega\,\theta(\Omega - \omega)$, cut off at the ultraviolet frequency $\Omega$, the kernel is the cosine transform of $\omega\operatorname{coth}(\hbar\omega/2k_BT)$:

$$
k_r(\tau) = \frac{2\eta}{\pi}\int_0^\Omega d\omega\;\omega\operatorname{coth}\!\Big(\frac{\hbar\omega}{2k_BT}\Big)\cos\omega\tau .
$$

Two scales sit inside this integral: the cutoff $\Omega$, which limits the bath's response time, and the thermal time $\tau_B = \hbar/k_BT$, which limits the range of the memory. The whole paper is about what happens to $k_r$ when $\tau_B$ is small but not zero.

Define the thermal weight and look at its small-argument behavior, since that is where the high-temperature limit lives:

```wl
Series[\[HBar] \[Omega] Coth[\[HBar] \[Omega]/(2 kB T)], {\[Omega], 0, 3}]
```

As one can see, the weight $\hbar\omega\operatorname{coth}(\hbar\omega/2k_BT)$ is $2k_BT$ at leading order, a constant independent of $\omega$, with the first correction quadratic in $\omega$. The leading constant is the classical equipartition value, and it is the term the standard high-temperature limit keeps. The quadratic correction is the term this essay is about.

## Part II: The Standard High-Temperature Limit: White Noise, and What It Discards

The textbook move is to send $\tau_B \to 0$ and keep only the leading term of the thermal weight. Concretely, replace $\operatorname{coth}(\hbar\beta\omega/2)$ by its small-argument value.

Expand the hyperbolic cotangent for small argument, which is the substitution the standard limit makes:

```wl
Series[Coth[x/2], {x, 0, 3}]
```

As one can see, $\operatorname{coth}(x/2) = 2/x + x/6 - \dots$. Keeping only the $2/x$ term is exactly the white-noise approximation $\operatorname{coth}(\hbar\beta\omega/2)\approx 2/\hbar\beta\omega$. Notice the very next term, $x/6$: it is small, it is usually thrown away, and it is the entire subject of the correction we will recover. The standard limit does not merely approximate the kernel; it discards a specific term whose fingerprint we can track.

With only the leading term, the kernel becomes $k_r(\tau)\approx (4\eta/\pi\tau_B)\int_0^\Omega\cos\omega\tau\, d\omega = (4\eta/\pi\tau_B)\sin(\Omega\tau)/\tau$. As the cutoff is removed, $\sin(\Omega\tau)/\tau$ is a nascent Dirac delta.

Confirm that $\sin(\Omega\tau)/\pi\tau$ carries unit area, so it converges to $\delta(\tau)$ as $\Omega\to\infty$:

```wl
Integrate[Sin[\[CapitalOmega] \[Tau]]/(Pi \[Tau]), {\[Tau], -Infinity, Infinity},
  Assumptions -> \[CapitalOmega] > 0]
```

As one can see, the area is exactly one, independent of $\Omega$, so the finite-cutoff kernel collapses to the singular white-noise kernel

$$
k_r^{(0)}(\tau) = \frac{4\eta}{\tau_B}\,\delta(\tau).
$$

This memoryless kernel is what produces the standard Caldeira-Leggett master equation. We will see in Part IV that it is precisely one term short of positivity. Recovering the missing term is a matter of not discarding the $x/6$ we just watched go by.

## Part III: The Kernel as a Tempered Distribution: The $\sinh^{-2}$ Form and the $\delta''$ Correction

The key methodological choice is to treat $k_r$ not as an ordinary function but as a tempered distribution, an object that only ever appears smeared against a smooth test function. The paper's central representation of the Ohmic kernel, valid up to corrections that vanish as $\Omega\to\infty$, is

$$
k_r(\tau) = \frac{2\eta}{\pi}\left[-\frac{\pi^2}{\tau_B^2}\,\frac{1}{\sinh^2(\pi\tau/\tau_B)} + \Omega^2\,O\!\Big(\tfrac{1}{(\Omega\tau)^3}\Big)\right],
$$

a hyperbolic-cosecant-squared profile whose width is set by the thermal time $\tau_B$, plus a remainder controlled entirely by the cutoff. Let us reconstruct the $\sinh^{-2}$ core from primitives, following Sec. I of the Supplemental Material, then extract the correction the white-noise limit dropped.

The kernel is a second derivative of a simpler thermal primitive. Compute the thermal part of that primitive, $G^{\rm th}(s) = \int_0^\infty (\operatorname{coth}\tfrac{x}{2} - 1)\sin sx\, dx$, where subtracting the $1$ isolates the thermal from the vacuum contribution:

```wl
Integrate[(Coth[x/2] - 1) Sin[s x], {x, 0, Infinity}, Assumptions -> s > 0]
```

As one can see, the thermal primitive is $G^{\rm th}(s) = \pi\operatorname{coth}\pi s - 1/s$, a closed form. The physics of the thermal weight has been packed into a single hyperbolic cotangent.

Differentiate the primitive once more to reach the kernel itself, since the Feynman-Vernon kernel is $dG/ds$:

```wl
D[Pi Coth[Pi s] - 1/s, s] // Simplify
```

As one can see, the derivative is $1/s^2 - \pi^2/\sinh^2(\pi s)$. The $1/s^2$ is a vacuum term that cancels against the vacuum part of the primitive; the surviving thermal piece is $-\pi^2/\sinh^2(\pi s)$, which is exactly the profile quoted above. So the kernel really is a $\sinh^{-2}$, and this is the function we must now integrate against the paths.

Here is the subtlety that makes the whole construction work, and the reason the kernel must be handled as a distribution. The memory kernel $k_r^{(m)}(\tau) = (\pi^2/\tau_B^2)/\sinh^2(\pi\tau/\tau_B)$ behaves like $1/\tau^2$ near the origin, so it is not integrable there. Cut the integral off at a lower limit $\epsilon$ and integrate the tail:

```wl
Integrate[Csch[Pi s]^2, {s, \[Epsilon], Infinity}, Assumptions -> \[Epsilon] > 0]
```

The tail integral is $(\operatorname{coth}\pi\epsilon - 1)/\pi$, finite for any $\epsilon > 0$. Now let the lower limit approach the origin:

```wl
Limit[(-1 + Coth[Pi \[Epsilon]])/Pi, \[Epsilon] -> 0, Direction -> "FromAbove"]
```

As one can see, it runs off to infinity: a naive numerical integration of the bare kernel would diverge here, and the object is only meaningful as a distribution.

The regularization is the one distribution theory prescribes: subtract the value of the test function at the origin. The decoherence functional, after the integration by parts of Sec. II of the Supplemental Material, reduces to $T_F[x] = 2\pi\big(x(0) - \pi\int_0^\infty ds\,\sinh^{-2}(\pi s)\,(x(s) - x(0))\big)$, where the subtraction $[x(s) - x(0)]$ tames the $1/s^2$ singularity because a smooth even test function has $x(s) - x(0) = O(s^2)$. Confirm that the subtracted integral is finite by evaluating it against a concrete even test function $x(s) = \cos ks$:

```wl
Integrate[Csch[Pi s]^2 (Cos[k s] - 1), {s, 0, Infinity}, Assumptions -> k > 0]
```

As one can see, the subtracted integral is the finite closed form $(2 - k\operatorname{coth}(k/2))/(2\pi)$. The lesson to carry forward is a rule of practice: integrate the subtracted integrand $\sinh^{-2}(\pi s)[x(s) - x(0)]$, never the bare kernel.

Now expand the subtraction to lowest order to read off the correction. For a smooth even test function, $x(s) - x(0) \approx \tfrac12 x''(0)\,s^2$, so the leading correction is weighted by the second moment of the kernel. Compute that moment, the single integral that fixes the size of the whole correction:

```wl
Integrate[x^2 Csch[x]^2, {x, 0, Infinity}]
```

As one can see, $\int_0^\infty x^2/\sinh^2 x\,dx = \pi^2/6$. This one number, the same $\pi^2/6$ that appears in the Basel problem, sets the weight of the term the white-noise limit discarded. In the variables of the kernel the same moment reads $\int_0^\infty \sinh^{-2}(\pi s)\,s^2\, ds = 1/(6\pi)$:

```wl
Integrate[Csch[Pi s]^2 s^2, {s, 0, Infinity}]
```

Verify that the two computations agree by expanding the subtracted test-function integral at small $k$, which isolates the same second moment:

```wl
Series[(2 - k Coth[k/2])/(2 Pi), {k, 0, 3}] // Normal
```

As one can see, the small-$k$ expansion is $-k^2/(12\pi)$, which is exactly $-\tfrac12 k^2$ times the second moment $1/(6\pi)$: the coefficient in front of $x''(0)$ is fixed, and it is nonzero. Collecting the leading singular term and this correction, the kernel expanded one step past white noise is

$$
k_r(\tau) \approx \frac{4\eta}{\tau_B}\,\delta(\tau) - \frac{\eta\tau_B}{3}\,\delta''(\tau).
$$

The first term is the standard white-noise kernel of Part II. The second term, proportional to the second derivative of the delta and to the same $\pi^2/6$ moment, is the correction. It vanishes as $\tau_B\to 0$, so it is invisible in the strict high-temperature limit, but at any finite temperature it is present, and it is what the next part turns into a $(\hat P,\hat P)$ term in the master equation.

One more closed form deserves its own cell, because it is the antiderivative that controls the boundary contributions to the action (Sec. III of the Supplemental Material). Compute $B(y) = \int_0^y x^2/\sinh^2 x\, dx$, the running version of the $\pi^2/6$ moment:

```wl
Integrate[x^2 Csch[x]^2, {x, 0, y}, Assumptions -> y > 0]
```

Confirm that this antiderivative saturates at the moment we already found, so the boundary terms reach a finite plateau:

```wl
Limit[y (y - y Coth[y] + 2 Log[1 - E^(-2 y)]) - PolyLog[2, E^(-2 y)] + Pi^2/6,
  y -> Infinity]
```

As one can see, $B(y)$ climbs from zero and levels off at $\pi^2/6$, the polylogarithm $\mathrm{Li}_2(e^{-2y})$ carrying the approach to the plateau. The rate of that approach is set by $\tau_B$, so as the temperature rises the boundary contributions collapse onto the endpoints of the time interval and drop out, which is why the bulk correction above is the whole story at high temperature.

Visualize how quickly $B(y)$ reaches its plateau for a few thermal times, reproducing the behavior of Fig. 3 of the paper:

```wl
With[{b = Function[{tt, tauB}, With[{y = Pi tt/tauB},
     y (y - y Coth[y] + 2 Log[1 - E^(-2 y)]) - PolyLog[2, E^(-2 y)] + Pi^2/6]]},
  Plot[Evaluate@Table[b[t, tauB], {tauB, {0.5, 1, 2}}], {t, 0, 3},
   PlotLegends -> (Row[{"\[Tau]B = ", #}] & /@ {0.5, 1, 2}),
   AxesLabel -> None, Frame -> True, FrameLabel -> {"t", "B(\[Pi] t/\[Tau]B)"},
   GridLines -> Automatic, PlotLabel -> "Boundary antiderivative saturates at \[Pi]^2/6",
   Epilog -> {Dashed, Gray, Line[{{0, Pi^2/6}, {3, Pi^2/6}}]}, ImageSize -> 460]]
```

As one can see, a smaller thermal time reaches the plateau sooner, so at high temperature the deviation of $B$ from $\pi^2/6$ is confined to a vanishing neighborhood of the endpoints. This is the computational content of the statement that the end terms are negligible.

## Part IV: From Kernel to Lindblad: The Kossakowski Matrix and Guaranteed Positivity

Now translate the kernel into a generator. Rewriting the influence functional in superoperator form turns the corrected kernel into a Markovian master equation for the reduced density operator, the Liouvillian

$$
\mathcal{L}\hat\rho = -\frac{i}{\hbar}[\hat H_S,\hat\rho] - i\frac{\gamma}{\hbar}[\hat X,\{\hat P,\hat\rho\}] - D_{PP}[\hat X,[\hat X,\hat\rho]] - D_{XX}[\hat P,[\hat P,\hat\rho]] - 2D_{XP}[\hat X,[\hat P,\hat\rho]],
$$

with friction coefficient $\gamma = \eta/M$ and Dekker diffusion coefficients

$$
D_{XX} = \frac{\gamma}{6Mk_BT}, \qquad D_{PP} = \frac{2\gamma Mk_BT}{\hbar^2}, \qquad D_{XP} = 0.
$$

These coefficients are not fresh inputs: they are what the two kernel terms of Part III become under the Feynman-Vernon to superoperator map (Aurell 2020), which sends a $\delta(\tau)$ term of weight $w$ to a decoherence coefficient $w/(2\hbar)$ on the $(\hat X,\hat X)$ channel, and a $\delta''(\tau)$ term of weight $w$ to $w/(2M^2\hbar)$ on the $(\hat P,\hat P)$ channel. Substitute the physical identifications $\eta = M\gamma$ and $\tau_B = \hbar/k_BT$ into the white-noise weight $4\eta/\tau_B$ and the $\delta''$ weight $\eta\tau_B/3$, and read off the two coefficients:

```wl
Simplify[{(4 \[Eta]/\[Tau]B)/(2 \[HBar]), (\[Eta] \[Tau]B/3)/(2 M^2 \[HBar])}
   /. {\[Eta] -> M \[Gamma], \[Tau]B -> \[HBar]/(kB T)}]
```

As one can see, the two weights land exactly $D_{PP} = 2\gamma Mk_BT/\hbar^2$ and $D_{XX} = \gamma/(6Mk_BT)$. The white-noise $\delta$ term is the standard position-decoherence coefficient of the old treatment; the $\delta''$ term, carrying the same $\pi^2/6$ moment, is the new momentum-diffusion coefficient, and the $1/M^2$ that sets it apart is the pair of time derivatives in $\delta''$ turning position differences into momentum.

The single new object is $D_{XX}$, the position-diffusion coefficient carried by the $(\hat P,\hat P)$ double commutator. It is exactly the $\delta''$ correction of Part III, and in the standard Caldeira-Leggett treatment it is zero. Two facts are worth stating before any computation: this generator holds for any potential $V(\hat X)$, because $V$ only shapes the free evolution and not the diffusion timescales, and every assembled coefficient is independent of the cutoff $\Omega$. That the finite-$\Omega$ remainder of the kernel (the $\Omega^2\,O((\Omega\tau)^{-3})$ tail of Part III) drops out in the scaling limit is the paper's distributional argument; what we can check directly, against Table I below, is that no cutoff survives in the assembled coefficients.

Whether such a generator describes a legitimate quantum evolution is decided by one matrix. In Lindblad (GKSL) form the dissipator is $\sum_{ij} a_{ij}(\hat L_i\hat\rho\hat L_j^\dagger - \tfrac12\{\hat L_j^\dagger\hat L_i,\hat\rho\})$, and the dynamics is completely positive if and only if the Kossakowski matrix $\boldsymbol a$ is positive semidefinite. With the two Hermitian Lindblad operators $\hat L_1 = \hat X$ and $\hat L_2 = \hat P$, the diffusion coefficients and the friction assemble into a $2\times 2$ matrix.

Define the Kossakowski matrix as a function of the three Dekker coefficients, with the friction fixing the off-diagonal:

```wl
ClearAll[kossakowski];
kossakowski[dxx_, dpp_, dxp_] := {{2 dpp, 2 dxp - I \[Gamma]/\[HBar]},
                                   {2 dxp + I \[Gamma]/\[HBar], 2 dxx}};
```

Build the matrix for the coefficients of this work and display it:

```wl
kossakowski[\[Gamma]/(6 M kB T), 2 \[Gamma] M kB T/\[HBar]^2, 0] // MatrixForm
```

The diagonal carries the two diffusion coefficients and the off-diagonal carries the friction as a pure imaginary $\mp i\gamma/\hbar$. The question of positivity is now the sign of a determinant.

Compute the determinant of the Kossakowski matrix for this work:

```wl
Det[kossakowski[\[Gamma]/(6 M kB T), 2 \[Gamma] M kB T/\[HBar]^2, 0]] // FullSimplify
```

As one can see, the determinant is $+\gamma^2/(3\hbar^2)$, strictly positive. Positivity of a Kossakowski matrix is a statement about a Hermitian matrix, so confirm both that the matrix is Hermitian and that its two leading principal minors are positive, which is Sylvester's criterion for positive-definiteness, assuming the physical parameters are positive:

```wl
With[{a = kossakowski[\[Gamma]/(6 M kB T), 2 \[Gamma] M kB T/\[HBar]^2, 0]},
  FullSimplify[{a == ConjugateTranspose[a], a[[1, 1]] > 0, Det[a] > 0},
    {\[Gamma] > 0, \[HBar] > 0, M > 0, kB > 0, T > 0}]]
```

As one can see, the matrix is Hermitian and both leading principal minors are positive, so it is positive definite and the generator is genuinely Lindblad. Now set the position-diffusion coefficient to zero, which is the standard Caldeira-Leggett master equation, and recompute the determinant:

```wl
Det[kossakowski[0, 2 \[Gamma] M kB T/\[HBar]^2, 0]] // FullSimplify
```

As one can see, the determinant is now $-\gamma^2/\hbar^2$, strictly negative. A negative determinant of a Hermitian matrix means one eigenvalue is negative, so the Caldeira-Leggett generator is not completely positive. This is the discrepancy the paper resolves: the missing $(\hat P,\hat P)$ term is exactly what the sign of this determinant is waiting for.

Assemble the comparison that the paper presents as Table I, computing each generator's determinant symbolically and reading its Lindblad verdict:

```wl
ClearAll[detA];
detA[dxx_, dxp_] := Det[kossakowski[dxx, 2 \[Gamma] M kB T/\[HBar]^2, dxp]] // FullSimplify;
Grid[{
   {"master equation", "det a", "Lindblad"},
   {"This work", detA[\[Gamma]/(6 M kB T), 0], "yes"},
   {"Di\[OAcute]si", detA[\[Gamma]/(6 M kB T), \[CapitalOmega] \[Gamma]/(6 Pi kB T)], "yes"},
   {"Breuer-Petruccione", detA[\[Gamma]/(8 M kB T), -\[Gamma] kB T/(\[HBar]^2 \[CapitalOmega])], "marginal"},
   {"Caldeira-Leggett", detA[0, 0], "no"}},
  Frame -> All, Alignment -> Left]
```

As one can see, the four rows separate cleanly by the sign of a single determinant. This work is positive. The Di\[OAcute]si determinant is positive in its stated regime $k_BT\ge\hbar\Omega$, since the $\Omega^2$ subtraction stays below the leading $3/\hbar^2$. The Breuer-Petruccione determinant is $-4k_B^2T^2\gamma^2/\Omega^2\hbar^4$, negative for finite $\Omega$ and only marginally positive as its cross-diffusion is dropped, which is the approximation its authors note. And Caldeira-Leggett is negative. Notice the columns for the other three all contain the cutoff $\Omega$; only this work's determinant is free of it.

Make that observation literal. Collect this work's two nonzero coefficients alongside the cross-diffusion terms that Di\[OAcute]si and Breuer-Petruccione carry, and test each for freedom from the cutoff:

```wl
FreeQ[#, \[CapitalOmega]] & /@ {\[Gamma]/(6 M kB T), 2 \[Gamma] M kB T/\[HBar]^2,
   \[CapitalOmega] \[Gamma]/(6 Pi kB T), -\[Gamma] kB T/(\[HBar]^2 \[CapitalOmega])}
```

As one can see, this work's coefficients are free of the cutoff while both other cross-diffusion terms retain it. Keeping the next term of the kernel is exactly what lets $\Omega\to\infty$ be taken at fixed temperature, so the cutoff leaves no trace in the assembled generator, which is the sense in which this is the genuinely cutoff-free Markovian limit.

Before leaving the matrix, verify that it really does encode the generator we wrote down, not just a matrix with a convenient determinant. This check and the covariance cross-check in Part V both need the same ingredients: the position and momentum operators of a truncated harmonic oscillator, built from the annihilation matrix, together with the commutator and anticommutator as plain functions (neither the built-in `Commutator` nor `Anticommutator` reduces on explicit matrices, both staying as an inert `NonCommutativeMultiply`, so we build `Dot`-based versions). Define them once, and confirm the operators satisfy the canonical commutator $[\hat X,\hat P] = i\hbar$ away from the top of the truncation:

```wl
ClearAll[comm, acomm, oscillatorOps];
comm[u_, v_] := u . v - v . u;
acomm[u_, v_] := u . v + v . u;
oscillatorOps[n_, hbar_, mm_, w0_] := With[
   {a = SparseArray[Band[{1, 2}] -> Sqrt[Range[n - 1]], {n, n}],
    x0 = Sqrt[hbar/(2 mm w0)], p0 = Sqrt[hbar mm w0/2]},
   {x0 (a + Transpose[a]), I p0 (Transpose[a] - a)}];
With[{xp = oscillatorOps[8, 1, 1, 1]}, Chop[(comm @@ xp)[[1 ;; 5, 1 ;; 5]] - I IdentityMatrix[5]]]
```

As one can see, on the retained subspace the commutator equals $i\hbar$ (here $\hbar = 1$): the zero matrix shown is its difference from $i\,\mathbb{1}$. The one place this fails is the last row of the truncation, where the coupling to the discarded level $n+1$ is missing, which is why every numerical check below reads its moments from random states supported well below the cutoff, where the operators are canonical.

Now assemble the GKSL dissipator from $\boldsymbol a$ with $\hat L_1 = \hat X$, $\hat L_2 = \hat P$, and compare it to the double-commutator Liouvillian acting on a random density matrix:

```wl
Block[{n = 50, hbar = 1, mm = 1, w0 = 1, gam = 1/10, kt = 5, xop, pop, ls, amat, dpp, dxx, gksl, dekker, lamb, rho, blk},
  {xop, pop} = oscillatorOps[n, hbar, mm, w0];
  dpp = 2 gam mm kt/hbar^2; dxx = gam/(6 mm kt);
  ls = {xop, pop}; amat = {{2 dpp, -I gam/hbar}, {I gam/hbar, 2 dxx}};
  gksl[r_] := Sum[amat[[i, j]] (ls[[i]] . r . ConjugateTranspose[ls[[j]]]
       - (1/2) acomm[ConjugateTranspose[ls[[j]]] . ls[[i]], r]), {i, 2}, {j, 2}];
  dekker[r_] := -I (gam/hbar) comm[xop, acomm[pop, r]] - dpp comm[xop, comm[xop, r]]
       - dxx comm[pop, comm[pop, r]];
  lamb[r_] := (I gam/(2 hbar)) comm[acomm[xop, pop], r];
  SeedRandom[7]; blk = RandomComplex[{-1 - I, 1 + I}, {12, 12}]; blk = blk . ConjugateTranspose[blk];
  rho = SparseArray[{}, {n, n}] + 0.; rho[[1 ;; 12, 1 ;; 12]] = blk; rho = rho/Tr[rho];
  Max@Abs@Flatten[(gksl[rho] - dekker[rho] - lamb[rho])[[1 ;; 20, 1 ;; 20]]]]
```

As one can see, the GKSL dissipator built from $\boldsymbol a$ reproduces the double-commutator generator to machine precision, once a single Hamiltonian piece $\tfrac{i\gamma}{2\hbar}[\{\hat X,\hat P\},\hat\rho]$ is included. That piece is not dissipation: it is a frequency-renormalization term, the Lamb shift that the paper's effective Hamiltonian absorbs, and it does not touch the positivity verdict. The Kossakowski matrix is faithful to the generator. And because the generator is now certified to be in GKSL form with a Hermitian, positive-definite $\boldsymbol a$, the other half of complete positivity comes for free: the map preserves the trace and the Hermiticity of $\hat\rho$ exactly, and Part V watches the completely-positive half play out as a covariance that never pierces its floor. Whether that positivity survives once more terms of the kernel are kept is the question Part VII takes up, and the answer is that it does not: keeping exactly two terms is special.

## Part V: The Physics Made Visible: Covariance Dynamics and the Uncertainty Floor

A determinant tells us the generator is completely positive; the physics of what that buys is clearest for a harmonic potential $V(\hat X) = \tfrac12 M\omega_0^2\hat X^2$, because there the master equation closes. The Liouvillian is quadratic in $\hat X$ and $\hat P$, so the first moments and the $2\times 2$ covariance matrix $\sigma$ form a closed linear system, a drift-plus-diffusion flow

$$
\dot\sigma = A\sigma + \sigma A^\top + D,
$$

with a drift matrix $A$ built from the free evolution and the friction, and a diffusion matrix $D$ built from the Dekker coefficients. Let us assemble both, verify them, and then watch positivity live or die.

The first moments obey the classical damped oscillator. Write the drift matrix and read its eigenvalues:

```wl
Amat = {{0, 1/M}, {-M \[Omega]0^2, -2 \[Gamma]}};
FullSimplify[Eigenvalues[Amat], {M > 0, \[Gamma] > 0, \[Omega]0 > 0}]
```

As one can see, the eigenvalues are $-\gamma\pm\sqrt{\gamma^2 - \omega_0^2}$, the two roots of a damped harmonic oscillator with damping rate $2\gamma$ and frequency $\omega_0$. Both have negative real part, so the mean position and momentum spiral to rest: the first moments know nothing about temperature or $\hbar$, only about friction.

The diffusion lives entirely in the second moments. For this work the diffusion matrix is diagonal, carrying position diffusion from the new $D_{XX}$ and momentum diffusion from $D_{PP}$. Write it:

```wl
Dmat[dxx_, dpp_] := 2 \[HBar]^2 {{dxx, 0}, {0, dpp}};
Dmat[\[Gamma]/(6 M kB T), 2 \[Gamma] M kB T/\[HBar]^2] // MatrixForm
```

The lower entry $4\gamma Mk_BT$ is the familiar thermal momentum diffusion that heats the particle; the upper entry $\gamma\hbar^2/(3Mk_BT)$ is the new position diffusion, small and proportional to $\hbar^2$, and identically zero in the standard treatment. These two matrices, $A$ and $D$, are the whole harmonic dynamics, so they had better be right.

Verify $A$ and $D$ against the generator itself, by building the Liouvillian from primitives on a truncated oscillator and reading the moment derivatives straight off a random state. If the drift and diffusion are correct, the residuals vanish:

```wl
Block[{n = 60, hbar = 1, mm = 1, w0 = 1, gam = 1/10, kt = 5, xop, pop, hs,
   ldot, rho, blk, sym, mean, cov, lr, dxdt, dpdt, dcov, ap, dp, sxp0, sxpdot},
  {xop, pop} = oscillatorOps[n, hbar, mm, w0];
  hs = pop . pop/(2 mm) + (1/2) mm w0^2 xop . xop;
  ldot[r_] := -(I/hbar) comm[hs, r] - I (gam/hbar) comm[xop, acomm[pop, r]]
     - (2 gam mm kt/hbar^2) comm[xop, comm[xop, r]] - (gam/(6 mm kt)) comm[pop, comm[pop, r]];
  SeedRandom[42]; blk = RandomComplex[{-1 - I, 1 + I}, {15, 15}]; blk = blk . ConjugateTranspose[blk];
  rho = SparseArray[{}, {n, n}] + 0.; rho[[1 ;; 15, 1 ;; 15]] = blk; rho = rho/Tr[rho];
  sym = (xop . pop + pop . xop)/2;
  mean = {Re@Tr[xop . rho], Re@Tr[pop . rho]};
  sxp0 = Re@Tr[sym . rho] - mean[[1]] mean[[2]];
  cov = {{Re@Tr[xop . xop . rho] - mean[[1]]^2, sxp0}, {sxp0, Re@Tr[pop . pop . rho] - mean[[2]]^2}};
  lr = ldot[rho]; dxdt = Re@Tr[xop . lr]; dpdt = Re@Tr[pop . lr];
  sxpdot = Re@Tr[sym . lr] - dxdt mean[[2]] - mean[[1]] dpdt;
  dcov = {{Re@Tr[xop . xop . lr] - 2 mean[[1]] dxdt, sxpdot}, {sxpdot, Re@Tr[pop . pop . lr] - 2 mean[[2]] dpdt}};
  ap = Amat /. {M -> mm, \[Omega]0 -> w0, \[Gamma] -> gam};
  dp = 2 hbar^2 {{gam/(6 mm kt), 0}, {0, 2 gam mm kt/hbar^2}};
  Chop[dcov - (ap . cov + cov . Transpose[ap] + dp), 10.^-6]]
```

As one can see, the covariance residual is the zero matrix: the hand-assembled drift and diffusion reproduce the exact generator, computed by an entirely independent route on truncated matrices. The Lyapunov flow is trustworthy.

Now solve for the steady state. Setting $\dot\sigma = 0$ gives a linear system for the three independent entries of $\sigma$; solve it symbolically with real positive parameters:

```wl
ClearAll[steady];
steady[dmat_] := Module[{sxx, sxp, spp, sig, sol},
   sig = {{sxx, sxp}, {sxp, spp}};
   sol = First@Solve[Thread[Flatten[Amat . sig + sig . Transpose[Amat] + dmat] == 0], {sxx, sxp, spp}];
   FullSimplify[sig /. sol, {\[HBar] > 0, M > 0, \[Omega]0 > 0, \[Gamma] > 0, kB > 0, T > 0}]];
sigTW = steady[Dmat[\[Gamma]/(6 M kB T), 2 \[Gamma] M kB T/\[HBar]^2]];
sigTW // MatrixForm
```

Read off the steady momentum variance and the steady correlation, the two entries where the quantum correction shows:

```wl
<|"\[Sigma]PP" -> sigTW[[2, 2]], "\[Sigma]XP" -> sigTW[[1, 2]]|>
```

As one can see, the steady momentum variance is $Mk_BT + M\omega_0^2\hbar^2/(12k_BT)$, the classical equipartition value $Mk_BT$ plus a small correction of order $\hbar^2$, and the position-momentum correlation is $-\gamma\hbar^2/(6k_BT)$, a purely quantum cross-correlation that is exactly zero in the standard treatment. The new position diffusion leaves its fingerprint on the steady state, not just on positivity.

Confirm that this steady state reduces to the classical thermal state when $\hbar\to 0$, which is the sanity check every high-temperature master equation must pass:

```wl
Limit[sigTW, \[HBar] -> 0] // MatrixForm
```

As one can see, the classical limit is $\mathrm{diag}(k_BT/M\omega_0^2,\, Mk_BT)$, the Boltzmann covariance of a classical oscillator in thermal equilibrium: position variance $k_BT/M\omega_0^2$ and momentum variance $Mk_BT$, with no correlation. The quantum corrections are the $\hbar^2$ terms that this limit discards.

Physical states must respect the Robertson-Schrödinger uncertainty relation, which for a $2\times 2$ covariance is the single inequality $\det\sigma \ge \hbar^2/4$. Check the margin by which the steady state clears this floor:

```wl
Det[sigTW] - \[HBar]^2/4 // FullSimplify
```

As one can see, the steady margin is $k_B^2T^2/\omega_0^2$ plus corrections, dominated by a large positive term at high temperature, so the steady state sits comfortably above the floor. The standard Caldeira-Leggett steady state clears it too. The failure of positivity is therefore not a steady-state effect; it is transient, and to see it we must watch a specific initial state evolve.

Here is where the missing term earns its place. Consider a minimum-uncertainty initial state that is squeezed in position, narrow in $\hat X$ and correspondingly broad in $\hat P$, sitting exactly on the floor $\det\sigma = \hbar^2/4$. On the floor with zero correlation, $d(\det\sigma)/dt = \sigma_{XX}'\,\sigma_{PP} + \sigma_{XX}\,\sigma_{PP}'$, with $\sigma_{XX}' = D_{11}$ and $\sigma_{PP}' = 4\gamma(Mk_BT - \sigma_{PP})$. Compute this rate for standard Caldeira-Leggett, where the position diffusion $D_{11}$ is absent:

```wl
ClearAll[ddetFloor, \[Sigma]XX, \[Sigma]PP];
ddetFloor[d11_] := Simplify[d11 \[Sigma]PP + \[Sigma]XX (4 \[Gamma] (M kB T - \[Sigma]PP))];
ddetFloor[0]
```

As one can see, the standard rate is $4\gamma\,\sigma_{XX}(Mk_BT - \sigma_{PP})$, which is negative precisely when the state is broad in momentum ($\sigma_{PP} > Mk_BT$): the determinant immediately drops below $\hbar^2/4$ and the state becomes unphysical. Now compute the term by which this work's rate differs, the rescue:

```wl
ddetFloor[2 \[HBar]^2 \[Gamma]/(6 M kB T)] - ddetFloor[0] // Simplify
```

As one can see, this work adds exactly $2\hbar^2 D_{XX}\,\sigma_{PP} = \gamma\hbar^2\sigma_{PP}/(3Mk_BT)$, proportional to the very momentum variance that caused the trouble, so for a strongly squeezed state that term dominates and keeps the rate positive. The rescue is not a coincidence of parameters; it is the same $D_{XX}$ whose sign fixed the Kossakowski determinant, acting now on the covariance.

Watch it happen. Solve the covariance flow in closed form from a position-squeezed minimum-uncertainty state under both generators, and plot how far $\det\sigma(t)$ sits above or below the uncertainty floor, zoomed on the early transient where the two part ways (the inset carries the full history). The determinant is real, so `Re` only strips the machine-precision residue of the complex-exponential solution:

```wl
Block[{hbar = 1, mm = 1, w0 = 1, gam = 1/10, kt = 5, eps = 1/40, a0, dtw, dcl, floor, run},
  a0 = Amat /. {M -> mm, \[Omega]0 -> w0, \[Gamma] -> gam};
  dtw = 2 hbar^2 {{gam/(6 mm kt), 0}, {0, 2 gam mm kt/hbar^2}};
  dcl = ReplacePart[dtw, {1, 1} -> 0];   (* Caldeira-Leggett = this work with D_XX -> 0 *)
  floor = hbar^2/4;
  run[dm_] := Module[{sxx, sxp, spp, s, d, sol, tt},
    s = {{sxx[tt], sxp[tt]}, {sxp[tt], spp[tt]}}; d = a0 . s + s . Transpose[a0] + dm;
    sol = DSolveValue[{sxx'[tt] == d[[1, 1]], sxp'[tt] == d[[1, 2]], spp'[tt] == d[[2, 2]],
       sxx[0] == eps/2, sxp[0] == 0, spp[0] == 1/(2 eps)}, {sxx, sxp, spp}, tt];
    Function[u, Re[sol[[1]][u] sol[[3]][u] - sol[[2]][u]^2] - floor]];
  With[{fcl = run[dcl], ftw = run[dtw]},
   Plot[{fcl[t], ftw[t]}, {t, 0, 0.15}, PlotRange -> {{0, 0.15}, {-0.0032, 0.006}},
    PlotStyle -> {{Red, Thick}, {Blue, Thick}}, Frame -> True, GridLines -> Automatic,
    FrameLabel -> {"t", "det \[Sigma](t) - \[HBar]^2/4"},
    PlotLabel -> "Early transient: only Caldeira-Leggett pierces the floor",
    PlotLegends -> Placed[{"Caldeira-Leggett (D_XX = 0)", "This work (D_XX > 0)"}, Below],
    Epilog -> {Dashed, Gray, Line[{{0, 0}, {0.15, 0}}],
      Inset[Plot[{fcl[t], ftw[t]}, {t, 0, 8}, PlotStyle -> {{Red}, {Blue}}, Frame -> True,
         PlotRange -> All, ImageSize -> 180, PlotLabel -> Style["full range: both thermalize", 8],
         Epilog -> {Dashed, Gray, Line[{{0, 0}, {8, 0}}]}], Scaled[{0.72, 0.34}]]},
    ImageSize -> 560]]]
```

As one can see, the standard Caldeira-Leggett curve dips below the dashed floor almost immediately: the squeezed state is driven into the forbidden region $\det\sigma < \hbar^2/4$, where the density operator has a negative eigenvalue and is no longer a state. The corrected curve never crosses the floor. Both eventually climb to the same thermal steady state, as the inset shows, so the distinction is a transient one, present exactly when it matters, at the moment of closest approach to the quantum limit. The same single term that lifted the Kossakowski determinant from negative to positive is what holds the covariance above its quantum limit, and the two facts are the same fact seen at the level of the matrix and at the level of the trajectory.

## Part VI: The Systematic Expansion of the Kernel

The paper stops one term past white noise, at the $\delta''$ correction, and names its own open problem in the Conclusions: the higher order terms in the expansion of the fluctuation paths. Nothing in Part III forced the stop. The regularized functional $T_F[x] = 2\pi\big(x(0) - \pi\int_0^\infty ds\,\sinh^{-2}(\pi s)\,(x(s) - x(0))\big)$ accepts the entire Taylor expansion of the fluctuation-path product, and every additional power of $s^2$ meets the same $\sinh^{-2}$ kernel and picks out one more of its moments. This part carries the expansion to all orders in closed form, and the tower it builds turns out to be the exact Ohmic kernel of Part I, resummed.

The object being expanded is the fluctuation-path product $x_{\bar t}(s) = \Delta X(\bar t + \tau_B s/2)\,\Delta X(\bar t - \tau_B s/2)$, even in $s$ because sending $s\to -s$ only swaps its two factors. Expand it in $s$ and collect the coefficient of each $(\tau_B s)^{2n}$, the objects $x^{(n)}$:

```wl
Clear[u, s, tb, \[Tau]B];
ser = Normal@Series[u[tb + \[Tau]B s/2] u[tb - \[Tau]B s/2], {s, 0, 4}];
Table[Coefficient[ser, s, 2 n]/\[Tau]B^(2 n) // Simplify, {n, 0, 2}]
```

The zeroth coefficient is $\Delta X^2$, the white-noise piece. The first is $\tfrac14(\Delta X\,\Delta X'' - \Delta X'^2)$, the single correction the paper keeps. The next brings in the fourth derivative $\Delta X\,\Delta X''''$ and the curvature $(\Delta X'')^2$, and the pattern holds at every order: $x^{(n)}$ is a bilinear in $\Delta X$ and its derivatives of total order $2n$, namely $x^{(n)} = \tfrac{1}{4^n}\sum_{j+k=2n}\tfrac{(-1)^k}{j!\,k!}\,\Delta X^{(j)}\Delta X^{(k)}$.

Inside the action each coefficient is integrated over $\bar t$, and integration by parts collapses the bilinear onto a single canonical term $\Delta X\,\partial_{\bar t}^{2n}\Delta X$. The $(-1)^k$ in the coefficient and the $(-1)^j$ from moving $j$ derivatives by parts multiply to $(-1)^{2n} = +1$, so every term adds with the same sign and the surviving weight is the pure combinatorial number $\tfrac{1}{4^n}\sum_{j+k=2n}\tfrac{1}{j!\,k!} = \tfrac{1}{(2n)!}$, independent of the path. Confirm the path-independence on two different even fluctuation paths, taking for each the ratio of $\int x^{(n)}\,d\bar t$ to $\int \Delta X\,\partial^{2n}\Delta X\,d\bar t$:

```wl
xcoef[n_, f_, t_] := (1/4^n) Sum[(-1)^(2 n - j)/(j! (2 n - j)!) D[f, {t, j}] D[f, {t, 2 n - j}], {j, 0, 2 n}];
ratio[f_] := Table[Integrate[xcoef[n, f, t], {t, -Infinity, Infinity}]/
      Integrate[f D[f, {t, 2 n}], {t, -Infinity, Infinity}] // Simplify, {n, 1, 3}];
ratio /@ {Exp[-t^2], Sech[t]}
```

The ratio is $1/(2n)!$ for both paths, so the weight is path-independent as the binomial sum promised: $\int d\bar t\,x^{(n)} = \tfrac{1}{(2n)!}\int d\bar t\,\Delta X\,\partial_{\bar t}^{2n}\Delta X$. At $n = 1$ this is the factor $\tfrac12$ sitting behind the paper's $S^{(1)}$; at every higher order it is the reciprocal factorial that keeps the assembled series under control.

Each order carries a moment of the kernel as its weight, the integral of $(\tau_B s)^{2n}$ against $\sinh^{-2}(\pi s)$. Part III already met the $n = 1$ case, $\int_0^\infty s^2\sinh^{-2}(\pi s)\,ds = 1/(6\pi)$, the Basel number $\pi^2/6$ in the $x = \pi s$ variable. Compute the family of moments the expansion needs:

```wl
Table[Integrate[x^(2 n) Csch[x]^2, {x, 0, Infinity}], {n, 1, 3}]
```

These are $\zeta(2)$, $\zeta(4)$, $\zeta(6)$ up to rational factors, and in closed form $\int_0^\infty x^{2n}\sinh^{-2}x\,dx = (2n)!\,\zeta(2n)/2^{2n-1} = (-1)^{n+1}B_{2n}\,\pi^{2n}$, with $B_{2n}$ the Bernoulli numbers. The $\pi^2/6$ that fixed the paper's single correction is the first member of an infinite family indexed by the Bernoulli numbers; in the $s$ variable the same family reads $\int_0^\infty s^{2n}\sinh^{-2}(\pi s)\,ds = (-1)^{n+1}B_{2n}/\pi$.

Multiply each coefficient $x^{(n)}$ by its moment and read the result back into the kernel. A term $\int d\bar t\,\Delta X\,\partial^{2n}\Delta X$ in the action is a $\delta^{(2n)}(\tau)$ in the kernel, so the systematic expansion is an even-derivative series $k_r(\tau)\approx\sum_n c_{2n}\,\delta^{(2n)}(\tau)$. Assemble its coefficients:

```wl
Clear[\[Eta], \[Tau]B];
c[n_] := 4 \[Eta] (-1)^n BernoulliB[2 n]/(2 n)! \[Tau]B^(2 n - 1);
Table[c[n], {n, 0, 3}]
```

The first two are exactly the white-noise and $\delta''$ coefficients of Part III. The third, $c_4 = -\eta\tau_B^3/180$, is the leading new term the paper leaves open, a $\delta''''(\tau)$ in the kernel, and its sign is negative. The general coefficient $c_{2n} = \tfrac{4\eta(-1)^n B_{2n}}{(2n)!}\tau_B^{2n-1}$ alternates in sign with the Bernoulli numbers, a fact that will decide positivity in the next part. Written out, the kernel three terms past white noise is

$$
k_r(\tau)\approx \frac{4\eta}{\tau_B}\delta(\tau) - \frac{\eta\tau_B}{3}\delta''(\tau) - \frac{\eta\tau_B^3}{180}\delta''''(\tau) - \dots
$$

The tower is not an arbitrary list of corrections. In Fourier space a $\delta^{(2n)}(\tau)$ is $(i\omega)^{2n}$, so the kernel's symbol is a power series in $\omega$ whose coefficients are the Bernoulli numbers, the signature of a single hyperbolic function. Sum the tower and compare it to the exact Ohmic thermal weight $2\eta\omega\operatorname{coth}(\hbar\omega/2k_BT)$ of Part I:

```wl
Clear[\[Omega]];
Normal@Series[
  Sum[c[n] (I \[Omega])^(2 n), {n, 0, 6}] - 2 \[Eta] \[Omega] Coth[\[Tau]B \[Omega]/2],
  {\[Omega], 0, 11}]
```

The difference vanishes to every order checked, so the delta-derivative tower is the term-by-term Fourier expansion of the exact kernel. The whole systematic expansion resums to the colored-noise spectrum we began with in Part I, with $\tau_B = \hbar/k_BT$, and the Bernoulli numbers are its Taylor coefficients. White noise is the zeroth term of $\tfrac{\tau_B\omega}{2}\operatorname{coth}\tfrac{\tau_B\omega}{2}$, the paper's $\delta''$ is the first, and keeping all of them is the same as not expanding the kernel at all.

## Part VII: Higher Orders in the Generator, Positivity, and Temperature Reach

The kernel expansion of Part VI is clean and closed. Turning it into a master equation is where the physics turns, because each $\delta^{(2n)}$ must be handed to an operator, and only the first two land in the tidy two-operator Lindblad form. This part maps the higher terms, asks whether the extended generator is still completely positive, and measures how far down in temperature each new term reaches.

The map of Part IV sent $\delta(\tau)$ to the position operator $\hat X$ and $\delta''(\tau)$ to $\hat P/M$, the velocity, because a white noise multiplying $\dot{\Delta X}$ couples to $\dot{\hat X} = \hat P/M$. The $\delta''''$ term multiplies $\ddot{\Delta X}$, so its operator is the acceleration $\ddot{\hat X}$ under the system Hamiltonian, fixed by the Heisenberg relation $\ddot{\hat X} = (i/\hbar)^2[\hat H_S,[\hat H_S,\hat X]]$. Evaluate that double commutator on a truncated oscillator, for a harmonic and a quartic potential, and compare it to $-V'(\hat X)/M$:

```wl
Block[{n = 40, hbar = 1, mm = 1, w0 = 1, lam = 1/7, xop, pop, hHarm, hQuart, top = 25},
  {xop, pop} = oscillatorOps[n, hbar, mm, w0];
  hHarm = pop . pop/(2 mm) + (1/2) mm w0^2 xop . xop;
  hQuart = pop . pop/(2 mm) + lam MatrixPower[xop, 4];
  {Max@Abs@Flatten@Normal[((I/hbar)^2 comm[hHarm, comm[hHarm, xop]] + w0^2 xop)[[1 ;; top, 1 ;; top]]],
   Max@Abs@Flatten@Normal[((I/hbar)^2 comm[hQuart, comm[hQuart, xop]]
       + 4 lam MatrixPower[xop, 3]/mm)[[1 ;; top, 1 ;; top]]]}]
```

Both residuals vanish, confirming $\ddot{\hat X} = -V'(\hat X)/M$, the force divided by the mass. The first time derivative $\dot{\hat X} = \hat P/M$ was kinematic, the same for every potential, which is why $\delta''$ produced a universal $(\hat P,\hat P)$ term. The second is not: it is the force, so the operator carried by $\delta''''$ depends on the potential. For a free particle $V' = 0$ and it vanishes, so the series truncates at $\delta''$ and the two-operator generator is exact to all orders. For the harmonic potential $V' = M\omega_0^2\hat X$ it is proportional to $\hat X$ again, so $\delta''''$ folds back into the position channel and only renormalizes its coefficient. For an anharmonic potential $V'\propto\hat X^3$ and beyond, it is a genuinely new Lindblad operator and the generator leaves the two-operator family. This is what the paper meant by higher orders being future work: past $\delta''$, the generator's operator content is potential-dependent.

The natural guess for a higher momentum-diffusion term is a quartic $[\hat P^2,[\hat P^2,\rho]]$, from reading each time derivative as one power of momentum. The computation says otherwise: two derivatives make an acceleration, not a squared momentum, so the operator is $V'(\hat X)$, not $\hat P^2$. A second surprise waits in the sign. Each channel carries a noise variance, related to the kernel coefficient by $\kappa_n = (-1)^n c_{2n} = 4\eta B_{2n}\tau_B^{2n-1}/(2n)!$: the $2n$ paired time derivatives put back the alternating sign that $c_{2n}$ carried, so the variance is left with the bare sign of $B_{2n}$, and complete positivity needs each $\kappa_n$ non-negative. Read those signs off the Bernoulli numbers:

```wl
Table[Sign[BernoulliB[2 n]], {n, 0, 5}]
```

The first two are positive, which is precisely why the two-term generator of Part IV is completely positive: the white-noise channel and the $\delta''$ momentum channel are both legitimate. But $B_4 < 0$, so the $\delta''''$ channel has a negative variance, an anti-diffusion. Complete positivity is therefore not a property the series carries term by term; it is a coincidence of the first two Bernoulli signs. Only the full resummation is safe, because the exact spectrum $2\eta\omega\operatorname{coth}(\hbar\omega/2k_BT)$ is non-negative at every $\omega$, so the exact non-Markovian influence functional is a genuine positive Gaussian even where its local truncations are not.

See the obstruction directly. A negative-variance channel is a double commutator $[\hat A,[\hat A,\rho]]$ entering with a positive coefficient, the wrong (anti-diffusion) sign for decoherence, and for the anharmonic case $\hat A\propto\hat X^3$. For a pure state $\rho = |\psi\rangle\langle\psi|$ and any $\phi\perp\psi$, the quadratic form $\langle\phi|[\hat A,[\hat A,\rho]]|\phi\rangle = -2|\langle\phi|\hat A|\psi\rangle|^2$ is negative, so one anti-diffusion step must push an eigenvalue of $\rho$ below zero. Apply it to the oscillator ground state and read the smallest eigenvalue:

```wl
Block[{n = 40, xop, pop, aop, psi, rho},
  {xop, pop} = oscillatorOps[n, 1, 1, 1];
  aop = MatrixPower[xop, 3];
  psi = Normal@SparseArray[{1} -> 1, {n}];
  rho = Outer[Times, psi, psi];
  Min@Eigenvalues[rho + (1/100) comm[aop, comm[aop, rho]]]]
```

The smallest eigenvalue is exactly $-3/80$, negative: after a single step under the anti-diffusion, $\rho$ has left the set of states. This is the honest verdict on the extended generator. Truncated at $\delta''''$ for an anharmonic potential it is not completely positive, because the $B_4 < 0$ channel is anti-diffusive and its operator $\hat X^3$ cannot be absorbed into the position and momentum channels. The paper's guaranteed Lindblad form is genuinely special to the two-term order, not the first step of a term-by-term positive series.

The higher terms are not useless, though. For the harmonic potential $\delta''''$ folds into the position channel, and there it sharpens the equilibrium. On the stationary paths $\ddot{\Delta X} = -\omega_0^2\Delta X$, so the $\delta''''$ term adds an effective white-noise weight $c_4\omega_0^4$ to the $(\hat X,\hat X)$ channel and shifts the momentum-diffusion coefficient by $\Delta D_{PP} = c_4\omega_0^4/(2\hbar)$, through the same weight-over-$2\hbar$ map as Part IV. Recompute the steady momentum variance with this shift and compare, order by order in $\hbar$, to the exact quantum Gibbs value $\tfrac{\hbar M\omega_0}{2}\operatorname{coth}\tfrac{\hbar\omega_0}{2k_BT}$. First the two-term generator:

```wl
sPPgibbs = \[HBar] M \[Omega]0/2 Coth[\[HBar] \[Omega]0/(2 kB T)];
DXX = \[Gamma]/(6 M kB T); DPP = 2 \[Gamma] M kB T/\[HBar]^2;
Normal@Series[steady[Dmat[DXX, DPP]][[2, 2]] - sPPgibbs, {\[HBar], 0, 5}]
```

The two-term generator reproduces the Gibbs momentum variance through order $\hbar^2$, and the first discrepancy is order $\hbar^4$. Now add the harmonic $\delta''''$ shift and recompute:

```wl
dppShift = (-\[Eta] \[Tau]B^3/180) \[Omega]0^4/(2 \[HBar]) /. {\[Eta] -> M \[Gamma], \[Tau]B -> \[HBar]/(kB T)};
Normal@Series[steady[Dmat[DXX, DPP + dppShift]][[2, 2]] - sPPgibbs, {\[HBar], 0, 7}]
```

The $\delta''''$ correction cancels the $\hbar^4$ discrepancy exactly, and the agreement with the Gibbs state now holds through order $\hbar^4$, the residual pushed down to $\hbar^6$. Each order of the fluctuation-path expansion buys one more power of $\hbar^2$ in the equilibrium match, so the systematic expansion is literally an expansion in $\hbar\omega_0/k_BT$: the higher terms carry the master equation's validity from high temperature down toward the quantum regime, one Bernoulli number at a time.

That downward reach has a floor. The expansion parameter is $z = \tau_B\omega = \hbar\omega/k_BT$, and the series is the Taylor expansion of $\tfrac{z}{2}\operatorname{coth}\tfrac{z}{2}$, which converges only up to the nearest singularity of the hyperbolic cotangent. Locate it:

```wl
z /. Solve[Sinh[z/2] == 0 && 0 < Im[z] < 4 Pi, z]
```

The nearest pole sits at $z = 2\pi i$, that is $\hbar\omega/k_BT = 2\pi$, the first bosonic Matsubara frequency $\omega = 2\pi k_BT/\hbar$. The delta-derivative tower converges only for frequencies below it, so the whole Markovian construction is a statement about modes with $\hbar\omega\lesssim 2\pi k_BT$. For the harmonic oscillator this is the condition $k_BT\gtrsim\hbar\omega_0/2\pi$: below it, the first Matsubara mode makes the bath memory genuinely non-Markovian and no finite number of delta-derivatives can follow. That $2\pi$ carries two meanings at once. It is the convergence radius of the harmonic observable, whose series inherited the $1/(2n)!$ of the smooth-path reduction of Part VI; and it is the boundary of useful truncation for the bare expansion, whose raw moments $(2n)!\,\zeta(2n)/2^{2n-1}$ grow factorially on their own, so the underlying series in $\tau_B$ is asymptotic rather than convergent, and the reciprocal factorials of a smooth path are exactly what rescue it. This is the precise temperature reach the paper's motivation invokes, and the precise place the systematic expansion ends.

## Part VIII: The Hu-Paz-Zhang Benchmark

The generator of Parts IV and V is Markovian: constant coefficients, and a $(\hat P,\hat P)$ term. The exact reduced dynamics of the same harmonic particle in a real bath is the Hu-Paz-Zhang dynamics, with time-dependent coefficients and no $(\hat P,\hat P)$ term at all, so the constant-coefficient generator is not a limit of it in any naive sense. The paper argues this is a feature: taking $\Omega\to\infty$ in the other master equations of Table I also forces $T\to\infty$, collapsing the kernel to a pure $\delta(\tau)$ and losing positivity, whereas keeping the next kernel term lets $\Omega\to\infty$ at fixed temperature. Whether the constant-coefficient generator is nonetheless the correct Markovian limit is a quantitative question, and for a harmonic system it is fully answerable, because the reduced dynamics of a linear bath is Gaussian and can be computed exactly from the microscopic model, with no master equation in between. The Hu-Paz-Zhang master equation is only the analytic packaging of that exact dynamics; here we compute the dynamics itself and let the constant-coefficient generator stand next to it.

The microscopic model is the Caldeira-Leggett Hamiltonian with the bath written out: the harmonic system coupled to $N$ oscillators whose frequencies and couplings discretize the Ohmic spectral density $J(\omega) = 2\eta\omega$, the counter-term keeping $\omega_0$ physical. The whole Hamiltonian is quadratic, so the full phase-space covariance evolves exactly under the symplectic flow $\Sigma(t) = e^{At}\Sigma(0)e^{A^\top t}$, and its $2\times 2$ system block is the exact reduced covariance $\sigma(t)$. Since $A = JG$ is a symplectic (Hamiltonian) generator, its spectrum is purely imaginary and the propagated covariance is real. Build the propagator, with the bath sampled thermally and the system in the squeezed state of Part V, and confirm the discretization has reached the continuum by comparing $N = 200$ against $N = 300$ modes over the whole window used below. The finite bath recurs at $t \approx 2\pi N/\Omega \approx 42$ for $N = 200$, so every window in this part stays before that revival:

```wl
ClearAll[exactPieces];
exactPieces[nn_, cutoff_] := Module[
   {mm = 1., w0 = 1., eta = 0.1, kT = 5., eps = 1/40., sys0, dw, wk, ck, dim, g, jm, amat, sig0, vals, vecs, vv, vinv, th},
   sys0 = {{eps/2, 0.}, {0., 1/(2 eps)}};
   dw = N[cutoff/nn]; wk = (Range[nn] - 1/2) dw; ck = wk Sqrt[4 eta dw/Pi];
   dim = 2 (nn + 1);
   g = SparseArray[Join[{{1, 1} -> mm w0^2 + Total[ck^2/wk^2], {2, 2} -> 1/mm},
       Flatten@Table[{{2 k + 1, 2 k + 1} -> wk[[k]]^2, {2 k + 2, 2 k + 2} -> 1.,
          {1, 2 k + 1} -> -ck[[k]], {2 k + 1, 1} -> -ck[[k]]}, {k, nn}]], {dim, dim}];
   jm = SparseArray[Band[{1, 1}] -> Table[{{0., 1.}, {-1., 0.}}, {nn + 1}]];
   amat = Normal[jm . g];
   th[w_] := {{Coth[w/(2 kT)]/(2 w), 0.}, {0., w Coth[w/(2 kT)]/2}};
   sig0 = Normal@SparseArray[Band[{1, 1}] -> Join[{sys0}, th /@ wk]];
   {vals, vecs} = Eigensystem[amat]; vv = Transpose[vecs]; vinv = Inverse[vv];
   {Function[t, Re@With[{p = vv[[1 ;; 2]] . (Exp[vals t] vinv)}, p . sig0 . Transpose[p]]],
    Function[t, Re@With[{p = vv[[1 ;; 2]] . (Exp[vals t] vinv), pd = vv[[1 ;; 2]] . (vals Exp[vals t] vinv)}, pd . sig0 . Transpose[p] + p . sig0 . Transpose[pd]]],
    Function[t, Re[vv[[1 ;; 2]] . (Exp[vals t] vinv[[All, 1 ;; 2]])]],
    Function[t, Re[vv[[1 ;; 2]] . (vals Exp[vals t] vinv[[All, 1 ;; 2]])]],
    Max@Abs@Re[vals]}];
pieces = exactPieces[200, 30]; p300 = First[exactPieces[300, 30]];
Max@Table[Max@Abs@Flatten[First[pieces][t] - p300[t]], {t, 0, 38, 0.5}]
```

The $N = 200$ versus $N = 300$ covariance difference is negligible across the whole window, so a few hundred modes reproduce the continuum Ohmic bath everywhere the dynamics is read below. Confirm too that the full drift is a symplectic generator, its spectrum purely imaginary, so the propagated covariance is real and reading the real part strips only machine residue:

```wl
Last[pieces]
```

The spectral real part is machine zero, so the construction is a genuine oscillatory, energy-conserving flow. This $\sigma(t)$ is the exact answer, the target the Markovian generator must hit.

Read off the drift and diffusion the exact dynamics carries. Written as a covariance flow $\dot\sigma = A(t)\sigma + \sigma A(t)^\top + D(t)$, the instantaneous drift is $A(t) = \dot G(t)\,G(t)^{-1}$ with $G$ the reduced first-moment propagator, and the diffusion is what remains. These coefficients ring persistently rather than settle, so read each one's time-average (its center) and its range over a window, and set the centers beside the Markovian values $-2\gamma$ and $4\gamma Mk_BT$:

```wl
{sig, sigd, gprop, gpropd} = pieces[[1 ;; 4]];
adAt[t_] := Module[{a = gpropd[t] . Inverse[gprop[t]], d},
   d = sigd[t] - a . sig[t] - sig[t] . Transpose[a]; {a, d}];
Module[{ts = Range[1., 20., 0.25], fr, dxx, dpp, dxp, stat},
  {fr, dxx, dpp, dxp} = Transpose@Table[
     With[{ad = adAt[t]}, {ad[[1, 2, 2]], ad[[2, 1, 1]], ad[[2, 2, 2]], ad[[2, 1, 2]]}], {t, ts}];
  stat[x_] := Chop[{Mean[x], Min[x], Max[x]}, 1.*^-10];
  Grid[{{"", "center", "min", "max"}, Prepend[stat[fr], "friction"],
     Prepend[stat[dxx], "D_XX"], Prepend[stat[dpp], "D_PP"], Prepend[stat[dxp], "D_XP"]},
   Frame -> All, Alignment -> Left]]
```

Against the constant drift of Part V, whose friction entry is $-2\gamma$, and its diagonal diffusion $D = 2\hbar^2\,\mathrm{diag}(D_{XX},D_{PP})$, three features stand out. The exact coefficients do not settle to constants: they oscillate persistently around centers that are exactly the Markovian values, the friction around $-2\gamma$ and the momentum diffusion around the classical Einstein value $4\gamma Mk_BT = 2\hbar^2 D_{PP}$. The propagator $G(t)$ rings slowly, at the system frequency $\omega_0$, but the extracted coefficient $\dot G G^{-1}$ divides out that decaying envelope and exposes the bath's fast memory underneath, so it rings far faster, at the cutoff $\Omega$; the oscillation does not damp, its Markovian center is reached almost at once on the short memory time $\sim 1/\Omega$, and the ringing is what the constant-coefficient generator smooths away. The position-diffusion entry, though, is not merely centered on zero, it is identically zero at every time, $D(t)_{XX} = 0$: Hu-Paz-Zhang has no $(\hat P,\hat P)$ term. So the dominant momentum diffusion agrees with this work, and the difference is in the small correction, where this work carries a constant position diffusion $2\hbar^2 D_{XX}$ and no cross, while the exact dynamics carries an anomalous cross $D_{XP}(t)$ centered on zero and no position diffusion: two memoryless stand-ins for the same memory, both completely positive.

Now the physical test, the uncertainty floor of Part V. A position-squeezed state starts on the floor $\det\sigma = \hbar^2/4$. Compute the smallest margin $\det\sigma(t) - \hbar^2/4$ reaches over the whole evolution, for the exact dynamics, this work, and standard Caldeira-Leggett:

```wl
Block[{mm = 1., w0 = 1., gam = 1/10, kT = 5., sys0 = {{1/80, 0.}, {0., 20.}}, am, dm, dcl, ssm, ssc, sm, sc},
  am = {{0., 1/mm}, {-mm w0^2, -2 gam}};
  dm = 2 {{gam/(6 mm kT), 0.}, {0., 2 gam mm kT}}; dcl = {{0., 0.}, {0., 4 gam mm kT}};
  ssm = LyapunovSolve[am, -dm]; ssc = LyapunovSolve[am, -dcl];
  sm[t_] := With[{e = MatrixExp[am t]}, ssm + e . (sys0 - ssm) . Transpose[e]];
  sc[t_] := With[{e = MatrixExp[am t]}, ssc + e . (sys0 - ssc) . Transpose[e]];
  Chop@{Min@Table[Det[sig[t]] - 1/4, {t, 0, 35, 0.05}],
        Min@Table[Det[sm[t]] - 1/4, {t, 0, 35, 0.05}],
        Min@Table[Det[sc[t]] - 1/4, {t, 0, 35, 0.05}]}]
```

The exact dynamics and this work never go negative over the entire evolution; only Caldeira-Leggett does, dipping below the floor into the region where the density operator has a negative eigenvalue. The exact curve cannot cross the floor, because the reduced state of a genuine unitary evolution is a state at every instant; that this work stays above it too is the covariance-level face of its guaranteed positivity. See it, zoomed on the early transient where the three part ways, with the full relaxation in the inset:

```wl
Block[{mm = 1., w0 = 1., gam = 1/10, kT = 5., sys0 = {{1/80, 0.}, {0., 20.}}, am, dm, dcl, ssm, ssc, sm, sc, dE, dM, dC},
  am = {{0., 1/mm}, {-mm w0^2, -2 gam}};
  dm = 2 {{gam/(6 mm kT), 0.}, {0., 2 gam mm kT}}; dcl = {{0., 0.}, {0., 4 gam mm kT}};
  ssm = LyapunovSolve[am, -dm]; ssc = LyapunovSolve[am, -dcl];
  sm[t_] := With[{e = MatrixExp[am t]}, ssm + e . (sys0 - ssm) . Transpose[e]];
  sc[t_] := With[{e = MatrixExp[am t]}, ssc + e . (sys0 - ssc) . Transpose[e]];
  dE[t_?NumericQ] := Det[sig[t]] - 1/4; dM[t_?NumericQ] := Det[sm[t]] - 1/4; dC[t_?NumericQ] := Det[sc[t]] - 1/4;
  Plot[{dC[t], dM[t], dE[t]}, {t, 0, 0.2}, PlotRange -> {{0, 0.2}, {-0.004, 0.03}},
   PlotStyle -> {{Red, Thick}, {Blue, Thick}, {Black, Dashed}}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"t", "det \[Sigma](t) - \[HBar]^2/4"},
   PlotLabel -> "Only Caldeira-Leggett pierces the floor; this work tracks the exact dynamics",
   PlotLegends -> Placed[{"Caldeira-Leggett (D_XX = 0)", "This work", "Exact (Hu-Paz-Zhang)"}, Below],
   Epilog -> {Dashed, Gray, Line[{{0, 0}, {0.2, 0}}],
     Inset[Plot[{dC[t], dM[t], dE[t]}, {t, 0, 35}, PlotStyle -> {{Red}, {Blue}, {Black, Dashed}},
        Frame -> True, PlotRange -> All, ImageSize -> 190,
        PlotLabel -> Style["full range: all thermalize", 8]], Scaled[{0.73, 0.34}]]},
   ImageSize -> 580]]
```

This work tracks the exact curve closely, the two differing most in the initial transient, where the constant-coefficient generator switches on its full diffusion at $t = 0$ while the exact dynamics ramps its up. Over the trajectory the constant-coefficient covariance follows the exact one to about a percent, the deviation living in that initial slip, the one regime a Markovian description is not meant to capture.

That leaves the equilibrium, and there the two genuinely differ. Average the exact covariance over a late window to strip the ringing, and compare the settled variances to the bare Gibbs values and to this work's steady state:

```wl
Block[{mm = 1., w0 = 1., gam = 1/10, kT = 5., gibbs, avg, am, dm, ssm},
  gibbs = Coth[w0/(2 kT)]/(2 w0);
  avg = Mean@Table[Diagonal[sig[t]], {t, 28, 38, 0.2}];
  am = {{0., 1/mm}, {-mm w0^2, -2 gam}}; dm = 2 {{gam/(6 mm kT), 0.}, {0., 2 gam mm kT}};
  ssm = LyapunovSolve[am, -dm];
  Grid[Chop@{{"", "\[Sigma]XX", "\[Sigma]PP"},
     {"exact (averaged)", avg[[1]], avg[[2]]},
     {"bare Gibbs", gibbs, gibbs},
     {"this work steady", ssm[[1, 1]], ssm[[2, 2]]}}, Frame -> All, Alignment -> Left]]
```

The exact position variance settles onto the bare Gibbs value, to the precision the averaging affords, which validates the microscopic construction against a known closed form. The exact momentum variance, though, sits about a percent above it, robustly across sub-windows: the finite-coupling mean-force enhancement, the shift of the true equilibrium away from the bare Gibbs state when the system is coupled to the bath at finite strength. The constant-coefficient generator, pinned to the bare Gibbs value plus its $O(\hbar^2)$ correction, does not carry this $O(\gamma)$ shift, so at equilibrium it and the exact dynamics agree only to the order of the coupling. That gap is the honest edge of the surrogate: the constant-coefficient generator is faithful to Hu-Paz-Zhang through the transient and positive where Caldeira-Leggett is not, but it is a memoryless stand-in, not the exact dynamics, and the mean-force shift is what a memoryless generator tied to the bare Gibbs state must leave on the table.

## Where This Leaves Us (and What Comes Next)

Before closing, let us summarize the most important points we computed and verified:

- The Ohmic Feynman-Vernon kernel is a hyperbolic-cosecant-squared profile of width $\tau_B$, obtained here as the second derivative of the thermal primitive $\pi\operatorname{coth}\pi s - 1/s$.
- White noise keeps only the $2/x$ term of $\operatorname{coth}(x/2)$; the discarded $x/6$ term is the whole correction, and its weight is fixed by the single moment $\int_0^\infty x^2/\sinh^2 x\, dx = \pi^2/6$.
- The kernel is only meaningful as a tempered distribution: the bare $\sinh^{-2}$ integral diverges at the origin, and the subtraction $[x(s) - x(0)]$ is the regularization that makes it finite.
- Keeping one term past white noise adds a $\delta''(\tau)$ piece to the kernel, hence a $(\hat P,\hat P)$ double commutator with position-diffusion coefficient $D_{XX} = \gamma/(6Mk_BT)$ to the master equation.
- With $\hat L_1 = \hat X$, $\hat L_2 = \hat P$, the Kossakowski determinant is $+\gamma^2/(3\hbar^2)$ for this work and $-\gamma^2/\hbar^2$ for standard Caldeira-Leggett: the correction is exactly what turns a non-positive generator into a Lindblad one, and the corrected generator contains no cutoff $\Omega$.
- For a harmonic potential the master equation is a Lyapunov flow whose steady state is the classical thermal covariance plus $\hbar^2$ corrections, and a position-squeezed state that pierces the Robertson-Schrödinger floor under Caldeira-Leggett stays above it under this work.
- Continued past $\delta''$, the fluctuation-path expansion has coefficients weighted by the Bernoulli numbers, $c_{2n} = 4\eta(-1)^n B_{2n}\tau_B^{2n-1}/(2n)!$, and the even-derivative tower $\sum_n c_{2n}\delta^{(2n)}(\tau)$ resums to the exact Ohmic weight $2\eta\omega\operatorname{coth}(\hbar\omega/2k_BT)$: the local expansion is the high-temperature expansion of the exact colored-noise spectrum.
- Past $\delta''$ the generator's operator content is potential-dependent, because the $\delta''''$ term carries the acceleration $-V'(\hat X)/M$ under the system Hamiltonian: it vanishes for a free particle, where the two-operator form is exact to all orders; it renormalizes the coefficients for a harmonic potential; and it becomes a genuinely new Lindblad operator for an anharmonic well.
- Complete positivity is special to the two-term order: the channel variances alternate in sign with the Bernoulli numbers, so the $\delta''''$ channel is anti-diffusive ($B_4 < 0$) and its truncation is not completely positive; positivity is restored only by the full resummation to a manifestly non-negative spectrum.
- For the harmonic oscillator each further term deepens the match to the Gibbs covariance by one power of $\hbar^2$ ($\delta''$ matches $\sigma_{PP}$ to $\hbar^2$, $\delta''''$ to $\hbar^4$), and the expansion converges only for $k_BT\gtrsim\hbar\omega_0/2\pi$, the first Matsubara frequency, below which the bath memory is genuinely non-Markovian.
- Benchmarked against the exact reduced dynamics of the microscopic model (the Hu-Paz-Zhang dynamics, computed by discretizing the bath), the constant-coefficient generator is a faithful memoryless surrogate: it tracks the exact covariance to about a percent through the transient and never pierces the uncertainty floor that Caldeira-Leggett does. It is a stand-in, not the exact dynamics, in two honest ways: it fixes positivity differently (a constant position diffusion, where the exact dynamics carries a time-dependent anomalous cross and no position diffusion), and at equilibrium it misses the $O(\gamma)$ mean-force enhancement of the momentum variance, because it is pinned to the bare Gibbs state.

You now have a complete, computation-first toolkit for the corrected high-temperature limit of quantum Brownian motion: the kernel algebra from primitives, the distributional regularization done honestly, the Kossakowski matrix and its spectrum as the positivity test, the covariance flow as the physical picture, the systematic expansion that settles the operator structure, positivity, and temperature reach of the higher orders in closed form, and the exact microscopic benchmark that weighs the whole construction against the dynamics it approximates. The natural continuations are the anharmonic generator on a position grid, where the new $V'(\hat X)$ operator and the loss of a finite Lindblad form show up directly in the density matrix, and the stochastic unraveling of the generator into a Schrodinger equation with a momentum-coupled noise, where the two-term generator's new piece appears as a $\hat P\phi(t)$ coupling. The two-term generator keeps the invariant we leaned on throughout, that nothing in it depends on the cutoff $\Omega$; the systematic expansion shows the price of going further is a return of the potential and, eventually, of the memory.

### Acknowledgment

This notebook follows the derivation of Pleasance, Aurell, and Petruccione, "High-temperature limit penalizing high-frequency quantum fluctuations" (arXiv:2508.14262), and reproduces its central results as live computations. Any errors in the reconstruction are mine.
