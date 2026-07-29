---
Template: Default
---

# Quantum in Continuous Space: Worked Solutions (Parts 1 to 3)

Pure Wolfram Language answers to [`Question-List.md`](Question-List.md): the wavefunction and the
state space $L^2$ (Part 1), the position and momentum operators and the Fourier bridge (Part 2),
and the time-independent Schrodinger equation as an ODE eigenvalue problem (Part 3). No quantum
package; natural units $\hbar=m=1$. Each answer is self-contained: it defines the state or operator
it uses, computes from that definition, and can be read or run on its own.

## Part 1. The wavefunction and the state space $L^2$

A particle on the line is a complex *wavefunction* $\psi(x)$ in the Hilbert space $L^2(\mathbb{R})$
of square-integrable functions, $\int_{-\infty}^{\infty}|\psi|^2\,dx<\infty$. Its *inner product*
$\langle\phi|\psi\rangle=\int_{-\infty}^{\infty}\phi^*\psi\,dx$ gives the norm, the orthogonality, and
the overlap whose modulus squared $|\langle\phi|\psi\rangle|^2$ is a transition probability. The
*Born rule* makes $|\psi(x)|^2$ a probability density, so $\int_a^b|\psi|^2\,dx$ is the probability of
finding the particle in $[a,b]$ and a physical state is normalized to $\langle\psi|\psi\rangle=1$.
*Expectation values* are inner products $\langle\hat A\rangle=\langle\psi|\hat A|\psi\rangle$, giving
the mean position $\langle x\rangle$ and the spread $\Delta x$. Probability is locally conserved by
the current $j$ through the continuity equation whenever the potential is real, and because only such
bilinear quantities are observable, an overall constant phase of $\psi$ carries no physics.
Throughout, $\hbar=m=1$.

### 1.1 [BSc] How do I represent a one-dimensional wavefunction $\psi(x)$ as a Wolfram Language function and plot its probability density $|\psi(x)|^2$?

A state on the line is a complex function $\psi(x)$ in $L^2(\mathbb{R})$. Measuring position returns a
single number, one point of the continuous spectrum of $\hat x$, and Born's rule fixes only the
statistics of that outcome: $|\psi(x)|^2\,dx$ is the probability that the result falls in $dx$, a
distribution an ensemble of identically prepared states reproduces and a single run never does. The
density is therefore not the state's whole physical content, only the outcome distribution of one
observable: the phase it discards carries the momentum distribution, the current, and the energy, so
two states with the same density can differ in every other measurement.

We take the ground state of the Morse oscillator, the standard model of a vibrating diatomic molecule.
With the bond measured from its equilibrium length as $x$, the potential is
$V(x)=D\big(e^{-2ax}-2e^{-ax}\big)$: a hard repulsive wall as the atoms are pushed together
($x\to-\infty$) and a soft tail as they are pulled apart ($x\to+\infty$), with well depth $D$ and range
set by $a$. The ground state is $\psi_0(x)\propto z^{\lambda-1/2}e^{-z/2}$ in the variable
$z=2\lambda\,e^{-ax}$, where $\lambda=\sqrt{2D}/a$ is the dimensionless depth. Unlike every symmetric
state in this Part, this one has no parity: it inherits the lopsidedness of the well, and that
asymmetry is the whole physical interest, so it is the feature to watch.

Substitute $\psi_0$ into the Morse stationary equation $-\tfrac12\psi_0''+V\psi_0=E_0\psi_0$ and confirm
the residual vanishes at $E_0=-\tfrac{a^2}{2}\big(\lambda-\tfrac12\big)^2$:

```wl
\[Psi]0[x_] := (2 \[Lambda] Exp[-a x])^(\[Lambda] - 1/2) Exp[-\[Lambda] Exp[-a x]];
Vmorse[x_] := (\[Lambda]^2 a^2/2) (Exp[-2 a x] - 2 Exp[-a x]);
FullSimplify[-(1/2) \[Psi]0''[x] + Vmorse[x] \[Psi]0[x] - (-(a^2/2) (\[Lambda] - 1/2)^2) \[Psi]0[x],
  a > 0 && \[Lambda] > 1/2 && x \[Element] Reals]
```

The residual is $0$ for every $a>0$ and $\lambda>\tfrac12$, so $\psi_0$ is the exact ground state with
energy $E_0=-\tfrac{a^2}{2}\big(\lambda-\tfrac12\big)^2$, a binding energy of
$\tfrac{a^2}{2}\big(\lambda-\tfrac12\big)^2$ below the dissociation threshold at $E=0$. The depth is
written $D=\tfrac{\lambda^2 a^2}{2}$ so that $\lambda$ is the only shape parameter. The state is not yet
normalized; the substitution $z=2\lambda e^{-ax}$ turns the norm into a Gamma integral,
$\int|\psi_0|^2\,dx=\tfrac1a\int_0^\infty z^{2\lambda-2}e^{-z}\,dz$. Compute it:

```wl
norm2 = Integrate[\[Psi]0[x]^2, {x, -Infinity, Infinity}, Assumptions -> a > 0 && \[Lambda] > 1/2]
```

It is $\Gamma(2\lambda-1)/a$, finite exactly when $\lambda>\tfrac12$: the Gamma pole at $2\lambda-1=0$ is
the edge of existence, the shallowest well that still binds a normalizable ground state, below which
$|\psi|^2$ is not a probability density at all. Divide by its square root and confirm the normalized
state carries unit probability:

```wl
\[Psi][x_] := \[Psi]0[x]/Sqrt[norm2];
Integrate[\[Psi][x]^2, {x, -Infinity, Infinity}, Assumptions -> a > 0 && \[Lambda] > 1/2]
```

Now read off the density
$|\psi(x)|^2=\tfrac{a}{\Gamma(2\lambda-1)}\big(2\lambda e^{-ax}\big)^{2\lambda-1}e^{-2\lambda e^{-ax}}$:

```wl
dens = FullSimplify[\[Psi][x]^2, a > 0 && \[Lambda] > 1/2 && x \[Element] Reals]
```

Plot it at $a=1$, $\lambda=\tfrac52$, a representative depth:

```wl
Plot[dens /. {a -> 1, \[Lambda] -> 5/2}, {x, -1, 7}, Filling -> Axis,
 PlotRange -> All, AxesLabel -> {"x", "|\[Psi]|^2"}]
```

The density is a single skewed hump, a steep rise on the compression side and a long tail on the
stretch side, nothing like a symmetric bell. That lopsidedness is what the plot is for. It shows first in where the density peaks: not at the bottom
of the well, but at $x_\star=\tfrac1a\log\!\frac{2\lambda}{2\lambda-1}>0$, displaced toward the soft
side. Confirm that this is where the density is stationary:

```wl
xstar = Log[2 \[Lambda]/(2 \[Lambda] - 1)]/a;
FullSimplify[D[dens, x] /. x -> xstar, a > 0 && \[Lambda] > 1/2]
```

The derivative vanishes at $x_\star>0$: the crest of the plotted density lies to the stretch side of the
force minimum at $x=0$. Cross-check the closed form against a direct numerical maximization at $a=1$,
$\lambda=\tfrac52$:

```wl
{N[xstar /. {a -> 1, \[Lambda] -> 5/2}], x /. Last@FindMaximum[dens /. {a -> 1, \[Lambda] -> 5/2}, {x, 0.2}]}
```

Both locate the peak at $x_\star\approx0.223$. Deepening the well washes the asymmetry out: as
$\lambda\to\infty$ the peak slides back to the origin and the hump tends to the symmetric Gaussian of the
harmonic oscillator.

```wl
Limit[xstar, \[Lambda] -> Infinity]
```

The limit is $0$. The single plot then carries the state's whole character: a normalized probability
density, manifestly lopsided at finite depth with its crest pushed to the soft side, and symmetric only
in the deep-well limit where the Morse well becomes harmonic.

### 1.2 [BSc] How do I normalize a wavefunction with $\int_{-\infty}^{\infty}|\psi|^2\,dx=1$, and recognize a state that cannot be normalized?

Normalizing a state means dividing it by the square root of its total probability, so that
$\int_{-\infty}^{\infty}|\psi|^2\,dx=1$. Equivalently, a wavefunction is admissible only when that
integral is finite, that is $\psi\in L^2$. Take the Lorentzian profile $\phi(x)=\dfrac{1}{x^2+a^2}$, a
heavy-tailed alternative to the usual bell curve. Its shape is the one familiar from Breit-Wigner
lineshapes, though there the Lorentzian is a function of energy rather than position; and note that
the density here, $|\phi|^2\propto(x^2+a^2)^{-2}$, is a *squared* Lorentzian, so the name describes
the amplitude and not the observable distribution.

Compute its total probability $\int|\phi|^2\,dx$:

```wl
\[Phi][x_] := 1/(x^2 + a^2);
norm2 = Integrate[Abs[\[Phi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> a > 0]
```

The integral is $\pi/(2a^3)$, finite, so the profile is normalizable. Before trusting that closed
form, cross-check it against a direct numerical integral at $a=1$, where $\pi/(2a^3)=\pi/2$:

```wl
{norm2 /. a -> 1, NIntegrate[Abs[\[Phi][x]]^2 /. a -> 1, {x, -Infinity, Infinity}]}
```

The symbolic value and the independent numeric quadrature agree, so the closed form is not a
`Simplify` artifact. Divide by its square root and confirm the normalized state
$\psi(x)=\sqrt{\tfrac{2a^3}{\pi}}\,\dfrac{1}{x^2+a^2}$ integrates to $1$:

```wl
\[Psi][x_] := \[Phi][x]/Sqrt[norm2];
Integrate[Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> a > 0]
```

As expected, the result is $1$: the normalized Lorentzian carries exactly unit probability, with a
width set by $a$. Now the opposite case. Compute the total probability of a plane wave $e^{ikx}$,
whose density $|\psi|^2=1$ is the same at every point:

```wl
Integrate[Abs[Exp[I k x]]^2, {x, -Infinity, Infinity}, Assumptions -> Element[{k, x}, Reals]]
```

The integral diverges to `Infinity`: a momentum eigenstate cannot be normalized to unity at all. It
lives outside $L^2$, and the most one can do is normalize it to a Dirac delta in the momentum label.

### 1.3 [BSc] How do I compute the Born-rule probability $\int_a^b|\psi|^2\,dx$ of finding the particle in a region?

The Born rule gives the probability of finding the particle in a window $[a,b]$ as the integral of
the density there, $\int_a^b|\psi|^2\,dx$. Take the sinc state $\psi(x)=\dfrac{\sin x}{x\sqrt\pi}$,
the position wavefunction of a particle whose momentum is spread uniformly over $[-1,1]$; it is the
Wannier function of the lowest free-electron band, built from Bloch waves of uniform weight across
the Brillouin zone, and equally the single-slit diffraction amplitude.

For a windowed integral to be a probability at all, the state must first carry unit total
probability. Confirm $\int_{-\infty}^{\infty}|\psi|^2\,dx=1$:

```wl
\[Psi][x_] := Sin[x]/(x Sqrt[Pi]);
Integrate[Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}]
```

The total norm is $1$; the $\sqrt\pi$ in the denominator is there precisely because
$\int(\sin x/x)^2\,dx=\pi$. Now compute the probability $P(|x|<b)=\int_{-b}^{b}|\psi|^2\,dx$ of
finding the particle within $b$ of the origin:

```wl
P = Integrate[Abs[\[Psi][x]]^2, {x, -b, b}, Assumptions -> b > 0]
```

The probability is $P(b)=\dfrac{-1+\cos 2b+2b\,\operatorname{Si}(2b)}{b\,\pi}$, where
$\operatorname{Si}$ is the sine integral. Cross-check this closed form against a direct numerical
integral at $b=3$:

```wl
{P /. b -> 3, NIntegrate[Abs[\[Psi][x]]^2, {x, -3, 3}]}
```

The symbolic and numeric values agree. A probability must swallow everything as the window opens and
nothing as it closes; verify $P\to1$ as $b\to\infty$ and $P\to0$ as $b\to0$:

```wl
{Limit[P, b -> Infinity], Limit[P, b -> 0]}
```

The limits are $1$ and $0$, exactly as a probability demands. Unlike a state with exponential tails,
the sinc decays only as $1/x$, so its density has long oscillatory tails and $P(b)$ approaches $1$
slowly (through the $\operatorname{Si}(2b)\to\pi/2$ term). That slow approach is the price of the
sharp momentum cutoff: confining momentum to a box smears position out into a $1/x$ tail.

### 1.4 [BSc] How do I compute the inner product $\langle\phi|\psi\rangle=\int\phi^*\psi\,dx$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?

The overlap of two states is $\int\phi^*\psi\,dx$; its squared modulus is the probability that a
system prepared as $\psi$ is found in $\phi$. We take the state $\psi(x)=2\kappa^{3/2}x\,e^{-\kappa x}$
on the half-line, which is the hydrogen $1s$ reduced radial $u(r)=rR(r)$. That claim is earned by
the radial equation: confirm $\psi$ solves the $l=0$ Coulomb radial equation
$-\tfrac12\psi''-\tfrac{\kappa}{x}\psi=E\psi$ with $E=-\kappa^2/2$:

```wl
\[Psi][x_] := 2 \[Kappa]^(3/2) x Exp[-\[Kappa] x];
FullSimplify[-(1/2) \[Psi]''[x] - (\[Kappa]/x) \[Psi][x] - (-\[Kappa]^2/2) \[Psi][x], \[Kappa] > 0 && x > 0]
```

The residual is $0$, so $\psi$ is the $1s$ radial with energy $E=-\kappa^2/2$, a binding energy of
$\kappa^2/2$, where $\kappa$ is the nuclear charge and $1/\kappa$ the Bohr radius. Being a bound state, it is normalized
on the half-line; confirm $\int_0^\infty|\psi|^2\,dx=1$:

```wl
Integrate[Abs[\[Psi][x]]^2, {x, 0, Infinity}, Assumptions -> \[Kappa] > 0]
```

The norm is $1$. Now take the *same* state after a sudden momentum kick $q$, $\phi=e^{iqx}\psi$ (the
survival amplitude in the sudden approximation), and compute the overlap
$\langle\phi|\psi\rangle=\int\phi^*\psi\,dx$:

```wl
\[Phi][x_] := Exp[I q x] \[Psi][x];
ovl = Integrate[Conjugate[\[Phi][x]] \[Psi][x], {x, 0, Infinity}, Assumptions -> \[Kappa] > 0 && q \[Element] Reals]
```

The overlap is $\langle\phi|\psi\rangle=8\kappa^3/(2\kappa+iq)^3$, genuinely complex: the phase
records the momentum kick. Cross-check it against a numerical integral at $\kappa=1,q=2$, where it
should equal $-\tfrac14-\tfrac{i}{4}$:

```wl
{ovl /. {\[Kappa] -> 1, q -> 2}, NIntegrate[Conjugate[\[Phi][x]] \[Psi][x] /. {\[Kappa] -> 1, q -> 2}, {x, 0, Infinity}]}
```

The symbolic overlap and the numeric quadrature agree. Read off its modulus and the survival
probability $|\langle\phi|\psi\rangle|^2$; `ComplexExpand` pulls the real modulus out of the complex
denominator (`Simplify` alone leaves an unreduced `Abs`):

```wl
Simplify[ComplexExpand[{Abs[ovl], Abs[ovl]^2}], \[Kappa] > 0]
```

The modulus is $8\kappa^3/(4\kappa^2+q^2)^{3/2}$ and the survival probability
$64\kappa^6/(4\kappa^2+q^2)^3$. Confirm the hard-kick falloff by expanding the survival probability
at large $q$:

```wl
Series[64 \[Kappa]^6/(4 \[Kappa]^2 + q^2)^3, {q, Infinity, 6}]
```

The leading term is $64\kappa^6/q^6$: a kick much larger than the orbital momentum scale $\kappa$
*leaves* the electron in its initial state with probability only $\sim(\kappa/q)^6$, so it is almost
certainly knocked out, with ejection probability $1-O\big((\kappa/q)^6\big)$. What this computes is a
one-dimensional overlap of a half-line state with its own momentum-shifted copy, the
sudden-approximation survival amplitude for a kick of momentum $q$. It is not the atomic form factor,
which lives in three dimensions and carries an angular integration that produces $\sin(qr)/(qr)$ in
place of the plane wave used here, and which falls off faster as a result.

### 1.5 [BSc] How do I compute the position expectation $\langle x\rangle$ and the spread $\Delta x=\sqrt{\langle x^2\rangle-\langle x\rangle^2}$?

The mean position is $\langle x\rangle=\int x\,|\psi|^2\,dx$ and the spread is the root-variance of
that distribution. Take the first excited state of the harmonic oscillator,
$\psi(x)=\left(\tfrac{4}{\pi}\right)^{1/4}(x-x_0)\,e^{-(x-x_0)^2/2}$, centered at $x_0$. It is a
Hermite function with a single node at $x_0$, not a Gaussian, and its density is a symmetric
double hump about $x_0$.

Every expectation is an average against $|\psi|^2$, so it is meaningful only for a normalized state.
Confirm $\int_{-\infty}^{\infty}|\psi|^2\,dx=1$:

```wl
\[Psi][x_] := (4/Pi)^(1/4) (x - x0) Exp[-(x - x0)^2/2];
Integrate[Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals]
```

The norm is $1$, so the averages need no denominator. Compute the mean position
$\langle x\rangle=\int x\,|\psi|^2\,dx$:

```wl
mx = Integrate[x Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals]
```

The mean is $x_0$: the node sits at the center and the two humps balance. Compute the second moment
$\langle x^2\rangle$:

```wl
mx2 = Integrate[x^2 Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals]
```

It is $x_0^2+\tfrac32$. Cross-check that closed form against a numerical integral at $x_0=0$, where
$\langle x^2\rangle=\tfrac32$:

```wl
{mx2 /. x0 -> 0, NIntegrate[x^2 Abs[\[Psi][x]]^2 /. x0 -> 0, {x, -Infinity, Infinity}]}
```

The symbolic and numeric second moments agree. The spread is
$\Delta x=\sqrt{\langle x^2\rangle-\langle x\rangle^2}$:

```wl
Simplify[Sqrt[mx2 - mx^2], x0 \[Element] Reals]
```

The spread is $\Delta x=\sqrt{3/2}$, larger than the ground state's $\sqrt{1/2}$ because the node
pushes probability outward into the two humps. This is the $n=1$ case of the oscillator law
$\langle x^2\rangle=n+\tfrac12$. Rather than assert that law, widen from this single state to the
whole Hermite ladder $\psi_n(x)=\left(2^n\,n!\,\sqrt{\pi}\right)^{-1/2}H_n(x)\,e^{-x^2/2}$ and read
off the second moment for the first few $n$ at $x_0=0$:

```wl
\[Psi]n[x_, n_] := (2^n n! Sqrt[Pi])^(-1/2) HermiteH[n, x] Exp[-x^2/2];
Table[Integrate[x^2 Abs[\[Psi]n[x, n]]^2, {x, -Infinity, Infinity}], {n, 0, 3}]
```

The moments come out $\left\{\tfrac12,\tfrac32,\tfrac52,\tfrac72\right\}$, exactly $n+\tfrac12$, so the
spread widens as $\sqrt{n+\tfrac12}$ with every node added. And it is independent of where the state
sits: shifting $x_0$ moves the mean but leaves the shape, and so $\Delta x$, untouched.

### 1.6 [BSc] How do I compute the probability current density $j(x,t)=\tfrac{1}{2i}\left(\psi^*\partial_x\psi-\psi\,\partial_x\psi^*\right)$ and verify the continuity equation $\partial_t|\psi|^2+\partial_x j=0$?

Probability is locally conserved: the density's rate of change is minus the divergence of the
current $j=\tfrac{1}{2i}(\psi^*\partial_x\psi-\psi\,\partial_x\psi^*)$ (with $\hbar=m=1$). Its two
terms are complex conjugates of each other, so equivalently $j=\operatorname{Im}(\psi^*\partial_x\psi)$.

The cleanest place to see the current at work is an ordinary normalizable state that actually moves.
Take a particle in an infinite well $[0,L]$ in an equal superposition of the ground and first
excited levels, $\psi(x,t)=\tfrac{1}{\sqrt2}\big(\chi_1 e^{-iE_1 t}+\chi_2 e^{-iE_2 t}\big)$ with
$\chi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$ and $E_n=n^2\pi^2/(2L^2)$: two stationary states beating
against each other. Confirm it is a bona fide state, $\int_0^L|\psi|^2\,dx=1$:

```wl
\[Chi][n_, x_] := Sqrt[2/L] Sin[n Pi x/L];
En[n_] := n^2 Pi^2/(2 L^2);
\[Psi][x_, t_] := (\[Chi][1, x] Exp[-I En[1] t] + \[Chi][2, x] Exp[-I En[2] t])/Sqrt[2];
Integrate[Abs[\[Psi][x, t]]^2, {x, 0, L}, Assumptions -> L > 0 && t \[Element] Reals]
```

The norm is $1$, with no caveats about tails or convergence. Read off the density $\rho=|\psi|^2$,
using `ComplexExpand` to turn the modulus of the complex superposition into real trigonometric form
(it treats $x,t,L$ as real):

```wl
dens = FullSimplify[ComplexExpand[Abs[\[Psi][x, t]]^2], L > 0]
```

The two stationary pieces interfere through a cross term oscillating at the beat frequency
$E_2-E_1=3\pi^2/(2L^2)$: the probability sloshes from one side of the well to the other. That
sloshing is exactly what the current reports. Compute $j=\operatorname{Im}(\psi^*\partial_x\psi)$:

```wl
j = FullSimplify[Im[Conjugate[\[Psi][x, t]] D[\[Psi][x, t], x]], {L > 0, Element[{x, t}, Reals]}]
```

The current is $j=\tfrac{2\pi}{L^2}\sin\!\big(\tfrac{3\pi^2 t}{2L^2}\big)\sin^3\!\big(\tfrac{\pi x}{L}\big)$:
it oscillates in time at the same beat frequency, positive while probability flows toward the right
wall and negative as it flows back, and it vanishes at both walls $x=0,L$, so no probability ever
leaves the box. Verify continuity $\partial_t\rho+\partial_x j=0$:

```wl
FullSimplify[D[dens, t] + D[j, x], {L > 0, Element[{x, t}, Reals]}]
```

The residual is $0$: on this everyday normalizable state, the local bookkeeping of density and
current balances exactly. A single stationary level (drop one term) would give a real wavefunction
with $j\equiv0$ and nothing flowing; it is the interference of two levels that sets the probability
in motion.

A more exotic solution carries the same identity into territory the well never reaches. The
Berry-Balazs *Airy beam* is the exact free-particle solution
$\psi(x,t)=\operatorname{Ai}\!\big(x-\tfrac{t^2}{4}\big)\,e^{i(tx/2-t^3/12)}$: a wave packet that
neither spreads nor decays but *accelerates*, its $\operatorname{Ai}^2$ profile sliding rigidly
along the parabola $x=t^2/4$ with no force acting. It is the quantum cousin of the optical Airy
beam and has been produced with free electrons.

```wl
\[Psi][x_, t_] := AiryAi[x - t^2/4] Exp[I (t x/2 - t^3/12)];
FullSimplify[I D[\[Psi][x, t], t] + (1/2) D[\[Psi][x, t], {x, 2}], Element[{x, t}, Reals]]
```

The residual is $0$: the Airy beam is an exact free-particle solution. Unlike every state above,
though, it is not square-integrable: its $\operatorname{Ai}^2$ tail decays only as $|x|^{-1/2}$ as
$x\to-\infty$, so the norm integral diverges. Confirm it does not converge:

```wl
Integrate[AiryAi[x]^2, {x, -Infinity, Infinity}]
```

The integral fails to converge (`Integrate::idiv`), so the Airy beam lives outside $L^2$, alongside
the plane wave. Its $|\psi|^2$ is therefore an improper density and $j$ a current of a
non-normalizable state, an idealized carrier of probability flux rather than a bona fide
probability distribution. Continuity, being a local identity, still holds pointwise regardless. Read
off the density $|\psi|^2$:

```wl
dens = FullSimplify[Abs[\[Psi][x, t]]^2, Element[{x, t}, Reals]]
```

The density is $\operatorname{Ai}(x-\tfrac{t^2}{4})^2$, a fixed Airy profile whose peak slides along
$x=t^2/4$: the beam accelerates without spreading. Compute the current
$j=\operatorname{Im}(\psi^*\partial_x\psi)$:

```wl
j = FullSimplify[Im[Conjugate[\[Psi][x, t]] D[\[Psi][x, t], x]], Element[{x, t}, Reals]]
```

The current is $j=\tfrac{t}{2}\operatorname{Ai}(x-\tfrac{t^2}{4})^2=\tfrac{t}{2}\,|\psi|^2$, the
density times the velocity $t/2$: every point of the profile moves at the same instantaneous speed,
so the shape is carried rigidly, which is precisely why the beam does not spread. Verify continuity
$\partial_t|\psi|^2+\partial_x j$:

```wl
FullSimplify[D[dens, t] + D[j, x], Element[{x, t}, Reals]]
```

The residual is $0$: probability is locally conserved. Continuity is not special to this state; it
follows from the equation of motion for any potential, provided that potential is real. Writing
$\psi$ and $\psi^*$ as independent functions $f,g$ and imposing
$\partial_t f=\tfrac{i}{2}\partial_x^2 f-iVf$ together with its conjugate makes the residual vanish
for every solution at once:

```wl
Simplify[(D[g[x, t] f[x, t], t]
     + D[(g[x, t] D[f[x, t], x] - f[x, t] D[g[x, t], x])/(2 I), x]) /.
   {Derivative[0, 1][f][x, t] -> (I/2) D[f[x, t], {x, 2}] - I V[x] f[x, t],
    Derivative[0, 1][g][x, t] -> -(I/2) D[g[x, t], {x, 2}] + I V[x] g[x, t]}]
```

The residual is $0$ for any $f,g$ and any $V$, so continuity holds for every state moving in a real
potential, the Airy beam (with $V=0$) being one concrete instance. The reality of $V$ is doing real
work here: the two potential terms cancel between $f$ and $g$ only because the same $V$ appears in
both. Let the potential be absorptive, $V=V_R-i\Gamma$, and the cancellation fails:

```wl
Simplify[(D[g[x, t] f[x, t], t]
     + D[(g[x, t] D[f[x, t], x] - f[x, t] D[g[x, t], x])/(2 I), x]) /.
   {Derivative[0, 1][f][x, t] -> (I/2) D[f[x, t], {x, 2}] - I (VR[x] - I \[CapitalGamma][x]) f[x, t],
    Derivative[0, 1][g][x, t] -> -(I/2) D[g[x, t], {x, 2}] + I (VR[x] + I \[CapitalGamma][x]) g[x, t]}]
```

The residual is now $-2\Gamma\,\rho$, so $\partial_t\rho+\partial_x j=-2\Gamma\rho$: probability
drains away at rate $2\Gamma$ instead of being conserved. That is not a defect but the design of an
optical potential, which models absorption out of the channel being described.

### 1.7 [BSc] How do I represent a state on a discrete spatial grid so that integrals become sums, building the bridge from $L^2$ to the finite vector a computer can hold?

Every numerical method replaces $\psi(x)$ by its samples and every integral by a Riemann sum
$\int|\psi|^2\,dx\approx\sum_n|\psi_n|^2\,\Delta x$. Take the tent (triangular) state
$\psi(x)=\sqrt{\tfrac{3}{2L}}\left(1-\tfrac{|x|}{L}\right)$ on $[-L,L]$ and zero outside, the classic
Rayleigh-Ritz variational trial wavefunction. It has kinks at the peak $x=0$ and at the two edges
$x=\pm L$, exactly the non-smooth features a grid must resolve, and its exact norm is $1$ by
construction.

Confirm the exact norm is $1$, then compare it to a Riemann sum over a grid of $4000$ points on
$[-1,1]$ (setting $L=1$):

```wl
\[Psi][x_] := Sqrt[3/(2 L)] (1 - Abs[x]/L);
exact = Integrate[Abs[\[Psi][x]]^2, {x, -L, L}, Assumptions -> L > 0];
grid = Subdivide[-1., 1., 4000]; dx = grid[[2]] - grid[[1]];
{exact, Total[Abs[\[Psi][grid] /. L -> 1.]^2] dx}
```

The exact norm is $1$, and the grid sum returns $1.0000\ldots$, converging to it as the grid is
refined. The kinks at $0$ and $\pm L$ cost a little accuracy (the density has corners there), which
is why a piecewise-linear trial state is a more honest grid test than a smooth bell curve: it
stresses exactly the corners a real grid struggles with. These tent functions are the building
blocks of the finite element method, where the wavefunction is a sum of them.

### 1.8 [BSc] How do I show that an overall constant phase is unobservable, while the local phase gradient $\partial_x\arg\psi$ carries the current (the velocity field)?

A global phase $e^{i\alpha}$ cancels in $|\psi|^2$ and every expectation value, so it is pure
gauge. A *position-dependent* phase is physical: its gradient is the local velocity $v=j/|\psi|^2$.
Take the chirped packet $\psi(x)=\tfrac{1}{\sqrt2}\operatorname{sech}(x)\,e^{i\beta x^2}$, whose phase
$S(x)=\beta x^2$ varies across the profile.

First, confirm a global phase $e^{i\alpha}$ leaves the density unchanged:

```wl
\[Psi][x_] := Sech[x] Exp[I \[Beta] x^2]/Sqrt[2];
Simplify[Abs[Exp[I \[Alpha]] \[Psi][x]]^2 - Abs[\[Psi][x]]^2, Element[{\[Alpha], \[Beta], x}, Reals]]
```

The difference is $0$: the global phase is unobservable. The $1/\sqrt2$ prefactor normalizes the
state; confirm $\int|\psi|^2\,dx=1$:

```wl
Integrate[Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Beta] \[Element] Reals]
```

The norm is $1$. Now compute the current, again the imaginary part
$j=\operatorname{Im}(\psi^*\partial_x\psi)$, which the position-dependent phase makes nonzero:

```wl
j = Simplify[Im[Conjugate[\[Psi][x]] D[\[Psi][x], x]], Element[{\[Beta], x}, Reals]]
```

The current is $j=\beta x\operatorname{sech}^2(x)$. Divide by the density to get the local velocity
field $v=j/|\psi|^2$:

```wl
Simplify[j/Abs[\[Psi][x]]^2, Element[{\[Beta], x}, Reals]]
```

The velocity field is $v(x)=2\beta x=\partial_x(\beta x^2)$, the slope of the local phase. Unlike a
plain plane-wave phase $kx$, which gives uniform flow, the chirp gives a position-dependent velocity
growing linearly with $x$, and its sign decides the fate of the packet. For $\beta>0$ the flow points
outward everywhere and the packet spreads; for $\beta<0$ it points inward, so the packet focuses,
contracting to a waist before spreading again, which is free spreading run backwards in time.

### 1.9 [BSc] How does a Galilean boost transform a wave packet (the boost phase factor), and why is $|\psi(x)|^2$ unchanged in shape?

A Galilean boost to velocity $v$ multiplies a state by the boost phase and shifts its argument,
$\psi_v(x,t)=e^{i(vx-v^2t/2)}\,\psi(x-vt,t)$ (with $m=1$). This is the *active* boost, which adds $v$
to the momentum of the state; viewing that same state from a frame moving at $+v$ would subtract it
instead.

The content of the statement is that this operation maps solutions to solutions. Rather than test
that on one packet, impose it on a generic free solution $f(x,t)$, replacing every time derivative of
$f$ by the free Schrodinger equation $\partial_t f=\tfrac{i}{2}\partial_x^2 f$:

```wl
\[Psi]B[x_, t_] := Exp[I (v x - v^2 t/2)] f[x - v t, t];
Simplify[(I D[\[Psi]B[x, t], t] + (1/2) D[\[Psi]B[x, t], {x, 2}]) /.
  Derivative[0, 1][f][a_, b_] :> (I/2) Derivative[2, 0][f][a, b]]
```

The residual is $0$ for every free solution $f$, so the boosted state solves the same free equation.
That is Galilean covariance, established once for all states rather than checked on one. Now for what
the phase does and does not do. Confirm it has unit modulus:

```wl
boost = Exp[I (v x - v^2 t/2)];
Simplify[Abs[boost]^2, Element[{v, x, t}, Reals]]
```

The modulus is $1$, so the boost cannot change a density: $|\psi_v(x,t)|^2=|\psi(x-vt,t)|^2$, the
density of the boosted solution being the original's carried along at speed $v$. A position
measurement sees that drift and nothing of the phase. The momentum is another matter, and that is
where the physics lives. Take a snapshot at $t=0$, where the boost is simply $e^{ivx}$, and apply it
to the harmonic oscillator's first excited state
$\psi(x)=\left(\tfrac{4}{\pi}\right)^{1/4}x\,e^{-x^2/2}$, a packet with internal structure (a node).
Compute the mean momentum $\langle p\rangle=\int\psi^*(-i\,\partial_x)\psi\,dx$ of the un-boosted
state:

```wl
\[Psi][x_] := (4/Pi)^(1/4) x Exp[-x^2/2];
Integrate[Conjugate[\[Psi][x]] (-I) D[\[Psi][x], x], {x, -Infinity, Infinity}]
```

The mean momentum is $0$: the stationary excited state carries no net current. Now compute
$\langle p\rangle$ for the boosted state $\psi_v=e^{ivx}\psi$:

```wl
\[Psi]v[x_] := Exp[I v x] \[Psi][x];
Integrate[Conjugate[\[Psi]v[x]] (-I) D[\[Psi]v[x], x], {x, -Infinity, Infinity}, Assumptions -> v \[Element] Reals]
```

The mean momentum is exactly $v$: the boost adds $v$ to every momentum component. It shifts the
energy too. Compute the free-particle energy $\langle H\rangle=\int\psi^*(-\tfrac12\partial_x^2)\psi\,dx$
of the un-boosted state:

```wl
Integrate[Conjugate[\[Psi][x]] (-(1/2)) D[\[Psi][x], {x, 2}], {x, -Infinity, Infinity}]
```

The un-boosted energy is $\tfrac34$, the kinetic energy of the first excited oscillator profile. Now
the boosted state:

```wl
Integrate[Conjugate[\[Psi]v[x]] (-(1/2)) D[\[Psi]v[x], {x, 2}], {x, -Infinity, Infinity}, Assumptions -> v \[Element] Reals]
```

The boosted energy is $\tfrac34+\tfrac12 v^2$: the boost raises it by exactly $\tfrac12 v^2$, the
kinetic energy of the added drift. This is the physical content that the unit-modulus phase alone
cannot show: at a fixed instant $e^{ivx}$ leaves the density identical while adding $v$ to every
momentum component and $\tfrac12 v^2$ to the energy. A position measurement made then sees nothing at
all; a momentum measurement sees the entire boost.

## Part 2. Operators, position and momentum, and the Fourier bridge

In the position representation the observables are operators on $L^2$: position $\hat x$ multiplies by
$x$, and momentum $\hat p=-i\,\partial_x$ differentiates. They do not commute, $[\hat x,\hat p]=i$, the
algebraic root of the uncertainty relation $\Delta x\,\Delta p\ge\tfrac12$. A physical observable must
be *self-adjoint*; formal Hermiticity, $\langle\phi|\hat A\psi\rangle=\langle\hat A\phi|\psi\rangle$, is
only half of that. For $\hat p$ the Hermitian half holds by integration by parts once a boundary term
vanishes, and on a finite domain the gap between the two conditions becomes the whole story. The Fourier transform
$\tilde\psi(p)=\tfrac{1}{\sqrt{2\pi}}\int\psi(x)e^{-ipx}\,dx$ carries a state to the momentum
representation, where $\langle p^n\rangle=\int p^n|\tilde\psi(p)|^2\,dp$ is an ordinary moment and the
momentum eigenfunctions $e^{ipx}$, not themselves in $L^2$, are normalized to a Dirac delta. Throughout,
$\hbar=m=1$.

### 2.1 [BSc] How do I represent the position operator $\hat x$ as multiplication and the momentum operator $\hat p=-i\,d/dx$ as differentiation, acting on a test function?

In the position representation the position operator acts by multiplication, $\hat x\psi=x\psi$, and the
momentum operator by differentiation, $\hat p\psi=-i\,\partial_x\psi$. The cleanest witness is a plane
wave $e^{ikx}$, which carries a single sharp momentum. Define the two operators as functions of the
wavefunction and apply them to $e^{ikx}$:

```wl
xop[f_] := x f;
pop[f_] := -I D[f, x];
{xop[Exp[I k x]], pop[Exp[I k x]]}
```

The results are $x\,e^{ikx}$ and $k\,e^{ikx}$: the plane wave is an eigenfunction of $\hat p$ with
eigenvalue $k$, the precise sense in which $e^{ikx}$ carries momentum $k$, while $\hat x$ only reweights
it by position. Every operator in this part is built from these two operations.

### 2.2 [BSc] How do I verify the canonical commutator $[\hat x,\hat p]=i$ by letting it act on an arbitrary function?

The canonical commutator $[\hat x,\hat p]=i$ is an operator identity, so it must hold acting on an
arbitrary function $f(x)$, not just on a special state. With $\hat x$ multiplying and
$\hat p=-i\,\partial_x$, form $[\hat x,\hat p]f=\hat x\hat p f-\hat p\hat x f$ and simplify:

```wl
xop[f_] := x f;
pop[f_] := -I D[f, x];
Simplify[xop[pop[f[x]]] - pop[xop[f[x]]]]
```

The result is $i\,f(x)$, so $[\hat x,\hat p]=i$. The surviving term is the one from differentiating
the product $x f$ inside $\hat p\hat x$: position and momentum fail to commute by exactly one unit,
and that single unit is what forbids a state from being sharp in both at once. One qualification the
manipulation quietly needs: $f$ must lie in a domain on which both orderings are defined, which on
the line the Schwartz functions supply. That is not decoration. On a box with a twisted boundary
condition, $\hat x\psi$ violates the condition and $\hat p\hat x\psi$ is undefined, which is why the
identity cannot be imposed on just any space.

### 2.3 [BSc] How do I compute the momentum expectation $\langle p\rangle$ and spread $\Delta p$ in the position representation?

Without leaving position space, $\langle p\rangle=\int\psi^*(-i\partial_x)\psi\,dx$, and using the
manifestly positive, boundary-free form $\langle p^2\rangle=\int|\partial_x\psi|^2\,dx$ (equal to
$\int\psi^*\hat p^2\psi\,dx$ because $\hat p$ is Hermitian), the spread is
$\Delta p=\sqrt{\langle p^2\rangle-\langle p\rangle^2}$. Take the first excited state of the
Poschl-Teller well $V(x)=-3\operatorname{sech}^2x$, namely
$\sqrt{\tfrac32}\,\operatorname{sech}(x)\tanh(x)$ at $E=-\tfrac12$ (that well's ground state being
$\operatorname{sech}^2x$ at $E=-2$), a state with a node, and give it a uniform drift with a boost
$e^{ikx}$:
$\psi(x)=\sqrt{\tfrac32}\,\operatorname{sech}(x)\tanh(x)\,e^{ikx}$. Confirm it is normalized:

```wl
\[Psi][x_] := Sqrt[3/2] Sech[x] Tanh[x] Exp[I k x];
Integrate[Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> k \[Element] Reals]
```

The norm is $1$. Compute the mean momentum $\langle p\rangle=\int\psi^*(-i\partial_x\psi)\,dx$:

```wl
mp = Integrate[Conjugate[\[Psi][x]] (-I D[\[Psi][x], x]), {x, -Infinity, Infinity}, Assumptions -> k \[Element] Reals]
```

The mean is $\langle p\rangle=k$, the drift stamped in by the boost; the real, odd envelope carries no
current of its own. Compute the second moment $\langle p^2\rangle=\int|\partial_x\psi|^2\,dx$:

```wl
p2 = Integrate[Abs[D[\[Psi][x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> k \[Element] Reals]
```

It is $\langle p^2\rangle=k^2+\tfrac75$: the drift energy $k^2$ plus an intrinsic piece set by the node.
The spread is $\Delta p=\sqrt{\langle p^2\rangle-\langle p\rangle^2}$:

```wl
Simplify[Sqrt[p2 - mp^2], k \[Element] Reals]
```

The spread is $\Delta p=\sqrt{7/5}$, independent of the boost $k$: shifting the mean momentum slides the
whole distribution rigidly without changing its width, exactly as a Galilean boost should. The width is
set by the node in the envelope, not by the drift.

### 2.4 [BSc] How do I show that $\hat p$ is Hermitian by integration by parts, and see that the boundary terms vanish for an $L^2$ state?

Hermiticity of $\hat p$ means $\langle\phi|\hat p\chi\rangle=\langle\hat p\phi|\chi\rangle$ for states
in its domain. Integration by parts turns the difference of the two matrix elements into a pure surface
term. Writing $\phi^*$ and $\chi$ as independent functions $g$ and $f$, so nothing is assumed of either,
confirm that the integrand of $\langle\phi|\hat p\chi\rangle-\langle\hat p\phi|\chi\rangle$ is the total
derivative $-i\,(gf)'$:

```wl
Simplify[g[x] (-I f'[x]) - (I g'[x]) f[x] + I D[g[x] f[x], x]]
```

The result is $0$, so the difference is $\int_{-\infty}^{\infty}-i\,(gf)'\,dx=-i\,[\phi^*\chi]_{-\infty}^{\infty}$,
a boundary term and nothing else. Hence $\hat p$ is Hermitian exactly when that term vanishes. For any
state that decays at infinity it does; confirm it for a concrete decaying pair
$\phi(x)=\tfrac{1}{1+x^2}$, $\chi(x)=\tfrac{e^{ix}}{1+(x-1)^2}$:

```wl
\[Phi][x_] := 1/(1 + x^2);
\[Chi][x_] := Exp[I x]/(1 + (x - 1)^2);
Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> Infinity] - Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> -Infinity]
```

The boundary term is $0$: both amplitudes vanish at infinity, so the two matrix elements agree. Note
what that does and does not establish. Belonging to $L^2$ does not by itself force a function to
decay pointwise (spikes of fixed height and shrinking width are square-integrable and do not tend to
zero), and $\hat p$ is unbounded, so it is not defined on all of $L^2$ but only on a dense domain. So
$\hat p$ is symmetric on the domain of functions for which this surface term vanishes, not on $L^2$
as a whole. That the whole content sits in a surface term is the crux: on a finite domain the term
need not vanish, and then everything depends on the boundary conditions imposed.

### 2.5 [BSc] How do I write the momentum eigenfunctions $e^{ipx}$, see that they are not square-integrable, and normalize them to a Dirac delta?

The momentum eigenfunctions $e^{ipx}$ have $|e^{ipx}|^2=1$ everywhere, so they are not
square-integrable and cannot be normalized to $1$; the continuum replacement is
$\langle p'|p\rangle=\delta(p-p')$ for $u_p(x)=\tfrac{1}{\sqrt{2\pi}}e^{ipx}$. The $2\pi$ that makes this
work is the Fourier integral of a constant. Compute $\int e^{iqx}\,dx$ as a Fourier transform:

```wl
FourierTransform[1, x, q]
```

The result is $\sqrt{2\pi}\,\delta(q)$, that is $\int e^{iqx}\,dx=2\pi\,\delta(q)$. Now transform the
normalized eigenstate $u_p=\tfrac{1}{\sqrt{2\pi}}e^{ipx}$ into momentum space and read off its overlap
with a momentum label $p'$:

```wl
FourierTransform[1/Sqrt[2 Pi] Exp[I p x], x, pp, FourierParameters -> {0, -1}]
```

The transform is exactly $\delta(p-p')$: a sharp momentum state is a delta spike at its eigenvalue, the
continuum analog of an orthonormal basis vector. The price of perfect momentum definiteness is a state
spread uniformly over all space, outside $L^2$ entirely.

### 2.6 [BSc] How do I obtain the momentum-space wavefunction $\tilde\psi(p)$ by Fourier transform and compute $\langle p^n\rangle$ as a moment in $p$?

The momentum-space wavefunction is $\tilde\psi(p)=\tfrac{1}{\sqrt{2\pi}}\int\psi(x)e^{-ipx}\,dx$ (the
option `FourierParameters -> {0, -1}` fixes the quantum-mechanics sign), after which
$\langle p^n\rangle=\int p^n|\tilde\psi(p)|^2\,dp$ is an ordinary moment. Take the ground state of the
infinite square well $[0,L]$, $\psi(x)=\sqrt{2/L}\,\sin(\pi x/L)$, a state with hard edges rather than
smooth tails. Transform it over its support:

```wl
\[Psi][x_] := Sqrt[2/L] Sin[Pi x/L];
\[Psi]T = FullSimplify[Integrate[\[Psi][x] Exp[-I p x]/Sqrt[2 Pi], {x, 0, L}, Assumptions -> L > 0 && p \[Element] Reals], L > 0]
```

The transform is $\tilde\psi(p)=-\dfrac{(1+e^{-iLp})\sqrt{L\pi}}{L^2p^2-\pi^2}$, peaked near $p=\pm\pi/L$,
the momenta of the two counter-propagating waves that make up the standing state. Confirm it carries unit
probability:

```wl
FullSimplify[Integrate[Abs[\[Psi]T]^2, {p, -Infinity, Infinity}, Assumptions -> L > 0], L > 0]
```

The norm is $1$. Compute the second moment $\langle p^2\rangle=\int p^2|\tilde\psi(p)|^2\,dp$:

```wl
FullSimplify[Integrate[p^2 Abs[\[Psi]T]^2, {p, -Infinity, Infinity}, Assumptions -> L > 0], L > 0]
```

The moment is $\langle p^2\rangle=(\pi/L)^2$, exactly twice the ground-state energy $E_1=\pi^2/(2L^2)$,
now read off in momentum space. Cross-check it against the position-space form
$\int|\partial_x\psi|^2\,dx$, which shares no machinery with the transform:

```wl
Integrate[Abs[D[\[Psi][x], x]]^2, {x, 0, L}, Assumptions -> L > 0]
```

The two routes agree at $(\pi/L)^2$, and $\langle p\rangle=0$ by the symmetry of $|\tilde\psi|^2$. The
hard walls put weight in tails of $\tilde\psi$ that fall only as $1/p^2$: sharp edges in position cost a
broad momentum spread, so the well ground state sits well above minimum uncertainty.

### 2.7 [BSc] How do I verify the position-momentum uncertainty relation $\Delta x\,\Delta p\ge 1/2$ and show that a Gaussian saturates it?

The bound $\Delta x\,\Delta p\ge\tfrac12$ (with $\hbar=1$) is saturated by a Gaussian and by nothing
else, so it is worth computing the product for a non-Gaussian state first and for the Gaussian second.
Take the hyperbolic-secant profile $\psi(x)=\tfrac{1}{\sqrt2}\operatorname{sech}(x)$, the ground state
and only bound state of the well $V(x)=-\operatorname{sech}^2x$, and form $\Delta x\,\Delta p$ from
its second moments:

```wl
\[Psi][x_] := Sech[x]/Sqrt[2];
dxS = Sqrt[Integrate[x^2 Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}]];
dpS = Sqrt[Integrate[Abs[D[\[Psi][x], x]]^2, {x, -Infinity, Infinity}]];
dxS dpS
```

The product is $\Delta x\,\Delta p=\pi/6\approx0.524$, just above the floor: a state can be beautifully
localized and still not minimal. Now the Gaussian $\psi(x)=\pi^{-1/4}\sigma^{-1/2}e^{-x^2/2\sigma^2}$, at
arbitrary width $\sigma$:

```wl
\[Psi]G[x_] := Pi^(-1/4) \[Sigma]^(-1/2) Exp[-x^2/(2 \[Sigma]^2)];
dxG = Sqrt[Integrate[x^2 Abs[\[Psi]G[x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Sigma] > 0]];
dpG = Sqrt[Integrate[Abs[D[\[Psi]G[x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Sigma] > 0]];
Simplify[dxG dpG, \[Sigma] > 0]
```

The Gaussian gives exactly $\tfrac12$ for every width $\sigma$: squeezing it in position inflates it in
momentum and vice versa, so it traces the very edge of the allowed region while every other state, the
sech included, lies strictly above. This is the one place the Gaussian earns its keep, as the unique
minimum-uncertainty state.

### 2.8 [MSc] How do I build a function $f(\hat x,\hat p)$ of noncommuting operators, and resolve the ordering ambiguity by Weyl symmetric ordering?

Because $\hat x$ and $\hat p$ do not commute, the classical product $xp$ has several inequivalent
quantizations. Of the combinations of $\hat x\hat p$ and $\hat p\hat x$, the Weyl (symmetric)
prescription $\tfrac12(\hat x\hat p+\hat p\hat x)$ is the only Hermitian one, up to an additive real
constant. Build the symmetric operator and let it act on a generic function $h(x)$:

```wl
xop[f_] := x f;
pop[f_] := -I D[f, x];
sym[h_] := (xop[pop[h]] + pop[xop[h]])/2;
Simplify[sym[h[x]]]
```

The symmetric operator acts as $-i\big(\tfrac12 h+x h'\big)$, the naive $\hat x\hat p=-i\,x h'$ plus half
the commutator. What the naive ordering lacks is Hermiticity; test the matrix-element difference for the
symmetric operator against the same difference for $\hat x\hat p$, using two convenient decaying test
functions (displaced Gaussians, used only as smooth probes, not as the object of study):

```wl
\[Phi][x_] := Pi^(-1/4) Exp[-(x - 1)^2/2];
\[Chi][x_] := Pi^(-1/4) Exp[-(x + 1)^2/2] Exp[I x];
{Integrate[Conjugate[\[Phi][x]] sym[\[Chi][x]], {x, -Infinity, Infinity}] - Integrate[Conjugate[sym[\[Phi][x]]] \[Chi][x], {x, -Infinity, Infinity}],
 Integrate[Conjugate[\[Phi][x]] xop[pop[\[Chi][x]]], {x, -Infinity, Infinity}] - Integrate[Conjugate[xop[pop[\[Phi][x]]]] \[Chi][x], {x, -Infinity, Infinity}]}
```

The first entry is $0$ and the second is nonzero: only $\tfrac12(\hat x\hat p+\hat p\hat x)$ is Hermitian,
so only it can represent an observable, while $\hat x\hat p$ cannot. The ambiguity between $\hat x\hat p$
and $\hat p\hat x$ is precisely the commutator $[\hat x,\hat p]=i$, now seen as an ordering choice with
physical consequences. Hermiticity settles the ordering here, but it does not do so in general: for
$x^2p^2$ the Weyl and Born-Jordan prescriptions are both Hermitian and differ from each other, so
picking a quantization takes more than requiring a real spectrum.

### 2.9 [MSc] How do I distinguish self-adjointness from formal Hermiticity through the operator domain and boundary conditions (the box versus the half-line)?

Formal Hermiticity is the bulk identity just proved: it needs the boundary term $-i[\phi^*\psi]$ to
vanish, which is a condition on the domain. Self-adjointness demands more, namely that the operator
and its adjoint have the *same* domain, $\mathcal{D}(\hat A^\dagger)=\mathcal{D}(\hat A)$. The two
requirements come apart differently on a box and on a half-line. On the box $[0,L]$, impose the
twisted boundary condition $\psi(L)=e^{i\theta}\psi(0)$ and evaluate the surface term:

```wl
bt = -I (Conjugate[pL] sL - Conjugate[p0] s0) /. {sL -> Exp[I \[Theta]] s0, pL -> Exp[I \[Theta]] p0};
Simplify[bt, \[Theta] \[Element] Reals]
```

The surface term is $0$ for every real $\theta$: each twist angle defines a genuine self-adjoint
momentum operator, a one-parameter family of extensions. Now the half-line $[0,\infty)$, where
$\psi(\infty)=0$ but $\psi(0)$ is free:

```wl
-I (0 - Conjugate[p0] s0)
```

The term reduces to $i\,\phi^*(0)\psi(0)$, and it *can* be made to vanish, by requiring $\psi(0)=0$.
So $\hat p$ is symmetric on the half-line, and the surface term alone does not distinguish this case
from the box. The obstruction lies deeper, in the deficiency indices: solving $\hat p^\dagger\psi=\pm
i\psi$ gives $\psi=e^{\mp x}$, and on $[0,\infty)$ only $e^{-x}$ is square-integrable, so the indices
are $(1,0)$. They are unequal, and von Neumann's theorem then forbids any self-adjoint extension
whatever. Unlike the box, the half-line has no second endpoint whose phase could be twisted to
restore the balance. Hermiticity is a property of a formula; self-adjointness is a property of a
formula together with its domain, and the difference is physical, since only a self-adjoint operator
has the real spectrum and unitary evolution that an observable requires.

## Part 3. The time-independent Schrodinger equation as an ODE eigenvalue problem

The stationary Schrodinger equation $-\tfrac12\psi''+V(x)\psi=E\psi$ is an eigenvalue problem, and it
is the boundary conditions, not the differential equation alone, that decide which energies are
allowed. Where $E$ lies below the potential's value at infinity, the solutions that stay normalizable
exist only at isolated energies, the discrete bound spectrum; above it every energy is permitted and
the solutions are non-normalizable scattering states. Because $\hat H$ is Hermitian, bound
eigenfunctions of different energies are orthogonal and, once normalized, satisfy
$\int\psi_m^*\psi_n\,dx=\delta_{mn}$, while together they are complete, so any state expands in them.
Two structural facts then organize the spectrum, both resting on the non-degeneracy of bound levels
in one dimension: in a symmetric potential $V(-x)=V(x)$ every bound eigenstate has definite parity,
and the $n$-th bound state has exactly $n$ nodes, so counting zeros counts levels. Closed forms exist
only for particular potentials; elsewhere the quantization condition becomes a transcendental
equation to be rooted, or the problem passes to a numerical eigensolver or to direct integration.
Throughout, $\hbar=m=1$.

### 3.1 [BSc] How do I write the stationary Schrodinger equation $-\tfrac12\psi''+V\psi=E\psi$, and what distinguishes a normalizable bound state from a scattering state?

The stationary equation asks which functions the Hamiltonian merely rescales. It is an
eigenvalue problem, and its two kinds of solution are separated not by the equation but by the
energy: below the potential's value at infinity a solution can decay and be normalized, above it
the solution must oscillate forever and cannot. To see both without changing potentials we take
the Poschl-Teller family $V_\lambda(x)=-\tfrac{\lambda(\lambda+1)}{2}\operatorname{sech}^2(x)$, a
whole one-parameter tower of wells rather than a single one, whose $\lambda=1$ member is the well
$V_1(x)=-\operatorname{sech}^2(x)$ with ground state $\operatorname{sech}(x)$.

```wl
Vpt[x_, \[Lambda]_] := -\[Lambda] (\[Lambda] + 1) Sech[x]^2/2;
asl = Assumptions -> \[Lambda] > 0 && Element[x, Reals];
{FullSimplify[-(1/2) D[Sech[x]^\[Lambda], {x, 2}] + Vpt[x, \[Lambda]] Sech[x]^\[Lambda]
   - (-\[Lambda]^2/2) Sech[x]^\[Lambda], asl],
 FullSimplify[Integrate[Sech[x]^(2 \[Lambda]), {x, -Infinity, Infinity}, asl], asl]}
```

The cell returns $\left\{0,\ \tfrac{\sqrt\pi\,\Gamma(\lambda)}{\Gamma(\lambda+\tfrac12)}\right\}$:
for **every** $\lambda>0$ the function $\operatorname{sech}^\lambda(x)$ is an exact bound state
with $E_0=-\lambda^2/2$, and its squared norm is a closed form in $\lambda$. At $\lambda=1$ that norm
is $\sqrt\pi\,\Gamma(1)/\Gamma(\tfrac32)=2$, the familiar $\int\operatorname{sech}^2=2$, which is one
member of a family rather than an isolated integral. The binding deepens as $\lambda^2$ while the
well deepens as $\lambda(\lambda+1)$, so a deeper well binds its ground state more tightly but not
proportionally.

One state is not a spectrum, and $\lambda=1$ is a poor advertisement for this family: it is the
*largest* coupling that binds only one state, since every $\lambda\in(0,1]$ binds exactly one. Push
$\lambda$ past $1$ and a tower appears, closing symbolically in **both** indices at once with no
integer pinned anywhere.

One caution decides which function to write down. A second-order equation has a two-dimensional
solution space, so a vanishing residual says only that a function solves the equation, never that it
is a state; normalizability is the extra condition, and it selects a branch. Here the normalizable
branch is the associated Legendre function of *negative* order, $P_\lambda^{-(\lambda-n)}(\tanh x)$,
because $P_\nu^{\mu}(z)\sim\frac{1}{\Gamma(1-\mu)}\big(\tfrac{1+z}{1-z}\big)^{\mu/2}$ as $z\to1$, so
positive order diverges at $x\to\infty$ except at positive integers, while negative order decays.
Check the two conditions together, the residual and the norm:

```wl
H[\[Lambda]_][f_] := -(1/2) D[f, {x, 2}] + Vpt[x, \[Lambda]] f;
\[Psi]t[\[Lambda]_, n_] := LegendreP[\[Lambda], -(\[Lambda] - n), Tanh[x]];
{FullSimplify[H[\[Lambda]][\[Psi]t[\[Lambda], n]] - (-(\[Lambda] - n)^2/2) \[Psi]t[\[Lambda], n], asl],
 NIntegrate[\[Psi]t[1.3, 0]^2, {x, -30, 30}]}
```

The residual is $0$ and the norm is finite: $\psi_n\propto P_\lambda^{-(\lambda-n)}(\tanh x)$ is an
exact **bound** state at $E_n=-\tfrac{(\lambda-n)^2}{2}$ for arbitrary $\lambda$ and arbitrary $n$,
the spectrum as a theorem rather than as a list. What ends the tower is normalizability and not the
sign of the energy: $E_n=-\tfrac{(\lambda-n)^2}{2}$ is negative for every $n\neq\lambda$, including
$n>\lambda$, but there the order $-(\lambda-n)$ turns positive and the function diverges instead of
decaying. The condition is $\lambda-n>0$, so the well holds $\lceil\lambda\rceil$ states, and at
$\lambda=1$ that count is $1$: a single level, and no tower in sight.

Nothing here is special to integer $\lambda$: $\lambda=\tfrac52$ carries
$\operatorname{sech}^{3/2}x\tanh x$ at $E_1=-\tfrac98$ just as happily. But the *elementary*
members are worth seeing explicitly, since $\lambda=2$ is where a second state first appears.

```wl
{FullSimplify[H[2][Sech[x]^2] - (-2) Sech[x]^2],
 FullSimplify[H[2][Sech[x] Tanh[x]] - (-1/2) Sech[x] Tanh[x]],
 FullSimplify[H[\[Lambda]][Sech[x]^(\[Lambda] - 1) Tanh[x]]
   - (-(\[Lambda] - 1)^2/2) Sech[x]^(\[Lambda] - 1) Tanh[x],
  Assumptions -> \[Lambda] > 1 && Element[x, Reals]]}
```

All three residuals vanish: the well $-3\operatorname{sech}^2x$ holds $\operatorname{sech}^2x$ at
$E_0=-2$ and $\operatorname{sech}x\tanh x$ at $E_1=-\tfrac12$, and the third cell shows the first
excited state is $\operatorname{sech}^{\lambda-1}x\tanh x$ at $-\tfrac{(\lambda-1)^2}{2}$ for every
$\lambda>1$, matching the assumption the cell carries: at $\lambda\le1$ the well holds only its
ground state and there is no excited level to find. That excited state is odd and has one node, which
is parity and the node theorem arriving early.

The coincidence $E_1(\lambda{=}2)=-\tfrac12=E_0(\lambda{=}1)$ is not a coincidence, and noticing it
is not the same as explaining it. The explanation is a factorization: define
$A_\lambda=\tfrac{1}{\sqrt2}\left(\partial_x+\lambda\tanh x\right)$ and its adjoint, and the whole
tower follows.

```wl
A[\[Lambda]_][f_] := (1/Sqrt[2]) (D[f, x] + \[Lambda] Tanh[x] f);
Ad[\[Lambda]_][f_] := (1/Sqrt[2]) (-D[f, x] + \[Lambda] Tanh[x] f);
{FullSimplify[Ad[\[Lambda]][A[\[Lambda]][f[x]]] - (H[\[Lambda]][f[x]] + \[Lambda]^2/2 f[x])],
 FullSimplify[A[\[Lambda]][Ad[\[Lambda]][f[x]]] - (H[\[Lambda] - 1][f[x]] + \[Lambda]^2/2 f[x])],
 FullSimplify[A[\[Lambda]][Sech[x]^\[Lambda]], Assumptions -> \[Lambda] > 0]}
```

Three zeros, each doing separate work. The first is the factorization
$H_\lambda=A_\lambda^\dagger A_\lambda-\tfrac{\lambda^2}{2}$, which alone proves no state can lie
below $-\lambda^2/2$, since $A^\dagger A$ is positive. The second is **shape invariance**: reversing
the product gives the *same family with $\lambda\to\lambda-1$*, so $A_\lambda$ maps the spectrum of
$V_\lambda$ onto that of $V_{\lambda-1}$ and the ladder is forced, not observed. The third says the
ground state is annihilated, $A_\lambda\operatorname{sech}^\lambda x=0$, which is why $-\lambda^2/2$
is attained rather than merely bounded. This is the machinery of 8.5 and 8.6, and the tower above is
its corollary.

Now the other half of the spectrum, at $\lambda=1$ where the scattering state is simplest.

```wl
\[Psi]k[x_] := Exp[I k x] (Tanh[x] - I k);
{FullSimplify[-(1/2) D[\[Psi]k[x], {x, 2}] + Vpt[x, 1] \[Psi]k[x] - (k^2/2) \[Psi]k[x],
  Assumptions -> Element[{k, x}, Reals]],
 {Limit[\[Psi]k[x] Exp[-I k x], x -> -Infinity, Assumptions -> Element[k, Reals]],
  Limit[\[Psi]k[x] Exp[-I k x], x -> Infinity, Assumptions -> Element[k, Reals]]}}
```

```wl
tk = FullSimplify[ComplexExpand[Arg[(1 - I k)/(-1 - I k)], TargetFunctions -> {Re, Im}], k > 0];
d0 = Limit[tk, k -> 0, Direction -> "FromAbove"]; dInf = Limit[tk, k -> Infinity];
{FullSimplify[Abs[(1 - I k)/(-1 - I k)]^2, Element[k, Reals]], tk, {d0, dInf, (d0 - dInf)/Pi}}
```

```wl
tAmp[\[Lambda]_, k_] := Gamma[1 + \[Lambda] - I k] Gamma[-\[Lambda] - I k]/(Gamma[1 - I k] Gamma[-I k]);
refl[\[Lambda]_, k_] := Sin[Pi \[Lambda]]^2/(Sin[Pi \[Lambda]]^2 + Sinh[Pi k]^2);
{FullSimplify[tAmp[1, k] - (1 - I k)/(-1 - I k), Assumptions -> k > 0],
 FullSimplify[1 - Abs[tAmp[\[Lambda], k]]^2 - refl[\[Lambda], k],
  Assumptions -> k > 0 && \[Lambda] > 0],
 N@Table[refl[\[Lambda], 1], {\[Lambda], {1, 2, 3, 1/2, 3/2, 13/10}}]}
```

```wl
{N[Arg[tAmp[13/10, 10.^-8]]/Pi], Ceiling[13/10],
 {Abs[tAmp[13/10, 10.^-8]]^2, Abs[tAmp[1, 10.^-8]]^2},
 {FullSimplify[LegendreP[1, 0, Tanh[x]]],
  FullSimplify[-(1/2) D[Tanh[x], {x, 2}] + Vpt[x, 1] Tanh[x]]}}
```

```wl
dens = FullSimplify[Abs[\[Psi]k[x]]^2, Assumptions -> Element[{k, x}, Reals]];
pn = Integrate[dens, {x, -a, a}, Assumptions -> a > 0 && Element[k, Reals]];
{dens, pn, Limit[pn, a -> Infinity, Assumptions -> Element[k, Reals]]}
```

The residual is $0$ for every real $k$: $\psi_k(x)=e^{ikx}(\tanh x-ik)$ solves the *same* equation
at $E=k^2/2>0$. The difference from the bound state is entirely in the norm. $|\psi_k|^2=k^2+\tanh^2x$
integrates over $[-a,a]$ to exactly $2\left(a+ak^2-\tanh a\right)$, growing linearly in the box and
diverging as $a\to\infty$. No constant can normalize it, so it must be delta-normalized in the sense
of a delta-normalized momentum eigenstate, and the linear growth is why: the natural measure of a continuum state is *per unit length*,
and $1+k^2$ is that density.

The reflectionless statement is done properly by stripping the carrier wave. $\psi_k e^{-ikx}$
tends to the constant $-1-ik$ as $x\to-\infty$ and to $1-ik$ as $x\to+\infty$. Constants, both: on
the left there is **no** $e^{-2ikx}$ piece, so the reflected amplitude is $r=0$ identically, not
merely small. (This is why the limit existing at all is the content: had there been a reflected
piece, $\lim_{x\to-\infty}e^{-2ikx}$ is `Indeterminate` and the cell would not have returned a
number.) Equal asymptotic densities on the two sides would have been necessary but not sufficient.

The transmission is then $t=(1-ik)/(-1-ik)$ with $|t|^2=1$ for every $k$: perfectly transparent at
every energy, all it does is imprint the phase $\delta(k)=\arg t=\arctan(k^2-1,2k)$, which is
$\pi-2\arctan k$. And at $\lambda=1$ that phase counts the bound states: the third entry returns
$\{\pi,\,0,\,1\}$, so $\left(\delta(0)-\delta(\infty)\right)/\pi=1$, matching the
$\lceil\lambda\rceil=1$ counted from the tower. That is **Levinson's theorem**, and the two halves of
this question are not independent facts.

**But the count comes out an integer here only because $\lambda=1$ is exceptional.** In one dimension
Levinson's theorem carries a half-unit: $\delta(0)=\left(n_b-\tfrac12\right)\pi$ generically. At
$\lambda=\tfrac{13}{10}$ the cell returns $\delta(0)/\pi=-\tfrac12$ against $n_b=\lceil1.3\rceil=2$,
and the two reconcile only modulo $2$: `Arg` is principal-valued on $(-\pi,\pi]$, so it pins the phase
only up to $2\pi$, and unwrapping gives $\tfrac32=n_b-\tfrac12$. That caveat is not cosmetic, because
it means this cell by itself cannot count bound states: $\lambda=\tfrac12$ with $n_b=1$ and
$\lambda=\tfrac{12}{5}$ with $n_b=3$ both return $+\tfrac12$. Counting requires following the phase
continuously in $k$ and unwrapping it, not reading one principal value. The integer $\lambda$ are the
exception, and the reason is in the same cell: $|t(0)|^2=1$ at $\lambda=1$ but $10^{-15}$ at
$\lambda=1.3$. Generic wells reflect *totally* at zero energy; integer-$\lambda$ wells transmit
perfectly even there, a zero-energy resonance, and it is that resonance that promotes the half-unit
to a whole one.

The resonance has a name and it has been hiding in our own formulas. At $k=0$ the scattering state
degenerates to $\psi_0=\tanh x$, which the cell confirms is exactly $P_\lambda^{\lambda-n}(\tanh x)$
at $n=\lambda=1$, the top member of the tower, sitting precisely at $E=0$. It solves $H_1\psi=0$
(residual $0$), it is bounded, and it is not normalizable: a **half-bound state**, neither in the
discrete spectrum nor properly in the continuum. That is also why the $k\to0$ limit needed
`Direction -> "FromAbove"`: at $\lambda=1$, $t(0)=-1$ exactly, and `Arg` sits on its branch cut
there. At generic $\lambda$ the limit fails for a different reason entirely, namely $t(0)=0$.

The reflection coefficient of the whole family says why transparency is quantized:
$R_\lambda(k)=\dfrac{\sin^2(\pi\lambda)}{\sin^2(\pi\lambda)+\sinh^2(\pi k)}$, which vanishes
**identically in $k$** exactly when $\sin(\pi\lambda)=0$. Numerically it is $\{0,0,0\}$ at
$\lambda=1,2,3$ against $\{0.0074,0.0074,0.0049\}$ at $\lambda=\tfrac12,\tfrac32,1.3$. That formula
must not be taken on faith, so the cell above ties it to the equation rather than to $\sin(\pi\lambda)$:
the transmission amplitude of the family is
$t=\frac{\Gamma(1+\lambda-ik)\Gamma(-\lambda-ik)}{\Gamma(1-ik)\Gamma(-ik)}$, which reduces to the
$(1-ik)/(-1-ik)$ derived above at $\lambda=1$ (first entry, $0$), and $1-|t|^2-R_\lambda$ vanishes
symbolically for arbitrary $\lambda$ and $k$ (second entry, $0$). Unitarity is the link: without it,
$R$ would be a formula asserted next to some numbers, and a wrong one with $\sinh\to\cosh$ would
sail through the same numerical spot-checks.

So reflectionlessness is not a property of $\operatorname{sech}^2$ wells in general, it is a property
of the ones whose coupling sits at an integer $\lambda$, where a half-bound state sits at threshold.
The $\lambda=1$ well is the first of them, and its single bound state and its transparency are the
same fact about $\operatorname{sech}$ seen from two sides.

### 3.2 [BSc] How do I solve the infinite square well by imposing boundary conditions on the ODE, and read off the spectrum and eigenfunctions?

The infinite well is the one bound-state problem that is exactly solvable with nothing but a
boundary condition, so it is the benchmark every later numerical route is measured against.
Routing it: the class is a second-order linear ODE eigenproblem on a bounded domain with
homogeneous Dirichlet conditions, which admits the exact symbolic eigenproblem (`DEigensystem`),
and exactness is the first criterion, so that is the primary route. It closes, and it closes with
$L$ still symbolic, which no numerical route could give.

```wl
{es, fs} = DEigensystem[{-(1/2) D[u[x], {x, 2}], DirichletCondition[u[x] == 0, True]},
   u[x], {x, 0, L}, 4, Assumptions -> L > 0, Method -> "Normalize"]
```

```wl
DEigenvalues[-(1/2) D[u[x], {x, 2}], u[x], {x, 0, L}, 4, Assumptions -> L > 0]
```

The first cell returns the spectrum $\left\{\tfrac{\pi^2}{2L^2},\tfrac{2\pi^2}{L^2},\tfrac{9\pi^2}{2L^2},\tfrac{8\pi^2}{L^2}\right\}$
with the normalized eigenfunctions $\sqrt{2/L}\,\sin(n\pi x/L)$: the energies scale as $1/L^2$, so
squeezing the box costs energy quadratically, which is the uncertainty principle wearing a
boundary condition.

The second cell is the trap, and it is worth running. Drop the `DirichletCondition` and
`DEigenvalues` does not complain: it returns $\left\{0,\tfrac{\pi^2}{2L^2},\tfrac{2\pi^2}{L^2},\tfrac{9\pi^2}{2L^2}\right\}$,
the *Neumann* spectrum, headed by a spurious zero-energy state (the constant function, which has
zero derivative at both walls and is not an infinite-well state at all). No boundary condition
does not mean no boundary condition; it means $\psi'=0$ at the walls. The physics of the infinite
well is the Dirichlet condition, and supplying it is not a formality.

Four levels are four data points, not a theorem, and `DEigensystem` cannot be asked for a symbolic
count (given $n$ in place of $4$ it returns unevaluated). So the general claim $E_n=n^2\pi^2/(2L^2)$
has to be made by supplying the closed form and verifying the eigenvalue equation as an identity
in $n$.

```wl
\[Psi]n[n_, x_] := Sqrt[2/L] Sin[n Pi x/L];
asn = Assumptions -> L > 0 && Element[n, Integers] && n >= 1;
{FullSimplify[-(1/2) D[\[Psi]n[n, x], {x, 2}] - (n^2 Pi^2/(2 L^2)) \[Psi]n[n, x], asn],
 FullSimplify[Integrate[\[Psi]n[n, x]^2, {x, 0, L}, asn], asn]}
```

The cell returns $\{0,1\}$: $H\psi_n=\tfrac{n^2\pi^2}{2L^2}\psi_n$ and $\langle\psi_n|\psi_n\rangle=1$
for *every* integer $n\ge1$ and every $L>0$, not merely the four the solver listed. The two checks
are independent and both are needed: the eigenvalue equation is invariant under $\psi\to c\psi$, so
it is blind to the $\sqrt{2/L}$ prefactor, and the norm is what pins it down.

### 3.3 [BSc] How do I verify the orthonormality and completeness of the well eigenfunctions?

Orthonormality says the eigenfunctions are a set of perpendicular unit vectors; completeness says
there are enough of them to build anything. The first is an integral, the second is a statement
about a limit, and only the second has real content: a basis can be orthonormal and still miss a
whole direction.

```wl
FullSimplify[Integrate[\[Psi]n[m, x] \[Psi]n[n, x], {x, 0, L},
   Assumptions -> L > 0 && Element[{m, n}, Integers] && m >= 1 && n >= 1 && m != n],
  Assumptions -> Element[{m, n}, Integers] && m != n]
```

```wl
fC[x_] := x^2 (L - x);
cn = FullSimplify[Integrate[\[Psi]n[n, x] fC[x], {x, 0, L}, asn], asn]
```

```wl
res = Sum[cn \[Psi]n[n, x], {n, 1, Infinity}];
Table[With[{lv = RandomReal[{0.5, 4}]},
   With[{xv = RandomReal[{0.01, 0.99}] lv},
    Chop[N[(res - fC[x]) /. {L -> lv, x -> xv}], 10^-10]]], {8}]
```

```wl
{FullSimplify[Sum[cn^2, {n, 1, Infinity}], Assumptions -> L > 0],
 Integrate[fC[x]^2, {x, 0, L}, Assumptions -> L > 0],
 FullSimplify[Sum[(cn /. n -> 2 j - 1)^2, {j, 1, Infinity}], Assumptions -> L > 0]}
```

The first cell returns $0$ for arbitrary integers $m\neq n$, which together with the unit norm of
3.2 is the statement $\langle\psi_m|\psi_n\rangle=\delta_{mn}$, proved symbolically rather than
sampled. The $m\neq n$ assumption is not removable: ask for $\langle\psi_m|\psi_n\rangle$ without
it and the kernel reasons about generic $m,n$ and returns a bare $0$, which would assert that the
states are orthogonal to *themselves* and contradict the unit norm of the well eigenfunctions. The delta is two statements
and they have to be made separately.

For completeness we expand a function that is *not* one of the modes. The obvious choice, the
symmetric parabola $x(L-x)$, is a trap: it is symmetric about the midpoint $L/2$ while the even
modes are antisymmetric there, so every even $c_n$ vanishes and the expansion never tests them.
Take the **asymmetric** cubic $f_C(x)=x^2(L-x)$ instead (named `fC` so it does not shadow the
generic test function $f$ used in the parity and delta-well arguments), whose coefficients do not vanish for any $n$.

Resumming the series is where the symbolic route runs out, and it is worth being
exact about how. `Sum` does return a closed form, but as a combination of trilogarithms, and
`FullSimplify` will not reduce it to $f_C(x)$. Rather than assert an identity the kernel declined
to prove, the spot-check cell takes the documented next rung and evaluates the difference at eight
random pairs $(L,x)$: all eight return $0$. The series does reproduce the cubic; only the
simplifier is beaten. Eight agreeing draws are not a proof, and the honest label for this is
"checked, not proved".

The Parseval cell is the statement that does close, and it is the better one anyway. The identity
$\sum_n|c_n|^2=\int_0^L f_C^2\,dx$ returns $L^7/105$ on both sides, exactly and symbolically in $L$.
Be careful about what that establishes: Parseval holding for one function says the whole weight of
$f_C$ is accounted for, so $f_C$ lies in the closed span of the basis. It does not by itself say the
basis spans everything, since a basis with a mode deleted still satisfies Parseval on any function
having no weight in that mode. Completeness is the statement that this holds for *every* $f\in L^2$,
and it comes from Sturm-Liouville theory rather than from a single test function.

The third entry is what makes it a test rather than an assertion. Summing over the **odd modes
only**, i.e. deliberately deleting every even eigenfunction from the basis, gives $L^7/120$, which
falls short of $L^7/105$ by exactly the weight the deleted modes carried. That is the failure mode
the identity is supposed to detect, exhibited rather than described. It also shows why the choice
of $f_C$ was not cosmetic: run the same deletion against the symmetric parabola $x(L-x)$ and the
odd-only sum still returns the full $L^5/30$, because its even coefficients were zero to begin with.
The check would have passed on a basis missing half its vectors. A completeness test is only as good
as the function you ask it to build, and a symmetric test function cannot see a missing antisymmetric mode.

### 3.4 [BSc] How do I find the finite-square-well bound states from the even/odd transcendental quantization conditions with `FindRoot`?

Make the walls finite and the exactness of the infinite well is gone: the wavefunction leaks into the
barrier, matching it across the edges produces a transcendental equation, and no closed form for the
levels exists. Routing it: the symbolic eigenproblem is the natural first try, and it is **blocked**
here, established by probe rather than by assumption.

```wl
Head@DEigenvalues[{-(1/2) D[u[x], {x, 2}] + Piecewise[{{-12.5, Abs[x] < 1.}}, 0.] u[x],
   DirichletCondition[u[x] == 0, True]}, u[x], {x, -8., 8.}, 4]
```

The head comes back as `DEigenvalues`, in about a third of a second: the expression is handed
straight back untouched. That is the cheapest possible verdict and the most honest kind, because
declining to evaluate is a *loud* failure. Contrast it with 3.5, where the numerical route answers
confidently and wrongly. A route that says nothing has told you more than a route that says
$\{-0.087,0.102,0.152,0.407\}$.

The refusal has a second mode worth knowing about, because it changes what "probe the gate" costs.
Write the same well with the depth and half-width left *symbolic*, which is what criterion 4 would
prefer, and `DEigensystem` does not decline at all: it grinds, still running after 95 seconds. So
the route is closed either way, but only the machine-real form says so quickly. Probing a gate is
cheap only when the answer is "no"; when the answer is "not in any reasonable time" the probe
costs exactly as much as you are willing to wait, and the discipline is to bound it with
`TimeConstrained` and record the bound rather than to sit and hope.

What remains is the transcendental condition solved by `FindRoot`, which is better than a
numerical PDE route anyway because it stays inspectable: the whole problem depends on the single
dimensionless combination $u_0=a\sqrt{2V_0}$, with $a$ the half-width and $V_0$ the depth.

Stay symbolic while the problem still is. Even the bound-state *count* has a closed form, because
even and odd roots interlace at spacing $\pi/2$ across $(0,u_0)$, so pin numbers only after the
count is known. The conditions are written pole-free so the root finder never meets the poles of
$\tan$ and $\cot$.

```wl
fEven[u0_][u_?NumericQ] := u Sin[u] - Sqrt[u0^2 - u^2] Cos[u];
fOdd[u0_][u_?NumericQ] := u Cos[u] + Sqrt[u0^2 - u^2] Sin[u];
bracket[f_, u0_] := With[{us = Subdivide[10.^-9, u0 (1 - 10.^-9), 4000]},
   With[{fs = f /@ us}, Pick[Partition[us, 2, 1], Sign[Most[fs] Rest[fs]], -1]]];
rootsOf[f_, u0_] := (u /. FindRoot[f[u], {u, #[[1]], #[[2]]}, Method -> "Brent"]) & /@ bracket[f, u0];
uAll[u0_] := Sort[Join[rootsOf[fEven[u0], u0], rootsOf[fOdd[u0], u0]]];
nBound[u0_] := Ceiling[2 u0/Pi];
{uAll[5.], nBound[5.]}
```

```wl
Table[Length[uAll[u0]] == nBound[u0], {u0, {0.05, 0.5, 1.5, 2.0, 3.2, 7.9, 12.8, 20.1, 50.3}}]
```

```wl
Esym[u0_, u_, a_] := -(u0^2 - u^2)/(2 a^2);
eFinite = Esym[5., uAll[5.], 1.];
{eFinite,
 Table[With[{u0 = a Sqrt[2 12.5]}, First@Esym[u0, uAll[u0], a]], {a, {1., 2., 4.}}],
 Table[With[{V0 = 12.5/a^2}, a^2 First@Esym[a Sqrt[2 V0], uAll[a Sqrt[2 V0]], a]],
  {a, {1., 2., 4.}}]}
```

```wl
{With[{us = uAll[200.]}, (us[[1 ;; 4]]^2)/(us[[1]]^2)], Length[uAll[0.01]]}
```

At $u_0=5$ the roots are $u\approx\{1.306,2.596,3.837,4.906\}$ and there are $\lceil 2u_0/\pi\rceil=4$
of them, so the closed-form count is right. One value is one data point, so the second cell sweeps
$u_0$ from $0.05$ to $50.3$ and the count formula holds at every one.

The energies keep the half-width symbolic: $E_n=-\left(u_0^2-u_n^2\right)/(2a^2)$. At $a=1$,
$V_0=12.5$ the levels are $\{-11.65,-9.13,-5.14,-0.46\}$, all lying strictly between $-V_0$ and
$0$, as a bound state of a well of depth $V_0$ must. Naming them `eFinite` also means 3.5 compares
against one computation of the roots rather than re-running the bracketing search each time.

The other two entries are there to kill a scaling law that looks obvious and is wrong. It is
tempting to read $E_n=-\left(u_0^2-u_n^2\right)/(2a^2)$ as "the spectrum scales as $1/a^2$", the
same squeeze 3.2 found for the infinite well. It does not. Widen the well at fixed depth $V_0=12.5$
and the ground level goes $\{-11.65,\,-12.25,\,-12.43\}$ for $a=1,2,4$: it *deepens* toward $-V_0$,
where $1/a^2$ would demand $\{-11.65,\,-2.91,\,-0.73\}$. The reason is that $u_0=a\sqrt{2V_0}$
carries its own $a$, so the roots $u_n$ move too, and a wider well simply holds its states more
securely. The document's own 3.9 already contradicts the naive reading: shrinking $a$ at fixed
$g=2aV_0$ sends $E\to-g^2/2$, a constant, not $1/a^2$.

The true statement is the third entry: $a^2E_n$ depends on $a$ and $V_0$ **only through $u_0$**.
Holding $u_0=5$ fixed while $a$ runs over $1,2,4$ (with $V_0$ compensating) returns
$-11.646607252$ every time, identical to ten digits. That is what "the problem depends on one
dimensionless parameter" means, and it is falsifiable, which the $1/a^2$ story was not.

The last cell closes the two limits that could fail. Deep well: at $u_0=200$ the ratios of
$u_n^2$ come back as $\{1.000000,3.999999,8.999993,15.999976\}$, which is $\{1,4,9,16\}$ to five
significant figures (the worst, $n=4$, is off by $1.5\times10^{-6}$ relative, since the highest
state is the least deeply bound and so the slowest to reach its limit). That is the infinite-well
infinite-well spectrum for a box of width $2a$, recovered as $V_0\to\infty$ exactly as the roots crowd
toward $u_n\to n\pi/2$. Shallow well: at $u_0=0.01$ there is still exactly one root. However weak the attraction, a one-dimensional
well always binds, which is a genuinely quantum statement with no classical counterpart and is the
same fact the delta well makes exact.

### 3.5 [BSc] How do I compute numerical bound states with `NDEigensystem` on a finite interval and compare them to the analytic answer?

This is the question the whole routing discipline exists for, so it is worth doing slowly. The
numerical route must (a) truncate the infinite line to a box and (b) discretize it, and both
truncations are lies that must be shown to stop mattering. Worse, the default is wrong here, and
it is wrong silently.

```wl
well[x_] := Piecewise[{{-12.5, Abs[x] < 1.}}, 0.];
NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + well[x] u[x], DirichletCondition[u[x] == 0., True]},
 u[x], {x, -8., 8.}, 4]
```

```wl
ndAll[box_, mcm_] := Sort@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + well[x] u[x],
    DirichletCondition[u[x] == 0., True]}, u[x], {x, -box, box}, 4,
   Method -> {"SpatialDiscretization" -> {"FiniteElement",
      {"MeshOptions" -> {MaxCellMeasure -> mcm}}},
     "Eigensystem" -> {"Arnoldi", "Shift" -> -12.5}}];
Table[ndAll[8., mcm], {mcm, {0.5, 0.1, 0.02, 0.005}}]
```

```wl
Table[ndAll[box, 0.005], {box, {4., 6., 8., 12.}}]
```

```wl
Max@Abs[(ndAll[12., 0.005] - eFinite)/eFinite]
```

The first cell returns $\{-0.087,\,0.102,\,0.152,\,0.407\}$ and **not one of those is a bound state
of this well**. The true levels are the $\{-11.65,-9.13,-5.14,-0.46\}$ found from the transcendental conditions. `NDEigensystem`
returns the eigenvalues smallest in *magnitude*, and for a spectrum that runs negative that is the
wrong end entirely: what came back are discretized continuum states of the surrounding box. There
is no error, no warning, and the numbers look perfectly reasonable. This is the documented
behavior, not a bug, and it is why the physics has to supply a spectral shift near the expected
eigenvalue.

With the shift in place the mesh sweep converges monotonically to the analytic answer, and the box
sweep tells a physical story worth reading rather than skipping: the three deep levels are already
converged at a box of $\pm4$, but the shallowest state at $-0.46$ keeps moving until the box
reaches $\pm12$. That state is the most weakly bound, so its exponential tail $e^{-\kappa|x|}$ has
the smallest $\kappa$ and reaches furthest; squeeze the box and you squeeze the state, raising its
energy. The last cell measures the final disagreement against 3.4 rather than eyeballing it: at a
box of $\pm12$ and a mesh of $0.005$ the two routes agree to $7\times10^{-9}$ relative, and they
share no machinery, one being a transcendental root of $\tan$ and $\cot$ and the other a finite
element discretization of the operator.

A message does fire during the mesh sweep, `Eigensystem::chnpdef`, warning that positive
definiteness cannot be assured for the shifted Arnoldi solve. It is reported rather than
suppressed: it is a warning about the method, and the agreement with 3.4 to seven digits is the
evidence that it did not spoil the answer here.

### 3.6 [BSc] How do I implement the shooting and Numerov methods from scratch to find a bound state?

A canned eigensolver hides the fact that a bound state is a *coincidence*: integrate the ODE
outward from one wall at an arbitrary energy and the solution generically fails to satisfy the far
boundary condition (on an unbounded domain it diverges outright in the classically forbidden region),
and only at an eigenvalue does it land on zero. Shooting makes the eigenvalue condition explicit
as a root of the endpoint value. Numerov is the sharpened version, a three-term recurrence exact
to fourth order in the step, obtained by using the ODE itself to eliminate the leading error term.
Both are written here with `NestList` rather than a loop.

The recurrence relates three consecutive points, so the natural WL object is not an index but the
sliding window `Partition[..., 3, 1]`, folded over. Written that way there is no counter, no `Part`
offset arithmetic, and no chance of an off-by-one. The window carries both $k^2$ and
$f=1+h^2k^2/12$ at each point, which matters: eliminating $k^2$ in favour of $f$ alone would force
the middle coefficient to be computed as $12-10f$, and subtracting two numbers near $12$ and $10$
to get a small one is exactly the cancellation that costs digits. Keeping $k^2$ lets it be written
as $2(1-5h^2k^2/12)$ and computed directly.

```wl
numerov[En_?NumericQ, v_, {x0_, x1_}, npts_Integer] := With[{h = (x1 - x0)/(npts - 1)},
   With[{ks = 2 (En - v[x0 + # h]) & /@ Range[0, npts - 1]},
    With[{w = Partition[Transpose[{ks, 1 + h^2 ks/12}], 3, 1]},
     Prepend[Last /@ FoldList[
        With[{fm = #2[[1, 2]], k0 = #2[[2, 1]], fp = #2[[3, 2]]},
          {#1[[2]], (2 (1 - 5 h^2 k0/12) #1[[2]] - fm #1[[1]])/fp}] &,
        {0., 1.*^-6}, w], 0.]]]];
shoot[En_?NumericQ] := Last[numerov[En, 0. &, {0., 1.}, 2001]];
{Length@numerov[4.9348, 0. &, {0., 1.}, 201],
 Table[{n, n^2 Pi^2/2., u /. FindRoot[shoot[u], {u, n^2 Pi^2/2. - 0.4, n^2 Pi^2/2. + 0.4},
     Method -> "Brent"]}, {n, 1, 4}]}
```

The `Prepend` restores $\psi_0=0$, the wall value the fold consumes as its seed rather than emits,
so the returned list is $\psi_0,\dots,\psi_{N-1}$ and pairs with the grid. Nothing above needs it,
since only `Last` is ever taken, but a reader who plots the eigenfunction would otherwise silently
misalign every point by one step. `Length` returning $201$ for a $201$-point grid is the contract.

```wl
shootND[En_?NumericQ] := NDSolveValue[{u''[x] == -2 En u[x], u[0] == 0, u'[0] == 1},
   u[1], {x, 0, 1}];
Table[{n, n^2 Pi^2/2., v /. FindRoot[shootND[v], {v, n^2 Pi^2/2. - 0.4, n^2 Pi^2/2. + 0.4},
    Method -> "Brent"]}, {n, 1, 4}]
```

```wl
shootN[En_?NumericQ, np_Integer] := Last[numerov[En, 0. &, {0., 1.}, np]];
Table[{np, Abs[(u /. FindRoot[shootN[u, np], {u, Pi^2/2. - 0.4, Pi^2/2. + 0.4},
     Method -> "Brent"]) - Pi^2/2.]}, {np, {51, 101, 201, 401, 801}}]
```

Run on the infinite well of width $1$, whose exact levels are $E_n=n^2\pi^2/2$, Numerov returns all
four to about $3.4\times10^{-10}$, evenly across the four. That evenness is not the scheme tracking
$h^4$ across the levels: the refinement sweep below shows this grid already sits on the roundoff
floor, where a *finer* grid does slightly worse, so the flat $3.4\times10^{-10}$ is arithmetic noise
rather than truncation error. `NDSolve` shooting, an
independent scheme sharing none of Numerov's machinery, is markedly worse and degrades with the
level: $5\times10^{-8}$ on the ground state but $2.2\times10^{-6}$ by $n=4$, since a more
oscillatory state is harder for a general-purpose stepper to track across the interval. The two
numerical routes therefore disagree with *each other* by about $2\times10^{-6}$, four orders more
than Numerov disagrees with the exact answer. That ordering is the reason to run both: Numerov is
fourth order by construction on exactly this equation, `NDSolve` is general, and the exact answer
is what reveals which of the two to believe. Had the analytic benchmark been unavailable, their
mutual agreement to $2\times10^{-6}$ would have been the only accuracy claim on offer, and it would
have understated Numerov by four orders.

The refinement sweep is the part that cannot be skipped. The ground-state error falls
$3.2\times10^{-7}$, $2.0\times10^{-8}$, $1.3\times10^{-9}$, $8.9\times10^{-11}$, $3.6\times10^{-11}$
as the grid doubles, which is ratios of $16.0$, $16.0$, $14.0$ and $2.5$. The first two are $h^4$
convergence measured rather than asserted: halve the step, cut the error sixteenfold. The third,
$14.0$, is already $12\%$ short of $16$, and the fourth abandons the pattern entirely. That decay
of the decay is the roundoff floor arriving: below about $10^{-10}$ the truncation error the scheme
controls has fallen beneath the arithmetic noise it does not, and refining the grid stops buying
accuracy. A scheme's order is a claim about a sequence, and a single run cannot express it; nor,
as the last two ratios show, can a sequence taken too far.

The infinite well is the right *benchmark* and the wrong *test*, because $V\equiv0$ makes $k^2$
constant and the recurrence degenerate: the weights $f=1+h^2k^2/12$ never vary, which is precisely
the thing Numerov exists to track. A scheme tuned to a varying $k^2(x)$ has not been exercised by a
problem with no potential. So run the same `numerov` on a real well, the Poschl-Teller $\lambda=2$
whose exact levels $\{-2,-\tfrac12\}$ we already know.

```wl
shootPT[En_?NumericQ] := Last[numerov[En, -3 Sech[#]^2 &, {-8., 8.}, 4001]];
{SameQ @@ (2 (0. - (-3 Sech[#]^2 &)[-8. + # 0.004]) & /@ Range[0, 5]),
 {u /. FindRoot[shootPT[u], {u, -2.3, -1.7}, Method -> "Brent"],
  u /. FindRoot[shootPT[u], {u, -0.8, -0.3}, Method -> "Brent"]}}
```

The first entry is `False`: on this well $k^2$ genuinely varies from point to point, unlike the
$V\equiv0$ case, so the $f$ weights are doing their job. The levels come back as
$\{-2.0000000000064,\ -0.4999986\}$ against the exact $\{-2,-\tfrac12\}$. Two things are
worth reading off. The ground state is recovered to $6\times10^{-12}$, better than the infinite
well managed, because the state is smooth and dies exponentially inside the box. The excited state
is only good to $1.4\times10^{-6}$ and `FindRoot` reports `brmp`, that it has bracketed the root as
tightly as machine precision allows: at $E=-\tfrac12$ the decay constant $\kappa=\sqrt{2|E|}=1$ is
the slowest of the two, so on $[-8,8]$ that state has not finished decaying, and shooting outward
through a long classically forbidden region amplifies the growing solution exponentially. That is
the real limitation of naive shooting for weakly bound states, and it is why 3.5's box sweep found
the same state to be the fussy one.

One trap here is worth more than the method. Guard `numerov` with `NumericQ` but hand
`FindRoot` the expression `Last[numerov[En, ...]]` directly, and the symbolic $En$ leaves the call
unevaluated, so `Last` returns the last *argument*, the integer `2001`. `FindRoot` then dutifully
searches for a root of the constant $2001$, reports a singular Jacobian, and returns its own
starting guess unchanged. Since the guess was near the answer, the output looks like a converged
eigenvalue and is off by exactly the offset you seeded. The fix is to guard the composite, so
`shoot[En]` stays a single unevaluated head until `FindRoot` supplies a number.

### 3.7 [BSc] How do I show that the eigenstates of a symmetric potential have definite parity?

Parity is the cleanest symmetry argument in quantum mechanics: if the potential is even then the
parity operator $P\psi(x)=\psi(-x)$ commutes with $H$, and a nondegenerate eigenstate of $H$ must
then be an eigenstate of $P$ too, with eigenvalue $\pm1$ since $P^2=1$. The whole statement is an
operator identity, so no solver is needed.

The potential here must be a genuinely *unknown* symbol, so it is called $W$: reusing the `V` of
3.1 would silently supply the Poschl-Teller well, which is already even, and the commutator would
collapse to zero for a reason having nothing to do with the argument being made.

```wl
comm = Simplify[(-(1/2) D[f[-x], {x, 2}] + W[x] f[-x])
    - ((-(1/2) D[f[x], {x, 2}] + W[x] f[x]) /. x -> -x)]
```

```wl
{Simplify[comm /. W[-x] -> W[x]],
 Simplify[(-(1/2) D[f[-x], {x, 2}] + Vpt[x, 1] f[-x])
   - ((-(1/2) D[f[x], {x, 2}] + Vpt[x, 1] f[x]) /. x -> -x)]}
```

The first cell returns $f(-x)\left(W(x)-W(-x)\right)$ for a *generic* $W$, which is the useful form:
the commutator is not zero in general, and it fails by exactly the odd part of the potential. This
matters, because a check that returns zero for every $W$ would prove nothing about symmetry, and it
is the easiest way to fool yourself here. The second cell then imposes $W(-x)=W(x)$ and gets $0$,
and repeats it on the concrete even potential $V_1(x)=-\operatorname{sech}^2(x)$, with no
hand substitution, again $0$.

So $[H,P]=0$ precisely when $V$ is even. The consequence is the physics: if $H\psi=E\psi$ then
$H(P\psi)=E(P\psi)$, so $P\psi$ is an eigenstate at the same energy; if that level is nondegenerate
then $P\psi$ can only be a multiple of $\psi$, and $P^2=1$ forces the multiple to be $\pm1$.

That conclusion is a claim about states, so it should be checked on states rather than left as an
inference. 3.4 already *assumed* it by splitting the finite well into separate even and odd
quantization conditions; here is the assumption tested, on the well's own numerical eigenfunctions.

```wl
{fvals, ffuns} = NDEigensystem[{-(1/2) D[u[x], {x, 2}] + well[x] u[x],
    DirichletCondition[u[x] == 0., True]}, u[x], {x, -12., 12.}, 4,
   Method -> {"SpatialDiscretization" -> {"FiniteElement",
      {"MeshOptions" -> {MaxCellMeasure -> 0.005}}},
     "Eigensystem" -> {"Arnoldi", "Shift" -> -12.5}}];
Table[With[{g = ffuns[[Ordering[fvals]]][[j]]},
   Round[NIntegrate[g (g /. x -> -x), {x, -12., 12.}]/NIntegrate[g^2, {x, -12., 12.}], 0.001]],
  {j, 1, 4}]
```

The expectation $\langle P\rangle$ comes back as $\{1,-1,1,-1\}$: the finite well's four bound states
well alternate strictly between even and odd, starting even, with no state of mixed parity and no
value between $\pm1$. So the even/odd split of the quantization conditions was not a convenience, it was the spectrum's
actual structure, and the numerical solver, which was told nothing about symmetry, discovered it
anyway.

The nondegeneracy clause is load-bearing rather than decorative: in one dimension bound states are
nondegenerate, so the conclusion is safe here, but drop it and the argument fails, since in a
degenerate level one can always form combinations of no definite parity. That is why the same
$\langle P\rangle=\pm1$ result cannot be expected of, say, the 2D box's degenerate levels (8.5 of
the routing reference makes exactly this point about eigenvectors in a degenerate subspace).

### 3.8 [MSc] How do I verify the node theorem, that the $n$-th bound state has exactly $n$ nodes?

The node theorem orders the spectrum by shape alone: count the interior zeros of a bound state and
you know its place in the spectrum, no energies needed. It is a theorem about a whole spectrum, so
checking it on one state is no check; it has to be checked as a statement in $n$, and it has to be
checked somewhere the nodes are not already obvious.

```wl
Solve[Sin[n Pi x/L] == 0 && 0 < x < L, x, Reals,
 Assumptions -> Element[n, Integers] && n >= 1 && L > 0]
```

```wl
FullSimplify[Sin[n Pi (kk L/n)/L], Assumptions -> Element[{n, kk}, Integers]]
```

For the infinite well the whole node family comes back at once, symbolically: `Solve` returns
$x=2LC_1/n$ subject to $C_1\in\mathbb{Z}$, $C_1\ge1$, $n>2C_1$, together with the odd companion.
Those constraints are the statement $x=kL/n$ for $k=1,\dots,n-1$, so the well's $n$-th eigenfunction
$\sin(n\pi x/L)$ has exactly $n-1$ interior nodes for **arbitrary** $n$ and **symbolic** $L$, not for
five values of $n$ at $L=1$. Counting levels from the ground state as $n=0$, as the oscillator and
the finite well below both do, that is the same law: level $n$ carries $n$ nodes. The second cell
confirms the zeros are where the constraints say, again for arbitrary integers.

```wl
Table[CountRoots[HermiteH[n, y], y], {n, 0, 6}]
```

```wl
nodesTol[f_, tol_] := With[{v = Table[f, {x, -11.9, 11.9, 0.005}]},
   With[{keep = Pick[v, UnitStep[Abs[v] - tol Max[Abs[v]]], 1]},
    Count[Sign[Most[keep] Rest[keep]], -1]]];
nodes[f_] := nodesTol[f, 1.*^-6];
{nodes /@ ffuns[[Ordering[fvals]]], nodes[Sin[6 Pi x/24]], nodes[1. + 0 x]}
```

```wl
Table[nodesTol[#, tol] & /@ ffuns[[Ordering[fvals]]], {tol, {0, 10.^-14, 10.^-6, 0.5}}]
```

The oscillator makes the same count with a different potential and a different special function:
$\psi_n\propto H_n(y)e^{-y^2/2}$ has zeros exactly where the Hermite polynomial does, and
`CountRoots` gives $\{0,1,2,3,4,5,6\}$ for $n=0,\dots,6$ (it counts *real* roots, so this can fail:
$y^6+1$ has degree $6$ and zero real roots).

But $\sin$ and $H_n$ both wear their nodes on their sleeve, and a theorem checked only where the
answer is visible by inspection has not been tested. The last cell counts nodes on the finite well's
numerical eigenfunctions, a potential whose **levels** have no closed form (its eigenfunctions are
elementary, but only once the transcendental $u_n$ is in hand, so the nodes still have to be counted
numerically), and gets $\{0,1,2,3\}$: the same law where the answer cannot simply be read off.

The mask in `nodes` is doing real work rather than tidying: outside the well the numerical state
decays to $10^{-16}$ and its *sign* there is pure roundoff, so a naive sign-change count returns
$153$ nodes for the nodeless ground state. Discarding samples below $10^{-6}$ of the peak keeps
only where the wavefunction actually is.

A tolerance introduced to make an answer come out right deserves suspicion, so the last cell
sweeps it. At $\text{tol}=0$ the count is the nonsense $\{153,16,2,3\}$; at every threshold from
$10^{-14}$ to $0.5$ it is $\{0,1,2,3\}$. That is a plateau spanning fourteen orders of magnitude,
which is what distinguishes a principled cut from a fudge factor: the answer does not depend on
where the cut is put, only on there being one. The two controls confirm the counter is not rigged to
return the level index: it gives $5$ on a five-zero function and $0$ on a constant.

The law being tested has hypotheses, and they are what set its reach. Sturm's oscillation theorem
applies to bound states of a one-dimensional Sturm-Liouville problem, ordered by increasing energy,
and needs those levels non-degenerate. Drop any of them and the count fails: in a two-dimensional box
the degenerate levels can be recombined at will, and their nodal patterns are lines whose count is
not fixed by the level index at all. So four confirming instances are not what gives the theorem its
scope; the hypotheses are.

Three exactly solvable problems and one whose levels are transcendental, sharing one counting law:
that is evidence the node theorem is about the structure of second-order eigenproblems, not about $\sin$
or $H_n$.

### 3.9 [BSc] How do I find the single bound state of an attractive delta-function potential?

The delta well is the sharpest bound-state problem there is: a potential with no width, one
parameter, and exactly one bound state no matter how weak. It is also where the symbolic route
runs out and has to be replaced by physics rather than by a bigger solver. `DSolve` on
$-\tfrac12\psi''-g\,\delta(x)\psi=E\psi$ returns unevaluated, so the route is closed; and
`NIntegrate` cannot help either, since it misses point measures entirely.

```wl
{Integrate[DiracDelta[x] f[x], {x, -e, e}, Assumptions -> e > 0],
 NIntegrate[DiracDelta[x] Exp[-x^2], {x, -1, 1}]}
```

The exact `Integrate` picks out $f(0)$; `NIntegrate` returns `0.` and says so with an `izero`
message. Any quadrature-based route to this problem is dead on arrival, which is the content of
the trap rather than a detail of this example.

So integrate the equation itself across the origin. Over $[-\epsilon,\epsilon]$ the $\psi''$ term
gives the jump $\psi'(0^+)-\psi'(0^-)$, the delta term gives $-2g\psi(0)$, and the $E\psi$ term
vanishes as $\epsilon\to0$ because $\psi$ is bounded. Away from the origin the potential is zero and a
bound state ($E<0$) must decay, so the general candidate is $Ae^{\kappa x}$ for $x<0$ and
$Be^{-\kappa x}$ for $x>0$ with $E=-\kappa^2/2$. The delta puts its jump in $\psi'$, not in $\psi$, so
$\psi$ stays continuous at the origin and $A=B$: every bound state of this potential is automatically
even, and there is no odd sector to check separately. That leaves
$\psi=\sqrt{\kappa}\,e^{-\kappa|x|}$, continuous by construction, and imposing the jump fixes
$\kappa$.

```wl
\[Psi]R[x_] := Sqrt[k] Exp[-k x];
\[Psi]L[x_] := Sqrt[k] Exp[+k x];
Solve[{\[Psi]R'[0] - \[Psi]L'[0] == -2 g \[Psi]R[0], k > 0}, k, Reals]
```

```wl
{-k^2/2 /. First@Solve[{\[Psi]R'[0] - \[Psi]L'[0] == -2 g \[Psi]R[0], k > 0}, k, Reals],
 Integrate[k Exp[-2 k Abs[x]], {x, -Infinity, Infinity}, Assumptions -> k > 0]}
```

The jump condition gives $\kappa=g$, so the bound state is $\psi(x)=\sqrt{g}\,e^{-g|x|}$ with
$E=-g^2/2$ and unit norm, the cusped exponential $\sqrt{\kappa}\,e^{-\kappa|x|}$ at $\kappa=g$. The
kink is now explained rather than merely observed: the cusp is the wavefunction's response to a
potential concentrated at a point, and the jump in $\psi'$ *is* the delta.

The constraint $k>0$ in the `Solve` is doing real work. Drop it and the equation
$-2\kappa\sqrt{\kappa}=-2g\sqrt{\kappa}$ also admits the root $\kappa=0$, which `Solve` returns
first; taking it silently yields $E=0$ and no bound state at all.

The physics is in the counting. There is exactly one root $\kappa=g$ for every $g>0$, so the delta
well binds exactly one state however weak the attraction, and its energy $-g^2/2$ is *quadratic*
in the coupling, so a weakly bound state is bound very weakly indeed.

This is the zero-width limit of a finite square well made exact, and that is a quantitative claim, so
it should be computed. Shrink the half-width $a$ while deepening the well at fixed $g=2aV_0$, and the
finite well's single surviving level must approach $-g^2/2$.

```wl
Table[{a, With[{u0 = a Sqrt[2 (1/(2 a))]}, -(u0^2 - First[uAll[u0]]^2)/(2 a^2)]},
 {a, {0.5, 0.1, 0.02, 0.004}}]
```

With $g=1$ fixed the finite well's ground level runs $-0.308$, $-0.442$, $-0.487$, $-0.497$ as the
half-width shrinks $0.5\to0.004$, converging to the delta well's exact $-g^2/2=-0.5$. The approach
is slow, and visibly linear in $a$ rather than quadratic, which is the finite well remembering it has
a width. So the delta well's cusped bound state is not merely analogous to the finite well's last
level, it *is* its zero-width limit, and two independent routes (a transcendental root and a jump
condition) agree in the limit where both apply.

### 3.10 [MSc] How do I solve the linear (triangular) potential in terms of Airy functions, the quantum bouncer?

A particle above a floor in a uniform field, $V(x)=Fx$ for $x>0$ with a hard wall at $x=0$, is the
quantum bouncer, and it is the cleanest demonstration that **the differential equation quantizes
nothing; the boundary condition does.** Routing it: the class is a second-order linear ODE with a
linear coefficient, whose exact symbolic solution is a classical special function, so `DSolve` is
the primary route and the quantization is imposed afterwards by hand. `NDEigensystem` on a
truncated domain is the independent cross-check.

```wl
DSolveValue[-(1/2) u''[x] + F x u[x] == En u[x], u[x], x]
```

```wl
FullSimplify[-(1/2) D[AiryAi[(2 F)^(1/3) (x - En/F)], {x, 2}]
   + F x AiryAi[(2 F)^(1/3) (x - En/F)] - En AiryAi[(2 F)^(1/3) (x - En/F)],
 Assumptions -> F > 0]
```

`DSolve` returns $C_1\operatorname{Ai}\!\left(\tfrac{2(Fx-E)}{(2F)^{2/3}}\right)+C_2\operatorname{Bi}\!\left(\cdots\right)$,
and the second cell confirms the residual is $0$ for **arbitrary** $E$ and $F$. That is the whole
point: every energy solves the equation. $\operatorname{Bi}$ is discarded because it grows without
bound as $x\to\infty$, which is normalizability, not the equation, doing the work. What is left is
the wall at $x=0$, and only there does the spectrum appear: $\psi(0)=0$ forces the argument of
$\operatorname{Ai}$ at the wall to be a zero of $\operatorname{Ai}$.

```wl
Ebounce[n_, F_] := -(F^2/2)^(1/3) AiryAiZero[n];
eAiry = Table[N[Ebounce[n, 1]], {n, 1, 4}];
{eAiry, Sort@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + x u[x], DirichletCondition[u[x] == 0, True]},
   u[x], {x, 0, 25}, 4]}
```

```wl
Table[{box, Max@Abs[(Sort@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + x u[x],
        DirichletCondition[u[x] == 0, True]}, u[x], {x, 0, box}, 4,
       Method -> {"SpatialDiscretization" -> {"FiniteElement",
          {"MeshOptions" -> {MaxCellMeasure -> 0.002}}}}] - eAiry)/eAiry]},
 {box, {10., 15., 20., 25.}}]
```

The quantization condition $\operatorname{Ai}\!\left(-(2F)^{1/3}E/F\right)=0$ gives the closed form
$E_n=-\left(F^2/2\right)^{1/3}a_n$ with $a_n$ the $n$-th (negative) zero of $\operatorname{Ai}$,
symbolic in the field strength $F$. At $F=1$ the first four levels are
$\{1.8558,3.2446,4.3817,5.3866\}$.

The comparison cell is worth reading rather than nodding at. `NDEigensystem` at its *default* mesh
returns $\{1.8601,3.2568,4.4064,5.4857\}$, which agrees with the closed form only to about two
significant figures, and the disagreement grows with the level: $0.2\%$ on the ground state but
$1.8\%$ by $n=4$. Reported alone, those numbers would look like confirmation to anyone reading them
to three digits. They are not; they are an unconverged mesh.

The sweep is what earns the agreement. With the mesh refined to $0.002$ the relative disagreement
drops to $2\times10^{-10}$ at a box of $10$ and to $3\times10^{-12}$ by a box of $20$, after which
it stops improving. The plateau is physical: the bouncer's potential rises forever, so unlike 3.5
there is no continuum to leak into and the truncation at $x=\text{box}$ is harmless *provided* the
box sits well beyond the classical turning point $x_n=E_n/F$ of the highest state wanted. Once
every state requested is comfortably inside its turning point, the remaining error is the mesh, not
the box, and growing the box further buys nothing. The $F^{2/3}$ scaling of the levels is visible in the closed form and is the
signature the experiment measures: it is why bouncing ultracold neutrons in the Earth's
gravitational field have level spacings of order peV.

So far only energies. A spectrum is not a state, and the brief's invariant is unit norm on every
bound state, so the eigenfunction has to be built and normalized too. It closes in closed form.

```wl
{NIntegrate[AiryAi[2^(1/3) x + AiryAiZero[1]]^2, {x, 0, Infinity}],
 N[AiryAiPrime[AiryAiZero[1]]^2/2^(1/3)],
 Table[AiryAi[AiryAiZero[n]], {n, 1, 3}]}
```

The normalization integral $\int_0^\infty\operatorname{Ai}^2$ evaluates to
$\operatorname{Ai}'(a_n)^2/(2F)^{1/3}$, and quadrature agrees at $0.39026$: the normalized bouncer
state is $\psi_n(x)=\operatorname{Ai}\!\left((2F)^{1/3}x+a_n\right)\big/\left[\operatorname{Ai}'(a_n)/(2F)^{1/6}\right]$,
with no numerical integration needed. The third entry is the wall condition itself,
$\psi_n(0)=\operatorname{Ai}(a_n)=0$ exactly, for every $n$ and every $F$: the boundary condition
that produced the spectrum is satisfied identically, which is the check that the quantization was
imposed rather than assumed.

Finally the limit that could fail. Restoring units,
$E_n=\left(\tfrac{\hbar^2F^2}{2m}\right)^{1/3}\left[\tfrac{3\pi(n-\tfrac14)}{2}\right]^{2/3}$, so
$\hbar$, $F$ and $m$ enter only through the prefactor and cancel entirely from the *ratio*
$E_n^{\mathrm{WKB}}/E_n$. The semiclassical parameter is the classical action measured in units of
$\hbar$, so shrinking $\hbar$ at fixed $n$ does not improve the relative WKB error at all; the limit
with content is $n\to\infty$, where the levels crowd and the WKB estimate must become exact. Bohr-Sommerfeld for a
hard wall plus a linear rise gives $E_n^{\mathrm{WKB}}=\left(\tfrac{3\pi(n-\tfrac14)}{2}\right)^{2/3}\left(\tfrac{F^2}{2}\right)^{1/3}$.

```wl
Table[{n, N@Ebounce[n, 1], N[(3 Pi (n - 1/4)/2)^(2/3) (1/2)^(1/3)],
   N@Abs[(3 Pi (n - 1/4)/2)^(2/3) (1/2)^(1/3)/Ebounce[n, 1] - 1]}, {n, {1, 5, 20, 100}}]
```

The relative error falls $7.6\times10^{-3}$, $2.1\times10^{-4}$, $1.2\times10^{-5}$,
$4.7\times10^{-7}$ at $n=1,5,20,100$: WKB is already good to a percent on the *ground* state and
becomes exact as $n\to\infty$, the correspondence principle with content. The ratio
$E_n^{\mathrm{WKB}}/E_n$ is independent of $F$ (both sides carry the same
$\left(F^2/2\right)^{1/3}$), so this is a statement in $n$ alone, and it is why the $-\tfrac14$
in the Maslov index matters: drop it and the $n=1$ error is an order of magnitude worse.

Note that the exact symbolic eigenproblem is not simply unavailable here. `DEigensystem` on the
bouncer in a box does return, as `Root` objects built from $\operatorname{Ai}$ and
$\operatorname{Bi}$ that agree with $E_n$ above, but it emits internal `ReplaceAll::reps` messages
while doing so. `DSolve` plus the Airy zero is the same physics obtained cleanly, and it keeps $F$
symbolic, which the eigensolver route does not.
