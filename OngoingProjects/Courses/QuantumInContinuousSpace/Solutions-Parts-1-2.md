---
Template: Default
---

# Quantum in Continuous Space: Worked Solutions (Parts 1 and 2)

Worked, pure Wolfram Language answers to the questions in
[`Question-List.md`](Question-List.md), Part 1 (the wavefunction and the state space $L^2$)
and Part 2 (operators, position and momentum, and the Fourier bridge). No quantum package is
used: everything is `Integrate`, `D`, `FourierTransform`, `Simplify`, and elementary algebra.

**The states are chosen to be physically real, not the path-of-least-resistance Gaussian.** A
Gaussian is the *unique* minimum-uncertainty state, so it is the natural answer only for the
saturation question (2.7), where it is the punchline. Everywhere else it makes the machinery
look trivial (its Fourier transform is itself, its current is featureless). So the solutions
run on states that actually occur in physics, each picked so its symbolic answer is interesting:

- the **cusp state** $\psi(x)=\sqrt{\kappa}\,e^{-\kappa|x|}$, the exact bound state of the
  attractive Dirac-delta potential and the one-dimensional analog of the hydrogen $1s$ orbital
  (Part 3); its momentum distribution is a heavy-tailed Lorentzian;
- the **sech soliton** $\psi(x)=\tfrac{1}{\sqrt2}\operatorname{sech}(x)$, the bright soliton of
  the nonlinear Schrodinger equation and the ground state of the Poschl-Teller well (Parts 8,
  21); its Fourier transform is again a sech, and it is nearly (but not exactly) minimum
  uncertainty;
- the **hydrogen $1s$ reduced radial** $u(x)=2\kappa^{3/2}x\,e^{-\kappa x}$ on the half-line
  (Part 12), which gives a clean complex recoil amplitude;
- the **Gaussian**, used once, in 2.7, as the state that saturates the uncertainty bound.

Each is used here as a given function: the computation needs only calculus, so introducing the
state before its home Part is motivation, not a dependency (the textbook habit of meeting the
Gaussian in chapter one). Natural units $\hbar=m=1$ throughout; $\kappa>0$ is an inverse length,
$k$ a wavenumber, $x_0$ a position.

**How each answer is checked.** A closed form is confirmed by substitution into the defining
relation, by agreement between two independent routes (a position-space integral against a
momentum-space integral), or by numerical convergence of a grid sum. The companion script
[`verify/part12b.wls`](verify/part12b.wls) reproduces every stated result with a `PASS`/`FAIL`.

A note on one convention trap. Wolfram's `FourierTransform` defaults to
`FourierParameters -> {0, 1}`, the kernel $\frac{1}{\sqrt{2\pi}}\int f(x)\,e^{+ipx}\,dx$. The
quantum-mechanics forward transform carries the opposite sign,
$\tilde\psi(p)=\frac{1}{\sqrt{2\pi}}\int\psi(x)\,e^{-ipx}\,dx$, the choice that makes
$\hat p=-i\,d/dx$ act as multiplication by $+p$. So the momentum transforms below use
`FourierParameters -> {0, -1}`, cross-checked against the position-space $\langle p^2\rangle$.

## Part 1. The wavefunction and the state space $L^2$

### 1.1 [BSc] How do I represent a one-dimensional wavefunction $\psi(x)$ as a Wolfram Language function and plot its probability density $|\psi(x)|^2$?

A state on the line is a complex function $\psi(x)$ whose observable content is the density
$|\psi(x)|^2$. We take the bound state of the attractive contact potential
$V(x)=-g\,\delta(x)$: a cusped exponential, the one-dimensional cousin of the hydrogen $1s$
orbital.

```wl
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]]
dens = Simplify[Abs[\[Psi][x]]^2, \[Kappa] > 0 && x \[Element] Reals]
```

```wl
Plot[dens /. \[Kappa] -> 1, {x, -5, 5}, Filling -> Axis, AxesLabel -> {"x", "|\[Psi]|^2"}]
```

The density is $|\psi(x)|^2=\kappa\,e^{-2\kappa|x|}$, with a sharp cusp at the origin where the
delta potential sits. Nothing about it is Gaussian: the exponential decay and the kink are the
signature of binding by a point interaction.

### 1.2 [BSc] How do I normalize a wavefunction with $\int_{-\infty}^{\infty}|\psi|^2\,dx=1$, and recognize a state that cannot be normalized?

Normalizing divides by the square root of the total probability; a state is admissible only if
that integral is finite ($\psi\in L^2$). We normalize the soliton profile $\operatorname{sech}(x)$.

```wl
\[Phi][x_] := Sech[x]
norm2 = Integrate[Abs[\[Phi][x]]^2, {x, -Infinity, Infinity}]
\[Psi]N = \[Phi][x]/Sqrt[norm2];
Integrate[Abs[\[Psi]N]^2, {x, -Infinity, Infinity}]
```

```wl
Integrate[Abs[Exp[I k x]]^2, {x, -Infinity, Infinity}, Assumptions -> Element[{k, x}, Reals]]
```

Since $\int\operatorname{sech}^2 x\,dx=2$, the normalized soliton is
$\psi(x)=\tfrac{1}{\sqrt2}\operatorname{sech}(x)$, which integrates to $1$. The plane wave
$e^{ikx}$ instead has $|\psi|^2=1$ everywhere and a divergent integral (`Infinity`): a momentum
eigenstate lives outside $L^2$ and must be delta-normalized (see 2.5).

### 1.3 [BSc] How do I compute the Born-rule probability $\int_a^b|\psi|^2\,dx$ of finding the particle in a region?

The probability of finding the particle in a window is the integral of the density over it. For
the soliton, the symmetric window $[-a,a]$ gives a strikingly simple answer.

```wl
\[Psi][x_] := Sech[x]/Sqrt[2];
Integrate[Abs[\[Psi][x]]^2, {x, -a, a}, Assumptions -> a > 0]
```

The probability is exactly $\tanh(a)$: it climbs from $0$ and saturates at $1$ as $a\to\infty$,
with $\tanh$ measuring how much of the soliton sits within $a$ of its center. The soliton's width
is order one, so a window of a few units already captures nearly all of it.

### 1.4 [BSc] How do I compute the inner product $\langle\phi|\psi\rangle=\int\phi^*\psi\,dx$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?

The overlap of two states is $\int\phi^*\psi\,dx$; its squared modulus is the probability that a
system prepared as $\psi$ is found in $\phi$. We take the hydrogen $1s$ radial state $u$ and the
*same* state after a sudden momentum kick $q$, $\phi=e^{iqx}u$ (the survival amplitude in the
sudden approximation).

```wl
u[x_] := 2 \[Kappa]^(3/2) x Exp[-\[Kappa] x];
\[Phi][x_] := Exp[I q x] u[x];
ovl = Integrate[Conjugate[\[Phi][x]] u[x], {x, 0, Infinity}, Assumptions -> \[Kappa] > 0 && q \[Element] Reals]
Simplify[{Abs[ovl], Abs[ovl]^2}, \[Kappa] > 0 && q \[Element] Reals]
```

The overlap is $\langle\phi|\psi\rangle=8\kappa^3/(2\kappa+iq)^3$, genuinely complex: the phase
records the momentum kick. Using $|2\kappa+iq|=\sqrt{4\kappa^2+q^2}$, its modulus is
$8\kappa^3/(4\kappa^2+q^2)^{3/2}$ and the survival probability is $64\kappa^6/(4\kappa^2+q^2)^3$,
both falling off as the kick $q$ grows past the orbital's momentum scale $\kappa$. This is the
recoil form factor probed when an atom is struck suddenly.

### 1.5 [BSc] How do I compute the position expectation $\langle x\rangle$ and the spread $\Delta x=\sqrt{\langle x^2\rangle-\langle x\rangle^2}$?

The mean position is $\int x\,|\psi|^2\,dx$ and the spread is the root-variance of that
distribution. We place the soliton at $x_0$.

```wl
\[Psi][x_] := Sech[x - x0]/Sqrt[2];
mx = Integrate[x Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals];
mx2 = Integrate[x^2 Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals];
{mx, Simplify[Sqrt[mx2 - mx^2], x0 \[Element] Reals]}
```

The soliton sits at $\langle x\rangle=x_0$ with spread $\Delta x=\pi/(2\sqrt3)$, a pure number
fixed by the $\operatorname{sech}^2$ profile (the $\pi$ comes from $\int x^2\operatorname{sech}^2 x\,dx=\pi^2/6$).
This $\Delta x$ is what enters the uncertainty product of 2.7.

### 1.6 [BSc] How do I compute the probability current density $j(x,t)$ and verify the continuity equation $\partial_t|\psi|^2+\partial_x j=0$?

Probability is locally conserved: the density's rate of change is minus the divergence of the
current $j=\tfrac{1}{2i}(\psi^*\partial_x\psi-\psi\,\partial_x\psi^*)$ (with $\hbar=m=1$). To
check continuity in full generality we treat $\psi$ and $\psi^*$ as independent functions $f,g$
and impose the free Schrodinger equation $\partial_t f=\tfrac{i}{2}\partial_x^2 f$ and its
conjugate.

```wl
\[Rho] = g[x, t] f[x, t];
j = (g[x, t] D[f[x, t], x] - f[x, t] D[g[x, t], x])/(2 I);
cont = D[\[Rho], t] + D[j, x] /.
   {Derivative[0, 1][f][x, t] -> (I/2) D[f[x, t], {x, 2}],
    Derivative[0, 1][g][x, t] -> -(I/2) D[g[x, t], {x, 2}]};
Simplify[cont]
```

The residual is exactly $0$: density and current satisfy $\partial_t|\psi|^2+\partial_x j=0$ for
*any* solution of the free Schrodinger equation. The argument used only the equation of motion,
so it holds for every state, the soliton and the cusp state included.

### 1.7 [BSc] How do I represent a state on a discrete spatial grid so that integrals become sums, building the bridge from $L^2$ to the finite vector a computer can hold?

Every numerical method replaces $\psi(x)$ by its samples and every integral by a Riemann sum
$\int|\psi|^2\,dx\approx\sum_n|\psi_n|^2\,\Delta x$. We test it on the cusp state, whose kink at
the origin is exactly the feature a grid must resolve.

```wl
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]];
grid = Subdivide[-12., 12., 6000]; dx = grid[[2]] - grid[[1]];
Total[Abs[(\[Psi][#] /. \[Kappa] -> 1.) & /@ grid]^2] dx
```

The sum returns $1.0000\ldots$, converging to the exact norm as the grid is refined. The cusp
costs a little accuracy near $x=0$ (the density is not smooth there), which is why a state with a
point-interaction kink is a more honest numerical test than a smooth Gaussian. Every later
numerical Part (eigenvalues by `NDEigensystem`, propagation by split-step) rests on this sampling.

### 1.8 [BSc] How do I show that an overall constant phase is unobservable, while the local phase gradient $\partial_x\arg\psi$ carries the current (the velocity field)?

A global phase $e^{i\alpha}$ cancels in $|\psi|^2$ and every expectation value, so it is pure
gauge. A *position-dependent* phase is physical: its gradient is the local velocity
$v=j/|\psi|^2$. We take a chirped soliton $\psi=\tfrac{1}{\sqrt2}\operatorname{sech}(x)\,e^{i\beta x^2}$,
whose phase $S(x)=\beta x^2$ varies across the packet.

```wl
\[Psi][x_] := Sech[x] Exp[I \[Beta] x^2]/Sqrt[2];
\[Psi]c[x_] := Sech[x] Exp[-I \[Beta] x^2]/Sqrt[2];
Simplify[Abs[Exp[I \[Alpha]] \[Psi][x]]^2 - Abs[\[Psi][x]]^2, Element[{\[Alpha], \[Beta], x}, Reals]]
```

```wl
j = Simplify[(\[Psi]c[x] D[\[Psi][x], x] - \[Psi][x] D[\[Psi]c[x], x])/(2 I), Element[{\[Beta], x}, Reals]];
Simplify[j/(\[Psi]c[x] \[Psi][x]), Element[{\[Beta], x}, Reals]]
```

The first cell returns $0$: the global phase changes nothing. The velocity field is
$v(x)=j/|\psi|^2=2\beta x=\partial_x(\beta x^2)$, the slope of the local phase. Unlike a plain
plane-wave phase $kx$ (uniform flow), the chirp gives a position-dependent velocity, an outflow
that grows linearly with $x$: this is exactly the phase profile that makes a free packet spread.

### 1.9 [BSc] How does a Galilean boost transform a wave packet (the boost phase factor), and why is $|\psi(x)|^2$ unchanged in shape?

Viewing a state from a frame moving at velocity $v$ multiplies it by the boost phase
$e^{i(vx-v^2t/2)}$ (with $m=1$) and translates its argument:
$\psi_v(x,t)=e^{i(vx-v^2t/2)}\,\psi(x-vt,t)$. We apply it to the soliton.

```wl
boost = Exp[I (v x - v^2 t/2)];
Simplify[Abs[boost]^2, Element[{v, x, t}, Reals]]
```

```wl
\[Psi]sol[x_] := Sech[x]/Sqrt[2];
Simplify[Abs[boost \[Psi]sol[x - v t]]^2 - Abs[\[Psi]sol[x - v t]]^2, Element[{v, x, t}, Reals]]
```

The boost factor has unit modulus, so $|\psi_v(x,t)|^2=|\psi(x-vt,t)|^2$: the soliton is carried
along rigidly at speed $v$, its $\operatorname{sech}^2$ shape untouched. The boost lives entirely
in the phase, invisible to a position measurement yet adding $v$ to every momentum component, the
hallmark of a soliton that translates without distortion.

## Part 2. Operators, position and momentum, and the Fourier bridge

### 2.1 [BSc] How do I represent the position operator $\hat x$ as multiplication and the momentum operator $\hat p=-i\,d/dx$ as differentiation, acting on a test function?

In the position representation $\hat x$ multiplies by $x$ and $\hat p$ differentiates. We define
them as operations and apply them to a plane wave.

```wl
xop[f_] := x f;
pop[f_] := -I D[f, x];
{xop[Exp[I k x]], pop[Exp[I k x]]}
```

$\hat x$ returns $x\,e^{ikx}$ and $\hat p$ returns $k\,e^{ikx}$: the plane wave is an eigenfunction
of momentum with eigenvalue $k$. That is the precise sense in which $e^{ikx}$ carries momentum $k$.

### 2.2 [BSc] How do I verify the canonical commutator $[\hat x,\hat p]=i$ by letting it act on an arbitrary function?

The defining relation $[\hat x,\hat p]=i$ is an operator identity, so we check it on a generic
test function $f(x)$.

```wl
Simplify[xop[pop[f[x]]] - pop[xop[f[x]]]]
```

The result is $i\,f(x)$ for arbitrary $f$, i.e. $[\hat x,\hat p]=i$. The surviving term comes from
differentiating the product $x f$: position and momentum fail to commute by exactly one unit.

### 2.3 [BSc] How do I compute the momentum expectation $\langle p\rangle$ and spread $\Delta p$ in the position representation?

Without leaving position space, $\langle p\rangle=\int\psi^*(-i\partial_x)\psi\,dx$. For
$\langle p^2\rangle$ we use the manifestly positive, boundary-term-free form
$\langle p^2\rangle=\int|\hat p\psi|^2\,dx=\int|\partial_x\psi|^2\,dx$ (equal because $\hat p$ is
Hermitian), which sidesteps the delta that the cusp would put into $\partial_x^2\psi$. We take the
moving cusp state $\psi=\sqrt{\kappa}\,e^{-\kappa|x|}e^{ikx}$.

```wl
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]] Exp[I k x];
mp = Integrate[Conjugate[\[Psi][x]] (-I D[\[Psi][x], x]), {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && k \[Element] Reals];
p2 = Integrate[Abs[D[\[Psi][x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && k \[Element] Reals];
{Simplify[mp, \[Kappa] > 0 && k \[Element] Reals], Simplify[Sqrt[p2 - mp^2], \[Kappa] > 0 && k \[Element] Reals]}
```

The mean momentum is $\langle p\rangle=k$ (the drift stamped into the phase), while
$\langle p^2\rangle=k^2+\kappa^2$, so the intrinsic spread is $\Delta p=\kappa$. The cusp sets the
momentum width: a more sharply bound state (larger $\kappa$) is broader in momentum, the
uncertainty relation at work.

### 2.4 [BSc] How do I show that $\hat p$ is Hermitian by integration by parts, and see that the boundary terms vanish for an $L^2$ state?

Hermiticity of $\hat p$ is $\langle\phi|\hat p\chi\rangle=\langle\hat p\phi|\chi\rangle$, which
integration by parts turns into the boundary term $-i\,[\phi^*\chi]_{-\infty}^{\infty}$: the
matrix-element difference *equals* that surface term. We take two smooth solitons; their momentum
matrix elements have no elementary closed form, so we read off the surface term symbolically and
confirm the equality by quadrature.

```wl
\[Phi][x_] := Sech[x]/Sqrt[2];
\[Chi][x_] := Sech[x - 1] Exp[I x]/Sqrt[2];
Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> Infinity] -
  Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> -Infinity]
```

```wl
Chop[NIntegrate[Conjugate[\[Phi][x]] (-I D[\[Chi][x], x]), {x, -Infinity, Infinity}, WorkingPrecision -> 30] -
     NIntegrate[Conjugate[-I D[\[Phi][x], x]] \[Chi][x], {x, -Infinity, Infinity}, WorkingPrecision -> 30]]
```

The boundary term is $0$ (the solitons decay at infinity), and the high-precision quadrature
returns numerical zero, so the two matrix elements agree. Because the difference *is* that surface
term, Hermiticity of $\hat p$ is exactly its vanishing, a point that turns subtle on a finite
domain (2.9).

### 2.5 [BSc] How do I write the momentum eigenfunctions $e^{ipx}$, see that they are not square-integrable, and normalize them to a Dirac delta?

The momentum eigenfunctions $e^{ipx}$ are not in $L^2$ (1.2), so they cannot be normalized to
$1$; instead $\langle p'|p\rangle=\delta(p-p')$ for $u_p(x)=\frac{1}{\sqrt{2\pi}}e^{ipx}$. The
$2\pi$ that makes this work is the Fourier integral of a constant.

```wl
FourierTransform[1, x, q]
```

```wl
FourierTransform[1/Sqrt[2 Pi] Exp[I p x], x, pp, FourierParameters -> {0, -1}]
```

The first cell gives $\sqrt{2\pi}\,\delta(q)$, i.e. $\int e^{iqx}\,dx=2\pi\,\delta(q)$. The second
transforms the normalized eigenstate $u_p$ into momentum space and returns exactly
$\delta(p-p')$: a sharp momentum state is a delta spike at its eigenvalue, the continuum analog of
an orthonormal basis vector.

### 2.6 [BSc] How do I obtain the momentum-space wavefunction $\tilde\psi(p)$ by Fourier transform and compute $\langle p^n\rangle$ as a moment in $p$?

The momentum-space wavefunction is
$\tilde\psi(p)=\frac{1}{\sqrt{2\pi}}\int\psi(x)e^{-ipx}\,dx$ (`FourierParameters -> {0, -1}`),
after which $\langle p^n\rangle=\int p^n|\tilde\psi(p)|^2\,dp$ is an ordinary moment. We transform
the soliton.

```wl
\[Psi][x_] := Sech[x]/Sqrt[2];
\[Psi]T = FullSimplify[FourierTransform[\[Psi][x], x, p, FourierParameters -> {0, -1}]]
```

```wl
{Integrate[Abs[\[Psi]T]^2, {p, -Infinity, Infinity}],
 Integrate[p^2 Abs[\[Psi]T]^2, {p, -Infinity, Infinity}]}
```

The transform is $\tilde\psi(p)=\tfrac{\sqrt\pi}{2}\operatorname{sech}\!\big(\tfrac{\pi p}{2}\big)$:
the Fourier transform of a sech is again a sech, a self-similarity special to this profile. Its
norm is $1$ and $\langle p^2\rangle=\tfrac13$, matching the position-space value
$\int|\partial_x\psi|^2\,dx$. (A boost $e^{ikx}$ would shift the whole momentum profile to center
$k$ without changing its width.)

### 2.7 [BSc] How do I verify the position-momentum uncertainty relation $\Delta x\,\Delta p\ge 1/2$ and show that a Gaussian saturates it?

Here the Gaussian earns its keep. The bound $\Delta x\,\Delta p\ge\tfrac12$ (with $\hbar=1$) is
saturated by a Gaussian and by nothing else. We compute the product for the soliton, then for the
Gaussian, side by side.

```wl
\[Psi]sol[x_] := Sech[x]/Sqrt[2];
dxSol = Sqrt[Integrate[x^2 Abs[\[Psi]sol[x]]^2, {x, -Infinity, Infinity}]];
dpSol = Sqrt[Integrate[Abs[D[\[Psi]sol[x], x]]^2, {x, -Infinity, Infinity}]];
dxSol dpSol
```

```wl
\[Psi]G[x_] := Pi^(-1/4) \[Sigma]^(-1/2) Exp[-x^2/(2 \[Sigma]^2)];
dxG = Sqrt[Integrate[x^2 Abs[\[Psi]G[x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Sigma] > 0]];
dpG = Sqrt[Integrate[Abs[D[\[Psi]G[x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Sigma] > 0]];
Simplify[dxG dpG, \[Sigma] > 0]
```

The soliton gives $\Delta x\,\Delta p=\pi/6\approx0.524$, just above the floor: a state can be
beautifully localized and still not minimal. The Gaussian gives exactly $\tfrac12$ for every
width $\sigma$, sitting right on the boundary. Squeezing it in position inflates it in momentum
and vice versa, tracing the edge of the allowed region; every non-Gaussian state, the soliton
included, lies strictly above.

### 2.8 [MSc] How do I build a function $f(\hat x,\hat p)$ of noncommuting operators, and resolve the ordering ambiguity by Weyl symmetric ordering?

Because $\hat x$ and $\hat p$ do not commute, the classical product $xp$ has several inequivalent
quantizations. The Weyl (symmetric) prescription, $\tfrac12(\hat x\hat p+\hat p\hat x)$, is the
unique Hermitian choice.

```wl
sym[h_] := (xop[pop[h]] + pop[xop[h]])/2;
Simplify[sym[h[x]]]
```

```wl
u[x_] := Pi^(-1/4) Exp[-(x - 1)^2/2];
w[x_] := Pi^(-1/4) Exp[-(x + 1)^2/2] Exp[I x];
{Integrate[Conjugate[u[x]] sym[w[x]], {x, -Infinity, Infinity}] -
   Integrate[Conjugate[sym[u[x]]] w[x], {x, -Infinity, Infinity}],
 Integrate[Conjugate[u[x]] xop[pop[w[x]]], {x, -Infinity, Infinity}] -
   Integrate[Conjugate[xop[pop[u[x]]]] w[x], {x, -Infinity, Infinity}]}
```

The symmetric operator acts as $-i\big(\tfrac12 f+x f'\big)$, and its matrix-element difference is
$0$ while the naive ordering $\hat x\hat p$ gives a nonzero difference. Only
$\tfrac12(\hat x\hat p+\hat p\hat x)$ is Hermitian, so only it can represent an observable; the
ambiguity $\hat x\hat p$ versus $\hat p\hat x$ is the commutator $[\hat x,\hat p]=i$ of 2.2. (The
two probe states here are displaced Gaussians, used only as convenient smooth test functions, not
as the object of study.)

### 2.9 [MSc] How do I distinguish self-adjointness from formal Hermiticity through the operator domain and boundary conditions (the box versus the half-line)?

Formal Hermiticity is the bulk identity of 2.4; self-adjointness additionally needs the boundary
term $-i[\phi^*\psi]$ to vanish on a domain the operator preserves. On a box this is achievable;
on the half-line (the natural home of the radial state $u$ of 1.4) it is not.

```wl
(* box [0,L] with twisted boundary condition \[Psi](L) = e^{i\[Theta]} \[Psi](0) *)
bt = -I (Conjugate[pL] sL - Conjugate[p0] s0) /. {sL -> Exp[I \[Theta]] s0, pL -> Exp[I \[Theta]] p0};
Simplify[bt, \[Theta] \[Element] Reals]
```

```wl
(* half-line [0,\[Infinity]): \[Psi](\[Infinity]) = 0 but \[Psi](0) is free *)
-I (0 - Conjugate[p0] s0)
```

On the box the twisted condition $\psi(L)=e^{i\theta}\psi(0)$ makes the surface term vanish for
every real $\theta$, so each $\theta$ defines a genuine self-adjoint momentum operator (a
one-parameter family of extensions). On the half-line the term reduces to $i\,\phi^*(0)\psi(0)$,
which cannot vanish for all admissible states: $\hat p=-i\,d/dx$ on $[0,\infty)$ is Hermitian but
has *no* self-adjoint extension. Hermiticity is a property of a formula; self-adjointness is a
property of a formula together with its domain.
