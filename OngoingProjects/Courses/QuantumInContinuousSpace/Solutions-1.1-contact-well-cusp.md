---
Template: Default
---

# Archived answer 1.1: the contact-well cusp state

This is the earlier worked answer to question 1.1 of `Question-List.md`, built around the cusped
exponential $\psi(x)=\sqrt{\kappa}\,e^{-\kappa|x|}$, the bound state of the attractive contact potential
$V(x)=-g\,\delta(x)$. It was archived on 2026-07-22 when 1.1 in `Solutions-Parts-1-3.md` was rewritten
around a different state (the Morse oscillator ground state), so that the contact-well derivation is
preserved without being lost from the record. The delta well itself still lives in the solution set as
answer 3.9 (found there from the jump condition). Pure Wolfram Language, natural units $\hbar=m=1$.

### 1.1 [BSc] How do I represent a one-dimensional wavefunction $\psi(x)$ as a Wolfram Language function and plot its probability density $|\psi(x)|^2$?

A state on the line is a complex function $\psi(x)$ in $L^2(\mathbb{R})$. Measuring position returns a
single number, one point of the continuous spectrum of $\hat x$, and Born's rule fixes only the
statistics of that outcome: $|\psi(x)|^2\,dx$ is the probability that the result falls in $dx$, a
distribution an ensemble of identically prepared states reproduces and a single run never does. The
density is therefore not the state's whole physical content, only the outcome distribution of one
observable: the phase it discards carries the momentum distribution, the current, and the energy, so
two states with the same density can differ in every other measurement.

We take the cusped exponential $\psi(x)=\sqrt{\kappa}\,e^{-\kappa|x|}$, the one-dimensional cousin of
the hydrogen $1s$ orbital: a normalizable amplitude with a sharp kink at the origin and exponential
tails. We will show it is the bound state of the attractive contact potential $V(x)=-g\,\delta(x)$ with
$g>0$: a normalizable solution of $-\tfrac12\psi''-g\,\delta(x)\psi=E\psi$ lying below the continuum edge
$E=0$. (Away from the origin the potential vanishes and the particle is free, so every $E\ge0$ is a
non-normalizable scattering state; a bound state is normalizable with $E<0$, too little energy to
escape to infinity.) The argument has two moves: away from the origin the delta sets the shape, and at
the origin it fixes the one remaining number.

Away from the origin the delta is silent and the equation is free, $-\tfrac12\psi''=E\psi$, whose
normalizable ($E<0$) solutions decay as $e^{-\kappa|x|}$ with $E=-\kappa^2/2$. Because $V$ is even the
bound state is even, so the two tails carry equal weight and the state is exactly
$\sqrt{\kappa}\,e^{-\kappa|x|}$, with only $\kappa$ left to find. Confirm the decaying branch solves the
free equation:

```wl
\[Psi]R[x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] x];
FullSimplify[-(1/2) \[Psi]R''[x] - (-\[Kappa]^2/2) \[Psi]R[x], \[Kappa] > 0]
```

The residual is $0$, so a state that decays at rate $\kappa$ has energy $E=-\kappa^2/2$. The exact
eigenstate therefore already lies in this one-parameter family $\sqrt{\kappa}\,e^{-\kappa|x|}$, with only
$\kappa$ left to pin down. Every member is normalized, so normalization is not what pins it:

```wl
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]];
Integrate[Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0]
```

The norm is $1$, confirming the family is admissible, so the decay rate is decided instead at the
origin, the one place the potential acts. A naive substitution cannot do the deciding, and it is worth
seeing why. Ask the kernel for the second derivative of $|x|$, the cusp's contribution to $\psi''$:

```wl
D[Abs[x], {x, 2}]
```

It comes back as a formal second derivative with no `DiracDelta`: the delta the kink should produce is
silently dropped, so substituting $\psi$ into $-\tfrac12\psi''-g\,\delta(x)\psi=E\psi$ would lose the
very term the potential is there to cancel. The cure is to use the delta the way it is defined, inside
an integral, $\int\delta(x)f(x)\,dx=f(0)$: pair the eigenvalue equation against an arbitrary smooth
probe $\varphi(x)=e^{-b x^2}$ and integrate over the line. One integration by parts moves a derivative
onto the probe, turning the kinetic term into $\tfrac12\int\psi'\varphi'\,dx$ with
$\psi'(x)=-\kappa\,\operatorname{sign}(x)\,\psi(x)$ the slope of the cusped exponential, while the delta
contributes exactly $-g\int\delta(x)\,\psi\varphi\,dx=-g\,\psi(0)\varphi(0)$. With the off-origin energy
$E=-\kappa^2/2$, the whole pairing must vanish for $\psi$ to be an eigenstate:

```wl
\[CurlyPhi][x_] := Exp[-b x^2];
weakResidual = FullSimplify[
   (1/2) Integrate[(-\[Kappa] Sign[x] \[Psi][x]) \[CurlyPhi]'[x], {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && b > 0]
    - g Integrate[DiracDelta[x] \[Psi][x] \[CurlyPhi][x], {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && b > 0]
    + (\[Kappa]^2/2) Integrate[\[Psi][x] \[CurlyPhi][x], {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && b > 0],
   \[Kappa] > 0 && b > 0]
```

The residual is $\sqrt{\kappa}\,(\kappa-g)$, with no trace of the probe width $b$: it vanishes for every
probe in this family exactly when $\kappa=g$. That is the point interaction doing its one job, fixing
the decay rate $\kappa=g$ and with it the energy $E=-\kappa^2/2=-g^2/2$, so $\psi(x)=\sqrt{g}\,e^{-g|x|}$
is the exact bound state of the contact well.

The same $\kappa=g$ falls out of energetics, and this second route says why the balance sits where it
does. Build $\langle H\rangle(\kappa)=\langle T\rangle+\langle V\rangle$ over the family and minimize
it. The potential acts at a single point, and a delta is defined by its action inside an integral,
$\int\delta(x)f(x)\,dx=f(0)$, so the potential energy samples the density at the origin,
$\langle V\rangle=-g\int\delta(x)|\psi(x)|^2\,dx=-g\,|\psi(0)|^2$. Compute it:

```wl
Vexp = -g Integrate[DiracDelta[x] Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0]
```

It is $\langle V\rangle=-g\kappa$: a deeper well (larger $g$) or a more concentrated state (larger
$\kappa$) both lower the potential energy, so the well on its own would drive $\kappa$ up without
limit. The kinetic term is what stops it, and here the cusp matters. Because $\psi$ has a kink at the
origin, $\psi''$ carries a delta there, and the form $-\tfrac12\int\psi\,\psi''\,dx$ would need that
singular piece handled separately. The equivalent gradient form
$\langle T\rangle=\tfrac12\int(\psi')^2\,dx$ has no such trouble: $\psi'$ is merely a step, so
$(\psi')^2$ is an ordinary integrable function. Its integrand is even, so the integral over the whole
line is twice the integral over the right half, and that factor $2$ cancels the $\tfrac12$ in front,
leaving $\langle T\rangle=\int_0^{\infty}(\psi_R')^2\,dx$ over the smooth branch. Compute it:

```wl
T = Integrate[D[\[Psi]R[x], x]^2, {x, 0, Infinity}, Assumptions -> \[Kappa] > 0]
```

It is $\langle T\rangle=\kappa^2/2$: squeezing the state costs kinetic energy, the price of a steeper
wavefunction. Setting the two against each other, $\langle H\rangle(\kappa)=\tfrac{\kappa^2}{2}-g\kappa$,
the well pulling $\kappa$ up and the kinetic cost pushing it down, with the bound state where they
balance. The minimum is where $d\langle H\rangle/d\kappa=0$:

```wl
Solve[D[T + Vexp, \[Kappa]] == 0, \[Kappa]]
```

So $\kappa=g$ again, a genuine minimum since $d^2\langle H\rangle/d\kappa^2=1>0$: the decay rate is set
by the coupling. Evaluate the energy there:

```wl
(T + Vexp) /. \[Kappa] -> g
```

The energy is $E=-g^2/2$, a binding energy of $g^2/2$, matching the $-\kappa^2/2$ branch at $\kappa=g$.
Because the weak-form pairing already proved $\sqrt{g}\,e^{-g|x|}$ solves the equation, this variational
minimizer is not merely a good trial state but the exact eigenstate. A weaker well ($g\to0$) unbinds it
as its energy rises to $0$; a deeper well binds it harder and squeezes it in. Now read off the
density:

```wl
dens = Simplify[Abs[\[Psi][x]]^2, \[Kappa] > 0 && x \[Element] Reals]
```

The density is $|\psi(x)|^2=\kappa\,e^{-2\kappa|x|}$. Now plot it:

```wl
Plot[dens /. \[Kappa] -> 1, {x, -5, 5}, Filling -> Axis, 
 PlotRange -> All, AxesLabel -> {"x", "|\[Psi]|^2"}]
```

It has a sharp cusp at the origin where the delta potential sits. Nothing about it is Gaussian: the
exponential decay and the kink are the signature of binding by a point interaction, and both fell
out of setting the delta's potential energy against the kinetic cost.
