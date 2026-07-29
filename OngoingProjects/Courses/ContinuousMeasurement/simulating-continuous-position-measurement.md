# Watching a Quantum Particle

**Simulating the stochastic Schrödinger equation of continuous position measurement with the Wolfram Language: where the equation comes from, what it tells you before you ever run it, how to integrate it correctly, and how to know that you did.**

Mads Bahrami (last updated: July 22, 2026)

### Setting the Stage: How This Document Flows

This is a computational essay on the continuous quantum measurement: how a weak Gaussian detector turns into a stochastic differential equation, why that equation is nonlinear but not violating no-signaling, what its conditional variances do before you simulate anything, and how to properly compute dynamical equations and relevant observables in Wolfram Mathematica.

This essay is meant to be understand computationally. It means the code should do the heavy lifting for grasping some complicated concepts. Although some codes may look complicated, but focus on the output and its meaning before worrying about every detail of the input. And remember that you are not locked into the code as given. You can (and should) modify it, and run your own numerical experiments.

I assume the reader is familiar with quantum theory at nonrelativistic Schrodinger equation, and enough Wolfram Language to read codes. No stochastic calculus is assumed: the two Itô facts we need are introduced where they are used. The pace steepens where we hand the kernel the canonical commutation relation and let it do the operator algebra. Throughout, $\hbar$, $m$ and the measurement strength $\lambda$ stay symbolic for as long as possible.

Let's start.

## Part I: From a Detector to a Differential Equation

### The Itô Equation: The Continuum Limit of Weak Measurement

Let $t_k=k\Delta t$ and write $|\psi_k\rangle=|\psi(t_k)\rangle$. During the interval $[t_k,t_{k+1})$, which we will call bin $k$, the detector directly returns one real number $\bar x_k$. The sequence $\{\bar x_k\}$ is the measurement record at finite time resolution.

The measurement during one interval is defined by the Gaussian measurement operators

$$
\hat M(\bar x_k)
=\left(\frac{2\lambda\Delta t}{\pi}\right)^{1/4}
e^{-\lambda\Delta t(\hat x-\bar x_k)^2}.
$$

Here $\lambda>0$ is the measurement strength. The prefactor ensures

$$
\int_{-\infty}^{\infty}
\hat M^\dagger(\bar x_k)\hat M(\bar x_k)\,d\bar x_k
=\hat I.
$$

Consequently, the same measurement operator determines both the probability density of the detector output,

$$
p(\bar x_k\mid\psi_k)
=\|\hat M(\bar x_k)|\psi_k\rangle\|^2,
$$

and the state conditioned on the measured value,

$$
|\psi_{k+1}\rangle
=\frac{\hat U\,\hat M(\bar x_k)\,|\psi_k\rangle}
{\bigl\|\hat U\,\hat M(\bar x_k)\,|\psi_k\rangle\bigr\|},
\qquad
\hat U=e^{-i\hat H\Delta t/\hbar},
\qquad
\hat H=\hat H^\dagger.
$$

To make the output density explicit, write $\psi_k(x)=\langle x|\psi_k\rangle$, where $\hat x|x\rangle=x|x\rangle$. Then

$$
p(\bar x_k\mid\psi_k)
=\sqrt{\frac{2\lambda\Delta t}{\pi}}
\int_{-\infty}^{\infty}
|\psi_k(x)|^2e^{-2\lambda\Delta t(\bar x_k-x)^2}\,dx.
$$

The Gaussian kernel is normalized and centered at $\bar x_k=x$, with variance $1/(4\lambda\Delta t)$. The detector-output density is therefore the Born position density blurred by the detector response, which is also a recipe for sampling one reading: draw $x$ from $|\psi_k(x)|^2$, then add independent detector noise of variance $1/(4\lambda\Delta t)$. Because the kernel is centered, the conditional mean and variance follow at once:

$$
\mathbb E[\bar x_k\mid\psi_k]=\langle\hat x\rangle_k,
\qquad
\operatorname{Var}(\bar x_k\mid\psi_k)
=\bigl\langle(\hat x-\langle\hat x\rangle_k)^2\bigr\rangle_k
+\frac{1}{4\lambda\Delta t},
$$

where $\langle\hat x\rangle_k\equiv\langle\psi_k|\hat x|\psi_k\rangle$ is evaluated before the reading in bin $k$.

The raw readings have no continuous-time limit, since that variance diverges as $1/\Delta t$. What does have one is the accumulated record, and the way it accumulates fixes both the power of $\Delta t$ and the constant in front. Scale the readings as $c\,\bar x_k\Delta t^{\alpha}$ and add up the $N=T/\Delta t$ bins in $[0,T]$; the detector noise is independent from bin to bin, so the variances add, giving $(c^2T/4\lambda)\,\Delta t^{2\alpha-2}$ to leading order. This is finite and nonzero only at $\alpha=1$; at $\alpha=1/2$ the accumulated noise diverges, and at $\alpha=3/2$ it vanishes, leaving a deterministic record. With $\alpha=1$ fixed, the choice $c=2\sqrt{\lambda}$ makes the accumulated variance exactly $T$, which is the statement $(dW_t)^2=dt$. Neither choice is free. Define

$$
Y_0=0,
\qquad
\Delta Y_k\equiv2\sqrt{\lambda}\,\bar x_k\Delta t,
\qquad
Y_{t_{k+1}}=Y_{t_k}+\Delta Y_k.
$$

The detector still measures $\bar x_k$ and $Y$ is computed from those readings; it is a repackaging, not an additional measured quantity, since $\bar x_k=(Y_{t_{k+1}}-Y_{t_k})/(2\sqrt{\lambda}\,\Delta t)$ recovers them exactly.

The mean above gives $\mathbb E[\Delta Y_k\mid\psi_k]=2\sqrt{\lambda}\,\langle\hat x\rangle_k\Delta t$. Remove it from the record increment:

$$
\Delta W_k
\equiv
\Delta Y_k-2\sqrt{\lambda}\,\langle\hat x\rangle_k\Delta t
=2\sqrt{\lambda}\,\Delta t
\bigl(\bar x_k-\langle\hat x\rangle_k\bigr).
$$

Nothing new is being assumed about $\Delta W_k$: this relation is exact and invertible, so its law is inherited from the measured-output law rather than imposed alongside it. Because the relation is affine, it carries the two-step sampling recipe through term by term:

$$
\Delta W_k
=\underbrace{2\sqrt{\lambda}\,\Delta t
\bigl(x-\langle\hat x\rangle_k\bigr)}_{\text{Born draw}}
+\underbrace{\xi_k}_{\text{detector noise}},
\qquad
x\sim|\psi_k(x)|^2,
\qquad
\xi_k\sim\mathcal N(0,\Delta t).
$$

The rescaled detector noise is exactly $\mathcal N(0,\Delta t)$ at every $\Delta t$, which is what the constant $2\sqrt{\lambda}$ was chosen to arrange. The Born draw has mean zero, so $\mathbb E[\Delta W_k\mid\psi_k]=0$ exactly, and it is the only state-dependent piece. Its variance relative to the detector's is

$$
\frac{
4\lambda\Delta t^2
\bigl\langle(\hat x-\langle\hat x\rangle_k)^2\bigr\rangle_k
}{\Delta t}
=
4\lambda\Delta t
\bigl\langle(\hat x-\langle\hat x\rangle_k)^2\bigr\rangle_k
\longrightarrow0.
$$

So provided the position variance stays finite, the state washes out of the noise and $\Delta W_k$ becomes Gaussian with mean zero and variance $\Delta t$. Since that limiting law no longer refers to $|\psi_k\rangle$, it does not refer to the earlier readings that determined $|\psi_k\rangle$ either; repeating the same detector rule in every bin therefore gives independent increments, whose cumulative sum converges to a Wiener process $W_t$.

We now derive the finite-bin state update. Every coefficient is evaluated from $|\psi_k\rangle$ before $\bar x_k$ is measured; this beginning-of-bin convention will give the Itô equation in the continuum limit. From the definition of $\Delta W_k$,

$$
\bar x_k
=\langle\hat x\rangle_k
+\frac{\Delta W_k}{2\sqrt{\lambda}\,\Delta t}.
$$

Define

$$
\hat A_k\equiv\hat x-\langle\hat x\rangle_k.
$$

Substitution into the exponent of $\hat M(\bar x_k)$ gives

$$
\begin{aligned}
-\lambda\Delta t(\hat x-\bar x_k)^2
&=-\lambda\Delta t
\left(
\hat A_k-\frac{\Delta W_k}
{2\sqrt{\lambda}\,\Delta t}
\right)^2\\
&=\sqrt{\lambda}\,\hat A_k\Delta W_k
-\lambda\hat A_k^2\Delta t
-\frac{(\Delta W_k)^2}{4\Delta t}.
\end{aligned}
$$

The last term is a scalar: it multiplies the entire state by the same positive number and cancels exactly in the normalized update. Only

$$
e^{\sqrt{\lambda}\hat A_k\Delta W_k-\lambda\hat A_k^2\Delta t}
$$

can change the normalized state. Since $\Delta W_k=O(\sqrt{\Delta t})$, its finite-bin expansion is

$$
e^{\sqrt{\lambda}\hat A_k\Delta W_k-\lambda\hat A_k^2\Delta t}
=\hat I+\sqrt{\lambda}\hat A_k\Delta W_k
-\lambda\hat A_k^2\Delta t
+\frac{\lambda}{2}\hat A_k^2(\Delta W_k)^2
+O(\Delta t^{3/2}).
$$

The finite-bin unitary factor is

$$
\hat U
=\hat I-\frac{i}{\hbar}\hat H\Delta t+O(\Delta t^2).
$$

All finite-bin quantities have now been expanded. Take the continuum limit once:

$$
\Delta t\longrightarrow dt,
\qquad
\Delta Y_k\longrightarrow dY_t,
\qquad
\Delta W_k\longrightarrow dW_t,
\qquad
(dW_t)^2=dt.
$$

The measurement factor becomes

$$
\hat I
+\sqrt{\lambda}\hat A_t\,dW_t
-\frac{\lambda}{2}\hat A_t^2dt,
\qquad
\hat A_t\equiv\hat x-\langle\hat x\rangle_t.
$$

Its squared norm differs from one by $2\sqrt{\lambda}\langle\hat A_t\rangle_t\,dW_t+\lambda\langle\hat A_t^2\rangle_t[(dW_t)^2-dt]=0$, so normalization adds no term. The record equation is therefore

$$
dY_t
=2\sqrt{\lambda}\,\langle\hat x\rangle_t\,dt+dW_t.
$$

Multiplying the measurement and unitary factors gives the state equation

$$
d|\psi_t\rangle
=\left[
-\frac{i}{\hbar}\hat H\,dt
-\frac{\lambda}{2}
\bigl(\hat x-\langle\hat x\rangle_t\bigr)^2dt
+\sqrt{\lambda}
\bigl(\hat x-\langle\hat x\rangle_t\bigr)dW_t
\right]|\psi_t\rangle.
$$

The same $dW_t$ appears here and in the record equation because it is the centered random part of the measured record.

### The Norm Cancellation Along a Conditioned Trajectory

Let us verify the norm claim rather than assert it. Suppress the time subscript and write $\hat A=\hat x-\langle\hat x\rangle$. For an Itô equation $d|\psi\rangle=(\hat D\,dt+\hat N\,dW)|\psi\rangle$, the product rule gives

$$d\langle\psi|\psi\rangle
=\langle\hat N+\hat N^\dagger\rangle\,dW
+\langle\hat D+\hat D^\dagger+\hat N^\dagger\hat N\rangle\,dt.$$

For the equation above, $\hat D=-i\hat H/\hbar-\lambda\hat A^2/2$ and $\hat N=\sqrt{\lambda}\hat A$. The noise coefficient vanishes because $\langle\hat A\rangle=0$. The remaining drift contains $\hat D+\hat D^\dagger$ and the Itô correction $\hat N^\dagger\hat N$. The constants $\hbar$ and $\lambda$ do not affect this cancellation, so set them to $1$ and define a function that tests the drift for a random Hamiltonian, a random measured observable and a random state in $n$ dimensions:

```wl
ClearAll[normDrift];
normDrift[n_] := With[
   {h = # + ConjugateTranspose[#] &@RandomComplex[{-1 - I, 1 + I}, {n, n}],
    xop = # + ConjugateTranspose[#] &@RandomComplex[{-1 - I, 1 + I}, {n, n}],
    ψ = Normalize[RandomComplex[{-1 - I, 1 + I}, n]]},
   {a = xop - (Conjugate[ψ] . xop . ψ) IdentityMatrix[n]},
   {dOp = -I h - a . a/2},
   Conjugate[ψ] . (dOp + ConjugateTranspose[dOp] + ConjugateTranspose[a] . a) . ψ];
```

The multi-argument `With` (successive binding lists `With[{h,...}, {a}, {dOp}, body]`, each list seeing the ones before it) keeps every intermediate a named constant rather than a reassigned variable, without a tower of nested brackets. Now run it across dimensions and read the result:

```wl
Chop[Table[normDrift[n], {n, 2, 8}]]
```

Here `Chop` sets tiny machine-precision residue to zero so the output reads cleanly. As one can see, the norm drift vanishes in every tested dimension for independently generated Hermitian operators and states. Together with $\langle\hat A\rangle=0$, this confirms norm preservation through Itô order. Note that $\hat H$ drops out entirely, so the cancellation is a property of the measurement terms.

This cancellation keeps every fully conditioned trajectory normalized and pure. It does not imply that the measurement causes no decoherence. If the detector record is ignored and the projectors $|\psi_t\rangle\langle\psi_t|$ are averaged over all possible records, the resulting density operator contains the measurement dissipator $-\frac{\lambda}{2}[\hat x,[\hat x,\rho]]$.

## Part II: What the Equation Says Before You Simulate It

### The Canonical Commutator: One Relation, Handed to the Kernel

Here is a habit worth acquiring: before writing a simulation, extract every exact statement you can, because those become the benchmarks that catch the simulation being wrong. For this equation a great deal is exactly computable, and all of it follows from one commutator, $[\hat x,\hat p] = i\hbar$.

Version 15.0 added a noncommutative-algebra system to the Wolfram Language, so we no longer hand-roll an operator ordering. We declare the two generators $\hat x$ and $\hat p$ and the single relation between them, and let `NonCommutativeAlgebra` carry the algebra. The commutation relation is entered leading-monomial-positive, $\hat p\hat x - \hat x\hat p + i\hbar = 0$, which says $\hat p\hat x = \hat x\hat p - i\hbar$. Define the algebra and two thin wrappers for expanding and for commutators:

```wl
ClearAll[xp, nc, comm];
xp = NonCommutativeAlgebra[<|
    "Generators" -> {{xo, po}},
    "CommutationRelations" -> {po ** xo - xo ** po + I ℏ}|>];
nc[e_] := NonCommutativeExpand[e, xp];
comm[a_, b_] := nc[Commutator[a, b, xp]];
xp["Unity"]
```

The algebra's unity is a plain `1`, not a stray formal object, which matters more than it looks: a hand-rolled `NonCommutativeMultiply` has no identity reduction and leaks empty products into everything downstream, and here that whole class of bug simply does not exist. Now check the one relation, the double commutator we will need for the heating law, and a normal-ordering:

```wl
{comm[xo, po], comm[xo, comm[xo, po ** po]], nc[xo ** po ** xo]}
```

As one can see, the canonical commutator is $i\hbar$, the double commutator $[\hat x,[\hat x,\hat p^2]] = -2\hbar^2$ is a c-number, and `NonCommutativeExpand` puts every product into canonical form with all $\hat x$ to the left of all $\hat p$, so $\hat x\hat p\hat x$ normal-orders to $\hat x^2\hat p - i\hbar\hat x$. That canonical form is exactly the normal ordering we need for the Gaussian closure below.

### The Heating Law: One Line, and It Is Exact

Averaging the stochastic equation over the noise gives back a Lindblad equation, $\dot\rho = -\frac{i}{\hbar}[\hat H,\rho] - \frac{\lambda}{2}[\hat x,[\hat x,\rho]]$. That means the ensemble-averaged rate of change of any observable $\hat O$ picks up a term $-\frac{\lambda}{2}\langle[\hat x,[\hat x,\hat O]]\rangle$. For $\hat O = \hat p^2$ we already have the double commutator in hand:

```wl
nc[-(λ/2) Commutator[xo, Commutator[xo, po ** po, xp], xp]]
```

Therefore $\frac{d}{dt}\mathbb{E}\langle \hat p^2\rangle = \lambda\hbar^2$, exactly, with no approximation. Look at what this does *not* depend on: not the state, not the mass, not the Hamiltonian (as long as $\hat H$ conserves $\langle p^2\rangle$ on its own, as a free particle does). Watching a particle's position heats it at a fixed rate forever, which is the price the uncertainty principle charges for position information. This will be our sharpest numerical benchmark, precisely because it is a straight line with a known slope and nothing to tune.

### Gaussian Closure: Wick's Theorem for Noncommuting Moments

To get the conditional variances we need expectation values of normal-ordered monomials in a Gaussian state. For Gaussian states all moments follow from the first two by Wick's theorem, pairing factors two at a time. The one subtlety is that the pairing is *order-dependent* here, because $\hat x$ and $\hat p$ do not commute: $\langle \delta x\,\delta p\rangle = C + i\hbar/2$ while $\langle \delta p\,\delta x\rangle = C - i\hbar/2$, where $C$ is the symmetrized covariance.

Set the two-point data (the covariances, ordered) and the pairing recursion that builds every higher even moment from them:

```wl
ClearAll[cov, wick];
cov["x", "x"] = Vx; cov["p", "p"] = Vp;
cov["x", "p"] = Cc + I ℏ/2; cov["p", "x"] = Cc - I ℏ/2;
wick[{}] := 1;
wick[l_ /; OddQ@Length@l] := 0;
wick[l_] := Sum[cov[First@l, l[[k]]] wick[Delete[l, {{1}, {k}}]], {k, 2, Length@l}];
```

In other words, `wick` pairs the first fluctuation with each later one, multiplies by that pair's covariance, and recurses on what is left. Now the map from a canonical operator monomial to its Gaussian expectation. The canonical form from `NonCommutativeExpand` writes repeated factors as `GeneralizedPower[NonCommutativeMultiply, g, k]`, so `letters` flattens a monomial into its list of $\hat x$'s and $\hat p$'s, `monExp` splits every letter into mean plus fluctuation and Wick-contracts the fluctuations, and `gExp` applies that to a whole expanded polynomial:

```wl
ClearAll[letters, monExp, gExp];
letters[HoldPattern[NonCommutativeMultiply[a__]]] := Join @@ (letters /@ {a});
letters[GeneralizedPower[NonCommutativeMultiply, g_, k_]] := ConstantArray[g, k];
letters[xo] := {xo}; letters[po] := {po};
monExp[ls_] := Total[(Times @@ MapThread[If[#2 == 0, #1 /. {xo -> mx, po -> mp}, 1] &, {ls, #}]) *
      wick[Pick[ls, #, 1] /. {xo -> "x", po -> "p"}] & /@ Tuples[{0, 1}, Length@ls]];
gExp[e_] := Expand[nc[e]] /. {
    m : (_NonCommutativeMultiply | _GeneralizedPower) :> monExp[letters@m],
    xo -> mx, po -> mp};
```

Test it on a mixed moment and a pure power, to see both the noncommuting piece and the Gaussian factorization:

```wl
{gExp[xo ** po], gExp[xo ** xo ** xo]}
```

As expected, $\langle \hat x\hat p\rangle = \mu_x\mu_p + C + i\hbar/2$, with the imaginary part exactly the commutator's contribution, and the third moment $\langle \hat x^3\rangle = \mu_x^3 + 3\mu_x V_x$ has no independent content, which is what "Gaussian" means. The pattern-dispatched `letters` is doing the whole job of walking the canonical monomials, and it never has to know how they were ordered, because the algebra already ordered them.

### The Conditional Variances: The Noise Cancels

Now the payoff. The general rule for how any observable evolves under the stochastic equation is

$$d\langle \hat O\rangle = \frac{i}{\hbar}\langle[\hat H,\hat O]\rangle\,dt - \frac{\lambda}{2}\langle[\hat x,[\hat x,\hat O]]\rangle\,dt + \sqrt{\lambda}\,\langle\{\hat A,\hat O\}\rangle\,dW, \qquad \hat A = \hat x - \langle\hat x\rangle.$$

In other words, an observable picks up a deterministic piece from the Hamiltonian, a second deterministic piece from the measurement, and a random piece proportional to how strongly it correlates with position. Encode that drift operator and that noise operator as two rules:

```wl
ClearAll[sseDrift, sseNoise];
sseDrift[o_, h_] := (I/ℏ) comm[h, o] - (λ/2) comm[xo, comm[xo, o]];
sseNoise[o_] := Sqrt[λ] (xo ** o + o ** xo);
```

Remember that Itô calculus adds a cross term to every product, so $d(\langle x\rangle^2) = 2\langle x\rangle\,d\langle x\rangle + (d\langle x\rangle)^2$. Apply the rules to the five moments that close on themselves for a free particle, take their Gaussian expectations, and assemble the three variance equations, showing drift and $dW$-coefficient side by side:

```wl
obs = {xo, po, xo ** xo, po ** po, (xo ** po + po ** xo)/2};
drift = gExp[sseDrift[#, po ** po/(2 m)]] & /@ obs;
noise = (gExp[sseNoise[#]] - 2 Sqrt[λ] mx gExp[#]) & /@ obs;
{dmx, dmp} = drift[[;; 2]]; {nmx, nmp} = noise[[;; 2]];
Simplify @ {
   {"dVx", drift[[3]] - 2 mx dmx - nmx^2, noise[[3]] - 2 mx nmx},
   {"dC ", drift[[5]] - mx dmp - mp dmx - nmx nmp, noise[[5]] - mx nmp - mp nmx},
   {"dVp", drift[[4]] - 2 mp dmp - nmp^2, noise[[4]] - 2 mp nmp}} // TableForm
```

Have you noticed something interesting? The third column, the coefficient of $dW$, is zero in all three rows. The means $\langle\hat x\rangle$ and $\langle\hat p\rangle$ wander randomly, but the *conditional variances are completely deterministic*: they obey ordinary differential equations,

$$\dot V_x = \frac{2C}{m} - 4\lambda V_x^2, \qquad \dot C = \frac{V_p}{m} - 4\lambda C V_x, \qquad \dot V_p = \lambda(\hbar^2 - 4C^2).$$

In other words, every experimental run produces a differently-placed wave packet, but all of them have the identical width at the identical time. That is a strong and very checkable prediction, and in Part V we will check it to fourteen digits. Note also that $d\langle \hat p\rangle$ has zero drift, so $\langle\hat p\rangle$ is a martingale: momentum information is not created by watching position, only shuffled between runs.

### The Steady State: Where Localization and Spreading Balance

The Riccati system has a fixed point. Find it, keeping everything symbolic:

```wl
fx = First@Select[
   Solve[{2 Cc/m - 4 λ Vx^2 == 0, Vp/m - 4 λ Cc Vx == 0, λ (ℏ^2 - 4 Cc^2) == 0},
    {Vx, Cc, Vp}],
   TrueQ@Simplify[(Vx /. #) > 0 && (Vp /. #) > 0, ℏ > 0 && m > 0 && λ > 0] &];
Simplify[{Sqrt[Vx], Cc, Sqrt[Vp], Vx Vp - Cc^2} /. fx, ℏ > 0 && m > 0 && λ > 0]
```

As one can see, a continuously watched free particle relaxes to a wave packet of fixed width

$$\sigma_x = \frac{1}{\sqrt2}\left(\frac{\hbar}{\lambda m}\right)^{1/4},$$

and the last entry says $V_xV_p - C^2 = \hbar^2/4$ exactly: the packet is a *minimum-uncertainty* state. Spreading and localization do not merely balance on average, they balance so precisely that the residual uncertainty is saturated. This is the quantitative content of the statement that measurement makes a particle classical, and it is the single most useful formula in this document, because $\sigma_x$ is what your spatial grid has to resolve.

Check the two limits that must hold. Watching infinitely hard should pin the particle; the classical limit should too, since $\hbar$ sets the scale of the residual spread:

```wl
{Limit[Sqrt[Vx /. fx], λ -> Infinity], Limit[Sqrt[Vx /. fx], ℏ -> 0],
 Vx /. First@DSolve[{Vx'[t] == 2 Cc[t]/m, Cc'[t] == Vp[t]/m, Vp'[t] == 0,
      Vx[0] == v0, Cc[0] == 0, Vp[0] == ℏ^2/(4 v0)}, {Vx, Cc, Vp}, t]}
```

Therefore $\sigma_x \to 0$ in both the Zeno limit $\lambda\to\infty$ and the classical limit $\hbar\to0$, and at $\lambda = 0$ the same Riccati system with the measurement switched off returns textbook free spreading, $V_x(t) = v_0 + \hbar^2 t^2/(4m^2v_0)$. Note that I extracted the solution by *name* rather than by position: `DSolve` returns its rules alphabetically, so here they arrive as `Cc, Vp, Vx` and indexing `[[1,1,2]]` would silently hand you $C(t)$ while you believe you are holding $V_x(t)$.

Before we move on to the simulation, let us summarize what we have established exactly, with no numerics at all:

- The Gaussian Kraus family $\hat M(\bar x)\propto e^{-\lambda \Delta t(\hat x-\bar x)^2}$ is complete, and its outcome law is $|\psi|^2$ blurred by a Gaussian of variance $1/(4\lambda\Delta t)$.
- The stochastic equation preserves the norm identically, for any Hamiltonian and any measured observable, in any dimension.
- The ensemble heating rate is exactly $\frac{d}{dt}\mathbb{E}\langle\hat p^2\rangle = \lambda\hbar^2$.
- The conditional variances are deterministic: all three $dW$ coefficients vanish identically.
- The steady state is a minimum-uncertainty packet, $V_xV_p - C^2 = \hbar^2/4$, of width $\sigma_x = (\hbar/\lambda m)^{1/4}/\sqrt2$, approached at rate $4\lambda V_x = 2\sqrt{\hbar\lambda/m}$.

Every one of those is a trap laid for the simulation we are about to write.

## Part III: The Right Shelf Is `ItoProcess`

### Motivation by Failure: The Solver That Cannot See the Noise

The reflex, faced with a differential equation in Wolfram Language, is `NDSolve`. Let us find out what happens, rather than assume:

```wl
NDSolveValue[{\[DifferentialD]x[t] == -x[t] \[DifferentialD]t + \[DifferentialD]w[t], x[0] == 1},
  x, {t, 0, 1}]
```

As one can see, this fails with `NDSolveValue::underdet`, and the message is more informative than a simple refusal: it reports "more dependent variables, `{w[t], x[t]}`, than equations." `NDSolve` has parsed $dw$ as a *second unknown function* to be solved for, because nothing in its model of the world says $w$ is a stochastic process. A stochastic differential equation is not a differential equation with an extra term; it is a different object, and it lives on a different shelf.

### Units and the Kinetic Operator: Setting the Numerical Stage

The right shelf is `ItoProcess`, which does know what $dw$ means. Everything from here on is numeric, so fix units $\hbar = m = 1$; the symbolic scalings of Part II now take numbers:

```wl
ℏ = 1.; mass = 1.;
```

To feed the equation to `ItoProcess` we discretize position on a grid and need the kinetic operator on that grid. A finite-difference second derivative with periodic wrap keeps it local, hence sparse. Define it as a function of the grid size and spacing:

```wl
ClearAll[fdKinetic];
fdKinetic[n_, dx_] := -ℏ^2/(2 mass dx^2) Normal@SparseArray[
    {{i_, i_} -> -2., {i_, j_} /; Abs[i - j] == 1 -> 1., {1, n} -> 1., {n, 1} -> 1.}, {n, n}];
```

### The Grid SSE as an Itô Process: A Real State With Scalar Noise

Now the model itself. We represent $|\psi\rangle$ by the real and imaginary parts of its grid values, stack them into one real state vector, and write the whole stochastic equation as a real Itô system driven by a single scalar Wiener process:

```wl
ClearAll[gridProcess];
gridProcess[n_, ll_, λ_, ψ0_] := Module[{dx, xs, kmat, us, vs, xav, a, drift, diff},
   dx = ll/n; xs = dx Range[-n/2, n/2 - 1]; kmat = fdKinetic[n, dx];
   us = Table[Unique["u"], n]; vs = Table[Unique["v"], n];
   xav = xs . (us^2 + vs^2); a = xs - xav;
   drift = Join[kmat . vs/ℏ - (λ/2) a^2 us, -kmat . us/ℏ - (λ/2) a^2 vs];
   diff = Transpose[{Sqrt[λ] Join[a us, a vs]}];
   ItoProcess[{drift, diff}, {Join[us, vs], Join[Re@ψ0, Im@ψ0]}, {t, 0}]];
```

Two idioms here earn their place. The state variables come from `Table[Unique["u"], n]`, which makes genuine plain symbols: `ItoProcess` insists on symbols, not indexed expressions like `u[1]`, and feeding it the latter triggers a `Function::flpar` complaint. And the diffusion is a $2n\times1$ matrix, a real state driven by a single scalar noise, which is exactly the shape that unlocks `RandomFunction`'s highest-order scalar-noise integrator.

### The Position Grid: Reading the Momentum Ordering Aloud

We also need a periodic position grid and a way to put a Gaussian wave packet on it. The grid carries both the positions and the FFT-ordered momenta the kinetic Fourier step will want later. Define both:

```wl
ClearAll[grid, gaussian];
grid[n_, ll_] := With[{dx = ll/n},
   <|"n" -> n, "dx" -> dx, "x" -> dx Range[-n/2, n/2 - 1],
     "p" -> (2 Pi ℏ/ll) RotateLeft[Range[-n/2, n/2 - 1], n/2]|>];
gaussian[g_, x0_, p0_, s_] := With[{w = Exp[-(g["x"] - x0)^2/(4 s^2) + I p0 g["x"]/ℏ]},
   w/Sqrt[w . Conjugate[w]]];
```

Look at the momentum list on a small grid before trusting it anywhere:

```wl
grid[8, 4.]["p"]
```

As one can see, the momenta are not in ascending order: they run $0, +1, +2, +3$ then wrap to $-4, -3, -2, -1$, in units of $2\pi\hbar/L$. That scrambled order is the convention `Fourier` expects, with the negative frequencies in the second half, and the `RotateLeft` is what places our physical momenta into it. Getting this ordering wrong is the classic way a spectral kinetic step runs without error and propagates the packet in the wrong direction.

### A Trajectory You Can Watch: One Run of the Built-In Integrator

With the process and the grid in hand, we need only to read a trajectory back. A raw `ItoProcess` sample is a real state vector; reassemble it into a wavefunction and read off the conditional mean position:

```wl
ClearAll[stateToψ, condX];
stateToψ[s_] := With[{n = Length[s]/2}, s[[;; n]] + I s[[n + 1 ;;]]];
condX[s_, xs_] := With[{p = Abs[stateToψ[s]]^2}, xs . p/Total[p]];
```

Now simulating is one call to `RandomFunction`, which returns a `TemporalData` object carrying every path. Evolve one free trajectory and watch the conditional mean $\langle\hat x\rangle$ random-walk as the detector localizes the packet:

```wl
With[{n = 32, ll = 16., λ = 1., tf = 3., dt = 0.002},
 {xs = (ll/n) Range[-n/2, n/2 - 1],
  td = RandomFunction[gridProcess[n, ll, λ, gaussian[grid[n, ll], 0., 0., 1.]],
     {0., tf, dt}, Method -> "StochasticRungeKuttaScalarNoise"]},
 ListLinePlot[Transpose[{td["Times"], condX[#, xs] & /@ td["ValueList"][[1]]}],
  Frame -> True, GridLines -> Automatic, ImageSize -> 480, AspectRatio -> 1/2,
  FrameLabel -> {"time t", "\[LeftAngleBracket]x\[RightAngleBracket]"}, PlotRange -> All,
  PlotLabel -> "One conditional trajectory: \[LeftAngleBracket]x\[RightAngleBracket] wanders as the packet localizes"]]
```

That plot is the single most important picture in continuous measurement: a pure state, conditioned on the noisy measurement record $Y_t$, whose center wanders while its width collapses to $\sigma_x$. The built-in stochastic machinery produced it with no integrator written by hand.

### The Ensemble Reproduces the Lindblad Equation: Two Routes That Must Agree

The defining property of the whole construction is that averaging the trajectories must reproduce the Lindblad equation. Test it against a completely different computational route, the matrix exponential of the Liouvillian, so that agreement is evidence and not tautology:

```wl
With[{n = 24, ll = 12., λ = 1., tf = 0.4, dt = 0.0005, ntraj = 400},
 {g = grid[n, ll]},
 {xm = DiagonalMatrix[g["x"]], kmat = fdKinetic[n, ll/n], id = IdentityMatrix[n],
  ψ0 = gaussian[g, 1., 1.2, 0.9]},
 {x2ito = Mean[Function[s, With[{p = Abs[stateToψ@s]^2}, g["x"]^2 . p/Total[p]]] /@
      RandomFunction[gridProcess[n, ll, λ, ψ0], {0., tf, dt}, ntraj,
         Method -> "StochasticRungeKuttaScalarNoise"]["ValueList"][[All, -1]]],
  rhoL = ArrayReshape[MatrixExp[tf (
        -(I/ℏ) (KroneckerProduct[kmat, id] - KroneckerProduct[id, Transpose@kmat])
        - (λ/2) (KroneckerProduct[xm . xm, id] - 2 KroneckerProduct[xm, xm]
           + KroneckerProduct[id, xm . xm]))] . Flatten@Outer[Times, ψ0, Conjugate@ψ0],
     {n, n}]},
 AssociationThread[{"<x^2> ensemble", "<x^2> Lindblad", "rel. diff", "Monte-Carlo scale"},
  {x2ito, Re@Tr[rhoL . xm . xm], Abs[x2ito/Re@Tr[rhoL . xm . xm] - 1], 1./Sqrt[ntraj]}]]
```

As one can see, the ensemble average agrees with the Lindblad reference up to Monte-Carlo scatter. Do not expect it to hit $1/\sqrt{\text{ntraj}}$ exactly: $\langle x^2\rangle$ is a variance-like estimator with a broader sampling distribution than the mean, so at $\text{ntraj} = 400$ a few percent is normal. Push `ntraj` up and the gap closes. The built-in route is correct.

### Why the Library Route Is Not the End: Two Visible Symptoms

So `ItoProcess` gives a correct, honest simulation with essentially no bespoke code. Why write anything more? Two reasons, and both are visible symptoms rather than opinions.

The first is that a general SDE integrator does not know our equation preserves the norm, so it does not. Read the norm-squared at the start and the end of a run:

```wl
With[{n = 32, ll = 16., λ = 3., tf = 1., dt = 0.002},
 {td = RandomFunction[gridProcess[n, ll, λ, gaussian[grid[n, ll], 0., 0., 1.]],
    {0., tf, dt}, Method -> "StochasticRungeKuttaScalarNoise"]},
 {Total[td["ValueList"][[1, 1]]^2], Total[td["ValueList"][[1, -1]]^2]}]
```

The norm-squared starts at $1$ and drifts, and the drift grows with $\lambda$, with the step size, and with time. You can renormalize by hand after each step, but that patches the symptom.

The deeper reason is the second one, and Part V will make it quantitative: the coefficients of this equation are *unbounded* in $x$, the drift carrying $(x-\langle x\rangle)^2$, and an explicit scheme applied to an unbounded generator is not unconditionally stable. The exponential of that generator is. That is the structural fact the next part exploits, and it is the difference between a simulation that runs and one you can trust to fourteen digits.

## Part IV: A Propagator That Exponentiates the Measurement

### Diagonal in Position: The One Fact That Picks the Integrator

Look again at the equation and ask what kind of operators appear. The Hamiltonian mixes positions and momenta. But *every measurement term*, both the $(\hat x - \langle\hat x\rangle)^2$ drift and the $(\hat x - \langle \hat x\rangle)$ noise, is a multiplication operator in $x$. On a position grid it is a diagonal matrix, or simply a list you multiply pointwise.

This is not a cosmetic observation, it is the whole method. Because the measurement terms commute among themselves, the measurement channel over a *finite* interval is exactly the Gaussian Kraus operator we started from: two measurements of strength $\lambda\,\Delta t/2$ compose exactly into one of strength $\lambda\,\Delta t$. So a scheme that applies the measurement as a Kraus multiplication carries **no discretization error at all from the measurement**, at any step size, and multiplies by the bounded factor $e^{-\lambda\,dt\,(x-\bar x)^2}$ rather than the unbounded $1 - \frac{\lambda}{2}(x-\langle x\rangle)^2 dt$. The only error left is the failure of $\hat H$ and $\hat x^2$ to commute.

Equivalently: the equation is telling you how to integrate it. Split off the part you can do exactly, and spend your step size only on the part you cannot. The lesson that transfers to any problem is the question, not the answer: *which part of my generator can I exponentiate exactly?*

### The Moment Readouts: One Line Each

We will need to read $\langle\hat x\rangle$, the conditional variance $V_x$, and $\langle\hat p^2\rangle$ off a grid wavefunction. Position moments are dot products against $|\psi|^2$; the momentum moment goes through the Fourier transform. Define the three:

```wl
ClearAll[xExp, vxOf, p2Exp];
xExp[g_, ψ_] := g["x"] . Abs[ψ]^2;
vxOf[g_, ψ_] := g["x"]^2 . Abs[ψ]^2 - xExp[g, ψ]^2;
p2Exp[g_, ψ_] := g["p"]^2 . Abs[Fourier[ψ, FourierParameters -> {0, -1}]]^2;
```

### The Unitary Substeps: Kinetic in Fourier, Potential Pointwise

We represent $\psi$ on a uniform periodic grid, so the kinetic term is a multiplication in the Fourier domain and the potential term is a pointwise phase. Each is a curried operator: give it the grid and the step size and it returns a function that maps a wavefunction to a wavefunction:

```wl
ClearAll[kinStep, potStep];
kinStep[g_, dt_][ψ_] := InverseFourier[
   Exp[-I g["p"]^2 dt/(2 mass ℏ)] Fourier[ψ, FourierParameters -> {0, -1}],
   FourierParameters -> {0, -1}];
potStep[g_, pot_, dt_][ψ_] := Exp[-I pot dt/ℏ] ψ;
```

That currying is what will let us compose the substeps without an explicit loop. I will not ask you to trust that the Fourier convention is right; propagate a free Gaussian and compare against the closed-form spreading $V_x(t) = s_0^2 + (\hbar t/2ms_0)^2$ from Part II:

```wl
With[{g = grid[512, 60.], s0 = 1., tf = 2., p0 = 1.5},
 {ψT = Nest[kinStep[g, tf/400], gaussian[g, -5., p0, s0], 400]},
 AssociationThread[{"Vx exact", "Vx grid", "<x> exact", "<x> grid"},
  {s0^2 + (ℏ tf/(2 mass s0))^2, vxOf[g, ψT], -5. + p0 tf/mass, xExp[g, ψT]}]]
```

The variance matches the closed form to thirteen digits and the packet has moved to the right, which is the check that actually matters: a sign error in the Fourier convention would send it left while leaving the variance perfectly correct. Whenever you validate a propagator, test something that the error you fear would actually break.

### The Measurement Substep: Exponentiate the Backaction

The measurement step is the whole quantum backaction of watching the particle. Draw an outcome $\bar x$ from the exact outcome law of Part I, then multiply pointwise by the Gaussian Kraus factor and renormalize:

```wl
ClearAll[drawExact, measStep];
drawExact[g_, λ_, dt_][ψ_] := RandomChoice[Abs[ψ]^2 -> g["x"]] +
   RandomVariate[NormalDistribution[0, 1/Sqrt[4 λ dt]]];
measStep[g_, λ_, dt_][ψ_] := With[{xb = drawExact[g, λ, dt][ψ]},
   {v = Exp[-λ dt (g["x"] - xb)^2] ψ}, v/Sqrt[v . Conjugate[v]]];
```

That Kraus factor can only ever damp the wavefunction, never amplify it. Keep that boundedness in mind: in Part V it is the reason this scheme is stable and the explicit one is not. Before composing anything, sample the detector readings $\bar x$ directly and check them against the outcome law of Part I, closing the factor-of-two worry by sampling rather than by algebra:

```wl
With[{g = grid[256, 40.], λ = 1., dt = 0.02},
 {ψ0 = gaussian[g, 0.7, 0., 1.]},
 {draws = Table[drawExact[g, λ, dt][ψ0], {200000}]},
 AssociationThread[{"mean predicted", "mean sampled", "variance predicted", "variance sampled"},
  {xExp[g, ψ0], Mean[draws], vxOf[g, ψ0] + 1/(4 λ dt), Variance[draws]}]]
```

As one can see, the sampled mean lands on $\langle\hat x\rangle$ and the sampled variance on $V_x + 1/(4\lambda\Delta t)$, the conditional spread of the packet plus the detector's own blur. This is the same $1/(4\lambda\Delta t)$ we carried symbolically in Part I, now confirmed by two hundred thousand draws.

### Assembling the Step: Composing Without a Loop

A full step is a composition of substeps. `RightComposition` reads left-to-right in application order, and `Nest` runs the trajectory, so there is no explicit iteration anywhere. Build the symmetric assembly, which splits the measurement in half and wraps it around the unitary step, together with a runner that iterates it:

```wl
ClearAll[stepSym, evolve];
stepSym[g_, pot_, λ_, dt_] := RightComposition[
   measStep[g, λ, dt/2], potStep[g, pot, dt/2], kinStep[g, dt],
   potStep[g, pot, dt/2], measStep[g, λ, dt/2]];
evolve[step_, tf_, dt_][ψ0_] := Nest[step, ψ0, Round[tf/dt]];
```

Two things have to be right about any experiment on this propagator, and I got both wrong on my first attempt. The first is *when* to measure the steady state: the Riccati system relaxes at rate $4\lambda V_x = 2\sqrt{\hbar\lambda/m}$, so comparing too early plateaus on the transient rather than on the step-size error. We measure at $t = 8$, many relaxation times in. The second is *how wide* the box must be: the conditional mean random-walks, and because $\langle\hat p\rangle$ random-walks too and $\langle\hat x\rangle$ integrates it, the spread grows like $\sqrt{t^3/3 + t}$, not $\sqrt t$. On a periodic box the packet must not wrap, because $e^{-\lambda dt(x-\bar x)^2}$ is not periodic and wrapping corrupts $V_x$ silently. So the box is sized for the $t^3$ wander, not for the packet width.

With both settled, measure the order by comparing the naive assembly (measurement last) against the symmetric one:

```wl
With[{g = grid[1024, 160.], λ = 1., tf = 8.},
 {ψ0 = gaussian[g, 0., 0., 1.]},
 TableForm[
  Table[{dt,
    Abs[vxOf[g, evolve[RightComposition[potStep[g, 0. g["x"], dt/2], kinStep[g, dt],
          potStep[g, 0. g["x"], dt/2], measStep[g, λ, dt]], tf, dt][ψ0]]/0.5 - 1],
    Abs[vxOf[g, evolve[stepSym[g, 0. g["x"], λ, dt], tf, dt][ψ0]]/0.5 - 1]},
   {dt, {0.008, 0.004, 0.002}}],
  TableHeadings -> {None, {"dt", "measurement last", "symmetrised"}}]]
```

Halving the step size halves the error in the first column but quarters it in the second: the naive assembly is first order and the symmetrized one is second order, for one extra random draw per step. The reason generalizes to every operator-splitting scheme you will ever write. Splitting error comes from the non-commutation of the pieces, and a *symmetric* arrangement makes the leading commutator cancel between the two halves. Putting the measurement last leaves it unpaired.

### The Record Is the Data, the State Is the Inference

Everything simulated so far is the conditional state, and no laboratory has one. What an experiment produces is the sequence of readings $\bar x_k$, and the state is what you compute from them. `measStep` already draws that number and discards it. Keep it instead, by splitting the draw from the deterministic update it implies:

```wl
ClearAll[measAt, stepAt, filter];
measAt[g_, λ_, dt_, xb_][ψ_] := With[{v = Exp[-λ dt (g["x"] - xb)^2] ψ}, v/Sqrt[v . Conjugate[v]]];
stepAt[g_, λ_, dt_, xb_][ψ_] := kinStep[g, dt][measAt[g, λ, dt, xb][ψ]];
filter[g_, λ_, dt_, rec_][ψ0_] := FoldList[stepAt[g, λ, dt, #2][#1] &, ψ0, rec];
```

This is the finite-bin update of Part I applied literally: one reading per bin, $\hat M(\bar x_k)$ and then $\hat U$. Once $\bar x_k$ is fixed, `stepAt` contains no randomness at all, since `RandomChoice` and `RandomVariate` live only in `drawExact`. Experiment and inference are therefore two different traversals of the same step. The experiment draws its own readings and emits data; the filter consumes data and emits states, and it is a plain `FoldList`, which is the whole content of quantum filtering.

Run one experiment. We keep its final state only so that the reconstruction can be checked against it, which is a luxury no experimentalist has:

```wl
SeedRandom[20260727];
gRec = grid[256, 40.]; λRec = 1.; dtRec = 0.005; ψRec = gaussian[gRec, 0., 0., 1.];
{ψEnd, {rec}} = Reap[Nest[Function[ψ,
    With[{xb = drawExact[gRec, λRec, dtRec][ψ]}, Sow[xb]; stepAt[gRec, λRec, dtRec, xb][ψ]]],
   ψRec, 800]];
```

Now reconstruct from `rec` alone, both from the true initial state and from a deliberately wrong guess placed three units away with twice the width:

```wl
tru = filter[gRec, λRec, dtRec, rec][ψRec];
gue = filter[gRec, λRec, dtRec, rec][gaussian[gRec, 3., 0., 2.]];
AssociationThread[
 {"detector noise sd", "packet width Sqrt[Vx]", "signal to noise per reading",
  "replay error max|dψ|", "overlap of wrong start with truth, initial", "overlap, final",
  "Vx, true trajectory", "Vx, wrong start", "predicted Vx"},
 {1/Sqrt[4 λRec dtRec], Sqrt[vxOf[gRec, Last@tru]], Sqrt[4 λRec dtRec vxOf[gRec, Last@tru]],
  Max[Abs[Last[tru] - ψEnd]], Abs[Conjugate[First@tru] . First@gue]^2,
  Abs[Conjugate[Last@tru] . Last@gue]^2, vxOf[gRec, Last@tru], vxOf[gRec, Last@gue],
  Sqrt[ℏ/(4 λRec mass)]}]
```

Three things to read, in the order that matters for an experiment.

**The record is not the trajectory.** The detector noise has standard deviation $1/\sqrt{4\lambda\,\Delta t}\approx7.1$, against a packet whose entire width is $\sqrt{V_x}\approx0.71$. One reading therefore carries a signal-to-noise ratio near $0.1$, and that figure is a property of the detector rather than of the run: the wander of $\langle\hat x\rangle$ accumulates with time, so any ratio formed against it would depend only on how long you watched. The smooth curve plotted earlier is an inference and never an observation.

**The state is a deterministic function of the record.** The replay error is exactly zero, not merely small: feeding the readings back through `filter` returns the original trajectory bit for bit. Every random number in the theory is spent producing the data, and none is left in the map from data to state.

**The filter forgets a wrong initial state.** Started from the wrong place and the wrong width, and fed the same readings, its overlap with the truth climbs from about $0.33$ to $0.9997$, and its conditional variance agrees with the true trajectory's to $7\times10^{-5}$. Both then sit about $0.4\%$ above the analytic $V_x=\sqrt{\hbar/(4\lambda m)}$, which is the splitting error of this assembly rather than any failure of the filter: `stepAt` applies the measurement last, so it is the first-order scheme measured just above, and the gap falls to $0.04\%$ at $\Delta t=0.001$. This is what makes the equation usable in a laboratory, where the initial state is exactly what you do not know.

The relation is one picture: the data, and the single trajectory it determines.

```wl
Show[
 ListPlot[Transpose[{dtRec Range[Length@rec], rec}],
  PlotStyle -> Directive[Gray, Opacity[0.35], PointSize[0.004]]],
 ListLinePlot[Transpose[{dtRec Range[0, Length@rec], xExp[gRec, #] & /@ tru}],
  PlotStyle -> Directive[Thick, Red]],
 Frame -> True, Axes -> False, ImageSize -> 520, AspectRatio -> 1/2,
 PlotRange -> {{0, dtRec Length@rec}, All},
 FrameLabel -> {"time t", "reading, and the inferred mean"},
 PlotLabel -> "800 detector readings and the one trajectory they determine"]
```

## Part V: Trusting the Answer

A number produced by a simulation is a hypothesis. This part is the ladder that turns it into a result, and I have ordered the rungs by how much they can catch.

### The Variances Are Deterministic: A Test With No Statistics In It

Part II predicted that all runs share the identical conditional variance. That is an unusually sharp claim, because it needs no averaging at all: two runs with different noise must agree to machine precision. Run six independent trajectories on a fine grid and compare their widths:

```wl
With[{g = grid[256, 40.], λ = 1., dt = 0.002},
 {runs = Table[vxOf[g, evolve[stepSym[g, 0. g["x"], λ, dt], 4., dt][gaussian[g, 0., 0., 1.]]], {6}]},
 AssociationThread[{"values", "spread across runs", "predicted"},
  {runs, StandardDeviation[runs], Sqrt[ℏ/(4 λ mass)]}]]
```

Six different noise realizations agree to roughly fourteen digits while landing on the predicted steady width. This one test simultaneously confirms the Riccati derivation, the sampling law, and the propagator, and it does so without a single statistical error bar.

### The Same Test Is Where the Library Route Shows Its Seams

That sharpness gives us a clean way to see what the exact-Kraus step buys over the general integrator of Part III. The split-step keeps the state *exactly* Gaussian, so its variance is deterministic. A general explicit integrator does not, so its variance both scatters across runs and drifts off the predicted value. Put the two side by side on the same grid:

```wl
With[{n = 32, ll = 16., λ = 1., tf = 2., dt = 0.002, nrun = 5},
 {g = grid[n, ll]},
 {ψ0 = gaussian[g, 0., 0., 1.]},
 {ss = Table[vxOf[g, evolve[stepSym[g, 0. g["x"], λ, dt], tf, dt][ψ0]], {nrun}],
  itψs = stateToψ /@ RandomFunction[gridProcess[n, ll, λ, ψ0], {0., tf, dt}, nrun,
      Method -> "StochasticRungeKuttaScalarNoise"]["ValueList"][[All, -1]]},
 AssociationThread[{"split-step spread", "ItoProcess spread", "split mean", "Ito mean", "predicted"},
  {StandardDeviation[ss], StandardDeviation[vxOf[g, #] & /@ itψs],
   Mean[ss], Mean[vxOf[g, #] & /@ itψs], Sqrt[ℏ/(4 λ mass)]}]]
```

On the same coarse grid, the split-step widths agree to a part in $10^{4}$ (limited only by the discreteness of `RandomChoice` on $32$ points) while the `ItoProcess` widths scatter by a thousand times more and sit biased below the predicted $0.5$. Both routes are simulating the same equation. Only the one that exponentiates the measurement keeps the exact Gaussianity the physics guarantees. This is the concrete payoff of asking which part of the generator you can exponentiate.

### Two Representations Agree: The Trajectory Ensemble Against the Lindblad Equation

We already saw the ensemble reproduce the Lindblad equation for the `ItoProcess` route in Part III. Confirm it holds for the split-step too, this time also watching purity fall as the ensemble decoheres while every single trajectory stays exactly pure. First the reference Hamiltonian as a matrix, built column by column from the spectral kinetic operator:

```wl
ClearAll[kinMat];
kinMat[g_] := Transpose[InverseFourier[g["p"]^2/(2 mass) Fourier[#, FourierParameters -> {0, -1}],
     FourierParameters -> {0, -1}] & /@ IdentityMatrix[g["n"]]];
```

Now build the Liouvillian, exponentiate it for the Lindblad reference, average three thousand split-step trajectories for the unravelling, and compare $\langle\hat x\rangle$, $\langle\hat x^2\rangle$, and the purity:

```wl
With[{nn = 32, ll = 12., λ = 1., tf = 0.4, dt = 0.0005, ntraj = 3000},
 {g = grid[nn, ll]},
 {id = IdentityMatrix[nn], xm = DiagonalMatrix[g["x"]], hm = kinMat[g], ψ0 = gaussian[g, -1., 1.2, 0.9]},
 Module[{liou, rhoL, ens, obs},
  liou = (
    -(I/ℏ) (KroneckerProduct[hm, id] - KroneckerProduct[id, Transpose[hm]])
    - (λ/2) (KroneckerProduct[xm . xm, id] - 2 KroneckerProduct[xm, xm]
       + KroneckerProduct[id, xm . xm]));
  rhoL = ArrayReshape[MatrixExp[liou tf] . Flatten@Outer[Times, ψ0, Conjugate@ψ0], {nn, nn}];
  ens = Mean[Table[Outer[Times, #, Conjugate@#] &@
       evolve[stepSym[g, 0. g["x"], λ, dt], tf, dt][ψ0], {ntraj}]];
  obs[r_] := {Re@Tr[r . xm], Re@Tr[r . xm . xm], Re@Tr[r . r]};
  TableForm[Transpose[{{"<x>", "<x^2>", "purity"}, obs[rhoL], obs[ens],
     Abs[obs[ens]/obs[rhoL] - 1]}],
   TableHeadings -> {None, {"observable", "Lindblad", "3000 trajectories", "rel. diff"}}]]]
```

The `hm` here is an entirely different construction from the split-step trajectories, so the reference and the ensemble share nothing but the physics. All three agree to within the Monte-Carlo scatter of $1/\sqrt{3000} \approx 1.8\%$, including the purity, which falls from $1$ to about $0.65$ while every individual trajectory stays exactly pure. That contrast is the whole conceptual content of an unravelling, and here you watch both sides of it in one cell.

### The Failure Mode: Heating Off the Edge of the Grid

Every simulation has a regime where it lies, and you want to meet yours deliberately. Recall the exact heating law: $\mathbb{E}\langle p^2\rangle$ grows as $\lambda\hbar^2 t$ without bound. But a grid of spacing $\Delta x$ can only represent momenta up to the Nyquist value $p_{\max} = \pi\hbar/\Delta x$. Sooner or later the packet heats past what the grid can hold. Watch a deliberately coarse grid meet that limit, against a fine one and against the exact law:

```wl
With[{λ = 40., dt = 0.0005},
 {coarse = grid[64, 40.], fine = grid[512, 40.]},
 TableForm[
  Table[{tt, 0.25 + λ ℏ^2 tt,
    Mean[Table[p2Exp[coarse, evolve[stepSym[coarse, 0. coarse["x"], λ, dt], tt, dt][
         gaussian[coarse, 0., 0., 1.]]], {40}]],
    Mean[Table[p2Exp[fine, evolve[stepSym[fine, 0. fine["x"], λ, dt], tt, dt][
         gaussian[fine, 0., 0., 1.]]], {40}]]}, {tt, {0.25, 0.5, 1.0}}],
  TableHeadings -> {None, {"t", "exact law", "N=64 (p_max^2=25)", "N=512"}}]]
```

The coarse grid does not diverge, does not warn, and does not stop. It saturates around $8$ and reports a stable, plausible, completely wrong number, while the fine grid tracks the linear growth up to Monte-Carlo scatter. Note that the saturation sits below the Nyquist ceiling $p_{\max}^2 = 25$, because aliasing folds probability back to low momenta rather than piling it up at the edge. This is the failure mode I find genuinely dangerous, because the wrong answer looks *more* stable than the right one. Refining the time step would never reveal it. Only the analytic law does. The practical rule that falls out: your grid must resolve $\sigma_x = (\hbar/\lambda m)^{1/4}/\sqrt2$ from below and accommodate $\sqrt{\lambda\hbar^2 t}$ from above, and those two requirements pull apart as $\lambda$ grows.

### Where This Leaves Us (and What Comes Next)

You now have a complete, computation-first toolkit for continuous measurement: a Gaussian detector built from its Kraus operator and validated by sampling, the built-in noncommutative algebra yielding the exact heating law and the deterministic Riccati system for the conditional variances, a closed-form minimum-uncertainty steady state to aim at, the built-in `ItoProcess` stack for an honest general-purpose trajectory simulation, a second-order split-step propagator that exponentiates the measurement channel exactly and beats the general integrator on the sharp benchmarks, and a verification ladder running from a statistics-free consistency test through a full Lindblad cross-check to a deliberate encounter with the scheme's own failure mode.

The natural continuations are all one modification away, and I would encourage you to try them before reading anything further. Put a double well in `pot` and watch a delocalized superposition collapse into one lobe, then check that it chooses each lobe with the Born probability. Set the detector efficiency below one, which breaks purity and forces you from a wavefunction to a stochastic master equation. Replace the sharp $\hat x$ by a spatially smeared version, which is what Diósi and the continuous-spontaneous-localization models do to keep the heating rate finite. Or measure two non-commuting observables at once and watch the conditional state stop being Gaussian.

The equation is short, and almost everything interesting about it is a consequence of one structural fact: the measurement is diagonal in position. Find the analogous fact in your problem, and the integrator will usually write itself.
