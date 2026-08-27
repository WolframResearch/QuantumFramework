---
Template: Default
Title: A Computational Guide to Open Quantum Systems: Lindblad Dynamics, Quantum Trajectories, and Feedback Control
Author: Mads Bahrami
---

# A Computational Guide to Open Quantum Systems: Lindblad Dynamics, Quantum Trajectories, and Feedback Control

<!-- #| style: Author -->
Mads Bahrami

<!-- #| style: Affiliation -->
Wolfram Research, Inc.

This computational essay is a hands-on tour of open quantum systems: systems that leak energy or information to their environment, watched by a detector on the output. We build a driven atom relaxing to its steady state, light in a leaky cavity, a superposition that decoheres almost the instant it forms, an atom whose fluorescence we read first as single photons and then as a homodyne current, a mechanical resonator cooled by measurement and feedback, and a dozen more. Each becomes a few lines of code you can run and change, and each ends in a picture.

This Wolfram notebook is a live animal, so evaluate the cells from top to bottom, change them, get errors, debug them, and practically do anything with it that comes to mind. This is the kind of document I wish I had when I was learning these topics as a graduate student. I remember reading many abstract stochastic dynamical equations and banging my head against the wall trying to truly understand what they meant, and now I am a little surprised (and perhaps disappointed) that I never tried exploring them the way I have done here in this notebook, or, more frankly, that nobody told me to learn by computing. Anyhow, as you read this notebook, you do not have to understand every tool to follow the physics: if you would rather learn the code from the examples, simply run the definitions in "The Toolkit: Everything the Twenty Examples Share" without studying them, and return to that part only when you want to understand the implementation. Look at each output and think about what it means before picking apart the code that produced it. And you are not locked into anything: change a number, watch a different system, reseed the randomness, and see what moves. That is the whole point.

I wrote this essay assuming you know three ideas. First, the state is a **density matrix** $\rho$: an $n\times n$ Hermitian, positive-semidefinite matrix ($n$ is number of levels, 2 for a qubit, more for an oscillator), with nonnegative level populations on its diagonal (summing to one) and the coherences off-diagonal. Second, the system has a Hamiltonian $\hat H$ that drives its unitary evolution when left alone. Third, one random quantity recurs below: the Wiener increment $dW$, a Gaussian random step with variance equal to the time step $dt$. Where the pace steepens I will say so.

We set $\hbar = 1$ throughout, so every rate and frequency is in units of inverse time.

Let's start!

## In Standard Notation: The Twenty Governing Equations

Before building anything, here is the catalog stated formally: the standard notation every example uses, then all twenty governing equations on one page. Skim it now and return to it later; nothing here needs memorizing on a first pass. The toolkit in the next section builds each of these objects as runnable code and the examples work them one at a time, so meeting the equations first means you already know the formal object each tool stands for.

### The machinery, in standard notation

Every example below opens with a one-line statement of its problem in standard notation, so the formal problem and the plain solution sit side by side. These are the symbols that statement uses. The state is the density matrix $\rho$; the averaged evolution is the **master equation**
$$\dot\rho = \mathcal{L}\rho = -i[\hat H,\rho] + \sum_k \mathcal{D}[\hat c_k]\rho, \qquad \hbar = 1,$$
with a Hamiltonian part $-i[\hat H,\rho]$ and one **dissipator** per leak channel,
$$\mathcal{D}[\hat c]\rho = \hat c\,\rho\,\hat c^\dagger - \tfrac12(\hat c^\dagger\hat c\,\rho + \rho\,\hat c^\dagger\hat c).$$
For a qubit the **Bloch vector** is $(x,y,z) = (\langle\hat\sigma_x\rangle,\langle\hat\sigma_y\rangle,\langle\hat\sigma_z\rangle)$, with $\rho = \tfrac12(\mathbb{1} + x\hat\sigma_x + y\hat\sigma_y + z\hat\sigma_z)$; the lowering operator is $\hat\sigma_- = |g\rangle\langle e|$, and for an oscillator the annihilation operator is $\hat a$ with $[\hat a,\hat a^\dagger] = 1$. When the record is kept, the state follows a **stochastic master equation** whose measurement (diffusion) term for a monitored $\hat c$ is $\mathcal{H}[\hat c]\rho = \hat c\rho + \rho\hat c^\dagger - \langle\hat c + \hat c^\dagger\rangle\rho$, driven by $dW$ with $dW^2 = dt$; photon counting instead applies the jump $\mathcal{G}[\hat c]\rho = \hat c\rho\hat c^\dagger/\mathrm{Tr}[\hat c\rho\hat c^\dagger] - \rho$ at click times.

### The Twenty Governing Equations, at a Glance

Here is the dynamics of every example in one place. All stochastic equations are Itô equations. For a Hermitian observable $\hat M$ monitored through $\hat c=\sqrt{k}\,\hat M$, the catalog uses
$$
k\mathcal D[\hat M]\rho=-\frac{k}{2}[\hat M,[\hat M,\rho]],\qquad
dJ=2\sqrt{k}\,\langle\hat M\rangle_c\,dt+dW,
$$
so $k$ is the coefficient carried by the monitored operator $\sqrt{k}\hat M$. Equivalently, the rescaled record $dy=dJ/(2\sqrt{k})$ obeys $dy=\langle\hat M\rangle_c\,dt+dW/(2\sqrt{k})$.

1. **Pure dephasing.**
   $$
   \dot\rho=\gamma\mathcal D[\hat\sigma_z]\rho,
   \qquad
   (\dot x,\dot y,\dot z)=(-2\gamma x,-2\gamma y,0).
   $$

2. **Driven, damped atom.**
   $$
   \dot\rho=-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho\right]
   +\gamma\mathcal D[\hat\sigma_-]\rho.
   $$

3. **Leaky cavity cat.**
   $$
   \dot\rho=\gamma\mathcal D[\hat a]\rho,
   \qquad
   \rho(0)=\frac{|\alpha\rangle\langle\alpha| + |\alpha\rangle\langle-\alpha| + |-\alpha\rangle\langle\alpha| + |-\alpha\rangle\langle-\alpha|}
   {2(1+e^{-2\left|\alpha\right|^2})}.
   $$

4. **Atom-made cat and its decoherence.** During the dispersive atom-cavity interaction and the subsequent storage interval, respectively,
   $$
   \dot\rho_{\mathrm{AC}}=-i[\chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\rho_{\mathrm{AC}}],
   \qquad
   \dot\rho_C=\gamma\mathcal D[\hat a]\rho_C.
   $$

5. **Dispersive qubit readout.**
   $$
   \dot\rho=-i[\epsilon(\hat a+\hat a^\dagger)
   +\chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\rho]
   +\gamma\mathcal D[\hat a]\rho,
   $$
   and, conditional on $\sigma_z=s=\pm1$, the coherent pointer obeys
   $$
   \dot\alpha_s=-i\epsilon-\left(is\chi+\frac{\gamma}{2}\right)\alpha_s.
   $$

6. **Quantum Brownian motion.** In the high-temperature Caldeira-Leggett limit,
   $$
   \dot\rho=-i\left[\frac{\hat p^2}{2M}+\frac12M\Omega^2\hat x^2,\rho\right]
   -i\gamma[\hat x,\{\hat p,\rho\}]
   -2\gamma Mk_BT[\hat x,[\hat x,\rho]].
   $$
   For the Gaussian simulation this closes as
   $$
   \dot{\bar r}=A\bar r,\qquad
   \dot\Sigma=A\Sigma+\Sigma A^T+D,\qquad
   A=\begin{pmatrix}0&1/M\\-M\Omega^2&-2\gamma\end{pmatrix},\quad
   D=\begin{pmatrix}0&0\\0&4\gamma Mk_BT\end{pmatrix}.
   $$
   The completely positive warm-oscillator comparison uses
   $$
   \dot\rho=\gamma(n_T+1)\mathcal D[\hat a]\rho
   +\gamma n_T\mathcal D[\hat a^\dagger]\rho.
   $$

7. **Photon counting.** With $\lambda_t=\mathrm{Tr}(\hat c^\dagger\hat c\rho_c)$ and $\mathbb E[dN_t|\rho_c]=\lambda_t\,dt$,
   $$
   d\rho_c=(-i[\hat H,\rho_c]+\mathcal D[\hat c]\rho_c)\,dt
   +\mathcal G[\hat c]\rho_c\,(dN_t-\lambda_t\,dt),
   $$
   where $\mathcal G[\hat c]\rho_c=\hat c\rho_c\hat c^\dagger/\lambda_t-\rho_c$.

8. **Homodyne detection.**
   $$
   d\rho_c=(-i[\hat H,\rho_c]+\mathcal D[\hat c]\rho_c)\,dt
   +\mathcal H[\hat c]\rho_c\,dW,
   \qquad
   dJ=\langle\hat c+\hat c^\dagger\rangle_c\,dt+dW.
   $$

9. **Heterodyne detection.** Put $\hat c_I=\hat c/\sqrt2$ and $\hat c_Q=i\hat c/\sqrt2$. Then
   $$
   d\rho_c=(-i[\hat H,\rho_c]+\mathcal D[\hat c]\rho_c)\,dt
   +\mathcal H[\hat c_I]\rho_c\,dW_I
   +\mathcal H[\hat c_Q]\rho_c\,dW_Q,
   $$
   with $dW_I^2=dW_Q^2=dt$ and $dW_IdW_Q=0$; the two records are
   $$
   dJ_j=\langle\hat c_j+\hat c_j^\dagger\rangle_c\,dt+dW_j,
   \qquad j\in\{I,Q\}.
   $$

10. **Quantum Zeno dynamics.** For the catalog's monitored operator $\sqrt{k}\,\hat\sigma_z$,
    $$
    d\rho_c=-i\frac{\Omega}{2}[\hat\sigma_x,\rho_c]\,dt
    +k\mathcal D[\hat\sigma_z]\rho_c\,dt
    +\sqrt{k}\,\mathcal H[\hat\sigma_z]\rho_c\,dW.
    $$
    Equivalently, the deterministic measurement term is $-\tfrac{k}{2}[\hat\sigma_z,[\hat\sigma_z,\rho_c]]\,dt$. This convention gives the Liouvillian pair $-k\pm i\sqrt{\Omega^2-k^2}$ and the exceptional point $k=\Omega$ used below.

11. **Charge qubit monitored by a QPC.** Since $\vec n_1=(\mathbb{1}-\hat\sigma_z)/2$,
    $$
    d\rho_c=-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho_c\right]\,dt
    +\frac{\kappa}{4}\mathcal D[\hat\sigma_z]\rho_c\,dt
    +\frac{\sqrt{\kappa}}{2}\mathcal H[\hat\sigma_z]\rho_c\,dW,
    \qquad
    dJ=\sqrt{\kappa}\,\langle\hat\sigma_z\rangle_c\,dt+dW.
    $$
    Dropping the innovation term gives the unconditional equation $\dot\rho=-i[(\Omega/2)\hat\sigma_x,\rho]+\kappa\mathcal D[\vec n_1]\rho$.

12. **Measurement-induced localization.**
    $$
    d\rho_c=k\mathcal D[\hat\sigma_z]\rho_c\,dt
    +\sqrt{k}\,\mathcal H[\hat\sigma_z]\rho_c\,dW,
    \qquad
    dy=\langle\hat\sigma_z\rangle_c\,dt+\frac{dW}{2\sqrt{k}}.
    $$

13. **Quantum Kalman filter.** With $\hat r=(\hat x,\hat p)^T$,
    $$
    \dot\Sigma=A\Sigma+\Sigma A^T+D-\Sigma C^TC\Sigma,
    \qquad
    d\bar r_c=A\bar r_c\,dt+\Sigma C^T(dJ-C\bar r_c\,dt),
    $$
    $$
    A=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
    D=\begin{pmatrix}0&0\\0&k\end{pmatrix},\qquad
    C=\begin{pmatrix}2\sqrt{k}&0\end{pmatrix},\qquad
    dJ=C\bar r_c\,dt+dW.
    $$

14. **Markovian measurement feedback.** For the feedback Hamiltonian $\hat H_{\rm fb}=\dot J\,\hat F$ with the current $\dot J=dJ/dt$,
    $$
    \dot\rho=-i[\hat H,\rho]+\mathcal D[\hat c]\rho
    -i[\hat F,\hat c\rho+\rho\hat c^\dagger]
    +\frac1\eta\mathcal D[\hat F]\rho.
    $$
    Equivalently,
    $$
    \dot\rho=-i\left[\hat H+\frac12(\hat c^\dagger\hat F+\hat F\hat c),\rho\right]
    +\mathcal D[\hat c-i\hat F]\rho+\frac{1-\eta}{\eta}\mathcal D[\hat F]\rho.
    $$
    The simulated ideal loop has $\hat c=\sqrt G\,\hat\sigma_z$, $\hat F=\sqrt G\,\hat\sigma_y$, and therefore $\dot\rho=G\mathcal D[\hat\sigma_z-i\hat\sigma_y]\rho$.

15. **Thermalization.**
    $$
    \dot\rho=-i[\omega\hat a^\dagger\hat a,\rho]
    +\gamma(n_T+1)\mathcal D[\hat a]\rho
    +\gamma n_T\mathcal D[\hat a^\dagger]\rho,
    \qquad
    n_T=(e^{\beta\omega}-1)^{-1}.
    $$

16. **Rapid purification by adaptive measurement.** For a monitored axis $\hat M$ at angle $\theta_t$ to the Bloch vector,
    $$
    d\rho_c=k\mathcal D[\hat M]\rho_c\,dt+\sqrt{k}\,\mathcal H[\hat M]\rho_c\,dW,
    $$
    and the impurity $\ell_t=(1-\left|\vec a_t\right|^2)/2$ drifts as
    $$
    \mathbb E[d\ell_t]=-4k\ell_t(\sin^2\theta_t+2\ell_t\cos^2\theta_t)\,dt.
    $$
    Ideal continuous-time feedback holds $\theta_t=\pi/2$, killing the noise and giving the pathwise deterministic $\dot\ell_t=-4k\ell_t$; the finite-step implementation below verifies convergence to it.

17. **Feedback cooling.** In this catalog's convention, with monitored operator $\sqrt{k}\,\hat x$,
    $$
    \begin{aligned}
    d\rho_c={}&-i[\omega\hat b^\dagger\hat b-f(t)\hat x,\rho_c]\,dt
    +k\mathcal D[\hat x]\rho_c\,dt+\sqrt{\eta k}\,\mathcal H[\hat x]\rho_c\,dW\\
    &+\gamma(n_T+1)\mathcal D[\hat b]\rho_c\,dt
    +\gamma n_T\mathcal D[\hat b^\dagger]\rho_c\,dt,
    \qquad f(t)=-G\langle\hat p\rangle_c.
    \end{aligned}
    $$

18. **Linear and nonlinear trajectory bookkeepings.** The unnormalized linear state obeys
    $$
    d\left|\tilde\psi\right\rangle=\left(-i\hat H-\frac12\sum_j\hat R_j^\dagger\hat R_j\right)\left|\tilde\psi\right\rangle dt
    +\sum_j\hat R_j\left|\tilde\psi\right\rangle dW_j.
    $$
    The normalized conditional state obeys
    $$
    d\rho_c=\mathcal L\rho_c\,dt+\sum_j\mathcal H[\hat R_j]\rho_c\,d\widehat W_j.
    $$

19. **Driven atom, averaged and watched.** The general unconditioned equation is
    $$
    \dot\rho=-i\left[\frac{\Delta\omega}{2}\hat\sigma_z+\frac{\Omega}{2}\hat\sigma_x,\rho\right]
    +\gamma(n_T+1)\mathcal D[\hat\sigma_-]\rho
    +\gamma n_T\mathcal D[\hat\sigma_+]\rho
    +\gamma k_d\mathcal D[\hat\sigma_z]\rho.
    $$
    In the resonant, zero-temperature, no-extra-dephasing corner watched by ideal homodyne detection,
    $$
    d\rho_c=\left(-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho_c\right]
    +\gamma\mathcal D[\hat\sigma_-]\rho_c\right)dt
    +\sqrt\gamma\,\mathcal H[\hat\sigma_-]\rho_c\,dW.
    $$

20. **Mollow spectrum.** The two-time correlations are generated by
    $$
    \dot\rho=-i\left[\frac{\Omega}{2}\hat\sigma_x,\rho\right]
    +\gamma\mathcal D[\hat\sigma_-]\rho,
    $$
    and quantum regression gives
    $$
    S_{\rm inel}(\mu)\propto\mathrm{Re}\int_0^\infty
    e^{-i\mu t}\,\langle\delta\hat\sigma_+(t)\delta\hat\sigma_-(0)\rangle_{\rm ss}\,dt.
    $$

## The Toolkit: Everything the Twenty Examples Share

Before any specific system, we build the small set of tools every example reuses. There are about ten of them, each a few lines. Once they are in hand, each of the twenty systems is just a short setup, a run, a picture, and a check.

### The density matrix, and quantities read from it

The smallest system is a qubit. Everything is built from the three Pauli matrices $X$, $Y$, $Z$ and the 2-by-2 identity:

```wl
{id2, X, Y, Z} = Table[PauliMatrix[j], {j, 0, 3}];
```

The basis states are the excited state $\left|e\right\rangle = \{1,0\}$ and the ground state $\left|g\right\rangle = \{0,1\}$; their equal superposition is $\left|+\right\rangle = (\left|e\right\rangle + \left|g\right\rangle)/\sqrt2$:

```wl
excited = {1, 0}; ground = {0, 1}; plus = {1, 1}/Sqrt[2];
```

The lowering operator $\hat\sigma_- = |g\rangle\langle e|$ takes the excited state to the ground state (an atom emitting a photon). Build it from the Pauli matrices:

```wl
lower = (X - I Y)/2;
```

Check that it sends $\left|e\right\rangle$ to $\left|g\right\rangle$:

```wl
lower . excited
```

It lands on the ground state. Its conjugate transpose is the raising operator $\hat\sigma_+$, written inline when needed.

A state vector describes only a pure state. The general state, pure or mixed, is the density matrix; for a pure state it is the outer product $\rho = |\psi\rangle\langle\psi|$:

```wl
densityMatrix[v_] := KroneckerProduct[v, Conjugate[v]];
```

Form the density matrix of the excited state:

```wl
densityMatrix[excited]
```

For $\left|e\right\rangle$ this is a single 1 in the top corner.

Two quantities summarize any qubit state. The first is the **Bloch vector** $(\langle\hat\sigma_x\rangle, \langle\hat\sigma_y\rangle, \langle\hat\sigma_z\rangle)$: it has length 1 for a pure state and shrinks as the state mixes. Define it:

```wl
blochVector[rho_] := Re[Tr[rho . #] & /@ {X, Y, Z}];
```

Check that $\left|+\right\rangle$ points along the $x$-axis:

```wl
blochVector[densityMatrix[plus]]
```

It points along the $x$-axis, length one.

The second quantity is the **purity** $\mathrm{Tr}(\rho^2)$: 1 for a pure state, one half for the maximally mixed qubit. Define it:

```wl
purity[rho_] := Re@Tr[rho . rho];
```

Read it for a pure state and for the maximally mixed one:

```wl
{purity[densityMatrix[plus]], purity[id2/2]}
```

The pure state gives one, the maximally mixed state one half.

Last, the **expectation value** $\langle A\rangle = \mathrm{Tr}(A\rho)$ of any operator. It is real when $A$ is Hermitian, as a physical meter observable must be, but it can be complex for a general operator:

```wl
expectation[op_, rho_] := Tr[op . rho];
```

### The two terms of the master equation

A state changes for two reasons. The first is unitary evolution under the Hamiltonian $\hat H$: the commutator term $-i[\hat H,\rho]$. Define it:

```wl
commutatorTerm[H_, rho_] := -I (H . rho - rho . H);
```

The second is a leak to the environment, one Lindblad (jump) operator $\hat c$ per channel (for a decaying atom, $\hat\sigma_-$). Its **dissipator** transfers population along the channel and damps coherences, $\mathcal{D}[\hat c]\rho = \hat c\rho\hat c^\dagger - \tfrac12(\hat c^\dagger\hat c\rho + \rho\hat c^\dagger\hat c)$:

```wl
dissipator[c_, rho_] := c . rho . ConjugateTranspose[c] -
   (ConjugateTranspose[c] . c . rho + rho . ConjugateTranspose[c] . c)/2;
```

A dissipator is trace-preserving, $\mathrm{Tr}(\mathcal{D}[\hat c]\rho) = 0$, so the total probability stays one. Confirm it for a general $\hat c$ and $\rho$:

```wl
FullSimplify[Tr[dissipator[{{c11, c12}, {c21, c22}}, {{r11, r12}, {r21, r22}}]]]
```

Zero with no assumptions.

Together the two terms give the **Lindblad master equation** $\dot\rho = \mathcal{L}\rho = -i[\hat H,\rho] + \sum_k\mathcal{D}[\hat c_k]\rho$, the evolution of an open system averaged over all measurement outcomes:

```wl
lindbladian[H_, leaks_List, rho_] := commutatorTerm[H, rho] + Total[dissipator[#, rho] & /@ leaks];
```

This is a system leaking to an environment whose record we set aside, whether no detector is present at all or a detector is present and we average over what it reports: smooth and deterministic, because averaging over every outcome removes the randomness. The first several examples are this equation with different $\hat H$ and leaks.

### The master equation, solved two independent ways

We solve the master equation two independent ways and check they agree. The first vectorizes $\rho$ (stack its rows) and applies the **Liouvillian**, the matrix whose action is $\mathcal{L}$. In this row-stacking convention (matching WL's `Flatten`; an overbar denotes entrywise conjugation) it is the Kronecker form
$$\mathcal{L} = -i(\hat H\otimes\mathbb{1} - \mathbb{1}\otimes\hat H^{T}) + \sum_k[\hat c_k\otimes\overline{\hat c_k} - \tfrac12(\hat c_k^\dagger\hat c_k\otimes\mathbb{1} + \mathbb{1}\otimes(\hat c_k^\dagger\hat c_k)^{T})],$$
which the code assembles directly, $\mathbb{1}$ the $d\times d$ identity (equivalently, the matrix whose columns are $\mathcal{L}$ applied to each basis matrix):

```wl
liouvillian[H_, leaks_List, d_] := With[{id = IdentityMatrix[d]},
   -I (KroneckerProduct[H, id] - KroneckerProduct[id, Transpose[H]]) +
    Total[Function[c, KroneckerProduct[c, Conjugate[c]] -
        (KroneckerProduct[ConjugateTranspose[c] . c, id] +
           KroneckerProduct[id, Transpose[ConjugateTranspose[c] . c]])/2] /@ leaks]];
```

The state at time $t$ is then a matrix exponential of the Liouvillian applied to the vectorized initial state, reshaped back to a matrix, $\mathrm{vec}\,\rho(t) = e^{\mathcal{L}t}\,\mathrm{vec}\,\rho_0$:

```wl
evolve[H_, leaks_List, rho0_, t_] := With[{d = Length[rho0]},
   ArrayReshape[MatrixExp[liouvillian[H, leaks, d] t] . Flatten[rho0], {d, d}]];
```

The Liouvillian also gives a **steady state** $\rho_{\mathrm{ss}}$, satisfying $\mathcal{L}\rho_{\mathrm{ss}} = 0$. A degenerate null space needs care: its arbitrary complex basis vectors are not themselves density matrices. `stationary` converts the complex null space into a complete real basis of Hermitian stationary operators, orthonormal in the Frobenius sense. When that basis is one-dimensional (the usual case, a single steady direction) it rescales the sole operator to a unit-trace Hermitian density matrix, so `First@stationary[...]` is the physical steady state. When the null space is degenerate, or the lone direction is traceless (a steady coherence, not a state), it returns the operator basis unchanged; the physical steady states are then the unit-trace, positive combinations of that basis, not the basis vectors separately.

```wl
stationary[H_, leaks_List, tol_: 10^-10] :=
  With[{d = Length[H]},
   With[{candidates =
       With[{raw = ArrayReshape[#, {d, d}] & /@ NullSpace[liouvillian[H, leaks, d]]},
        Join[(# + ConjugateTranspose[#])/2 & /@ raw,
         (# - ConjugateTranspose[#])/(2 I) & /@ raw]]},
    With[{basis = ArrayReshape[#, {d, d}] & /@
        Select[Orthogonalize[Flatten /@ candidates], Norm[N@#] > tol &]},
     If[Length[basis] =!= 1, basis,
      With[{weight = Tr[First@basis]},
       If[Abs[N@weight] <= tol, basis,
        With[{rho = First@basis/weight}, {(rho + ConjugateTranspose[rho])/2}]]]]]]];
```

The second way integrates the master equation directly with an ODE solver, sharing no code with the Liouvillian route:

```wl
evolveODE[H_, leaks_List, rho0_, t1_] :=
  NDSolveValue[{s'[t] == lindbladian[H, leaks, s[t]], s[0] == N@rho0}, s, {t, 0, t1}];
```

Now the cross-check. Fix a test system: a qubit driven by $X$, dephased through $Z$, starting in the ground state:

```wl
Hx = 1.0 X; oneLeak = {Sqrt[0.4] Z}; start = densityMatrix[ground];
```

Integrate its master equation once with the ODE route:

```wl
solved = evolveODE[Hx, oneLeak, start, 10.0];
```

Measure the largest disagreement between the two solvers over the whole span:

```wl
Max@Table[Norm[evolve[Hx, oneLeak, start, tt] - solved[tt], "Frobenius"], {tt, 0, 10, 0.5}]
```

The largest disagreement is at the ODE solver's tolerance: the matrix exponential and the ODE march give the same state. Either is a trustworthy reference.

### The positivity-preserving measurement step

Now put a detector on a leak and condition the state on its output. The naive update, nudging $\rho$ by the record, can push $\rho$ off the set of valid states and give a negative eigenvalue, a negative probability. See it break: take $\left|+\right\rangle$, measure $Z$ with no Hamiltonian, and feed one large noise kick into the naive update:

```wl
naiveStep = densityMatrix[plus] + dissipator[Z, densityMatrix[plus]] 0.1 +
   (Z . densityMatrix[plus] + densityMatrix[plus] . Z -
      Tr[(2 Z) . densityMatrix[plus]] densityMatrix[plus]) 0.5;
```

Read the smallest eigenvalue of the updated state:

```wl
Min@Eigenvalues[naiveStep]
```

A negative eigenvalue: the naive update has left the set of valid states. Shrinking $dt$ makes it rarer but never impossible.

The fix is the structure-preserving filter of [Rouchon and Ralph (2015)](https://arxiv.org/abs/1410.5345): write the update as a sum of Kraus terms $A\rho A^\dagger$, an operator times $\rho$ times its conjugate transpose, so the unnormalized result is manifestly positive whenever the supplied efficiencies satisfy $0\leq\eta_k\leq1$. With monitored operators $\{\hat c_k\}$ at detection efficiencies $\{\eta_k\}$, unmonitored leaks $\{\hat l_j\}$, record increments $dJ_k$, the no-jump generator $K = i\hat H + \tfrac12(\sum_k \hat c_k^\dagger \hat c_k + \sum_j \hat l_j^\dagger \hat l_j)$, and the measurement kick $s = \sum_k \sqrt{\eta_k}\,dJ_k\,\hat c_k$, the corrected Kraus operator is
$$M = \mathbb{1} - K\,dt + s + \tfrac12 s^2 - \tfrac{dt}{2}\sum_k \eta_k\,\hat c_k^2,$$
and the stepped state is the renormalized sum
$$\tilde\rho = M\rho M^\dagger + dt\sum_j \hat l_j\,\rho\,\hat l_j^\dagger + dt\sum_k (1-\eta_k)\,\hat c_k\,\rho\,\hat c_k^\dagger, \qquad \rho' = \frac{\tilde\rho}{\mathrm{Tr}\,\tilde\rho}.$$
The three terms of $\tilde\rho$ are the measured record, the unmonitored leaks, and the unrecorded fraction $(1-\eta_k)$ of each watched channel. Each has Kraus form, so a positive input gives a positive numerator; the normalized update is defined when that numerator has positive trace. This structural positivity is distinct from time-discretization accuracy: at finite $dt$ the update is an approximation to the continuous stochastic master equation, and its ensemble average has a systematic step-size bias in addition to Monte-Carlo scatter. Only one precondition matters physically: every detection efficiency must lie in $[0,1]$, since the upper bound $\eta_k\le1$ is exactly what keeps the Kraus-form update positive. The code asserts that and a positive time step once, up front; the returned stepper then preserves Hermiticity, positivity, and trace by construction, so it trusts its input and guards only the final normalization against a non-positive trace:

```wl
measurementStep[H_, watched_, effs_, unwatched_, dt_] := Enclose[
  Module[{d = Length[H], id = IdentityMatrix[Length[H]], channels = Join[watched, unwatched], drift, corr},
   ConfirmAssert[VectorQ[effs, 0 <= # <= 1 &], "efficiencies must lie in [0,1]"];
   ConfirmAssert[TrueQ[dt > 0], "time step must be positive"];
   drift = I H + Fold[Plus, 0 H, ConjugateTranspose[#] . # & /@ channels]/2;
   corr = Fold[Plus, 0 H, MapThread[#1 (#2 . #2) &, {effs, watched}]]/2;
   Function[{rho, dJ}, Module[{sig, M, top, nrm},
     sig = Fold[Plus, 0 H, MapThread[Sqrt[#2] #1 #3 &, {dJ, effs, watched}]];
     M = id - drift dt + sig + sig . sig/2 - corr dt;
     top = M . rho . ConjugateTranspose[M] +
       dt Fold[Plus, 0 H, # . rho . ConjugateTranspose[#] & /@ unwatched] +
       dt Fold[Plus, 0 H, MapThread[(1 - #2) #1 . rho . ConjugateTranspose[#1] &, {watched, effs}]];
     nrm = Re@Tr[top];
     If[TrueQ[nrm > 0], top/nrm,
       Failure["NonPositiveNormalization", <|"Normalization" -> nrm|>]]]]]];
```

Here `watched` are the monitored channels, `effs` their detection efficiencies (1 for perfect, less if the detector misses some output), and `unwatched` the undetected leaks. Revisit the kick that broke the naive step, now through the positivity-preserving step:

```wl
fixedStep = measurementStep[0 id2, {Z}, {1.}, {}, 0.1][densityMatrix[plus], {0.5}];
```

Read the smallest eigenvalue and the trace of the stepped state:

```wl
{Min@Re@Eigenvalues[fixedStep], Re@Tr[fixedStep]}
```

The smallest eigenvalue is now zero to numerical precision and the trace is one: the step lands on the boundary of the state space instead of leaving it. Positivity of the normalized Kraus-form update no longer depends on taking an infinitesimal step, although agreement with the continuous-time filter still does.

### The measurement record, and one trajectory

The homodyne record at each step is the signal $\langle\hat c + \hat c^\dagger\rangle\,dt$ plus fresh Gaussian noise. Write it:

```wl
measurementRecord[rho_, watched_List, effs_List, dt_, kick_List] :=
  MapThread[Sqrt[#3] Re@Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &,
   {watched, kick, effs}];
```

`kick` is the Wiener increment $dW$, drawn fresh each step. One **quantum trajectory** is: draw $dW$, form the record, advance with the positivity step. This is the only place randomness enters. It returns the times, the states, and the record:

```wl
trajectory[rho0_, H_, watched_List, effs_List, unwatched_List, dt_, tf_, seed_] :=
  BlockRandom[SeedRandom[seed];
   Module[{n = Round[tf/dt], step, kicks, states, record},
    step = measurementStep[H, watched, effs, unwatched, dt];
    kicks = RandomVariate[NormalDistribution[0, Sqrt[dt]], {n, Length[watched]}];
    states = FoldList[
      Function[{r, dw}, step[r, measurementRecord[r, watched, effs, dt, dw]]], rho0, kicks];
    record = MapThread[measurementRecord[#1, watched, effs, dt, #2] &, {Most[states], kicks}];
    <|"times" -> dt Range[0, n], "states" -> states, "record" -> record|>]];
```

### The check that ties the catalog together: trajectories average to the master equation

The central fact, used as a check almost everywhere: averaging many trajectories recovers the master-equation solution. Averaging over the record is what tracing out the detector means. Fix the step and the horizon:

```wl
dt = 0.005; tf = 8.0;
```

Run 200 trajectories of the test qubit:

```wl
crowd = Table[
   trajectory[start, Hx, oneLeak, {1.}, {}, dt, tf, k]["states"], {k, 1, 200}];
```

Integrate the master equation as the reference:

```wl
smoothRef = evolveODE[Hx, oneLeak, start, tf];
```

Measure the largest gap between the ensemble mean and that reference:

```wl
Max@MapThread[Norm[#1 - #2, "Frobenius"] &,
   {Mean[crowd], smoothRef /@ (dt Range[0, Round[tf/dt]])}]
```

The gap contains two errors: Monte-Carlo scatter, which shrinks as $1/\sqrt N$, and the finite-$dt$ bias of the discrete trajectory update, which shrinks only when the time step is refined. The trajectory ensemble and the master equation share no evolution code; agreement under both ensemble enlargement and time-step refinement is the numerical validation repeated throughout.

Two independent errors keep the ensemble mean from the master equation: too few trajectories, which leaves sampling scatter, and too coarse a step, which leaves discretization bias. `gapDiagnosis` separates them in a single run. It averages an ensemble built from the same `trajectory` used above, measures the gap between the mean and the master equation, and reads the ensemble's own scatter for how large a gap sampling alone would leave at this many trajectories; gap and scatter add in quadrature, so $\sqrt{\text{gap}^2 - \text{scatter}^2}$ is the systematic remainder. Beside it sits the per-step bias computed with no sampling at all: a quadrature over one full Gaussian record increment that averages one discrete `measurementStep` against the exact Lindblad step, one integration axis per watched channel, cheap for one channel or the two of a heterodyne pair and impractical past three. The gap is read at one instant, the final time by default, a chosen time, or over the whole grid with `All`:

```wl
gapDiagnosis[rho0_, H_, watched_List, effs_List, unwatched_List, dt_, tf_, nTraj_, atTime_ : Automatic] :=
  Module[{steps = Round[tf/dt], leaks = Join[watched, unwatched], d = Length[rho0],
     grid, runs, means, gaps, it, gap, sem, sysGap, bias, stepper, entry},
   stepper = measurementStep[H, watched, effs, unwatched, dt];
   entry[dw_?(VectorQ[#, NumericQ] &), j_, k_] :=
     (Times @@ (PDF[NormalDistribution[0, Sqrt[dt]], #] & /@ dw))
       stepper[rho0, measurementRecord[rho0, watched, effs, dt, dw]][[j, k]];
   bias = With[{vars = Table[Unique["dw"], {Length[watched]}]},
     Norm[Table[NIntegrate[entry[vars, j, k],
        Evaluate[Sequence @@ ({#, -Infinity, Infinity} & /@ vars)], Method -> "GlobalAdaptive"],
       {j, d}, {k, d}] - evolve[H, leaks, rho0, dt], "Frobenius"]];
   runs = Table[trajectory[rho0, H, watched, effs, unwatched, dt, tf, k]["states"], {k, nTraj}];
   grid = dt Range[0, steps];
   means = Mean[runs];
   gaps = MapThread[Norm[#1 - #2, "Frobenius"] &, {means, evolve[H, leaks, rho0, #] & /@ grid}];
   it = Which[atTime === All, First@Ordering[gaps, -1], atTime === Automatic, steps + 1,
      True, Clip[Round[atTime/dt] + 1, {1, steps + 1}]];
   gap = gaps[[it]];
   sem = Sqrt[Total[Variance[runs[[All, it]]], 2]/nTraj];
   sysGap = Sqrt[Max[gap^2 - sem^2, 0]];
   <|"Time" -> grid[[it]], "EnsembleGap" -> gap, "SamplingError" -> sem,
     "SystematicGap" -> sysGap, "StepBias" -> bias,
     "Verdict" -> If[gap <= 2 sem, "within sampling noise: add trajectories to tighten",
        "discretization bias dominates: refine dt"]|>];
```

At a fine step the gap sits inside the sampling scatter, the systematic remainder is zero, and only more trajectories would tighten it:

```wl
gapDiagnosis[start, Hx, oneLeak, {1.}, {}, 0.005, 8.0, 200]
```

At a coarse step the same run splits the other way: the gap outgrows the scatter, a systematic remainder survives, and the per-step bias, carrying no sampling, points to the step itself as the cause:

```wl
gapDiagnosis[start, Hx, oneLeak, {1.}, {}, 0.1, 8.0, 2000]
```

The first run is sampling-limited, the second step-limited. These are the two knobs behind every convergence check that follows: raise the trajectory count until the scatter is small, refine the step until the systematic remainder is small, and the ensemble mean meets the master equation. The homodyne and heterodyne examples apply `gapDiagnosis` directly to their diffusive ensembles; the jump unravellings, which this diffusive diagnostic does not cover, run the same two-knob check in their own idiom. The surviving gap is read from the ensemble rather than from the per-step bias times the step count, because the master equation is contractive: each step's error decays over the steps that follow, so what remains at a given time is smaller than the running sum of per-step biases.

One more fact from the 200-trajectory ensemble: each **conditional** trajectory stays pure, while the ensemble average is mixed. Compare the purity of one trajectory's final state with that of the ensemble mean:

```wl
{purity[Last@trajectory[start, Hx, oneLeak, {1.}, {}, dt, tf, 1]["states"]],
 purity[Last@Mean[crowd]]}
```

The single trajectory ends pure; the average ends mixed. This is the heart of the subject: under ideal monitoring (a pure initial state, every output channel watched, unit detector efficiency, no extra unrecorded noise) decoherence does not happen to any one conditioned state, which stays pure, and mixing is entirely the result of averaging over the unrecorded outcomes. Relax any of those conditions, add inefficiency, an unmonitored channel, thermal drive, or dephasing, and the conditioned state itself becomes mixed; the pure-trajectory picture is the ideal limit, and where the examples below invoke it they are in that limit.

### The Bloch-sphere plot

For a qubit the natural picture is the Bloch sphere: the Bloch vector's tip lies on the surface for a pure state, inside for a mixed one. Draw it with the `BlochSpherePlot` resource function, whose empty call returns the shaded ball with its axes and labels as a `Graphics3D` ready to overlay each trajectory on:

```wl
bloch[] := ResourceFunction["BlochSpherePlot"][]
```

A small octahedron will mark where a trajectory ends; build it from its six vertices:

```wl
diamond[p_, r_] := GraphicsComplex[p + # & /@ (r {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}),
   Polygon[{{1, 3, 5}, {3, 2, 5}, {2, 4, 5}, {4, 1, 5}, {3, 1, 6}, {2, 3, 6}, {4, 2, 6}, {1, 4, 6}}]];
```

The plotter uses one visual grammar for every trajectory: the path is a curve, a small ball marks where it starts, a diamond marks where it ends, and all three share one color per trajectory. A label names the picture, a legend names the trajectories when they mean different things, and extra marks (fixed points, targets) ride along as plain graphics primitives:

```wl
blochPlot[paths_, lbl_ : None, names_ : None, extras_ : {}] := With[{g = Show[bloch[],
     Graphics3D[MapIndexed[{ColorData[97, First@#2], Thick, Line[#1],
         Sphere[First[#1], 0.035], diamond[Last[#1], 0.06]} &, paths]],
     Graphics3D[extras], PlotLabel -> lbl, ImageSize -> 340]},
   If[names === None, g, Legended[g, LineLegend[ColorData[97, #] & /@ Range[Length[names]], names]]]];
```

Run one monitored trajectory of the test qubit:

```wl
demoRun = trajectory[{{1/2, 1/2}, {1/2, 1/2}}, Hx, oneLeak, {1.}, {}, 0.01, 3.0, 4];
```

Now visualize the monitored trajectory on the Bloch sphere, from its starting ball to its final diamond:

```wl
blochPlot[{blochVector /@ demoRun["states"]}, "one monitored trajectory, ball to diamond"]
```

The Bloch vector wanders as the record comes in; one seed is one run. Change the seed and it wanders differently; averaging many recovers the smooth relaxation checked above.

The sphere shows the shape of the path; the same run also reads as three numbers in time, one curve for each Bloch component:

```wl
ListLinePlot[Transpose[blochVector /@ demoRun["states"]],
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2], ColorData[97, 3]},
 PlotLegends -> {"\[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]", "\[LeftAngleBracket]\[Sigma]y\[RightAngleBracket]", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05}, ImageSize -> 520,
 FrameLabel -> {"time", "Bloch component"}, DataRange -> MinMax@demoRun["times"],
 PlotLabel -> "the same run, component by component"]
```

All three components move. $\langle\hat\sigma_x\rangle$ starts at one and decays toward zero, since the drive turns the state about the $x$-axis and cannot rebuild what watching $\hat\sigma_z$ erases; meanwhile $\langle\hat\sigma_y\rangle$ and $\langle\hat\sigma_z\rangle$ wander as that same drive rotates them and the record tips the state toward a pole. The sphere and the curves are one trajectory seen two ways: the sphere gives the shape of the path, the curves give each direction changing in time. This component-by-component reading is how a watched qubit is usually plotted, and most later examples show just the one component the measurement makes interesting, $\langle\hat\sigma_z\rangle$ for a $\hat\sigma_z$ measurement. With the qubit toolkit in hand, on to the first systems.

## Part One: Systems We Can Solve Symbolically

The first two systems have exact closed-form answers we can check against the code. Both are qubits, so first a convenience: read the master equation as three **Bloch equations** for $(\dot x, \dot y, \dot z)$. Write the general Bloch state $\rho = \tfrac12(\mathbb{1} + x X + y Y + z Z)$:

```wl
blochState[x_, y_, z_] := (id2 + x X + y Y + z Z)/2;
```

Feed it to the Lindbladian and extract each component's rate:

```wl
rates[H_, leaks_List] := ComplexExpand@Table[
   Tr[d . lindbladian[H, leaks, blochState[x, y, z]]], {d, {X, Y, Z}}];
```

### Pure Dephasing

**The problem.** A two-level system dephased along $\hat\sigma_z$, with no Hamiltonian:
$$\dot\rho = \gamma\,\mathcal{D}[\hat\sigma_z]\rho.$$
Since $\hat\sigma_z^\dagger\hat\sigma_z = \mathbb{1}$ this is pure dephasing, so in Bloch components $\dot x = -2\gamma x$, $\dot y = -2\gamma y$, $\dot z = 0$, giving $x(t) = x_0 e^{-2\gamma t}$, $y(t) = y_0 e^{-2\gamma t}$, $z(t) = z_0$. Both poles $\left|\hat\sigma_z = \pm 1\right\rangle$ (the $\hat\sigma_z$ eigenstates $\left|e\right\rangle$ and $\left|g\right\rangle$, not to be confused with the $\hat\sigma_x$ eigenstate the code calls `plus`) are fixed points, so there is no unique steady state: any mixture $\rho_{\mathrm{ss}} = p_+|e\rangle\langle e| + p_-|g\rangle\langle g|$ sits still forever.

Pure dephasing continuously measures $\hat\sigma_z$: it damps the off-diagonal coherences $\rho_{eg} = \rho_{ge}^*$ (the $x$ and $y$ Bloch components, perpendicular to the measured $z$ axis) and leaves the diagonal populations, set by $z$, unchanged. Turning off the Hamiltonian and with one dissipator, then the Bloch equations read as:

```wl
Clear[\[Gamma]]; Simplify[rates[0 id2, {Sqrt[\[Gamma]] Z}], \[Gamma] > 0]
```

The $x$ and $y$ components decay at rate $2\gamma$; $z$ is fixed. Hand those rates, exactly as the Lindbladian produced them, to the differential-equation solver as the Bloch equations, and solve from any initial Bloch vector:

```wl
Clear[x, y, z, t, x0, y0, z0];
closedDephasing[t_] = {x[t], y[t], z[t]} /. First@DSolve[Join[
     Thread[{x'[t], y'[t], z'[t]} ==
       (Simplify[rates[0 id2, {Sqrt[\[Gamma]] Z}], \[Gamma] > 0] /. {x -> x[t], y -> y[t], z -> z[t]})],
     {x[0] == x0, y[0] == y0, z[0] == z0}], {x, y, z}, t]
```

The $x$ and $y$ components decay exponentially, $z$ is frozen. Check the closed form against numerical results from the master equation over time:

```wl
With[{gg = 0.5, a0 = {0.8, 0.4, 0.3}},
 Max@Table[Norm[(closedDephasing[tt] /. {\[Gamma] -> gg, x0 -> a0[[1]], y0 -> a0[[2]], z0 -> a0[[3]]}) -
    blochVector[evolve[0 id2, {Sqrt[gg] Z}, blochState @@ a0, tt]]], {tt, 0, 3, 0.3}]]
```

The maximum discrepancy is at numerical solver precision.

Because $z$ never changes, there is no unique steady state: any mixture of $\left|e\right\rangle$ and $\left|g\right\rangle$ is stationary. The Hermitian stationary-operator basis is two-dimensional:

```wl
Chop@stationary[0 id2, {Sqrt[0.5] Z}]
```

The returned basis spans the identity and $Z$ directions. It is an operator basis, not a list of physical states: every physical member of the stationary family is $\rho_{\mathrm{ss}}(z)=\tfrac12(\mathbb{1}+zZ)$ with $-1\leq z\leq1$. Verify the whole family directly:

```wl
FullSimplify[lindbladian[0 id2, {Sqrt[1/2] Z}, blochState[0, 0, z]]]
```

It vanishes identically. The endpoints $z=\pm1$ are the two poles, and every convex mixture between them is stationary.

Now visualize the dephasing behavior above: start two Bloch vectors off-axis and watch each slide straight toward the $z$-axis (ball to diamond), keeping its height:

```wl
With[{gg = 0.6, starts = {{0.85, 0.25, 0.4}, {-0.6, -0.55, -0.35}}},
 blochPlot[
  Table[Table[closedDephasing[tt] /. {\[Gamma] -> gg, x0 -> s[[1]], y0 -> s[[2]], z0 -> s[[3]]}, {tt, 0, 2.5, 0.05}], {s, starts}],
  "dephasing slides each state to the z-axis", {"start (0.85, 0.25, 0.4)", "start (-0.6, -0.55, -0.35)"}]]
```

Each tip drifts to the $z$-axis and stops at its own height: the off-diagonal coherence is gone, the populations untouched. This is dephasing in its purest form, decoherence with nothing else attached.

Two quantities summarize the dephasing. The $x$ and $y$ components decay at rate $2\gamma$, so the coherence time is $T_2 = 1/(2\gamma)$, while the populations never move, so there is no energy relaxation, $T_1 = \infty$. The purity follows from the closed-form Bloch vector, since a state with Bloch vector $\vec r$ has purity $\tfrac12(1 + \left|\vec r\right|^2)$:

```wl
FullSimplify[(1 + closedDephasing[t] . closedDephasing[t])/2, \[Gamma] > 0]
```

The purity is $\tfrac12(1 + z_0^2 + e^{-4\gamma t}(x_0^2 + y_0^2))$: it relaxes from its start to $\tfrac12(1 + z_0^2)$, the purity of the leftover diagonal state, at twice the coherence rate. The environment has learned which $\hat\sigma_z$ eigenstate the atom is in and taken none of its energy.

### A Driven, Damped Atom Relaxing to Steady State

**The problem.** A two-level atom driven by a laser (a drive about $\hat\sigma_x$) while it decays to its ground state:
$$\mathcal{L}\rho = -i\,\frac{\Omega}{2}[\hat\sigma_x,\rho] + \gamma\,\mathcal{D}[\hat\sigma_-]\rho.$$
In Bloch components the drift is affine, $\dot{\vec r} = A\,\vec r + \vec b$ with $\vec r = (x,y,z)$,
$$\frac{d}{dt}\begin{pmatrix}x\\y\\z\end{pmatrix} = \begin{pmatrix}-\gamma/2 & 0 & 0\\ 0 & -\gamma/2 & -\Omega\\ 0 & \Omega & -\gamma\end{pmatrix}\begin{pmatrix}x\\y\\z\end{pmatrix} + \begin{pmatrix}0\\0\\-\gamma\end{pmatrix},$$
which spirals in to the steady state
$$x_{\mathrm{ss}} = 0,\qquad y_{\mathrm{ss}} = \frac{2\Omega\gamma}{\gamma^2 + 2\Omega^2},\qquad z_{\mathrm{ss}} = -\frac{\gamma^2}{\gamma^2 + 2\Omega^2},$$
at the rate set by the eigenvalues $\lambda_\pm = -\frac{3\gamma}{4} \pm i\tilde\Omega$, with $\tilde\Omega = \sqrt{\Omega^2 - (\gamma/4)^2}$ (a damped oscillation once $\Omega > \gamma/4$). Solve it and watch the spiral.

Given a laser drive about $\hat\sigma_x$ and spontaneous emission through $\hat\sigma_-$, read the Bloch equations:

```wl
Clear[\[Gamma], \[CapitalOmega]]; Simplify[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}], \[Gamma] > 0]
```

The component $x$ decays; the drive couples $y$ and $z$ while emission pulls $z$ toward the ground state. This has a unique steady state. Find it by setting all three rates to zero:

```wl
Clear[x, y, z];
steady = First@Solve[Thread[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}] == 0], {x, y, z}];
```

Read the steady state in closed form:

```wl
Simplify[{x, y, z} /. steady, \[Gamma] > 0]
```

The steady state in closed form: the $x$-component vanishes, so the Bloch vector lies in the $y$-$z$ plane, pulled up out of the ground state ($z=-1$) toward the equator by an amount set by the drive-to-emission ratio $\Omega/\gamma$.

Two physical numbers fall straight out of this steady state. The first is the excited-state population $\rho_{\mathrm{ee}}^{\mathrm{ss}} = \tfrac12(1 + z_{\mathrm{ss}}) = \Omega^2/(\gamma^2 + 2\Omega^2)$: it grows as $\Omega^2/\gamma^2$ when the drive is weak and saturates at one half when the drive is strong, so the excited population never passes one half, and no drive can invert the atom (make it more likely excited than in its ground state). The second is the steady rate at which the atom scatters photons, the fluorescence flux $\gamma\,\rho_{\mathrm{ee}}^{\mathrm{ss}}$: the spontaneous-emission rate $\gamma$ times that population. Read both off the closed form:

```wl
Grid[
 Transpose@{
   {"the excited-state population", "the fluorescence flux"},
   FullSimplify /@ With[{zss = Simplify[z /. steady, \[Gamma] > 0]}, {(1 + zss)/2, \[Gamma] (1 + zss)/2}]},
 Alignment -> Left, Frame -> All]
```

How does the atom approach the stationary state? Consider the affine drift $\dot{\vec r} = A\,\vec r + \vec b$ and subtract the steady state. Since $A\,\vec r_{\mathrm{ss}} + \vec b = 0$, the constant term cancels and the displacement $\delta\vec r = \vec r - \vec r_{\mathrm{ss}}$ obeys a purely linear law,
$$\frac{d}{dt}\,\delta\vec r = A\,\delta\vec r,$$
so the motion is a sum of eigenmodes of $A$, each of the form $e^{\lambda t}\vec v$: the real part of $\lambda$ shrinks a mode, the imaginary part rotates it. The eigenvalues of $A$ (the Jacobian of the rates) therefore set the entire approach:

```wl
Simplify[Eigenvalues[D[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}], {{x, y, z}}]], \[Gamma] > 0]
```

Two of the three are a complex-conjugate pair, so that mode is a **damped spiral**, a ring-down whenever the drive exceeds a quarter of the emission rate, $\Omega > \gamma/4$; the third is real and belongs to the $x$-direction, a mode that only shrinks. The eigenvalues fix these rates, but which of them a given trajectory excites depends on where it starts, which is the opening for the straight-line approach just below.

But a special pair of initial states skips the spiral entirely and slides straight in, and the reason sits in the drift matrix itself. Check that $(1,0,0)$ is one of its eigenvectors, with eigenvalue $-\gamma/2$:

```wl
Simplify[D[rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}], {{x, y, z}}] . {1, 0, 0}, \[Gamma] > 0]
```

The result is $-\tfrac{\gamma}{2}(1,0,0)$: a state displaced from the steady state purely along $x$ stays purely $x$-displaced as it evolves (the drift never turns that displacement into $y$ or $z$) and simply shrinks at rate $\gamma/2$. So the line through the steady state parallel to the $x$-axis is invariant: a state starting on it slides straight down it into the steady state, with no spiral. That line meets the Bloch sphere at two pure states, and since $x_{\mathrm{ss}} = 0$ they sit at $x = \pm\sqrt{1 - y_{\mathrm{ss}}^2 - z_{\mathrm{ss}}^2}$ with the steady $y_{\mathrm{ss}}, z_{\mathrm{ss}}$:

```wl
flatPair = Simplify[With[{ss = {x, y, z} /. steady},
    {# Sqrt[1 - ss[[2]]^2 - ss[[3]]^2], ss[[2]], ss[[3]]} & /@ {1, -1}], \[CapitalOmega] > 0 && \[Gamma] > 0]
```

These two states are one way to write the mixed steady state as a blend of two pure states, and a mixed state has infinitely many such blends. This one is special for a reason about measurement. If the atom's emitted light is continuously monitored, its conditioned state stays pure and wanders as the record arrives, and averaging over the record returns $\rho_{\mathrm{ss}}$; for certain blends a real detector can keep that wandering confined to the blend's own states, and such a blend is a **physically realizable ensemble** (Wiseman and Vaccaro). This two-state one is the simplest kind, a **physically realizable (PR) pair**. Not every blend qualifies, and the averaged, unmonitored drift $\dot{\vec r}$ already tells the two apart. Read it at a point off the same `rates`, kept symbolic in $\Omega,\gamma$:

```wl
driftAt[r_] := rates[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}] /. Thread[{x, y, z} -> r]
```

The PR pair is `flatPair`, the two states where the invariant line pierces the sphere. Its rival is the *orthogonal* ensemble, the steady-state density matrix's own two eigenstates, the antipodal pure states in the $\pm\vec r_{\mathrm{ss}}$ directions:

```wl
orthPair = Simplify[{#, -#} &@Normalize[{x, y, z} /. steady], \[CapitalOmega] > 0 && \[Gamma] > 0]
```

With both pairs in hand, watch how they sit on the sphere as the drive-to-emission ratio varies. Drag $\Omega/\gamma$ (at $\gamma = 1$): the steady state is green, the PR pair blue joined by its dashed chord, the orthogonal pair orange. The PR chord stays parallel to the $x$-axis at every setting, while the orthogonal chord runs through the origin along the steady direction, the geometric split the drifts below make precise.

```wl
Manipulate[
 With[{ss = {x, y, z} /. steady /. {\[Gamma] -> 1, \[CapitalOmega] -> OO},
   flat = flatPair /. {\[Gamma] -> 1, \[CapitalOmega] -> OO},
   orth = orthPair /. {\[Gamma] -> 1, \[CapitalOmega] -> OO}},
  Legended[
   Show[bloch[],
    Graphics3D[{{PointSize[.02], Green, Point@ss},
      {PointSize[.02], ColorData[97][1], Point@flat},
      {PointSize[.02], ColorData[97][2], Point@orth},
      {Dashed, Line[flat], Line[orth]}}]],
   SwatchLegend[{Green, ColorData[97][1], ColorData[97][2]}, {"Steady state", "PR pair", "Ortho pair"}]]],
 {{OO, 1, "\[CapitalOmega]/\[Gamma]"}, 0, 10, AppearanceElements -> All}]
```

The drift at a PR state, in closed form:

```wl
Simplify[driftAt[flatPair[[1]]], \[CapitalOmega] > 0 && \[Gamma] > 0]
```

and at an eigenstate:

```wl
Simplify[driftAt[orthPair[[1]]], \[CapitalOmega] > 0 && \[Gamma] > 0]
```

For every $\Omega$ and $\gamma$ the PR drift is purely along $x$, exactly along the segment joining the two PR states (they differ only in $x$), while the eigenstate drift is purely along $z$, which cuts across the segment joining the two eigenstates (tilted in the $y$-$z$ plane). The PR pair sits on a line the averaged flow leaves invariant, the eigenstate pair does not, and that is the intuition for why one blend is physically realizable and the other is not. The drift is only that intuition, though: it describes the averaged, unmonitored dynamics, and it neither constructs the measurement that realizes the PR pair nor shows the conditional state hopping between its two ends. That waits for the trajectory examples below, where the stochastic machinery is actually run.

Now visualize all five paths: the spiral from a generic pure state, the two flat slides from the PR pair, and the two spirals from the orthogonal pair. This is the *unconditional* picture, the record averaged away, so every start relaxes into the mixed steady state:

```wl
With[{gg = 1., om = 3.},
 Module[{startsFlat, startsOrtho, path},
  startsFlat = N[flatPair /. {\[Gamma] -> gg, \[CapitalOmega] -> om}];
  startsOrtho = N[orthPair /. {\[Gamma] -> gg, \[CapitalOmega] -> om}];
  path[r0_] := With[{run = evolveODE[(om/2) X, {Sqrt[gg] lower}, r0, 14.]},
    blochVector[run[#]] & /@ Range[0, 14, 0.02]];
  blochPlot[{path[densityMatrix[{Cos[\[Pi]/8], E^(I \[Pi]/4) Sin[\[Pi]/8]}]], path[blochState @@ startsFlat[[1]]],
    path[blochState @@ startsFlat[[2]]], path[blochState @@ startsOrtho[[1]]],
    path[blochState @@ startsOrtho[[2]]]},
   "every start lands on the steady state",
   {"from a generic pure state", "physically realizable start +", "physically realizable start -",
    "orthogonal (eigenstate) start +", "orthogonal (eigenstate) start -"}]]]
```

The two physically realizable starts slide straight in along the invariant line, while the generic start and the two orthogonal starts, all off that line, ring down in a damped spiral. This is the drift geometry of the previous cells made visible: the averaged flow keeps the PR pair on its own chord, so those starts slide, and carries the orthogonal pair off theirs, so those must spiral. All five diamonds land on the same Bloch vector inside the ball, a driven, damped atom's *unconditional* steady state: one mixed point set by the balance of drive and emission, whichever pure state you start from.

## Part Two: Light in a Box, and States That Live in Two Places

The next systems are a single mode of light in a cavity. Its states span the Fock ladder $\left|0\right\rangle, \left|1\right\rangle, \left|2\right\rangle, \dots$, an infinite-dimensional space, so we **truncate** it: a matrix of dimension $n$ keeps $\left|0\right\rangle$ through $\left|n-1\right\rangle$ and discards $\left|n\right\rangle$ upward, and that $n$ is the **cutoff**. We choose the cutoff high enough that the retained top level has negligible population, then check by raising it.

The **annihilation** and **creation** operators step down and up the ladder, $\hat a\left|n\right\rangle = \sqrt{n}\left|n-1\right\rangle$ and $\hat a^\dagger\left|n\right\rangle = \sqrt{n+1}\left|n+1\right\rangle$, so $\hat a$ carries $\sqrt{n}$ on its superdiagonal:

```wl
annihilation[n_] := SparseArray[Band[{1, 2}] -> Sqrt[Range[n - 1]], {n, n}];
creation[n_] := ConjugateTranspose[annihilation[n]];
```

Finite truncation changes the canonical commutator only at the top retained level:

```wl
With[{n = 7, aa = annihilation[7]},
 aa . ConjugateTranspose[aa] - ConjugateTranspose[aa] . aa ==
  IdentityMatrix[n] - n SparseArray[{{n, n} -> 1}, {n, n}]]
```

It returns `True`. Thus $[a,a^\dagger]=\mathbb 1$ on all states below the boundary, with the required finite-dimensional correction on $\left|n-1\right\rangle$; cutoff checks ensure that correction is dynamically irrelevant.

The **coherent state** $\left|\alpha\right\rangle$ (the most classical state, an eigenstate of $\hat a$) is the Fock-space sum
$$\left|\alpha\right\rangle = e^{-\left|\alpha\right|^2/2}\sum_{n=0}^{\infty}\frac{\alpha^n}{\sqrt{n!}}\left|n\right\rangle.$$
Build it on the coefficients $\alpha^n/\sqrt{n!}$ through the highest retained level, normalized in the truncated basis (that finite-sum normalization replaces the infinite-series prefactor $e^{-\left|\alpha\right|^2/2}$):

```wl
coherentState[n_, a_] := With[{v = Prepend[Table[a^k/Sqrt[k!], {k, n - 1}], 1]}, v/Norm[v]];
```

The **displacement operator** $\hat D(\alpha) = e^{\alpha\hat a^\dagger - \alpha^*\hat a}$, which shifts the vacuum to $\left|\alpha\right\rangle$:

```wl
displacement[n_, a_] := MatrixExp[a creation[n] - Conjugate[a] annihilation[n]];
```

The **Wigner function** $W(x,p)$ pictures a state on the phase-space plane: positive where the state is concentrated and, for a nonclassical state, negative in places. In the canonical $(x,p)$ normalization its displaced-parity definition, with $\hat{Π} = e^{i\pi\hat a^\dagger\hat a}$ the photon-number parity operator, is
$$W(x,p) = \frac{1}{\pi}\,\mathrm{Tr}\left[\hat D^\dagger(\alpha)\rho\hat D(\alpha)\hat{Π}\right], \qquad \alpha = \frac{x + i p}{\sqrt{2}}.$$

A word on the coordinate convention, kept fixed for every phase-space picture below. The axes $x$ and $p$ are the dimensionless quadratures $\hat x = (\hat a + \hat a^\dagger)/\sqrt2$, $\hat p = (\hat a - \hat a^\dagger)/(i\sqrt2)$ (so $[\hat x,\hat p] = i$ and the vacuum $\left|0\right\rangle$ has variances $\Sigma_{\mathrm{xx}} = \Sigma_{\mathrm{pp}} = 1/2$), and $\alpha = (x + ip)/\sqrt2$ is the complex amplitude. The prefactor $1/\pi$ is the textbook normalization on these axes, $\int W\,dx\,dp = 1$, so the vacuum is $W(x,p) = \tfrac1\pi e^{-x^2-p^2}$.

Rather than exponentiate the displacement operator at every grid point, evaluate this trace directly in the truncated Fock basis, writing $n$ for the cutoff: $\rho$ is then the $n\times n$ matrix on the states $\left|0\right\rangle, \dots, \left|n-1\right\rangle$, and every sum below runs over its $n$ superdiagonals $\ell = 0, \dots, n-1$. Define
$$z = 2\alpha = \sqrt{2}(x + i p), \qquad b = \left|z\right|^2 = 4\left|\alpha\right|^2.$$
For the $\ell$th superdiagonal of $\rho$, define the radial Laguerre sum
$$F_\ell(b) = \sum_{m=0}^{n-\ell-1}\rho_{m,m+\ell}(-1)^m\sqrt{\frac{m!\,\ell!}{(m+\ell)!}}\,L_m^{(\ell)}(b),$$
where $L_m^{(\ell)}$ is a generalized Laguerre polynomial. The complete Wigner function is then
$$W(x,p) = \frac{1}{\pi}e^{-b/2}\,\mathrm{Re}\left[F_0(b) + 2\sum_{\ell=1}^{n-1}\frac{z^\ell}{\sqrt{\ell!}}F_\ell(b)\right].$$
In other words, the recurrence computes the radial function $F_\ell$ for one diagonal of $\rho$; the powers $z^\ell$ supply the angular dependence, and the Gaussian $e^{-b/2}$ supplies the envelope. Build the whole evaluation as one function: the Laguerre recurrence lives inside as a local helper, the input may be a state vector or a density matrix (a vector becomes its density matrix first), and the cutoff is read off the input's size. Each grid value evaluates the finite-cutoff expression; the returned interpolating function is a numerical approximation between those grid points, with a wide enough frame that its far corners stay clean:

```wl
wignerRepresentation[psi_?VectorQ, reach_, pts_] := wignerRepresentation[densityMatrix[psi], reach, pts];
wignerRepresentation[rho_?MatrixQ, reach_, pts_] := Module[{n = Length[rho], grid, a2, b, lag, w},
   grid = N@Subdivide[-reach, reach, pts - 1];
   a2 = Sqrt[2.] Outer[#2 + I #1 &, grid, grid]; b = Abs[a2]^2;
   lag[L_, c_] := With[{m = Length[c]},
     Which[
      m == 1, c[[1]] + 0. b,
      m == 2, c[[1]] - c[[2]] Sqrt[1./(L + 1)] (L + 1. - b),
      True, With[{yy = Fold[Function[{y, cj}, With[{j = cj[[2]]},
           {cj[[1]] - y[[2]] Sqrt[((j - 1.) (L + j - 1.))/((L + j) j)],
            y[[1]] - y[[2]] (L + 2. j - 1 - b) Sqrt[1./((L + j) j)]}]],
         {c[[-2]] + 0. b, c[[-1]] + 0. b}, Transpose[{Reverse[c[[;; -3]]], N@Range[m - 1, 2, -1]}]]},
        yy[[1]] - yy[[2]] Sqrt[1./(L + 1)] (L + 1. - b)]]];
   w = Re[Fold[lag[#2 - 1, If[#2 == 1, 1, 2] Diagonal[rho, #2 - 1]] + #1 a2/Sqrt[#2] &,
       ConstantArray[2. rho[[1, -1]], {pts, pts}], Range[n - 1, 1, -1]]] Exp[-b/2]/Pi;
   ListInterpolation[Transpose[w], {{-reach, reach}, {-reach, reach}}]];
```

Check it on the vacuum, whose Wigner function is the closed form $\tfrac{1}{\pi}e^{-x^2-p^2}$, feeding the state vector directly:

```wl
With[{w = wignerRepresentation[coherentState[8, 0], 3., 25]},
 Max@Table[Abs[w[xx, pp] - (1/Pi) Exp[-xx^2 - pp^2]], {xx, -2, 2, 0.5}, {pp, -2, 2, 0.5}]]
```

The sampled discrepancy is at floating-point precision on the grid's own points; between grid nodes, accuracy is limited by the interpolation resolution.

Wigner pictures are compared, across times or across states, and the comparison is only honest on one shared color scale. The plotter takes a list of Wigner interpolating functions and a label per member, finds the largest $\left|W\right|$ over the whole list, and draws every panel and the bar legend with the same diverging map: a light neutral at zero so empty phase space reads as empty, warm (orange) for positive, cold (blue) for negative, with the $W=0$ contour drawn as a thin gray line so the sign change is unmistakable. An optional color map replaces the default; it receives the signed fraction $W/\max\left|W\right|$ in $[-1,1]$:

```wl
ClearAll[wignerPlot]
wignerPlot[ws : {__InterpolatingFunction}, labels_List, color_ : Automatic] :=
  With[{m = Max[Abs[#["ValuesOnGrid"]] & /@ ws], reach = Max[Abs[ws[[1]]["Domain"]]],
    cf = If[color === Automatic,
      Blend[{{-1, RGBColor[0.15, 0.4, 0.8]}, {0, GrayLevel[0.97]}, {1, RGBColor[0.9, 0.45, 0.1]}}, #] &, color]},
   With[{raw = cf[#/m] &},
    Legended[GraphicsRow[MapThread[
       Show[DensityPlot[#1[xx, pp], {xx, -reach, reach}, {pp, -reach, reach}, PlotPoints -> 55,
          PlotRange -> All, ColorFunction -> raw, ColorFunctionScaling -> False,
          FrameLabel -> {"x", "p"}, PlotLabel -> #2],
         ContourPlot[#1[xx, pp] == 0, {xx, -reach, reach}, {pp, -reach, reach},
          ContourStyle -> Directive[Thin, Gray], PlotPoints -> 40]] &,
       {ws, labels}]],
     BarLegend[{raw, {-m, m}}, LegendMarkerSize -> 130]]]];
wignerPlot[w_InterpolatingFunction, label_, color_ : Automatic] := wignerPlot[{w}, {label}, color];
```

### A Cavity Cat, and Why It Decoheres Almost at Once

**The problem.** A single leaky mode of light (a cavity dribbling photons), with no drive:
$$\dot\rho = \gamma\,\mathcal{D}[\hat a]\rho.$$
A lone coherent state only shrinks, staying pure, $\left|\alpha\right\rangle \to \left|\alpha e^{-\gamma t/2}\right\rangle$, but a *cat* $\left|\alpha\right\rangle + \left|\beta\right\rangle$ of two well-separated amplitudes keeps its cross-term only through a coherence factor
$$\left|C(t)\right| = \exp\left(-\tfrac12\,\left|\alpha-\beta\right|^2\,(1 - e^{-\gamma t})\right) \;\approx\; \exp\left(-\tfrac12\,\left|\alpha-\beta\right|^2\,\gamma t\right)\ \text{at short time,}$$
so the coherence decays at a rate $\propto \left|\alpha-\beta\right|^2$, quadratic in the phase-space separation. Build it in a truncated Fock basis and watch the fringes vanish.

Set a high cutoff and a cavity decay rate:

```wl
topCat = 30; blankCat = IdentityMatrix[topCat]; downCat = annihilation[topCat]; \[Gamma]cat = 0.7;
```

First a single coherent state: under damping it stays coherent, its amplitude decaying as $\alpha \to \alpha e^{-\gamma t/2}$. Check that the damped state equals a fresh, smaller coherent state:

```wl
With[{later = evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[coherentState[topCat, 2.]], 1.0]},
 Re[Conjugate[coherentState[topCat, 2. Exp[-\[Gamma]cat/2]]] . later . coherentState[topCat, 2. Exp[-\[Gamma]cat/2]]]]
```

Overlap one: a coherent state stays coherent.

Now the interesting state, a **cat**, the superposition of two coherent states $\left|\alpha\right\rangle + \left|{-}\alpha\right\rangle$:

```wl
\[Alpha] = 2.; cat = Normalize[coherentState[topCat, \[Alpha]] + coherentState[topCat, -\[Alpha]]];
```

First a rigorous check on the solver. Amplitude damping has an exact operator-sum (Kraus) solution, summed over the number of photons lost $m$,
$$\rho(t) = \sum_{m\ge 0} E_m\,\rho_0\,E_m^\dagger, \qquad E_m = \sqrt{\frac{(1-e^{-\gamma t})^{m}}{m!}}\;e^{-\gamma t\,\hat a^\dagger\hat a/2}\,\hat a^{m},$$
with $\hat a^m$ removing $m$ photons and $e^{-\gamma t\,\hat a^\dagger\hat a/2}$ shrinking the survivors; the set is trace-preserving, so it needs no renormalization. Write it:

```wl
photonLoss[rho0_, t_] := With[{shrink = DiagonalMatrix[Exp[-\[Gamma]cat t Range[0, topCat - 1]/2]], keep = 1 - Exp[-\[Gamma]cat t]},
   Sum[With[{M = keep^(m/2)/Sqrt[m!] shrink . If[m == 0, blankCat, MatrixPower[downCat, m]]},
      M . rho0 . ConjugateTranspose[M]], {m, 0, topCat - 1}]];
```

Confirm the master equation and this Kraus sum agree on the cat:

```wl
Max@Abs@Flatten[evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], 0.6] - photonLoss[densityMatrix[cat], 0.6]]
```

The discrepancy is at the numerical integration tolerance.

Now visualize the cat's decoherence. Build the Wigner film over a short window:

```wl
catTimes = Range[0, .4, .1];
catFilm = Table[wignerRepresentation[evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], t], 4.2, 55], {t, catTimes}];
```

Draw it on one shared color scale:

```wl
wignerPlot[catFilm, "time " <> ToString[#] & /@ catTimes]
```

The thin gray lines in each panel are the $W = 0$ contour, where the Wigner function changes sign: the negative (blue) patches of quantum interference lie inside them.

At first there are two Gaussian lobes (the two coherent states) with interference fringes between them: the fringes are the coherence, the off-diagonal cross-term. Across this short window the fringes fade fast while the lobes barely shrink: the cat is collapsing toward a classical mixture of the two coherent states while both are still nearly full size. Decoherence is far faster than energy loss.

How fast? The coherence $\left|C(t)\right|$ is the closed-form envelope from the problem statement, and the master equation can be made to report the same quantity on its own. Evolve the cat once over the plotted span:

```wl
catRun = evolveODE[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], 2.];
```

Under amplitude damping the state stays in the span of its two shrinking lobes, so it keeps the form
$$\rho(t) = p\,(|\alpha_t\rangle\langle\alpha_t| + |-\alpha_t\rangle\langle-\alpha_t|) + c\,(|\alpha_t\rangle\langle-\alpha_t| + |-\alpha_t\rangle\langle\alpha_t|), \qquad \alpha_t = \alpha e^{-\gamma t/2},$$
and the coherence is the weight ratio $\left|C(t)\right| = c/p$. The two lobes are not orthogonal (they overlap by $o = \langle\alpha_t|{-\alpha_t}\rangle$), so read $p$ and $c$ off two matrix elements of $\rho(t)$ by undoing that overlap with a two-by-two linear solve:

```wl
catCoherence[t_?NumericQ] := Module[{va, vb, o, m1, m2, p, cc},
   va = coherentState[topCat, \[Alpha] Exp[-\[Gamma]cat t/2]]; vb = coherentState[topCat, -\[Alpha] Exp[-\[Gamma]cat t/2]];
   o = Re[Conjugate[va] . vb];
   m1 = Re[Conjugate[va] . catRun[t] . va]; m2 = Re[Conjugate[va] . catRun[t] . vb];
   {p, cc} = LinearSolve[{{1 + o^2, 2 o}, {2 o, 1 + o^2}}, {m1, m2}];
   cc/p];
```

Now visualize both, the closed-form envelope and the master-equation extraction, for our chosen separation:

```wl
Plot[{Exp[-2 \[Alpha]^2 (1 - Exp[-\[Gamma]cat t])], catCoherence[t]}, {t, 0, 2},
 PlotStyle -> {Automatic, Dashed},
 PlotLegends -> {"closed form", "master equation"}, Frame -> True, GridLines -> Automatic,
 FrameLabel -> {"time", "coherence |C(t)|"}, PlotRange -> All,
 PlotLabel -> "master-equation coherence lies on the closed-form envelope"]
```

The dashed master-equation extraction lies on the closed form within the displayed numerical tolerance. The envelope is the exact infinite-dimensional amplitude-damping law; the finite Fock simulation approaches it as its cutoff and solver accuracy are increased. The decoherence rate grows as the *square* of the separation, $\left|\alpha - (-\alpha)\right|^2 = 4\left|\alpha\right|^2$.

Every Fock-space example truncates the infinite number-state ladder at a finite cutoff, so each needs a check that the cutoff is high enough not to change the answer. Recompute the mean photon number of a decaying coherent state at cutoffs 20, 30, and 40:

```wl
With[{observable = Function[nn, expectation[creation[nn] . annihilation[nn], evolve[0 IdentityMatrix[nn],
      {Sqrt[\[Gamma]cat] annihilation[nn]}, densityMatrix[coherentState[nn, 2.]], 0.5]]]},
 observable /@ {20, 30, 40}]
```

The three values are numerically indistinguishable at the shown precision. This bounds, rather than eliminates, the truncation error: for this observable and state, raising the dimension beyond 20 changes the result by less than the displayed tolerance.

### Making a Cat With an Atom, and Reading Its Decoherence

**The problem.** A *far-detuned* atom, its transition frequency sitting far from the cavity's so it can trade no energy with the light and only turns its phase, imprints a state-dependent phase shift under
$$\hat H_C = \chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\qquad \left|\alpha\right\rangle \to \left|\alpha e^{\pm i\phi}\right\rangle,\quad \phi = \chi\tau.$$
An atom prepared in the equal superposition $(\left|e\right\rangle + \left|g\right\rangle)/\sqrt2$ therefore leaves the field in a two-phase **cat**, $\left|\alpha e^{-i\phi}\right\rangle + \left|\alpha e^{+i\phi}\right\rangle$, entangled with the atom. After the cat decoheres for a delay (the leak of the last example), a second atom probes the same field, and the two-atom correlation $\eta = p(e_2|e_1) - p(e_2|g_1)$ (the second atom's chance of coming out excited when the first did, minus when the first came out in the ground state) is nonzero only while the field keeps its coherence, so $\eta(\tau)$ tracks the surviving coherence $C(\tau)$ and falls to zero as the cat blurs into a mixture. This is the [Haroche experiment](https://link.aps.org/doi/10.1103/PhysRevLett.77.4887); build the atom-made cat and measure how fast the correlation dies.

In plainer terms, the atom is a phase stamp, not an absorber: far off resonance it takes no photon from the field, it only rotates the field's phase by $+\phi$ or $-\phi$ according to its state, so an atom in $(\left|e\right\rangle+\left|g\right\rangle)/\sqrt2$ leaves the two phase-rotated fields superposed and entangled with it. One probe atom cannot reveal this coherence, since its own outcome averages over the two phases; only the correlation $\eta$ between two atoms reading the same field exposes it, and only while it survives.

The protocol runs as a fixed sequence: the first atom writes the cat onto the field, a delay lets it decohere, and the second atom reads it back. The grid below is that arc as a chain of states, one line per operation, running from the initial coherent field through the written cat and its decoherence to the final two-atom readout. Reading the first atom is what conditions the cat into being; discard its record instead and the field is left a diagonal mixture of the two phase-pointers, no cross term and so no fringes.

```wl
#| file: grid_states.wxf
#| style: Output
```

Set the cutoff, the cavity decay rate, and the number operator:

```wl
topHar = 28; blankHar = IdentityMatrix[topHar]; downHar = annihilation[topHar]; numHar = creation[topHar] . downHar; \[Gamma]har = 0.5;
```

The **trace distance** $\tfrac12\mathrm{Tr}\left|\rho_1 - \rho_2\right|$ will measure how far a field sits from a classical mixture:

```wl
traceDistance[a_, b_] := Total[SingularValueList[a - b]]/2;
```

Because the singular values of any matrix sum to its trace norm, this is exactly $\tfrac12\mathrm{Tr}\left|a-b\right|$; for the Hermitian difference of two states those singular values are the magnitudes of its eigenvalues.

Now the machinery that makes the cat with an atom, three pieces on the joint atom-cavity space (atom first, cavity second). A **Ramsey $\pi/2$ pulse** rotates the atom by $\pi/2$ about $\hat\sigma_y$, turning $\left|g\right\rangle$ into an equal superposition. On the joint atom-cavity space it is the operator $e^{-i(\pi/4)\hat\sigma_y}\otimes\mathbb{1}_{\mathrm{cavity}}$, the $\pi/2$ rotation appearing as $\pi/4$ in the exponent since the generator turns the state by twice its angle:

```wl
ramsey = KroneckerProduct[MatrixExp[-I (Pi/4) Y], blankHar];
```

The **dispersive interaction** $e^{-i\phi\,\hat a^\dagger\hat a\,\hat\sigma_z}$ rotates the field's phase by $-\phi$ when the atom is $\left|e\right\rangle$ and by $+\phi$ when it is $\left|g\right\rangle$, the phase stamp of the problem statement. Since $\hat\sigma_z$ is diagonal the exponential splits into one field rotation per atom state, $|e\rangle\langle e|\otimes e^{-i\phi\,\hat a^\dagger\hat a} + |g\rangle\langle g|\otimes e^{+i\phi\,\hat a^\dagger\hat a}$:

```wl
disperse[\[Phi]_] := KroneckerProduct[{{1, 0}, {0, 0}}, MatrixExp[-I \[Phi] numHar]] + KroneckerProduct[{{0, 0}, {0, 1}}, MatrixExp[I \[Phi] numHar]];
```

And a function returning the cavity conditioned on finding the atom in $\left|e\right\rangle$ (block 1) or $\left|g\right\rangle$ (block 2): the partial matrix element $\langle s|\rho|s\rangle$ over the atom state $|s\rangle$, which brackets only the atom and leaves the atom's diagonal block of the joint $\rho$, a cavity operator whose trace is the probability of that outcome:

```wl
cavityGiven[joint_, s_] := With[{r = ArrayReshape[joint, {2, topHar, 2, topHar}]}, r[[s, All, s, All]]];
```

Start the atom in $\left|g\right\rangle$ and the cavity in the coherent state of amplitude $\alpha = 2$:

```wl
jointStart = densityMatrix[Flatten@KroneckerProduct[ground, coherentState[topHar, 2.]]];
```

Ramsey the atom into a superposition, then let the dispersive interaction entangle them at phase $\phi$. Applied to $|g\rangle|\alpha\rangle$ the two together produce the pure state $|\psi_\phi\rangle\langle\psi_\phi|$ with
$$|\psi_\phi\rangle = \tfrac{1}{\sqrt2}(|g\rangle|\alpha_+\rangle - |e\rangle|\alpha_-\rangle), \qquad \alpha_\pm = \alpha\,e^{\pm i\phi},$$
each atom state now entangled with its own phase-rotated field:

```wl
entangle[\[Phi]_] := With[{U = disperse[\[Phi]]}, U . (ramsey . jointStart . ConjugateTranspose[ramsey]) . ConjugateTranspose[U]];
```

Here is the step that is easy to skip. If we now *ignore* the atom and trace it out, the cavity is not a cat. Tracing out the atom sums its two diagonal blocks, $\mathrm{Tr}_{\mathrm{atom}}\,\rho = \langle e|\rho|e\rangle + \langle g|\rho|g\rangle$, and for the entangled state that collapses to the equal mixture $\tfrac12(|\alpha_+\rangle\langle\alpha_+| + |\alpha_-\rangle\langle\alpha_-|)$, no cross term:

```wl
reducedCavity[\[Phi]_] := cavityGiven[entangle[\[Phi]], 1] + cavityGiven[entangle[\[Phi]], 2];
```

Read its purity:

```wl
purity[reducedCavity[1.3]]
```

against the pure cat $\propto |\alpha_-\rangle + |\alpha_+\rangle$ the naive picture would insert by hand:

```wl
purity[densityMatrix@Normalize[coherentState[topHar, 2. Exp[-I 1.3]] + coherentState[topHar, 2. Exp[I 1.3]]]]
```

A half against a one: the unread atom leaves the cavity a **mixture** of the two phase-rotated coherent states, not a pure cat. Its Wigner function says the same, two lobes with nothing between them:

```wl
wignerPlot[wignerRepresentation[reducedCavity[1.3], 4.0, 55], "atom ignored: a fringeless mixture"]
```

Two clean lobes, no fringes: the which-phase information has leaked into the atom, and averaging over it erases the coherence, exactly as tracing out any unwatched detector does.

The cat appears only when the atom is *read*, and the basis it is read in decides what comes out. Because $\left|e\right\rangle$ rode with $\left|\alpha e^{-i\phi}\right\rangle$ and $\left|g\right\rangle$ with $\left|\alpha e^{+i\phi}\right\rangle$, reading the atom directly in $\{\left|e\right\rangle, \left|g\right\rangle\}$ would reveal which phase the field took and collapse it to that one coherent state, a pointer and not a cat. The second Ramsey pulse, the analysis pulse, first rotates the atom into the $\left|\pm\right\rangle = (\left|e\right\rangle \pm \left|g\right\rangle)/\sqrt2$ basis, so its outcome cannot say which phase; keep the runs where it comes out $\left|e\right\rangle$ and that measurement projects the cavity into the even superposition $\left|\alpha e^{-i\phi}\right\rangle + \left|\alpha e^{+i\phi}\right\rangle$ of the two pointers, the cat; the $\left|g\right\rangle$ outcome gives the odd superposition $\left|\alpha e^{-i\phi}\right\rangle - \left|\alpha e^{+i\phi}\right\rangle$. Condition on the atom after the analysis pulse. In symbols $\mathrm{conditionalCat}(\phi,s) = \langle s|R\,\rho_\phi\,R^\dagger|s\rangle\,/\,\mathrm{Tr}[\cdots]$: the analysis pulse $R=$ `ramsey` on the entangled $\rho_\phi$, then the atom-$s$ block, renormalized:

```wl
conditionalCat[\[Phi]_, s_] := With[{r = ramsey . entangle[\[Phi]] . ConjugateTranspose[ramsey]}, With[{blk = cavityGiven[r, s]}, blk/Tr[blk]]];
```

Read its purity:

```wl
purity[conditionalCat[1.3, 1]]
```

Purity one: reading the atom has turned the mixture into a pure cat. Its Wigner function carries the fringes the mixture lacked:

```wl
wignerPlot[wignerRepresentation[conditionalCat[1.3, 1], 4.0, 55], "atom read as |e\[RightAngleBracket]: a cat, fringes and all"]
```

The interference fringes are back. Making the cat and reading the atom are one act: no cat is inserted, it is conditioned into being by the atomic measurement. The phase $\phi$ still sets the separation: small $\phi$ leaves the two pointers close, large $\phi$ far apart.

Now watch the conditional cat decohere under the leak. Sample three delays:

```wl
harTimes = {0., 0.5, 1.0};
```

Evolve the small-phase cat under the cavity leak across that span:

```wl
harSmall = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, conditionalCat[0.5, 1], Max[harTimes]];
```

and the large-phase cat the same way:

```wl
harLarge = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, conditionalCat[1.3, 1], Max[harTimes]];
```

Turn each into a Wigner film, one frame per delay, the small-phase cat first:

```wl
harFilmSmall = Table[wignerRepresentation[harSmall[t], 4.0, 45], {t, harTimes}];
```

and the large-phase cat:

```wl
harFilmLarge = Table[wignerRepresentation[harLarge[t], 4.0, 45], {t, harTimes}];
```

Draw both films on one shared color scale so the two separations compare directly, the small-phase triple then the large-phase triple:

```wl
wignerPlot[Join[harFilmSmall, harFilmLarge],
 Join["small \[Phi], \[Tau]=" <> ToString[#] & /@ harTimes, "large \[Phi], \[Tau]=" <> ToString[#] & /@ harTimes]]
```

Read the two triples against each other: the small-phase cat keeps its fringes across the window, while the large-phase cat, its two states farther apart, has collapsed to a fringeless mixture by the end. Same cavity, same leak, same clock; only the separation differs, and the fringes die at a rate set by its square.

Turn that into a number. The mixture the cat decays toward is the equal blend of its two pointer states, each amplitude shrinking with the cavity,
$$\rho_{\mathrm{mix}}(\phi,t) = \tfrac12[|\alpha_+(t)\rangle\langle\alpha_+(t)| + |\alpha_-(t)\rangle\langle\alpha_-(t)|], \qquad \alpha_\pm(t) = 2\,e^{-\gamma t/2}\,e^{\pm i\phi},$$
with $\gamma = \gamma_{\mathrm{har}}$:

```wl
mixture[\[Phi]_, t_] := With[{sh = Exp[-\[Gamma]har t/2]},
   (densityMatrix[coherentState[topHar, 2. sh Exp[-I \[Phi]]]] + densityMatrix[coherentState[topHar, 2. sh Exp[I \[Phi]]]])/2];
```

The surviving coherence is the trace distance between the leak-evolved conditional cat and that mixture, over a list of delays, one master-equation run per phase:

```wl
coherenceLeft[\[Phi]_, delays_] := With[{run = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, conditionalCat[\[Phi], 1], Max[delays]]},
   Table[traceDistance[run[t], mixture[\[Phi], t]], {t, delays}]];
delays = Range[0, 3, 0.25];
{close, far} = {coherenceLeft[0.5, delays], coherenceLeft[1.3, delays]};
```

Visualize the two normalized coherences against the delay:

```wl
ListLinePlot[{Transpose[{delays, close/First@close}], Transpose[{delays, far/First@far}]},
 Frame -> True, GridLines -> Automatic, PlotMarkers -> Automatic,
 PlotLegends -> {"small phase shift", "large phase shift"}, PlotRange -> All,
 FrameLabel -> {"delay before the second atom reads", "surviving coherence"},
 PlotLabel -> "cat coherence decays faster for a larger phase shift"]
```

Both fall, but the large-$\phi$ curve plunges: the separation-squared decoherence law, the separation now set by the atom's phase shift.

That trace distance is the theorist's coherence. The experiment reads it out with a *second atom*, and that two-atom correlation is what the section is named for. Send a second probe (in $\left|g\right\rangle$, Ramsey'd) through the decohered field, close its Ramsey, and read its excited chance; the correlation is $\eta = p(e_2|e_1) - p(e_2|g_1)$, the second atom's excited probability conditioned on the first atom's two outcomes. The probe's excited chance on a given field $\rho$ is $\mathrm{Tr}\,\langle e|\,V(|g\rangle\langle g|\otimes\rho)\,V^\dagger\,|e\rangle$ with $V = R\,D(\phi)\,R$ the probe's own Ramsey-stamp-Ramsey, the $|e\rangle$-block trace reading its excited probability:

```wl
probe[\[Phi]_, cav_] := With[{U = disperse[\[Phi]], sig = KroneckerProduct[{{0, 0}, {0, 1}}, cav]},
   With[{s3 = ramsey . (U . (ramsey . sig . ConjugateTranspose[ramsey]) . ConjugateTranspose[U]) . ConjugateTranspose[ramsey]},
    Re@Tr[cavityGiven[s3, 1]]]];
```

The correlation differences that probe on the two conditioned fields, each decohered for the delay, $\eta(\tau) = \mathrm{probe}(\rho_{e_1}(\tau)) - \mathrm{probe}(\rho_{g_1}(\tau)) = p(e_2|e_1) - p(e_2|g_1)$:

```wl
correlation[\[Phi]_, delays_] := With[{runE = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, conditionalCat[\[Phi], 1], Max[delays]], runG = evolveODE[0 blankHar, {Sqrt[\[Gamma]har] downHar}, conditionalCat[\[Phi], 2], Max[delays]]},
   Table[probe[\[Phi], runE[t]] - probe[\[Phi], runG[t]], {t, delays}]];
```

For the well-separated cat, plot the two-atom correlation against the trace-distance coherence, each normalized to its start:

```wl
etaFar = correlation[1.3, delays];
ListLinePlot[{Transpose[{delays, far/First@far}], Transpose[{delays, etaFar/First@etaFar}]},
 Frame -> True, GridLines -> Automatic, PlotMarkers -> Automatic,
 PlotLegends -> {"coherence (trace distance)", "two-atom correlation \[Eta]"}, PlotRange -> All,
 FrameLabel -> {"delay before the second atom reads", "normalized signal"},
 PlotLabel -> "the two-atom correlation traces the coherence"]
```

The two curves fall together: the two-atom correlation, the thing a laboratory actually records, is a faithful readout of the field coherence the trace distance measures, and it dies as the cat blurs into the mixture. (At small phase the raw $\eta$ additionally carries a slow Ramsey fringe, since the shrinking pointers keep rotating as they fade; for the well-separated cat here the coherence is gone before that fringe turns, so $\eta$ traces it cleanly.) This is the Haroche experiment: the cat is written on the field by one atom and read by the next, and its life is the coherence between two places.

### Dispersive Readout: A Qubit Measured by a Cavity

**The problem.** A qubit read out by a driven, damped cavity whose frequency it shifts. The coupling is dispersive, so the qubit's $\hat\sigma_z$ pushes the cavity resonance one way or the other, turning the cavity into a meter for it:
$$\dot\rho = -i[\epsilon(\hat a + \hat a^\dagger),\rho] - i[\chi\,\hat a^\dagger\hat a\,\hat\sigma_z,\rho] + \gamma\,\mathcal{D}[\hat a]\rho.$$
Two idealizations are folded into that compactness. It is written in the frame rotating at the drive frequency for the cavity and the qubit frequency for the qubit, with the cavity driven on resonance ($\Delta_c = \omega_c - \omega_d = 0$), so the bare $\omega_c\hat a^\dagger\hat a$ and $\tfrac{\omega_q}{2}\hat\sigma_z$ are rotated away and $\hat\sigma_z$ survives only as a static label that shifts the cavity resonance by $\pm\chi$. And the only bath is the cavity's own leak: the qubit is given no intrinsic relaxation or dephasing, so all of its decoherence is the leak it induces on the cavity. The drive $\epsilon$ and damping $\gamma$ settle the cavity into one of two qubit-conditioned coherent states, the meter's two readings, while the cavity's own leak decoheres any superposition of the two qubit states. Build the joint qubit-plus-cavity system and watch the two readings pull apart as the qubit's coherence blurs.

A qubit is read out by a driven, damped cavity whose frequency it shifts (dispersive coupling). The cavity settles into one of two qubit-conditioned coherent states, each a dot sitting in one of two places in the phase plane: the cavity is the meter, and these two settled amplitudes are its two **pointer states** (Zurek's name for the apparatus states that record a measurement outcome). Set the drive, coupling, and decay:

```wl
\[Epsilon] = 1.2; \[Chi] = 0.6; \[Gamma]disp = 1.0;
```

Conditional on the qubit in a $\hat\sigma_z$ eigenstate $s$, the cavity sees $H_s = \epsilon(\hat a + \hat a^\dagger) + s\chi\,\hat a^\dagger\hat a$ and leaks through $\hat a$, a linear harmonic problem that keeps a coherent state coherent; its mean $\alpha_s = \mathrm{Tr}(\hat a\rho_s)$ therefore closes on itself, and a single amplitude carries the whole conditional cavity (confirmed against the full master equation below).

The equation that amplitude obeys is the generator's own action on $\hat a$, not a coefficient put in by hand. The adjoint (Heisenberg) generator $\mathcal L_s^{\dagger}[\hat a] = i[H_s,\hat a] + \gamma\,\mathcal D^{\dagger}[\hat a](\hat a)$ is, in the Fock basis, a constant plus a multiple of $\hat a$, and for each qubit state the pair is $\{-i\epsilon,\,-(\tfrac{\gamma}{2}+is\chi)\}$: the adjoint dissipator $\mathcal D^{\dagger}[\hat a](\hat a) = -\tfrac12\hat a$ supplies the $\gamma/2$, and $[\hat a^\dagger\hat a,\hat a] = -\hat a$ the $is\chi$. So the amplitude obeys $\dot\alpha = -i\epsilon - (is\chi + \tfrac{\gamma}{2})\alpha$ with $\alpha(0)=0$, one per qubit state ($s = +1$ for $\left|e\right\rangle$, $s = -1$ for $\left|g\right\rangle$), giving in closed form
$$\alpha_s(t) = \frac{-i\epsilon}{is\chi + \gamma/2}(1 - e^{-(is\chi + \gamma/2)\,t}),$$
settling to the mirror-image steady amplitudes $\alpha_s^{\mathrm{ss}} = -i\epsilon/(is\chi + \gamma/2)$. Solve it:

```wl
pointer[which_] := First[amp /. DSolve[{amp'[t] == -I \[Epsilon] - (I which \[Chi] + \[Gamma]disp/2) amp[t], amp[0] == 0}, amp, t]];
```

Evaluate the two pointers, one per qubit state:

```wl
{\[Alpha]e, \[Alpha]g} = {pointer[1], pointer[-1]};
```

Read the two late-time amplitudes:

```wl
Chop[{\[Alpha]e[9.], \[Alpha]g[9.]}]
```

The two steady amplitudes are mirror images, one per qubit state.

Now visualize the pointer-state separation by drawing both trajectories from the origin to their steady points:

```wl
ListLinePlot[{Table[Sqrt[2] {Re@\[Alpha]e[t], Im@\[Alpha]e[t]}, {t, 0, 9, 0.03}], Table[Sqrt[2] {Re@\[Alpha]g[t], Im@\[Alpha]g[t]}, {t, 0, 9, 0.03}]},
 PlotLegends -> {"pointer if qubit |e\[RightAngleBracket]", "pointer if qubit |g\[RightAngleBracket]"}, Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"x", "p"}, PlotLabel -> "the meter pointer settles to two places, one per qubit state"]
```

What a real readout does with those two pointers is discriminate them in the IQ plane. Detecting the cavity output (heterodyne) turns each qubit-conditioned coherent pointer into a Gaussian blob of measurement outcomes, centered at its amplitude with the shot-noise width of a coherent state. A snapshot draws $(I,Q) = \sqrt2\,(\mathrm{Re}\,\alpha_s(t),\,\mathrm{Im}\,\alpha_s(t)) + (\xi_I,\xi_Q)$ with independent $\xi_I,\xi_Q\sim\mathcal N(0,1)$, Gaussian shot noise added to the pointer amplitude; because the conditional cavity is the coherent state $|\alpha_s(t)\rangle$, that draw is precisely its Husimi-$Q$ distribution. As the pointers separate the two blobs pull apart and the assignment becomes reliable. Sample the two outcome clouds at three readout times:

```wl
iqCloud[ptr_, t_, seed_] := BlockRandom[SeedRandom[seed];
   Table[Sqrt[2] {Re@ptr[t], Im@ptr[t]} + RandomVariate[NormalDistribution[0, 1], 2], 250]];
iqTimes = {0.5, 1.5, 4.};
```

Now visualize the readout separating in time, the two clouds (qubit $\left|e\right\rangle$ in blue, $\left|g\right\rangle$ in orange) drifting apart across the IQ plane:

```wl
GraphicsRow[Table[
   ListPlot[{iqCloud[\[Alpha]e, t, 1], iqCloud[\[Alpha]g, t, 2]}, PlotStyle -> {Directive[Opacity[0.5], PointSize[0.015]], Directive[Opacity[0.5], PointSize[0.015]]},
    PlotRange -> {{-6, 6}, {-6, 6}}, AspectRatio -> 1, Frame -> True, GridLines -> Automatic, ImageSize -> 260,
    FrameLabel -> {"I", "Q"}, PlotLabel -> "t = " <> ToString[t]], {t, iqTimes}], ImageSize -> 800]
```

Early the clouds overlap and the qubit state is ambiguous; by the last frame they are cleanly resolved. Turn that into numbers: discriminating the two pointer states from a single snapshot of the cavity, at one instant, has error $\tfrac12\mathrm{erfc}(\left|\alpha_e - \alpha_g\right|/2)$, the overlap of the two Gaussians: the nearest-centroid error of the two heterodyne clouds, not the ultimate quantum limit over all measurements. Define the pointer separation, the distance between the two coherent amplitudes:

```wl
sepDisp[t_] := Abs[\[Alpha]e[t] - \[Alpha]g[t]];
```

The wider that separation, the smaller the overlap error the discrimination leaves. Plot both against readout time:

```wl
ListLinePlot[
 Transpose@Table[{{t, sepDisp[t]}, {t, 0.5 Erfc[sepDisp[t]/2]}}, {t, 0, 6, 0.05}],
 Frame -> True, GridLines -> Automatic, AspectRatio -> 1/2, PlotRange -> All,
 PlotLegends -> Placed[{"pointer separation |\!\(\*SubscriptBox[\(\[Alpha]\), \(e\)]\) - \!\(\*SubscriptBox[\(\[Alpha]\), \(g\)]\)|", "assignment error \!\(\*FractionBox[\(1\), \(2\)]\) erfc(|\!\(\*SubscriptBox[\(\[Alpha]\), \(e\)]\) - \!\(\*SubscriptBox[\(\[Alpha]\), \(g\)]\)|/2)"}, Right],
 FrameLabel -> {"readout time", "separation / error"},
 PlotLabel -> "instantaneous discrimination: \nseparation and error saturate once the cavity settles"]
```

The separation rises to its steady value and the instantaneous error falls with it, then both flatten: once the cavity has settled, a snapshot of it says the same thing at every later instant. But that snapshot is not what a real readout uses. The cavity leaks the whole time, and a heterodyne receiver integrates that output record, so its information keeps accumulating even after the cavity is stationary: the integrated signal-to-noise grows as $\mathrm{SNR}^2(t) \propto \gamma\int_0^t\left|\alpha_e(s) - \alpha_g(s)\right|^2\,ds$. Code that running integral of the squared pointer separation:

```wl
snrAccum[t_?NumericQ] := \[Gamma]disp NIntegrate[Abs[\[Alpha]e[s] - \[Alpha]g[s]]^2, {s, 0, t}];
```

Plot the saturating snapshot separation-squared against this accumulating integral:

```wl
ListLinePlot[{Table[{t, sepDisp[t]^2}, {t, 0, 6, 0.05}], Table[{t, snrAccum[t]}, {t, 0, 6, 0.1}]},
 Frame -> True, GridLines -> Automatic, AspectRatio -> 1/2, PlotRange -> All,
 PlotLegends -> Placed[{"pointer separation squared |alpha_e - alpha_g|^2 (saturates)", "integrated SNR^2 (grows)"}, Right],
 FrameLabel -> {"readout time", "distinguishability"},
 PlotLabel -> "the pointer separation saturates;\n the integrated output record keeps accumulating"]
```

The snapshot separation-squared levels off as the pointers reach their steady amplitudes, but the integrated $\mathrm{SNR}^2$ climbs without bound, linearly once the cavity is stationary. That growing line, not the saturating one, is the single-shot readout a circuit-QED experiment optimizes, and it is the flip side of the qubit coherence the same measurement destroys. See [Gambetta et al.](https://arxiv.org/abs/cond-mat/0602322).

Measurement has a cost: as the pointer separates and reveals the qubit state, the qubit's coherence is destroyed. To see it we need the joint qubit-cavity density matrix and two functions. First, embed a cavity operator on the joint space as $\mathbb{1}_{\mathrm{qubit}}\otimes M$, identity on the qubit:

```wl
onCavity[metering_] := KroneckerProduct[IdentityMatrix[2], metering];
```

Second, the partial trace over the cavity, recovering the qubit's reduced state $\rho_{\mathrm{qubit}} = \mathrm{Tr}_{\mathrm{cavity}}[\rho]$:

```wl
reducedQubit[joint_, rungs_] := TensorContract[ArrayReshape[joint, {2, rungs, 2, rungs}], {{2, 4}}];
```

Build a function for the qubit coherence as the pointer resolves the qubit: it returns $2\left|\rho_{\mathrm{eg}}(t)\right| = 2\left|\langle e|\rho_{\mathrm{qubit}}(t)|g\rangle\right|$, the length of the transverse Bloch vector:

```wl
qubitCoherence[rungs_] := Module[{count = creation[rungs] . annihilation[rungs], H, run},
   H = \[Epsilon] onCavity[annihilation[rungs] + creation[rungs]] + \[Chi] KroneckerProduct[Z, count];
   run = evolveODE[H, {onCavity[Sqrt[\[Gamma]disp] annihilation[rungs]]},
     densityMatrix[Flatten@KroneckerProduct[plus, coherentState[rungs, 0]]], 4.];
   Function[t, 2 Abs[reducedQubit[run[t], rungs][[1, 2]]]]];
```

Confirm the cutoff by comparing two truncations (cutoff at $n = 9$ and $n = 12$) at the same time:

```wl
reading = qubitCoherence[9];
{reading[2.], qubitCoherence[12][2.]}
```

The two cutoffs agree, so 9 is enough.

Now visualize the qubit coherence draining as the pointer states separate:

```wl
Plot[reading[t], {t, 0, 4}, Frame -> True, GridLines -> Automatic, PlotRange -> All, AspectRatio -> 1/3,
 FrameLabel -> {"time", "qubit coherence"},
 PlotLabel -> "the pointer measures away the qubit coherence"]
```

The coherence falls to near zero as the pointer states separate: a measurement that distinguishes the qubit states must destroy their superposition.

The cavity, the meter itself, tells the same story, and what it shows depends on how you read the qubit. Trace the qubit out, $\rho_{\mathrm{cavity}} = \mathrm{Tr}_{\mathrm{qubit}}[\rho]$, and the cavity is a fringeless mixture of the two pointer states as soon as they separate, since the which-way record now lives in the qubit. But keep only the records where the qubit is found in $\left|+\right\rangle$ (the basis where its coherence lives), and the conditional cavity is a cat of the two pointers, its off-diagonal coherence carrying the qubit's. A function for the cavity conditioned on finding the qubit in a state $\left|q\right\rangle$ takes the qubit matrix element and renormalizes it, $\rho_{\mathrm{cavity};q} = \langle q|\rho|q\rangle/\mathrm{Tr}\langle q|\rho|q\rangle$, where $\langle q|\rho|q\rangle$ leaves a cavity operator and its trace is the probability $\langle q|\rho_{\mathrm{qubit}}|q\rangle$ of that outcome:

```wl
cavityGiven[joint_, rungs_, q_] := With[{r = ArrayReshape[joint, {2, rungs, 2, rungs}]},
   With[{m = Sum[Conjugate[q[[a]]] q[[b]] r[[a, All, b, All]], {a, 2}, {b, 2}]}, m/Tr[m]]];
```

Run the joint master equation once, at a cutoff wide enough for the pointer blobs:

```wl
dispRungs = 16;
dispRun = evolveODE[\[Epsilon] onCavity[annihilation[dispRungs] + creation[dispRungs]] + \[Chi] KroneckerProduct[Z, creation[dispRungs] . annihilation[dispRungs]],
   {onCavity[Sqrt[\[Gamma]disp] annihilation[dispRungs]]},
   densityMatrix[Flatten@KroneckerProduct[plus, coherentState[dispRungs, 0]]], 3.];
```

In the untruncated driven-damped model, the conditional cavity has a closed two-coherent-state form. Its two diagonal blocks are the pointer coherent states $\left|\alpha_e\right\rangle$ and $\left|\alpha_g\right\rangle$; the qubit coherence carried into the cavity as the off-diagonal weight is
$$T(t) = \tfrac12\exp\left(-2i\chi\int_0^t \alpha_e\,\alpha_g^{*}\,dt'\right), \qquad \alpha_e\,\alpha_g^{*} = \frac{\epsilon^2}{\Lambda^2}(1 - e^{-\Lambda t'})^2,\quad \Lambda = i\chi + \tfrac{\gamma}{2},$$
an elementary integral. Dividing $T$ by the pointers' overlap gives the coherence coefficient $c = T/\langle \alpha_g|\alpha_e\rangle$, so
$$\rho_{\mathrm{cavity};+}(t) \;\propto\; \tfrac12|\alpha_e\rangle\langle\alpha_e| + \tfrac12|\alpha_g\rangle\langle\alpha_g| + c|\alpha_e\rangle\langle\alpha_g| + c^{*}|\alpha_g\rangle\langle\alpha_e|.$$
Build it directly:

```wl
\[CapitalLambda] = I \[Chi] + \[Gamma]disp/2;
cohT[t_] := (1/2) Exp[-2 I \[Chi] (\[Epsilon]^2/\[CapitalLambda]^2) (t + (2/\[CapitalLambda]) (Exp[-\[CapitalLambda] t] - 1) - (1/(2 \[CapitalLambda])) (Exp[-2 \[CapitalLambda] t] - 1))];
cavityCat[t_, rungs_] := With[{ke = coherentState[rungs, \[Alpha]e[t]], kg = coherentState[rungs, \[Alpha]g[t]]},
   With[{cc = cohT[t]/(Conjugate[kg] . ke)},
    With[{m = densityMatrix[ke]/2 + densityMatrix[kg]/2 + cc KroneckerProduct[ke, Conjugate[kg]] + Conjugate[cc] KroneckerProduct[kg, Conjugate[ke]]},
     m/Tr[m]]]];
```

Confirm the closed form against the simulated conditional cavity across the whole readout:

```wl
Max@Table[traceDistance[cavityCat[t, dispRungs], cavityGiven[dispRun[t], dispRungs, plus]], {t, 0.1, 3, 0.1}]
```

The maximum trace-distance discrepancy is small but nonzero because the ODE integration and Fock truncation are numerical; decreasing the solver tolerance and increasing the cutoff reduce it. Within that controlled error, the two-coherent-state formula reproduces the simulated conditional cavity. Now the question a phase-space picture cannot answer cleanly: cat or mixture? Measure the conditional cavity against both references, the pure cat of the current pointers $(\left|\alpha_e\right\rangle + \left|\alpha_g\right\rangle)$ and their even mixture, as the readout runs:

```wl
ListLinePlot[{
   Table[{t, traceDistance[cavityGiven[dispRun[t], dispRungs, plus], densityMatrix@Normalize[coherentState[dispRungs, \[Alpha]e[t]] + coherentState[dispRungs, \[Alpha]g[t]]]]}, {t, 0.1, 3, 0.05}],
   Table[{t, traceDistance[cavityGiven[dispRun[t], dispRungs, plus], densityMatrix[coherentState[dispRungs, \[Alpha]e[t]]]/2 + densityMatrix[coherentState[dispRungs, \[Alpha]g[t]]]/2]}, {t, 0.1, 3, 0.05}]},
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"distance to the pure cat", "distance to the classical mixture"},
 FrameLabel -> {"readout time", "trace distance"},
 PlotLabel -> "the conditional cavity starts a cat and decoheres into the mixture"]
```

Early the conditional cavity hugs the pure cat: the two pointers still coincide, so there is little to tell apart and little coherence to lose. As they separate its distance to the cat climbs while its distance to the mixture falls, the two crossing near $t \approx 2$ for these parameters; past there the state is closer to the classical mixture, and by the time the pointers are cleanly resolved it is close to that mixture. The cat is fleeting because distinguishability and dephasing grow together. Contrast the Haroche cat of the earlier section: there the dispersive atom imprinted the two phases *before* the cavity leaked, so separation came first and decoherence after; here they arise from the same readout process. This is dispersive qubit readout in circuit QED, where the pointer is a microwave tone reflected off a resonator.

## Part Three: A Heavy Particle, and the Two Faces of a Warm Bath

**The problem.** A particle of mass $M$ in a gentle trap $\tfrac12 M\Omega^2\hat x^2$, coupled to a warm bath (the Caldeira-Leggett model). In the high-temperature limit,
$$\dot\rho = -i[\frac{\hat p^2}{2M} + \tfrac12 M\Omega^2\hat x^2,\ \rho] - i\gamma[\hat x,\{\hat p,\rho\}] - 2\gamma M k_B T\,[\hat x,[\hat x,\rho]].$$
Friction damps the momentum at rate $2\gamma$; the double commutator $[\hat x,[\hat x,\rho]]$ decoheres a superposition of two positions a distance $\Delta x$ apart at a rate $\propto M k_B T\,\Delta x^2$, quadratic in the separation, the same law as the light-cat. For a smooth (Gaussian) blob the first and second moments of $(\hat x,\hat p)$ close into a small set of rate equations, exactly the five below.

A heavy particle in a harmonic trap $\tfrac12 M\Omega^2\hat x^2$, coupled to a warm bath (Caldeira-Leggett). Friction damps the momentum; the bath's fluctuations grow the position variance to thermal equilibrium. For a Gaussian state the first and second moments of $(\hat x,\hat p)$ close into five equations: the mean $(\langle\hat x\rangle, \langle\hat p\rangle)$ and the covariance matrix $\Sigma$ with entries $(\Sigma_{\mathrm{xx}}, \Sigma_{\mathrm{xp}}, \Sigma_{\mathrm{pp}})$. The high-temperature master equation above holds only when $k_B T > \hbar\Omega$, so fix the mass, trap frequency, friction, and a temperature well above the level spacing $\hbar\Omega = 1$:

```wl
M = 1.; \[CapitalOmega]qbm = 1.; \[Gamma]qbm = 0.15; kT = 5.; Ddiff = 4 \[Gamma]qbm M kT;
```

With $k_B T = 5$ against a level spacing of $1$, the bath is genuinely warm. Written out, those five equations are
$$\dot{\langle\hat x\rangle} = \frac{\langle\hat p\rangle}{M}, \qquad \dot{\langle\hat p\rangle} = -M\Omega^2\langle\hat x\rangle - 2\gamma\langle\hat p\rangle,$$
$$\dot\Sigma_{\mathrm{xx}} = \frac{2\Sigma_{\mathrm{xp}}}{M}, \qquad \dot\Sigma_{\mathrm{xp}} = \frac{\Sigma_{\mathrm{pp}}}{M} - M\Omega^2\Sigma_{\mathrm{xx}} - 2\gamma\Sigma_{\mathrm{xp}}, \qquad \dot\Sigma_{\mathrm{pp}} = -2M\Omega^2\Sigma_{\mathrm{xp}} - 4\gamma\Sigma_{\mathrm{pp}} + D,$$
with diffusion coefficient $D = 4\gamma M k_B T$. The initial covariance must itself be a valid quantum state, so start from the trap's ground state, the minimum-uncertainty blob with $\Sigma_{\mathrm{xx}} = \Sigma_{\mathrm{pp}} = \hbar/2 = 1/2$ (so $\det\Sigma = \hbar^2/4 = 1/4$), displaced off the origin. Integrate them from there:

```wl
particle = NDSolveValue[{
   cx'[t] == cp[t]/M, cp'[t] == -M \[CapitalOmega]qbm^2 cx[t] - 2 \[Gamma]qbm cp[t],
   vx'[t] == 2 vc[t]/M, vc'[t] == vp[t]/M - M \[CapitalOmega]qbm^2 vx[t] - 2 \[Gamma]qbm vc[t],
   vp'[t] == -2 M \[CapitalOmega]qbm^2 vc[t] - 4 \[Gamma]qbm vp[t] + Ddiff,
   cx[0] == 2., cp[0] == 0., vx[0] == 0.5, vc[0] == 0., vp[0] == 0.5},
  {cx, cp, vx, vc, vp}, {t, 0, 30}];
```

Two functions returning the mean and the covariance matrix:

```wl
center[t_] := {particle[[1]][t], particle[[2]][t]};
spread[t_] := {{particle[[3]][t], particle[[4]][t]}, {particle[[4]][t], particle[[5]][t]}};
```

Two checks. First, the mean spirals in from its displaced start to rest:

```wl
{center[0], center[25]}
```

It starts displaced and ends near the origin. Second, the position variance $\Sigma_{\mathrm{xx}}$ grows from its ground-state value $\hbar/2 = 1/2$ to the thermal value $k_B T/(M\Omega^2)$, equipartition for a harmonic trap:

```wl
{particle[[3]][0], particle[[3]][25], kT/(M \[CapitalOmega]qbm^2)}
```

It lands on the thermal value, ten times the ground-state spread. Check it a second way: the steady covariance solves the Lyapunov equation, where diffusion balances friction:

```wl
drift = {{0, 1/M}, {-M \[CapitalOmega]qbm^2, -2 \[Gamma]qbm}};
LyapunovSolve[drift, -{{0, 0}, {0, Ddiff}}]
```

The Lyapunov steady state matches the integrated one.

The same covariance bookkeeping returns under monitoring in Part Five. The Kalman filter's watched oscillator obeys $\dot\Sigma = A\Sigma + \Sigma A^{T} + D - \Sigma C^{T} C\Sigma$, the Lyapunov form above plus one term the bath cannot supply, the measurement's information gain $-\Sigma C^{T} C\Sigma$: diffusion widens the blob to thermal here, watching narrows it below the ground-state spread there.

One more check the model must pass to be believed: a real quantum state obeys the Robertson-Schrodinger uncertainty relation $\det\Sigma = \Sigma_{\mathrm{xx}}\Sigma_{\mathrm{pp}} - \Sigma_{\mathrm{xp}}^2 \ge \hbar^2/4 = 1/4$, so the margin $\det\Sigma - 1/4$ must stay non-negative for the whole run. The high-temperature Caldeira-Leggett generator is not of Lindblad form, so this is not automatic; watch it directly:

```wl
Plot[Det@spread[t] - 1/4, {t, 0, 25}, Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "uncertainty margin  det \[CapitalSigma] - 1/4"},
 PlotLabel -> "the state stays above the uncertainty floor"]
```

The margin starts at zero, where the ground state saturates the bound, and only grows for this valid high-temperature run: the blob swells above the Heisenberg floor and never dips below it. Starting instead from $\Sigma_{\mathrm{xx}}=\Sigma_{\mathrm{pp}}=0.1$ would already violate the bound at $t=0$; the moment equations can later carry that covariance into the valid region, but doing so does not retroactively make the initial state physical. Conversely, at temperatures outside the model's high-temperature regime this non-Lindblad generator can drive an initially valid covariance below the floor. The point of the plotted check is therefore local and explicit: the chosen parameters and initial state remain physical throughout the displayed interval.

Those checks read the mean and the variance as numbers; the phase-space picture puts them in one object. That object is the state's **Wigner function** $W(x,p)$, the quasiprobability the cat used in Part Two, except that for a Gaussian state no Fock-basis sum is needed: $W$ is exactly the bivariate normal set by the five moments, centered at the mean $\bar r = (\langle\hat x\rangle, \langle\hat p\rangle)$ and shaped by the covariance $\Sigma$,
$$W(x,p) = \frac{1}{2\pi\sqrt{\det\Sigma}}\,\exp\left(-\tfrac12\,(r - \bar r)^{T}\,\Sigma^{-1}\,(r - \bar r)\right), \qquad r = (x,p),$$
normalized so $\int W\,dx\,dp = 1$. Transcribe it as a function of time, reading $\bar r$ and $\Sigma$ off the solution at each $t$:

```wl
wignerGaussian[t_?NumericQ] := With[
   {r0 = center[t], iS = Inverse[spread[t]], nrm = 1/(2 Pi Sqrt[Det[spread[t]]])},
   Function[{x, p}, nrm Exp[-(({x, p} - r0) . iS . ({x, p} - r0))/2]]];
```

Trace the mean's whole path once, a faint guide to lay under the moving blob:

```wl
centerLine = ListLinePlot[Table[center[t], {t, 0, 25, 0.05}], PlotStyle -> Directive[Thin, Gray]];
```

Now set it in motion from the displaced start: the Wigner bump rides the gray path in toward the origin as friction drains the mean, and fills out to the thermal size as diffusion widens it, the red dot marking where the mean sits at each instant. Press play:

```wl
Animate[Show[
   DensityPlot[wignerGaussian[t][x, p], {x, -6, 6}, {p, -6, 6}, PlotPoints -> 60, PlotRange -> All,
    ColorFunction -> "SunsetColors", AspectRatio -> 1, Frame -> True, ImageSize -> 400,
    FrameLabel -> {"x", "p"}, PlotLabel -> "Wigner W(x,p) at t = " <> ToString[t],
    Epilog -> {Red, Point[center[t]]}], centerLine],
 {t, 0, 25, .1}, AnimationRunning -> False]
```

It stays a single Gaussian bump the whole way, never developing finer structure: the state never leaves the five-number family, the closure that let five equations describe it in the first place. Friction and diffusion are the same bath's two effects, tied by the fluctuation-dissipation relation: the friction rate $\gamma$ is set by the coupling, while the diffusion (and so the thermal spread the blob swells to) grows with temperature.

The sharper effect, and why this sits beside the cavity cat: a superposition of two positions decoheres, and a warm bath does it faster than a cold one. Here the model changes from the Caldeira-Leggett moment equations above to the thermal-oscillator Lindblad master equation (a decay channel $\hat a$ at rate $\gamma(n_T+1)$ and an excitation channel $\hat a^\dagger$ at rate $\gamma n_T$, the same warm-bath form used later for thermalization), integrated in a truncated Fock basis. Build the same coherent-state cat and evolve it under that bath. Set the oscillator cutoff and the bath's rate:

```wl
twoRung = 22; \[Gamma]warm = 0.5;
```

The surviving coherence is the trace distance between the coherently evolved cat and the classical mixture its two lobes would form, both carried through the same thermal bath $\mathcal{L}\rho = \gamma(n_T+1)\,\mathcal{D}[\hat a]\rho + \gamma n_T\,\mathcal{D}[\hat a^\dagger]\rho$:
$$\mathrm{warmCoherence}(n_T,\alpha,t) = \tfrac12\,\mathrm{Tr}\left|\,\rho_{\mathrm{cat}}(t) - \tfrac12[\rho_{+\alpha}(t) + \rho_{-\alpha}(t)]\,\right|,\quad \rho_{\mathrm{cat}}(t) = e^{\mathcal{L}t}[|\mathrm{cat}\rangle\langle\mathrm{cat}|],\ \ \rho_{\pm\alpha}(t) = e^{\mathcal{L}t}[|{\pm\alpha}\rangle\langle{\pm\alpha}|],$$
with $\alpha = \text{gap}$. Here `gap` is the coherent amplitude of either lobe and therefore half their phase-space amplitude separation $2\left|\alpha\right|$. It is largest while the cat is still distinct from that mixture and falls to zero as the cross-terms decohere (both thermalize to the same steady state). The bath's $\hat a^\dagger$ channel heats a coherent state, so the two lobes will not stay coherent and are evolved rather than written down; compute it over a list of delays, one run each for the cat and its lobes:

```wl
warmCoherence[nT_, gap_, delays_] := With[
   {fall = Sqrt[\[Gamma]warm (nT + 1)] annihilation[twoRung], climb = Sqrt[\[Gamma]warm nT] creation[twoRung],
    pair = densityMatrix@Normalize[coherentState[twoRung, gap] + coherentState[twoRung, -gap]],
    blobA = densityMatrix[coherentState[twoRung, gap]], blobB = densityMatrix[coherentState[twoRung, -gap]]},
   With[{runs = evolveODE[0 IdentityMatrix[twoRung], {fall, climb}, #, Max[delays]] & /@ {pair, blobA, blobB}},
    Table[Total@SingularValueList[runs[[1]][t] - (runs[[2]][t] + runs[[3]][t])/2]/2, {t, delays}]]];
```

Now visualize the temperature dependence by comparing a cold and a warm bath at the same separation:

```wl
warmTimes = Range[0, 1.5, 0.06];
ListLinePlot[{Transpose[{warmTimes, warmCoherence[0.15, 1.5, warmTimes]}], Transpose[{warmTimes, warmCoherence[1.5, 1.5, warmTimes]}]},
 PlotMarkers -> Automatic, Frame -> True, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"cold bath", "warm bath"}, FrameLabel -> {"time", "surviving coherence"},
 PlotLabel -> "a warm bath decoheres the cat faster than a cold one"]
```

The warm bath decoheres the cat far faster: the rate grows with temperature. As with the cavity cat, the rate also grows with the separation: together the decoherence rate scales as $(2\bar n + 1)$ times the *square* of the separation. The $+1$ is the zero-point (vacuum) term that decoheres even a cold cavity cat; temperature adds to it, so heating speeds decoherence but is not what starts it. For macroscopic separations the rate is enormous at any temperature: the fuller answer to why superpositions of large objects become so difficult to preserve as their separation grows. Sweep the cat size continuously at the fixed short delay, and the two effects trade off:

```wl
ListPlot[With[{sep = Range[0, 3, 0.05]}, {#, warmCoherence[1.5, #, {0.15}][[1]]} & /@ sep],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"coherent amplitude |\[Alpha]| (half-separation)", "surviving coherence at fixed short delay"},
 PlotLabel -> "surviving coherence peaks at a modest cat size, then falls"]
```

The curve is not monotonic: a bigger cat starts with more coherence but decoheres faster, at a rate set by the separation squared, so at a fixed short delay the two trade off and the coherence peaks at a modest size before the quadratic decoherence overruns it. (Stop at separation 3: past there the Fock cutoff clips the coherent states.)

## Part Four: One Leak, Three Ways to Watch It

So far we mostly averaged over the record. Now we keep it. The key fact: the same master equation can be **unravelled** into trajectories in genuinely different ways, each giving different single histories that all average to the same master equation. We watch one driven, emitting atom three ways. Fix the atom, its drive, its decay channel, and the time grid:

```wl
\[CapitalOmega]atom = 3.0; \[Gamma]atom = 1.0; Hdrive = (\[CapitalOmega]atom/2) X; cAtom = Sqrt[\[Gamma]atom] lower;
dtAtom = 0.01; when = dtAtom Range[0, 600];
```

Given a density matrix, define a function returning the excited-state population:

```wl
excitedPop[rho_] := (1 + blochVector[rho][[3]])/2;
```

Integrate the master equation for the population every unravelling must average back to:

```wl
smoothAtom = evolveODE[Hdrive, {cAtom}, densityMatrix[excited], 6.0];
averageChance = excitedPop[smoothAtom[#]] & /@ when;
```

### Counting the Flashes, One at a Time

**The problem.** Any leak $\dot\rho = -i[\hat H,\rho] + \mathcal{D}[\hat c]\rho$ can be *unravelled* by watching its output, each way of detecting the emitted light giving one unravelling into pure-state histories. Counting photons one click at a time gives the jump trajectory
$$d\left|\psi\right\rangle = \left[dN\left(\tfrac{\hat c}{\sqrt{\langle\hat c^\dagger\hat c\rangle}} - \mathbb{1}\right) + dt\left(\tfrac{\langle\hat c^\dagger\hat c\rangle}{2} - \tfrac{\hat c^\dagger\hat c}{2} - i\hat H\right)\right]\left|\psi\right\rangle,$$
where the click $dN \in \{0,1\}$ arrives at rate $\langle\hat c^\dagger\hat c\rangle\,dt$: between clicks the norm decays under the non-Hermitian $\hat H - \tfrac{i}{2}\hat c^\dagger\hat c$, and at a click the state jumps by $\hat c$. In continuous time, averaging $\left|\psi\right\rangle\left\langle\psi\right|$ over histories returns the master equation. The finite-step simulation below also carries Bernoulli time-discretization error, so its comparison contains both that bias and ensemble scatter.

The first way, **photon counting**: detect each emitted photon. Between clicks the state evolves under the non-Hermitian $\hat H - \tfrac{i}{2}\hat c^\dagger\hat c$; at each click it jumps by $\hat c$ to the ground state. A single run is one sweep that does both jobs at once: `FoldList` advances and keeps the conditioned state at every step (a jump with probability $\langle\hat c^\dagger\hat c\rangle\,dt$, else the no-click drift), while `Sow` writes down the time whenever the jump fires. It returns the inferred state history and the raw click record together, because in a real experiment they are one run, not two. Everything that follows, one trajectory, the ensemble average, and the waiting-time statistic, reads off realizations of this single generator:

```wl
countingRun[seed_, steps_] := BlockRandom[SeedRandom[seed];
   With[{noClick = MatrixExp[-I (Hdrive - I/2 ConjugateTranspose[cAtom] . cAtom) dtAtom]},
    With[{r = Reap@FoldList[
         Function[{psi, n},
          If[RandomReal[] < dtAtom Re[Conjugate[cAtom . psi] . (cAtom . psi)],
            (Sow[n dtAtom]; Normalize[cAtom . psi]),
            Normalize[noClick . psi]]],
         excited, Range[steps]]},
     <|"states" -> First[r], "clicks" -> Flatten[Last[r]]|>]]];
```

The same population, read from a state vector:

```wl
excitedPopK[psi_] := (1 + Re[Conjugate[psi] . Z . psi])/2;
```

Record one counted history over the displayed window, and read its two parts, the conditioned states and the click times, straight off the one run:

```wl
oneRun = countingRun[3, 600]; oneCount = oneRun["states"]; oneClicks = oneRun["clicks"];
```

Now visualize one counted trajectory's excited-state population against the master equation: the trajectory drives up and down and plunges at each emission, while the average remains smooth:

```wl
ListLinePlot[{Transpose[{when, excitedPopK /@ oneCount}], Transpose[{when, averageChance}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"one trajectory (jumps)", "master equation"}, PlotRange -> All,
 FrameLabel -> {"time", "excited-state population"}, PlotLabel -> "photon counting: one trajectory jumps, the average is smooth"]
```

The plunges are the emitted photons, one at a time. No single run is smooth; smoothness is in the average.

The plot shows the inferred *state*; the same run also carries the *record*, the list of click times a photon counter actually writes down, already returned as `oneRun["clicks"]` and bound to `oneClicks` above. Every plunge is one entry in it, because the jump that plunges the state is the very step that sowed the time. Count the clicks in this record:

```wl
Length[oneClicks]
```

A handful of photons over the run. Now put the record and the state on one time axis with two vertical scales: the conditioned population on the left, each click marked by a red line, and the cumulative count $N(t)$ on the right, a staircase that rises by one at every click:

```wl
Legended[
 Labeled[
  ListLinePlot[{Transpose[{when, excitedPopK /@ oneCount}],
     Transpose[{Riffle[oneClicks, oneClicks], Riffle[Range@Length@oneClicks - 1, Range@Length@oneClicks]}]},
   MultiaxisArrangement -> All, PlotRange -> All, Frame -> True,
   GridLines -> {{#, Red} & /@ oneClicks, None}],
  {Rotate["excited population", 90 Degree], Rotate["cumulative count N(t)", 270 Degree], "time",
   "clicks (red lines) fall at the population plunges; the staircase counts them"}, {Left, Right, Bottom, Top}],
 LineLegend[{ColorData[97, 1], ColorData[97, 2]}, {"excited population", "cumulative count N(t)"}]]
```

Every red line sits on a plunge: the click is the emission, and the plunge is the state collapsing to the ground state the instant the photon is registered. The staircase (right scale) is the raw experimental data; the population trajectory (left scale) is what the theorist infers from it.

Now run an ensemble of these same trajectories, carried well past the display window so each click record is long enough for the timing statistic below; our single run above is seed 3 of this same generator, the ensemble's third member seen over just the plotted window:

```wl
nCount = 200; nWaitSteps = 4000;
countRuns = countingRun[#, nWaitSteps] & /@ Range[nCount];
```

Averaging the conditioned states across the ensemble must recover the master equation. Compare the ensemble-mean excited population, over the displayed window, against the integrated reference:

```wl
Max@Abs[(excitedPop /@ Mean[Table[densityMatrix /@ Take[countRuns[[k]]["states"], Length[when]], {k, nCount}]]) - averageChance]
```

At fixed `dtAtom`, the sampling component of the gap shrinks as $1/\sqrt N$; the residual time-step bias requires refinement of `dtAtom`. Together those limits recover the master equation.

The click record carries a statistic no smooth curve can: the distribution of *waiting times* between clicks. A click leaves the atom in $\left|g\right\rangle$, so the time $\tau$ to the next click has the density
$$w(\tau) = \langle \psi(\tau),\hat c^\dagger\hat c\,\psi(\tau)\rangle, \qquad \left|\psi(\tau)\right\rangle = e^{-i\hat H_{\mathrm{eff}}\tau}\left|g\right\rangle,\quad \hat H_{\mathrm{eff}} = \hat H - \tfrac{i}{2}\hat c^\dagger\hat c,$$
the unnormalized no-click state evolved from the ground state. Since it starts in $\left|g\right\rangle$, $w(0) = \gamma\left|\langle e|g\rangle\right|^2 = 0$: the atom cannot emit twice at once, it must be re-driven to the excited state first, so very short gaps are forbidden. This is antibunching. Read the density directly off the no-click evolution:

```wl
waitDensity[t_] := With[{psi = MatrixExp[-I (Hdrive - I/2 ConjugateTranspose[cAtom] . cAtom) t] . ground},
   Re[Conjugate[psi] . ConjugateTranspose[cAtom] . cAtom . psi]];
```

No separate simulation is needed. Because every click resets the atom exactly to $\left|g\right\rangle$, the gap between two consecutive clicks in any trajectory is already a draw from $w(\tau)$: it begins where $w$ begins. Pool those gaps straight from the ensemble's click records, one set of differences per run (the leading interval, which starts from the excited state rather than from $\left|g\right\rangle$, drops out on its own, since a difference needs two clicks):

```wl
waitTimes = Flatten[Differences[#["clicks"]] & /@ countRuns];
```

Overlay their histogram on the closed form:

```wl
Show[Histogram[waitTimes, {0, 6, 0.2}, "PDF", ChartStyle -> ColorData[97, 1]],
 Plot[waitDensity[t], {t, 0, 6}, PlotStyle -> Directive[Thick, Red]],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"waiting time between clicks", "density"},
 PlotLabel -> "the trajectory's own inter-click gaps land on the closed-form density"]
```

The bars land on the closed-form density: it vanishes at zero (no two clicks at once), rises to a first hump near half a Rabi period, then oscillates and decays, the re-driven atom's Rabi cycle written into the timing of the clicks. This antibunching is the hallmark of single-photon (nonclassical) light, and it lives entirely in the click timing, invisible in the smooth population.

### Homodyne: Reading One Quadrature

**The problem.** Mixing the output with a strong reference beam before detecting it (homodyne) turns the clicks into a continuous noisy current and gives a *diffusive* stochastic master equation
$$d\rho = -i[\hat H,\rho]\,dt + \mathcal{D}[\hat c]\rho\,dt + \mathcal{H}[\hat c]\rho\,dW,$$
with a real Wiener increment $dW$ ($dW^2 = dt$). The conditioned state now diffuses continuously instead of jumping: starting pure in the ideal example here it stays pure, and a mixed initial state would purify. The point of this section is the bookkeeping: what a single run gives us, which quantities are *observed* and which are *inferred*, and how every homodyne object is one of those two or is built from them.

The second way, **homodyne** detection: mix the output with a strong reference beam and measure one quadrature as a continuous noisy current. Run one homodyne trajectory of the same atom:

```wl
smoothHistory = trajectory[densityMatrix[excited], Hdrive, {cAtom}, {1.}, {}, dtAtom, 6.0, 3];
```

`trajectory` hands back exactly two things, and everything else on this page is one of them or is computed from them. The **conditional state** $\rho_c(t)$ at each time node (`smoothHistory["states"]`) is the observer's running estimate of the atom given the record so far. The **record increment** $dJ$ at each step (`smoothHistory["record"]`) is the raw output of the detector. The model tying them is
$$dJ = \sqrt\gamma\,\langle\hat\sigma_x\rangle\,dt + dW,$$
the signal read from $\rho_c$, the noise $dW$ a Wiener increment. Hold the two apart: $dJ$ is **observed**, and $\rho_c$ is **inferred**, reconstructed from $dJ$ by the filter (the positivity step from the toolkit).

Take the conditional state first. Its excited population is one smooth diffusing history, no jumps:

```wl
ListLinePlot[{Transpose[{when, excitedPop /@ smoothHistory["states"]}], Transpose[{when, averageChance}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"one conditional history", "master equation"}, PlotRange -> All,
 FrameLabel -> {"time", "excited-state population"}, PlotLabel -> "the conditional state diffuses smoothly, no jumps"]
```

Smooth diffusion, not the plunges of a counted trajectory, yet averaging many conditional states recovers the same master equation. Hand the ensemble to the toolkit's `gapDiagnosis`, which averages the trajectories, compares the mean to the master equation, and splits sampling scatter from step bias in one call:

```wl
gapDiagnosis[densityMatrix[excited], Hdrive, {cAtom}, {1.}, {}, dtAtom, 6.0, 150]
```

The ensemble gap sits inside the sampling scatter it reports, so the systematic remainder is zero and the per-step bias beside it is negligible: at this step and this many trajectories the record-conditioned states average to the unconditioned master equation within Monte-Carlo noise, and the verdict says it is more trajectories, not a finer step, that would tighten the agreement.

Now the record itself, and the four quantities read from it, each labeled by how it is obtained and whether it is observed or inferred.

**The record increment $dJ$** (observed). *How obtained:* the detector's output each step, here generated together with the state. *Status:* observed, the one thing that enters from the world; everything else is downstream of it. Plot it:

```wl
dJ = Flatten@smoothHistory["record"];
```

```wl
ListLinePlot[Transpose[{Most@when, dJ}],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "record increment  dJ"}, PlotLabel -> "dJ, the observed record: almost pure noise per step"]
```

It looks like noise because it is almost all noise: per step the signal $\sqrt\gamma\langle\hat\sigma_x\rangle\,dt \approx 0.005$ sits under a $dW$ of size $\sqrt{dt} = 0.1$, twenty times larger. A single increment carries almost no information about the atom.

**The integrated record $J(t) = \int_0^t dJ$** (observed). *How obtained:* sum the observed increments. *Status:* observed, the same data accumulated. It is built from three running arrays. First the inferred signal $\sqrt\gamma\langle\hat\sigma_x\rangle$, read off each conditional state (it gets its own panel below):

```wl
signal = Sqrt[\[Gamma]atom] (blochVector[#][[1]] & /@ Most[smoothHistory["states"]]);
```

The integrated record itself, the running sum of the observed increments:

```wl
Jt = Prepend[Accumulate[dJ], 0.];
```

and its drift, the inferred signal accumulated the same way:

```wl
driftJ = Prepend[Accumulate[signal dtAtom], 0.];
```

Plot the record against its drift:

```wl
ListLinePlot[{Transpose[{when, Jt}], Transpose[{when, driftJ}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"integrated record  J(t)", "inferred drift  \[Integral]signal dt"},
 FrameLabel -> {"time", "integrated record  J(t)"}, PlotLabel -> "the observed record integrates to a drifting path; its drift is the signal"]
```

Integrating is itself a filter: the signal accumulates coherently into a drift while the noise only grows as $\sqrt t$, so in $J(t)$ the signal is comparable to the noise instead of twenty times below it. This is the honest, window-free picture of the data. (The `signal` used for the drift is inferred, and gets its own panel below.)

**The current $dJ/dt$** (observed, singular). *How obtained:* divide the observed increment by $dt$. *Status:* observed but not a proper quantity: $dJ/dt = \text{signal} + dW/dt$, and $dW/dt$ is white noise of variance $1/dt$, which diverges as $dt\to0$. On the grid it is a $\pm 1/\sqrt{dt}$ band around a $\pm 1$ signal:

```wl
ListLinePlot[Transpose[{Most@when, dJ/dtAtom}],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "current  dJ/dt"}, PlotLabel -> "the current dJ/dt: a white-noise band with no dt\[Rule]0 limit"]
```

There is no current *curve*. The instantaneous homodyne current is a distribution, not a function, which is why any smooth current needs an arbitrary averaging window and why $J(t)$ above is the cleaner observed object.

**The signal $\sqrt\gamma\langle\hat\sigma_x\rangle$** (inferred). *How obtained:* read off the conditional state. *Status:* inferred, the conditional mean of the current, never measured directly. This is the clean thing buried in the record:

```wl
ListLinePlot[Transpose[{Most@when, signal}],
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "signal  Sqrt[gamma] <sigma_x>"},
 PlotLabel -> "the signal, inferred from the conditional state, not measured"]
```

No detector emits this curve. It is what the filter reconstructs, smooth because the conditional state is smooth, and it is the conditional drift about which the observed $J(t)$ fluctuates.

**The innovation $dW$** (inferred). *How obtained:* $dW = dJ - \text{signal}\,dt$, the observed increment minus the inferred prediction. *Status:* inferred, a residual, not a second measurement. Form the residual:

```wl
dW = dJ - signal dtAtom;
```

and plot it:

```wl
ListLinePlot[Transpose[{Most@when, dW}],
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "innovation  dW"}, PlotLabel -> "the innovation dW = dJ - signal dt: the residual"]
```

If the filter is optimal, the residual is zero-mean white noise. Check its autocorrelation, with the $\pm 2/\sqrt N$ band white noise would sit inside:

```wl
With[{innov = dW/Sqrt[dtAtom], ci = 2/Sqrt[Length[dW]]},
 ListPlot[CorrelationFunction[innov, {1, 10}], Filling -> Axis, PlotRange -> {-0.2, 0.2},
  Frame -> True, GridLines -> {None, {{ci, Directive[Gray, Dashed]}, {-ci, Directive[Gray, Dashed]}}},
  PlotMarkers -> Automatic, AspectRatio -> 1/3,
  FrameLabel -> {"lag", "innovation autocorrelation"}, PlotLabel -> "the innovation is white within confidence limits"]]
```

Every nonzero lag sits in the band: the filter has extracted every predictable part of the record, and what remains is irreducible measurement noise. This whiteness is the diffusive counterpart of the counting record's antibunched clicks, and the same optimal-filter signature the Kalman tracker will show.

So the section is one ledger. A single observed record $dJ$ (equivalently its integral $J$) enters from the world; from it the conditional state, the signal, and the innovation $dW$ are all computed. Only $dJ$ is data. Everything that looks like the atom, the smooth population, the signal curve, is inferred from it.

### Heterodyne: Reading Both Quadratures

**The problem.** Reading *both* shadows at once (heterodyne) splits the one monitored leak into two half-strength homodyne channels a quarter-turn apart, $\hat c_I = \hat c/\sqrt2$ and $\hat c_Q = i\hat c/\sqrt2$, each read as its own real current:
$$d\rho = -i[\hat H,\rho]\,dt + \mathcal{D}[\hat c]\rho\,dt + \mathcal{H}[\hat c_I]\rho\,dW_I + \mathcal{H}[\hat c_Q]\rho\,dW_Q,\qquad dW_I^2 = dW_Q^2 = dt,\quad dW_I\,dW_Q = 0.$$
It is one point in a whole *general-dyne* family running from one sharp shadow (homodyne) to two blurry ones, every member averaging to the identical master equation while giving genuinely different single histories. Show that the two channels reproduce the one leak, then watch.

The third way, **heterodyne**: measure both quadratures at once. It is equivalent to two half-strength homodyne channels a quarter-turn apart, $\hat c/\sqrt2$ and $i\hat c/\sqrt2$. First prove those two dissipators sum to the original, for any $\hat c$ and $\rho$:

```wl
With[{c = Array[Subscript[\[FormalC], ##] &, {2, 2}],
   r = Array[Subscript[\[FormalR], ##] &, {2, 2}]},
  dissipator[c/Sqrt[2], r] + dissipator[I c/Sqrt[2], r] ==
   dissipator[c, r]] // FullSimplify
```

It returns True, with no assumptions: the two channels reproduce the same dissipator, so the master equation is unchanged, and as in the counting and homodyne runs an ensemble of these trajectories averages back to it.

The trajectory explores differently, so hold apart what heterodyne *observes* from what it *infers*, exactly as with homodyne. `trajectory` still returns two objects, the conditional state $\rho_c$ (inferred) and the record (observed), but the record now has *two* columns, one real current per quadrature:
$$dJ_I = \sqrt{\gamma/2}\,\langle\hat\sigma_x\rangle\,dt + dW_I,\qquad dJ_Q = \sqrt{\gamma/2}\,\langle\hat\sigma_y\rangle\,dt + dW_Q.$$
Equivalently the two real records are the real and imaginary parts of a single *complex* record $dJ = dJ_I + i\,dJ_Q$, whose increment satisfies $dJ^2 = 0$ (the two shot noises are independent) and $dJ\,dJ^* = 2dt$ (each quadrature contributes its own variance $dt$): two independent real noises carried in one complex number. The single new fact over homodyne is that both field quadratures are recorded at once, so the two records carry the two transverse signals $\langle\hat\sigma_x\rangle$ and $\langle\hat\sigma_y\rangle$ directly, each at half the power an optimally phased homodyne puts in its one. Run one history and keep its record:

```wl
hetOne = trajectory[densityMatrix[excited], Hdrive, {cAtom/Sqrt[2], I cAtom/Sqrt[2]}, {1., 1.}, {}, dtAtom, 6.0, 3];
```

**The conditional state $\rho_c$** (inferred). *How obtained:* the filter's running estimate of the atom given the record so far. *Status:* inferred. Where homodyne resolved one transverse coordinate, heterodyne resolves both; plot the run's $\langle\hat\sigma_x\rangle$ and $\langle\hat\sigma_y\rangle$ together:

```wl
ListLinePlot[{Transpose[{when, blochVector[#][[1]] & /@ hetOne["states"]}], Transpose[{when, blochVector[#][[2]] & /@ hetOne["states"]}]},
 PlotLegends -> {"\[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]", "\[LeftAngleBracket]\[Sigma]y\[RightAngleBracket]"},
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "conditional transverse Bloch vector"}, PlotLabel -> "one heterodyne run resolves both quadratures"]
```

Both components diffuse smoothly, no jumps. The conditional state is a full Bloch vector under homodyne too, tracked by the filter from its one record; what heterodyne adds is a record that measures both transverse signals directly, where homodyne's single record measured one and left the other to the conditional dynamics.

**The record increments $dJ_I, dJ_Q$** (observed). *How obtained:* the detector's two outputs each step. *Status:* observed, the only data that enters from the world; everything else is downstream. Read the two columns:

```wl
dJI = hetOne["record"][[All, 1]];
dJQ = hetOne["record"][[All, 2]];
```

```wl
ListLinePlot[{Transpose[{Most@when, dJI}], Transpose[{Most@when, dJQ}]},
 PlotLegends -> {"I record", "Q record"},
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "record increments"}, PlotLabel -> "both quadrature records: almost pure noise per step"]
```

Each is nearly pure noise, and for a sharper reason than homodyne's: each channel carries only $\sqrt{\gamma/2}$ of signal, half the power homodyne put into its one quadrature, against the same shot noise.

**The integrated record $J(t) = \int_0^t dJ$** (observed). *How obtained:* accumulate each column. *Status:* observed, the same data summed. A complex record integrates to a path in the $IQ$ plane, an observed walk that drifts along the inferred signal.

That signal is the per-channel coupling $\sqrt{\gamma/2}$ times each conditional quadrature, read off $\rho_c$ rather than measured, one channel for $\langle\hat\sigma_x\rangle$ and one for $\langle\hat\sigma_y\rangle$:

```wl
hetSigI = Sqrt[\[Gamma]atom/2] (blochVector[#][[1]] & /@ Most[hetOne["states"]]);
hetSigQ = Sqrt[\[Gamma]atom/2] (blochVector[#][[2]] & /@ Most[hetOne["states"]]);
```

The observed walk is the running sum of the two record columns, each starting from the origin:

```wl
JI = Prepend[Accumulate[dJI], 0.]; JQ = Prepend[Accumulate[dJQ], 0.];
```

The drift is that same running sum on the inferred signal, the smooth curve the walk fluctuates around:

```wl
driftI = Prepend[Accumulate[hetSigI dtAtom], 0.]; driftQ = Prepend[Accumulate[hetSigQ dtAtom], 0.];
```

Now trace the observed walk against that inferred drift:

```wl
ListLinePlot[{Transpose[{JI, JQ}], Transpose[{driftI, driftQ}]},
 PlotLegends -> {"integrated record", "inferred drift"},
 Frame -> True, GridLines -> Automatic, AspectRatio -> 1,
 FrameLabel -> {"I", "Q"}, PlotLabel -> "the observed record is a drifting walk in the IQ plane"]
```

The jagged path is the raw data; the smooth curve it drifts along is inferred. Integrating both quadratures at once is a 2D filter: each signal accumulates coherently into the drift while its noise only spreads as $\sqrt t$, so the walk's net displacement carries the record's information where a single increment held almost none.

**The innovations $dW_I, dW_Q$** (inferred). *How obtained:* $dW_I = dJ_I - \text{signal}_I\,dt$ and $dW_Q = dJ_Q - \text{signal}_Q\,dt$, each observed increment minus its inferred prediction. *Status:* inferred residuals. If the filter is optimal each is zero-mean white noise, and because $dJ^2 = 0$ the two are mutually independent. Check the whiteness of each channel first, against the $\pm 2/\sqrt N$ band:

```wl
dWI = dJI - hetSigI dtAtom;
dWQ = dJQ - hetSigQ dtAtom;
```

```wl
With[{ci = 2/Sqrt[Length[dWI]]},
 ListPlot[{CorrelationFunction[dWI/Sqrt[dtAtom], {1, 10}], CorrelationFunction[dWQ/Sqrt[dtAtom], {1, 10}]}, AspectRatio -> 1/3,
  PlotRange -> {-0.2, 0.2}, PlotMarkers -> Automatic, PlotLegends -> {"I channel", "Q channel"},
  Frame -> True, GridLines -> {None, {{ci, Directive[Gray, Dashed]}, {-ci, Directive[Gray, Dashed]}}},
  FrameLabel -> {"lag", "innovation autocorrelation"}, PlotLabel -> "each innovation channel is white"]]
```

Both correlograms hug zero within the band, the signature of white residuals: the filter has extracted every predictable part of each record.

Now the heterodyne-specific check with no homodyne counterpart, whether the two residuals are independent of each other. Read off the 2-by-2 correlation matrix of the innovation channels:

```wl
Correlation[Transpose[{dWI, dWQ}]] // MatrixForm
```

The off-diagonal sits near zero, within the same $\pm 2/\sqrt N$ band: the two residuals are not just each white but mutually independent, the simulation realizing the $dJ_I\,dJ_Q = 0$ half of $dJ^2 = 0$. Reading both conjugate quadratures buys two independent white innovations, one per shadow.

Keeping the record splits purity into two quantities that need not agree: the purity of each conditioned state, and the purity of the state you get by averaging them. Run ninety homodyne histories of the same atom:

```wl
homRuns = Table[trajectory[densityMatrix[excited], Hdrive, {cAtom}, {1.}, {}, dtAtom, 6.0, s]["states"], {s, 90}];
```

Follow both purities in time:

```wl
ListLinePlot[{Transpose[{when, Mean[purity /@ #] & /@ Transpose[homRuns]}],
   Transpose[{when, purity[Mean[#]] & /@ Transpose[homRuns]}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"mean purity of conditioned states", "purity of the mean state"},
 PlotRange -> All, FrameLabel -> {"time", "purity"},
 PlotLabel -> "every run stays pure; only the record-blind average mixes"]
```

The conditioned purity holds at one for the whole run; the averaged state's purity falls and settles at the mixed steady value. Decoherence sits entirely in the second curve: it is not an event in any watched atom's history but the cost of discarding the record.

The two watched channels average to the same master equation the rest of the catalog shares. Hand the heterodyne ensemble to `gapDiagnosis`, its two integration axes the pair the diagnostic was built to handle, splitting sampling scatter from step bias in one call:

```wl
gapDiagnosis[densityMatrix[excited], Hdrive, {cAtom/Sqrt[2], I cAtom/Sqrt[2]}, {1., 1.}, {}, dtAtom, 6.0, 150]
```

Counting, homodyne, heterodyne: three different measurements of one atom, three different trajectory ensembles, one shared master equation. The measurement you choose sets what you learn moment to moment; the unconditioned dynamics is the same. Between homodyne and heterodyne lies the continuous general-dyne family, all on this same average.

## Part Five: What Watching Does, and How to Act on It

Keeping the record does more than reveal the state: measurement reshapes the dynamics, and the record can be fed back to control the system. Five examples follow: a watched-frozen qubit, a charge qubit read by a detector, measurement-induced collapse, a Kalman filter tracking an oscillator, and measurement feedback (all qubits except the oscillator). First, two functions returning the $z$ and $x$ components of a list of states:

```wl
zOf[states_] := blochVector[#][[3]] & /@ states; xOf[states_] := blochVector[#][[1]] & /@ states;
```

### Watching Harder Freezes It

**The problem.** A driven qubit continuously measured along $\hat\sigma_z$, the quantum Zeno setup. The conditioned state obeys the stochastic master equation
$$d\rho = -i\,\tfrac{\Omega}{2}[\hat\sigma_x,\rho]\,dt - \tfrac{k}{2}\,[\hat\sigma_z,[\hat\sigma_z,\rho]]\,dt + \sqrt{k}\,(\hat\sigma_z\rho + \rho\hat\sigma_z - 2\langle\hat\sigma_z\rangle\rho)\,dW,$$
with drive $\Omega$ and measurement strength $k$. When $k < \Omega$ the qubit swings almost freely; when $k > \Omega$ the watching pins it near a pole and $z(t)$ becomes a random telegraph, the "watched pot" that never boils. Watch a gently-watched and a hard-watched history.

A driven qubit measured along $\hat\sigma_z$. Measure weakly ($k < \Omega$) and it swings almost freely; measure strongly ($k > \Omega$) and it is pinned near a $\hat\sigma_z$ eigenstate with occasional jumps, the quantum Zeno effect. Set the drive and the fine time grid:

```wl
\[CapitalOmega]zeno = 2.0; fastZeno = 0.004; clockZeno = fastZeno Range[0, 3000];
```

Run one weakly and one strongly measured trajectory, same start and seed. The gently watched history takes measurement strength $k = 0.15$, well below the drive $\Omega = 2$:

```wl
gentle = trajectory[densityMatrix[excited], \[CapitalOmega]zeno/2 X, {Sqrt[0.15] Z}, {1.}, {}, fastZeno, 12., 7];
```

and the hard-watched one takes $k = 6$, well above it:

```wl
hard = trajectory[densityMatrix[excited], \[CapitalOmega]zeno/2 X, {Sqrt[6.] Z}, {1.}, {}, fastZeno, 12., 7];
```

Compare the magnitude of the time-averaged $\langle\hat\sigma_z\rangle$:

```wl
{Abs@Mean[zOf@gentle["states"]], Abs@Mean[zOf@hard["states"]]}
```

The weakly measured $\langle\hat\sigma_z\rangle$ averages to near zero (it swings freely); the strongly measured one keeps a large mean (pinned).

Now visualize the weakly and strongly measured trajectories together:

```wl
ListLinePlot[{Transpose[{clockZeno, zOf@gentle["states"]}], Transpose[{clockZeno, zOf@hard["states"]}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"weak measurement (free)", "strong measurement (pinned)"}, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
 PlotLabel -> "strong measurement freezes the qubit (Zeno)"]
```

The strongly measured qubit clings to an eigenstate for long stretches, then jumps. Strong measurement freezes the dynamics, the Zeno effect: the more often you measure, the less the state can evolve. The jumps are the drive occasionally winning.

The single trace shows the freezing but cannot pin its rate. That rate is exact and sits in the master equation's own spectrum, and following it across the measurement strength tells a two-stage story: watching harder first kills the coherent Rabi oscillation, then freezes the relaxation left behind. Both stages are eigenvalues of the Liouvillian, and at measurement strength $k$ against drive $\Omega$ they come out in closed form:

```wl
Eigenvalues@FullSimplify[liouvillian[(\[CapitalOmega]/2) X, {Sqrt[k] Z}, 2], \[CapitalOmega] > 0 && k > 0]
```

Four eigenvalues: $0$, the steady state; $-2k$, the transverse dephasing, off to the side and taking no part in the freezing; and the Rabi pair $-k \pm \sqrt{k^2 - \Omega^2}$. Below the critical strength $k < \Omega$ the radical is imaginary, $-k \pm i\sqrt{\Omega^2 - k^2}$, a damped oscillation the measurement only slows; above it, $k > \Omega$, the radical turns real, two overdamped rates; the two meet and coalesce at the exceptional point $k = \Omega$.

For the sweep, build them numerically at a chosen fraction $r = k/\Omega$ of the drive, keeping the three that carry a rate: the steady-state $0$ is dropped, so the slowest decay reads as a genuine rate rather than that zero, and the exceptional point falls at $r = 1$:

```wl
zenoEigs[ratio_] := Select[Eigenvalues[liouvillian[(\[CapitalOmega]zeno/2) X, {Sqrt[ratio \[CapitalOmega]zeno] Z}, 2]], Abs[#] > 10^-9 &];
```

The oscillation lives in the pair's imaginary part and the relaxation in the least-negative real part, so read one frequency and one rate off the spectrum at each ratio:

```wl
zenoFreq[ratio_] := Max[Abs[Im[zenoEigs[ratio]]]];
zenoRate[ratio_] := -Max[Re[zenoEigs[ratio]]];
```

Follow both across the sweep, the oscillation frequency and the slowest relaxation rate against the ratio $k/\Omega$:

```wl
ListLinePlot[
 {Table[{r, zenoFreq[r]}, {r, 0.01, 3., 0.01}], Table[{r, zenoRate[r]}, {r, 0.01, 3., 0.01}]},
 Frame -> True, GridLines -> {{{1, Directive[Gray, Dashed]}}, Automatic},
 PlotLegends -> Placed[{"Rabi oscillation  |Im \[Lambda]|", "relaxation rate  |Re \[Lambda]|"}, {Right, Top}],
 FrameLabel -> {"k/\[CapitalOmega] ratio", "Liouvillian eigenvalue  (1/time)"},
 PlotLabel -> "watching harder: the oscillation dies at k = \[CapitalOmega], the relaxation rate peaks there"]
```

The oscillation frequency slides to zero at $k = \Omega$ and stays there, the coherent swinging gone; at that same strength the relaxation rate reaches its peak and turns over, so watching harder past $\Omega$ slows the decay rather than quickening it. The pair has met on the real axis at an exceptional point, where the eigenvalues and their eigenvectors coincide and the oscillation can no longer be told apart from the decay.

Past that peak the surviving rate falls as one over the strength. Read it across a wide range on log-log, with the $\Omega^2/(2k)$ law drawn on the strong side:

```wl
zenoKs = \[CapitalOmega]zeno Table[2.^e, {e, -3, 4, 0.08}];
ListLogLogPlot[
 {Table[{k, zenoRate[k/\[CapitalOmega]zeno]}, {k, zenoKs}], Table[{k, \[CapitalOmega]zeno^2/(2 k)}, {k, Select[zenoKs, # >= \[CapitalOmega]zeno &]}]},
 Joined -> {False, True},
 PlotLegends -> {"measured rate  |Re \[Lambda]|", "\!\(\*FractionBox[SuperscriptBox[\(\[CapitalOmega]\), \(2\)], \(2 k\)]\)  (the Zeno freeze)"},
 Frame -> True, GridLines -> {{{\[CapitalOmega]zeno, Directive[Gray, Dashed]}}, None},
 FrameLabel -> {"measurement strength k", "slowest relaxation rate  |Re \[Lambda]|"},
 PlotLabel -> "the rate follows k up to k = \[CapitalOmega], then freezes as \!\(\*FractionBox[SuperscriptBox[\(\[CapitalOmega]\), \(2\)], \(2 k\)]\)"]
```

Below $k = \Omega$ the rate climbs along $k$; above it the points settle onto the $\Omega^2/(2k)$ line, the leftover relaxation frozen as $1/k$, the small gap near the peak closing as the watching strengthens. Watching past the optimal strength does not pin the qubit harder, it slows what motion is left, the quantum Zeno suppression.

The whole crossover is exact, read straight off the deterministic Liouvillian with no gate, no time step, and no ensemble. Its slowest eigenvalue returns at once in the next section as the width of the central peak the quantum point contact reads off its own current, so the Zeno rate here and that spectrum there are one number seen twice. Coherent Rabi giving way to overdamped freezing through a Liouvillian exceptional point is the modern reading of when watching stops helping and starts freezing ([Snizhko, Kumar, and Romito, 2020](https://arxiv.org/abs/2003.10476)), seen on a driven transmon by [Slichter and Siddiqi, 2016](https://arxiv.org/abs/1512.04006).

### A Charge Qubit Read by a Quantum Point Contact

**The problem.** The last two examples watched an abstract qubit; here is the same watching realized in a solid-state device. A charge qubit (one electron shared between two quantum dots, "dot 1" or "dot 2") is monitored by a nearby quantum point contact (QPC) whose current depends on which dot holds the electron. With interdot tunnelling $\Omega$ and no detuning, the unconditional master equation is
$$\dot\rho = -i[\tfrac{\Omega}{2}\hat\sigma_x,\rho] + \kappa\,\mathcal{D}[\vec n_1]\rho,\qquad \vec n_1 = \tfrac12(\mathbb{1} - \hat\sigma_z),$$
with $\hat\sigma_z = \vec n_2 - \vec n_1$ the which-dot operator. Here $\kappa$ is the coefficient of $\mathcal D[n_1]$; because $n_1$ is a projector, the off-diagonal density-matrix element decays at rate $\kappa/2$. The QPC current is a continuous weak measurement of $\hat\sigma_z$: when $\kappa < \Omega$ the electron tunnels coherently and the current carries an oscillation at the Rabi frequency; when $\kappa > \Omega$ the qubit is pinned and the current is a random telegraph. Build the conditioned dynamics, the detector current, and the current's power spectrum.

The dephasing $\vec n_1$ is a projector, and $\mathcal{D}[\vec n_1] = \tfrac14\mathcal{D}[\hat\sigma_z]$ (the identity part drops from any dissipator), so watching the QPC is exactly watching $\hat\sigma_z$ at strength $\kappa/4$: the machinery of the last two examples applies unchanged. Fix the tunnelling and the two regimes:

```wl
OmQpc = 2.0; kSlow = 0.6; kFast = 12.0; qpcLeak[\[Kappa]_] := Sqrt[\[Kappa]/4] Z;
```

The conditioned trajectory uses the toolkit's positivity-preserving step. After subtracting the state-independent mean QPC current and choosing the sign convention $Z=n_2-n_1$, the rescaled record is $\dot J = \sqrt{\kappa}\,\langle\hat\sigma_z\rangle + \text{shot noise}$. An experimental current can differ by an overall gain and sign; neither changes the conditioned dynamics. Fix the step and horizon:

```wl
dtQpc = 0.004; tfQpc = 10.;
```

the weakly measured, coherent history ($\kappa = 0.6$):

```wl
qpcSlow = trajectory[densityMatrix[ground], (OmQpc/2) X, {qpcLeak[kSlow]}, {1.}, {}, dtQpc, tfQpc, 4];
```

the strongly measured, Zeno history ($\kappa = 12$):

```wl
qpcFast = trajectory[densityMatrix[ground], (OmQpc/2) X, {qpcLeak[kFast]}, {1.}, {}, dtQpc, tfQpc, 4];
```

and their shared time axis:

```wl
whenQpc = qpcSlow["times"];
```

The QPC hands back a single real record, $dJ = \sqrt{\kappa}\,\langle\hat\sigma_z\rangle\,dt + dW$, one continuous measurement of $\hat\sigma_z$: it is the homodyne ledger again in a solid-state skin, the record the only data and the coherent which-dot swing inferred from it, buried per step under shot noise many times larger.

Contrast the two regimes directly, the coherent swing against the strongly-watched telegraph:

```wl
ListLinePlot[{Transpose[{whenQpc, zOf@qpcSlow["states"]}], Transpose[{whenQpc, zOf@qpcFast["states"]}]},
 Frame -> True, GridLines -> Automatic,
 PlotLegends -> {"weak (coherent tunnelling)", "strong (Zeno telegraph)"}, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"}, PlotLabel -> "measurement strength selects coherent oscillation or telegraph"]
```

Weak measurement leaves the electron tunnelling freely; strong measurement pins it in one dot with rare jumps.

The distinctive QPC signature is not in any single record but in the *spectrum*, because a coherent oscillation is a *frequency*: it cancels in the time-integral of the record and shows itself only in the spectrum, painting a peak at the Rabi frequency in the fluctuations of $\hat\sigma_z$ that the detector current inherits. The power spectrum is the Fourier transform of the steady-state autocorrelation $\langle\delta\hat\sigma_z(\tau)\,\delta\hat\sigma_z(0)\rangle_{\mathrm{ss}}$, where $\delta\hat\sigma_z = \hat\sigma_z - \langle\hat\sigma_z\rangle_{\mathrm{ss}}$ is the fluctuation about the steady mean. Quantum regression is the statement that this correlation runs on the *same* Liouvillian as the state: plant the seed $\delta\hat\sigma_z\,\rho_{\mathrm{ss}}$, propagate it by $e^{\mathcal L\tau}$, and read $\delta\hat\sigma_z$ back off it, so $C(\tau) = \mathrm{Tr}[\delta\hat\sigma_z\,e^{\mathcal L\tau}(\delta\hat\sigma_z\,\rho_{\mathrm{ss}})]$.

With the convention $S(\omega)=\int_{-\infty}^{\infty}e^{-i\omega\tau}C(\tau)\,d\tau$, Fourier-transforming in $\tau$ replaces the propagator $e^{\mathcal L\tau}$ by a single matrix inverse, the resolvent $(i\omega - \mathcal L)^{-1}$. If the Liouvillian is diagonalizable, the inverse splits into simple poles, one per eigenvalue $\lambda_j$. An eigenvalue $\lambda_j = -\gamma_j + i\nu_j$ supplies a pole at $\nu_j$ with width $\gamma_j$, but a complex residue can mix absorptive (Lorentzian) and dispersive line shapes; it is not generally a positive Lorentzian by itself. The coefficient $c_j$ is the seed's coordinate along mode $j$, and $r_j$ is the readout overlap. Summing the poles, with the steady mode $\lambda_j = 0$ dropped because $\delta\hat\sigma_z$ has zero steady mean,
$$S(\omega) = 2\,\mathrm{Re}\sum_{\lambda_j\neq 0}\frac{c_j\,r_j}{i\omega - \lambda_j}.$$
This is where the slowest eigenvalue of the Zeno crossover returns: it is the pole nearest the imaginary axis and controls the narrowest long-time spectral scale. When its residue is real and it dominates the other modes, that scale appears as a Lorentzian half-width equal to the relaxation rate.

Build the general spectrum for any Liouvillian, seed, and readout: eigendecompose the Liouvillian, find the seed's coordinates $c_j$ by solving for them in the eigenbasis, read each mode's overlap $r_j$ with the readout, and sum the poles at every requested frequency:

```wl
regressionSpectrum[big_, seedMat_, readout_, tones_] := Module[{d = Length[readout], n = Length[big], vals, vecs, basis, singular, coeffs, resu},
   {vals, vecs} = Eigensystem[big];
   basis = Transpose[vecs]; singular = SingularValueList[N@basis];
   If[MatrixRank[basis] < n || Last[singular]/First[singular] < 10^-7,
    Return@Failure["DefectiveEigenbasis", <|"NumericalRank" -> MatrixRank[N@basis], "ConditionNumber" -> First[singular]/Last[singular]|>]];
   coeffs = LinearSolve[basis, Flatten[seedMat]];
   resu = Table[Tr[readout . ArrayReshape[vecs[[j]], {d, d}]], {j, d^2}];
   Table[2 Re@Total@Table[If[Abs[vals[[j]]] < 10^-9, 0, coeffs[[j]] resu[[j]]/(I w - vals[[j]])], {j, d^2}], {w, tones}]];
```

`vals` are the eigenvalues $\lambda_j$; `coeffs` are the $c_j$, produced by `LinearSolve` writing the seed as a combination of the eigenmodes; `resu` are the $r_j$, each mode's trace against the readout; the last line sums $c_j r_j/(i\omega - \lambda_j)$ over the nonzero eigenvalues.

The eigenmode sum is one way to the spectrum. There is a second that never diagonalizes at all: Fourier-transforming $C(\tau)$ turns $e^{\mathcal L\tau}$ straight into the resolvent $(i\omega - \mathcal L)^{-1}$, so $S(\omega) = 2\,\mathrm{Re}\,\mathrm{Tr}[\hat R\,(i\omega - \mathcal L)^{-1}\,\hat S]$ for readout $\hat R$ and seed $\hat S$. In other words, solve one linear system per frequency instead of splitting the inverse into poles. At $\omega = 0$ the propagated seed stays traceless, so we pin the physical branch by appending the trace condition and solving in the least-squares sense:

```wl
resolventSpectrum[big_, seedMat_, readout_, tones_] := With[{n = Length[big], d = Length[readout]},
   Table[2 Re@Tr[readout . ArrayReshape[If[Abs[w] < 10^-9,
       LeastSquares[Join[N[-big], {1. Flatten[IdentityMatrix[d]]}], Join[Flatten[seedMat], {0.}]],
       LinearSolve[I w IdentityMatrix[n] - big, Flatten[seedMat]]], {d, d}]], {w, tones}]];
```

Now specialize to the QPC. Build its Liouvillian and steady state, center $\hat\sigma_z$ on its steady value to form $\delta\hat\sigma_z$, and hand the symmetrized seed $(\delta\hat\sigma_z\,\rho_{\mathrm{ss}} + \rho_{\mathrm{ss}}\,\delta\hat\sigma_z)/2$ together with that same $\delta\hat\sigma_z$ as the readout to the resolvent form:

```wl
qpcSpectrum[\[Kappa]_, tones_] := With[{big = liouvillian[(OmQpc/2) X, {qpcLeak[\[Kappa]]}, 2], rss = First@stationary[(OmQpc/2) X, {qpcLeak[\[Kappa]]}]},
   With[{dz = Z - Re[Tr[Z . rss]] id2}, resolventSpectrum[big, (dz . rss + rss . dz)/2, dz, tones]]];
```

Here `dz` is $\delta\hat\sigma_z$, and symmetrizing the seed selects the real, even correlation measured by the detector. `qpcSpectrum` deliberately uses the mode-free resolvent, so the parameter sweep remains defined at the exceptional point; the modal representation is used only as a cross-check away from coalescence.

Confirm the normalization is right: the spectrum integrates to the steady variance of $\hat\sigma_z$, which is one for the maximally mixed steady state:

```wl
With[{wide = Range[-40, 40, 0.02]}, 0.02 Total[qpcSpectrum[1., wide]]/(2 Pi)]
```

It returns one, so the spectral weight is accounted for.

Confirm the two ways give the same QPC spectrum, at a ratio far from any coalescence:

```wl
With[{big = liouvillian[(OmQpc/2) X, {qpcLeak[4.]}, 2], rss = First@stationary[(OmQpc/2) X, {qpcLeak[4.]}], grid = Range[-6, 6, 0.1]},
 With[{dz = Z - Re[Tr[Z . rss]] id2}, Max@Abs[regressionSpectrum[big, (dz . rss + rss . dz)/2, dz, grid] - resolventSpectrum[big, (dz . rss + rss . dz)/2, dz, grid]]]]
```

Away from coalescence, the maximum difference is at numerical precision: two representations of one spectrum. They part company at one place only. Where two Liouvillian eigenvalues and their eigenvectors coalesce, an **exceptional point** (met in the Zeno crossover above, and reached again in the squeezing spectrum of the last section), the eigenbasis the mode sum needs turns singular while the resolvent never forms it. The QPC sweep therefore uses the resolvent for every ratio, including the exact exceptional point.

Now sweep the ratio $r = k/\Omega$ of the effective measurement strength $k = \kappa/4$ to the drive, from weak, past the point where the two spectral peaks merge, into the Zeno regime. Take five ratios along that path:

```wl
ratioSweep = {0.25, 0.5, 0.75, 1., 1.5};
```

Resolve each spectrum on a fine frequency grid wide enough to hold both Rabi peaks near $\pm\Omega$:

```wl
tonesQpc = Range[-6, 6, 0.04];
```

Colour the curves warm to cool across the sweep, so the peaks' inward march reads at a glance:

```wl
sweepCols = ColorData["TemperatureMap"] /@ Subdivide[Length[ratioSweep] - 1];
```

Overlay the five power spectra, one per ratio (recovering $\kappa = 4\Omega r$ for each), with gridlines at $\pm\Omega$:

```wl
ListLinePlot[Table[Transpose[{tonesQpc, qpcSpectrum[4 OmQpc r, tonesQpc]}], {r, ratioSweep}],
 PlotStyle -> (Directive[#, Thickness[0.006]] & /@ sweepCols),
 PlotLegends -> LineLegend[sweepCols, {"k/\[CapitalOmega] = 0.25  (weak)", "k/\[CapitalOmega] = 0.5", "k/\[CapitalOmega] = 0.75", "k/\[CapitalOmega] = 1  (exceptional point)", "k/\[CapitalOmega] = 1.5  (Zeno)"}],
 Frame -> True, GridLines -> {{-OmQpc, 0, OmQpc}, None}, PlotRange -> All, AspectRatio -> 1/2,
 FrameLabel -> {"frequency", "\!\(\*SubscriptBox[\(S\), \(zz\)]\)(omega)"},
 PlotLabel -> "watching harder marches the spectral peaks inward,\n merging them before the exceptional point"]
```

At weak watching the spectrum is split, twin peaks near $\pm\Omega$, the coherent Rabi oscillation written into the fluctuations of $\hat\sigma_z$, the solid-state echo of the Mollow triplet's sidebands. Watching harder marches those peaks inward, but not to the pole frequency: the damping pushes each maximum below it, to $\pm\sqrt{\Omega^2 - 2k^2}$, so the peaks reach zero and merge already at $k = \Omega/\sqrt2$ (ratio $r \approx 0.71$). That merger is one event; the exceptional point is another. The Liouvillian poles sit at $-k \pm i\sqrt{\Omega^2 - k^2}$ and keep oscillating until $k = \Omega$ (ratio $r = 1$, that is $\kappa = 4\Omega$), where the two eigenvalues *and their eigenvectors* coalesce and the Liouvillian turns defective. Spectral merger comes first, eigenvector coalescence later: at the plotted $r = 0.75$ the spectrum already shows a single central maximum even though the poles are still a damped-oscillating pair. Deep in the Zeno regime the narrow slow pole dominates the central feature and its width scales as $\Omega^2/(2k)$, so the rate tracked in the previous section reappears spectrally.

The two thresholds are easiest to see if we track the peak location itself. Plot the spectral maximum $\sqrt{\Omega^2 - 2k^2}$ and, beside it, the pole frequency $\sqrt{\Omega^2 - k^2}$ against the ratio, marking where each reaches zero:

```wl
With[{Om = OmQpc},
 ListLinePlot[{Table[{r, Re@Sqrt[Om^2 - 2 (r Om)^2]}, {r, 0, 1.5, 0.005}], Table[{r, Re@Sqrt[Om^2 - (r Om)^2]}, {r, 0, 1.5, 0.005}]},
  PlotStyle -> {Automatic, Dashed}, PlotLegends -> {"spectral maximum", "pole frequency"},
  Frame -> True, GridLines -> {{{1/Sqrt[2], Directive[Gray, Dashed]}, {1, Directive[Gray, Dashed]}}, Automatic},
  PlotRange -> All, FrameLabel -> {"k/\[CapitalOmega] ratio", "frequency"},
  PlotLabel -> "spectral peak reaches zero at k/\[CapitalOmega] = \!\(\*FractionBox[\(1\), SqrtBox[\(2\)]]\), poles only at k/\[CapitalOmega] = 1"]]
```

The solid curve, the spectral maximum, hits zero first at $k = \Omega/\sqrt2$; the dashed curve, the pole frequency, survives to $k = \Omega$. Between the two gray markers the spectrum is already single-peaked while the trajectory still rings, so the five-ratio sweep above steps over both thresholds without landing on either.

What a lab actually records is the detector *current* spectrum, the flat shot-noise floor with this $S_{\mathrm{zz}}(\omega)$ riding on it (scaled by the measurement coupling); a fundamental result, the Korotkov-Averin bound, caps that current peak at four times the shot-noise background for a quantum-limited detector, the solid-state statement of the measurement limit that runs through this whole catalog.

Two different knobs move two different crossovers, and they should not be confused: the *measurement strength* $\kappa$ sets coherent-oscillation versus Zeno-telegraph behaviour (shown above), while the detector *transparency* sets whether the record is read as discrete electron jumps (a counting unravelling like the photon counter) or as this diffusive current. Both unravellings average to the one master equation above.

### Measurement-Induced Localization

**The problem.** The general continuous measurement of an observable (here $\hat\sigma_z$) on an *undriven* qubit, written in the catalog's measurement-strength convention:
$$d\rho = k\,\mathcal D[\hat\sigma_z]\rho\,dt + \sqrt{k}\,(\hat\sigma_z\rho + \rho\hat\sigma_z - 2\langle\hat\sigma_z\rangle\rho)\,dW,\qquad dy = \langle\hat\sigma_z\rangle\,dt + \frac{dW}{2\sqrt{k}}.$$
With nothing pushing, the state's knowledge sharpens toward an eigenstate of $\hat\sigma_z$ at the rate set by $k$: each history collapses to one pole, the choice fixed only by the noise, while the average over histories stays balanced. Watch many histories make their choice.

Turn off the drive and measure $\hat\sigma_z$ on a state leaning toward $+z$, $\langle\hat\sigma_z\rangle = z_0 = 0.5$. With no Hamiltonian, measurement alone drives each run to a $\hat\sigma_z$ eigenstate, $+1$ or $-1$, and the record that drives the collapse is also what reads out which pole it chose. Fix the biased start, the measured operator, and a fine grid:

```wl
z0Loc = 0.5; x0Loc = Sqrt[1 - z0Loc^2]; measZ = Sqrt[1.5] Z;
dtLoc = 0.002; tfLoc = 2.; gridLoc = dtLoc Range[0, Round[tfLoc/dtLoc]];
```

One run reads the $z$ component throughout; run a fan of them, and beside them solve the master equation for the record-blind mean:

```wl
locZ[seed_] := zOf@trajectory[blochState[x0Loc, 0, z0Loc], 0 id2, {measZ}, {1.}, {}, dtLoc, tfLoc, seed]["states"];
fanLoc = Table[locZ[s], {s, 60}];
locMeanRun = evolveODE[0 id2, {measZ}, blochState[x0Loc, 0, z0Loc], tfLoc];
meanLoc = zOf[locMeanRun /@ gridLoc];
```

Now visualize the fan against that conserved mean: every displayed path approaches a pole while the thick master-equation curve holds flat at $z_0$:

```wl
Show[ListLinePlot[Transpose[{gridLoc, #}] & /@ fanLoc, PlotStyle -> Directive[ColorData[97, 1], Opacity[0.12]]],
 ListLinePlot[Transpose[{gridLoc, meanLoc}], PlotStyle -> Directive[Thick, ColorData[97, 2]]],
 Frame -> True, GridLines -> Automatic, PlotRange -> {-1.05, 1.05},
 FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket]"},
 PlotLabel -> "each run approaches a pole;\n the master-equation mean holds at \!\(\*SubscriptBox[\(z\), \(0\)]\)"]
```

By the finite endpoint every displayed run is numerically close to a pole, while the thick curve never leaves $z_0$: the dephasing $\mathcal{D}[\hat\sigma_z]$ leaves the populations untouched, so $\langle\hat\sigma_z\rangle$ is conserved on average. In continuous time, $z(t)$ is a bounded martingale and converges asymptotically to $\pm1$; it does not hit a pole at a generic finite time. The limiting split follows from martingale convergence (equivalently, optional-stopping arguments applied to bounded approximating stopping times): $\mathbb E[z(\infty)]=z_0$, so $p_+ - p_- = z_0$ and $p_+ + p_- = 1$ give $p_\pm = (1 \pm z_0)/2$. Born is a theorem about the limiting martingale, not a fraction the finite run count has to fit. The symmetric $z_0 = 0$ start, the state the earlier examples called $\left|+\right\rangle$, is the equal-weight special case.

The record is what a detector actually delivers, and it both drives the collapse and reads it out. Accumulate each run's record into the integrated readout $J(t) = \int_0^t dJ$ and keep the pole it lands on:

```wl
readOut[seed_] := With[{run = trajectory[blochState[x0Loc, 0, z0Loc], 0 id2, {measZ}, {1.}, {}, dtLoc, tfLoc, seed]},
   {Prepend[Accumulate[Flatten@run["record"]], 0.], Sign[Last@zOf@run["states"]]}];
readData = Table[readOut[s], {s, 400}];
```

Now visualize the integrated record, one path per run, colored by the pole the run reached: each single increment is almost pure noise, but the accumulated record fans into two well-separated groups:

```wl
With[{show = readData[[;; 50]]},
 ListLinePlot[Transpose[{gridLoc, #[[1]]}] & /@ show,
  PlotStyle -> (If[#[[2]] > 0, Directive[ColorData[97, 1], Opacity[0.45]], Directive[ColorData[97, 2], Opacity[0.45]]] & /@ show),
  Frame -> True, GridLines -> Automatic, PlotRange -> All,
  FrameLabel -> {"time", "integrated record  J(t)"}, PlotLabel -> "the integrated record splits by the pole each run chose"]]
```

The sign of $J(T)$ is the readout's verdict on the eigenstate; measure how often it is right as the integration window grows:

```wl
fidLoc[T_] := With[{iT = Round[T/dtLoc] + 1}, Mean[Boole[Sign[#[[1, iT]]] == #[[2]]] & /@ readData]];
ListLinePlot[Table[{T, fidLoc[T]}, {T, dtLoc, tfLoc, 0.01}], PlotRange -> All,
 Frame -> True, GridLines -> Automatic,
 FrameLabel -> {"integration time T", "assignment fidelity"}, PlotLabel -> "longer integration reads the eigenstate out with certainty"]
```

For these 400 seeded runs the empirical fidelity reaches one by the displayed endpoint. That is a finite-sample observation, not a finite-time certainty: at any fixed duration a Gaussian-noise record retains a nonzero error probability, while the assignment error tends to zero asymptotically as more record accumulates. Two things keep this honest rather than a lookup. The limiting eigenstate is not fixed in advance; it is selected by the same record used to infer it. And the flat-weight integral $J$ is not the sharpest estimator: the conditional $\langle\hat\sigma_z\rangle(t)$ of the fan above is the Bayesian-optimal filter, which a uniform integral only approximates. This trajectory-read-from-its-own-record is what [Murch, Weber, Macklin, and Siddiqi](https://arxiv.org/abs/1305.7270) tracked on a superconducting qubit, reconstructing the conditioned state from the readout and confirming it by tomography.

### Quantum Kalman Filter: Tracking an Oscillator

**The problem.** Continuously measure the position of a harmonic oscillator; a Gaussian state stays Gaussian and the stochastic master equation collapses to a quantum Kalman-Bucy filter. The conditional covariance $\Sigma$ obeys a *deterministic* matrix Riccati equation, independent of the record,
$$\dot \Sigma = A \Sigma + \Sigma A^{T} + D - \Sigma\,C^{T} C\,\Sigma,$$
with the free-rotation drift, the position readout, and the measurement back-action
$$A = \begin{pmatrix}0 & 1\\-1 & 0\end{pmatrix},\qquad C = \begin{pmatrix}2\sqrt{k} & 0\end{pmatrix},\qquad D = \begin{pmatrix}0 & 0\\0 & k\end{pmatrix}.$$
In the three independent entries $(\Sigma_{\mathrm{xx}}, \Sigma_{\mathrm{xp}}, \Sigma_{\mathrm{pp}})$ of $\Sigma$ this reads
$$\frac{d}{dt}\begin{pmatrix}\Sigma_{\mathrm{xx}}\\ \Sigma_{\mathrm{xp}}\\ \Sigma_{\mathrm{pp}}\end{pmatrix} = \begin{pmatrix}2\Sigma_{\mathrm{xp}} - 4k\,\Sigma_{\mathrm{xx}}^{2}\\ \Sigma_{\mathrm{pp}} - \Sigma_{\mathrm{xx}} - 4k\,\Sigma_{\mathrm{xx}}\Sigma_{\mathrm{xp}}\\ k - 2\Sigma_{\mathrm{xp}} - 4k\,\Sigma_{\mathrm{xp}}^{2}\end{pmatrix},$$
the quadratic $-4k(\cdots)$ terms being the information gain $-\Sigma C^{T} C \Sigma$, which is what makes this a Riccati and not a linear Lyapunov equation. Balanced against the back-action ($+D$), that gain relaxes the uncertainty to a steady value fixed only by how hard you watch, squeezing it below the oscillator's own quietest spread; the running best-guess then follows a linear equation driven by the record. Track the shrinking uncertainty and check its standstill.

Continuously measure the position of a harmonic oscillator and you can track it with a Gaussian estimate: a conditional mean and covariance. The covariance obeys a *deterministic* Riccati equation, independent of the record, set only by the measurement strength. It is the warm particle's covariance law from Part Three, $\dot\Sigma = A\Sigma + \Sigma A^{T} + D$, now carrying the measurement's information gain $-\Sigma C^{T} C\Sigma$: diffusion widened that blob to thermal, watching squeezes this one below the ground state. Set the strength and integrate the Riccati equations from a broad start:

```wl
kKal = 1;
riccati = NDSolveValue[{vx'[t] == 2 vc[t] - 4 kKal vx[t]^2, vc'[t] == vp[t] - vx[t] - 4 kKal vx[t] vc[t],
   vp'[t] == -2 vc[t] + kKal - 4 kKal vc[t]^2, vx[0] == 3., vc[0] == 0., vp[0] == 3.}, {vx, vc, vp}, {t, 0, 12}];
```

Assemble the covariance matrix at time $t$ from the Riccati solution:

```wl
covariance[t_] := {{riccati[[1]][t], riccati[[2]][t]}, {riccati[[2]][t], riccati[[3]][t]}};
```

Compare the position variance at the start and at the end:

```wl
{riccati[[1]][0], riccati[[1]][12]}
```

The position variance has shrunk from its broad start below the ground-state variance of one half: continuous measurement conditionally squeezes it below the zero-point level. Two checks. First, the conditional state is pure, its uncertainty product at the Heisenberg floor:

```wl
Det@covariance[12.]
```

The determinant is within the ODE tolerance of the Heisenberg floor; the continuum conditional state is pure, and the numerical solution converges toward that value. Second, the steady covariance is independently the fixed point of the Riccati equation:

```wl
{N[vx /. Solve[{2 vc - 4 kKal vx^2 == 0, vp - vx - 4 kKal vx vc == 0,
      -2 vc + kKal - 4 kKal vc^2 == 0, vx > 0}, {vx, vc, vp}, Reals][[1]]], riccati[[1]][12.]}
```

The fixed point matches the late value.

That fixed point has a clean closed form. With $s = \sqrt{1 + 4k^2}$,
$$\Sigma_{\mathrm{xx}} = \frac{1}{\sqrt{2(1+s)}}, \qquad \Sigma_{\mathrm{xp}} = \frac{s-1}{4k}, \qquad \Sigma_{\mathrm{pp}} = s\,\Sigma_{\mathrm{xx}},$$
and it sits exactly on the Heisenberg floor $\det\Sigma = 1/4$ for every measurement strength. Put the closed form beside the numerical fixed point:

```wl
With[{s = Sqrt[1 + 4 kKal^2]},
 TableForm[{{1/Sqrt[2 (1 + s)], (s - 1)/(4 kKal), s/Sqrt[2 (1 + s)]}, Through[riccati[12.]]},
  TableHeadings -> {{"closed form", "Riccati at t = 12"}, {"\!\(\*SubscriptBox[\(\[CapitalSigma]\), \(xx\)]\)", "\!\(\*SubscriptBox[\(\[CapitalSigma]\), \(xp\)]\)", "\!\(\*SubscriptBox[\(\[CapitalSigma]\), \(pp\)]\)"}}]]
```

The two rows agree: position is conditionally squeezed below one half, momentum broadens to hold the product at the floor, and the tilt $\Sigma_{\mathrm{xp}}$ carries the correlation the measurement builds. The conditional state is pure at any watching strength, not only the one plotted.

The conditional state is Gaussian, so picture it as a one-sigma ellipse: the covariance $\Sigma$ around a center $c$ traces $c + \Sigma^{1/2}\{\cos s, \sin s\}$ as $s$ runs the circle. Define it:

```wl
covarianceEllipse[c_, \[CapitalSigma]_] := With[{pts = Table[c + Re[MatrixPower[\[CapitalSigma], 1/2]] . {Cos[s], Sin[s]}, {s, 0, 2 Pi, 0.12}]}, Append[pts, First[pts]]];
```

Now visualize the covariance ellipse over time, broad at first and then tightening and tilting to a small steady one:

```wl
trackTimes = {0, 0.1, 0.25, 0.5, 1, 2, 4, 8, 12};
Legended[Graphics[{Thick, MapIndexed[{ColorData["SunsetColors"][First[#2]/Length[trackTimes]], Line@covarianceEllipse[{0, 0}, covariance[#1]]} &, trackTimes]},
  Frame -> True, GridLines -> Automatic, AspectRatio -> 1, ImageSize -> 400, FrameLabel -> {"position", "momentum"},
  PlotLabel -> "measurement shrinks the covariance to a tight, tilted, steady one"],
 SwatchLegend[ColorData["SunsetColors"][#/Length[trackTimes]] & /@ Range[Length[trackTimes]],
  "t = " <> ToString[#] & /@ trackTimes]]
```

This conditional estimate is exactly what a lab computes in real time to track a trapped particle or membrane, more precisely than its own zero-point spread.

The covariance is only half the filter. It says how sharp the estimate is; the other half is the estimate itself, the conditional *mean*, threading through the noisy record in real time. Watch the filter track. Set a fine grid and draw one record's worth of Wiener increments:

```wl
kalDt = 0.004; kalN = 1500; kalGrid = kalDt Range[0, kalN];
kalDW = BlockRandom[SeedRandom[7]; RandomVariate[NormalDistribution[0, Sqrt[kalDt]], kalN]];
```

A note on the constant $k$. In the catalog's measurement-strength convention the monitored leak is $\hat c = \sqrt{k}\,\hat x$, so the informative terms carry $4k$ and the mean gain is $2\sqrt{k}\,V$; the record noise and the gain must then be consistent, which pins the record normalization below. What the code calls a "reference conditional mean" is the physical conditional mean of one simulated record, driven by its own innovation $dW$ with that gain. Roll it out from a displaced start:

```wl
kalTrue = FoldList[Function[{xp, jw}, With[{j = jw[[1]], dw = jw[[2]]},
     {xp[[1]] + xp[[2]] kalDt + 2 Sqrt[kKal] riccati[[1]][kalGrid[[j]]] dw,
      xp[[2]] - xp[[1]] kalDt + 2 Sqrt[kKal] riccati[[2]][kalGrid[[j]]] dw}]],
   {3., 0.}, Transpose[{Range[kalN], kalDW}]];
```

All the observer ever sees is the record, the signal $\langle\hat x\rangle\,dt$ buried in white noise. Consistency with the $4k$-Riccati and the $2\sqrt k$ gain fixes its noise level at $dy = \langle\hat x\rangle\,dt + dW/(2\sqrt k)$:

```wl
kalRec = MapThread[#1[[1]] kalDt + #2/(2 Sqrt[kKal]) &, {Most@kalTrue, kalDW}];
```

Now the filter. Write its running estimates of the conditional means as $\langle x\rangle$ and $\langle p\rangle$, without the operator hat that marks the true expectations $\langle\hat x\rangle, \langle\hat p\rangle$ the record carries, so the innovation below is exactly their difference. Start it from a deliberately *wrong* guess for the position and feed it only the record: at each step it forms its own innovation, the record minus its own prediction, $d\widehat W = 2\sqrt{k}\,(dy - \langle x\rangle\,dt)$, and corrects the estimate with the Riccati gain,
$$d\langle x\rangle = \langle p\rangle\,dt + 2\sqrt{k}\,\Sigma_{\mathrm{xx}}(t)\,d\widehat W, \qquad d\langle p\rangle = -\langle x\rangle\,dt + 2\sqrt{k}\,\Sigma_{\mathrm{xp}}(t)\,d\widehat W,$$
the deterministic Riccati covariance $\Sigma(t)$ setting the gain:

```wl
kalFilt = FoldList[Function[{xp, jdy}, With[{j = jdy[[1]], dwf = 2 Sqrt[kKal] (jdy[[2]] - xp[[1]] kalDt)},
     {xp[[1]] + xp[[2]] kalDt + 2 Sqrt[kKal] riccati[[1]][kalGrid[[j]]] dwf,
      xp[[2]] - xp[[1]] kalDt + 2 Sqrt[kKal] riccati[[2]][kalGrid[[j]]] dwf}]],
   {-3., 0.}, Transpose[{Range[kalN], kalRec}]];
```

Read the $\pm 1\sigma$ band off the covariance:

```wl
kalBand = Sqrt[riccati[[1]][#]] & /@ kalGrid;
```

Now visualize the tracking: the reference conditional mean and the estimate from the incorrect prior. The deliberately wrong initial estimate is outside its own $\pm1\sigma$ band; the band quantifies conditional covariance under the assumed filter model, not uncertainty over an arbitrarily misspecified initial mean. The record then pulls the estimate toward the reference path:

```wl
Legended[Show[
  ListLinePlot[{Transpose[{kalGrid, kalFilt[[All, 1]] + kalBand}], Transpose[{kalGrid, kalFilt[[All, 1]] - kalBand}]},
   Filling -> {1 -> {2}}, PlotStyle -> Directive[ColorData[97, 1], Opacity[0.15]]],
  ListLinePlot[{Transpose[{kalGrid, kalTrue[[All, 1]]}], Transpose[{kalGrid, kalFilt[[All, 1]]}]},
   PlotStyle -> {Directive[Thick, ColorData[97, 2]], Directive[Thick, ColorData[97, 1]]}],
  Frame -> True, GridLines -> Automatic, PlotRange -> All,
  FrameLabel -> {"time", "position"}, PlotLabel -> "the filter locks onto the reference conditional mean from the record alone"],
 LineLegend[{ColorData[97, 2], ColorData[97, 1]}, {"reference conditional mean", "estimate from wrong prior \[PlusMinus]\[Sigma]"}]]
```

The estimate begins wildly wrong (opposite sign and outside its stated band), then the record pulls it close to the reference conditional mean while the band approaches the steady conditional width: this observable model forgets the deliberately wrong mean. The reference mean is not a hidden classical position; it is the conditional mean generated by the simulated innovation.

Confirm the lock, comparing the two once the filter has settled:

```wl
Abs[Last@kalTrue[[All, 1]] - Last@kalFilt[[All, 1]]]
```

The endpoint difference is about $3\times10^{-3}$ for this seeded Euler discretization. It is small on the unit position scale but not zero; refining `kalDt` is the relevant convergence check.

The engine of the filter is its own innovation, the record minus the *filter's* prediction (not the reference mean's). For an optimal filter this must be white noise, zero-mean and uncorrelated; measure its autocorrelation past the transient, with the $\pm 2/\sqrt N$ band that white noise would sit inside:

```wl
kalInnov = MapThread[2 Sqrt[kKal] (#2 - #1[[1]] kalDt) &, {Most@kalFilt, kalRec}][[Round[kalN/3] ;;]]/Sqrt[kalDt];
With[{ci = 2/Sqrt[Length[kalInnov]]},
 ListPlot[CorrelationFunction[kalInnov, {1, 10}], Filling -> Axis, PlotRange -> {-0.2, 0.2},
  Frame -> True, GridLines -> {None, {{ci, Directive[Gray, Dashed]}, {-ci, Directive[Gray, Dashed]}}},
  ImageSize -> 440, PlotMarkers -> Automatic,
  FrameLabel -> {"lag", "innovation autocorrelation"}, PlotLabel -> "the filter's innovation is white within confidence limits"]]
```

For this finite record the displayed nonzero-lag correlations lie within the approximate white-noise bands. Innovation whiteness is a necessary and useful model diagnostic, not by itself a proof of optimality: a wrong model can pass a short autocorrelation check, and a correct one can show excursions in finite data.

### Measurement Feedback: Steering With the Record

**The problem.** Feed the homodyne record back as a Hamiltonian $\hat H_{\mathrm{fb}} = \dot J\,\hat F$, the feedback signal being the homodyne current $\dot J = dJ/dt$ read off the record $dJ = \langle\hat c + \hat c^\dagger\rangle\,dt + dW$, so over a step the loop applies $\hat F\,dJ$; averaged over the noise, the conditioned-plus-feedback evolution becomes the Wiseman-Milburn feedback master equation
$$\dot\rho = -i[\hat H,\rho] + \mathcal{D}[\hat c]\rho - i[\hat F,\ \hat c\rho + \rho\hat c^\dagger] + \tfrac{1}{\eta}\,\mathcal{D}[\hat F]\rho,$$
with $\eta$ the detection efficiency, this is just a plain leak with a shifted operator, $\mathcal{D}[\hat c - i\hat F]\rho$, plus a Hamiltonian correction. The feedback engineers the effective damping, here steering the maximally mixed qubit onto a chosen state. Drive $\langle\hat\sigma_x\rangle$ to one with the record itself and check it against this equation.

Close the loop: feed the measurement record back as a Hamiltonian in proportion, so the record steers the state. Concretely the measured operator is $\hat c = \sqrt{G}\,\hat\sigma_z$ and the feedback generator is $\hat F = \sqrt{G}\,\hat\sigma_y$, so the loop measures $\hat\sigma_z$ and rotates about $Y$ by an amount set by the record. The one constant $G$ does two jobs here: it is the measurement rate, $\mathcal{D}[\hat c] = G\,\mathcal{D}[\hat\sigma_z]$, and it is the feedback gain, entering through $\hat F = \sqrt{G}\,\hat\sigma_y$; the two are matched on purpose, which is exactly what makes the target a dark state of the shifted leak below. The shifted operator is then $\hat c - i\hat F = \sqrt{G}\,(\hat\sigma_z - i\hat\sigma_y)$, and the Hamiltonian correction $\tfrac12(\hat c^\dagger\hat F + \hat F\hat c)$ vanishes because $\hat\sigma_z$ and $\hat\sigma_y$ anticommute, so the ideal loop is the single leak $G\,\mathcal{D}[\hat\sigma_z - i\hat\sigma_y]$. Set the gain and the time grid:

```wl
Gfb = 0.8; tickFb = 0.01; track = tickFb Range[0, 200];
```

A rotation about $Y$ has a closed form, $R_y(\phi) = e^{-i\phi\hat\sigma_y} = \cos\phi\,\mathbb{1} - i\sin\phi\,\hat\sigma_y$, and it is what the generator $\hat F$ produces: the feedback unitary applied each step is $e^{-i\hat F\,dJ} = R_y(\sqrt{G}\,dJ)$, so $\hat F$ is the generator and $R_y$ the finite rotation, not one and the same. The step is light:

```wl
rotateY[\[Phi]_] := Cos[\[Phi]] id2 - I Sin[\[Phi]] Y;
```

One steered run, one step at a time. Write $\rho^{\mathrm{meas}}$ for the positivity-preserving measurement update of $\rho$ under $\hat c = \sqrt{G}\,\hat\sigma_z$, given the record $dJ = 2\sqrt{G}\,\langle\hat\sigma_z\rangle_\rho\,dt + dW$. One step measures, then rotates by the commanded angle $\phi = \sqrt{G}\,dJ$:
$$\rho' = R_y(\phi)\,\rho^{\mathrm{meas}}\,R_y(\phi)^\dagger,\qquad \phi = \sqrt{G}\,dJ,$$
run $N = 200$ times from $\rho_0 = \tfrac12\mathbb{1}$. Keeping the pair $(\rho,\phi)$ at each step, the run is the whole causal chain at once,
$$((\rho_0,\,0),\ (\rho_1,\phi_1),\ \dots,\ (\rho_N,\phi_N)),$$
each state beside the command that produced it, the seed angle a placeholder $0$. One function then feeds every plot below, the response curves from its states and the command panel from its angles:

```wl
steer[seed_] := BlockRandom[SeedRandom[seed];
   With[{measure = measurementStep[0 id2, {Sqrt[Gfb] Z}, {1.}, {}, tickFb]},
    FoldList[
     Function[{prev, dw},
      With[{rho = prev[[1]]},
       With[{dJ = measurementRecord[rho, {Sqrt[Gfb] Z}, {1.}, tickFb, {dw}]},
        With[{u = rotateY[Sqrt[Gfb] dJ[[1]]]}, {u . measure[rho, dJ] . ConjugateTranspose[u], Sqrt[Gfb] dJ[[1]]}]]]],
     {id2/2, 0.},
     RandomVariate[NormalDistribution[0, Sqrt[tickFb]], 200]]]];
```

This feedback has a target: the $+X$ state is a fixed point the loop leaves alone, so from the maximally mixed state the loop should build it. Confirm $+X$ is the dark state (the shifted operator annihilates it):

```wl
Chop[(Sqrt[Gfb] (Z - I Y)) . plus]
```

The zero vector confirms it. Integrate the feedback master equation for the predicted $\langle\hat\sigma_x\rangle$:

```wl
predicted = blochVector[evolveODE[0 id2, {Sqrt[Gfb] (Z - I Y)}, id2/2, 2.][#]][[1]] & /@ track;
```

Average sixty steered runs:

```wl
steeredMean = Mean[Table[xOf[steer[s][[All, 1]]], {s, 60}]];
```

Now visualize one run, the ensemble mean, and the master-equation prediction climbing together:

```wl
ListLinePlot[{Transpose[{track, xOf[steer[2][[All, 1]]]}], Transpose[{track, steeredMean}], Transpose[{track, predicted}]},
 PlotStyle -> {Directive[ColorData[97, 1], Opacity[0.45]], Directive[Thick, ColorData[97, 2]], Directive[Dashed, Red]},
 PlotLegends -> {"one steered run", "average of many", "predicted"}, Frame -> True, GridLines -> Automatic,
 PlotRange -> {-0.35, 1.05}, ImageSize -> 520, FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]"},
 PlotLabel -> "feedback drives the maximally mixed state to a pure state"]
```

$\langle\hat\sigma_x\rangle$ climbs from zero to one: the maximally mixed state is driven to a pure state by its own record. This is measurement-based feedback control, the basis of feedback cooling and qubit stabilization.

The loop has a causal chain worth seeing whole: the noisy record drives a control command, which drives the state's response. That same run already carries both, so take seed 2 and read its per-step command against its response:

```wl
feedRun = steer[2];
```

Now visualize the loop in time, two panels side by side: the feedback rotation the raw record commands at each step, a noise band, and the state's response $\langle\hat\sigma_x\rangle$ climbing to the target:

```wl
Row[{
  ListLinePlot[Transpose[{track[[2 ;;]], feedRun[[2 ;;, 2]]}], Frame -> True, GridLines -> Automatic,
   ImageSize -> Medium, PlotRange -> All, FrameLabel -> {"time", "feedback rotation"},
   PlotLabel -> "the command the raw record drives"], "  ",
  ListLinePlot[Transpose[{track[[2 ;;]], blochVector[#[[1]]][[1]] & /@ feedRun[[2 ;;]]}], Frame -> True, GridLines -> Automatic,
   ImageSize -> Medium, PlotRange -> All, FrameLabel -> {"time", "\[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]"},
   PlotLabel -> "the response it produces"]}]
```

Command to response, left to right: the loop reads the noisy record as a rotation, a noise band, and the rotations walk the state smoothly to the target.

But a real loop never has a perfect detector, and inefficiency caps how pure the steered state can get. The Wiseman-Milburn feedback master equation carries an efficiency $\eta$ that adds an irreducible feedback-noise term $\tfrac{1-\eta}{\eta}\mathcal{D}[\hat F]$; solve its steady state against $\eta$:

```wl
steadyXeff[\[Eta]_] := blochVector[First@stationary[0 id2, {Sqrt[Gfb] (Z - I Y), Sqrt[(1 - \[Eta])/\[Eta]] Sqrt[Gfb] Y}]][[1]];
etasFb = Range[0.3, 1, 0.05];
ListLinePlot[Transpose[{etasFb, steadyXeff /@ etasFb}], PlotMarkers -> Automatic, PlotStyle -> ColorData[97, 2], Frame -> True,
 GridLines -> Automatic, PlotRange -> {0, 1.05}, ImageSize -> 460, FrameLabel -> {"detection efficiency \[Eta]", "steady \[LeftAngleBracket]\[Sigma]x\[RightAngleBracket]"},
 PlotLabel -> "a leaky detector cannot fully purify: feedback noise sets the floor"]
```

Steady $\langle\hat\sigma_x\rangle$ climbs to one only as $\eta \to 1$; a lossy detector leaves the loop feeding its own missing information back as noise, and the target is reached only in part. Perfect stabilization needs a perfect record, the same efficiency limit that capped the feedback cooling.

## Part Six: Settling Warm, Purifying Fast, Cooling Cold

### Thermalization to a Bath

**The problem.** A warm bath both takes and gives energy, settling a system at its own temperature. For an oscillator,
$$\dot\rho = -i[\omega\hat a^\dagger\hat a,\rho] + \gamma(n_T+1)\,\mathcal{D}[\hat a]\rho + \gamma\,n_T\,\mathcal{D}[\hat a^\dagger]\rho,\qquad n_T = \frac{1}{e^{\beta\omega} - 1},$$
where the down-rate $\gamma(n_T+1)$ beats the up-rate $\gamma n_T$ by exactly the Boltzmann factor $e^{\beta\omega}$. This *detailed balance* forces the steady state to the thermal distribution $\rho_{\mathrm{ss}}\propto e^{-\beta\omega\,\hat a^\dagger\hat a}$ with mean occupation $n_T$. Settle a cold and a hot start and confirm both land on the thermal distribution.

A warm bath both absorbs and emits, settling a system at its temperature. Model it as two dissipators: decay $\hat a$ at rate $\gamma(n_T+1)$ and excitation $\hat a^\dagger$ at rate $\gamma n_T$. Set the cutoff, the ladder, the bath occupation and its rate:

```wl
topTh = 24; a = annihilation[topTh]; count = creation[topTh] . a; nT = 1.2; \[Gamma]th = 0.5;
```

The two channels, decay and excitation, with their thermal weights:

```wl
fall = Sqrt[\[Gamma]th (nT + 1)] a; climb = Sqrt[\[Gamma]th nT] creation[topTh];
```

The mean occupation of a state:

```wl
units[rho_] := Re@expectation[count, rho];
```

Settle a cold (vacuum) and a hot start:

```wl
fromCold = evolveODE[0 IdentityMatrix[topTh], {fall, climb}, densityMatrix[coherentState[topTh, 0]], 16.];
fromHot = evolveODE[0 IdentityMatrix[topTh], {fall, climb}, densityMatrix[UnitVector[topTh, 9]], 16.];
```

Read as a dynamical law, this equation has three consequences worth drawing out one at a time, each sharper than the one before.

The first is a **single attractor**: the evolution forgets the initial state, so every start relaxes to the same steady state. Test that from opposite extremes, the cold vacuum and the hot $n=9$ state, by reading their late mean occupations against $n_T$:

```wl
{units[fromCold[16.]], units[fromHot[16.]], nT}
```

Both land on $n_T$. Watch them get there, the cold climbing and the hot relaxing until they meet:

```wl
ListLinePlot[{Transpose[{Range[0, 16, 0.05], units[fromCold[#]] & /@ Range[0, 16, 0.05]}],
   Transpose[{Range[0, 16, 0.05], units[fromHot[#]] & /@ Range[0, 16, 0.05]}], {{0, nT}, {16, nT}}},
 PlotLegends -> {"from cold", "from hot", "warm equilibrium"}, Frame -> True, GridLines -> Automatic,
 FrameLabel -> {"time", "mean occupation"}, PlotLabel -> "both starts reach the thermal mean", PlotRange -> All]
```

They settle on the same $n_T$ from opposite sides: the bath drains the hot oscillator and feeds the cold one, so wherever the state began, it ends where decay and excitation balance. This two-way traffic is what makes the attractor unique; the pure-dephasing bath of Part One, which only froze populations, left a whole line of steady states instead.

In the infinite oscillator model, the mean occupation follows an exact closed path. Tracing the master equation against $\hat a^\dagger\hat a$ gives $\dot{\langle n\rangle} = -\gamma(\langle n\rangle - n_T)$, so $\langle n(t)\rangle = n_T + (\langle n(0)\rangle - n_T)e^{-\gamma t}$. The finite-cutoff simulation should approach that law; check the cold start ($\langle n(0)\rangle = 0$):

```wl
Max@Table[Abs[units[fromCold[t]] - nT (1 - Exp[-\[Gamma]th t])], {t, 0, 16, 1}]
```

The reported nonzero residual is the combined cutoff and ODE error; increasing the dimension and tightening solver tolerances reduce it. In the untruncated law one relaxation rate $\gamma$ carries the mean from its start to the bath value. The first field moment obeys the same kind of law but half as fast and while rotating, $\dot{\langle\hat a\rangle} = -(i\omega + \gamma/2)\langle\hat a\rangle$; the code sets $\omega = 0$, so no rotation is drawn here: the frequency appears in coherence rather than number populations.

The second consequence is stronger: the equation fixes **the whole state**, not just its mean. A mean is weak evidence, since a coherent, a squeezed, and a thermal state can all share the occupation $n_T$. So build the settled populations and the geometric thermal shape $p_n = (1-r)\,r^n$ with $r = n_T/(n_T+1) = e^{-\beta\omega}$:

```wl
pops = Re@Diagonal[fromCold[16.]][[;; 8]];
warmShape = With[{ratio = nT/(nT + 1)}, (1 - ratio) ratio^Range[0, 7]];
```

Read the largest gap between them:

```wl
Max@Abs[pops - warmShape]
```

They agree to a fraction of a percent. Level by level, the settled populations sit on the thermal bars rung for rung:

```wl
BarChart[Transpose@{pops, warmShape}, ChartLegends -> {"settled", "thermal"}, ChartLabels -> {Range[0, 7], None},
 Frame -> True, FrameLabel -> {"Fock level", "probability"}, PlotLabel -> "the settled distribution is thermal"]
```

The attractor is a complete density matrix; reaching the thermal mean was only its shadow.

The third consequence is the deepest, the point where a dynamical equation reproduces equilibrium statistical mechanics. Stationarity alone forces the population ratio on adjacent rungs to equal the ratio of the two rates, excitation $\gamma n_T$ over decay $\gamma(n_T+1)$, and that flux balance is mere bookkeeping. The physics is that the ratio is not free: through **detailed balance** the bath locks the downward rate above the upward one by exactly $e^{\beta\omega}$, its own thermal signature, so an oscillator with no temperature of its own is driven into the Boltzmann distribution at the bath's temperature. The canonical ensemble is derived here, not assumed:

```wl
{(\[Gamma]th nT)/(\[Gamma]th (nT + 1)), nT/(nT + 1)}
```

The two are equal, both $e^{-\beta\omega}$. The settled state shows the same thing directly: a geometric distribution is a straight line on a log axis, its constant slope the Boltzmann factor. Plot the settled populations against the thermal shape:

```wl
With[{fock = Range[0, 6]},
 ListLogPlot[{Transpose[{fock, pops[[;; Length@fock]]}], Transpose[{fock, warmShape[[;; Length@fock]]}]},
  Joined -> True, PlotMarkers -> Automatic, PlotStyle -> {Automatic, Dashed},
  PlotLegends -> {"settled", "thermal"}, Frame -> True, FrameTicks -> {fock, Automatic}, GridLines -> {fock, Automatic},
  PlotRange -> All, FrameLabel -> {"Fock level n", "population (log)"},
  PlotLabel -> "geometric tail: a straight line on a log axis"]]
```

The settled populations fall on a straight line, parallel to the thermal one and with the same constant slope: every rung is less likely than the one below by the bath's Boltzmann factor $e^{-\beta\omega}$, the visible fingerprint of detailed balance. The ratio of rates in the equation has become the slope of the ladder in the state. Damping and thermalization are then one process: decay $\hat a$ and excitation $\hat a^\dagger$ are the same bath, and the ground-state relaxation of the earlier examples is its $n_T \to 0$ corner, an oscillator settling to its bath temperature with its whole distribution fixed by two rates.

### Rapid Purification by Feedback

**The problem.** Continuously watching an observable extracts information and sharpens the state toward one of its eigenstates, so a monitored qubit purifies with time. Measure $\hat M = \sin\phi\,\hat\sigma_x + \cos\phi\,\hat\sigma_z$ at strength $k$ with no drive (the leak is $\hat c = \sqrt{k}\,\hat M$), a stationary axis tilted by $\phi$ from $z$; the conditioned state obeys the same measurement equation as the localization example, now along that tilted axis,
$$d\rho_c = k\,\mathcal D[\hat M]\rho_c\,dt + \sqrt{k}\,(\hat M\rho_c + \rho_c\hat M - 2\langle\hat M\rangle\rho_c)\,dW.$$

Watch it happen before dissecting the rate: fix one tilted axis, $\phi = 0.7$, measure it at unit strength on the maximally mixed state, and follow the purity $\mathrm{Tr}(\rho^2)$ of five conditioned runs.

```wl
mAxis = Sin[0.7] X + Cos[0.7] Z;
purities = Table[purity /@ trajectory[id2/2, 0 id2, {mAxis}, {1.}, {}, 0.01, 1.5, k]["states"], {k, 5}];
ListLinePlot[Transpose[{0.01 Range[0, Length[#] - 1], #}] & /@ purities,
 Frame -> True, FrameLabel -> {"time", "purity"}, PlotRange -> {Automatic, {0.4, 1.02}},
 GridLines -> {None, {0.5, 1.}}, PlotStyle -> Thick,
 PlotLabel -> "monitoring purifies every conditioned run"]
```

All five runs climb from $\tfrac12$ toward $1$: monitoring alone has turned the maximally mixed qubit nearly pure, the raw content of "purifies with time." Each climb is jagged because these are *conditioned* histories, each driven by its own record's noise; averaging the conditioned states over many records leaves the qubit fixed at $\tfrac12$, since the maximally mixed state is a steady state of the unconditional evolution. Purification is a fact about each watched trajectory, not the record-blind average. Measuring one stationary axis, held put for the whole run, is the unsteered case; the steered strategy below re-aims the axis each step instead.

How fast it purifies depends on where the state sits relative to that axis, an angle that itself moves: let $\theta_t$ be the instantaneous angle between the Bloch vector $\vec a_t$ and $\hat M$. ($\theta_t$ is not $\phi$: $\phi$ fixes where $\hat M$ points and stays constant, while $\theta_t$ measures where the state sits relative to it and moves with the record.) Writing the impurity as the random variable $\ell_t = \tfrac12(1 - \left|\vec a_t\right|^2)$, Itô's rule turns the measurement equation above into a full stochastic increment,
$$d\ell_t = -4k\,\ell_t\,(\sin^2\theta_t + 2\ell_t\cos^2\theta_t)\,dt \;-\; 4\sqrt{k}\,\ell_t\sqrt{1 - 2\ell_t}\,\cos\theta_t\,dW_t.$$
The noise carries the factor $\cos\theta_t$, and that factor is the whole story. If feedback holds the axis *perpendicular* to the Bloch vector, $\theta_t = \pi/2$, the noise coefficient vanishes and the impurity is deterministic for *every* record, $d\ell_t = -4k\,\ell_t\,dt$, so $\ell_t = \ell_0\,e^{-4kt}$: the fast case. If instead the axis lies *along* the Bloch vector, $\theta_t = 0$, the noise is maximal and the impurity keeps fluctuating,
$$d\ell_t = -8k\,\ell_t^2\,dt - 4\sqrt{k}\,\ell_t\sqrt{1 - 2\ell_t}\,dW_t,$$
with no deterministic pathwise solution; only its ensemble mean $\mu(t) = \mathbb{E}[\ell_t]$ can be tracked, through the nonclosed $\dot\mu(t) = -8k\,\mathbb{E}[\ell_t^2]$. A stationary measurement axis lets the state drift into alignment ($\theta_t \to 0$), where it learns slowly and stochastically: the back-action (the measurement's kick back on the state) pulls the Bloch vector toward whatever axis is measured, so an axis that never moves is one the state keeps falling into. Holding $\theta_t = \pi/2$ therefore takes constant re-aiming, and that re-aiming is the feedback: it keeps the axis perpendicular to the moving state and with it the fast, deterministic rate. Because the feedback resets the axis to perpendicular only at discrete steps, the state drifts off $\theta_t = \pi/2$ within each step, and the $\cos\theta_t$ noise term in $d\ell_t$ no longer cancels exactly; the residue is a record-dependent path spread that vanishes as $h\to0$. This is [Jacobs' rapid-purification protocol](https://arxiv.org/abs/quant-ph/0301056).

There is a catch in what "fast" means. Keeping the measurement perpendicular maximizes the *mean purity at a fixed time*, yet it *lengthens* the mean time to reach a chosen purity, because its determinism forbids the lucky fast runs an aligned measurement occasionally enjoys; the two figures of merit disagree by a factor of two ([Wiseman and Ralph, 2006](https://arxiv.org/abs/quant-ph/0603062)). We build both — an **unsteered** $\hat\sigma_z$ measurement and a **steered** axis re-aimed perpendicular to the state at every step — and let the distributions settle it.

As loops, both share a measure-observe-update core; the steered one wraps it in two extra steps that recompute the axis each pass.

*Steered:* estimate $\rho_j$ from the past record $\to$ compute $\vec a_j = \langle\vec\sigma\rangle$ $\to$ choose $\vec n_j$ perpendicular to $\vec a_j$ $\to$ measure $\vec n_j\cdot\vec\sigma$ $\to$ observe $\Delta J_j$ $\to$ update $\rho_{j+1}$ $\to$ repeat with a new $\vec n_{j+1}$.

*Unsteered:* estimate $\rho_j$ from the past record $\to$ measure the fixed $\vec n\cdot\vec\sigma$ $\to$ observe $\Delta J_j$ $\to$ update $\rho_{j+1}$ $\to$ repeat with the same $\vec n$.

Both strategies start from the maximally mixed state and are scored by the same impurity $\ell_t = \tfrac12(1 - \left|\vec a_t\right|^2)$, which is $\tfrac12$ when mixed and $0$ when pure. Given a density matrix, define a function to compute the impurity:

```wl
impurity[rho_] := (1 - Norm[blochVector[rho]]^2)/2;
```

Set the measurement strength and the time grid:

```wl
kPur = 1.; tickPur = 0.01; span = tickPur Range[0, 120];
```

The **unsteered** strategy measures one axis, $\hat\sigma_z$, for the whole run. Each step feeds the record $dJ = 2\sqrt{k}\,\langle\hat\sigma_z\rangle\,dt + dW$ into the positivity-preserving update, and the conditioned state diffuses and lengthens; but as its Bloch vector swings toward the $z$-axis the measurement falls into alignment and the purification slows. This is exactly the toolkit's `trajectory` with a fixed $\hat\sigma_z$ leak, so it needs no function of its own; the ensembles below call `trajectory` for it directly and score the result with `impurity`.

The **steered** strategy uses the record differently. Instead of a stationary axis it reads the current Bloch vector off the conditioned state and measures along an axis *perpendicular* to it, holding the fast $\theta_t = \pi/2$ angle at every step.

In the plainest words: look at where the qubit points now, aim the measurement axis perpendicular to it, make one weak measurement, update the state with the outcome, and repeat. The feedback is just this repeated re-aiming of the measurement direction, the steered chain above. The code never physically rotates the state with a feedback unitary (a control pulse applied conditioned on the record); it only changes which observable is measured next.

One function carries the whole steered loop, with the same signature as `trajectory`: an initial state $\rho_0$ and Hamiltonian $H$, the measurement strength $k$, the step $dt$, the horizon $tf$, and a seed. Each step reads the Bloch vector, aims the axis perpendicular to it (the horizontal $\vec n = \vec a \times \hat z$), and advances the conditioned state by the positivity-preserving update. It returns the observer's data in `trajectory`'s shape, the times, states, and record, plus the axes it chose:

```wl
steeredBlur[rho0_, H_, k_, dt_, tf_, seed_] := BlockRandom[SeedRandom[seed];
  Module[{n = Round[tf/dt], perp, dws, states, axes, records},
   perp[lean_] := With[{c = Cross[lean, {0, 0, 1.}]}, If[Norm[c] < 1.*^-6, {1., 0., 0.}, Normalize[c]]];
   dws = RandomVariate[NormalDistribution[0, Sqrt[dt]], n];
   states = FoldList[
     Function[{rho, dw}, With[{leak = Sqrt[k] (perp[blochVector[rho]] . {X, Y, Z})},
        measurementStep[H, {leak}, {1.}, {}, dt][rho, measurementRecord[rho, {leak}, {1.}, dt, {dw}]]]],
     rho0, dws];
   axes = perp[blochVector[#]] & /@ Most[states];
   records = MapThread[measurementRecord[#1, {Sqrt[k] (#2 . {X, Y, Z})}, {1.}, dt, {#3}] &, {Most[states], axes, dws}];
   <|"times" -> dt Range[0, n], "states" -> states, "axes" -> axes, "record" -> records|>]];
```

In formulas: with $\vec a = \langle\vec\sigma\rangle_\rho$ the current Bloch vector, each step measures the leak $\hat c = \sqrt{k}\,(\vec n\cdot\vec\sigma)$ along a unit axis $\vec n$ perpendicular to $\vec a$. Because that axis is perpendicular to the state the signal vanishes and the record is pure noise,
$$dJ = 2\sqrt{k}\,(\vec n\cdot\vec a)\,dt + dW = dW,$$
which is nonetheless the fast case: a state perpendicular to the axis is an equal superposition of the measured operator's eigenstates, where the back-action is strongest. The step is then the positivity-preserving measurement update of $\rho$ under $\hat c$ with that record, and nothing else, with no feedback unitary; the steering is the re-choice of $\vec n$ at the next node. Any direction perpendicular to $\vec a$ gives $\theta_t = \pi/2$ and the same rate; `steeredBlur` takes the horizontal one, $\vec a \times \hat z$, so the conditioned state runs in the horizontal plane while $\hat M$ stays fixed in its tilted $xz$-plane, unmeasured: only the angle $\theta_t$, never the plane the axis lies in, sets the rate. In continuous time this holds $\theta_t=\pi/2$ pathwise; at finite $h$ the residual spread is measured below. Folding from $\rho_0 = \tfrac12\mathbb{1}$ over $N$ Wiener increments, `steeredBlur` returns the states $\rho_0,\dots,\rho_N$ along with the record it saw and the axes it chose.

Look at what the observer actually holds on one run. The record it reads is featureless noise, while the axes it picks swing around to stay perpendicular to the sharpening state:

```wl
obs = steeredBlur[id2/2, 0 id2, kPur, tickPur, 1.2, 1];
Row[{
   ListLinePlot[Flatten[obs["record"]], DataRange -> MinMax@Most[obs["times"]], Frame -> True,
    FrameLabel -> {"time", "record"}, ImageSize -> Medium, PlotLabel -> "what the observer sees: pure noise"],
   ListLinePlot[Transpose@obs["axes"], DataRange -> MinMax@Most[obs["times"]], Frame -> True, ImageSize -> Medium,
    PlotLegends -> {"\!\(\*SubscriptBox[\(n\), \(x\)]\)", "\!\(\*SubscriptBox[\(n\), \(y\)]\)", "\!\(\*SubscriptBox[\(n\), \(z\)]\)"},
    FrameLabel -> {"time", "axis"}, PlotLabel -> "what the observer chooses: the re-aimed axis"]}]
```

The record is a driftless noise stream: with the axis held perpendicular to the state, the signal $\vec n\cdot\vec a$ vanishes and the meter reads pure $dW$, so it shows no drift toward either eigenstate. The striking part is that the purification *rate* does not depend on this record at all: the impurity falls as $\ell_0 e^{-4kt}$ on every history, the deterministic case. The geometry of the axis choices sets how fast the state purifies; the record's fluctuations only decide which of the two perpendicular eigenstates it lands on. The axis stays horizontal throughout ($n_z$ flat at zero) and swings in the plane to track the sharpening Bloch vector.

The steering pins $\theta_t = \pi/2$ at every re-aiming node by construction, so the quantity worth reading is how far the state slips within a step before the axis resets. Pair each held axis with the Bloch vector at the step's *end* and take the angle (`VectorAngle` normalizes both):

```wl
thetaDrift = MapThread[VectorAngle, {obs["axes"], Rest[blochVector /@ obs["states"]]}];
ListLinePlot[Transpose[{Rest[obs["times"]], thetaDrift}], Frame -> True,
 FrameLabel -> {"time", "\[Theta] between axis and state"},
 GridLines -> {None, {{\[Pi]/2, Directive[Red, Dashed]}}},
 PlotRange -> {\[Pi]/4, 3 \[Pi]/4},
 FrameTicks -> {{Range[0, 2 \[Pi], \[Pi]/4], None}, Automatic},
 PlotLabel -> "\[Theta] scatters about \[Pi]/2 within each step (the finite-h residue)"]
```

On average $\theta$ sits on $\pi/2$, but it scatters about it from step to step: within each interval the state tilts toward or away from the frozen axis before the next re-aiming pulls it back. That symmetric scatter, not any systematic bias, is the finite-$h$ residue behind the run-to-run spread measured below. The ragged opening is the near-mixed transient, where the Bloch vector is still near zero and points nowhere definite; the range clips its excursions below $\pi/4$.

Average eighty runs of each strategy:

```wl
unsteeredAvg = Mean[Table[impurity /@ trajectory[id2/2, 0 id2, {Sqrt[kPur] Z}, {1.}, {}, tickPur, 1.2, s]["states"], {s, 80}]];
steeredAvg = Mean[Table[impurity /@ steeredBlur[id2/2, 0 id2, kPur, tickPur, 1.2, s]["states"], {s, 80}]];
```

Now visualize both purification strategies on a log scale, where a straight line represents exponential decay:

```wl
ListLogPlot[{Transpose[{span, unsteeredAvg}], Transpose[{span, steeredAvg}]}, Joined -> True,
 Frame -> True, GridLines -> Automatic, PlotLegends -> {"unsteered measurement", "steered measurement"},
 FrameLabel -> {"time", "impurity (log scale)"}, PlotLabel -> "steering purifies faster; finite steps approach the deterministic limit"]
```

The steered mean is steeper and smoother, tracking the continuous-time exponential more closely than the unsteered mean. But determinism is a strong claim: theory says the steered impurity follows $\ell_0 e^{-4kt}$ on *every* record, with zero run-to-run spread, yet the simulated runs show a small spread and a small gap from that exponential. Is the residue real, or an artifact of re-aiming the axis only once per discrete step? The step is the suspect: within it the state drifts off $\theta_t = \pi/2$, so the noise no longer cancels. If that is the cause, refining $h$ must shrink both the spread and the gap to the deterministic target. Test at the same physical time $t = 0.6$ across three step sizes (each uses $n = 0.6/h$ steps to reach it):

```wl
TableForm[
  Table[With[{finals = Last[impurity /@ steeredBlur[id2/2, 0 id2, kPur, h, 0.6, #]["states"]] & /@ Range[80], target = 0.5 Exp[-4 kPur 0.6]},
     {h, Mean[finals], StandardDeviation[finals], Abs[Mean[finals] - target]}],
    {h, {0.02, 0.01, 0.005}}],
  TableHeadings -> {None, {"step h", "mean impurity", "spread", "|mean - target|"}}]
```

Both the spread and the bias contract as $h$ shrinks, so the residue is discretization (the finite step), not physics: determinism is a genuine property of the continuum limit, not something to read off the $h = 0.01$ grid.

But the mean hides the caveat, so look at the distributions. Set a four-unit window (400 steps) long enough that every one of the unsteered seeded runs reaches the target, as the explicit check below will confirm, and a count of 200 runs:

```wl
purGrid = tickPur Range[0, 400]; purN = 200;
```

Run the unsteered strategy over that ensemble, keeping each run's full impurity history:

```wl
unsteeredRunsPur = Table[impurity /@ trajectory[id2/2, 0 id2, {Sqrt[kPur] Z}, {1.}, {}, tickPur, 4., s]["states"], {s, purN}];
```

and the steered strategy the same way:

```wl
steeredRunsPur = Table[impurity /@ steeredBlur[id2/2, 0 id2, kPur, tickPur, 4., s]["states"], {s, purN}];
```

First, compare the finite-step steered runs with the continuous-time law $\ell_t = \ell_0 e^{-4kt}$, $\ell_0 = 1/2$. Overlay thirty paths on the analytic exponential, on a log scale:

```wl
ListLogPlot[Append[Transpose[{purGrid, #}] & /@ steeredRunsPur[[;; 30]], Transpose[{purGrid, 0.5 Exp[-4 kPur purGrid]}]],
 Joined -> True, PlotRange -> {{0, 2}, {10^-4, 1}}, Frame -> True, GridLines -> Automatic, ImageSize -> 480,
 PlotStyle -> Append[ConstantArray[Directive[ColorData[97, 2], Opacity[0.2]], 30], Directive[Thick, ColorData[97, 1], Dashed]],
 PlotLegends -> LineLegend[{ColorData[97, 1], ColorData[97, 2]}, {"continuous-time law", "finite-step steered runs"}],
 FrameLabel -> {"time", "impurity (log)"}, PlotLabel -> "finite-step paths fluctuate around the continuum exponential (dashed)"]
```

The paths cluster around, but do not coincide with, the dashed line at this step size. The refinement table above is the evidence that this cluster collapses onto the pathwise deterministic continuum law.

Now the objective-dependence. At a *fixed time* the steered strategy wins on the *mean*, but the distributions show why that is subtle. Animate both ensembles as they purify: the steered runs stay a tight cluster sliding toward zero, while the unsteered runs scatter, many already purer than the steered cluster and a slow tail lagging behind. It is that tail, not a uniformly worse spread, that lifts the unsteered mean above the steered's tight, deterministic value:

```wl
Animate[
  Histogram[{unsteeredRunsPur[[All, t]], steeredRunsPur[[All, t]]}, {0, .5, 0.01}, "Probability",
   PlotLegends -> {"unsteered", "steered"}, PlotRange -> {{0, .5}, All}, Frame -> True, AspectRatio -> 1/3,
   FrameLabel -> {"impurity at t = " <> ToString[purGrid[[t]]], "fraction"},
   PlotLabel -> "steered clusters; the unsteered runs scatter"],
  {t, 2, Floor[Length[purGrid]/2], 1}, AnimationRunning -> False]
```

But the *first-passage time* to a chosen purity tells the opposite story. For each run, record the time it first drops below $\ell = 0.05$; a run that does not cross within the window is left explicitly *right-censored*, logged as still short of the target at the end rather than handed the endpoint as if it had crossed there:

```wl
fptPur[run_] := With[{pos = FirstPosition[run, l_ /; l < 0.05]},
   If[MissingQ[pos], Missing["RightCensored", (Length[run] - 1) tickPur], (First[pos] - 1) tickPur]];
```

Collect that time for every run of each strategy:

```wl
{fptUnsteered, fptSteer} = {fptPur /@ unsteeredRunsPur, fptPur /@ steeredRunsPur};
```

Check the censor counts before computing a histogram or mean:

```wl
{Count[fptUnsteered, _Missing], Count[fptSteer, _Missing]}
```

The extended four-unit window gives zero censored runs for these seeds. If either count were nonzero, a plain (uncensored) average would be biased, ignoring the runs still short of the target, and the window would need extending or a survival-analysis estimate (the statistics of partly-observed times).

Histogram the two distributions:

```wl
Histogram[{DeleteMissing[fptUnsteered], DeleteMissing[fptSteer]}, {0, 1.2, 0.06}, "PDF", PlotLegends -> {"unsteered", "steered"},
 Frame -> True, AspectRatio -> 1/3,
 FrameLabel -> {"first-passage time to impurity < 0.05", "density"}, PlotLabel -> "yet the unsteered measurement often reaches the target sooner"]
```

The unsteered strategy reaches down to shorter times: its runs vary, and the lucky ones cross the target before the steered runs. At finite step the steered crossing times still have a narrow distribution; under refinement that distribution collapses toward the continuum strategy's single deterministic hitting time.

Compare the two means directly, fixed-time impurity against mean hitting time:

```wl
Grid[{{"", Style["unsteered", Bold], Style["steered", Bold]},
   {Style["impurity at t = 0.6", Bold], Mean[unsteeredRunsPur[[All, 61]]], Mean[steeredRunsPur[[All, 61]]]},
   {Style["first-passage time to impurity < 0.05", Bold], Mean@DeleteMissing[fptUnsteered], Mean@DeleteMissing[fptSteer]}},
  Frame -> All, Alignment -> Left, Spacings -> {2, 1}]
```

The steered strategy has the lower fixed-time impurity but the *longer* mean first-passage: its determinism forbids the lucky fast records that let an unsteered measurement occasionally reach the target early. Which protocol is "optimal" depends on the question, sharp average purity by a deadline or shortest typical time to a goal, and no single curve can answer both. This is the [Wiseman-Ralph](https://arxiv.org/abs/quant-ph/0603062) point, made by the distributions rather than asserted.

### Feedback Cooling of a Mechanical Oscillator

**The problem.** Continuously measure a mechanical oscillator's position through a cavity and push against its estimated motion to cool it, against a thermal bath. Its annihilation operator is $\hat b$, and position and momentum are the dimensionless quadratures
$$\hat x=\frac{\hat b+\hat b^\dagger}{\sqrt2}, \qquad \hat p=\frac{\hat b-\hat b^\dagger}{i\sqrt2}, \qquad \left[\hat x,\hat p\right]=i,$$
so the vacuum has $\Sigma_{\mathrm{xx}}=\Sigma_{\mathrm{pp}}=\tfrac12$ and a thermal state $\Sigma_{\mathrm{xx}}=\Sigma_{\mathrm{pp}}=n_T+\tfrac12$, the thermal start the covariance equations take below. The conditioned state obeys
$$d\rho = -i[\omega\hat b^\dagger\hat b - f(t)\hat x,\ \rho]\,dt + k\,\mathcal D[\hat x]\rho\,dt + \sqrt{\eta k}\,(\hat x\rho + \rho\hat x - 2\langle\hat x\rangle\rho)\,dW + \gamma(n_T+1)\mathcal{D}[\hat b]\rho\,dt + \gamma n_T\mathcal{D}[\hat b^\dagger]\rho\,dt.$$
Watching sharpens the estimate, and the feedback force $f \propto -G\langle\hat p\rangle$ drains the energy by opposing the motion: it takes the estimated momentum, which way the oscillator is heading, and pushes back, a velocity-proportional drag that bleeds off kinetic energy. What is left is a balance of measurement rate, efficiency $\eta$, gain $G$, and heating $\gamma n_T$. Watching alone only heats; feedback cools, but only up to an optimal gain beyond which it feeds the detector's noise back onto the oscillator. Cool a warm start, find the best gain, and sweep the detector quality.

The same measure-and-feedback loop cools a mechanical oscillator. The loop is the closed cycle of control: each step measures the oscillator's position, updates the running estimate of its motion, turns that estimate into a force, and applies it back, then repeats. Cooling takes more than watching, though: the loop must act. Three pieces enter now that the Kalman tracker did not have. A thermal bath (damping $\gamma$, occupation $n_T$) keeps trying to warm the oscillator; the continuous position measurement (strength $k$, efficiency $\eta$) tracks it; and a feedback force pushes back against the *estimated* velocity, $f = -G\,\langle\hat p\rangle_c$ (the conditional mean), draining energy. The state is Gaussian throughout, so, as in the Kalman filter, the conditional covariance is deterministic and the conditional mean is stochastic. This example stands alone, so fix its own constants:

```wl
kCool = 1.; etaCool = 1.; gammaCool = 0.1; nCool = 10.;
```

The conditional covariance is the tracking Riccati of Part Five again, now carrying the bath (a drift toward the thermal spread $n_T + \tfrac12$ at rate $\gamma$) and the efficiency $\eta$ on the information-gain terms, while the back-action term $+k$ (the momentum kick that measuring position unavoidably costs) stays full. In the three entries $(\Sigma_{\mathrm{xx}}, \Sigma_{\mathrm{xp}}, \Sigma_{\mathrm{pp}})$,
$$\frac{d}{dt}\begin{pmatrix}\Sigma_{\mathrm{xx}}\\ \Sigma_{\mathrm{xp}}\\ \Sigma_{\mathrm{pp}}\end{pmatrix} = \begin{pmatrix}2\Sigma_{\mathrm{xp}} - 4\eta k\,\Sigma_{\mathrm{xx}}^{2} - \gamma\Sigma_{\mathrm{xx}} + \gamma(n_T + \tfrac12)\\ \Sigma_{\mathrm{pp}} - \Sigma_{\mathrm{xx}} - 4\eta k\,\Sigma_{\mathrm{xx}}\Sigma_{\mathrm{xp}} - \gamma\Sigma_{\mathrm{xp}}\\ k - 2\Sigma_{\mathrm{xp}} - 4\eta k\,\Sigma_{\mathrm{xp}}^{2} - \gamma\Sigma_{\mathrm{pp}} + \gamma(n_T + \tfrac12)\end{pmatrix},$$
integrated from the thermal start $\Sigma_{\mathrm{xx}} = \Sigma_{\mathrm{pp}} = n_T + \tfrac12$, $\Sigma_{\mathrm{xp}} = 0$:

```wl
riccatiCool[\[Eta]_] := NDSolveValue[{
    vx'[t] == 2 vc[t] - 4 \[Eta] kCool vx[t]^2 - gammaCool vx[t] + gammaCool (nCool + 1/2),
    vc'[t] == vp[t] - vx[t] - 4 \[Eta] kCool vx[t] vc[t] - gammaCool vc[t],
    vp'[t] == -2 vc[t] + kCool - 4 \[Eta] kCool vc[t]^2 - gammaCool vp[t] + gammaCool (nCool + 1/2),
    vx[0] == nCool + 1/2, vc[0] == 0., vp[0] == nCool + 1/2}, {vx, vc, vp}, {t, 0, 80}];
```

The steady covariance, read off the late-time solution, is the sharpest the observer can ever know the oscillator:

```wl
steadyCovCool[\[Eta]_] := Through[riccatiCool[\[Eta]][80.]];
```

Take its three entries at perfect efficiency, $\eta = 1$:

```wl
{vxCool, vcCool, vpCool} = steadyCovCool[etaCool];
```

Its mean phonon number (the count of vibrational quanta) is the *conditional floor*, the lowest occupation the record can pin down:

```wl
condFloorCool = (vxCool + vpCool)/2 - 0.5;
```

The Kalman gains $2\sqrt{k}\,\Sigma_{\mathrm{xx}}$ and $2\sqrt{k}\,\Sigma_{\mathrm{xp}}$ set how hard each record increment steers the estimate:

```wl
kxCool = 2 Sqrt[kCool] vxCool; kpCool = 2 Sqrt[kCool] vcCool;
```

Here is the key point, and the reason watching alone is not cooling: averaged over all records, the unconditional occupation is the conditional floor plus the spread of the estimate about it.

The conditional mean $\vec m = (\langle\hat x\rangle, \langle\hat p\rangle)$ follows a linear stochastic differential equation (SDE) driven by the innovation $dW$ (the unpredicted part of the record),
$$d\vec m = A_G\,\vec m\,dt + \vec K\,dW, \qquad A_G = \begin{pmatrix} -\gamma/2 & \omega \\ -\omega & -\gamma/2 - G \end{pmatrix}, \qquad \vec K = \begin{pmatrix} k_x \\ k_p \end{pmatrix},$$
component by component (with $\omega = 1$ here),
$$d\langle\hat x\rangle = (-\tfrac{\gamma}{2}\langle\hat x\rangle + \omega\langle\hat p\rangle)\,dt + k_x\,dW, \qquad d\langle\hat p\rangle = (-\omega\langle\hat x\rangle - (\tfrac{\gamma}{2}+G)\langle\hat p\rangle)\,dt + k_p\,dW.$$
The drift $A_G$ is the closed loop written out: the off-diagonal $\pm\omega$ is the free rotation, the $-\gamma/2$ on the diagonal is bath damping, and the lone $-G$ on the $\langle\hat p\rangle$ row is the feedback; the gains $(k_x, k_p) = (2\sqrt{k}\,\Sigma_{\mathrm{xx}}, 2\sqrt{k}\,\Sigma_{\mathrm{xp}})$ are the ones just defined. Averaged over records, the covariance $S$ of $\vec m$ obeys
$$\dot S = A_G S + S A_G^{T} + \vec K\vec K^{T},$$
and at steady state,
$$A_G S + S A_G^{T} + \vec K\vec K^{T} = 0,$$
the continuous Lyapunov equation, with the rank-one Kalman-gain noise $\vec K\vec K^{T} = \begin{pmatrix} k_x^2 & k_x k_p \\ k_x k_p & k_p^2 \end{pmatrix}$, solved by:

```wl
nUncCool[gg_] := With[{s = LyapunovSolve[{{-gammaCool/2, 1}, {-1, -gammaCool/2 - gg}},
     -{{kxCool^2, kxCool kpCool}, {kxCool kpCool, kpCool^2}}]},
   (s[[1, 1]] + vxCool + s[[2, 2]] + vpCool)/2 - 0.5];
```

With no feedback, $G = 0$, watching does not cool at all: it *heats*, because measurement back-action kicks the momentum and nothing removes the energy. Confirm this two ways, the conditional-plus-spread bookkeeping against the plain unconditional master equation (back-action on $\Sigma_{\mathrm{pp}}$, bath, no information gain):

```wl
{nUncCool[0.],
 With[{u = NDSolveValue[{vx'[t] == 2 vc[t] - gammaCool vx[t] + gammaCool (nCool + 1/2),
      vc'[t] == vp[t] - vx[t] - gammaCool vc[t], vp'[t] == -2 vc[t] + kCool - gammaCool vp[t] + gammaCool (nCool + 1/2),
      vx[0] == nCool + 0.5, vc[0] == 0., vp[0] == nCool + 0.5}, {vx, vp}, {t, 0, 200}]},
  (u[[1]][200.] + u[[2]][200.])/2 - 0.5]}
```

Both give the same number, above the bath's $n_T = 10$: the meter alone leaves the oscillator hotter than it found it.

Now close the loop. The conditional mean is a stochastic trajectory driven by the innovation, with the feedback force damping the momentum estimate:

```wl
Gcool = 2.; dtCool = 0.01; nStepCool = 2500; timesCool = dtCool Range[0, nStepCool];
```

At $t=0$ the oscillator is thermal, so the total covariance is $(n_T + 1/2)\mathbb{1}$. That total splits into the conditional floor the filter always carries plus the spread of the estimate about it, so the means must be drawn with the full *difference* $(n_T+\tfrac12)\mathbb{1} - \Sigma$, including its off-diagonal $-\Sigma_{\mathrm{xp}}$, to avoid double-counting the initial uncertainty and to make the initial unconditional state exactly thermal (the conditional covariance $\Sigma$ is held at its steady value throughout, the filter assumed already settled at its steady covariance, so the cooling transient sits entirely in the mean-spread $S$):

```wl
meanRunCool[gg_, seed_] := BlockRandom[SeedRandom[seed];
   Module[{init = RandomVariate[MultinormalDistribution[{0., 0.},
        {{nCool + 0.5 - vxCool, -vcCool}, {-vcCool, nCool + 0.5 - vpCool}}]],
     kicks = RandomVariate[NormalDistribution[0, Sqrt[dtCool]], nStepCool]},
    FoldList[Function[{xp, dw},
      {xp[[1]] + (xp[[2]] - gammaCool/2 xp[[1]]) dtCool + kxCool dw,
       xp[[2]] + (-xp[[1]] - gammaCool/2 xp[[2]] - gg xp[[2]]) dtCool + kpCool dw}],
     init, kicks]]];
```

The unconditional occupation in time is the ensemble spread of these means plus the conditional floor. Average over an ensemble of runs:

```wl
nEnsCool = 120;
```

For a feedback gain, `nDrainCool` runs that ensemble and returns the occupation at each sampled time, $n(t) = \tfrac12\big[(S_{\mathrm{xx}} + \Sigma_{\mathrm{xx}}) + (S_{\mathrm{pp}} + \Sigma_{\mathrm{pp}})\big] - \tfrac12$, the mean-spread $S$ (read off as the ensemble variance of each conditional mean) added to the conditional floor $\Sigma$:

```wl
nDrainCool[gg_] := Module[{ens = Table[meanRunCool[gg, s], {s, nEnsCool}]},
   Table[(Variance[ens[[All, j, 1]]] + vxCool + Variance[ens[[All, j, 2]]] + vpCool)/2 - 0.5,
    {j, 1, nStepCool + 1, 25}]];
```

Track it with the loop open (gain $0$) and closed (gain $G$):

```wl
{coolOff, coolOn} = {nDrainCool[0.], nDrainCool[Gcool]};
```

Now visualize the cooling: the occupation with the loop off (heating above the bath) and on (draining toward, but not reaching, the conditional floor):

```wl
ListLinePlot[{Transpose[{timesCool[[1 ;; ;; 25]], coolOff}], Transpose[{timesCool[[1 ;; ;; 25]], coolOn}]},
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"loop open (watching only)", "loop closed (feedback)"},
 FrameLabel -> {"time", "phonon number"}, PlotLabel -> "feedback cools; watching alone heats",
 Epilog -> {Directive[ColorData[97][3], Dashed], Line[{{0, nCool}, {Last@timesCool, nCool}}],
   Directive[ColorData[97][4], Dashed], Line[{{0, condFloorCool}, {Last@timesCool, condFloorCool}}]}]
```

The mean phonon number $n=\langle\hat b^\dagger\hat b\rangle$ is the oscillator's stored energy in quanta of $\omega$, so a higher curve is a hotter oscillator and a lower one colder; cooling drives $n$ down toward the ground state. The open loop climbs above the thermal line (green); the closed loop drains fast and settles well below it, above the conditional floor (red).

The occupation plot shows that the loop cools, but not how. Take it apart, from a single record to the whole ensemble.

Start with one record seen twice. The quantity that separates watching from cooling is the energy carried by the center of the estimate, $E_{\mathrm{center}} = \tfrac12(\langle\hat x\rangle^2 + \langle\hat p\rangle^2)$, formed step by step from a run:

```wl
centerEnergy[run_] := ((#[[1]]^2 + #[[2]]^2)/2 &) /@ run;
```

A fixed seed draws the same starting point and the same noise kicks at any gain, so running one seed with the loop open and then closed holds the disturbance fixed and leaves the feedback as the only difference between the two histories:

```wl
{recOff, recOn} = {meanRunCool[0., 4], meanRunCool[Gcool, 4]};
```

Lay the two center energies on one axis, the open record against the closed one:

```wl
ListLinePlot[{Transpose[{timesCool, centerEnergy[recOff]}], Transpose[{timesCool, centerEnergy[recOn]}]},
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, Frame -> True, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"loop open (same noise)", "loop closed (same noise)"},
 FrameLabel -> {"time", "center energy"}, PlotLabel -> "same record: feedback drains the estimate's energy"]
```

Both traces start together. The open one wanders and stays high, since the bath damps the mean only at $\gamma/2$, far too slowly to remove the energy the noise pours in; the closed one, under the identical noise, is pulled down within a couple of periods and holds low. Same disturbance, opposite fate: the loop turns the record into a force that carries energy out.

How much energy? The feedback force $-G\langle\hat p\rangle$ opposes the velocity, so it does work at rate $G\langle\hat p\rangle^2$. Accumulate that rate along the closed record for the energy removed by time $t$, $E_{\mathrm{removed}}(t) = G\int_0^t \langle\hat p\rangle^2\,ds$:

```wl
extractedWork[run_, gg_] := 
  Prepend[gg dtCool Accumulate[(#[[2]]^2 &) /@ Most[run]], 0.];
```

Plot that running total along the closed record:

```wl
ListLinePlot[Transpose[{timesCool, extractedWork[recOn, Gcool]}], AspectRatio -> 1/2,
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"time", "energy removed by feedback"}, PlotLabel -> "cooling is a steady drain, not a one-time act"]
```

The curve climbs and never flattens, and the work it reports far outruns the few units the estimate's energy actually fell. That gap is the mechanism: measurement back-action and the bath keep pouring energy in, and the feedback keeps taking it back out. Once the oscillator is cold the loop does not stop; it runs on to hold it there, the shock absorber built from a measurement, doing work forever against the noise.

One record is one draw; cooling is a statement about the whole ensemble. From each of many records at a fixed gain, take the endpoint of the estimate's center, $(\langle\hat x\rangle, \langle\hat p\rangle)$:

```wl
cloudCool[gg_] := (meanRunCool[gg, #][[-1]] &) /@ Range[nEnsCool];
```

Scatter that cloud of endpoints, the loop open against closed:

```wl
ListPlot[{cloudCool[0.], cloudCool[Gcool]}, PlotStyle -> Opacity[0.5],
 PlotLegends -> {"loop open", "loop closed"}, AspectRatio -> 1, Frame -> True, GridLines -> Automatic,
 PlotRange -> All, FrameLabel -> {"<x>", "<p>"}, PlotLabel -> "feedback pulls the cloud of estimates in"]
```

The open-loop centers sprawl across phase space, still spreading because nothing damps them; the closed-loop centers collapse into a tight knot at the origin. That knot is the mean-spread $S$, the covariance of where the estimate sits across records, and feedback is what pulls it in.

The knot is only part of the uncertainty. What an outside observer sees is that spread $S$ blurred by the conditional floor $\Sigma$, how sharply the filter knows each single history, the two adding to the total covariance $\Sigma + S$. Assemble the conditional floor $\Sigma$ from the three steady entries found above:

```wl
condCovCool = {{vxCool, vcCool}, {vcCool, vpCool}};
```

The mean-spread $S$ solves the Lyapunov balance of the previous block, so it is exact, with no sampling:

```wl
spreadCool[gg_] := 
  LyapunovSolve[{{-gammaCool/2, 
     1}, {-1, -gammaCool/2 - gg}}, -{{kxCool^2, 
      kxCool kpCool}, {kxCool kpCool, kpCool^2}}];
```

Draw the one-sigma ellipse of the floor $\Sigma$ alone and of the total $\Sigma + S$ with the loop open and closed:

```wl
ListLinePlot[{covarianceEllipse[{0, 0}, condCovCool + spreadCool[0.]],
   covarianceEllipse[{0, 0}, condCovCool + spreadCool[Gcool]], 
  covarianceEllipse[{0, 0}, condCovCool]}, Frame -> True, 
 GridLines -> Automatic, AspectRatio -> 1, PlotRange -> All,
   FrameLabel -> {"<x>", "<p>"}, 
 PlotLabel -> "feedback shrinks the total toward the fixed floor",
 PlotLegends -> {"total \[CapitalSigma]+S, open-loop", 
   "total \[CapitalSigma]+S, closed-loop", 
   "conditional floor \[CapitalSigma]"}]
```

The two total ellipses are the whole state the world sees; the third is the conditional floor $\Sigma$ by itself, the uncertainty the filter keeps on every single history. Closing the loop shrinks the total from the wide open-loop ellipse to the much smaller closed-loop one, but even that halts at the floor and cannot cross it. The floor itself is identical open and closed: feedback is a Hamiltonian, so it moves where the estimate sits, never how sharply the filter knows it, shrinking $S$ alone and leaving $\Sigma$ fixed. The gap between the closed-loop ellipse and the floor is the residual the next plot minimizes over the gain; that floor is the cold that no gain can beat.

How hard should the loop push? Sweep the gain and read the steady occupation off the exact Lyapunov balance:

```wl
gainsCool = Range[0, 12, 0.1];
ListLinePlot[Transpose[{gainsCool, nUncCool /@ gainsCool}], Frame -> True,
 GridLines -> Automatic, PlotRange -> All, FrameLabel -> {"feedback gain G", "steady phonon number"},
 PlotLabel -> "an optimal gain: too little leaves it warm, too much feeds noise back",
 Epilog -> {Directive[ColorData[97][3], Dashed], Line[{{0, nCool}, {Last@gainsCool, nCool}}],
   Directive[ColorData[97][4], Dashed], Line[{{0, condFloorCool}, {Last@gainsCool, condFloorCool}}]}]
```

The occupation falls steeply from the heated open-loop value, reaches a minimum near the gain chosen above, then rises again. Beyond the optimum the loop feeds the estimate's own measurement noise back onto the oscillator: this is the imprecision-back-action tradeoff, between the shot-noise blur on the record (the *imprecision*) and the momentum kick that reading the position costs (the back-action), and it caps every real feedback-cooling experiment. It never reaches the conditional floor (red): that would need a noiseless actuator, the element that applies the feedback force, acting on a perfectly known state.

Finally, a better detector lowers both the floor and the achievable cold. Repeat the floor-plus-spread Lyapunov balance at efficiency $\eta$, where the conditional covariance sharpens and the gain becomes $2\sqrt{\eta k}\,\Sigma$, holding $G$ fixed at its ideal-detector value (not re-optimized per $\eta$, so the cooled curve is conservative at low $\eta$):

```wl
nUncCoolEta[\[Eta]_, gg_] := Module[{c = steadyCovCool[\[Eta]], kx, kp, s},
   kx = 2 Sqrt[\[Eta] kCool] c[[1]]; kp = 2 Sqrt[\[Eta] kCool] c[[2]];
   s = LyapunovSolve[{{-gammaCool/2, 1}, {-1, -gammaCool/2 - gg}}, -{{kx^2, kx kp}, {kx kp, kp^2}}];
   (s[[1, 1]] + c[[1]] + s[[2, 2]] + c[[3]])/2 - 0.5];
```

Sweep the efficiency from a poor detector, $\eta = 0.2$, to a perfect one:

```wl
etasCool = Range[0.2, 1, 0.01];
```

Plot the conditional floor and the fixed-gain cooled occupation against efficiency:

```wl
ListLinePlot[{Transpose[{etasCool, (#[[1]] + #[[3]])/2 - 0.5 & /@ (steadyCovCool /@ etasCool)}],
   Transpose[{etasCool, nUncCoolEta[#, Gcool] & /@ etasCool}]},
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 PlotLegends -> {"conditional floor", "cooled occupation (fixed gain)"},
 FrameLabel -> {"detector efficiency", "phonon number"}, PlotLabel -> "a better detector cools colder"]
```

Both curves fall as the detector improves, and the gap between them is the price of noisy actuation. This is measurement-based cooling of a nanomechanical drum or a levitated nanoparticle toward its ground state, not by contact with something colder but by watching precisely and pushing back: the same conditional covariance that told the Kalman filter how well it could track now sets how cold the feedback can pull.

## Part Seven: Two Bookkeepings, and the Colour of the Light

### Linear vs Nonlinear: One Trajectory, Two Descriptions

**The problem.** One diffusive history, written under two probability measures. Under a reference Wiener measure $Q$, the *linear* stochastic Schrodinger equation runs an unnormalized state under plain noise,
$$d\tilde\psi = K\,\tilde\psi\,dt + \sum_j \hat R_j\,\tilde\psi\,dW_j^{(0)},\qquad K = -i\hat H - \tfrac12\sum_j \hat R_j^\dagger \hat R_j.$$
Its squared norm is the Radon-Nikodym likelihood density $dP/dQ$ of the physical record measure relative to $Q$ up to that time, not an absolute probability for one continuous record. Under the physical measure $P$, the normalized master equation is driven by the innovation $d\widehat W_j=dJ_j-v_jdt$ of the observed increment $dJ_j$,
$$d\rho = \mathcal{L}\rho\,dt + \sum_j(\hat R_j\rho + \rho \hat R_j^\dagger - v_j\rho)\,d\widehat W_j,\qquad v_j = 2\,\mathrm{Re}\,\mathrm{Tr}[\hat R_j\rho].$$
The useful fact is that the reference-measure average of the unnormalized projectors is the master-equation state. Build the linear version and check.

One diffusive trajectory has two equivalent descriptions. The **nonlinear** SME (used so far) carries the normalized state and renormalizes each step. The **linear** SSE carries an unnormalized state under the zero-drift reference noise $dW^{(0)}$, the noise of a reference measure $Q$ that reads every record as signal-free noise, and never renormalizes,
$$d\left|\tilde\psi\right\rangle = (-i\hat H - \tfrac12\hat c^\dagger\hat c)\left|\tilde\psi\right\rangle dt + \hat c\left|\tilde\psi\right\rangle dJ,\qquad dJ = dW^{(0)}\ \text{under}\ Q.$$
Two facts turn this single equation into the section. Its squared norm $w=\langle\tilde\psi|\tilde\psi\rangle$ is the likelihood ratio $dP/dQ$ that reweights $Q$ into the physical measure $P$; and the plain $Q$-average of the unnormalized projectors $|\tilde\psi\rangle\langle\tilde\psi|$ is the master-equation state $\bar\rho$, carried with no reweighting. Take the average first, the weight second.

The average is two lines of Ito calculus. With $\tilde\rho=|\tilde\psi\rangle\langle\tilde\psi|$ and $dJ^2=dt$, the drift terms collect into the Lindbladian and the noise into one trace-free term,
$$d\tilde\rho = \mathcal L\tilde\rho\,dt + (\hat c\tilde\rho + \tilde\rho\hat c^\dagger)\,dW^{(0)}.$$
Under $Q$ the increment $dW^{(0)}$ has zero mean, so averaging erases the stochastic term and the mean unnormalized state obeys the master equation on its own, $\tfrac{d}{dt}\mathbb E_Q[\tilde\rho]=\mathcal L\,\mathbb E_Q[\tilde\rho]$, hence $\mathbb E_Q[\tilde\rho(t)]=\bar\rho(t)$. Fix the atom and the non-Hermitian drift $K = -i\hat H - \tfrac12\hat c^\dagger\hat c$:

```wl
\[CapitalOmega]lin = 2.; cLin = lower; Knh = -I (\[CapitalOmega]lin/2) X - ConjugateTranspose[cLin] . cLin/2;
step = 0.005; steps = 600;
```

One unnormalized run under plain noise, never renormalizing:

```wl
linearRun[seed_] := BlockRandom[SeedRandom[seed];
    Fold[#1 + Knh . #1 step + cLin . #1 #2 &, excited, RandomVariate[NormalDistribution[0, Sqrt[step]], steps]]];
```

Run four hundred of them:

```wl
linearStates = linearRun /@ Range[400];
```

Average the unnormalized projectors, each already carrying its weight $w=\|\tilde\psi\|^2$ in its trace, so no separate reweighting is applied:

```wl
linearAverage = Mean[Outer[Times, #, Conjugate[#]] & /@ linearStates];
```

Measure the gap to the master equation:

```wl
Max@Abs@Flatten[linearAverage - evolve[(\[CapitalOmega]lin/2) X, {cLin}, densityMatrix[excited], steps step]]
```

The unweighted reference-measure average lands on the master equation within Monte-Carlo scatter and Euler time-step error, the numerical face of $\mathbb E_Q[\tilde\rho]=\bar\rho$, reached with no pathwise renormalization.

Now the weight. The squared norm $w=\langle\tilde\psi|\tilde\psi\rangle$ is not conserved, and its own increment is where the likelihood sits. The non-Hermitian drift is built so that $K^\dagger+K+\hat c^\dagger\hat c=0$, which cancels every $dt$ term in $dw$ and leaves the norm driven purely by the record,
$$dw = w\,v\,dJ,\qquad v = \langle\hat c+\hat c^\dagger\rangle_{\rho_c} = 2\,\mathrm{Re}\,\mathrm{Tr}[\hat c\rho_c],$$
with $\rho_c=\tilde\rho/w$ the normalized state. Under $Q$ this drift-free increment makes $w$ a mean-one martingale, $\mathbb E_Q[w]=1$: exactly the statement that $w=dP/dQ$ is a probability reweighting. Its logarithm accumulates the running log-likelihood,
$$d\log w = v\,dJ - \tfrac12 v^2\,dt,$$
climbing when a record increment aligns with the state's predicted signal $v$ and falling when it contradicts it. The finite seeded Euler sample of $w$ should sit near one:

```wl
Mean[Re[Conjugate[#] . #] & /@ linearStates]
```

Different reference-noise histories carry different likelihood-ratio weights.

Now visualize the distribution of their likelihood weights:

```wl
Histogram[Re[Conjugate[#] . #] & /@ linearStates, 30, "PDF", Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"likelihood-ratio weight dP/dQ", "density"},
 PlotLabel -> "reference histories carry likelihood-ratio weights"]
```

These weights organize the section into a single identity. The unnormalized projector is the weighted normalized state, $\tilde\rho = w\rho_c$, so the unweighted $Q$-average of $\tilde\rho$, the weighted $Q$-average of $\rho_c$, and the physical $P$-average of $\rho_c$ are one and the same master-equation state,
$$\mathbb E_Q[\tilde\rho] \;=\; \mathbb E_Q[w\,\rho_c] \;=\; \mathbb E_P[\rho_c] \;=\; \bar\rho.$$
The middle equality is the change of measure, $\mathbb E_Q[w\,X]=\mathbb E_P[X]$ when $w=dP/dQ$, and it is where the weight earns its keep. Drop it, average the same normalized states under $Q$ with equal weight, and the result is not $\bar\rho$, because $Q$ draws records at the wrong physical frequencies. Put the three averages side by side at the final time:

```wl
With[{ws = Re[Conjugate[#] . #] & /@ linearStates,
   tildes = Outer[Times, #, Conjugate[#]] & /@ linearStates,
   ref = evolve[(\[CapitalOmega]lin/2) X, {cLin}, densityMatrix[excited], steps step]},
 {rows = {{"unnormalized projectors (Q)", Mean[tildes]},
     {"likelihood-weighted (Q)", Total[tildes]/Total[ws]},
     {"equal-weight normalized (Q)", Mean[MapThread[#1/#2 &, {tildes, ws}]]}}},
 Grid[Prepend[{#[[1]], Re[#[[2, 1, 1]]], Norm[#[[2]] - ref, "Frobenius"]} & /@ rows,
    {"reference-measure average", "excited population", "gap to master equation"}],
  Frame -> All, Alignment -> Left]]
```

The two weighted rows sit on the master equation; the equal-weight row misses it, and by a bias no larger ensemble will close, since the error is in the sampling law and not the sample count. That gap is the whole content of the weight: it is what makes the reference measure physical. Two readings of one diffusive process follow from it: normalize the linear state and drive it with the physical innovation for the nonlinear state an experimenter tracks, or carry the unnormalized state under reference noise and keep its norm as the $dP/dQ$ factor.

One record carries both bookkeepings at once. Roll out a single linear run keeping its whole history:

```wl
linPathA5[seed_] := BlockRandom[SeedRandom[seed];
    FoldList[#1 + Knh . #1 step + cLin . #1 #2 &, excited, RandomVariate[NormalDistribution[0, Sqrt[step]], steps]]];
```

Take one such run for a fixed seed:

```wl
oneA5 = linPathA5[7];
```

and its time grid:

```wl
gridA5 = step Range[0, steps];
```

Read the two descriptions off that one run: the normalized state (the state direction the physical trajectory takes for this record) and the log-weight it accumulates (its running log-likelihood), side by side:

```wl
Row[{
  ListLinePlot[Transpose[{gridA5, (1 + Re[Conjugate[#] . Z . #])/2 &[Normalize[#]] & /@ oneA5}],
   Frame -> True, GridLines -> Automatic, PlotRange -> {{0, steps step}, {0, 1.05}}, ImageSize -> Medium,
   FrameLabel -> {"time", "excited population (normalized)"}, PlotLabel -> "the state you track"], "  ",
  ListLinePlot[Transpose[{gridA5, Log[Re[Conjugate[#] . #]] & /@ oneA5}],
   Frame -> True, GridLines -> Automatic, PlotRange -> {{0, steps step}, All}, ImageSize -> Medium,
   FrameLabel -> {"time", "log weight  log ||psi||^2"}, PlotLabel -> "the likelihood it carries"]}]
```

The same record has a normalized-state reading and a running weight. The normalized state is the physical trajectory an experimenter tracks; the weight is the extra likelihood factor, needed only to turn averages taken under the reference measure $Q$ into physical ones.

The catch of the linear method is that the weights spread, so a handful of records come to dominate the average, the importance-sampling degeneracy. Measure it with the effective sample size $N_{\mathrm{eff}} = (\sum w)^2/\sum w^2$. Run an ensemble of linear paths:

```wl
essPathsA5 = linPathA5 /@ Range[400];
```

Read $N_{\mathrm{eff}}/N$ from the squared-norm weights at sampled times:

```wl
essA5 = Table[With[{w = (Re[Conjugate[#] . #] &@#[[j]]) & /@ essPathsA5}, (Total[w]^2/Total[w^2])/400.], {j, 1, steps + 1, 20}];
```

and plot its decay over time:

```wl
ListLinePlot[Transpose[{gridA5[[1 ;; ;; 20]], essA5}], Frame -> True, GridLines -> Automatic,
 PlotRange -> All, FrameLabel -> {"time", "effective sample size / N"},
 PlotLabel -> "the linear weights degenerate: fewer records carry the average"]
```

The effective fraction falls from one toward a fraction of the ensemble: as time goes on, the unweighted linear average is carried by ever fewer high-weight records, and its Monte-Carlo error grows even though it remains unbiased. That is the price of the linear picture's convenience, and the reason the nonlinear (renormalized) equation, which keeps every trajectory equally weighted, is what a long real-time filter uses.

### The Driven Atom, Averaged and Watched

**The problem.** One driven, damped, emitting atom, carried in full generality and then watched three ways. Its unconditional master equation keeps detuning, a thermal drive, and pure dephasing together,
$$\mathcal{L}\rho = -i[\hat H,\rho] + \gamma(n_T+1)\,\mathcal{D}[\hat\sigma_-]\rho + \gamma\,n_T\,\mathcal{D}[\hat\sigma_+]\rho + \gamma\,k_d\,\mathcal{D}[\hat\sigma_z]\rho,\qquad \hat H = \tfrac{\Delta\omega}{2}\hat\sigma_z + \tfrac{\Omega}{2}\hat\sigma_x,$$
and relaxes to a mixed steady state inside the ball. Three separate environmental channels carry information away, one dissipator each: the downward channel $\hat R_\downarrow = \sqrt{\gamma(n_T+1)}\,\hat\sigma_-$ is excited-to-ground emission, the upward channel $\hat R_\uparrow = \sqrt{\gamma n_T}\,\hat\sigma_+$ is thermal absorption from the bath, and the dephasing channel $\hat R_\phi = \sqrt{\gamma k_d}\,\hat\sigma_z$ scrambles the phase while moving no population. Detuning is different in kind: $\tfrac{\Delta\omega}{2}\hat\sigma_z$ is Hamiltonian, so it turns the state without opening a fourth channel and without mixing it. We hold all three channels active throughout and vary only how much of their output the observer keeps: nothing, the fluorescence alone, or all three ideally. These are three levels of knowledge about one atom, not three physical models, running from the mixed interior point (no record) through mixed conditioned states (partial record) to pure conditioned states (complete record), every one sharing the single unconditional average $\mathbb{E}[\rho_c(t)] = \rho(t)$. Heterodyne split one emitted field into two quadratures; here three distinct mechanisms carry three independent records, and the question is how many of them the observer keeps.

Start with the level we never record, the unconditional steady state, the mixed point every watched trajectory below will circle. Two rates organize it: the longitudinal rate $\Gamma_1 = \gamma(2n_T + 1)$ at which populations relax, and the transverse rate $\Gamma_2 = \tfrac12\Gamma_1 + 2\gamma k_d$ at which coherences decay, faster than $\Gamma_1/2$ by the extra dephasing. Define the steady Bloch vector in closed form:

```wl
steadyGeneral[\[CapitalDelta]_, \[CapitalOmega]_, \[Gamma]_, nT_, kd_] := With[{G1 = \[Gamma] (2 nT + 1), G2 = \[Gamma] (2 nT + 1)/2 + 2 \[Gamma] kd},
   With[{den = G1 (G2^2 + \[CapitalDelta]^2) + \[CapitalOmega]^2 G2}, {-\[Gamma] \[CapitalOmega] \[CapitalDelta], \[Gamma] \[CapitalOmega] G2, -\[Gamma] (G2^2 + \[CapitalDelta]^2)}/den]];
```

Check it against the toolkit's steady-state finder with detuning, temperature, and dephasing all switched on:

```wl
Max@Abs[steadyGeneral[1.3, 3., 1., 0.4, 0.25] - blochVector[First@stationary[(1.3/2) Z + (3./2) X, {Sqrt[1. (0.4 + 1)] lower, Sqrt[1. 0.4] ConjugateTranspose[lower], Sqrt[1. 0.25] Z}]]]
```

The discrepancy is at numerical precision: the closed form is exact for the stated finite-dimensional Bloch equations. The detuning enters through $\Gamma_2^2 + \Delta\omega^2$ in the denominator, so the resonance has a Lorentzian denominator of width $\Gamma_2$ that the drive power-broadens through the $\Omega^2\Gamma_2$ term; temperature lifts the ground-state bias, since at $\Omega = 0$ the state rests at $z_{\mathrm{eq}} = -1/(2n_T + 1)$ rather than the south pole. Plot the steady excited population $\tfrac12(1 + z)$ against detuning, weak drive against strong:

```wl
ListLinePlot[Table[Table[{d, (1 + steadyGeneral[d, om, 1., 0., 0.][[3]])/2}, {d, -8, 8, 0.05}], {om, {1., 4.}}],
 PlotStyle -> {ColorData[97, 1], ColorData[97, 2]}, PlotLegends -> {"\[CapitalOmega] = 1", "\[CapitalOmega] = 4"},
 Frame -> True, GridLines -> Automatic, PlotRange -> All,
 FrameLabel -> {"detuning \[CapitalDelta]\[Omega]", "steady excited population"},
 PlotLabel -> "resonance is a power-broadened Lorentzian in detuning"]
```

The excited population peaks on resonance and falls off as a Lorentzian to either side; the stronger drive both lifts the peak toward saturation at one half and widens it, the power broadening. That is the general averaged atom, the no-record level: with the detector ignored, the state relaxes deterministically to one mixed interior point. To watch it we must say what a detector does with each of the three outputs, and that choice, not the Liouvillian, is what picks a trajectory out of the average.

Now the same atom, watched. Fix the drive $\Omega$, the decay rate $\gamma$, the detuning $\Delta\omega$, the bath occupation $n_T$, and the dephasing $k_d$:

```wl
\[CapitalOmega]dr = 3.; \[Gamma]dr = 1.; \[CapitalDelta]dr = 0.8; nTdr = 0.3; kdDr = 0.2;
```

Build the Hamiltonian, a detuning about $\hat\sigma_z$ plus a drive about $\hat\sigma_x$:

```wl
Hdr = (\[CapitalDelta]dr/2) Z + (\[CapitalOmega]dr/2) X;
```

Build the three channel operators from the toolkit's `lower`: emission down, thermal absorption up, and dephasing about $\hat\sigma_z$:

```wl
Rdown = Sqrt[\[Gamma]dr (nTdr + 1)] lower; Rup = Sqrt[\[Gamma]dr nTdr] ConjugateTranspose[lower]; Rphi = Sqrt[\[Gamma]dr kdDr] Z;
```

Collect the three as the leak list:

```wl
drLeaks = {Rdown, Rup, Rphi};
```

With the detector ignored, the state relaxes to the mixed interior point $\vec a_{\mathrm{ss}}$ of the closed form above; read its length, which a pure state would carry as one:

```wl
settled = steadyGeneral[\[CapitalDelta]dr, \[CapitalOmega]dr, \[Gamma]dr, nTdr, kdDr];
Norm[settled]
```

Well short of one: with no record kept, the atom is mixed, sitting inside the ball where drive, emission, thermal absorption, and dephasing balance.

To watch it, we must put a detector on each output, and that choice, not $\mathcal{L}$, is what selects a trajectory. Read all three channels by phase-zero homodyne, each carrying its own independent noise $dW_\alpha$ (with $dW_\alpha\,dW_\beta = \delta_{\alpha\beta}\,dt$), so each record is $dJ_\alpha = \langle\hat R_\alpha + \hat R_\alpha^\dagger\rangle_c\,dt + dW_\alpha$. Which system observable does each carry? Since $\hat\sigma_- + \hat\sigma_+ = \hat\sigma_x$, both emission records track $\langle\hat\sigma_x\rangle$, differing only by their thermal weights, while the dephasing operator is Hermitian, so $\hat R_\phi + \hat R_\phi^\dagger = 2\sqrt{\gamma k_d}\,\hat\sigma_z$ and its record tracks $\langle\hat\sigma_z\rangle$. Confirm the three record drifts are exactly those combinations:

```wl
With[{probe = densityMatrix[Normalize[{1., 0.5 + 0.5 I}]]}, {a = blochVector[probe]},
 Chop[(Re@Tr[(# + ConjugateTranspose[#]) . probe] & /@ {Rdown, Rup, Rphi}) -
   {Sqrt[\[Gamma]dr (nTdr + 1)] a[[1]], Sqrt[\[Gamma]dr nTdr] a[[1]], 2 Sqrt[\[Gamma]dr kdDr] a[[3]]}]]
```

Zero in every slot: three independent records collapse onto two system observables. Two different baths both report $\langle\hat\sigma_x\rangle$, one reports $\langle\hat\sigma_z\rangle$, and $\langle\hat\sigma_y\rangle$, generally nonzero, appears in no record at all; it reaches the conditioned state only through the Hamiltonian rotation. Channel identity and measured observable are not the same thing.

Integrate the record-blind master equation once, the baseline every ensemble must return to:

```wl
refRun = evolveODE[Hdr, drLeaks, densityMatrix[excited], 8.];
```

Run one conditioned history that keeps only the fluorescence, watching $\hat R_\downarrow$ while the thermal-absorption and dephasing outputs escape:

```wl
watchedOne = trajectory[densityMatrix[excited], Hdr, {Rdown}, {1.}, {Rup, Rphi}, 0.01, 8., 5];
```

Run another that keeps all three channels ideally:

```wl
watchedAll = trajectory[densityMatrix[excited], Hdr, drLeaks, {1., 1., 1.}, {}, 0.01, 8., 5];
```

Take the shared time grid off that run:

```wl
whenDr = watchedAll["times"];
```

The fluorescence-only history obeys $d\rho_c = \mathcal{L}\rho_c\,dt + \mathcal{H}[\hat R_\downarrow]\rho_c\,dW_\downarrow$: the unwatched channels stay in $\mathcal{L}$ but condition nothing. Read the smallest purity each run reaches:

```wl
{Min[purity /@ watchedOne["states"]], Min[purity /@ watchedAll["states"]]}
```

The fluorescence-only run drops well below one, genuinely mixed; the fully watched run holds at one. Temperature and dephasing are equally present in both, active in the same $\mathcal{L}$; the only difference is whether their information is collected. Missing two of the three channels mixes the conditioned state; catching all three keeps it pure, thermal drive and dephasing notwithstanding. That is the sharper form of "watching keeps purity": not watching anything, but holding a complete ideal record of every channel that carries information away.

See it as three curves in time, the record-blind purity against the two conditioned runs:

```wl
ListLinePlot[{Transpose[{whenDr, purity[refRun[#]] & /@ whenDr}],
   Transpose[{whenDr, purity /@ watchedOne["states"]}],
   Transpose[{whenDr, purity /@ watchedAll["states"]}]},
 PlotLegends -> {"no record (unconditional)", "fluorescence only", "all three channels"},
 Frame -> True, GridLines -> Automatic, PlotRange -> {0, 1.1},
 FrameLabel -> {"time", "purity"}, PlotLabel -> "only the complete record keeps a conditioned run pure"]
```

The unconditional purity slides smoothly from one to its mixed steady value; the fluorescence-only run jitters below it, genuinely mixed; the fully watched run rides the top at one. Purity is a property of the record kept, not of the atom. The same three levels read on the Bloch sphere, the pure surface path and the mixed interior path circling the one record-blind steady point:

```wl
blochPlot[{blochVector /@ watchedAll["states"], blochVector /@ watchedOne["states"]},
 "complete record stays on the surface; fluorescence only wanders inside",
 {"all three channels (pure)", "fluorescence only (mixed)"},
 {Black, PointSize[0.035], Point[settled]}]
```

The fully watched trajectory stays on the surface, pure, circling the black interior point; the fluorescence-only trajectory wanders through the interior, mixed. Both must average back to that same black point, the atom with no record at all. Check it: build an ensemble of each monitoring and compare its mean to the unconditional reference:

```wl
crowdOne = Table[trajectory[densityMatrix[excited], Hdr, {Rdown}, {1.}, {Rup, Rphi}, 0.01, 8., k]["states"], {k, 100}];
crowdAll = Table[trajectory[densityMatrix[excited], Hdr, drLeaks, {1., 1., 1.}, {}, 0.01, 8., k]["states"], {k, 100}];
With[{ref = refRun /@ whenDr},
 {Max@MapThread[Norm[#1 - #2, "Frobenius"] &, {Mean[crowdOne], ref}],
  Max@MapThread[Norm[#1 - #2, "Frobenius"] &, {Mean[crowdAll], ref}]}]
```

Both gaps sit at the hundred-trajectory sampling level: averaging the record away returns the same $\rho(t)$ whether the conditioned states were mixed or pure. The two knobs from the toolkit still set the residue, more trajectories against the scatter and a finer step against the bias.

One honest caveat on the pure fully watched run. Writing the thermal bath as the two operators $\sqrt{\gamma(n_T+1)}\,\hat\sigma_-$ and $\sqrt{\gamma n_T}\,\hat\sigma_+$ dilates it into two effective vacuum channels, and monitoring both of those is what keeps the trajectory pure. A real thermal reservoir is itself mixed, and ordinary access to its physical output need not deliver the two clean records this idealization assumes. Complete monitoring of all three effective channels is an information-complete unravelling: it shows what full environmental knowledge would imply, not necessarily what a laboratory detector on a warm bath can reach.

So one Liouvillian carries three answers. The black interior point is the atom with the record thrown away; the interior stochastic path is the same atom watched incompletely; the surface stochastic path is the atom watched in full. None is more true than the others, and all three average to the first: they answer three different questions about one open system.

### The Mollow Triplet: Resonance-Fluorescence Spectrum

**The problem.** The spectrum (power spectrum) of the light a driven atom emits, computed by quantum regression from the two-time correlations of the master equation. The atom is the resonant emitter of the previous section,
$$\dot\rho = -i[\tfrac{\Omega}{2}\hat\sigma_x,\rho] + \gamma\,\mathcal{D}[\hat\sigma_-]\rho.$$
Its light divides into two parts: an *elastic* (coherent) line at the drive frequency, a delta function carrying the weight $\left|\langle\hat\sigma_-\rangle_{\mathrm{ss}}\right|^2$, and the *inelastic* (fluorescence) spectrum of the fluctuations $\delta\hat\sigma_- = \hat\sigma_- - \langle\hat\sigma_-\rangle_{\mathrm{ss}}$, given by quantum regression as the steady-state correlation
$$S_{\mathrm{inel}}(\mu) \propto \mathrm{Re}\int_0^\infty e^{-i\mu t}\,\langle\delta\hat\sigma_+(t)\,\delta\hat\sigma_-(0)\rangle_{\mathrm{ss}}\,dt.$$
Regression runs this correlation on the same Liouvillian $\mathcal{L}$ that evolves the state, so Fourier transform turns the propagator $e^{\mathcal{L}t}$ into the resolvent $(i\mu-\mathcal{L})^{-1}$. For weak driving the inelastic part is a single line; for strong driving $\Omega > \gamma$ it splits into the three-peaked *Mollow triplet*, a central line and two sidebands at $\pm\Omega$, the dynamical Stark effect.

Homodyne detection instead reads a chosen quadrature, and its spectrum can dip *below* the shot-noise floor: squeezed, nonclassical light. Build both from the same big matrix.

One last set of quantities from the master equation, all through the eigen-based spectrum tool built for the quantum point contact (`regressionSpectrum`). It needs no shifted frequency grid, because the steady mode is dropped explicitly. The inelastic fluorescence spectrum is its correlation of $\delta\hat\sigma_-$:

```wl
mollowInelastic[\[CapitalOmega]_, \[Gamma]_, tones_] := With[{big = liouvillian[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}, 2],
    rss = First@stationary[(\[CapitalOmega]/2) X, {Sqrt[\[Gamma]] lower}]},
   With[{dm = lower - Tr[lower . rss] id2}, regressionSpectrum[big, dm . rss, ConjugateTranspose[dm], tones]]];
```

Sweep frequencies (centered on zero, no offset needed) for a weak and a strong drive:

```wl
tones = Range[-9, 9, 0.05];
{weak, strong} = {mollowInelastic[0.6, 1., tones], mollowInelastic[6., 1., tones]};
```

Weakly driven, the inelastic spectrum is a single line at the drive frequency. Strongly driven, it splits into the **Mollow triplet**: a central line and two sidebands at $\pm\Omega$.

Now visualize the weak- and strong-drive inelastic spectra together. Their absolute scales differ by orders of magnitude (the weak line would vanish on the strong axis), so plot each normalized to its own peak to compare the *shapes*, a single line against the triplet:

```wl
ListLinePlot[{Transpose[{tones, weak/Max[weak]}], Transpose[{tones, strong/Max[strong]}]},
 Frame -> True, GridLines -> {{-6, 0, 6}, None},
 PlotLegends -> {"weak drive", "strong drive"}, PlotRange -> All,
 FrameLabel -> {"frequency (relative to drive)", "inelastic spectrum (each / its peak)"},
 PlotLabel -> "strong driving splits the inelastic spectrum into the Mollow triplet"]
```

The sidebands sit near $\pm\Omega$: under strong driving the atom's levels are dressed into a four-level ladder, whose transitions produce the three lines. Confirm the sidebands land at the drive strength:

```wl
Select[tones[[#]] & /@ FindPeaks[strong][[All, 1]], Abs[Abs[#] - 6] < 0.6 &]
```

The peaks sit at $\pm\Omega$. Confirm too that the inelastic spectrum is nowhere negative, as a genuine power spectrum must be:

```wl
Min[strong]
```

Nonnegative within numerical precision.

The elastic line lives separately, a delta at zero frequency whose weight is the coherent dipole; it dominates at weak drive and all but vanishes at strong drive as the light turns fully inelastic:

```wl
elasticWeight[\[CapitalOmega]_] := With[{rss = First@stationary[(\[CapitalOmega]/2) X, {Sqrt[1.] lower}]}, Abs[Tr[lower . rss]]^2];
```

At weak drive the coherent weight is substantial:

```wl
elasticWeight[0.6]
```

and at strong drive it all but vanishes:

```wl
elasticWeight[6.]
```

The coherent weight falls by more than an order of magnitude from weak to strong drive. To watch the triplet being born, sweep the drive continuously:

```wl
mollowDrives = Range[0.3, 8, 0.25];
```

Stack the inelastic spectra into a heatmap, each row scaled to its own peak so the sidebands stay visible against the central peak that dwarfs them:

```wl
mollowHeat = (#/Max[#]) & /@ (mollowInelastic[#, 1., Range[-11, 11, 0.1]] & /@ mollowDrives);
```

Plot it as a heatmap, weak drive at the bottom rising to strong at the top:

```wl
ArrayPlot[mollowHeat, DataReversed -> True,
 ColorFunction -> "SunsetColors", Frame -> True, AspectRatio -> 1/2,
 FrameLabel -> {"frequency", "drive \[CapitalOmega]"},
 PlotLabel -> "the Mollow triplet is born as the drive grows (rows scaled to peak)"]
```

The bright ridge splits into three, the outer two fanning apart as the drive grows: the dressed-state sidebands at $\mu = \pm\Omega$. The triplet is what an ordinary spectrometer sees.

A homodyne detector sees something an intensity spectrometer cannot: it measures one quadrature $\hat\sigma_- e^{-i\theta} + \hat\sigma_+ e^{i\theta}$, and its spectrum, normalized so the shot-noise floor is one, can fall below that floor. The weak drive we examine, $\Omega = \gamma/4$, is exactly an exceptional point of this Liouvillian, the same coalescence the quantum point contact met: its two dressed modes merge there, so we read this spectrum with the mode-free `resolventSpectrum` rather than the eigenmode sum, which would divide by a singular eigenbasis. Build it from the measurement current fluctuations:

```wl
mollowHom[\[CapitalOmega]_, \[Gamma]_, \[Theta]_, tones_] := Module[{c = Sqrt[\[Gamma]] lower, big, rss, cth, mq, mean},
   big = liouvillian[(\[CapitalOmega]/2) X, {c}, 2]; rss = First@stationary[(\[CapitalOmega]/2) X, {c}];
   cth = c Exp[-I \[Theta]]; mq = cth + ConjugateTranspose[cth]; mean = Re@Tr[mq . rss];
   1 + resolventSpectrum[big, cth . rss + rss . ConjugateTranspose[cth] - mean rss, mq - mean id2, tones]];
```

At weak drive, one quadrature is squeezed and the perpendicular one is not. Compare the two at $\Omega = 0.25\gamma$:

```wl
ListLinePlot[{Transpose[{tones, mollowHom[0.25, 1., Pi/2, tones]}], Transpose[{tones, mollowHom[0.25, 1., 0., tones]}]},
 Frame -> True, GridLines -> {None, {{1, Directive[Gray, Dashed]}}},
 PlotLegends -> {"squeezed quadrature (\[Theta] = \[Pi]/2)", "anti-squeezed (\[Theta] = 0)"}, PlotRange -> {0.7, All},
 FrameLabel -> {"frequency", "homodyne spectrum / shot noise"}, PlotLabel -> "resonance fluorescence dips below the shot-noise floor"]
```

The $\theta = \pi/2$ quadrature dips below the dashed shot-noise line while its partner rises above it: the atom's fluorescence is squeezed, its noise pushed below the vacuum in one quadrature at the cost of excess in the other. Confirm the dip is real:

```wl
Min[mollowHom[0.25, 1., Pi/2, tones]]
```

Below one: sub-shot-noise, a signature of nonclassical light that only phase-sensitive (homodyne) detection can reveal, and that the necessarily-nonnegative intensity spectrum of the triplet can never show.

How far the squeezing reaches, and where it is deepest, is a question about the drive, and it is cleanest at zero frequency, the bottom of the dip we just drew. There the homodyne spectrum has the closed form $S(0) = 1 + 8\Omega^2(4\Omega^2 - 1)/(1 + 2\Omega^2)^3$ in units of $\gamma$, which drops below one exactly when $\Omega < \gamma/2$. Read $S(0)$ off the same machinery:

```wl
zeroFreqSqueeze[\[CapitalOmega]_] := mollowHom[\[CapitalOmega], 1., Pi/2, {0.}][[1]];
```

and sweep it across the drive:

```wl
ListLinePlot[Table[{om, zeroFreqSqueeze[om]}, {om, 0.05, 1.2, 0.02}],
 Frame -> True, GridLines -> {{0.5}, {{1, Directive[Gray, Dashed]}}}, PlotRange -> All,
 FrameLabel -> {"drive \[CapitalOmega]/\[Gamma]", "zero-frequency spectrum S(0)"},
 PlotLabel -> "the in-phase quadrature is squeezed only for \[CapitalOmega] < \[Gamma]/2"]
```

$S(0)$ dips below the shot-noise line on $\Omega < \gamma/2$ (the vertical gridline) and rises above it past that drive, deepest at $\Omega = \tfrac12\sqrt{3 - \sqrt7}\,\gamma$, where it reaches $(323 - 56\sqrt7)/243$. Push harder and the anti-squeezed quadrature's excess noise floods over the floor: the in-phase squeezing is a weak-drive effect, and the shallow spectral squeezing that lingers to higher drive has slid off zero frequency. Confirm the plotted curve is that closed form and not an artifact of the resolvent:

```wl
Max@Table[Abs[zeroFreqSqueeze[om] - (1 + 8 om^2 (4 om^2 - 1)/(1 + 2 om^2)^3)], {om, 0.05, 1.2, 0.05}]
```

The maximum discrepancy is at numerical precision. The whole spectroscopy, elastic line, inelastic triplet, and squeezing dip, comes from the two-time correlations of the same master equation built on the first page.

## Where This Leaves Us

From plain matrices and a few rules, we built a small laboratory for open quantum systems: the master equation for the unconditioned average, a positivity-preserving finite-step filter for a kept record, and plotters for the Bloch sphere, phase space, covariance ellipses, and spectra. On that base, twenty systems: pure dephasing, a driven damped atom, a cavity cat decohering at a rate set by the square of its separation, the Haroche experiment, dispersive qubit readout, quantum Brownian motion, the jump/homodyne/heterodyne unravellings, the quantum Zeno effect, a charge qubit read by a quantum point contact, measurement-induced collapse, a quantum Kalman filter, measurement feedback, thermalization, rapid purification, feedback cooling, the linear-vs-nonlinear trajectory, and the Mollow triplet. The recurring theme is conditionality: a complete ideal record can keep an initially pure state on a pure stochastic trajectory, whereas missed channels, inefficiency, or mixed initial data produce mixed conditional states. In every unravelling, ensemble recovery of the master equation is a joint large-sample and small-step limit. Change a number, rerun the convergence checks, and watch the physics move.