---
Template: Default
Title: "How fast does photon loss erase the Heisenberg advantage of a NOON state?"
---

# How fast does photon loss erase the Heisenberg advantage of a NOON state?

To measure an optical phase $\theta$, we send light through the arm that carries it and read the phase off the output. The precision of that reading is set by the state we send in. With $n$ independent photons, each an independent estimate of $\theta$, the error averages down as $1/\sqrt n$: this is the shot-noise limit. A NOON state, the two-mode superposition $(\lvert n,0\rangle+\lvert0,n\rangle)/\sqrt2$ of all $n$ photons in one arm or all $n$ in the other, does better. Because the two branches accumulate phase at opposite ends, the superposition winds $n$ times faster than a single photon, and the error falls as $1/n$: the Heisenberg limit.

The figure of merit behind both scalings is the quantum Fisher information $F$: the Cramer-Rao bound says the error of any unbiased estimator of $\theta$ obeys $\delta\theta\geq1/\sqrt F$ for one shot. Shot noise is $F=n$; the NOON state reaches $F=n^2$. The quadratic gain is the whole point of the state.

That gain is also fragile. A NOON state is a maximally delocalized superposition, and losing even one of its $n$ photons tells the environment which arm the light was in, which is exactly the information the interferometer was hiding. So the question of this essay is: **as photons leak out of the interferometer, how fast does the Fisher information fall, and when does the Heisenberg advantage drop below plain shot noise?** The answer, computed below, is $F=n^2\eta^n$ with $\eta$ the per-photon transmission: the $n^2$ sensitivity survives only on the branch where every photon makes it through, and that branch has weight $\eta^n$.

The lossless half of this story is the representation worked out in Chapter 6 (SU(2) interferometry in the joint Fock-polarization basis) and Chapter 3 (the two-mode NOON state and its Wigner function). The lossy half is the extension.

Load the checked-out QF paclet, then its main and bosonic constructors:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"]
```

Load the main context:

```wl
Needs["Wolfram`QuantumFramework`"]
```

Expose the second-quantization context, which carries `AnnihilationOperator` and the Fock constructors:

```wl
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

## The NOON state and its Fisher information: a two-mode cat that winds $n$ times faster

We work with two bosonic modes, the two arms of the interferometer. A NOON state of order $n$ puts all $n$ photons in the first arm or all $n$ in the second, in equal superposition. Since loss can only lower the photon number and no arm ever holds more than $n$, a Fock cutoff of $d=n+1$ levels per mode holds every state we will meet: the truncation is exact here, not an approximation.

Build the NOON state as a vector in the $d\otimes d$ two-mode space, with $\lvert n,0\rangle$ and $\lvert0,n\rangle$ the two populated basis kets:

```wl
noonState[n_] := With[{d = n + 1},
    QuantumState[(UnitVector[d^2, n d + 1] + UnitVector[d^2, n + 1])/Sqrt[2.], {d, d}]]
```

Look at the order-two state, $(\lvert2,0\rangle+\lvert0,2\rangle)/\sqrt2$:

```wl
noonState[2]
```

Confirm it is a legitimate normalized state:

```wl
noonState[2]["Norm"]
```

The phase is imprinted on one arm: a phase shift $\theta$ on the second mode acts by $e^{-i\theta\hat n_2}$, so the generator of the estimation is the number operator $\hat n_2=\hat a_2^\dagger\hat a_2$ of that mode. Build it on the two-mode space. The identity on the untouched mode is written with an explicit dimension, `QuantumOperator[IdentityMatrix[d], {k}, d]`, so the two-mode tensor product never has to guess a cutoff:

```wl
numberTwo[n_] := With[{d = n + 1},
    With[{a2 = QuantumTensorProduct[QuantumOperator[IdentityMatrix[d], {1}, d], AnnihilationOperator[d, {2}]]},
        a2["Dagger"] @ a2]]
```

For a pure state the quantum Fisher information is four times the variance of the generator, $F=4\,(\langle\hat n_2^2\rangle-\langle\hat n_2\rangle^2)$: this is exactly the $4\,\mathrm{Var}(\hat J_z)$ of the SU(2) treatment in Chapter 6, with $\hat J_z$ here specialized to the single-arm number operator. Take expectations by tracing the operator against the density matrix:

```wl
expval[op_, qs_] := Re @ Tr[op["Matrix"] . qs["DensityMatrix"]]
```

Define the pure-state Fisher information and evaluate it across $n$:

```wl
pureQFI[n_] := With[{g = numberTwo[n], qs = noonState[n]},
    4 (expval[g @ g, qs] - expval[g, qs]^2)]
```

Read $F$ for the first several NOON orders:

```wl
Table[{n, Chop @ pureQFI[n]}, {n, 1, 6}]
```

As one can see, $F=n^2$ exactly: the two branches are eigenstates of $\hat n_2$ with eigenvalues $0$ and $n$, so the variance is $n^2/4$ and four times it is $n^2$. The sensitivity $\delta\theta=1/\sqrt F=1/n$ is the Heisenberg limit, a factor $\sqrt n$ below the shot-noise $1/\sqrt n$ of $n$ independent photons. This is the lossless advantage we now expose to loss.

## Photon loss as a Lindblad channel: evolving each arm toward vacuum

Loss is each arm leaking photons into an unwatched environment. For a mode with annihilation operator $\hat a$, pure loss is the Lindblad equation $\dot{\hat\rho}=\kappa\,(\hat a\hat\rho\hat a^\dagger-\tfrac12\{\hat a^\dagger\hat a,\hat\rho\})$, under which the mean photon number decays as $\langle\hat n\rangle(t)=\langle\hat n\rangle(0)\,e^{-\kappa t}$. So a channel of transmission $\eta$, meaning each photon survives with probability $\eta$, is loss at rate $\kappa=1$ run for a time $t=-\log\eta$. In other words, we do not need a separate loss parameter: the transmission $\eta$ fixes the evolution time.

`QuantumEvolve` takes the generator, the jump operators with their rates, the initial state, and the time: with a zero Hamiltonian and unit rate on each arm's annihilation operator, evolving to $t=-\log\eta$ applies transmission $\eta$ to both modes. Build the annihilation operators on the two-mode space with the same explicit-dimension identity, and return the lossy density matrix:

```wl
lossyNOON[n_, eta_] := With[{d = n + 1},
    With[{
        a1 = QuantumTensorProduct[AnnihilationOperator[d, {1}], QuantumOperator[IdentityMatrix[d], {2}, d]],
        a2 = QuantumTensorProduct[QuantumOperator[IdentityMatrix[d], {1}, d], AnnihilationOperator[d, {2}]]
        },
        With[{evolution = QuantumEvolve[0. (a1["Dagger"] @ a1), {a1, a2} -> {1., 1.},
                noonState[n], {t, 0, -Log[eta]}]},
            (evolution /. t -> -Log[eta])["DensityMatrix"]]]]
```

Loss turns the pure NOON state into a mixture. Check that the order-two state at $\eta=0.9$ is no longer pure: its purity $\mathrm{Tr}(\hat\rho^2)$ drops below one:

```wl
With[{r = Normal @ lossyNOON[2, 0.9]}, Re @ Tr[r . r]]
```

As expected, the purity is below one: a fraction of the population has fallen out of the two-photon NOON manifold into states with a missing photon, which carry which-arm information and cannot interfere.

## The Fisher information under loss: $F=n^2\eta^n$

For a mixed state the Fisher information is no longer four times a variance. The general form sums over the eigenpairs $(\lambda_i,\lvert i\rangle)$ of $\hat\rho$,
$$
F=2\sum_{i,j}\frac{(\lambda_i-\lambda_j)^2}{\lambda_i+\lambda_j}\,\lvert\langle i\rvert\hat G\lvert j\rangle\rvert^2,
$$
with $\hat G$ the phase generator and the sum restricted to pairs with $\lambda_i+\lambda_j>0$. In words, only coherences between eigenstates of different population contribute, weighted by how distinguishable the two populations are. Implement it, symmetrizing $\hat\rho$ to exactly Hermitian first so the eigensolver returns real eigenvalues:

```wl
qfi[rho_, g_] := Module[{l, v, m},
    {l, v} = Eigensystem[(rho + ConjugateTranspose[rho])/2];
    l = Chop[Re[l]];
    m = Conjugate[v] . g . Transpose[v];
    2 Total[Flatten @ Table[
        If[l[[i]] + l[[j]] > 1*^-12, (l[[i]] - l[[j]])^2/(l[[i]] + l[[j]]) Abs[m[[i, j]]]^2, 0],
        {i, Length[l]}, {j, Length[l]}]]]
```

Compare the computed Fisher information with the closed form $n^2\eta^n$ at $\eta=0.9$, for NOON orders up to ten:

```wl
Table[With[{d = n + 1},
    {n, qfi[Normal @ lossyNOON[n, 0.9], Normal[numberTwo[n]["Matrix"]]], N[n^2 0.9^n]}],
    {n, 2, 10, 2}]
```

As one can see, the two columns agree to six or seven digits: the Fisher information of the lossy NOON state is $F=n^2\eta^n$. The mechanism is clean. The $n^2$ sensitivity lives entirely on the coherence between $\lvert n,0\rangle$ and $\lvert0,n\rangle$, the only pair whose populations differ by the full $n$. That coherence survives only when all $n$ photons are transmitted, an event of probability $\eta^n$. Every term in which a photon is lost maps $\lvert n,0\rangle\to\lvert n-1,0\rangle$, a state orthogonal to anything the other branch collapses to, so it is which-arm information and carries no phase sensitivity at all. The surviving coherence keeps its full $n^2$, scaled by its weight $\eta^n$.

Check the lossless limit as a consistency test: at $\eta=1$ the evolution time is zero and the Fisher information must return to the pure-state $n^2$:

```wl
Table[qfi[Normal @ lossyNOON[n, 1.], Normal[numberTwo[n]["Matrix"]]], {n, 2, 4}]
```

As expected, $F=n^2$ is recovered: with nothing lost, $\eta^n=1$ and the mixed-state formula reduces to the pure-state result.

## When the advantage collapses: loss turns $n^2$ into $n^2\eta^n$

Having checked $F=n^2\eta^n$ against QF, we can read the sensitivity straight off the closed form. The lossy NOON error is $\delta\theta=1/\sqrt{n^2\eta^n}=1/(n\,\eta^{n/2})$, to be set against the lossless Heisenberg $1/n$ and the shot-noise $1/\sqrt n$. Define the three sensitivities:

```wl
sensitivities[eta_] := {
    "NOON, lossless" -> Function[n, 1/n],
    "shot noise" -> Function[n, 1/Sqrt[n]],
    "NOON, loss " <> ToString[eta] -> Function[n, 1/Sqrt[n^2 eta^n]]}
```

For fixed $\eta<1$ the factor $\eta^{n/2}$ shrinks faster than $n$ grows, so the lossy NOON error stops falling, turns around, and eventually climbs past the shot-noise line. Solve for the photon number where the lossy NOON Fisher information drops back to shot noise, $n^2\eta^n=n$:

```wl
n /. FindRoot[n 0.9^n == 1, {n, 30}]
```

At $\eta=0.9$ the crossover sits near this photon number: below it the NOON state still beats independent photons, above it the surviving coherence is too rare to pay for the delocalization, and more photons make the estimate worse. Plot the three sensitivities on a log scale to see the lossy curve bottom out and cross:

```wl
ListLogPlot[
    Table[sensitivities[0.9][[k, 2]] /@ Range[40], {k, 3}],
    Joined -> True, Mesh -> All, Frame -> True, GridLines -> Automatic,
    PlotLegends -> (First /@ sensitivities[0.9]),
    FrameLabel -> {"photon number n", "\[Delta]\[Theta] = 1/\!\(\*SqrtBox[\(F\)]\)"},
    PlotLabel -> "Phase sensitivity: loss turns the Heisenberg curve around", ImageSize -> 560]
```

The lossless NOON line falls as $1/n$ forever, always below shot noise; the lossy line tracks it at small $n$, reaches a minimum near the photon number where $n^2\eta^n$ peaks, and then rises through the shot-noise line at the crossover just computed. The Heisenberg scaling is not merely degraded by loss, it is reversed: past the optimum, adding photons to a lossy NOON state costs precision rather than buying it.

This is the qualitative reason NOON states are not how large-$n$ optical metrology is actually done. The loss-limited bound of [Demkowicz-Dobrzanski et al.](https://arxiv.org/abs/1201.3940) says the best achievable error under transmission $\eta$ scales as $\delta\theta\gtrsim\sqrt{(1-\eta)/(\eta N)}$, a shot-noise-like $1/\sqrt N$ with a loss-dependent prefactor; coherent-plus-squeezed light approaches it, while the NOON state, with its $\eta^n$ cliff, falls short for all but the smallest $n$.

## The answer

The quantum Fisher information of a NOON state of order $n$ is $F=n^2$ in the ideal case and $F=n^2\eta^n$ under per-photon transmission $\eta$, a result the mixed-state Fisher-information formula confirms against the closed form to six digits for $n$ up to ten. The $n^2$ factor is the sensitivity of the single coherence between the two extremal Fock states; the $\eta^n$ factor is the probability that this coherence survives, since every lost photon is a which-arm witness that kills it. For any $\eta<1$ the product $n^2\eta^n$ rises, peaks, and falls, so the phase error $1/\sqrt{n^2\eta^n}$ has a minimum at a finite $n$ and climbs back through the shot-noise line beyond it. The Heisenberg advantage is real only while the whole state stays intact, and photon loss makes that a small-$n$ luxury.

<!-- Cells smoke-tested in Wolfram Language with the checked-out QuantumFramework paclet (2.1.0), 2026-10-12: pure QFI = n^2 exactly; lossy QFI matches n^2 eta^n to 6-7 digits for n up to 10 at eta = 0.9; eta = 1 recovers n^2. -->
