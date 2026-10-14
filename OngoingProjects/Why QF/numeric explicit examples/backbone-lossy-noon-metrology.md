---
Template: Default
Title: "Backbone: How Fast Photon Loss Erases the Heisenberg Advantage of a NOON State"
Author: Mads Bahrami
---

# Backbone: How Fast Photon Loss Erases the Heisenberg Advantage of a NOON State

This page is the statement layer of the essay *Lossy-NOON-Metrology.md*: a graph of assumptions, choices, and derived statements, each given in words and formulas only. Every node names its parents, so the dependency graph is recoverable from the page, and the layers order the reading from the setup to the headline and then to what holds when an assumption is weakened. Three conventions hold throughout. The page has two layers: a construction layer, where the physical roots are translated into the loss equation S0 by cited standard theory, and an analysis layer, where every consequence consumes S0 as a contract and therefore survives any construction that yields the same channel; the equation is never a root here, because the weakening below (set the loss to zero) reaches beneath it. Every assumption must be consumed by at least one statement that is not its own weakening; a root failing that test is not an assumption of the theory but only a direction along which the model can be varied. And every statement's tail separates provenance (derived in place, cited, or conjectured) from verification (exact decision, analytic cross-check, numerical check, or none), naming the clause scope. Node labels: R root, C choice, S derived statement with S0 reserved for the equation of motion, H headline, N recorded non-assumption, prime a weakened variant. Nothing here executes; the verifications all live in the essay, whose cells are the only computation, and are numerical throughout, as befits a "numeric explicit examples" showcase.

## The setup, and the choice

**R1 (assumption).** Two bosonic modes are the two arms of an interferometer, carrying the NOON state of order $n$: all $n$ photons in the first arm or all $n$ in the second, in equal superposition, $(|n,0\rangle + |0,n\rangle)/\sqrt2$. This maximally delocalized two-mode superposition is the probe whose fate under loss the page tracks.

<!-- box: two arms, order-n NOON: all n | in one arm or the other -->

**R2 (assumption).** The phase $\theta$ to be estimated is imprinted on one arm, $e^{-i\theta\hat n_2}$, so the generator of the estimation is that arm's number operator $\hat n_2 = \hat a_2^\dagger\hat a_2$. The two branches wind at opposite ends, so the superposition winds $n$ times faster than a single photon.

<!-- box: phase on one arm; | generator n̂₂ -->

**R3 (assumption).** Each arm leaks photons into an unwatched environment: pure loss as the Lindblad channel with jump operator $\hat a$ at rate $\kappa$, zero temperature. Transmission $\eta$, the probability a photon survives, is loss at unit rate run for time $t=-\log\eta$, so $\eta$ alone fixes the evolution, no separate loss parameter needed.

<!-- box: per-arm pure loss; | transmission η = e^{−t} -->

**C1 (choice).** Precision is measured by the quantum Fisher information $F$: the quantum Cramer-Rao bound $\delta\theta \ge 1/\sqrt F$ for one shot is the best error over all measurements, not the variance of one particular estimator. Shot noise is $F=n$, the Heisenberg limit $F=n^2$. Choosing the state-intrinsic bound rather than a fixed readout is what lets the page speak of the advantage itself, not of one instrument's share of it. *Derived in place (the quantum Fisher information upper-bounds every measurement's classical Fisher information). Verified by numerical check, whole node: $F$ computed from the state and matched to the closed form; no shared code.*

<!-- box: figure of merit: quantum | Fisher information F, δθ ≥ 1/√F -->

## The equation, derived from the setup

**S0 (from R3; cited).** Standard open-system theory turns per-arm loss into the Lindblad channel
$$\dot{\hat\rho} = \kappa\left(\hat a\hat\rho\hat a^\dagger - \tfrac12\hat a^\dagger\hat a\,\hat\rho - \tfrac12\hat\rho\,\hat a^\dagger\hat a\right), \qquad \langle\hat n\rangle(t) = \langle\hat n\rangle(0)\,e^{-\kappa t},$$
run on each mode, with the mean photon number decaying exponentially; transmission $\eta$ is this channel evolved to $t=-\log\eta$. Root R3 is visible as the single zero-temperature dissipator on $\hat a$ with no upward counterpart. *Cited: pure photon loss as a Lindblad channel, Walls and Milburn, Quantum Optics; the amplitude-damping channel. Verification: none; the translation is trusted to its source, and the essay evolves with exactly this generator.*

<!-- box: ρ̇ = κ D[a]ρ per arm; | ⟨n⟩(t) = ⟨n⟩(0) e^{−κt} -->

## Consequences

**S1 (from R1, R2, C1).** With no loss the NOON state is pure, so C1's Fisher information is four times the variance of the generator R2, $F = 4(\langle\hat n_2^2\rangle - \langle\hat n_2\rangle^2)$. The two branches of R1 are eigenstates of $\hat n_2$ with eigenvalues $0$ and $n$, so the variance is $n^2/4$ and $F=n^2$ exactly: the Heisenberg limit $\delta\theta = 1/n$, a factor $\sqrt n$ below the shot-noise $1/\sqrt n$ of $n$ independent photons. *Derived in place (the two branches are $\hat n_2$-eigenstates with eigenvalues $0$ and $n$, variance $n^2/4$). Verified by numerical check ($F=n^2$ exactly for $n=1$ through $6$), upgradeable to exact decision; no shared code.*

<!-- box: lossless: F = 4 Var(n̂₂) = n² | (Heisenberg) -->

**S2 (from S0, R1).** The loss channel S0 turns the pure NOON state of R1 into a mixture. Every term in which a photon is lost maps $|n,0\rangle \to |n-1,0\rangle$, a state orthogonal to every descendant of the other branch $|0,n\rangle$, so a lost photon is a which-arm witness that cannot interfere. The coherence between $|n,0\rangle$ and $|0,n\rangle$, the only part carrying phase sensitivity, survives only when all $n$ photons are transmitted, an event of weight $\eta^n$. *Derived in place (loss sends $|n,0\rangle$ to a state orthogonal to the other branch; the intact coherence keeps weight $\eta^n$). Verified by numerical check, whole node: the purity of the lossy order-two state at $\eta=0.9$ drops below one; no shared code.*

<!-- box: loss → mixture; intact coherence | survives with weight η^n -->

**S3 (from S2, R2, C1).** For the mixture of S2, C1's Fisher information is no longer four times a variance but the eigenpair sum $F = 2\sum_{i,j}(\lambda_i-\lambda_j)^2/(\lambda_i+\lambda_j)\,|\langle i|\hat G|j\rangle|^2$ over the generator R2. Only the surviving coherence between $|n,0\rangle$ and $|0,n\rangle$ carries the full population gap $n$, so it alone contributes the $n^2$, scaled by its weight $\eta^n$:
$$F = n^2\eta^n.$$
*Derived in place (in the mixed-state Fisher sum only the surviving coherence spans the full gap $n$, contributing $n^2$ at weight $\eta^n$). Cited for the closed form in the lossy-metrology literature: Demkowicz-Dobrzanski et al., arXiv:1201.3940. Verified by analytic cross-check and numerical check: the mixed-state formula agrees with $n^2\eta^n$ to six or seven digits for $n$ up to ten at $\eta=0.9$; single transmission, others untested; no shared code.*

<!-- box: mixed-state QFI: F = n²η^n -->

## The headline

**H (from S3, S1).** Photon loss does not merely degrade the Heisenberg advantage, it reverses it. The Fisher information is $F=n^2$ ideally (S1) and $F=n^2\eta^n$ under per-photon transmission $\eta$ (S3): the $n^2$ sensitivity is the single coherence between the two extremal Fock states, and the $\eta^n$ is the probability that this coherence survives, since every lost photon is a which-arm witness that kills it. So the phase error $\delta\theta = 1/\sqrt{n^2\eta^n} = 1/(n\,\eta^{n/2})$ falls, bottoms out at a finite $n$, and climbs back through the shot-noise line where $n^2\eta^n = n$: beyond the optimum, adding photons costs precision. The Heisenberg advantage is real only while the whole state stays intact, a small-$n$ luxury under any loss. Audit: H rests on the loss contract S0 together with R1, R2, C1 and the derived S3, S1; S0 is realized from R3 by the cited translation, so the analysis survives any construction that yields the same loss channel. No conjecture anywhere in the closure. *Derived in place from S3 and S1 ($\delta\theta = 1/(n\eta^{n/2})$ has an interior minimum; $n^2\eta^n=n$ at the shot-noise crossover). Verified by numerical check, whole node, through S3's closed form and S1's lossless value, the crossover located by FindRoot; no shared code.*

<!-- box: F = n²η^n: advantage peaks | then falls past shot noise -->

## When an assumption is weakened

**S3′ (from R3 weakened; S3).** Set the loss to zero, $\eta=1$: the evolution time $t=-\log\eta$ is zero, the mixture of S2 collapses back to the pure NOON state, and S3's formula returns $F=n^2\eta^n = n^2$, S1's lossless value. This is the consistency corner that proves the $\eta^n$ factor is loss and nothing else: with nothing lost the surviving coherence has full weight and the Heisenberg limit is recovered exactly. *Derived in place ($\eta=1$ sends the mixed-state formula to the pure-state $n^2$). Verified by numerical check: $F=n^2$ recovered at $\eta=1$ for $n=2$ through $4$.*

<!-- box: η = 1: F = n², | Heisenberg recovered -->

**N (a recorded non-assumption on the truncation in R1).** No cutoff-error assumption is needed. Loss only lowers photon number and no arm ever holds more than $n$, so a Fock cutoff of $d=n+1$ levels per mode holds every reachable state exactly: the truncation is exact here, not an approximation to be converged. The boundary: this is special to pure loss from a bounded-number probe; a gain or drive channel would populate levels above $n$ and the cutoff would again be an approximation. *Derived in place (loss lowers number and no arm exceeds $n$, so $d=n+1$ is closed under the dynamics). Verified by numerical check, whole node: the truncation carries the exact state rather than a converged approximation.*

<!-- box: cutoff d = n+1 is exact, | not a convergence parameter -->

## The abstract, read off the graph

Two bosonic arms carry a NOON state of order $n$, all $n$ photons in one arm or the other, and a phase is imprinted on one arm so the superposition winds $n$ times faster than a single photon. With precision measured by the quantum Fisher information, the state-intrinsic best over all measurements, the lossless NOON state reaches $F=n^2$, the Heisenberg limit, a factor $\sqrt n$ below shot noise. Pure photon loss on each arm, standard Lindblad amplitude damping with transmission $\eta$ fixing the evolution time, turns that pure state into a mixture: every lost photon is a which-arm witness, so the phase-carrying coherence between the two extremal Fock states survives only with weight $\eta^n$, and the mixed-state Fisher information is $F=n^2\eta^n$. The phase error $1/(n\eta^{n/2})$ therefore falls, reaches a minimum at a finite photon number, and climbs back through the shot-noise line, so loss does not merely degrade the Heisenberg advantage but reverses it beyond the optimum. Setting the loss to zero recovers $F=n^2$ exactly, proving the $\eta^n$ factor is the loss alone; and because loss only lowers the photon number, a cutoff of $n+1$ levels per mode is exact rather than an approximation.
