---
Template: Default
Title: "Backbone: Rapid Purification by Choosing What to Measure"
Author: Mads Bahrami
---

# Backbone: Rapid Purification by Choosing What to Measure

This page is the statement layer of the essay: a graph of assumptions, choices, and derived statements, each given in words and formulas only. Every node names its parents, so the dependency graph is recoverable from the page, and the layers order the reading from the setup to the headline and then to what holds when an assumption is weakened. Three conventions hold throughout. The page has two layers: a construction layer, where the physical roots are translated into the equation S0 by cited standard theory, and an analysis layer, where every consequence consumes S0 as a contract and therefore survives any construction that yields the same equation; the equation is never a root here, because the weakenings below (lose part of the record, open a channel) reach beneath it. Every assumption must be consumed by at least one statement that is not its own weakening; a root failing that test is not an assumption of the theory but only a direction along which the model can be varied. And every statement's tail separates provenance (derived in place, cited, or conjectured) from verification (exact decision, analytic cross-check, numerical check, or none), naming the clause scope and what the checking route shares with the finding route. Node labels: R is a root, a setup assumption; C is a choice, a meaning we fixed rather than a fact; S is a derived statement, with S0 reserved for the equation of motion (E would name it if posited as a contract without its construction, which this page does not do); H is the headline, the terminal statement; N is a recorded non-assumption; a prime marks the same statement under a weakened root. Nothing here executes; the verifications live in the companion file *two-examples-from-minimal-assumptions.md*, downstream of this page, whose computations share no code with the derivations here.

## The setup, and the choice

**R1 (assumption).** One observable of the qubit, $\hat M = \vec n\cdot\vec\sigma$, is watched continuously and weakly, with diffusive readout: the detector's output is a noisy current with Gaussian noise, not a train of discrete clicks, and the watching has a single strength $k$. This fixes which unravelling the story lives in.

**R2 (assumption).** The record is complete: everything the measurement emits is collected, detection efficiency $\eta = 1$.

**R3 (assumption).** The watching is the only coupling: no other channel touches the qubit. This is a closure statement, a different kind of claim from R1: R1 asserts that a particular coupling exists and has a particular character, R3 asserts that the list of couplings ends there, and closure never follows from existence. The model of S4, which keeps R1 and violates R3, separates the two.

**C1 (choice).** The working variable: the impurity $\ell = \tfrac12\left(1 - \left|\vec a\right|^2\right) = 1 - \mathrm{Tr}\,\rho^2$, one half for the maximally mixed state, zero for a pure one. A definition, so nothing here is true or false; what can be wrong is only whether ranking by it is fair, which is C2's business.

**C2 (choice; from C1).** The score: protocols are ranked by $\ell$. For a qubit the von Neumann entropy is a strictly increasing function of $\ell$, so every ranking and every passage-time statement below is score-independent: the choice is provably harmless. *Derived in place (monotonicity of entropy in impurity for a qubit). Verified by numerical check, whole node: monotonicity on a grid of states.*

No "no drive" assumption is listed; node N records why.

## The equation, derived from the setup

**S0 (from R1, R2, R3; cited).** Quantum trajectory theory turns the setup into the conditioned-state equation
$$d\rho = k\,\mathcal D[\hat M]\rho\,dt + \sqrt{k}\left(\hat M\rho + \rho\hat M - 2\langle\hat M\rangle\rho\right)dW, \qquad dy = \langle\hat M\rangle\,dt + \frac{dW}{2\sqrt{k}},$$
the record normalization being a declared convention. Each root is visible in the equation's shape: R1 fixes the diffusive form and the single strength; R2 is the noise coefficient being $\sqrt{k}$ rather than $\sqrt{\eta k}$, equivalently the fact that this equation keeps pure states exactly pure: nothing is being discarded; R3 is the absence of any further generator term. *Cited: Wiseman and Milburn, Quantum Measurement and Control (Cambridge 2010), and Jacobs and Steck, Contemporary Physics 47, 279 (2006). Verification: none; the translation is trusted to its sources, and the companion computes with exactly this equation.*

## Consequences

**S1 (from S0, C1).** The exact increment of C1's impurity, with $\theta$ the angle between the Bloch vector and the measured axis:
$$d\ell = -4k\,\ell\left(\sin^2\theta + 2\ell\cos^2\theta\right)dt \;-\; 4\sqrt{k}\,\ell\sqrt{1 - 2\ell}\,\cos\theta\,dW.$$
In words: the drift is never positive, so every kept run purifies; how fast is pure geometry, and the randomness rides entirely on $\cos\theta$. Averaged over records the state stays maximally mixed, so "purifies" already means three inequivalent things: every run, the mean, or the time to a target. The headline is their disagreement. *Cited: Jacobs 2003, quant-ph/0301056; rederivable by Ito calculus on S0. Verified by numerical check, through its consequences only: the increment itself was never tested directly, but S2's law and H's reversal, which follow from it, were; no shared code.*

**S2 (from S1).** Two protocols, one dial. *Fixed axis*: the sharpening state drifts into alignment, $\theta \to 0$, where the drift is weakest (of order $\ell^2$) and the noise largest: slow on average, lucky on occasion. *Steered*: re-aim the axis perpendicular to the state at every instant, $\theta \equiv \pi/2$: the noise coefficient vanishes and
$$\ell(t) = \ell_0\,e^{-4kt} \ \ \text{on every record.}$$
Steering is a choice of what to measure next; no Hamiltonian is ever applied. *Derived in place from S1. Verified by numerical check, whole node: thirty steered runs clustering on the law, with a step refinement contracting both the spread and the gap; no shared code.*

## The headline

**H (from S2, C2).** "Faster" is two claims, and they disagree. At any fixed deadline the steered mean impurity is the lower one; the mean time to reach a fixed target is shorter for the fixed axis, whose lucky records cross early, a luck that determinism forbids. Mechanism: determinism buys the best average at a deadline by selling off the fast tail of the passage-time distribution. Any sentence of the form "this protocol purifies faster," with the score unnamed, has not yet said anything. Audit: H rests on the contract S0 together with C2 and, through it, C1; S0 is realized from R1, R2, R3 by the cited translation, so the analysis survives any construction that yields the same equation. No conjecture anywhere in the closure. *Cited: Wiseman and Ralph 2006, quant-ph/0603062. Verified by numerical check, whole node: both scores computed on the same two ensembles, the reversal observed directly; no shared code.*

## When an assumption is weakened

**N (a recorded non-assumption; from S2).** A known rotation between measurements is absorbed by the next re-aim, so the steered law survives a drive unchanged; for the fixed axis a drive even helps, pulling the state off alignment, a partial substitute for steering. Stillness was propping up the baseline's slowness, not the headline. *Derived in place from S2 (a known rotation is absorbed by the next re-aim). Verified by numerical check at a single drive strength, $\Omega = 1.5$: the driven and undriven steered ensembles coincide, and the driven fixed axis falls between the two; other strengths untested.*

**S3 (from R2 weakened; S0, C1).** Collect only a fraction $\eta < 1$ and trajectory theory changes S0 in exactly one place: the noise coefficient becomes $\sqrt{\eta k}$ while the dissipator keeps its full $k$, because the measurement disturbs at full strength and informs only in proportion to what is kept; the unrecorded fraction acts as dephasing about the measured axis. The verdict then reverses at long times. The steered protocol holds the state perpendicular to that axis, precisely where the unrecorded dephasing bites hardest; at the held angle the drift balance
$$2k\left(1 - 2\ell\right) - 2\eta k = 0 \quad\Longrightarrow\quad \ell_\ast = \frac{1 - \eta}{2}$$
is a floor steering cannot cross. The fixed axis has no floor: its destinations are the measured eigenstates, invisible to that dephasing, so its runs still end pure, later, and its mean passes below the steered floor and keeps going. Optimality bought against an ideal assumption pays rent on it forever. *Derived in place (the drift balance). Verified by numerical check at a single efficiency, $\eta = 0.7$: the floor lands on the prediction and the fixed axis crosses below it there; other efficiencies untested.*

**S4 (from R3 weakened; S0).** Open an additional channel the record does not see, and purity is capped exactly when that channel's own stationary state is mixed: an unwatched warm pair imposes a lasting impurity plateau that no choice of measured axis removes, while an unwatched decay does not, since its attractor is the pure ground state and every run still ends pure by that other route. Weakening an assumption is itself a physical act: it requires reading the new channel's fixed point first. *Derived in place (reading each channel's attractor). Verified by numerical check at one warm pair, $\gamma = 0.4$ and $\bar n = 1$, against one decay channel; the computation corrected this statement's first draft, which had claimed decay would cap purity.*

## The abstract, read off the graph

A qubit watched along one axis purifies on every kept record: the measurement gives the impurity a never-positive drift whose size depends only on the angle between the state and the measured axis, while the record-averaged state never purifies at all. Holding the axis fixed lets the state align with it, where learning is slowest and luckiest; re-aiming the axis perpendicular to the state at every instant removes the randomness and makes the impurity fall deterministically at $4k$ on every record. Which protocol is faster is then two different questions: steering wins any fixed deadline on average, the fixed axis wins the race to a target, because determinism forbids lucky early crossings. A leaky detector reverses the long-time verdict, capping steering at impurity $(1-\eta)/2$ while the fixed axis drives toward eigenstates immune to the lost record; an extra unwatched channel caps purity only if its own attractor is mixed; and a drive is a partial substitute for steering, never its defeat.
