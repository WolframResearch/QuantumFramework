---
Template: Default
Title: "Backbone: A Cat State, and Why Its Wigner Function Goes Negative"
Author: Mads Bahrami
---

# Backbone: A Cat State, and Why Its Wigner Function Goes Negative

This page is the statement layer of the essay *Cat-State-Wigner.md*: a graph of assumptions, choices, and derived statements, each in words and formulas only. Every node names its parents, so the dependency graph is recoverable from the page. Conventions. This essay is kinematic, not dynamical: there is no equation of motion and no time evolution, so the central mathematical object the analysis consumes is the Weyl-Wigner map itself, posited as a contract E with a scope statement, construction (the derivation of the map from the Weyl correspondence) declared out of scope. Rule 1 permits this because no weakening on the page reaches beneath the map: the two weakenings act at the Hilbert-space level (the cutoff) and the state level (the amplitude), never on the correspondence. Every assumption must be consumed by at least one statement that is not its own weakening. Every statement's tail separates provenance (derived in place, cited, or conjectured) from verification (exact decision, analytic cross-check, numerical check, or none), naming the clause scope. Node labels: R root, C choice, E the posited contract, S derived statement, H headline, N recorded non-assumption, prime a weakened variant. Nothing here executes; the verifications all live in the essay, whose cells are the only computation, and are numerical throughout, as befits a "numeric explicit examples" showcase.

## The setup, and the choices

**R1 (assumption).** A single bosonic mode, its infinite Fock ladder $|n\rangle$ truncated at a finite cutoff $N$ (here $N=40$). In words: the state lives on a finite-dimensional space large enough that the physical state's weight above the cutoff is negligible. R1 fixes the arena; it is a numerical convenience whose irrelevance is recorded below (N2), not a physical claim.

<!-- box: single mode, Fock ladder | cut at N = 40 -->

**R2 (assumption).** The state is the normalized superposition of two coherent states, $|\psi\rangle = \mathcal N(|\alpha\rangle + e^{i\phi}|{-\alpha}\rangle)$, with $\alpha$ real and large enough that the lobes are well separated: at $\alpha=2.8$ the mean occupation is $|\alpha|^2\approx 8$ and the overlap $|\langle\alpha|{-\alpha}\rangle| = e^{-2\alpha^2}$ is negligible. Real $\alpha$ places both lobes symmetrically on the position axis, at positions $\pm\sqrt2\,\alpha$.

<!-- box: two coherent lobes, α real; | overlap e^{−2α²} ≪ 1 -->

**C1 (choice).** The relative phase $\phi$ selects the parity sector: $\phi=0$ is the even cat (support on even $n$), $\phi=\pi$ the odd cat (odd $n$), $\phi=\pi/2$ the Yurke-Stoler state (neither a parity eigenstate). The essay fixes $\phi\in\{0,\pi\}$; any other phase breaks the parity eigenstate property that S4 rests on, so this choice and the parity statements travel together. *Derived in place (the $(-1)^n$ phase of $|{-\alpha}\rangle$ makes $\phi=0,\pi$ the two parity eigenstates). Verified by numerical check, whole node: the even and odd number combs computed directly; no shared code.*

<!-- box: φ = 0 even, φ = π odd; | π/2 Yurke-Stoler -->

**C2 (choice).** The phase-space normalization. `WignerRepresentation` works in the $(x,p)$ plane and normalizes the vacuum to $W(0,0)=1/\pi$ with $|W|\le 1/\pi$; `WignerFunction` works in the complex amplitude $\alpha=(x+ip)/\sqrt2$ and normalizes the vacuum to $2/\pi$. The two differ only by the Jacobian between $d^2\alpha$ and $dx\,dp$; the essay reads the origin value in both and reports the same physics. *Derived in place (the factor of two is the change of area element). Verified by numerical check, whole node: the two routes evaluated at the origin, ratio exactly two; no shared code.*

<!-- box: convention: W(0,0) = 1/π (xp) | vs 2/π (complex α) -->

## The map, posited as a contract

**E (contract; cited).** The Wigner function is the phase-space representation of the density matrix $\rho$: a real, normalized quasiprobability whose marginals are the true $x$ and $p$ distributions but which may take negative values. It is bounded, $|W|\le 1/\pi$ in the $(x,p)$ convention, and its value at any point is a displaced photon-number parity expectation,
$$W(x,p) = \tfrac{2}{\pi}\,\big\langle \hat D^\dagger(x,p)\,\hat\Pi\,\hat D(x,p)\big\rangle, \qquad \hat\Pi = (-1)^{\hat a^\dagger\hat a},$$
so at the origin $W(0,0)=\tfrac{2}{\pi}\langle\hat\Pi\rangle$ (and $\tfrac1\pi\langle\hat\Pi\rangle$ in the $(x,p)$ normalization). Scope: a single mode, the standard Weyl-Wigner correspondence, read on the truncated space of R1 to the tail weight recorded in N2. *Cited: Royer, Phys. Rev. A 15, 449 (1977) for the displaced-parity form; Cahill and Glauber, Phys. Rev. 177, 1882 (1969); Haroche and Raimond, Exploring the Quantum (Oxford 2006). Verification: none; the map is trusted to its source, and the essay computes with exactly this object.*

<!-- box: Wigner map: real, |W| ≤ 1/π, | W(0,0) = (2/π)⟨Π⟩ -->

## Consequences

**S1 (from R2, C1).** Expanding $|{\pm\alpha}\rangle$ on Fock states, the $|{-\alpha}\rangle$ amplitudes carry $(-1)^n$, so the even cat ($\phi=0$) keeps only even $n$ and doubles them while the odd cat ($\phi=\pi$) keeps only odd $n$. The two number-probability combs are strictly interleaved: the even-or-odd alternation is the interference of the two coherent branches written in the number basis. *Cited (coherent-state Fock expansion, standard); the parity selection derived in place. Verified by numerical check, whole node: the first ten number probabilities of each cat, even-only and odd-only after chopping machine residue; no shared code.*

<!-- box: number comb: even cat on even n, | odd cat on odd n -->

**S2 (from S1).** Each lobe $|{\pm\alpha}\rangle$ is an eigenstate of $\hat a$ with eigenvalue $\pm\alpha$, so each is an eigenstate of $\hat a^2$ with the common eigenvalue $\alpha^2$, and so is any superposition of them: the cat is an eigenstate of $\hat a^2$ with eigenvalue $\alpha^2$. Equivalently, removing two photons preserves the parity sector of S1. *Derived in place. Verified by numerical check, whole node: $\hat a^2|\psi\rangle$ normalized is similarity one with $|\psi\rangle$, and $\langle\psi|\hat a^2|\psi\rangle=\alpha^2$; a decidable claim computed numerically here, upgradeable to exact decision.*

<!-- box: â² eigenstate, | eigenvalue α² -->

**S3 (from E, R2).** In phase space the state is two Gaussian lobes at $x=\pm\sqrt2\,\alpha$, the images of $|\alpha\rangle$ and $|{-\alpha}\rangle$, with a ridge of fringes between them that oscillate along $p$ and dip below zero. Two separated Gaussians alone give a non-negative function; the negativity is created entirely by the cross term $|\alpha\rangle\langle{-\alpha}| + \text{h.c.}$ and lives between the lobes, not on them. *Cited (cat Wigner function, Haroche and Raimond); the non-negativity of two bare Gaussians derived in place. Verified by numerical check, whole node: the grid minimum of the even cat's Wigner function is negative; no shared code.*

<!-- box: two lobes + fringe ridge; | grid minimum < 0 -->

**S4 (from E, S1, C2).** The origin sits halfway between the lobes, and its value is fixed, not free: by E it is $\tfrac{2}{\pi}\langle\hat\Pi\rangle$, and by S1 the even cat is a $+1$ parity eigenstate and the odd cat a $-1$ one. So the even cat is pinned to the top of the allowed range and the odd cat to the bottom; in C2's two conventions these extremes read $\pm1/\pi$ and $\pm2/\pi$, and the deepest and the highest fringe values both sit at the origin. *Derived in place from E and S1's parity eigenvalue. Verified by numerical check, whole node: the origin value in both conventions ($1/\pi$ and $2/\pi$) and the parity trace $\pm1$ for the two cats; no shared code.*

<!-- box: origin pinned to parity: | even +1/π, odd −1/π -->

## The headline

**H (from S3, S4, C1).** The negativity of the cat's Wigner function is the phase-space fingerprint of the superposition. It is made entirely by the interference cross term between the two coherent lobes, which are classical non-negative Gaussians on their own; it lives in the region between them; and its value at the origin is locked to the photon-number parity, $+1/\pi$ for the even cat and $-1/\pi$ for the odd one, the extreme values the map allows. The same even-or-odd parity that forbids half the Fock comb (S1, S2) is what pins the origin sign (S4): one parity eigenvalue read two ways, in the number basis and in phase space. A coherent state of the same amplitude is a single non-negative Gaussian with none of it. Audit: H rests on the contract E together with R2, C1, C2, and the derived S3 and S4; E is posited with the scope stated in its node (a single mode, the standard Weyl-Wigner correspondence, read on a truncated space whose tail weight is negligible), so the analysis holds for any realization of that map. No conjecture anywhere in the closure. *Derived in place from S3 and S4. Verified by numerical check, whole node, through the grid minimum and origin values of S3 and S4; no shared code.*

<!-- box: negativity = interference; | origin locked to parity -->

## When an assumption is weakened

**S3′ (from R2 superposition weakened: one lobe).** Drop the superposition and keep a single coherent state $|\alpha\rangle$ of the same amplitude: its Wigner function is one non-negative Gaussian, grid minimum zero to machine precision, no fringes. This is the corner model that proves the negativity of S3 and H is the cross term and nothing else, not the lobes. *Derived in place (a coherent state is Gaussian and non-negative). Verified by numerical check, whole node: the coherent-state grid minimum is zero at $\alpha=2.8$; no shared code.*

<!-- box: one lobe: non-negative | Gaussian, min 0 -->

**N1 (a recorded non-assumption on the separation in R2).** The headline's *geometry* needs well-separated lobes, but its origin value does not. For any $\alpha\neq0$ the even and odd cats remain exact parity eigenstates, so $W(0,0)=\pm1/\pi$ exactly, independent of $\alpha$; shrinking $\alpha$ only merges the fringe ridge into the lobes until the "two lobes with interference between them" picture dissolves. Separation is load-bearing for the phase-space image and for reading the negativity as a *between-the-lobes* fringe, never for the origin-parity lock. *Derived in place (parity is an exact eigenvalue at any amplitude). Verified by numerical check: the origin value stays $\pm1/\pi$ as the amplitude is reduced; the single-point manifest lives beside the check if run.*

<!-- box: separation sets the geometry, | not the origin value ±1/π -->

**N2 (a recorded non-assumption on the cutoff in R1).** The cutoff $N$ is a numerical convenience, not a physical assumption: the cat's weight above level 30 is vanishingly small, so the number comb, the fringes, and the origin value are unchanged by raising $N$. The boundary of the deletion: $N$ must stay comfortably above the mean occupation $|\alpha|^2$; at a cutoff near the occupation the truncation would distort the lobes. *Derived in place (negligible tail population). Verified by numerical check, whole node: the total probability above level 30 is vanishingly small at $\alpha=2.8$.*

<!-- box: cutoff irrelevant above n̄; | tail weight ≈ 0 -->

## The abstract, read off the graph

A single bosonic mode carries the normalized superposition of two well-separated coherent states, even or odd by the relative phase. In the number basis this superposition is a strict even-or-odd comb, which makes the cat an eigenstate of $\hat a^2$ with eigenvalue $\alpha^2$. The Wigner map represents the state in phase space as two Gaussian lobes with a ridge of interference fringes between them, and those fringes go negative: the negativity is the cross term alone, since two bare Gaussians are non-negative and a single coherent state of the same amplitude shows none of it. The value at the origin is not free but fixed by photon-number parity, the same parity that empties half the number comb, pinning the even cat to $+1/\pi$ and the odd cat to $-1/\pi$, the extremes the map allows. Separation sets the geometry but not that origin value, which is exact at any amplitude; the Fock cutoff is irrelevant once it clears the occupation. The dip below zero is the phase-space fingerprint of the superposition, the signature of a nonclassicality no classical field can supply.
