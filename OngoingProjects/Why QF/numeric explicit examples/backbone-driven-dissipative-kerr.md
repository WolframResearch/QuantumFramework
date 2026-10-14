---
Template: Default
Title: "Backbone: A Nonclassical Steady State from Balanced Drive and Loss"
Author: Mads Bahrami
---

# Backbone: A Nonclassical Steady State from Balanced Drive and Loss

This page is the statement layer of the essay *Driven-Dissipative-Kerr.md*: a graph of assumptions, choices, and derived statements, each given in words and formulas only. Every node names its parents, so the dependency graph is recoverable from the page, and the layers order the reading from the setup to the headline and then to what holds when an assumption is weakened. Three conventions hold throughout. The page has two layers: a construction layer, where the physical roots are translated into the equation S0 by cited standard theory, and an analysis layer, where every consequence consumes S0 as a contract and therefore survives any construction that yields the same equation; the equation is never a root here, because the weakening below (turn off the Kerr term) reaches beneath it. Every assumption must be consumed by at least one statement that is not its own weakening; a root failing that test is not an assumption of the theory but only a direction along which the model can be varied. And every statement's tail separates provenance (derived in place, cited, or conjectured) from verification (exact decision, analytic cross-check, numerical check, or none), naming the clause scope. Node labels: R root, C choice, S derived statement with S0 reserved for the equation of motion, H headline, N recorded non-assumption, prime a weakened variant. Nothing here executes; the verifications all live in the essay, whose cells are the only computation, and are numerical throughout, as befits a "numeric explicit examples" showcase.

## The setup, and the choice

**R1 (assumption).** A coherent laser drives the mode, adding the term $F(\hat a + \hat a^\dagger)$ to the Hamiltonian in the frame rotating at the drive frequency, so the drive is time independent. In words: photons are pumped in at a fixed amplitude $F$, the only source term in the dynamics.

<!-- box: coherent drive F(a + a†) -->

**R2 (assumption).** The mode is a Kerr medium: the energy to add the next photon grows with the number already present, $\tfrac{K}{2}\,\hat a^\dagger\hat a^\dagger\hat a\hat a = \tfrac{K}{2}\,\hat n(\hat n-1)$, with a detuning term $-\Delta\,\hat a^\dagger\hat a$ that only sets the working point. The intensity-dependent resonance is the one nonlinearity in the model, and the whole nonclassical story rests on it.

<!-- box: Kerr: energy per photon | grows with photon number -->

**R3 (assumption).** The mode leaks through an imperfect mirror into a zero-temperature, memoryless bath: a single jump operator $\hat a$ at rate $\kappa$, absorbing photons and never returning them. This is the only dissipation.

<!-- box: zero-T loss, jump a at rate κ -->

**R4 (assumption).** The infinite Fock ladder is truncated at a finite cutoff $d$ (here $d=24$), large enough to hold the steady state's weight to the tail recorded in N. This is a numerical convenience that fixes the arena, not a physical claim.

<!-- box: Fock cutoff d (numerical) -->

**C1 (choice).** Nonclassicality is read from the photon-number statistics: the Mandel $Q = (\langle\Delta\hat n^2\rangle - \langle\hat n\rangle)/\langle\hat n\rangle$, negative exactly when the distribution is narrower than the Poisson of equal mean, together with the second-order coherence $g^{(2)}(0) = \langle\hat n(\hat n-1)\rangle/\langle\hat n\rangle^2$. The reasonable rival witness is phase-space negativity of the Wigner function; the two are inequivalent, and S3 shows they disagree here, so the choice of witness is load-bearing. *Derived in place (both $Q$ and $g^{(2)}(0)$ measure the departure of the number variance from Poisson). Verified by numerical check, whole node: both computed from the steady-state populations; no shared code.*

<!-- box: witness: Mandel Q, g^{(2)}(0) | (number statistics) -->

## The equation, derived from the setup

**S0 (from R1, R2, R3; cited).** Standard open-system theory turns the setup into the Lindblad master equation
$$\dot{\hat\rho} = -i[\hat H,\hat\rho] + \kappa\left(\hat a\hat\rho\hat a^\dagger - \tfrac12\hat a^\dagger\hat a\,\hat\rho - \tfrac12\hat\rho\,\hat a^\dagger\hat a\right), \qquad \hat H = -\Delta\,\hat a^\dagger\hat a + \tfrac{K}{2}\,\hat a^\dagger\hat a^\dagger\hat a\hat a + F\left(\hat a + \hat a^\dagger\right),$$
with each root visible in its shape: the drive R1 is the source term $F(\hat a + \hat a^\dagger)$, the Kerr nonlinearity R2 is the quartic term, and the zero-temperature memoryless loss R3 is the single dissipator at one constant rate $\kappa$ with no upward counterpart. *Cited: Drummond and Walls, J. Phys. A 13, 725 (1980), the driven-dissipative Kerr (anharmonic-oscillator) model; the Lindblad translation standard, Breuer and Petruccione, The Theory of Open Quantum Systems (Oxford 2002), chapter 3. Verification: none; the translation is trusted to its sources, and the essay computes with exactly this equation.*

<!-- box: dρ/dt = −i[H,ρ] + κ D[a]ρ -->

## Consequences

**S1 (from S0, R4).** At steady state $\dot{\hat\rho}=0$, so writing S0 as $\dot{\hat\rho} = \mathcal L[\hat\rho]$ the steady state is whatever the Liouvillian $\mathcal L$ annihilates: its null space, found as a single linear-algebra step on the truncated space of R4 rather than by waiting out an integration. Two facts are read off $\mathcal L$'s spectrum. The null space is one dimensional at these parameters, so the steady state is unique; uniqueness is a computed property, not announced by the model, and a Kerr cavity can carry metastable branches elsewhere. And the smallest nonzero $|\mathrm{Re}\,\lambda|$, the Liouvillian gap, sets the relaxation time $1/\mathrm{gap}$, which here is a sizeable fraction of the evolved window, so a time trajectory is still drifting when the null space already has the exact answer. *Derived in place (steady state is equivalent to a null vector of $\mathcal L$; every other mode decays as $e^{-|\mathrm{Re}\,\lambda|\,t}$). Verified by numerical check, by clause: the null-space dimension one and the gap value at the working parameters, the cutoff sweep of N confirming both; single parameter point, others untested; no shared code.*

<!-- box: ℒ spectrum: null space 1-dim | (unique); gap sets 1/gap approach -->

**S2 (from S1, R2, C1).** The unique steady state of S1 is sub-Poissonian: in C1's witnesses $Q<0$ and $g^{(2)}(0)<1$, a narrowing of the number distribution and antibunching of the emitted light that no coherent field has. The mechanism is the Kerr term of R2: once a photon is in the cavity it shifts the resonance away from the drive, so the drive is less able to add the next one, and the number spread is squeezed below the Poissonian value the drive alone would give. *Cited for the sub-Poissonian statistics of the driven Kerr cavity: Walls and Milburn, Quantum Optics; the witness signs derived in place from the populations. Verified by numerical check, whole node: $Q$ negative and $g^{(2)}(0)$ below one at the working parameters, with the number comb below the Poisson of equal mean; single point; no shared code.*

<!-- box: sub-Poissonian: Q < 0, g^{(2)}(0) < 1 -->

**S3 (from S1).** The same steady state has a non-negative Wigner function at these parameters: its grid minimum is zero to machine precision. So the nonclassicality lives entirely in the number statistics, not in phase-space negativity; the two witnesses are different instruments and a state can trip one without the other. *Derived in place (the Wigner minimum read off the steady state). Verified by numerical check, whole node: the grid minimum of the steady-state Wigner function is zero to machine precision at the working parameters; single point; no shared code.*

<!-- box: Wigner ≥ 0 here: | nonclassical in number only -->

## The headline

**H (from S2, S3, C1).** Drive and loss do not cancel into nothing: they settle into a unique steady state, and the Kerr nonlinearity makes that state nonclassical. The nonclassicality is sub-Poissonian number statistics, $Q<0$ and $g^{(2)}(0)<1$, the intensity-dependent resonance resisting each next photon and squeezing the number spread below Poisson; and here it is a number-statistics effect only, since the Wigner function stays non-negative, so which witness one reads (C1) decides whether the state even looks nonclassical. Kerr is the whole cause: turn it off and the state is the coherent, classical steady state of the linear driven cavity. Audit: H rests on the contract S0 together with C1 and the derived S2, S3; S0 is realized from R1, R2, R3 (and read on the cutoff R4) by the cited translation, so the analysis survives any construction that yields the same equation. No conjecture anywhere in the closure. *Derived in place from S2 and S3. Verified by numerical check, whole node, through S2's witness signs and S3's Wigner minimum, with the Kerr-off corner S2′ supplying the contrast; no shared code.*

<!-- box: nonclassical steady state; | Kerr is the whole cause -->

## When an assumption is weakened

**S2′ (from R2 weakened; S1).** Set $K=0$ and the Hamiltonian is linear: the model is a driven, leaking cavity whose unique null state, by S1's null-space construction, is the coherent state $|\beta\rangle$ with $\beta = -2F/\kappa\,i$, Poissonian with $Q=0$. This is the corner model that proves R2 is the whole source of the sub-Poissonian statistics in S2: drive and loss alone balance into a state as classical as the laser that makes it, and the nonlinearity is the entire difference. *Derived in place ($K=0$ makes the Hamiltonian linear; drive and loss balance to a coherent state). Cited for the linear-cavity steady state: standard input-output theory, Walls and Milburn. Verified by numerical check, whole node: the $K=0$ null state has $Q=0$ to rounding and QuantumDistance zero from $|{-2F/\kappa\,i}\rangle$; no shared code.*

<!-- box: K = 0: coherent state β = −2F/κ i, | Poissonian -->

**N (a recorded non-assumption on the cutoff in R4).** The cutoff $d$ is a numerical control, not a physical assumption: once $d$ clears the mean occupation the steady-state photon number stops moving, agreeing digit for digit, so the uniqueness, the gap, and the witnesses are unchanged by raising it. The boundary of the deletion: below a cutoff near the occupation the box clips the distribution's tail and the numbers still drift. *Derived in place (negligible steady-state weight above the cutoff once $d$ clears the occupation). Verified by numerical check, whole node: the steady-state photon number is unchanged from $d=12$ to $d=24$ and still moving at $d=8$.*

<!-- box: cutoff irrelevant above n̄; | sweep converges -->

## The abstract, read off the graph

A single cavity mode is driven by a laser, made nonlinear by a Kerr medium whose photon energy grows with intensity, and leaks into a zero-temperature memoryless bath. Standard open-system theory gives the Lindblad master equation, and its steady state is the null space of the Liouvillian, found in one linear-algebra step: at these parameters that null space is one dimensional, so the steady state is unique, and the Liouvillian gap sets a relaxation time long enough that watching a trajectory settle is the slow way to the same answer. That steady state is sub-Poissonian, with Mandel $Q<0$ and $g^{(2)}(0)<1$, because the intensity-dependent resonance resists adding each next photon and squeezes the number spread below Poisson; yet its Wigner function stays non-negative, so the nonclassicality is a number-statistics effect only and which witness one chooses decides whether the state looks nonclassical at all. The Kerr term is the whole cause: removing it returns the coherent, classical steady state of the linear driven cavity. The Fock cutoff is irrelevant once it clears the mean occupation.
