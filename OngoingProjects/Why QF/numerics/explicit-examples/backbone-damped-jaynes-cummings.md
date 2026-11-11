---
Template: Default
Title: "Backbone: When Does the Photon Come Back? Dissipative Jaynes-Cummings"
Author: Mads Bahrami
---

# Backbone: When Does the Photon Come Back? Dissipative Jaynes-Cummings

This page is the statement layer of the essay *Damped-Jaynes-Cummings.md*: a graph of assumptions, choices, and derived statements, each given in words and formulas only. Every node names its parents, so the dependency graph is recoverable from the page, and the layers order the reading from the setup to the headline and then to what holds when an assumption is weakened. Three conventions hold throughout. The page has two layers: a construction layer, where the physical roots are translated into the equation S0 by cited standard theory, and an analysis layer, where every consequence consumes S0 as a contract and therefore survives any construction that yields the same equation; the equation is never a root here, because the weakening below (switch off atomic decay) reaches beneath it. Every assumption must be consumed by at least one statement that is not its own weakening; a root failing that test is not an assumption of the theory but only a direction along which the model can be varied. And every statement's tail separates provenance (derived in place, cited, or conjectured) from verification (exact decision, analytic cross-check, numerical check, or none), naming the clause scope. Node labels: R root, C choice, S derived statement with S0 reserved for the equation of motion, H headline, N recorded non-assumption, prime a weakened variant. Nothing here executes; the verifications all live in the essay, whose cells are the only computation, and are numerical throughout, as befits a "numeric explicit examples" showcase.

## The setup, and the choice

**R1 (assumption).** A two-level atom and one cavity mode exchange a quantum through the Jaynes-Cummings coupling in the rotating-wave approximation, $\hbar g(\hat\sigma_-\hat a^\dagger + \hat\sigma_+\hat a)$, on resonance $\nu=\omega$. The coupling conserves the total excitation number $\hat a^\dagger\hat a + \hat\sigma_+\hat\sigma_-$: every photon the field gains is a quantum the atom has lost. On resonance the free terms contribute only a phase within each excitation manifold, so the excited-state population is set by the coupling and the losses alone.

<!-- box: JC exchange g, RWA, | resonant; excitation conserved -->

**R2 (assumption).** The atom decays into modes outside the cavity: a jump operator $\hat\sigma_-$ at rate $\gamma$, lowering the atom without touching the field, into a zero-temperature Markovian bath.

<!-- box: atomic decay, jump σ₋ at rate γ -->

**R3 (assumption).** The cavity mirror leaks: a jump operator $\hat a$ at rate $\kappa$, removing a photon without touching the atom, into a zero-temperature Markovian bath.

<!-- box: cavity loss, jump a at rate κ -->

**C1 (choice).** The initial excitation to track. A single quantum, $|e,0\rangle$, closes the problem on three states and makes the excited amplitude follow a $2\times2$ generator exactly; a multi-photon start, $|e,n\rangle$, sets the exchange frequency by the occupied rung. The essay computes both, and reads the same competition of coupling against loss off each. *Derived in place (the one-excitation sector is closed under the coupling and the losses; a higher rung fixes the Rabi frequency). Verified by numerical check, whole node: both initial states evolved under the same master equation; no shared code.*

<!-- box: initial excitation: one quantum | vs a multi-photon start -->

## The equation, derived from the setup

**S0 (from R1, R2, R3; cited).** Standard open-system theory turns the setup into the Lindblad master equation with two jump operators,
$$\frac{d\hat\rho}{dt} = -\frac{i}{\hbar}[\hat H,\hat\rho] + \mathcal D[\sqrt{\gamma}\,\hat\sigma_-]\hat\rho + \mathcal D[\sqrt{\kappa}\,\hat a]\hat\rho, \qquad \hat H = \hbar\nu\,\hat a^\dagger\hat a + \tfrac12\hbar\omega\,\hat\sigma_z + \hbar g(\hat\sigma_-\hat a^\dagger + \hat\sigma_+\hat a),$$
with each root visible: the JC coupling R1 is the exchange term in $\hat H$, atomic decay R2 is the dissipator on $\hat\sigma_-$, and cavity loss R3 is the dissipator on $\hat a$. *Cited: the Jaynes-Cummings model, Jaynes and Cummings, Proc. IEEE 51, 89 (1963); the dissipative extension and its two jump operators, Carmichael, Statistical Methods in Quantum Optics, and Breuer and Petruccione, chapter 3. Verification: none; the translation is trusted to its sources, and the essay computes with exactly this equation.*

<!-- box: dρ/dt = −(i/ħ)[H,ρ] | + D[√γ σ₋]ρ + D[√κ a]ρ -->

## Consequences

**S1 (from S0, C1).** In the single-excitation sector of C1 the coupling connects only $|e,0\rangle$ and $|g,1\rangle$, and either loss drains to $|g,0\rangle$, which is dark: it has no excitation and the coupling cannot lift it. Because the recycling term of S0 feeds $|g,0\rangle$ but never $|e,0\rangle$, the excited amplitude follows the non-Hermitian no-jump generator alone,
$$\frac{d}{dt}\begin{pmatrix}c_e\\c_1\end{pmatrix} = \begin{pmatrix}-\gamma/2 & -ig\\ -ig & -\kappa/2\end{pmatrix}\begin{pmatrix}c_e\\c_1\end{pmatrix}, \qquad P_e = |c_e|^2.$$
The full open-system dynamics of the excited population is captured exactly by this $2\times2$ picture. *Derived in place ($|g,0\rangle$ is dark and the recycling term never feeds $|e,0\rangle$). Verified by numerical check, whole node: the full master-equation excited population agrees with the $2\times2$ result to the integrator's accuracy across all three regimes; no shared code.*

<!-- box: one excitation: excited amplitude | obeys a 2×2 no-jump generator -->

**S2 (from S1).** The eigenvalues of S1's generator are $\tfrac14\!\left(-\gamma-\kappa \pm \sqrt{(\gamma-\kappa)^2 - 16g^2}\right)$, and the sign of the discriminant $(\gamma-\kappa)^2 - 16g^2$ splits the excited population into three regimes. Strong coupling $g>|\kappa-\gamma|/4$: a complex pair, so $P_e$ is a vacuum Rabi oscillation under a decaying envelope, the photon emitted and reabsorbed several times before it escapes. Crossover $g=|\kappa-\gamma|/4$: the eigenvalues merge, the oscillation disappears, critically damped. Purcell regime $\kappa\gg g$: the cavity empties before reabsorption, so the atom decays as a pure exponential at the enhanced rate $\gamma + 4g^2/\kappa$ rather than its bare $\gamma$. *Derived in place (the eigenvalues of the $2\times2$ generator and the sign of the discriminant). Verified by analytic cross-check (the closed-form eigenvalues) and numerical check (the three regimes reproduced by the full evolution, the Purcell curve matching $e^{-(\gamma + 4g^2/\kappa)t}$); one point per regime; no shared code.*

<!-- box: discriminant (γ−κ)² − 16g² | splits into three regimes -->

**S3 (from S0, C1).** From the multi-photon start $|e,n\rangle$ of C1 the excited population is a Rabi oscillation at $2g\sqrt{n+1}$, the frequency of the $|e,n\rangle\leftrightarrow|g,n+1\rangle$ pair, under an envelope that decays toward the joint ground state $|g,0\rangle$ that neither loss can leave. The oscillation is the coherent exchange, the envelope is the loss, and loss sets how long the exchange survives, not its frequency. *Derived in place (the occupied rung oscillates at $2g\sqrt{n+1}$; the two losses set the envelope). Verified by numerical check, whole node: the damped oscillation at three loss strengths, collapsing sooner as the rate grows; no shared code.*

<!-- box: multi-photon: damped Rabi | at 2g√(n+1) -->

## The headline

**H (from S2, S3).** When does the photon come back? One dimensionless comparison decides it: the coupling against the loss, $g$ versus $|\kappa-\gamma|/4$. When the coupling wins, the excitation is shared coherently and $P_e$ oscillates under decay, the photon reabsorbed before it escapes; at the crossover the oscillation is critically damped away; when the loss wins, the excitation leaves before it can return and the atom simply decays, enhanced above its bare rate by the cavity (the Purcell effect). The same comparison governs the multi-photon ladder: a decaying envelope over a Rabi oscillation whose frequency is set by the occupied rung, not by the loss. Audit: H rests on the contract S0 together with C1 and the derived S2, S3; S0 is realized from R1, R2, R3 by the cited translation, so the analysis survives any construction that yields the same equation. No conjecture anywhere in the closure. *Derived in place from S2 and S3 (the competition $g$ versus $|\kappa-\gamma|/4$). Verified by numerical check, whole node, through S2's three regimes and S3's damped oscillation; no shared code.*

<!-- box: one comparison g vs |κ−γ|/4 | decides if the exchange survives -->

## When an assumption is weakened

**S2′ (from R2 weakened; S2).** Switch off atomic decay, $\gamma=0$: S2's eigenvalues send the Purcell rate $\gamma + 4g^2/\kappa$ to the cavity-limited $4g^2/\kappa$, while the three-way split by $g$ versus $\kappa/4$ survives intact. This separates the enhanced decay into its two pieces, the bare atomic rate owned by R2 and the cavity-induced $4g^2/\kappa$ owned jointly by the coupling R1 and the cavity loss R3, and is the corner model that holds R2 apart from R3. *Derived in place (setting $\gamma=0$ in S2's eigenvalues). Verified by analytic cross-check against S2's closed form at $\gamma=0$.*

<!-- box: γ = 0: Purcell rate → 4g²/κ; | split survives -->

**N (a recorded non-assumption on the equal rates).** Nothing on this page needs $\gamma=\kappa$. The multi-photon plot of S3 fixes equal rates only to read one dial, but S3's frequency $2g\sqrt{n+1}$ and its decaying envelope, and S2's three regimes, all hold for independent $\gamma$ and $\kappa$; the headline's comparison is written in $|\kappa-\gamma|$ precisely because the two rates are free. *Derived in place (S2 and S3 carry $\gamma$ and $\kappa$ independently throughout). Verified by numerical check: the three-regime comparison uses unequal rates, the damped-Rabi plot equal ones, with no statement resting on the equality.*

<!-- box: equal rates γ = κ not needed; | both rates free -->

## The abstract, read off the graph

A two-level atom and one cavity mode exchange a quantum through the resonant Jaynes-Cummings coupling, while the atom decays to the outside and the cavity mirror leaks, two zero-temperature Markovian channels draining the exchange. Standard theory gives the Lindblad master equation with those two jump operators. Started with a single quantum, the problem closes on three states with a dark joint ground state, so the excited amplitude obeys a non-Hermitian $2\times2$ no-jump generator exactly, and the sign of its discriminant $(\gamma-\kappa)^2 - 16g^2$ splits the behavior into three regimes: strong coupling with a vacuum Rabi oscillation under decay, a critically damped crossover, and the Purcell regime where the atom decays as a pure exponential at the cavity-enhanced rate $\gamma + 4g^2/\kappa$. Started higher up the ladder, the excited population is a damped Rabi oscillation at $2g\sqrt{n+1}$ under a decaying envelope. One dimensionless comparison, the coupling $g$ against the loss $|\kappa-\gamma|/4$, decides in every case whether the photon comes back. Switching off atomic decay leaves the three-way split intact and reduces the Purcell rate to the cavity-limited $4g^2/\kappa$; the equality of the two rates is never needed.
