# Part 18 Plan: Open quantum systems in continuous space

5 questions. Class census per the class-census table in `Route-Table.md`: C5 (truncated-basis
ODE-IVP systems): all of 18.1 through 18.5; 18.1 is the class's probed representative, so entries
cite the C5 gates directly.

## Common ground

Everything in this part rests on the Lindblad master equation
$\dot\rho=-i[H,\rho]+\sum_k\left(L_k\rho L_k^\dagger-\tfrac12\{L_k^\dagger L_k,\rho\}\right)$ for
the reduced density operator of a mode in contact with an environment, and four of its
consequences: trace preservation and complete positivity, which hold for every Lindblad generator
and fail for non-Lindblad approximations such as Caldeira-Leggett; the exact damped-coherent
benchmark, that under $H=\hat n$, $L=\sqrt\gamma\,a$ a coherent state stays coherent with
$\alpha(t)=\alpha e^{-(\gamma/2+i)t}$ and Poisson populations of mean $|\alpha|^2e^{-\gamma t}$;
the equivalence, for quadratic systems, of the operator master equation to a phase-space
drift-diffusion (Fokker-Planck) flow with exactly closing Gaussian moments; and the unravelling
identity, that the deterministic $\rho(t)$ is the ensemble average of stochastic pure-state
trajectories punctuated by quantum jumps.

All five entries share the certified contract of the C5 verdict in `Route-Table.md`: the Fock
space is truncated at $N$, the ladder operator is a `SparseArray` with $a_{n,n+1}=\sqrt{n+1}$, and
the vectorized Liouvillian is assembled as a `SparseArray` of `KroneckerProduct` terms under the
row-major convention $\mathrm{vec}(A\rho B)=(A\otimes B^{T})\,\mathrm{vec}\,\rho$; evolution goes
through the `MatrixExp[L t, v]` action form, exact in time, so truncation is the only error and it
decays with the Poisson tail mass; the independent cross-check is `NDSolveValue` on the
matrix-valued ODE $\rho'(t)=\mathcal L[\rho(t)]$, with the probed caveat that route agreement is
blind to truncation (agreement $1.1\times10^{-7}$ while both routes were wrong by
$5.6\times10^{-3}$ at $N=8$), so the exact benchmark and the $N$-sweep are mandatory wherever they
exist; default-tolerance norm and trace drift grows with integration length, and the positivity
floor at default tolerances is $\sim-10^{-6}$, so every positivity assertion uses a
tolerance-aware threshold.

### 18.1 [MSc] How do I solve the Lindblad master equation for a damped oscillator?

Damp the coherent state $\alpha=2$ under $H=\hat n$, $L=\sqrt\gamma\,a$, with $\gamma$ symbolic in
the assembly and pinned to $\gamma=\tfrac25$ for the runs: this is the C5 verdict's own probed
representative, and the coherent state is the pointer-like carrier whose stability under damping
is the physics. Assemble the vectorized Liouvillian per the shared contract, evolve by
`MatrixExp[L t, v]`, and cross-check with `NDSolveValue` on the matrix ODE. The refuting check is
the exact benchmark, reusing the stored $a$ and $\rho(t)$: $\mathrm{Tr}[a\,\rho(t)]$ must equal
$\alpha e^{-(\gamma/2+i)t}$ and the diagonal populations must stay Poisson with mean
$|\alpha|^2e^{-\gamma t}$, so the energy $\mathrm{Tr}[\hat n\,\rho(t)]$ decays at exactly
$\gamma$. Because route agreement is blind to truncation, run the $N$-sweep from $N=8$ to $40$
against the benchmark and confirm the super-exponential tail-mass decay the verdict records
(population error $5.4\times10^{-3}$ at $N=8$ falling to $10^{-13}$ by $N=40$), while trace error
and the minimum eigenvalue of $\rho$ (`Eigenvalues`) are monitored at tolerance-aware thresholds.
Close on the pointer reading: damping sends $|\alpha\rangle$ to another coherent state with purity
intact, the stability that 18.3 shows failing spectacularly for superpositions of two of them.

### 18.2 [MSc] How do I model quantum Brownian motion (Caldeira-Leggett) on a truncated basis?

Evolve a displaced position-squeezed packet ($x_0=2$, squeezing $r=1$, so $\Delta x$ sits well
below the thermal length; the ledger pins "Caldeira-Leggett packet", and this specification sits
inside that class) in the trap $H=\tfrac12(\hat p^2+\hat x^2)$ under the high-temperature
Caldeira-Leggett generator
$\dot\rho=-i[H,\rho]-i\gamma[\hat x,\{\hat p,\rho\}]-2\gamma T[\hat x,[\hat x,\rho]]$ with
$\gamma=\tfrac1{10}$, $T=10$, $N=30$ (prefactor conventions verify at authoring), building
$\hat x=(a+a^\dagger)/\sqrt2$ and $\hat p=i(a^\dagger-a)/\sqrt2$ from the shared `SparseArray`
ladder, vectorizing each commutator term by the row-major convention, evolving by `MatrixExp`
action with the `NDSolveValue` cross-check. The class trap inverts here, per the C5 members-sanity
note: this generator is NOT completely positive, so the transient negativity of $\rho$ is physics,
not solver error, and the positivity monitor must not be used as a correctness gate; instead
exhibit the early-time minimum eigenvalue of $\rho(t)$ dipping far below the $-10^{-6}$ tolerance
floor and let the honest failure regime carry the lesson. The refuting check moves to the moments:
for quadratic $H$ the first and second moments close as a small linear ODE system that
`DSolveValue` solves exactly (damped $\langle p\rangle$, diffusion feeding $\Delta p^2$ toward its
thermal value), and the truncated-basis moments $\mathrm{Tr}[\hat x\,\rho(t)]$,
$\mathrm{Tr}[\hat p\,\rho(t)]$, reusing the stored operators, must land on those closed forms with
the trace pinned at 1 and the $N$-sweep guarding truncation. Close with the repair: adding the
minimal $[\hat p,[\hat p,\rho]]$ diffusion term (coefficient verify at authoring) restores
Lindblad form, and the same monitor that indicted Caldeira-Leggett now certifies positivity,
separating approximation from artifact.

### 18.3 [MSc] How does dephasing select a pointer basis and decohere a spatial superposition (einselection)?

Damp the superposition $\mathcal N(|\alpha_1\rangle+|\alpha_2\rangle)$ of two coherent states with
$\alpha_1=2$, $\alpha_2=-1$ (an asymmetric pair, distinct from 16.4's even cat) under $H=\hat n$,
$L=\sqrt\gamma\,a$, $\gamma=\tfrac15$, $N\geq30$, built and evolved per the shared contract with
the `NDSolveValue` cross-check and $N$-sweep: the channel that drains energy through $a$ acts on
this state as dephasing in the coherent-state basis it selects. Track two rates on the same run,
reusing the stored $\rho(t)$: the energy $\mathrm{Tr}[\hat n\,\rho(t)]$ relaxing at $\gamma$, and
the off-diagonal coherence weight $|\langle\beta_2(t)|\rho(t)|\beta_1(t)\rangle|$ between the
drifting pointer components $\beta_i(t)=\alpha_i e^{-(\gamma/2+i)t}$, whose suppression factor is
$\exp[-\tfrac12|\alpha_1-\alpha_2|^2(1-e^{-\gamma t})]$ with short-time rate
$\Gamma_{\mathrm{dec}}=\tfrac\gamma2|\alpha_1-\alpha_2|^2$ (closed forms verify at authoring). The
refuting check is the scaling law itself: sweep the separation $d=|\alpha_1-\alpha_2|$ over
$\{\tfrac32,2,3,4\}$ at fixed midpoint with `Table`, extract each rate from the initial slope of
the log-coherence, and the rates must fall on $\tfrac\gamma2 d^2$ while the energy rate stays
pinned at $\gamma$; that quadratic-over-constant contrast is einselection, and a wrong jump of
machinery anywhere would bend the fitted exponent off 2. Close by cashing to the lab: at
macroscopic separation the same law kills a cat's coherence orders of magnitude before its energy
moves, which is why the world looks classical in the basis the damping selects.

### 18.4 [MSc] How do I cast the optical master equation as a Fokker-Planck equation in the Wigner representation?

Take the optical master equation for a mode in thermal contact, $H=\hat n$ with
$L_1=\sqrt{\gamma(\bar n+1)}\,a$ (emission) and $L_2=\sqrt{\gamma\bar n}\,a^\dagger$ (absorption),
at $\gamma=\tfrac15$, $\bar n=2$. The derivation is symbolic algebra, since per the C5
members-sanity note this member is probed only through the small-system `DSolve` gate and the
Fokker-Planck framing itself is unprobed: apply the Wigner correspondence rules
($a\rho\leftrightarrow(\alpha+\tfrac12\partial_{\bar\alpha})W$ and companions, verify at
authoring) with `D` on a generic $W(x,p,t)$ and `Simplify`, and read off the drift and diffusion
coefficients of the resulting form: drift is the classical rotation from $H$ plus an isotropic
contraction at $\gamma/2$, diffusion is isotropic with $D=\tfrac\gamma2(\bar n+\tfrac12)$ (exact
coefficients verify at authoring). Then close the Gaussian moment ODEs exactly, the probed gate:
$\langle x\rangle$, $\langle p\rangle$ and the three second moments form a small linear system
`DSolveValue` solves in closed form, relaxing to the thermal Gaussian of variance
$\bar n+\tfrac12$, and residual-check that this Gaussian annihilates the derived right-hand side
(`D` plus `Simplify` to 0, reusing the derived drift and diffusion, so a wrong coefficient cannot
survive). The refuting cross-check is the truncated-basis route: the same moments as
$\mathrm{Tr}[\hat x\,\rho(t)]$ and companions from the shared `SparseArray` plus `MatrixExp`
machinery must land on the `DSolveValue` closed forms, with the $N$-sweep guarding truncation.
Close on the reading: the mode performs a classical Ornstein-Uhlenbeck relaxation whose noise
floor $\bar n+\tfrac12$ never reaches zero, the $\tfrac12$ being vacuum fluctuation that survives
at $\bar n=0$.

### 18.5 [MSc] How do I unravel a damped mode into quantum trajectories (the quantum-jump picture)?

Unravel the damped Fock state $|2\rangle$ under $H=\hat n$, $L=\sqrt\gamma\,a$, $\gamma=\tfrac12$:
the dynamics is confined to $\mathrm{span}\{|0\rangle,|1\rangle,|2\rangle\}$, so truncation is
exact and any ensemble discrepancy indicts the jump logic itself. Run the probed functional loop:
precompute the non-Hermitian step `MatrixExp[-I Heff dt]` with
$H_{\mathrm{eff}}=\hat n-\tfrac{i\gamma}2\,a^\dagger a$, then `FoldList` the unnormalized state
over the time grid, firing a jump (apply $a$, renormalize, redraw the `RandomReal` threshold)
whenever the squared norm crosses the threshold. The probed warning binds: a coherent-state
ensemble test is VACUOUS, $|\alpha\rangle$ being an eigenstate of the jump operator (probed
trajectory spread $3\times10^{-6}$), which is exactly why the ledger pins the Fock start. The
ensemble-vs-master-equation check has teeth here: the mean of $\langle n\rangle(t)$ over roughly
300 trajectories must match $2e^{-\gamma t}$ from the shared vectorized-Liouvillian `MatrixExp`
route within the Monte Carlo error bar; per trajectory, the integer staircase $2\to1\to0$ (a Fock
state stays Fock between jumps; probed at $4.4\times10^{-16}$) and the waiting-time statistics,
the first jump exponential at rate $2\gamma$ because the no-jump norm of $|n\rangle$ decays as
$e^{-n\gamma t}$, with mean jump count $2(1-e^{-\gamma t})$ by photon bookkeeping. Cross-check a
single trajectory on the verified `WhenEvent` plus `DiscreteVariables` machinery inside `NDSolve`
(probed clean), and sweep $dt$ at fixed ensemble since jump-time resolution enters at $O(dt)$, a
flagged open risk. Close on the punchline: each trajectory is what a photodetector records, clicks
and staircases, and the smooth exponential of 18.1 is nothing but their average.
