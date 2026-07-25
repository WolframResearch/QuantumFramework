# Part 14 Plan: Spin coupled to spatial motion

Six questions, all MSc. Class census per `Route-Table.md`: C0: 14.1, 14.2, 14.3, 14.6 (no
differential equation; WL machinery named per entry, C0 kernel facts binding); C3: 14.4
(time-dependent Schrodinger PDE, the probed spinor gate); C1: 14.5 (exactly solvable stationary
ODE, reached only after a hand reduction).

Everything in this part lives in the tensor product $L^2(\mathbb R)\otimes\mathbb C^2$: a state is
a two-component spinor $\Psi=(\psi_\uparrow,\psi_\downarrow)$ with norm
$\sum_\sigma\int\vert\psi_\sigma\vert^2\,dx=1$, and tracing out position leaves the spin density
matrix $(\rho_s)_{\sigma\sigma'}=\int\psi_\sigma\overline{\psi_{\sigma'}}\,dx$, whose purity
$\operatorname{Tr}\rho_s^2$ quantifies position-spin entanglement. Spin couples to motion through
three recurring structures: the identity $\vec L\cdot\vec S=\tfrac12(J^2-L^2-S^2)$, which converts
spin-orbit coupling into $j$-dependent numbers; the Clebsch-Gordan change of basis
$\vert j,m_j\rangle=\sum_{m_l,m_s}C^{j m_j}_{l m_l; s m_s}\vert l,m_l\rangle\vert s,m_s\rangle$;
and the Pauli Hamiltonian $H=\tfrac12(\hat p-A)^2-\tfrac{g}{2}\vec B\cdot\vec S$ with $g=2$, whose
gradient-field limit is Stern-Gerlach and whose uniform-field spectrum is Landau plus Zeeman. The
Zeeman sign is fixed by the family convention of a unit positive charge (`Shared-Brief.md`,
"Charge and sign conventions"): writing the kinetic term as $(\hat p-A)$ forces
$-\tfrac{g}{2}\vec B\cdot\vec S$. The C0 row's kernel facts bind throughout (`E` is Euler's number,
so energies and field strengths need other names; `Integrate` and `Simplify` need explicit reality
and positivity assumptions to close norms and residuals); 14.4 evolves a PDE and cites the C3
route, 14.5 cites the C1 route.

### 14.1 [MSc] How do I build a two-component spinor wavefunction $\psi_\sigma(x)$ by tensoring a spin-1/2 with the spatial state?

Build the pinned entangled spinor as a one-parameter family: normalize each branch on its own,
$\hat\psi_\uparrow=\operatorname{sech}(x-a)/\sqrt2$ (an even nodeless profile, the single bound
state of $V=-\operatorname{sech}^2x$) and
$\hat\psi_\downarrow=\sqrt{3/2}\,\operatorname{sech}(x+a)\tanh(x+a)$ (an odd one-node profile, the
excited state of the deeper well $V=-3\operatorname{sech}^2x$), then mix them with a symbolic angle,
$\Psi=\cos\theta\,\hat\psi_\uparrow\vert\uparrow\rangle+\sin\theta\,\hat\psi_\downarrow\vert\downarrow\rangle$,
whose $\theta=\pi/6$ member is the ledger's pinned state. Every spin observable follows from
`Integrate` under `Assumptions -> {a > 0, Element[x, Reals]}`: the populations are
$p_\uparrow=\cos^2\theta$, $p_\downarrow=\sin^2\theta$ with no $a$ dependence at all (both
$\int\operatorname{sech}^2u\,du=2$ and $\int\operatorname{sech}^2u\tanh^2u\,du=2/3$ are
translation invariant), so the whole geometry sits in the single real overlap
$s(a)=\int\hat\psi_\uparrow\hat\psi_\downarrow\,dx$ (closed form in $a$: verify that `Integrate`
closes it) and the purity is
$\operatorname{Tr}\rho_s^2=1-\tfrac12\sin^2(2\theta)\,[1-s(a)^2]$. The refuting contrast is the
product state $\phi(x)(\alpha\vert\uparrow\rangle+\beta\vert\downarrow\rangle)$, where the identical
pipeline must return purity exactly 1, with $\operatorname{Tr}\rho_s=1$ standing on both and an
`NIntegrate` spot check of $s$ at $a=1$ reusing the definitions. Two limits carry the physics:
$s$ vanishes both at $a=0$ (odd integrand, exact parity zero) and as $a\to\infty$ (disjoint
supports), so the entanglement is maximal at zero displacement and at infinite displacement alike
and is weakest in between, and the floor is $\tfrac12$ only at $\theta=\pi/4$, while the pinned
$\theta=\pi/6$ member floors at $\tfrac58$. Close on the measurement reading: an apparatus probing
spin alone sees the mixed $\rho_s$, and the missing purity is exactly the which-branch information
stored in position, which parity can supply as completely as distance does.

### 14.2 [MSc] How do I add the spin-orbit term $\propto \vec L\cdot\vec S$ and compute the hydrogen fine structure?

Assemble the pinned hydrogen $2p$ fine structure from its radial and angular factors. Radially,
define $R_{nl}(r)$ from `LaguerreL` (the associated-Laguerre argument convention and the resulting
normalization constant: verify at authoring) and earn the norms by `Integrate` with the $r^2$
measure; then compute $\langle r^{-3}\rangle_{nl}$ from those same definitions over the family
$(n,l)\in\{(2,1),(3,1),(3,2),(4,1)\}$ and check that all four reproduce
$1/[n^3\,l(l+\tfrac12)(l+1)]$, a definition-reusing test that exposes a wrong radial function or a
wrong normalization at once and hands over the $n^{-3}$ scaling of the whole fine structure;
$(2,1)$ gives $\langle r^{-3}\rangle_{2p}=\tfrac1{24}$. Angularly, build the $6\times6$
$\vec L\cdot\vec S=\tfrac12(J^2-L^2-S^2)$ from `KroneckerProduct` of the $l=1$ matrices and
`PauliMatrix[i]/2`, and let `Eigenvalues` return $\tfrac12$ four times and $-1$ twice, matching
$\tfrac12[j(j+1)-l(l+1)-s(s+1)]$ at $j=\tfrac32,\tfrac12$ with multiplicities $2j+1$; the two
trace identities $\operatorname{Tr}(\vec L\cdot\vec S)=0$ and
$\operatorname{Tr}[(\vec L\cdot\vec S)^2]=3$ pin the assembly independently of the eigensolver.
With $\xi(r)=\alpha^2/(2r^3)$ the splitting is
$\tfrac32\cdot\tfrac{\alpha^2}{2}\langle r^{-3}\rangle_{2p}=\alpha^2/32$ hartree. Note the edge the
family exposes: at $l=0$ the radial factor diverges while $\langle\vec L\cdot\vec S\rangle$
vanishes, so the product is indeterminate and the $2s$ shift is carried by the Darwin term, not by
spin-orbit. Close in the lab: $\alpha^2/32$ hartree is about $10.9\,\mathrm{GHz}$, the measured
fine-structure doubling of the Lyman-$\alpha$ line.

### 14.3 [MSc] How do I form the total angular momentum $\vec J=\vec L+\vec S$ and change to the coupled basis?

Construct $J_i=L_i\otimes\mathbb 1_2+\mathbb 1_3\otimes S_i$ with `KroneckerProduct`, the $l=1$
matrices built from the ladder elements $\sqrt{l(l+1)-m(m\pm1)}$ via `SparseArray` and
$S_i=$ `PauliMatrix[i]/2`, and certify the algebra first: `Simplify` of $[J_x,J_y]-iJ_z$ to the
zero matrix. Then form the $6\times6$ change-of-basis matrix $U$ whose rows are
`ClebschGordan[{1, ml}, {1/2, ms}, {j, mj}]` coefficients (argument order and Condon-Shortley
phase: verify at authoring) and let exact linear algebra deliver the payoff: $U.U^\dagger=\mathbb 1$,
and $U.J^2.U^\dagger$, $U.J_z.U^\dagger$ exactly diagonal, symbolically, with entries $\tfrac{15}4$
four times and $\tfrac34$ twice beside the matching $m_j$, every off-diagonal an exact 0. Two
independent refutations: build the same basis by ladder algebra, starting from the stretched state
$\vert\tfrac32,\tfrac32\rangle=\vert1,1\rangle\vert\uparrow\rangle$, stepping down with $J_-$ and
renormalizing, then taking the orthogonal complement inside each $m_j$ block for the $j=\tfrac12$
pair, and compare row by row (a phase-convention slip shows up as a sign mismatch); and check the
general $s=\tfrac12$ closed form
$\vert l\pm\tfrac12,m_j\rangle=\pm\sqrt{\tfrac{l\pm m_j+1/2}{2l+1}}\vert m_j-\tfrac12\rangle\vert\uparrow\rangle+\sqrt{\tfrac{l\mp m_j+1/2}{2l+1}}\vert m_j+\tfrac12\rangle\vert\downarrow\rangle$
(overall sign convention: verify at authoring) against `ClebschGordan` over a `Table` of
$l\in\{1,2,3\}$ and every $m_j$, since `ClebschGordan` at symbolic $l,m_j$ may simply not evaluate.
Close with the forward pointer: these coupled columns are the weak-field eigenvectors 14.6
recovers, and their multiplet structure is what 14.2's eigenvalue count already displayed.

### 14.4 [MSc] How do I model Stern-Gerlach as a spin-dependent spatial deflection of a spinor wave packet in a field gradient?

Evolve the pinned two-component sech packet $\Psi(x,0)=\tfrac{1}{\sqrt2}\psi_0(x)(1,1)$, with
$\psi_0=\operatorname{sech}(x)/\sqrt2$, under
$i\partial_t\Psi=[-\tfrac12\partial_x^2+\lambda x\,\sigma_z]\Psi$ at the probed gradient
$\lambda=0.2$ to $T=8$, using the C3 spec in full: `NDSolveValue` on the coupled pair with
`Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> n, "MaxPoints" -> n, "DifferenceOrder" -> "Pseudospectral"}}`
and `AccuracyGoal -> 10, PrecisionGoal -> 10`, the point pins giving a fixed $dx$ for the
$dx$-weighted grid sums that must replace `NIntegrate` on oscillatory interpolants, and tight goals
being load-bearing (defaults are silently wrong at the $10^{-3}$ level and refinement at default
goals is non-monotone). Size the box in both position and momentum: on $[-40,40]$ with $n=513$ the
walls sit far beyond the deflected centers $\mp\lambda t^2/2=\mp6.4$ plus the dispersive spread,
while $k_{\max}=\pi/dx\approx20$ clears the momentum the tilt pumps in, $\vert\langle p\rangle_\sigma\vert=\lambda T=1.6$,
plus several $k$-widths of $\tilde\psi_0\propto\operatorname{sech}(\pi k/2)$; the linear potential
is not periodic, so the periodic identification is honest only while $\psi$ is negligible at both
walls, and the C3 wrap-around contamination is structural and invisible to norm conservation, so
carry an $n$-doubling and $T$-lengthening sweep (the probe evidence reaches only $t\le2$ at this
$\lambda$). The deliverable is the measurement record, the coherence
$c(t)=\int\psi_\uparrow\overline{\psi_\downarrow}\,dx$ and the spin purity falling from 1 while
both populations hold at $\tfrac12$, so the check must bite on shape, not on centroids: the tilt is
a momentum translation, giving the exact momentum-space form
$c(t)=\int\tilde\psi_0(k-\lambda t)\overline{\tilde\psi_0(k+\lambda t)}\,e^{ik\lambda t^2}\,dk$
(closed form for the sech pair: verify that `Integrate` closes it), and the C3 cross-check route,
a hand-built split-step Fourier evolution of the same coupled pair, must reproduce $\vert c(t)\vert$
itself. The probed $\Omega=0$ Ehrenfest anchor
($\langle x\rangle_\uparrow(2)=-0.39999$ against the exact $-0.4$, total norm drift
$5\times10^{-8}$) stays as a cheap standing sanity, but it is exact for any initial packet in a
linear potential and therefore cannot see a wrong width, wrong dispersion, or an aliased tail.
Close on measurement as dynamics: the coherence dies mainly through momentum separation, so once
$\vert c\vert$ is negligible the beams are which-path records and no downstream spin rotation
revives the interference, the projection postulate emerging from unitary evolution plus a
traced-out position.

### 14.5 [MSc] How do I solve the Pauli equation for a nonrelativistic spin in an electromagnetic field?

Take uniform $\vec B=B\hat z$ in the Landau gauge $A=(0,Bx,0)$, where the Pauli Hamiltonian
$H=\tfrac12(\hat p-A)^2-\tfrac{g}{2}BS_z$ (the family sign convention for a unit positive charge)
is diagonal in $S_z$: split into the two spin blocks and translate along $y$ with the ansatz
$e^{iky}f(x)$. Per the C1 verdict this member is covered only after the hand reduction: complete the
square to a shifted oscillator of frequency $B$ centered at $x_0=k/B$ as hand algebra, certify the
reduction with a `Simplify` residual on the operator identity, and only then hand the oscillator to
the C1 route, `DSolve` with symbolic energy followed by the manual series-termination read-off and
`HermiteH` eigenfunctions with a `FullSimplify` residual of zero; the read-off is not a stylistic
choice, `DSolve` handed the full decay conditions hangs silently rather than failing (the C1 gate
timed out at 60 s with no message), and the whole-line `DEigensystem` domain form is a free exact
confirmation only if it accepts a symbolic $B$ (the gate probed it at numeric coefficients: verify
at authoring). The spectrum $E_{n,m_s}=B(n+\tfrac12)-\tfrac{g}{2}Bm_s$ becomes at $g=2$ the two
interleaved ladders $E_{n,\uparrow}=Bn$ and $E_{n,\downarrow}=B(n+1)$, so `Simplify` of
$E_{n,\downarrow}-E_{n+1,\uparrow}$ is 0 symbolically in $B$ and $n$, and the discriminating
refutation is that `Solve` on $E_{n,\downarrow}=E_{n+1,\uparrow}$ for $g$ returns exactly $g=2$:
the coincidence is a property of the gyromagnetic ratio, not of magnetism. Add one numeric
confirmation at $B=1$, an `NDEigenvalues` run per spin block on a truncated box with explicit
`DirichletCondition` (omitting boundary conditions is a silent Neumann-zero answer, wrong at the
percent level on a ground state, per the probed gates), whose sorted merge is
$\{0,1,1,2,2,\dots\}$: the coincidence appears as numerical degeneracy. The infinite per-level
Landau degeneracy was not probed and its flux counting belongs to 13.2; here only the $k$
independence of $E_{n,m_s}$, read directly off the closed form, records it. Close at the bottom of
the ladder: at $g=2$ the state $(n=0,m_s=+\tfrac12)$ sits at exactly zero energy, the zero-point
cost cancelled by the Zeeman gain, the nonrelativistic shadow of the Dirac value of $g$ that 22.3
derives.

### 14.6 [MSc] How do I compute the anomalous Zeeman effect and the Paschen-Back crossover, where spin-orbit and the magnetic field compete?

On the hydrogen $2p$ manifold build the exact $6\times6$
$H=\xi\,\vec L\cdot\vec S+\tfrac{B}{2}(L_z+2S_z)$ from `KroneckerProduct` with $\xi$ and $B$
symbolic (rebuild the $l=1$ and Pauli matrices in-entry; the answer stays self-contained). Certify
the block structure first, `Simplify` of $[H,J_z]$ to the zero matrix, then make the blocks explicit
rather than hoping the eigensolver finds them: order the uncoupled basis by $m_j=m_l+m_s$ with
`Ordering` (or conjugate by the coupled-basis matrix of 14.3) so $H$ becomes block diagonal with
sizes 1, 2, 2, 1, and take `Eigenvalues` on the $1\times1$ and $2\times2$ blocks, where the
Breit-Rabi square roots are guaranteed closed forms in $\xi$ and $B$; a direct `Eigenvalues` on the
full two-parameter $6\times6$ may return `Root` objects instead (verify at authoring). The two
limits are the refuting checks because they rest on independent physics: `Series` in $B$ at fixed
$\xi$ must reproduce the anomalous pattern $E_{j,m_j}\approx E_j+g_j m_j B/2$ with Lande factors
$g_{3/2}=\tfrac43$ and $g_{1/2}=\tfrac23$, matching the coupled-basis diagonal elements of
$(L_z+2S_z)/2$, while the opposite `Series` in $\xi$ at fixed $B$ must land on the Paschen-Back
pattern $B(m_l+2m_s)/2$ with its characteristic degeneracies; and the coupling-resolved identity
$\operatorname{Tr}(H^2)=3\xi^2+\tfrac52B^2$, exact and symbolic in both, catches an assembly error
that a traceless sum rule would let through. `Plot` the six branches against $B$ at fixed $\xi$ to
exhibit the crossover near $B\sim\xi$: branches of different $m_j$ cross freely while the paired
same-$m_j$ branches avoid crossing, two-level repulsion inside each block. Close in the
spectroscope: weak-field lines splitting by $g_j$ rather than by 1 is the historical anomalous
Zeeman effect, and the field where the pattern reorganizes, Zeeman energy comparable to the
$\alpha^2/32$ splitting of 14.2, is about $0.8\,\mathrm{T}$ for hydrogen.
