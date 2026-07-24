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
C0 row's kernel facts bind throughout (`E` is Euler's number, so energies and field strengths need
other names; `Integrate` and `Simplify` need explicit reality and positivity assumptions to close
norms and residuals); 14.4 evolves a PDE and cites the C3 route, 14.5 cites the C1 route.

### 14.1 [MSc] How do I build a two-component spinor wavefunction $\psi_\sigma(x)$ by tensoring a spin-1/2 with the spatial state?

Build the pinned entangled spinor
$\Psi=N[\operatorname{sech}(x-a)\vert\uparrow\rangle+\operatorname{sech}(x+a)\tanh(x+a)\vert\downarrow\rangle]$,
the Poschl-Teller ground-state shape riding spin up at $+a$ and its nodal first-excited shape
riding spin down at $-a$, as an explicit two-component vector; every spin observable then follows
from the overlap integrals $(\rho_s)_{\sigma\sigma'}=\int\psi_\sigma\overline{\psi_{\sigma'}}\,dx$
by `Integrate` under `Assumptions -> {a > 0, Element[x, Reals]}`. Normalize, assemble $\rho_s$ from
the branch weights $p_\sigma$ and the cross overlap $c$, and compute the purity
$\operatorname{Tr}\rho_s^2=p_\uparrow^2+p_\downarrow^2+2\vert c\vert^2<1$, keeping $a$ symbolic as
far as `Integrate` closes the overlaps (closed form in $a$: verify at authoring). The contrast is
the product state $\phi(x)(\alpha\vert\uparrow\rangle+\beta\vert\downarrow\rangle)$: the identical
pipeline must return purity exactly 1 there, the check that could refute the entanglement claim,
with $\operatorname{Tr}\rho_s=1$ as the standing sanity on both and an `NIntegrate` spot check of
$c$ at $a=1$ reusing the definitions. Read off two regimes: as $a\to\infty$ the purity falls to
$p_\uparrow^2+p_\downarrow^2$ (the $\tfrac12$ floor at equal weights), while even at $a=0$ the
cross overlap vanishes by parity (odd integrand), so the spin stays mixed at zero displacement:
entanglement tracks spatial distinguishability, not distance. Close on the measurement reading: an
apparatus probing spin alone sees the mixed $\rho_s$, and the missing purity is exactly the
which-branch information stored in position.

### 14.2 [MSc] How do I add the spin-orbit term $\propto \vec L\cdot\vec S$ and compute the hydrogen fine structure?

Assemble the pinned hydrogen $2p$ fine structure from its radial and angular factors. Radially,
define $R_{21}(r)=r\,e^{-r/2}/\sqrt{24}$ in atomic units, earn its norm, and compute
$\langle 1/r^3\rangle_{2p}$ by `Integrate` with the $r^2$ measure, expecting $\tfrac1{24}$ (verify
at authoring), cross-checked against the general
$\langle 1/r^3\rangle_{nl}=1/[n^3\,l(l+\tfrac12)(l+1)]$ at $(n,l)=(2,1)$ (verify at authoring), a
pairing that would expose a wrong radial function. Angularly, build the $6\times6$
$\vec L\cdot\vec S=\tfrac12(J^2-L^2-S^2)$ from `KroneckerProduct` of the $l=1$ matrices and
`PauliMatrix[i]/2`, and let `Eigenvalues` return $\tfrac12$ four times and $-1$ twice, matching
$\tfrac12[j(j+1)-l(l+1)-s(s+1)]$ at $j=\tfrac32,\tfrac12$ with multiplicities $2j+1$: matrix
assembly and operator identity refute each other if either is off, and the degeneracy-weighted
mean $4\cdot\tfrac12+2\cdot(-1)=0$ (the trace rule) is a further refutable identity. With
$\xi(r)=\alpha^2/(2r^3)$ the splitting comes out
$\tfrac32\cdot\tfrac{\alpha^2}{2}\langle1/r^3\rangle_{2p}=\alpha^2/32$ hartree (verify at
authoring). Close in the lab: that is about $10.9\,\mathrm{GHz}$, the measured fine-structure
doubling of the Lyman-$\alpha$ line.

### 14.3 [MSc] How do I form the total angular momentum $\vec J=\vec L+\vec S$ and change to the coupled basis?

Construct $J_i=L_i\otimes\mathbb 1_2+\mathbb 1_3\otimes S_i$ with `KroneckerProduct`, the $l=1$
matrices built from the ladder elements $\sqrt{l(l+1)-m(m\pm1)}$ via `SparseArray` and
$S_i=$ `PauliMatrix[i]/2`, and certify the algebra first: `Simplify` of $[J_x,J_y]-iJ_z$ to the
zero matrix. Then form the $6\times6$ change-of-basis matrix $U$ whose rows are
`ClebschGordan[{1, ml}, {1/2, ms}, {j, mj}]` coefficients (argument order and Condon-Shortley
phase: verify at authoring) and let exact linear algebra deliver the payoff: $U.U^\dagger=\mathbb 1$,
and $U.J^2.U^\dagger$, $U.J_z.U^\dagger$ exactly diagonal, symbolically, with entries $\tfrac{15}4$
four times and $\tfrac34$ twice beside the matching $m_j$, every off-diagonal an exact 0. The
refuting cross-check builds the same basis with independent machinery: start from the stretched
state $\vert\tfrac32,\tfrac32\rangle=\vert1,1\rangle\vert\uparrow\rangle$, ladder down with $J_-$
(normalizing each step), take the orthogonal complement within each $m_j$ block for the
$j=\tfrac12$ pair, and compare row by row against the `ClebschGordan` output; a phase-convention
slip appears as a sign mismatch. Close with the forward pointer: these coupled columns are the
weak-field eigenvectors 14.6 recovers, and their multiplet structure is what 14.2's eigenvalue
count already displayed.

### 14.4 [MSc] How do I model Stern-Gerlach as a spin-dependent spatial deflection of a spinor wave packet in a field gradient?

Evolve the pinned two-component sech packet
$\Psi(x,0)=\tfrac{1}{\sqrt2}\,\psi_0(x)\,(1,1)$ with $\psi_0$ a normalized $\operatorname{sech}$,
under $i\partial_t\Psi=[-\tfrac12\partial_x^2+\lambda x\,\sigma_z]\Psi$ with $\lambda=0.2$, the
value the C3 spinor gate probed; this class's certified spec applies verbatim: `NDSolveValue` on
the coupled pair with
`Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "DifferenceOrder" -> "Pseudospectral"}}`,
periodic boundary identification, and `AccuracyGoal -> 10, PrecisionGoal -> 10` (tight goals are
load-bearing: defaults are silently wrong at the $10^{-3}$ level and refinement at default goals is
non-monotone), with the box sized so the split packets never reach the walls over the run (centers
at $\mp\lambda t^2/2$ plus the dispersive width; walls reflect at norm 1, so norm conservation
cannot detect the contamination). The refuting anchor is the probed $\Omega=0$ limit, which here is
the model itself: each component sees a linear potential where Ehrenfest is exact, so the computed
$\langle x\rangle_\sigma(t)$ must track $\mp\lambda t^2/2$ (the probe reproduced
$\langle x\rangle_\uparrow(2)=-0.4$ to $9\times10^{-6}$ with total norm drift $5\times10^{-8}$),
with all moments and overlaps taken as $dx$-weighted grid sums, never `NIntegrate` on oscillatory
interpolants (C3 trap). The deliverable is the measurement record: the coherence
$c(t)=\int\psi_\uparrow\overline{\psi_\downarrow}\,dx$ decays as the packets separate, the reduced
spin purity falls from 1 toward $\tfrac12$ while both populations hold at $\tfrac12$. Close on
measurement as dynamics: once $\vert c\vert$ is negligible the beams are which-path records and no
downstream spin rotation revives the interference, the projection postulate emerging from unitary
evolution plus a traced-out position.

### 14.5 [MSc] How do I solve the Pauli equation for a nonrelativistic spin in an electromagnetic field?

Take uniform $\vec B=B\hat z$ in the Landau gauge $A=(0,Bx,0)$, where the Pauli Hamiltonian
$H=\tfrac12(\hat p-A)^2+\tfrac{g}{2}B S_z$ (sign set by the charge convention: verify at authoring
so the pinned coincidence lands as stated) is diagonal in $S_z$: split into the two spin blocks and
translate along $y$ with the ansatz $e^{iky}f(x)$. Per the C1 verdict this member is covered only
after the hand reduction: complete the square to a shifted oscillator of frequency $B$ centered at
$x_0=k/B$ as hand algebra, certify the reduction with a `Simplify` residual on the operator
identity, and only then hand the oscillator to the C1 route (`DSolve` with symbolic energy, the
manual termination read-off, `HermiteH` eigenfunctions with a `FullSimplify` residual of zero; the
whole-line `DEigensystem` domain form is the probed free exact confirmation). The spectrum
$E_{n,m_s}=B(n+\tfrac12)+B g\,m_s/2$ at $g=2$ exhibits the pinned coincidence: `Simplify` of
$E_{n,+1/2}-E_{n+1,-1/2}$ to 0 symbolically in $B$ and $n$, with the discriminating refutation that
`Solve` on $E_{n,+1/2}=E_{n+1,-1/2}$ for $g$ returns exactly $g=2$, so the doubled ladder splits
the moment $g$ moves: the degeneracy is a property of $g=2$, not of magnetism. The infinite
per-level Landau degeneracy was not probed and its flux counting belongs to 13.2; here only the $k$
independence of $E_{n,m_s}$, read directly off the closed form, records it. Close at the bottom of
the ladder: at $g=2$ the state $(n=0,m_s=-\tfrac12)$ sits at exactly zero energy, zero-point cost
cancelled by Zeeman gain, the nonrelativistic shadow of the Dirac value of $g$ that 22.3 derives.

### 14.6 [MSc] How do I compute the anomalous Zeeman effect and the Paschen-Back crossover, where spin-orbit and the magnetic field compete?

On the hydrogen $2p$ manifold build the exact $6\times6$
$H=\xi\,\vec L\cdot\vec S+\tfrac{B}{2}(L_z+2S_z)$ from `KroneckerProduct` with $\xi$ and $B$
symbolic (rebuild the $l=1$ and Pauli matrices in-entry; the answer stays self-contained). Certify
the block structure first, `Simplify` of $[H,J_z]$ to the zero matrix, so $H$ splits into $m_j$
blocks of sizes 1, 2, 2, 1 and `Eigenvalues` returns closed forms in $\xi,B$, the $2\times2$
blocks giving Breit-Rabi-style square roots. The two limits are the refuting checks because they
rest on independent physics: `Series` in $B$ at fixed $\xi$ must reproduce the anomalous pattern
$E_{j,m_j}\approx E_j+g_j m_j B/2$ with Lande factors $g_{3/2}=\tfrac43$, $g_{1/2}=\tfrac23$
(verify at authoring), matching the coupled-basis diagonal elements of $(L_z+2S_z)/2$; the
opposite `Series` in $\xi$ at fixed $B$ must land on the Paschen-Back pattern $B(m_l+2m_s)/2$ with
its characteristic degeneracies; and the eigenvalue sum must equal `Tr` of $H$ exactly, symbolic
in both couplings. `Plot` the six branches against $B$ at fixed $\xi$ to exhibit the crossover
near $B\sim\xi$: branches of different $m_j$ cross freely while the paired same-$m_j$ branches
avoid crossing, two-level repulsion inside each block. Close in the spectroscope: weak-field lines
splitting by $g_j$ rather than by 1 is the historical anomalous Zeeman effect, and the field where
the pattern reorganizes, Zeeman energy comparable to the $\alpha^2/32$ splitting of 14.2, is of
order one tesla for hydrogen (verify at authoring).
