# Part 23 Plan: From one particle to fields: the second-quantization bridge

Questions: 6 (23.1 through 23.6), all MSc. Class census: C0 for 23.1 through 23.5 (no differential
equation, per the `Route-Table.md` C0 row); C5 for 23.6 (time-dependent generator, per the C5
verdict block).

## Common ground

Everything in this part is one move repeated: promote the mode amplitudes of a one-particle problem
to ladder operators. The single ladder obeys $[a,a^{\dagger}]=1$ with number operator
$\hat n=a^{\dagger}a$; the multimode Fock space is the tensor product of per-mode ladders, with
occupation basis $\vert n_1,n_2,\ldots\rangle$ and operators of different modes commuting; the
field operator is the mode sum $\hat\psi(x)=\sum_n a_n\varphi_n(x)$ over an orthonormal set, with
$[\hat\psi(x),\hat\psi^{\dagger}(x')]=\delta(x-x')$ holding only in the untruncated limit; the
second-quantized Hamiltonian is
$\hat H=\sum_{ij}h_{ij}\,a_i^{\dagger}a_j+\tfrac12\sum_{ijkl}V_{ijkl}\,a_i^{\dagger}a_j^{\dagger}a_l a_k$
with $h_{ij}$ and $V_{ijkl}$ first-quantized matrix elements; the vacuum energy is the mode sum
$E_0=\tfrac12\sum_n\omega_n$; and a time-dependent quadratic Hamiltonian mixes $a$ with
$a^{\dagger}$ through a Bogoliubov transformation $b=\alpha a+\beta^{*}a^{\dagger}$ constrained by
$\vert\alpha\vert^{2}-\vert\beta\vert^{2}=1$, whose $\vert\beta\vert^{2}$ counts particles created
from vacuum. Every construction lives on an explicitly truncated space (a per-mode quantum cap, a
mode-number cap $N$, a site count, a frequency cutoff), and each answer states what its truncation
costs and where the cap is exact. Standing C0 inheritances (`Route-Table.md`, C0 row): truncated
matrices via SparseArray and KroneckerProduct with Eigenvalues/Eigensystem; delta-function matrix
elements go through Integrate only, never NIntegrate; `E` is Euler's number, energies are `en`;
Simplify needs explicit positivity, reality, and integer assumptions to close sums and residuals.

## Per-question entries

### 23.1 [MSc] How do I build multimode occupation-number (Fock) space from copies of the oscillator ladder?

Build one truncated ladder $a$ as a SparseArray with $\sqrt{1},\ldots,\sqrt{n_{\max}}$ on the
superdiagonal at $n_{\max}=4$ (five levels per mode), then embed three copies as
$a_1,a_2,a_3$ by KroneckerProduct with `IdentityMatrix[5, SparseArray]` into the $125$-dimensional
three-mode space, whose product basis is the occupation basis $\vert n_1,n_2,n_3\rangle$. Verify
the algebra where it can refute the construction: each per-mode commutator equals the truncated
identity with its corner defect, $[a_j,a_j^{\dagger}]=\mathbb 1-(n_{\max}{+}1)\,P_j$ with $P_j$
projecting mode $j$ onto its top level, exactly (the defect 4.4 exhibits for $[\hat x,\hat p]$;
demanding the ideal $\mathbb 1$ would be a false test on any finite matrix), while every cross-mode
commutator $[a_i,a_j^{\dagger}]$ with $i\ne j$ is exactly the zero matrix, defect-free. Then do
occupation arithmetic: the total number operator $\hat N=\sum_j a_j^{\dagger}a_j$ must come out
diagonal with entries $n_1+n_2+n_3$ read against the product-basis index arithmetic, and
$a_2^{\dagger}a_1$ applied to the explicit unit vector for $\vert 1,0,2\rangle$ must return exactly
the unit vector for $\vert 0,1,2\rangle$, a check a wrong tensor-slot order fails loudly. Close on
what the checks just earned: cross-mode commutativity is exact even when truncated because it is
tensor-product structure, while each mode's own ladder carries the finite-dimensional corner defect.

### 23.2 [MSc] How do I define the field operators $\hat\psi(x)$ as mode sums and impose the (anti)commutation relations?

Assemble $\hat\psi(x)=\sum_{n=1}^{N}a_n\varphi_n(x)$ on the first $N$ infinite-well modes
$\varphi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$, and impose $[a_n,a_m^{\dagger}]=\delta_{nm}$
symbolically with KroneckerDelta so that Sum reduces the operator commutator to the c-number kernel
$K_N(x,x')=\sum_{n=1}^{N}\varphi_n(x)\varphi_n(x')$: the content of the question is what this
kernel does as $N$ grows. Get the exact partial sum in closed form, a Dirichlet-type kernel: the
product-to-sum split writes $K_N$ as a difference of Dirichlet kernels in $x-x'$ and $x+x'$, and
Sum with Simplify under $0<x<L$, $0<x'<L$ should close it (verify at authoring). The refuting
check is weak convergence, never a pointwise delta claim: smear against the smooth test function
$f(x')=x'(L-x')$, whose exact Integrate against $K_N$ at the pinned interior point $x=3L/10$ must
converge to $f(3L/10)$ as $N$ doubles through $4,8,16,32$ (Fourier coefficients of $f$ decay as
$n^{-3}$, so the error should fall roughly as $N^{-2}$), while the diagonal $K_N(x,x)$ grows
without bound like $N/L$, exhibited beside it. For fermions the identical c-number kernel follows
from $\{c_n,c_m^{\dagger}\}=\delta_{nm}$, one KroneckerDelta line, so the statistics are invisible
at this level. Close on the statement the pair of computations earns: $\hat\psi(x)$ is an
operator-valued distribution, $K_N\to\delta(x-x')$ only against test functions, and the divergent
diagonal is that same fact seen pointwise.

### 23.3 [MSc] How do I represent a free scalar field as infinitely many oscillators and read off the vacuum energy?

Massless scalar field in a box of length $L$ with Dirichlet walls: mode frequencies
$\omega_n=n\pi/L$, vacuum energy $E_0=\tfrac12\sum_n\omega_n$. Exhibit the divergence with an
explicit sharp cutoff first: Sum gives $E_0(n_{\max})=\tfrac{\pi}{4L}n_{\max}(n_{\max}{+}1)$,
growing as the square of the cutoff. Then regulate smoothly with $e^{-\omega_n/\Lambda}$: Sum
closes geometrically to $\tfrac{\pi}{2L}\,e^{\pi/(L\Lambda)}/(e^{\pi/(L\Lambda)}-1)^{2}$ (verify at
authoring), and Series in $1/\Lambda$ splits it as $\tfrac{L\Lambda^{2}}{2\pi}-\tfrac{\pi}{24L}
+O(\Lambda^{-2})$: a bulk divergence proportional to $L$ plus a finite geometry-dependent term. The
Casimir-flavored difference is two geometries at one cutoff: partition a box of fixed total length
$D$ at position $L$ and form $\Delta E(L)=E_0(L)+E_0(D{-}L)-2E_0(D/2)$, reusing the same regulated
summand; the bulk terms cancel identically and Limit as $\Lambda\to\infty$ leaves
$-\tfrac{\pi}{24}\left(\tfrac1L+\tfrac1{D-L}-\tfrac4D\right)$. State plainly what is and is not
claimed: this is mode-sum regularization of a one-dimensional scalar toy, not the electromagnetic
Casimir experiment; the claim is only that the difference of two divergent sums has a finite,
cutoff-independent limit. The refuting check is exactly that independence: sweep $\Lambda$ over a
decade numerically with the same summand and watch $\Delta E$ settle on the Series value; residual
$\Lambda$ dependence refutes the regularization. Close with the force: $-\partial_L\Delta E$ pulls
the partition toward the nearer wall, approaching the attraction $-\pi/(24L^{2})$ as $D\to\infty$.

### 23.4 [MSc] How do I show the second-quantized many-body Hamiltonian equals the symmetrized first-quantized one?

Take the first $M=3$ well modes on $(0,1)$, one-body part diagonal with $E_i=i^{2}\pi^{2}/2$, and
the contact pair interaction $g\,\delta(x_1-x_2)$, whose tensor
$V_{ijkl}=g\int_0^1\varphi_i\varphi_j\varphi_k\varphi_l\,dx$ is exact rationals times $g$ by
Integrate only (NIntegrate returns 0. on a point measure, the C0 binding fact). First-quantized
side: the six symmetrized pair states
$\vert ij\rangle_S=(\vert ij\rangle+\vert ji\rangle)/\sqrt{2(1+\delta_{ij})}$ and the $6\times6$
matrix of $h\otimes\mathbb 1+\mathbb 1\otimes h+g\,\delta(x_1-x_2)$ by exact Integrate, the delta
collapsing the double integral to one dimension. Second-quantized side: rebuild per-mode ladders
capped at two quanta (a $27$-dimensional space) via SparseArray and KroneckerProduct, assemble
$\hat H=\sum_i E_i\,a_i^{\dagger}a_i+\tfrac12\sum_{ijkl}V_{ijkl}\,a_i^{\dagger}a_j^{\dagger}a_l a_k$,
and restrict to the total-number-two eigenspace of $\hat N$ in the normalized pair basis
$a_i^{\dagger}a_j^{\dagger}\vert 0\rangle/\sqrt{1+\delta_{ij}}$; the two-quanta cap represents the
two-boson sector exactly, since normal-ordered $\hat H$ never crosses the truncation edge inside
the sector, so equality can be demanded exactly rather than approximately. The interacting term is
the content: the free parts agree trivially, so the refuting check is Simplify of the entrywise
difference of the two $6\times6$ matrices to the zero matrix with $g$ symbolic, where a wrong
$\tfrac12$, a dropped exchange term, or an index-order slip in $V_{ijkl}$ shifts the pair-changing
off-diagonal entries. Close by reading one physical number off the shared matrix: the
$\vert 11\rangle$ diagonal entry $\pi^{2}+\tfrac{3g}{2}$ (verify at authoring), the first-order
shift of two bosons in one orbital, contrasted with the pair-changing entries a mean-field
treatment discards.

### 23.5 [MSc] How do I put a bosonic field on a finite lattice (truncated) and compute its dispersion relation?

Harmonic chain of $N=24$ sites, unit masses and springs, periodic boundary conditions: the
dynamical matrix is the circulant SparseArray with $2$ on the diagonal, $-1$ on the neighbors, and
$-1$ in the corners. Take the exact route first: substitute the Fourier mode $u_j=e^{ikj}$ with
$k=2\pi m/N$ into the eigenvalue equation and Simplify to $2-2\cos k=4\sin^{2}(k/2)$, so
$\omega(k)=2\vert\sin(k/2)\vert$ exactly, doubly degenerate in $\pm m$. The refuting check reuses
the matrix, never a re-typed formula: Sqrt of the numeric Eigenvalues of the SparseArray, sorted,
against the sorted `Table` of $2\vert\sin(\pi m/N)\vert$ over $m=0,\ldots,N{-}1$, agreeing to
machine zero; an open-chain slip or sign error yields the shifted sine family and the wrong
degeneracy pattern, failing loudly. Then quantize: in normal coordinates
$\hat H=\sum_k\omega(k)\,(a_k^{\dagger}a_k+\tfrac12)$, the lattice field is $N$ independent
ladders, the multimode structure of 23.1 realized by a mechanical system. Quantify the continuum
limit with Series: $2\sin(k/2)=k-k^{3}/24+O(k^{5})$, the linear phonon branch $\omega\to\vert
k\vert$ with unit sound speed, and at the zone edge $\omega(\pi)=2$ against the linear
extrapolation $\pi$, a deficit of $1-2/\pi\approx36\%$, with the group velocity $\cos(k/2)$ (by D)
vanishing there. Close on the physics: the chain emulates a massless continuum field only at
wavelengths long against the spacing, and the zone-edge mode is a standing wave that transports
nothing.

### 23.6 [MSc] How does a time-dependent boundary or parameter create particles from the vacuum (a dynamical-Casimir-type capstone, with the truncation made explicit)?

C5 member (`Route-Table.md`, C5 verdict): the time-dependent generator goes through NDSolveValue at
tight goals, norm drift grows with integration length, and parametric-resonance stiffness was
probed only by proxy (Landau-Zener), flagged. Single mode with modulated frequency
$\omega(t)=\omega_0(1+\epsilon\sin 2\omega_0 t)$, $\omega_0=1$, $\epsilon=1/20$: the mode function
obeys $\ddot f+\omega(t)^{2}f=0$ with vacuum data $f(0)=1/\sqrt{2\omega_0}$,
$\dot f(0)=-i\sqrt{\omega_0/2}$; integrate the equivalent first-order complex system with
NDSolveValue at PrecisionGoal and AccuracyGoal 10 or higher (WorkingPrecision above machine only
with exact rationalized input, C5 trap). The truncation is explicit and exact here: the field is
reduced to one mode, legitimate because a parametric frequency drive does not couple modes, whereas
a genuinely moving boundary would mix them and demand the multimode machinery of 23.1. The
adiabatic particle number is phase-free,
$n(t)=\bigl(\vert\dot f\vert^{2}+\omega(t)^{2}\vert f\vert^{2}\bigr)/\bigl(2\omega(t)\bigr)
-\tfrac12=\vert\beta(t)\vert^{2}$, and the refuting invariant is the Wronskian
$f\dot f^{*}-f^{*}\dot f=i$, equivalently $\vert\alpha\vert^{2}-\vert\beta\vert^{2}=1$, monitored
along the whole trajectory: drift beyond the goal tolerance refutes the run before any physics is
read. On the resonant drive the number grows exponentially, $\vert\beta\vert\sim\sinh(\lambda t)$
with parametric rate $\lambda=\epsilon\omega_0/2$ for this modulation depth (verify at authoring):
fit the late-time slope of $\log n(t)$ out to $t=200$ with LinearModelFit and compare rates in the
asymptotic window only, never pointwise at early times (the C5 wrong-benchmark trap at finite
time); because the stiffness was probed only by proxy, a goals sweep at authoring must show the
fitted rate solver-independent. The contrast is the detuned drive
$\omega(t)=\omega_0(1+\epsilon\sin(2.3\,\omega_0 t))$, outside the resonance tongue of width of
order $\epsilon\omega_0$: $n(t)$ stays bounded and oscillatory. Close on the capstone reading:
modulating a cavity parameter at $2\omega_0$ converts vacuum into pairs, the dynamical Casimir
mechanism, and it is resonance, not modulation as such, that creates particles secularly.
