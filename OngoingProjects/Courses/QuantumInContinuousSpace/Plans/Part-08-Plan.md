# Part 8 Plan: General theorems and structural methods

Six questions. Class census (per the class census table in `Route-Table.md`): C0, no differential
equation, for 8.1, 8.2, 8.3, 8.5, 8.6; C9, WKB / semiclassical asymptotics, for 8.4.

## Common ground

The part rests on identities every stationary state must obey, and on operator constructions that
generate spectra without solving a differential equation. On a normalized eigenstate $H\psi=E\psi$
with $\langle O\rangle=\int\bar\psi\,O\psi$, the virial theorem fixes
$2\langle T\rangle=\langle\vec r\cdot\nabla V\rangle$; a parameter in $H_\lambda$ moves the energy
by the Hellmann-Feynman theorem $\partial E/\partial\lambda=\langle\partial H_\lambda/\partial\lambda\rangle$;
the double commutator $[\hat x,[\hat H,\hat x]]=1$ (natural units) has the Thomas-Reiche-Kuhn sum
rule $\sum_n f_{kn}=1$, $f_{kn}=2(E_n-E_k)\vert\langle n\vert\hat x\vert k\rangle\vert^2$, as its
diagonal matrix element; at large quantum number the locally averaged stationary density approaches
the classical dwell-time density $\rho_{\mathrm{cl}}(x)=2/(T_{\mathrm{cl}}\,v(x))$, for the
oscillator the arcsine law $1/(\pi\sqrt{A^2-x^2})$; the factorization $H=A^\dagger A+E_0$ with
$A=\tfrac{1}{\sqrt2}(\tfrac{d}{dx}+W)$ and $W=-\psi_0'/\psi_0$ builds partner potentials
$V_\pm=\tfrac12(W^2\pm W')+E_0$ sharing every level except the annihilated ground state; and when
the partner reproduces its own shape, $V_+(x;a_0)=V_-(x;a_1)+R(a_1)$, the spectrum telescopes to
$E_n=\sum_{k=1}^{n}R(a_k)$ with no equation solved.

### 8.1 [BSc] How do I verify the virial theorem $2\langle T\rangle=\langle \vec r\cdot\nabla V\rangle$ on a stationary state?

Compute $\langle T\rangle$, $\langle V\rangle$, and $\langle\vec r\cdot\nabla V\rangle$ as three
independent computations on two power-law systems, so the theorem's exponent dependence is what
gets tested (C0 per `Route-Table.md`): the hydrogen $2p$ reduced radial
$u_{21}(r)=\tfrac{1}{2\sqrt6}\,r^2e^{-r/2}$ via exact Integrate,
$\langle T\rangle=\int_0^\infty(\tfrac12u'^2+\tfrac{1}{r^2}u^2)\,dr=\tfrac18$ and
$\langle V\rangle=-\tfrac14$ with $\langle T\rangle+\langle V\rangle=E_2=-\tfrac18$ as a second
identity, and the ground state of the quartic well $V=x^4/4$, which has no closed form but stays
inside C0 as a truncated harmonic number basis: ladder matrices $\hat x=(a+a^\dagger)/\sqrt2$,
$\hat p=i(a^\dagger-a)/\sqrt2$ through SparseArray and Eigensystem, the three expectations as
separate quadratic forms swept in the truncation $N$. The ratios $2\langle T\rangle=-\langle V\rangle$
and $2\langle T\rangle=4\langle V\rangle$ must then emerge, never be imposed, with an NIntegrate
spot-check on the Coulomb integrals. Close by evaluating the same combination on a 50/50
superposition of two hydrogen levels, where it fails: stationarity is load-bearing.

### 8.2 [MSc] How do I apply the Hellmann-Feynman theorem to get $\partial E/\partial\lambda$ from an expectation value?

Treat $l$ as a continuous Hellmann-Feynman parameter in the radial family
$H_l=-\tfrac12\partial_r^2+\tfrac{l(l+1)}{2r^2}-\tfrac1r$ (C0 per `Route-Table.md`), stating the
continuation honestly: at fixed node count $n_r$ the normalizable eigenfunctions exist for real
$l>-\tfrac12$ with $E=-\tfrac{1}{2(n_r+l+1)^2}$, the physical integer-$l$ spectrum is the
restriction of this analytic family, and the derivative is taken at fixed $n_r$, never fixed $n$.
D in $l$ then gives $\langle\tfrac{2l+1}{2r^2}\rangle=\tfrac{1}{(n_r+l+1)^3}$, hence
$\langle 1/r^2\rangle=\tfrac{1}{n^3(l+1/2)}$ with no integral computed; close the loop by an
independent Integrate of $\int_0^\infty u_{n_r,l}^2\,r^{-2}\,dr$ on the C1-verdict radial forms
$u_{n,l}=r^{l+1}e^{-r/n}L_{n-l-1}^{2l+1}(2r/n)$ (LaguerreL, FullSimplify under $r>0$,
$l>-\tfrac12$): fully symbolic in $l$ on the nodeless family $u_{0,l}\propto r^{l+1}e^{-r/(l+1)}$,
where both sides come out $(l+1)^{-3}$, and exact at the nodal instance $(n_r,l)=(1,1)$, a fixed
instance because the symbolic-$(n,l)$ norm is recorded unprobed in the C1 verdict. Spot-check one
case with NIntegrate and close with the large-$l$ circular-orbit limit,
$\langle 1/r^2\rangle\to n^{-4}$ (verify at authoring).

### 8.3 [MSc] How do I verify the Thomas-Reiche-Kuhn oscillator-strength sum rule?

Verify the rule where convergence is genuine content: the infinite well on $(0,L)$ from the ground
state (C0 per `Route-Table.md`), with dipole elements $x_{1n}$ from Integrate, symbolic in $n$,
nonzero only for even $n$ by parity and falling as $n^{-3}$. Build
$f_{1n}=2(E_n-E_1)\vert x_{1n}\vert^2$ reusing those elements and watch the partial sums
$S_N=\sum_{n\le N}f_{1n}$ from Sum converge to 1, led by $f_{12}=\tfrac{256}{27\pi^2}\approx0.961$,
with the tail $f_{1n}=O(n^{-4})$ and the exact closure $\lim_{N\to\infty}S_N=1$ attempted
symbolically (verify at authoring). Contrast with the oscillator, where saturation must fall out
rather than be asserted: Hermite-Gaussian Integrate gives $f_{0n}=0$ for every $n\neq1$ and
$f_{01}=1$ exactly, and from an excited state $k$ the negative emission strength $f_{k,k-1}=-k$
joins $f_{k,k+1}=k+1$ to keep the sum at 1. Close on that pair: one system spreads the strength
over an infinite even-parity ladder, the other concentrates it in a single line.

### 8.4 [MSc] How do I exhibit the correspondence principle, the classical limit of stationary states at large quantum number?

Make the limit quantitative on the oscillator $n=30$ (C9 per `Route-Table.md`): the exact
eigenfunction via HermiteH with $E_{30}=\tfrac{61}{2}$, against the classical arcsine density
$\rho_{\mathrm{cl}}(x)=\tfrac{1}{\pi\sqrt{A^2-x^2}}$ with $A=\sqrt{2E_{30}}=\sqrt{61}$, normalized
exactly on $(-A,A)$; both legs are certified machinery in the C9 verdict, whose cross-check, the
C2-certified FD recipe (SparseArray Band tridiagonal, Eigenvalues with Arnoldi shift, Richardson),
independently reproduces the $n=30$ eigenpair. The C9 verdict flags the local-averaging window as
an open authoring decision with no kernel evidence: candidates are the integral mean
$\bar\rho(x)=\tfrac{1}{2w}\int_{x-w}^{x+w}\vert\psi_{30}\vert^2$ or a grid MovingAverage, with $w$
tied to the local de Broglie wavelength $2\pi/p_{\mathrm{cl}}(x)$; record the choice and a
$w$-sensitivity check, never silently. The principle itself is the refutable object: the windowed
error $\max\vert\bar\rho-\rho_{\mathrm{cl}}\vert$ must shrink from $n=10$ to $n=30$ with the same
definitions reused at both $n$, while the raw density does not converge, its oscillation staying
comparable to the density itself. Close at the turning points, where $\rho_{\mathrm{cl}}$ diverges
while the quantum density stays finite with tunneling tails beyond: the edge where the
correspondence stops.

### 8.5 [MSc] How do I factorize a Hamiltonian as $H=A^\dagger A+E_0$ and build its supersymmetric partner potential?

Factorize the infinite well on $(0,L)$ by pure algebra (C0 per `Route-Table.md`):
$W=-\psi_1'/\psi_1=-\tfrac{\pi}{L}\cot\tfrac{\pi x}{L}$, $A=\tfrac{1}{\sqrt2}(\partial_x+W)$,
$A^\dagger=\tfrac{1}{\sqrt2}(-\partial_x+W)$. Expand $A^\dagger A$ on a generic $f(x)$ and
FullSimplify: $-\tfrac12f''-E_1f$ with $E_1=\tfrac{\pi^2}{2L^2}$ falls out rather than being
imposed, so $H_1=A^\dagger A+E_1$, and the reversed product $AA^\dagger$ yields the partner
$V_2=\tfrac{\pi^2}{L^2}\csc^2\tfrac{\pi x}{L}$, which looks nothing like a box yet must carry the
box spectrum minus its ground level. Check isospectrality against exact solvability, not against
$H_1$ itself: $V_2$ is the trigonometric Poschl-Teller well
$\tfrac{\lambda(\lambda-1)}{2}\tfrac{\pi^2}{L^2}\csc^2\tfrac{\pi x}{L}$ at $\lambda=2$ with exact
spectrum $\tfrac{\pi^2(n+2)^2}{2L^2}$ (verify at authoring), equal to the well levels for $n\ge2$;
certify the mapped eigenfunctions by the residual $H_2(A\psi_n)-E_n(A\psi_n)=0$ with
$\psi_n=\sqrt{2/L}\sin(n\pi x/L)$ at symbolic integer $n$ (fixed-$n$ instances if the kernel
stalls), with the C2-certified FD recipe on $V_2$ under Dirichlet walls as the numeric fallback.
Close on the deleted level: $A\psi_1=0$, the annihilated ground state, is exactly the state the
partner cannot have.

### 8.6 [MSc] How do I use shape invariance to obtain the oscillator and Poschl-Teller spectra purely algebraically?

Obtain both spectra from the recursion alone, no differential equation anywhere (C0 per
`Route-Table.md`): oscillator $W=x$ with constant shift $R=1$, and Poschl-Teller
$W=\lambda\tanh x$ with $V_-=\tfrac12\bigl(\lambda^2-\lambda(\lambda+1)\operatorname{sech}^2x\bigr)$
and partner parameter $\lambda\to\lambda-1$. Verify the shape-invariance identity
$V_+(x;a_0)=V_-(x;a_1)+R(a_1)$ as a zero FullSimplify residual at symbolic $x$ and $\lambda$, then
telescope $E_n=\sum_{k=1}^{n}R(a_k)$ with Sum: $E_n=n$ above the oscillator ground energy, and
$E_n=\tfrac12\bigl(\lambda^2-(\lambda-n)^2\bigr)$ on the $\lambda$ ladder. The recursion output
must match independently known exact spectra: the oscillator anchor
$\{\tfrac12,\tfrac32,\tfrac52,\tfrac72\}$ is kernel-evidenced in the C1 verdict (DEigensystem
domain form), and the $\lambda=2$ ladder must reproduce the exactly solvable pair $E_0=-2$,
$E_1=-\tfrac12$ of $V=-3\operatorname{sech}^2x$ after the asymptote shift, with the C2-certified
FD recipe on $V_-$ as the numeric fallback. Close on termination: the ladder stops when
$\lambda-n\le0$, so at $\lambda=2$ it predicts exactly two bound states, a count the algebra must
get right, not just spacings; pairing an infinite ladder with a terminating one shows the
recursion carries the bound-state count too.
