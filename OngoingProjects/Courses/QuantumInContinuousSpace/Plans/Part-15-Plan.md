# Part 15 Plan: Identical particles in continuous space

Five questions, all MSc. Class census (`Route-Table.md`): C0 for 15.1, 15.2, 15.3, 15.4 (no
differential equation, symbolic `Integrate` and matrix algebra); C6 for 15.5 (mean-field SCF).

## Common ground

Two identical particles admit only states of definite exchange symmetry,
$\psi(x_2,x_1)=\pm\psi(x_1,x_2)$: from distinct orbitals $\phi_a,\phi_b$ the normalized pair states
are $\psi_{\pm}=(\phi_a(x_1)\phi_b(x_2)\pm\phi_b(x_1)\phi_a(x_2))/\sqrt{2}$, and for $N$ fermions the
antisymmetrizer is the Slater determinant $\Psi=\det[\phi_i(x_j)]/\sqrt{N!}$. Statistics becomes
observable through the two-body density $\rho_2(x_1,x_2)=|\psi(x_1,x_2)|^2$ at coincidence, and
through the energy via the direct and exchange integrals
$J=\iint|\phi_a(\mathbf r_1)|^2\,v_{12}\,|\phi_b(\mathbf r_2)|^2$ and
$K=\iint\phi_a^{*}(\mathbf r_1)\phi_b(\mathbf r_1)\,v_{12}\,\phi_b^{*}(\mathbf r_2)\phi_a(\mathbf r_2)$,
which split spatially symmetric from antisymmetric pairs as $J\pm K$; when the pair interaction is too
hard to diagonalize exactly, a mean field closes it self-consistently, $F[\phi]\,\phi=\varepsilon\phi$.

### 15.1 [MSc] How do I build a two-particle wavefunction $\psi(x_1,x_2)$ and symmetrize or antisymmetrize it?

Work in the infinite well of width $L$ with orbitals $\chi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$, $n=1,2$,
distinct so exchange has something to act on. Define the product $\chi_1(x_1)\chi_2(x_2)$ and the pair
$\psi_{\pm}=(\chi_1(x_1)\chi_2(x_2)\pm\chi_2(x_1)\chi_1(x_2))/\sqrt{2}$; earn the normalization with
`Integrate` of $|\psi_{\pm}|^2$ over $[0,L]^2$ under $L>0$ (exactly 1, riding on orbital
orthonormality), and earn the symmetry by applying the swap rule `{x1 -> x2, x2 -> x1}` and
confirming with `Simplify` that $\psi_{+}$ is invariant and $\psi_{-}$ flips sign. The fermion node
line is `Simplify` of $\psi_{-}$ at $x_2\to x_1$ returning 0 identically; the boson coincidence
density $|\psi_{+}(x,x)|^2=2\,\chi_1(x)^2\chi_2(x)^2$, exactly twice the distinguishable product's
coincidence value, quantifies the enhancement. All three checks reuse the defined $\psi_{\pm}$ and can
refute: a wrong prefactor fails the norm, a dropped sign fails the swap or the node. Close on the
equal-orbital edge: with both orbitals $\chi_1$ the antisymmetric combination vanishes identically
(exclusion falls out of the algebra), while the naive symmetric $1/\sqrt{2}$ state has norm
$\sqrt{2}$, so the general prefactor is $1/\sqrt{2(1+\delta_{ab})}$.

### 15.2 [MSc] How do I exhibit the exchange hole for fermions and the bunching for bosons (Pauli exclusion versus enhancement)?

Use the oscillator pair $\phi_0(x)=\pi^{-1/4}e^{-x^2/2}$ and
$\phi_1(x)=\sqrt{2}\,\pi^{-1/4}\,x\,e^{-x^2/2}$ (built from `HermiteH`), the same orbitals for both
statistics so only the sign differs. Form $\rho_2^{\pm}=|\psi_{\pm}(x_1,x_2)|^2$ and reduce with
`Simplify` under `Element[{x1, x2}, Reals]`: the fermion density collapses to the closed form
$\rho_2^{-}=(x_1-x_2)^2\,e^{-x_1^2-x_2^2}/\pi$, whose quadratic zero along $x_1=x_2$ is the exchange
hole (exhibit the $s^2$ opening with `Series` in $s=x_1-x_2$), while the boson coincidence value
$\rho_2^{+}(x,x)=2\,\phi_0(x)^2\phi_1(x)^2$ doubles the distinguishable reference
$\rho_2^{\mathrm{d}}=\tfrac12\left(|\phi_0(x_1)\phi_1(x_2)|^2+|\phi_1(x_1)\phi_0(x_2)|^2\right)$:
factor 0 against factor 2, both closed forms. For the pair correlation compute the marginal
$\rho(x)=\int\rho_2^{\pm}\,dx_2$ by `Integrate` and form
$g(x_1,x_2)=\rho_2^{\pm}/(\rho(x_1)\rho(x_2))$; the refuting check is that the marginal must equal
$\tfrac12(\phi_0(x)^2+\phi_1(x)^2)$ for both signs (orthogonality kills the cross term), so any
normalization slip in $\psi_{\pm}$ surfaces immediately. Close in the lab: coincidence counters on
cold-atom pairs read exactly these two numbers, a doubled accidental rate for bosons and an
antibunching null for fermions, the atomic Hanbury Brown-Twiss experiment.

### 15.3 [MSc] How do I build a Slater determinant for $N$ fermions in a well or oscillator?

Fill the well orbitals $\chi_n(x)=\sqrt{2/L}\,\sin(n\pi x/L)$, $n=1,2,3$, and build
$\Psi(x_1,x_2,x_3)=\det[\chi_n(x_j)]/\sqrt{3!}$ with `Det` on the `Table`-built $3\times3$ orbital
matrix, a genuine $N=3$ object rather than the two-body shortcut. Verify antisymmetry under one
explicit transposition, `Simplify` of $(\Psi$ at `{x1 -> x2, x2 -> x1}`$)+\Psi$ to 0, and the
coincidence node, `Simplify` of $\Psi$ at $x_2\to x_1$ to 0. The load-bearing check is the
determinant identity for orthonormal orbitals: the one-body density
$\rho(x)=3\iint|\Psi(x,x_2,x_3)|^2\,dx_2\,dx_3$ must equal $\sum_{n=1}^{3}\chi_n(x)^2$; compute the
left side by nested `Integrate` under $0<x<L$ and `Simplify` the difference against
`Sum[chi[n, x]^2, {n, 3}]` to 0, which refutes a wrong $1/\sqrt{3!}$ or a non-orthonormal orbital
choice while stating the physics: filled levels add in density with no interference. Close by
occupying $\chi_1$ twice in the matrix: `Det` with two equal rows returns 0 identically, Pauli
exclusion as a determinant fact rather than a postulate.

### 15.4 [MSc] How do I compute the exchange energy and the para/ortho splitting of a two-electron (helium-like) model?

Use hydrogenic orbitals with symbolic nuclear charge $Z$: $\phi_{1s}=Z^{3/2}e^{-Zr}/\sqrt{\pi}$ and
$\phi_{2s}=Z^{3/2}(2-Zr)\,e^{-Zr/2}/(4\sqrt{2\pi})$, orthonormality earned first by `Integrate` with
$Z>0$. Route both Coulomb integrals through the multipole expansion
$1/r_{12}=\sum_{l}r_{<}^{l}\,r_{>}^{-l-1}\,P_l(\cos\theta_{12})$: the angular integrals of s orbitals
kill every $l>0$ term, so $J$ and $K$ collapse to nested radial quadratures with kernel $1/r_{>}$,
done as two explicit `Integrate` pieces (inner $\int_0^{r_1}dr_2$ weighted by $1/r_1$ plus
$\int_{r_1}^{\infty}dr_2$ weighted by $1/r_2$, then the outer $r_1$ integral, assumptions $Z>0$,
$r_1>0$); expected closed forms $J=\tfrac{17}{81}Z$ and $K=\tfrac{16}{729}Z$ (verify at authoring),
splitting $2K=\tfrac{32}{729}Z$. Cross-check with `NIntegrate` on the same nested radial integrand at
$Z=2$, reusing the orbital definitions; the physics check is the sign chain $J>K>0$ ($K$ is the
Coulomb self-energy of the overlap density $\phi_{1s}\phi_{2s}$, hence positive), which places the
ortho triplet at $J-K$ below the para singlet at $J+K$ from a spin-independent Hamiltonian. Close
with the number: at $Z=2$, $2K\approx0.088$ Hartree $\approx2.4$ eV, overshooting the measured helium
$1s2s$ splitting near $0.8$ eV because the unscreened $2s$ orbital sits too deep; screening is
exactly what a self-consistent field supplies.

### 15.5 [MSc] How do I solve the Hartree-Fock self-consistent-field equations for a model atom numerically?

This is the C6 SCF member (`Route-Table.md`, C6 verdict, probe p6): softened-Coulomb one-dimensional
helium on a grid, one-body Hamiltonian by finite differences via `SparseArray` with `Band`, the
electron-electron kernel a dense matrix on the $dx$-weighted grid; the softening parameters and the
closed-shell reduction are re-derived in the part, showing that for the doubly occupied orbital the
exchange term equals the direct self-term, so RHF collapses to the Hartree operator
$F[\phi]=h+J[\phi]$. Iterate with `NestList`, dense `Eigensystem` per step, linear mixing: the probed
pattern converges geometrically with contraction $\approx0.135$ at plain mixing, but convergence is
declared only on the orbital change norm $\|\Delta\phi\|$ (`Norm`), never from an eigenvalue plateau
(the C6 refutation: a plateau can mask a limit cycle); the converged orbital must pass the stationary
residual $\|(h+J[\phi])\phi-\varepsilon\phi\|$, two mixings $\alpha=1$ and $\alpha=0.5$ must land on
one fixed point, and the probe anchors $\varepsilon\approx-0.75032$,
$E_{\mathrm{HF}}=2\varepsilon-U\approx-2.22450$ with a parity-symmetric orbital (energies named `en`,
never `E`). Grade against the exact benchmark: the two-electron grid Hamiltonian
$h\otimes\mathbb{1}+\mathbb{1}\otimes h+V_{\mathrm{ee}}$ assembled with `KroneckerProduct` on
`SparseArray`, ground state by the C2-certified recipe (`Eigenvalues` with Arnoldi and `"Shift"`
below the ground energy, since unshifted magnitude ordering silently misses the negative spectrum)
plus a grid-doubling sweep. The honest close is the correlation energy
$E_{\mathrm{corr}}=E_{\mathrm{exact}}-E_{\mathrm{HF}}<0$, the part of the physics no mean field can
reach, pinned to a number by the exact grid.
