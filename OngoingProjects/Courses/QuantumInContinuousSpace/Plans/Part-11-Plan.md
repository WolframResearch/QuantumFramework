# Part 11 Plan: Orbital angular momentum in continuous space

Six questions. Class census (per the class census table in `Route-Table.md`): C0 for all of 11.1
through 11.6, differential-operator algebra and quadrature on the sphere, no solver anywhere.

## Common ground

The part rests on one operator family and its representation theory. Orbital angular momentum is
$\vec L=\vec r\times\vec p$ with $\vec p=-i\nabla$, which in spherical coordinates
($x=r\sin\theta\cos\phi$, $y=r\sin\theta\sin\phi$, $z=r\cos\theta$) closes into first-order
differential operators in the angles alone; the algebra $[L_i,L_j]=i\epsilon_{ijk}L_k$ with
Casimir $L^2=L_x^2+L_y^2+L_z^2$ and ladders $L_\pm=L_x\pm i\,L_y$; the sphere inner product
$\langle f,g\rangle=\int_0^{2\pi}\!\!\int_0^\pi\bar f\,g\,\sin\theta\,d\theta\,d\phi$ under which
the joint eigenfunctions obey $L^2Y_{lm}=l(l+1)\,Y_{lm}$ and $L_zY_{lm}=m\,Y_{lm}$; rotations act
on a multiplet by $R\,Y_{lm}=\sum_{m'}D^l_{m'm}\,Y_{lm'}$ without leaving it; and products couple
by Clebsch-Gordan coefficients $C^{LM}_{l_1m_1,l_2m_2}$. Every entry runs on D, Integrate,
Simplify and FullSimplify under the assumptions $r>0$, $0<\theta<\pi$, $\phi$ real (the C0
assumption discipline), with the built-ins SphericalHarmonicY, WignerD, and ClebschGordan serving
as anchors that constructions are checked against, never as sources the derivations quote.

### 11.1 [BSc] How do I write the orbital angular-momentum operators $L_x,L_y,L_z$ as differential operators on the angles?

Earn all three forms by explicit change of variables on the pinned generic $f(\theta,\phi)$, never
by quotation (C0 per `Route-Table.md`): form the composite $f(\theta(x,y,z),\phi(x,y,z))$ with
$\theta=\arccos\bigl(z/\sqrt{x^2+y^2+z^2}\bigr)$ and the two-argument ArcTan for $\phi$, let D
chain-rule it through $-i(x\,\partial_y-y\,\partial_x)$ and the two cyclic partners of
$\vec r\times\vec p$, then substitute $x=r\sin\theta\cos\phi$, $y=r\sin\theta\sin\phi$,
$z=r\cos\theta$ and Simplify under $r>0$, $0<\theta<\pi$, $\phi$ real, reading off
$L_z=-i\,\partial_\phi$, $L_x=i(\sin\phi\,\partial_\theta+\cot\theta\cos\phi\,\partial_\phi)$,
$L_y=i(-\cos\phi\,\partial_\theta+\cot\theta\sin\phi\,\partial_\phi)$ with every power of $r$
cancelling in the output. The read-off from the simplified generic result is where a transcription
slip hides, so verify refutably there: apply the Cartesian and the angular representation of each
component to the concrete carrier $(x+iy)^2z/r^3$ (a $Y_{32}$ shape, alive under all three
operators) and Simplify the three differences to zero, reusing both definitions. Note the
$\cot\theta$ poles at $\theta=0,\pi$ as chart artifacts of operators that are regular in Cartesian
form: the edge regime of the coordinate change. Close on what cancelled: no $L_i$ contains $r$ or
$\partial_r$, so orbital angular momentum generates pure rotations and the rest of the part lives
on the unit sphere.

### 11.2 [BSc] How do I verify $[L_i,L_j]=i\epsilon_{ijk}L_k$ and that $L^2$ is the Casimir with eigenvalue $l(l+1)$?

Redefine the three angular operators inside the answer as pure functions acting on expressions in
$(\theta,\phi)$ and run the algebra on the pinned generic $f(\theta,\phi)$ (C0 per
`Route-Table.md`): Simplify $[L_x,L_y]f-i\,L_zf$ and its two cyclic partners to zero in a single
Table over the index triples, so the whole $\epsilon_{ijk}$ structure is exercised rather than one
instance; assemble $L^2f=L_x(L_xf)+L_y(L_yf)+L_z(L_zf)$ and Simplify $[L^2,L_i]f$ to zero for all
three $i$, the Casimir property proved on a generic function rather than a lucky state. Then the
eigenvalue on the pinned concrete family: over a Table of SphericalHarmonicY spanning the full
$l=2$ multiplet $m=-2,\dots,2$, FullSimplify $L^2Y_{2m}-6\,Y_{2m}$ to zero for all five members,
so $l(l+1)=6$ is earned on the whole block. Discriminate rather than illustrate: on the
cross-multiplet mixture $Y_{11}+Y_{22}$ the same $L^2$ fails to act as a scalar (the ratio
$L^2\psi/\psi$ does not Simplify to a constant), which refutes the reading of $l(l+1)$ as a
property of the operator instead of the multiplet. Close: commutators plus Casimir spectrum are
the complete rotation fingerprint that every later construction in the part must carry.

### 11.3 [BSc] How do I build the spherical harmonics $Y_{lm}$ as the simultaneous eigenfunctions of $L^2$ and $L_z$ and check their orthonormality?

Construct the pinned $l=2$ family from its highest weight with $L_z$ and $L_\pm=L_x\pm i\,L_y$
redefined in their angular forms (C0 per `Route-Table.md`): $L_zY_{22}=2\,Y_{22}$ is a first-order
equation in $\phi$ whose DSolve solution forces $Y_{22}=g(\theta)\,e^{2i\phi}$, and $L_+Y_{22}=0$
then collapses to the one-variable ODE $g'(\theta)=2\cot\theta\,g(\theta)$, DSolve giving
$g\propto\sin^2\theta$; normalize by Integrate of ComplexExpand of $Y\bar Y$ against the sphere
measure $\sin\theta\,d\theta\,d\phi$. Ladder down with NestList of $L_-$, normalizing each rung
the same way, and confirm a fifth application annihilates the bottom state, the multiplet closing
from within. The refuting check reuses the constructed family: the $5\times5$ Gram matrix by
Outer and Integrate must equal IdentityMatrix, and any wrong rung breaks its row. Only then
compare with the built-in: each ratio of constructed $Y_{2m}$ to SphericalHarmonicY must Simplify
to a $(\theta,\phi)$-independent unimodular constant, the comparison holding up to a fixed phase
per $m$ because the ladder fixes relative phases while Condon-Shortley fixes the built-in's
absolute ones (verify at authoring which $m$ land on phase exactly $1$). Close: one first-order
PDE plus four ladder steps produce the entire multiplet; the built-in is a lookup table for what
the algebra just constructed.

### 11.4 [MSc] How do the raising and lowering operators $L_\pm$ act on the $Y_{lm}$?

Take the built-in $l=3$ family as the pinned carrier, ladder operators redefined in their angular
forms (C0 per `Route-Table.md`): over a Table in $m=-3,\dots,3$, FullSimplify
$L_+Y_{3m}-\sqrt{12-m(m+1)}\,Y_{3,m+1}$ and $L_-Y_{3m}-\sqrt{12-m(m-1)}\,Y_{3,m-1}$ to zero, the
whole multiplet at once with SphericalHarmonicY supplying every carrier, and the edge
annihilations $L_+Y_{33}=0$, $L_-Y_{3,-3}=0$ computed as the exact zeros of the factor where
$m(m\pm1)=l(l+1)$. A single-rung identity can hide a convention error, so triangulate magnitude
and phase separately, reusing the definitions: the two-step diagonal
$L_-L_+Y_{3m}=(12-m(m+1))\,Y_{3m}$ is phase-insensitive, and the norm
$\int\vert L_+Y_{3m}\vert^2\,d\Omega=12-m(m+1)$ by Integrate and ComplexExpand checks magnitude
alone; the positive real $\sqrt{l(l+1)-m(m\pm1)}$ one-step factors then hold exactly in the
Condon-Shortley phase the built-in carries, and an $m$-dependent sign in the one-step residuals is
precisely what a convention mismatch would look like. Close on termination: the factor vanishes
only at the edges, so the ladder stops itself after $2l+1=7$ states; the multiplet dimension is a
computed consequence, not an input.

### 11.5 [MSc] How do I rotate a wavefunction and represent the rotation on the $Y_{lm}$ by the Wigner $D$-matrix?

Rotate the pinned $l=2$ multiplet through its argument (C0 per `Route-Table.md`): each $Y_{2m}$ is
a homogeneous quadratic in the unit vector, so substitute the inverse rotation of $(x,y,z)$ built
from RotationMatrix in $zyz$ Euler order with $\beta$ symbolic about the $y$ axis; the flanking
$z$ rotations act as computed $\phi$ shifts, pure $e^{-im\alpha}$ and $e^{-im'\gamma}$ phases, so
the middle $\beta$ block is the honest content. Rewrite the rotated functions in angles and
project back on the basis: coefficients $c_{m'm}(\beta)$ by Integrate of $\bar Y_{2m'}$ times the
rotated $Y_{2m}$ over the sphere measure, symbolic $\beta$ surviving the quadrature. The
$5\times5$ coefficient matrix must match the WignerD entries $D^2_{m'm}$: FullSimplify each
entrywise difference to zero at symbolic $\beta$; WignerD's Euler-angle ordering, its
active-versus-passive sense, and its conjugation convention vary across references, so fix the
dictionary on one nonvanishing entry and hold it for the whole matrix (verify at authoring).
Refuting checks that reuse the machinery: the matrix is unitary at symbolic $\beta$ (ComplexExpand
under $\beta$ real), collapses to IdentityMatrix at $\beta=0$, and leaves the addition-theorem
invariant $\sum_{m}\vert Y_{2m}\vert^2=5/(4\pi)$ untouched. Close: no $l=1$ or $l=3$ component
ever appears in the expansion; the five-dimensional block is rotationally closed, which is what
irreducible means in practice.

### 11.6 [MSc] How do I couple two angular momenta and compute the Clebsch-Gordan coefficients (with the three-$Y_{lm}$ Gaunt integral as the orbital instance), the change of basis reused later for adding spin?

Compute the pinned $1\otimes1\to0,1,2$ decomposition twice on independent machinery (C0 per
`Route-Table.md`). First, ClebschGordan tabulated with Table over the physical grid only,
$M=m_1+m_2$: off-shell argument triples emit a message while returning zero, so never generate
them (the no-Quiet rule). Second, the Gaunt integrals
$\int Y_{1m_1}Y_{1m_2}\overline{Y_{LM}}\,d\Omega$ by Integrate over the sphere measure with
SphericalHarmonicY and Conjugate, compared entry by entry against the closed Gaunt form
$\sqrt{9/(4\pi(2L+1))}\,C^{L0}_{10,10}\,C^{LM}_{1m_1,1m_2}$ assembled from the same ClebschGordan
values: FullSimplify each difference to zero over all $(L,M,m_1,m_2)$. The $L=0,2$ channels pin
the coefficients up to the computed $m$-independent reduced factor; the $L=1$ channel is the edge,
both sides vanishing identically since $C^{10}_{10,10}=0$ by parity ($1+1+1$ odd), equivalently
the $L=1$ combination is antisymmetric under $m_1\leftrightarrow m_2$ while a one-sphere product
$Y_{1m_1}Y_{1m_2}$ is symmetric. Certify the odd channel, and every other, with the part's own
operators on two spheres: assemble
$\Psi_{LM}=\sum C^{LM}_{1m_1,1m_2}\,Y_{1m_1}(\theta_1,\phi_1)\,Y_{1m_2}(\theta_2,\phi_2)$ and
verify $L_z^{\mathrm{tot}}\Psi_{LM}=M\,\Psi_{LM}$ and
$(L^{\mathrm{tot}})^2\Psi_{LM}=L(L+1)\,\Psi_{LM}$ for all nine states, the total operators built
from the per-sphere angular forms; one wrong coefficient in any channel breaks it. Close in the
lab: the same Gaunt zero is the parity selection rule forbidding electric-dipole transitions
between two $p$ states, and the $9=1+3+5$ change of basis is exactly the recoupling matrix reused
when spin is added.
