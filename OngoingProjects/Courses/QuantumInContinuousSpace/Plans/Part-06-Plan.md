# Part 6 Plan: Scattering in one dimension

Questions: 6 (6.1 through 6.6). Class census: C4 for all six questions, per the `Route-Table.md`
class census row C4 and the full C4 verdict block; no other class appears in this part.

## Common ground

Every question here is a fixed-energy scattering boundary-value problem: the stationary Schrodinger
equation at a continuum energy $E$, with connection amplitudes, not eigenvalues, as the unknowns.
The probability current $j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$ defines the physical
coefficients as current ratios, $R=|r|^2$ and $T=(k'/k)\,|t|^2$, with $k=\sqrt{2E}$ the incident and
$k'$ the far-side asymptotic wavevector; $|t|^2$ alone is not $T$ whenever the two media differ. All
dynamics enters through matching: continuity of $\psi$ and $\psi'$ at every interface of a
piecewise-constant potential, replaced for a point interaction $\lambda\,\delta(x)$ by the jump
condition $\psi'(0^{+})-\psi'(0^{-})=2\lambda\,\psi(0)$ with $\psi$ continuous. The scattering
matrix collects the amplitudes, outgoing $=S\,\cdot$ incoming with
$S=\begin{pmatrix} r & t' \\ t & r' \end{pmatrix}$; current conservation is unitarity
$S^{\dagger}S=\mathbb{1}$, and the continuation $k\to i\kappa$ turns poles of $S$ into bound states.
Standing C4 inheritances (Route-Table.md, C4 verdict): answers are phrased region by region because
DSolve on a Piecewise coefficient silently echoes its input (probe 3); the primary route is symbolic
matching, Solve on the continuity system plus ComplexExpand plus FullSimplify under positivity
assumptions with wavevectors substituted last; the transfer matrix is the independent cross-check;
and the standing verification is current-based unitarity, $R+T=1$ derived from currents, never from
amplitude moduli alone.

## Per-question entries

### 6.1 [BSc] How do I compute the reflection and transmission coefficients at a potential step?

Set up the step $V=V_0\,\theta(x)$ (a heterojunction) region by region, because DSolve refuses the
piecewise potential silently: for $E>V_0$ match $\psi$ and $\psi'$ at $x=0$ with transmitted
wavevector $k'=\sqrt{2(E-V_0)}$, Solve for $r,t$, and reduce with ComplexExpand and FullSimplify
under $k,k'>0$. The point of the pinned two-regime example is that only currents survive as
probabilities: near threshold $|t|^2=4k^2/(k+k')^2\to4$ while the current-based
$T=(k'/k)|t|^2\to0$, so build $R$ and $T$ from $j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$ on
the solved amplitudes and verify $R+T=1$ there, a check any amplitude-only definition fails for
every $E>V_0$. Take the limits $E\to V_0^{+}$ ($T\to0$) and $E\gg V_0$ ($T\to1$), spot-check
numerically with the single-interface transfer matrix, then rerun the matching for $E<V_0$: the
evanescent branch $e^{-\kappa x}$ with $\kappa=\sqrt{2(V_0-E)}$ carries zero current, so $R=1$
exactly with penetration depth $1/\kappa$. Close on the contrast the ledger pins: at a step,
amplitude moduli mislead and currents do not.

### 6.2 [BSc] How do I compute the tunneling transmission $T(E)$ through a rectangular barrier?

Set up plane waves in the three regions of a barrier of height $V_0$ and width $L$, match $\psi$
and $\psi'$ at both walls, and Solve the linear system for the transmission amplitude; $T=|t|^2$
comes out in closed form via ComplexExpand and FullSimplify under $k,\kappa,L>0$ (region by
region, because DSolve refuses the piecewise potential silently), with $k=\sqrt{2E}$ and
$\kappa=\sqrt{2(V_0-E)}$ substituted last, and the result provably equal to the textbook
$[1+\tfrac{V_0^2\sinh^2\kappa L}{4E(V_0-E)}]^{-1}$ (probe 1: difference 0). Verify $R+T=1$ from
the probability currents, recover the over-barrier resonances by the continuation $\kappa\to iq$
(with $T=1$ exactly at $\sin qL=0$), and cross-check the whole thing against a transfer-matrix
product (probe 2: difference 0, $\det M=1$). Close with the opaque limit
$T\approx16\,\tfrac{E(V_0-E)}{V_0^2}\,e^{-2\kappa L}$ by Series, and the threshold edge
$E\to V_0$, where Limit gives $T\to(1+V_0L^2/2)^{-1}$.

### 6.3 [MSc] How do I exhibit transmission resonances in a well (the Ramsauer-Townsend effect)?

Take the rectangular well of depth $V_0$ and width $a$ and continue the barrier matching by
$\kappa\to ik'$ with inside wavevector $k'=\sqrt{2(E+V_0)}$ (probe 1 of the C4 verdict: $T=1$
exactly at $\sin k'a=0$). Predict the resonance ladder first, Solve/Reduce under positivity giving
$E_n=\tfrac{n^2\pi^2}{2a^2}-V_0$, with $V_0$ and $a$ placed so the $n=1$ resonance sits at low
energy (Ramsauer-Townsend); then substitute $E_n$ back into the closed-form $T(E)$, reused and
never re-typed, and FullSimplify $T-1$ to 0, with strict $T<1$ off resonance. Verify $R+T=1$ from
the currents throughout and bracket the predicted $E_1$ with transfer-matrix numeric values of $T$.
Close at the two ends, the strong-scattering dip as $E\to0$ between resonances and the deep-well
limit where the resonances sharpen, against the barrier of 6.2, whose closed form has no $T=1$
point below $V_0$.

### 6.4 [MSc] How do I assemble the transfer-matrix method for a piecewise-constant potential?

Build the interface and free-propagation $2\times2$ matrices for a symmetric double barrier
(heights $V_0$, widths $a$, gap $b$: a resonant-tunneling diode) and assemble the total $M$ with
Inverse and Dot; the transfer matrix is primary here per the C4 verdict, because the cell structure
is the physics, and the symbolic product stays tractable (probe 2: LeafCount 9775 simplifies to 186
in under a minute). Read $T=1/|M_{11}|^2$ off the product and locate the inter-barrier resonance,
invisible in either factor: on resonance $T=1$ for the symmetric structure, off resonance
$T\sim T_1^2$, far sharper than either barrier alone. Verify $\det M=1$ per factor and for the
product (flux conservation, refutable factor by factor), confirm the single-barrier factor
reproduces 6.2's closed form and that $R+T=1$ holds from the currents, and pin the numerics against
the direct 8-unknown Solve matching on the same geometry (probe 2: both give
$T=0.2965025061851854$). Close with the opaque limit, where the resonance width collapses toward a
bound-state-like line.

### 6.5 [MSc] How do I build the one-dimensional scattering matrix, verify its unitarity, and identify bound states as poles?

For the finite square well of depth $V_0$ and width $a$, run the symbolic matching for both
incidence directions (Solve, ComplexExpand, FullSimplify under $k,k'>0$, region by region as
everywhere in this part) and assemble $S=\begin{pmatrix} r & t \\ t & r \end{pmatrix}$; then
FullSimplify $S^{\dagger}S-\mathbb{1}$ to the zero matrix, which is the current-conservation
statement itself, its diagonal the standing $R+T=1$. Continue $k\to i\kappa$ and show the poles of
$S$ land on the even and odd conditions $\kappa=k'\tan(k'a/2)$ and $\kappa=-k'\cot(k'a/2)$ with
$k'=\sqrt{2(E+V_0)}$; this pole step is algebra on the probed closed form, and the C4 verdict's
open risks flag it as not itself probed, so the continuation must be shown explicitly in a cell,
never cited. The refuting check is independence: match the bound sector directly
(evanescent-outside ansatz, Solve) and locate both conditions with FindRoot at pinned parameters;
the two routes must give the same energies, and bound states enter as poles of the continued
denominator, never through an eigensolver. Close with the limits: $V_0\to0$ gives
$S\to\sigma_x=\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$, pure transmission with $r=0$ and
$t=1$, since on this convention the identity would mean $r=1$, a hard wall rather than free
space; and a pole arriving at threshold $\kappa\to0$ announces a new bound state.

### 6.6 [MSc] How do I check Levinson's theorem relating the phase shift to the number of bound states?

Give the attractive delta well $V=\lambda\,\delta(x)$, $\lambda=-g<0$ (a contact interaction),
entirely to the symbolic jump condition $\psi'(0^{+})-\psi'(0^{-})=2\lambda\,\psi(0)$ with $\psi$
continuous: nothing numeric may touch the potential, since NIntegrate returns 0. on a point measure
and `D[Abs[x],{x,2}]` yields no delta (C4 traps (c) and (d)). Solve the jump system for $r,t$,
extract the even-channel phase $\delta_{+}(k)=\arctan(g/k)$ via Arg/ArcTan, and take Limit at both
ends, $\delta\to0$ as $k\to\infty$ and the Levinson value as $k\to0$; the odd channel is exactly
unscattered ($\delta_{-}\equiv0$), the discriminating contrast the zero-range interaction gives for
free. The bound-state count comes from exact delta-well algebra, not a C4 route (the verdict's
members-sanity note): the continued $t=ik/(ik-\lambda)$ has its single pole at $k=ig$, hence
exactly one bound state at $E=-g^2/2$. Cross-check the phase as the eigenphase of $S$ (the
symmetric eigenvector $(1,1)$ gives $r+t$, the antisymmetric $(1,-1)$ gives $r-t$, which is the
free value $-1$ here, matching $\delta_{-}\equiv0$) against the jump-condition $\delta_{+}$ at a
pinned $k$, verify $R+T=1$ from the
currents, and close by cashing the theorem: the phase winding at threshold counts the spectrum.
Concern: the pinned $\delta(0)-\delta(\infty)=\pi$ is the winding of the even eigenphase
$2\delta_{+}$ of $S$; the phase shift itself obeys the half-corrected 1D Levinson
$\delta_{+}(0)-\delta_{+}(\infty)=\pi\,(n_{+}-\tfrac12)=\pi/2$ for one bound state and no
zero-energy resonance, so authoring must pin the convention before quoting $\pi$.
