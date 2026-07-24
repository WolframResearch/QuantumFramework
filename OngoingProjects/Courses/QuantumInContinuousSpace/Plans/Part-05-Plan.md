# Part 5 Plan: The time-dependent Schrodinger equation as a PDE

10 questions. Class census per the class-census table in `Route-Table.md`: C3 (time-dependent
linear Schrodinger PDE): 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.9, 5.10; C0 (no differential equation):
5.1, 5.8.

## Common ground

Everything in this part rests on the time-dependent Schrodinger equation
$i\,\partial_t\psi = -\tfrac12\,\partial_{xx}\psi + V(x)\,\psi$ and four of its consequences: a
stationary state evolves only by the phase $e^{-iEt}$; a general state propagates in the energy
eigenbasis as $\psi(x,t)=\sum_n c_n e^{-iE_n t}\psi_n(x)$ with $c_n=\int\bar\psi_n\psi_0\,dx$;
the same evolution is carried by the propagator, $\psi(x,t)=\int K(x,t;x')\,\psi(x',0)\,dx'$;
the norm $\int|\psi|^2\,dx$ is conserved through the continuity equation
$\partial_t\rho+\partial_x j=0$ with $j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$; and
expectation values obey Ehrenfest's relations $\frac{d\langle x\rangle}{dt}=\langle p\rangle$ and
$\frac{d\langle p\rangle}{dt}=-\langle V'(x)\rangle$, which reduce to classical mechanics only
when $\langle V'(x)\rangle=V'(\langle x\rangle)$.

The eight C3 entries share the certified contract of the C3 verdict in `Route-Table.md`: primary
`NDSolveValue` with `Method -> {"MethodOfLines", "SpatialDiscretization" ->
{"TensorProductGrid", "DifferenceOrder" -> "Pseudospectral"}}`, periodic identification
$\psi(t,L_1)=\psi(t,L_2)$, and `AccuracyGoal -> 10, PrecisionGoal -> 10` (tight goals are
load-bearing: defaults are silently wrong at the $0.25\%$ level with no message, and grid
refinement at default goals is non-monotone); independent cross-check by a hand-built Strang
split-step Fourier propagator on `Fourier`/`InverseFourier` (`FourierParameters -> {1,-1}`,
packed arrays); a box sized so the exact tail at the edge is negligible, edge
$\gtrsim k_0T+8\sigma_T$, because norm conservation does not detect wall reflection; and moments
or probabilities taken as $dx$-weighted grid sums, never `NIntegrate` on an interpolant.

### 5.1 [BSc] How do I show that a stationary state evolves only by the phase $e^{-iEt}$, while a superposition acquires genuine time dependence?

Work in the Poschl-Teller well $V=-3\operatorname{sech}^2x$, whose two bound states
$\psi_0\propto\operatorname{sech}^2x$ at $E_0=-2$ and $\psi_1\propto\operatorname{sech}x\tanh x$
at $E_1=-\tfrac12$ are closed forms in a non-box system (the infinite well is reserved for 5.5).
This is C0 machinery per the `Route-Table.md` C0 row: earn both levels first with `D` +
`FullSimplify` residuals $\hat H\psi_n-E_n\psi_n\to0$, then form $e^{-iE_0t}\psi_0$ and the
equal-weight superposition and reduce both densities with `ComplexExpand` under
`Element[{x, t}, Reals]`: the stationary density is $t$-free identically, while the superposition
beats as $\cos(\tfrac32 t)$ at $\Delta E=\tfrac32$. Refute or confirm with a numeric spot value
of the density at one $(x,t)$ against the period $4\pi/3$. Close on the edge that defines
"genuine": a superposition of degenerate levels would still have a static density, so time
dependence lives in energy differences, not in superposition itself.

### 5.2 [BSc] How do I integrate the time-dependent Schrodinger equation directly as a PDE with `NDSolve` for a given initial $\psi(x,0)$?

Release the cusp packet $\psi_0\propto e^{-|x-x_0|}$, $x_0=3$, in the harmonic trap $V=x^2/2$:
its momentum amplitude falls only as $(1+k^2)^{-1}$, a genuine grid stress a smooth default
packet would never apply. Integrate under the shared C3 contract (running a default-goals twin
alongside exhibits the silent-error trap live), with the box beyond the classical turning point
of the highest resolved momentum and the split-step propagator as cross-check. The trap makes
the run refutable through exact anchors that reuse $\psi_0$: the density at $t=\pi$ must equal
$\rho_0(-x)$ (mirror), the recurrence overlap $|\langle\psi_0|\psi(2\pi)\rangle|$ must return to
1, and the grid-sum drift of $\langle H\rangle$ checks independently. Close with the convergence
law under grid doubling, tail-limited by the cusp, against the spectral convergence a smooth
packet would show.

### 5.3 [BSc] How do I follow a free Gaussian wave packet's spreading and group velocity, both analytically and by `NDSolve`?

Take the question's own carrier, the free Gaussian with symbolic width $\sigma_0$ and drift
$k_0$ (named-trivial policy), and get the closed form analytically first: `DSolve` on the free
IVP, cross-derived by the `FourierTransform` route, both certified in the C3 verdict; this
closed form is the class's permanent benchmark. Then run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$ and compare pointwise, plus the moment-law residuals
$\langle x\rangle(t)-k_0t$ and $\sigma^2(t)-[\sigma_0^2+t^2/(4\sigma_0^2)]$ from grid sums, and
exhibit the norm-blind boundary trap deliberately: an undersized box reflects at norm 1. The
discriminating contrast is a same-width sech packet, whose rescaled density fails to collapse
onto its initial shape: self-similar spreading is a Gaussian privilege, not a free-particle law.
Close with the ballistic limit $\sigma(t)\to t/(2\sigma_0)$.

### 5.4 [MSc] How do I implement the split-step Fourier propagator from scratch and benchmark it against `NDSolve`?

Build the Strang map $e^{-iV\,dt/2}e^{-iK\,dt}e^{-iV\,dt/2}$ from scratch on
`Fourier`/`InverseFourier` (`FourierParameters -> {1,-1}`, packed arrays, the $k$-grid
convention checked before use) and benchmark it against the shared `NDSolveValue` contract on a
packet oscillating in the symmetric quartic double well $V=(x^2-a^2)^2/2$: a
local-harmonic-width Gaussian displaced off the left minimum, $a$ fixed at authoring so
inter-well leakage is visible on the run time. This question is the C3 verdict's own certified
cross-check pair, and the tunneling tail is where phase errors make the two integrators
measurably disagree. Verify that the pointwise disagreement shrinks under $dt$ halving and goal
tightening, and measure the Strang order in a coarse-$dt$ window, since the verdict flags the
order ratio as contaminated at the $10^{-13}$ floor; the discriminator is the per-step norm,
split-step holding $|\,\lVert\psi\rVert^2-1|$ near $10^{-13}$ while `NDSolveValue` drifts near
$10^{-5}$. Close with the $V\to0$ limit, which must reproduce the free-Gaussian closed form.

### 5.5 [MSc] How do I evolve a state by expanding it in energy eigenstates and propagating term by term, and observe wave-packet revivals?

Expand the off-center triangular packet peaked at $x_0=L/3$ in the infinite well of width $L$:
the kink populates many modes with $c_n$ decaying only as $1/n^2$, and the quadratic spectrum
$E_n=n^2\pi^2/(2L^2)$ makes the revival exact, $E_nT_{\mathrm{rev}}=2\pi n^2$ at
$T_{\mathrm{rev}}=4L^2/\pi$. The primary is the exact eigenbasis per the `Route-Table.md`
cross-class hand-off (5.5 leans on the box's exact spectrum): $c_n$ by `Integrate` in closed
form with `Assumptions -> Element[n, Integers]`, then a truncated `Sum` bounded by the tail mass
$\sum_{n>N}|c_n|^2$. The revival is exact by construction there, so the refutable content sits
in the mirror identity $\psi(x,T_{\mathrm{rev}}/2)=-\psi_0(L-x)$, compared pointwise between the
truncated series and the reflected initial function, and in an independent `NDSolveValue`
evolution whose overlap $|\langle\psi_0|\psi(T_{\mathrm{rev}})\rangle|$ must return to 1; that
run takes explicit Dirichlet conditions $\psi(t,0)=\psi(t,L)=0$ at tight goals rather than the
periodic-pseudospectral primary, because the walls are physical and the C3 verdict certifies
only the evolution cross-check for this member. As a long-time member it inherits the flagged
risks: default-goal drift grows with integration length and pseudospectral cost above $n=513$ is
unprobed, so the numeric leg carries a goal-and-resolution sweep. Close with the quarter-revival
two-image cat at $T_{\mathrm{rev}}/4$.

### 5.6 [MSc] How does an oscillator coherent state evolve, so that $\langle x\rangle(t)$ traces the classical oscillation while the packet stays minimal and does not spread?

Evolve the coherent state $\alpha=2$ ($x_0=2\sqrt2$, $p_0=0$, width
$\sigma_{\mathrm{gs}}^2=\tfrac12$) beside a squeezed-displaced Gaussian with the same $(x_0,p_0)$
but $\sigma_0^2=\tfrac18$ ($r=\ln 2$) in the trap $V=x^2/2$: the same-centroid partner isolates
what "coherent" adds to "displaced", turning "does not spread" into a statement about the width
alone. Get the exact side first: for quadratic $V$ the moment system
$\langle x\rangle,\langle p\rangle,\langle x^2\rangle,\langle xp{+}px\rangle,\langle p^2\rangle$
closes and `DSolve` solves it, giving $\langle x\rangle(t)=2\sqrt2\cos t$ for both states and
$\sigma^2(t)=\sigma_0^2\cos^2t+\tfrac{1}{4\sigma_0^2}\sin^2t$. Run both under the shared C3
contract with moments as $dx$-weighted grid sums, and refute or confirm on: the coherent width
holding $\tfrac12$ for all $t$ while the squeezed width breathes between $\tfrac18$ and $2$ at
frequency 2, the two centroids agreeing pointwise, and the exact recurrence
$|\langle\psi_0|\psi(2\pi)\rangle|=1$. Close with the $r\to0$ limit collapsing the pair onto one
state, noting the edge that a strong squeeze under-resolves the narrow phase of the breath.

### 5.7 [MSc] How do I watch a wave packet scatter off a barrier in real time and see tunneling and reflection split the packet?

Boost the compact-support raised-cosine packet
$\psi_0\propto\left(1+\cos\frac{\pi(x-x_c)}{w}\right)e^{ik_0x}$ on $|x-x_c|<w$ at the Eckart
barrier $V_0\operatorname{sech}^2(x/a)$: the initial barrier overlap is exactly zero, so the
reflected and transmitted fragments are unambiguous, and the barrier carries the exact anchor
$T(E)=\sinh^2(\pi ka)\,/\,[\sinh^2(\pi ka)+\cosh^2(\tfrac\pi2\sqrt{8V_0a^2-1})]$ with
$k=\sqrt{2E}$, valid for $8V_0a^2>1$. Run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$ on both sides (both fragments travel) and the split-step cross-check.
The refuting check: the late-time transmitted mass (grid sum beyond the barrier) against the
momentum-resolved prediction $\int T(k)\,|\tilde\psi_0(k)|^2\,dk$ built from the same $\psi_0$
via `FourierTransform`, with $R+T=1$ from the same sums. The C3 verdict flags 5.7's ringing as
its unmeasured open risk: the smooth Eckart barrier retires the potential-side Gibbs ringing
while keeping an exact anchor, but the packet's $C^1$ support seams still seed high-$k$ content,
so a resolution sweep of the seam ringing is mandatory at authoring. Close with $T\to1$ at
$E\gg V_0$ set against the deep-tunneling suppression.

### 5.8 [MSc] How do I construct the propagator (kernel) $K(x,t;x',0)$ for the free particle and for the oscillator (the Mehler kernel)?

Construct the free kernel $K_0=(2\pi it)^{-1/2}e^{i(x-x')^2/(2t)}$ from the momentum integral
$\frac{1}{2\pi}\int e^{ik(x-x')-ik^2t/2}\,dk$ by `Integrate` with a $t\to t-i\epsilon$ regulator
and `Limit`, and the Mehler kernel
$K_{\mathrm{osc}}=(2\pi i\sin t)^{-1/2}\exp\!\left\{\tfrac{i}{2\sin t}\left[(x^2+x'^2)\cos t-2xx'\right]\right\}$
by summing the Wick-rotated eigen-sum $\sum_n e^{-(n+1/2)\tau}\psi_n(x)\psi_n(x')$ (`HermiteH`)
into its closed Gaussian form and continuing $\tau\to it$; at real $t$ the eigen-sum oscillates
without converging pointwise, so the Wick-rotated comparison is the honest second
representation (C0 per the `Route-Table.md` C0 row). Each kernel earns its name three ways, any
failure refuting it: the TDSE residual in $(x,t)$ is identically $0$ under `D` + `FullSimplify`;
the composition law $\int K(x,t_1;x'')\,K(x'',t_2;x')\,dx''=K(x,t_1{+}t_2;x')$ closes; and the
$t\to0^+$ action on a test packet returns the packet, the $\delta(x-x')$ limit taken through a
test function, never as a bare limit. Then let the free kernel act on the Berry-Balazs Airy beam
(`AiryAi`, transformation-object reuse allowed by the ledger) and land on the accelerating form
$\operatorname{Ai}(x-t^2/4)\,e^{it(x-t^2/6)/2}$. Close at the caustics $t=n\pi$, where the
Mehler prefactor diverges and the kernel degenerates to mirrored delta transport, with the
$\sin t\to t$ free-kernel limit as the bridge between the two.

### 5.9 [BSc] How do I verify Ehrenfest's theorem numerically, that $\langle x\rangle$ and $\langle p\rangle$ obey the classical equations of motion?

Evolve the asymmetric two-lobe packet, the normalized unequal-weight sum
$e^{-(x-1.2)^2}+\tfrac12 e^{-(x+0.6)^2}$, in the quartic well $V=x^4/4$ (the ledger pins
"packet in the quartic well"; this skewed specification sits inside that class): its skewness
makes $\langle x^3\rangle-\langle x\rangle^3$ order one from $t=0$, so Ehrenfest's actual
content, $\frac{d\langle p\rangle}{dt}=-\langle V'(x)\rangle=-\langle x^3\rangle$ against the
classical $-V'(\langle x\rangle)=-\langle x\rangle^3$, is the largest feature on the plot rather
than a subtlety. Run the shared C3 contract with the box beyond the packet's classical turning
points, moments as $dx$-weighted grid sums and $\langle p\rangle$ from the spectral derivative
on the split-step $k$-grid. The machinery's refuter is the harmonic control: the same packet in
$V=x^2/2$, where quadratic $V$ closes Ehrenfest on $\langle x\rangle$ exactly, so any gap there
is an error. On the quartic run both residuals $\frac{d\langle x\rangle}{dt}-\langle p\rangle$
and $\frac{d\langle p\rangle}{dt}+\langle x^3\rangle$ must vanish to tolerance (the $d/dt$ of a
sampled moment series has a finite-difference noise floor, stated as such) while the classical
curve departs. Close with the narrow-packet limit shrinking the gap: classicality is a property
of the state's spread, not of the theorem.

### 5.10 [BSc] How do I confirm that norm and probability current are conserved along an `NDSolve` evolution (a numerical-fidelity check)?

Launch the boosted supergaussian $\psi_0\propto e^{-(x-x_0)^4}e^{ik_0x}$ from $x_0=-8$ across
the reflectionless well $V=-\operatorname{sech}^2x$, the $\lambda=1$ Poschl-Teller well with the
single bound state $E_0=-\tfrac12$: $|t(k)|=1$ for every momentum component, so any reflected
lobe above tolerance is guaranteed spurious and the fidelity checks ride a striking exact
property. Derive the current symbolically before the run, from the TDSE for generic $\psi$ and
real $V$, writing $\rho=\psi\bar\psi$ (never differentiate $|\psi|^2$ literally: `D[Abs[..]]`
yields no delta, a C3-cited trap) to get $\partial_t\rho+\partial_x j=0$ with
$j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$. Run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$ and grade three checks of increasing strength on the same run: the
norm grid-sum drift; the continuity residual on the grid, shrinking under refinement; and
station-flux bookkeeping, $\int_0^T j(x_s,t)\,dt$ equal to the probability mass transferred past
$x_s$, with the transmitted mass $\to1$ as the reflectionless anchor. Exhibit the norm-blindness
lesson live, an undersized box reflecting at norm 1, so norm alone certifies nothing about
walls. Close on that ranking: a fidelity check that cannot fail teaches nothing.
