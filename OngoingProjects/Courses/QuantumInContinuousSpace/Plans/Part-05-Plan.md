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

The solver contract below is authoring guidance for this plan, not material for the Part opening.
Seven of the eight C3 entries, all but 5.5, share the certified contract of the C3 verdict in
`Route-Table.md`: primary `NDSolveValue` with `Method -> {"MethodOfLines",
"SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> n, "MaxPoints" -> n,
"DifferenceOrder" -> "Pseudospectral"}}`, periodic identification $\psi(t,L_1)=\psi(t,L_2)$, and
`AccuracyGoal -> 10, PrecisionGoal -> 10`. Pinning `MinPoints` and `MaxPoints` to the same $n$ is
what makes $n$ the swept knob rather than a solver choice, and every resolution sweep below moves
that one number over the verdict's measured ladder $n=65,129,257,513$; the tight goals are
load-bearing, since defaults are silently wrong at the $0.25\%$ level with no message and
refinement at default goals is non-monotone. Cross-check independently with a hand-built Strang
split-step Fourier propagator on `Fourier`/`InverseFourier` (`FourierParameters -> {1,-1}`,
packed arrays). Size the box from the state's own energy content and hold it fixed across every
sweep, because norm conservation does not detect wall reflection. Take moments and probabilities
as $dx$-weighted grid sums, never `NIntegrate` on an interpolant. Entry 5.5 leaves this contract
and states its own discretization.

### 5.1 [BSc] How do I show that a stationary state evolves only by the phase $e^{-iEt}$, while a superposition acquires genuine time dependence?

Work in the Poschl-Teller family $V=-\tfrac{\lambda(\lambda+1)}{2}\operatorname{sech}^2x$, whose
bound levels are $E_n=-\tfrac{(\lambda-n)^2}{2}$, and specialize to $\lambda=2$, that is
$V=-3\operatorname{sech}^2x$ with exactly two states, $\psi_0\propto\operatorname{sech}^2x$ at
$E_0=-2$ and $\psi_1\propto\operatorname{sech}x\tanh x$ at $E_1=-\tfrac12$: a finite spectrum
makes the beat a single closed-form frequency with no truncation question. This is C0 machinery
per the `Route-Table.md` C0 row: earn both levels with `D` + `FullSimplify` residuals
$\hat H\psi_n-E_n\psi_n\to0$, then reduce the densities of $e^{-iE_0t}\psi_0$ and of the
equal-weight superposition with `ComplexExpand` under `Element[{x, t}, Reals]`. The stationary
density is $t$-free identically; the superposition carries a cross term $\cos(\Delta E\,t)$ with
$\Delta E=\lambda-\tfrac12$, so the whole family beats, and $\lambda=2$ gives $\Delta E=\tfrac32$
with period $4\pi/3$, refutable by a numeric spot value of the density at one $(x,t)$. Close by
taking `Limit` of that same cross term as $\Delta E\to0$: the density goes static, so genuine
time dependence lives in energy differences, not in superposition itself.

### 5.2 [BSc] How do I integrate the time-dependent Schrodinger equation directly as a PDE with `NDSolve` for a given initial $\psi(x,0)$?

Release the cusp packet $\psi_0\propto e^{-|x-x_0|}$, $x_0=3$, in the harmonic trap $V=x^2/2$:
its momentum amplitude falls only as $(1+k^2)^{-1}$, a genuine grid stress a smooth default
packet would never apply. Integrate under the shared C3 contract, with a default-goals twin run
alongside to exhibit the silent-error trap live. Size the box from the state's own energy, edge
$\gtrsim\sqrt{2(\langle H\rangle+4\sigma_H)}$ plus several cusp decay lengths (the tail falls as
$e^{-|x-x_0|}$, so a handful of units suffices), and hold that box fixed while $n$ doubles, so
the convergence law is measured on one geometry rather than on a box that moves with the grid.
The trap supplies exact anchors that reuse $\psi_0$: the density at $t=\pi$ must equal
$\rho_0(-x)$, the recurrence overlap $|\langle\psi_0|\psi(2\pi)\rangle|$ must return to 1, and
the grid-sum drift of $\langle H\rangle$ refutes independently, with the split-step propagator as
the cross-check. Close with the measured convergence law in $n$, tail-limited by the cusp, set
against the spectral convergence a smooth packet would show.

### 5.3 [BSc] How do I follow a free Gaussian wave packet's spreading and group velocity, both analytically and by `NDSolve`?

The question names the free Gaussian, so it is the carrier, with width $\sigma_0$ and drift $k_0$
kept symbolic. Get the closed form analytically first: the C3 verdict certifies `DSolve` on the
free IVP only for the fixed packet whose $t=0$ form is $e^{-x^2/4}e^{2ix}$ ($\sigma_0=1$,
$k_0=2$), so try symbolic-parameter `DSolve` and tag it (verify at authoring), with the transform
route as the named fallback: `FourierTransform` of $\psi_0$, multiply by $e^{-ik^2t/2}$,
transform back under `Assumptions -> {sigma0 > 0, Element[{x, t}, Reals]}`, then residual-verify
the resulting closed form against the TDSE. Run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$, compare pointwise, and check the moment laws $\langle x\rangle(t)=k_0t$
and $\sigma^2(t)=\sigma_0^2+t^2/(4\sigma_0^2)$ from grid sums; then exhibit the norm-blind
boundary trap deliberately, an undersized box reflecting at norm 1. The discriminating contrast
is a same-width sech packet whose rescaled density fails to collapse onto its initial shape:
self-similar spreading is a Gaussian privilege, not a free-particle law. Close with the ballistic
limit $\sigma(t)\to t/(2\sigma_0)$.

### 5.4 [MSc] How do I implement the split-step Fourier propagator from scratch and benchmark it against `NDSolve`?

Build the Strang map $e^{-iV\,dt/2}e^{-iK\,dt}e^{-iV\,dt/2}$ from scratch on
`Fourier`/`InverseFourier` (`FourierParameters -> {1,-1}`, packed arrays, the $k$-grid convention
checked before use) and benchmark it against the shared `NDSolveValue` contract on a packet in
the symmetric quartic double well $V=(x^2-a^2)^2/2$, started as a local-harmonic-width Gaussian
displaced onto the left minimum. Pin $a=1.5$, where the barrier $a^4/2\approx2.5$ is comparable
to the local zero-point energy $\omega/2=a=1.5$: the doublet splitting is then of order $10^{-1}$
and the tunneling period $\pi/\Delta E$ of order $10^2$, inside reach of the integrator, whereas
the C2-measured deep cases $a=2$ (splitting $7.58\times10^{-4}$, period about $8\times10^3$) and
$a=2.2$ ($2.9\times10^{-5}$, about $2\times10^5$) sit three to five orders beyond the probed
windows and would be dominated by drift. Confirm the splitting with a quick eigenvalue check (the
C2 recipe) before committing run time. Verify that the pointwise disagreement between the two
integrators shrinks under $dt$ halving and goal tightening at fixed $n$, and measure the Strang
order in a coarse-$dt$ window, since the verdict flags the order ratio as contaminated at the
$10^{-13}$ floor. The discriminator is per-step norm, measured for this run rather than
transplanted: the verdict's $10^{-13}$ split-step and $10^{-5}$ `NDSolveValue` drifts are
free-packet references at $T\lesssim6$, not predictions for a longer double-well run. Close with
the $V\to0$ limit, which must reproduce the free-Gaussian closed form.

### 5.5 [MSc] How do I evolve a state by expanding it in energy eigenstates and propagating term by term, and observe wave-packet revivals?

Expand the off-center triangular packet peaked at $x_0=L/3$ in the infinite well of width $L$:
the kink populates many modes with $c_n$ decaying only as $1/n^2$. Certify the spectrum before
leaning on it, since the whole revival argument needs $E_n$ exactly quadratic in $n$: substitute
$\psi_n=\sqrt{2/L}\sin(n\pi x/L)$ with symbolic $n$ into the stationary equation, `FullSimplify`
the residual to $0$ under `Element[n, Integers]`, and check $\psi_n(0)=\psi_n(L)=0$; then
$E_n=n^2\pi^2/(2L^2)$ gives $E_nT_{\mathrm{rev}}=2\pi n^2$ at $T_{\mathrm{rev}}=4L^2/\pi$. The
primary is the exact eigenbasis per the `Route-Table.md` cross-class hand-off: $c_n$ by
`Integrate` in closed form with `Assumptions -> Element[n, Integers]`, a truncated `Sum` bounded
by the tail mass $\sum_{n>N}|c_n|^2$. The revival is then exact by construction, so the refutable
content is the mirror identity $\psi(x,T_{\mathrm{rev}}/2)=-\psi_0(L-x)$, compared pointwise
between the truncated series and the reflected initial function, plus an independent numeric
evolution whose overlap $|\langle\psi_0|\psi(T_{\mathrm{rev}})\rangle|$ must return to 1. That
run leaves the periodic-pseudospectral primary because these walls are physical, and the
verdict's only Dirichlet evidence is a free packet striking an artificial wall, which does not
transfer; pseudospectral on a non-periodic grid with Dirichlet conditions is gated in no verdict,
so take `"DifferenceOrder" -> 4` with `MinPoints` and `MaxPoints` pinned and swept, and tag a
pseudospectral variant (verify at authoring) if the kink's high modes demand it. As a long-time
member it inherits the flagged risks, default-goal drift growing with integration length and
pseudospectral cost above $n=513$ unprobed, so the numeric leg carries a goal-and-resolution
sweep. Close with the quarter-revival two-image cat at $T_{\mathrm{rev}}/4$.

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
$\psi_0\propto\left(1+\cos\frac{\pi(x-x_c)}{w}\right)e^{ik_0x}$ on $|x-x_c|<w$, with $x_c=-20$,
$w=5$, $k_0=2$, at the Eckart barrier $V_0\operatorname{sech}^2(x/a)$ with $V_0=2$ and $a=1$.
That places the mean energy $E_0=k_0^2/2=2$ mid-ramp on the exact anchor
$T(E)=\sinh^2(\pi ka)\,/\,[\sinh^2(\pi ka)+\cosh^2(\tfrac\pi2\sqrt{8V_0a^2-1})]$ with
$k=\sqrt{2E}$, valid for $8V_0a^2>1$: here $\sinh^2(2\pi)=7.17\times10^4$ against
$\cosh^2(\tfrac\pi2\sqrt{15})=4.83\times10^4$, so $T=0.60$ and both fragments are large enough to
see and to weigh. The barrier vanishes nowhere, so compact support does not give zero initial
overlap; it gives a quantitatively negligible one, $V/V_0=\operatorname{sech}^2(d/a)\approx
4e^{-2d/a}=3.7\times10^{-13}$ at the nearest support edge $d=15$, an energy shift far below
tolerance, while the real payoff is an unambiguous late-time mass split and a clean initial
momentum amplitude. Run the shared C3 contract with box edge $\gtrsim k_0T+8\sigma_T$ on both
sides, since both fragments travel, and the split-step cross-check. The refuting check is the
late-time transmitted mass (grid sum beyond the barrier) against
$\int T(k)\,|\tilde\psi_0(k)|^2\,dk$ built from the same $\psi_0$ via `FourierTransform`, with
$R+T=1$ from the same sums. The C3 verdict flags 5.7's ringing as its unmeasured open risk: the
smooth Eckart barrier retires the potential-side Gibbs ringing while keeping the exact anchor,
but the packet's $C^1$ support seams still seed high-$k$ content, so a resolution sweep in $n$ of
the seam ringing is mandatory at authoring. Close with $T\to1$ at $E\gg V_0$ set against the
deep-tunneling suppression.

### 5.8 [MSc] How do I construct the propagator (kernel) $K(x,t;x',0)$ for the free particle and for the oscillator (the Mehler kernel)?

Construct the free kernel $K_0=(2\pi it)^{-1/2}e^{i(x-x')^2/(2t)}$ from the momentum integral
$\frac{1}{2\pi}\int e^{ik(x-x')-ik^2t/2}\,dk$ by `Integrate` with a $t\to t-i\epsilon$ regulator
and `Limit`, and the Mehler kernel
$K_{\mathrm{osc}}=(2\pi i\sin t)^{-1/2}\exp\!\left\{\tfrac{i}{2\sin t}\left[(x^2+x'^2)\cos t-2xx'\right]\right\}$
by closing the Wick-rotated eigen-sum $\sum_n e^{-(n+1/2)\tau}\psi_n(x)\psi_n(x')$ with `Sum`
(`HermiteH`) and continuing $\tau\to it$; if `Sum` refuses the bilinear Hermite series, posit the
closed form and certify it instead by the imaginary-time residual
$\partial_\tau K=\tfrac12\partial_{xx}K-\tfrac12x^2K$ and the $\tau\to0^+$ delta limit. At real
$t$ the eigen-sum oscillates without converging pointwise, which is why the comparison runs
through imaginary time (C0 per the `Route-Table.md` C0 row). Each kernel earns its name three
ways, any failure refuting it: the TDSE residual in $(x,t)$ is identically $0$ under `D` +
`FullSimplify`; the composition law $\int K(x,t_1;x'')\,K(x'',t_2;x')\,dx''=K(x,t_1{+}t_2;x')$
closes; and the $t\to0^+$ action on a test packet returns the packet, the $\delta(x-x')$ limit
taken through a test function, never as a bare limit. Add a numeric leg: propagate one packet by
explicit kernel quadrature at fixed $t$ and compare against a locally rebuilt split-step Fourier
map. Then let the free kernel act on the Berry-Balazs Airy beam (`AiryAi`), which propagates as
$\operatorname{Ai}(x-t^2/4)\,e^{it(x-t^2/6)/2}$, an accelerating profile out of a free kernel.
Close at the caustics $t=n\pi$, where the Mehler prefactor diverges and the kernel degenerates to
mirrored delta transport, with the $\sin t\to t$ free-kernel limit as the bridge between the two.

### 5.9 [BSc] How do I verify Ehrenfest's theorem numerically, that $\langle x\rangle$ and $\langle p\rangle$ obey the classical equations of motion?

Evolve the asymmetric two-lobe packet, the normalized unequal-weight sum
$e^{-(x-1.2)^2}+\tfrac12 e^{-(x+0.6)^2}$, in the quartic well $V=x^4/4$, where Ehrenfest's actual
content is the difference between $\frac{d\langle p\rangle}{dt}=-\langle V'(x)\rangle=
-\langle x^3\rangle$ and the classical $-V'(\langle x\rangle)=-\langle x\rangle^3$. That gap is
$\langle x^3\rangle-\langle x\rangle^3=3\langle x\rangle\operatorname{Var}(x)+\mu_3$, so what the
example must supply is a large variance at nonzero $\langle x\rangle$ together with a nonzero
third central moment $\mu_3$; skewness alone is not the mechanism, and even a symmetric packet
displaced off the origin already opens the first term. Run the shared C3 contract with the box
beyond the packet's classical turning points, moments as $dx$-weighted grid sums and
$\langle p\rangle$ from the spectral derivative on the $k$-grid. The machinery's refuter is the
harmonic control, the same packet in $V=x^2/2$, where quadratic $V$ closes Ehrenfest on
$\langle x\rangle$ exactly, so any gap there is an error. On the quartic run both residuals
$\frac{d\langle x\rangle}{dt}-\langle p\rangle$ and
$\frac{d\langle p\rangle}{dt}+\langle x^3\rangle$ must vanish to tolerance, allowing for the
finite-difference noise floor of $d/dt$ on a sampled moment series, while the classical curve
departs. Close with the narrow-packet limit shrinking both terms of the gap: classicality is a
property of the state's spread, not of the theorem.

### 5.10 [BSc] How do I confirm that norm and probability current are conserved along an `NDSolve` evolution (a numerical-fidelity check)?

Launch the boosted supergaussian $\psi_0\propto e^{-(x-x_0)^4}e^{ik_0x}$ from $x_0=-8$ across the
reflectionless well $V=-\operatorname{sech}^2x$, the $\lambda=1$ Poschl-Teller well with the
single bound state $\psi_b=\operatorname{sech}(x)/\sqrt2$ at $E_b=-\tfrac12$: $|t(k)|=1$ for
every momentum component, so any reflected lobe above tolerance is guaranteed spurious. Derive
the current symbolically before the run, from the TDSE for generic $\psi$ and real $V$, writing
$\rho=\psi\bar\psi$ (never differentiate $|\psi|^2$ literally, since `D[Abs[..]]` yields no
delta, a C3-cited trap) to get $\partial_t\rho+\partial_x j=0$ with
$j=\operatorname{Im}(\bar\psi\,\partial_x\psi)$. Run the shared C3 contract with box edge
$\gtrsim k_0T+8\sigma_T$ and grade three checks of increasing strength on the same run: the norm
grid-sum drift; the continuity residual on the grid, shrinking under refinement in $n$; and
station-flux bookkeeping, $\int_0^T j(x_s,t)\,dt$ against the probability mass transferred past
$x_s$. Reflectionless does not mean everything transmits: the well captures its bound component
permanently, so predict $T=1-|\langle\psi_b|\psi_0\rangle|^2$ by `Integrate` from the same
definitions, a truthful and strictly stronger anchor than $T\to1$. Exhibit the norm-blindness
lesson live, an undersized box reflecting at norm 1, so norm alone certifies nothing about walls.
Close on that ranking: a fidelity check that cannot fail teaches nothing.
