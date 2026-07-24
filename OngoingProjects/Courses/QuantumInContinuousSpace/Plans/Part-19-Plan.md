# Part 19 Plan: Path integrals

Part 19, "Path integrals": 6 questions, all MSc.
Class census (`Route-Table.md`): C0 for 19.1, 19.2, 19.3, 19.5 (no differential equation, C0 row); C9 for 19.4, 19.6 (WKB / semiclassical asymptotics, C9 verdict).

Everything in this part rests on the time-sliced propagator
$K(x,x';t)=\langle x\vert e^{-i\hat{H}t}\vert x'\rangle=\lim_{N\to\infty}\int\prod_{j=1}^{N-1}dx_{j}\,(2\pi i\epsilon)^{-N/2}e^{iS_{N}}$
with step $\epsilon=t/N$ and discretized action $S_{N}=\sum_{j}\bigl[(x_{j+1}-x_{j})^{2}/(2\epsilon)-\epsilon V(x_{j})\bigr]$;
the Wick rotation $t\to-i\beta$ turns the weight into $e^{-S_{E}}$ and the trace over periodic paths into the
partition function $Z(\beta)=\operatorname{Tr}e^{-\beta\hat{H}}$; stationary phase $\delta S=0$ selects the classical
path with action $S_{\mathrm{cl}}$, and the Gaussian fluctuations around it give the Van Vleck form
$K_{\mathrm{sc}}=\sqrt{-\partial^{2}S_{\mathrm{cl}}/(\partial x\,\partial x')/(2\pi i)}\;e^{iS_{\mathrm{cl}}}$, exact
whenever the action is quadratic. Every kernel-checkable claim in the part reduces to three facts: a propagator
solves the time-dependent Schrodinger equation (TDSE), collapses to $\delta(x-x')$ as $t\to0^{+}$, and composes,
$\int K(x,y;t_{1})K(y,x';t_{2})\,dy=K(x,x';t_{1}+t_{2})$.

### 19.1 [MSc] How do I derive the free-particle propagator from the Gaussian path integral?

Slice the free evolution into $N$ steps of width $\epsilon=t/N$, each carrying the short-time kernel
$K_{\epsilon}(x,y)=e^{i(x-y)^{2}/(2\epsilon)}/\sqrt{2\pi i\epsilon}$, and let a single Integrate do the two-slice
composition $\int K_{\epsilon_{1}}(x,y)K_{\epsilon_{2}}(y,x')\,dy$ with symbolic widths under assumptions
$\epsilon_{1}>0$, $\epsilon_{2}>0$, $x,x'$ real (the Fresnel Gaussian stalls or picks a wrong branch without them;
if the oscillatory form still refuses, evaluate the Wick-rotated Gaussian and continue back, recording the swap):
the output must be $K_{\epsilon_{1}+\epsilon_{2}}(x,x')$, which is the induction step. A symbolic Nest of the same
composition at fixed $N=4$ collapses the sliced integral to $K_{4\epsilon}$, and with the induction established the
general $N$-slice form is $K_{N\epsilon}$, so Limit under $\epsilon=t/N$, $N\to\infty$ lands on
$K=e^{i(x-x')^{2}/(2t)}/\sqrt{2\pi i t}$, the limit taken, not asserted. Verify with checks that can refute and
reuse the derived $K$: FullSimplify of the TDSE residual $i\partial_{t}K+\tfrac{1}{2}\partial_{x}^{2}K$ to zero
under $t>0$, and the $t\to0^{+}$ delta limit taken as $\lim_{t\to0^{+}}\int K(x,x';t)f(x')\,dx'=f(x)$ on an explicit
test function through Integrate then Limit, never NIntegrate, which silently returns $0.$ on an emerging point
measure (C0 row). Close by reading the exponent as the classical action of the straight-line path,
$S_{\mathrm{cl}}=(x-x')^{2}/(2t)$: the free kernel is its own stationary-phase form.

### 19.2 [MSc] How do I get the oscillator propagator (Mehler) from the path integral?

Run the same $N$-slice discretization on the oscillator ($\omega=1$): integrating the $N-1$ interior points of
$S_{N}=\sum_{j}\bigl[(x_{j+1}-x_{j})^{2}/(2\epsilon)-\epsilon x_{j}^{2}/2\bigr]$ leaves a Gaussian whose prefactor
carries the determinant $d_{N-1}$ of the tridiagonal fluctuation matrix, obeying the Gelfand-Yaglom recursion
$d_{j+1}=(2-\epsilon^{2})d_{j}-d_{j-1}$, $d_{0}=1$, $d_{1}=2-\epsilon^{2}$; solve it with RSolve (verify at
authoring; the fallback is the explicit closed form $d_{j}=\sin((j+1)\theta)/\sin\theta$ with
$\cos\theta=1-\epsilon^{2}/2$, certified by FullSimplify of the recursion residual), cross-check the recursion at
fixed $N=3$ against the brute-force sliced Integrate, and take Limit of $\epsilon\,d_{N-1}\to\sin t$ under
$\epsilon=t/N$, $N\to\infty$; the boundary quadratic form closes through cofactors of the same tridiagonal matrix,
landing on the Mehler kernel $K=e^{i[(x^{2}+x'^{2})\cos t-2xx']/(2\sin t)}/\sqrt{2\pi i\sin t}$. Verify by the
refutable pair reusing the derived $K$: the TDSE residual with $V=x^{2}/2$ FullSimplified to zero under $0<t<\pi$,
and the $t\to0^{+}$ delta limit against a test function via Integrate and Limit. The honest edge is the caustic at
$t=\pi$: $\sin t\to0$, the prefactor diverges, and the kernel collapses onto $\delta(x+x')$; state the derived
branch as $0<t<\pi$ and name the Maslov phase beyond it as out of scope. Close on the caustic as the time where
every classical path refocuses regardless of $x'$.

### 19.3 [MSc] How do I use the imaginary-time path integral to compute the partition function of a finite system?

Compute the oscillator $Z(\beta)$ three ways and make them collide: the spectrum sum
$\sum_{n\geq0}e^{-\beta(n+1/2)}$ closed by Sum; the exact $1/(2\sinh(\beta/2))$, their identity certified by
FullSimplify under $\beta>0$; and the discretized imaginary-time trace as a finite matrix product: grid $[-L,L]$ by
Subdivide, symmetrized Euclidean short-time matrix
$T_{jk}=\sqrt{1/(2\pi\epsilon)}\,e^{-(x_{j}-x_{k})^{2}/(2\epsilon)-\epsilon(V(x_{j})+V(x_{k}))/2}\,\delta x$ with
$V=x^{2}/2$ built by Outer on packed machine arrays (apply N to the grid first so MatrixPower stays packed), then
$Z_{N}=\operatorname{Tr}T^{N}$ via MatrixPower and Tr at $\epsilon=\beta/N$. The refuting check: the deviation of
$Z_{N}$ from the closed form at fixed $\beta$ must fall fourfold per doubling of $N$ (the symmetrized Trotter error
is $O(\epsilon^{2})$); a plateau above that law exposes grid truncation, so sweep $L$ as well. Anchor the classical
limit: Series of $1/(2\sinh(\beta/2))$ at $\beta\to0$ gives $1/\beta$, and the phase-space integral
$\int dx\,dp\,e^{-\beta(p^{2}+x^{2})/2}/(2\pi)$ closed by Integrate gives the same $1/\beta$, equipartition
regained. Close on the two ends of $\beta$: quantum ground-state dominance $e^{-\beta/2}$ against the classical
$1/\beta$.

### 19.4 [MSc] How do I apply the stationary-phase (semiclassical) approximation and identify the classical action?

State the Van Vleck stationary-phase propagator
$K_{\mathrm{sc}}=\sqrt{-\partial^{2}S_{\mathrm{cl}}/(\partial x\,\partial x')/(2\pi i)}\,e^{iS_{\mathrm{cl}}}$ and
identify the classical action by computing it: DSolve the classical boundary problem $\ddot{x}=-x$, $x(0)=x'$,
$x(t)=x$, then $S_{\mathrm{cl}}=\int_{0}^{t}(\dot{x}^{2}-x^{2})/2\,ds$ by Integrate, closing to
$[(x^{2}+x'^{2})\cos t-2xx']/(2\sin t)$; D gives the mixed derivative $-1/\sin t$, and the assembled kernel is
exact for the oscillator, confirmed by the refutable pair reusing the assembled $K_{\mathrm{sc}}$: FullSimplify of
the TDSE residual to zero under $0<t<\pi$ and the $t\to0^{+}$ delta limit (the C9 verdict rates this member
low-risk pure symbolic Gaussian quadrature). Then locate the boundary of stationary phase: with
$V=x^{2}/2+\lambda x^{4}$, expand the path-integral weight about the classical path; the quartic fluctuation term
contributes at $O(\lambda)$ through the Gaussian moments $\int_{0}^{t}x_{\mathrm{cl}}(s)^{2}G(s,s)\,ds$ and
$\int_{0}^{t}G(s,s)^{2}\,ds$ with $G(s,s)=\sin s\,\sin(t-s)/\sin t$ the Dirichlet fluctuation propagator, each
closed by Integrate and neither vanishing, so the correction enters at exactly $O(\lambda)$ and at no lower order.
Close: for a quadratic action stationary phase is not an approximation at all, and $\lambda x^{4}$ is precisely
where it becomes one.

### 19.5 [MSc] How do I evaluate a discretized path integral numerically for a simple potential?

Get the quartic-well ground energy from pure matrix multiplication: on a Subdivide grid over $[-L,L]$ build the
symmetrized Euclidean transfer matrix
$T_{jk}=\sqrt{1/(2\pi\epsilon)}\,e^{-(x_{j}-x_{k})^{2}/(2\epsilon)-\epsilon(V(x_{j})+V(x_{k}))/2}\,\delta x$ with
$V=x^{4}/4$ by Outer on packed machine arrays, and read the ground energy from the large-$\beta$ decay of the
trace: $-\epsilon^{-1}\ln\bigl(\operatorname{Tr}T^{N+1}/\operatorname{Tr}T^{N}\bigr)$, computed with MatrixPower
and Tr (repeated Dot of $T$ onto a vector doubles as an internal power-iteration consistency check), flattens to a
plateau once $e^{-\beta(E_{1}-E_{0})}$ is negligible; store the estimate as en0, never E, which is Euler's number
(C0 row). The Trotter step $\epsilon$ is the swept knob: halve it at fixed $\beta$ and the plateau must converge as
$O(\epsilon^{2})$; sweep $L$ too, since a short box lifts $E_{0}$ silently. The refuting benchmark is independent
machinery: the C2-certified FD recipe (SparseArray Band tridiagonal Hamiltonian, Eigenvalues with the Arnoldi
shift, Richardson in the grid step) pins $E_{0}$ far below the transfer-matrix error, and the extrapolated plateau
must land on it within its own $O(\epsilon^{2})$ bar or the discretization is wrong. Close spectrally:
$Z(\beta)=e^{-\beta E_{0}}(1+e^{-\beta(E_{1}-E_{0})}+\dots)$, so the pre-plateau bend of the same data measures the
gap $E_{1}-E_{0}$ as the leading correction.

### 19.6 [MSc] How do I compute the double-well tunneling splitting from an instanton (imaginary-time bounce)?

Deep-barrier quartic double well $V=(x^{2}-a^{2})^{2}/2$: compute the instanton action
$S_{0}=\int_{-a}^{a}\sqrt{2V(x)}\,dx$ by the C9 quadrature discipline, scaling the turning points out with $x=au$
before Integrate (raw Integrate with symbolic turning-point limits stalls, and the antiderivative-plus-Limit
Newton-Leibniz route silently returns $0$ when the limit direction lands in the forbidden region: the C9 traps
that bite exactly this integral), which closes to $S_{0}=a^{3}\int_{-1}^{1}(1-u^{2})\,du=4a^{3}/3$. The prediction
$\Delta E\propto e^{-S_{0}}$ is asymptotic-only per the C9 verdict (the one-loop prefactor has no kernel
evidence), so never demand equality: measure the splitting $\Delta E=E_{1}-E_{0}$ with the C2-certified FD recipe
(SparseArray Band + Eigenvalues with the Arnoldi shift + Richardson; certified at the $3\times10^{-5}$ splitting
for $a=2.2$; below roughly $10^{-7}$ switch to the parity-resolved half-domain solves the C2 verdict prescribes)
over a sweep of $a$ from about $1.6$ to $2.2$, and compare $\ln\Delta E$ trends: the local slope
$d\ln\Delta E/da$ extracted with Differences must track $-dS_{0}/da=-4a^{2}$ increasingly well as $a$ grows, the
trend that refutes or confirms the exponent while leaving the prefactor alone. Close: $S_{0}$ is the Euclidean
action of the path that crosses the forbidden barrier in imaginary time, a number invisible to every order of
perturbation theory around either minimum.
