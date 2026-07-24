# Part 16 Plan: Density operators, mixed states, and the Wigner function

Six questions, all MSc. Class census (`Route-Table.md`): C0 for all of 16.1 through 16.6 (kernels,
quadrature, and integral transforms; no differential-equation solver).

## Common ground

A mixed state enters as a kernel $\rho(x,x')=\sum_{k}p_k\,\psi_k(x)\,\overline{\psi_k(x')}$ with
unit trace $\int\rho(x,x)\,dx=1$ and purity
$\operatorname{Tr}\rho^2=\iint|\rho(x,x')|^2\,dx\,dx'$, equal to 1 exactly when the state is pure;
a subsystem of a pure pair state is the partial trace
$\rho_A(x,x')=\int\Psi(x,y)\,\overline{\Psi(x',y)}\,dy$. Phase space sees $\rho$ through the
Wigner-Weyl transform $W(x,p)=\frac{1}{2\pi}\int e^{-ipy}\,\rho(x+\tfrac{y}{2},x-\tfrac{y}{2})\,dy$,
whose marginals are the position and momentum densities and which obeys $|W|\le 1/\pi$; smoothing
orders the representations, with $Q(\alpha)=\frac{1}{\pi}\langle\alpha|\rho|\alpha\rangle\ge0$ at
the regular end and $P$, defined by $\rho=\int P(\alpha)\,|\alpha\rangle\langle\alpha|\,d^{2}\alpha$,
at the singular end. Evolution is the Moyal equation
$\partial_{t}W=-p\,\partial_{x}W+V'(x)\,\partial_{p}W-\tfrac{1}{24}V'''(x)\,\partial_{p}^{3}W+\ldots$,
which truncates to the classical Liouville flow whenever $V$ is quadratic. Oscillator entries set
$\omega=1$, so the coherent state $|\alpha\rangle$ is a vacuum-width Gaussian centered at
$(x_0,p_0)=(\sqrt{2}\,\operatorname{Re}\alpha,\sqrt{2}\,\operatorname{Im}\alpha)$.

### 16.1 [MSc] How do I represent a density operator as a kernel $\rho(x,x')$ and compute its purity $\iint|\rho(x,x')|^2\,dx\,dx'$?

Build the two levels of the Poschl-Teller well $V=-3\operatorname{sech}^2x$ in closed form,
$\psi_0(x)=\tfrac{\sqrt{3}}{2}\operatorname{sech}^2x$ at $E_0=-2$ and
$\psi_1(x)=\sqrt{3/2}\,\operatorname{sech}x\tanh x$ at $E_1=-\tfrac{1}{2}$, earning both by a
`FullSimplify` residual of the stationary equation and unit norms by `Integrate` (energies named
`en`, never `E`). Form the discriminating pair: the equal-population mixture
$\rho_{\mathrm{mix}}(x,x')=\tfrac{1}{2}\psi_0(x)\psi_0(x')+\tfrac{1}{2}\psi_1(x)\psi_1(x')$ against
the coherent superposition $\psi=(\psi_0+i\psi_1)/\sqrt{2}$ with kernel
$\rho_{\mathrm{sup}}(x,x')=\psi(x)\overline{\psi(x')}$; the relative phase $i$ makes the two
diagonals coincide pointwise (`FullSimplify` of the difference to 0), so no position histogram
separates the pair. Purity does: `Integrate` the double integral $\iint|\rho|^2\,dx\,dx'$ under
reality assumptions to exactly $\tfrac{1}{2}$ for the mixture and $1$ for the superposition, and
cross-check both numbers as $\sum_{k}p_k^2$ from the `Integrate`-verified orthonormality of the
levels, reusing the same definitions, so a wrong kernel fails one route or the other. Close on
where the difference lives: `FullSimplify` $\rho_{\mathrm{sup}}-\rho_{\mathrm{mix}}$ down to the
pure off-diagonal cross term $\propto i\,(\psi_1(x)\psi_0(x')-\psi_0(x)\psi_1(x'))$, the coherence
a position measurement never sees and decoherence destroys.

### 16.2 [MSc] How do I obtain a reduced density matrix by tracing out one particle (integrating over its coordinate)?

Take the coupled pair $H=\tfrac{1}{2}(p_1^2+p_2^2)+\tfrac{1}{2}(x_1^2+x_2^2)+\tfrac{1}{2}g(x_1-x_2)^2$
with normal modes $x_{\pm}=(x_1\pm x_2)/\sqrt{2}$ at $\omega_{+}=1$, $\omega_{-}=\sqrt{1+2g}$, and
write the exactly Gaussian ground state
$\Psi(x_1,x_2)=(\omega_{+}\omega_{-}/\pi^2)^{1/4}\,e^{-\omega_{+}x_{+}^2/2-\omega_{-}x_{-}^2/2}$,
earned by the 2D stationary residual (`D` plus `FullSimplify`, energy
$\tfrac{1}{2}(\omega_{+}+\omega_{-})$). Trace out the partner by `Integrate` over $y$ with
`Assumptions` $\omega_{\pm}>0$: $\rho_A(x,x')=\int\Psi(x,y)\,\overline{\Psi(x',y)}\,dy$ is again
Gaussian, $\propto e^{-\gamma(x^2+x'^2)/2+\beta xx'}$ with $\gamma,\beta$ read off by `Coefficient`.
Purity comes from the Gaussian double `Integrate` of $|\rho_A|^2$, a
$2\sqrt{\omega_{+}\omega_{-}}/(\omega_{+}+\omega_{-})$-type closed form (verify at authoring); the
independent route that could refute it computes $\langle x^2\rangle$ and $\langle p^2\rangle$ from
the same kernel (`D` across the diagonal before setting $x'=x$, then `Integrate`), forms the
symplectic eigenvalue $\nu=\sqrt{\langle x^2\rangle\langle p^2\rangle}$, and demands purity
$=1/(2\nu)$, after which the entanglement entropy follows as
$S=(\nu+\tfrac{1}{2})\ln(\nu+\tfrac{1}{2})-(\nu-\tfrac{1}{2})\ln(\nu-\tfrac{1}{2})$ (verify at
authoring). Anchor with `Limit` at $g\to0$: purity 1 and $S=0$ exactly; a large-$g$ `Series` shows
the purity falling like $(\omega_{+}/\omega_{-})^{1/2}$ toward zero. Close: the global state stays
pure at every $g$, yet each oscillator alone is strictly mixed once $g\neq0$; mixedness of a
subsystem is the operational face of entanglement.

### 16.3 [MSc] How do I compute the Wigner quasi-probability function by the Wigner-Weyl transform?

Define the $n=1$ Fock state $\psi_1(x)=\sqrt{2}\,\pi^{-1/4}\,x\,e^{-x^2/2}$ and the transform
$W(x,p)=\frac{1}{2\pi}\int e^{-ipy}\,\psi_1(x+\tfrac{y}{2})\,\overline{\psi_1(x-\tfrac{y}{2})}\,dy$,
then evaluate by one direct `Integrate` under `Element[{x, p}, Reals]`: the closed form
$W_1(x,p)=\frac{1}{\pi}\,(2x^2+2p^2-1)\,e^{-x^2-p^2}$ lands with $W_1(0,0)=-1/\pi$ exactly, the
extremal value permitted by the bound $|W|\le1/\pi$. Three checks reuse the definitions and could
each refute the result: `Integrate` over $p$ must return $|\psi_1(x)|^2$, `Integrate` over $x$ must
return $|\tilde{\psi}_1(p)|^2$ with $\tilde{\psi}_1$ obtained independently by `FourierTransform`,
and the full phase-plane `Integrate` must give 1. `Solve` locates the negative region as the disk
$x^2+p^2<\tfrac{1}{2}$ around the origin. Close: a genuine probability density is nowhere negative,
so the clean $-1/\pi$ certifies the Fock state as non-classical in the strongest phase-space sense,
and no state can dip deeper.

### 16.4 [MSc] How do I compute the Wigner function of coherent, number, and cat states and exhibit its negativity?

Normalize the even cat at $\alpha=2$: with $x_0=\sqrt{2}\,\alpha$,
$\psi_{\mathrm{cat}}(x)=N\,(e^{-(x-x_0)^2/2}+e^{-(x+x_0)^2/2})$ with $N$ fixed by `Integrate`
($N^{-2}=2\sqrt{\pi}\,(1+e^{-x_0^2})$), and push it through the Wigner kernel of the common ground
beside the coherent state $|\alpha\rangle$ (a displaced vacuum Gaussian, positive everywhere) and
the $n=1$ number state (the negative disk): one `Integrate` under reality assumptions produces the
two blobs at $(\pm x_0,0)$ plus the interference term $\propto e^{-x^2-p^2}\cos(2x_0p)$ at the
origin, fringes of period $\pi/x_0$ in $p$ set by the blob separation $2x_0$. The discriminating
contrast is the incoherent mixture
$\tfrac{1}{2}(|\alpha\rangle\langle\alpha|+|{-\alpha}\rangle\langle{-\alpha}|)$: its kernel through
the same transform must give exactly the two blobs with no fringe term, so the fringe visibility
(interference amplitude over blob amplitude at the origin, unity for the pure cat up to
$e^{-2\alpha^2}$ normalization corrections) is the check that separates coherence from classical
ignorance, with the marginal and normalization checks of the transform repeated on the cat. Close
in the lab: the momentum marginal $\propto e^{-p^2}\cos^{2}(x_0p)$ shows the fringes to a momentum
measurement while the position histogram sees two dead bumps, and the fringe period shrinking as
$1/x_0$ is the seed of why large cats are the first to decohere.

### 16.5 [MSc] How do I compute the Husimi-$Q$ and Glauber-Sudarshan-$P$ representations?

Run the same trio, the coherent state $\alpha_0=2$, the $n=1$ number state, and the even cat
$\alpha=2$ (each redefined here), through both representations. For Husimi, take
$Q(\alpha)=\tfrac{1}{\pi}|\langle\alpha|\psi\rangle|^2$ with the coherent-state wavefunction
$\phi_{\alpha}(x)=\pi^{-1/4}\,e^{-(x-\sqrt{2}\operatorname{Re}\alpha)^2/2+i\sqrt{2}\operatorname{Im}\alpha\,x}$
(phase convention fixed in-entry) and overlaps by `Integrate` under real
$\operatorname{Re}\alpha,\operatorname{Im}\alpha$: the coherent state gives
$\tfrac{1}{\pi}e^{-|\alpha-\alpha_0|^2}$, the number state $\tfrac{|\alpha|^2}{\pi}e^{-|\alpha|^2}$,
the cat two bumps with exponentially suppressed interference, all manifestly nonnegative. The check
that could refute: $Q$ is the Wigner function smoothed by the vacuum Gaussian, so convolving each
state's independently computed $W$ with $\tfrac{1}{\pi}e^{-x^2-p^2}$ by `Integrate` must reproduce
each $Q$ exactly (convolution normalization to be pinned, verify at authoring). For $P$: the
coherent state is $P=\delta^{(2)}(\alpha-\alpha_0)$, exercised only through `Integrate` against
`DiracDelta` (never `NIntegrate`, which silently returns `0.` on a point measure); the thermal
state at mean occupation $\bar{n}$ joins as the regular anchor,
$P=\tfrac{1}{\pi\bar{n}}\,e^{-|\alpha|^2/\bar{n}}$, a true density; the number and cat states admit
no density, so exhibit them honestly through a regularized representation, the $n=1$ case as the
delta-derivative form $e^{|\alpha|^2}\partial_{\alpha}\partial_{\bar{\alpha}}\delta^{(2)}(\alpha)$
(verify at authoring), stating that the limit exists only as a distribution. Verify all of them by
optical equivalence: `Integrate` each $P$ against $|\alpha|^2$ must return the independently
computed $\langle a^{\dagger}a\rangle$ ($|\alpha_0|^2$, $\bar{n}$, $1$, and the cat's closed form).
Close on the ladder: one Gaussian smoothing separates $P$ from $W$ and another separates $W$ from
$Q$, so singular $P$, negative $W$, positive $Q$ are three verdicts on the same state, and
"classical" means the state whose $P$ already sits at the top rung as a density. Concern: the
ledger's $P$ triple (delta, Gaussian, singular) does not map onto the pinned trio, since the number
state's $P$ is a delta derivative rather than a Gaussian; the Gaussian $P$ belongs to the thermal
state, kept here as the regular anchor beside the trio.

### 16.6 [MSc] How do I build the thermal (Gibbs) oscillator state, compute its Wigner function, and evolve a Wigner function by the Moyal equation?

State the Mehler-form Gibbs kernel
$\rho_{\beta}(x,x')=\sqrt{\tfrac{\tanh(\beta/2)}{\pi}}\,\exp\left(-\tfrac{\tanh(\beta/2)}{4}(x+x')^2-\tfrac{\coth(\beta/2)}{4}(x-x')^2\right)$
with $Z=1/(2\sinh(\beta/2))$ and earn it rather than assert it: the Bloch-equation residual
$-\partial_{\beta}\rho_{\beta}=(-\tfrac{1}{2}\partial_x^2+\tfrac{1}{2}x^2)\rho_{\beta}$ must
`FullSimplify` to 0 under $\beta>0$, with unit trace by `Integrate` (whether `Sum` closes the
Fock-ladder Mehler sum symbolically: verify at authoring). Displace by $\alpha$: the kernel picks
up the shift $x\to x-x_0$, $x'\to x'-x_0$ and the phase $e^{ip_0(x-x')}$ with
$(x_0,p_0)=(\sqrt{2}\operatorname{Re}\alpha,\sqrt{2}\operatorname{Im}\alpha)$. One Wigner
`Integrate` under $\beta>0$ then gives the displaced Gaussian
$W(x,p)=\tfrac{1}{\pi\coth(\beta/2)}\,\exp\left(-\tfrac{(x-x_0)^2+(p-p_0)^2}{\coth(\beta/2)}\right)$
of width $\tfrac{1}{2}\coth(\beta/2)$: `Limit` at $\beta\to\infty$ recovers the coherent state's
vacuum width $\tfrac{1}{2}$, a small-$\beta$ `Series` the classical equipartition width $1/\beta$,
and the refuting moment check
$\iint\tfrac{1}{2}(x^2+p^2)\,W\,dx\,dp=\tfrac{1}{2}\coth(\beta/2)+\tfrac{1}{2}(x_0^2+p_0^2)$ must
match the independently known $\langle H\rangle$, reusing the kernel. For the evolution write the
Moyal equation through its first correction,
$\partial_{t}W=-p\,\partial_{x}W+V'(x)\,\partial_{p}W-\tfrac{1}{24}V'''(x)\,\partial_{p}^{3}W$: for
$V=\tfrac{1}{2}x^2$ the correction vanishes identically, so the flow is exactly classical
Liouville, and substituting the rigidly rotating Gaussian
$W_0(x\cos t-p\sin t,\;p\cos t+x\sin t)$ must `FullSimplify` the residual to 0; give the check
teeth by repeating the substitution with $V=x^4/4$, where the retained $V'''$ term leaves a
visibly nonzero residual. Close: a displaced thermal blob orbits the trap rigidly at the
oscillator frequency with no spreading at any temperature; what shears a phase-space distribution
is anharmonicity, never thermal width, and that is the exact content of Moyal reducing to
Liouville for quadratic $H$.
