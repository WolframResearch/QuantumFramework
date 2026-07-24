# Part 9 Plan: Periodic potentials and band structure

Part 9, "Periodic potentials and band structure", 6 questions. Class census (`Route-Table.md`):
C1: 9.1; C4: 9.2, 9.3; C0: 9.4, 9.5; C3: 9.6.

A potential with period $a$, $V(x+a)=V(x)$, commutes with the lattice translation, so energy
eigenstates take the Bloch form $\psi_k(x)=e^{ikx}u_k(x)$ with periodic $u_k(x+a)=u_k(x)$ and
quasimomentum $k$ in the Brillouin zone $(-\pi/a,\pi/a]$; the spectrum organizes into bands
$E_n(k)$ separated by gaps. For a piecewise cell the one-cell transfer matrix $M(E)$ with
$\det M=1$ carries the whole spectrum through $\cos ka=\tfrac12\operatorname{Tr}M(E)$, allowed
bands being $\lvert\operatorname{Tr}M\rvert\le2$. The density of states per site is
$g(E)=\frac{a}{\pi}\lvert dE/dk\rvert^{-1}$ with $\int_{\mathrm{band}}g\,dE=1$ and van Hove
$1/\sqrt{\lvert E-E_{\mathrm{edge}}\rvert}$ edges; Wannier functions
$w_m(x)=\frac{a}{2\pi}\int_{\mathrm{BZ}}e^{-ikma}\psi_k(x)\,dk$ localize when the Bloch gauge is
smooth; a deep lattice gives the tight-binding band $E(k)=\varepsilon-2J\cos ka$, and a static
tilt $Fx$ replaces transport by the Wannier-Stark ladder with spacing and Bloch frequency
$\omega_B=Fa$.

### 9.1 [BSc] How do I state and use Bloch's theorem, writing eigenstates as $e^{ikx}u_k(x)$ with periodic $u_k$?

Pinned example: the cosine lattice $V(x)=V_0\cos2x$ (period $\pi$, the optical-lattice standard)
at $V_0=1$; the substitution $\psi''+(2E-2V_0\cos2x)\psi=0$ maps it onto the Mathieu equation
with $a=2E$, $q=V_0$, so the entry rides the probed C1 Mathieu gate: band energies
$E_n(k)=\tfrac12\,$`MathieuCharacteristicA`/`MathieuCharacteristicB` (certified against a
periodic-cell `NDEigenvalues` solve with `PeriodicBoundaryCondition` at the $3\times10^{-6}$
level) and Bloch states from the `MathieuC`/`MathieuS` Floquet pair (certified against direct
`NDSolve` integration near $10^{-11}$). State the theorem, then verify the form instead of
asserting it: assemble $\psi_k$ at an interior $k$ (say $k=0.3$ in zone-edge units), extract
$u_k(x)=e^{-ikx}\psi_k(x)$, and check $u_k(x+\pi)-u_k(x)$ at sample points, a check that fails
outright if the Floquet assembly is wrong. Exact anchor and refuting limit reusing the same band
function: at $q=0$ the characteristic reduces to $\nu^2$, so `FullSimplify`/`Limit` must return
the free dispersion $E(k)=k^2/2$ exactly. Close at the zone edge: a small-$q$ `Series` of the
first characteristic gap gives $\Delta E\approx V_0$, the nearly-free-electron gap opening
linear in lattice depth.

### 9.2 [MSc] How do I solve the Kronig-Penney model and plot the allowed and forbidden energy bands?

Pinned example: the rectangular lattice, wells of width $a=1$ between barriers of height $V_0=8$
and width $b=1/4$, period $d=a+b$. Work region by region with $e^{\pm iqx}$ in wells and
$e^{\pm\kappa x}$ under barriers ($q=\sqrt{2E}$, $\kappa=\sqrt{2(V_0-E)}$ substituted last),
because `DSolve` on the `Piecewise` potential silently echoes unevaluated (C4 trap): build the
one-cell matrix $M(E)$ by `Dot` from the probed interface matrices, confirm `Det`$[M]=1$, and
read off $\cos kd=\tfrac12\operatorname{Tr}M(E)$, simplified with `ComplexExpand` and
`FullSimplify` under $q,\kappa,a,b>0$. The Bloch trace step is flagged unprobed algebra in the
C4 verdict, so certify it independently: the identical condition must emerge from the direct
continuity-plus-Bloch-phase `Solve` on one cell, and the `FullSimplify` of the difference to $0$
can refute either derivation. `Plot` the trace function with the $[-1,1]$ strip marking allowed
bands, invert to the band diagram $E_n(k)$ with `FindRoot` on a $k$-grid, and pin the first two
bands and their gap from the edge roots $\lvert\tfrac12\operatorname{Tr}M\rvert=1$ (the
$E>V_0$ branch via the continuation $\kappa\to iq'$). Refuting limits reusing the trace:
$V_0\to0$ collapses onto the folded free parabola $E=k^2/2$, and $\kappa b\gg1$ flattens each
band onto an isolated-well level with width $\sim e^{-2\kappa b}$. Close with that opaque-barrier
scaling: bandwidth is tunneling.

### 9.3 [MSc] How do I treat a Dirac comb (a lattice of delta potentials) and find its band edges?

Pinned example: the attractive Dirac comb $V(x)=-g\sum_{n}\delta(x-nd)$ with $g=1$, $d=2$. The
delta never touches numerics (`NIntegrate` silently returns $0.$ on a point measure, and `D` of
`Abs` yields no `DiracDelta`, both binding here per C4): it enters only through the symbolic
jump condition $\psi'(nd^{+})-\psi'(nd^{-})=-2g\,\psi(nd)$, derived once by `Integrate` across
the point. Plane waves between sites plus continuity, jump, and Bloch phase go through `Solve`
to the band condition $\cos kd=\cos qd-\frac{g}{q}\sin qd$ for $E=q^2/2>0$, continued by
$q\to i\kappa$ to $\cos kd=\cosh\kappa d-\frac{g}{\kappa}\sinh\kappa d$ for the negative band;
band edges are the `FindRoot`/`Reduce` roots at $\cos kd=\pm1$. Cross-checks that could refute:
the delta-interface transfer matrix (`Dot` per cell, `Det` $=1$) must reproduce the same trace,
and the ledger's connecting limit, a rectangular-cell trace re-derived inside this entry for
attractive cells and collapsed by `Limit` $b\to0$ at fixed area $V_0b=-g$, must land exactly on
the comb condition. Anchor: as $d\to\infty$ the negative band pinches onto the isolated
delta-well level $E=-g^2/2$ with exponentially small width. Close there: a band is bound-state
hybridization, its width the tunneling overlap of neighboring wells.

### 9.4 [MSc] How do I compute the density of states within a band?

Pinned example: the lowest band of the cosine lattice $V_0\cos2x$ at $V_0=1$, dispersion
$E(k)=\tfrac12\,$`MathieuCharacteristicA`$[k,V_0]$ for $k\in(0,1)$ in zone-edge units, a C0
entry riding already-probed Mathieu evaluations with no new equation solved. The density of
states per site is $g(E)=\frac{a}{\pi}\lvert dE/dk\rvert^{-1}$, normalized so
$\int_{\mathrm{band}}g\,dE=1$; the symbolic $\nu$-derivative of `MathieuCharacteristicA` is
uncertain (verify at authoring), with the fallback a fine `Table` $k$-grid and centered
`Differences`, keeping energies named `en` since `E` is Euler's number (C0 row). Van Hove
content: `Series` of $E(k)$ at $k=0$ and $k=1$ is quadratic in $k-k_{\mathrm{edge}}$, so
$g\sim1/\sqrt{\lvert E-E_{\mathrm{edge}}\rvert}$; exhibit the $-\tfrac12$ slope of $\log g$
against $\log\lvert E-E_{\mathrm{edge}}\rvert$ near both edges. Verification by state counting
reusing the dispersion, able to refute a wrong $g$: sample $E(k)$ on a uniform $N$-point
$k$-grid, bin with `BinCounts`, and compare bin fractions to $\int g\,dE$ over the same bins
(`NIntegrate`, safe here, no point measure), plus the sum-rule total $\int g\,dE=1$. Close on
the contrast between the divergent edges and the flat band interior: the $1/\sqrt{}$ spikes are
what a tunneling or absorption measurement actually resolves.

### 9.5 [MSc] How do I take the tight-binding limit and build Wannier functions from Bloch states?

Pinned example: the deep cosine lattice $V_0\cos2x$ at $V_0=10$ (cold-atom regime), lowest band.
Assemble Bloch states $\psi_k$ on a uniform $k$-grid from the `MathieuC`/`MathieuS` Floquet pair
at `MathieuCharacteristicA` (C1-probed evaluations), and fix the gauge per $k$, say $u_k(0)$
real and positive, before summing; without a smooth gauge the Wannier sum delocalizes, which is
the one nontrivial choice in the construction. Build
$w_0(x)=\frac{a}{2\pi}\int_{\mathrm{BZ}}\psi_k(x)\,dk$ as a uniform $k$-sum (`Table`, `Total`;
the integrand is periodic in $k$ so the uniform sum converges fast) and $w_1(x)=w_0(x-a)$ by
translation. Extract $J$ two independent ways: from the band,
$J=\big(E(\pi/a)-E(0)\big)/4$ off `MathieuCharacteristicA`/`MathieuCharacteristicB`, with the
residual of $E(k)$ against $\varepsilon-2J\cos ka$ quantifying the neglected next-nearest
hopping; and from the states, $J=-\int w_0\,\hat H\,w_1\,dx$ with
$\hat Hw=-\tfrac12w''+V_0\cos(2x)\,w$ via `NIntegrate` (delta-free, safe). Agreement of the two
$J$ values is the refuting check; disagreement convicts the gauge or the $k$-sum. Show
exponential localization as linear $\log\lvert w_0\rvert$ tails steepening across
$V_0\in\{5,10,15\}$. Close with the experimental knob: $J$ falls roughly exponentially in
$\sqrt{V_0}$, which is how a lattice-depth ramp freezes tunneling in cold-atom experiments.

### 9.6 [MSc] How do I exhibit Bloch oscillations (a Wannier-Stark ladder) under a static field?

Pinned example: the tilted deep lattice in its tight-binding frame,
$H=-J\sum_{n}(\vert n{+}1\rangle\langle n\vert+\mathrm{h.c.})+Fa\sum_{n}n\,\vert n\rangle\langle n\vert$
on $N=101$ sites via `SparseArray` with `Band`, pinned at $J=1$, $Fa=0.4$ so $2J/(Fa)=5$ spans
several sites; the ladder is primary because the C3 verdict flags the continuum run as an open
conflict (the unbounded tilt $Fx$ is incompatible with the periodic pseudospectral grid as
written, and the long-time cost is unprobed). Sorted `Eigenvalues` show constant bulk spacing
$Fa$ to machine precision, the Wannier-Stark ladder itself, with chain-edge deviations as the
finite-size failure edge; the ladder amplitudes are `BesselJ`$[n-m,2J/(Fa)]$ (verify at
authoring). Evolve a broad $k_0=0$ packet with the `MatrixExp` action form
`MatrixExp[-I H t, v]`: $\langle x\rangle(t)=a\sum_{n}n\lvert\psi_n\rvert^2$ oscillates at
$\omega_B=Fa$ with peak-to-peak excursion $4J/F$, the semiclassical
$x(t)=x_0+\frac{2J}{F}\big(\cos(Fat)-1\big)$, and no DC transport, the counterintuitive point.
Refuting checks reusing the evolution: the equally spaced spectrum forces exact periodicity,
$1-\lvert\langle\psi(0)\vert\psi(T_B)\rangle\rvert$ at machine zero for $T_B=2\pi/(Fa)$, and the
measured excursion must equal $4J/F$. Flagged cross-check: one Bloch period of the continuum C3
recipe (`NDSolveValue`, tight goals 10) on the $V_0=10$ lattice tilted by $Fx$ in a large box
with explicit Dirichlet walls sized several times $4J/F$ (a time-dependent gauge sweeping
$k\to k-Ft$ would restore periodicity but is unprobed); norm conservation cannot detect wall
reflection (walls reflect at norm 1), so certify by doubling the box and demanding
$\langle x\rangle(t)$ unchanged, with $dx$-weighted grid sums, never `NIntegrate` on the
oscillatory interpolant. Keep $Fa$ well below the band gap or Landau-Zener leakage breaks the
single-band picture. Close with the experiment: cold atoms in a vertical optical lattice read
the local gravitational force off $\omega_B=Fa$.
