# Part 22 Plan: Relativistic wave equations

Questions: 5 (22.1 through 22.5). Class census per `Route-Table.md`: C0 (no differential
equation) for 22.1, 22.2, 22.3; C8 (coupled singular radial system, sole member) for 22.4, per
the full C8 verdict as amended by revision R1 in its log; C4 (fixed-energy scattering BVP) for
22.5, per the C4 verdict and its members-sanity note.

## Common ground

The part works in $\hbar=c=1$ with the mass $m$ kept explicit, because the rest energy $m$ is the
protagonist: every pathology and every limit in these five questions happens on the scale
$E\sim m$. The Klein-Gordon equation $(\partial_t^2-\partial_x^2+m^2)\phi=0$ quantizes the
relativistic dispersion $E^2=k^2+m^2$, so both roots $E=\pm\sqrt{k^2+m^2}$ solve it, and its
conserved density $\rho=\tfrac{i}{2m}(\bar\phi\,\partial_t\phi-\phi\,\partial_t\bar\phi)$ is not
sign-definite. The Dirac equation $(i\gamma^\mu\partial_\mu-m)\psi=0$ rests on the Clifford
algebra $\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}$ with
$\eta=\operatorname{diag}(1,-1,-1,-1)$; it buys the positive density $\rho=\psi^\dagger\psi$ and
current $j^\mu=\bar\psi\gamma^\mu\psi$ ($\bar\psi=\psi^\dagger\gamma^0$) at the price of keeping
the negative-energy branch, which returns as the Klein paradox and zitterbewegung. Its exact
hydrogen spectrum $E_{nj}$ is the fine structure that 14.2 obtains perturbatively. The notation
is crowded, so code names carry the distinctions the TeX cannot: the gamma matrices are
`gmat[mu]` with `alphaD` and `betaD` for the $1+1$ pair, while the fine-structure combination is
`zalpha` ($Z\alpha$) and the Dirac-Coulomb exponent is `gam` ($\gamma$); energies are `en`
(`E` is Euler's number); the step potential is `HeavisideTheta` while 22.1's mixing angle is
`th`. Standing kernel facts binding all five entries (Route-Table.md, C0 row and C8 verdict):
densities and currents are written as explicit $\phi\bar\phi$ products before differentiation,
never through Abs; FullSimplify closes residuals only under explicit assumptions ($m>0$, reality
of $x,t,k$, $0<Z\alpha<1$); and every DSolve on a coupled or piecewise system is time-boxed,
because the failure modes in this part are silent hangs and silent echoes, not messages.
Concern: PIPELINE fixes $\hbar=m=1$, and this part departs from it, keeping $\hbar=c=1$ with $m$
symbolic, because $m=1$ would erase the rest energy that every question here measures against
($V_0>2m$, $g=2$, $E_{nj}\to m$, the $1/(2m)$ tremor); the departure is deliberate and part-wide,
and the coordinator should confirm it before authoring.

## Per-question entries

### 22.1 [MSc] How do I solve the Klein-Gordon equation, write its plane-wave modes, and confront the negative-energy and indefinite-density problem?

Write both branches $\phi_{\pm}=e^{i(kx\mp\omega t)}$ with $\omega=\sqrt{k^2+m^2}$ and earn them
twice: the residual of $(\partial_t^2-\partial_x^2+m^2)\phi_{\pm}$ FullSimplifies to 0, and Solve
on the dispersion polynomial returns both roots $\pm\omega$, the negative branch arriving
uninvited. The pinned example is a two-mode superposition with negative-energy admixture, a
positive-energy rest mode plus a negative-energy moving mode,
$\phi=\cos\theta\,e^{-imt}+\sin\theta\,e^{i(\sqrt3\,mx+2mt)}$ at $\theta=\pi/6$. Build
$\rho=\tfrac{i}{2m}(\bar\phi\,\partial_t\phi-\phi\,\partial_t\bar\phi)$ and
$j=\tfrac{1}{2im}(\bar\phi\,\partial_x\phi-\phi\,\partial_x\bar\phi)$ with the conjugate written
out and reduced by ComplexExpand under reality of $x,t,\theta$ and $m>0$ (never differentiate an
Abs square), then FullSimplify the continuity $\partial_t\rho+\partial_x j$ to 0, a check any
sign slip refutes. Compute the indefiniteness instead of narrating it: $\rho$ reduces to
$\cos^2\theta-2\sin^2\theta-\sin\theta\cos\theta\,\cos\Theta$ with travelling phase
$\Theta=3mt+\sqrt3\,mx$, so Minimize over one period gives the closed-form minimum
$(1-\sqrt3)/4<0$ at the pinned angle, and one Solve of
$\cos^2\theta-2\sin^2\theta-\sin\theta\cos\theta=0$ widens that single case into the general
threshold $\tan\theta_c=\tfrac12$, above which $\rho$ goes negative somewhere, while each branch
alone gives the constant $\pm\omega/m$. Two plane waves are not square integrable, so this is a
pointwise statement; the normalizable version, worth one Integrate on a mixed-branch $L^2$
profile, is that the total $\int\rho\,dx$ itself can be negative. Anchor with the nonrelativistic
limit: substitute $\phi=e^{-imt}\psi(x,t)$ with generic $\psi$ and Series in $1/m$ to recover
$\rho=|\psi|^2+O(1/m)$, the Schrodinger density. Close on the reinterpretation: $\rho$ survives
as a charge density, and the sign it refuses to fix belongs to the antiparticle.

### 22.2 [MSc] How do I build the Dirac equation from gamma matrices and write the free spinor solutions?

Construct the Dirac-basis gammas as explicit $4\times4$ matrices, $\gamma^0=\sigma_3\otimes\mathbb 1$
and $\gamma^i=i\sigma_2\otimes\sigma_i$ via KroneckerProduct and PauliMatrix with
$\eta=$ DiagonalMatrix of $\{1,-1,-1,-1\}$, and verify the Clifford algebra exhaustively: Table
over all sixteen $(\mu,\nu)$ pairs of $\gamma^\mu\gamma^\nu+\gamma^\nu\gamma^\mu-2\eta^{\mu\nu}\mathbb 1$,
Flatten and Union collapsing to $\{0\}$, an exact matrix computation with nothing sampled. The
pinned example keeps the momentum nontrivial: $\vec p=p\,(\sin\vartheta,0,\cos\vartheta)$
symbolic in the $x$-$z$ plane, helicity two-spinors $\chi_{\pm}$ in half-angle form checked as
eigenvectors of $\vec\sigma\cdot\hat p$, then
$u_{\pm}=\sqrt{E+m}\,\bigl(\chi_{\pm},\pm\tfrac{p}{E+m}\chi_{\pm}\bigr)^T$ and $v_{\pm}$ with the
blocks swapped, all simplification under $p>0$, $m>0$, $0<\vartheta<\pi$ with
$E=\sqrt{p^2+m^2}$ (`en` in code). Verify per spinor, reusing the defined gammas: the Dirac
residuals $(\gamma^\mu p_\mu-m)u_{\pm}=0$ and $(\gamma^\mu p_\mu+m)v_{\pm}=0$ by FullSimplify,
the normalizations $\bar uu=2m$ and $\bar vv=-2m$ with $\bar u=u^\dagger\gamma^0$ via
ConjugateTranspose, and helicity $\vec\Sigma\cdot\hat p\,u_{\pm}=\pm u_{\pm}$ with
$\Sigma_i=\mathbb 1\otimes\sigma_i$; each equation fails under any wrong sign or normalization.
Take both limits with Series on those same expressions: at $p\ll m$, $u_{\pm}\to\sqrt{2m}
(\chi_{\pm},0)^T$ with the lower block $O(p/m)$, which is precisely the small component 22.3
eliminates, and at $m\to0$ the normalization $\bar uu=2m$ degenerates to 0 while
$\vec\Sigma\cdot\hat p$ and $\gamma^5$ agree on $u_{\pm}$, helicity collapsing onto chirality.
The genuine edge is $p\to0$ exactly, where $\hat p$ is undefined and with it the helicity
spinors: at rest only the spin projection survives, and the $\pm$ labels stop meaning anything.
Close with completeness: $\sum_s u_s\bar u_s-(\gamma^\mu p_\mu+m)$ FullSimplifies to the zero
matrix, tying the spinors back to the algebra that built them.

### 22.3 [MSc] How do I take the nonrelativistic limit to the Pauli equation and recover the $g=2$ prediction?

Couple the Dirac equation minimally, $\vec p\to\vec\pi=\vec p-q\vec A$ with a fully symbolic
$\vec A(x,y,z)$ and scalar $\Phi$, charge-sign convention tagged (verify at authoring), and let
every operator act on a generic two-spinor test function $f$ through D, so the noncommutativity
of $\vec p$ and $\vec A$ is handled by differentiation, never by symbol shuffling. First earn the
key identity on generic $f$: $(\vec\sigma\cdot\vec\pi)^2f-(\vec\pi^2-q\,\vec\sigma\cdot\vec B)f$
FullSimplifies to $\{0,0\}$ with $\vec B=$ Curl of $\vec A$ and PauliMatrix carrying the algebra.
The identity fixes the coefficient of $\vec\sigma\cdot\vec B$, nothing more; $g$ enters one step
later, when $-\tfrac{q}{2m}\vec\sigma\cdot\vec B$ is read against $-\vec\mu\cdot\vec B$ with
$\vec\mu=g\,\tfrac{q}{2m}\vec S$ and $\vec S=\tfrac12\vec\sigma$, giving $g=2$. Then reduce:
split $\psi=e^{-imt}(\varphi,\chi)^T$, stripping the rest phase before any expansion (Series in
$1/m$ chokes on the oscillatory factor otherwise), solve the lower component as
$\chi=\tfrac{\vec\sigma\cdot\vec\pi}{2m}\varphi+O(1/m^2)$, substitute into the upper equation,
and a Series at $m\to\infty$ lands on
$i\partial_t\varphi=[\tfrac{\vec\pi^2}{2m}-\tfrac{q}{2m}\vec\sigma\cdot\vec B+q\Phi]\varphi$.
Coefficient on that expansion confirms the $\vec\sigma\cdot\vec B$ term, but it and the identity
share the elimination $\chi=\tfrac{\vec\sigma\cdot\vec\pi}{2m}\varphi$, so a factor-of-2 slip
there moves $g$ and both readings still pass; the independent refuter is the exact Dirac Landau
spectrum in uniform $\vec B=B\hat z$, obtained by squaring the Dirac operator with the same
$\vec A$ and Solve, $E^2=m^2+p_z^2+qB(2n+1-\sigma)$ with $\sigma=\pm1$, whose Series at
$m\to\infty$ returns $E\approx m+\tfrac{p_z^2}{2m}+\tfrac{qB}{2m}(2n+1-\sigma)$ with no
elimination assumed anywhere, and whose $(n,\sigma=+1)$ and $(n+1,\sigma=-1)$ coincidence (the
degeneracy 14.5 exhibits) happens only at $g=2$. Close in the lab: the measured electron value
2.00232 sits a tenth of a percent above the Dirac prediction, and that excess belongs to
radiative corrections outside this equation.

### 22.4 [MSc] How do I obtain the Dirac fine-structure spectrum of hydrogen?

State the coupled radial system $F'=-\tfrac{\kappa}{r}F+(E+m+\tfrac{Z\alpha}{r})G$,
$G'=\tfrac{\kappa}{r}G-(E-m+\tfrac{Z\alpha}{r})F$ (the convention the C8 probes
residual-verified) and never hand it to DSolve: the coupled system hangs even without boundary
conditions, and naive elimination returns only DifferentialRoot (C8 gates; time-box any DSolve at
authoring). Fix the quantum-number dictionary first, because it carries the structure of the
answer: $\kappa=-(l+1)$ for $j=l+\tfrac12$ and $\kappa=+l$ for $j=l-\tfrac12$, hence
$|\kappa|=j+\tfrac12$ and $\gamma=\sqrt{(j+\tfrac12)^2-(Z\alpha)^2}$ (code `gam`, `zalpha`), the
identity from which the spectrum's independence of $l$ follows rather than being asserted. The
probed backbone is the squared Biedenharn route: DSolve the effective Coulomb ODE
$-u''+(\tfrac{s(s+1)}{r^2}-\tfrac{2EZ\alpha}{r})u=(E^2-m^2)u$ with non-integer
$s\in\{\gamma-1,\gamma\}$, which returns WhittakerM with symbolic energy (`en`, and
$0<Z\alpha<1$, $m>0$ explicit in every Simplify); quantization is the manual termination
read-off $EZ\alpha/\sqrt{m^2-E^2}=n_r+s+1$, and Solve gives the exact
$E_{nj}=m\,[1+(\tfrac{Z\alpha}{n_r+\gamma})^2]^{-1/2}$ with $n=n_r+|\kappa|$, the $s=\gamma-1$
and $s=\gamma$ branches paired by hand and $n_r=0$ existing only for $\kappa<0$. Certify every
exhibited state by substituting the assembled $(F,G)$ into the first-order system and
FullSimplifying the residual to $\{0,0\}$ (the probes certify general-$n$ closed forms only for
the ground state, so the residual runs per state), and assert the normalization
$\int_0^\infty(F^2+G^2)\,dr=1$ by Gamma integrals. Cross-check against 14.2's perturbative fine
structure: Series of $E_{nj}$ in $Z\alpha$ reproduces
$m\,[1-\tfrac{(Z\alpha)^2}{2n^2}-\tfrac{(Z\alpha)^4}{2n^4}(\tfrac{n}{j+1/2}-\tfrac34)]$ exactly,
remainder scaling as $(Z\alpha)^6$ (probed). The numeric refuter is 32-digit NDSolve shooting on
the coupled system from $r_0=10^{-8}$ with the $E$-independent indicial ratio
$G/F=(\gamma+\kappa)/(Z\alpha)$, the balance of the $r^{\gamma-1}$ terms, which reduces to the
probed $(\gamma-1)/(Z\alpha)$ only at $\kappa=-1$; seeding a $\kappa=+1$ state with that special
value picks the irregular branch and silently disarms the discriminator, whereas the general
seed blows the far tail up by a factor near $10^7$ for a $10^{-2}$ energy error (probed), with
node counting and $r_{\max}\sim1/\sqrt{m^2-E^2}$ for excited states. Exhibit the degeneracy at
fixed $n$ as $\kappa\to-\kappa$: $2s_{1/2}$ ($\kappa=-1$) and $2p_{1/2}$ ($\kappa=+1$) coincide
exactly, since $E_{nj}$ sees $\kappa$ only through $|\kappa|=j+\tfrac12$, while $2p_{3/2}$
($\kappa=-2$) stays split off by $\alpha^4m/32$ in the same Series. The edges are where the
formula ends: as $Z\alpha\to|\kappa|^{-}$ the exponent $\gamma\to0$ and Limit on
$E_{1,1/2}=m\gamma$ sends the ground level to 0 with divergent slope near $Z\approx137$, the fall
to the center; and at the origin the components go as $r^\gamma$ with $\gamma<1$, finite in
amplitude but with divergent derivative, still normalizable since $2\gamma>-1$, against the
Schrodinger $1s$ reduced radial that leaves the origin linearly as $r$. Close on that cusp and
its missing partner: the relativistic origin behavior is qualitative, and the splitting the
Dirac equation refuses to produce, $2s_{1/2}$ against $2p_{1/2}$, is the Lamb shift.

### 22.5 [MSc] How do I exhibit the Klein paradox and zitterbewegung numerically?

Work in the $1+1$ Dirac representation $\alpha=\sigma_1$, $\beta=\sigma_3$ (code `alphaD`,
`betaD`; representation choice verify at authoring), so spinors are two-component. Klein step:
the pinned supercritical step $V=V_0\,$HeavisideTheta$(x)$ with $V_0>E+m$ (hence $V_0>2m$),
matched region by region because DSolve on a Piecewise coefficient silently echoes (C4 verdict;
the probes covered the Schrodinger case only, so time-box the Dirac matching at authoring): the
equation is first order, so only $\psi$ itself is continuous at $x=0$, a two-component Solve, and
the transmitted momentum $q=\sqrt{(E-V_0)^2-m^2}$, real again in the Klein zone, has its sign
fixed by the group velocity $dE/dq$ (verify at authoring). Transmission comes from the current
ratio $j=\psi^\dagger\alpha\psi$, never from $|t|^2$ (the C4 members-sanity note); FullSimplify
$R+T-1$ to 0 from the currents, then sweep $T(V_0)$ at fixed $E$ across all three windows:
ordinary transmission for $V_0<E-m$, a closed window of width exactly $2m$ for
$E-m<V_0<E+m$ where $q$ is imaginary and $T=0$, and the Klein zone $V_0>E+m$ where $T$ reopens
and stays finite however high the step is raised. The same sweep on the Schrodinger step decays
monotonically and never reopens, so the discriminator is the reopening itself, and the closed
window being exactly $2m$ wide is the pair-creation threshold in disguise. Zitterbewegung is a
C0 mode sum (Route-Table C0 row), no PDE solver: pin
$\psi_0(x)=\operatorname{sech}^2(x/a)\,e^{ik_0x}\,(1,0)^T$ with $k_0=m/2$ and $a=4/m$, the flat
spinor chosen because it is not a $\pm$ eigenvector, and report $\|\Lambda_{-}(k_0)\psi_0\|^2$
as the number that sets the tremor amplitude (of order $1/(2m)$, half the Compton wavelength);
starting from $u_{+}(k_0)$ instead would send it to zero and leave the checks with nothing to
measure. The momentum amplitudes stay symbolic through
$\int\operatorname{sech}^2(x)e^{-ikx}dx=\pi k/\sinh(\pi k/2)$ by Integrate, rescaled for width
$a$, so the only discretization is the $k$-quadrature (Subdivide); split each mode with the
projectors $\Lambda_{\pm}(k)=\tfrac12(\mathbb 1\pm H(k)/E_k)$, $H(k)=\alpha k+\beta m$, rebuild
$\psi(x,t)$ by Total over the grid, and take $\langle x\rangle(t)$ as a $dx$-weighted grid sum:
a drift with a superposed tremor. Three refuters reuse those same definitions: reconstructing
$\psi(x,0)$ from the discrete mode sum must reproduce the pinned profile and keep improving
under a $k$-grid doubling, or the quadrature is inventing the signal; the tremor period must
match $\pi/E_{k_0}=2.81/m$; and rerunning the identical sum on $\Lambda_{+}$-projected data must
flatten the tremor while leaving the drift. State the validity window with the period check,
since it is the reason $a$ is pinned wide: $a=4/m$ gives $\Delta k\approx m/4$, hence
$\Delta(2E_k)\approx0.22m$ and dephasing near $4.5/m$, about one and a half clean oscillations
(and $a=8/m$ doubles that), while at the literal Compton width $a=1/m$ the tremor dephases in
$1.1/m$ and the period is unmeasurable; the drift over one period, $\approx1.26/m$ at
$k_0=m/2$, stays comparable to the tremor amplitude, so both live on one plot. Close where the
two computations meet: both effects are $\pm E$ interference, and both dissolve only in the
many-particle reinterpretation the supercritical step forces.
