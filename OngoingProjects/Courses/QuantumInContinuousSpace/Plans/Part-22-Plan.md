# Part 22 Plan: Relativistic wave equations

Questions: 5 (22.1 through 22.5). Class census per `Route-Table.md`: C0 (no differential
equation) for 22.1, 22.2, 22.3; C8 (coupled singular radial system, sole member) for 22.4, per
the full C8 verdict; C4 (fixed-energy scattering BVP) for 22.5, per the C4 verdict and its
members-sanity note.

## Common ground

The part works in $\hbar=c=1$ with the mass $m$ kept explicit, because the rest energy $m$ is the
protagonist: every pathology and every limit in these five questions happens on the scale
$E\sim m$ (WL code writes energies as `en`, since `E` is Euler's number). The Klein-Gordon
equation $(\partial_t^2-\partial_x^2+m^2)\phi=0$ quantizes the relativistic dispersion
$E^2=k^2+m^2$, so both roots $E=\pm\sqrt{k^2+m^2}$ solve it, and its conserved density
$\rho=\tfrac{i}{2m}(\bar\phi\,\partial_t\phi-\phi\,\partial_t\bar\phi)$ is not sign-definite. The
Dirac equation $(i\gamma^\mu\partial_\mu-m)\psi=0$ rests on the Clifford algebra
$\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}$ with $\eta=\operatorname{diag}(1,-1,-1,-1)$; it buys
the positive density $\rho=\psi^\dagger\psi$ and current $j^\mu=\bar\psi\gamma^\mu\psi$
($\bar\psi=\psi^\dagger\gamma^0$) at the price of keeping the negative-energy branch, which
returns as the Klein paradox and zitterbewegung. Its exact hydrogen spectrum $E_{nj}$, a function
of $n$ and $j$ only, is the fine structure that 14.2 obtains perturbatively. Standing kernel
facts binding all five entries (Route-Table.md, C0 row and C8 verdict): densities and currents
are written as explicit $\phi\bar\phi$ products before differentiation, never through Abs;
FullSimplify closes residuals only under explicit assumptions ($m>0$, reality of $x,t,k$,
$0<Z\alpha<1$); and every DSolve on a coupled or piecewise system is time-boxed, because the
failure modes in this part are silent hangs and silent echoes, not messages.

## Per-question entries

### 22.1 [MSc] How do I solve the Klein-Gordon equation, write its plane-wave modes, and confront the negative-energy and indefinite-density problem?

Write both branches $\phi_{\pm}=e^{i(kx\mp\omega t)}$ with $\omega=\sqrt{k^2+m^2}$ and earn them
twice: the residual of $(\partial_t^2-\partial_x^2+m^2)\phi_{\pm}$ FullSimplifies to 0, and Solve
on the dispersion polynomial returns both roots $\pm\omega$, the negative branch arriving
uninvited. The pinned example is a packet with negative-energy admixture, a rest mode plus a
moving negative-energy mode, $\phi=\cos\theta\,e^{-imt}+\sin\theta\,e^{i(\sqrt3\,mx+2mt)}$ at
$\theta=\pi/6$. Build $\rho=\tfrac{i}{2m}(\bar\phi\,\partial_t\phi-\phi\,\partial_t\bar\phi)$ and
$j=\tfrac{1}{2im}(\bar\phi\,\partial_x\phi-\phi\,\partial_x\bar\phi)$ with the conjugate written
out and reduced by ComplexExpand under reality of $x,t,\theta$ and $m>0$ (never differentiate an
Abs square), then FullSimplify the continuity $\partial_t\rho+\partial_x j$ to 0, a check any
sign slip refutes. Compute the indefiniteness instead of narrating it: $\rho$ reduces to
$\cos^2\theta-2\sin^2\theta-\sin\theta\cos\theta\,\cos\Theta$ with a travelling interference
phase $\Theta(x,t)$, and Minimize over one period (or Reduce on $\rho<0$) exhibits the
closed-form minimum $(1-\sqrt3)/4<0$, while each branch alone gives the constant $\pm\omega/m$:
positivity is exactly what the admixture destroys. Anchor with the nonrelativistic limit:
substitute $\phi=e^{-imt}\psi(x,t)$ with generic $\psi$ and Series in $1/m$ to recover
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
Close with completeness: $\sum_s u_s\bar u_s-(\gamma^\mu p_\mu+m)$ FullSimplifies to the zero
matrix, tying the spinors back to the algebra that built them.

### 22.3 [MSc] How do I take the nonrelativistic limit to the Pauli equation and recover the $g=2$ prediction?

Couple the Dirac equation minimally, $\vec p\to\vec\pi=\vec p-q\vec A$ with a fully symbolic
$\vec A(x,y,z)$ and scalar $\Phi$, charge-sign convention tagged (verify at authoring), and let
every operator act on a generic two-spinor test function $f$ through D, so the noncommutativity
of $\vec p$ and $\vec A$ is handled by differentiation, never by symbol shuffling. First earn the
key identity on generic $f$: $(\vec\sigma\cdot\vec\pi)^2f-(\vec\pi^2-q\,\vec\sigma\cdot\vec B)f$
FullSimplifies to $\{0,0\}$ with $\vec B=$ Curl of $\vec A$ and PauliMatrix carrying the algebra;
this identity is where $g=2$ lives, and a $g=1$ guess fails it. Then reduce: split
$\psi=e^{-imt}(\varphi,\chi)^T$, stripping the rest phase before any expansion (Series in $1/m$
chokes on the oscillatory factor otherwise), solve the lower component as
$\chi=\tfrac{\vec\sigma\cdot\vec\pi}{2m}\varphi+O(1/m^2)$, substitute into the upper equation,
and a Series at $m\to\infty$ lands on the Pauli equation
$i\partial_t\varphi=[\tfrac{\vec\pi^2}{2m}-\tfrac{q}{2m}\vec\sigma\cdot\vec B+q\Phi]\varphi$. The
refuting check is Coefficient on the expansion, reusing the same definitions: the
$\vec\sigma\cdot\vec B$ coefficient must equal $-q/(2m)$, the value the identity fixes, which is
the moment $g\,\tfrac{q}{2m}\vec S$ with $g=2$ exactly. Close in the lab: the measured electron
value 2.00232 sits a tenth of a percent above the Dirac prediction, and that excess belongs to
radiative corrections outside this equation.

### 22.4 [MSc] How do I obtain the Dirac fine-structure spectrum of hydrogen?

State the coupled radial system $F'=-\tfrac{\kappa}{r}F+(E+m+\tfrac{Z\alpha}{r})G$,
$G'=\tfrac{\kappa}{r}G-(E-m+\tfrac{Z\alpha}{r})F$ (the convention the C8 probes
residual-verified) and never hand it to DSolve: the coupled system hangs even without boundary
conditions, and naive elimination returns only DifferentialRoot (C8 gates; time-box any DSolve at
authoring). The probed backbone is the squared Biedenharn route: DSolve the effective Coulomb ODE
$-u''+(\tfrac{s(s+1)}{r^2}-\tfrac{2EZ\alpha}{r})u=(E^2-m^2)u$ with non-integer
$s\in\{\gamma-1,\gamma\}$, $\gamma=\sqrt{\kappa^2-Z^2\alpha^2}$, which returns WhittakerM with
symbolic energy (`en`, and $0<Z\alpha<1$, $m>0$ explicit in every Simplify); quantization is the
manual termination read-off $EZ\alpha/\sqrt{m^2-E^2}=n_r+s+1$, and Solve gives the exact
$E_{nj}=m\,[1+(\tfrac{Z\alpha}{n_r+\gamma})^2]^{-1/2}$ with $n=n_r+|\kappa|$, the $s=\gamma-1$
and $s=\gamma$ branches paired by hand and $n_r=0$ existing only for $\kappa<0$. Certify every
exhibited state by substituting the assembled $(F,G)$ into the first-order system and
FullSimplifying the residual to $\{0,0\}$ (the probes certify general-$n$ closed forms only for
the ground state, so the residual runs per state), and assert the normalization
$\int_0^\infty(F^2+G^2)\,dr=1$ by Gamma integrals. Cross-check against 14.2's perturbative fine
structure: Series of $E_{nj}$ in $Z\alpha$ reproduces
$m\,[1-\tfrac{(Z\alpha)^2}{2n^2}-\tfrac{(Z\alpha)^4}{2n^4}(\tfrac{n}{j+1/2}-\tfrac34)]$ exactly,
remainder scaling as $(Z\alpha)^6$ (probed); the numeric refuter is 32-digit NDSolve shooting on
the coupled system from $r_0=10^{-8}$ with the $E$-independent indicial ratio
$G/F=(\gamma-1)/(Z\alpha)$, where a $10^{-2}$ energy error blows the far tail up by a factor near
$10^7$ (probed). Exhibit the $\kappa\to-\kappa-1$ degeneracy: $2s_{1/2}$ and $2p_{1/2}$ coincide
exactly because $E_{nj}$ depends only on $n$ and $j$, and close there: the splitting the Dirac
equation refuses to produce is the Lamb shift, the exit toward field theory.

### 22.5 [MSc] How do I exhibit the Klein paradox and zitterbewegung numerically?

Work in the 1+1 Dirac representation $\alpha=\sigma_1$, $\beta=\sigma_3$ (representation choice
verify at authoring), so spinors are two-component. Klein step: the pinned supercritical step
$V=V_0\,\theta(x)$ with $V_0>E+m$ (hence $V_0>2m$), matched region by region because DSolve on a
Piecewise coefficient silently echoes (C4 verdict; the probes covered the Schrodinger case only,
so time-box the Dirac matching at authoring): the equation is first order, so only $\psi$ itself
is continuous at $x=0$, a two-component Solve, and the transmitted momentum
$q=\sqrt{(E-V_0)^2-m^2}$, real again in the Klein zone, has its sign fixed by the group-velocity
direction (verify at authoring). Transmission comes from the current ratio
$j=\psi^\dagger\alpha\psi$, never from $|t|^2$ (the C4 members-sanity note); FullSimplify
$R+T-1$ to 0 from the currents, and Limit as $V_0\to\infty$ leaves $T$ finite and positive, the
paradox itself, against the Schrodinger step where the same limit kills transmission.
Zitterbewegung is a C0-style mode sum, no PDE solver: boost a $\operatorname{sech}^2x$ profile to
$k_0=\tfrac12$ (in units of $m$), decompose it on a symmetric momentum grid (Subdivide), split
every mode into the $\pm E_k$ branches with the projectors
$\Lambda_{\pm}(k)=\tfrac12(\mathbb 1\pm H(k)/E_k)$, $H(k)=\alpha k+\beta m$, rebuild $\psi(x,t)$
by Total over the grid, and take $\langle x\rangle(t)$ as a $dx$-weighted grid sum: a drift with
a superposed tremor. Two refuters reuse the same definitions: the tremor period must match
$\pi/E_{k_0}$ (interference at $2E_k\approx2m$), and rerunning the identical sum on the
$\Lambda_{+}$-projected data must kill the tremor while keeping the drift. Close where the two
computations meet: both effects are $\pm E$ interference, and both dissolve only in the
many-particle reinterpretation the supercritical step forces.
