# Part 20 Plan: Three-dimensional scattering theory

Questions: 6 (20.1 through 20.6). Class census: C4 for 20.1, 20.4, 20.5 (per the `Route-Table.md`
C4 verdict, radial members); C0 for 20.2, 20.3 (the C0 row); C1 for 20.6 (the C1 verdict,
continuum member).

## Common ground

Every question lives in the stationary continuum at energy $E=k^2/2$, where the scattering state
obeys the asymptotics $\psi(\vec r)\to e^{ikz}+f(\theta)\,e^{ikr}/r$ and the amplitude carries all
observables through $d\sigma/d\Omega=|f(\theta)|^2$. A central potential decouples into partial
waves, each radial channel ending in one number, the phase shift defined by
$u_l(r)\to\sin(kr-l\pi/2+\delta_l)$, so that
$f(\theta)=\tfrac1k\sum_l(2l+1)e^{i\delta_l}\sin\delta_l\,P_l(\cos\theta)$ and
$\sigma_{tot}=\tfrac{4\pi}{k^2}\sum_l(2l+1)\sin^2\delta_l$; unitarity compresses into the optical
theorem $\sigma_{tot}=\tfrac{4\pi}{k}\operatorname{Im}f(0)$. Weak scatterers admit the Born
amplitude $f_B(q)=-\tfrac2q\int_0^\infty r\,V(r)\sin(qr)\,dr$ at momentum transfer
$q=2k\sin(\theta/2)$, and the Coulomb $1/r$ tail is the exception that never reaches free
asymptotics: its radial waves are the Coulomb wave functions with logarithmically distorted phases.
Standing inheritances: the C4 radial members use the certified fixed-$E$ recipe (NDSolveValue
outward integration, log-derivative matching to SphericalBesselJ/SphericalBesselY, hard-sphere
anchor $\delta_0=-ka$ extracted to $10^{-9}$ in a clean sweep), matching away from nodes of $u$
where $\beta=u'/u$ diverges (or matching the pair $(u,u')$ linearly) and evaluating interpolants
only at explicit numbers, since trap (e) silently echoes symbolic arguments and extrapolates
out-of-domain; the C0 members name their machinery directly and keep the symbol `E` out of code
(`en` for energies); the C1 member leans only on the probed fact that CoulombF/CoulombG evaluate
in WL 15, with the Rutherford phase and normalization conventions unprobed.

## Per-question entries

### 20.1 [MSc] How do I expand a scattering state in partial waves and extract the phase shifts $\delta_l$?

Expand the plane wave as $e^{ikz}=\sum_l(2l+1)\,i^l j_l(kr)\,P_l(\cos\theta)$ and let the potential
act channel by channel; the pinned pair is the hard sphere of radius $a$ as the exact anchor with a
finite well of depth $V_0$ and the same radius beside it as the numeric generalization. For the
hard sphere, $u_l(a)=0$ forces the exterior combination $\cos\delta_l\,j_l-\sin\delta_l\,y_l$ to
vanish at the surface, so Solve and FullSimplify give $\tan\delta_l=j_l(ka)/y_l(ka)$ exactly
(SphericalBesselJ, SphericalBesselY), the probed C4 anchor. For the well, run the certified
fixed-$E$ recipe: NDSolveValue outward from the regular start $u_l\sim r^{l+1}$, form the
log-derivative $\beta$ at a matching radius $r_m$ outside the well, and extract
$\tan\delta_l=\frac{k\,j_l'(kr_m)-\beta\,j_l(kr_m)}{k\,y_l'(kr_m)-\beta\,y_l(kr_m)}$; place $r_m$
away from zeros of the interpolant (the log-derivative diverges at nodes of $u$, and matching the
pair $(u,u')$ linearly is the fallback), and evaluate the interpolant only at explicit numeric
points, since trap (e) bites exactly at computed matching points. Verify per the probe: the
PrecisionGoal sweep $\{4,6,8,10\}$ must drop the hard-sphere error monotonically toward $10^{-9}$
against $\delta_0=-ka$, the extraction must not move when $r_m$ does, and the well's numeric
$\delta_0$ must land on its exact region-matching closed form from Solve on the continuity system
(a disagreement refutes the recipe). Close on the low-energy reading: $\delta_0=-ka$ says the
hard sphere's scattering length is its geometric radius, and $\delta_l$ collapses for $l\gg ka$,
the centrifugal barrier screening the distant waves.

### 20.2 [MSc] How do I assemble the scattering amplitude and the differential cross section?

Recompute the hard-sphere shifts in-entry from $\tan\delta_l=j_l(ka)/y_l(ka)$ (SphericalBesselJ,
SphericalBesselY; each answer is self-contained, nothing imported from 20.1), then assemble
$f(\theta)=\tfrac1k\sum_{l=0}^{l_{\max}}(2l+1)e^{i\delta_l}\sin\delta_l\,P_l(\cos\theta)$ with
LegendreP over a Table of shifts and read off $d\sigma/d\Omega=|f(\theta)|^2$ via ComplexExpand.
Sweep the truncation $l_{\max}$ at $ka\sim1$ until the amplitude stops moving (the shifts die once
$l\gtrsim ka$, and the sweep is the honest convergence exhibit). The refuting check computes
$\sigma_{tot}$ two ways from the same shift table: Integrate of $|f(\theta)|^2$ over the sphere,
which Legendre orthogonality collapses, against the partial-wave sum
$\tfrac{4\pi}{k^2}\sum(2l+1)\sin^2\delta_l$, with the $ka\to0$ end anchored at
$\sigma_{tot}\to4\pi a^2$ (pure $l=0$, four times the geometric cross section). Then push $ka$
high with $l_{\max}$ scaled up and watch the forward peak sharpen while $\sigma_{tot}\to2\pi a^2$;
close on that shadow doubling: the classical $\pi a^2$ is doubled by forward diffraction that a
detector at any finite angle misses.

### 20.3 [MSc] How do I compute a cross section in the Born approximation?

Take the Yukawa potential $V(r)=-g\,e^{-\mu r}/r$ (screened Coulomb, the Born staple) and evaluate
the central-potential Born amplitude in one Integrate:
$f_B(q)=-\tfrac2q\int_0^\infty r\,V(r)\sin(qr)\,dr=\tfrac{2g}{q^2+\mu^2}$ under Assumptions
$\mu>0$, $q>0$, $g>0$ (the assumptions are load-bearing for the clean closed form, per the C0
row). Substitute $q=2k\sin(\theta/2)$ for $d\sigma/d\Omega=|f_B|^2$, and let a second Integrate
close the total Born cross section $\sigma_B=\tfrac{16\pi g^2}{\mu^2(\mu^2+4k^2)}$ over the solid
angle. The refuting check is the screening-off limit: Limit $\mu\to0$ must turn $d\sigma/d\Omega$
into exactly the Rutherford form $(\eta/2k)^2\csc^4(\theta/2)$ with $\eta=g/k$ (any stray factor
refutes the amplitude convention), while $\sigma_B$ itself diverges in the same limit, the forward
signature of an unscreened $1/r$ tail. Close on the bridge to 20.6: for the Coulomb potential the
Born modulus is not merely a first-order estimate, the exact treatment with Coulomb wave functions
returns the same Rutherford cross section.

### 20.4 [MSc] How do I extract the low-energy scattering length and verify the optical theorem?

Define $a_s=-\lim_{k\to0}\delta_0/k$ and work the finite well of radius $a$ as a function of depth
$V_0$, class C4. Get $\delta_0(k)$ at a ladder of small $k$ by the certified fixed-$E$ recipe
(NDSolveValue outward, log-derivative matched to SphericalBesselJ/SphericalBesselY away from
interpolant zeros, or the pair $(u,u')$ matched linearly, arguments numeric per trap (e)), and
beside it the exact route: Solve the region-matching continuity system for $\tan\delta_0$ in
closed form and Limit $k\to0$ to $a_s(V_0)=a\left(1-\tan(\kappa_0 a)/(\kappa_0 a)\right)$ with
$\kappa_0=\sqrt{2V_0}$. The refuting structure is the depth sweep: $a_s$ must diverge exactly at
$\kappa_0 a=(2n+1)\pi/2$, the well's known thresholds for a new $l=0$ bound state, and the numeric
extraction must reproduce each pole location (a shifted pole refutes the small-$k$ extraction).
Verify the optical theorem at a moderate $k$ on the same shift set: assemble
$f(0)=\tfrac1k\sum(2l+1)e^{i\delta_l}\sin\delta_l$ and check
$\sigma_{tot}=\tfrac{4\pi}{k}\operatorname{Im}f(0)$ against the partial-wave sum
$\tfrac{4\pi}{k^2}\sum(2l+1)\sin^2\delta_l$; the identity holds term by term because
$\operatorname{Im}(e^{i\delta_l}\sin\delta_l)=\sin^2\delta_l$, so any residual exposes a wrongly
assembled amplitude rather than roundoff. Close on the ultracold reading: near each threshold
$a_s$ swings through $\pm\infty$ and the low-energy cross section $4\pi a_s^2$ dwarfs the
geometric size, the zero-range universality behind Feshbach-tuned atomic gases.

### 20.5 [MSc] How do I fit a resonance to the Breit-Wigner form and check Levinson's theorem in three dimensions?

Build the well-plus-barrier $V(r)=-V_0$ for $r<a$, $V_b>0$ for $a<r<b$, zero beyond, and choose
$(V_0,V_b,a,b)$ so a single $l=0$ shape resonance sits inside the barrier window $0<E_R<V_b$. The
potential is piecewise constant, so the C4 primary applies region by region (DSolve silently
echoes a Piecewise coefficient): Solve the continuity system for $\tan\delta_0(E)$ in closed form,
tabulate $\delta_0$ on an energy grid through the window with the arctan branch unwrapped so the
phase is continuous, and cross-check a few grid points with the fixed-$E$ radial recipe, matching
away from nodes of $u$ and keeping interpolant arguments numeric (trap (e) bites at computed
matching points). Fit $\delta_0(E)=\delta_{bg}+\arctan\tfrac{\Gamma/2}{E_R-E}$ with
NonlinearModelFit (verify at authoring; flagged unprobed in the C4 verdict), seeded at the
steepest phase rise; the refuting checks are that the fitted $E_R$ must sit on the near-$\pi$ jump
of the tabulated curve, and that a $b$-sweep must shrink $\Gamma$ as the barrier thickens, or the
shape-resonance reading is wrong. Re-run Levinson in three dimensions on the same closed form:
$\delta_0(0)-\delta_0(\infty)=n_b\pi$ with $n_b$ the well's actual bound-state count from its
exact quantization conditions, the resonance contributing its rapid $\pi$ rise at $E_R$ without
changing the count. Close on the lifetime $\tau=1/\Gamma$ growing exponentially with barrier
thickness: the fit parameters carry the decay physics.

### 20.6 [MSc] How do I treat Coulomb scattering (Rutherford) with the Coulomb wave functions?

Take the repulsive Coulomb potential $V=\alpha/r$ with Sommerfeld parameter $\eta=\alpha/k$, class
C1: continuum physics with no quantization, and the verdict certifies exactly that CoulombF and
CoulombG evaluate in WL 15, while every phase-shift and normalization convention beyond that is
unprobed, so each convention-dependent equality below carries (verify at authoring). Exhibit the
$1/r$ obstruction first: the radial waves never reach free asymptotics, CoulombF with arguments
$l$, $\eta$, $\rho=kr$ oscillating as $\sin(\rho-\eta\ln2\rho-l\pi/2+\sigma_l)$ with the Coulomb
phase $\sigma_l=\arg\Gamma(l+1+i\eta)$ computed from Gamma and Arg. The refuting check reuses
those definitions: evaluate CoulombF at large $\rho$ numerically and compare against the predicted
distorted sinusoid built from the computed $\sigma_l$ (a constant phase offset refutes the
convention, which is precisely the unprobed item). Then assemble the closed Coulomb amplitude
$f_C(\theta)=-\tfrac{\eta}{2k}\csc^2(\theta/2)\,e^{-i\eta\ln\sin^2(\theta/2)+2i\sigma_0}$ (verify
at authoring) and close on the exact Rutherford
$d\sigma/d\Omega=|f_C(\theta)|^2=(\eta/2k)^2\csc^4(\theta/2)$: the logarithmic phases drop out of
the modulus, the quantum result coincides with the classical formula at every angle, and the
forward $\csc^4$ divergence is the unscreened long-range tail that 20.3's Yukawa regulator cuts
off.
