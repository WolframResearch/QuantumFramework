---
Template: Default
Title: Caldeira-Leggett High-T Lindblad, Symbolic Exploration
Author: Mads Bahrami
---

# Caldeira-Leggett High-T Lindblad: Symbolic Exploration

Working document. It develops the symbolic (closed-form) track of `exploration-plan.md`: given the paper's main outcome, the master equation of arXiv:2508.14262, which of its features can be extracted in closed form that are physically important and carry useful impact. Each feature below is stated with why it matters, then explored; every symbolic result was reproduced in a kernel. Companions: `caldeira-leggett-high-T-lindblad.md` (Layers 0-1 essay), `related-papers/tex/2508.14262/` (paper source).

## The generator, and what "symbolic" means here

The object of study is the general three-coefficient Dekker generator, the most general Markovian quantum-Brownian-motion master equation:

$$
\mathcal{L}\rho = -\frac{i}{\hbar}[H_S,\rho] - i\frac{\gamma}{\hbar}[X,\{P,\rho\}] - D_{PP}[X,[X,\rho]] - D_{XX}[P,[P,\rho]] - 2D_{XP}[X,[P,\rho]],
$$

with $H_S = P^2/2M + V(X)$ and $\gamma=\eta/M$. Its quadratic structure is inherited from the model (a particle coupled bilinearly to a harmonic bath); reducing to a Markovian, constant-coefficient generator at all is the Born-Markov assumption, the exact dynamics being non-Markovian (Hu-Paz-Zhang for a harmonic potential). The paper's result is one distinguished point in this family,

$$
D_{XX} = \frac{\gamma}{6Mk_BT}, \qquad D_{PP} = \frac{2\gamma Mk_BT}{\hbar^2}, \qquad D_{XP}=0,
$$

pinned by four assumptions: an Ohmic bath, weak coupling $\hbar\gamma\ll k_BT$, the cutoff removed at fixed temperature ($k_BT\ll\hbar\Omega$, $\Omega\to\infty$), and keeping only the leading $\delta''$ correction past white noise. The first three make the coefficients cutoff-free and give $D_{XP}=0$; relax one and they move, a different bath rescales them, the other scale ordering gives Diósi or Breuer-Petruccione with a cutoff-dependent $D_{XP}$, and higher order adds lower-temperature corrections. Everything symbolic below explores this one point.

Every term is at most quadratic in $(X,P)$ except the Hamiltonian commutator $[V(X),\cdot]$: for quadratic $V$ the dynamics closes on the first two moments and everything is closed form, while for anharmonic $V$ the hierarchy opens and simulation takes over. The high-impact symbolic features live in the quadratic (harmonic and free) case, plus one basis-level fact that holds for any $V$.

## The features worth exploring, and why

Two questions must be kept apart: which features close in symbolic form (many do), and which are genuinely new (few are). Take the audit first, then the results that clear the bar.

**The four features below reproduce or reframe. By themselves they do not add to the literature.**

- **A. Thermalization pins $D_{XX}$.** Reading $D_{XX}=\gamma/6Mk_BT$ off the equilibrium momentum variance, with complete positivity following, is a clean *reframing* of the paper's positivity argument, and the sharpest of the four. It is not new physics. That a completely positive Markovian generator cannot sit exactly at the Gibbs state, and that the position-diffusion term restoring positivity is precisely what sustains a steady phase-space current with nonzero entropy production, was established for the general Dekker family by Artini, Lo Monaco, Imparato, Paternostro and Donadi (arXiv:2507.23322); the inverse route, fixing $D_{XX}$ *from* the entropy production, is Bernád, Homa and Csirik (Eur. Phys. J. D **72**, 212 (2018)). Feature A is their statement, specialized to the paper's point, and section A below adds the exact version of the tension.
- **B. Decoherence and classicalization.** Quantifies the paper's own title claim by reading the two double-commutator rates off the generator. Reproduction.
- **C. The Wigner picture.** The textbook Wigner transform of the generator into a Klein-Kramers equation. Reproduction.
- **Engine. The moment hierarchy.** The standard Ehrenfest-plus-diffusion moment ladder from the adjoint generator. Reproduction.

So of the four, only A carries a useful idea, and even that idea's deep half, the thermalization-positivity tension, is already known.

**The bar, and where the frontier actually is.** A referee credits what neither the paper nor the literature already holds. Two directions looked promising, and testing them was the point: the systematic low-temperature extension of the kernel, and a sharp thermalization error against Gibbs. Both turn out to be settled, and settled correctly, in places this document must defer to rather than duplicate. Reporting that plainly is the honest result here.

- **The systematic expansion is done, and it converges.** The companion essay carries the fluctuation-path expansion to all orders (`caldeira-leggett-high-T-lindblad.md`, Parts VI-VII). The $\delta^{(2n)}$ coefficients are $c_{2n}=4\eta(-1)^nB_{2n}\tau_B^{2n-1}/(2n)!$; the even-derivative tower resums to the exact Ohmic spectrum $2\eta\omega\coth(\hbar\omega/2k_BT)$; complete positivity is special to the two-term order, since the $\delta''''$ channel variance is negative ($B_4<0$); the $\delta''''$ operator is the free acceleration $-V'(\hat X)/M$, so it stays in the harmonic (X,X) channel and becomes a new Lindblad operator only for anharmonic $V$; and each order deepens the Gibbs match by one power of $\hbar^2$, the $\delta''''$ term cancelling the $\hbar^4$ miss exactly. Decisively, the series is *convergent*, radius set by the first Matsubara frequency $\hbar\omega=2\pi k_BT$. The raw moments $\pi^{2n}|B_{2n}|$ do grow factorially, but the smooth-path reduction supplies a $1/(2n)!$ that turns them into the decaying coefficients above, so the tower converges down to $T\sim\hbar\omega_0/2\pi k_B$. Any "asymptotic, keep-one-term-is-optimal" reading of this hierarchy is wrong: the $1/(2n)!$ is the whole point, and it lives in the coefficient, not the bare moment.
- **The equilibrium error is a known kind.** That a completely positive Markovian steady state departs from Gibbs beyond leading order in $\hbar$ is established (arXiv:2305.06312), and the broken-detailed-balance side, the steady correlation $\sigma_{XP}=-M\hbar^2D_{XX}$ and its entropy production, is Artini, Lo Monaco, Imparato, Paternostro and Donadi (arXiv:2507.23322) and Bernád, Homa and Csirik (Eur. Phys. J. D **72**, 212 (2018)). Section A below makes the completely-positive-versus-thermalization tension exact for this family, but as an illustration of that known result, not a new one.

So neither promising direction is new. The genuinely open frontier, settled by neither the essay nor the broken-detailed-balance literature, was two questions. The first is now answered, in closed form, in its own section below.

- **Mean force, not bare Gibbs: answered.** A weakly coupled open system equilibrates not to the bare Gibbs state but to the mean-force Gibbs state, whose reduced covariance differs at $O(\gamma)$. Are the completely positive steady state's deviations from bare Gibbs those mean-force corrections in disguise? The section "Mean force, not bare Gibbs" below answers no in two lines (the generator's diagonal has no $O(\gamma)$ deviation at all; its steady correlation sits against a time-reversal zero), then computes what holding the mean-force covariance would take. On the Caldeira-Leggett friction slice, exact stationarity of any diagonal covariance costs $\det a\le-\gamma^2/\hbar^2$ (the Dekker-Valsakumar floor; the slice impossibility is in Isar's thermal-generator formulas), no matter how the Hamiltonian is renormalized; the $O(\gamma)$ mean-force corrections are cutoff-dependent in momentum and odd in $\hbar$ in position, the latter structurally beyond the local kernel expansion's even, diagonal sources (the parity scoping is the one genuinely unasked question here); and the slice-versus-wider-family dichotomy is Artini and co-workers' translation-covariance theorem, rendered as a single closed-form slack whose completely positive equilibrating window $|\mu|\le\lambda\,\mathrm{sech}(\hbar\omega_0/2k_BT)$ shows the obstruction to be exactly critical. The exact Hu-Paz-Zhang dynamics shares the slice signature (no $(\hat P,\hat P)$ term, a time-dependent cross coefficient, mean-force endpoint), which yields one falsifiable prediction, $\det a_{\rm HPZ}^{\,t\to\infty}\le-\gamma_\infty^2/\hbar^2$, posed conventions-pinned in Next steps.
- **The Hu-Paz-Zhang benchmark, quantitatively.** Both the two-term generator and the exact non-Markovian dynamics are Gaussian for harmonic $V$, so their covariance difference along the transient is a closed-form contest between truncation error and memory. The equilibrium endpoints are settled by the mean-force section below (the exact late-time state is the mean-force Gibbs state); what remains open is the transient, which the essay sets up (Part VIII).

Sections A and B carry the established closed-form results, and section A makes the equilibrium tension exact. The mean-force section is the new work. The Wigner picture (C) and the moment hierarchy (Engine) remain the general-$V$ tools, set up at the end.

## A. Thermalization pins $D_{XX}$; positivity is a corollary

Take the harmonic case $V=\tfrac12 M\omega_0^2 X^2$, where the steady state is Gaussian and exact. The physical question is whether the master equation relaxes to the quantum canonical (Gibbs) state $\propto e^{-H_S/k_BT}$, whose covariance is known in closed form (the positivity assumptions are declared once here and every algebraic step below inherits them):

```wl
$Assumptions = {hbar > 0, M > 0, w0 > 0, gam > 0, kB > 0, T > 0, Om > 0,
   aa > 0, bb > 0, lam > 0, Ms > 0, ws > 0};
sXXgibbs = hbar/(2 M w0) Coth[hbar w0/(2 kB T)];
sPPgibbs = hbar M w0/2 Coth[hbar w0/(2 kB T)];
Normal @ Series[sPPgibbs, {hbar, 0, 2}]
```

The high-temperature expansion of the Gibbs momentum variance is $Mk_BT + \hbar^2 M\omega_0^2/(12k_BT) + \dots$: the classical equipartition value plus the leading quantum thermal correction of order $\hbar^2$.

Now the master equation's steady covariance. The Lyapunov equation solves in one stroke for the *full* three-coefficient family, the diffusion matrix carrying the Dekker coefficients as $D = 2\hbar^2\{\{D_{XX},-D_{XP}\},\{-D_{XP},D_{PP}\}\}$:

```wl
Amat = {{0, 1/M}, {-M w0^2, -2 gam}};
diffusion[dxx_, dpp_, dxp_] := 2 hbar^2 {{dxx, -dxp}, {-dxp, dpp}};
sigmaGen = FullSimplify[LyapunovSolve[Amat, -diffusion[dxx, dpp, dxp]],
   Element[{dxx, dpp, dxp}, Reals]];
sigmaGen
```

The general steady covariance of the family in closed form, and its structure is already the physics: $\sigma_{XP}=-\hbar^2MD_{XX}$ depends on the position-diffusion coefficient alone, $\sigma_{PP}$ sees only $(D_{PP},D_{XX})$, and the cross coefficient enters $\sigma_{XX}$ only. It is the genuine stationary solution, residual identically zero:

```wl
FullSimplify[Amat . sigmaGen + sigmaGen . Transpose[Amat] + diffusion[dxx, dpp, dxp]]
```

And the steady state exists and is unique because the drift is Hurwitz:

```wl
FullSimplify @ Eigenvalues[Amat]
```

The damped-oscillator poles $-\gamma\pm\sqrt{\gamma^2-\omega_0^2}$; that both stay in the open left half-plane for all positive parameters is one quantifier elimination via the Routh-Hurwitz form:

```wl
Resolve[ForAll[{M, w0, gam}, M > 0 && w0 > 0 && gam > 0,
   Tr[Amat] < 0 && Det[Amat] > 0], Reals]
```

Relaxation happens on the $1/2\gamma$ timescale, independent of temperature and of $\hbar$. Now specialize to the paper's channel structure, $D_{PP}$ at the Einstein value and $D_{XP}=0$, keeping $D_{XX}$ free so we can see what value the physics wants:

```wl
steady[d_] := steady[d] = Simplify[sigmaGen /. {dxx -> d, dpp -> 2 gam M kB T/hbar^2, dxp -> 0}];
steady[dxx][[2, 2]]
```

The steady momentum variance is $\sigma_{PP} = Mk_BT + D_{XX}\,\hbar^2 M^2\omega_0^2/(2\gamma)$, linear in the new coefficient. Setting it equal to the Gibbs value fixes $D_{XX}$:

```wl
Solve[steady[dxx][[2, 2]] == Normal@Series[sPPgibbs, {hbar, 0, 2}], dxx]
```

The unique solution is $D_{XX} = \gamma/(6Mk_BT)$, exactly the paper's coefficient. So $D_{XX}$ is not a free positivity-restoring parameter: it is the value, and the only value, that makes the equilibrium momentum fluctuation quantum-mechanically correct to order $\hbar^2$.

Complete positivity is then automatic, and it is worth building the Kossakowski matrix as an object rather than quoting its determinant from memory. With $L_1=X$, $L_2=P$ the determinant is $\det a = 4D_{XX}D_{PP}-4D_{XP}^2-\gamma^2/\hbar^2$, and its zero over the whole family:

```wl
kossakowski[dxx_, dpp_, dxp_] := {{2 dpp, 2 dxp - I gam/hbar}, {2 dxp + I gam/hbar, 2 dxx}};
Solve[Det[kossakowski[dxx, 2 gam M kB T/hbar^2, dxp]] == 0, dxx]
```

The positivity floor is $D_{XX}\ge(\gamma^2+4D_{XP}^2\hbar^2)/8\gamma Mk_BT$, minimized at $D_{XP}=0$ where it reads $\gamma/8Mk_BT$: cross diffusion only raises it, so every bound derived below at $D_{XP}=0$ holds a fortiori across the family. Two facts make this matrix trustworthy rather than quoted. First, it is faithful to the generator: the GKSL dissipator built from it with $L_1=X$, $L_2=P$ equals the Dekker double-commutator form once a frequency-renormalization Hamiltonian $\tfrac{\gamma}{2}\{X,P\}$ is split off into $H_{\rm eff}$ (the Lamb shift the paper's footnote absorbs), an operator identity on a truncated Fock carrier with everything symbolic:

```wl
Block[{nn = 6, a, x, p, comm, acomm, rho, aK, ls, gksl, dekker, lamb},
  a = SparseArray[Band[{1, 2}] -> Sqrt[Range[nn - 1]], {nn, nn}];
  x = Sqrt[hbar/(2 M w0)] (a + Transpose[a]); p = I Sqrt[hbar M w0/2] (Transpose[a] - a);
  comm[u_, v_] := u . v - v . u; acomm[u_, v_] := u . v + v . u;
  rho = Array[r, {nn, nn}]; ls = {x, p}; aK = kossakowski[dxx, dpp, dxp];
  gksl = Sum[aK[[i, j]] (ls[[i]] . rho . ls[[j]] - (1/2) acomm[ls[[j]] . ls[[i]], rho]), {i, 2}, {j, 2}];
  dekker = -I (gam/hbar) comm[x, acomm[p, rho]] - dpp comm[x, comm[x, rho]] -
     dxx comm[p, comm[p, rho]] - 2 dxp comm[x, comm[p, rho]];
  lamb = (I gam/(2 hbar)) comm[acomm[x, p], rho];
  Union @ Flatten @ Simplify @ Normal[(gksl - dekker - lamb)[[1 ;; 4, 1 ;; 4]]]]
```

The interior block vanishes for arbitrary symbolic $\hbar$, $M$, $\omega_0$, coefficients, and state, so the Kossakowski matrix is the generator, up to a Hamiltonian shift that touches no positivity statement. Second, the determinant is the binding constraint but not the whole criterion; positive semidefiniteness needs the diagonal too, and at the paper's point the full spectrum clears:

```wl
Resolve[ForAll[{hbar, M, gam, kB, T}, hbar > 0 && M > 0 && gam > 0 && kB > 0 && T > 0,
   And @@ Thread[Eigenvalues[kossakowski[gam/(6 M kB T), 2 gam M kB T/hbar^2, 0]] >= 0]]]
```

Both eigenvalues are nonnegative at every temperature, mass, and coupling. The three coefficients line up as $\{\,\gamma/(6Mk_BT)\ \text{(thermal, this work)},\ \gamma/(8Mk_BT)\ \text{(minimal positivity, Breuer-Petruccione)},\ 0\ \text{(Caldeira-Leggett)}\,\}$. The thermal value exceeds the positivity threshold by a factor $4/3$, so getting the equilibrium fluctuation right guarantees a completely positive generator with margin to spare. Breuer-Petruccione sits exactly on the positivity boundary but undershoots the thermal correction ($\hbar^2 M\omega_0^2/16k_BT$ against the correct $\hbar^2 M\omega_0^2/12k_BT$); Caldeira-Leggett has neither.

How closely the full covariance matches Gibbs, term by term in $\hbar$:

```wl
{Normal@Series[steady[gam/(6 M kB T)][[2, 2]] - sPPgibbs, {hbar, 0, 4}],
 Normal@Series[steady[gam/(6 M kB T)][[1, 1]] - sXXgibbs, {hbar, 0, 4}],
 Normal@Series[steady[gam/(6 M kB T)][[1, 2]], {hbar, 0, 4}]}
```

The momentum variance matches the Gibbs state through order $\hbar^3$, its first surviving difference the $\gamma$-free $\hbar^4$ term now visible. The position variance matches to order $\hbar^2$ with a residual $\gamma^2\hbar^2/(3Mk_BT\omega_0^2)$, and there is a steady position-momentum correlation $-\gamma\hbar^2/(6k_BT)$ that the Gibbs state does not have. The $\hbar^2$ residuals are simultaneously quantum and dissipative, vanishing in the weak-coupling limit $\gamma\to 0$, exactly the regime where a Markovian generator is meant to hold; the $\hbar^4$ pieces are $\gamma$-free truncation residues that survive even that limit. Standard Caldeira-Leggett, by contrast, misses the quantum correction entirely: its steady state is the classical equipartition covariance $\mathrm{diag}(k_BT/M\omega_0^2,\,Mk_BT)$, off from Gibbs by order $\hbar^2$ in both entries.

Whatever its distance from Gibbs, the steady state is a legal quantum state, and not just at the paper's point: for *every* completely positive member of the family the Robertson-Schrödinger bound holds, provable outright as one quantifier elimination over the coefficients themselves:

```wl
Resolve[ForAll[{hbar, M, w0, gam, dxx, dpp, dxp},
   hbar > 0 && M > 0 && w0 > 0 && gam > 0 && dxx >= 0 && dpp >= 0 &&
    4 dxx dpp - 4 dxp^2 - gam^2/hbar^2 >= 0,
   Det[sigmaGen] >= hbar^2/4], Reals]
```

Complete positivity of the generator implies the uncertainty floor for its steady covariance, at every temperature, mass, frequency, coupling, and admissible coefficient triple: the advertised job, done at the Gaussian level, family-wide.

The section's closed forms deserve one reference computed by an entirely different route: the full density-matrix Liouvillian on a truncated Fock space, evolved to stationarity, against the Lyapunov closed form at the same point:

```wl
Block[{nn = 40, a, x, p, id, kp, h, lmat, v0, v, rho, tr, mom, gamN = 1/10},
  a = N @ SparseArray[Band[{1, 2}] -> Sqrt[Range[nn - 1]], {nn, nn}];
  x = (a + Transpose[a])/Sqrt[2.]; p = I (Transpose[a] - a)/Sqrt[2.];
  id = N @ IdentityMatrix[nn, SparseArray]; kp = KroneckerProduct;
  h = p . p/2 + x . x/2;
  lmat = -I (kp[h, id] - kp[id, Transpose[h]]) -
    I gamN (kp[x . p, id] + kp[x, Transpose[p]] - kp[p, Transpose[x]] - kp[id, Transpose[p . x]]) -
    2 gamN (kp[x . x, id] - 2 kp[x, Transpose[x]] + kp[id, Transpose[x . x]]) -
    (gamN/6) (kp[p . p, id] - 2 kp[p, Transpose[p]] + kp[id, Transpose[p . p]]);
  v0 = Flatten @ SparseArray[{{1, 1} -> 1. + 0. I}, {nn, nn}];
  v = MatrixExp[lmat 120., v0];
  rho = Partition[v, nn]; tr = Tr[rho];
  mom = Re @ ({Tr[x . x . rho], Tr[(x . p + p . x) . rho]/2, Tr[p . p . rho]}/tr);
  {mom,
   N @ Flatten[{#[[1, 1]], #[[1, 2]], #[[2, 2]]}] &[
     steady[gam/(6 M kB T)] /. {hbar -> 1, M -> 1, w0 -> 1, gam -> 1/10, kB -> 1, T -> 1}],
   {Abs[tr - 1], Max @ Abs[rho - ConjugateTranspose[rho]],
    Min @ Re @ Eigenvalues[(rho + ConjugateTranspose[rho])/(2 tr)],
    Total @ Abs @ Diagonal[rho][[-4 ;;]]}}]
```

The first two rows agree: the density operator relaxed from the ground state under the full generator reaches the covariance the Lyapunov algebra predicts, an end-to-end check that touches the operators, not the moment equations. The third row holds the run's own legality certificates, trace drift, Hermiticity defect, smallest eigenvalue, and the population of the top truncation levels; all sit at numerical noise, so the evolution neither leaked probability nor left the physical state space, and the agreement above is not laundered by the normalization.

**Why this matters.** The paper motivates $D_{XX}$ by positivity. The symbolic reading is cleaner: $D_{XX}$ is fixed by thermalization, it is the coefficient that recovers the leading quantum correction to the equilibrium state, and positivity follows because that thermal value happens to exceed the positivity threshold. The new term repairs the equilibrium, not merely the spectrum of the generator, and it upgrades the master equation from correct-to-$O(\hbar^0)$ to correct-to-$O(\hbar^2)$ at fixed temperature.

**The tension, made exact.** The match to Gibbs is only to $O(\hbar^2)$, and it must be. Ask instead for the generator whose steady state is *exactly* the Gibbs covariance, at all orders in $\hbar$. The Lyapunov equation inverts, $D=-(A\,\sigma_{\rm Gibbs}+\sigma_{\rm Gibbs}A^\top)$, and returns a unique diffusion matrix:

```wl
lyapInv[a_, s_] := -(a . s + s . Transpose[a]);
sigGibbs = {{sXXgibbs, 0}, {0, sPPgibbs}};
Dreq = FullSimplify @ lyapInv[Amat, sigGibbs];
Dreq // MatrixForm
```

The exact-thermalizing diffusion is diagonal with a vanishing position slot. Read its Dekker coefficients back through the same `diffusion` map:

```wl
FullSimplify @ Solve[diffusion[dxx, dpp, dxp] == Dreq, {dxx, dpp, dxp}]
```

The generator that thermalizes exactly has $D_{XX}=0$: it is the position-coupling-only generator carrying the full Bose-Einstein momentum diffusion $(\gamma M\omega_0/\hbar)\coth(\hbar\omega_0/2k_BT)$, the quantum-corrected Caldeira-Leggett equation. Its Kossakowski determinant is negative at every temperature:

```wl
FullSimplify @ Det[kossakowski[0, gam M w0 Coth[hbar w0/(2 kB T)]/hbar, 0]]
```

So the exact-thermalizing generator is not completely positive: $\det a=-\gamma^2/\hbar^2<0$, the diffusion-constraint floor of Dekker and Valsakumar (Phys. Lett. A **104**, 67 (1984)). Complete positivity needs $D_{XX}\ge\gamma/8Mk_BT>0$, and any $D_{XX}>0$ pushes the steady correlation $\sigma_{XP}=-D_{XX}\hbar^2M$ off its equilibrium value of zero, opening the phase-space current. Exact thermalization and complete positivity are therefore mutually exclusive *on the Caldeira-Leggett friction slice and for a confining potential* (for a free particle the momentum marginal is exactly Maxwellian, though the steady correlation $-\hbar^2MD_{XX}$ survives even there); the wider two-parameter friction family escapes (the quantum-optical damped oscillator is completely positive with an exact Gibbs steady state, see the mean-force section). Both halves of that dichotomy are the theorem of Artini and co-workers on translation-covariant completely positive quadratic master equations, with the entropy production as its thermodynamic reading. The paper's $D_{XX}=\gamma/6Mk_BT$ and Breuer-Petruccione's $\gamma/8Mk_BT$ are then two named compromises inside the completely positive cone: the former is the point that gets the equilibrium *energy* right (the momentum variance exact to $O(\hbar^2)$), the latter the *minimal* completely positive point (the $\det a=0$ boundary), their ratio $4/3$ being the paper's positivity margin.

## B. Decoherence and classicalization, quantified

The two double commutators are diagonal in position and in momentum respectively, so their decohering action is read off directly. The load-bearing identity, $\langle x|[X,[X,\rho]]|x'\rangle = (x-x')^2\langle x|\rho|x'\rangle$, is one cell on a position grid:

```wl
With[{xs = Array[xg, 5]},
  With[{x = DiagonalMatrix[xs], rho = Array[r, {5, 5}]},
   Simplify[x . x . rho - 2 x . rho . x + rho . x . x == Outer[Subtract, xs, xs]^2 rho]]]
```

The double commutator multiplies each matrix element by the squared separation of its grid points, and by Fourier symmetry $[P,[P,\rho]]$ does the same in momentum with weight $(p-p')^2$. Hence an off-diagonal element decays exponentially at a rate quadratic in the separation:

- position superposition of separation $\Delta x$: $\Gamma_x = D_{PP}\,\Delta x^2 = 2\gamma Mk_BT\,\Delta x^2/\hbar^2$;
- momentum superposition of separation $\Delta p$: $\Gamma_p = D_{XX}\,\Delta p^2 = \gamma\,\Delta p^2/(6Mk_BT)$.

The second rate is the new channel and is exactly the paper's "penalizing high-frequency fluctuations": a superposition that is broad in momentum (fine, high-spatial-frequency interference structure) is damped in proportion to $\Delta p^2$. The two rates carry opposite temperature dependence, seen in their $T$-derivatives:

```wl
{D[2 gam M kB T dx^2/hbar^2, T], D[gam dp^2/(6 M kB T), T]}
```

Position decoherence strengthens with temperature (the familiar thermal-wavelength suppression of spatial coherence), while momentum decoherence weakens with temperature. Their ratio,

$$
\frac{\Gamma_p}{\Gamma_x} = \frac{\hbar^2\,\Delta p^2}{12\,M^2 k_B^2 T^2\,\Delta x^2},
$$

says the new channel dominates for momentum-broad superpositions and at lower temperature. This is a concrete, testable fingerprint of the corrected generator that the standard Caldeira-Leggett equation ($D_{XX}=0$) does not produce at all: it decoheres position but never momentum.

## Mean force, not bare Gibbs: where the equilibrium statement ends

The frontier question above has a short answer and a longer lesson. The comparison target is exactly known: the late-time state of the Caldeira-Leggett model is the *mean-force Gibbs state*, the reduced state of the global Gibbs state of system plus bath (Subaşı, Fleming, Taylor and Hu, Phys. Rev. E **86**, 061132 (2012)); for the Ohmic damped oscillator its covariance is a Matsubara sum with a closed form (Grabert, Schramm and Ingold, Phys. Rep. **168**, 115 (1988); the modern mean-force program is Cresser and Anders, Phys. Rev. Lett. **127**, 250601 (2021) and the review of Trushechkin, Merkli, Cresser and Anders, AVS Quantum Sci. **4**, 012301 (2022)). At $\gamma\to 0$ the mean-force state is the bare Gibbs state, so feature A's target was correct at its order; the two part ways at $O(\gamma)$. The short answer is then already on the page: by feature A, the generator's diagonal deviations from bare Gibbs are $O(\gamma^2\hbar^2)$ with no $O(\gamma)$ term at all, while the mean-force corrections computed below are $O(\gamma)$; and its steady correlation $-\gamma\hbar^2/6k_BT$ sits against an exact zero. The deviations are not mean-force corrections in disguise; nothing was hiding. The longer lesson is what holding the mean-force covariance stationary would actually take, with the Caldeira-Leggett friction structure and with the bare Hamiltonian, and what each of those two qualifiers is worth. That requires the mean-force corrections in closed form first.

The exact equilibrium position variance is the Matsubara sum over $\nu_n = 2\pi n k_BT/\hbar$ with the damping entering through $\nu_n^2 + 2\gamma\nu_n + \omega_0^2$, in closed form:

```wl
nu[n_] := 2 Pi n kB T/hbar;
sxxMF = (kB T/M) (1/w0^2 + 2 Sum[1/(w0^2 + nu[n]^2 + 2 gam nu[n]), {n, 1, Infinity}]) //
   FullSimplify
```

A digamma pair at the two relaxation frequencies $\gamma\pm\sqrt{\gamma^2-\omega_0^2}$ of the damped oscillator (the Euler constants of the two harmonic numbers cancel in the difference), one expression for both branches: overdamped the arguments are real, underdamped they are complex conjugates and the imaginary parts cancel. Three independent routes must give one number: this closed form, the raw Matsubara sum, and the fluctuation-dissipation integral $\sigma_{XX}=(\hbar/\pi)\int_0^\infty \coth(\hbar\omega/2k_BT)\,\mathrm{Im}\,\chi(\omega)\,d\omega$ with $\chi=[M(\omega_0^2-\omega^2-2i\gamma\omega)]^{-1}$. At one benign underdamped point:

```wl
Table[Block[{hbar = 1, M = 1, kB = 1, T = 7/10, w0 = 1, gam = g},
   N[{Chop @ N[sxxMF, 20],
     kB T/(M w0^2) + 2 kB T/M NSum[1/(w0^2 + nu[n]^2 + 2 gam nu[n]), {n, 1, Infinity},
        NSumTerms -> 200, WorkingPrecision -> 30],
     NIntegrate[(hbar/Pi) Coth[hbar w/(2 kB T)] Im[1/(M (w0^2 - w^2 - 2 I gam w))],
      {w, 0, Infinity}, WorkingPrecision -> 20]}, 12]], {g, {3/10, 5/2}}]
```

One number three ways in each row, the first row underdamped and the second overdamped, so both branches of the closed form are exercised; a convention slip in any route, a factor of two in the friction, a misplaced $2\pi$ in the Matsubara frequency, would break the agreement loudly. Two sanity limits pin the closed form to known physics. The uncoupled limit must return the Gibbs variance:

```wl
FullSimplify[
  Limit[sxxMF, gam -> 0, Assumptions -> {hbar > 0, M > 0, kB > 0, T > 0, w0 > 0}] - sXXgibbs,
  {hbar > 0, M > 0, kB > 0, T > 0, w0 > 0}]
```

The difference vanishes: at zeroth order in the coupling, mean force and bare Gibbs coincide. And the classical limit must lose the coupling entirely, at every $\gamma$:

```wl
FullSimplify[
  Limit[sxxMF, hbar -> 0, Assumptions -> {M > 0, kB > 0, T > 0, w0 > 0, gam > 0}],
  {M > 0, kB > 0, T > 0, w0 > 0, gam > 0}]
```

Classical equipartition, $\gamma$-free: for the harmonic Caldeira-Leggett model with its counterterm, the classical mean-force correction is identically zero, so every mean-force correction below is a quantum effect. The one place the closed form degenerates is critical damping, where the two relaxation frequencies collide and the square root vanishes; the limit closes:

```wl
Limit[sxxMF, gam -> w0, Assumptions -> {hbar > 0, M > 0, kB > 0, T > 0, w0 > 0}]
```

Finite, a single trigamma at $1+\hbar\omega_0/2\pi k_BT$: the degeneracy is removable. The free-particle boundary is not, and on the generator's side all three entries of the story sit in one limit of the covariance already computed:

```wl
Limit[{steady[dxx][[2, 2]], steady[dxx][[1, 2]], steady[dxx][[1, 1]]}, w0 -> 0,
  Assumptions -> {hbar > 0, M > 0, gam > 0, kB > 0, T > 0, dxx > 0}]
```

The momentum marginal is exactly Maxwellian ($Mk_BT$, no $\hbar$ at any order), the steady correlation $-\hbar^2MD_{XX}$ survives untouched, and the position variance has no stationary value at all; this is the computation behind the free-particle statements in section A and below. Now the object of interest, the $O(\gamma)$ coefficient of $\sigma_{XX}^{\rm MF}-\sigma_{XX}^{\rm Gibbs}$:

```wl
c1 = FullSimplify[(2 kB T/M) Sum[Evaluate[D[1/(w0^2 + nu[n]^2 + 2 gam nu[n]), gam] /. gam -> 0],
    {n, 1, Infinity}]]
```

A trigamma pair at $1\pm i\hbar\omega_0/2\pi k_BT$. Its high-temperature expansion exposes the structure that matters:

```wl
Normal @ Series[c1, {hbar, 0, 7}]
```

The series through order $\hbar^7$ has only odd terms, $\hbar^3$, $\hbar^5$, $\hbar^7$, with no even power anywhere. The leading coefficient is a zeta value:

```wl
FullSimplify[PolyGamma[2, 1] == -2 Zeta[3]]
```

So the true equilibrium position variance carries the correction $-\zeta(3)\,\gamma\hbar^3/2\pi^3M(k_BT)^2$ at leading order: finite, cutoff-free, negative, and *odd in $\hbar$* at every order shown. Neither the number nor the physics is a discovery. The coefficient is the next term of the standard large-$\nu_n$ expansion of the Grabert-Schramm-Ingold Matsubara sum, the first $\gamma$-dependent one (the preceding $\zeta(2)$ term is the $\gamma$-free Wigner-Kirkwood correction, where Hilt, Thomas and Lutz stop, Phys. Rev. E **84**, 031110 (2011), whose reading of the sign, the bath as a continuous position measurement squeezing the variance, carries over). What matters for this section is not the constant but the parity: the correction is odd in $\hbar$. Against it the generator offers nothing to compare: at the paper's point the steady $\sigma_{XX}$ is even in $\gamma$, so its $O(\gamma)$ mean-force correction is not badly approximated, it is absent.

Odd is one obstruction, with a sharp scope. On the Markovian side, the Lyapunov diffusion source contributed by the $n$-th term of the local kernel expansion (the $\delta^{(2n)}$ tower of the essay) is $2\hbar^2$ times the coefficient over $2\hbar$, with $\hbar$-free channel factors; that weight-to-coefficient map is quoted from the essay's Part IV, and given it, the $\hbar$-content is even at every order:

```wl
c2n[n_] := 4 eta (-1)^n BernoulliB[2 n]/(2 n)! (hbar/(kB T))^(2 n - 1);
Simplify[Exponent[2 hbar^2 c2n[n]/(2 hbar), hbar], n \[Element] NonNegativeIntegers]
```

The exponent is $2n$ at every order, even by construction, and, just as important, the sources are diagonal only: on stationary paths each $\delta^{(2n)}$ folds into the $(X,X)$ or $(P,P)$ channel (the essay's Part VII), so the local truncation generates no cross channel at any order. The conclusion is then itself a one-line theorem on the general covariance: feed the Lyapunov solution any even, diagonal source and the steady covariance is even in $\hbar$,

```wl
With[{par = sigmaGen /. {dxp -> 0, dxx -> u hbar^(2 n - 2), dpp -> v hbar^(2 n - 2)}},
  FullSimplify[par - (par /. hbar -> -hbar), n \[Element] Integers]]
```

identically, for arbitrary weights $u,v$ and every order $n$. An even function cannot contain the odd $\zeta(3)$ term. What the local expansion cannot produce, the exact theory does: the missing object is a cross coefficient $D_{XP}$, which enters $\sigma_{XX}$ but not $\sigma_{XP}$, and the inversion below shows it is precisely what holding the mean-force covariance demands.

The momentum entry fails harder. The momentum variance has its own Matsubara representation (Grabert-Schramm-Ingold), $\sigma_{PP}^{\rm MF}=(Mk_BT)[1+2\sum_n(\omega_0^2+\nu_n\hat\gamma(\nu_n))/(\nu_n^2+\nu_n\hat\gamma(\nu_n)+\omega_0^2)]$ with the Drude damping $\hat\gamma(\nu)=2\gamma\Omega/(\nu+\Omega)$; its $O(\gamma)$ coefficient comes by the same derivative-under-the-sum route as $c_1$, and its large-cutoff expansion:

```wl
p1 = (M kB T) Sum[Evaluate[
     2 D[(w0^2 + nu[n] (2 gam Om/(nu[n] + Om)))/(w0^2 + nu[n]^2 + nu[n] (2 gam Om/(nu[n] + Om))),
        gam] /. gam -> 0], {n, 1, Infinity}];
p1asy = FullSimplify[Normal @ Series[p1, {Om, Infinity, 0}]]
```

The expansion carries $\log(\hbar\Omega/2\pi k_BT)$: the momentum correction has no cutoff-free limit. Its logarithmic slope is exact:

```wl
FullSimplify[D[p1asy, Om] Om]
```

The representation itself deserves the same independent check the position variance got: the Matsubara sum at finite $\gamma$ against the fluctuation-dissipation integral with the Drude response $\chi(\omega)=[M(\omega_0^2-\omega^2-i\omega\,2\gamma\Omega/(\Omega-i\omega))]^{-1}$, two genuinely different routes to $\sigma_{PP}^{\rm MF}$:

```wl
Block[{hbar = 1, M = 1, kB = 1, T = 7/10, w0 = 1, Om = 50, gam = 1/10},
  N[{M kB T (1 + 2 NSum[
       (w0^2 + nu[n] 2 gam Om/(nu[n] + Om))/(w0^2 + nu[n]^2 + nu[n] 2 gam Om/(nu[n] + Om)),
       {n, 1, Infinity}, NSumTerms -> 3000, WorkingPrecision -> 25]),
    M^2 hbar/Pi NIntegrate[w^2 Coth[hbar w/(2 kB T)] Im[1/(M (w0^2 - w^2 - I w 2 gam Om/(Om - I w)))],
      {w, 0, Infinity}, WorkingPrecision -> 20]}, 10]]
```

The two routes coincide, so the momentum representation, and with it the log and its slope, stands on its own numbers rather than on the citation. The true equilibrium momentum variance grows as $(2M\hbar\gamma/\pi)\ln(\hbar\Omega/2\pi k_BT)$, divergent in exactly the scaling limit the paper takes, $\Omega\to\infty$ at fixed $T$. The divergence is a classical result, the ultraviolet sensitivity of the Ohmic momentum dispersion of Grabert, Weiss and Talkner (Z. Phys. B **55**, 87 (1984)); what is presented here is its form in the paper's own ordering of scales. The ordering matters and must be stated: the logarithm, and its oddness in $\hbar$, belong to $\hbar\Omega\gg k_BT$, the regime the paper works in. In the opposite ordering $\hbar\Omega\ll k_BT$ the same sum gives $\gamma M\hbar^2\Omega/6k_BT$, even in $\hbar$ and linear in the cutoff, which is the form in Hilt, Thomas and Lutz's high-temperature appendix. Either way the conclusion for the family is the same exclusion: the $O(\gamma)$ momentum correction is cutoff-dependent in every ordering, so a generator whose coefficients are free of $\Omega$ cannot track it, and matching the *bare Gibbs* variance at $O(\gamma^0)$, which is what the paper's $D_{XX}$ achieves, is the only cutoff-free equilibrium statement available. There is an irony worth stating: the generator built to penalize high-frequency quantum fluctuations necessarily omits the physical high-frequency dressing of the momentum, because that dressing is exactly what diverges when the cutoff is removed.

The correlation entry needs no computation at all. The global Gibbs state is invariant under time reversal, which flips $P$ and every bath momentum while fixing $X$, and $\{X,P\}$ is odd under it, so $\sigma_{XP}^{\rm MF}=0$ exactly, at every order in $\gamma$ and $\hbar$. The completely positive steady state cannot comply once $D_{PP}$ is pinned near the Einstein value $2\gamma Mk_BT/\hbar^2$ by classical equipartition, which any candidate thermal generator must respect: then complete positivity requires $D_{XX}\ge\gamma/8Mk_BT$ (feature A), and the steady correlation is $\sigma_{XP}=-\hbar^2MD_{XX}$, so every such completely positive member carries $|\sigma_{XP}|\ge\gamma\hbar^2/8k_BT$, bounded away from the true value by the positivity floor itself (the constraint of Artini and co-workers, read as an equilibrium error bound; the paper's point sits at $\gamma\hbar^2/6k_BT$, above the floor by the same $4/3$ as in feature A).

**The inversion, and what it costs on the Caldeira-Leggett slice.** Ask directly for the constant-coefficient generator whose steady covariance is a prescribed equilibrium-grade $\mathrm{diag}(a,b)$, diagonal because time reversal forbids a steady correlation. Two restrictions define the class being asked about, and both must be said out loud. The friction is the *translation-covariant*, position-coupled one of the paper's generator, acting on the momentum row alone; this is fidelity to the generator as stated, not a microscopic necessity (a position-coupled bath at finite cutoff produces frequency-dependent friction, which is where the section ends up), and the general quadratic Lindbladian has a wider, two-parameter friction (Isar, Sandulescu and Scheid, Int. J. Mod. Phys. E **3**, 635 (1994)), revisited below; "slice" abbreviates the translation-covariant class throughout. And the Hamiltonian is the bare one, no $O(\gamma)$ frequency renormalization, matching the generator as the paper states it. Within the slice the Lyapunov map inverts uniquely:

```wl
sigT = {{aa, 0}, {0, bb}};
solT = First @ Solve[diffusion[dxx, dpp, dxp] == lyapInv[Amat, sigT], {dxx, dpp, dxp}];
solT
```

One generator per target, holding its second moments stationary (stationary moments, not equilibrium in the strong sense: quantum detailed balance is a stricter requirement, the one imposed in the general characterization of thermal Gaussian semigroups by Toscano and Nicacio, Phys. Rev. A **106**, 062207 (2022), and not tested here). The anatomy is forced: $D_{XX}=0$ identically, which is what keeps $\sigma_{XP}=0$; the momentum channel carries the full target variance; the cross coefficient carries the anisotropy, $D_{XP}=(b-aM^2\omega_0^2)/2\hbar^2M$. This is the Caldeira-Leggett-slice specialization of thermal-generator formulas that exist in the literature (Isar-Sandulescu-Scheid give exactly $D_{qq}\propto\lambda-\mu$, vanishing on this slice, and $D_{pq}=0$ at the Gibbs target), so the inversion itself is not new; what it is being used for is. For the bare Gibbs target the cross coefficient vanishes, by the isotropy of the Gibbs covariance:

```wl
FullSimplify[sPPgibbs - M^2 w0^2 sXXgibbs]
```

So at $O(\gamma^0)$ the unique stationary-moment generator is the quantum-corrected Caldeira-Leggett equation of section A, cross-free. At $O(\gamma)$ the mean-force target breaks isotropy, and because the inversion's cross coefficient is linear in the target, substituting $a=\sigma_{XX}^{\rm Gibbs}+\gamma c_1$, $b=\sigma_{PP}^{\rm Gibbs}+\gamma p_1$ into it reduces, by the isotropy just shown, to $D_{XP}=\gamma\,(p_1-M^2\omega_0^2c_1)/2\hbar^2M+O(\gamma^2)$: the cutoff dependence and the $\zeta(3)$ term land in the one coefficient the local kernel expansion cannot supply. The positivity ledger follows in one line from the general covariance's $\sigma_{XP}=-\hbar^2MD_{XX}$: a diagonal target forces $D_{XX}=0$, hence $\det a=-\gamma^2/\hbar^2-4D_{XP}^2$, and the same line holds verbatim under any drift renormalization $M\to M^*$, $\omega_0\to\omega^*$. The cells spell it out. The slice inversion gives

```wl
FullSimplify[Det[kossakowski[dxx, dpp, dxp]] /. solT]
```

$\det a = -\gamma^2/\hbar^2 - (b-aM^2\omega_0^2)^2/\hbar^4M^2$. The quadratic term is *not* an equilibration cost; it is the price of insisting on the bare Hamiltonian. Every diagonal covariance above the vacuum floor ($2\sqrt{ab}>\hbar$) is the Gibbs covariance of a unique effective oscillator $(M^*,\omega^*)$ at the same temperature (the mean-force effective Hamiltonian of Hilt, Thomas and Lutz; the mean-force literature's own prescription is to renormalize the system Hamiltonian, Timofeev and Trushechkin, arXiv:2204.00599), and re-running the inversion with the renormalized drift removes the quadratic term exactly:

```wl
AmatR = {{0, 1/Ms}, {-Ms ws^2, -2 gam}};
solR = First @ Solve[diffusion[dxx, dpp, dxp] == lyapInv[AmatR, sigT], {dxx, dpp, dxp}];
FullSimplify[Det[kossakowski[dxx, dpp, dxp]] /. solR /. bb -> aa Ms^2 ws^2]
```

What no renormalization removes: $\det a=-\gamma^2/\hbar^2$ exactly, the Caldeira-Leggett value, still strictly negative. On the slice, holding *any* diagonal covariance exactly stationary costs $\det a\le-\gamma^2/\hbar^2$, Hamiltonian renormalization or not. This statement is not new either: the floor is the diffusion constraint of Dekker and Valsakumar (Phys. Lett. A **104**, 67 (1984)), and its failure on the slice at *every* frequency and mass, which is exactly the "any target, any renormalization" content, is stated below Eq. (3.26) of Isar's review of the thermal-generator formulas (quant-ph/0602149): at $\mu=\lambda$ the thermal-stationarity constraint reads $0\ge\lambda^2$. What the cells add is only the packaging as a target-independent $\det a$ ledger. Complete positivity therefore buys its Lindblad form by giving up exactness, and the paper's point is the thermal-optimal completely positive compromise at $O(\gamma^0)$.

The slice matters, because off it the story inverts. Widen the friction to the two-parameter family, drift $A=\{\{-(\lambda-\mu),1/M\},\{-M\omega_0^2,-(\lambda+\mu)\}\}$, and evaluate the Dekker-Valsakumar slack of the generator that holds the exact Gibbs covariance stationary:

```wl
A2 = {{-(lam - mu), 1/M}, {-M w0^2, -(lam + mu)}};
slack = With[{d = lyapInv[A2, sigGibbs]/2}, FullSimplify[Det[d] - lam^2 hbar^2/4]];
slack
```

The slack is $\hbar^2[(\lambda^2-\mu^2)\coth^2(\hbar\omega_0/2k_BT)-\lambda^2]/4$. On the Caldeira-Leggett slice $\mu=\lambda$ it is $-\hbar^2\lambda^2/4$, negative at every temperature, the statement above. At $\mu=0$ it is $+\hbar^2\lambda^2\operatorname{csch}^2(\hbar\omega_0/2k_BT)/4$, *positive at every temperature*: the textbook quantum-optical damped oscillator, completely positive, constant coefficients, exact Gibbs stationary state to all orders in $\hbar$. This dichotomy is the published theorem of Artini and co-workers (arXiv:2507.23322): every translation-covariant completely positive quadratic master equation with a confining potential ends in a non-equilibrium steady state, and only breaking translation covariance, by fine-tuning against the system Hamiltonian, restores an equilibrium; the cell above is their theorem as a single closed form. Two refinements are worth reading off it. The tension needs the potential as well as the friction: for a free particle the momentum marginal is exactly Maxwellian on the slice (the equilibration statement Artini and co-workers make there; the steady correlation $-\hbar^2MD_{XX}$ still survives, the limit cell above), which is consistent with the essay's Part VII, where the free-particle kernel tower truncates exactly at two operators. And the boundary of the completely positive equilibrating region is exactly critical, solvable outright:

```wl
FullSimplify @ Reduce[slack >= 0, mu, Reals]
```

The window is $|\mu|\le\lambda\,\mathrm{sech}(\hbar\omega_0/2k_BT)$, and its distance from the Caldeira-Leggett point at high temperature:

```wl
Normal @ Series[lam (1 - Sech[hbar w0/(2 kB T)]), {T, Infinity, 2}]
```

The excluded band around $\mu=\lambda$ has width $\lambda\hbar^2\omega_0^2/8(k_BT)^2$ at high temperature: the obstruction sits on a knife edge, not in a broad forbidden zone. The physical price of $\mu\ne\lambda$ stays visible: symmetric damping of both quadratures is an energy-basis (rotating-wave) coupling, not the position-coupled bath of the microscopic model, and it gives up the translation-covariant free-particle limit that the Caldeira-Leggett friction keeps.

The exact dynamics shares the slice's signature. The Hu-Paz-Zhang master equation, derived from the same Feynman-Vernon influence functional without the locality truncation, carries no $(\hat P,\hat P)$ term and does carry a time-dependent cross-diffusion coefficient (the $\Gamma(t)f(t)$ term in Halliwell and Yu, quant-ph/9508004), and its late-time state is the mean-force state (Subaşı and co-workers). The ledger then makes one falsifiable prediction about that published object: since the asymptotic Hu-Paz-Zhang generator sits on the slice, ends on a diagonal covariance, and so must have a vanishing $(\hat P,\hat P)$ coefficient, its asymptotic Kossakowski data must obey $\det a\le-\gamma_\infty^2/\hbar^2$, exact quantum Brownian motion is asymptotically non-completely-positive with a quantified deficit. The quantitative identification of the asymptotic $\Gamma(t)f(t)$ with the slice inversion's $D_{XP}$ is untested here, and is not even well posed until the frequency convention is pinned, since Hu-Paz-Zhang also renormalizes the frequency and $D_{XP}$ shifts linearly under $\omega_0^2\to\omega^{*2}$ at fixed target. Next steps 3 states both checks with the convention fixed.

**Verdict, and what is new here.** The completely positive steady state's deviations from bare Gibbs are not mean-force corrections in disguise (the short answer needed only feature A and time reversal). Nearly everything else in this section is the literature, assembled: the mean-force state as the true late-time state (Subaşı et al.), its closed-form covariance and high-temperature expansion (Grabert-Schramm-Ingold; the $\zeta(3)$ term is simply that expansion's first $\gamma$-dependent order), the momentum-dispersion cutoff sensitivity in both orderings (Grabert-Weiss-Talkner; Hilt-Thomas-Lutz), the effective mean-force oscillator (Hilt-Thomas-Lutz), the thermal-generator inversion formulas and the slice's thermal-stationarity impossibility at every frequency and mass (Isar-Sandulescu-Scheid; Isar quant-ph/0602149, Eq. (3.26)), the positivity floor (Dekker-Valsakumar), the translation-covariance dichotomy including the fine-tuned escape (Artini and co-workers), the general second-order mismatch for completely positive steady states (Lee and Yeo, discrete spectra) and Redfield's escape through frequency-dependent rates (Thingna-Wang-Hänggi). What this section adds is small and specific: the parity scoping (no one had asked whether the $\delta^{(2n)}$ tower's even, diagonal Lyapunov sources can produce the odd-in-$\hbar$ mean-force corrections; they cannot), the mean-force target run through the slice inversion with the resulting $D_{XP}$ in closed form, the exactly-critical boundary $|\mu|\le\lambda\,\mathrm{sech}(\hbar\omega_0/2k_BT)$ read off the slack, and the prediction that the asymptotic Hu-Paz-Zhang generator is non-completely-positive with $\det a\le-\gamma_\infty^2/\hbar^2$. The practical audience is first of all this document's own feature A, which motivates $D_{XX}$ by thermalization and thereby invites over-reading: the equilibrium content of the paper's generator is the bare Gibbs state at $O(\gamma^0)$ for quadratic $V$, and no $O(\gamma)$ equilibrium physics should be read off it (for anharmonic $V$ even the $O(\gamma^0)$ verdict is open, feature Engine's territory, since the position diffusion then sustains entropy production, Artini). In the generator's own validity window the missed position correction is doubly small, suppressed by $\hbar\gamma/k_BT$ and $(\hbar\omega_0/k_BT)^2$ together; the momentum entry, cutoff-divergent, is the one with teeth. A reader who needs the $O(\gamma)$ equilibrium has published routes: the mean-force Hamiltonian renormalization (Timofeev-Trushechkin), Redfield-class frequency-dependent rates (Thingna-Wang-Hänggi), or the canonically consistent master equation of Becker, Schnell and Thingna (Phys. Rev. Lett. **129**, 200403 (2022)), whose steady state is the mean-force state to second order and which is, consistently with everything above, not completely positive. Two reconciliations keep the documents consistent: the essay's convergent tower refines the $O(\gamma^0)$ covariance toward bare Gibbs, exactly the right target at that order, and says nothing about $O(\gamma)$; and the paper makes no equilibrium claim at all, so nothing here contradicts it.

## C. The Wigner Fokker-Planck picture (set up)

TO DEVELOP. Wigner-transforming the operator equation gives, for quadratic $V$, an exact Klein-Kramers / Fokker-Planck equation

$$
\partial_t W = -\frac{p}{M}\partial_x W + V'(x)\,\partial_p W + 2\gamma\,\partial_p(pW) + \underbrace{\hbar^2 D_{PP}}_{=\,2\gamma Mk_BT}\,\partial_p^2 W + \underbrace{\hbar^2 D_{XX}}_{=\,\hbar^2\gamma/(6Mk_BT)}\,\partial_x^2 W - \underbrace{2\hbar^2 D_{XP}}_{=\,0\ \text{here}}\,\partial_x\partial_p W ,
$$

written for the general three-coefficient family and specialized by the underbraces to the paper's point, where the cross term vanishes (it is the very term the mean-force section shows correct equilibration would need). The momentum diffusion is the classical Einstein term; the position diffusion, from the new $(P,P)$ term, is of order $\hbar^2$ and is a purely quantum smoothing in position that erases sub-Planck (interference) structure. This is the classicalization of feature B seen in phase space, and it makes the classical limit transparent: as $\hbar\to 0$ the position diffusion vanishes and the equation becomes the classical Kramers equation. Tasks: derive the transform and pin the coefficients symbolically (confirming the $O(\hbar^2)$ position diffusion), solve the Gaussian propagator, and hand the anharmonic PDE (with its Moyal $V'''\partial_p^3$ corrections) to the simulation track.

## Engine. The moment hierarchy for general $V$ (the tool for A, B, and beyond)

Features A and B are the harmonic-case shadow of a general structure. The adjoint generator gives, for any observable $O$, $d\langle O\rangle/dt=\langle\mathcal{L}^\dagger[O]\rangle$ with

$$
\mathcal{L}^\dagger[O] = \tfrac{i}{\hbar}[H_S,O] + i\tfrac{\gamma}{\hbar}\{P,[X,O]\} - D_{PP}[X,[X,O]] - D_{XX}[P,[P,O]] - 2D_{XP}[P,[X,O]],
$$

evaluated with $[X,P]=i\hbar$. This is a Weyl-algebra computation, so the tool is `/nca-route` (encode $X,P$ as generators, apply $\mathcal{L}^\dagger$, reduce to a normal form, cross-check against a Fock-matrix carrier). The first two moment layers, for any $V$:

$$
\frac{d\langle X\rangle}{dt} = \frac{\langle P\rangle}{M}, \qquad
\frac{d\langle P\rangle}{dt} = -\langle V'(X)\rangle - 2\gamma\langle P\rangle,
$$
$$
\frac{d\langle X^2\rangle}{dt} = \frac{\langle\{X,P\}\rangle}{M} + 2\hbar^2 D_{XX}, \qquad
\frac{d\langle P^2\rangle}{dt} = -\langle\{P,V'(X)\}\rangle - 4\gamma\langle P^2\rangle + 2\hbar^2 D_{PP},
$$
$$
\frac{d\langle\{X,P\}\rangle}{dt} = \frac{2\langle P^2\rangle}{M} - \langle\{X,V'(X)\}\rangle - 2\gamma\langle\{X,P\}\rangle - 4\hbar^2 D_{XP} .
$$

The diffusion enters only as three constant sources: $2\hbar^2 D_{XX}$ (position, the new term), $2\hbar^2 D_{PP}$ (momentum, the heating), and $-4\hbar^2 D_{XP}$ in the correlation equation (zero at the paper's point, but the very coefficient the mean-force section turns on). All five displayed equations are operator identities, checkable at once on a truncated Fock carrier with $\hbar$, $M$, $\omega_0$, the anharmonicity, and every coefficient symbolic (keeping $\hbar$ symbolic is what pins the $\hbar^2$ powers of the sources; a carrier at $\hbar=1$ would pass a wrong power silently), on an interior block clear of the truncation edge:

```wl
Block[{nn = 18, a, x, p, comm, acomm, vp, hs, ldag, obs, rhs},
  a = SparseArray[Band[{1, 2}] -> Sqrt[Range[nn - 1]], {nn, nn}];
  x = Sqrt[hbar/(2 M w0)] (a + Transpose[a]); p = I Sqrt[hbar M w0/2] (Transpose[a] - a);
  comm[u_, v_] := u . v - v . u; acomm[u_, v_] := u . v + v . u;
  vp = M w0^2 x + g3 x . x + g4 x . x . x;
  hs = p . p/(2 M) + M w0^2 x . x/2 + g3 x . x . x/3 + g4 x . x . x . x/4;
  ldag[o_] := (I/hbar) comm[hs, o] + (I gam/hbar) acomm[p, comm[x, o]] -
     dpp comm[x, comm[x, o]] - dxx comm[p, comm[p, o]] - 2 dxp comm[p, comm[x, o]];
  obs = {x, p, x . x, p . p, acomm[x, p]};
  rhs = {p/M, -vp - 2 gam p,
    acomm[x, p]/M + 2 hbar^2 dxx IdentityMatrix[nn],
    -acomm[p, vp] - 4 gam p . p + 2 hbar^2 dpp IdentityMatrix[nn],
    2 p . p/M - acomm[x, vp] - 2 gam acomm[x, p] - 4 hbar^2 dxp IdentityMatrix[nn]};
  Union @ Flatten @ Simplify @ Normal[(ldag /@ obs - rhs)[[All, 1 ;; 8, 1 ;; 8]]]]
```

All five ladder lines vanish identically on the interior block, for a genuinely anharmonic force $V'=M\omega_0^2X+g_3X^2+g_4X^3$ with two independent anharmonicities, pinning every constant source with its $\hbar^2$ power and the sign of the cross term. The closure-breakers are exactly $\langle V'(X)\rangle$, $\langle\{P,V'(X)\}\rangle$, $\langle\{X,V'(X)\}\rangle$: for quadratic $V$ they are linear in the second moments and the system closes into the Lyapunov flow behind feature A; for anharmonic $V$ they climb the hierarchy and name precisely what a cumulant closure must model. TO DEVELOP: run `/nca-route` to produce these symbolically for general $V$, verify against the Fock carrier, and derive the energy balance $d\langle H_S\rangle/dt$ (heating rate versus friction, the fluctuation-dissipation content).

## Further features (identified)

- **Relaxation spectrum.** Now computed in section A: the drift eigenvalues are $-\gamma\pm\sqrt{\gamma^2-\omega_0^2}$ for the first moments, their pairwise sums governing the covariance, so the thermalization time is $\sim 1/2\gamma$, independent of temperature and $\hbar$.
- **Entropy production / second law.** The completely positive steady state is a non-equilibrium steady state, so its entropy production rate is nonzero. The closed form for the Gaussian case, $\Pi_{\rm ss}\propto D_{XX}$, and its reading as the thermodynamic signature of the $(P,P)$ term, are given by Artini and co-workers (arXiv:2507.23322) and by Bernád, Homa and Csirik (Eur. Phys. J. D **72**, 212 (2018)); the tension is made exact for the paper's point in section A. What is open *here* is the transient entropy production from a non-Gaussian initial state under anharmonic $V$, which the closed-form track cannot reach and the simulation track must.

## Next steps

1. Fill feature C (the Wigner Fokker-Planck) symbolically: the phase-space face of feature B and the cleanest classical limit.
2. Run `/nca-route` on the engine for general $V$, promoting features A and B off the harmonic case and setting up the anharmonic simulation track.
3. Test the prediction the ledger makes about exact quantum Brownian motion: the asymptotic Hu-Paz-Zhang generator should be non-completely-positive with $\det a\le-\gamma_\infty^2/\hbar^2$, and its asymptotic cross coefficient should match the slice inversion's $D_{XP}$ once the renormalized frequency is used in the drift on both sides. A caution on sources: the closed-form Hu-Paz-Zhang coefficients of Homa, Bernád and Csordás (arXiv:2211.15722) are evaluated at $T=0$, where the high-temperature generator does not exist, so they support only a $T\to 0$ consistency check; the finite-$T$ comparison needs the finite-$T$ coefficients (Hu-Paz-Zhang; Halliwell-Yu). Secondary: the Gaussian relative entropy $S(\rho_{\rm ss}\Vert\rho_{\rm MF})$ as the operational size of the equilibrium error, the $O(\gamma^2)$ anatomy, and the regularization dependence of the finite part beside the $\log\Omega$.
4. Carry out the transient half of the Hu-Paz-Zhang benchmark the essay sets up (`caldeira-leggett-high-T-lindblad.md`, Part VIII). The equilibrium endpoints are now settled (the exact dynamics ends at the mean-force covariance, the generator at its even-in-$\hbar$ steady state); what remains is the approach, a closed-form Gaussian contest between truncation error and memory.
