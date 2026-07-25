# Plans: Parts 4 through 23

Outline plans for the worked solutions to the remaining 130 questions of
`../Question-List.md`. Parts 1 through 3 already have full solutions and are not planned here.

Each `Part-NN-Plan.md` gives, for every question in that Part: the question verbatim as a heading,
then one paragraph naming the explicit example with its parameters, the Wolfram Language route with
its named functions, the traps that bite it, the check that could refute the result, and the closing
move. `All-Parts-Plan.md` is the same content concatenated for straight-through reading.

## Index

| plan | theme | questions |
|---|---|---|
| [Part-04-Plan.md](Part-04-Plan.md) | The harmonic oscillator | 8 |
| [Part-05-Plan.md](Part-05-Plan.md) | The time-dependent Schrodinger equation as a PDE | 10 |
| [Part-06-Plan.md](Part-06-Plan.md) | Scattering in one dimension | 6 |
| [Part-07-Plan.md](Part-07-Plan.md) | Approximation methods | 10 |
| [Part-08-Plan.md](Part-08-Plan.md) | General theorems and structural methods | 6 |
| [Part-09-Plan.md](Part-09-Plan.md) | Periodic potentials and band structure | 6 |
| [Part-10-Plan.md](Part-10-Plan.md) | Two and three dimensions: separation of variables | 7 |
| [Part-11-Plan.md](Part-11-Plan.md) | Orbital angular momentum in continuous space | 6 |
| [Part-12-Plan.md](Part-12-Plan.md) | Central potentials and the hydrogen atom | 9 |
| [Part-13-Plan.md](Part-13-Plan.md) | Electromagnetic coupling | 6 |
| [Part-14-Plan.md](Part-14-Plan.md) | Spin coupled to spatial motion | 6 |
| [Part-15-Plan.md](Part-15-Plan.md) | Identical particles in continuous space | 5 |
| [Part-16-Plan.md](Part-16-Plan.md) | Density operators, mixed states, and the Wigner function | 6 |
| [Part-17-Plan.md](Part-17-Plan.md) | Continuous-variable quantum optics and information | 7 |
| [Part-18-Plan.md](Part-18-Plan.md) | Open quantum systems in continuous space | 5 |
| [Part-19-Plan.md](Part-19-Plan.md) | Path integrals | 6 |
| [Part-20-Plan.md](Part-20-Plan.md) | Three-dimensional scattering theory | 6 |
| [Part-21-Plan.md](Part-21-Plan.md) | Nonlinear and mean-field wave mechanics | 4 |
| [Part-22-Plan.md](Part-22-Plan.md) | Relativistic wave equations | 5 |
| [Part-23-Plan.md](Part-23-Plan.md) | From one particle to fields: the second-quantization bridge | 6 |
| | **total** | **130** |

Shared files: [Shared-Brief.md](Shared-Brief.md) is the contract every plan follows,
[Example-Ledger.md](Example-Ledger.md) pins one example per question (collision-checked against the
`PIPELINE.md` Appendix A ledger for Parts 1 through 3), and [Route-Table.md](Route-Table.md) holds
the nine kernel-probed solver verdicts the plans cite.

## Status, honestly

**Drafting: complete.** All 20 plans exist, 130 entries, one per question. Mechanically verified:
every entry heading matches its question in `../Question-List.md` verbatim including the
`[BSc]`/`[MSc]` tag, no em or en dashes anywhere, no leftover scaffolding sections.

**Critique: 7 of 20 parts** (4, 5, 6, 12, 14, 17, 22), all graded by fresh-context reviewers and all
findings repaired. Every closed form the reviewers hand-derived was correct; the failures were one
level above the algebra, in signs, conventions, and the transfer of probe evidence into prose, and
they were systematically invisible to the checks the entries proposed for themselves.

The four that would have inverted the physics:

- Part 17 wrote the two-mode squeeze generator with the sign that anti-squeezes the pair the entry
  claims to squeeze, turning the Duan criterion from a violation into no violation at any squeezing.
  Both of that entry's own checks are provably sign-blind.
- Part 14 paired a positive-charge kinetic term with a negative-charge Zeeman term, mirroring the
  Landau degeneracy labels; and quoted a purity floor of $\tfrac12$ for a state whose branch norms
  pin the populations at $\tfrac34$ and $\tfrac14$ for every separation, so equal weights are
  unreachable and the floor is $\tfrac58$.
- Part 6 claimed $S\to\mathbb{1}$ for free space, where its own convention gives $\sigma_x$;
  $\mathbb{1}$ would describe a hard wall.

**The other 13 plans are drafted but ungraded** (7, 8, 9, 10, 11, 13, 15, 16, 18, 19, 20, 21, 23).
Their physics is stated but unverified by a fresh context, so treat their formulas as claims to check
at authoring time rather than as settled. The graded parts averaged about five findings each, with a
physics-inverting error in three of the five in the first batch, which is the calibration for how
much to trust an ungraded one.

**Route revisions taken: one.** `Route-Table.md` revision R1 corrects the C8 (Dirac hydrogen)
verdict: the recorded degeneracy $\kappa\to-\kappa-1$ is the $\kappa\!\leftrightarrow\!l$ dictionary
misread as a degeneracy, mapping $2p_{1/2}$ onto the state the fine structure splits from it. The
correct map is $\kappa\to-\kappa$, since the energy depends on $\kappa$ only through $|\kappa|$. The
shooting boundary condition was likewise the $\kappa=-1$ instance quoted as general. Both were
claims the verdict wrote around its probes rather than probe outputs, which is the useful signal:
what the kernel measured held up, what got written around the measurements did not.

**Concerns raised by drafters** (recorded in the entries, each needing an authoring decision):

- 6.6: the pinned Levinson value $\delta(0)-\delta(\infty)=\pi$ is the winding of the even
  eigenphase $2\delta_{+}$; the phase shift itself obeys the half-integer 1D form
  $\delta_{+}(0)-\delta_{+}(\infty)=\pi(n_{+}-\tfrac12)=\pi/2$ for one bound state. The reviewer
  confirmed this analysis, so the convention must be pinned before $\pi$ is quoted.
- 10.6: the attractive Coulomb effective potential has no local maximum, so "barrier height" names
  no finite number; the entry reads the centrifugal barrier as the zero crossing, the well depth,
  and the inner turning point instead.
- 16.5: the number state's Glauber $P$ is a delta derivative, not a Gaussian, so the Gaussian $P$
  anchor is the thermal state rather than a member of the pinned trio.

## What authoring inherits

Every route in `Route-Table.md` was chosen by probing a live kernel (WL 15.0.0), not by recall, and
the traps below are the ones that bite silently, each now paired with a loud detector:

- `DSolve` never quantizes on its own. Imposing decay at both ends hangs it; the Dirac coupled
  system hangs even with no boundary conditions at all (use the squared Coulomb form), and a
  `Piecewise` coefficient makes it echo its input back unevaluated with no message.
- `NDEigensystem` returns eigenvalues smallest in magnitude, so an unshifted call silently returns
  none of the bound spectrum; without an explicit `DirichletCondition` it silently solves the
  Neumann problem instead. An $\varepsilon$ cutoff at the origin is actively harmful for $l\ge1$.
- Refining a time-dependent grid at default tolerances makes the answer worse, non-monotonically,
  and norm conservation cannot detect wall reflection, because a unitary wall reflects at norm 1.
- Agreement between two routes sharing a truncation or a box size certifies neither: measured
  agreement to $10^{-7}$ while both were wrong by $5\times10^{-3}$. Only an exact benchmark or a
  parameter sweep detects it.
- Self-consistent iteration on the Gross-Pitaevskii equation limit-cycles while its eigenvalue
  plateau imitates convergence; imaginary-time split-step is first order in $dt$, not second.
- A converged solver can fail a textbook benchmark because the benchmark is the wrong one: the
  Landau-Zener asymptote is invalid at finite sweep time once the transition probability drops
  below $(\Delta/2vT)^2$.
- Radial WKB with the Langer correction reproduces $E_n=-1/(2n^2)$ exactly in closed form; without
  it the answer is off by 6 percent. Turning-point integrals stall unless the turning point is
  scaled out, and evaluating an antiderivative across one can silently return zero.
