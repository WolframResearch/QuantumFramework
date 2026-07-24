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

**Critique: 2 of 20 parts.** The adversarial pass was interrupted by an account spend limit that
terminated the remaining reviewer agents mid-run. What completed:

- Part 4: `OPEN ISSUES: 0`. The reviewer hand-derived every closed form in the plan (ladder
  normalization, the truncated commutator defect and its tracelessness, the truncated top eigenvalue
  $N/2$, the displacement splitting, the squeezed quadrature variance, the oscillator-length lab
  scales) and confirmed each.
- Part 6: `OPEN ISSUES: 1`, since fixed. The closing limit of 6.5 claimed $S\to\mathbb{1}$ as
  $V_0\to0$; on this document's own convention, with reflection on the diagonal, free space gives
  $S\to\sigma_x$, and $\mathbb{1}$ would describe a hard wall. Two minor items were fixed with it
  (the antisymmetric $S$ eigenvalue is $r-t$, and one entry defined its state by cross-reference).

**The other 18 plans are drafted but ungraded.** Their physics is stated but unverified by a fresh
context, so treat their formulas as claims to check at authoring time rather than as settled. The
Part 6 finding is the calibration: a wrong closing limit survived drafting and was caught only by
the adversarial read.

**Route revisions taken: none.** No completed critique refuted a `Route-Table.md` verdict, so its
revision log is empty and every plan cites the routes as originally probed.

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
