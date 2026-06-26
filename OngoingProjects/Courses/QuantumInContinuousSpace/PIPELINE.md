# Pipeline & Brief: "Quantum in Continuous Space" question list

The continuous-space (wave-mechanics) companion to
`../QuantumInDiscreteSpace/Question-List.md`. Deliverable so far: `Question-List.md` (the
curriculum). This file is the reusable brief and the build tracker.

## DRAFT BRIEF (adversarial-draft, Compose track, research-rigor mode)

- **Goal:** A dependency-ordered question list that teaches continuous-space quantum mechanics
  by *computing* each item in pure Wolfram Language; a reader works straight down and builds or
  solves each task. The Schrodinger equation read two ways: as an ODE eigenvalue problem (the
  time-independent spectrum) and as a PDE (time-dependent propagation).
- **Reader:** Advanced BSc / early MSc physicist; strong calculus, linear algebra, ODE/PDE,
  Fourier analysis; wants "how do I compute X", not a proof course.
- **Claim (spine):** Every concept of single/few-particle quantum mechanics in $L^2$ has a
  short executable WL realization, either a closed form via special functions or a convergent
  numerical scheme.
- **Genre/length:** One Markdown file, curriculum-as-questions, tagged `[BSc]`/`[MSc]`, with a
  coverage map. Target ~150-200 items; landed at **158 across 23 parts (54 BSc / 104 MSc).**

## Decisions (AskUserQuestion)

1. **Scope ceiling = comprehensive + relativistic/QFT bridge.** 1D potentials, the harmonic
   oscillator, the TDSE/PDE, 1D and 3D scattering, approximation methods, periodic potentials,
   3D/central/hydrogen, angular momentum, EM coupling, spin coupled to space, identical
   particles, density operators/Wigner, CV quantum optics, open systems, path integrals,
   nonlinear/mean-field, the Klein-Gordon and Dirac equations, and a second-quantization
   bridge. OUT: full interacting QFT, lattice gauge theory, hardware.
2. **Tooling = strictly pure Wolfram Language, no QuantumFramework anywhere.** `DSolve`,
   `NDEigensystem`, `NDSolve`, `FourierTransform`/`Fourier`, `NIntegrate`, `FindRoot`, and
   hand-built schemes (Numerov, split-step Fourier, transfer matrix, imaginary-time descent).
   Continuous objects that must become finite arrays (truncated number basis, spatial grid) are
   truncated explicitly with a convergence check.
3. **Units = natural, $\hbar=m=1$** (and $\omega=1$ where natural). A few items restore
   constants to show a physical scale (oscillator length, Bohr radius).
4. **Size = comparable to the companion** (~150-200; companion is 194/25/83-111).

## The structural difference from the finite-dimensional companion (the disanalogy, named)

The discrete list could settle a claim by *exact* equality of two finite arrays
(`WL===QF`). Continuous space cannot: there is no exact finite matrix in general, the spectrum
can be continuous, an eigenstate can be non-normalizable (delta-normalized). So the
verification idiom shifts to, in order of preference: (1) a closed form when the special
function exists (Hermite, Laguerre, Airy, spherical harmonics, parabolic-cylinder, Coulomb
wave functions), checked by back-substitution into the differential equation; (2) agreement
between two independent methods on the same quantity (analytic spectrum vs `NDEigensystem`,
`DSolve` propagation vs split-step Fourier); (3) numerical convergence under grid/basis
refinement. The representation and discretization machinery is therefore early *content* (1.7
grid discretization, 2.5 delta-normalization, 4.4 honest truncation), not a footnote. This is
stated in the file's intro so the metaphor "companion volume" never papers over it.

## Pressure-test (skeptic's pass, with resolutions baked into the list)

- *Weakest assumption:* "everything is computable in finite closed form." False; resolved by
  making the toolkit early content and shifting verification to method-agreement (above).
- *Strongest objection:* "this is just a Mathematica tutorial on `NDEigensystem`." Answer: each
  item leads with a physics question; continuous QM has genuine conceptual content with no
  finite-dim analog (continuous spectra, non-normalizable states, self-adjointness vs
  Hermiticity, boundary conditions) that is the real subject.
- *Dependency hazards specific to wave mechanics:* Fourier/momentum rep before momentum-space
  work; the harmonic oscillator before anything CV (Wigner, optics, Fock, second quantization);
  angular momentum and central potentials before 3D scattering; spin and EM coupling before the
  Dirac equation; the density operator before reduced states and Wigner; the propagator before
  path integrals. Enforced in the ordering.

## Dependency-ordering rule (hard constraint, mirrors companion DNA 13)

No question may need a concept introduced by a later-numbered question. Tag every item
`[BSc]`/`[MSc]`. One concept per question; inseparable facets stay together.

## Review (adversarial-draft Step 6: two fresh-context reviewers, parallel)

Reviewer A (physics correctness + completeness, quantum-review lens) and Reviewer B (strict
forward-reference audit) ran in parallel on the draft.

- **Reviewer B: 1 genuine forward-reference** = old 4.8 (coherent-state *time evolution*) sat in
  Part 4 but needed time evolution, first introduced in Part 5. **Fixed:** moved to Part 5 as
  5.6 (next to wave-packet revivals, the eigenstate-expansion dynamics); Part 4 keeps only the
  static coherent/squeezed-state construction. Confirmed no `[BSc]`-on-`[MSc]` dependency.
- **Reviewer A findings, applied:** (1) virial typeset $\langle x\cdot\nabla V\rangle\to
  \langle \vec r\cdot\nabla V\rangle$ (8.1); (2) displacement/squeeze operators were named
  formally only at 17.2 but used at 4.6/4.7: self-contained 4.6/4.7 by writing the operators
  inline ($e^{\alpha a^\dagger-\alpha^* a}$, $e^{\frac12(\xi^* a^2-\xi a^{\dagger2})}$) and
  reworded 17.2 to "revisit ... as phase-space transformations"; (3) added the two-body
  center-of-mass / reduced-mass separation (new 10.4, a real BSc gap before hydrogen); (4)
  added homodyne/heterodyne detection (new 17.7, the missing CV-measurement item, motivates the
  Husimi-$Q$ operationally); (5) added the radial-WKB Langer correction (new 12.9); (6)
  generalized 11.6 Clebsch-Gordan beyond the orbital Gaunt instance so 14.3 ($\vec J=\vec
  L+\vec S$) has an in-list foundation; (7) disambiguated 13.4 to the *normal* (orbital) Zeeman
  effect and added the anomalous Zeeman / Paschen-Back item (new 14.6), removing the latent
  spin forward-reference.
- **Soft edges judged non-violations (kept, not "fixed"):** 1.9 Galilean boost uses a c-number
  velocity, not the $\hat p$ operator; 7.7 Fermi-golden-rule density of states is the
  free/scattering continuum (from plane-wave delta-normalization, 2.5), distinct from the band
  DOS at 9.4; 11.5 "Wigner $D$-matrix" is the rotation matrix, not the Wigner quasi-probability
  (16.3).

Post-edit checks: per-part item counts match the coverage table; numbering sequential 1..n in
all 23 parts; tag totals 54 BSc / 104 MSc / 158 total; no em or en dashes; all math in
delimited TeX.

## Status

| Stage | State |
|-------|-------|
| Curriculum question list (`Question-List.md`) | DONE, 158 items / 23 parts, two-reviewer pass applied |
| Worked solutions, Parts 1-2 (`Solutions-Parts-1-2.md`) | DONE, 18 Q, pure WL, kernel-verified + /wl-verify 0 |
| Notebook build (`Solutions-Parts-1-2.nb`, md2nb) | DONE, Template Default, Evaluate->False, 115 KB, structure + code-fidelity verified |
| Worked solutions, Parts 3-23 | NOT STARTED |

## Solutions work log

**Parts 1-2 worked solutions (`Solutions-Parts-1-2.md`), pure WL, hbar=m=1.** Per the user's
recommended order: first fetched the WL doc pages (`wl-docs/`, 17 functions: NDSolve(+Value),
DSolve(+Value), Solve, DEigensystem(+values), NDEigensystem(+values), FourierTransform, Fourier,
Integrate, NIntegrate, FindRoot, Resolve, Reduce, SchrodingerPDEComponent) via the nb-reader
`doc2md.py`; read the Parts 1-2 relevant ones (Integrate definite/generalized-function sections,
FourierTransform in full). The ODE/PDE/eigensystem pages are fetched and will be read when Parts
3+ use them.

**Key correctness point (from the FourierTransform doc):** WL's default
`FourierParameters -> {0, 1}` is the kernel $\frac{1}{\sqrt{2\pi}}\int f\,e^{+i\omega t}dt$, the
OPPOSITE sign from the QM forward transform $\tilde\psi(p)=\frac{1}{\sqrt{2\pi}}\int\psi\,e^{-ipx}dx$.
Using the default would silently flip the sign of $\langle p\rangle$. So 2.6 uses
`FourierParameters -> {0, -1}` and cross-checks $\langle p\rangle=k$ against the position-space
value of 2.3 (both $+k$). Verified independently by the wl-verifier.

**Verification.** Harness `verify/part12.wls` (22/22 PASS) plus a full extract-and-run of the
exact `.md` cells (zero error messages, runs to completion). Caught and fixed: (a) 2.5 plane-wave
delta-normalization is best shown via `FourierTransform[1,x,q] = Sqrt[2Pi] DiracDelta[q]` and the
momentum-eigenstate image `DiracDelta[p-pp]` (the bare `Integrate[e^{i(p-pp)x},...]` does not
converge in this version); (b) a real cell-ordering defect: cells 1.8 and 2.3 had assigned
`\[Psi] = ...` as a bare OwnValue (a product), after which the later `\[Psi][x_] :=` definitions in
2.4/2.6/2.7 failed with `SetDelayed::write Tag Times` (the symbol head was a product). Fixed by
using the function form `\[Psi][x_] :=` consistently. Minimal-code adversary removed 3 inert
`Simplify`/`ComplexExpand`/`FullSimplify` wrappers (cells 1.1, 1.4, 2.6). `/wl-verify` fresh-context
verifier: **OPEN ISSUES 0**, all 18 items reproduced, Fourier convention judged consistent.

**Notebook build (md2nb).** Built `Solutions-Parts-1-2.nb` from the `.md` twin via
`MarkdownToNotebook["Solutions-Parts-1-2.md", "Solutions-Parts-1-2.nb", "Evaluate" -> False]`
(load `Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/MarkdownToNotebook.wl"]` first;
the local `wolframscript` is the standard one, so `-file`/`-code` work, NOT the M2N-sandbox
custom-wrapper flags). `Template: Default` -> StyleDefinitions Default.nb; Title from the body H1.
`Evaluate -> False` keeps Input cells only (no Output), matching the reader-runs-cells design.
Verified by reading back with `Get` (build scripts `verify/build-nb.wls`, `verify/check-nb.wls`,
`verify/check-nb-inputs.wls`): 1 Title, 2 Sections, 18 Subsections, 27 Input, 40 Text, 195
InlineFormula; **0 empty/dropped mathboxes**; `\hbar` rendered to the `\[HBar]` glyph (6x, the M2N
WolframParser path knows it; only the ImportString fallback drops it), `\tfrac`/`\frac` -> 14
FractionBoxes, `\operatorname{erf}` -> visible "erf"; and **all 27 Input cells reparse + run clean**
from the built nb (no `ToExpression::sntx`, no `SetDelayed::write Tag Times`), so md2nb did not
mangle the WL (the `_\[Name]` blank-head trap was avoided by using `\[Psi][x_]` form). Rebuild the
nb whenever the `.md` changes (it is a build output, never hand-edited).

**Non-trivial-states revision (user feedback: "you went with the most trivial answers e.g.
Gaussian; that is not good").** Reworked Parts 1-2 to run on physically important non-Gaussian
states, each chosen so its symbolic answer is interesting, rather than defaulting to the Gaussian
(which is degenerate: it self-Fourier-transforms, saturates uncertainty, has a featureless
current). States and their verified results (`verify/part12b.wls`, all PASS):
- **Cusp state** $\sqrt{\kappa}e^{-\kappa|x|}$ (Dirac-delta-potential bound state / 1D hydrogen 1s):
  density $\kappa e^{-2\kappa|x|}$ (1.1), grid sum resolves the cusp (1.7), moving-state
  $\langle p\rangle=k$, $\langle p^2\rangle=k^2+\kappa^2$, $\Delta p=\kappa$ via the boundary-term-free
  $\int|\partial_x\psi|^2$ form (2.3, sidesteps the cusp's $\delta$ in $\psi''$); Lorentzian
  momentum wavefunction $\sqrt{2/\pi}\kappa^{3/2}/(\kappa^2+p^2)$.
- **Sech soliton** $\tfrac{1}{\sqrt2}\operatorname{sech}(x)$ (NLS bright soliton / Poschl-Teller
  ground): $\int\operatorname{sech}^2=2$ (1.2), $P(|x|<a)=\tanh a$ (1.3), displaced $\Delta x=\pi/(2\sqrt3)$
  (1.5), Hermiticity probes (2.4), FT is a sech $\tfrac{\sqrt\pi}{2}\operatorname{sech}(\pi p/2)$ with
  $\langle p^2\rangle=1/3$ (2.6), and the uncertainty contrast $\Delta x\Delta p=\pi/6\approx0.524>1/2$
  (2.7, vs the Gaussian's exact $1/2$).
- **Hydrogen 1s reduced radial** $2\kappa^{3/2}x e^{-\kappa x}$ on $[0,\infty)$: complex recoil/survival
  amplitude $8\kappa^3/(2\kappa+iq)^3$ (1.4); motivates the half-line in 2.9.
- **Chirped soliton** $\operatorname{sech}(x)e^{i\beta x^2}/\sqrt2$: position-dependent velocity field
  $v(x)=2\beta x$ (1.8).
- **Gaussian**: kept ONLY in 2.7 as the unique saturating state (the punchline).
Caveat: these states live in later Parts (delta potential P3, soliton/PT P8/P21, hydrogen P12);
used here as GIVEN functions (the computation needs only calculus), with a one-line physical-origin
note, so no dependency violation (textbook habit of meeting a non-trivial state early). Full
extract-run of the 28 rewritten cells: zero error messages, reached end. TRAP: the shifted+boosted
sech Hermiticity integral in 2.4 is SLOW (~2 min) but correct.

**Solution-style settled (apply to Parts 3-23):** symbolic inputs where the kernel allows; no helper
functions beyond operator defs (code is the explanation); compute from the definition; result stated
inline in prose as TeX; each cell verified by closed-form / two-method-agreement / numerical
convergence. Always use the function form `name[x_] :=` for states (never bare `name = ...`, which
collides with later `name[x_] :=`). Operators `xop`/`pop` defined once (2.1), reused.

## Next (if the worked manual is wanted)

Author a sibling deliverable giving each question a worked **pure-WL** answer (no QF), to the
companion's answer-style DNA: symbolic wherever the solver allows, no helper functions (code is
the explanation), compute from the definition, honest truncation/convergence caveats stated.
Verification per the three-way idiom above (closed form / two-method agreement / convergence),
then `/wl-verify` for any `.wls` harness. Build the `.nb` from the `.md` twin via
`MarkdownToNotebook[src, out, "Evaluate" -> False]` only if a notebook is requested.
