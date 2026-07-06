# Answer Key: Quantum in Finite Dimensions (WL and QuantumFramework)

Per-part answer key. Each `Part-0N.md` is a self-contained source that builds to `Part-0N.nb` via md2nb.
Every question carries **two** answers that agree by exact equality (WL `===` QF), kernel-verified.

## Files
- `Part-0N.md` / `Part-0N.nb`: one Part each (Parts 1-10 authored; Part 11+ pending).
- `_check.sh` / `_check-one.wls`: verify + rebuild driver.
- Shared spine (parent folder): `Question-List.md` (curriculum, 25 parts, coverage map) and
  `PIPELINE.md` (the design DNA: read this before revising).
- `../_archive/`: the pre-split monolith `Quantum-in-Finite-Dimensions-WL-and-QF-Answers.md`/`.nb`
  (kept for reference; safe to delete once the split is trusted).

## Working on one part
Edit `Part-0N.md`, then run the driver:

```
./_check.sh Part-06.md     # fresh-kernel verify every wl cell + rebuild Part-06.nb
./_check.sh                # all parts
```

Per part it runs every `wl` cell in a fresh kernel (flags any that error or emit a message), then
rebuilds the sibling `.nb`. A clean run prints `errored/messaged={}`.

## Conventions (full list in `../PIPELINE.md`)
- **WL** cells are native Wolfram Language only (no `QuantumState`/`QuantumOperator`/... in a WL cell);
  **QF** cells use the paclet's objects and property downvalues.
- The two answers must agree by exact equality where the prose claims it.
- Default state = the generic Bloch **mixed** state
  `1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])`, $|\vec r|\le1$; use a **pure** state
  $\{\cos\tfrac\theta2, e^{i\varphi}\sin\tfrac\theta2\}$ only when the concept requires purity. Write
  $\vec\sigma$ as `PauliMatrix[{1, 2, 3}]` in every $\vec r\cdot\vec\sigma$ / $\hat n\cdot\vec\sigma$ dot
  product; keep bare `PauliMatrix[k]` for a single named observable.
- Answers are **independent**: no answer reuses a variable defined in another (a one-time `Setup` load is
  the only shared thing). Cross-references are part-qualified ("Part 3, 3.4").
- No em dash or en dash; valid delimited TeX for all math; complex conjugate as `^*` not overline.

## Deferred work (pre-verified in-kernel, ready to apply)
- **Part 6**: none; the whole observable-generalization batch is applied. 6.1 and 6.2 use the generic pair
  $A=\hat n_1\cdot\vec\sigma$, $B=\hat n_2\cdot\vec\sigma$ via `FromSphericalCoordinates` (6.1 slack
  $|\hat n_1\times\hat n_2|^2(1-|\vec r|^2)$ with the inequality confirmed by a `Reduce` violation search;
  6.2 gap $=$ squared covariance, saturation $(\hat m\cdot\hat n_1)(\hat m\cdot\hat n_2)=\hat n_1\cdot\hat n_2$);
  6.4 derives the unbiased-basis form and proves maximality; 6.4b is the collapse-to-$\{1/d\}$ test; 6.5
  defines $c$ precisely and verifies by global `NMinimize`; 6.6 derives the saturating states.
- **Part 3**: 3.5 matrix-entry state $\to$ Bloch form (optional; the entries already are the 3 parameters).
- **Parts 7-10** authored and verified (`/wl-verify` OPEN ISSUES 0, one round). Part 7 (20 Q, unitary
  dynamics: split 7.1/7.5/7.9/7.10 into a/b), Part 8 (7 Q, spin), Part 9 (7 Q, SU(2): split 9.1 into a/b),
  Part 10 (4 Q, circuits: split 10.2 into a/b). Heavy MSc items (Fermi golden rule, Landau-Zener, Berry
  phase, Solovay-Kitaev) drop to a discriminating concrete instance where a closed symbolic form is
  intractable, with the caveat stated. Spin operators are built from the ladder relation, not the QF
  named `"JX"`/`"JY"` operators (which return the diagonal `Jz` at HEAD; a kernel-fix task was spun off).
- **Part 11+** (entanglement, mixed states, channels, ...): not yet authored.
