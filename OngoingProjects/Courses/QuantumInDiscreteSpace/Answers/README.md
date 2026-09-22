# Answer Key: Quantum in Finite Dimensions (WL and QuantumFramework)

Per-part answer key to the curriculum in `../Question-List.md`. Each `Part-NN.md` is a self-contained
source that builds to `Part-NN.nb` via md2nb. Every question carries **two** answers, one in native
Wolfram Language (**WL**) and one through the QuantumFramework paclet (**QF**), and where the prose
claims agreement the two agree by exact equality, kernel-verified.

## Files
- `Part-NN.md` / `Part-NN.nb`: one Part each, numbered as in `../Question-List.md`.
- `IDIOMS.md`: the minimal-idiom catalog (preferred WL and QF construction per task, anti-patterns to
  flag in review) and its audit log.
- `AI-Actionable-Audit.md`: a review backlog for Parts 1-10 dated 2026-07-24. It is an audit only; its
  items have not been applied.
- `_check.sh` / `_check-one.wls`: verify + rebuild driver.
- Shared spine (parent folder): `../Question-List.md` (curriculum, 25 parts, coverage map) and
  `../PIPELINE.md` (the design DNA: read this before revising).
- `../_archive/`: the pre-split monolith `Quantum-in-Finite-Dimensions-WL-and-QF-Answers.md`/`.nb`
  (kept for reference; safe to delete once the split is trusted).

## Status

| Part | Theme | Questions answered | Author |
|------|-------|--------------------|--------|
| 1 | Pure states and the Born rule | 10 of 10 | Mads Bahrami |
| 2 | Observables, spectra, approximation | 16 of 17 (2.14 open) | Mads Bahrami |
| 3 | States as density operators | 10 of 10 | Mads Bahrami |
| 4 | Composite systems: tensor product, partial trace | 2 of 2 | Mads Bahrami |
| 5 | Projective measurement | 4 of 4 | Mads Bahrami |
| 6 | Uncertainty and incompatibility | 7 of 7 | Mads Bahrami |
| 7 | Unitary dynamics and pictures | 20 of 25 (7.17 to 7.21 open) | Mads Bahrami |
| 8 | Spin and angular momentum | 7 of 7 | Mads Bahrami |
| 9 | Single-qubit operations, SU(2) | 7 of 7 | Mads Bahrami |
| 10 | Elementary multi-qubit gates and circuits | 4 of 4 | Mads Bahrami |
| 11 | Composite systems and entanglement | 14 of 14 | Bruno Tenorio |
| 12 | Mixed states: distinguishability, thermal | 9 of 9 | Bruno Tenorio |
| 14 | Quantum information and entropy | 8 of 8 | Bruno Tenorio |
| 17 | Tomography, estimation, metrology | 4 of 4 | Bruno Tenorio |

Parts 13, 15, 16, and 18 to 25 have no answers yet. The questions still open inside answered Parts
(2.14 and 7.17 to 7.21) were added to the curriculum after those Parts were written.

The notebooks of Parts 1 to 10 are unevaluated builds (Input cells only, `"Evaluate" -> False`); the
notebooks of Parts 11, 12, 14, and 17 were saved with their evaluated outputs. Running the driver on a
Part rebuilds its notebook without outputs, so the driver skips a Part whose notebook carries Output
cells unless it is given `--force`.

## Working on one part
Edit `Part-NN.md`, then run the driver:

```
./_check.sh Part-06.md           # fresh-kernel verify every wl cell + rebuild Part-06.nb
./_check.sh                      # all parts whose notebook carries no evaluated outputs
./_check.sh --force Part-11.md   # rebuild even though Part-11.nb carries evaluated outputs
```

Per part it runs every `wl` cell in a fresh kernel (flags any that error or emit a message), then
rebuilds the sibling `.nb`. A clean run prints `errored/messaged={}` and exits 0; a flagged cell makes
the driver exit 1 and leaves the existing notebook untouched. The driver evaluates cells with
`ToExpression`, so a cell that refers to the previous output through `%` is not checked: `%` there is
the empty `Out[0]`, and the cell passes without computing anything. Bind a result to a name instead.
The driver also checks only that cells run clean, not that they return the claimed values; a cell that
returns an inert expression passes. Read the outputs of any cell you change. The driver loads the paclet
with `Needs` (the Parts carry no Setup cell); without that load every framework cell evaluates inert and
passes unchecked.

## Conventions (full list in `../PIPELINE.md`)
- **WL** cells are native Wolfram Language only (no `QuantumState`/`QuantumOperator`/... in a WL cell);
  **QF** cells use the paclet's objects and property downvalues. Each answer labels its two routes
  **WL** : and **QF** : before the first cell of each, and closes with one paragraph on what both routes
  give and what it means.
- The two answers must agree by exact equality where the prose claims it. Compare a framework matrix
  with a hand-built one as `Normal[obj["Matrix"]] == m`; the postfix form `obj["Matrix"] // Normal == m`
  parses as `obj["Matrix"] // (Normal == m)` and returns an inert expression, not a Boolean.
- Default state = the generic Bloch **mixed** state
  `1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])`, $|\vec r|\le1$; use a **pure** state
  $\{\cos\tfrac\theta2, e^{i\varphi}\sin\tfrac\theta2\}$ only when the concept requires purity. Write
  $\vec\sigma$ as `PauliMatrix[{1, 2, 3}]` in every $\vec r\cdot\vec\sigma$ / $\hat n\cdot\vec\sigma$ dot
  product; keep bare `PauliMatrix[k]` for a single named observable.
- Answers are **independent**: no answer reuses a variable defined in another. Cross-references are
  part-qualified ("Part 3, 3.4") and point backward only.
- One computation per cell, a bridge sentence before each cell, objects shown (assign without `;` so the
  summary box displays) with a matrix read from the object as one representation.
- No em dash or en dash; valid delimited TeX for all math; complex conjugate as `^*` not overline.
