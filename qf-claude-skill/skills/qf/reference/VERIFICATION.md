# Verification record

The skill's behavioral claims and examples were independently re-derived against the live kernel,
not taken on trust from the prior audit. This file records what was checked, how, and what was
corrected.

| Field | Value |
|---|---|
| Date | 2026-06-19 (delta below); 2026-06-17 (base record) |
| Paclet | `Wolfram/QuantumFramework 2.0.0` |
| Kernel | working tree at commit `b4fcdb31` (loaded via `PacletDirectoryLoad`); base record was `a601a187` |
| Method | fresh `wolframscript` batteries (no `Quiet` hiding results), direct reads of the cited `Kernel/*.m` source, and a cross-check of this `reference/` set against the maintainer's `audit/` source |

## 2026-06-19 delta: multi-gate bracket dispatch (kernel `a601a187` → `b4fcdb31`)

One kernel-touching commit landed since the base record, `c3f21a43` (multi-gate bracket forms on
`PauliStabilizer` / `StabilizerFrame`; `Stabilizer/PauliStabilizer.m`, `Stabilizer/StabilizerFrame.m`,
`Stabilizer/Compiled.m`). No `PacletInfo` `Symbols` change. The new claims were re-derived against the
working-tree kernel and the cites read from source:

- **List = variadic = per-gate fold.** `PauliStabilizer[2][{"H","CNOT"}]`, `PauliStabilizer[2]["H","CNOT"]`,
  and `Fold[#1[#2]&, PauliStabilizer[2], {"H"->1,"CNOT"->{1,2}}]` give identical `"Tableau"` and `"Signs"`.
- **Bare-name default targets** (one-qubit → qubit 1, two-qubit → `{1,2}`) equal the explicit-target list.
- **`ps[{"X"}]` applies X; `ps["X"]` is the property accessor** (`{_List,_List}`); `ps[{1,2}]` stays a
  measurement (`AssociationQ` True); a single bare `ps["H"]` stays inert; an unrecognized two-token call
  (`ps["Wobble","Splork"]`) is rejected with **no** message (`PauliStabilizer::badgate` is unreachable
  from the public path).
- **Catalogs** `$stabilizerCliffordNames` = `{H,S,X,Y,Z,V,T,CNOT,CX,CZ,SWAP}`,
  `$stabilizerTwoQubitNames` = `{CNOT,CX,CZ,SWAP}` (`Compiled.m:57-58`).
- **Frame parity + the `args : ___Integer` fix.** `PauliStabilizer[2]["T",1]` is a `StabilizerFrame`;
  both `sf[{"H","S"->2}]` and `sf["H","S"->2]` reproduce `sf["H",1]["S",2]["StateVector"]` to amplitude 0.
  The single-gate rules at `StabilizerFrame.m:201, 212` carry `args : ___Integer` (a trailing gate name is
  no longer mis-read as a qubit index); helpers at `PauliStabilizer.m:141-168`, frame multi-gate at `:263-266`.
- **Test suites** re-run on the working tree: `Tests/Stabilizer/AuditMatrix.wlt` **216/216**,
  `Tests/Stabilizer/PauliStabilizer.wlt` **336/336**, all green. A fresh adversarial reviewer reproduced
  all of the above and the file:line cites in a clean kernel: **OPEN ISSUES: 0**.

The base 2026-06-17 record (everything below) is unchanged and still applies; the line numbers it cites are
point-in-time for `a601a187` (only the stabilizer files above moved).

## Method

Every falsifiable claim in `SKILL.md` was turned into a `wolframscript` check run against the
loaded paclet, and every cited behavior was confirmed by reading the kernel `.m` source. The
skill's own rule applies to the skill itself: line numbers are point-in-time (see
[VERSION.md](VERSION.md)); the *behavior* notes are the durable part.

## Coverage

- **`SKILL.md`: fully behaviorally verified.** Every §2 pitfall, every §3 construction example
  (states, operators, circuits, stabilizers), the §4 / §5 / §11 kernel cites, the §6 Wigner
  parity dimensions, the §8 nonexistent-symbol claims, the §10 property counts, the §12
  expectation-value-four-ways and composition idioms, and the catalog counts.
- **`reference/` set: cross-checked + spot-verified.** Confirmed in sync with the `audit/` source
  (the same content, minus the maintainer-only banners and notes). The QuantumDistance and
  property-count entries in `mistakes.md` / `function-index.md` were already correct; the
  VonNeumann-entropy entry carried the same defect as `SKILL.md` and was corrected in both. The
  ~90 individual `mistakes.md` entries were not each re-run end to end.
- **Not reachable offline:** claims that need live hardware (IBM QPU) or the Python/qiskit bridge
  (§7 hardware patterns, parts of the interop hub) were not executed; they are read-confirmed
  against the kernel source only.

## Confirmed correct (representative)

| Claim | Check |
|---|---|
| 2.1 list-form gate specs misparse silently | `QuantumOperator[{"R",θ,"XY"}]` builds a `d=3` operator; `{"RX",θ}` a `d=2` one |
| 2.2 `"SX"` is 2-qubit, `"V"` is √X | arity 2 vs 1; `V·V − X = 0` |
| 2.2 lowercase alias dead; `"Z"[d]` inverse clock | `QuantumOperator["x"]` → `Failure`; `Z[3]` diagonal `{1, ω⁻¹, ω⁻²}` |
| 2.7 `MatrixInverse` silent pseudo-inverse | `{{1,1},{1,1}}` → `{{1/4,1/4},{1/4,1/4}}`, no message |
| 2.8 Magic injects `i` on `\|10⟩` | returns `{0, i/√2, i/√2, 0}` |
| 3.x all named constructors evaluate | states 9/9, operators, circuits, stabilizers all build |
| catalog counts | `$QuantumOperatorNames` 108, `$QuantumCircuitOperatorNames` 44, `$PauliStabilizerNames` 9 |
| 6 Wigner parity | odd `d=3` → `3×3`; even `d=2` → `4×4` |
| 8 nonexistent symbols | `QuantumSchmidtDecomposition` / `…Spectral…` / `QuantumCompile` stay unevaluated; `qs["SchmidtDecompose"]`, `qs["SpectralBasis"]`, `qco["Compile"]` work |
| 10 property split | `QuantumOperator["X"]`: 67 `"Properties"` ⊂ 209 `"AllProperties"` |
| 12 expectation four ways | all four agree to `0` on a random `⟨Z⟩` |

## Corrections applied (this verification)

1. **§2.3 VonNeumann entropy / purity guard: narrowed.** The blanket "`Indeterminate` for ≥10
   qubits" was wrong: both handlers fast-path pure states (`If[qs["PureStateQ"], 0, …]`;
   `If[StateType==="Vector", 1, …]`), so a *pure* state of any size returns `0` / `1`. The
   `ConfirmAssert[qs["Dimension"] < 2^10]` guard fires only on the *mixed-state* branch.
   `QuantumState["GHZ"[10]]` returns entropy 0 / purity 1, not `Indeterminate`. Cites corrected
   `:461,:481` (and `:455,:475` in `mistakes.md`) → `:460, :480`. Fixed in `SKILL.md`,
   `reference/mistakes.md`, and the `audit/mistakes.md` source.
2. **§2.4 / §12 QuantumDistance "Fidelity": corrected.** `"Fidelity"` **is** the default measure
   (`QuantumDistance.m:13`) and returns the *distance* `1 − Re Tr(ρ₁ρ₂)^(1/2)` (`= 1 − |⟨a|b⟩|`
   for pure), not the fidelity. The old "for raw fidelity call `QuantumDistance[…,"Fidelity"]`"
   was self-contradictory. Use `QuantumSimilarity[a,b,"Fidelity"] = 1 − d = |⟨a|b⟩|` for the
   overlap, or `|⟨a|b⟩|²` for the squared fidelity. (Verified: `|0⟩,|+⟩` → distance `1−1/√2`,
   similarity `1/√2`, true `F = 1/2`.) The fuller `mistakes.md` / `function-index.md` entries were
   already correct; only `SKILL.md` carried the error (§2.4 and the §12 idiom), now fixed.
3. **Cite drift: corrected.** `Utilities.m` references in `SKILL.md` were off by one:
   `MatrixPartialTrace` `:135 → :136`, `eigensystem` `:154 → :155`, `MatrixInverse`
   `:309-312 → :310-313`.
4. **Lindblad caveat: added (verified).** `QuantumOperator["Hamiltonian"[H, {Ls}, {γs}]]` /
   `"Liouvillian"` wrap `H` automatically but require the **jump operators to be
   `QuantumOperator`s**; bare matrices throw `QuantumOperator::dim` ("…must have the same
   dimensions") because the dimension assert calls `["OutputDimension"]` on the unwrapped jump ops
   (`NamedOperators.m:900-901`). Pass raw matrices to `QuantumEvolve` instead.

## Note for users on a different paclet version

These results are exact for `2.0.0 @ a601a187`. On another build, re-run the behavioral checks
(the `SKILL.md` examples are copy-pasteable) before trusting an exact line number; the qualitative
behavior is stable across patch releases.
