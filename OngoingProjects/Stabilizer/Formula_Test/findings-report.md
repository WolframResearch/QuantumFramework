# Stabilizer Subsystem — Findings Report

**Date:** 2026-05-07
**Reference:** [`stabilizer-formulas.md`](stabilizer-formulas.md) (extracted from 28 papers under `OngoingProjects/Stabilizer/tex/`)
**Test battery:** [`stabilizer-formulas-test.wls`](stabilizer-formulas-test.wls)
**Paclet under test:** `QuantumFramework` at HEAD of `master`

---

## Executive summary

Running 132 strict `VerificationTest`s mapped 1:1 to the formulas in
`stabilizer-formulas.md` produced **126 passing, 6 failing**. The 6 failures
collapse into **3 distinct findings**, all confirmed by direct kernel probes
that bypass the test harness:

| # | Severity | Where | What |
|---|---|---|---|
| **F1** | High (correctness) | `PauliStabilizer["5QubitCode" \| "SteaneCode" \| "9QubitCode"]` | Returns `|+_L⟩`, not `|0_L⟩` as the API doc claims |
| **F2** | High (correctness) | `ps["M", pauli_string]` | Pauli-string measurement does **not** collapse the state to a +1 eigenspace of the measured Pauli when the operator anti-commutes with the stabilizer; only the sign of an existing stabilizer is flipped |
| **F3** | Medium (contract gap) | `CliffordChannel[ps][ps]` | `nA == 0` (state-prep) channel applied to a same-dim PauliStabilizer emits `CliffordChannel::stateevol` and falls back to dense; API doc says this dispatch case is supported |

The remaining 126 tests covering Pauli encoding, the full Heisenberg gate
table, Bell / GHZ / cluster states, MacWilliams identities, the sum-of-
stabilizers projection, inner products, the symplectic ↔ Clifford isomorphism,
local complementation, the Choi-tableau composition `S² = Z`, and the 13-item
"Stabilizer Olympics" all pass.

---

## F1 — Named codes return `|+_L⟩`, not `|0_L⟩`

**Severity:** High (correctness — affects every downstream amplitude check).

**Failing tests (4):**

| TestID | Expected | Actual |
|---|---|---|
| `S3-5QubitCode-Got97-16termSum-uptoPhase` | `True` | `False` |
| `S3-SteaneCode-LogicalZero-8amps-STRICT` | `8` | `16` |
| `S7-Steane-CSS-Codeword-8terms-STRICT` | `8` | `16` |
| `S24-SteaneCode-LogicalZero-8Amps-STRICT` | `8` | `16` |

**Root cause.** The named-code constructors in
[`Kernel/Stabilizer/NamedCodes.m`](../../../QuantumFramework/Kernel/Stabilizer/NamedCodes.m)
list the **logical X̄** as the n-th stabilizer instead of the **logical Z̄**.
For `[[5,1,3]]`:

```
PauliStabilizer["5QubitCode"]["Stabilizers"]
= {XZZXI, IXZZX, XIXZZ, ZXIXZ, XXXXX}      ← XXXXX is logical X̄
```

The **+1 eigenstate** of the *full* listed group (4 stabilizers ∪ {X̄}) is by
construction `|+_L⟩ = (|0_L⟩ + |1_L⟩)/√2`, not `|0_L⟩`. To return `|0_L⟩` the
generator set must include logical Z̄ = ZZZZZ instead.

**Direct kernel reproduction:**

```wolfram
Needs["Wolfram`QuantumFramework`"];
ps = PauliStabilizer["5QubitCode"];
v  = Normal[ps["State"]["StateVector"]];

(* Logical Z̄ is NOT a stabilizer of the named-code state *)
zzz = KroneckerProduct @@ ConstantArray[{{1, 0}, {0, -1}}, 5];
Chop[zzz . v - v]               (* = nonzero -- v is NOT in +1 eigenspace of Z̄ *)

(* The state IS in +1 eigenspace of XXXXX *)
xxx = KroneckerProduct @@ ConstantArray[{{0, 1}, {1, 0}}, 5];
Chop[xxx . v - v]               (* = 0 *)

(* And |⟨5QubitCode | 5QubitCode1⟩| = 0 (5QubitCode1 lists -XXXXX) *)
v1 = Normal[PauliStabilizer["5QubitCode1"]["State"]["StateVector"]];
Chop[Conjugate[v] . v1]         (* = 0 *)
```

The two named states span the orthogonal `±1` eigenspaces of X̄, i.e.
`|+_L⟩` and `|-_L⟩` — not `|0_L⟩` and `|1_L⟩`.

**Verifying the underlying math is correct.** Projecting `|+_L⟩` onto the +1
eigenspace of Z̄ (i.e., applying `(I + Z̄)/2`) recovers the textbook Got97 §3.7
16-term `|0_L⟩` exactly:

```wolfram
proj0 = (IdentityMatrix[32] + zzz) / 2;
v0 = proj0 . v;  v0 = v0 / Sqrt[Conjugate[v0] . v0];
Position[v0, x_ /; Abs[x] > 1*^-10] // Flatten
(* {1, 4, 6, 7, 10, 11, 13, 16, 18, 19, 21, 24, 25, 28, 30, 31}  -- exactly Got97 §3.7 *)
```

So the **state machinery is correct**; only the choice of n-th generator in the
named code is wrong.

**Affected codes (and their *Code1* siblings):** `"5QubitCode"`,
`"7QubitCode"` / `"SteaneCode"`, `"9QubitCode"`.

**Suggested fix.** In `NamedCodes.m`, replace the logical-X̄ generator with
logical-Z̄ for each of these named codes, *or* update the API doc and
`*Code1` semantics:

| Current behavior | Documented intent | Fix option A | Fix option B |
|---|---|---|---|
| `5QubitCode` ⇒ `|+_L⟩` | `|0_L⟩` | replace `XXXXX` with `ZZZZZ` | rename to `"5QubitCode_PlusL"` |
| `5QubitCode1` ⇒ `|-_L⟩` | `|1_L⟩` | replace `-XXXXX` with `-ZZZZZ` | rename to `"5QubitCode_MinusL"` |

(I lean toward Fix A — the API docstring "|0_L⟩ of [[5,1,3]] (Got97 §3.5)" is
the principled contract, and Got97 §3.5 itself defines `|0_L⟩` as the +1
eigenstate of the four real stabilizers + Z̄.)

---

## F2 — `ps["M", pauli_string]` does not collapse the state on random outcomes

**Severity:** High (correctness — multi-qubit Pauli measurement is broken).

**Failing test:**

| TestID | Expected | Actual |
|---|---|---|
| `S5-PostMeas-State-IsEigenstate-STRICT` | `True` | `False` |

**Symptom.** After `ps["M", pauli_string]` with a *random* outcome (the
measured operator anti-commutes with one or more stabilizers), the
post-measurement state should be the +1 (or -1) eigenstate of `pauli_string`
per AG §3 Case I. The simulator instead returns a state that:

1. Has the same support as the *input*, not the projected support.
2. Has expectation value **0** for the measured Pauli, not ±1.
3. Leaves the anti-commuting stabilizer in place with only its sign flipped,
   instead of *replacing* it with `±pauli_string`.

The single-qubit form `ps["M", q_integer]` works correctly — only the
Pauli-string path is broken.

**Direct kernel reproduction (Bell |Φ+⟩, measure ZI):**

```wolfram
ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];

(* CORRECT path: single-qubit form *)
ps["M", 1]
(* <|0 -> ps with stabs={ZI, ZZ}, state={1,0,0,0},
    1 -> ps with stabs={-ZI, ZZ}, state={0,0,0,1}|>  -- collapsed to |00⟩ or |11⟩ *)

(* BUGGY path: Pauli-string form *)
ps["M", "ZI"]
(* <|0 -> ps with stabs={XX, ZZ}, state={1/√2, 0, 0, 1/√2},
    1 -> ps with stabs={-XX, ZZ}, state={-1/√2, 0, 0, 1/√2}|>  -- still |Φ+⟩ / |Φ-⟩ *)

(* Direct expectation values on the buggy path *)
v0 = Normal[ps["M", "ZI"][0]["State"]["StateVector"]];
zi = KroneckerProduct[{{1, 0}, {0, -1}}, IdentityMatrix[2]];
Conjugate[v0] . zi . v0          (* = 0  -- should be +1 *)
```

The same bug shows for `XI`, `IX`, `IZ`, etc. (any single-Z / single-X
embedded in a multi-qubit Pauli string), and is reproducible across every
2-qubit stabilizer state we tried.

**Root cause analysis.** The Pauli-string measurement path in
[`Kernel/Stabilizer/PauliMeasure.m`](../../../QuantumFramework/Kernel/Stabilizer/PauliMeasure.m)
appears to perform only step (4) of AG §3 Case I (toggle the sign of
non-commuting generators), without step (3) (replace the chosen
anti-commuting generator with the measurement operator). The single-qubit
path in [`Measurement.m`](../../../QuantumFramework/Kernel/Stabilizer/Measurement.m)
performs both steps and is correct.

**Suggested fix.** In `PauliMeasure.m`, after rowsumming all but one
anti-commuting stabilizer, the chosen anti-commuting generator should be
overwritten with the bit pattern of `pauli_string` and a random sign
(matching what `Measurement.m` does for `Z_q`).

**Workaround.** For users who need a multi-qubit Pauli measurement, conjugate
to the Z basis first:

```wolfram
(* Measuring ZI ≡ Z₁ — already works via the integer path *)
ps["M", 1]

(* Measuring XI = (H_1)(Z_1)(H_1) — push through Heisenberg conjugation *)
ps["H", 1]["M", 1]   (* outcome 0 ↔ XI = +1; outcome 1 ↔ XI = -1 *)

(* Measuring an arbitrary Pauli P: rotate to Z, measure, rotate back *)
```

---

## F3 — `CliffordChannel[ps][ps']` falls back to dense for state-prep + dim-matched arg

**Severity:** Medium (contract gap; documented case unimplemented).

**Failing test:**

| TestID | Expected | Actual |
|---|---|---|
| `S16-CliffordChannel-StatePrep-Yields-Zeros-STRICT` | `{"ZI", "IZ"}` | wrapped `QuantumChannel[…]` (dense fallback) |

**Symptom.** A state-prep channel built from a `PauliStabilizer` (so
`InputQubits = 0`, `OutputQubits = n`) applied to a same-dim `PauliStabilizer`
should — per the API doc — *return the prepared state*. Instead the kernel
emits

```
CliffordChannel::stateevol: cc[ps]: only identity / matching-input cases
  are implemented in Phase 8.2.
```

and falls back to wrapping the result in a dense `QuantumChannel[…][…]`
expression.

**Direct kernel reproduction:**

```wolfram
cc = CliffordChannel[PauliStabilizer[2]];
{cc["UA"], cc["UB"], cc["c"], cc["InputQubits"], cc["OutputQubits"]}
(* {{}, {{0,0,1,0}, {0,0,0,1}}, {0,0}, 0, 2}  -- a valid nA=0 prep channel *)

cc[PauliStabilizer[2]]   (* fires CliffordChannel::stateevol, returns dense *)
```

**API contract that's currently violated** (from `Documentation/api.md`):

> ### State evolution: `cc[ps]`
> Apply the channel to a `PauliStabilizer` state. Three recognized cases:
>
> - **Identity channel** (`Source -> "Identity"`, dimensions match): return ps unchanged.
> - **State-prep channel** (`nA == 0`): return the state encoded by cc.
> - **Dim-matched channel**: build `CliffordChannel[ps]`, compose, convert back to `PauliStabilizer`.

The middle row ("State-prep channel" with `nA == 0`) doesn't fire in the
dispatch tree — falling through to the `CliffordChannel::stateevol` fallback.

**Suggested fix.** In
[`Kernel/Stabilizer/CliffordChannel.m`](../../../QuantumFramework/Kernel/Stabilizer/CliffordChannel.m)
(around line 451 per the API doc), add a dispatch arm:

```wolfram
cc_CliffordChannel[ps_PauliStabilizer] /;
    cc["InputQubits"] === 0 :=
        cliffordChannelToPauliStabilizer[cc]
```

or equivalent — the `cliffordChannelToPauliStabilizer[cc]` helper is already
listed in the api.md as a PackageScope.

---

## What passed (the 126 successful tests)

For completeness, the following formulas all reproduce correctly in the
current QF kernel:

- §0 Pauli matrix identities, group cardinalities, eigenstate table (10/10)
- §1 xz / GF(4) encodings, symplectic predicate vs. matrix commutator (5/5)
- §2 Heisenberg conjugation table for H, S, S†, X/Y/Z, CNOT, CZ, SWAP (11/11)
- §3 Bell, GHZ, cluster, code stabilizer *generator-set* outputs (13/15) — only
       the |0_L⟩ amplitude-level checks fail (F1 above)
- §4 Density-matrix and projector identities (3/3)
- §5 Single-qubit measurement (`["M", q]`), all 5-qubit / Steane stabilizer
       measurements deterministic (7/8) — only the Pauli-string path fails (F2)
- §6 Inner products, AarGot04 |⟨ψ|φ⟩|² = 2⁻ˢ (7/7)
- §7 Sum-of-stabilizers projection (Bell + 5-qubit) (2/3) — Steane amplitude
       count downstream of F1 (1/3)
- §8 CSS structure, transversal-CNOT closure on Steane (2/2)
- §9 MacWilliams `A_d = (1, 0, 0, 0, 15, 0)` for [[5,1,3]] (1/1)
- §10 Hamming, Singleton, GV bounds (4/4)
- §11 T-gate magic states, StabilizerFrame component count (5/5)
- §12 Symplectic ↔ Clifford isomorphism, MJM^T = J, M⁻¹ = JM^T J (3/3)
- §13 BSM, teleportation 4-outcome, two-Bell fusion (3/3)
- §14 Local complementation (3/3)
- §15 Random-Clifford ↔ canonical-form round-trips (2/2)
- §16 Identity channel, state-prep tableau structure, S² = Z composition (3/4)
       — only state-prep evolution fails (F3)
- §17 Qudit sanity (1/1)
- §18 Pauli tracking through CNOT (3/3)
- §19 Universal one-line identities — H², S², S⁴, (HS)³, (CNOT)², (CZ)²,
        H_b·CNOT·H_b = CZ, Bell expectations, Toffoli sanity (10/10)
- §20 Tableau ↔ State ↔ Operator round-trips, Y-canary (Phase 5c) (6/6)
- §21 Stim-style identities (3/3)
- §22 Stim format normalization (1/1)
- §23 Caveats: Y = iXZ, equivalent stabilizer sets, W ∉ stabilizer (3/3)
- §24 Concrete numeric fixtures: |+⟩, GHZ-3 amplitudes (2/3) — Steane
        amplitude count downstream of F1 (1/3)
- §25 Stabilizer Olympics 13-item conformance (13/13)

---

## Recommended remediation order

1. **F2 first** — it's the most user-visible (`ps["M", "ZI"]` on a Bell
   state silently returns wrong post-states; this *will* break any test or
   protocol that does multi-qubit Pauli measurement).
2. **F1 second** — easy fix in `NamedCodes.m`; high impact on every test
   that compares against textbook codeword amplitudes (Got97 §3.7,
   Got00 §4, Nielsen & Chuang 10.5).
3. **F3 last** — the most contained; users have an obvious workaround
   (`ps_state["State"]` directly, or just don't compose state-prep
   channels with same-dim arguments).

---

## Reproduction

From a fresh kernel:

```bash
wolframscript -f OngoingProjects/Stabilizer/Formula_Test/stabilizer-formulas-test.wls
```

Expected output: `132 tests, 126 pass, 6 fail`. The 6 failures cluster into
the 3 findings above. Each `TestID` matches its row in the formula reference
[`stabilizer-formulas.md`](stabilizer-formulas.md), so a fix can be verified
by re-running the affected tier alone:

```bash
# After applying a fix, narrow the run via grep on TestID strings:
wolframscript -f OngoingProjects/Stabilizer/Formula_Test/stabilizer-formulas-test.wls 2>&1 \
    | grep -E "S3-|S5-PostMeas|S7-Steane|S16-Clifford|S24-Steane"
```
