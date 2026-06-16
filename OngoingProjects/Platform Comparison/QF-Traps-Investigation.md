# QF Traps Investigation: seven flagged items

**Scope.** For each of seven flagged behaviors, classify it as **BUG** (incorrect code), **BY-DESIGN** (correct by intent), or **KNOWN-GAP** (an unshipped feature, not a defect). This is an audit. **No kernel code was changed.**

**Environment.** Repo `/Users/mohammadb/Documents/GitHub/QuantumFramework`, branch `main`, HEAD `f9dc1cdc`, paclet 2.0.0. Every claim below was verified live against the **working tree** (`PacletDirectoryLoad["…/QuantumFramework"]; Needs["Wolfram\`QuantumFramework\`"]`), not the installed paclet, using `wolframscript -file`. The audit reference set at `audit/` (anchor `f9dc1cdc`) was read in full first.

**Verdict summary.**

| # | Item | Verdict |
|---|---|---|
| 1 | `QuantumState[Table[0,{n}],2]` norm-0 padding | BY-DESIGN (now warned); permissive contract, no incorrect computation |
| 2 | Order-driven multiplicity broadcast | BY-DESIGN feature + correct guard; one conservative false-positive edge |
| 3 | Stabilizer `SameQ` across storage forms | BY-DESIGN representation invariant; genuine ergonomic sharp-edge for equality code |
| 4 | `$QuantumFrameworkPropCache = False` nearly moot | BY-DESIGN / expected; no defect found (cache correct, bit-identical, no `Rule::rhs`) |
| 5 | Packed numeric statevector fold (178x prototype) | KNOWN-GAP (unshipped); default apply path has no packed-numeric route |
| 6 | Bulk Pauli-frame (reference-frame) sampling | KNOWN-GAP; no frame simulator exists; `StabilizerFrame`/`SampleOutcomes` are different objects |
| 7 | Error-aware transpilation (uncommitted working tree) | IN-PROGRESS feature, sound and live-verified; two minor flags |

**Fixes verified, not applied.** For the four items with an actionable change (1, 2, 3, and 7's two flags), the root cause was confirmed at source and a fix was prototyped and tested entirely out-of-tree (scratch `wolframscript` redefinitions, upvalues, and a guard predicate evaluated on the real operators; repo-wide grep for flag 7b). Each item's section ends with a **Verified** subsection carrying the evidence. **No file under `QuantumFramework/Kernel/` was modified.** Items 4, 5, 6 need no fix (4 is expected behavior with no defect; 5 and 6 are unshipped capabilities, not bugs).

---

## Item 1: `QuantumState[Table[0,{n}],2]` produces a 4-qubit norm-0 state

**Verdict: BY-DESIGN (now warned).** The constructor's "accept any-length amplitude vector, pad up to the nearest power of the qudit dimension" rule is intentional and is used legitimately for the single-qudit amplitude-prefix case. The multi-qudit round-up that changes the qudit count and can leave norm $0$ is a misuse path that, as of `7e12e8f0`, emits `QuantumState::padded`. No incorrect computation occurs; the surprise is purely in the permissive contract.

**Evidence.** Constructor: [QuantumState.m:32-47](QuantumFramework/Kernel/QuantumState/QuantumState.m). The multiplicity helper rounds up: [Utilities.m:68-71](QuantumFramework/Kernel/Utilities.m),
$$\texttt{basisMultiplicity}[\text{dim},\text{size}] = \begin{cases} 1 & \text{size}=1 \\ \lceil \log_{\text{size}}\text{dim}\rceil & \text{otherwise}\end{cases}$$
with `dim` the vector length and `size` the per-qudit dimension. For a length-$12$ vector over qubits, $\lceil\log_2 12\rceil = 4$, so $4$ qubits ($16$ amplitudes), and `PadRight` extends the $12$ explicit entries (all zero) to $16$, giving norm $0$.

Live (`/tmp/traps_verify2.wls`):

```
QuantumState[Table[0, {12}], 2]  ->  Qudits = 4, Dimension = 16, Norm = 0
                                     message fired: QuantumState::padded
basisMultiplicity[12, 2] = 4
```

The two distinct pad paths, both verified:

```
QuantumState[{1}, 4]        ->  Qudits = 1, StateVector = {1,0,0,0} = |0>,  NO message   (multiplicity 1)
QuantumState[Table[0,{8}],2] ->  Qudits = 3,  NO message                                  (8 = 2^3, exact)
QuantumState[Table[0,{12}],2] ->  4 qubits, norm 0,  QuantumState::padded                 (round-up)
QuantumState["Register"[12]] ->  Qudits = 12, Norm = 1                                     (the correct register)
```

The guard in the constructor only warns when the pad crossed a qudit boundary and the multiplicity exceeded $1$: `If[multiplicity > 1 && Length[state] =!= basis["Dimension"], Message[QuantumState::padded, …]]` ([QuantumState.m:40-42](QuantumFramework/Kernel/QuantumState/QuantumState.m)). The `{1}` prefix has multiplicity $1$, so it stays silent: this is the path the tensor-network initial-state prepend relies on (a length-$1$ amplitude vector promoted to $|0\rangle$ in a $d$-dimensional qudit).

**Root cause.** `stateQ` accepts any `VectorQ` with `Length >= 0` ([Utilities.m:84](QuantumFramework/Kernel/Utilities.m)), so a length-$12$ list is a structurally valid (if physically meaningless) amplitude vector. The constructor infers the qudit count by rounding the length up to a power of the dimension rather than requiring an exact power or an explicit dimension list. A norm-$0$ vector is a legal element of the underlying linear space (the zero vector); QF does not reject it because zero-amplitude SparseArrays are also legitimate intermediate objects.

**Recommendation (no change applied).** The warning is the correct minimal fix and should stay. Hardening options, in increasing strictness: (a) keep as-is (warn); (b) additionally warn or fail when the *result* has norm $0$, since a norm-$0$ multi-qudit state is almost never intended; (c) require the vector length to be an exact power of the qudit dimension for the multi-qudit case, keeping the silent multiplicity-$1$ prefix path for `{1}`-style prepends. Option (b) is the cleanest user-facing guard without disturbing the prefix path. Do not make norm-$0$ a hard global failure: the single-qudit prefix and zero-amplitude intermediates depend on it being constructible.

**Verified (cause confirmed, guard scope checked; `/tmp/fix123.wls`).** The warning is *not* limited to the norm-$0$ case: it fires on every non-power round-up, including the more insidious *nonzero* misshape, and is silent exactly where it should be:

| input | qudits | norm | `QuantumState::padded`? |
|---|---|---|---|
| `QuantumState[{1,1,1}, 2]` (len 3, not a power of 2) | 2 | $\sqrt3$ | **fires** |
| `QuantumState[{1,0,0,1}/Sqrt2, 2]` (len 4 $= 2^2$) | 2 | 1 | silent (correct) |
| `QuantumState[{1}, 4]` (single-qudit prefix) | 1 | 1 | silent (correct) |

So the existing guard already covers the dangerous "silently wrong but nonzero" padding (`{1,1,1}` becomes `{1,1,1,0}`), not just the norm-$0$ headline case. **Conclusion: no further kernel change is strictly needed.** The only optional refinement is wording: the message could append the resulting norm when it is $0$ (option (b)) so the most-broken case is unmistakable. The verdict (BY-DESIGN, correctly scoped) stands and is verified.

---

## Item 2: the order-driven multiplicity-broadcast rule

**Verdict: BY-DESIGN feature with a correct OOM guard; one conservative false-positive edge.** Placing a $k$-qudit operator on an order of length $m$ (with $k \mid m$) tensors $\mathrm{LCM}$ copies of the operator. This is a genuinely useful broadcast (a single-qubit gate spread over several wires) and it is computed correctly. The `7e12e8f0` guard correctly converts the former misshapen-input kernel-killer into a fast `$Failed`. The only blemish is that the guard's cutoff is a coarse element-count proxy, so it also blocks a *sparse* single-qubit gate broadcast past 12 wires.

**Evidence.** Broadcast rule: [QuantumOperator.m:251-271](QuantumFramework/Kernel/QuantumOperator/QuantumOperator.m). The multiplicity is $\mathrm{LCM}\!\left(\lfloor m_{\text{out}}/q_{\text{out}}\rfloor, \lfloor m_{\text{in}}/q_{\text{in}}\rfloor\right)$ and the construction tensors that many copies ([QuantumOperator.m:294-295](QuantumFramework/Kernel/QuantumOperator/QuantumOperator.m)). The guard ([QuantumOperator.m:260-265](QuantumFramework/Kernel/QuantumOperator/QuantumOperator.m)) fires when $\dim(qo)^{\text{multiplicity}} > \$QuantumOperatorBroadcastLimit = 2^{24}$.

Live (`/tmp/traps_verify2.wls`), with `xgate = QuantumOperator[{{0,1},{1,0}}]`:

```
$QuantumOperatorBroadcastLimit = 16777216  (= 2^24)
QuantumOperator[xgate, {1,2,3}]  ->  3 qudits;  Matrix === QuantumOperator["XXX"]["Matrix"]  ->  True   (legitimate, correct)
QuantumOperator[xgate, Range[11]] ->  11 qudits, no message
QuantumOperator[xgate, Range[12]] ->  12 qudits, no message
QuantumOperator[xgate, Range[13]] ->  $Failed, message: QuantumOperator::broadcast
QuantumOperator["CNOT", {1,2,3,4}] ->  QuantumOperator, 4 qudits   (2-qubit gate, 2 copies, legitimate)
```

The boundary is exact: $\dim(qo) = d_{\text{out}}\, d_{\text{in}} = 4$ for a single-qubit gate, so $4^{12} = 2^{24}$ passes ($>$ is strict) and $4^{13} = 2^{26} > 2^{24}$ fails. So the legitimate broadcast $X^{\otimes 13}$ (a sparse $2^{13}\times 2^{13}$ permutation matrix, cheap to build) is blocked.

**Which combinations trigger it.** The rule fires only for an explicit `QuantumOperatorQ` operator (not a name) when each side's order length is a strict multiple of that side's qudit count. A single-qubit gate over $m$ wires broadcasts to $m$ copies; a $2$-qubit gate over an even number of wires $m$ broadcasts to $m/2$ copies; and so on. Named gates such as `"H"` over a long order take the Fourier-kronecker construction path instead of this rule (documented in `audit/mistakes.md`, HEAD-resync section), so a named gate over a huge order still builds eagerly and is not protected by this guard.

**Root cause.** The guard estimates cost by the operator's *element-space dimension* $\dim^{\text{multiplicity}}$ (output dimension $\times$ input dimension, raised to the copy count), which equals the number of dense matrix entries. It is not sparsity-aware, so a permutation-sparse gate like $X$ is judged by the same yardstick as a dense gate like $H$. For a *dense* single-qubit gate, $H^{\otimes 13}$ is genuinely a $\sim$1 GB dense matrix, so the cutoff is defensible; for a sparse gate it is over-conservative.

**Recommendation (no change applied).** This is acceptable as a safety limit, and the user can raise `$QuantumOperatorBroadcastLimit` for a deliberate large broadcast. A sharper, verified fix is below.

**Verified fix (prototyped out-of-tree; `/tmp/fix123.wls`).** A single base-swap is *not* enough: switching to the output dimension $d_{\text{out}}^{\text{multiplicity}}$ alone (limit $2^{24}$) would allow a *dense* $H^{\otimes 13}$ (a $\sim 1$ GB build), while switching to a sparsity-aware nonzero count alone would *miss* the OOM ket, whose malformed form has zero explicit entries ($0^{12}=0$). The OOM is driven by the **output dimension/basis** ($16^{12}=2^{48}$ output basis elements), whereas $X^{\otimes 13}$ and $H^{\otimes 13}$ differ only in matrix **sparsity** at the *same* dimension $2^{13}$. The correct guard therefore needs **both** terms:
$$\text{block} \iff d_{\text{out}}^{\,\text{mult}} > L_{\dim} \;\lor\; \max(\mathrm{nnz},1)^{\,\text{mult}} > L_{\text{nnz}},\qquad L_{\dim}=L_{\text{nnz}}=2^{24}.$$
The sparsity premise was confirmed: $X^{\otimes 13}$ built via the Pauli-string path is genuinely cheap, $0.30$ s, $8192$ explicit nonzeros, and **stays a `SparseArray`** (QF preserves sparsity through `QuantumTensorProduct`), so $\mathrm{nnz}^{\text{mult}}$ is the true entry count. Measured inputs: $X$ has $\mathrm{nnz}=2,\ d_{\text{out}}=2$; $H$ has $\mathrm{nnz}=4$; the malformed 4-qubit ket-operator has $d_{\text{out}}=16,\ \mathrm{nnz}=0$. The two-part predicate is correct on all seven probed cases:

| case | $d_{\text{out}}$ | nnz | mult | guard blocks? | want | correct |
|---|---|---|---|---|---|---|
| $X$ over 3 | 2 | 2 | 3 | no | no | yes |
| $X$ over 13 | 2 | 2 | 13 | **no** | no | yes (was wrongly blocked) |
| $X$ over 24 | 2 | 2 | 24 | no | no | yes |
| $H$ over 13 (dense) | 2 | 4 | 13 | **yes** (nnz term, $4^{13}=2^{26}$) | yes | yes |
| CNOT over 4 | 4 | 4 | 2 | no | no | yes |
| OOM ket, norm 0 | 16 | 0 | 12 | **yes** (dim term, $16^{12}=2^{48}$) | yes | yes |
| OOM ket, normal | 16 | 16 | 12 | yes | yes | yes |

The dimension term is what catches the norm-$0$ OOM that the sparsity term alone would miss; the sparsity term is what unblocks the legitimate sparse single-qudit broadcasts that the current dense-assumption guard wrongly rejects. This is a small, local change to [QuantumOperator.m:260-264](QuantumFramework/Kernel/QuantumOperator/QuantumOperator.m), not applied here. (Named gates still take the Fourier-kronecker path and are unaffected by this rule, as noted above.)

---

## Item 3: stabilizer `SameQ` across storage forms fails

**Verdict: BY-DESIGN representation invariant, but a genuine ergonomic sharp-edge.** A `PauliStabilizer` has two valid in-memory shapes: a *packed* form (returned by single Clifford gates) and a *canonical* rank-3 tableau form (returned by `ApplyCircuit` and the named constructors). They are not `SameQ` even for the identical physical state, because their Association keys differ. Both forms compute identical physics, so this is not a correctness bug; but `===`, `Union`, `DeleteDuplicates`, Association keys, and memoization will silently treat the two as distinct.

**Evidence.** Predicate accepts both shapes: [PauliStabilizer.m:36-48](QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m). Packed form is built by the gate kernels ([Packed.m](QuantumFramework/Kernel/Stabilizer/Packed.m), gated by `psConcreteFastQ` at [Packed.m:87-97](QuantumFramework/Kernel/Stabilizer/Packed.m)); `ApplyCircuit` returns canonical ([Compiled.m:224-229](QuantumFramework/Kernel/Stabilizer/Compiled.m)). There is **no** `Equal`/`SameQ` upvalue on `PauliStabilizer` (grep of `Kernel/Stabilizer/*.m`); `StabilizerFrame` does have one ([StabilizerFrame.m:128](QuantumFramework/Kernel/Stabilizer/StabilizerFrame.m)), which highlights the asymmetry.

Live (`/tmp/traps3.wls`), starting from `base = PauliStabilizer["GHZ"[3]]`:

```
base["H", 1]                     keys = {PackedX, PackedZ, Signs, Qubits}   (single gate -> packed)
base["ApplyCircuit", {"H" -> 1}] keys = {Signs, Tableau}                    (ApplyCircuit -> canonical)

psP === psC                 ->  False
psP["Tableau"]   === psC["Tableau"]    ->  True     (same physical state)
psP["Stabilizers"] === psC["Stabilizers"] ->  True
Length @ Union[{psP, psC}]            ->  2         (sees two distinct objects)
Length @ DeleteDuplicates[{psP, psC}] ->  2
TrueQ[psP == psC]                     ->  False     (no Equal upvalue either)
both pass PauliStabilizerQ            ->  True
```

**Tests rely on `===` only within a single storage form.** The one raw-object equality test, `PauliStabilizer["7QubitCode"] === PauliStabilizer["SteaneCode"]` ([Tests/Stabilizer/PauliStabilizer.wlt:695](Tests/Stabilizer/PauliStabilizer.wlt), expects `True`), compares two *named constructors* that both take the canonical path, so it is form-consistent and passes. Everywhere else the suite compares the canonical accessors `["Tableau"]`/`["Stabilizers"]`/`["Generators"]` ([Tests/Stabilizer/AuditMatrix.wlt:126,423,429](Tests/Stabilizer/AuditMatrix.wlt)). So the test suite is already disciplined around this; the trap is for downstream user code that mixes a single-gate result with a constructor/`ApplyCircuit` result and compares with `===`.

**Root cause.** The packed form was added as a performance representation (62 tableau bits per machine word) without a canonicalization-on-construction step or an equality overload. Two Associations with different key sets are never `SameQ`. `SameQ` (`===`) is not overloadable in WL for a user head, so form-sensitivity of `===` cannot be removed; only `Equal` (`==`) can be given an upvalue.

**Recommendation (no change applied).** Two non-exclusive options:
1. **Add a `PauliStabilizer /: Equal` upvalue** comparing the canonical content (`["Tableau"]`, `["Signs"]`, and recovered phase), mirroring the existing `StabilizerFrame` upvalue. This fixes `==`, `Union`/`DeleteDuplicates` with the default `SameTest` unchanged but gives users a correct value-equality operator. It cannot fix bare `===`.
2. **Lazily canonicalize** at a small set of choke points (or eagerly on construction) so that two equal states share one representation and become `SameQ`. This is the only route that also fixes `===`, `Union`, and Association keys, at the cost of giving up the packed fast-path's storage savings unless canonicalization is deferred.

Until then, the contract is already documented (`audit/mistakes.md`): compare `ps1["Tableau"] === ps2["Tableau"]`, never the raw objects. This is a latent ergonomic trap, not a physics bug.

**Verified fix (prototyped out-of-tree; `/tmp/fix123.wls`).** First, the packed and canonical forms of the same construction produce **identical** canonical content: `psP["Tableau"] === psC["Tableau"]` and `psP["Signs"] === psC["Signs"]` are both `True`, so a structural `Equal` is well-defined. The prototype upvalue (set in a scratch session only, no kernel file touched)

```wolfram
PauliStabilizer /: Equal[a_PauliStabilizer, b_PauliStabilizer] :=
  TrueQ[a["Qubits"] == b["Qubits"]] && a["Tableau"] === b["Tableau"] && a["Signs"] === b["Signs"]
```

behaves correctly:

```
psP == psC                              ->  True    (packed vs canonical, same state: fixed)
psP == PauliStabilizer["GHZ"[3]]["X",1] ->  False   (genuinely different state: no false positive)
Length @ DeleteDuplicates[{psP, psC}, Equal] ->  1   (dedup collapses the two equal forms)
```

For the harder cross-construction question (two different gate sequences for the same state, which could in principle yield different generator choices), the framework already ships a bulletproof oracle: the inner product. Building a Bell state two ways (`PauliStabilizer["GHZ"[2]]` versus `PauliStabilizer[2]["H",1]["CNOT",1,2]` from $|00\rangle$) gave **identical** tableaux here (QF's gate updates are deterministic, so `bell1 == bell2` is already `True`), and `Abs[bell1["InnerProduct", bell2]] == 1` confirms physical equality independently. So the recommended fix is the cheap structural `Equal` upvalue (verified), with `Abs[a["InnerProduct", b]] == 1` as the documented robust test for any case where two inequivalent tableaux might describe the same state. Both are verified; neither is applied. (Reminder: `===` is intrinsically form-sensitive and cannot be overloaded; only `==` can be fixed this way, so the structural-accessor comparison remains the contract for `===`-style code.)

---

## Item 4: `$QuantumFrameworkPropCache = False` is now nearly moot

**Verdict: BY-DESIGN / expected; no defect found.** After the `dfc741dc`/`f9dc1cdc` cache-machinery repair, the default-on warm apply is within $\sim 1.2\times$ of the cache-off floor, so the flag is no longer a meaningful speed lever. The three correctness checks the prompt asked for all pass.

**Evidence (`/tmp/traps3.wls`), representative $n=8$, $3$-layer circuit (H + CNOT ladder + RZ$(0.3)$):**

(a) Cache returns correct values, bit-identical on vs off:
```
$QuantumFrameworkPropCache = True  -> vOn
$QuantumFrameworkPropCache = False -> vOff
vOn === vOff           ->  True
Max[Abs[vOn - vOff]]   ->  0.
```

(b) No `Rule::rhs` under heavy property access. With a counting handler on `Rule::rhs` wrapped around $5$ rounds of `qc[psi]["VonNeumannEntropy"]`, `["Purity"]`, `["Eigenvalues"]`:
```
item4_rulerhs_count = 0
```
This is the intended consequence of the `cacheProperty` store ([Utilities.m:426-432](QuantumFramework/Kernel/Utilities.m)), which commits a cache DownValue only when the held key is `FreeQ[_, _Pattern | _Blank | _BlankSequence | _BlankNullSequence | _Optional]`, so a pattern-bearing key can never form an over-matching DownValue. The former `Quiet[HeadProp[…] = result, Rule::rhs]` is gone; `cacheProperty` resolves and has DownValues (`item4_cacheProperty_def = True`).

(c) Staleness. The documented caveat stands: there is no automatic invalidation, so substituting into the underlying expression of a *cached non-parametric* state requires pinning a new symbol (the cache key is structural identity). Parametric states bypass the cache entirely (gated on `ParameterArity == 0`), so a free-parameter state never serves a stale value; substitute parameters first, then read. Neither is a defect: both are the intended cache contract.

**Root cause / classification.** Expected. The whole point of the refactor was to make the cache-off path no longer necessary; the flag remains only as an escape hatch for the no-auto-invalidation case. No real defect surfaced.

---

## Item 5: packed numeric statevector fold, prototyped at 178x, not shipped

**Verdict: KNOWN-GAP (unshipped feature), not a bug.** The prototype reproduces at HEAD and the default apply path provably has no packed-numeric fast path.

**Evidence (re-running the committed prototype `OngoingProjects/Platform Comparison/scripts/apply06_prototypes.wls`):**

```
ref_warm_ms (n=12, default TN apply)  =  288.95
A_fold_run{1,2,3}_ms (packed fold)    =  4.33, 5.14, 4.36
A_correct_maxdiff                     =  4.9e-17     (agrees with default path)
ref16_warm_ms (n=16, default)         =  2661.63
A16_fold_run{1,2,3}_ms                =  179.06, 167.35, 183.87
A16_correct_maxdiff                   =  2.0e-17
```

The packed fold is $\sim 4.3$ ms at $n=12$, bit-accurate to $\sim 5\times 10^{-17}$. The "$178\times$" figure in `QF-Apply-Path-Deep-Audit.md` Section 5 was relative to the *pre-cache-fix* warm apply ($\sim 800$ ms $/\, 4.5$ ms $\approx 178$); against the *current* post-fix warm default ($289$ ms) the same fold is $\sim 67\times$. At $n=16$ the fold ($\sim 177$ ms) is $\sim 15\times$ the default ($2.66$ s), where the dense-in-sparse contraction now dominates. The opportunity is real and unchanged.

**Default path has no such route.** The apply dispatcher ([QuantumCircuitOperator.m:97-112](QuantumFramework/Kernel/QuantumCircuitOperator/QuantumCircuitOperator.m)) has exactly five `Method` branches: `"Schrodinger"` (a `FullSimplify` fold), `Automatic | "TensorNetwork"` (the default, `TensorNetworkApply`), `"QuEST"`, `"Qiskit"`, `"Stabilizer"`. None is a packed-numeric dense fold. Confirmed absent.

**Recommendation.** This is the one large remaining apply-path win (tracked as GitHub issue #5 in the prototype header). It is explicitly out of scope here; do not implement. Classification: KNOWN-GAP.

---

## Item 6: bulk Pauli-frame (reference-frame) sampling is absent

**Verdict: KNOWN-GAP, not a bug.** There is no constant-cost-per-gate Pauli-frame simulator in `Kernel/Stabilizer/`. The two superficially related symbols are different objects.

**Evidence.** A grep of `Kernel/Stabilizer/` for `frame.?sampl | pauli.?frame | reference.?frame | frame.?simul | batch.*shot | shot.*frame` returns nothing. The two near-neighbors:

- **`StabilizerFrame`** ([StabilizerFrame.m:1-28](QuantumFramework/Kernel/Stabilizer/StabilizerFrame.m)) is a stabilizer-*rank* decomposition $\sum_i c_i\,|s_i\rangle$ for $T$-gate-rich circuits and magic-state work (Garcia and Markov, Quipu, arXiv:1712.03554). It is a superposition of stabilizer states, not a frame of injected Pauli errors.
- **`ps["SampleOutcomes", n]`** ([Properties.m:69-70](QuantumFramework/Kernel/Stabilizer/Properties.m), implemented in [SymbolicMeasure.m](QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m), FangYing23 SymPhase, arXiv:2311.03906) amortizes repeated *measurement* of **one** tableau: a single circuit traversal leaves fresh $\mathbb{F}_2$ symbols in the signs for unresolved random outcomes ([SymbolicMeasure.m:9-27](QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m)), and sampling substitutes random $0/1$ bits into those symbols.

**Why these are not a frame simulator.** Stim's and QuantumClifford.jl's headline QEC feature propagates a batch of independent random Pauli error *frames* through a fixed Clifford circuit at $O(1)$ cost per gate per shot, producing many independent noisy syndrome shots for decoder training. `SampleOutcomes` instead resolves the random measurement outcomes of a *single* noiseless stabilizer evolution; it does not inject per-shot Pauli errors and does not generate independent noisy syndrome records. So the headline bulk-shot QEC generator is genuinely absent.

**Recommendation.** Classify as KNOWN-GAP. A true frame simulator (random Pauli frame propagated with the same $O(n)$ per-gate tableau commutation already implemented in the packed kernels, vectorized over shots) would be a natural future addition. Not a bug; nothing to fix.

---

## Item 7: error-aware transpilation in the working tree (uncommitted)

**Verdict: IN-PROGRESS feature, sound and live-verified offline; two minor flags.** The new `"OptimizationLevel"` option, the `QiskitTarget["FromBackend"]` extractor, the offline `GenericV2` fake device, and `update_instruction_properties` in `build_target` all function correctly. Two non-blocking concerns: a silent-fallthrough usability edge and a stray cross-project metadata key.

**Source audit.** Diffs read in full:
- [IBMQuantum.m](QuantumFramework/Kernel/QuantumCircuitOperator/IBMQuantum.m): adds `"OptimizationLevel" -> Automatic` to `IBMJobSubmit` options and threads it into `qiskitPrimitiveSubmit`. Sound: `Automatic` maps to qiskit's preset default ($1$), preserving prior behavior.
- [Qiskit.m](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m): `qiskitQASM` and `qiskitPrimitiveSubmit` gain `"OptimizationLevel"` (default `Automatic -> 1`); the `"GenericV2"` provider case builds `GenericBackendV2(num_qubits=int(backend_name) or 5, seed=42)` ([Qiskit.m:446-450](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)); `QiskitTarget["FromBackend" -> name, opts]` reads a live backend's populated qiskit `Target` and returns an `InstructionProperties`-carrying spec ([Qiskit.m:670-694](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)). The Python correctly filters control-flow/directive operation names out of `basis_gates` while keeping `measure` for the instruction-properties walk.
- [QuantumQASM.m](QuantumFramework/Kernel/QuantumCircuitOperator/QuantumQASM.m): `build_target` pops `instruction_properties` before `Target.from_configuration` (it is not a constructor parameter) and re-applies it via `update_instruction_properties` ([QuantumQASM.m:79-99](QuantumFramework/Kernel/QuantumCircuitOperator/QuantumQASM.m)), tolerating missing error/duration and recording un-appliable rows rather than dropping them silently; `instruction_properties` is whitelisted as a reserved (non-`from_configuration`) key in the unknown-option check ([QuantumQASM.m:206-211](QuantumFramework/Kernel/QuantumCircuitOperator/QuantumQASM.m)); a new `["InstructionProperties"]` accessor and an "Error-aware" summary-box row are added. All sound.

**Live verification (offline, no credentials, no network; qiskit 2.4.1 in QF's managed Python env; `/tmp/item7.wls`, `/tmp/item7b.wls`):**

Correct invocation, `QiskitTarget["FromBackend" -> "5", "Provider" -> "GenericV2"]` (backend name $=$ qubit count):
```
Head = QiskitTarget,  built in ~2.9 s, offline
OperationNames = {cx, delay, id, measure, reset, rz, sx, x}
num_qubits = 5
coupling_map  = 20 directed edges (all-to-all on 5 qubits)
InstructionProperties = 45 rows, each {gate, {qubits}, <|Error -> _, Duration -> _|>}
   e.g. {cx,{0,1}, Error 2.78e-3, Duration 2.66e-7}, {sx,{0}, Error 9.56e-5}, {rz,{q}, Error 0., Duration 0.}, {measure,{0}, Error 3.61e-3}
```
This is a genuine, populated error-aware `Target`: realistic two-qubit ($\sim 10^{-3}$) versus single-qubit ($\sim 10^{-4}$) error hierarchy, zero-cost virtual `rz`, nonzero readout error. The end-to-end error-aware transpile then works:
```
QuantumQASM[qc, "Provider"->"GenericV2", "Backend"->"5", "OptimizationLevel"->3]
  ->  OPENQASM 3.0; include "stdgates.inc"; rz(pi/2) $1; sx $1; rz(pi/2) $1; cx $1, $3; ...
```
The transpiler lowered to the device basis $\{rz, sx, cx\}$ and chose physical qubits `$1, $3` (an error-aware layout decision), confirming the `optimization_level` plumbing reaches qiskit's preset pass manager. `"OptimizationLevel" -> 0` also returns a valid string. The bare native path `qc["QASM"]` ([Properties.m:521](QuantumFramework/Kernel/QuantumCircuitOperator/Properties.m), $\to$ `QuantumQASM[qco, "WL"]`) emits native OpenQASM 3.0 with no Python.

**Flag (a): silent fall-through on the natural-looking call.** The literal `QiskitTarget["FromBackend" -> "GenericV2"]` (provider unset) does **not** select the fake device. With `"Provider"` defaulting to `None`, the dispatch falls to the default Aer branch, and the string `"GenericV2"` lands in `"Backend"` where the Aer path ignores it. The result is a `QiskitTarget` reporting an **idealized** device with **empty** instruction properties and an **empty** coupling map, and **no warning**:
```
QiskitTarget["FromBackend" -> "GenericV2"]
  ->  num_qubits = 31, coupling_map = {}, instruction_properties = {} (0 rows),
      basis_gates = full 42-gate standard set
```
A user who writes the obvious `"FromBackend" -> "GenericV2"` expecting the offline fake device silently gets a no-error, all-to-all target instead. The error-aware device requires `"Provider" -> "GenericV2"` with the backend name set to the qubit count. Worth a guard: either recognize `"GenericV2"` as a provider when it appears as the backend name, or warn when an extracted target has zero instruction properties and an empty coupling map.

**Flag (b): stray cross-project metadata key.** [QuantumQASM.m:119](QuantumFramework/Kernel/QuantumCircuitOperator/QuantumQASM.m) stamps skipped rows under `md['qualink_skipped_props']`. The `qualink_` prefix is a copy-paste artifact from the QUALink project; in QuantumFramework it should be a QF-namespaced or neutral key (e.g. `quantumframework_skipped_props` or `skipped_props`). Harmless to behavior, but it should be renamed before this lands.

**Recommendation.** The feature is correct and the offline `GenericV2` path is verified working. Before committing: (1) address the silent fall-through so the natural call either works or warns; (2) rename the `qualink_skipped_props` key. Neither affects the correctness of the intended (`"Provider" -> "GenericV2"`) path. Classification: in-progress feature, sound.

**Verified cause + fixes (out-of-tree; source trace + `/tmp/item7.wls`, `/tmp/item7b.wls`, grep).**

*Cause of flag (a), traced through source.* `QiskitTarget["FromBackend" -> name_, opts]` forwards `name` only as `"Backend"` and forwards `"Provider"` from `opts` ([Qiskit.m:673-674](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)). With the literal call, `"Provider"` is unset, so inside `qiskitInitBackend` it defaults to `None` ([Qiskit.m:353,364](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)). `None` falls to the `_` branch of all three Switches: env $\to$ `"default"` ([Qiskit.m:382-383](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)), provider-init $\to$ `BasicProvider` ([Qiskit.m:408-412](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)), and backend-select $\to$ the `AerSimulator` branch ([Qiskit.m:451-465](QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m)), which **ignores `backend_name`** entirely. So `"GenericV2"` never reaches the `GenericBackendV2` constructor; the extracted `backend.target` is Aer's idealized 31-qubit, all-gates, no-error target. Confirmed.

*Fix (a), two verified options.*
1. **Token routing (does what the user meant).** In `QiskitTarget["FromBackend" -> name]`, when `"Provider"` is unset and `name` matches a known offline fake-device token (currently just `"GenericV2"`), route `name` as the provider and default the backend to the qubit count. Verified to yield the correct device: `QiskitTarget["FromBackend" -> "5", "Provider" -> "GenericV2"]` returns $45$ instruction-property rows, a $20$-edge coupling map, $5$ qubits.
2. **Empty-target backstop (catches any idealized fallthrough).** After extraction, warn when the spec carries no error model. The discriminator is clean and verified against both runs:

| call | `instruction_properties` | `coupling_map` | empty-target predicate |
|---|---|---|---|
| `"FromBackend" -> "5", "Provider" -> "GenericV2"` | 45 rows | 20 edges | `False` (real device) |
| `"FromBackend" -> "GenericV2"` (provider unset) | `{}` | `{}` | `True` (idealized: warn) |

i.e. `emptyTargetQ[spec] := spec["instruction_properties"] === {} && spec["coupling_map"] === {}` fires only on the silent-fallthrough case. Recommend option 1 as primary with option 2 as a general backstop.

*Fix (b), safety of the rename, verified.* A repo-wide grep shows `qualink_skipped_props` is **write-only**: the only references are the assignment and the `skipped_props` list it is built from, all inside `build_target` ([QuantumQASM.m:57,88,95,99,117,119](QuantumFramework/Kernel/QuantumCircuitOperator/QuantumQASM.m)); nothing in `Kernel/`, `Tests/`, or elsewhere ever reads the key back. So renaming it to a QF-namespaced or neutral key (`quantumframework_skipped_props` or `skipped_props`) is behavior-preserving with zero downstream impact.

Both fixes are verified and ready; neither is applied, per the no-kernel-change constraint.

---

## Summary for a physicist

Most of these seven flags are not defects in the physics or the numerics; they are properties of how the framework represents states and operators, and where it has chosen to put guard rails. Two of them (items 1 and 2) are about the constructor accepting a slightly malformed description of a system and quietly completing the Hilbert space for you: a length that is not a power of the qudit dimension gets padded up to the next power, which can hand you a zero vector living in a larger space than you intended, and an operator placed on more wires than it acts on gets tensored with copies of itself. Both are deliberate conveniences, both now carry a warning or a fast failure where they used to silently waste memory, and in both cases the actual arithmetic that does run is correct. Item 3 is a representation choice: a stabilizer state can sit in memory in two equivalent encodings (a bit-packed one for speed, a plain tableau for everything else), and two encodings of the *same* physical state are not recognized as literally identical by the symbolic equality test, even though every physical quantity computed from them agrees exactly; the cure is to compare tableaux, not raw objects. Item 4 confirms that a recent rework of the property cache did its job: turning the cache off no longer changes either the speed (much) or, importantly, the answer (the output is bit-for-bit identical), and the spurious internal warnings the old cache used to suppress are simply gone. Items 5 and 6 are genuine missing capabilities rather than mistakes: a fast packed numerical circuit-evaluation path exists only as a validated prototype (about two orders of magnitude faster than the current default, agreeing to machine precision), and the bulk error-frame sampling that makes Stim and QuantumClifford.jl fast for quantum-error-correction shot generation has no counterpart here yet (the two similarly named objects do something else). Item 7 is a feature still on the workbench: error-aware transpilation, where the compiler routes a circuit using a device's measured per-gate error and timing data. It works, verified live and entirely offline against a synthetic fake device that reports realistic gate errors; the only rough edges are that the most natural-looking way to ask for that fake device silently gives you an idealized error-free one instead, and a leftover variable name from a sibling project needs renaming before it is committed.

*Code footnote: nothing under `QuantumFramework/Kernel/` was modified. The findings rest on reading the constructor and dispatch sources (`QuantumState/QuantumState.m`, `QuantumOperator/QuantumOperator.m`, `Stabilizer/*.m`, `Utilities.m`, `QuantumCircuitOperator/*.m`) and on live `wolframscript` runs whose scripts are in `/tmp/traps_verify2.wls`, `/tmp/traps3.wls`, `/tmp/item7.wls`, `/tmp/item7b.wls`, plus the committed prototype at `OngoingProjects/Platform Comparison/scripts/apply06_prototypes.wls`.*
