# QuantumFramework Apply-Path Deep Audit

**Scope.** A Wolfram-Language-native performance audit of the default circuit-evaluation hot path of the Wolfram QuantumFramework paclet (`qc[psi]`, TensorNetwork method), with a function-level cost decomposition, the WL-mechanism explanation for each cost, ablations separating hypotheses, and prototype-validated fix assessments. This supersedes the first-pass profile in section 12 of `QF-Master-Platform-Comparison.md` and materially corrects it.

**Environment.** Apple M2 Pro, Mathematica 15.0, paclet 2.0.0 loaded from the working tree (`PacletDirectoryLoad`). The original audit was at repo HEAD `6f953ad5`; kernel citations `QuantumFramework/Kernel/<path>:line` are at that SHA unless a status note says otherwise. All scripts live in `scripts/apply0*.wls` next to this document and emit `RESULT|label|value` lines; timings are `AbsoluteTiming` wall times in ms. Standard benchmark: $n = 12$ qubits, 3 layers of $H^{\otimes n}$, a CNOT chain, and $R_Z(0.3)^{\otimes n}$, 105 gates total, applied to the register $|0\rangle^{\otimes 12}$.

---

## Status update (2026-06-14, HEAD `f9dc1cdc`): the dominant finding has been fixed in the kernel

This audit was written as a diagnosis-and-recommendation document. Between 2026-06-13 and 2026-06-14 the kernel team **acted on its top finding**. A four-commit property-cache refactor (`7e12e8f0`, `5b696969`, `dfc741dc`, `f9dc1cdc`) implemented what this report calls **Fix D** (Section 5), and the OOM guard from Section 3.7 was committed. The numbers below are the original audit's; this banner records what is now shipped, re-measured at `f9dc1cdc` on the same machine.

| Audit item | Recommendation | Status | Measured outcome |
|---|---|---|---|
| **Fix D** (§5, §3.1): stop paying the property-cache machinery per hit | direct `"Basis"` getter + exclude $O(1)$ structural getters from the cache-eligibility guard; replace `Quiet[...=result, Rule::rhs]` with a pattern-safe store | **SHIPPED** (`dfc741dc`, `f9dc1cdc`) | warm apply $832 \to 273\,\mathrm{ms}$; cold $2956 \to 783\,\mathrm{ms}$; `op["Basis"]`$\times105$ $32 \to 2\,\mathrm{ms}$; output bit-identical |
| §3.1 mechanism: `Quiet`-masked `Rule::rhs` over-match in the cache store | replace with a key-pattern-checked store (`cacheProperty`) | **SHIPPED** (`f9dc1cdc`) | no `Rule::rhs` under heavy access (verified); pattern-bearing keys left uncached |
| §3.5 / Memoize: `AppendTo`-DownValue cache trips `Rule::rhs` on patterned args | rewrite `Memoize` to a Hash-keyed `Association` + a `"CacheQ"` filter | **SHIPPED** (`5b696969`) | `Memoize` now usable for the `*Prop` accessors; integer cache key, never a pattern |
| §3.7 / §0.1: OOM multiplicity-broadcast on a misshapen state; silent zero-pad | guard the broadcast rule with a dimension message; warn on the vector-pad constructor | **SHIPPED** (`7e12e8f0`) | `QuantumOperator::broadcast` + `$Failed` over `$QuantumOperatorBroadcastLimit` ($2^{24}$); `QuantumState::padded` warning on the misparse. Verified live |
| **Fix A** (§5): packed numeric statevector fold | strict-predicate packed fold in `TensorNetworkApply`, fallback otherwise | **NOT SHIPPED** | prototype stands at $178\times$ ($n{=}12$, $4.5\,\mathrm{ms}$); the one large win still open |
| **Fix C** (§5): reusable compiled form across applies | cache the extracted gate list / network on the circuit | **NOT SHIPPED** | prototype $22\times$ for repeat-apply at $n{=}8$ |
| **Fix B** (§5): memoized parametric-gate matrices | (rescoped) closed-form named rotations | **NOT SHIPPED** | validated negative for the apply path; helps construction-heavy loops only |

**Net effect on the headline.** The warm $n=12$ apply that this report decomposes as $832\,\mathrm{ms}$ is now $\sim 273\,\mathrm{ms}$ by default, because the $\sim 581\,\mathrm{ms}$ "cache machinery" term (Section 3.1, the single largest finding) is essentially gone. The `$QuantumFrameworkPropCache = False` lever this report recommended as an interim user-side workaround is now nearly moot: the default-on path is within $\sim 1.2\times$ of the cache-off floor. Everything in Sections 1 through 6 below is the *pre-fix* analysis; read it as "why the apply was slow and how the dominant cost was removed," and read Section 5 Fix A as the remaining opportunity. The mechanism detail is preserved because it is the rationale for the shipped change (direct `"Basis"` getters + screened cache exclusions) and because the dense-in-sparse contraction cost (Section 3.4) and per-apply re-assembly (Section 3.3) are untouched and now the leading terms.

---

## 0. Two corrections that invalidate parts of the first pass

### 0.1 The first-pass benchmark simulated the zero vector

The first-pass input state `QuantumState[Table[0, {12}], 2]` is **not** the 12-qubit register $|0\rangle^{\otimes 12}$. The constructor `QuantumState[state_ ? stateQ, basisArgs__]` (`QuantumState/QuantumState.m:31`) treats the 12-element list as a *state vector* of length 12, pads it into the nearest qubit space, and produces a **4-qubit state with norm 0** (16-dim, zero explicit entries; `apply00e_state.wls`):

| check | value |
|---|---|
| `psi["Dimensions"]` | $\{2,2,2,2\}$ |
| `psi["Norm"]` | $0$ |
| explicit SparseArray entries | $0$ |

The apply wrapper (`QuantumCircuitOperator/QuantumCircuitOperator.m:120-137`) then silently pads it with eight fresh $|0\rangle$ qubits, so `qc[psi]` "works", but every tensor contraction runs over SparseArrays containing **no explicit values**. Amplitude arithmetic was therefore almost free in the first pass, which is what produced the (wrong) conclusions "contraction is 1% of the time" and "cost is flat in qubit number". With the correct register `QuantumState["Register"[12]]` the output has 4096 explicit amplitudes, norm 1, and the contraction costs $\sim 130$ ms at $n = 12$ (16% of the warm total) and $\sim 2.2$ s at $n = 16$ (where it **dominates**).

The bottom-line diagnosis of the first pass ("the time is per-gate object overhead, not arithmetic") survives at $n = 12$, but for the wrong quantitative reasons, and it fails at $n = 16$.

### 0.2 `H["Normal"]` never measured a normalization

`"Normal"` is not a `QuantumOperator` property. `QuantumOperator["H", {1}]["Normal"]` takes this path: the property dispatcher rejects it with `QuantumOperator::undefprop` (`QuantumOperator/Properties.m:70-82`), and the expression then matches the *application fallback* `(qo_QuantumOperator)[opts___] := QuantumCircuitOperator[qo][opts]` (`QuantumOperator/QuantumOperator.m:605`), so it evaluates `QuantumCircuitOperator[{H}]["Normal"]` and returns a one-gate **circuit**, after $6.3$ ms of message emission, failed dispatch, and circuit construction (`apply03_micro.wls`). The first-pass figure "`H["Normal"]` costs 3.8 ms vs 0.39 ms for the matrix" timed this accident, not a per-gate normalization. The genuine circuit-level `qco["Normal"]` costs **0.45 ms warm** for the whole 105-gate circuit (it is memoized, and `"NormalQ"` short-circuits), not 50% of the apply.

---

## 1. Reproduction and headline numbers

`apply00e_state.wls`, `apply01_stages.wls`, `apply06_prototypes.wls`:

| measurement | time |
|---|---:|
| cold apply, fresh kernel, correct register | $2.8$ - $2.9$ s |
| warm apply (identical circuit and state re-applied) | $800$ - $832$ ms ($7.7$ ms/gate) |
| warm apply with `$QuantumFrameworkPropCache = False` (memos already populated) | $223$ - $226$ ms |
| TN contraction alone, same network, real state | $129$ - $162$ ms |
| packed-array fold of the same 105 gate matrices (prototype A) | $4.3$ - $4.9$ ms |
| Qiskit Aer / Cirq on the same circuit (prior measurement, context) | $2$ - $12$ ms |

The first-pass "760 ms" is the warm number. The cold first apply is $3.6\times$ that; the difference ($\sim 2.1$ s) is one-time gate construction and first-compute of memoized properties (section 4.5).

Correctness anchors: the cache-off apply returns a state with max amplitude deviation $0.$ from the default path; prototype A agrees with the default path to $5 \times 10^{-17}$.

---

## 2. Where the warm 800 ms goes: location decomposition

Stage replication of `quantumCircuitApply` $\to$ `TensorNetworkApply` $\to$ `TensorNetworkCompile`, fully warm (run 3 of `apply01_stages.wls`; stages sum to 849 ms against a measured 811 ms total, $\sim 5\%$ double-count from shared memo hits):

| stage | code | ms | share |
|---|---|---:|---:|
| s0 apply wrapper + dimension guard | `QuantumCircuitOperator.m:97,120-137` | $101.2$ | $12\%$ |
| s1 `Flatten` | `QuantumCircuitOperator.m:106` | $0.9$ | |
| s2 circuit `Sort` | `TensorNetwork.m:252` | $0.5$ | |
| s3 prepend state as a gate | `TensorNetwork.m:255` | $148.7$ | $18\%$ |
| s4 `qco["Normal"]` | `TensorNetwork.m:185` | $0.5$ | $0.06\%$ |
| s5 `TensorNetworkBasis` (+`Info`) | `TensorNetwork.m:189` | $0.9$ | |
| s6 computational-basis test | `TensorNetwork.m:191-199` | $49.7$ | $6\%$ |
| s7 trace / eigen / order scans | `TensorNetwork.m:200-203` | $42.2$ | $5\%$ |
| s8 `QuantumTensorNetwork` build | `TensorNetwork.m:213` (body `:92-167`) | $368.7$ | $45\%$ |
| s9 `TensorNetworkContract` | `TensorNetwork.m:215` | $132.5$ | $16\%$ |
| s10 result wrap + unwrap | `TensorNetwork.m:216-246,258-265` | $3.4$ | $0.4\%$ |

Per gate: network build $3.5$ ms, state-prepend $1.4$ ms, contraction $1.3$ ms, wrapper/guard $1.0$ ms, scans $0.9$ ms.

Notes on individual stages:

- **s0** is dominated by `qco["InputDimensions"]` (`QuantumCircuitOperator/Properties.m:231-235`), which runs a per-qudit `SelectFirst` over `NormalOperators` and calls `ResourceFunction["LookupPart"]` per input qudit. Its *first* evaluation on a fresh circuit costs $365$ ms (resource-function resolution plus first `NormalOrders` fold); warm it is $0.5$ ms, and the remaining warm s0 cost is the composition wrapper that rebuilds the circuit (`qco /* QuantumCircuitOperator[{...identities...}]`).
- **s3** pays three `QuantumCircuitOperator` constructions (the rule list, the `RightComposition` UpValue flatten, the relabel), each re-running `FromCircuitOperatorShorthand` over all 106 elements and re-deriving the label composition.
- **s8** re-walks all 106 operators every apply: per operator it checks `#["Order"] === #["FullOrder"]`, tests input and output `ComputationalQ`, and fetches `#["Tensor"]`; the per-op property dispatches (at the warm-hit cost of section 3.1) are the bulk of the 369 ms. Nothing in s8 is cached at the network level: `"Flatten"`, `"Elements"`, `"Operators"` are on `$QuantumCircuitPreventCache` (`QuantumCircuitOperator/Properties.m:30-33`), and `QuantumTensorNetwork` itself is a plain function.
- **s9** is genuine contraction work, but over SparseArray tensors (section 3.4).

---

## 3. The WL mechanisms: why each millisecond burns

### 3.1 Mechanism 1 (dominant): the property-cache machinery taxes every dispatch

This is the single largest finding. The warm apply splits as

$$ 804 \text{ ms} \approx \underbrace{581 \text{ ms}}_{\text{cache machinery}} + \underbrace{94 \text{ ms}}_{\text{dispatch floor}} + \underbrace{129 \text{ ms}}_{\text{contraction}} , $$

established by the ablation: with `$QuantumFrameworkPropCache = False` but the memoized values still present as DownValues (they keep matching; the flag only disables the *machinery*), the identical apply returns a byte-identical state in $223$ ms (`apply09_cacheoff_validate.wls`).

The mechanism, from `QuantumOperator/Properties.m:70-82` (the same shape repeats for `QuantumState`, `QuantumBasis`, and the other heads):

```wolfram
(qo_QuantumOperator[prop_ ? propQ, args___]) /; QuantumOperatorQ[qo] := With[{result = QuantumOperatorProp[qo, prop, args]},
    If[ TrueQ[$QuantumFrameworkPropCache] &&
        ! MemberQ[{"Properties", "AllProperties", "State", "Basis"}, prop] &&
        QuantumOperatorProp[qo, "Basis"]["ParameterArity"] == 0,
        Quiet[QuantumOperatorProp[qo, prop, args] = result, Rule::rhs],
        result
    ] /; ...
]
```

Three compounding costs, all paid **on every access including cache hits**:

1. **The guard chain.** `QuantumOperatorProp[qo, "Basis"]["ParameterArity"]` is evaluated per access. `"Basis"` and `"State"` are deliberately excluded from caching, so the guard is never memoized; it walks operator $\to$ state $\to$ basis and dispatches `qb["ParameterArity"]` (`QuantumBasis/Properties.m:310`) through the basis head's own cached dispatcher. Measured: $341$ - $383\,\mu$s per evaluation (`apply03_micro.wls`).
2. **Write-on-hit.** The `Quiet[... = result, Rule::rhs]` re-executes `Set` even when `result` came from the memo. The `Set` itself is only $4.4\,\mu$s (`apply07_verify.wls`), but it is pure waste, and each assignment re-sorts the rule into a growing DownValues table.
3. **Net effect per access.** A warm *hit* on the trivial getter `op["Order"]` costs $367$ - $410\,\mu$s; with the cache flag off, the same access (recomputing from scratch) costs $12.6\,\mu$s. The cache makes trivial getters $\sim 30\times$ slower. A bare memo lookup (`QuantumOperatorProp[qo,"Order"]` directly) costs $1.6\,\mu$s, so $> 99\%$ of the hit cost is guard plus write-on-hit.

Call volume (`apply02_callcounts.wls`, counting rules prepended to the internal `Prop` tables; one fresh-circuit apply, cache on): $32{,}687$ `QuantumOperatorProp` plus $32{,}834$ `QuantumStateProp` plus $\sim 1{,}200$ basis/circuit computes, i.e. $\sim 620$ internal property evaluations **per gate**. The top entries are pure overhead, not physics:

| compute | count | role |
|---|---:|---|
| `QuantumStateProp[qs, "Basis"]` | $19{,}777$ | cache-guard chain (uncacheable by design) |
| `QuantumOperatorProp[qo, "State"]` | $19{,}373$ | getter hop in guard and delegation |
| `QuantumStateProp[qs, "AllProperties"]` | $12{,}890$ | delegation-rule condition (uncacheable) |
| `QuantumOperatorProp[qo, "AllProperties"]` | $6{,}445$ | delegation-rule condition |

With the cache flag off, the guard short-circuits at `TrueQ[...]` and total computes *drop* to $\sim 22$k. The cache does more work than it saves on this workload.

**Why the memo values still matter.** A fully cold apply with the cache disabled (no memos anywhere) costs $\sim 2.05$ s at 105 gates ($\sim 85$ ms/gate, `apply05_cacheoff_scaling.wls` and `apply08_scalingcounts.wls`), because the expensive per-operator properties really do get recomputed: with the cache off, `op["Tensor"]` recomputes at $2.8$ ms (H) to $16.8$ ms (CNOT), `op["MatrixRepresentation"]` at $5.2$ - $23.7$ ms, since every `"Tensor"`/`"Matrix"` access re-runs `qo["Sort"]` (which **always rebuilds** the operator and re-processes its label through `collectLabel`/`sortLabel`, even when already sorted; `QuantumOperator/Properties.m:275-287`, measured $2.3$ - $14.1$ ms). So the right model is: *memoized values good, per-hit machinery bad.*

### 3.2 Mechanism 2: delegation recomputes `AllProperties` three times per delegated access

Any property an operator answers via its state (`"Label"`, `"Picture"`, `"ComputationalQ"`, `"ParameterSpec"`, `"AmplitudesList"`, ...) goes through the last DownValue (`QuantumOperator/Properties.m:1072`):

```wolfram
QuantumOperatorProp[qo_, args : ...] /; MemberQ[Intersection[qo["State"]["AllProperties"], qo["AllProperties"]], prop] := qo["State"][args]
```

`"AllProperties"` is on every head's cache-prevention list (it must reflect the wrapped object), so the condition recomputes two $\sim 200$-element string unions plus an `Intersection` and `MemberQ` on **every** delegated access: measured $176\,\mu$s (operator) $+ 87\,\mu$s (state) per evaluation, $427\,\mu$s for one delegated `qo["AmplitudesList"]` warm. The counts above show $\sim 6{,}445$ delegated lookups per apply, so this condition alone accounts for $\mathcal{O}(100)$ ms. The same pattern sits on `QuantumCircuitOperator` (`Properties.m:560-562`).

### 3.3 Mechanism 3: nothing network-level is reused across applies

`TensorNetworkApply` rebuilds from scratch each call: the composition wrapper (s0), the state-prepend circuit (s3), and the whole `TensorNetwork` object (s8) including re-fetching all 106 tensors. The properties that would make reuse possible (`"Flatten"`, `"Operators"`, `"Elements"`) are on `$QuantumCircuitPreventCache`, and the network construction is a function, not a property. Hence the warm per-apply floor of $\sim 94$ ms of pure re-assembly even with all memos hot, paid identically on the thousandth application of the same circuit to the same state.

### 3.4 Mechanism 4: dense data lives in sparse containers

There are no packed arrays anywhere in the apply path. `quantumStateQ` *requires* the state to be a `SparseArray` (`QuantumState/QuantumState.m:14`); gate tensors come out of `SparseArray`-producing constructors. After one layer of Hadamards the state is **fully dense**: the output at $n = 12$ is a `SparseArray` with 4096 explicit entries (density 1.0), which costs $1.5\times$ the memory of a packed vector and $\sim 3.4\times$ on a dot product ($7.1$ vs $2.1\,\mu$s at $2^{12}$; `apply05_cacheoff_scaling.wls`). `On["Packing"]` during an apply shows only trivial unpacking (dimension-$\{1\}$ lists and one $\{4096, 11\}$ index array); the cost is not *unpacking*, it is *never being packed*. Consequences:

- contraction at $n = 12$: $129$ - $162$ ms vs $4.5$ ms for the identical arithmetic on packed arrays ($\sim 29\times$);
- contraction at $n = 16$: $2.2$ s vs $128$ - $156$ ms packed ($\sim 15\times$), and at $n=16$ contraction is $\sim 70\%$ of the QF total, so the "overhead, not arithmetic" diagnosis flips with scale.

### 3.5 Mechanism 5: cold construction (the other 2.1 s)

Named-gate construction costs, cold (`apply03_micro.wls`):

| construction | cold | mechanism |
|---|---:|---|
| `QuantumOperator["H", {q}]` | $2.7$ - $10.5$ ms | constructor cascade (dozens of conditional `DownValues` with `?QuantumStateQ`/`?orderQ` tests), `QuantumBasis` assembly ($1.5$ ms), label plumbing |
| `QuantumOperator["CNOT", {q, q{+}1}]` | $\sim 10$ ms | same plus controlled-label processing |
| `QuantumOperator["RZ"[0.3], {1}]`, fresh numeric angle | $44.5$ ms | see below |
| same angle repeated (system caches warm) | $16.7$ ms | |
| `QuantumOperator["RZ"[θ]]`, symbolic | $76.5$ ms | |

**Question 3 answered: where the RZ milliseconds go.** The rule (`QuantumOperator/NamedOperators.m:172`) is

```wolfram
QuantumOperator[("ZRotation" | "RZ")[angle_, dim_ : 2], opts___] :=
    QuantumOperator[FullSimplify @ QuantumOperator[Exp[- I angle / 2 QuantumOperator["PauliZ"[dim]]], "Label" -> ...], opts]
```

Measured decomposition for the numeric angle (repeated-arg regime, sums to $\sim 16.5$ ms): `QuantumOperator["PauliZ"[2]]` eigenbasis-operator construction $2.0$ ms; the `Exp` UpValue (a `MatrixExp` on a $2 \times 2$ SparseArray plus full operator reconstruction) $8.5$ ms; `FullSimplify` on the already-numeric operator only $0.37$ ms; the label re-wrap $2.6$ ms; the outer order wrap $2.0$ ms. So the first-pass hypothesis "symbolic simplification on numeric data" is **wrong for numeric angles**: `FullSimplify` is cheap there (it is the $\sim 30$ ms extra in the symbolic case); the cost is the *generic operator-exponential route* (build a `"PauliZ"` eigenbasis operator object, exponentiate through object algebra, rebuild, relabel, reorder) for a gate whose closed form is $\mathrm{diag}(e^{-i\theta/2}, e^{i\theta/2})$. A memoized or closed-form construction is $\sim 10^5 \times$ cheaper on repeat (memo hit $0.1\,\mu$s).

A note on `Memoize`: contrary to a plausible static reading, `Memoize[FromOperatorShorthand]` (`Utilities.m:382`, applied in `Init.m:7`) **does** work; WL's DownValues auto-ordering places appended literal entries before the catch-all trampoline, hits cost $< 1\,\mu$s, and the table does not grow on repeated identical calls (`apply07b_memoize.wls`). Its only pathology is unbounded growth keyed by full argument expressions (a memoized `psi -> Range[12]` entry retains the entire state).

### 3.6 Mechanism 6: superlinear growth in circuit length, caused by the cache, not the algorithm

Warm per-gate cost grows with gate count $G$ (cache on, `apply04_packing_ablations.wls`): $8.0$ ms/gate at $G = 35$, $9.5$ at $280$, $14.9$ at $560$, $22.5$ at $840$; a fit gives roughly $t \approx (7.3\,\text{ms}) \cdot G + (18\,\mu\text{s}) \cdot G^2$. Two pieces of evidence localize the quadratic term in the cache machinery:

1. With caching fully disabled and everything cold, per-gate cost is flat ($91.6 / 80.6 / 86.1$ ms/gate at $G = 35 / 140 / 560$; `apply08_scalingcounts.wls`): the *algorithm* is linear.
2. Circuit-level property call counts are **constant** in $G$ (`Elements` 248, `Operators` 169, `Flatten` 71, `NormalOrders` 60 per apply, identical at $L = 1, 4, 16$), ruling out call-count growth.

What grows is the cache work per operation: memo keys for circuit-level properties are $\mathcal{O}(G)$-sized expressions that must be hashed per lookup and per write-on-hit, and DownValues insertion re-sorts ever-larger rule tables ($2{,}532$ rules / 69 MB after one cold apply; $2{,}948$ / 118 MB after a second distinct circuit; `apply01_stages.wls`). Caveat: the sweeps ran ascending in one kernel, so table accumulation across sweep points is partially confounded with $G$; isolating per-$G$ in fresh kernels is a cheap follow-up. At the benchmark size the quadratic term contributes $\mathcal{O}(100)$ ms; beyond $G \approx 500$ it dominates.

### 3.7 Mechanism 7 (question 7): the OOM trap is a shape-broadcast rule, fed by the malformed state

`QuantumCircuitOperator[{psi -> Range[12]} /* qc]["Normal"]` killed kernels in the first pass because `psi` was the malformed 4-qubit object of section 0.1. The wrap routes through `FromOperatorShorthand[psi -> Range[12]]` (`NamedOperators.m:99`) $\to$ `QuantumOperator[QuantumOperator[psi], Range[12]]`, and since the inner operator has 4 output qudits while 12 are requested, the **multiplicity broadcast rule** (`QuantumOperator/QuantumOperator.m:245-258`, then `:279-281`) fires: it computes `LCM[Quotient[12,4], Quotient[12,1]] = 12` and builds `QuantumTensorProduct` of **twelve copies** of the 16-dim ket, a $16^{12} \approx 2.8 \times 10^{14}$-amplitude object. Verified: the construction aborts a 6 GB `MemoryConstrained` at $n = 12$ and is flat $\sim 250$ ms at $n \le 10$ (`apply00c_bisect.wls`, `apply00d_bisect2.wls`). With a correctly-shaped 12-qubit state the rule never fires (output count matches the order length) and the framework's own apply peaks at 171 MB. Two fixes: the silent broadcast should be guarded by a dimension message, and the `QuantumState[vector, dim]` constructor should reject or warn on a non-power-of-$d$ vector length rather than zero-padding.

---

## 4. Ablations (question 6)

`apply04_packing_ablations.wls`, warm standard benchmark unless noted:

| ablation | time | vs default |
|---|---:|---|
| default warm | $832$ ms | |
| `$QuantumFrameworkPropCache = False` (memos populated) | $226$ ms | $3.7\times$ faster, identical output |
| `$QuditBasisCaching = False` | $809$ ms | no effect (basis cache is off the warm path) |
| pre-normalized circuit `qc["Normal"]` | $814$ / $842$ ms | no effect (`Normal` was never the cost) |
| circuit from raw numeric matrices | $828$ ms warm | no effect (construction is not the warm cost) |
| parametric circuit, `qcP[<|θ -> v|>][psi]` per point | $3{,}160$ - $3{,}316$ ms | $4\times$ **slower** (substitution rebuilds every gate and parametric objects bypass the cache) |
| width $n = 8 / 12 / 16$ (warm, per gate) | $7.3$ / $7.6$ / $22.6$ ms | **not** flat in $n$ |
| depth $G = 35 \to 840$ (warm, per gate) | $8.0 \to 22.5$ ms | superlinear (section 3.6) |

The parametric-substitution result matters for variational workflows: for numeric sweeps, rebuilding the circuit ($\sim 180$ ms cold build) plus applying beats the `"Parameters"` substitution idiom by $\sim 4\times$ per point on this circuit.

---

## 5. Prototype-validated fixes (question 8)

`apply06_prototypes.wls`. All prototypes are session-local monkey-patches or user-level code; nothing was committed to the kernel.

### Fix A: numeric fast path (packed statevector fold). PAYS, by far the most.

Extract each gate's `MatrixRepresentation` once (dedup by label), convert to packed complex arrays, and fold over the state with the standard axes-permutation kernel: reshape the state to rank $n$, `Transpose` the target axes to the front, one matrix-matrix multiply $\left(2^k \times 2^k\right) \cdot \left(2^k \times 2^{n-k}\right)$, undo.

| measurement | value |
|---|---:|
| correctness vs default path | max deviation $4.9 \times 10^{-17}$ ($n{=}12$), $2.0 \times 10^{-17}$ ($n{=}16$) |
| fold, $n = 12$, 105 gates | $4.3$ - $4.9$ ms (**178$\times$** vs 800 ms warm; $41\,\mu$s/gate) |
| fold, $n = 16$, 141 gates | $128$ - $156$ ms (**21-25$\times$** vs 3.3 s warm) |
| one-time matrix extraction (3 distinct labels) | $233$ ms |
| end-to-end repeat including per-gate label lookups through the property layer | $123$ ms |

The $n = 12$ fold is inside Aer's measured envelope ($2.1$ - $8.6$ ms) on the same machine. The end-to-end repeat shows the integration risk: even *looking up* 105 labels through the property layer costs $25\times$ the simulation itself, so the fast path must capture its gate list without per-gate property dispatch (or with the cache machinery fixed first).

### Fix B: memoized parametric-gate matrices. Does NOT pay for circuit building; pays only for distinct-parameter churn.

Building the 105-gate circuit costs $182$ ms cold / $50$ ms repeat with plain named gates, and $178$ / $49$ ms with fully memoized gate objects: **no win**, because `FromOperatorShorthand` is already memoized and the build cost is the `QuantumCircuitOperator` constructor itself (per-element shorthand parse, parameter and label plumbing). The memoization wins only where many *distinct* numeric angles are constructed ($44.5$ ms each cold vs $0.1\,\mu$s memo hit), i.e. inside variational loops that rebuild gates; it does nothing for the apply path.

### Fix C: normalized-circuit / compile-once reuse. PAYS for repeated application.

At $n = 8$: `qop = qc["QuantumOperator"]` costs $1.1$ s once; subsequent `qop[psi]` cost $57$ ms vs $1{,}256$ ms per default apply, a **22$\times$** amortized win (break-even after one reuse). The materialized-operator route is exponential in memory ($4^n$), so it is only for small $n$ or few qubits; the right general form of this fix is reusing the *tensor network* (or the fast-path gate list) across applies, which the current code cannot do because `TensorNetworkApply` re-derives everything per call (section 3.3). For perspective, the raw matrix-vector product after extraction is $0.08$ ms, i.e. the QF wrapper on a precompiled operator still costs $700\times$ the arithmetic.

### Fix D (new, found by this audit): stop paying the cache machinery per hit. PAYS 3.7$\times$, smallest diff. **[SHIPPED 2026-06-14, `dfc741dc`/`f9dc1cdc`; see the status banner at the top.]**

Not in issue #5's list, but the cheapest large win: keep the memoized values, eliminate the per-hit guard and write-on-hit (section 3.1). The ablation bounds the win at $832 \to 226$ ms with byte-identical output. *(As shipped, the implementation was a direct `"Basis"` getter on every head plus widening the cache-exclusion list to the $O(1)$ structural getters so they skip the eligibility guard, and a `cacheProperty` store that commits only pattern-free keys; the live result is $832 \to 273\,\mathrm{ms}$ warm, $32 \to 2\,\mathrm{ms}$ for `op["Basis"]`$\times105$.)* Implementation sketch: compute the parametric-arity flag once per object (it is structural; cache it with the validity bit or read `ParameterSpec` by direct association access instead of three dispatches), and skip the `Set` when the value came from an existing memo (e.g. probe the memo with a cheap head-test before computing, or write through an explicit `Association` cache with single-write semantics instead of DownValues).

---

## 6. Answers to the audit questions, compressed

1. **Decomposition of the 7.7 ms/gate (warm):** network re-build $3.5$ ms/gate, state-prepend wrap $1.4$, contraction $1.3$, apply wrapper/dim guard $1.0$, basis/order scans $0.9$, everything else $< 0.1$ (section 2). By *mechanism*: cache machinery $5.5$ ms/gate, sparse contraction $1.2$, machinery floor $0.9$, arithmetic $0.04$ (section 3.1). Per apply there are $\sim 66$k internal property evaluations, $\sim 620$ per gate.
2. **Why "`Normal`" looked expensive:** it never existed; the measured 4 ms was an `undefprop` failure path plus an accidental circuit construction (section 0.2). The real `qco["Normal"]` is 0.45 ms warm.
3. **Why RZ costs 17-45 ms cold:** generic object-exponential construction route; `MatrixExp` UpValue $8.5$ ms, `PauliZ` operator object $2.0$ ms, label and order rewraps $4.6$ ms; `FullSimplify` is only $0.37$ ms for numeric angles (it is the cost only for symbolic ones) (section 3.5).
4. **Packed arrays:** none exist on the path; states are `SparseArray` by constructor contract, dense-as-sparse costs $3.4\times$ on kernels and $15$ - $29\times$ on contraction vs packed (section 3.4).
5. **Is the cache helping?** Its *values* yes (cold-no-memo is $2.05$ s vs $832$ ms), its *machinery* no: $367$ - $410\,\mu$s per hit vs $12.6\,\mu$s uncached recompute for getters; $581$ ms/apply of pure tax; plus $\sim 25$ - $50$ MB of memo growth per distinct circuit with no eviction. `Options`/`FilterRules` plumbing is real but secondary (inside the $\sim 100$ ms constructor stages).
6. **Ablations:** table in section 4.
7. **OOM:** multiplicity-broadcast rule tensoring 12 copies of a misshapen state; the framework's own path peaks at 171 MB (section 3.7).
8. **Fixes:** A $178\times$ (validated), C $22\times$ for reuse (validated), B no effect on apply (validated negative), D (new) $3.7\times$ (validated by ablation).
9. **Missed by the first pass:** zero-state benchmark artifact; contraction is 16% at $n{=}12$ and dominant at $n{=}16$; cost is not flat in $n$; superlinear depth scaling from cache-key hashing; write-on-hit and guard-chain dispatch tax; `AllProperties` delegation tax; parametric-substitution idiom slower than rebuilding; `ResourceFunction["LookupPart"]` in the dimension guard (365 ms first call, a latent network dependency in the hot path).

---

## 7. Prioritized implementation plan

| # | change | effort | risk | expected win (measured basis) |
|---|---|---|---|---|
| 1 ✅ **DONE** | **Cache machinery repair** (fix D): per-object parametric flag, no write-on-hit, optionally Association-backed memo with eviction | small (the 8 dispatch SubValues) | low: output verified identical; keep the flag as escape hatch | warm apply $3.7\times$ ($832 \to 226$ ms); also fixes the $G^2$ depth scaling and the 25-50 MB/circuit leak. **Shipped `dfc741dc`/`f9dc1cdc`; live $832 \to 273$ ms** |
| 2 | **Numeric fast path in `TensorNetworkApply`** (fix A): strict predicate (all ops numeric vector-form unitaries, computational bases, qubit dims, no trace/eigen/bend, no measurements), packed fold, fallback to current path otherwise | medium | medium: predicate must be exact; mitigate with property-based tests vs default path | $178\times$ at $n{=}12$ ($800 \to 4.5$ ms fold), $21\times$ at $n{=}16$; lands QF in Aer's bracket for statevector circuits |
| 3 | **Reusable compiled form on the circuit** (fix C generalized): cache the extracted gate list / network keyed on the circuit, reuse across applies and parameter values | medium | low-medium | $22\times$ measured for repeat-apply at $n{=}8$; unlocks variational loops (currently $3.2$ s/point) |
| 4 | **Guards on the broadcast rule and the vector-pad constructor** (section 3.7) | tiny | none | turns a kernel-killing OOM into an error message |
| 5 | **Closed-form / memoized named rotations** (fix B, rescoped): construct $R_{X,Y,Z}$ from closed-form matrices, memoize on (name, angle, dim) | small | low | $44.5 \to \sim 0.1$ ms per distinct angle; helps construction-heavy loops only |
| 6 | **Delegation condition cheapening** (section 3.2): precompute each head's static property-name set once per session instead of `AllProperties` per access | small | low | $\mathcal{O}(100)$ ms/apply, included in fix 1's ceiling |

Items 1 and 4 are safe immediate wins; item 2 is the headline competitive fix; item 3 falls out of item 2's compile step almost for free. **Update 2026-06-14:** items 1 and 4 shipped (`dfc741dc`/`f9dc1cdc`/`7e12e8f0`), so the warm path is now $\sim 273\,\mathrm{ms}$ by default and the misshapen-state OOM is guarded; item 2 (the packed numeric fold) is the remaining large win. Interim user-side guidance for the *pre-fix* kernel was to set `$QuantumFrameworkPropCache = False` after a warm-up; on the current kernel that lever is nearly moot (the default already pays no per-hit tax). Still useful: prefer rebuilding circuits over the `"Parameters"` substitution idiom in numeric sweeps, and never pass `QuantumState[Table[0, {n}], 2]` as a register (use `QuantumState["Register"[n]]`; it now warns).

---

## Appendix A: draft follow-up comment for issue #5 (not posted)

> Follow-up with a deeper decomposition; three corrections to the numbers in this issue and one new finding.
>
> 1. The benchmark state in the original measurements, `QuantumState[Table[0,{n}],2]`, is parsed as a length-$n$ *vector* (then zero-padded), not as an $n$-qubit register, so the circuits were applied to a norm-0 state and the contraction step was contracting empty SparseArrays. With a real register (`QuantumState["Register"[n]]`) the contraction is ~16% of the warm time at n=12 and ~70% at n=16, so "amplitude arithmetic is 1%" should be read as "1% at n=12 on a zero state". The per-gate-overhead diagnosis itself stands at n=12.
> 2. `qco["Normal"]` is not 50% of the apply: warm it is ~0.5 ms (memoized, and `NormalQ` short-circuits). The per-operator figure `H["Normal"]` ~4 ms in the issue text is an artifact: `"Normal"` is not an operator property; the call emits `undefprop` and then falls through to the application SubValue, returning a one-gate circuit.
> 3. The dominant warm cost is the property-dispatch cache machinery itself: every access, including cache hits, evaluates `qo["Basis"]["ParameterArity"]` (uncacheable by design, ~0.35 ms) and re-executes the memo `Set`. Measured: a warm hit on `op["Order"]` costs ~400 µs vs 12.6 µs with `$QuantumFrameworkPropCache = False`; disabling the cache flag (with memos already populated) makes the standard 105-gate apply 3.7x faster (832 -> 226 ms) with byte-identical output. This "fix the cache guard / skip write-on-hit" item is cheaper than any of the three fixes proposed above and should go first.
> 4. Of the three proposed fixes, prototypes validate: numeric fast path 178x at n=12 (4.5 ms for the full fold, agreeing to 5e-17); compiled-circuit reuse 22x for repeated application; parametric-gate memoization does **not** speed up circuit construction (FromOperatorShorthand is already memoized) and only pays for many distinct numeric angles.
> 5. New: warm cost is superlinear in gate count (~8 ms/gate at 35 gates -> ~22 ms/gate at 840) and the quadratic term disappears when caching is disabled; circuit-level property call counts are constant in depth, pointing at O(G)-sized memo keys and growing DownValues tables. Also, `QuantumCircuitOperator[{state -> order} /* qc]` with a state whose qudit count does not match the order length hits the multiplicity-broadcast rule (`QuantumOperator.m:245`) and builds an LCM-fold tensor power of the state (a 16^12-amplitude object in the reported OOM); a dimension guard would turn that kernel kill into a message.

## Appendix B: scripts

All under `scripts/`, runnable with `wolframscript -file <name>` against the working tree:

- `apply00_diag.wls`, `apply00b_diag.wls`, `apply00c_bisect.wls`, `apply00d_bisect2.wls`: OOM bisection (section 3.7).
- `apply00e_state.wls`: benchmark-state misparse and zero-state artifact (section 0.1).
- `apply01_stages.wls`: stage decomposition, cold/perturbed/warm, memo-table growth (section 2).
- `apply02_callcounts.wls`: per-property compute counts via prepended counting rules, cache on/off (section 3.1).
- `apply03_micro.wls`: dispatch, guard, Sort/Tensor/Matrix, constructor, and RZ micro-costs (sections 3.1, 3.5).
- `apply04_packing_ablations.wls`: packing audit, cache/basis-cache/pre-normal/raw-matrix/parametric ablations, depth and width sweeps (sections 3.4, 4).
- `apply05_cacheoff_scaling.wls`: cache-off sweeps, contraction-only scaling, sparse-vs-packed micro (sections 3.4, 3.6).
- `apply06_prototypes.wls`: fix prototypes A, B, C with correctness checks (section 5).
- `apply07_verify.wls`, `apply07b_memoize.wls`: write-on-hit cost, `InputDimensions` first-call cost, `Memoize` behavior (sections 3.1, 3.5).
- `apply08_scalingcounts.wls`: circuit-property call counts vs depth, cold cache-off regime (section 3.6).
- `apply09_cacheoff_validate.wls`: cache-off output equality and the 223 ms reconciliation (section 3.1).

Timing discipline: stage and sweep numbers are single-shot `AbsoluteTiming` in stated cache regimes (cold / perturbed / warm are distinct rows, not minima across regimes); micro-costs are min-of-batches with `ClearSystemCache[]` between batches and `HoldFirst` timing helpers. The known first-pass pitfall (timing a pre-evaluated expression) is avoided throughout; correctness anchors accompany every prototype.
