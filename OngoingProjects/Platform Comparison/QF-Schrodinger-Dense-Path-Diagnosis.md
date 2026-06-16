# QuantumFramework Dense State-Vector Diagnosis: Is `Method -> "Schrodinger"` Good?

**Scope.** A focused performance diagnosis of QuantumFramework's dense (full $2^n$ amplitude
vector) circuit simulation, with special attention to the non-default `Method -> "Schrodinger"`
path. Tests the hypothesis "the non-TN / Schrodinger method of QF is good." Companion to
`QF-Apply-Path-Deep-Audit.md` (which decomposes the *default* TensorNetwork apply) and Section
3.1 / 12.2 of `QF-Master-Platform-Comparison.md`.

**Environment.** Apple M2 Pro, Mathematica 15.0, paclet 2.0.0 from the working tree
(`PacletDirectoryLoad`), HEAD `f9dc1cdc`. Timing discipline: `HoldFirst` helper,
`ClearSystemCache[]` before every repetition, min over repetitions. Scripts:
`scripts/schrodinger_diag.wls`, `scripts/schrodinger_symbolic.wls`, plus the existing
`scripts/qf_bench4_register.wl` and `scripts/apply06_prototypes.wls`. Standard circuit: $H$ on
all qubits, a CNOT chain, $R_Z(0.3)$ on all qubits, repeated for $L$ layers, applied to the true
register `QuantumState["Register"[n]]`.

---

## 0. Verdict up front

**The non-TN `Method -> "Schrodinger"` path is not "good" as a performance choice. It is the
slowest dense option QF offers, at every size measured, with no speed crossover.** The qf-skill
note that Schrodinger is "the right choice for small circuits and fidelity benchmarking, where TN
overhead dominates" is **not supported by the data as a speed claim**: even a 1-gate circuit on 1
qubit runs faster through the default TensorNetwork path. Schrodinger has exactly one legitimate
niche, and it is about **correctness/memory, not speed**: circuits whose TensorNetwork contraction
order produces an intermediate larger than the $2^n$ statevector (QPE / HHL at precision $\geq 5$),
where the default path exhausts memory and gate-by-gate folding is the only thing that finishes.

The single dominant bottleneck of the Schrodinger path is structural and surprising: it is **not**
the `FullSimplify` wrapper (that is $0.3\%$ to $5\%$ of the time), and **not** the $2^n$ amplitude
arithmetic. It is that Schrodinger issues **one full framework circuit-application per gate**, each
carrying a fixed $\sim 40\,\mathrm{ms}$ re-assembly overhead that is essentially flat in $n$. With
$G$ gates that is $G \times 40\,\mathrm{ms}$; the default path pays that fixed overhead roughly
once for the whole circuit.

---

## 1. Reproduction and the crossover hunt (Task 1)

### 1.1 Headline reproduction at `f9dc1cdc`

`scripts/qf_bench4_register.wl`, 105-gate circuit ($L=3$), statevector readout, min-of-3. The
clean documented numbers and this re-measurement (taken under some concurrent load, so absolute
values run high, but the ratio is stable):

| $n$ | QF default (TN) | QF `Schrodinger` | ratio Sch/TN | Aer (documented) |
|---:|---:|---:|---:|---:|
| 12 | $283\,\mathrm{ms}$ (repro $323$) | $6564\,\mathrm{ms}$ (repro $7286$) | $\sim 23\times$ | $2.1\,\mathrm{ms}$ |
| 16 | $2716\,\mathrm{ms}$ (repro $3799$) | $62{,}545\,\mathrm{ms}$ (repro $69{,}038$) | $\sim 23\times$ | $8.6\,\mathrm{ms}$ |

Confirmed: Schrodinger is $\sim 23\times$ slower than the TN default, which is itself $\sim 130\times$
($n=12$) to $\sim 320\times$ ($n=16$) slower than Aer. Schrodinger is $\sim 3000\times$ to
$\sim 7000\times$ Aer.

### 1.2 There is no speed crossover (numeric, shallow)

`scripts/schrodinger_diag.wls`, 1 layer, varying $n$, min-of-3/4:

| $n$ | gates | TN default | Schrodinger | ratio Sch/TN |
|---:|---:|---:|---:|---:|
| 1 | 2 | $32\,\mathrm{ms}$ | $77\,\mathrm{ms}$ | $2.4\times$ |
| 2 | 5 | $38\,\mathrm{ms}$ | $182\,\mathrm{ms}$ | $4.8\times$ |
| 3 | 8 | $39\,\mathrm{ms}$ | $297\,\mathrm{ms}$ | $7.7\times$ |
| 4 | 11 | $43\,\mathrm{ms}$ | $431\,\mathrm{ms}$ | $10.0\times$ |
| 5 | 14 | $46\,\mathrm{ms}$ | $559\,\mathrm{ms}$ | $12.1\times$ |
| 6 | 17 | $55\,\mathrm{ms}$ | $703\,\mathrm{ms}$ | $12.8\times$ |
| 8 | 23 | $61\,\mathrm{ms}$ | $1030\,\mathrm{ms}$ | $16.8\times$ |

Even the smallest cases lose. A bare $H$ on 1 qubit: TN $29\,\mathrm{ms}$, Schrodinger
$41\,\mathrm{ms}$. A 2-gate Bell circuit: TN $30\,\mathrm{ms}$, Schrodinger $78\,\mathrm{ms}$. The
ratio only *grows* with size, because Schrodinger's cost is $\propto G$ while the TN path amortizes
its fixed overhead. **The claimed "small/shallow circuits favor Schrodinger" crossover does not
exist for speed.** Depth sweep at $n=4$ confirms the same shape: $L=1$ TN $45$ / Sch $420$; $L=16$
($176$ gates) TN $433$ / Sch $6598$, a flat $\sim 10$ to $15\times$ gap.

### 1.3 No symbolic crossover either

`scripts/schrodinger_symbolic.wls`, symbolic $R_Z(\theta)$, 1 layer. The plausible defense of
Schrodinger is that its per-gate `FullSimplify` yields more compact symbolic amplitudes. It does
not: the default TN path already returns clean closed forms, and the two outputs have **identical**
`LeafCount` ($n=1$: $31$ vs $31$; $n=2$: $29$ vs $29$; $n=3$: $145$ vs $145$), agree after
`Simplify`, and Schrodinger is still $2$ to $7\times$ slower. For example both paths return
$\{e^{-i\theta/2}/\sqrt{2},\, e^{i\theta/2}/\sqrt{2}\}$ for `{H, RZ[θ]}`. The per-gate
`FullSimplify` is redundant work, not a feature, on these circuits.

### 1.4 The one genuine niche: TN contraction blow-up (memory, not speed)

The default path converts the circuit to a tensor network and contracts it
(`qf-tn-integration.md`). For local circuits like the benchmark the contraction stays at
$2^n$ and is fast. For QPE / HHL structures with many controlled $U^{2^k}$ powers, the
contraction order bloats the intermediate qudit count: documented at
`project_qf_2_0_0_hhl_tn_crash` and the qf-skill, precision 3 gives $2^{11}$, precision 4 gives
$2^{19}$, precision 6 gives $2^{35}$ and exhausts memory. Schrodinger never builds the network: it
folds gate by gate, holding the state at $2^n$ throughout, so it finishes where the default path
dies. **This is the real and only reason to reach for `Method -> "Schrodinger"`: as a
memory-safety fallback for circuits whose TN contraction explodes, not as a fast path.**

---

## 2. Where the Schrodinger time actually goes (Task 2)

The kernel definition (`QuantumCircuitOperator/QuantumCircuitOperator.m:100-105`, `f9dc1cdc`):

```wolfram
"Schrodinger" :> Fold[
   (FullSimplify @ #2[#1, Method -> {"TensorNetwork", "Computational" -> False}]) &,
   qs, qco["NormalOperators"]]
```

So each gate `#2` is applied to the running state `#1` via a **single-gate TensorNetwork apply**,
and the whole resulting state object is wrapped in `FullSimplify`. Profiling
(`scripts/schrodinger_diag.wls`, Part 3, min-of-2/3, $L=3$):

| $n$ | gates | full Sch | fold without `FullSimplify` | `FullSimplify` tax | one single-gate apply | bulk TN |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 51 | $2213\,\mathrm{ms}$ | $2129\,\mathrm{ms}$ | $85\,\mathrm{ms}$ ($3.8\%$) | $40\,\mathrm{ms}$ | $104\,\mathrm{ms}$ |
| 8 | 69 | $3077\,\mathrm{ms}$ | $3067\,\mathrm{ms}$ | $10\,\mathrm{ms}$ ($0.3\%$) | $44\,\mathrm{ms}$ | $146\,\mathrm{ms}$ |
| 10 | 87 | $4441\,\mathrm{ms}$ | $4228\,\mathrm{ms}$ | $213\,\mathrm{ms}$ ($4.8\%$) | $46\,\mathrm{ms}$ | $195\,\mathrm{ms}$ |

Three readings, the first counterintuitive:

1. **`FullSimplify` is not the bottleneck.** It is $0.3\%$ to $5\%$ of the runtime. `FullSimplify`
   on the *final* numeric state object costs $0.05$ to $0.06\,\mathrm{ms}$ (it threads to the
   machine-number amplitudes and finds nothing to do). The intuition that "FullSimplify on a $2^n$
   numeric vector is the killer" is **wrong** here. It is wasteful and should be gated behind a
   `SymbolicQ` check, but removing it would only recover a few percent.

2. **The dominant cost is the per-gate apply itself.** The fold without `FullSimplify` is
   essentially the whole runtime, and it equals (gate count) $\times$ (one single-gate apply):
   $87 \times 46\,\mathrm{ms} = 4002\,\mathrm{ms}$ vs measured $4228\,\mathrm{ms}$ at $n=10$.
   Schrodinger re-enters the entire framework apply pipeline once per gate.

3. **That per-gate apply is $\sim 40\,\mathrm{ms}$ and nearly flat in $n$** ($40 / 44 / 46\,\mathrm{ms}$
   at $n = 6 / 8 / 10$). If the $2^n$ amplitude arithmetic dominated it would quadruple every $+2$
   qubits; it does not. The cost is the fixed framework re-assembly per apply: the apply wrapper +
   dimension guard, prepending the state as a gate, building a one-gate `TensorNetwork` object,
   contraction setup, and reconstructing a `QuantumState` from the contracted tensor (the same
   per-apply stages `QF-Apply-Path-Deep-Audit.md` Section 2 decomposes for the bulk path: s0, s3,
   s8, s9, s10). Schrodinger pays that fixed cost $G$ times; the bulk TN path pays it roughly once
   and amortizes the network build across all gates ($195\,\mathrm{ms}/87 \approx 2.2\,\mathrm{ms}$
   effective per gate vs $46\,\mathrm{ms}$ for Schrodinger, the $\sim 20\times$ gap).

**Ruled in / ruled out** against the hypotheses in the task brief:
- Gates applied as objects, re-deriving structure per apply: **confirmed dominant** (the per-gate
  framework re-assembly).
- Dense amplitudes in `SparseArray` not packed `Developer`PackedArray`: **confirmed present** (Part
  4: the $n=10$ output is a `SparseArray` with all $1024$ entries explicit, `PackedArrayQ` is
  `False`); this is the $15$ to $29\times$ contraction penalty from the deep audit and it taxes
  *both* dense paths, but it is secondary to the per-gate re-assembly for Schrodinger specifically.
- Per-apply tensor-network re-assembly: **confirmed, and multiplied by $G$**.
- Pattern-matched dispatch in inner loops: present inside each apply (the property dispatch), now
  much cheaper after the June cache refactor, but still paid $G$ times.
- `FullSimplify` on numeric data: **present but negligible** ($\leq 5\%$).
- `Dot` / `KroneckerProduct` materializing a full $2^n \times 2^n$ operator: **not** what
  Schrodinger does (it applies each gate on its own axes); that trap is `qc["Matrix"]` /
  `qc["QuantumOperator"]`, a different path.

---

## 3. The structural gap vs Aer (Task 3)

| | per gate | data layout | arithmetic |
|---|---|---|---|
| **Aer** | in-place update of the statevector on the gate's axes; gate fusion | one contiguous packed complex array, SIMD/BLAS, multi-threaded | $O(2^n)$ per gate, no object churn |
| **QF default TN** | one network build + one contraction sweep for the whole circuit | dense data in `SparseArray` containers ($15$ to $29\times$ penalty) | $O(2^n)$, but wrapped in framework objects |
| **QF Schrodinger** | one **full framework apply** per gate (network build + contract + reconstruct), $\sim 40\,\mathrm{ms}$ fixed each | same `SparseArray` containers | $O(2^n)$ swamped by per-apply overhead |

The structural difference is stark: Aer mutates one packed array $G$ times with no per-gate
allocation; QF Schrodinger constructs and tears down a tensor-network apply (objects, bases,
labels, a fresh `QuantumState`) $G$ times. Even with the heavy arithmetic on BLAS, the orchestration
overhead per gate ($\sim 40\,\mathrm{ms}$) is four orders of magnitude above Aer's per-gate cost.

---

## 4. The packed numeric fold fast path (Task 4)

`scripts/apply06_prototypes.wls`, re-measured at `f9dc1cdc`. The prototype extracts each gate's
numeric matrix once (dedup by label), converts to packed complex arrays, and folds over the
statevector with the textbook axes-permutation kernel (reshape to rank $n$, move the target axes to
the front, one $2^k \times 2^k$ by $2^k \times 2^{n-k}$ matrix multiply, undo):

| measurement | $n=12$ | $n=16$ |
|---|---:|---:|
| correctness vs default path (max amplitude deviation) | $4.9\times 10^{-17}$ | $2.0\times 10^{-17}$ |
| packed fold (min-of-3) | $6.7\,\mathrm{ms}$ | $201\,\mathrm{ms}$ |
| default TN warm (this run) | $337\,\mathrm{ms}$ | $4226\,\mathrm{ms}$ |
| speedup vs default | $\sim 50\times$ | $\sim 21\times$ |
| output `PackedArrayQ` | `True` | `True` |
| one-time matrix extraction | $77\,\mathrm{ms}$ | |
| end-to-end incl. per-gate label lookups through property layer | $22\,\mathrm{ms}$ | |

At $n=12$ the fold ($6.7\,\mathrm{ms}$) is in Aer's measured envelope ($2$ to $9\,\mathrm{ms}$); the
quoted "$178\times$" is against the pre-cache-fix $800\,\mathrm{ms}$ warm path, $\sim 42$ to
$50\times$ against today's faster default. At $n=16$ the fold is $\sim 24\times$ Aer: the residual
gap is the per-gate `Transpose`/`ArrayReshape` copies (each $O(2^n)$, no true in-place update) and
the absence of gate fusion. **The end-to-end number is the integration warning:** merely *looking
up* 105 gate labels through the property layer costs $\sim 22\,\mathrm{ms}$, more than $3\times$ the
fold, so the fast path must capture its gate list without per-gate property dispatch.

**What shipping it takes, and what it does not cover.** It must be a strict-predicate fast path
inside `TensorNetworkApply` with a fallback to the current path. The predicate has to hold: all
operators are numeric, vector-form, unitary gates on qubit ($d=2$) computational bases, with no
trace / eigen / bend, no density-matrix or channel operands, and no mid-circuit measurements. It
**does not cover**, and must fall back for: symbolic amplitudes (packed arrays do not apply),
qudits $d>2$ (the reshape is hard-coded to `ConstantArray[2,n]`; mixed-radix needs a different
kernel), density-matrix / channel inputs (need the doubled representation), and mid-circuit
measurement / classical control (need branching). Risk is predicate exactness (endianness, basis
and order conventions must match the default path bit-for-bit); it is mitigable, since the
prototype already agrees with the default path to $5\times 10^{-17}$ via property-based testing.

---

## 5. Verdict, regime map, and ceiling (Task 5)

**Is `Method -> "Schrodinger"` good?** No, not for performance. Precisely:

- **Where it loses (everywhere, for speed):** every $n$ and every depth measured. It is
  $\sim 23\times$ the default TN path on the 105-gate benchmark and $2.4$ to $17\times$ on shallow
  circuits, with the gap growing in both width and depth. There is no small-circuit speed
  crossover.
- **Where it wins (correctness/memory, narrow):** circuits whose default TN contraction order
  builds an intermediate larger than the $2^n$ statevector (QPE / HHL at precision $\geq 5$), where
  the default path runs out of memory. Schrodinger's gate-by-gate fold keeps the state at $2^n$ and
  finishes. Use it as a memory-safety fallback, not a fast path. (For symbolic work the default TN
  path is both faster and equally compact, so Schrodinger has no symbolic advantage either.)
- **Single dominant bottleneck of the dense path:** for Schrodinger, one full framework apply per
  gate ($\sim 40\,\mathrm{ms}$ fixed re-assembly, flat in $n$), times $G$. For the dense path in
  general, the framework's per-apply object overhead plus dense amplitudes living in `SparseArray`
  rather than packed arrays. Both are framework cost, not physics and not algorithm.
- **Realistic ceiling:** QF can get *close* to Aer at moderate $n$ but not to parity. The packed
  fold prototype lands in Aer's bracket at $n=12$ ($6.7$ vs $2$ to $9\,\mathrm{ms}$) and is
  $\sim 24\times$ at $n=16$. The remaining floor is single-threaded WL orchestration plus per-gate
  array copies (WL has no in-place statevector mutation and no gate fusion); the heavy arithmetic
  itself rides BLAS, so the floor is a small constant factor at moderate $n$ widening with $n$, not
  a fundamental $100\times$ wall. Object overhead and copy-per-gate are the hard part.
- **Escape hatches that already exist, in order of preference for fast dense work:**
  1. **Default TN** (`Automatic`): the right in-product choice for this workload, $\sim 23\times$
     faster than Schrodinger.
  2. **`Method -> "QuEST"`** (qubits only, `QuEST.m:58`): a C-backed in-place statevector engine,
     the closest thing to Aer inside QF. It requires the QuEST backend to be compiled/installed; it
     did not load in this session, so no local number, but architecturally it is the in-product
     near-Aer path for pure qubit circuits.
  3. **The qiskit bridge** (`Method -> {"Qiskit", ...}`): delegates to qiskit/Aer through Python,
     i.e. literally Aer speed minus the round-trip serialization, when a Python qiskit stack is
     present.
  4. **The unshipped packed fold:** the one large native win still on the table
     (`QF-Apply-Path-Deep-Audit.md` Section 7, item 2).

**One-line correction to propagate:** the qf-skill bullet "`Method -> "Schrodinger"` is also the
right choice for small circuits and fidelity benchmarking, where TN overhead dominates" should be
rewritten to "`Method -> "Schrodinger"` is a memory-safety fallback for circuits whose TensorNetwork
contraction blows up (QPE / HHL high precision); it is never the faster dense path."

---

## Appendix: scripts

- `scripts/schrodinger_diag.wls`: crossover sweep, depth sweep, the per-gate / `FullSimplify`
  profile decomposition, and the packed-vs-sparse container check (Sections 1.2, 2).
- `scripts/schrodinger_symbolic.wls`: symbolic-regime compactness and timing (Section 1.3).
- `scripts/qf_bench4_register.wl`: headline $n = 12 / 16$ TN-vs-Schrodinger reproduction
  (Section 1.1).
- `scripts/apply06_prototypes.wls`: the packed numeric fold prototype A (Section 4).

All run with `wolframscript -file <abs-path>` against the working tree and emit `RESULT|label|value`
lines. Timing helper is `HoldFirst` with `ClearSystemCache[]` per repetition and min-over-reps.
