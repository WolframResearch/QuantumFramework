# Computing probabilities and expectations in QuantumFramework without paying the measurement-dilation tax

**Status:** audit + benchmark. Written 2026-06-15.
**Companion to:** [`measurement-hilbert-space-blowup.md`](measurement-hilbert-space-blowup.md) (the design note on *why* QF dilates and the menu of escape routes). This document is the empirical, decision-oriented half: given that dilation is expensive, *what should you actually call* to get a probability distribution, a single outcome probability, or an expectation value, and how much does each route cost.

**Anchor.** Live working-tree paclet at HEAD `f9dc1cdc` (the property-cache-refactor anchor; see `audit/last-synced.md`), loaded via `PacletDirectoryLoad`. All timings on the same machine as the platform-comparison docs (M2 Pro), `ClearSystemCache[]` before every measurement, no `Quiet`. Repro scripts: `/private/tmp/qf_meas_probe{,2,3,4,5}.wls` (also reproduced inline in the appendix).

---

## 0. The one-paragraph answer

QF realizes a projective measurement as a *unitary* that appends a pointer qudit (von Neumann / Stinespring dilation), so the state vector grows from $2^n$ to $2^{n+k}$ after measuring $k$ qubits. That is physically honest and indispensable when you need the collapsed state (feedforward, teleportation, error correction), but it is the wrong tool for the routine job of *reading out numbers*. To get the **full Born distribution**, build the state once and call `state["ProbabilitiesList"]`: cost $O(2^n)$, no dilation. To get a **single expectation value** $\langle\psi|O|\psi\rangle$ or a **single overlap** $\langle\phi|\psi\rangle$, contract a bra-operator-ket inner product (your `QuantumCircuitOperator[{ket, op, bra}]` form, or equivalently `(\[Psi]^\[Dagger] @ O @ \[Psi])["Scalar"]`): it stays bounded where the measurement-operator route blows up. The two options you proposed are both correct and both avoid dilation; they answer different questions (whole distribution vs single scalar). The route to avoid for large systems is `QuantumMeasurementOperator[...]["Mean"]` / applying a measurement operator just to harvest a number.

---

## 1. The question in physics terms

You have a state $|\psi\rangle = U|0\rangle^{\otimes n}$ produced by some circuit, and you want one of three things:

1. the Born distribution $p_x = |\langle x|\psi\rangle|^2$ over a basis $\{|x\rangle\}$ (a full readout, $2^n$ numbers);
2. the probability of one specific outcome, $p_{x_0} = |\langle x_0|\psi\rangle|^2$ (one number);
3. the expectation of an observable, $\langle O\rangle = \langle\psi|O|\psi\rangle$ (one number).

There are three families of QF calls that compute these, and they differ enormously in cost:

- **Route M (measurement operator).** `QuantumMeasurementOperator[...]` applied to the state, or appended to the circuit. This is the only route that *dilates*: it builds the pointer register.
- **Route S (state then probabilities).** Compute $|\psi\rangle$ once with `qc[]`, then read `state["Probabilities"]` / `["ProbabilitiesList"]`.
- **Route I (inner product).** Contract a bra-operator-ket sandwich to a scalar: $\langle\psi|O|\psi\rangle$ or $\langle x_0|\psi\rangle$, with no pointer qudit. This is the `QuantumCircuitOperator[{QuantumState[...]["Dagger"], ..., QuantumState[...]}]` form you asked about.

Routes S and I are what you want almost always. Route M is for when you genuinely need the *measurement object* (the collapsed state, a classical record, shot counts), not a number.

---

## 2. How the dilation cost arises (kernel anchors)

`QuantumMeasurementOperator/Properties.m:206-276` ("SuperOperator") implements the dilation: it diagonalizes the target observable (`:222`), builds a pointer `eigenBasis` (`:226-232`), and tensors it onto the output (`:241-259`). Each measured register adds a qudit. Companion accessors: `"ExtraQudits"` (`:81`), `"Eigenqudits"`/`"Eigendimensions"` (`:83-86`), the cheaper alternatives `"DiscardExtraQudits"` (`:300-315`), `"POVM"` (`:278`), and the channel form `"Bend"` (`:286-294`).

Two separate cost mechanisms make Route M expensive:

- **Dimension squaring.** Measuring all $n$ qubits appends an $n$-qubit pointer, so the dilated state lives on $2^{n+n}=2^{2n}$. Verified exactly (Table A).
- **Generic eigendecomposition.** To measure an arbitrary observable $O$, the operator must be diagonalized. QF's `eigensystem` (`Utilities.m:151`) is generic: it wraps `Eigensystem` with `Simplify` and `Chop`, with **no Hermitian fast path**. For an $n$-qubit observable this is an $O((2^n)^3)=O(8^n)$ dense solve plus `Simplify`. This, not the dilation, is what makes Route M for expectation values catastrophic (Table D).

By contrast, `state["ProbabilitiesList"]` is just $|\text{amplitude}|^2$, an $O(2^n)$ elementwise operation with no dimension cap (unlike `VonNeumannEntropy`/`Purity`, which silently cap at $2^{10}$). The inner product $\langle\psi|O|\psi\rangle$ is a sparse matrix-vector product plus a contraction.

---

## 3. Benchmarks

Circuit: $n$ qubits, a Hadamard wall, a CNOT ladder, and a layer of `RY` rotations (about $3n$ gates), applied to $|0\rangle^{\otimes n}$ by the default tensor-network path. All times in ms.

### Table A. Dilation squares the dimension (measure all $n$ qubits)

| $n$ | system $2^n$ | dilated `["SuperOperator"]` output | $2^{2n}$ | applied dilated state `["Dimension"]` |
|----:|----:|----:|----:|----:|
| 4 | 16 | 256 | 256 | 256 |
| 6 | 64 | 4096 | 4096 | 4096 |
| 8 | 256 | 65536 | 65536 | 65536 |

Exact $2^{n+k}$ growth. At $n=13$ the dilated vector for a full readout is $2^{26}\approx 6.7\times10^7$ complex amplitudes ($\sim 1$ GB) before any contraction; the squaring is the memory wall the design note flags.

### Table B. Shared cost: building the state once (`qc[]`)

| $n$ | 6 | 8 | 10 | 12 |
|---|---:|---:|---:|---:|
| state prep (ms) | 289 | 320 | 344 | 383 |

Roughly flat over this range (the default TN apply path, post cache-refactor). This is the unavoidable cost that every route pays; the question is what each route adds *on top*.

### Table C. Full distribution: Route S vs Route M

| $n$ | Route S: `state["ProbabilitiesList"]` (given $\psi$) | Route M, one joint $n$-qubit measurement (end-to-end) | Route M, $n$ separate 1-qubit measurements (end-to-end) |
|----:|----:|----:|----:|
| 4 | 0.6 | 207 | 384 |
| 6 | 0.7 | 291 | 344 |
| 8 | 0.7 | 562 | 380 |
| 10 | 0.8 | 2059 | 443 |
| 12 | 1.2 | (grows) | (grows) |

All routes return the identical distribution (max abs difference $\le 4\times10^{-17}$). Reading: extracting the distribution from an already-built state is essentially free ($\sim1$ ms); the measurement-operator route's end-to-end cost is dominated by the dilated contraction. A *joint* $n$-qubit measurement builds one dim-$2^n$ pointer and grows fast (2.06 s at $n=10$); $n$ *separate* single-qubit measurements distribute the pointer and the tensor network handles them better (it stays near the bare-prep cost, $\sim 0.4$ s, through $n=10$). Either way Route S wins, and the gap widens exponentially because Route M is contracting toward a $2^{2n}$ object.

### Table D. Expectation value $\langle Z^{\otimes n}\rangle$: Route I is the only one that stays bounded

| $n$ | `Tr[O . \[Rho]]` (dense $4^n$) | Route I: $(\psi^\dagger\,O\,\psi)["Scalar"]$ | Route M: `QuantumMeasurementOperator[O][\[Psi]]["Mean"]` |
|----:|----:|----:|----:|
| 4 | 7.4 | 127 | 584 |
| 6 | 8.2 | 130 | **21162** |
| 8 | (OOM-prone) | 135 | **kernel killed (OOM)** |
| 10 | - | 135 | - |
| 12 | - | 148 | - |

All agree to $3\times10^{-17}$ where they complete. The measurement-operator route is **the worst possible choice for an expectation value**: it diagonalizes the observable with the generic `Simplify`-laden eigensolver *and* dilates, exploding from 0.58 s ($n=4$) to 21 s ($n=6$) to an out-of-memory kernel death ($n=8$). The bra-operator-ket inner product is flat at $\sim 135$ ms out to $n=12$. `Tr[O\cdot\rho]` is fastest of all for small $n$ but densifies the $2^n\times2^n$ density matrix, so it hits the memory wall around $n=12$-$13$; the inner product avoids that by staying in vector form.

### Table E. Single-outcome probability $p_{0\ldots0}$ (state already in hand)

| $n$ | 4 | 6 | 8 | 10 | 12 |
|---|---:|---:|---:|---:|---:|
| `First[state["ProbabilitiesList"]]` | 0.6 | 0.7 | 0.7 | 0.8 | 1.2 |
| projector overlap `Abs[(QuantumState["0.."]^\[Dagger] \[Psi])["Scalar"]]^2` | 887* | 36 | 40 | 42 | 49 |

Agreement to $10^{-18}$. (*$n=4$ is first-call warmup.) **Subtlety worth stating plainly:** once you already hold $|\psi\rangle$, indexing the probability list beats building a projector overlap, because the overlap pays a tensor-network contraction ($\sim 40$ ms) while the index is a lookup. The inner-product form (Route I) for a single outcome is the right tool only when you want to *avoid materializing or holding the full state*, or want the result symbolic/parametric. For expectation values of a non-trivial observable (Table D) the inner product is the clear winner; for a single computational-basis outcome from an existing state, just index `ProbabilitiesList`.

### Table F. Non-computational basis: rotate with gates, never rebase the state

To get probabilities in, say, the $X$ eigenbasis:

| $n$ | Route via gates: append $H$ wall, read computational probs | Route via dense rebase: `QuantumState[\[Psi], QuantumBasis["PauliX", n]]["ProbabilitiesList"]` |
|----:|----:|----:|
| 4 | 146 | 24 |
| 6 | 237 | 114 |
| 8 | 197 | 3829 |
| 10 | 225 | **~207000** |
| 12 | 261 | (OOM-prone) |

Agreement to $1.5\times10^{-16}$. Rebasing a state into a non-computational basis goes through `MatrixInverse` on a dense `ReducedMatrix` (`QuantumState/QuantumState.m:145-177`), an $O(4^n)$ dense change of basis that explodes to **207 seconds at $n=10$**. The correct pattern is to rotate the measurement basis into the computational one with *gates* inside the circuit (here, an $H$ on each qubit), then read computational-basis probabilities; that stays flat (TN apply) and beats the rebase by $\sim900\times$ at $n=10$. This is the same trick the QASM emitter uses for non-$Z$ measurements (`Inverse[B]` via `QuantumShortcut`).

---

## 4. Verdict on your two proposed options

> *Is it better to compute the state, then probabilities from the state? Or to put an inner product in a QCO, `QCO[{QS[...]["Dagger"], ..., QS[...]}]`?*

**Both are correct, both avoid the dilation, and they are not competitors: they answer different questions.**

| You want | Use | QF call | Cost on top of state prep |
|---|---|---|---|
| the whole distribution over a basis | **Route S** | `qc[]["ProbabilitiesList"]` (or `["Probabilities"]` for an outcome-keyed association) | $O(2^n)$, $\sim1$ ms |
| one expectation value $\langle O\rangle$ | **Route I** | `(\[Psi]^\[Dagger] @ O @ \[Psi])["Scalar"]`, or your sandwich `QuantumCircuitOperator[{\[Psi], O, \[Psi]^\[Dagger]}][]["Scalar"]` | bounded, $\sim135$ ms flat to $n=12$ |
| one overlap $\langle\phi|\psi\rangle$ / a single amplitude, without holding $\psi$ | **Route I** | `(qc /* {\[Phi]^\[Dagger]})[]["Number"]` | one contraction, $\sim40$ ms |
| one computational-basis outcome, state already built | **Route S** | `state["ProbabilitiesList"][[k]]` | lookup, $\sim1$ ms |
| probabilities in a non-computational basis | **Route S + gates** | append basis-change gates, then `["ProbabilitiesList"]` | flat (TN), not the dense rebase |

The mental model: the **inner product is "the framework computing one number for you in one contraction,"** while **state-then-probabilities is "the framework handing you the whole vector once, after which every readout is free."** Your `QCO[{bra, ..., ket}]` sandwich is exactly Route I, and it is verified to equal $\mathrm{Tr}[O\rho]$ to machine precision (appendix, probe 5). Use it for scalars; use the state for the whole picture.

Both routes do still build the $2^n$ amplitudes at some point: QF is a dense / tensor-network simulator, so neither beats the exponential, they only avoid the *extra* penalties (the $2^{2n}$ dilation, the $8^n$ eigensolve, the $4^n$ rebase). To go genuinely beyond $2^n$, leave this regime entirely: stabilizer tableau for Clifford + Pauli measurements (`PauliStabilizer`, `ps["Expectation", "ZZI"]`), light-cone-localized observables via the tensor-network engine (n-independent for local $\langle Z_c\rangle$ on a 1D chain; see the TN light-cone draft), or an MPS apply with bounded bond dimension. The escape routes are catalogued in the design note's Section 3.

The stabilizer escape route has a hard boundary that is exactly the `RY[theta]` layer in the benchmark circuit above: the stabilizer path is *Clifford only*. A verified study of what happens at that boundary is in [`../Stabilizer/StabilizerFrame-NonClifford-Blowup.md`](../Stabilizer/magic-and-pathsums/StabilizerFrame-NonClifford-Blowup.md). In short, a non-diagonal rotation `RY[theta]` returns `$Failed` (one `PauliStabilizer::nonclifford`), so this circuit cannot use the stabilizer route at all; a diagonal phase rotation `P[theta]` or a `T` gate is instead lifted to a `StabilizerFrame` whose component count **doubles per non-Clifford gate** (measured: exactly $2^k$ out to $k=18$), so it is usable only for a handful of $T$ gates. That report also found, and now fixes in the kernel, a materialization defect: the frame's dense `["StateVector"]` was phase-incoherent for genuine superpositions and is now phase-correct (verified amplitude-exact against the dense state). The remaining limit is scalability, not correctness: the frame does no stabilizer-rank compression, so its dense readout is a correct numeric output only for a handful of $T$ gates.

---

## 5. When you *should* use the measurement operator

Route M is not a mistake; it is the right object when you need the *physics of measurement*, not a statistic:

- **Collapsed / post-measurement states.** Mid-circuit measurement with classically-controlled follow-up: teleportation, error-correction syndrome extraction then correction, deferred-measurement experiments. The pointer register *is* the classical record you condition on. (See the feedforward showcase, `OngoingProjects/Why QF/QF-Feedforward-Pipeline.md`.)
- **The `QuantumMeasurement` outcome object's API.** `["MultivariateDistribution"]` for use inside `Probability[]`/`NProbability[]`, `["MixedStates"]`, `["Mean"]`/`["Entropy"]` when $n$ is small.
- **Shot statistics.** `qm["SimulatedMeasurement", shots]` / `Counts[...]` samples the exact distribution; it does *not* re-dilate per shot, so this is cheap and is the way to mimic QPU counts.

Even here, cap the memory: if you want the post-measurement state but not the explicit pointer, prefer `qmo["DiscardExtraQudits"]` (folds the pointer back out, cost $O(4^n)$ independent of how many measurements), the `"POVM"` form, or the `QuantumChannel` (Kraus) form. These keep you at $4^n$ for the whole register instead of letting pointers stack to $2^{2n}$.

---

## 6. Verified recipes (copy-paste, checked against `Tr[O\[Rho]]` / the full distribution)

```wolfram
Needs["Wolfram`QuantumFramework`"];
qc  = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "RY"[0.5] -> 2}];
psi = qc[];                              (* build the state once *)

(* --- full distribution (Route S) --- *)
psi["ProbabilitiesList"]                 (* {p_00, p_01, p_10, p_11} *)
psi["Probabilities"]                     (* <|"00" -> p_00, ...|>     *)

(* --- single computational-basis outcome (Route S, cheapest once psi exists) --- *)
psi["ProbabilitiesList"][[1]]            (* P(|00>) *)

(* --- single overlap / amplitude without holding psi (Route I) --- *)
amp00 = (qc /* {QuantumState["00"]["Dagger"]})[]["Number"];   (* <00|psi> *)
Abs[amp00]^2                                                   (* P(|00>)  *)

(* --- expectation value <O> (Route I, the route to use for large n) --- *)
O = QuantumOperator["ZZ"];
(psi["Dagger"][O[psi]])["Scalar"]                          (* psi^dagger . O . psi *)
QuantumCircuitOperator[{psi, O, psi["Dagger"]}][]["Scalar"]  (* your sandwich form, identical *)

(* --- probabilities in a non-computational (X) basis: rotate with gates --- *)
(qc /* QuantumCircuitOperator[{"H" -> 1, "H" -> 2}])[]["ProbabilitiesList"]
(* NOT QuantumState[psi, QuantumBasis["PauliX", n]] -- that is a dense O(4^n) rebase *)
```

Avoid, for anything but small $n$ or when you specifically need the measurement object:

```wolfram
QuantumMeasurementOperator[O][psi]["Mean"]                   (* expectation: 21 s at n=6, OOM at n=8 *)
QuantumState[psi, QuantumBasis["PauliX", n]]["ProbabilitiesList"]  (* X-basis: 207 s at n=10 *)
```

---

## 7. Caveats and provenance

- All numbers are single-run wall times on one M2 Pro after `ClearSystemCache[]`; treat them as order-of-magnitude, not benchmark-grade. The qualitative hierarchy (Route I/S flat; Route M-for-expectation super-exponential; dense rebase $\sim4^n$) is robust and reproduced across $n$.
- `eigensystem` being generic (no Hermitian path, with `Simplify`) is the documented root cause of the expectation-route blow-up (`Utilities.m:151`; `audit/architecture.md` §3, `audit/idioms-and-performance.md` "Eigendecomposition"). A Hermitian/diagonal-aware path in the measurement constructor would shrink Route M's penalty but not remove the dilation.
- The kernel was not modified. This is a usage/cost audit; the only "fix" is choosing the right call. The standing kernel-side opportunities (a `Method -> "Trajectories"` apply path, automatic pointer garbage-collection, a stabilizer/MPS dispatch for measurement) are in the design note's Sections 3.B-C.
- File:line anchors are against HEAD `f9dc1cdc`; re-verify with the live kernel if the repo has moved (per the QF audit ritual).

---

## Appendix: probe scripts

Saved under [`scripts/`](scripts/) next to this report: `qf_meas_probe.wls` (dimension table, full-distribution Route S vs M, agreement), `qf_meas_probe2.wls` (expectation three ways, single outcome, dilated-state dimension), `qf_meas_probe3.wls` (single outcome + X-basis rebase blow-up), `qf_meas_probe4.wls` (X-basis via gates vs rebase, bounded bra-op-ket), `qf_meas_probe5.wls` (idiom verification against `Tr[O\[Rho]]`). Each loads the working tree with `PacletDirectoryLoad` and times with `ClearSystemCache[]; AbsoluteTiming`.
