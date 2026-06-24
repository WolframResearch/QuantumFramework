# StabilizerFrame Observables: Design and Verified Out-of-Tree Prototype

**Status.** Investigation and out-of-tree prototype only. No kernel `.m` file is modified in this pass. Every numerical claim below was produced by `wolframscript -file` against the working-tree paclet (`PacletDirectoryLoad`), with no `Quiet` and no `Print`. Kernel anchor: live `HEAD` `b81f3332`. The three cited files (`InnerProduct.m`, `StabilizerFrame.m`, `GateUpdates.m`) are unchanged since the audit anchor `b4fcdb31`; the only intervening `Kernel/Stabilizer/` change is an unrelated `SymbolicMeasure.m` update, so the cited `file:line` numbers are current.

Scripts (all under `OngoingProjects/Stabilizer/stabilizerframe-observables-scripts/`):

| File | Role |
|---|---|
| `framealgebra.wls` | the prototype library (symplectic Pauli algebra + the frame observables) |
| `verify_core.wls` | expectation, norm, shared-reference inner product, Born vs dense |
| `verify_amplitude.wls` | reference amplitudes, full complex amplitude, computational Born vs dense |
| `verify_crossref.wls` | distinct-reference inner product, shortcoming 3 demonstration |
| `bench_scaling.wls`, `bench_densify.wls` | rank path vs densification, large-$n$ correctness anchor |

---

## 1. The object, and the three gaps

A `StabilizerFrame` is the Quantum Framework representation of a state past the Clifford boundary: a superposition
$$
|\psi\rangle \;=\; \sum_{i=1}^{\chi} c_i\,|s_i\rangle
$$
of stabilizer states $|s_i\rangle$ with complex coefficients $c_i$. The component count $\chi$ is the stabilizer rank; each non-Clifford $T$ or $P[\theta]$ gate at most doubles it, so a few-$T$ state stays small ($\chi = 2^{t}$ for $t$ such gates) independent of the qubit number $n$. The internal form is `StabilizerFrame[<|"Components" -> {{c_i, ps_i}, ...}, "Paulis" -> {...}|>]`.

The frame is the right data structure for $T$-rich circuits, but today you cannot read a single number off it without paying $2^n$. The three gaps, each verified live:

**Gap 1: no native Pauli expectation or measurement.** The property contract `_StabilizerFrame["Properties"]` (`Kernel/Stabilizer/StabilizerFrame.m:47`) lists only `Components, Coefficients, Stabilizers, Length, Qubits, GeneratorCount, StateVector, State, InnerProduct`. There is no `"Expectation"` and no `"M"`. Calling `frame["Expectation", "X"]` or `frame["M", "X"]` returns **inert**: the only string-keyed multi-argument rule is the gate-update `f_StabilizerFrame[gate_String, args : ___Integer]` (`StabilizerFrame.m:201`), whose `args : ___Integer` guard rejects a string second argument such as `"X"`, so nothing fires. The only way to get an observable is to densify to the $2^n$ state vector.

**Gap 2: the frame inner product densifies.** `stabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame]` (`Kernel/Stabilizer/InnerProduct.m:111`) is literally
```
vA = fA["StateVector"];   (* builds a 2^n vector *)
vB = fB["StateVector"];   (* builds a 2^n vector *)
Conjugate[vA] . vB
```
The comment one line above (`InnerProduct.m:110`) states the efficient formula $\sum_{ij}\overline{c_i}\,c_j\,\langle s_i|s_j\rangle$, but the code does not use it. A closed-form pairwise stabilizer overlap **already exists**: `stabilizerInnerProductClosedForm[psi, phi]` (`InnerProduct.m:50`), $O(n^3)$ per pair via Aaronson-Gottesman row-sum phases. Its one limitation is documented in place (`InnerProduct.m:41`): it returns the magnitude (and the $\pm 1$ sign for perfectly aligned states) but not the complex phase, which "requires Gaussian-sum bookkeeping over the intersection generators (TODO)."

**Gap 3: the naive $\sum_i c_i\,|s_i\rangle$ from `Components` does not reconstruct the state.** The stored component stabilizer states each carry their own arbitrary global phase, so summing them with the stored coefficients gives the wrong state. The relating Paulis in the `"Paulis"` key fix this. We demonstrate the failure quantitatively in §2.

A single-state Pauli expectation already exists and is exact, with phase: `stabilizerExpectation[ps_PauliStabilizer, P_String]` (`InnerProduct.m:156`) returns $\langle s|P|s\rangle \in \{0,+1,-1\}$ for a single stabilizer state $|s\rangle$. This primitive, together with the relating Paulis, is all that is needed to close gaps 1 and 3 exactly, and gap 2 for the common case, in $\chi^2\,\mathrm{poly}(n)$ time. That is the content of this report.

---

## 2. The relating-Pauli reconstruction (gap 3), reproduced

A gate-built frame carries the `"Paulis"` key: one relating Pauli per component, encoded $\{x\text{-bits}, z\text{-bits}, \mathrm{coeff}\}$ with $\mathrm{coeff}\in\{1,-1,i,-i\}$ (`StabilizerFrame.m:73`). Component 1 is the **reference** and its relating Pauli is the identity. The contract is
$$
|s_i\rangle \;=\; \mathrm{matrix}[P_i]\,|\mathrm{ref}\rangle ,
\qquad |\mathrm{ref}\rangle = \text{component 1's stabilizer state}.
$$
The decisive structural fact: a $Z$ doubling leaves the stabilizer tableau unchanged, and a Clifford update transforms every component identically, so **all $\chi$ components share one tableau** and differ only by the relating Pauli applied to the single reference. The non-Clifford boundary that builds the frame (`Kernel/Stabilizer/GateUpdates.m:182`) sets component 1 to the original stabilizer $|\mathrm{ref}\rangle$ and component 2 to $Z_j|\mathrm{ref}\rangle$, with relating Paulis $\{0,0,1\}$ (identity) and $\{0, e_j, 1\}$ ($Z_j$); subsequent gates conjugate or $Z$-multiply these.

`frame["StateVector"]` (`StabilizerFrame.m:174`) therefore materializes the reference **once** and relates the rest:
$$
|\psi\rangle \;=\; \sum_{i} c_i\,\mathrm{matrix}[P_i]\,|\mathrm{ref}\rangle
\;=\;\Big(\sum_i c_i\,P_i\Big)\,|\mathrm{ref}\rangle .
$$
Here each relating Pauli, via `sfApplyPauli` (`StabilizerFrame.m:147`), denotes the operator $\mathrm{coeff}\cdot i^{\,x\cdot z}\,X^{x}Z^{z}$, a genuine signed-and-phased Pauli. So the whole frame is **one Pauli-sum operator applied to one reference stabilizer state**.

We reproduce this exactly in the prototype with `pauliApplyVec` (the `sfApplyPauli` convention) and confirm both halves of gap 3 (`verify_crossref.wls`, 32 random Clifford+$T$ frames, $n\le 4$):

| quantity | worst $|\,\cdot\,|$ error vs `frame["StateVector"]` |
|---|---|
| naive $\sum_i c_i\,|s_i\rangle$ over stored components | $1.41$ (i.e. $\sqrt{2}$: the trap is real) |
| relating-Pauli reconstruction $\sum_i c_i\,P_i\,|\mathrm{ref}\rangle$ | $0$ (exact) |

The naive sum is off by a full $\sqrt 2$ in amplitude; the relating-Pauli reconstruction is bit-exact. Every observable below threads the relating Paulis.

---

## 3. The core reduction: Pauli expectation on a frame (gap 1)

With $|\psi\rangle = \sum_i c_i P_i|\mathrm{ref}\rangle$ and a Hermitian Pauli observable $P$,
$$
\langle\psi|P|\psi\rangle
\;=\;\sum_{ij}\overline{c_i}\,c_j\,\big\langle \mathrm{ref}\big|\,P_i^{\dagger}\,P\,P_j\,\big|\mathrm{ref}\big\rangle .
$$
Every factor is a Pauli, so $Q_{ij} \equiv P_i^{\dagger}P P_j$ is a Pauli times a scalar phase in $\{1,i,-1,-i\}$. Writing $Q_{ij} = \sigma_{ij}\,\hat Q_{ij}$ with $\hat Q_{ij}$ a Hermitian Pauli string,
$$
\big\langle \mathrm{ref}\big|Q_{ij}\big|\mathrm{ref}\big\rangle
\;=\;\sigma_{ij}\,\underbrace{\big\langle \mathrm{ref}\big|\hat Q_{ij}\big|\mathrm{ref}\big\rangle}_{\in\{0,+1,-1\}} ,
$$
and the underbraced quantity is exactly what `stabilizerExpectation[ref, ...]` computes, with its full Aaronson-Gottesman sign. **No densification, no Gaussian-sum phase, no gate application.** The phases come entirely from exact Pauli multiplication (the $i$ factors of $XZ=-iY$ and the $(-1)$ of $ZX$); the magnitudes and signs come from the existing exact single-state primitive.

The bookkeeping is a small symplectic Pauli algebra (operator $s\,X^{x}Z^{z}$, scalar $s$):
$$
g(x,z,s)\,g(x',z',s') = g\!\big(x\!\oplus\! x',\,z\!\oplus\! z',\,s s'(-1)^{z\cdot x'}\big),
\qquad
g(x,z,s)^{\dagger} = g\!\big(x,z,\ \overline{s}\,(-1)^{x\cdot z}\big).
$$
The relating Pauli $\{x,z,\mathrm{coeff}\}$ enters as $s = \mathrm{coeff}\cdot i^{x\cdot z}$. This is `frameSandwich` in `framealgebra.wls`; `frameExpectation[f, P] = frameSandwich[f, P, f]`.

**Cost.** Exactly $\chi^2$ Pauli products (each $O(n)$) and $\chi^2$ calls to `stabilizerExpectation` (each $O(n^{2\text{-}3})$). Total $\chi^2\,\mathrm{poly}(n)$, independent of $2^n$.

**Verification** (`verify_core.wls`, seed 2026): 72 random Clifford+$T$ frames, all Pauli strings for $n\le 3$ and 12 sampled strings for $n=4$, against the dense $\langle P\rangle = \overline{v}\cdot(P v)$ with $v = $ `frame["StateVector"]`:

| check | cases | worst absolute error |
|---|---|---|
| $\langle\psi|P|\psi\rangle$ vs dense | 72 | $5.1\times10^{-16}$ |
| $\langle\psi|\psi\rangle$ (norm) vs dense | 32 | $3.3\times10^{-16}$ |

A representative single point: on $T|+\rangle$ the prototype returns $\langle X\rangle = \langle Y\rangle = 1/\sqrt2$, $\langle Z\rangle = 0$, the correct Bloch position, matching dense to $5\times10^{-16}$.

---

## 4. Inner product without densifying (gap 2)

The same sandwich generalizes to two frames:
$$
\langle\psi_A|\psi_B\rangle = \sum_{ij}\overline{a_i}\,b_j\,\big\langle \mathrm{ref}_A\big|\,P^{A\,\dagger}_i\,P^{B}_j\,\big|\mathrm{ref}_B\big\rangle .
$$

**Shared reference (the common case).** A single frame's norm, and any two frames produced from a common reference (for example a frame and the same frame with one extra $T$, or a scalar multiple), have $\mathrm{ref}_A = \mathrm{ref}_B$. Then each term collapses to a single-state Pauli expectation exactly as in §3: **exact, phase-correct, $\chi^2\,\mathrm{poly}(n)$.** This is `frameInnerProduct`.

| check | cases | worst absolute error |
|---|---|---|
| shared-reference $\langle\psi_A|\psi_B\rangle$ vs dense (`verify_core`) | 32 | $4.2\times10^{-16}$ |
| general reduction agrees with shared-ref specialization (`verify_crossref`) | 32 | $0$ |

**Distinct references.** When $\mathrm{ref}_A \neq \mathrm{ref}_B$, each term $\langle \mathrm{ref}_A|Q_{ij}|\mathrm{ref}_B\rangle$ is a Pauli transition amplitude between two different stabilizer states. Here a subtlety bites, and it is worth stating because the brief's hypothesis ("just sum the existing closed-form overlaps") is not quite enough: applying a Pauli **gate** to a `PauliStabilizer` is phase-blind ($Z|1\rangle$ returns $|1\rangle$, not $-|1\rangle$), so one cannot read $\langle \mathrm{ref}_A|Q|\mathrm{ref}_B\rangle$ off a gate-applied tableau. Two honest routes:

- **Exact (oracle).** Apply $Q$ to $\mathrm{ref}_B$'s vector with explicit phase, then dot with $\mathrm{ref}_A$'s vector. Phase-correct but $2^n$ per pair. This is a correctness oracle, not a speed win, and we use it to confirm the reduction **formula** is right. `verify_crossref.wls`, 32 random distinct-reference frame pairs: worst error vs the densifying `frame["InnerProduct"]` is $1.6\times10^{-16}$.
- **Poly magnitude.** $\big|\langle \mathrm{ref}_A|Q|\mathrm{ref}_B\rangle\big|$ depends only on the tableaux, not the dropped global phase, so it is recovered exactly by `stabilizerInnerProductClosedForm` on the (phase-blind) gate-applied tableau, $O(n^3)$ per pair. Verified per-pair against the dense magnitude: worst error $0$ over the same sweep.

The one genuinely missing piece for a poly **and** phase-correct distinct-reference overlap is the per-pair complex phase, which is exactly the Gaussian-sum TODO already named at `InnerProduct.m:41`. Completing that one closed-form phase turns the distinct-reference inner product poly as well. Crucially, the **shared-reference** path (norm, expectation, same-family overlaps) needs none of it: it routes through `stabilizerExpectation`, which already tracks phase exactly.

---

## 5. Born probabilities and measurement (gap 1, derived)

Everything measurement-related follows from the two primitives above, with no new machinery.

**Pauli measurement.** For a Hermitian Pauli $P$ (eigenvalues $\pm1$) with projectors $(\mathbb{1}\pm P)/2$,
$$
\Pr(\pm) \;=\; \frac{\langle\psi|\psi\rangle \pm \langle\psi|P|\psi\rangle}{2\,\langle\psi|\psi\rangle}.
$$
**Single-qubit computational marginal.** $\Pr(x_q) = \big(\langle\psi|\psi\rangle + (-1)^{x_q}\langle Z_q\rangle\big)/(2\langle\psi|\psi\rangle)$.

`verify_core.wls`, 20 random frames:

| check | result |
|---|---|
| $\Pr(+)+\Pr(-) = 1$ (and marginal sums) | worst deviation $0$ |
| probabilities real | worst imaginary residue $0$ |
| probabilities $\ge 0$ | minimum $0$ |

**Full computational string.** $\Pr(x) = |\langle x|\psi\rangle|^2$ (see §6) is exact and $2^n$-free.

---

## 6. Amplitudes and computational Born

The relating-Pauli form gives a computational amplitude as a $\chi$-term sum over reference amplitudes:
$$
\langle x|\psi\rangle = \sum_i c_i\,\langle x|P_i|\mathrm{ref}\rangle
= \sum_i c_i\,s_i\,(-1)^{\,z_i\cdot(x\oplus x_i)}\,\langle x\oplus x_i\,|\,\mathrm{ref}\rangle ,
$$
so only single amplitudes of the **reference** stabilizer state are needed. Those are obtained from one anchor $x_0$ in the support plus a poly phase ratio: for $y, x_0$ in the support, picking any group element $g = s_g\,X^{a}Z^{b}\in S$ with $a = y\oplus x_0$,
$$
\frac{\langle y|\mathrm{ref}\rangle}{\langle x_0|\mathrm{ref}\rangle} = \langle y|g|\mathrm{ref}\rangle/\langle x_0|\mathrm{ref}\rangle = s_g\,(-1)^{\,b\cdot x_0},
$$
since $g|\mathrm{ref}\rangle = |\mathrm{ref}\rangle$. The anchor magnitude is $|\langle x_0|\mathrm{ref}\rangle| = 2^{-k/2}$ with $k$ the $X$-rank of the reference tableau. Therefore:

- **Computational Born $\Pr(x) = |\langle x|\psi\rangle|^2$ is fully poly and gauge-free** (the anchor magnitude $2^{-k/2}$ is all that is needed).
- The **full complex amplitude** $\langle x|\psi\rangle$ additionally needs the anchor's phase, which is the reference state's global-phase convention. In a kernel implementation this is the standard CH-form global phase; in the prototype we read the single anchor amplitude from the reference vector (the one place that touches $2^n$, purely to pin the convention) and confirm the rest is reproduced by the poly ratios.

`verify_amplitude.wls`, seed 7, against `frame["StateVector"]` entrywise:

| check | cases | worst absolute error |
|---|---|---|
| reference single-amplitude ratio vs dense reference vector | 20 | $0$ |
| full complex $\langle x|\psi\rangle$ vs `frame["StateVector"]` | 36 | $1.2\times10^{-16}$ |
| computational $\Pr(x)$ vs $|\langle x|\psi\rangle|^2$ | 24 | $1.1\times10^{-16}$ |
| $\sum_x \Pr(x)$ vs $\langle\psi|\psi\rangle$ | 24 | $3.1\times10^{-16}$ |
| $\Pr(x)\ge 0$ | 24 | minimum $0$ |

---

## 7. Scaling: the rank path beats densification

The reductions in §§3-5 never build a $2^n$ object. `bench_scaling.wls` (seed 11) times `frameExpectation` on few-$T$ frames ($\chi = 8$, three $T$ gates, a depth-$n$ random Clifford layer), with `ClearSystemCache[]` before each timing:

| $n$ | 4 | 6 | 8 | 10 | 12 | 16 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|
| $\chi$ | 8 | 8 | 8 | 8 | 8 | 8 | 8 | 8 | 8 |
| rank `frameExpectation` (ms) | 14.7 | 14.5 | 15.2 | 15.6 | 19.4 | 23.4 | 18.7 | 22.3 | 29.4 |

The rank path is essentially flat from $n=4$ to $n=40$: the cost is $\chi^2\,\mathrm{poly}(n)$, with the mild growth coming from the $\mathrm{poly}(n)$ inside `stabilizerExpectation`, not from $2^n$.

For contrast, today the only way to get $\langle P\rangle$ is to densify the frame to its $2^n$ state vector first (`frame["StateVector"]`). `bench_densify.wls` times that materialization against `frameExpectation` on the same few-$T$ frames:

| $n$ | 4 | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|---|
| rank `frameExpectation` (ms) | 14.1 | 14.4 | 15.0 | 15.1 | 19.4 | 21.5 |
| densify `frame["StateVector"]` (ms) | 2.0 | 4.2 | 26.2 | 586 | 19722 | 208453 |

Densification crosses the rank path near $n=8$ and is already four orders of magnitude slower by $n=14$ ($208$ s vs $21$ ms). The rank path stays in the tens of milliseconds and runs comfortably at $n=40$ (§7 table), where building a $2^{40}$ vector is out of the question. The current `frame["StateVector"]` cost above is inflated further by a non-vectorized per-amplitude loop in the kernel, but even an optimal densification is $\Theta(2^n)$ and cannot follow the rank path past a few tens of qubits.

**Large-$n$ correctness anchor (no dense reference).** To check correctness where densification is impossible, `bench_scaling.wls` builds qubits $1,2,3$ as $T|+\rangle$ and the rest as $|0\rangle$, an $n$-independent state with $\langle Y_1Y_2Y_3\rangle = 2^{-3/2}$, $\langle Z_1\rangle = 0$, $\langle Z_4\rangle = +1$, norm $1$. The prototype reproduces all four to absolute error $0$ at $n = 6, 10, 20, 40, 60$ ($\chi = 8$).

---

## 8. Verification summary

| gap | prototype | what was checked | worst error | scaling |
|---|---|---|---|---|
| 3 | relating-Pauli reconstruction | vs `frame["StateVector"]`, 32 frames | $0$ (naive sum: $\sqrt2$ off) | $\chi$ |
| 1 | `frameExpectation` | $\langle P\rangle$ vs dense, 72 frames | $5.1\times10^{-16}$ | $\chi^2\,\mathrm{poly}(n)$ |
| 2 | `frameNorm2` / shared-ref `frameInnerProduct` | vs dense, 32+32 cases | $4.2\times10^{-16}$ | $\chi^2\,\mathrm{poly}(n)$ |
| 2 | distinct-ref reduction (exact oracle) | formula vs densifying IP, 32 pairs | $1.6\times10^{-16}$ | $\chi^2\,2^n$ per pair |
| 2 | distinct-ref per-pair magnitude (closed form) | vs dense magnitude | $0$ | $\chi^2\,\mathrm{poly}(n)$, magnitude only |
| 1 | Pauli + marginal Born | sum-to-one, real, $\ge 0$, 20 frames | $0$ | $\chi^2\,\mathrm{poly}(n)$ |
| 1 | `frameAmplitude` / computational $\Pr(x)$ | vs dense, 36+24 cases | $1.2\times10^{-16}$ | $\chi\,\mathrm{poly}(n)$ per $x$ |
| 1 | large-$n$ anchor | $n$ up to 60, no dense | $0$ | $n$-independent |

---

## 9. Proposed kernel methods (deferred, do not implement this pass)

The prototype is a thin layer over machinery the kernel already has. The eventual change is small and localized.

**`Kernel/Stabilizer/StabilizerFrame.m`.**
1. Add the symplectic Pauli helpers (multiply, dagger, "relating Pauli to operator", "Pauli string to symplectic", "symplectic to Hermitian string") as `PackageScope` functions, or factor them next to `sfApplyPauli` (`StabilizerFrame.m:147`), which already uses the same encoding.
2. Add, beside the existing `f_StabilizerFrame["InnerProduct", other_]` dispatch (`StabilizerFrame.m:55`):
   - `f_StabilizerFrame["Expectation", obs_String]` $\to$ the §3 sandwich, requiring the `"Paulis"` key (gate-built frame) and a concrete reference; reuse `stabilizerExpectation` (`InnerProduct.m:156`) for each $\langle\mathrm{ref}|\hat Q_{ij}|\mathrm{ref}\rangle$.
   - `f_StabilizerFrame["M", obs_String]` and `f_StabilizerFrame["Probability", obs_String]` $\to$ the §5 Born formulas, derived from `"Expectation"` and the norm.
   - optionally `f_StabilizerFrame["Amplitude", x_]` and a computational `"Probabilities"` from §6.
3. Add `"Expectation"`, `"M"`, `"Probability"`, `"Amplitude"` to the property contract list (`StabilizerFrame.m:47`).

**`Kernel/Stabilizer/InnerProduct.m`.** Replace the densifying body of `stabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame]` (`InnerProduct.m:111`):
- if both frames carry `"Paulis"` and share the reference tableau, use the §4 shared-reference reduction (exact, poly) via `stabilizerExpectation`;
- otherwise fall back to the per-component double sum $\sum_{ij}\overline{a_i}b_j\langle\mathrm{ref}_A|Q_{ij}|\mathrm{ref}_B\rangle$, threading an optional `Method`: the default per-pair `stabilizerInnerProduct` stays exact (at $2^n$ per pair, the present cost), and `Method -> "ClosedForm"` gives the poly per-pair magnitude. The existing densify can remain as the final fallback for hand-built frames with no `"Paulis"` key, where there is no rank structure to exploit and densifying once is already optimal.

The single piece that would let the distinct-reference inner product be poly **and** phase-correct is the Gaussian-sum phase already flagged at `InnerProduct.m:41`. It is independent of, and not required by, the shared-reference methods above.

No new public symbols; no `PacletInfo` change; the additions are property methods on an existing head plus one rewritten internal branch.

---

## 10. Out of scope: rank compression

This pass addresses **reading observables off a frame of whatever rank it currently has**. It does not address a separate, larger gap: the frame keeps all $2^t$ terms after $t$ non-Clifford gates (each $T$ doubles $\chi$, with no merging or deduplication), whereas the stabilizer rank of a magic state is bounded by the Bravyi-Gosset estimate $\chi \lesssim 2^{\,0.4\,t}$. Compressing a frame toward that bound is a stabilizer-decomposition search (matching pursuit over stabilizer states), orthogonal to the observable machinery here. Our reductions cost $\chi^2\,\mathrm{poly}(n)$ on the actual $\chi$; halving $\chi$ by compression would quarter that cost but is not needed for correctness, and is left for a future pass.
