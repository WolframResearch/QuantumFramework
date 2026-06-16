# Measurement by Hilbert-space dilation in QuantumFramework: the blow-up and ways around it

**Status:** exploration / design note. Started 2026-06-10.
**Scope:** why QF realizes measurement by enlarging the Hilbert space, why that is physically faithful, where it becomes computationally catastrophic, and a menu of escape routes (some already in the kernel, some standard QIS techniques to adopt, some open).
**Empirical companion:** [`Computing-Probabilities-And-Expectations.md`](Computing-Probabilities-And-Expectations.md) runs the Section 4 experiments and turns them into a decision guide (state-then-probabilities vs inner-product vs measurement-operator; verified at HEAD `f9dc1cdc`).

Kernel anchors below are against paclet 2.0.0 (audit anchor SHA `ac32bda8`, 2026-06-09). Re-verify `file:line` against the live kernel before relying on any single line.

---

## 1. What QF does, and why it is correct

QF implements a projective measurement the way von Neumann's measurement scheme prescribes it: not as an instantaneous collapse postulate bolted onto the dynamics, but as a *unitary* interaction between the measured system and an added pointer system (the "measuring device").

Given a system in state $|\psi\rangle$ on a Hilbert space $\mathcal{H}_S$, and a target observable with spectral projectors $\{P_k\}$, QF builds the Stinespring/von Neumann dilation

$$
U_{\text{meas}} \;:\; |\psi\rangle_S \otimes |0\rangle_A \;\longmapsto\; \sum_k \bigl(P_k|\psi\rangle\bigr)_S \otimes |k\rangle_A ,
$$

where $\mathcal{H}_A$ is a freshly appended **pointer register** (one extra qudit per measurement) whose basis $\{|k\rangle\}$ is labelled by the eigenvalues. After this unitary, the reduced state of $S$ is exactly the post-measurement mixture $\sum_k P_k\rho P_k$, and reading the pointer reproduces the Born rule. This is the physically honest picture: a measurement *is* an entangling interaction with another quantum system, and "collapse" is what the system half looks like once the pointer is traced out.

This is genuinely close to the core of quantum theory, and it is what makes QF's measurement compose cleanly inside circuits: a measurement is just another operator, so a mid-circuit measurement, a deferred measurement, and an erased measurement are all the same object in different orderings.

### Where it lives in the kernel

- `QuantumMeasurementOperator/Properties.m:206-276` — the `"SuperOperator"` property. It partial-traces away the non-target qudits (`:220`), diagonalizes the target observable (`:222`), builds the pointer `eigenBasis` (`:226-232`), and tensors it on the output: `QuantumTensorProduct[eigenBasis, ...]` at `:241-245` and `:255-259`. The `eigenBasis` qudit is literally the added measuring device.
- `QuantumMeasurementOperator/Properties.m:81` — `"ExtraQudits"` counts the pointer registers (output orders with non-positive index).
- `QuantumMeasurementOperator/Properties.m:83` — `"Eigenqudits"`; `:86` — `"Eigendimensions"`.
- `audit/mistakes.md:635` — confirms each separate single-qubit measurement gets *its own* eigenvalue register: `{1}, {2}` in a circuit is two pointer qudits, not one joint register.

---

## 2. The blow-up, quantified

The faithfulness has a price: every measurement *adds a qudit instead of removing one*. In a classical statevector simulation, a projective collapse should let you keep tracking $2^n$ amplitudes (or fewer, if you commit to one branch). The dilation does the opposite.

Let the system have $n$ qubits and let us measure $k$ of them in the computational basis (each pointer of dimension $2$).

| representation | no measurement | after $k$ dilated measurements |
|---|---|---|
| pure statevector length | $2^{n}$ | $2^{\,n+k}$ |
| density matrix entries | $4^{n}$ | $4^{\,n+k}$ |
| measure-everything worst case ($k=n$), statevector | $2^{n}$ | $2^{2n}=4^{n}$ |

So measuring all $n$ qubits and keeping the dilated pure state *squares the length of the statevector*. In a circuit with mid-circuit measurements the registers accumulate: $m$ rounds of single-qubit measurement on the same line append $m$ pointer qudits, and because the object stays a single pure unitary evolution (no branch is ever discarded), there is no point at which the dimension comes back down. This is the "computational catastrophe very fast" you hit: a few mid-circuit measurements on a modest register push the dense statevector past memory.

The HHL precision crash already on record (`project_qf_2_0_0_hhl_tn_crash`, default `TensorNetwork.m:215` apply path bloating output qudits) is the same disease in a different organ: an apply path that keeps appending qudits it never contracts away.

---

## 3. Ways around it

Grouped by how ready each one is. **(A)** already exists in the kernel and just needs to be the default / documented; **(B)** standard QIS techniques QF could adopt; **(C)** open research.

### A. Already in QF — use or default to these

1. **`"DiscardExtraQudits"` — trace the pointer back out.**
   `QuantumMeasurementOperator/Properties.m:300-315`. After the dilation it folds a `"Trace"` over each eigenbasis register, returning a `QuantumOperator` on the original target order. This is the channel $\rho \mapsto \sum_k P_k\rho P_k$ with no surviving pointer. Cost is $O(4^{n})$ (density matrix), independent of $k$: it never lets the registers stack. **This is the single most important escape hatch** and the natural default for "I want the post-measurement state, not a record."

2. **Apply as `QuantumMeasurement` (the collapse/outcomes object), not as a dilated operator.**
   `QuantumMeasurement.m`. Measuring a *state* returns a `QuantumMeasurement` carrying outcome probabilities and the post-measurement (mixed) states keyed by eigenvalue (`"States"` `:141`, `"MixedStates"` `:127`, `"Mean"` `:223`). Here the pointer is summarized as a classical distribution rather than tensored on. The blow-up is specifically a property of keeping the measurement as a *unitary inside a circuit*; pulled out to the state level it collapses to the outcome object.

3. **POVM / `QuantumChannel` (Kraus) form — stay in $O(d^2)$, never dilate the statevector.**
   `"POVM"` (`:278`) and the `QuantumChannel` path (`"Bend"` at `:286-294`). Working in the density-matrix / Kraus picture caps the cost at $4^{n}$ for the whole register regardless of how many measurements occur, because measurement is a channel on a fixed space, not a dimension-adding unitary. This is the standard answer for "many measurements, I don't need the explicit pointer."

4. **`Method -> "Schrodinger"` to dodge the TN apply-path bloat** when the dilation interacts with the tensor-network default that itself adds qudits (`project_qf_2_0_0_hhl_tn_crash`).

### B. Standard QIS techniques to adopt / surface

5. **Principle of deferred measurement.** Push all measurements to the end of the circuit. Mathematically the dilation is unchanged, but with measurements last you can read the pointer once and never propagate the enlarged state through subsequent gates. QF already has the ingredients (`"Shift"` `:283`, `"Bend"` `:286`, circuit reordering); the missing piece is an automatic rewrite that recognizes a measurement followed by classically-controlled gates and defers it.

6. **Quantum trajectories / Monte-Carlo wavefunction.** Instead of one dilated pure state holding *all* branches at dimension $2^{n+k}$, sample: at each measurement draw an outcome $k$ with probability $\langle\psi|P_k|\psi\rangle$, project, renormalize, and continue with a *single* $2^{n}$ statevector. Average over shots. This trades the $4^n$ memory wall for $O(\text{shots}\times 2^n)$ time and is the textbook way real simulators handle mid-circuit measurement. A `Method -> "Trajectories"` (or shot-based) apply path is the highest-value addition here.

7. **Stabilizer / Clifford fast path.** For Clifford circuits with Pauli measurements the whole thing is polynomial (Gottesman-Knill). QF already has `PauliStabilizer` / `StabilizerStateQ` / `CliffordChannel` (see `reference_paulistabilizer_deep_audit`). Routing measurement through the stabilizer tableau when the circuit is Clifford sidesteps the dilation entirely. The connective tissue (detect Clifford + Pauli-measurement, dispatch to the tableau) is the work.

8. **Tensor-network / MPS apply with bounded bond dimension.** Mid-circuit measurement on an MPS is local: project one site, renormalize, optionally truncate. Cost stays $O(n\chi^3)$ rather than $4^n$, provided the apply path *contracts* rather than appends (contrast `TensorNetwork.m:215`). This is the same fix the HHL crash needs.

### C. Open / research directions

9. **Lazy pointer registers.** Keep the pointer symbolic (a labelled classical register) and only materialize the tensor factor if and when a later gate genuinely acts on it (classical control, erasure experiments). Most pointers are read and discarded; they never need a Hilbert dimension.

10. **Garbage collection of dead registers in circuits.** A static pass that detects pointer qudits which are never subsequently coupled and rewrites them to `"DiscardExtraQudits"` automatically, so the circuit object self-limits its dimension.

11. **Cost-model-driven `Method -> Automatic`.** Like NDSolve's nested method selection: pick statevector-dilation vs. channel vs. trajectories vs. stabilizer vs. MPS from the circuit's measurement count, Cliffordness, and entanglement, rather than always dilating. (Mirrors the `qf_tn_contraction_path` recommendation of an NDSolve-style nested `Method`.)

---

## 4. Suggested experiments (to ground the doc with real numbers)

Run on the live kernel, `ClearSystemCache[]` before timing:

1. Build an $n$-qubit register, measure $k$ qubits mid-circuit as a dilated `QuantumCircuitOperator`, and record statevector dimension and wall-time vs. $k$ for $n = 4,\dots,10$. Confirm the $2^{n+k}$ growth and locate the practical wall.
2. Same circuit via `"DiscardExtraQudits"` and via the `QuantumChannel` path; confirm the cost flattens to $4^n$.
3. Same circuit sampled as a trajectory by hand (project + renormalize per measurement); confirm $O(2^n)$ memory.
4. A Clifford version through `PauliStabilizer`; confirm polynomial scaling.

These four numbers turn this note into a benchmark figure and tell us which escape route to make the default.

---

## 5. Open questions to verify against the kernel

- Does applying a `QuantumMeasurementOperator` *inside* a `QuantumCircuitOperator` always keep the dilated pure form, or does it ever auto-collapse? (Check the circuit apply dispatch, not just `QuantumMeasurement.m`.)
- Is there already an automatic deferral or register-discard anywhere in the circuit compiler? (grep the circuit `Properties.m` / `IBMQuantum.m` transpile paths.)
- Does the `QuantumChannel` path's `"Bend"` truly avoid the pointer tensor factor, or does it dilate internally?

---

## Appendix: file-line index

- `QuantumMeasurementOperator/Properties.m:206-276` — `"SuperOperator"` dilation (pointer construction).
- `QuantumMeasurementOperator/Properties.m:81-86` — `"ExtraQudits"`, `"Eigenqudits"`, `"Eigendimensions"`.
- `QuantumMeasurementOperator/Properties.m:278` — `"POVM"`.
- `QuantumMeasurementOperator/Properties.m:286-294` — `"Bend"` (channel form).
- `QuantumMeasurementOperator/Properties.m:300-315` — `"DiscardExtraQudits"`.
- `QuantumMeasurement.m:127,141,223` — outcome / mixed-state accessors.
- `audit/mistakes.md:635` — sequential vs. joint measurement registers.
- `audit/PROGRESS.md:33` — Stinespring dilation construction note.
- `TensorNetwork.m:215` — the TN apply path that similarly appends qudits (HHL crash).
