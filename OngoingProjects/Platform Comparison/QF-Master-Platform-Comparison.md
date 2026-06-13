# QuantumFramework vs the Quantum Software Landscape: A Verified, Benchmarked Comparison

**Date:** 2026-06-12 (full revision: Section 3.2 stabilizer rows re-measured after the June 2026 optimization shipped; Section 12 rewritten as a plain-language account of that work, with before/after examples and run-it-yourself commands)
**QuantumFramework version:** 2.0.0, repo working tree at HEAD `02cc92d9` (verified live, not the lagging installed paclet). Feature claims were verified at `77b9ccba`; the commits since then touch only the stabilizer subsystem and the Liouvillian basis convention, both re-verified here.
**Author's note:** This document merges and supersedes three earlier, now-stale comparison sources (see Section 11). Every capability claim about QuantumFramework here was re-verified against the live 2.0.0 kernel; every performance number was measured head-to-head on one machine. Nothing in this file is carried over from the prior documents without re-checking.

---

## Evidence tiers (read this first)

Every claim carries one of three tags so you can tell measurement from citation:

- **[BENCH]** : I ran it head-to-head on this machine (Section 3, Section 10). Numbers are wall-clock minimums.
- **[SRC]** : verified by executing it live or reading the source. For QuantumFramework this means the 2.0.0 kernel (loaded from the repo working tree) or its `Kernel/*.m` files; for another package it means its installed code or repository source.
- **[DOC]** : taken from official documentation or a second-hand audit, **not** executed locally. Treat as a claim, not a measurement.

When a cell in a table is unlabeled it inherits the tag stated in that table's caption.

---

## 0. Executive summary

QuantumFramework is a **symbolic-first, object-algebra** quantum framework embedded in the Wolfram Language. The packages it is usually compared against are **numerically-specialized engines**: Qiskit/Aer and Cirq for gate-model state vectors, Stim and QuantumClifford.jl for Clifford/stabilizer circuits, QuTiP for open-system dynamics, PennyLane for differentiable quantum machine learning. The honest one-line summary is that these solve different problems, and the benchmark numbers below show exactly where that boundary falls.

Three head-to-head benchmarks on identical tasks (Apple M2 Pro, single machine, Section 3):

1. **State-vector simulation** of a fixed 105-gate, $n$-qubit circuit. QuantumFramework's default engine is **about $200\times$ to $400\times$ slower** than Qiskit Aer or Cirq (e.g. $n=16$: QF $3.4\,\mathrm{s}$ vs Aer $8.6\,\mathrm{ms}$). **[BENCH]**
2. **Stabilizer simulation** of an identical Clifford gate stream. After the June 2026 rework of the stabilizer subsystem (Section 12.1), QuantumFramework is **within $2.5\times$ to $5.3\times$ of Stim** and **beats QuantumClifford.jl at $n=1000$** ($n=1000$, $2\times10^4$ gates: QF $58\,\mathrm{ms}$, QC.jl $108\,\mathrm{ms}$, Stim $11\,\mathrm{ms}$). The pre-rework implementation was $660\times$ to $5800\times$ slower than Stim on the same task; the cumulative speedup is $266\times$ to $1095\times$. QuantumFramework's stabilizer results are **numerically cross-validated against both Stim and QuantumClifford.jl** by its own test suite. **[BENCH]/[SRC]**
3. **Open-system (Lindblad) dynamics** of a driven, damped qubit. Here QuantumFramework is only **about $9\times$ slower** than QuTiP ($17\,\mathrm{ms}$ vs $1.9\,\mathrm{ms}$) and the two agree on the physics to $4\times10^{-4}$. This is the regime where QuantumFramework is genuinely competitive, because it is doing the same thing QuTiP does (numerical ODE integration) on a small system. **[BENCH]**

The shape after June 2026: QuantumFramework's value remains the symbolic object model (exact algebra, free parameters, qudits, ZX, second quantization, phase space, basis reconciliation) and seamless Wolfram-Language integration, but the performance story is no longer uniform. The **stabilizer engine is now genuinely competitive** with the dedicated tools (single-digit factors from Stim, ahead of QC.jl at large $n$). The **state-vector path remains two orders of magnitude behind** Aer/Cirq; a completed audit (Section 12.2) explains exactly why, has prototype-validated fixes that would close most of the gap, and gives user-side levers available today. Open-system work is competitive on small systems already.

A second headline: **both prior comparison documents are substantially out of date.** The stabilizer audit (2026-04-30, against paclet 1.6.5) listed a set of "Tier 1 gaps" (Pauli-string measurement, Pauli expectation values, inner products, graph states, a test suite) that have **since all shipped** in the 2.0.0 stabilizer rewrite, and its Tier-3 performance items (packed representation, native kernel) shipped in June 2026. The Qiskit/QuTiP gap analysis (2026-02) predates the 2.0.0 transpilation, IBM-job, KAK-decomposition, and parametrized-layer surfaces. Section 11 lists the specific corrections.

---

## 1. Scope and method

### 1.1 Packages compared, with versions

| Package | Version | Role | Runnable here? |
|---|---|---|---|
| **QuantumFramework** (Wolfram) | 2.0.0 @ `02cc92d9` | Symbolic object-algebra quantum framework | Yes (wolframscript, Mathematica 15.0.0) |
| **Qiskit** (IBM) | 2.4.1 | Gate-model SDK + hardware | Yes (venv) |
| **Qiskit Aer** | 0.17.2 | High-performance simulator backend | Yes (venv) |
| **Cirq** (Google) | cirq-core 1.6.1 | Gate-model SDK + simulator | Yes (venv) |
| **QuTiP** | 5.3.0 | Open quantum systems / dynamics | Yes (venv) |
| **PennyLane** (Xanadu) | 0.45.0 | Differentiable QML | Yes (venv), not benchmarked |
| **Stim** (Google) | 1.16.0 | Stabilizer / QEC sampler (C++) | Yes (venv) |
| **QuantumClifford.jl** | 0.11.4 | Stabilizer simulator (Julia) | Yes (Julia 1.12.6) |

Wolfram-ecosystem packages (DiraqQ, SNEG, MulAtoLEG, etc.) and non-WL reference tools (QuantumOptics.jl, Strawberry Fields, Bosonic Qiskit, Perceval, CUDA-Q) are surveyed in Section 8 at the **[DOC]** level only.

### 1.2 What "verified" means for QuantumFramework

The separately-installed QuantumFramework paclet lags the repo. All QF claims and benchmarks here load the working tree explicitly:

```wolfram
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
```

This matters: a bare `Needs` of the installed paclet would, for example, lack the `"Path"` contraction option and report a different feature set. Sub-package surfaces (`SecondQuantization`, `QuantumOptimization`) require their own `Needs` (Section 4.4).

### 1.3 Benchmark protocol

- **Machine:** Apple M2 Pro, 12 cores, macOS (Darwin 23.2.0, arm64). Single machine, all runs.
- **Timing:** QuantumFramework via `AbsoluteTiming` with `ClearSystemCache[]` before each run, minimum of 2 to 3 repetitions (the `HoldFirst` discipline matters: timing an already-evaluated WL expression measures nothing). Python via `time.perf_counter`, minimum of 5. Julia via `BenchmarkTools.@belapsed` (minimum). All numbers are **minimums**, so they understate variance but are the standard "best case" comparison.
- **Threading:** all runs single-process. Stim, Aer, and QuantumClifford.jl can exploit SIMD and (Aer) multithreading; QuantumFramework's kernels are single-threaded WL. This is part of the honest picture, not a controlled variable.
- **Fairness:** each benchmark uses an **identical task** across frameworks. The stabilizer benchmark replays the exact same deterministic gate stream (shared op list); the state-vector benchmark builds the identical gate structure; the Lindblad benchmark integrates the identical physical model.

---

## 2. The contenders, neutrally

- **QuantumFramework.** A quantum object lives as a first-class Wolfram-Language expression (`QuantumState`, `QuantumOperator`, `QuantumChannel`, `QuantumCircuitOperator`, `QuantumMeasurementOperator`, `PauliStabilizer`). Everything is exact and symbolic by default, supports qudits ($d>2$), and flows into `Plot`, `Reduce`, `Dataset`, `DSolve`, etc. Strength: expressiveness and integration. Weakness: raw numeric throughput.
- **Qiskit + Aer.** Industry-standard gate-model SDK with the strongest transpiler (routing, layout, optimization levels 0 to 3), a full noise/error-mitigation suite, and a high-performance C++ simulator. Strength: hardware path and performance. Weakness: not symbolic, qubit-only.
- **Cirq.** Google's gate-model SDK and simulator. Similar niche to Qiskit, lighter transpiler, strong for NISQ experiments and Google hardware.
- **QuTiP.** The reference for open quantum systems: Lindblad, Monte-Carlo trajectories, Bloch-Redfield, Floquet, steady-state, and optimal control (GRAPE/CRAB). Strength: dynamics and control. Not a gate-model/QEC tool.
- **PennyLane.** Differentiable programming: autodiff through circuits, PyTorch/JAX/TensorFlow integration, quantum-kernel and QNN primitives. Strength: QML and gradients. 
- **Stim.** A specialized stabilizer/QEC sampler in C++; its headline is reference-frame bulk sampling at "constant cost per gate," billions of Pauli operations per second. Strength: QEC at scale. Not a general simulator.
- **QuantumClifford.jl.** Pure-Julia stabilizer simulator with four state types (`Stabilizer`/`Destabilizer`/`MixedStabilizer`/`MixedDestabilizer`), bit-packed kernels, Pauli-frame trajectories, and a large QEC code library. Strength: stabilizer breadth in a high-level language.

---

## 3. Benchmarks **[BENCH]**

All numbers measured on the machine and with the protocol in Section 1.3. Raw outputs are in Section 10.

### 3.1 State-vector simulation

**Task.** Build the circuit (apply $H$ to every qubit; a CNOT chain $0\to1\to\cdots\to n{-}1$; then $R_Z(0.3)$ on every qubit), repeated for 3 layers. For $n=12$ this is 105 gates. Prepare $|0\rangle^{\otimes n}$, run, and read out the full statevector (all frameworks).

> **Correction (2026-06-12).** The first edition of this table prepared the QF input as `QuantumState[Table[0,{n}],2]`, which the constructor parses as a length-$n$ state *vector* zero-padded into a 4-qubit, norm-0 state, not the $n$-qubit register $|0\rangle^{\otimes n}$; QF was therefore simulating the zero vector while Aer/Cirq simulated real states (see `QF-Apply-Path-Deep-Audit.md`, section 0.1, evidence script `scripts/apply00e_state.wls`). **[SRC]** The QF rows below are re-measured at repo SHA `02cc92d9` with the correct register `QuantumState["Register"[n]]` (input norm 1; output norm 1 with all $2^n$ amplitudes explicit) and a statevector readout, matching what Aer/Cirq read; the earlier density-matrix readout was free on the zero state but would be a dense $2^{16}\times2^{16}$ object on a real state. Script: `scripts/qf_bench4_register.wl`. **[BENCH]** The Aer/Cirq columns are unchanged: they were measured on real states.

| $n$ | QF, default (TensorNetwork) | QF, `Method->"Schrodinger"` | Qiskit Aer (statevector) | Cirq |
|---:|---:|---:|---:|---:|
| 12 | $800\,\mathrm{ms}$ | $14{,}298\,\mathrm{ms}$ | $2.1\,\mathrm{ms}$ | $3.7\,\mathrm{ms}$ |
| 16 | $3{,}378\,\mathrm{ms}$ | $86{,}438\,\mathrm{ms}$ | $8.6\,\mathrm{ms}$ | $12.1\,\mathrm{ms}$ |

**Reading it.**
- QuantumFramework's **default TensorNetwork engine is the one to quote** ($\approx 215\times$ to $395\times$ slower than Aer/Cirq), not the Schrödinger fold, which is another $18\times$ to $26\times$ slower still. Choice of `Method` is the single biggest QF performance lever; the default is the right one for this workload.
- The $n$-scaling is genuine: $800\to3378\,\mathrm{ms}$ ($4.2\times$) for a $16\times$ larger Hilbert space, comparable in shape to Aer's $2.1\to8.6\,\mathrm{ms}$ ($\approx4\times$). The apply-path deep audit decomposes where the time goes: at $n=12$ the cost is dominated by per-gate framework overhead (property-dispatch machinery and per-apply network re-assembly; the tensor contraction itself is $\sim16\%$), while at $n=16$ the contraction dominates ($\sim2.2\,\mathrm{s}$, $\sim65\%$ to $70\%$), slowed $15\times$ to $29\times$ by fully dense amplitudes living in `SparseArray` containers instead of packed arrays. Section 12.2 explains this in plain language, with the user-side levers; `QF-Apply-Path-Deep-Audit.md` sections 2 to 3.4 carry the full decomposition. (The first edition claimed near-flat $n$-scaling and "contraction is 1% of the time"; both were artifacts of contracting the zero state.)
- Practical takeaway: for textbook and few-qubit symbolic work QuantumFramework is fine; past $\sim18$ to $20$ qubits, or for tight inner loops, use a dedicated simulator (including, from inside QF, the qiskit bridge or `Method->"QuEST"`).

### 3.2 Stabilizer / Clifford simulation

**Task.** Replay an identical deterministic stream of single- and two-qubit Cliffords (a fixed $H$/$S$/CNOT pattern) on $n$ qubits, starting from $|0\rangle^{\otimes n}$. The exact same op list is fed to all three simulators (shared file for Stim, replicated deterministically in WL and Julia).

> **Revision (2026-06-12).** The QuantumFramework stabilizer subsystem was reworked for performance in three rounds during June 2026 (shipped through `730abcec`; Section 12.1 explains the work in plain language). The "QF original" column below is the pre-rework implementation that earlier editions of this table reported as the QF number; the "QF current" column is what ships at HEAD. Stim and QC.jl were re-measured fresh on the revision date. Total wall-clock ms including the constructor, min-of-3, `ClearSystemCache[]` per repetition.

| $n$ / gates | QF original | **QF current** (`ps["ApplyCircuit"]`) | Stim `TableauSimulator` | QuantumClifford.jl `MixedDestabilizer` |
|---|---:|---:|---:|---:|
| 100 / 2000 | $681\,\mathrm{ms}$ | $2.6\,\mathrm{ms}$ | $1.03\,\mathrm{ms}$ | $1.08\,\mathrm{ms}$ |
| 500 / 10000 | $15{,}836\,\mathrm{ms}$ | $26.2\,\mathrm{ms}$ | $5.14\,\mathrm{ms}$ | $23.1\,\mathrm{ms}$ |
| 1000 / 20000 | $63{,}469\,\mathrm{ms}$ | $58.0\,\mathrm{ms}$ | $10.9\,\mathrm{ms}$ | $108.4\,\mathrm{ms}$ |

(Measured at `6f953ad5`; re-verified at HEAD `02cc92d9` on 2026-06-12: $3.6$ / $25.5$ / $57.0\,\mathrm{ms}$, within run-to-run noise.)

**Reading it.**
- The cumulative speedup over the pre-rework implementation is **$266\times$ / $606\times$ / $1095\times$**. The gap to Stim is now **$2.5\times$ / $5.1\times$ / $5.3\times$** (it was $660\times$ to $5800\times$), and QuantumFramework **beats QuantumClifford.jl at $n=1000$** by $1.9\times$ while roughly matching it at $n=500$. Per-gate cost is $1.3$ to $2.9\,\mu\mathrm{s}$, in the same regime as the specialists.
- **The route matters.** The numbers above use the bulk entry point `ps["ApplyCircuit", specs]` (equivalently `PauliStabilizerApply`, which `qc[ps, Method -> "Stabilizer"]` reaches automatically when every gate is a recognized Clifford and the state has numeric signs). Applying gates one at a time through the property interface (`ps["H", 1]["CNOT", 1, 2]...`) still works and is itself $\sim66\times$ faster than the original ($951\,\mathrm{ms}$ at $n=1000$), but pays a fixed bookkeeping cost per gate that the bulk route does not. Building a `QuantumCircuitOperator` object out of $10^4$ individual gates is also slow in the host framework; feed `ApplyCircuit` the gate-spec list directly for long streams.
- **Correctness is not in question, and was actively strengthened by this work.** QuantumFramework ships cross-package test suites that validate its stabilizer engine against these exact tools: `Tests/Stabilizer/CrossPackage_Stim.wlt` (44 tests) and `Tests/Stabilizer/CrossPackage_QuantumClifford.wlt` (15 tests), part of a 965-test stabilizer suite (all passing at HEAD). The dense-simulation and Stim cross-checks built for the rework also surfaced and fixed two pre-existing measurement bugs that predate it (Section 12.1). **[SRC]**
- Where QF still genuinely lacks a feature: **bulk reference-frame (Pauli-frame) sampling**, the trick that makes Stim and QC.jl fast for QEC shot generation. QF's `StabilizerFrame` represents stabilizer mixtures and supports `"SampleOutcomes"`, but there is no constant-cost-per-gate frame sampler. For distance-$\geq7$ surface-code shot counts this is the operative gap.

### 3.3 Open-system (Lindblad) dynamics

**Task.** A driven, damped qubit: $H=\tfrac{\omega}{2}\sigma_x$ with $\omega=2\pi$, one collapse operator $\sqrt{\gamma}\,\sigma_-$ with $\gamma=0.3$, initial state $|1\rangle$ (excited), integrated over $t\in[0,2]$ (201 points). Compare wall time and the final excited-state population $\langle n\rangle(t{=}2)$.

| Quantity | QF `QuantumEvolve` | QuTiP `mesolve` |
|---|---:|---:|
| Wall time | $17.0\,\mathrm{ms}$ | $1.9\,\mathrm{ms}$ |
| Final $\langle n\rangle$ | $0.81862$ | $0.81901$ |

**Reading it.**
- QuantumFramework is only **$\approx 9\times$ slower** here, the closest of the three regimes, because both tools are doing the same thing: numerical integration of a (Liouvillian) ODE on a tiny system. QuTiP's edge is its specialized C-backed solvers; QF rides `NDSolve`. For small systems the constant factors are comparable.
- **The physics agrees:** $|0.81862 - 0.81901| = 3.9\times10^{-4}$, consistent with two different adaptive ODE integrators at default tolerances. This is the correctness cross-check for QF's open-system path. **[BENCH]** Re-verified at HEAD `02cc92d9` (2026-06-12): $17.6\,\mathrm{ms}$, final $\langle n\rangle = 0.818625$, unchanged. (The June Liouvillian fix `3225a899` made the four equivalent Lindblad call forms report their states in a common basis; it did not change the dynamics, which were always physically equal across routes.)
- Caveat: QuTiP's advantage grows with Hilbert-space dimension and with features QF lacks (Monte-Carlo trajectories, Bloch-Redfield, Floquet, steady-state, see Section 6). This benchmark is the small-system best case for QF; it is not a claim of parity at scale.

### 3.4 A correctness spot-check beyond timing: HHL / linear solve

`QuantumLinearSolve[{{2,1},{1,2}}, {1,0}]` returns $\{0.6667, -0.3333\}$, the exact solution $x=[\tfrac{2}{3}, -\tfrac{1}{3}]$ of $\begin{pmatrix}2&1\\1&2\end{pmatrix}x=\begin{pmatrix}1\\0\end{pmatrix}$. **[BENCH]** QuantumFramework's quantum linear-solve (HHL-family) path is numerically correct on this small instance. (It lives in the `QuantumOptimization` sub-package; see Section 4.4.)

---

## 4. Feature comparison matrix, refreshed against the 2.0.0 kernel

Legend: **QF** QuantumFramework 2.0.0, **QK** Qiskit, **QT** QuTiP, **PL** PennyLane. QF cells are **[SRC]** (verified live against the working tree). Competitor cells are **[DOC]** unless noted. "via bridge" means QuantumFramework reaches the capability by delegating to qiskit through its Python interop, not natively. Cells that **correct the prior stale document** are marked $\Delta$.

### 4.1 Core simulation

| Feature | QF | QK | QT | PL |
|---|:--:|:--:|:--:|:--:|
| State-vector simulation | Yes | Yes | Yes | Yes |
| Density-matrix simulation | Yes | Yes | Yes | Yes |
| Tensor-network simulation | Yes (default engine) | Yes | Yes | Yes |
| Stabilizer / Clifford | Yes (`PauliStabilizer`) | Yes | partial | No |
| MPS / MPO | Yes (`QuantumMPS`, Experimental) $\Delta$ | Yes | No | Yes |
| Qudit ($d>2$) native | **Yes** | No | partial | No |
| GPU acceleration | No | Yes (cuQuantum) | Yes (cuQuantum) | Yes (JAX) |
| Performance at $n\gtrsim18$ | weak (Section 3.1) | strong | n/a | strong |

### 4.2 Noise, mitigation, transpilation

| Feature | QF | QK | QT | PL |
|---|:--:|:--:|:--:|:--:|
| Named channels (depolarizing, amplitude/phase damping, bit/phase flip, reset, generalized amplitude damping) | Yes (8 named) | Yes | Yes | Yes |
| Thermal-relaxation $T_1/T_2$ named channel | No (compose from damping + reset) | Yes | Yes | No |
| Error mitigation (ZNE, PEC, measurement) | No | Yes | No | Yes |
| Dynamical decoupling / twirling | No | Yes | partial | partial |
| Transpiler: routing / layout / optimization | via bridge (qiskit) $\Delta$ | Yes (native, L0 to L3) | n/a | Yes |
| Hardware-aware target spec | Yes (`QiskitTarget`) $\Delta$ | Yes | n/a | Yes |
| Native gate-set transpile | via bridge (qiskit) $\Delta$ | Yes | n/a | Yes |

Note on transpilation: QuantumFramework 2.0.0 has a real, current path to a device basis and topology, but it **delegates to qiskit** through the Python bridge (`QuantumQASM[circuit, gateSet | QiskitTarget | options]`, `IBMJobSubmit`). It has no native SABRE-style router. So the honest cell is "via bridge," which is a capability upgrade over the prior document's flat "No" but is not native-parity with Qiskit.

### 4.3 Dynamics and open systems

| Feature | QF | QK | QT | PL |
|---|:--:|:--:|:--:|:--:|
| Lindblad solver | Yes (`QuantumEvolve`, Section 3.3) | Yes | Yes | No |
| Hamiltonian / unitary time evolution, symbolic | **Yes (closed form)** | No | No | No |
| Monte-Carlo trajectories | No | Yes | Yes | No |
| Steady-state solver | No | No | Yes | No |
| Bloch-Redfield / Floquet / non-Markovian | No | No | Yes | No |
| Adiabatic evolution | Yes (`QuantumAdiabaticEvolve`) $\Delta$ | partial | Yes | No |

### 4.4 Algorithms, variational, QML

| Feature | QF | QK | QT | PL |
|---|:--:|:--:|:--:|:--:|
| Parametrized circuit layers | Yes (`ParametrizedLayer`, `EntanglementLayer`) $\Delta$ | Yes | n/a | Yes |
| Parameter-shift gradients | Yes (`SPSRGradientValues`, `ASPSRGradientValues`) | Yes | No | Yes |
| Quantum natural gradient | Yes (`QuantumNaturalGradientDescent`, `FubiniStudyMetricTensor`) | partial | n/a | Yes |
| VQE building blocks | Yes (layers + gradients) | Yes | n/a | Yes |
| QAOA as a named routine | No | Yes | n/a | Yes |
| HHL / quantum linear solve | Yes (`QuantumLinearSolve`, Section 3.4) | partial | n/a | No |
| Autodiff through circuits | No (uses parameter-shift) | No | No | **Yes** |
| ML-framework bridge (PyTorch/JAX/TF) | No | partial | No | **Yes** |

**Important load-context caveat [SRC].** The variational/QML and linear-solve symbols (`ParametrizedLayer`, `EntanglementLayer`, `GenerateParameters`, `GradientDescent`, `QuantumNaturalGradientDescent`, `SPSRGradientValues`, `ASPSRGradientValues`, `FubiniStudyMetricTensor`, `QuantumLinearSolve`, `QuantumAdiabaticEvolve`) live in the **`Wolfram`QuantumFramework`QuantumOptimization`` sub-package**, and the Fock/bosonic symbols (`FockState`, `AnnihilationOperator`, `CoherentState`) in `Wolfram`QuantumFramework`SecondQuantization``. They are **not** exposed by a bare `Needs["Wolfram`QuantumFramework`"]` and resolve to `` Global` `` (inert) until you also load the sub-package:

```wolfram
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];   (* then ParametrizedLayer, QuantumLinearSolve, ... work *)
```

I verified each of the above is callable and correct once the sub-package is loaded. This is a real discoverability wrinkle, not a missing feature: the prior document's "QML/QNN layers: No" is now wrong, but the capability is one `Needs` away rather than on the main surface.

### 4.5 Hardware, visualization, interop

| Feature | QF | QK | QT | PL |
|---|:--:|:--:|:--:|:--:|
| IBM Quantum submission | Yes (`IBMJobSubmit`/`IBMJob`) $\Delta$ | Yes | No | Yes |
| OpenQASM 3 import/export, native (no Python) | **Yes** (`QuantumQASM`) $\Delta$ | partial | No | partial |
| AWS Braket | partial | No | No | Yes |
| Bloch sphere / circuit diagrams / Wigner | Yes | Yes | Yes | partial |
| ZX-calculus | Yes (`ZXTensorNetwork`, `ZXExpression`) | No | No | partial (via plugins) |
| Multiway / causal-graph circuit views | **Yes** (unique) | No | No | No |

---

## 5. Where QuantumFramework is genuinely unique **[SRC]**

These are capabilities the numeric specialists structurally do not have, verified live against the 2.0.0 kernel:

- **Exact symbolic object algebra.** States, operators, channels, and measurements can carry free symbolic parameters and combine in closed form. A Hamiltonian with a symbolic coupling, its symbolic time-evolution operator $e^{-iHt}$, and a parametric circuit are all first-class and composable. Qiskit/Cirq/QuTiP are numeric throughout.
- **Native qudits.** $d>2$ systems are first-class (the tableau kernel even carries a `d` parameter). The others are qubit-centric.
- **One object model across all of quantum theory.** Pure states, density matrices, channels (CPTP), POVMs, continuous-variable/Fock states, and phase-space quasi-probabilities are the same kind of expression, not separate libraries.
- **Basis as a first-class, reconciled attribute.** Every object carries its basis; QF reconciles mismatched bases automatically (Pauli decomposition, mutually unbiased bases, custom bases). This underwrites the QASM measurement-basis rebasing and the KAK magic-basis machinery.
- **Second quantization** (`FockState`, `AnnihilationOperator`, `CoherentState`) and **phase space** (`QuantumWignerTransform`, Husimi) in the same framework as the gate model.
- **Multiway and causal-graph circuit views** (`QuantumCircuitMultiwayGraph`, `...CausalGraph`, `...TokenEventGraph`), a Wolfram-Physics-flavored perspective no other package offers.
- **Wolfram-Language integration.** A quantum object flows directly into `Plot`, `Reduce`, `DSolve`, `Dataset`, `NDSolve`, the visualization stack, and the rest of the language. This is the real moat, and it is not a feature any standalone library can match.

---

## 6. Where QuantumFramework lags (verified gaps, 2.0.0) **[SRC]**

Confirmed absent in the current kernel (no matching source, verified by search across `Kernel/`):

- **State-vector performance at scale.** Two orders of magnitude behind Aer/Cirq (Section 3.1). Single-threaded; no packed numeric fast path yet (prototyped and validated, Section 12.2, but not shipped); no GPU. The *stabilizer* workload no longer belongs on this list (Section 3.2).
- **Error mitigation.** No ZNE, PEC, measurement-error mitigation, dynamical decoupling, or twirling.
- **Native transpiler.** No native routing/layout/optimization passes; the device path is delegated to qiskit.
- **Open-system breadth.** No Monte-Carlo trajectory solver, no steady-state solver, no Bloch-Redfield/Floquet/non-Markovian machinery (all QuTiP strengths).
- **Optimal control.** No GRAPE/CRAB/Krotov, no pulse-level control.
- **Named algorithm routines.** No QAOA as a packaged routine (the variational primitives exist; the named algorithm does not).
- **Differentiable programming.** No autodiff through circuits and no PyTorch/JAX/TF bridge (QF uses parameter-shift, not backprop). This is PennyLane's core niche.
- **Thermal-relaxation $T_1/T_2$ as a single named channel** (must be composed from amplitude/phase damping + reset).
- **Bulk Pauli-frame sampling** for QEC shot generation (Section 3.2).

---

## 7. Stabilizer subsystem: QuantumFramework vs Stim vs QuantumClifford.jl, refreshed

The 2026-04-30 audit compared the old `PauliStabilizer.m` (494 lines, paclet 1.6.5) against Stim and QuantumClifford.jl and produced a prioritized upgrade list. The subsystem was then **rewritten** for 2.0.0 (a 17-file `Kernel/Stabilizer/` subdirectory, six new public symbols). The design comparison still holds; the gap list largely does not.

### 7.1 Design philosophy (still accurate) **[DOC]/[SRC]**

| Axis | Stim | QuantumClifford.jl | QuantumFramework |
|---|---|---|---|
| State object | one `TableauSimulator` over the **inverse** tableau | four types (`Stabilizer`/`Destabilizer`/`MixedStabilizer`/`MixedDestabilizer`) | one `PauliStabilizer` head, always destabilizer-augmented ($\equiv$ `MixedDestabilizer`, rank $=n$) |
| Tableau storage | packed `simd_bits`, 256-bit lanes | bit-packed `xzs` chunks | canonical rank-3 array `{2, n, 2n}` **plus, since June 2026, a packed machine-word fast representation** (62 tableau bits per word, generator-major) used automatically for concrete states |
| Bulk gate execution | native C++ loop | native Julia loop | **one compiled-to-C kernel** (`ps["ApplyCircuit"]`) running the whole encoded circuit at $0.22$ to $0.27\,\mu\mathrm{s}$/gate, flat in $n$ |
| Phase | one `bool` (Hermitian, $\pm1$ only) | `UInt8`, two bits ($\pm1,\pm i$) | integer signs $\in\{-1,+1\}$ (symbolic $\mathbb{F}_2$-polynomial signs also supported, on the canonical path) |
| Bulk sampling | `FrameSimulator` (reference-frame) | `PauliFrame` + multithreaded trajectories | `StabilizerFrame` (mixtures), `"SampleOutcomes"`; **no** constant-cost frame sampler |
| Random Clifford | uniform symplectic | uniform symplectic | **Bravyi-Maslov / Mallows** (the mathematically uniform-on-Cliffords sampler; ahead of both) |
| Lines of core | $\sim20$k C++ | $\sim8.4$k Julia | 20-file WL subdir incl. `Packed.m`/`Compiled.m` (was a 494-line monolith) |

The two structural conclusions from the audit remain correct: QuantumFramework is closest to QC.jl's `MixedDestabilizer` semantics, and it is **ahead of both** on random-Clifford sampling (Mallows distribution, not symplectic-uniform).

### 7.2 The "Tier 1 gaps" the audit flagged: now closed **[SRC]**

Verified present and callable in the 2.0.0 kernel (each was listed as missing in 2026-04):

| Audit Tier-1 gap (2026-04, v1.6.5) | Status in 2.0.0 | Evidence |
|---|---|---|
| Pauli-string measurement `ps["M","XZZXI"]` | **Implemented** | `Kernel/Stabilizer/PauliMeasure.m`; `ps["Measure","XXX"]` returns the deterministic outcome live |
| Pauli expectation $\langle P\rangle$ | **Implemented** | `ps["Expectation","ZZI"]` returns $1$ for GHZ |
| Inner product $\langle\psi|\phi\rangle$ | **Implemented** | `Kernel/Stabilizer/InnerProduct.m`; `ps["InnerProduct",ps] = 1` |
| Bulk sampling primitive | **Partial** | `StabilizerFrame` + `"SampleOutcomes"` exist; constant-cost frame sampler still absent |
| Graph-state extraction | **Implemented** | public `GraphState`, `LocalComplement` |
| Test suite (`Tests/PauliStabilizer.wlt` "does not exist") | **Implemented** | 965-test suite at HEAD incl. `CrossPackage_Stim.wlt` (44), `CrossPackage_QuantumClifford.wlt` (15), `AuditMatrix.wlt` (203) |
| Logicals / entropy / Clifford channels | **Implemented** | `CliffordChannel`, stabilizer `"Entropy"`, `StabilizerStateQ` |

So the prior audit's headline recommendations have been executed. The audit's Tier-3 performance items (a packed machine-word representation and a native-code kernel) **have now also shipped**, in the June 2026 rework (Section 12.1): the tableau bits live in packed 62-bit words, and a single compiled-to-C routine executes whole Clifford circuits. Section 3.2 shows the result: single-digit factors from Stim, ahead of QC.jl at $n=1000$. The remaining genuine throughput item is constant-cost frame sampling (next section).

### 7.3 What both references have that QF still does not

- **Reference-frame bulk sampling** at constant per-gate cost (Stim's and QC.jl's reason to exist).
- **Detector error models** (Stim) for downstream decoders.
- **QEC code-family generators and decoders** (QC.jl's `QECCore`, $\sim100$ files).
- ~~Bit-packed performance kernels~~ : **closed** by the June 2026 rework (Section 12.1). The residual $2.5$ to $5.3\times$ to Stim is SIMD lane width and circuit-encoding overhead, not a representation gap.

---

## 8. The Wolfram-ecosystem packages (survey) **[DOC]**

From the in-repo notebook `OngoingProjects/OtherMathematicaQuantumPackages/quantum-packages.nb`. These are mostly **symbolic / second-quantization / Liouville** tools, not gate-model frameworks, and they complement rather than compete with QuantumFramework:

| Package | Niche | Note |
|---|---|---|
| **DiracQ** | symbolic many-body commutator algebra | tested |
| **SNEG** | second-quantization operator algebra (non-commuting) | tested |
| **Quantum Algebra** (Cesar) | older symbolic package, some useful pieces | |
| **secondQuant** | fermionic second quantization | fermionic only |
| **MMA-Quantum-Computing** | Dirac bra-ket, quantum algebra, QHD Heisenberg EOM | |
| **Quantum** | quantum computing + non-commutative algebra | tested |
| **MulAtoLEG** | Liouville superoperators for multi-level atoms | good for symbolic equations; no field in Hilbert space |

Non-WL tools listed there "for reference": **QuantumOptics.jl**, **QuTiP** (used as the reference for QF's Wigner/Husimi numerics), **Strawberry Fields** (photonic), **Bosonic Qiskit** (bosonic modes on Qiskit), **Perceval** (photonic, Quandela), **CUDA-Q** (GPU). These mark directions (photonics, bosonic, GPU) where QuantumFramework is thin and an interop or reference implementation could be borrowed.

---

## 9. How to choose (decision guide)

- **Symbolic derivation, teaching, qudits, exact algebra, Wolfram integration** : QuantumFramework, clearly.
- **Largest gate-model simulations, real IBM hardware, transpilation, error mitigation** : Qiskit + Aer (QuantumFramework can reach IBM and qiskit transpilation through its bridge, but Qiskit is native).
- **NISQ experiments, Google hardware** : Cirq.
- **Open quantum systems at scale, trajectories, steady state, optimal control** : QuTiP (QuantumFramework is fine for small Lindblad problems and is uniquely able to keep them symbolic).
- **Quantum machine learning, autodiff, hybrid training** : PennyLane (QuantumFramework has parametrized layers and parameter-shift gradients, but no autodiff/ML-framework bridge).
- **QEC / stabilizer simulation** : since June 2026 QuantumFramework's `PauliStabilizer` is a real option: cross-validated against Stim and QC.jl, within $2.5$ to $5.3\times$ of Stim and ahead of QC.jl at $n=1000$ on raw Clifford streams, with symbolic signs, qudit support, and the Wolfram environment around it. Choose Stim when you need **bulk shot sampling** (Pauli frames) or **detector error models** for decoders; choose QC.jl for its QEC code library. **[BENCH]**

---

## 10. Reproducibility

### 10.1 Environment

| Component | Version |
|---|---|
| Machine | Apple M2 Pro, 12 cores, macOS (Darwin 23.2.0, arm64) |
| Mathematica / WL | 15.0.0 (Feb 2026) |
| QuantumFramework | 2.0.0, repo working tree @ HEAD `02cc92d9` (via `PacletDirectoryLoad`); feature matrix verified @ `77b9ccba`, Sections 3.1/3.2 measured @ `02cc92d9`/`6f953ad5` |
| Python | 3.13.1 (venv `/tmp/qcbench`) |
| Qiskit / Aer / IBM-runtime | 2.4.1 / 0.17.2 / 0.47.0 |
| Cirq | cirq-core 1.6.1 |
| QuTiP | 5.3.0 |
| PennyLane | 0.45.0 |
| Stim | 1.16.0 |
| Julia | 1.12.6 |
| QuantumClifford.jl | 0.11.4 |

### 10.2 Raw benchmark output

```
# State vector (105-gate circuit at n=12), milliseconds
# QF rows re-measured 2026-06-12 at SHA 02cc92d9 with the correct register
# input QuantumState["Register"[n]] and statevector readout (Section 3.1)
QF default (TensorNetwork)   n=12: 799.706     n=16: 3378.029
QF Method->"Schrodinger"     n=12: 14298.303   n=16: 86437.715
Qiskit Aer (statevector)     n=12: 2.137       n=16: 8.575
Cirq                         n=12: 3.698       n=16: 12.061

# Stabilizer (identical deterministic Clifford stream), milliseconds
# QF current = ps["ApplyCircuit"] at 6f953ad5 (2026-06-12); externals re-measured same day;
# "QF original" = pre-rework implementation (the number earlier editions reported)
QF original                  n=100: 681       n=500: 15836     n=1000: 63469
QF current (ApplyCircuit)    n=100: 2.56      n=500: 26.2      n=1000: 58.0
  re-verified @ 02cc92d9     n=100: 3.58      n=500: 25.5      n=1000: 57.0
Stim TableauSimulator        n=100: 1.03      n=500: 5.14      n=1000: 10.9
QuantumClifford.jl           n=100: 1.08      n=500: 23.1      n=1000: 108.4

# Stabilizer hot paths, before -> after the June 2026 rework (ms; Section 12.1)
PauliStabilizer[n] ctor      n=1000: 355 -> 4.9
ps["M", 1] (single qubit)    n=500:  16.5 -> 0.19  (re-verified 0.187 @ 02cc92d9)
ps["M", "Z..Z"] (string)     n=500:  42300 -> 58.7 (re-verified 50.3 @ 02cc92d9)
ps["Stabilizers"] display    n=100:  128.6 -> 2.0  (re-verified 2.2 @ 02cc92d9)
composition a[b]             n=100:  2577 -> 46.1  (re-verified 45.4 @ 02cc92d9)

# Statevector apply-path levers, re-verified @ 02cc92d9 working tree (Section 12.2)
qc[psi] cold (fresh kernel)          2956 ms
qc[psi] warm                          836-839 ms
warm + $QuantumFrameworkPropCache=False  232-237 ms  (max amplitude deviation 0.)

# Lindblad (driven damped qubit)
QF QuantumEvolve   time: 16.987 ms   final <n>(t=2): 0.818625
QuTiP mesolve      time:  1.876 ms   final <n>(t=2): 0.819010
|delta <n>| = 3.85e-4

# Linear solve correctness
QuantumLinearSolve[{{2,1},{1,2}},{1,0}] = {0.6667, -0.3333}  (exact {2/3, -1/3})
```

### 10.3 Notes on fair comparison

- QuantumFramework numbers come from the **repo working tree**, not the installed paclet (which lags and reports a different feature set). Timing uses `HoldFirst` + `ClearSystemCache[]`; without `HoldFirst` the computation runs during argument evaluation and the timer measures nothing (an easy and silent mistake).
- Stim, Aer, and QuantumClifford.jl exploit SIMD/threads; QuantumFramework is single-threaded WL. The gaps in Section 3 are "as shipped, single machine," not a controlled-thread comparison.
- The stabilizer task is byte-identical across the three engines (shared op list, replicated deterministically). The state-vector task is identical in gate structure and, after the register-state correction of Section 3.1, in input state and readout. The Lindblad task is the identical physical model, and the two engines' results agree to $4\times10^{-4}$.
- Benchmarks are **minimums** over repetitions, so they favor every framework's best case equally and understate variance.

---

## 11. Corrections to the prior comparison documents

This master file supersedes three sources. The substantive corrections:

**A. `Notebooks/Function report/QuantumFramework_Comparison_and_Roadmap.md` (2026-02, pre-2.0):**
- "Transpiler: routing/mapping: No", "Topology-aware: No", "Native gate-set conversion: No" : now **via the qiskit bridge** (`QuantumQASM` with a gate set / `QiskitTarget`), and `QiskitTarget` is a real device-spec carrier. Still not native routing.
- "QML/QNN layers: No" : **wrong for 2.0**. `ParametrizedLayer`, `EntanglementLayer`, parameter-shift gradients (`SPSRGradientValues`), and `QuantumNaturalGradientDescent` exist (in the `QuantumOptimization` sub-package).
- "MPS/MPO: Basic" : `QuantumMPS` exists (Experimental sub-package) and converts to circuits.
- IBM Quantum: the path is now `IBMJobSubmit`/`IBMJob` (async SamplerV2/EstimatorV2), not the older interface the document assumed.
- Adiabatic evolution: `QuantumAdiabaticEvolve` now exists.
- Still accurate: no GPU, no ZNE/PEC/measurement mitigation, no GRAPE/CRAB, no QAOA routine, no autodiff/ML-bridge, no Monte-Carlo trajectories, no steady-state solver, no thermal-relaxation named channel. Its phased roadmap remains a reasonable wish-list for these.

**B. `OngoingProjects/Stabilizer/external-packages-audit.md` (2026-04-30, against v1.6.5):**
- Its "Tier 1 functional gaps blocking real use" (Pauli-string measurement, expectation values, inner product, graph states, a test suite) have **shipped** in the 2.0.0 stabilizer rewrite (Section 7.2). The design analysis (Section 7.1) and the performance gap (Section 3.2, Section 7.3) remain accurate.
- Its line-number citations (`PauliStabilizer.m:NNN`) are against the deleted 494-line monolith; the code now lives in `Kernel/Stabilizer/*.m`.

**C. `OngoingProjects/OtherMathematicaQuantumPackages/quantum-packages.nb`:** folded into Section 8 unchanged; it was a survey, not a claim set, so nothing needed correcting.

---

## 12. The June 2026 performance and correctness work, explained

This section is the plain-language account of two efforts that changed the numbers in this document, written for a physicist rather than a software engineer. Each claim links to a runnable check; Section 12.3 collects the commands. Detailed technical reports sit next to this file: `QF-Stabilizer-Optimization-Report.md` (stabilizer, shipped) and `QF-Apply-Path-Deep-Audit.md` (state-vector path, diagnosed and prototyped, **kernel deliberately untouched**).

### 12.1 Stabilizer engine: from $1000\times$ slower than Stim to single-digit factors (shipped)

**The physics.** A stabilizer state on $n$ qubits is specified not by its $2^n$ amplitudes but by the abelian group of Pauli strings that fix it. Concretely: $2n$ generator rows (stabilizers plus the conjugate "destabilizers"), each row a string of $n$ Paulis encoded as $2n$ bits (an X bit and a Z bit per qubit) plus a $\pm1$ sign. This is the binary symplectic tableau. A Clifford gate acts by conjugation, $P \mapsto U P U^\dagger$, which rewrites only the one or two qubit columns the gate touches: cost $O(n)$ per gate, independent of $2^n$. Measurement is the Aaronson-Gottesman update. None of this changed; the representation and the algorithms were already right.

**Why it was slow anyway.** The tableau was stored as an ordinary Wolfram array of integers, and every gate update walked it element by element through the general Wolfram evaluator. That evaluator is what makes everything in the language symbolic and composable, but it charges a fixed overhead of roughly a microsecond per elementary step, and a single gate update was taking dozens of such steps plus the cost of unwrapping and rewrapping the `PauliStabilizer` object. Stim pays none of this: it keeps 256 tableau bits in each processor register and updates them with single machine instructions.

**What was done (three rounds, all shipped).**
1. *Pack the bits into machine words.* Each qubit's column of X (or Z) bits across all $2n$ generators is stored as $\lceil 2n/62 \rceil$ machine integers, so one gate update is a handful of word-level XOR/AND operations instead of $2n$ array writes. This fixed the per-gate scaling.
2. *Run the whole circuit in one native loop.* The remaining cost was the per-step evaluator overhead itself. The entire gate loop was translated once into machine code (via WL `Compile` to C, built lazily and reused), so `ps["ApplyCircuit", specs]` hands the packed tableau and the whole encoded circuit to a single native routine: $0.22$ to $0.27\,\mu\mathrm{s}$ per gate, flat in $n$, which is Stim's own regime. The same round gave the register constructor its closed form (destabilizers $X_i$, stabilizers $Z_i$, signs $+1$, assembled directly), moved single-qubit measurement onto the packed words, and vectorized the Pauli-letter display (a letter is just $x + 2z$ indexed into $\{I, X, Z, Y\}$).
3. *Same packed treatment for Pauli-string measurement.* Measuring an $n$-qubit string like $Z_1 Z_2 \cdots Z_n$ needs an anticommutation scan of all $2n$ rows against the string and then row-by-row group updates; both now run on the packed words (the scan is a bit-count parity per 62-qubit chunk).

**Before and after, on identical inputs** (Apple M2 Pro; "before" is the pre-June-2026 kernel, "after" is HEAD; every row re-verified at `02cc92d9` on 2026-06-12):

| Operation | Before | After | Factor |
|---|---:|---:|---:|
| 20000-gate Clifford stream, $n=1000$ (Section 3.2) | $63.5\,\mathrm{s}$ | $58\,\mathrm{ms}$ | $1095\times$ |
| `PauliStabilizer[1000]` (the $\lvert 0\rangle^{\otimes 1000}$ register) | $355\,\mathrm{ms}$ | $4.9\,\mathrm{ms}$ | $\sim70\times$ |
| `ps["M", 1]` single-qubit measurement, $n=500$ | $16.5\,\mathrm{ms}$ | $0.19\,\mathrm{ms}$ | $\sim85\times$ |
| `ps["M", "Z..Z"]` string measurement, $n=500$ | $42.3\,\mathrm{s}$ | $59\,\mathrm{ms}$ | $720\times$ |
| `ps["Stabilizers"]` (display the group), $n=100$ | $129\,\mathrm{ms}$ | $2\,\mathrm{ms}$ | $63\times$ |
| Composition `a[b]`, $n=100$ | $2.58\,\mathrm{s}$ | $46\,\mathrm{ms}$ | $55\times$ |

A minimal before/after you can paste into a notebook (the "before" timing is what the same input cost on the pre-rework kernel; on the current kernel you will see only the "after"):

```wolfram
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
n = 1000;
specs = Catenate @ Table[
    {"H" -> q, "CNOT" -> {q, Mod[q, n] + 1}, "S" -> q}, {q, n}];  (* 3000 Clifford gates *)
AbsoluteTiming[PauliStabilizer[n]["ApplyCircuit", specs];]
(* first call ~0.5 s (one-time native compilation of the gate kernel),
   every repeat ~44 ms; the pre-rework kernel took ~10 s for this input *)
```

**Correctness work the speed work forced (and the bugs it caught).** Every packed path is gated on *exact* equality with the canonical slow path, not on statistics, and the cross-checks were run against dense state-vector simulation and against Stim. That discipline caught four real bugs, two of them **pre-existing** (present before any of this work):

- *Measurement update omitted the destabilizers* (pre-existing). Aaronson-Gottesman requires the row-sum update for **every** row that anticommutes with the measured operator, destabilizers included; the kernel updated only the stabilizer rows above the pivot. Single measurements, gate tests, and textbook states all pass under this bug; what breaks is the *correlation structure of later measurements* on some branches. Found by enumerating sequential full-register measurements of random Clifford circuits against dense simulation (3 of 15 circuits disagreed; Stim arbitrated that dense was right). Fixed in `1c524432`, with 15 permanent regression tests.
- *Pauli-string measurement left a non-abelian "stabilizer" group* (pre-existing, same disease in the string branch). Fixed in `52d3ea89`, cross-validated against dense simulation and Stim.
- *`PauliStabilizer["Random"[n]]` died near $n \approx 500$.* Two separate integer blow-ups: a structured-array wrapper pushed a matrix product off the fast packed route onto a generic path with tens-of-GB intermediates, and a matrix inverse was computed over $\mathbb{Z}$ (where unit-triangular inverses grow exponentially) instead of over $\mathbb{F}_2$, where the tableau actually lives. Fixed in `52d3ea89`; a residual machine-underflow warning in the Mallows sampler for $n \gtrsim 538$ was fixed at `02cc92d9` (exact rational arithmetic in the acceptance ratio; seeded output bit-identical).
- *`PauliStabilizer["GHZ"[5]]` returned unevaluated* (the constructor shortcut never matched the v2.0 `"Name"[args]` call form, and bare `"Random"` could return a corrupt object depending on definition load order). Fixed in `6f953ad5` with nine regression tests.

Verification status at HEAD: stabilizer suite **965/965**, full framework suite **1493/1493**, cross-package suites vs Stim (44) and QC.jl (15) all green, plus exact packed-vs-canonical equality pins at $n \in \{5, 20, 40\}$ and budgeted performance pins.

**Caveats that remain.** `ApplyCircuit` returns the canonical (unpacked) object while single gates return packed ones; both are accepted everywhere, but compare states via `["Tableau"]`/`["Stabilizers"]`, not `SameQ` on raw objects. Non-Clifford gates ($T$, $P[\theta]$) still route to `StabilizerFrame` exactly as before (term count doubles per non-Clifford gate). Composition and circuit synthesis remain $O(n^3)$ algorithms (correct and now fast enough at these sizes).

### 12.2 State-vector apply path: diagnosis complete, fixes validated, kernel untouched

The first-pass profile in earlier editions of this section was wrong, and understanding *why* it was wrong is half the story. Nothing in this subsection changed the kernel; it is a measurement campaign with session-local prototypes (`QF-Apply-Path-Deep-Audit.md` has the full decomposition; scripts in `scripts/apply0*.wls`).

**Correction 1: the old benchmark was simulating the zero vector.** The input state had been built as `QuantumState[Table[0, {12}], 2]`, on the assumption that this means "12 qubits in $\lvert 0\rangle$". It does not: the constructor reads a 12-component list as a *state vector* (amplitudes), zero-pads it to the nearest qubit dimension, and returns a 4-qubit state of norm $0$. Every "circuit application" was therefore propagating an empty array, amplitude arithmetic was free, and two first-pass conclusions ("contraction is 1% of the cost", "cost is flat in qubit number") were artifacts. See for yourself:

```wolfram
psiWrong = QuantumState[Table[0, {12}], 2];
{psiWrong["Qudits"], psiWrong["Norm"]}      (* {4, 0}  : 4 qubits, norm 0 *)
psiRight = QuantumState["Register"[12]];
{psiRight["Qudits"], psiRight["Norm"]}      (* {12, 1} : the real |0...0> *)
```

The corrected Section 3.1 numbers use the real register. The practical rule: **always build registers as `QuantumState["Register"[n]]`** (a warning guard for the misparse now exists in the working tree, uncommitted).

**Correction 2: `H["Normal"]` never measured a normalization.** `"Normal"` is not an operator property; the call fails property lookup and falls through to "apply the operator as a circuit", returning a one-gate circuit after several milliseconds of message handling. The first-pass claim that normalization was half the apply cost timed this accident.

**Where the time actually goes.** With the real register, applying the 105-gate benchmark circuit at $n=12$ costs $\sim 2.9\,\mathrm{s}$ the first time and $\sim 830\,\mathrm{ms}$ on every repeat. The repeat cost splits, by measurement and ablation, into three mechanisms:

$$
830\,\mathrm{ms} \;\approx\; \underbrace{580\,\mathrm{ms}}_{\text{memo bookkeeping}} \;+\; \underbrace{130\,\mathrm{ms}}_{\text{contraction in sparse containers}} \;+\; \underbrace{95\,\mathrm{ms}}_{\text{per-apply re-assembly}} .
$$

- *Memo bookkeeping (the dominant cost, and the surprise).* QF objects answer questions ("your matrix", "your qubit order") through a property system, and the answers are remembered so they are computed once. The remembered values are doing their job: with no memory at all the same apply costs $2.05\,\mathrm{s}$. But the bookkeeping that decides "should I remember this?" runs on **every** access, including the $\sim 66{,}000$ accesses per apply whose answer is already stored, and each such check costs $\sim 0.4\,\mathrm{ms}$ where the stored answer itself costs $1.6\,\mu\mathrm{s}$ to fetch. Disabling just the bookkeeping after a warm-up run (the stored values remain) gives the same state, byte-identical, $3.7\times$ faster. This also explains a superlinear growth with circuit depth: the bookkeeping keys grow with circuit size.
- *Dense data in sparse containers.* By constructor contract every state and gate tensor lives in a `SparseArray`. After one layer of Hadamards the state is fully dense (all $2^n$ amplitudes nonzero), and contracting dense data held in sparse containers costs $15$ to $29\times$ the identical arithmetic on packed numeric arrays. At $n=12$ contraction is 16% of the apply; at $n=16$ it dominates ($\sim 70\%$).
- *Per-apply re-assembly.* The tensor network is rebuilt from scratch on every application of the same circuit, even the thousandth.

**What you can do today, with no kernel change** (each verified, output identical):

```wolfram
(* 1. Repeated numeric application of a fixed circuit: warm up once,
      then switch the memo bookkeeping off. 836 ms -> 232 ms here. *)
qc[psi];   (* warm-up populates the stored values *)
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
qc[psi]    (* 3.7x faster, byte-identical state *)

(* 2. Registers: always "Register"[n], never Table[0, {n}]. *)

(* 3. Parameter sweeps: rebuilding the circuit at each value is ~4x faster
      than the substitution idiom qcP[<|theta -> v|>][psi]. *)

(* 4. Few qubits, many applications: compile once, reuse.
      22x amortized at n = 8 (memory is 4^n, so small n only). *)
qop = qc["QuantumOperator"];
qop[psi]
```

**What a kernel fix would buy** (prototyped and validated against the default path to $5\times10^{-17}$, but **not shipped**, per the decision to stop kernel changes): extracting each gate's numeric matrix once and folding it over a packed state vector runs the $n=12$ benchmark in $4.5\,\mathrm{ms}$, i.e. $178\times$ faster and inside Aer's measured $2$ to $9\,\mathrm{ms}$ envelope on this machine; the same prototype is $21$ to $25\times$ at $n=16$. Together with repairing the memo bookkeeping (cheap, $3.7\times$, fixes the depth scaling too) this is the prioritized plan in `QF-Apply-Path-Deep-Audit.md` Section 7. The point of the audit: the state-vector gap is **framework bookkeeping, not physics or algorithm**, and it is closable.

One more trap the audit ran into, worth knowing even though it is now guarded in the working tree: applying a state whose qudit count does not match the wires you assign it to (`{psi -> Range[12]}` with a 4-qubit `psi`) used to fire a silent "broadcast" rule that builds the tensor power of the state needed to cover the mismatch, here $16^{12} \approx 2.8\times10^{14}$ amplitudes, killing the kernel. With correctly-shaped states the rule never fires.

### 12.3 Run it yourself

Everything below assumes the repo working tree (`PacletDirectoryLoad` as in Section 1.2) and the environment of Section 10.1. Expected numbers are M2 Pro; expect $\pm20\%$.

```bash
# 1. Stabilizer headline benchmark (Section 3.2 QF column + the hot-path table):
#    expects the shared op streams in /tmp
cp "OngoingProjects/Platform Comparison/scripts/"stab_ops_*.json /tmp/
wolframscript -file "OngoingProjects/Platform Comparison/scripts/qf_final_all.wl"
#    -> gates_n1000_AC_ms ~ 57-58, ctor_n1000_ms ~ 4.9, measPauli_n500_ms ~ 50-59

# 2. The same op streams through the competitors:
/tmp/qcbench/bin/python "OngoingProjects/Platform Comparison/scripts/bench_stim_full.py"
JULIA_DEPOT_PATH=/tmp/qc_depot julia "OngoingProjects/Platform Comparison/scripts/bench_qc_full.jl"

# 3. Statevector: correct-register benchmark (Section 3.1 QF rows)
wolframscript -file "OngoingProjects/Platform Comparison/scripts/qf_bench4_register.wl"

# 4. Statevector: the cache-off lever and ablations (Section 12.2)
wolframscript -file "OngoingProjects/Platform Comparison/scripts/apply04_packing_ablations.wls"
wolframscript -file "OngoingProjects/Platform Comparison/scripts/apply09_cacheoff_validate.wls"

# 5. The test suites that gate all of this (from the repo root)
wolframscript -file Tests/RunTests.wls            # full framework, 1493 tests
```

Quick interactive spot-checks (paste into a notebook after loading the working tree):

```wolfram
(* stabilizer speed: 3000 gates at n=1000; ~0.5 s on the first call
   (one-time kernel compilation), ~44 ms on every repeat *)
n = 1000;
specs = Catenate @ Table[{"H" -> q, "CNOT" -> {q, Mod[q, n] + 1}, "S" -> q}, {q, n}];
AbsoluteTiming[PauliStabilizer[n]["ApplyCircuit", specs];]

(* measurement returns an Association: outcome bit b (result (-1)^b) -> post-state *)
ps = PauliStabilizer["GHZ"[3]];   (* this call form is itself a June 2026 fix *)
Keys @ ps["M", 1]                 (* {0, 1}: random outcome, two branches *)
Keys @ ps["M", "ZZI"]             (* {0}: the ZZ correlation is deterministic, +1 *)
Keys @ ps["M", "ZZI"][0]["M", "IZZ"]   (* {0}: and so is IZZ on the post-state *)
ps["Expectation", "ZZI"]          (* 1, the closed-form route to the same fact *)

(* the register trap *)
QuantumState[Table[0, {12}], 2]["Norm"]    (* 0  : not a register! *)
QuantumState["Register"[12]]["Norm"]       (* 1 *)

(* the cache lever, end to end *)
qc = QuantumCircuitOperator[Catenate @ Table[Join[Table["H" -> q, {q, 12}],
      Table["CNOT" -> {q, q + 1}, {q, 11}], Table["RZ"[0.3] -> q, {q, 12}]], {3}]];
psi = QuantumState["Register"[12]];
AbsoluteTiming[out1 = qc[psi];]   (* cold ~3 s, warm ~0.84 s on repeat *)
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
AbsoluteTiming[out2 = qc[psi];]   (* ~0.23 s *)
Max @ Abs[out1["AmplitudesList"] - out2["AmplitudesList"]]   (* 0. *)
```

---

## Appendix: scripts

The exact benchmark scripts are persisted next to this document in `scripts/`:
- QuantumFramework statevector (current, Section 3.1): `scripts/qf_bench4_register.wl` (correct register input `QuantumState["Register"[n]]`, statevector readout, `HoldFirst` timing).
- QuantumFramework stabilizer (current, Sections 3.2 and 12.1): `scripts/qf_final_all.wl` (`ApplyCircuit` + per-gate routes, constructor, measurements, formatting, composition). Companions on the same op streams: `scripts/bench_stim_full.py`, `scripts/bench_qc_full.jl`. Stabilizer test sweep: `scripts/run_stab_tests.wl`.
- Apply-path audit (Section 12.2): `scripts/apply00_diag.wls` through `scripts/apply09_cacheoff_validate.wls`; the per-script map is in `QF-Apply-Path-Deep-Audit.md`, Appendix B.
- Superseded provenance: `scripts/qf_bench2.wl` (pre-rework stabilizer numbers and the misparsed statevector input), `scripts/qf_bench3.wl` (Lindblad section still current). Do not re-run them for statevector or stabilizer numbers.
- Python (Aer, Cirq, QuTiP, Stim): `scripts/bench_py.py`, `scripts/bench_stim_shared.py`.
- Julia (QuantumClifford.jl): `scripts/bench_qc.jl`, `scripts/bench_qc_full.jl` (`@belapsed`).
- Shared stabilizer op lists: `scripts/stab_ops_100.json`, `scripts/stab_ops_500.json`, `scripts/stab_ops_1000.json` (generated by `scripts/gen_ops.py`, replicated deterministically in WL and Julia).

Re-running requires the environment in Section 10.1: a Python venv with the listed packages, Julia with `QuantumClifford` + `BenchmarkTools`, and Mathematica 15 with the QuantumFramework working tree loaded via `PacletDirectoryLoad`.
