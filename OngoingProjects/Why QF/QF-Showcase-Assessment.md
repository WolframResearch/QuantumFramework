# Assessment of `QF-Showcase.md`

**Audit date:** 4 August 2026

**Source audited:** `QF-Showcase.md`, 1,051 lines, 169 Wolfram Language code blocks

**Source location:** `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Why QF/QF-Showcase.md`

**Live environment:** Wolfram Language 15.0.0; current QuantumFramework working tree 2.1.0; Arrays 1.3.2; TensorNetworks 1.0.11
**Audit type:** complete static reading, live code execution, independent mathematical checks, package-source inspection, reproducibility review, and primary-source checks for external claims

## Executive assessment

The showcase is ambitious, unusually broad, and built around a genuinely strong set of examples. Much of its quantum mechanics is correct. The Rabi solution, time-dependent Lindblad relaxation law, SIC-POVM probabilities, oscillator identities, spin-frame invariance, zero-field splitting spectrum, tomography reconstruction, Grover probability, cluster-state stabilizers, operator algebra, Werner-state entanglement boundary, CHSH winning probability, and local GHZ construction all survived either direct execution or an independent derivation.

The document nevertheless **does not presently qualify as a verified, reproducible showcase**. Its central idea—that exact symbolic and executable quantum objects can connect many parts of quantum theory—is defensible. Several stronger formulations are not. The most serious problems are:

1. Four source examples hang in the current package because multi-qubit string shorthands such as `"0000"`, `"01"`, and `"00"` regress in QuantumFramework 2.1.0. Two additional state-arithmetic/application examples hang, and a probability-plot path stalled separately.
2. Several central prose claims are mathematically false or materially overstated, especially the proposed “radicals” boundary for exact symbolic computation, the assertion that Schrödinger simulation is a literal dense matrix-vector product, automatic stabilizer fallback at the first non-Clifford gate, and the claim that Qiskit-target transpilation runs in Wolfram Language.
3. The OpenQASM example is accepted by QuantumFramework's permissive parser but is not standard OpenQASM 3 as written. The native emitter also produces constructs that do not conform to the published grammar.
4. The measurement sections conflate reversible unitary premeasurement with a completed projective measurement. The teleportation section removes a classical wire only by retaining a coherent quantum record, which is a change of resource rather than elimination of communication.
5. The benchmark, hardware run, package versions, and installation procedure are not captured with enough provenance to reproduce the document's results.
6. The manuscript repeatedly shifts from “this example produced an exact expression” to universal claims such as “every object,” “all of quantum theory,” or “the exact object that numerical simulation estimates.” Those universal statements do not survive scrutiny.

**Publication verdict: major revision required.** The document should not be published with “Every cell below was run” or “the exact object behind every numerical approximation” as its warranty until the current failures, versioning, external provenance, and overclaims are corrected. After those changes, the verified core would make a persuasive technical showcase.

## Severity and evidence standard

- **Critical:** invalidates a headline claim, interoperability claim, or the document's verification warranty.
- **Major:** source code currently fails, a scientific claim is false or substantially misleading, or the result cannot be reproduced as presented.
- **Moderate:** a local statement needs qualification, an example is numerically rather than exactly verified, or provenance is incomplete.
- **Minor:** wording, notation, structure, or maintainability issue that does not change the result.

“Verified” below means that the source expression executed in the stated environment and its relevant output was checked. “Corrected verification” means that the intended mathematics worked after replacing a broken shorthand without otherwise changing the calculation. “Independent check” means that the claimed result was derived or computed outside the showcased QF call. “Static only” means that the syntax or implementation was inspected but the operation was not safely executable in this audit, usually because it required credentials, Python/Qiskit, graphics, or installation side effects.

## Audit method and limits

The audit read the complete Markdown source, inventoried all 169 Wolfram code fences, inspected the adjacent `verify-showcase.wls`, loaded the live sibling source trees, and evaluated the document in order with bounded timeouts. Ninety-eight blocks completed in the main ordered run. Additional targeted runs covered examples hidden behind earlier failures, including QAOA, the Heisenberg pair, quantum error correction, unmeasurement, stabilizers, Werner states, CHSH, and local GHZ creation.

The audit also inspected the relevant package implementation where the behavior of a claim depended on internals. This was essential for the circuit simulator, stabilizer handling, OpenQASM parser/emitter, and Qiskit target. Algebraic claims were checked independently where a QF call stalled. Hardware submissions and package installation were not performed: they would require credentials or mutate the user's environment. The IBM device's present metadata was checked separately, but the historical run asserted in the manuscript could not be reconstructed from the document.

This is a current-tree audit, not a historical recreation of QuantumFramework 2.0.0. The Markdown says it was run with 2.0.0, but it was committed after the repository version had advanced, the companion verifier describes 2.0.2, and the live tree is 2.1.0. Because no immutable environment or commit is pinned, the historical execution warranty cannot be certified.

### Evidence manifest for this assessment

| Item | Audited value |
|---|---|
| QuantumFramework repository commit | `f77a664025e2b1e7a2c4cbc05cf8d399001efb36` |
| QF paclet | 2.1.0, loaded from the working source directory |
| Arrays source | 1.3.2 at commit `87c6332afe36e22bdf7172c0e4e8556580d52189` |
| TensorNetworks source | 1.0.11 at commit `437a2d0c58a2f014102b45b7f9681c2a7f5ac725` |
| IBMQuantumPlatform source | 0.0.5 in the QF repository |
| Kernel | Wolfram Language 15.0.0, launched with `-noinit` |
| Platform | MacBook Pro, Apple M2 Pro (12 cores), 32 GB RAM, arm64, macOS 14.2 |
| Main timeout policy | 45 seconds per source block; known graphics, external, setup, and dependency-blocked blocks skipped explicitly |
| Targeted timeout policy | Up to 60 seconds for isolated reproducers and corrected paths |
| Raw audit logs used | `/tmp/qf_showcase_eval_1_142.log` and `/tmp/qf_showcase_eval_146_161.log`; these are ephemeral audit aids, not a release artifact |
| Source-tree state | The Markdown source itself had no audit edits. The repository was already dirty in unrelated course files and `QF-Showcase.nb`; this assessment is a new file. |

The kernel load order was Arrays, TensorNetworks, IBMQuantumPlatform, and QuantumFramework via `PacletDirectoryLoad`, followed by `Needs` for QF, SecondQuantization, and TensorNetworks. The minimal current-tree reproducers are the source calls `QuantumState["0000"]`, `QuantumState["01"]`, `QuantumState["00"]`, the block-132 state equality, and the block-143 `RY` application. Each exceeded the bounded probe rather than returning a value or `Failure`. Explicit `"Register"[n, integer]` state specifications bypassed the first family.

This manifest makes the current assessment traceable, but it does not turn the manuscript into a locked release test. A publication bundle should retain the runner, machine-readable results, expected outputs, and dependency archives in the repository rather than `/tmp`.

## Priority findings

| ID | Severity | Location | Finding | Required correction |
|---|---|---|---|---|
| F1 | Critical | lines 1–3, 1,019–1,023, 1,031 | “QF returns the exact object that numerical simulation or finite measurement estimates” is a category error. Hardware probabilities, noise, and finite-shot parameters are not made exact by a symbolic representation; several examples themselves use machine numbers, truncation, sampling, or numerical solvers. | Restrict the thesis to exact inputs and symbolically tractable models. Explicitly distinguish exact algebra, controlled numerical approximation, sampling, inference, and experimental data. |
| F2 | Critical | lines 1,019–1,023 | The proposed exactness boundary at matrices of dimension four and solvability in radicals is false. `Root` is an exact algebraic representation, not a failure; structured matrices above dimension four can have elementary spectra; and symbolic time evolution depends on more than radical solvability of a characteristic polynomial. | Replace the entire boundary paragraph with a discussion of expression growth, algebraic numbers represented by `Root`, special functions, model structure, and numerical fallback. Do not identify Abel–Ruffini as QF's symbolic boundary. |
| F3 | Major | blocks 21, 31, 88, and 96 | `QuantumState["0000"]`, `QuantumState["01"]`, and `QuantumState["00"]` hang in the live 2.1.0 tree. These failures block QAOA, Heisenberg evolution, teleportation, and QEC as written. | Replace them with explicit register states such as `QuantumState["Register"[n, 0]]`, `QuantumState["Register"[2, 1]]`, and `QuantumState["Register"[2, 0]]`; add regression tests for multi-digit computational-basis strings. |
| F4 | Major | blocks 132 and 143 | A simple arithmetic-equality example and a direct `RY`-on-state application hang in the current tree. The overloaded-expression thesis is therefore not demonstrated reliably by the current source. | Isolate and fix the package regressions, pin a passing commit, and make the verifier fail fast with a timeout rather than hanging. |
| F5 | Critical | lines 919–921, blocks 159–161 | `"QiskitTarget"` transpilation does not run entirely in Wolfram Language. Current QF source routes it through Qiskit using the Python bridge. Round-tripped objects are also not unchanged in the general case. | Say explicitly that this target requires a compatible Python/Qiskit environment and that only a supported, concrete qubit-circuit subset round-trips. |
| F6 | Critical | blocks 156–158 | The text calls the shown string standard OpenQASM 3, but it uses missing operand commas, uppercase `Pi`, and zero-count modifiers. QF accepts this as a permissive extension. The native emitter produces the same nonstandard forms. | Conform both emitter and example to the OpenQASM grammar; distinguish permissive import extensions from standard syntax; test output with an independent OpenQASM parser. |
| F7 | Major | lines 591–617 | The manuscript conflates reversible coherent premeasurement with a completed measurement. A coherently retained record can be uncomputed; a read or decohered record is not generally deterministically reversible. The cited superconducting-qubit experiment concerns probabilistic weak-measurement reversal. | Rename the construction “coherent premeasurement/unmeasurement,” state its resource assumptions, and sharply separate it from measurement reversal after amplification or environmental decoherence. |
| F8 | Major | lines 498–535 | “Classical communication is not required” is misleading. The circuit retains a coherent quantum record and uses controlled operations, replacing a classical message with a quantum-control/communication resource. The later paragraph partly acknowledges this, but the opening and conclusion still overstate the result. | Say consistently that deferred measurement yields a coherent equivalent of teleportation before readout; it does not remove the causal feed-forward resource. |
| F9 | Major | lines 621–663 | Schrödinger simulation is described as a literal dense `2^n × 2^n` matrix-vector product, but QF applies circuit gates sequentially and may use tensor-network machinery. The stabilizer-frame path does not automatically accept every first non-Clifford gate. | Describe the implemented algorithms, their actual asymptotic costs, and the supported phase-gate/frame decompositions. Show the failure mode for an unsupported arbitrary rotation. |
| F10 | Major | blocks 119–125 and related prose | The source constructs 19,998 gates, not 20,000. A cold live run took about 1.55 s; a warm run took about 0.295 s. The document times the first call yet calls it “a fraction of a second.” “Per-gate cost is flat in n” is also not a defensible asymptotic statement. | Report exact gate count, warm-up policy, hardware/OS/kernel/QF commit, repetitions, distribution, and memory. Replace “flat in n” with the measured range and implementation complexity. |
| F11 | Major | block 31 and lines 207–224 | The displayed two-qubit entropy formula has a removable `0 Log[0]` singularity. Direct substitution at a product-state endpoint is `Indeterminate`; zero is obtained by a limit. | Use `Limit`, `PiecewiseExpand`, or an entropy routine that defines the continuous extension. State the endpoint convention. |
| F12 | Major | lines 392–398 | The basis convention is internally inconsistent: “rows are basis vectors” is paired with the column-basis formula `B^-1 A B`. | Choose a row- or column-vector convention and make the matrix, prose, and formula agree, including complex conjugation. |
| F13 | Major | blocks 74–77 | “Returns exactly” is false because the code explicitly calls `N` and `Chop`. It extracts the diagonal but never verifies that the off-diagonal elements vanish. | Use exact assumptions and `FullSimplify`, or say “numerically diagonalizes”; assert `DiagonalMatrixQ` or compare the full transformed matrix. |
| F14 | Major | lines 488–493 | The document infers a fidelity-deviation scaling law from one seeded tomography sample. The observed infidelity is about `1.5×10^-5`, but one realization cannot establish or disprove an asymptotic rate; that rate depends on estimator, state rank, and metric. | Separate estimator-element standard errors from fidelity/infidelity behavior and support any scaling claim with repeated trials across several shot counts or an applicable theorem. |
| F15 | Moderate | blocks 52–57 | The text says the truncated oscillator examples give `g^(2)(0)=0,1,2`; live results were `0`, `0.999999999968...`, and `1.999675737326...`. | Say “approaches” or “is numerically approximately” one and two at the chosen truncation, and report the cutoff. |
| F16 | Major | lines 23–24 | “SymPy and every Julia package cannot represent [this] at all” is an unbounded and unsupported competitor claim. The Mathematica competitor is unnamed and unversioned. | Name every compared package and version, provide a runnable benchmark, define “represent” and “solve,” and publish failures as artifacts. Otherwise remove the comparison. |
| F17 | Major | blocks 13–16 | The finite-time Landau–Zener section calls machine-precision evaluation “exact survival,” uses an unguarded immediate definition that stalls in the current run, and labels a single-pass plot “Stückelberg oscillations.” | Use `survival[t_?NumericQ] := ...`, distinguish the exact symbolic propagator from its numerical plot, and use “finite-time Landau–Zener oscillations” unless a multipassage interferometer is modeled. |
| F18 | Moderate | blocks 38–43 | The time-dependent Lindblad solution is described as valid for “every state,” but `Array[rho,{2,2}]` is an arbitrary complex matrix without Hermiticity, positivity, or trace-one constraints. A global `$Assumptions` assignment is never restored. | Say “generic matrix entries” and impose physical-state constraints when relevant. Replace the global assignment with `Assuming` or `Block`. |
| F19 | Major | blocks 162–168 | The hardware story has no immutable job ID, submission date, backend calibration, captured result object, or raw counts. The prose assertion that the job was queued and returned 4,096 shots is not auditable from the Markdown. | Record a sanitized job identifier, timestamps, backend metadata, transpiled circuit, status history, counts, and an exported result artifact. Keep live availability separate from historical provenance. |
| F20 | Critical | blocks 169 and closing note | The installation procedure is unpinned and order-dependent. It loads QF before installing the development TensorNetworks build, can mix published and development dependencies in one kernel, omits several dependency versions, and mutates the environment with forced installs. | Pin repository commits or paclet archives, install all dependencies before `Needs`, start a fresh kernel, capture `$MaxExtraPrecision`/system data where relevant, and provide a lock manifest plus a non-mutating verifier. |

## Claim-by-claim assessment

The source's claim numbers are retained below. For technical dependency, Claims 4 and 5 are discussed before Claim 3; this does not create an additional claim. The basis material is a subsection of source Claim 2, and source Claim 6 receives its own visualization verdict.

### Source Claim 1 — “Symbolic-first, exact algebra by default”

**Verdict: the examples demonstrate real symbolic capability, but the competitive and universal claims are overstated.**

The generic qubit-state representation, amplitude access, Bloch vector, Rabi Hamiltonian, commutator, propagator, and transition probability all executed and agreed with independent algebra. The Landau–Zener Hamiltonian and symbolic propagator also evaluated; the appearance of special functions is a legitimate strength.

The defects are in the framing and the numerical tail of the example:

- The manuscript cannot establish that an unnamed Mathematica package, SymPy, and “every Julia package” are categorically unable to represent the solution. That requires named versions, exact inputs, outputs, timeouts, and a definition of success.
- A library can represent Kummer or generalized hypergeometric functions without automatically deriving this particular solution. “Cannot solve this model automatically under this benchmark” would be defensible; “cannot represent it at all” is not.
- The survival-probability plot uses machine numbers and `N`, so the curve is a numerical evaluation of an exact symbolic expression, not itself an exact output.
- The source's immediate numeric definition stalled. A plot-facing function should be numeric-guarded and delayed.
- “Landau–Zener–Stückelberg” normally signals interference between multiple passages. The shown model is a finite-time single sweep; the wording should not imply a multipassage interferometer.
- The asymptotic Landau–Zener result is a limit in the appropriate preparation/scattering regime, not simply a “long-time average” without further qualification.

The Heisenberg exchange example is mathematically sound after the broken input state is expressed explicitly. At `t = π/(8J)`, the evolved state is maximally entangled and the one-qubit entropy is one bit. At product-state endpoints, however, the displayed logarithmic formula must be interpreted by continuity. Direct substitution gives `Indeterminate`, so the text should show `Limit` rather than say the formula itself directly returns zero.

The complete four-vertex graph (`K4`) QAOA formula was independently checked despite the source-state hang. It gives three when either angle is zero, is periodic in `θ` with period `π` and in `γ` with period `2π`, and reaches approximately `3.69751609925`. Its optimum value satisfies

```text
32 z^3 - 393 z^2 + 1242 z - 837 = 0,
```

so the cubic-algebraic optimum claim is supported. The source should nevertheless stop calling the currently hanging program a live verified example until the input regression and runtime are addressed.

### Source Claim 2 — “One object model spans all of quantum theory”

**Verdict: mostly supported, with important exactness and state-domain qualifications.**

The amplitude-damping channel reduces the purity of the tested input as stated. The time-dependent Lindblad relaxation matrix and `T2 <= 2 T1` reduction, SIC-POVM probabilities, measurement probabilities, oscillator operator identity, and coherent-state eigenvalue relation worked in substance. Point samples established one positive and one negative Wigner value and a small positive Husimi value; they did **not** establish global normalization, global Husimi nonnegativity, extrema, lobe heights, or grid/cutoff convergence.

The north-pole SIC-POVM probabilities are exactly `{1/2, 1/6, 1/6, 1/6}`. The coherent and thermal `g^(2)(0)` examples are finite-cutoff numerical states, so their live values are only approximately one and two. The Wigner value at the origin was approximately `0.318237`, close to `1/π ≈ 0.318310`; the manuscript already gestures at truncation and should apply that qualification consistently.

“For every state” is too strong for the time-dependent calculation as coded. `Array[rho,{2,2}]` creates unconstrained entries, not a physical density matrix domain. The result is a formal solution for generic matrix entries; applying it to physical states additionally assumes Hermiticity, positivity, and unit trace. Those invariants were not independently asserted in the showcase.

The section also sets `$Assumptions` globally and never restores it. In a long notebook this can silently influence unrelated later simplifications. The showcase should use locally scoped assumptions throughout.

#### Basis and dimension subsection of source Claim 2

**Verdict: the executable examples work, but the written basis convention must be repaired.**

The `X`, spin-`Jx`, computational, Hamiltonian-eigenvector, and Fourier-basis examples evaluated. The spin-one expectation value of `Jx` is one before and after the frame change. The zero-field-splitting Hamiltonian has exact eigenvalues `{0, 2 d - e, 2 d + e}` as stated.

The prose says the rows of `B` are basis vectors and then writes `A_B = B^-1 A B`, which is the standard column-basis change formula. With basis vectors stored as rows, the conjugation order changes. This is not a cosmetic issue: the section is explicitly teaching a translation convention. Pick one convention and state it once.

For the five-site ring, the code transforms numerically and extracts a diagonal. The values agree with the expected spectrum, but “returns exactly” is false because `N` is called explicitly. Moreover, taking `Diagonal[m]` does not prove that `m` is diagonal. The example needs a full-matrix assertion.

The opening phrase “every object can live in any finite dimension” also needs narrowing. QF's core state/operator/channel abstractions support broad finite-dimensional systems, but stabilizer methods, OpenQASM, hardware backends, and several named constructs are qubit-specific; continuous-variable examples are truncated. “The core linear-algebra object model supports arbitrary finite dimensions” is a safer claim.

### Source Claim 4 — “One object model, three computational engines”

**Verdict: the representation-switching idea is valuable; several implementation descriptions and the benchmark narrative are inaccurate.**

The Grover example produced the exact success probability and the numerical value `0.9999470421032737`. Tensor-network conversion and the most-probable marked state were supported by the companion verification. The cluster-state graph, stabilizer generators, selected expectation values, stabilizer-to-state-vector conversion, and structural equality check also worked.

The manuscript's description of Schrödinger simulation as “a literal `2^n` matrix-vector product” should be removed. A full dense operator of size `2^n × 2^n` would be much more expensive than sequential local-gate application. Current QF circuit code folds over gates and can use tensor-network application paths; it does not establish the claimed literal dense multiplication.

The stabilizer-frame implementation supports particular non-Clifford decompositions, notably phase-like gates. It is not a universal automatic fallback at the first arbitrary non-Clifford gate. An unsupported rotation can return a non-Clifford failure. The text should specify the supported gate family and the approximation/decomposition policy, if any.

The benchmark requires a rewrite:

- The constructed circuit has 19,998 gates, not 20,000.
- In the live environment, a cold first evaluation was about 1.55 seconds and a warm same-session evaluation about 0.295 seconds.
- The displayed `AbsoluteTiming` is a first call, so a warmed-cache number cannot silently support its caption.
- “Per-gate cost is flat in `n`” is neither proved by the single example nor consistent with the state-size dependence of tableau/frame data structures.
- The entropy value of 496 across a 500-qubit cut is striking and near maximal, but “Page-like” should be presented as an analogy. Random stabilizer states are not the Haar ensemble used in Page's theorem.
- The referenced companion filename `QF-Stabilizer-vs-Packages.md` was not found. A nearby optimization report exists, but its benchmark conditions and dated commit must be linked accurately.

The current official Aer documentation does list several simulation methods, so the comparison can be made, but it needs versioned citations and like-for-like tasks rather than a rhetorical list. Stim is an exceptionally optimized stabilizer simulator; comparisons should report interface overhead separately from compiled-kernel time. See the official [Qiskit Aer simulator reference](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html) and [Stim repository](https://github.com/quantumlib/Stim).

### Source Claim 5 — “A quantum object is an ordinary Wolfram Language expression”

**Verdict: conceptually and mathematically strong, but the current execution path has regressions.**

The exponential of a symbolic Pauli Hamiltonian, operator sum/decomposition, commutator identity, time-dependent Hamiltonian reconstruction, and second-order Schrödinger-equation check all executed successfully. The recovered Hamiltonian coefficients were

```text
<|X -> Cos[t^2]/2, Y -> Sin[t^2]/2, Z -> t, I -> 0|>,
```

and the equation residual was zero.

The Werner-state realignment criterion simplified to `(-1 + Abs[3 - 4 p])/2`, and the witnessed entangled region `0 < p < 1/2` is correct for the convention used in this document.

However, block 132—a simple structural equality in this narrative—and block 143—the application of a symbolic `RY` rotation to a state—hang in the current tree. The intended variational energy is independently `Cos[θ] + Sin[θ]`, with maximum `Sqrt[2]`, but the showcased path did not finish. This does not disprove finite operator overloading; it does defeat the claim that these cells are presently verified. The section's own finite-overload caveat is good and should be retained.

### Source Claim 3 — “Measurement records are quantum wires, and feedforward is a controlled gate”

**Verdict: the coherent circuits are useful; the interpretation overreaches.**

The teleportation circuit and input Bloch vector evaluated. After replacing the broken `"00"` shorthand, Bob's final Bloch vector exactly matched the input:

```text
{Cos[φ] Sin[θ], Sin[φ] Sin[θ], Cos[θ]}.
```

The three-qubit correction-code syndrome and recovery also worked after the same state-input correction, returning unit fidelity. The unitary “unmeasurement” circuit returned the original state exactly in the tested ideal model.

These are coherent dilations or deferred-measurement constructions. They are not proof that an already amplified, decohered, projective laboratory measurement can generally be deterministically undone. Once the record remains coherent, it is still part of the quantum system and can be uncomputed. Once it has irreversibly leaked to an environment, the conditions change. The review literature describes this distinction explicitly; see the unitary-premeasurement discussion in [Bassi et al., arXiv:1412.7862](https://arxiv.org/abs/1412.7862).

The cited superconducting-qubit experiment is about conditional, probabilistic reversal of a weak measurement, not deterministic reversal of a projective measurement. See [Katz et al., arXiv:0806.3547](https://arxiv.org/abs/0806.3547). Calling the showcase's unitary construction an “idealized deterministic limit” does not bridge that physical distinction.

Similarly, coherent teleportation has not abolished the feed-forward resource. Alice's record qubits coherently control Bob's corrections. If Alice and Bob are separated, those controls still require a causal quantum interaction or transmission. The correct claim is that a deferred-measurement circuit can represent the feed-forward coherently until terminal readout.

The stated dimension growth from `2^n` to `2^(n+k)` assumes qubits and one binary record qubit per outcome bit. For general qudits or multi-outcome measurements, the enlarged dimension is the product of the system and record dimensions. The QEC discussion is also specifically about coherent syndrome extraction, not fault-tolerant measurement-free error correction under realistic noise. The symbolic fidelity-one proof uses real `α,β` assumptions; the prose extension to arbitrary complex amplitudes was not separately demonstrated and needs a normalized phase-bearing test.

Deferred-measurement transforms are not unique to QF. Cirq and PennyLane both document program transforms of this kind; their exact supported semantics differ. See [Cirq `defer_measurements`](https://quantumai.google/reference/python/cirq/defer_measurements) and [PennyLane `defer_measurements`](https://docs.pennylane.ai/en/stable/code/api/pennylane.defer_measurements.html). A comparison should focus on symbolic integration and object interoperability rather than imply absence elsewhere.

### Source Claim 6 — “Every object draws itself”

**Verdict: the examples show a broad visualization surface, but this claim was only statically audited and is not fully verified.**

The showcase calls state, amplitude, measurement, stabilizer-tableau, circuit, tensor-network, and hypergraph views. Most of those blocks return graphics or formatted notebook objects; this audit did not perform image-level comparison. The report therefore does not certify the claimed arrow angles, interference-comb null, phase colors, histogram labels, tableau layout, graph topology, option forwarding, or publication quality. Several underlying mathematical objects were checked elsewhere, but a correct underlying state does not prove that its renderer depicts every promised feature correctly.

The universal wording “every object draws itself” is broader than the evidence. A bounded, testable formulation would list the object heads and properties supported in the pinned version. Verification should assert returned heads and underlying plotted data, render representative outputs, compare them visually, and test that documented graphics options are forwarded. The phase-estimation example must also use an explicit or seeded unitary; its bare random default makes the diagram and network non-reproducible.

The source's assertions about roughly three dozen named bases, Bell-derived KAK compilation, generalized POVM/Naimark/qudit visualization, tomography internals, and the complete plot-property inventory were not established by the worked cells. They should be linked to versioned API documentation or moved to a separately tested capability matrix.

### Source Claim 7 — “An interoperability hub, down to real hardware”

**Verdict: local symbolic-to-circuit examples are promising; the standard-compliance, Python-free, and hardware-provenance claims fail as written.**

The native QASM exporter executed. The manuscript's QASM string was also accepted by QF and the internal CHSH calculation returned the exact optimal winning probability `(2 + Sqrt[2])/4`. Bare `PhaseEstimation[3]` constructed a circuit, but it defaults to a random unitary and is not seeded, so its rendered output is not reproducible.

Acceptance by QF does not establish OpenQASM compliance. The shown source uses constructs such as space-separated gate operands, uppercase `Pi`, and `ctrl(1) @ negctrl(0)`. The published grammar specifies comma-separated operand lists; predefined mathematical constants include lowercase `pi`; and `ctrl`/`negctrl` counts are positive when supplied. QF's parser explicitly accepts whitespace as a permissive extension. The document must call this a QF dialect until the emitter is fixed. Primary references are the [OpenQASM grammar](https://openqasm.com/grammar/index.html), [type and constant specification](https://openqasm.com/versions/3.0/language/types.html), and [gate-modifier specification](https://openqasm.com/versions/3.1/language/gates.html).

The `QiskitTarget` path requires Qiskit through QF's Python bridge. In this audit environment Qiskit was absent, and the call did not complete. The statement that transpilation “runs in Wolfram Language” is directly contradicted by the package implementation. “Same symbolic objects leave Wolfram Language and re-enter unchanged” is also too broad: QASM targets qubit circuits, concrete angles, and a supported instruction subset; information can be lowered, canonicalized, or lost.

The IBM backend name remains current. At audit time, IBM's live systems page listed `ibm_fez` online with 156 programmable qubits, but that current fact does not verify the manuscript's historical job. See [IBM Quantum Platform system information](https://quantum.cloud.ibm.com/computers?order=eplg+DESC&search=fez&system=ibm_fez&view=table). A real hardware showcase needs the original result artifact and submission metadata.

MATLAB has supported IBM hardware workflows since R2023b according to its version history, so that comparison is broadly grounded. The manuscript's stronger “Python-free” characterization should be sourced precisely or omitted. See MathWorks' [IBM hardware workflow](https://www.mathworks.com/help/matlab/math/run-quantum-circuit-on-hardware-using-IBM.html) and [`run` documentation](https://www.mathworks.com/help/matlab/ref/quantumcircuit.run.html).

## Regimes and limits matrix

| Family | General symbolic case | Exact special case demonstrated | Limit/asymptotic case | Independent numeric reference | Failure or edge case |
|---|---|---|---|---|---|
| Rabi / Landau–Zener | Symbolic qubit Hamiltonians and propagators | Rabi transition probability; special-function LZ propagator | LZ scattering/asymptotic regime needs explicit preparation and limit | Finite-time survival plot uses machine evaluation | Unguarded plot function stalled; “Stückelberg” terminology overreaches a single sweep |
| Lindblad dynamics | Generic matrix entries and symbolic rates | Closed relaxation/dephasing law and `T2 <= 2 T1` | Long-time ground-state limit follows for physical normalized inputs and positive rates | No independent ODE integration was retained | Input matrix is not constrained to a density matrix; trace/positivity preservation not asserted; assumptions leak globally |
| Oscillator / phase space | Symbolic truncated coherent-state algebra | Coherent-state eigenvalue/evolution identities away from truncation edge | Infinite-cutoff values are approached numerically | Selected `g^(2)`, Wigner, and Husimi point values checked | Finite-cutoff error; no grid/cutoff convergence, quadrature normalization, or global positivity test |
| Tomography | State reconstruction object from seeded counts | One 2,000-shot-per-basis synthetic run | Large-shot scaling is not established by one sample | Reconstructed matrix and similarity evaluated | Boundary/rank/metric effects; trace, positivity, bias, and coverage should be tested across repetitions |
| Tensor/stabilizer engines | Tensor-network and tableau representations | Grover probability, cluster generators, representation equality | Near-maximal 1,000-qubit half-cut entropy in one seeded circuit | Cold and warm timings measured | 19,998 rather than 20,000 gates; unsupported non-Clifford rotations; `%` history dependence; no asymptotic benchmark sweep |
| Measurement dilation | Coherent system-record unitary | Teleportation, repetition-code syndrome, and uncomputation after state-input correction | Completed/decohered measurement lies outside the reversible model | Symbolic Bloch vector, syndrome, and fidelity checks | Quantum communication resource remains; complex-amplitude QEC extension and realistic noisy/fault-tolerant regime are open |

## Open invariants and API assertions

The following statements were read and assessed but remain **unverified**, not silently accepted:

- positivity, completeness, and equal pairwise overlaps of the returned SIC-POVM effects;
- trace preservation and physical-state preservation of the symbolic Lindblad map over a constrained density-matrix domain;
- Wigner/Husimi normalization, global Husimi nonnegativity, extrema, lobe heights, and convergence with cutoff/grid size;
- unitarity of every named basis, full off-diagonal cancellation in the Fourier ring example, and the claimed named-basis count;
- trace one, positive semidefiniteness, bias, and confidence behavior of the tomography estimate across repeated samples;
- semantic equivalence of QASM import/export under an independent standards-compliant parser and simulator;
- exact current signatures and documented support for every plot/property name, Naimark dilation, arbitrary-observable/degenerate-projector handling, KAK basis conversion, tensor-network graph, and hypergraph claim; and
- pixel-level correctness and option forwarding for the visualization examples.

Each should either gain a small executable assertion or be presented as a linked, versioned API capability rather than as an output demonstrated by the current cells.

## Complete code-block ledger

Every Wolfram code fence in the source is accounted for below. Block numbers follow source order.

| Blocks | Source lines | Status | Assessment |
|---|---:|---|---|
| 1–3 | 27–39 | Verified | Generic symbolic qubit state, amplitudes, and Bloch vector worked. |
| 4–5 | 48–56 | Verified | Rabi Hamiltonian and commutator worked. |
| 6–9 | 61–84 | Verified | Propagator, evolved state, and transition probability agree with the standard Rabi solution. |
| 10–13 | 89–106 | Verified | Landau–Zener Hamiltonian and symbolic propagator evaluated; special-function structure is genuine. |
| 14–16 | 110–131 | Stalled / numerical | Plot-facing evaluation did not finish under bounded execution. The code uses machine precision and `N`; it is not an exact plot. |
| 17–20 | 136–152 | Verified | `K4`, its maximum cut, register size, and diagonal cost Hamiltonian worked. |
| 21 | 156–158 | Current failure | Multi-digit string state preparation hangs in QF 2.1.0. |
| 22–29 | 160–198 | Dependency-blocked; independently checked | QAOA depends on block 21. The landscape formula, periods, maximum, and cubic optimum were checked independently. |
| 30 | 203–205 | Verified | Heisenberg Hamiltonian constructed. |
| 31 | 207–209 | Current failure | `QuantumState["01"]` hangs. Explicit register form works. |
| 32–34 | 211–224 | Corrected verification | Evolved state, entropy, and Bell time are correct after fixing block 31. Entropy endpoints need `Limit`. |
| 35–43 | 234–282 | Verified with caveats | Channel purity, time-dependent Lindblad relaxation, and `T2 <= 2 T1` worked. The initial matrix is unconstrained and assumptions leak globally; trace/positivity preservation was not separately asserted. |
| 44–50 | 287–316 | Verified with bounded scope | SIC probabilities and circuit construction worked, including the north-pole values. Full POVM positivity/completeness, pairwise SIC overlaps, and plot rendering were not independently checked. |
| 51 | 318–320 | Stalled | `circ[]["ProbabilitiesPlot"]` did not finish in the bounded live run. |
| 52–58 | 325–357 | Verified with truncation | Oscillator algebra worked. Coherent/thermal `g^(2)` values are approximate at finite cutoff; block 58 is the bar chart and was not visually certified. |
| 59–62 | 362–378 | Point values verified | Wigner/Husimi representations constructed and three point values were checked. Global normalization, extrema, positivity, and convergence were not established. |
| 63 | 380–388 | Static/visual | Side-by-side phase-space surface plot was not image-verified. |
| 64–73 | 393–439 | Verified with caveat | `X`, spin, computational, eigenvector, and Fourier-basis examples worked. The prose row/column convention is inconsistent; named-basis inventory and KAK/export claims remain static. |
| 74–77 | 443–461 | Numerically verified | Ring spectrum is correct numerically. The source neither returns it exactly nor proves the transformed matrix diagonal. |
| 78–83 | 466–491 | Verified with statistical caveat | Counts, reconstructed density matrix, and similarity worked. One sample does not establish a fidelity scaling law; trace/positivity should be asserted explicitly. |
| 84 | 504–511 | Verified | Teleportation circuit constructed. |
| 85 | 513–515 | Static/visual | Teleportation diagram not image-verified. |
| 86–87 | 519–525 | Verified | Symbolic input state and Bloch vector worked. |
| 88 | 527–529 | Current failure | `QuantumState["00"]` hangs. Explicit register form verifies the intended teleportation result. |
| 89 | 531–533 | Corrected verification | Bob's Bloch vector matches the input after correcting block 88. |
| 90–94 | 540–566 | Verified | Coherent QEC circuit components constructed. |
| 95 | 568–570 | Static/visual | Circuit diagram not used as proof. |
| 96 | 572–574 | Current failure | Another `QuantumState["00"]` hang. |
| 97–98 | 576–587 | Corrected but scope-limited | Syndrome probabilities and unit fidelity verified under the stated real-amplitude assumptions. The claimed complex-amplitude extension remains open. |
| 99–102 | 592–606 | Verified | Unmeasurement circuit components constructed. |
| 103 | 608–610 | Static/visual | Diagram-only block. |
| 104 | 612–614 | Verified | Final unitary output equals the input state in the ideal coherent model. |
| 105–110 | 629–651 | Verified | Grover setup and iteration count worked. |
| 111 | 653–655 | Static/visual | Flattened circuit diagram not image-verified. |
| 112–115 | 657–673 | Verified with maintainability caveat | Tensor-network application, exact/numeric success probability, and network extraction worked. Block 114 uses `%`, so it depends on notebook evaluation history. |
| 116 | 675–680 | Static/visual/API | Contraction-tree construction inspected; graphic and option forwarding were not image-verified. |
| 117 | 685–689 | Verified | Five-qubit cluster circuit constructed. |
| 118 | 691–693 | Static/visual | Cluster diagram not image-verified. |
| 119–125 | 697–730 | Verified with benchmark correction | Stabilizer state, generators, expectations, scale setup, timing, and entropy worked. Gate count/timing/asymptotic prose is inaccurate. |
| 126 | 735–737 | Verified | Three-qubit cluster state constructed. |
| 127 | 739–741 | Static/visual | State plot only. |
| 128 | 745–747 | Verified | Stabilizer representation constructed. |
| 129 | 749–751 | Static/formatted | Formatted state output not visually certified. |
| 130 | 753–755 | Verified | Converted state is structurally equal to the original. |
| 131 | 759–761 | Static/visual | Stabilizer circuit diagram not image-verified. |
| 132 | 772–774 | Current failure | The overloaded state arithmetic/equality expression hangs in the live tree. |
| 133–142 | 779–827 | Verified | Symbolic exponential, algebra, decomposition, commutator, Hamiltonian recovery, and state-calculus residual worked. |
| 143 | 832–834 | Current failure | Direct symbolic `RY` application hangs. |
| 144–145 | 836–842 | Independently checked | Depend on block 143. Energy `Cos[θ]+Sin[θ]` and maximum `Sqrt[2]` are correct. |
| 146–147 | 844–850 | Verified | Werner realignment criterion and witnessed `0<p<1/2` region are correct under the chosen convention. |
| 148–151 | 861–889 | Static/visual/API | Bloch, amplitude, probability, and tableau views were not image-verified; property heads and option forwarding remain only statically inspected. |
| 152 | 894–896 | Verified but nondeterministic | Bare phase estimation uses a random default unitary and is not seeded. |
| 153–155 | 898–910 | Static/visual/API | Circuit, tensor-network graph, and hypergraph views were not visually or topologically certified. |
| 156 | 923–925 | Executed; nonconforming output | Native export runs, but its output should not be labeled standard OpenQASM 3. |
| 157–158 | 930–961 | Internally verified; nonstandard input | QF's permissive importer accepts the string; CHSH value is correct. Independent standards compliance fails. |
| 159–161 | 966–977 | Dependency-blocked | Qiskit/Python bridge unavailable. Source inspection confirms the dependency and contradicts the prose. |
| 162 | 982–986 | External/authenticated | Runtime service query not executed. |
| 163 | 988–990 | Verified locally | GHZ circuit object constructed. |
| 164–165 | 992–998 | External/authenticated | Hardware submission and status were not run. |
| 166 | 1,000–1,002 | Verified locally | Ideal GHZ measurement gives probability one-half at `000` and `111`. |
| 167–168 | 1,004–1,012 | External/authenticated | Historical result retrieval and histogram cannot be audited without captured job data. |
| 169 | 1,039–1,047 | Static; deliberately not run | Installation block mutates the environment and is unpinned/order-dependent. |

## Prose-claim and editorial assessment

Every sentence was read. The findings below dispose of the recurring scientific and editorial claims; explanatory sentences inherit the status of the code or formula they describe in the ledger. The document's overall progression is good: each claim has a concrete model, and the prose usually explains why the output matters. The tone is energetic and accessible. The strongest sentences connect a formula to a physical invariant or computational consequence. The weakest sentences substitute universal quantifiers for evidence.

### Prose/API claims not demonstrated by their adjacent cells

| Source assertion | Disposition |
|---|---|
| “About three dozen named bases are provided” | Plausible from implementation, but not counted or versioned in the manuscript; open until a generated capability list is attached. |
| The same basis mechanism handles noncomputational QASM measurements and Bell-derived KAK compilation | Static implementation/API claim; no adjacent execution establishes either path. |
| Arbitrary observables, degenerate rank-`k` projectors, POVMs, Naimark dilations, and qudits all share the reversible record model | Much broader than the qubit examples. Needs one assertion per distinct construction and a statement of dimensional/outcome limits. |
| The tomography object performs the stated maximum-likelihood procedure | The returned object/matrix was checked; optimizer, likelihood, constraints, convergence, and uncertainty method were not audited from the example. |
| Phase-space surfaces have the stated normalization, ceiling, lobe heights, and global sign properties | Only selected point values were checked. The integral and global claims remain open pending converged quadrature/grid tests. |
| Every graphics property forwards Wolfram graphics options | Not demonstrated. It requires documented option lists and representative forwarding tests. |
| Tensor-network graph/hypergraph structure depicts the actual contraction computation | Objects were constructed, but topology and semantic correspondence were not independently checked. |
| All three computational representations convert exactly in general | One three-qubit stabilizer/state-vector equality worked. That example does not prove unrestricted bidirectional conversion across symbols, qudits, arbitrary channels, or non-Clifford circuits. |
| Foreign circuits enter and leave unchanged | False as a universal statement: lowering, supported-subset restrictions, concrete-parameter requirements, parser extensions, and transpilation can alter or discard structure. |
| The shown hardware histogram has the described error pattern | Not auditable without the historical counts/result artifact. |

### Opening and thesis

- **Line 1, “What Makes Wolfram QuantumFramework Unique: A Verified Showcase”:** “verified” is not presently supported by the versioning and execution evidence. The broader phrase “all of quantum theory” appears in Claim 2 and is not literal: the document covers a broad finite-dimensional slice plus truncated oscillator examples, not relativistic QFT, infinite-dimensional domain issues, path integrals, quantum gravity, or every algorithmic regime.
- **Line 3:** the first two sentences are clear. The third overclaims exactness across numerical simulation and finite measurement. Replace it with a layered thesis: one object system can retain exact structure when available, carry controlled numerical approximations when necessary, and represent sampled/experimental data without pretending those data are exact.
- **Lines 5–16:** the list is useful but mixes a software philosophy, API facts, benchmark claims, and hardware claims. Add version and evidence links at first mention. “Runs on actual IBM quantum hardware” should mean “can submit supported circuits through configured IBM credentials,” not that every object in the showcase can run there.

### Claim headings and transitions

- The seven-claim structure is readable and worth retaining.
- Repeated constructions such as “one object,” “the same object,” and “without translation” become misleading where the code actually converts between representations, lowers a circuit to QASM, uses a Python bridge, or reconstructs a state statistically. “Coordinated object model” is more accurate than “no translation layer.”
- Sentences saying competitors “cannot” do something must include a reproducible comparison or be reframed as a feature statement about QF.
- “Exact,” “ordinary,” “automatic,” “unchanged,” and “every” should be reserved for statements with a test that actually enforces the scope.

### Mathematical exposition

- The Rabi exposition is concise and correct.
- The QAOA section should specify the objective convention and angle domains before discussing periods and maxima.
- The entropy section should teach the continuous `0 Log[0] = 0` convention explicitly.
- The basis section needs a single declared row/column convention.
- The tomography section should distinguish generated multinomial counts, element-wise estimator error, state similarity, and statistical confidence.
- The stabilizer section should separate asymptotic algorithm class from one-machine timing.
- The “boundary of exactness” section should be replaced rather than lightly edited; its organizing mathematical idea is wrong.

### Measurement language

- Use **premeasurement** for a unitary system-record correlation before terminal readout.
- Use **deferred measurement** for replacing intermediate measurement/feed-forward with coherent controls under stated conditions.
- Use **uncomputation** for reversing a still-coherent premeasurement.
- Reserve **measurement reversal** or **uncollapse** for protocols whose probability, strength, conditioning, and physical record are stated.
- Avoid saying a classical channel “is not required” when a coherent control wire supplies the same causal dependence.

### Reproducibility language

- “Every cell below was run” is a provenance claim, not a stylistic flourish. It must be backed by a dated, machine-readable run log against an immutable commit.
- “The only external requirements” is false as written. Depending on the cells, the document also needs compatible Arrays, TensorNetworks, IBMQuantumPlatform, Python, Qiskit, network access, an IBM account/token, and a backend that still accepts the circuit.
- A random default object must be seeded or made explicit.
- A graph or histogram shown as evidence should have its generating data captured alongside the manuscript.

### Style and maintainability

- The prose has few copy-editing defects. The main problem is evidentiary precision, not grammar.
- Long paragraphs near the end combine mathematical exactness, software architecture, benchmarking, and philosophy. Split these into independently falsifiable sentences.
- Keep the candid finite-overload caveat in Claim 5; it is one of the document's most credible passages.
- Avoid global state such as `$Assumptions`; a showcase should be rerunnable cell-by-cell and in fresh kernels.
- Replace block 114's `%` with an explicit named result. Reusing global names such as `psi`, `n`, and `true` across sections makes isolated or out-of-order execution fragile; scope examples or give section-specific names.
- Add expected-output assertions after every nonvisual example. A displayed expression is not a test.

## Reproducibility and verification assessment

The adjacent `verify-showcase.wls` is valuable: it contains many explicit assertions and confirms that a substantial mathematical core worked in an earlier 2.0.2-era environment. It does not establish the present document's warranty because:

1. it is not tied to an immutable QF commit or complete dependency lock;
2. its own header identifies a different QF version from the Markdown;
3. it uses the same now-broken string state shorthands;
4. it does not independently validate OpenQASM conformance;
5. external hardware results are not stored as auditable fixtures;
6. graphics are not checked against numerical invariants; and
7. no timeout wrapper prevents a regression from hanging the full suite.

The installation block should be replaced with a reproducible manifest. At minimum it should record:

- exact QF, Arrays, TensorNetworks, and IBMQuantumPlatform commits or paclet hashes;
- Wolfram Language version, operating system, and architecture;
- Python and Qiskit versions for Qiskit-dependent cells;
- random seeds;
- hardware backend, job ID, submission/result timestamps, and sanitized calibration metadata;
- whether each timing is cold, warm, compiled, or inclusive of bridge/interface overhead;
- expected numerical tolerances and truncation dimensions; and
- a fresh-kernel command that exits nonzero on any failure or timeout.

The current order is unsafe for a verification recipe: it installs and loads QF before installing the development TensorNetworks build, so the active kernel can retain an already-loaded published dependency. Forced version installation also changes the user's paclet environment. Install into an isolated location, load all pinned directories before `Needs`, and run the suite in a fresh kernel.

## Confirmed strengths worth preserving

The necessary revisions do not erase the document's strongest contribution. The following parts are good evidence for QF when presented with bounded wording and passing assertions:

- exact symbolic Rabi dynamics and the associated invariant checks;
- special-function Landau–Zener propagation, separated from its numerical visualization;
- the analytic `K4` QAOA landscape and cubic optimum;
- open-system channel/Lindblad calculations and SIC-POVM probabilities;
- oscillator identities and unified Wigner/Husimi access with explicit truncation;
- spin-frame invariance and exact zero-field-splitting spectrum;
- simulated tomography with transparent shot counts;
- Grover success probability and conversion to tensor-network form;
- cluster-state stabilizer generators, entropy, and exact representation conversion;
- symbolic operator arithmetic, commutators, decomposition, and Hamiltonian reconstruction;
- Werner-state entanglement criterion and CHSH value; and
- ideal local GHZ construction and measurement probabilities.

These examples collectively support a more modest and stronger thesis: **QF provides a coherent Wolfram Language object model that often preserves symbolic structure while allowing selected numerical, sampled, tensor-network, stabilizer, interoperability, and hardware workflows.** That sentence is impressive and testable.

## Required revision plan

### P0 — required before calling the document verified

1. Fix or bypass the multi-qubit string-state regression and blocks 132/143; add timeout-protected regression tests.
2. Remove or rewrite the radicals/Abel–Ruffini boundary paragraph.
3. Correct the QiskitTarget/Python statement and OpenQASM syntax/emitter.
4. Recast teleportation, QEC, and unmeasurement as coherent premeasurement/deferred-measurement constructions with explicit resource assumptions.
5. Pin the complete software environment and replace the installation block with a fresh-kernel, isolated verifier.
6. Attach real hardware provenance or remove the historical output claim.
7. Rerun all 169 blocks, store the log, and make every nonvisual claim an assertion.

### P1 — scientific and performance corrections

1. Repair the row/column basis convention.
2. Correct the ring “exact” statement and assert full diagonalization.
3. Separate tomography element error from fidelity/infidelity scaling.
4. Qualify finite-cutoff `g^(2)` and phase-space numbers.
5. Fix Landau–Zener terminology and numeric function definition.
6. Scope the Lindblad “every state” statement and localize assumptions.
7. Replace the simulator/stabilizer implementation narrative with source-accurate algorithms.
8. Rebuild the stabilizer benchmark with exact gate count, cold/warm separation, repetitions, machine metadata, and a valid companion link.

### P2 — editorial strengthening

1. Replace universal quantifiers with versioned, bounded statements.
2. Seed or explicitly specify the phase-estimation unitary.
3. Name and version every competitor or remove the comparison.
4. Make external dependencies visible at the start of the relevant section.
5. Split philosophy from empirical claims and link each empirical claim to a check.

## Final judgment

The document is **scientifically promising, technically substantial, and presently overclaimed**. Its verified mathematical core is larger than its failures, but the failures occur exactly where the manuscript makes its strongest warranties: universal exactness, automatic representation choice, standard interchange, Python-free transpilation, reversible measurement, and reproducible execution.

The right response is not to discard the showcase. It is to narrow the thesis, fix the current regressions, correct the standards and measurement language, and turn the examples into a pinned test artifact. With those revisions, this can become a credible demonstration of QuantumFramework's distinctive strength: symbolic quantum objects that remain connected to multiple computational representations without hiding where exactness ends and approximation, inference, or external tooling begins.
