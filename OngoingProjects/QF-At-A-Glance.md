# Wolfram QuantumFramework: A Physicist's Overview

A compact answer to four questions, addressed to a physicist. Grounded in the QF-authored
reports under `OngoingProjects/` and the kernel audit under `audit/`, all verified against the
2.0.0 working tree at HEAD `f9dc1cdc` (2026-06-14). Numbers are wall-clock minimums on one
machine (Apple M2 Pro), so they understate variance but are the standard best-case comparison.

---

## 1. What is QuantumFramework?

QuantumFramework (QF) is a quantum-theory laboratory built into the Wolfram Language. The objects
of quantum theory are first-class **symbolic expressions**: a `QuantumState`, a `QuantumOperator`,
a `QuantumChannel`, a `QuantumMeasurementOperator`, a `QuantumCircuitOperator`, a `PauliStabilizer`.
They are exact and symbolic by default, they carry free parameters, they compose in closed form, and
each answers questions about itself through one uniform interface, `obj["Property"]`.

The defining consequence is that a computation typically ends in a **formula**, not a table of
numbers: the propagator $U(t)$ of a time-dependent Hamiltonian, an entanglement entropy $S(t)$ you
can differentiate, an exact variational optimum that is a `Root` of a polynomial. Where a numerical
run or a finite set of shots returns an estimate, QF returns the exact object that estimate
converges to.

Two generalities run underneath the whole model:

- **Native qudits.** Every object lives in any finite dimension $d$, qubits, qutrits, spin-$j$
  multiplets, or $N$-level systems alike. The others in the field are essentially qubit-centric.
- **Basis as a first-class attribute.** Every object carries the basis it is written in, and QF
  reconciles mismatched bases automatically through the computational frame, $A_{\mathcal B} =
  B^{-1} A B$. About three dozen named bases ship (Bell, Fourier, Gell-Mann, MUB, SIC, Wigner, ...).

The scope is the whole theory under one roof: the gate model, open-system Lindblad dynamics
(`QuantumEvolve`), generalized measurements (POVMs), continuous-variable / Fock space
(`AnnihilationOperator`, `CoherentState`), phase-space quasi-probabilities ($W$, Husimi $Q$), the
stabilizer formalism (`PauliStabilizer`), and tomography (`QuantumStateEstimate`). All of it is
ordinary Wolfram Language, so `Exp`, `D`, `Plus`, `Commutator`, `Maximize`, `Reduce`, `Plot`, and
`DSolve` act directly on the objects. That last point is the real moat: it is not a feature any
standalone library can match.

---

## 2. Where it stands versus other quantum packages

The honest framing is that QF and the usual comparison targets solve **different problems**. The
others are numerically specialized engines; QF is a symbolic object algebra. The boundary shows up
cleanly in three head-to-head benchmarks on identical tasks.

| Regime | Task | QF vs the specialist | Verdict |
|---|---|---|---|
| State vector | $105$-gate circuit, full readout | $\approx 130\times$ slower than Qiskit Aer at $n=12$ ($283\,\mathrm{ms}$ vs $2.1\,\mathrm{ms}$), $\approx 320\times$ at $n=16$ | behind by $\sim 10^2\times$ |
| Stabilizer | identical Clifford stream | within $2.5\times$ to $5.3\times$ of Stim; **beats** QuantumClifford.jl at $n=1000$ ($58$ vs $108\,\mathrm{ms}$) | genuinely competitive |
| Open system | driven, damped qubit (Lindblad) | $\approx 9\times$ slower than QuTiP ($17$ vs $1.9\,\mathrm{ms}$); physics agrees to $4\times10^{-4}$ | competitive at small $n$ |

Two of these moved a lot in June 2026. The **stabilizer engine** was rewritten (bit-packed tableau
words plus a compile-to-C bulk gate fold, `ps["ApplyCircuit"]`), a cumulative speedup of $266\times$
to $1095\times$ over the old path, taking it from $\sim 1000\times$ slower than Stim into single-digit
factors. The **state-vector path** got a property-cache refactor that cut the warm $n=12$ apply from
$836\,\mathrm{ms}$ to $273\,\mathrm{ms}$; the residual gap to Aer is now genuine numeric work (dense
amplitudes held in sparse containers) plus per-application network re-assembly, not framework
bookkeeping. A packed numeric fold prototyped at $178\times$ (into Aer's envelope) is the one large
win still on the table, validated but not yet shipped.

**Where QF is structurally unique** (the numeric specialists cannot do these): exact symbolic object
algebra with free parameters; native qudits; one object model across pure states, density matrices,
CPTP channels, POVMs, Fock space, and phase space; basis as a reconciled attribute; second
quantization and phase space in the same framework as gates; multiway / causal-graph circuit views;
and full Wolfram-Language integration.

**Where QF genuinely lags** (verified absent in 2.0.0): error mitigation (no ZNE, PEC, measurement
mitigation, dynamical decoupling); a native transpiler (the device path delegates to qiskit); open-
system breadth (no Monte-Carlo trajectories, steady state, Bloch-Redfield, or Floquet, all QuTiP
strengths); optimal control (no GRAPE / CRAB / Krotov); QAOA as a packaged routine; differentiable
programming (no autodiff or PyTorch / JAX bridge, QF uses parameter-shift); GPU; and bulk Pauli-frame
sampling for QEC shot generation (Stim's and QC.jl's reason to exist).

Rule of thumb: QF for symbolic derivation, teaching, qudits, exact algebra, and anything that lives
inside Wolfram; Aer/Cirq for the largest gate-model runs; QuTiP for open systems at scale; PennyLane
for QML / autodiff; Stim for QEC shot sampling and detector error models. QF's `PauliStabilizer` is
now a real option for raw Clifford simulation cross-validated against both Stim and QC.jl.

---

## 3. Where it stands with respect to QPU hardware

QF is an **interoperability hub** that reaches real devices, with one important honesty caveat.

- **OpenQASM, native.** `QuantumQASM` imports and exports OpenQASM 2/3 in pure Wolfram Language, no
  Python and no external account. Import keeps angles symbolic (`Pi/2` stays the exact constant), so
  a foreign circuit becomes an exact QF object. Export of OpenQASM 3.0 is native; a 2.0 dump or a
  faithful qiskit dump routes through the Python bridge.
- **Device-aware transpilation.** `QiskitTarget` carries a validated device spec (qubit count, native
  gate set, coupling map); `QuantumQASM[circuit, "Target" -> target]` transpiles to that backend's
  native gates and connectivity. This is real and current, but it **delegates to qiskit** through the
  Python interop. QF has no native SABRE-style router, so the honest label is "via bridge." A
  working-tree refinement (not yet committed) makes the bridge **error-aware**: an
  `"OptimizationLevel"` option, transpilation against the backend's populated qiskit `Target`
  (per-instruction error and duration), and an offline `GenericV2` fake device for credential-free
  testing.
- **Submission to real hardware.** `IBMJobSubmit[circuit, backend, "Shots" -> ...]` sends a circuit
  asynchronously through qiskit's SamplerV2 / EstimatorV2 and returns an `IBMJob` handle; the
  completed result comes back as the **same kind of `QuantumMeasurement` object** the exact simulation
  produces, decoded into ascending-qubit order, so hardware noise is one chart away from the ideal
  distribution. A measurement in QF is just a list of qubit indices placed in the circuit.

Caveats worth stating to a user: the device path is bridge-based, not native; `IBMJob["ExecutedCircuit"]`
(the server-transpiled circuit) is unreliable for QF's ISA+V2 jobs because IBM stores no such object
for them and the download URL is VPC-private (a guard now distinguishes the two failure modes); and
AWS Braket support is only partial. The dependency surface is heavier than a pure-WL tool: the
TensorNetworks paclet and `IBMQuantumPlatform` are auto-installed at load, and qiskit / pyzx / classiq
integrations each need their own Python environment.

---

## 4. How it is designed fundamentally

QF is an **object algebra**. Each head is a thin symbolic wrapper around a smaller one, so the whole
system is a few primitives composed:

- `QuantumState` $=$ a `SparseArray` (a vector for pure, a density matrix for mixed) plus a
  `QuantumBasis`. The basis itself wraps two `QuditBasis` factors (output and input), and a
  `QuditBasis` maps each (name, index) pair to its tensor representation. Dimension semantics are
  multiplicative across qudits.
- `QuantumOperator` $=$ a `QuantumState` plus a 2-tuple of qudit orders $\{\text{out}, \text{in}\}$.
  An operator is literally "a state with wires," which is why states can sit inside circuits and a
  bra-ket sandwich $\langle\psi|H|\psi\rangle$ is a single composable object.
- `QuantumChannel` $=$ a `QuantumOperator` in **Stinespring dilation form**: extra output qudits
  (non-positive order indices) are the environment to be traced. A CPTP map is an isometry plus a
  deferred partial trace.
- `QuantumMeasurementOperator` realizes measurement by **von Neumann unitary dilation**, not by a
  collapse postulate: it appends a pointer qudit and entangles it, $V|\psi\rangle = \sum_k
  (P_k|\psi\rangle)\otimes|k\rangle$. The measurement record is therefore itself a quantum wire
  (numbered $0, -1, -2, \dots$), feedforward is just a gate controlled on a record, and because $V$
  is an isometry ($V^\dagger V = I$) a measurement is reversible until the record is discarded. The
  price is explicit: $k$ measurements on $n$ qubits live in dimension $2^{n+k}$, so for many
  mid-circuit measurements one drops the pointer with `"DiscardExtraQudits"` or works in the
  channel / POVM picture to stay at $O(4^n)$.
- `QuantumCircuitOperator` $=$ an ordered list of operators and barriers.

**Three interchangeable engines** compute a circuit, selected by `Method`, all sharing the one object
model so their representations convert exactly:

1. **TensorNetwork (default).** The circuit is converted to a tensor network and contracted by the
   separate `Wolfram/TensorNetworks` paclet; an index is a (vertex, wire) pair, and the contraction
   `"Path"` (greedy / optimal / explicit) is a user-exposed sub-option that sets the cost via the
   largest intermediate tensor. This is the engine that wins on locality (light-cone observables,
   shallow circuits) and the one that can blow up memory on the wrong contraction order (HHL / QPE at
   precision $\geq 5$).
2. **Schrödinger.** A literal gate-by-gate fold over the dense $2^n$ amplitude vector. Despite the
   name it is the **slowest** dense option ($\sim 23\times$ the TN default, with no speed crossover at
   any size), because it pays one full framework apply per gate. Its only real niche is memory safety
   on circuits whose TN contraction over-bloats.
3. **Stabilizer.** A compiled Aaronson-Gottesman tableau (`PauliStabilizer`), $O(n^2)$ to $O(n^3)$
   per Clifford gate, carrying a thousand qubits with ease. It integrates by UpValues on the
   stabilizer heads, deliberately bypassing the basis machinery that would re-impose an $O(2^n)$ cost.

Two cross-cutting design facts shape performance and behavior. Property access goes through a
**dispatch-and-cache** layer: $O(1)$ structural getters (`"Basis"`, `"Order"`, `"State"`) now have
direct rules and skip the cache, and the cache store only commits keys free of pattern symbols.
And **everything is symbolic by default**: matrices are `SparseArray`, the generic eigensolver is not
Hermitian-specialized, and free parameters propagate untouched, which is exactly what buys the
closed-form results and exactly what costs the dense numeric throughput.

---

### Scope note on sources

I read the QF-authored reports in `OngoingProjects/` and the active `audit/` set fully, focusing on
the most recently updated (the master platform comparison, the showcase, the apply-path and
Schrödinger diagnoses, the stabilizer reports, the measurement-dilation note, and the architecture /
TN-integration / function-index / activity-status audit files). The `OngoingProjects/Stabilizer/
External Packages/` tree contains roughly sixty vendored copies of **other** tools' documentation
(Stim, QuantumClifford.jl); those are third-party reference, not QF content, and are not summarized
here.
