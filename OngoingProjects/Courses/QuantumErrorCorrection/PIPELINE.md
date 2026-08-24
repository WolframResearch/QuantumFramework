# PIPELINE: Authoring and Verification Contract for Quantum Error Correction

This file governs `Question-List.md`, every plan, every benchmark, and all future worked answers.
Markdown is canonical. Generated notebooks, if added later, are build artifacts.

## 0. Governing principle

The reader learns QEC by following the causal chain

$$
\text{physical noise} \longrightarrow \text{syndrome information} \longrightarrow
\text{decoder inference} \longrightarrow \text{recovery or frame update}
\longrightarrow \text{logical failure rate}.
$$

A code name by itself teaches almost nothing. Every computational answer must make the code, noise,
available information, decoder or recovery, and logical success criterion explicit. Every theorem
answer must make its assumptions, conclusion, and failure when an assumption is removed explicit.

## 1. Scope and reader

The course begins with a compact bridge from finite-dimensional quantum mechanics and classical
coding, develops the stable QEC core, and ends with a dated research frontier. It is not a catalogue
of named codes and not a survey whose conclusions depend on brand or platform claims.

The reader is intelligent and technically mature but may be new to error correction. No answer may
silently assume binary symplectic notation, a decoder convention, a noise convention, or a meaning
of "threshold" that has not been defined locally.

## 2. Question contract

Each question has one central task. Its prompt or answer plan must pin the following fields:

1. **Object:** code, code space, circuit, channel, decoder, theorem, or experiment.
2. **Noise:** an explicit error set or stochastic/coherent channel when noise is relevant.
3. **Information:** syndrome bits, analogue readout, erasure locations, final data measurement, or
   no classical side information.
4. **Action:** derive, construct, simulate, decode, compare, prove, or reconstruct.
5. **Output:** a matrix, circuit, logical operator, recovery, bound, curve, confidence interval, or
   sharply stated conclusion.
6. **Success criterion:** an invariant or benchmark capable of rejecting a wrong result.
7. **Scope boundary:** what the question deliberately does not claim.

If those fields cannot be stated without adding a second independent task, split the question.

## 3. Evidence routes

Every answer declares one primary route from `Plans/Evidence-Route-Table.md`:

- `Q0`: exact finite-field or Pauli algebra.
- `Q1`: exact finite-dimensional channel or recovery calculation.
- `Q2`: circuit fault propagation and exhaustive bounded enumeration.
- `Q3`: decoder construction on a fixed syndrome model.
- `Q4`: stochastic simulation, threshold estimation, or resource analysis.
- `Q5`: proof or conceptual argument.
- `Q6`: reconstruction of a dated primary result.

Code is evidence only when the literal code is run and its output tests the claim. A plot without a
quantitative observable, a syndrome without a logical-equivalence check, or two algorithms sharing
the same hidden approximation do not constitute independent verification.

## 4. Conventions

- Arithmetic for binary check matrices is over $\mathrm{GF}(2)$.
- A Pauli operator modulo global phase is represented by
  $(\mathbf{x}\mid\mathbf{z})\in\mathrm{GF}(2)^{2n}$.
- The symplectic product is
  $\mathbf{x}\cdot\mathbf{z}' + \mathbf{z}\cdot\mathbf{x}' \pmod{2}$.
- A stabilizer group never contains $-I$. Generator counts refer to independent generators.
- For a stabilizer code, distance is the minimum weight of an element of the normalizer that is not
  a stabilizer. Subsystem codes must state the corresponding dressed or bare convention.
- A code of distance $d$ corrects arbitrary errors on at most
  $t = \lfloor(d-1)/2\rfloor$ unknown locations and erasures on at most $d-1$ known locations.
- `physical error rate`, `logical error rate`, `per round`, `per cycle`, `per gate`, and `per unit
  time` are never interchanged.
- `below threshold`, `beyond break-even`, and `fault-tolerant` are distinct claims.
- A decoder may return a recovery, a logical class, or a Pauli-frame update; the answer states which.

## 5. Answer structure

Each finished answer follows this order:

1. **Physical question.** One short paragraph stating what must be protected or inferred.
2. **Defining objects.** Introduce the smallest sufficient code, channel, circuit, and conventions.
3. **Derivation or computation.** Make the central reasoning visible. One result per code cell.
4. **Adversarial check.** Test a limit, enumerate the bounded fault set, verify a logical invariant,
   compare with an exact anchor, or run an independently formulated route.
5. **Failure edge.** Show where the method stops working or which assumption carries the result.
6. **Meaning.** State what was learned about protection, inference, or overhead.

No answer refers to another answer for a definition needed to understand it. Cross-references may
guide sequence, but each answer is locally readable.

## 6. Verification requirements

### Exact algebra

- Check dimensions and ranks over the correct field.
- Verify all proposed stabilizers commute.
- Verify logical Paulis commute with stabilizers and have the required mutual anticommutation.
- Distinguish equality from equivalence modulo stabilizers or gauge operators.
- Exhaustively enumerate the claimed bounded error set when its size permits.

### Channels and recovery

- Verify complete positivity and trace preservation when a map is constructed.
- Test the recovery on a symbolic logical state or on an entangled logical-reference state.
- Report entanglement fidelity or a specified logical observable, not only state-vector overlap for
  one basis state.
- Keep coherent, stochastic, leakage, loss, and erasure errors distinct.

### Circuits

- Enumerate every single fault when claiming distance-three fault tolerance.
- Propagate faults through the literal gate order.
- Include preparation, measurement, reset, idle, and leakage locations when the noise model claims
  to be circuit-level.
- Check both residual data weight and resulting logical class.

### Decoders and simulations

- Fix the code instance, noise model, number of rounds, boundary conditions, and decoding weights.
- Store seeds or sufficient replay data.
- Report sample count and uncertainty.
- Sweep code distance before claiming suppression or a threshold.
- Compare decoders on identical generated events and state whether weights use privileged knowledge
  of the data-generating model.
- Search explicitly for error floors and finite-size artifacts.

### Proofs

- State all hypotheses before the conclusion.
- Identify whether the result is exact, asymptotic, approximate, or architecture-dependent.
- Give a counterexample or failure mode when a hypothesis is dropped whenever one is available.

## 7. Wolfram Language and QuantumFramework policy

Use Wolfram Language and QuantumFramework where they expose the physics cleanly: Pauli algebra,
finite-dimensional channels, stabilizer transformations, exact linear algebra, circuit simulation,
and small exhaustive searches. Do not invent framework symbols. A symbol not confirmed in the live
installation is marked `verify at authoring` in a plan and is not printed as working code.

Small transparent implementations of binary symplectic algebra or a decoder are preferable when the
implementation is itself the lesson. External libraries may be used for serious decoder benchmarks,
but their version, configuration, and data interface become part of the benchmark record.

## 8. Frontier policy

Part 16 and `Frontier-Snapshot-2026.md` are pinned to 2026-08-23. Each frontier claim carries one
status:

- `theorem`
- `simulation`
- `architectural proposal`
- `experimental demonstration`
- `open problem`

Every frontier answer cites a primary source, reproduces one bounded result, distinguishes what was
measured from what was inferred, and names at least one limitation. Frontier entries expire after
twelve months or when a cited result is superseded, whichever comes first. Updating the frontier
must not reorder the stable prerequisite spine without a separate pedagogical argument.

## 9. Review loop

Each Part receives at most three passes:

1. **Independent physics pass:** rederive equations and recompute examples without following the
   draft's route.
2. **Adversarial pedagogy pass:** look for bundled tasks, missing prerequisites, false mental models,
   and checks that cannot reject the proposed answer.
3. **Literal execution pass:** run all code and benchmarks from a clean state.

A Part ships only when no fatal or misleading defect remains. Unresolved limitations are reported,
not hidden by prose.

## 10. Release order

1. Author and audit the vertical pilot in `Plans/Pilot-Vertical-Slice.md`.
2. Freeze the stable core, Parts 0 through 10.
3. Author topological codes and decoding, Parts 11 and 12.
4. Author the qLDPC, bosonic, and hardware-adapted branches, Parts 13 through 15.
5. Re-audit primary sources and author the dated frontier, Part 16.
6. Build any notebooks only after their Markdown answers pass verification.

## 11. Pre-ship checklist

- [ ] The question is atomic and prerequisite-safe.
- [ ] The object, noise, information, output, and success criterion are explicit.
- [ ] The evidence route is declared and its known traps are addressed.
- [ ] Every equation and example was independently reproduced.
- [ ] Logical equivalence is checked, not guessed from the syndrome.
- [ ] Simulation settings and uncertainty are recorded.
- [ ] `below threshold`, `break-even`, and `fault-tolerant` are used precisely.
- [ ] Frontier claims have status, primary source, snapshot date, and limitation.
- [ ] No answer leaks unexplained framework machinery into the pedagogy.
- [ ] The closing paragraph states physical meaning rather than narrating the computation.
