# Proposal to the QF side: conjugating a Pauli through a Clifford, with its sign

*From the QEC layer (`OngoingProjects/QEC`), step 6 of `QEC-API-Redesign-Plan.md`. A proposal
only: nothing in `QuantumFramework/Kernel/` was changed.*

## What the QEC layer needs

One operation: given a Clifford unitary `U` on `m` qubits and a Pauli `P` on the same qubits,
return `U P U†` **as a Pauli with its phase** — `±Q`, or `i^e Q` for a non-Hermitian input.

The QEC layer uses it for every transversal-gate question (does `U^⊗n` map the stabilizer
group to itself, *with signs*; what logical gate does it perform) and for the register's
transversal CNOT. The sign is the whole point: on the Steane code `S^⊗7` sends `X̄` to `−Ȳ`,
which is why transversal `S` is logical `S†` (Got26 eq. 11.22), and an image equal to minus a
generator preserves the stabilizer *as a set of Paulis* while destroying the code space.

## What the layer does today, and why it is a workaround

`transversalAction` (in `QECCore/Transversal.wl`) builds the dense matrix of the gate,
conjugates each one- or two-qubit Pauli densely, and reads the image and its phase off traces:
`Tr[σ · U P U†] / 2^m`. It is correct and pinned by 24 tests, but it is a dense computation
standing in for a tableau one, and it only works because the gates are one or two qubits wide.
The frame propagator in the engine (`framePropagate`) is not usable for this: it tracks Pauli
*frames*, which by design forget signs.

## What QF already has

`CliffordChannel` (Kernel/Stabilizer/CliffordChannel.m, after Yashin25) already carries the
right data. For a Clifford unitary with `|A| = |B| = n`, its tableau has `2n` rows
`[u_A | u_B | c]`: each row pairs a Pauli generator with its conjugate, and `c` is the sign.
That table *is* the conjugation map on generators; images of products follow by the row-sum
phase rule the file already implements (`stabilizerRowSumAGPhase`).

What is missing is only the two ends:

1. **A constructor from a Clifford gate.** Today `CliffordChannel` is built from a
   `PauliStabilizer`, as the identity, or from a `QuantumChannel` whose label is a single Pauli
   string (the "deterministic Pauli" case); a Clifford `QuantumOperator` (`"H"`, `"S"`,
   `"CNOT"`, a `QuantumCircuitOperator` of Clifford gates) has no route in. Proposed:
   `CliffordChannel[op_QuantumOperator]` for Clifford `op`, filling the `2n` rows by
   conjugating `X_j` and `Z_j` — once, per gate, by the engine's own stabilizer gate updates
   (`GateUpdates.m` already has H, S, CNOT on tableaux), so no dense matrix is needed.
2. **A query.** `cc["Conjugate", P]` (or `cc[P]` on a Pauli string) returning the image as a
   phased Pauli: decompose `P` over the input generators, multiply the images with the row-sum
   phase, carry `i^e` through for a non-Hermitian `P`.

Both are small, separable, and testable against the dense definition on every one- and
two-qubit Clifford (24 + 11520), which is exactly the check the QEC tests already run.

## What it would retire

`transversalAction`, `transversalPairConjugate`'s dense path and the trace read-out in
`Transversal.wl`, about 80 lines, and the restriction to one- and two-qubit gates: with a
tableau route, a transversal *circuit* per qubit (not just a single gate) becomes as cheap as
a gate.

## What it would not change

The QEC layer's own Pauli rows (`{x | z | e}` with `e ∈ ℤ₄`) stay as they are; the proposal
asks the engine to answer in its own sign convention, and the layer converts (`e = 2c` for a
Hermitian image). Nothing else in the layer depends on it.
