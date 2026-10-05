# Brief 3: building a Liouvillian (investigate first)

Read `README.md` in this folder first: it has the workflow, the house rules and the traps every brief shares.

## Goal

Since `191f1ac3` and `f6299b43`, the exponential of a pure-dephasing Liouvillian on 7 qubits takes about 0.02 s. Building that Liouvillian takes about 10 s, so construction is now nearly the whole cost of an open-system evolution of this kind. The task: find where the construction time goes, then make it scale with the size of the superoperator rather than with repeated operator arithmetic.

## What QF does today

Measured on main at `f6299b43` by `baseline.wls liouvillian` (raw output in `baseline-out/liouvillian.txt`), for H a diagonal Ising operator on n qubits and jump operators Z₁, ..., Zₙ at rate 0.2:

| n | full Liouvillian | Hamiltonian part alone (`"Liouvillian"[H]`) | one jump operator alone (`"Liouvillian"[None, {Z₁}, {0.2}]`) | building the n operators Zₖ |
|---|---|---|---|---|
| 5 | 1.8 s | 0.04 s | 0.07 s | 0.13 s |
| 6 | 3.4 s | 0.08 s | 0.15 s | 0.20 s |
| 7 | 10.3 s | 0.28 s | 0.51 s | 0.28 s |

At n = 7 the parts measured alone add up to about 0.28 + 7 × 0.51 + 0.28 ≈ 4.1 s, so about 6 s of the 10.3 s is spent somewhere the table does not isolate.

The construction is `QuantumOperator["Liouvillian"[H, Ls, Gammas]]` in `QuantumFramework/Kernel/QuantumOperator/NamedOperators.m:908-935`. It:

1. pads jump operators that act on a subsystem to the full register (`padJumpOperators`, in `QuantumOperator/QuantumOperator.m`);
2. rewrites H and every jump operator in one common basis (`QuantumOperator[#, basis]`);
3. builds the superoperator of each term, `HamiltonianMixedOperator` (`NamedOperators.m:889-892`) and `LindbladMixedOperator` (`:894-906`), each through QF operator products, transposes and a reshape;
4. adds the terms with `PadRight[gammas, Length[ls], 1] . (LindbladMixedOperator /@ ls)` for a vector of rates, or with a double `MapThread` over pairs for a matrix of rates. Each `+` between two `QuantumOperator`s goes through `addQuantumOperators` and `padQuantumOperators` in `QuantumOperator.m`.

## Investigation (do this first, report before changing code)

1. At n = 7, time each of the four stages above separately, and each `+` of the sum, with `ClearSystemCache[]` before each. State which stage holds the unexplained 6 s.
2. For the stage that dominates, find what it repeats: re-sorting orders, rebuilding bases, densifying a `SparseArray`, or simplifying.
3. Measure how each stage grows from n = 5 to 7, and give the cost the dominant stage would have at n = 8.

## Likely change, to confirm with the investigation

Build the superoperator as one sparse matrix and wrap it in a `QuantumOperator` once. For a jump operator L with rate γ, the dissipator on vec(ρ) is γ (L ⊗ L* − ½ (L†L ⊗ 1 + 1 ⊗ (L†L)ᵀ)), and the Hamiltonian part is −i (H ⊗ 1 − 1 ⊗ Hᵀ). Sum those as `SparseArray`s, in the same index convention `HamiltonianMixedOperator` and `LindbladMixedOperator` use now, which the investigation must pin down by comparing on small cases. That replaces n + 1 operator additions by sparse matrix sums.

## What must not change

- The matrix, basis, order and label of every Liouvillian, compared exactly (symbolic rates and symbolic H) on small cases and to a stated bound on numeric ones.
- Jump operators on a subsystem, padded by `padJumpOperators` (since `236c7db5`).
- A matrix of rates (the Kossakowski form), whose `If` skips terms with a provably zero rate.
- `None` for H, an empty list of jump operators, and named operators carrying their own basis (the comment at `NamedOperators.m:912-915`).
- `QuantumOperator["Hamiltonian"[...]]` (`:937`), which is i times the Liouvillian.

## Tests to write

`Tests/Liouvillian.wlt` already covers construction; read it first and keep all of it passing. Add:

- The 7-qubit dephasing Liouvillian builds inside a time limit chosen from the investigation, and its exponential's diagonal entries are e^(t Lₖₖ).
- Equality with today's construction, computed on the unchanged copy and stored as the expected value, for symbolic 2-qubit cases: dephasing, amplitude damping, a jump operator on one qubit of two, and a 2 × 2 matrix of rates.

## How to verify and land

Follow the README. Report the stage timings of the investigation, then the before and after times of the table above.

## Traps

- `"Liouvillian"` stores a matrix-type operator. Since `f6299b43` its `"Matrix"` stays a `SparseArray`; keep it sparse to the end, or the 7-qubit exponential goes back to holding a dense 16384 × 16384 list.
- Machine rates combined with exact operators give machine superoperators; compare numeric results against a bound.
