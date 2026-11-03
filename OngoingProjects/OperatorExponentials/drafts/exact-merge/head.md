# Operator Exponentials: Formula Audit and Exact QuantumFramework Code

Each entry below is one distinct formula from the atlas [Bosonic-Fermionic-Operator-Exponentials.md](Bosonic-Fermionic-Operator-Exponentials.md). It gives a short statement of what the formula says, the formula, every atlas ID where it appears, an audit verdict with how it was checked, and the QuantumFramework (QF) code that produces it with the kernel's output. When several books state the same identity, the entry shows one version first and each other book's version under "Also stated in". All 428 atlas IDs are covered: {N} formula entries, of which {QF} have an exact QF cell, plus {NF} records that carry no formula.

Everything here is exact. No truncated Fock matrix is used anywhere, neither as QF code nor as a check. A QF cell uses only the symbolic second-quantization algebra (`**`, `FieldVariables`, `BosonicNormalOrder`, `BosonicExpOrder`, `BosonicExpSimplify`, `BosonicExpBraid`, `BosonicMatrixElement`, `BosonicVEV`, `NormalOrdered`, and `DisplacementOperator` / `SqueezeOperator` with size `Infinity`). Where QF has no exact route, the entry says "QF code (exact): none" and names where QF stops; the formula is then checked by exact Wolfram computation (symbolic algebra, series to a stated order, exact integrals, Bargmann or position realizations acting on symbolic functions, exact 2×2 or 2N×2N group matrices, exact Jordan-Wigner matrices for a few fermion modes) or by a hand derivation. In a few places two exact closed forms are compared at random points to high precision because `FullSimplify` stalls; the entry says so.

## Loading the QF code

Every cell assumes these two lines have been run. Without the last `Needs` the second-quantization symbols stay unevaluated and raise no error.

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

## Conventions to fix before comparing with QF

- `SqueezeOperator[ξ, Infinity]` is exp((ξ̄ a² − ξ a†²)/2), the sign of the atlas's common convention (X-002) and of Gerry and Knight. Combescure and Robert's S(ξ) in C-025, C-027 and C-048 is QF's `SqueezeOperator[-ξ, Infinity]`.
- `DisplacementOperator[α, Infinity]` satisfies D(α)†aD(α) = a + α. Berezin's implementer, oriented as U a U⁻¹ = a + f, is QF's D(−f). Without `Infinity`, `DisplacementOperator[α]` is a truncated matrix; it stays exact only when held inside `BosonicMatrixElement` with integer indices.
- `DisplacementOperator[α, Infinity]` and `SqueezeOperator[ξ, Infinity]` are written in `\[FormalA]1` (from `FieldVariables[{1}]`), while `FieldVariables[]` returns `\[FormalA]`. QF treats every distinct formal symbol as its own mode, so mixing the two silently makes a two-mode expression. Build everything from `FieldVariables[{1}]`.

## QF problems found (exact layer)

- Fixed in commit fbebae87: `BosonicMatrixElement[{0, 0}, NormalOrdered[Exp[-a† ** a]]]` returned `Indeterminate` (a 0⁰ in `numberDiagonal`); it now returns 1 at the vacuum and 0 elsewhere.
- Open, minor: `BosonicMatrixElement[{0, 0}, 0^(-a† ** a)]` and `{1, 1}` of `0^(a† ** a - 1)` still give `Indeterminate`, because the fix covers a zero base only when the exponent is provably nonnegative.
- Open: `BosonicNormalOrder[{e1, e2}, vars]` with a list as first argument returns a wrong value with no message (for example `{{0 &}, a† ** a}` for `{a ** a†, a† ** a}`). The usage documents a single expression; mapping over the list gives the right answer.
- Missing rules (the call returns unevaluated or unchanged, never a wrong value):
  - `BosonicMatrixElement` in the coherent basis does not take exponentials of a or a† (`BosonicAlgebra.m:406`); with several modes it does not take exponentials in the Fock basis; with symbolic m it does not reduce {m, m} of D or {2m, 0} of S.
  - `BosonicExpOrder` leaves exp(x a² + y a†² + w a†a), exp(c a†a + linear), a general quadratic with linear terms, the Kerr exponential exp(λ (a†a)²), and exp(−(w + μ) a†a) with a sum as coefficient unchanged.
  - `BosonicExpSimplify` leaves conjugation of a sum of monomials, of a + f by a squeeze, by a mixed quadratic and linear generator, and two-mode conjugation by a general complex-symmetric squeeze unevaluated; it leaves `S(ξ) ** S(−ξ)` as `E^(X/2 − X/2)` uncombined.
  - `BosonicVEV` does not take a product of exponentials.

These gaps, and minimal ways to close them, are worked out in [Minimal-Bosonic-Design.md](Minimal-Bosonic-Design.md).

## Atlas issues found

No formula in the atlas has a wrong sign, factor or condition beyond the source errors it already flags. The corrections the atlas makes were rechecked exactly, and in each case the printed form fails while the corrected one holds: B-036, B-071, B-075, B-077, B-081, P-019, P-029, P2-013, P2-017, P2-023, P2-038, P2-054, C-050, C-085, G-028, G-077, G-079. These statements are ambiguous or loose and should be tightened:

- B-099 says "the first derivative is the right derivative". Read as the derivative written first (in a), it gives the wrong fermionic contraction sign and makes B-101 and B-102 fail; read as the right derivative in a* and the left one in a, both hold. The text should name the variable, not the position.
- X-006 refers to a matrix U that it never defines; it means the canonical block u of X-009.
- X-007 calls 0 ≤ M < I equivalent to ‖M‖ < 1. That holds for positive M; for a general M the trace converges when the spectral radius is below 1, and ‖M‖ < 1 is only sufficient.
- P-045, P-039 and P2-012 are right only under an index or basis reading the atlas leaves implicit: P-045's A is ⟨c_j c_i⟩, P-039's B is the transpose of the natural reading, and P2-012's K₂ pairs with a K₁ of flipped sign (the relabeling a → ia).
- G-050's fringe sign depends on an interferometer layout the atlas does not reproduce; it stays unchecked.

Every cell and check is stored in [drafts/exact](drafts/exact) (one `<section>-work` folder per book section, with the cells as run and their raw output).
