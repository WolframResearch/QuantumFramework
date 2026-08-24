# Exact, Symbolic Stabilizer-Code Analysis, Cross-Checked Against Stim

**Mads Bahrami**, Wolfram Research Inc, USA

## Abstract

We give an exact, symbolic layer for stabilizer-code analysis in the Wolfram Language and verify it two independent ways. For the standard small codes (the three-qubit bit-flip and phase-flip codes, the nine-qubit Shor code, the five-qubit perfect code, and the seven-qubit Steane code) the layer returns the exact code distance and a minimum-weight witness, the logical Pauli operators with their full commutation algebra, and a sign-exact +1 Clifford encoder. Each result is confirmed by a from-scratch GF(2)/symplectic re-derivation and by a differential comparison against Stim, the field-standard stabilizer simulator: the syndromes agree across every weight-one and weight-two error, and Stim's own tableau confirms every encoder to +1. The layer then does one thing a numeric simulator structurally cannot: it carries a symbolic parameter, certifying the distance of an entire code family indexed by $n$ in a single evaluation, and returning a coherent error's syndrome as an exact function of a continuous angle. The value is exactness and symbolic generality, not scale or speed; exact minimum distance is NP-hard, so this is a small-code instrument by design, complementary to Stim rather than competing with it.

## 1. Positioning

Stim is the reference tool for stabilizer circuits: it samples large circuits and detector-error models at enormous scale, which is what decoder benchmarking and threshold estimation need. It does not return the exact minimum distance of a code, and for a structural reason: computing the exact minimum distance of a stabilizer code is NP-hard (Kapshikar and Kundu, arXiv:2203.04262). Any tool that returns the *exact* answer is therefore confined to small codes. The textbook codes are small, and small is exactly where an exact, symbolic treatment is most useful: for teaching, for certifying a construction, and for stating results about whole families rather than single instances.

This note describes such a layer (the `Wolfram`QuantumFramework`QEC`` package, following Gottesman's thesis) and, more importantly, the discipline used to verify it. Nothing here is asserted without an independent check.

## 2. Method

**Representation.** Each Pauli string on $n$ qubits maps to a binary symplectic vector $(x \mid z)$ of length $2n$: an $X$ sets an $x$-bit, a $Z$ a $z$-bit, a $Y$ both. Two Paulis commute iff the symplectic product $\omega(u, v) = x_u \cdot z_v + z_u \cdot x_v$ vanishes mod 2. A code is a set of $m$ independent commuting generators; it encodes $k = n - m$ logical qubits.

**Independent re-derivation (in-language).** We reimplement the symplectic layer from scratch and let it judge the package, without calling the package's own analysis routines to check themselves.

- *Distance.* A Pauli is a logical operator iff it commutes with every generator (zero syndrome) *and* its symplectic vector is not in the GF(2) span of the generators. The span test is the one to get right: membership means appending the vector leaves the rank unchanged, so "not in the group" is the rank strictly *rising*. The distance is the least weight at which such a logical exists, found by an exhaustive weight-increasing search.
- *Logical algebra.* The package's $\bar X_i$, $\bar Z_i$ are checked to satisfy the full algebra: $\bar X_i$ anticommutes with $\bar Z_i$, commutes with every other logical, all $\bar X$ mutually commute, all $\bar Z$ mutually commute, and every logical commutes with every stabilizer while lying outside the group.
- *Encoder.* The real encoding circuit is applied to $\lvert 0\ldots 0\rangle$ and, using independently built Pauli matrices, each generator is checked to fix the resulting state vector *exactly* (a symbolic zero test on $G\lvert\psi\rangle - \lvert\psi\rangle$, not a fidelity), so a wrong $-1$ sign cannot hide.

**Stim differential (out-of-language).** A small converter maps a code to Stim Pauli strings and its encoder to a Stim circuit. Three differentials are run: the package's `Syndrome` is compared to Stim's commutation-based syndrome over all weight-one and weight-two errors; Stim's `peek_observable_expectation` is required to be +1 for every generator on the encoded state (a sign-exact encoder check through an independent tableau engine); and each distance witness is confirmed by Stim to be both undetectable (commutes with all generators) and nontrivial (anticommutes with a logical), which places it in the normalizer but outside the stabilizer group.

## 3. Results

All checks pass. For the five codes:

| code | $[[n,k,d]]$ | witness | independent $d$ | Stim syndromes | Stim encoder |
|---|---|---|---|---|---|
| bit-flip | $[[3,1,1]]$ | `ZII` | 1 | agree | +1 |
| phase-flip | $[[3,1,1]]$ | `XII` | 1 | agree | +1 |
| Shor | $[[9,1,3]]$ | `XXXIIIIII` | 3 | agree | +1 |
| five-qubit | $[[5,1,3]]$ | `XYXII` | 3 | agree | +1 |
| Steane | $[[7,1,3]]$ | `XXXIIII` | 3 | agree | +1 |

The independent distance equals the package distance for every code; the logical operators pass the full algebra check; the encoders are sign-exact by both the dense-matrix and the Stim-tableau route; and the package syndrome equals Stim's on the entire weight-one and weight-two error set (36, 36, 351, 105, and 210 errors respectively).

One result is worth stating plainly because it is often misremembered: the three-qubit bit-flip and phase-flip codes have quantum distance one, not three. A single $Z$ on one qubit of the bit-flip code commutes with both $Z$-type stabilizers (zero syndrome) yet is not in the stabilizer group, so it is an undetectable weight-one logical. The classical repetition code has distance three against bit flips; as a quantum code it leaves phase errors unguarded. The exact layer reports the mathematics, not the folklore.

**Symbolic generality.** The symplectic structure lets a symbolic parameter be carried through the analysis. For the $[[n, n-2, 2]]$ family (generators $X^{\otimes n}$ and $Z^{\otimes n}$), the generator commutator is $\omega(X^{\otimes n}, Z^{\otimes n}) = n$, so the code is defined exactly when $n$ is even; every weight-one error overlaps a generator on one qubit and is detected; and a weight-two logical is lighter than every nontrivial (weight-$n$) group element as soon as $n > 2$. Together these certify distance exactly two for the entire even-$n$ family in a single evaluation, with $n$ symbolic. Likewise the repetition family $[[n,1,1]]$ has distance one for all $n$ (a lone $Z_1$ is undetectable and, by a parity argument on the $Z$-support, never in the group). And for a coherent (non-Clifford) error $R_x(\theta)$ on the bit-flip codeword, the syndrome expectation $\langle Z_1 Z_2 \rangle = \cos\theta$ is returned as an exact function of the continuous angle. None of these family-level or continuous-parameter statements is expressible in a single object by a numeric stabilizer sampler, which must instantiate each $n$ and cannot represent $\theta$ symbolically at all.

## 4. Limits

This is deliberately a small-code, exact instrument. The distance routine is an exhaustive weight-increasing search, feasible only because $n$ is small; exact minimum distance is NP-hard, so no reorganization of the code changes that ceiling. It does not sample noisy circuits, build detector-error models, decode at scale, or estimate thresholds; those are Stim's domain and Stim is far faster and larger there. The two tools are complementary: reach for this layer when you want the exact answer for a textbook code or a statement about a whole family, and for Stim when you want a threshold for a distance-25 surface code. The point of the cross-check is that, wherever the two overlap, they agree, which is what makes the exact answers trustworthy and the symbolic ones worth having.

## 5. Reproduce

The full computational essay (`exact-stabilizer-codes-cross-checked.md`, built to a notebook by MarkdownToNotebook) carries every cell above. To run it you need a QuantumFramework checkout and, for the Stim differential, `stim` installed in the Python that `ExternalEvaluate` uses (`pip install stim`; without it the independent re-derivation still runs). See the accompanying `README.md` for exact commands.

## References

- D. Gottesman, *Stabilizer Codes and Quantum Error Correction*, Ph.D. thesis, Caltech (1997), arXiv:quant-ph/9705052.
- U. Kapshikar and S. Kundu, *On the hardness of the minimum distance problem of quantum codes*, arXiv:2203.04262 (2022).
- C. Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum **5**, 497 (2021), arXiv:2103.02202.
