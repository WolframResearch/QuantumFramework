# Stabilizer Codes in the Wolfram Language, Cross-Checked Against Stim

**Mads Bahrami**, Wolfram Research Inc, USA

## Abstract

We verify the stabilizer-code layer that ships in the Wolfram QuantumFramework paclet against [Stim](https://github.com/quantumlib/Stim), the field-standard stabilizer simulator, and add an exact, symbolic code analysis built on top of it in the Wolfram Language. The work runs in two parts. Part A takes the shipped `PauliStabilizer` engine on the standard small codes (the three-qubit bit-flip and phase-flip codes, the nine-qubit Shor code, the five-qubit perfect code, and the seven-qubit Steane code) and checks, against Stim, that its stabilizer generators form a consistent group and that its syndromes agree across every one- and two-qubit error, and checks against the exact algebraic amplitudes that each codeword carries the correct sign; it also exercises a symbolic measurement that Stim structurally cannot represent. Part B builds a from-scratch GF(2)/symplectic analysis that returns the exact code distance and a minimum-weight witness, the logical operators with their full commutation algebra, and then carries a symbolic parameter: it certifies the distance of a whole code family indexed by $n$ in a single evaluation, and returns a coherent error's syndrome as an exact function of a continuous angle. The value is exactness and symbolic generality, not scale or speed; exact minimum distance is NP-hard, so this is a small-code instrument by design, complementary to Stim rather than competing with it.

## 1. Positioning

Stim is the reference tool for stabilizer circuits: it samples large circuits and detector-error models at enormous scale, which is what decoder benchmarking and threshold estimation need. It does not return the exact minimum distance of a code, and for a structural reason: computing the exact minimum distance of a stabilizer code is NP-hard (Kapshikar and Kundu, arXiv:2203.04262). Any tool that returns the *exact* answer is therefore confined to small codes. The textbook codes are small, and small is exactly where an exact, symbolic treatment is most useful: for teaching, for certifying a construction, and for stating results about whole families rather than single instances.

Two claims are separated deliberately. The first is that QuantumFramework's shipped stabilizer engine, `PauliStabilizer` (an Aaronson-Gottesman tableau), is *correct* on the standard codes, checked against an independent engine. The second is that the Wolfram Language around it can do things a numeric stabilizer simulator cannot: carry a measurement outcome or a code parameter as a symbol. Nothing here is asserted without an independent check, and a companion harness (`verify.wls`) re-runs every check below and exits nonzero on any failure.

## 2. Method

**Representation.** Each Pauli string on $n$ qubits maps to a binary symplectic vector $(x \mid z)$ of length $2n$: an $X$ sets an $x$-bit, a $Z$ a $z$-bit, a $Y$ both. Two Paulis commute iff the symplectic product $\omega(u, v) = x_u \cdot z_v + z_u \cdot x_v$ vanishes mod 2. A code is a set of $m$ independent commuting parity checks; it encodes $k = n - m$ logical qubits.

**Part A, the shipped engine against Stim.** The five codes are `PauliStabilizer` states (the five-qubit, Steane, and Shor codes are named states; the two three-qubit codes are built from their generators). Three checks run per code. First, sign-exactness: the shipped state has exact algebraic amplitudes, and applying each generator (with its sign) to the state vector returns it unchanged under a symbolic zero test on $G\lvert\psi\rangle - \lvert\psi\rangle$, so a wrong $-1$ cannot hide. Second, Stim accepts the same generators through `stim.Tableau.from_stabilizers`, which succeeds only for an independent, mutually commuting set and returns Stim's own canonical form of the group. Third, the syndrome differential: the engine's syndrome, read as $\tfrac12(1 - \langle g\rangle)$ for each check $g$ on the errored codeword, is compared to Stim's commutation-based syndrome over every weight-one and weight-two error. (Sign-exactness is established by the exact dense test, not by Stim; `from_stabilizers` checks group validity, and the syndrome comparison checks commutation, neither of which is sign-sensitive.)

**Part A, a symbolic measurement.** `PauliStabilizer`'s `"SymbolicMeasure"` carries a random outcome as a fresh formal bit $s_k$ written into a stabilizer sign. Substituting a value reproduces each post-measurement branch exactly, and when a later outcome is forced by the ones already taken, it is returned as an explicit function of the free bits. On a three-qubit parity circuit the third $Z$-outcome comes back as the sum of the first two, which is the syndrome of a single $Z$-check written as a polynomial in the measured bits. No numeric stabilizer sampler can hold that symbol.

**Part B, an independent exact analysis.** We reimplement the symplectic layer from scratch and let it judge the codes, without calling any packaged distance routine to check itself.

- *Distance.* A Pauli is a logical operator iff it commutes with every check (zero syndrome) *and* its symplectic vector is not in the GF(2) span of the checks. The span test is the one to get right: membership means appending the vector leaves the rank unchanged, so "not in the group" is the rank strictly *rising*. The distance is the least weight at which such a logical exists, found by an exhaustive weight-increasing search.
- *Logical algebra.* The logical pair $\bar X$, $\bar Z$ is checked to satisfy the full algebra: $\bar X$ anticommutes with $\bar Z$, both commute with every check, and neither lies in the group.
- *Witness against Stim.* Each distance witness is confirmed by Stim to be both undetectable (it commutes with every check) and nontrivial (it anticommutes with a logical), which places it in the normalizer but outside the stabilizer group. Stim confirms the witness is a genuine logical; its *minimality*, that nothing lighter is logical, is what the exhaustive search establishes and is the part no sampler provides.

## 3. Results

All checks pass, with Stim available (version 1.16.0 at time of writing). For the five codes:

| code | $[[n,k,d]]$ | witness | independent $d$ | Stim syndromes | sign-exact | Stim witness |
|---|---|---|---|---|---|---|
| bit-flip | $[[3,1,1]]$ | `ZII` | 1 | agree | +1 | logical |
| phase-flip | $[[3,1,1]]$ | `XII` | 1 | agree | +1 | logical |
| Shor | $[[9,1,3]]$ | `XXXIIIIII` | 3 | agree | +1 | logical |
| five-qubit | $[[5,1,3]]$ | `XYXII` | 3 | agree | +1 | logical |
| Steane | $[[7,1,3]]$ | `XXXIIII` | 3 | agree | +1 | logical |

The shipped engine's syndrome equals Stim's on the entire weight-one and weight-two error set (36, 36, 351, 105, and 210 errors respectively); Stim accepts every generator set as a consistent stabilizer group; each codeword is a sign-exact $+1$ state by the dense exact test; the independent distance is found without reference to any packaged routine; the logical operators pass the full algebra check; and Stim confirms each witness is a genuine logical of the stated weight.

One result is worth stating plainly because it is often misremembered: the three-qubit bit-flip and phase-flip codes have quantum distance one, not three. A single $Z$ on one qubit of the bit-flip code commutes with both $Z$-type checks (zero syndrome) yet is not in the stabilizer group, so it is an undetectable weight-one logical. The classical repetition code has distance three against bit flips; as a quantum code it leaves phase errors unguarded. The analysis reports the mathematics, not the folklore.

**Symbolic generality.** The symplectic structure lets a symbolic parameter be carried through the analysis. For the $[[n, n-2, 2]]$ family (generators $X^{\otimes n}$ and $Z^{\otimes n}$), the generator commutator is $\omega(X^{\otimes n}, Z^{\otimes n}) = n$, so the code is defined exactly when $n$ is even; every weight-one error overlaps a generator on one qubit and is detected; and a weight-two logical is lighter than every nontrivial (weight-$n$) group element as soon as $n > 2$. Together these certify distance exactly two for the entire even-$n$ family in a single evaluation, with $n$ symbolic, and the from-scratch search confirms distance two at $n = 4, 6, 8, 10$. The repetition family $[[n,1,1]]$ has distance one for all $n$ (a lone $Z_1$ is undetectable and, by a parity argument on the $Z$-support, never in the group). And for a coherent (non-Clifford) error $R_x(\theta)$ on the bit-flip codeword, the syndrome expectation $\langle Z_1 Z_2 \rangle = \cos\theta$ is returned as an exact function of the continuous angle. None of these family-level or continuous-parameter statements is expressible in a single object by a numeric stabilizer sampler, which must instantiate each $n$ and cannot represent $\theta$ symbolically at all.

## 4. Limits

This is deliberately a small-code, exact instrument. The distance routine is an exhaustive weight-increasing search, feasible only because $n$ is small; exact minimum distance is NP-hard, so no reorganization of the code changes that ceiling. It does not sample noisy circuits, build detector-error models, decode at scale, or estimate thresholds; those are Stim's domain, and Stim is far faster and larger there. The two tools are complementary: reach for this layer when you want the exact answer for a textbook code, a sign-checked stabilizer state, or a statement about a whole family, and for Stim when you want a threshold for a distance-25 surface code. The point of the cross-check is that, wherever the two overlap, they agree, which is what makes the exact answers trustworthy and the symbolic ones worth having.

## 5. Reproduce

The full computational essay (`exact-stabilizer-codes-cross-checked.md`, built to a notebook by MarkdownToNotebook) carries every cell above, and `verify.wls` re-runs every check and prints PASS/FAIL. To run either you need a QuantumFramework checkout and, for the Stim differential, `stim` installed in the Python that `ExternalEvaluate` uses (`pip install stim`; without it the sign-exact, distance, logical-algebra, and symbolic-family checks still run). See the accompanying `README.md` for exact commands.

## References

- D. Gottesman, *Stabilizer Codes and Quantum Error Correction*, Ph.D. thesis, Caltech (1997), arXiv:quant-ph/9705052.
- S. Aaronson and D. Gottesman, *Improved simulation of stabilizer circuits*, Phys. Rev. A **70**, 052328 (2004), arXiv:quant-ph/0406196.
- U. Kapshikar and S. Kundu, *On the hardness of the minimum distance problem of quantum codes*, arXiv:2203.04262 (2022).
- C. Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum **5**, 497 (2021), arXiv:2103.02202.
