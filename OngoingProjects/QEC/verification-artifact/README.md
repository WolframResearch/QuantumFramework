# Stabilizer codes, cross-checked against Stim

A public, reproducible verification artifact for the stabilizer-code layer that
ships in QuantumFramework (`PauliStabilizer`), plus an exact, symbolic code
analysis built on top of it in the Wolfram Language. It shows that the shipped
stabilizer engine agrees with the field-standard simulator
[Stim](https://github.com/quantumlib/Stim) on the standard textbook codes, and
that the Wolfram Language around it does one thing Stim structurally cannot:
carry a symbolic parameter (a whole code family, or a continuous error angle) and
a symbolic measurement outcome.

The value is **exactness and symbolic generality, not scale or speed**. Exact
minimum distance is NP-hard (Kapshikar and Kundu, arXiv:2203.04262), so the
distance analysis is a small-code instrument by mathematical necessity,
complementary to Stim rather than competing with it. Stim owns scale, sampling,
and decoding; this layer owns the exact, small, and symbolic corner.

## Contents

- `exact-stabilizer-codes-cross-checked.md` - the computational essay, authored as
  literate markdown (`Template: Default`) for MarkdownToNotebook (md2nb). **Source of truth.**
- `exact-stabilizer-codes-cross-checked.nb` - the built notebook, every cell
  evaluated and cached. **A build output; do not hand-edit** (regenerate with `build.wls`).
- `methods-note.md` - a short methods note (arXiv note / Wolfram Community / Example
  Repository style): the method, the agreement result, the symbolic-family result, the limits.
- `verify.wls` - a standalone harness that re-runs every check and prints PASS/FAIL
  (exits nonzero on any failure). The fastest way to reproduce the agreement.
- `build.wls` - rebuilds the notebook from the essay markdown.

## What is verified

The five codes are shipped `PauliStabilizer` states (the five-qubit, Steane, and
Shor codes are named; the three-qubit bit-flip `[[3,1,1]]` and phase-flip
`[[3,1,1]]` codes are built from their generators). The work runs in two parts.

**Part A: the shipped stabilizer engine, cross-checked against Stim.**

1. Each codeword is a **sign-exact +1 state**: applying every generator (sign
   included) to the exact algebraic state vector leaves it fixed, under a symbolic
   zero test, so a wrong sign cannot hide.
2. Stim **accepts every generator set** as a consistent stabilizer group
   (`stim.Tableau.from_stabilizers`).
3. The engine's **syndrome equals Stim's** across every weight-one and weight-two
   error, code by code (36, 36, 351, 105, and 210 errors for bit-flip, phase-flip,
   Shor, five-qubit, and Steane).
4. A **symbolic measurement** (`"SymbolicMeasure"`) carries a random outcome as a
   formal bit and returns a forced outcome as an explicit function of the free
   bits, which a numeric stabilizer sampler cannot represent.

**Part B: an in-language GF(2)/symplectic analysis, cross-checked against Stim.**
Written from scratch (it calls no packaged distance routine to check itself), it
reproduces the exact code **distance** and a minimum-weight **witness**, verifies
the full **logical commutation algebra**, and has **Stim confirm each witness** is
undetectable and nontrivial. It then certifies the `[[n, n-2, 2]]` and repetition
families at symbolic `n` in one evaluation, and returns a coherent error's
syndrome as the exact function `Cos[θ]`.

Note (correct, and worth stating): the three-qubit bit-flip and phase-flip codes
have **quantum distance one**, not three. A single `Z` (or `X`) is an undetectable
logical. The classical repetition code has distance three; as a quantum code it
leaves the conjugate error type unguarded.

The exactness/symbolic content is Wolfram Language built on the shipped stabilizer
primitives; the exact **distance** and **logical** operators are computed in the
essay, not returned by a packaged function. What Stim independently cross-checks is
the generator groups, the syndromes, and the witness validity; sign-exactness is the
exact dense test on the Wolfram side.

## Requirements

- A QuantumFramework checkout (this repository). The package is **not** installed as a
  paclet; it is loaded from the checkout. The scripts locate the checkout automatically
  by walking up from their own location; the essay's first cell uses an explicit path
  that you adjust to your own checkout.
- `wolframscript` on the PATH.
- For the Stim differential only: `stim` importable from the Python that
  `ExternalEvaluate` uses. Without it, `verify.wls` reports the Stim sections as SKIP
  (not FAIL) and the sign-exact, distance, logical-algebra, and symbolic-family checks
  still run.

## Reproduce

From this directory:

**1. (Stim differential only) install stim into the ExternalEvaluate Python.** Find that
interpreter, then install:

```bash
wolframscript -code 'ExternalEvaluate["Python", "import sys; sys.executable"]'
```

```bash
/path/to/that/python -m pip install stim
```

**2. Run the full verification harness** (sign-exact states + Stim group acceptance +
syndrome differential + symbolic measurement + in-language distance/logical/witness +
symbolic families):

```bash
wolframscript -file verify.wls
```

It prints a PASS/FAIL line per check and exits `0` only if all pass.

**3. Rebuild the notebook** from the essay markdown (evaluates and caches every cell):

```bash
wolframscript -file build.wls
```

MarkdownToNotebook (md2nb) is expected at `~/Documents/GitHub/MarkdownToNotebook`; edit
the `md2nb` path in `build.wls` if yours differs, or use the public cloud deployment
linked from that project.

## References

- D. Gottesman, *Stabilizer Codes and Quantum Error Correction*, Caltech Ph.D. thesis
  (1997), arXiv:quant-ph/9705052.
- S. Aaronson and D. Gottesman, *Improved simulation of stabilizer circuits*,
  Phys. Rev. A **70**, 052328 (2004), arXiv:quant-ph/0406196. The tableau engine
  `PauliStabilizer` follows.
- U. Kapshikar and S. Kundu, *On the hardness of the minimum distance problem of quantum
  codes*, arXiv:2203.04262 (2022). Exact quantum minimum distance is NP-hard.
- C. Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum **5**, 497 (2021),
  arXiv:2103.02202.
