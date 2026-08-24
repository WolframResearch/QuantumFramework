# Exact stabilizer codes, cross-checked against Stim

A public, reproducible verification artifact for the exact stabilizer-code analysis
layer in QuantumFramework (`Wolfram`QuantumFramework`QEC``). It shows that the
symbolic, exact code analysis agrees with the field-standard simulator
[Stim](https://github.com/quantumlib/Stim) on the standard textbook codes, and does
one thing Stim structurally cannot: carries a symbolic parameter and certifies a
whole family (or a continuous error angle) in a single evaluation.

The value is **exactness and symbolic generality, not scale or speed**. Exact
minimum distance is NP-hard (Kapshikar and Kundu, arXiv:2203.04262), so this is a
small-code instrument by mathematical necessity, complementary to Stim rather than
competing with it. Stim owns scale, sampling, and decoding; this layer owns the
exact, small, and symbolic corner.

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

For the bit-flip `[[3,1,1]]`, phase-flip `[[3,1,1]]`, Shor `[[9,1,3]]`, five-qubit
`[[5,1,3]]`, and Steane `[[7,1,3]]` codes, the exact code distance and witness, the
logical operators, and a sign-exact +1 encoder are each confirmed two independent ways:

1. **An in-language GF(2)/symplectic re-derivation** written from scratch, which does
   not call the package's own analysis routines to check themselves. It reproduces the
   distance, verifies the full logical commutation algebra, and confirms each encoder
   fixes every generator with an exact (sign-aware) state-vector test.
2. **A Stim differential.** A small converter maps each code to Stim, and the package's
   `Syndrome` is diffed against Stim's over every weight-one and weight-two error; Stim's
   tableau confirms each encoder to +1; and each distance witness is confirmed by Stim to
   be undetectable and nontrivial.

The symbolic differentiator then certifies the `[[n, n-2, 2]]` and repetition families
at symbolic `n` in one evaluation, and returns a coherent error's syndrome as the exact
function `Cos[θ]`.

Note (correct, and worth stating): the three-qubit bit-flip and phase-flip codes have
**quantum distance one**, not three. A single `Z` (or `X`) is an undetectable logical.
The classical repetition code has distance three; as a quantum code it leaves the
conjugate error type unguarded.

## Requirements

- A QuantumFramework checkout (this repository). The package is **not** installed as a
  paclet; it is loaded from `OngoingProjects/QEC/QEC/QEC.wl` in the checkout. The scripts
  locate the checkout automatically by walking up from their own location.
- `wolframscript` on the PATH.
- For the Stim differential only: `stim` importable from the Python that
  `ExternalEvaluate` uses. Without it, `verify.wls` reports the Stim section as SKIP (not
  FAIL) and the in-language re-derivation still runs.

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

**2. Run the full verification harness** (package tests + independent re-derivation +
Stim differential + symbolic families):

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
  (1997), arXiv:quant-ph/9705052. The constructions the QEC package follows.
- U. Kapshikar and S. Kundu, *On the hardness of the minimum distance problem of quantum
  codes*, arXiv:2203.04262 (2022). Exact quantum minimum distance is NP-hard.
- C. Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum **5**, 497 (2021),
  arXiv:2103.02202.
