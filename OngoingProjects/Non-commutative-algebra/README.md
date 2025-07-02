# Zassenhaus Formula and Baker–Campbell–Hausdorff Expansion

**Some Applications of Non-Commutative Algebra**  
**Mathematica Version:** 14.3.0 (June 12, 2025)

---

## Overview

This folder contains Mathematica notebooks that provide implementations of two fundamental tools in non-commutative algebra:

1. **Zassenhaus Formula**  
   Decomposes $e^{X_1 + X_2 + \cdots + X_n}$ into an (infinite) product of exponentials of nested commutators (i.e. the degree-n term of Zassenhaus Formula)

2. **Baker–Campbell–Hausdorff (BCH) Expansion**  
   Computes the logarithm $\log(e^{X_1}e^{X_2}...e^{X_n})$ as a formal series in nested commutators (i.e the degree-n term of BCH Expansion)

In addition to the core algorithms, the notebooks include:

- Symbolic tests verifying known identities from the literature.  
- Support for multi-operator generalizations ($n > 2$).  
- Both “raw” expansion form and a more compact commutator-only form.

The non-commutative operators $X_j$ can be matrices, arrays, operators or [symbolic arrays](https://reference.wolfram.com/language/guide/SymbolicArrays.html) introduced in V.14.1 of Mathematica; also the operation of non-commutative operators on each other, $X_i(X_j)$, can be [Dot](https://reference.wolfram.com/language/ref/Dot.html), [Composition](https://reference.wolfram.com/language/ref/Composition.html), [NonCommutativeMultiply](https://reference.wolfram.com/language/ref/NonCommutativeMultiply.html) or a customized operation.