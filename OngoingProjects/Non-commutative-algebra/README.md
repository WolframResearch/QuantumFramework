# Zassenhaus Formula and Baker–Campbell–Hausdorff Expansion

**Some Applications of Non-Commutative Algebra**  
**Mathematica Version:** 14.3.0 (June 12, 2025)

---

## Overview

This folder contains Mathematica notebooks that provides implementation of two fundamental tools in non-commutative algebra:

1. **Zassenhaus Formula**  
   Decomposes \(e^{X_1+X_2+\cdots+X_n}\) into an (infinite) product of exponentials of nested commutators.
2. **Baker–Campbell–Hausdorff (BCH) Expansion**  
   Computes the logarithm \(\log(e^{X_1}e^{X_2}\cdots e^{X_n})\) as a formal series in nested commutators.

In addition to the core algorithms, the notebook includes:
- Symbolic tests verifying known identities from the literature.
- Support for multi-operator generalizations (\(n>2\)).
- Both “raw” expansion form and a more compact commutator-only form.

