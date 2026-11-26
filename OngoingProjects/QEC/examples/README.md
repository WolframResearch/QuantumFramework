# Examples

Two studies done with the QEC layer alone, evaluated, with their outputs saved. Open them in
Mathematica from this folder: the first input cell finds the repository from the notebook's
location and loads the paclet and the QEC core.

| notebook | what it is | evaluation time |
|---|---|---|
| `SurfaceCodeCircuitLevel.nb` | The distance-3 surface code under CNOT and readout errors: the code, its extraction circuit, the 352 faults, and the exact logical error rate under two decoders (every single fault, and maximum likelihood) as polynomials in pC and pM. Then more rounds, distance 5 exactly, and sampled memory experiments up to distance 7. | about 50 s |
| `ToricSquidSymbolic.nb` | The review of the driven-SQUID toric-code proposal (arXiv:2611.03329), complete. The code part (Question 4 and the appendix) runs on the current QEC layer. | about 4.5 min, 3 of them in one NDSolve |

Both redo, with the current layer, the two notebooks of November 2026 that started the speed work
(`surface-code.nb` and `toric-squid-symbolic.nb`).

- **Surface code.** The original built the faults with a stabilizer simulation per fault and wrote
  the Walsh–Hadamard sum by hand for this one circuit. Here the same numbers come from three calls:
  `QECDetectorModel` (0.01 s), and `QECLogicalErrorRate` with the default decoder (0.16 s) and with
  maximum likelihood (0.22 s). The single-fault polynomial is the original's, symbol for symbol.
  Maximum likelihood replaces the decoder that knew every pair of faults. It has the same leading
  terms, and its rate is lower: 5.41e-5 against 5.55e-5 at p = 1/1000. Because the same calls take
  any code and any number of rounds, the notebook then goes where the original could not.
- **Toric code.** Every result is the original's; each output was compared with a run of the
  original. Two things changed:
  - `QECCode["CSS"[hx, hz]]` is the current spelling of the CSS construction.
  - Pauli products become matrices through `QECPauli[...]["Matrix"]`, which builds the sparse matrix
    directly, about 100 times faster than through `QuantumOperator`. The second-order sums of
    Question 4 drop from 43 s to 1.5 s.

  The rest of the time is the physics: the Floquet check (NDSolve, 177 s), a `Reduce` over the
  tone phases (37 s) and an `Asymptotic` (15 s).
