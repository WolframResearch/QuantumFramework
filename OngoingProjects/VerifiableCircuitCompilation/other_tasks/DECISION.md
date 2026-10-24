# Reviewed feasibility decision

Completed 2026-10-24. Advance the GF(8) multiplier to a bounded DQI decoder
integration experiment. Retain orthogonal pair mixing as a verified, modest
DFTHC improvement. No production decoder integration, commit, push or new
submission packaging was performed in this stage.

## Measured baselines and candidate

| Task / design | Official cost, lower is better | Evidence |
|---|---:|---|
| DQI unchanged baseline | 24,439,458,623,294.45 | Full PASS, 102,400 shots, zero value/phase errors |
| DFTHC unchanged baseline | 3,355,599,588.3155365 | Full PASS, 254,464 shots |
| DFTHC pair mixing | 3,313,336,907.312215 | Full PASS, 254,464 shots; 1.012755x improvement |

Both full baselines and the mixed candidate were independently rerun by the
correctness reviewer. Every deterministic report field reproduced; durations
are observations rather than speed benchmarks.

## DQI: choose the GF(8) route for integration

The exact quartic interpolation construction over GF(8) lowers a native
degree-12 multiplication from 63 to 42 Toffolis, controlled multiplication
from 75 to 54, and inverse phase from 315 to 210. Its fragment peak grows by
10 qubits. It dominates the tested GF(16) route on qubits, Toffolis and total
operations. The GF(64) route offers no Toffoli saving and adds more scratch.

Native field multiplications account for 64.05% of static DQI Toffolis.
Replacing their static counts and adding 10 qubits to the baseline peak gives
a **conditional** cost of 19,426,257,704,104.727, or **1.258063x** improvement.
This assumes unchanged repetition penalty and the modeled scratch peak.
There is no modified full-decoder result yet.

The next integration experiment should preserve the native bit basis and
classical phase masks, then measure the actual global peak and executed
Toffoli count. Acceptance requires the unchanged full official verifier,
including value and phase errors and its confidence penalty, and an actual
whole-task cost improvement. Fragment savings alone do not meet that gate.

## DFTHC: keep pair mixing; reject the spectral control

Rotating weights and identity offsets together preserves the represented
Hamiltonian and SOS energy. The tested pair search reduces actual official
cost by **1.259467%**, with the original angle words, dimensions, shifts and
circuit resources. The spectral control increases cost and is rejected.
The bounded coordinate search has no global-optimality certificate.

This is a small gain. The profile identifies a larger remaining circuit
target: angle loading contributes 79,480 of 93,237 Toffolis. Total workspace
is instead dominated by preparation's 544 temporary qubits. A future study
must price both effects in the official full objective.

## Aggregate implication and review gates

Using the existing verified inverse-square-root factor 3.9385957568764147,
the verified mixing factor and baseline DQI gives a geometric-mean factor
of **1.5859225744**. Adding the modeled DQI factor would give **1.7120487472**.
These are arithmetic projections from the recorded factors, not a new Kaggle
submission score. The latter depends on successful decoder integration;
these changes alone do not approach the requested 8.43991 target.

- [Correctness round 1](reviews/correctness-1.md): zero open issues. Real-kernel
  certificates, all fragments, both full baselines and full mixed candidate
  reproduced, plus independent polynomial arithmetic.
- [Quality round 1](reviews/quality-1.md): axes A-G pass, zero open issues.
  Current APIs, functional WL, symbolic derivations, nontrivial invariants
  and all five regimes checked. An independent noncommuting Hermitian SOS
  example preserved the shifted operator exactly.

The local suite passes 25 WL checks and 8 Python tests. See [README](README.md)
for reproduction, assumptions and full comparisons, and
[evidence-manifest.json](evidence-manifest.json) for the final file hashes.
Competition inputs remain external; the generated candidate is ignored.
