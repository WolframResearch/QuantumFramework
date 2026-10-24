# DQI and DFTHC feasibility graph

Scope: baseline profiling and reviewed feasibility only. No production decoder,
new competition submission, commit or push is authorized by this stage.
The verified inverse-square-root incumbent remains fixed.

## Nodes and signals

- INIT: pin downloaded verifier/reference hashes, read implementation and WL
  documentation. Signal `CONTRACTS_PINNED`.
- BASE: run unchanged baseline verifiers; collect actual objectives and native
  stage profiles. Signal `BASELINES_PROFILED`, or explicit incomplete coverage.
- DQI: prove exact subfield multiplication tensors and mask/basis conventions;
  compile bounded arithmetic fragments including cleanup and scratch.
  Signal `IDENTITY_AND_RESOURCE_PASS` or `REJECTED` for each design.
- DFTHC: prove orthogonal ensemble mixing invariants, derive the piecewise
  pair-rotation objective and screen fixed-angle component mixing.
  Signal `INVARIANTS_AND_METRICS_PASS` or `REJECTED`.
- V: fresh wl-verify reviewer, exit `OPEN ISSUES: n`.
- Q: fresh wl-quality reviewer, exit `OPEN ISSUES: n`.
- DECIDE: compare full objectives and state a next implementation recommendation.
  Signal `RECOMMENDATION_ONLY`.

INIT precedes BASE, DQI and DFTHC. Both feasibility branches and BASE precede
V, then Q, then DECIDE. A failed identity rejects that design. A resource or
metric failure permits one redesign. Correctness findings return to the
responsible source. Quality fixes affecting executable WL return to V.
Each review has at most three fresh rounds; correctness has at most three
visits. This is a code deliverable; no notebook build is required.

## Concrete seeds

The downloaded source is external and is never vendored.

- `dqi/dqi_reference_circuit.py:522`: Chien/Forney uses field inversions and
  phase-by-product at every codeword position.
- `dqi/dqi_reference_circuit.py:613`: default EEA iteration count is
  `3*(n+1)+9`; benchmark n=70.
- Pinned kickmix `gen_gf_mul.cc`: recursive Karatsuba is already implemented.
  `gen_gf_inverse.cc`: the explicit Itoh-Tsujii scaffold already permits
  zero-Toffoli cleanup. Neither is a new improvement.
- `dfthcsos/verify.py:389,419,472`: whole-block costs share temporary workspace;
  the objective is not just the submitted ANGLE circuit count.
- `dfthcsos/scoreinfra/dfthc_math.py:79,298`: tensor reconstruction and effective
  normalization depend on factorization weights, enabling a test of exact
  orthogonal component mixing while retaining the angle-word table.

Initial DQI comparison: native degree-12 Karatsuba; a quadratic extension of
GF(64); cubic interpolation over GF(16); quartic interpolation over GF(8).
All converted linear maps and nonlinear products must be charged in compiled
fragments. An optimized GF(8) six-product formula is a hypothesis until proved.

Initial DFTHC comparison: unchanged components; orthogonal pair coordinate
mixing; a bounded whole-component spectral rotation as a control. All use the
same factorization dimensions and decoded angle vectors. Changes are checked
against the unmodified official metrics; fixed symmetry shifts are retained.
