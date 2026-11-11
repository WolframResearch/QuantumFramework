# Cross-checking the circuit-level machinery against Stim

The `.wlt` suite must run on a bare Wolfram install, so it cannot call Stim.
This directory holds the cross-check that does, and it is the independent
oracle behind `Tests/CircuitNoise.wlt`.

```
pip install stim pymatching
wolframscript -file emit.wls
python3 check.py
```

`emit.wls` writes, per case, the memory experiment in Stim syntax and a JSON of
this package's **exact** detector firing rates and observable flip rates.
`check.py` samples the same circuit with Stim and compares. The comparison is
closed form against Monte Carlo, so a genuine disagreement shows up as a large
sigma and not as noise.

It also hands each exported circuit to PyMatching. That is a second, structural
check: a malformed circuit does not produce a detector error model that
decomposes into a matching graph.

## Result of the run on 2026-08-31

Five cases -- three at circuit level, two at phenomenological -- two million
shots each. `max sigma` is over every detector in the case.

| case | detectors | max sigma | observable (ours, exact) | observable (Stim) | sigma | MWPM |
|---|---|---|---|---|---|---|
| `bf-circ` | 6 | 2.25 | 0.00532 | 0.00524 | 1.61 | 0.00013 |
| `bf-phen` | 6 | 1.01 | 0.01324 | 0.01323 | 0.22 | 0.00055 |
| `shor-phen` | 16 | 2.26 | 0.03894 | 0.03912 | 1.29 | 0.00159 |
| `steane-circ` | 18 | 1.83 | 0.04094 | 0.04132 | 2.65 | 0.03840 |
| `5q-circ` | 12 | 1.94 | 0.07127 | 0.07096 | 1.75 | 0.02767 |

58 detector rates and 5 observable rates, all within 2.7 sigma. PyMatching
decoded every exported circuit unchanged.

`5q-circ` is the one that matters most for the emitter: the 5-qubit code is not
CSS, its generators carry `Y`, and so it is the only case that exercises the
sqrt(X) basis change and its round trip through Stim's `SQRT_X`.

## Result of the run on 2026-11-03: the scalable families

Two more cases, added with `QECCode["Surface", d]` and `QECCode["Toric", L]`, two million shots
each (Stim 1.16.0, PyMatching 2.4.0).

| case | detectors | max sigma | observable (ours, exact) | observable (Stim) | sigma | MWPM |
|---|---|---|---|---|---|---|
| `surface3-circ` | 32 | 2.49 | 0.06574 | 0.06549 | 1.42 | 0.00443 |
| `toric3-phen` | 48 | 2.23 | 0.03869 | 0.03874 | 0.39 | 0.00129 |

80 detector rates and 2 observable rates, all within 2.5 sigma; PyMatching decodes both. The
`surface3-circ` MWPM column is the first number in this table that a native matching decoder in the
layer will have to reproduce.

## Result of the run on 2026-11-04: the surface-code check order and the native matching decoder

`surface3-circ` re-emitted after `QECCode["Surface", d]` started carrying its CNOT order (`code["CheckOrder"]`: X checks
row by row, Z checks column by column, so no hook error lies along a logical operator), two million
shots:

| case | detectors | max sigma | observable (ours, exact) | observable (Stim) | sigma | MWPM |
|---|---|---|---|---|---|---|
| `surface3-circ` | 32 | 1.61 | 0.06574 | 0.06582 | 0.43 | 0.00443 |

The `Z-bar` statistics do not move with the Z-check order -- that order only places Z errors -- which
is what the unchanged exact rate shows. The layer's own `QECDecoder[..., "Method" -> "Matching"]`,
on the graph it builds from its own detector model, fails on `Z-bar` in 0.439 % of 200 000 shots of
the same circuit, against PyMatching's 0.443 %. On Stim's own decomposed graph the two decoders
agree shot for shot (100 000 shots, 100 %). The native decoder also reports `X-bar`, which this
single-basis experiment cannot: 0.40 %, and its scaling from p = 1/500 to 1/2000 (a factor 13)
is quadratic. With the old row-by-row Z checks it was 0.98 % and fell by only a factor 5 -- the
linear term of a hook error along logical `Z`.

## What this does and does not establish

It establishes that the emitted circuit, the fault enumeration, the frame
propagation and the detector definitions agree with Stim's, which is the part
that is easy to get subtly wrong.

It does not compare logical error rates directly, and deliberately. Stim's
memory experiment is single-basis by construction -- prepare `|0_L>`, keep it,
measure `Z-bar` -- while `QECLogicalErrorRate` asks whether the decoder got the
residual's whole class right, covering logical `X` and logical `Z` damage at
once. The two are different questions, and the `MWPM` column is answering
Stim's. Comparing the built-in minimum-weight-over-faults decoder against
PyMatching on equal terms belongs with the scalable families of roadmap item 4.
