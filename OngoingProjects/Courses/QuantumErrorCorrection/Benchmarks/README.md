# Benchmark Contract

This directory will contain reproducible data, configurations, and results for decoder, threshold,
logical-memory, and resource questions. No benchmark data are committed in the initial curriculum
draft.

## Required record for each benchmark

- Question identifier and benchmark name.
- Code family and literal instance, including check matrices or a generator script.
- Verified $[[n,k,d]]$ parameters and the convention for subsystem distance if applicable.
- Noise model, units, circuit locations, correlations, leakage/erasure semantics, and round count.
- Syndrome-extraction schedule and boundary conditions.
- Decoder name, version, commit, weights, stopping conditions, and hardware used for timing.
- Random seed or replayable event stream.
- Number of trials, logical failures, estimator, and uncertainty interval.
- Wall-clock and per-shot time, excluding and including preprocessing and model training as separate
  quantities.
- Exact or exhaustively enumerable small-instance anchor.
- Known limitations and the parameter range over which the benchmark is claimed.

## Directory pattern

```text
Benchmarks/
  <question-id>-<short-name>/
    README.md
    config.json
    code-instance/
    events/
    results/
```

Large generated event streams should not be committed blindly. Store a generator, checksums,
summary statistics, and a durable data location when the data exceed the repository's policy.

## Fair-comparison rules

1. Competing decoders receive the same events.
2. A decoder using exact generating-model probabilities is labelled oracle-weighted.
3. Training, tuning, and test samples are disjoint and their split is recorded.
4. Accuracy and latency are reported separately before any combined figure of merit.
5. Timeouts and decoder failures count under a declared policy; they are not discarded.
6. Threshold estimates use several distances and report finite-size drift.
7. Resource estimates include routing, syndrome extraction, decoding, reset, and magic-state cost
   when those resources are in scope.
