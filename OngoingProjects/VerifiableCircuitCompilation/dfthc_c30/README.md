# Standalone fixed-factor C30 loader

`generator.py` builds the public (65,60,30) factorization's 25 ten-bit,
three-word XOR circuits. It depends on NumPy, kickmix and the official
competition verifier package, with no imports from downloaded public notebook
modules. Public W and H2_I stay exactly fixed. The official shifts and tree
stay fixed. The signed-angle gauge changes only row signs, with measured
component drift below 1e-11; row signs preserve the physical projectors.
See [NOTICE.md](NOTICE.md) for source attribution and license.

Each circuit implements `|a,z> -> |a,z XOR delta_k(a)>` on arbitrary target
words. Four low-address cases share a base and three nine-bit differences.
The redundant sign bit restores the tenth bit. A five-wire pebbled AND
traversal and split-ANF measurement phase repair keep the FeMoCo circuits
at Q128. Sequential deltas load consecutive batches; reverse traversal clears
them. This port retains the public factors and circuit design; it does not
claim a new factorization or loader invention.

## Reproduction

Use the project's build environment with Python, NumPy, SciPy and kickmix.
Set `DATA` to the downloaded competition bundle and `PUBLIC_FACTORS` to the
version-2 `femoco_public_words.npz` named in NOTICE.md. Its hash is enforced.
Use a new or empty output directory for a new build.

```sh
python generator.py --data "$DATA" --public-factors "$PUBLIC_FACTORS"
python verify.py --data "$DATA"
python probe_loader.py
WolframKernel -noinit -noprompt -script run.wls
```

`build.json` pins code, fixed arrays, word table, factorization and each circuit.
`verified.json` records the unmodified official verifier at its default 64
shots per address. All 3976 addresses and 25 circuits are exercised in both
directions; `verified_shots=254464` counts shots, not circuit applications.
`probes.json` checks independently randomized dirty targets, two applications
of each oracle, phase repair, selection preservation and clean workspace.
Its 795200 released and 4704 synthetic circuit applications supplement the
official run. No sampled test proves every possible measurement branch.

`wl-checks.json` records exact symbolic checks and all 298200 released angle
words. [DOCUMENTATION.md](DOCUMENTATION.md) states the mathematics, API
contracts and required regimes. [RESULTS.md](RESULTS.md) records the final
cost, combined local estimate and review verdicts.

Raw factors, transformed words and generated binaries stay in ignored
`local/`. They are not included in the public artifact set. No training or
Kaggle submission is performed.
