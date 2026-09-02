#!/usr/bin/env python3
"""Compare the exported circuits against Stim.

For every case emitted by emit.wls this reads the .stim source and the .json of
exact rates, samples the circuit with Stim, and checks that each detector's
firing rate and the observable's flip rate agree.

The comparison is exact-versus-sampled: the rates in the JSON are closed-form
values computed by QECDetectorModel, not estimates, so any real disagreement
shows up as a large sigma rather than as noise.

Also runs PyMatching on each exported circuit, which is a second check: a
malformed circuit will not produce a decomposable detector error model.

    pip install stim pymatching
    wolframscript -file emit.wls
    python3 check.py [outputDirectory] [shots]
"""
import glob
import json
import os
import sys

import numpy as np
import stim
import pymatching

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "out")
SHOTS = int(sys.argv[2]) if len(sys.argv) > 2 else 2_000_000
TOLERANCE_SIGMA = 5.0

failures = 0
print(f"{'case':14s} {'dets':>5s} {'max sigma':>10s} {'obs ours':>10s} {'obs stim':>10s} "
      f"{'obs sigma':>10s} {'MWPM':>10s}")
print("-" * 78)

for path in sorted(glob.glob(os.path.join(OUT, "*.json"))):
    meta = json.load(open(path))
    name = meta["name"]
    circuit = stim.Circuit(open(os.path.join(OUT, name + ".stim")).read())

    assert circuit.num_detectors == meta["detectors"], (
        name, circuit.num_detectors, meta["detectors"])

    dets, obs = circuit.compile_detector_sampler().sample(
        SHOTS, separate_observables=True)

    ours = np.asarray(meta["detectorRates"], dtype=float)
    theirs = dets.mean(axis=0)
    se = np.sqrt(np.maximum(theirs * (1 - theirs), 1e-12) / SHOTS)
    sigma = np.abs(ours - theirs) / se

    # The Stim circuit is a single-basis memory experiment, so its one observable
    # corresponds to the second half of our observable columns: the symplectic
    # products with the logical Z operators.
    our_obs = np.asarray(meta["observableRates"], dtype=float)
    k = len(our_obs) // 2
    our_z = our_obs[k]
    their_obs = obs.mean(axis=0)[0]
    obs_se = np.sqrt(max(their_obs * (1 - their_obs), 1e-12) / SHOTS)
    obs_sigma = abs(our_z - their_obs) / obs_se

    try:
        matcher = pymatching.Matching.from_detector_error_model(
            circuit.detector_error_model(decompose_errors=True))
        mwpm = (matcher.decode_batch(dets)[:, 0] != obs[:, 0]).mean()
        mwpm_text = f"{mwpm:10.5f}"
    except Exception as exc:                       # noqa: BLE001
        mwpm_text = f"{type(exc).__name__[:10]:>10s}"

    print(f"{name:14s} {circuit.num_detectors:5d} {sigma.max():10.2f} {our_z:10.5f} "
          f"{their_obs:10.5f} {obs_sigma:10.2f} {mwpm_text}")

    if sigma.max() > TOLERANCE_SIGMA or obs_sigma > TOLERANCE_SIGMA:
        failures += 1
        for i, (a, b, s) in enumerate(zip(ours, theirs, sigma)):
            if s > TOLERANCE_SIGMA:
                print(f"    detector {i:3d}: ours {a:.6f}  stim {b:.6f}  sigma {s:.2f}")

print()
print("STATUS", "RED" if failures else "GREEN")
sys.exit(1 if failures else 0)
