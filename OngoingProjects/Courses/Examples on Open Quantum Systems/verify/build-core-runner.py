#!/usr/bin/env python3
"""Extract the catalog's shared WL definitions and append regression assertions."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path


MARKERS = (
    "{id2, X, Y, Z} =",
    "excited =",
    "lower =",
    "densityMatrix[v_]",
    "blochVector[rho_]",
    "purity[rho_]",
    "expectation[op_",
    "commutatorTerm[H_",
    "dissipator[c_",
    "lindbladian[H_",
    "liouvillian[H_",
    "evolve[H_",
    "steadyState[H_",
    "evolveODE[H_",
    "measurementStep[H_",
    "measurementRecord[rho_",
    "trajectory[rho0_",
    "Hx =",
    "oneStepBias[h_",
    "blochState[x_",
    "annihilation[n_]",
    "creation[n_]",
    "coherentState[n_",
    "OmQpc =",
    "regressionSpectrum[big_",
    "qpcSpectrum[\\[Kappa]_, tones_] :=",
    "resolventSpectrum[big_",
    "linearEntropy[rho_]",
    "kPur =",
    "fixedBlur[seed_",
    "crosswise[lean_]",
    "steerStep[rho_",
    "steeredBlur[seed_",
    "fptPur[run_]",
)


def wolfram_blocks(markdown: str) -> list[str]:
    return re.findall(r"^```wl\s*$\n(.*?)^```\s*$", markdown, flags=re.M | re.S)


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: build-core-runner.py CATALOG ASSERTIONS OUTPUT", file=sys.stderr)
        return 2

    catalog, assertions, output = map(Path, sys.argv[1:])
    source = catalog.read_text(encoding="utf-8")
    blocks = wolfram_blocks(source)
    chosen = [block for block in blocks if any(marker in block for marker in MARKERS)]
    missing = [marker for marker in MARKERS if not any(marker in block for block in blocks)]
    if missing:
        print("missing required catalog definitions: " + ", ".join(missing), file=sys.stderr)
        return 2

    digest = hashlib.sha256(source.encode("utf-8")).hexdigest()
    preamble = (
        f'catalogSourceHash = "{digest}";\n'
        f"catalogExtractedDefinitionCells = {len(chosen)};\n"
    )
    runner = preamble + "\n\n".join(chosen) + f'\n\nGet["{assertions}"];\n'
    output.write_text(runner, encoding="utf-8")
    print(f"extracted {len(chosen)} definition cells from {catalog.name}; sha256={digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
