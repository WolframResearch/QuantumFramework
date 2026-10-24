# Independent WL quality review, round 2/3

Scope: the thirteen supplied DFTHC feasibility/audit artifacts. All target artifacts were read; complete JSON objects were parsed and the regenerated numerical reports were compared in full. No prior reviews were read, and no target was edited. The empty applicable parent AGENTS.md imposed no additional repository rules. This review follows the wl-quality critic mandate and the wl-native-style/wl-style-guide rubrics. The physical object is a real shifted two-electron coefficient tensor, not a many-body state.

Review is bound to the SHA256 manifest at the end. Executions took place in a temporary copy, including ignored local matrix/spectrum caches and the sibling incumbent summary. The established Python preparation, Python released-loader audit, and Wolfram run.wls commands all exited 0. No training, submission export, full circuit verifier, commit, or push was run.

## A — Documented symbols: PASS

- RankCost.wl:13–16 explicitly extracts real parts, selects negative and positive modes separately, and algebraically sorts the positive modes. It does not misuse the absolute-value ordering of numeric Eigenvalues/Eigensystem.
- run.wls:9–15 uses BinaryReadList[file,"Real64",ByteOrdering->-1], row-oriented Eigensystem eigenvectors, and explicit matrix Frobenius norms. Installed Wolfram 15.0 documentation confirms the BinaryReadList file/type/options form and file closure; Eigensystem Details confirms {values,vectors}, normalized numeric vectors, correspondence to rows, and the exact matrix eigen-equation used at line 12; Norm Details explicitly distinguishes spectral and Frobenius matrix norms. PolynomialMod Details confirms coefficient reduction modulo 2, as used in RankCost.wl:71.
- Live Options[BinaryReadList] returns AllowIncomplete, ByteOrdering, and Path; Options[PolynomialMod] returns CoefficientDomain, Modulus, and Trig. All used WL functions are core; no QF or paclet API is inferred from memory.
- The official costs_prepare.py implementation was read, including its power-of-two B handling and separate qubit scratch exponent (lines 84–99, 115–147, 169–190). verify.py:389–516 confirms sharing temporary workspace and charging both angle/Givens passes. Every one of the 482 generated rows agrees with the WL implementation.
- All three released module hashes are enforced by audit_loader.py:19–33 and reproduced successfully. The starter-kit inspected-file hashes and public NPZ hash match the provenance manifest. The original nine Git blob IDs match the fetched source content after removing the audit cache's one extra trailing newline; the notebook source hash matches the fetched JSON response's source/sourceNullable field, rather than the enclosing API response.

No unsupported signature, shape mismatch, or memory-only API claim found. Smallest fix: none.

## B — Idiomatic functional WL: PASS

- The search for Do/For/While/Table/AppendTo finds none in either WL artifact. PSDTail uses Select/Sort/Drop/Total; qrom uses MinimalBy; the scratch fit uses SelectFirst; four-case interpolation uses MapThread; run.wls uses AllTrue, Map, and MinimalBy.
- RankCost.wl:26–52 uses a Module for named, dependent resource intermediates; RankCost.wl:64–104 additionally needs local symbolic variables. These are not trivial constant-only Modules. The two If expressions are failure propagation, not finite-type dispatch cascades.
- run.wls:33 uses Which for continuous threshold/guard classification. Pattern dispatch would obscure this scalar screen. The dense official 2926 matrix and complete eigenvectors are appropriate for the full dense spectrum; sparse conversions would not eliminate the dense eigensystem or its residual matrices.

No meaningful idiom defect found. Smallest fix: none.

## C — Terse physics argument: PASS

- DOCUMENTATION.md:7–27 gives the needed chain: orthonormal symmetric-pair coordinates -> Frobenius isometry -> sum of outer products -> PSD -> shared-projector span -> R min(B,C) cap -> spectral PSD-tail minimum. The trace-inequality argument supplies arbitrary-dimensional reasoning, not a numerical slogan.
- The code names physical intermediates: pairMatrix, spectrum, eigenResidual, orthogonality, rank, bound, persistent/temporary qubits. The state-preparation decomposition maps directly to official model terms.
- DOCUMENTATION.md:29–46, README.md:3–9 and 26, and AUDIT.md:80–91 keep numerical screening, conditional resource estimates, metric acceptance, and unperformed full circuit verification distinct. The shared-projector constraint can only make fitting harder, as stated.

No dead or redundant physics scaffolding needing a substantive fix. Smallest fix: none.

## D — Applicable physics invariants: PASS

- RankCost.wl:74–83 checks pair inner products, the tensor norm, projector Gram identity, row-sign invariance, all four Boolean address cases, all 512 signed words, the exact mixed-sign PSD residual, and the decoupling limit. All are True in a fresh Wolfram 15.0 kernel.
- prepare.py:68–83 rejects missing within-pair/exchange symmetry and compares the official N^4 reconstruction residual with its packed residual. Reproduction returned both symmetry errors 0; full residual 0.018639570630973906 Ha; packed residual 0.0186395706309741 Ha.
- run.wls:7–18 checks the binary SHA256, real symmetric matrix, eigen-equation and eigenvector orthogonality, and an independent SciPy eigenvalue reference. Reproduction gave dimension 2926, eigen residual 3.624148053680244e-14, orthogonality 2.5268021680905354e-12, and SciPy spectrum l2 difference 1.3048146537550791e-14, all below the explicit 1e-10 acceptance tolerance.
- audit_loader.py reran the official metric gates after signed gauge canonicalization: maximum component drift up to row sign 7.771561172376096e-16; both residual gates and strict normalization/gap gates pass. All 25 static circuits have Q128; table and ANF-table identities are PASS. Both-pass angle T56450 plus heuristic T8794 gives T65244; shared scratch gives Q661; conditional cost 2051023203.3721037 and gain 1.6154555940004633.
- State normalization, density-matrix trace, unitarity, and a many-body eigenproblem are not claims of these artifacts. Full emitted-circuit measurement/phase behavior remains explicitly outside this audit, so it is not silently graded as verified.

No failing coefficient-tensor invariant found. Smallest fix: none. The missing resource-cap rejection *coverage* is counted once under G below.

## E — Nontrivial examples: PASS

The trivial case would be a scalar or an already PSD diagonal with no discarded mode. The actual artifacts include mixed-sign exact spectra, degeneracy and zero-rank cases, symbolic off-diagonal packing and projector coordinates, exhaustive Boolean identities, the complete nontrivial 2926 target and twelve DFTHC shapes. The fixed-shift paper shape (15,57,19) has cap 285 and error lower bound 0.021548084729907786 Ha, which exceeds 0.020 after subtracting the 1e-8 guard. Thus the work exercises a genuine physically relevant rejection, not just a feasible toy case.

Smallest fix: none.

## F — Symbolic-first derivation: PASS

- The unrestricted PSD-tail formula and rank-cap argument retain arbitrary dimension and parameters in DOCUMENTATION.md:13–27.
- RankCost.wl:66–83 retains symbolic matrix entries, projectors, positive/negative magnitudes, and the decoupling parameter. Exact tests include Sqrt[2], Pi, and rational-vs-radical ordering (lines 89–91); the implementation does not force N or machine precision.
- An independent probe PSDTail[{N[Sqrt[2],60],N[3/2,60]},1] retained precision 60. The 2926 numerical calculation is appropriate because the official input consists of fixed floating-point integral data and is independently crosschecked. It instantiates the spectral result rather than replacing its derivation.

Smallest fix: none.

## G — Required regimes: FINDING (one issue)

All rows except resource-cap failure coverage have executable evidence:

| Regime | Evidence |
|---|---|
| General symbolic | DOC:7–27; symbolic pair/Gram/Boolean checks, all True |
| Exact special | DOC:71; RankCost.wl:80–91, mixed signs, exact constants, degeneracy |
| Limit | RankCost.wl:83–86; discarded positive mode decouples, negative floor remains at unrestricted rank |
| Independent numeric | Official full vs packed norm; full WL vs SciPy spectrum; 482 cost rows; released module audit |
| Failure/edge | Empty, malformed, negative-rank, invalid-shape, exact/high-precision real-valued complex inputs, impossible scratch budget and propagated failure; all 28 published checks True |

**G1 — The 128-qubit cap's rejection branch and WL's inclusive boundary are absent from the published executable checks.**

Evidence:

- prepare.py:109–111 discards every row with angle_qubits >128 before writing inputs.json.
- Complete parsing of all 482 inputs.json rows confirms `max(angle_qubits)=127`, with zero rows at 128 and zero above 128.
- run.wls:21–26 compares toffoli, qubits, angle_qubits, and conditional_cost, but never compares LookupCost's `within_cap` flag.
- RankCost.wl:74–103 contains 28 checks, including scratch failures, but no cap assertion; the flag itself is defined at line 60.
- DOCUMENTATION.md:74 cites "the cap boundary in the sweep", but that sweep does not reach equality and cannot exercise failure. The released loader's Q128 is generated and checked by a separate Python module; it does not exercise LookupCost's flag.
- Independent fresh-kernel calls show the correct expected boundary/failure behavior:
  `LookupCost[{65,60,30},7,15,1]` gives angle_qubits 128, within_cap True.
  `LookupCost[{65,60,30},10,11,1]` gives angle_qubits 133, within_cap False.
  All 28 existing checks and all 482 comparisons could still pass if line 60 always returned True, or if its <= became <.

Smallest meaningful fix: add one paired built-in assertion (or two named checks) covering exactly 128 and >128, for example compare the {angle_qubits,within_cap} pairs above with {{128,True},{133,False}}. Also assert model["within_cap"] in the existing legal-row crosscheck, update the published check count, regenerate comparison.json, and point the regimes table to these actual checks. No loader simulation or factor optimization is needed.

The implementation currently returns the correct flags; the finding is missing required failure/boundary evidence, not a claim that the cost formula is wrong. This is one issue, counted once.

## Source binding

```json
{
  "RankCost.wl": "8a37bc53ec545a9009b600616d3510a2e1ec3c6bf68a542c61a715fa33adceed",
  "run.wls": "a9e295e743bd84abf44a642376706f2a542d823bceb0e89a195c97810b922757",
  "prepare.py": "646b4d94e58fec92bd6bd0173676f4a94dce047625750fd98e2f42960a780d50",
  "audit_loader.py": "b7947facc9d2eda9e93038002156fa598eda35ed63b5a467907aaf3bdc843184",
  "DOCUMENTATION.md": "58a1a5d9fd581588a9d12b6b89a399ddda4451e5286c43be2c678075a9ead530",
  "README.md": "7ceacf02fd8cfab49c5eb4758e6149c83e51abdb93ee8fde7d792422682c5c09",
  "AUDIT.md": "77baa2a09dc91e4d3f010e82d2f694a344b049b045876f019f2e413d5a378195",
  "PLAN.md": "8a9e8571aa4cb99c0c45efbd6bd5534555dda04e1219f4ead4b091249b001bbf",
  "inputs.json": "e6d37f4bbf685562b9e837f9abd48de91e3ca5c2666b57e9140e99aa41e6906f",
  "comparison.json": "0aea18740443a27bb4124ec4d23bb86fa0dbae635b3914dc474470c3decbd72d",
  "released-loader.json": "a045dcccfecfe6bd58746aebeb2ac0a6cffb0444045d52e2016dcd2c9dbbef0f",
  "released-sources.json": "7463fd40fa460055507186c7e989e63ffe9ebc424c9bf5e8a99eda0ef09ccca9",
  "public-release-provenance.json": "34eb48a8887bceda16847b21ff4186dce2b3a5cc5ee64306129375748aa6bde2"
}
```

The reproduced inputs.json, comparison.json, and released-loader.json are byte-for-byte equal to the supplied target reports. Independent probes and the snapshot manifest remain beside this review in temporary storage.

OPEN ISSUES: 1
