# ROUND: 3/3 — independent adversarial correctness verification

## VERDICT

PASS for the current audit/feasibility scope. No target artifacts were edited. All three established scripts were executed from a temporary copy, including the ignored matrix/spectrum cache and the necessary incumbent-summary dependency. The updated RankCost.wl and current generated report were refreshed and rerun during this round. No prior reviews were read.

## REPRODUCED

### Claims and methods

1. Pair packing is an orthonormal coordinate map for real symmetric matrices. The diagonal weight is 1 and off-diagonal weight is sqrt(2); applying it on each symmetric tensor pair preserves the full N^4 Frobenius norm. The independent symbolic identities and the official full-tensor baseline residual check both pass.
2. The official reconstruction is a sum of outer products of packed symmetric X factors, hence PSD. At fixed rho, C such vectors lie in the span of B packed rank-one projectors. Its rank is therefore at most Min[2926,R Min[B,C]]. The unrestricted PSD rank-k minimizer retains the k largest positive eigenvalues and discards all negative eigenvalues. This establishes a necessary bound, not representability.
3. All 23 built-in symbolic/edge checks passed, including exact radicals, exact constants, mixed exact ordering, the impossible scratch budget, and propagation of that precise Failure value through LookupCost.
4. Independently reran prepare.py against official fixed shifts, reference integrals, sample parameters, and public factor arrays. Its entire freshly generated inputs.json is identical to the target JSON, including all 482 resource rows, source hashes, public metric gates, and five-array research-baseline identity.
5. Independently reran run.wls using the freshly prepared matrix. Its entire generated comparison.json is identical to the current target JSON. All 12 shape comparisons, both selected reference-loader variants per shape, all hashes, and every aggregate numeric value reproduce.
6. Independently reran audit_loader.py with the three pinned released modules and public word-table hash. Its entire released-loader.json is identical to the target JSON, including every circuit profile. All 25 in-memory builds have 128 qubits and 1129 static Toffolis each.
7. Extra kernel tests independently checked 108 official state-preparation cost rows across R={1,3,65}, B={1,2,4,16,32,64}, C={1,4}, beta={1,10,16}. All 108 match every returned WL cost field. Additional checks cover exact ordering of Pi/E/radicals/rationals, negative-only and zero spectra, oversized rank, malformed/symbolic/nonfinite/nested spectra, invalid rank/shape/precision/batch/lambda, exact ANGLE Q128 equality, Q143 overflow rejection through within_cap, and precise scratch Failure propagation. All 23 additional check flags passed.
8. Installed Wolfram 15.0 documentation confirms numerical eigenvalues sort by magnitude, eigenvectors are rows, approximate eigenvectors are normalized, matrix default Norm is spectral, explicit Frobenius is the flattened-vector 2-norm, BinaryReadList supports ByteOrdering and closes a file it opens, and PolynomialMod reduces integer coefficients modulo 2. The numerical and symbolic code uses these documented signatures correctly.
9. Source statements were checked as literature/provenance statements: original source READMEs contain the reported 1140/Q142, 442397/Q1837, and 14978/Q126 counts, modulus 4331, sampling contract, and 42/66 multiplication/inversion counts. All nine pinned Git blob IDs match the retrieved texts after removing the single newline appended by their retrieval. The nine starter-kit file hashes match, and the public notebook JSON source payload hash matches the claimed hash. The research paper explicitly reports FeMoCo (15,57,19), with a different reported Frobenius error/selection contract. These historical counts were not locally recertified.
10. Public target artifacts contain no personal/private filesystem paths or credentials under the checked strings. No performance benchmark was asserted, so no performance conclusion is drawn.

### Reproduction commands

Commands below use externally supplied DATA, AUDIT_DIR, SOURCE_DIR and WORDS_NPZ; execution used a temporary project copy.

```sh
python prepare.py --data DATA --released AUDIT_DIR
WolframKernel -noinit -noprompt -script run.wls
python audit_loader.py --data DATA --source SOURCE_DIR --words WORDS_NPZ
WolframKernel -noinit -noprompt -script extra.wls
```

All four commands exited 0. Kernel output used Echo; no Quiet or message suppression was used. There were no Wolfram diagnostic messages.

Environment: Python 3.13.1; NumPy 2.5.4; SciPy 1.18.1; kickmix 0.1.dev0; Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026.

### Numeric results

| Quantity | Reproduced result |
|---|---:|
| Pair dimension | 2926 |
| Full official baseline residual / Ha | 0.018639570630973906 |
| Packed baseline residual / Ha | 0.0186395706309741 |
| Eigenpair Frobenius residual | 3.6241480536802438e-14 |
| Orthogonality Frobenius residual | 2.5268021680905354e-12 |
| WL versus SciPy spectrum L2 difference | 1.3048146537550791e-14 |
| Minimum unrestricted PSD rank at 20 mHa | 290 |
| Tail at rank 289 / Ha | 0.020187308109326476 |
| Tail at rank 290 / Ha | 0.019858425813203292 |
| Negative spectrum floor / Ha | 0.00029764978429244179 |
| Released angle Toffolis, both passes | 56450 |
| Released total Toffolis / total qubits | 65244 / 661 |
| Released conditional cost | 2051023203.3721037 |
| Released gain over incumbent | 1.6154555940004633 |

The (15,57,19) cap is 285; its tail is 0.021548084729907786 Ha. It is numerically rejected even after subtracting the explicit 1e-8 guard. The released public C30/10-bit metric gates pass before and after the signed gauge. The strict gap and two-body error margins reproduce.

### Binding SHA256 hashes

| Artifact | SHA256 |
|---|---|
| RankCost.wl | 46f9eef74d6ff8a6066a9aa2fad77cd0e6c19d73d376eb22ec4ab1849ecd052a |
| run.wls | a9e295e743bd84abf44a642376706f2a542d823bceb0e89a195c97810b922757 |
| prepare.py | 646b4d94e58fec92bd6bd0173676f4a94dce047625750fd98e2f42960a780d50 |
| audit_loader.py | b7947facc9d2eda9e93038002156fa598eda35ed63b5a467907aaf3bdc843184 |
| DOCUMENTATION.md | 58a1a5d9fd581588a9d12b6b89a399ddda4451e5286c43be2c678075a9ead530 |
| README.md | 7be43f16dc13c409ee2abda9351fa83cfd765e36f7c8c80c44a39af73aacc21b |
| AUDIT.md | 77baa2a09dc91e4d3f010e82d2f694a344b049b045876f019f2e413d5a378195 |
| PLAN.md | 8a9e8571aa4cb99c0c45efbd6bd5534555dda04e1219f4ead4b091249b001bbf |
| inputs.json | e6d37f4bbf685562b9e837f9abd48de91e3ca5c2666b57e9140e99aa41e6906f |
| comparison.json | 1813ba6a6aa84c821cce8e97f7409fcf96ad6d5caa35ee9c5a15bd4f350b332d |
| released-loader.json | a045dcccfecfe6bd58746aebeb2ac0a6cffb0444045d52e2016dcd2c9dbbef0f |
| released-sources.json | 7463fd40fa460055507186c7e989e63ffe9ebc424c9bf5e8a99eda0ef09ccca9 |
| public-release-provenance.json | 34eb48a8887bceda16847b21ff4186dce2b3a5cc5ee64306129375748aa6bde2 |

## FINDINGS

None within the specified audit/feasibility scope. The machine-eigenspectrum bound remains a numerically crosschecked result with an explicit guard, not an interval certificate. The released-loader profile remains a static projection, not a verified circuit score. Both limits are accurately stated in the target documentation.

The full released circuit simulation, measurement-phase repair checks, dirty output targets, chained batches, official forward/reverse shot coverage, training, new factor fitting, submission export, and original published large-shot experiments were not run. They are outside this task and are not claimed as verified by this report.

## RECOMMENDED FIXES

None. Preserve the current scope labels and hash binding when publishing these results.

Reference checked: [Low et al., arXiv:2502.15882v1](https://arxiv.org/html/2502.15882v1), Tables 4–5 and optimization/selection discussion.

OPEN ISSUES: 0
