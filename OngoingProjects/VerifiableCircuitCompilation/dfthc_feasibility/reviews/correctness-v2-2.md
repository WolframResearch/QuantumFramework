# VERDICT

PASS for the stated audit/feasibility scope, round 2/3. No unresolved correctness findings. Target files were never edited. All 13 target SHA256s remained unchanged after verification; regenerated inputs.json, comparison.json and released-loader.json are each byte-for-byte identical to the targets. No prior reviews were read.

# REPRODUCED

## Environment and execution

Wolfram kernel: 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. Python 3.13 environment with NumPy 2.5.4 and SciPy 1.18.1; the supplied kickmix environment built the released circuits successfully. There were no dependency, approval or timeout barriers.

The complete target files were read, including all JSON records programmatically. Directly used official cost/evaluator/math implementations, public loader/gauge/ANF dependencies and relevant converted installed documentation were inspected. Applicable ancestor AGENTS.md was empty; no more specific AGENTS.md applied. Verification used a temporary copy, including ignored matrix/spectrum caches and the sibling incumbent summary. Established command forms were:

```sh
"$PYTHON" "$COPY/prepare.py" --data "$OFFICIAL_DATA" --released "$RELEASE_AUDIT"
"$KERNEL" -noinit -noprompt -script "$COPY/run.wls"
"$PYTHON" "$COPY/audit_loader.py" --data "$OFFICIAL_DATA" --source "$PUBLIC_LOADER_SOURCE" --words "$PUBLIC_WORDS"
"$KERNEL" -noinit -noprompt -script "$PROBES/independent.wls"
"$KERNEL" -noinit -noprompt -script "$PROBES/cost-edges.wls"
```

All five commands exited 0. WL probes used Echo, with no Quiet, message suppression or Print. Kernel logs contain no messages or warnings. No speed benchmark is asserted by the artifacts, so repeated timing trials are not applicable.

## Claim inventory and results

1. **Weighted full-tensor Frobenius isometry** (DOCUMENTATION.md mathematical contract; prepare.py packing). The diagonal weight 1 and off-diagonal sqrt(2) form an orthonormal basis for real symmetric matrices. Its tensor-product basis preserves the full N^4 Frobenius norm. Independently checked a general nonseparable symbolic pair coefficient matrix, rather than only a separable tensor. Official N=76 baseline full residual 0.018639570630973906 versus packed residual 0.0186395706309741; pair dimension 2926; target norm 12.184184977999532. Target pair symmetries had exactly zero measured error. Rebuilt binary matrix SHA256 is 4c7e70360e876c9a13a7ffc3482e5cb7979b18c9cef5fa597923465563ae1efd.

2. **PSD structure and rank cap** (RankCost.wl:18; DOCUMENTATION.md). Each factor contributes x x^T; for fixed rho, C vectors lie in a span of at most B projectors. Summing rho gives rank <= Min[d,R Min[B,C]]. The argument is valid for arbitrary real weights. An independent exact constructive example with R=2,B=2,C=4,d=6 achieved rank 4 and matched the cap. All twelve reported caps agree with this formula.

3. **PSD spectral tail** (RankCost.wl:13-17). The documented trace-inequality minimization correctly keeps the k largest positive eigenvalues and retains all negative squares in the residual. Built-in mixed, negative-floor, degenerate, zero-rank, empty, exact radical/constant and limit checks passed. Independent orthogonally rotated exact matrix gave the same rank-one residual as explicit truncation. Additional exact constants {Pi,E,sqrt(2),1/3,-7/5} gave squared rank-two tail 916/225. Rank beyond length retained only the negative floor; all-negative, zero and empty spectra behaved correctly.

4. **Numerical rank screen** (run.wls:9-19,28-35; README table). Fresh WL Eigensystem, independently regenerated SciPy spectrum and the exact matrix hash reproduced: eigen residual 3.624148053680244e-14; orthogonality 2.5268021680905354e-12; spectrum L2 difference 1.3048146537550791e-14. The documented guard is 1e-8 and is explicitly numerical, not an interval proof. Minimum unrestricted PSD rank at 0.020 Ha is 290: rank 289 gives 0.020187308109326476 Ha, rank 290 gives 0.019858425813203292 Ha. Negative floor is 0.0002976497842924418 Ha. Shape (15,57,19) has cap 285 and bound 0.021548084729907786 Ha; subtracting the guard still exceeds 0.020, so REJECT is correct. All remaining eleven shapes and every README table value were reproduced, with NOT_RULED_OUT correctly scoped as necessary-only.

5. **28 built-in symbolic/edge checks** (RankCost.wl:64-104). All 28 passed in a fresh kernel. Nineteen additional independent check groups passed, including genuine 80-digit values and zero-imaginary complex representations, preservation of numerical precision, exact constants, arbitrary nonseparable packing, rotated truncation, invalid/nonfinite/symbolic/nested spectra, negative/noninteger ranks, malformed shapes and dimensions, scratch Failure and its propagation. The cap boundary probe gave angle_qubits=128/within_cap=True, while the next batch gave 133/False. Documented numeric-only signatures were respected.

6. **Official resource model and all 482 rows** (RankCost.wl:21-62; run.wls:20-27; prepare.py). Regeneration called the unmodified official Python heuristic model, including every full stored resource breakdown. The regenerated 482 complete rows, all twelve shapes, five precisions and legal batch widths are byte-identical to inputs.json. WL comparisons of all retained cost fields and conditional costs passed; comparison.json is also byte-identical. An additional 700 (R,B,C,beta) cases, covering R={1,2,3,65,120}, B={1,2,3,16,32,64,65}, C={1,2,3,30,60}, beta={1,2,10,16}, all matched the official model, including power-of-two B and minimal shapes. Lambda is held fixed only for this conditional comparison, as the documents state.

7. **Research baseline and public C30 data audit** (prepare.py; inputs.json; AUDIT.md). All five stored FeMoCo research arrays exactly equal the sample parameters. Public word table SHA256, shape (3976,75), range 0..1023 and weight shape (65,60,30) reproduced. Official evaluator returned error_1=0.01391173198501434, error_2=0.019980679640672218, Lambda=290.7542054993138, E_gap=3.9159356215095613 and lambda_eff=47.55854233515931; official metric gates pass. The stated error/gap margins follow directly from these numbers and official thresholds.

8. **Released signed-gauge/static loader projection** (audit_loader.py:31-64; released-loader.json; AUDIT.md). All three source hashes and the word-table hash passed before execution. Row-sign component drift was 7.771561172376096e-16; official gates remained true after the gauge. All 25 in-memory circuits compiled with 128 qubits. Every profile showed static_toffoli=1129 and table/ANF identities PASS. Angle Toffolis are 56450 for both passes; official heuristic T=8794, total T=65244, total block-encoding Q=661. Conditional cost is 2051023203.3721037 and incumbent/projected ratio is 1.6154555940004633. All released-loader.json bytes were reproduced. This is a static resource projection, not a circuit-execution score.

9. **API contracts** (DOCUMENTATION.md APIs; live kernel). Converted installed Eigenvalues/Eigensystem documentation confirms numeric ordering by absolute magnitude and eigenvectors in rows; the implementation explicitly sorts positive eigenvalues algebraically and uses the documented matrix eigen-equation. Norm documentation confirms matrix default operator norm and explicit Frobenius residuals. BinaryReadList options and little-endian Real64 usage are valid; live options were echoed. PolynomialMod's integer-modulus coefficient reduction supports the four-address identity. No QuantumFramework API is used here.

10. **Provenance/source-attributed claims** (released-sources.json; public-release-provenance.json; AUDIT.md). All nine upstream blob IDs were confirmed by the GitHub connector at the recorded exact commit ref. All nine downloaded starter-file SHA256 values, the public word SHA256 and the public notebook embedded JSON source SHA256 match. Eight research baseline entries are recorded, with one FeMoCo76o113e file. Historical Coulomb, DQI and chemistry resource counts are accurately treated as source-attributed literature statements, not local score verification. Their large experiments were not rerun. Their stated contract differences are material to transfer. The complete external public notebook, other molecules, training claims and original released sampling experiments are outside this audit's stated scope.

11. **Scope and omissions** (README.md, AUDIT.md, PLAN.md, output statuses). Full official circuit verification, dirty-target/measurement-phase/chained-batch coverage, factor training, a new verified score and a submission are explicitly not claimed. No training, submission export, commit or push was performed. The lack of full circuit verification therefore does not contradict the delivered audit/feasibility claim; it remains required before adoption.

# FINDINGS

None. No unresolved or unverifiable in-scope claims. Numerical results carry exactly the documented guard limitation, and released resource numbers retain their static/conditional status.

# RECOMMENDED FIXES

None required for correctness in the stated scope.

# SHA256 binding

- RankCost.wl: `8a37bc53ec545a9009b600616d3510a2e1ec3c6bf68a542c61a715fa33adceed`
- run.wls: `a9e295e743bd84abf44a642376706f2a542d823bceb0e89a195c97810b922757`
- prepare.py: `646b4d94e58fec92bd6bd0173676f4a94dce047625750fd98e2f42960a780d50`
- audit_loader.py: `b7947facc9d2eda9e93038002156fa598eda35ed63b5a467907aaf3bdc843184`
- DOCUMENTATION.md: `58a1a5d9fd581588a9d12b6b89a399ddda4451e5286c43be2c678075a9ead530`
- README.md: `7ceacf02fd8cfab49c5eb4758e6149c83e51abdb93ee8fde7d792422682c5c09`
- AUDIT.md: `77baa2a09dc91e4d3f010e82d2f694a344b049b045876f019f2e413d5a378195`
- PLAN.md: `8a9e8571aa4cb99c0c45efbd6bd5534555dda04e1219f4ead4b091249b001bbf`
- inputs.json: `e6d37f4bbf685562b9e837f9abd48de91e3ca5c2666b57e9140e99aa41e6906f`
- comparison.json: `0aea18740443a27bb4124ec4d23bb86fa0dbae635b3914dc474470c3decbd72d`
- released-loader.json: `a045dcccfecfe6bd58746aebeb2ac0a6cffb0444045d52e2016dcd2c9dbbef0f`
- released-sources.json: `7463fd40fa460055507186c7e989e63ffe9ebc424c9bf5e8a99eda0ef09ccca9`
- public-release-provenance.json: `34eb48a8887bceda16847b21ff4186dce2b3a5cc5ee64306129375748aa6bde2`

Evidence files: target-hashes.json, prepare.log, kernel.log, loader.log, independent.wls, independent.log, independent-checks.json, cost-edge-fixtures.json, cost-edges.wls, cost-edges.log, provenance-checks.json. Commands and results were retained in the temporary verification directory.

OPEN ISSUES: 0
