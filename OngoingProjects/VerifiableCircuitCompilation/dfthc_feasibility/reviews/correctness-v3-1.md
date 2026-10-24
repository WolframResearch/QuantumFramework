# Independent WL correctness report — correctness visit 3, round 1

## VERDICT

PASS within the declared audit and feasibility scope. No unresolved correctness findings. Targets were copied with ignored local caches and the sibling incumbent summary into an isolated temporary directory. Target files were not edited, and all thirteen target SHA256 values were rechecked unchanged at completion. No prior reviews were read. No training, circuit submission export, commit, or push was performed.

## REPRODUCED

### Commands and versions

External directories are represented by variable names to avoid publishing personal paths. Commands were executed against the temporary copy.

- `python prepare.py --data DATA --released RELEASED`: exit 0.
- `WolframKernel -noinit -noprompt -script run.wls`: exit 0.
- `python audit_loader.py --data DATA --source SOURCE --words PUBLIC_WORDS`: exit 0.
- `WolframKernel -noinit -noprompt -script probe.wls`: exit 0.
- Deep JSON equality of original versus regenerated `inputs.json`, `comparison.json`, and `released-loader.json`: all fields equal, including every one of the 482 cost rows, all twelve shape comparisons, and all 25 circuit profiles.

Wolfram kernel: `15.0.0 for Mac OS X ARM (64-bit) (February 17, 2026)`.
Python: 3.13.1; NumPy: 2.5.4; SciPy: 1.18.1; kickmix: 0.1.dev0.
No WL messages occurred in the established run or independent probes. No Quiet, message suppression, or Print was used in the independent WL probes.

### Claim inventory and evidence

1. **Packing and PSD/span theorem** (`DOCUMENTATION.md:7–26`, `RankCost.wl:13–19,74–83`). Independently checked the general pair-tensor Frobenius identity on a symbolic 3^4 tensor with 36 independent packed-pair entries, rather than relying only on the built-in decomposable X tensor X example. Its polynomial difference was exactly zero. An exact random R=2, B=2, C=3 projector construction had PSD pair Gram matrix rank 4, matching R Min(B,C). The general argument is sound: each rho contributes at most Min(B,C) dimensions, and the ambient cap is 2926. Spectral truncation of the positive eigenvalues, with the full negative-spectrum floor retained, is the unrestricted PSD rank bound, not a DFTHC existence proof.

2. **Full official N^4 residual and reference provenance** (`prepare.py:62–84`, `inputs.json`). Reconstructed the official shifted target from the hash-checked float32-rounded FCIDUMP and fixed shifts. Both exact target pair-symmetry errors were 0. Full baseline Frobenius residual: 0.018639570630973906 Ha; packed residual: 0.0186395706309741 Ha. Matrix dimension 2926 and SHA256 `4c7e70360e876c9a13a7ffc3482e5cb7979b18c9cef5fa597923465563ae1efd` reproduced. The FeMoCo research baseline's five parameter arrays are each exactly equal to the competition sample.

3. **Spectrum, minimum rank, guards, and twelve shape screens** (`run.wls:9–19,28–35`; `README.md:5,13–24`). Fresh full WL eigensystem versus separately recomputed SciPy eigvalsh:
   - eigenpair Frobenius residual: 3.624148053680244e-14;
   - orthogonality Frobenius residual: 2.5268021680905354e-12;
   - SciPy spectrum L2 difference: 1.3048146537550791e-14;
   - negative-spectrum floor: 0.0002976497842924418 Ha;
   - minimum unconstrained PSD rank at 0.020 Ha: 290;
   - rank-289 bound: 0.020187308109326476 Ha;
   - rank-290 bound: 0.019858425813203292 Ha.
   `(15,57,19)` has cap 285 and bound 0.021548084729907786 Ha, so its numerical rejection exceeds the stated 1e-8 guard comfortably. All twelve saved comparisons and rounded README table entries agree with regenerated outputs. The guard is correctly described as numerical evidence rather than an interval certificate.

4. **All 32 built-in checks** (`RankCost.wl:64–109`). Every returned check was True. Coverage includes pair/projector identities, four-address GF(2) interpolation, all 512 signed-word redundant-bit cases, exact mixed-sign and degenerate spectra, rank zero, empty input, negative floor, exact radicals/Pi/mixed ordering, real-valued complex machine and high-precision representations, Q128 acceptance, Q133 rejection, negative normalization rejection, zero normalization limit, ScratchBudget Failure and propagation, span cap and power-of-two B.

5. **Additional domain and signature probes** (`RankCost.wl:9–62`). Live DownValues and exported names were inspected. Independent checks all passed: negative-only spectrum, excessive allowed rank, exact algebraic ordering, 80-digit spectrum accuracy, a positive exact mode of size 1e-100, symbolic/nested/complex/nonfinite spectra, invalid ranks/shapes/precision, symbolic/negative/complex/nonfinite normalization, and exact Pi normalization. Valid cost outputs are Associations with integer resource fields and Boolean within_cap; invalid cases return Failure. Independently constructed 192 additional official resource comparisons cover B=1,2,3,5,16,32,64; beta=1,2,7,10,16,40; and J=1,3,10,15. Every emitted WL heuristic key and every total/angle/cap field matched the official implementation.

6. **All 482 published cost rows** (`prepare.py:104–107`, `run.wls:20–27`). Regenerated all 482 complete rows from the unmodified official Python model. Every JSON field matched, including the thirteen official heuristic subfields. The fresh WL run matched its emitted heuristic fields, total Toffolis/qubits/angle width, fixed-normalization conditional cost, and within_cap=True for every retained row. Twelve shapes and five precisions were reproduced; power-of-two B handling and distinct inner scratch exponents agree with official source. All twelve displayed best loaders use beta=10, J=10, and gains use incumbent 3313336907.312215 and reference lambda_eff=48.68043460482008.

7. **Public C30/10-bit metric gates** (`prepare.py:93–103`, `AUDIT.md`). Hash-pinned public weights and words reproduced error_1=0.01391173198501434, error_2=0.019980679640672218, Lambda=290.7542054993138, E_gap=3.9159356215095613, lambda_eff=47.55854233515931, and acceptance through the unmodified official evaluator. Both stated small acceptance margins agree arithmetically.

8. **Released loader static projection** (`audit_loader.py:31–65`, `released-loader.json`, `AUDIT.md`). The three external module hashes and public word hash passed. Gauge drift 7.771561172376096e-16, official metric gates, all 25 circuits at Q128, 1129 static Toffoli-like gates per circuit, angle total 56450, heuristic total 8794, full static projection T=65244 and Q=661 reproduced exactly. Conditional cost 2051023203.3721037 and improvement 1.6154555940004633 reproduced exactly. All stored circuit-profile fields matched.

9. **API contracts** (`DOCUMENTATION.md:48–62`). Read provided converted installed documentation and live definitions/source. Confirmed numerical eigensystem eigenvalue ordering by absolute magnitude and eigenvectors as rows; Frobenius versus default matrix norm; BinaryReadList little-endian Real64 shape and automatic closing when opened by the call; PolynomialMod integer coefficient reduction. Live options are BinaryReadList `{AllowIncomplete->False, ByteOrdering:>$ByteOrdering, Path:>$Path}` and PolynomialMod `{CoefficientDomain->Rationals, Modulus->0, Trig->False}`. Direct official source confirmed fixed chi=10, N=76, strict metric floors, shared scratch, ANGLE cap, and resource score formula.

10. **Source-attributed historical statements and hashes** (`AUDIT.md`, `released-sources.json`, `public-release-provenance.json`). Historical Coulomb/DQI/Fe2S2 counts and contracts, including DQI 42-product multiplication and 66-product inversion, are supported by the pinned released sources and remain literature/source statements. They are not represented as local reproduced competition scores. All nine original-release Git blob IDs match the supplied cached text after excluding its one additional trailing newline. All nine starter-kit inspected-file SHA256 values match. Public notebook metadata identifies version 8 / kernel ID 140416287, and its `blob.source` SHA256 matches the recorded JSON-source hash. Public word and research baseline hashes match.

### Explicit exclusions, not failures of the declared scope

Full circuit measurement-phase correctness, arbitrary dirty targets, chained batches, official forward/reverse shot coverage, a smaller-rank accurate fit, new-shape normalization and error acceptance, and a new official verified score are explicitly disclaimed and were not tested. No speed benchmark is asserted. Historical large sampling experiments and the entire released notebook were not rerun. These exclusions do not invalidate the reported necessary conditions, metric audits, or static resource projections.

## FINDINGS

None. No missing dependencies, kernel limitation, timeout, or approval barrier reduced the requested coverage.

## RECOMMENDED FIXES

None required for the declared audit and feasibility scope. Retain the existing explicit qualification that released-loader resources are a static projection until full circuit verification is performed.

## SHA256 binding

- `RankCost.wl`: `016dbf7e9884d762d08bede6b2a4b45fd1aa6853b988dc8b1b031a226f5794e0`
- `run.wls`: `454bd174ca688f604f00e923e3c840c937c7f6ef9fd53423395377ae94affec1`
- `prepare.py`: `646b4d94e58fec92bd6bd0173676f4a94dce047625750fd98e2f42960a780d50`
- `audit_loader.py`: `b7947facc9d2eda9e93038002156fa598eda35ed63b5a467907aaf3bdc843184`
- `DOCUMENTATION.md`: `e35682d7d6c20b6d85d3729c6ed4a9a7d584dad605d8905d9d9b15626e7e5514`
- `README.md`: `f09b6826ee3bdcc27c0620797508fa9f2538eff11e632657aace4fcf2e7fd5cf`
- `AUDIT.md`: `77baa2a09dc91e4d3f010e82d2f694a344b049b045876f019f2e413d5a378195`
- `PLAN.md`: `8a9e8571aa4cb99c0c45efbd6bd5534555dda04e1219f4ead4b091249b001bbf`
- `inputs.json`: `e6d37f4bbf685562b9e837f9abd48de91e3ca5c2666b57e9140e99aa41e6906f`
- `comparison.json`: `47a7118f2c635fe82113e867d6a79f3cd4f3456e0505f8306bd8dd6e419553ac`
- `released-loader.json`: `a045dcccfecfe6bd58746aebeb2ac0a6cffb0444045d52e2016dcd2c9dbbef0f`
- `released-sources.json`: `7463fd40fa460055507186c7e989e63ffe9ebc424c9bf5e8a99eda0ef09ccca9`
- `public-release-provenance.json`: `34eb48a8887bceda16847b21ff4186dce2b3a5cc5ee64306129375748aa6bde2`

OPEN ISSUES: 0
