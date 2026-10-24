ROUND: 1/3

Independent fresh-context wl-verify review. Personal workspace prefixes in file
links were removed for public storage; the findings are preserved.

VERDICT

Core numerical and audit claims reproduce. One API contract issue remains.

REPRODUCED

Ran prepare.py, run.wls and audit_loader.py in a temporary copy with ignored
caches. All exited 0 without Wolfram messages. Wolfram 15.0.0 ARM, NumPy 2.5.4,
SciPy 1.18.1. All three result JSONs regenerated identically.

- Full/packed residuals 0.018639570630973906 and 0.0186395706309741; matrix hash matches.
- Rank cap Min[2926,R Min[B,C]] follows from the official shared projector span and PSD reconstruction.
- Minimum rank290, rank289 error0.020187308109326476, rank290 error0.019858425813203292 Ha; negative floor0.0002976497842924418 Ha.
- Shape(15,57,19) cap285 has error0.021548085 Ha and fails the guarded screen.
- Eigen residual3.624148053680244e-14, orthogonality2.5268021680905354e-12, SciPy difference1.3048146537550791e-14.
- All18 original symbolic/edge checks and all482 cost rows pass. Additional120 comparisons cover B=1,2,16,32,64 and beta=1,10,16,63.
- Independent mixed-sign/random PSD probes match direct spectral reconstruction. Tested invalid, empty, degenerate, negative-rank, malformed-shape and cap-exceeding cases.
- Public C30/10bit metric gates reproduce.
- Released loader25 circuits eachQ128/static1129 Toffolis; angle56450, total65244, Q661, conditional cost2051023203.3721037 and gain1.6154555940004633.
- Literature resource claims distinguished from local verification; applicability contracts supported. Git blobs match after removing the retrieval-added trailing newline.
- Installed documentation/live checks support eigenvalue order, eigenvector orientation, Frobenius norms, binary options and polynomial reduction.

Full circuit correctness/projected-score validation remain outside scope. No
runtime speed benchmark is asserted.

FINDINGS

[P2] RankCost.wl:2 promises real numeric scalar spectra, but line11 uses
NumberQ, excluding exact radicals/constants. Fresh-kernel PSDTail[{Sqrt[2]},0]
and PSDTail[{Pi},0] both return InvalidSpectrum, rather than Sqrt[2] and Pi.
This does not change the reported machine-spectrum results.

RECOMMENDED FIXES

Use a suitable finite-real NumericQ guard and add exact-radical/constant tests.

OPEN ISSUES: 1
