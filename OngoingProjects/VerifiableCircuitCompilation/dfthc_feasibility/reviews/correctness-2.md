ROUND: 2/3

VERDICT
The feasibility and static-resource claims reproduce exactly. One low-priority API boundary finding remains; none of the twelve reported shapes is affected.

REPRODUCED
- Ran prepare.py, WolframKernel -noinit -noprompt -script run.wls, and audit_loader.py from a temporary copy containing the ignored matrix/spectrum cache. All exited 0. Regenerated inputs.json, comparison.json, and released-loader.json are each exactly equal as complete parsed JSON objects to their target counterparts.
- Kernel: Wolfram 15.0.0 ARM64 (February 17, 2026). NumPy 2.5.4; SciPy 1.18.1. No Wolfram messages were emitted; none were suppressed.
- Weighted pair packing: dimension 2926; exact target pair symmetries; full N^4 residual 0.018639570630973906 versus packed residual 0.0186395706309741. Independently checked the isometry on an arbitrary exact symmetric-pair tensor rather than only a rank-one tensor product. The span argument gives Min[2926,R Min[B,C]].
- All 20 published symbolic/edge checks pass, including exact Sqrt[2] and Pi inputs. Further probes covered symbolic/nested/nonreal/nonfinite spectra, invalid rank and shape, exact mixed constants, 60-digit values, rank exceeding dimension, and an exact orthogonally rotated mixed-sign matrix. Invalid inputs returned Failure; spectral truncations match direct Frobenius residuals.
- Minimum unrestricted PSD rank is 290: rank 289 bound 0.020187308109326476 Ha, rank 290 bound 0.019858425813203292 Ha; negative floor 0.0002976497842924418 Ha. Shape (15,57,19) has rank cap 285 and bound 0.0215480846 Ha, hence REJECT under the explicit guard. Eigensystem residual 3.62415e-14, orthogonality 2.52680e-12, SciPy difference 1.30481e-14.
- All 482 complete generated cost rows and twelve shape summaries reproduce. The legal J=10 / illegal J=11 cap boundary for shape (65,60,30), beta=10 also agrees. Official power-of-two B treatment and distinct Toffoli/scratch exponents agree with the implementation source.
- Public (65,60,30), ten-bit metric gates reproduce: error_1=0.01391173198501434, error_2=0.019980679640672218, Lambda=290.7542054993138, E_gap=3.9159356215095613, lambda_eff=47.55854233515931.
- Released gauge drift 7.77156e-16; 25 circuits build within Q128. Static projection reproduces angle Toffolis 56450, total 65244, Q661, conditional cost 2051023203.3721037, gain 1.6154555940004633. No circuit execution correctness or new verified score was claimed or tested.
- Official research baseline's five factor arrays exactly match the sample. Starter-file hashes match. The original-release blob IDs match local source bytes after removing the download wrapper's one trailing newline; the public notebook JSON source string matches its stated SHA256. Attributed original design counts/contracts are supported by the pinned sources. The research paper uses correlation-energy selection and reports FeMoCo (15,57,19) Frobenius error 0.0825 Ha, supporting the documented non-transferability to this fixed-shift 0.020-Ha gate: https://arxiv.org/html/2502.15882v1 . Local GF(8)-based fragment compilation reproduces 42 multiplication and 210 inversion-phase Toffolis.
- Installed documentation contracts and live probes agree for Eigensystem row vectors, absolute-magnitude eigenvalue ordering, explicit Frobenius norms, little-endian BinaryReadList options, and PolynomialMod coefficient reduction. No notebook deliverable or speed claim is present.

FINDINGS
[P3] Handle an empty admissible scratch-exponent search before arithmetic. RankCost.wl:32-34 accepts PreparationCost[{2^75,1,1},1], but no exponent satisfies the strict scratch budget. SelectFirst returns Missing["NotFound"], which then contaminates an apparently successful Association with unevaluated costs. The official heuristic_costs(76,2**75,1,1,1) instead raises ValueError("Could not find k2_opt that lowers temp inner-qroam costs below rotation persistent ancilla"). This is outside the current sweep, so the stored numerical results remain valid. LookupCost at lines 53-59 should also propagate the preparation failure rather than form resource arithmetic around it.

RECOMMENDED FIXES
Return a structured Failure immediately when the exponent search finds no integer, and propagate that Failure from LookupCost. Add one failure-regime check for this no-admissible-exponent input. No other reproducible correctness finding remains.

OPEN ISSUES: 1
