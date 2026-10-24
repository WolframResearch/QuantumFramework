# Independent DFTHC quality review — round 3 of 3

Read all thirteen named target artifacts, including every JSON field and row, without consulting prior reviews or editing targets. No applicable AGENTS.md was present in the repository or ancestor directories. This review is bound to the SHA256 manifest below. Live work used an isolated copy with ignored matrix/spectrum caches and the sibling incumbent summary. The copied prepare.py, run.wls and audit_loader.py all exited 0. Parsed inputs.json, comparison.json and released-loader.json reproduced exactly, field for field. The only independent-probe serialization adjustment was converting Options rules to strings before RawJSON export; it did not touch target code. The final probe run produced no messages.

A. PASS — Documented symbols.

RankCost.wl uses core WL symbols, with no paclet API. Read installed Wolfram 15.0 documentation contracts for Eigenvalues, Eigensystem, Norm, BinaryReadList and PolynomialMod. Eigensystem returns eigenvectors as rows, normalized for approximate matrices; numeric eigenvalues sort by magnitude. run.wls:11–15 uses the documented eigen-equation, explicitly takes Frobenius matrix norms, and algebraically sorts both spectra for comparison. RankCost.wl:13–16 independently sorts positive eigenvalues before tail truncation, avoiding magnitude-order mistakes for mixed signs. BinaryReadList[...,"Real64",ByteOrdering->-1] matches the installed option contract. Live Options[Eigensystem] and Options[BinaryReadList], exact Norm and PolynomialMod probes passed. The official Python costs_prepare.py and verify.py implementations directly support the state-preparation, both Givens passes, shared-scratch and scoring formulas used in RankCost.wl:21–61. The three released module SHA256s and public NPZ are pinned and passed the loader audit checks. All nine starter-kit source/data SHA256s match. All nine original Git blob IDs match the inspected source after removing the temporary fetch files' single extra terminal blank line. No unsupported signature or memory-only API claim remains. Smallest fix: none.

B. PASS — Idiomatic functional WL.

The source search for Do/For/While/AppendTo/Table/Switch/Quiet/Print found none in either WL target. PSDTail uses Select/Sort/Drop/Total with listable Re and powers (RankCost.wl:13–16). The official decrement search becomes SelectFirst over a finite reversed range (:33–34). qrom and LookupCost use With for dependent constants; the nontrivial preparation and symbolic-test scopes use Module for their local symbols. run.wls:21–35 uses AllTrue, Map, Select and MinimalBy. Guards are boundary checks rather than a repeated type-dispatch cascade. The 2926 reference and its full eigenbasis are genuinely dense; their residual calculations produce dense matrices. No meaningful imperative or storage defect found. Smallest fix: none.

C. PASS — Terse physics reasoning.

DOCUMENTATION.md:3–27 identifies the real coefficient tensor, orthonormal symmetric-pair packing, PSD reconstruction, shared projector span cap, and arbitrary-dimensional spectral minimization in that order. RankCost.wl names its costs and scratch resources for their physical roles; run.wls names its eigensystem invariants. DOCUMENTATION.md:29–46 and README explicitly separate necessary screens, conditional fixed-normalization costs, static released-loader projections and circuit verification. No dead physics scaffolding or misleading success claim found. Smallest fix: none.

D. PASS — Applicable invariants.

All 32 built-in checks returned True in Wolfram 15.0.0, including pair inner product/isometry, projector Gram and sign invariance, all four Boolean address cases and all 512 signed words (RankCost.wl:74–109). Independently reproduced full N^4 baseline Frobenius residual 0.018639570630973906 and packed residual 0.0186395706309741. Full 2926 eigensystem: eigen-equation Frobenius residual 3.624148053680244e-14, orthogonality residual 2.5268021680905354e-12, independent SciPy spectrum L2 difference 1.3048146537550791e-14. All 482 cost rows match with model["within_cap"] explicitly asserted True (run.wls:24).

Additional independent probe: an exact non-diagonal orthogonal conjugation of diag(5,3,0,-2) has eigenvalues {-2,0,3,5}; its retained rank-one PSD matrix has {0,0,0,5}; the exact squared residual is 13 and matches PSDTail. A separate N=3, R=2, B=3, C=4 construction has six independent projector coordinates, factor rank 6, reconstruction rank 6 and RankCap=6, with symmetric Gram reconstruction and minimum numerical eigenvalue 1.1619183780607638. This tests saturation of the B<C shared-span cap rather than only a scalar bound.

Public C30 metrics and signed-gauge metrics pass the official gates; all 25 released circuit profiles have width 128 and table_identity PASS. Static projection remains angle T=56450, total T=65244, total Q=661, cost=2051023203.3721037, gain=1.6154555940004633. These results do not verify circuit execution or measurement phase repair. State normalization/unitarity/many-body conservation tests are inapplicable because the target is a real coefficient-tensor/resource audit, not a propagated state or verified circuit. Smallest fix: none.

E. PASS — Nontrivial examples.

The hiding place would be a scalar or positive diagonal spectrum. The target goes beyond it with exact mixed signs and degeneracy (RankCost.wl:80–91), the arbitrary-dimensional spectral argument (DOCUMENTATION.md:18–27), the 2926 physical reference, twelve distinct shapes, a real rejection at {15,57,19}, and fixed public C30 factors. Independent exact non-diagonal and rank-six saturation probes above further confirm structure unavailable in the scalar case. Smallest fix: none.

F. PASS — Symbolic first.

Packing, projector identities, four-address interpolation and the signed word identity use exact arithmetic. The PSD special case retains a,b,c with a>=b>=0,c>=0, and the decoupling parameter epsilon is taken to zero symbolically (RankCost.wl:80–83). The arbitrary-dimensional PSD-tail argument derives the form before applying it to the large official machine-precision reference. Exact radicals, Pi and exact mixed ordering remain exact. The independent rotated example also remains exact. Numerical diagonalization belongs to the supplied 2926 reference rather than replacing a tractable symbolic derivation. Smallest fix: none.

G. PASS — Required regimes.

The explicit table at DOCUMENTATION.md:66–74 is supported by execution:

- General symbolic: packing/projector/Boolean identities plus arbitrary-dimensional spectral argument.
- Exact special: a,b,-c symbolic truncation and {5,3,-2}; independently rotated mixed-sign matrix.
- Limit: epsilon->0 leaves c; allowing full rank retains the negative-spectrum floor.
- Numeric reference: full N^4/pair residual, independent SciPy/WL eigenspectra and all 482 official cost comparisons. Minimum unrestricted PSD rank at 20 mHa is 290; ranks 289/290 give 0.020187308109326476/0.019858425813203292 Ha. Shape {15,57,19} has cap 285 and bound 0.021548084729907786 Ha, so REJECT.
- Edge/failure: zero/empty/degenerate spectra, exact radical/constant ordering, machine/high-precision real-valued complex inputs, malformed complex spectra, negative ranks, invalid shape, scratch-budget failure and its propagation, negative/zero normalization, power-of-two B, exact Q128 acceptance and Q133 rejection. Independent probes additionally reject infinity and a genuinely nonreal imaginary part of 10^-40, while a 60-digit mixed-sign tail agrees to better than 10^-50.

No required regime is missing. Full circuit verification is explicitly excluded from this audit and remains required before adoption; the work makes no new verified score claim. Smallest fix: none.

Reproduction evidence is retained in the isolated copy, run.log, independent.wls, independent.log and independent.json. All original target SHA256s were unchanged after review.

SHA256 binding:

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
