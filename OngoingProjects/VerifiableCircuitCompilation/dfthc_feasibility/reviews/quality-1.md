# Independent WL quality review, round 1

Reviewed the complete target code/prose and loaded every field of the JSON artifacts. No applicable AGENTS.md was present in the project ancestry. No prior review report was read. All execution used an isolated temporary copy, including the local matrix and spectrum cache. Target files were not edited.

## A. Documented symbols — PASS

The WL targets use core language symbols; no paclet or QuantumFramework API needs an unsupported-signature exception. Installed Wolfram 15.0 Eigensystem documentation states that eigenvectors occupy rows and satisfy m.Transpose[vectors] == Transpose[vectors].DiagonalMatrix[values]; run.wls:11–13 uses and checks precisely that convention. Installed Norm documentation distinguishes the matrix operator norm from the explicit Frobenius norm used at run.wls:12–13. Installed BinaryReadList documentation supports the Real64/ByteOrdering option contract at run.wls:9. Installed PolynomialMod documentation supports coefficient reduction modulo 2 at RankCost.wl:70. Documentation of numerical eigenvalue ordering supports the explicit descending positive-value sort at RankCost.wl:15.

Read the official costs_prepare.py implementation and verify.py:419–518: the separate scratch exponent, power-of-two B condition, persistent workspace, and shared temporary maximum match the implementation. The isolated run crosschecked all 482 rows. The released loader's hash-pinned source contracts and official resource_costs/evaluate_words/score_from_metrics calls were exercised successfully.

## B. Idiomatic functional WL — FINDING

**One concrete defect:** RankCost.wl:78 maps a pure composition of listable bit operations over all 512 words. A Wolfram 15.0 kernel returned these attributes:

- BitXor: {Flat,Listable,OneIdentity,Orderless,Protected}
- BitShiftRight: {Listable,Protected}
- BitAnd: {Flat,Listable,OneIdentity,Orderless,Protected}

The smallest native rewrite is:

```wolfram
BitXor[BitShiftRight[signedWords,9],
  BitAnd[BitShiftRight[signedWords,8],1]] === ConstantArray[0,512]
```

The independent kernel probe returned True for this direct array expression, so removing Map preserves the full exhaustive invariant and lets the language thread the bit operations. This is one low-impact idiom issue, not a failed physics result. No imperative list construction, needless index-Table, or branching cascade was found. The many named dependent intermediates in PreparationCost serve the official resource arithmetic.

## C. Terse physics reasoning — PASS

DOCUMENTATION.md:7–27 develops the physical argument in order: orthonormal symmetric-pair packing, shared-projector span, PSD reconstruction, then the unrestricted spectral lower bound. It explicitly distinguishes that necessary condition from DFTHC representability. RankCost.wl:25–50 names the state-preparation and shared-workspace quantities directly. The public-source audit clearly labels static resource projection and required future circuit verification (AUDIT.md:84–92), and does not claim a new score.

## D. Physics invariants — PASS

The isolated run.wls execution completed successfully in Wolfram 15.0.0 and reproduced all 23 symbolic/edge checks. It obtained eigenpair Frobenius residual 3.624148053680244e-14, orthogonality residual 2.5268021680905354e-12, and independent SciPy spectrum difference 1.3048146537550791e-14; all are below the explicitly enforced guard/100.

The independent prepare.py run reproduced exact target pair symmetries and full N^4 residual 0.018639570630973906 versus packed residual 0.0186395706309741. Every field of its inputs.json matched the checked-in artifact. The independent loader run reproduced maximum vector drift up to row sign 7.771561172376096e-16, both official metric gates, all 25 Q128 circuit profiles, total T65244/Q661, cost 2051023203.3721037, and gain 1.6154555940004633. Every field of released-loader.json matched. No wavefunction normalization or many-body conservation law is implied by this real coefficient-tensor audit.

Additional exact kernel probes gave True for Gram quadratic-form positivity, its rank-two construction, an exact orthogonal change of basis, and the equality of Frobenius and spectral-tail error for a nondiagonal target with spectrum {5,3,-2}. The unrestricted PSD rank threshold reproduced as 290, with rank-289 error 0.020187308109326476 and rank-290 error 0.019858425813203292. The paper shape's rank cap 285 remains rejected.

## E. Nontrivial examples — PASS

The potential trivial case is a positive one-dimensional or diagonal-only target. The artifacts extend beyond it with symbolic mixed signs, degeneracy, arbitrary coefficient/projector identities, the physical 2926-dimensional shifted tensor, a genuine rejected shape, and a real Q128 source-built loader. An independent exact nondiagonal orthogonal conjugate of diag(5,3,-2) reproduced sqrt(13) for both rank-one truncation error and PSDTail, confirming the argument survives basis rotation.

## F. Symbolic-first derivations — PASS

RankCost.wl:65–82 keeps symbolic matrix entries, vector entries, Boolean table values, and the decoupling parameter. Its exact special spectra include radicals and Pi. DOCUMENTATION.md:18–27 supplies the arbitrary-dimension spectral argument. Machine diagonalization is reserved for the actual float32-rounded integral reference and checked independently; fixed orbital/precision constants belong to the expressly task-specific N76/chi10 resource model.

## G. Regimes and limits — PASS

The required table is present at DOCUMENTATION.md:66–74. Executed evidence covers the arbitrary-symbol identities and general spectral proof; exact mixed-sign case; decoupling/negative floor; official N^4 residual, full independent 2926 spectrum, and all 482 cost rows; and invalid, empty, degenerate, exact-radical/constant spectra plus impossible scratch/failure propagation. The independent cap probes additionally returned angle_qubits=128/within_cap=True for LookupCost[{65,60,30},15,7,1] and angle_qubits=143/within_cap=False for the same arguments with batch 8. The real released loader also exercises the Q128 boundary. Full circuit verification is explicitly outside the claimed result, so its deferred measurement/dirty-target coverage is not an unverified claim of this audit.

OPEN ISSUES: 1
