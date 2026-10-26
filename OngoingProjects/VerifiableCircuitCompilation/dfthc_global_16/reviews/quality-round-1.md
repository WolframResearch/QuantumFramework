# Independent WL quality review — round 1

The mathematical study passes the seven quality axes below. This is a quality verdict on a bounded, accuracy-rejected feasibility study: it does not establish a feasible rank-(16,32,20) candidate or a circuit-verified improvement.

Reviewed all named source, method, input and result artifacts. Parsed the complete results object, including all 3,600 history values; every history was finite and decreasing. Checked the exact applicable ancestor instruction paths; the existing ancestor AGENTS.md was empty. No prior reviews or other study material were consulted. Evidence-generating commands ran against isolated copies of the named sources and dependencies, with the named large arrays linked read-only in practice. No optimization, circuit execution or incumbent mutation was performed.

References below use repository-relative paths. Set `D=OngoingProjects/VerifiableCircuitCompilation/dfthc_global_16` when reading the abbreviated `D/...` references.

Reproduced commands in the isolated copy:

- `WolframKernel -noinit -noprompt -script symbolic.wls`: exit 0, PASS, all **53** checks true.
- `python -m unittest test_refit -v`: exit 0, all **9** tests pass.
- `WolframKernel -noinit -noprompt -script verify-data.wls`: exit 0, PASS, all **11** checks true; five complete resource rows checked.

The runtime was Wolfram Language 15.0.0 for Mac OS X ARM (February 17, 2026). Original source hashes and the results digest still matched the recorded manifests at the end of the review.

## A. Documented symbols — PASS

The non-core WL APIs are defined in the supplied packages, rather than recalled from another library. `D/GramCores.wl:2–21` supplies usage and definitions for `AugmentedCore[w,v]`, `ProjectorGram[u,v]` and `ProjectorGradient[u,w,u0,w0]`; lines 23–120 define the three nullary check functions. `D/GlobalRefit.wl:2–82` defines its nullary check function. The resource dependency `OngoingProjects/VerifiableCircuitCompilation/dfthc_feasibility/RankCost.wl` supplies the usage and definition of `LookupCost[shape,beta,J,lambda]` and its delegated `PreparationCost`.

Live `Options` for the eight exported functions called directly by the target returned eight empty lists. Executable representative calls returned augmented-core dimensions `{3,3}`, rectangular projector-Gram dimensions `{2,3}`, projector-gradient dimensions `{2,2}`, and an Association for `LookupCost[{16,32,20},16,6,1]`. The reproduced symbolic script confirmed Association-valued check functions and 53 Boolean results; the reproduced data script exercised the actual resource shapes and legal loader batches. No QF APIs occur. No unsupported signature, return-shape claim or unnecessary replacement for a documented built-in was found.

## B. Idiomatic functional WL — PASS

The substantive WL uses `Map`, `MapThread`, `AllTrue`, `Outer`, array differentiation, `GroebnerBasis`, `PolynomialReduce`, `Reduce` and `Minimize`. There are no `Do`, `For`, `While` or `AppendTo` constructs in the target WL. The `Table` in `D/GramCores.wl:32` enumerates principal-minor sizes and is a natural indexed construction, not a list-accumulation loop. The `Module` scopes also isolate the symbolic atoms and array heads needed for repeatable algebra; they are not wrappers around constant numeric results. Dense two- and three-dimensional symbolic fixtures are appropriate here. The numerical coordinate sweeps live in Python and carry sequential global-residual state; their structure does not imply an avoidable imperative WL implementation.

## C. Terse physics reasoning — PASS

`D/GlobalRefit.wl:13–58` follows the physical chain from Gram residual to cross-block weight/projector derivatives, constrained frame, polar adjoint, core projection, allocation and LCU chamber gradient. Intermediate names such as `residual`, `tangent`, `compressed`, `coreLoss`, `probability` and `chamberNorm` make those objects identifiable. `D/README.md:23–115` states the assumptions and equations needed to interpret the code. The added numerical machinery computes the full residual and actual decoded representation, rather than repeating independent factor fits. No material dead definition, duplicated calculation or prose/code mismatch requiring a fix was found.

## D. Physics invariants — PASS

The exact run independently confirmed core block reconstruction, offset contraction/norm, PSD quadratic form, orthogonal gauge, Schur range condition, frame isometry, polar adjoint and tensor/core projection identities (`D/GramCores.wl:49–71,93–120`; `D/GlobalRefit.wl:62–82`). A separate real exact fixture formed `A=Join[{{1,-2},{3,4}},{{-1,2}}]` and returned `True` for both symmetry and `PositiveSemidefiniteMatrixQ[A.Transpose[A]]`. The factor shape also gives augmented-core rank at most 20 directly.

Separately reconstructed every saved floating and quantized factor set, rather than accepting the report fields. With `x=einsum('rbc,rbi,rbj->rcij',w,u,u)`, `a=concatenate([w,v[:,None,:]],axis=1)` and `k=a@a.swapaxes(-1,-2)`, checked unit rows, real symmetric generators, symmetric PSD augmented cores, packed tensor reconstruction, global offset contraction and offset norm. All seven NPZ digests matched. Across those arrays:

- Maximum row-norm error: `4.45e-16`; maximum generator symmetry error: `6.94e-17`.
- Core symmetry error: zero; minimum computed core eigenvalue: `-5.72e-14`, consistent with rounding of the exact `A A^T` form.
- Maximum packed reconstruction discrepancy: `1.32e-15`; direct pair-residual errors agreed with saved values within `4.2e-15` Ha.
- Floating offset-matrix and norm drifts stayed below `2.29e-13`.
- Quantized offset drift ranged from `0.23976531874861534` to `0.0037703714889413824` Ha; the norm drift stayed below `2.28e-13`. Thus the report correctly separates floating preservation from quantized drift (`D/README.md:101–108`; `D/RESULTS.md:31–33`).
- Comparing the saved beta-16 and gauged beta-16 arrays gave pair-tensor drift `1.0041e-14` Ha, zero native-offset difference, and identical two-body words. This supports the claimed gauge preservation (`D/refit.py:348–363`).

These checks cover the invariants applicable to this real factorization. Evolution unitarity and cleanup are outside the explicitly stated mathematical scope.

## E. Nontrivial examples — PASS

The base case would be one commuting projector or a scalar core. The artifacts go beyond it: the fully symbolic rectangular gradient has independent weights and unequal source/target row counts (`D/GramCores.wl:35–40,61–62`); the joint symbolic example has two interacting blocks (`D/GlobalRefit.wl:17–29,63–64`); and the mixed-sign PSD boundary returns exact minimum 5 (`D/GlobalRefit.wl:78–80`). The noncommuting projectors in `D/GramCores.wl:67` have exact commutator `{{0,12/25},{-12/25,0}}`, reproduced separately. The normalization fixture distinguishes equal Gram cores from equal LCU normalization (`D/GramCores.wl:45–46,66`). Full-map finite differences use a nontrivial 3-orbital, two-block random case, followed by actual 76-orbital candidate checks. No extra circuit complexity is needed to expose the claimed structure.

## F. Symbolic-first derivations — PASS

All WL coefficients of the repaired rectangular gradient remain independent symbols. The global Gram and chain-rule gradients are differentiated exactly before numeric fitting. Frame and polar relations are reduced algebraically under explicit orthogonality constraints; allocation retains symbolic logits and positive symbolic tau; LCU differentiation retains arbitrary sign-chamber parameters (`D/GlobalRefit.wl:13–77`). Rational projectors, exact minimization and exact limits supply special cases. Floating arithmetic is confined to the large reference problem, finite-difference checks, numerical optimization and word decoding. It does not substitute for those symbolic derivations.

## G. Regimes and limits — PASS

The required regimes table exists at `D/README.md:134–140`, with applicable assumptions at lines 28–55. Its evidence is supported as follows:

| Regime | Reproduced evidence |
|---|---|
| General symbolic | Rectangular independent weights, cross-block polynomial gradients, tight-frame identity, polar/Sylvester adjoints, fixed-frame core projection, allocation and LCU derivatives: all true in the 53-check run. |
| Exact special | Nonzero rational projector commutator, mixed-sign rank boundary and same-core/different-normalization fixture: exact results reproduced. |
| Limits | Identity rotation, vanishing projection error and radial polar checks pass in the target. Independently specialized the general rescaling identity at zero covariance scale: with `o={1,0,0}`, `v={0,1,0}`, symbolic `w={{a,b,c},{d,e,f}}`, `q=I-o o^T`, `pv=v v^T` and `A=pv+s(I-pv)`, substituting `s->0` returned exactly `{w.pv.Transpose[w],w.v,v.v}` for the core, offset and norm expressions. |
| Numerical reference | All nine tests pass, including complete-map finite differences, direct global loss and core monotonicity; direct saved-array reconstruction checks both floating fits and five 76-orbital candidates. The WL data check independently matches five complete actual-normalization resource estimates and recorded score arithmetic. |
| Failure/boundary | Existing domain tests pass. Additional calls using the implemented guards reject rank-deficient polar input, zero projector rows, singular dictionaries, positive offsets and indefinite offsets. The saved two 1200-iteration fits, four 300-iteration fits and one 150-iteration gauge fit all honestly report the cap and lack of optimizer success; every actual two-body metric fails the 20 mHa gate. |

For reproducibility of the additional edge probes, use `small=GlobalModel(eye(6),-eye(3),7.,(2,3,2))` in Python. Form each raw parameter vector with `concatenate([frame.ravel(),zeros(small.w_size+small.r)])`. A frame made from `vstack([eye(3),zeros((3,3))])` raises `ValueError: Zero projector row.` in `small.materialize(raw)`. The frame `eye(3)[[0,0,1,1,2,2]]/sqrt(2)` raises `ValueError: Initial projector dictionary is singular.` in `core_initialize(small,frame,1)`. `polar(zeros((6,3)))` raises the rank-deficiency error. Replacing the negative offset with `eye(3)` or `diag([-1,1,-1])` raises `LinAlgError` at Cholesky. These checks perform no new optimization.

The failed bounded fits remain negative feasibility evidence for the stated restricted family; neither the prose nor the machine-readable result turns them into an impossibility theorem, a valid candidate or a new competition score.

OPEN ISSUES: 0
