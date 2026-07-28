# ArrayUtilities Container-Candidates Survey

> **REVISED 2026-07-26 (section F):** the admission criterion changed to
> shape-introspection + materialization (compute-nativeness is a capability flag, not
> a gate). Several section A/B rejections are superseded - QuantityArray and
> TabularColumn turn out to be COMPUTE-NATIVE, Tabular/Dataset/ByteArray/EventSeries/
> DataStructure["DynamicArray"] admit as storage-only. See section F for the
> live-kernel re-probe evidence.

Synthesis of three live-kernel probe reports. All claims kernel-verified on
**WL 15.0.0 for Mac OS X ARM (64-bit), May 19 2026** unless noted. Tiers refer
to the ArrayUtilities three-tier model: EXPLICIT / LAZY-PARAMETRIC / SYMBOLIC.

---

## A. Verdict table

| Candidate | Tier admitted | Recognition pattern | Shape introspection | One-line rationale |
|---|---|---|---|---|
| InterpolatingFunction[...][t] | LAZY (reference impl) | `HoldPattern[_InterpolatingFunction[_]]` (+ `Derivative[__][_InterpolatingFunction][_]`) | `Head[expr]["OutputDimensions"]` - free, no materialization | Inert on symbolic t, one packed whole-array eval per substitution, derivatives free |
| ParametricFunction | LAZY (flagship two-stage) | `_ParametricFunction`, `HoldPattern[_ParametricFunction[__]]`, `HoldPattern[_ParametricFunction[__][_]]` | probe-solve once: `pfv[RandomReal[]]["OutputDimensions"]`, then cache | Two-stage laziness (params, then t), param-solves auto-cached, packed complex output |
| BSplineFunction / BezierFunction | LAZY (second-class) | `HoldPattern[(_BSplineFunction\|_BezierFunction)[args__]]` with arg count == `obj["Rank"]` | `Dimensions[obj["ControlPoints"]]` (last = output dim), `obj["Rank"]` = arity | Real vectors only - no complex, no matrix output; silent out-of-domain inertness |
| Pure Function (unapplied) | LAZY (generic fallback) | `_Function`, arity from arg spec | formal-symbol probe `Dimensions[f[\[FormalT]]]`, else numeric probe, else caller-declared | Covers analytic gates exactly; NOT inert when applied - must store unapplied; unpacked output |
| Array-valued Piecewise | LAZY (fallback) | array-valued `_Piecewise` in one symbolic var | **CORRECTED**: shape comes from the BRANCH VALUES, not from `Dimensions` (see correction note below) | Shape is free without materializing, but the claimed `Dimensions` privilege is false in 15.0 |
| TimeSeries / TemporalData | LAZY (second-class) + interop | `_TimeSeries \| _TemporalData`; applied: `HoldPattern[(ts_TimeSeries)[t_]]` | `ts["ValueDimensions"]` (unlisted but works) | Correct lazy semantics but v15 TabularColumn leakage and broken Map make it accept-and-convert, not store |
| LinearSolveFunction | LAZY-OPERATOR (only if a lazy-map tier exists) | `_LinearSolveFunction` | `ls["MatrixDimensions"]` | An operator (lazy matrix inverse), not an array - has shape but no elements |
| FunctionInterpolation | REJECT (constructor) | n/a | n/a | Silently corrupts array-valued input in 15.0; use `Interpolation[Table[...]]` instead |
| PermutationMatrix / LowerTriangularMatrix / UpperTriangularMatrix / BlockDiagonalMatrix | EXPLICIT | head whitelist (atoms in 15.0; `AtomQ` and `ArrayQ` both True) | `Dimensions` works without materializing | Dot/Transpose/Part native and structure-preserving; but `/.` does not penetrate (Normal first) |
| SymmetrizedArray | EXPLICIT (flagged) | `_SymmetrizedArray` | `Dimensions` works | Dot and Map silently densify to List; Normal is NOT packed; only TensorContract preserves structure |
| ToeplitzMatrix / HankelMatrix | (not containers) | n/a | n/a | Plain constructors returning packed Lists in 15.0 - handled by the ordinary List path |
| NumericArray | EXPLICIT (storage-only) | `NumericArrayQ` | `Dimensions`, `ArrayDepth` work | Transpose/Part/Reshape native, but ALL arithmetic and Dot silently stay unevaluated - Normal (free) before compute |
| QuantityArray | REJECT (strip-adapter) | `_QuantityArray` | `Dimensions` works | Elements are Quantity with `NumericQ -> False` - poisons amplitude numericity; strip via QuantityMagnitude |
| ByteArray | REJECT (lift-adapter) | `ByteArrayQ` | `Dimensions` works | Integer-only, no arithmetic; lift via `NumericArray[ba]` on ingest |
| Around / CenteredInterval arrays | ELEMENT CLASS inside List | per-element head test | container is plain List | Not a container tier - Dot/Norm propagate correctly, but NumericQ/packed probes wrongly reject |
| Compress string | REJECT | `_String` (danger) | `Dimensions[str] -> {}` - misclassification trap | No compressed-array object in 15.0; a naive Dimensions probe reads it as a scalar |
| Tabular | REJECT / interop | `_Tabular` | `Dimensions`, RowCount/ColumnCount | No arithmetic, no Dot/Transpose, no laziness in v15; ingest via Normal + pack |
| Dataset | REJECT / interop | `_Dataset` | `Dimensions`; free type via ``Dataset`GetType`` | No arithmetic/Dot; Normal returns packed - cheap unwrap on ingest |
| DataStructure vectors | REJECT | n/a | n/a | Method-call-only access, mutable reference semantics, rank-1, unpacked extraction |
| EventSeries | REJECT (fold into TimeSeries ingest) | n/a | `"ValueDimensions"` works | Exact-time lookup only, no interpolation - not a lazy array |
| TimeSeriesModel | REJECT | n/a | n/a | Scalar forecasts only, not array-valued |
| MatrixSymbol / VectorSymbol / $Assumptions symbols | SYMBOLIC | `_MatrixSymbol` etc.; or `TensorDimensions` succeeds | `TensorDimensions` / `TensorRank` (NEVER `Dimensions`) | Shape flows through Dot symbolically; TensorExpand/Reduce give modest canonicalization |

---

## B. Per-candidate detail

### B.1 InterpolatingFunction (WL 15.0.0) - LAZY reference implementation

- Construction: NDSolveValue vector/matrix ODEs, or `Interpolation[Table[{t, matrix}, ...]]` including complex matrices.
- Shape: `ifn["OutputDimensions"]` -> `{2}` / `{2,2}` / `{}`; `ifn["Domain"]` -> `{{0.,1.}}`. Caveat: in 15.0 these are UNLISTED in `ifn["Properties"]` (which shows only ParameterAddress/Properties/Region/Unpack) but all still work. For the applied form use `Head[expr]["OutputDimensions"]`.
- Inertness: `ifn[tt]` stays inert on a symbolic argument. Eval trigger is NumericQ - `ifn[1/2]` and `ifn[Pi/6]` DO evaluate; only non-numeric symbols stay lazy.
- Evaluation: one whole-array call; real case packed. Packability of complex output depends on provenance: NDSolve-produced complex IF returns packed Complex, Interpolation-from-complex-data returned UNPACKED.
- Batch fast path: `ifn[{0.1, 0.2}]` threads over a time list in one call.
- Footguns (verified): `ifn[tt][[1]]` returns `tt` (Part hits the expression tree); `Map[f, ifn[tt]]` gives `ifn[f[tt]]`. `Dimensions[ifn[tt]]` -> `{1}` (arg count, garbage). Container must intercept Part/Map/Dimensions.
- Derivatives: `ifn'` is itself an InterpolatingFunction, array-valued, inert on symbolic arg.
- `Normal[ifn]` is a no-op.

### B.2 ParametricFunction (WL 15.0.0) - LAZY flagship, two-stage

- Construction: `ParametricNDSolveValue[..., {a}]` (multi-param and complex/Schrodinger RHS verified).
- Introspection: `"Properties"` -> Creator, DependentVariables, Expression, IndependentVariables, Parameters, Properties, TooltipTable. NO shape property on the object itself; `Dimensions[pfv]` -> `{7}` is garbage (internal expression parts).
- Shape protocol: one cheap probe solve `pfv[RandomReal[]]["OutputDimensions"]`; solves are cached internally ("Cache"->True confirmed; repeat call at same param took 3e-6 s), so cache the shape in the container at construction.
- Inertness: `pfv[aa]` and `pfv[aa][tt]` fully inert (LeafCount 214, cheap to store); `N[...]` does not force. Warning: exact numerics like Pi trigger a real solve (NumericQ trigger).
- Evaluation: `pfv[aa][tt] /. {aa -> 1., tt -> 0.5}` is ONE whole-array eval returning a packed array (packed Complex verified). `pfv[1.]` yields an InterpolatingFunction (list-of-functions form yields a list of IFs). `pfv[1.]'[0.5]` works.
- Same expression-tree footguns as IF on the inert forms; Transpose/Dot wrap symbolically without threading.

### B.3 BSplineFunction / BezierFunction (WL 15.0.0) - LAZY second-class

- Shape without evaluation: `Dimensions[bs["ControlPoints"]]` (output dim = last element); `bs["Rank"]` = input arity; `bs["Domain"]` works though UNLISTED.
- Hard limits (verified): complex control points NEVER evaluate (stay inert silently); a list of matrices is reinterpreted as a rank-2 spline surface (arity 2), so matrix-valued output is impossible. Output is only ever a real vector.
- Footgun: out-of-domain numeric args also stay inert - silent, no message. Guard with a domain check.
- Derivatives work numerically (`bs'[0.3]` verified). Exact args (1/3, Pi/12) evaluate.
- Do not route amplitudes here; usable for real parametrized curves (Bloch paths, pulse envelopes).

### B.4 Pure Function and array-valued Piecewise (WL 15.0.0) - LAZY fallback

- Function is NOT inert under symbolic application: `f[zz]` immediately expands to the elementwise tree. Store the UNAPPLIED Function; arity is introspectable from the arg spec (`f[[1]]`).
- Shape protocol (in order): `Dimensions[f[\[FormalT]]]` works for dimension-transparent bodies (elementwise math), fails for If-bodies; numeric probe at a domain-interior point is hazardous (`ConstantArray[1/u,{2,2}]` at 0. -> ComplexInfinity matrix); else require caller-declared shape.
- Evaluation is exactly one closure call (side-effect counter verified) but the result is UNPACKED - repack with Developer`ToPackedArray.
- Piecewise privilege (verified): `Dimensions` works on the unevaluated array-valued Piecewise (WL threads through branches). `If` does NOT get this (`Dimensions` -> `{3}`, arg count).

### B.5 TimeSeries / TemporalData (WL 15.0.0) - LAZY second-class + interop

- Vector-, matrix-, and complex-valued all construct and interpolate whole-array (`mts[0.5]` -> interpolated matrix).
- Shape: `ts["ValueDimensions"]` (works though unlisted in v15 Properties). `"PathLengths"` does NOT exist - returns silently inert `obj[PathLengths]`; unknown properties never error (footgun). `Dimensions[ts]` is the time length, not value shape.
- `ts[t]` fully inert on symbolic t; `ts[0.25]` returns a packed vector in one call.
- v15 hazards (all verified): `ts[{t1,t2}]` returns a raw TabularColumn, not a packed matrix; `ts["Values"]`, `ts["Path"]`, and `ts[[1]]` leak TabularColumn objects (unpacked); `Map[Total, ts]` is a silent no-op; symbolic values fall back to previous-value semantics (NOT linear interpolation); `ts["PathFunction"]` returns a TimeSeries, not an InterpolatingFunction - no cheap unwrap.
- Arithmetic (`ts+ts`, `2 ts`) is container-preserving and correct. TimeSeriesResample works.
- Verdict rationale: accept on input (convert), offer on output; keep InterpolatingFunction as the canonical internal lazy carrier.

### B.6 LinearSolveFunction (WL 15.0.0) - LAZY-OPERATOR

- `ls["MatrixDimensions"]` is the shape API (dense and sparse verified); `ls["Dimensions"]` stays inert and `Dimensions[ls]` -> `{2}` is garbage.
- `ls[vv]` (symbolic atom) inert, but `ls[{x1,x2}]` evaluates symbolically elementwise - a List argument is always consumed.
- Matrix RHS works (`ls[IdentityMatrix[2]]` = inverse); `Normal` is a no-op; getL/getU expose factors.
- Has shape but no elements: admissible only if ArrayUtilities adds a "lazy linear map" notion alongside "lazy array".

### B.7 FunctionInterpolation (WL 15.0.0) - REJECT

- Verified broken for array-valued expressions: vector input produced a SCALAR IF (value = Cos branch only, ncvb messages); 2x2 matrix input produced OutputDimensions `{2}` equal to just the first row. Silent truncation.
- Replacement: `Interpolation[Table[{t, f[t]}, ...]]` - verified correct for complex matrices.

### B.8 StructuredArray family (WL 15.0.0) - EXPLICIT

- Version-critical: in 15.0 these are atoms whose Head is the constructor symbol itself (LowerTriangularMatrix, SymmetrizedArray, PermutationMatrix, BlockDiagonalMatrix, QuantityArray), `AtomQ -> True`, internal form `head[StructuredArray`StructuredData[dims, data]]`. No `StructuredArrayQ` exists. Recognition = head whitelist.
- Critical trap: ReplaceAll does NOT penetrate (atomic) - `LowerTriangularMatrix[{{t,0},{2t,3}}] /. t -> 5` leaves t inside. Substitution must be Normal-first (repacks fine). So these cannot serve as lazy-parametric carriers even though symbolic entries are allowed inside.
- Structure preservation under Dot (verified): PermutationMatrix, LT/UT (Transpose converts LT <-> UT), BlockDiagonalMatrix all preserve; SymmetrizedArray is silently DESTROYED to a plain List by Dot and Map, and its Normal is not packed; TensorContract[TensorProduct[sa,sa],...] DOES stay SymmetrizedArray.
- Toeplitz/Hankel return plain packed Lists in 15.0 - nothing to recognize.
- QuantityArray: Dot works (`qa.qa` -> Quantity) but elements are non-NumericQ - reject for amplitudes, strip units on ingest.

### B.9 NumericArray (WL 15.0.0) - EXPLICIT storage-only

- Weaker than folklore: `ArrayQ`/`MatrixQ` -> False, TensorRank unevaluated; `Dimensions`/`ArrayDepth`/`Part` work.
- Native (stays NumericArray): Transpose, ConjugateTranspose (correct on ComplexReal64), Flatten, ArrayReshape, Part, Total.
- Silently UNEVALUATED (worst failure mode - no error, symbolic residue): Dot in every combination, Plus, scalar Times, Sin, Map, KroneckerProduct, Norm. Dot does NOT auto-convert.
- No symbolic story: construction with a symbol fails.
- Normal is essentially free (0.0004-0.0005 s round trip on 1000x1000 Real64, ByteCount identical to packed). Route all compute through Normal, re-wrap on request; never let a bare NumericArray reach Dot/Plus.

### B.10 Around / CenteredInterval element arrays (WL 15.0.0) - element class, not a tier

- Container is a plain (necessarily unpacked) List: `ArrayQ -> True`, `PackedArrayQ -> False`, elements `NumericQ -> False`.
- Propagation genuinely works: `av.av` -> Around, `Norm[av]` -> Around, `Total[Abs[av]^2]` -> Around, matrix.vec per component; CenteredInterval same.
- Design consequence: ArrayUtilities needs a per-element numericity notion ("uncertain-numeric") orthogonal to the container tier, or NumericQ/packability probes will wrongly reject these.

### B.11 Tabular / Dataset / DataStructure / EventSeries / TimeSeriesModel / ByteArray / Compress (WL 15.0.0) - REJECT

- Tabular: constructs, Dimensions works, Normal packs cheaply, but arithmetic/Dot/Transpose all stay inert; `tab[[All,1]]` returns TabularColumn; no lazy or out-of-core backing exists (100k x 5 reals fully in-core). Interop only.
- Dataset: free type introspection via ``Dataset`GetType`` (TypeSystem shape without traversal); Normal returns PACKED; Transpose/Map/queries work natively; but Plus/Dot inert. Interop only.
- DataStructure: $DataStructures has DynamicArray, FixedArray, BitVector, ExtensibleVector, ImmutableVector, StringVector (no "Array", no "SparseVector"); method-call access only, unpacked Normal, mutable reference semantics (aliasing hazard for value-semantics QuantumState). Reject.
- EventSeries: exact-time lookup only (`es[0.5]` -> Missing), no interpolation. Fold into TimeSeries ingest.
- TimeSeriesModel: scalar forecasts. Reject.
- ByteArray: Dimensions/Part/Normal work, no arithmetic, integer-only; lift via `NumericArray[ba]` (verified). Reject as container.
- Compress: Head is String; `Dimensions[string] -> {}` would misclassify as a scalar under a naive probe. At most at-rest serialization behind explicit Uncompress.

### B.12 Symbolic tier: MatrixSymbol / VectorSymbol / $Assumptions (WL 15.0.0)

- With `A,B in Matrices[{n,n}]`: `TensorDimensions[A.B]` -> `{n,n}`, `TensorRank[A.B]` -> 2; `TensorExpand[A.(B+A)]` -> `A.B + MatrixPower[A,2]`; `TensorReduce[Transpose[A.B]]` -> `Transpose[B].Transpose[A]`. Limitation: `TensorReduce[Tr[A.B] - Tr[B.A]]` does NOT reduce to 0.
- MatrixSymbol: TensorDimensions correct; Transpose inert; Normal does NOT materialize; NumericQ/ArrayQ False. Trap verified: `Dimensions[MatrixSymbol[...]] -> {2}` measures the expression - the symbolic tier MUST use TensorDimensions, never Dimensions.

---

## C. Recommendations for ArrayUtilities

### C.1 Admit to the LAZY tier NOW (v1)

1. **InterpolatingFunction applied form** - canonical lazy container. Recognition `HoldPattern[_InterpolatingFunction[_]] | HoldPattern[Derivative[__][_InterpolatingFunction][_]]`; shape via `Head[expr]["OutputDimensions"]`.
2. **ParametricFunction** (all three arities: object, `pf[params]`, `pf[params][t]`) - the two-stage flagship. Cache shape from one probe solve at construction.
3. **Unapplied Function** as the generic analytic fallback, with the three-step shape protocol (formal-symbol probe, numeric probe, caller-declared) and repack-after-eval.
4. **Array-valued Piecewise** in one symbolic variable (free Dimensions on the unevaluated form).

Cross-cutting v1 requirements (all verified as universal):
- Intercept Part/Map/Dimensions on EVERY inert applied form - they operate on the expression tree, not the array. Either force evaluation or rewrite the parameter.
- Document the NumericQ evaluation trigger: exact symbolic values (Pi, 1/2) do NOT stay lazy; only non-numeric symbols do. Inertness cannot represent "exact symbolic time".
- One-substitution semantics: a single `/. {params}` on the stored inert form is exactly one whole-array evaluation. Packed for IF/ParametricFunction; unpacked for Function and data-Interpolation-complex - repack downstream.
- Dispatch by head-of-head: a table on {InterpolatingFunction, ParametricFunction, BSplineFunction, BezierFunction, LinearSolveFunction, Function, Piecewise} covers all lazy survivors.

### C.2 Behind an extensible registry LATER (v2)

- **BSplineFunction / BezierFunction** - real-vector-only; needs the arity check (`Length[{args}] == obj["Rank"]`) and a domain guard against silent out-of-domain inertness.
- **TimeSeries / TemporalData** - accept-on-input via conversion to InterpolatingFunction, offer-on-output; do not store internally until the v15 TabularColumn leakage settles.
- **LinearSolveFunction** - only if/when a "lazy linear map" concept is added beside "lazy array"; recognition `_LinearSolveFunction`, shape `"MatrixDimensions"`.
- **EXPLICIT structured heads** (PermutationMatrix, LT/UT, BlockDiagonalMatrix) via head whitelist, with a Normal-first substitution rule; SymmetrizedArray with a densify-on-Dot flag.
- **NumericArray** as storage-only EXPLICIT: Normal-before-compute, re-wrap on request.
- **Element-class registry**: Around / CenteredInterval as "uncertain-numeric" elements inside plain Lists - a numericity axis orthogonal to the container tier.
- **SYMBOLIC tier**: probe with TensorDimensions/TensorRank (never Dimensions/NumericQ); canonicalize with TensorExpand/TensorReduce, knowing their limits.

### C.3 Reject

- FunctionInterpolation (data corruption on array input), Tabular, Dataset, DataStructure vectors, EventSeries, TimeSeriesModel, ByteArray, QuantityArray (for amplitudes), Compress strings, Toeplitz/Hankel (nothing to recognize).
- Provide cheap ingest adapters where verified: Dataset/Tabular -> Normal + pack; ByteArray -> `NumericArray[ba]`; QuantityArray -> QuantityMagnitude; EventSeries -> TimeSeries path.

### C.4 Regression tests to pin

The following work in 15.0 but are UNLISTED/undocumented API - pin tests so a future version break is caught:
- IF `"OutputDimensions"`, `"Domain"`, `"ValuesOnGrid"` (absent from `"Properties"`).
- BSpline/Bezier `"Domain"`.
- TimeSeries `"ValueDimensions"`.
- Piecewise Dimensions-threading on the unevaluated form.
- Structured-array atomicity (Head = constructor symbol) and the ReplaceAll opacity that follows from it.
- Dataset Normal returning packed; NumericArray Normal round-trip cost.

---

## D. QuantumFramework-specific opportunities

1. **QuantumEvolve with symbolic Hamiltonian parameters**: store `ParametricFunction[...][a][t]` whole in the state container. Parameter binding is a single ReplaceAll doing one whole-array solve + interpolation eval; per-parameter solves are cached by the kernel automatically. This is the natural backend for parameter-sweep states and variational-circuit scans.
2. **Trajectory sampling fast path**: `psi[{t1, t2, ...}]` on an InterpolatingFunction returns the matrix of states in one call - use for QuantumEvolve trajectory output.
3. **Derivative containers for free**: `ifn'` is itself an array-valued InterpolatingFunction - d(psi)/dt, Heisenberg-picture operator derivatives, etc., need no new machinery.
4. **Analytic gates without interpolation error**: unapplied Functions cover `t |-> {{Cos t, -Sin t}, {Sin t, Cos t}}` exactly; document repacking and possible undiscoverable shape.
5. **BlockDiagonalMatrix** as a carrier for direct-sum operators and superselection-sector density matrices; **PermutationMatrix** for qudit-permutation gates (O(n) storage, Dot preserves structure).
6. **LinearSolveFunction** as a QuantumOperator backend representing H^-1 or solving (I - K)x = b in steady-state/Lindblad problems without forming an inverse.
7. **TimeSeries as interop export**: it carries time-domain metadata (FirstTime/LastTime/ResamplingType) that the basis currently stores as ParameterSpec - offer QuantumEvolve results as TimeSeries on output, accept on input, convert internally.
8. **NumericArray "ComplexReal64"** for compact at-rest amplitude storage and GPU/LibraryLink boundaries (conversion is essentially free), never as a compute representation.
9. **Around/CenteredInterval amplitudes** for noise studies: Dot/Norm/Abs^2 propagate uncertainty out of the box at unpacked-List speed, once the numericity probe is fixed.
10. **BSpline curves** for Bloch-sphere paths and real pulse envelopes (not amplitudes).

---

## E. Open questions

1. **TabularColumn leakage scope (v15)**: is the TimeSeries rebuild on Tabular internals a transitional regression (Values/Part/list-of-times leaking TabularColumn, Map a no-op) or the intended long-term surface? Determines whether TimeSeries can ever be promoted from interop to internal storage.
2. **Complex-IF packability by provenance**: NDSolve-produced complex IFs return packed, Interpolation-from-data ones do not. Is there a construction option or post-hoc repack that normalizes this, or must the container track provenance?
3. **Lazy-operator tier**: should ArrayUtilities grow a first-class "lazy linear map" alongside "lazy array" (LinearSolveFunction today; potentially MatrixFunction/affine maps later), or does that belong in QuantumFramework's operator layer?
4. **Shape declaration protocol**: for Function bodies where both formal-symbol and numeric probes fail (If-bodies, singular points), what is the API for caller-declared shape, and is it required or optional at construction?
5. **Exact-time inertness**: the NumericQ trigger means Pi cannot stay lazy. Is a wrapper (e.g. holding the argument) worth the interception cost to represent exact symbolic times, or is "non-numeric symbols only" an acceptable documented contract?
6. **Structured-array substitution**: should the EXPLICIT tier auto-Normal structured atoms when a ReplaceAll carries parameter bindings, or refuse and require explicit materialization? (Silent no-substitution is the current kernel behavior and is a correctness hazard.)
7. **SymmetrizedArray routing**: is TensorContract-heavy code common enough in QuantumFramework (symmetric state spaces, bosonic sectors) to justify keeping SymmetrizedArray structured, versus normalizing to packed on ingest?
8. **Unknown-property silence**: TimeSeries (and structured objects generally) return inert expressions for unknown properties instead of erroring. Should the recognition layer wrap property access with a validity check?
9. **Version pinning**: all verdicts are 15.0.0-specific (atomic structured arrays, TimeSeries internals, unlisted properties). What is the minimum supported WL version, and do the 13.x/14.x behaviors (e.g. `StructuredArray` head) need parallel recognition paths?

---

## F. Revision under the shape-based admission criterion (WL 15.0, re-probed)

### F.1 Vector wrappers

MINI-REPORTS — WL 15.0.0 (Mac ARM), wolframscript live-kernel, all snippets verified this session.

---

**1. QuantityArray — REVISED VERDICT: compute-native explicit**

1. Shape without materializing: `Dimensions[qa]` → `{3}`, `Length`, `ArrayQ[qa]` → True; rank-3 `Dimensions[QuantityArray[..., {4,5,6}, "Kelvins"]]` → `{4,5,6}`. `AtomQ[qa]` → True (atomic wrapper, but Dimensions still works).
2. Materialization: `QuantityMagnitude[qa]` — effectively free (1e6 elements: 1e-5 s; it hands back the internal packed array) and IS packed. WARNING: `Normal[qa]` is the wrong path — 1e5 elements took 0.314 s, returns unpacked list of `Quantity` objects. Adapter must use `QuantityMagnitude`, never `Normal`.
3. Reconstruction: `QuantityArray[mags, "Meters"]` — 3.5e-5 s at 1e6 elements. Unit metadata for rebuild: `qa["UnitBlock"]` → `Meters` (scalar for uniform, `{Meters, Seconds}` for per-column mixed); mixed round-trip `QuantityArray[newMags, mixed["UnitBlock"]]` verified → QuantityArray. (`QuantityUnit[qa]` also works but materializes a full unit array — use `"UnitBlock"`.)
4. Arithmetic (all verified, old survey wrong to reject): `qa+qa` → QuantityArray `{2.,4.,6.}` m; `2 qa` → QuantityArray; `qa.qa` → `Quantity[14., Meters^2]`; `Total[qa]` → `Quantity[6., Meters]`; `Mean` works; `qa*qa` → Meters² QuantityArray. Capability flag: full native arithmetic incl. unit algebra.
5. Numericity: numeric at module level (numeric-with-units). Complex magnitudes verified: `QuantityArray[{1.+2.I, 3.-I}, "Meters"]` constructs, Dimensions `{2}`, magnitudes round-trip exactly.
6. Recognition: `MatchQ[x, _QuantityArray]` → True. No `QuantityArrayQ` exists (`NameQ` → False). Pattern: `_QuantityArray` (plus optional `ArrayQ` guard).
7. Adapter recipe: shape = `Dimensions`; extract = `QuantityMagnitude[#]&` (packed, O(1)); rebuild = `QuantityArray[newData, #["UnitBlock"]]&`; compute natively when both operands share the wrapper, else compute on magnitudes and rebuild.

---

**2. TabularColumn — REVISED VERDICT: compute-native explicit** (biggest reversal)

1. Constructor: `TabularColumn[{1.,2.,3.}]` works directly → `TabularColumn[<|Data->{{1.,2.,3.},{},None}, ElementType->Real64|>]`; typed form `TabularColumn[{1.,2.,3.}, "Real64"]` also works. Leak reproduced: `TimeSeries[{1.,2.,3.},{1,3}]["Values"]` returns a TabularColumn (identical structure). (Aside: `Tabular[rows, {"Columns"->...}]` with that option syntax failed to construct; irrelevant to the column type.)
2. Shape: `Dimensions[tc]` → `{3}`, `Length[tc]` → 3, without materializing.
3. Materialization: `Normal[tc]` → `{1.,2.,3.}`, PACKED (`Developer`PackedArrayQ` → True) for clean numeric columns.
4. Arithmetic — the old survey never tested it; it is native: `tc+tc` and `2 tc` return TabularColumn (Data `{2.,4.,6.}`); `Total[tc]` → 6.; `tc.tc` → 14.; `Mean` → 2.; `Map[#^2&, tc]` → materializes to `{1.,4.,9.}` (List, not wrapper).
5. Missing-data semantics: `TabularColumn[{1., Missing[], 3.}]` constructs, ElementType stays Real64, `Dimensions` → `{3}`, `Normal` → `{1., Missing[], 3.}` (NOT packed), `Total` → 4. (Missing skipped — NaN-ignoring semantics). Numericity: numeric at module level; Missing-tolerant column is still numeric per revised criterion, with a "may contain Missing / unpacked Normal" caveat flag from `#["ElementType"]` + a Missing scan. Non-numeric columns exist too (`TabularColumn[{"a","b"}]` constructs; Normal → `{"a","b"}`) — numericity must be checked per instance via ElementType, not per head.
6. Recognition: `MatchQ[x, _TabularColumn]` → True; `AtomQ` → True. Properties: `tc["ElementType"]` → `Real64`.
7. Reconstruction: `TabularColumn[newList]` or `TabularColumn[newList, type]` — verified. Adapter: shape=`Dimensions`, extract=`Normal`, rebuild=`TabularColumn[#, tc["ElementType"]]&`, numericity gate on ElementType.

---

**3. ByteArray — REVISED VERDICT: storage-only explicit (integer vector)**

1. Shape: `Dimensions[ba]` → `{4}`, `Length` works — no materialization.
2. Materialization: `Normal[ba]` at 1e7 bytes: 7.4e-4 s, packed. Lift: `NumericArray[ba]` at 1e7: 6.3e-4 s → type `UnsignedInteger8`. String round-trip: `BaseEncode[ba]` → `"AQID/w=="`, `ByteArray["AQID/w=="] === ba` → True.
3. Reconstruction: `ByteArray[list]` (or via BaseEncode string) — trivial.
4. Arithmetic flag: NOT native — `ba+ba` and `2 ba` stay symbolic (`2*ByteArray[<4>]`), `ba.ba` unevaluated, `Mean[ba]` unevaluated. Exception: `Total[ba]` → 261 works natively. Flag: reductions-partial (Total only); materialize-before-compute otherwise.
5. Numericity: numeric (uint8 vector), integer-typed.
6. Recognition: `ByteArrayQ[ba]` → True; pattern `_ByteArray?ByteArrayQ`.
7. Adapter: shape=`Dimensions`, extract=`Normal` (packed) or lift `NumericArray[ba]`, rebuild=`ByteArray[ints]` (consumer must clamp to 0–255). Admit as storage-only under revised criterion.

---

**4. NumericArray (reconfirm anchor)**

`NumericArray[data, "ComplexReal64"]` + `Normal` round trip at 2^16 elements: **9.0e-5 s average** (100 iterations); Normal returns packed, equality holds, `Dimensions` → `{65536}`. Confirms the materialize-before-compute contract is ~free at this scale. Side note: `na+na` stays symbolic (`Times[2, na]` unevaluated wrapper), `na.na` unevaluated, but `Total[na]` computes natively — same "reductions-partial" flag as ByteArray. `NumericArrayQ` → True. Storage-only contract anchored.

---

**5. Association (index→value) — REVISED VERDICT: genuinely reject (no shape)**

Nuance the old survey missed: `Dimensions[<|1->1.5, 2->2.5, 5->-1.0|>]` does NOT fail — it returns `{3}`. But that is the entry count, not the logical vector shape (max index 5): shape introspection reports the wrong thing for a sparse-vector reading, and `ArrayQ` → False. `Normal` returns a list of RULES, not a dense vector — so there is no materialization path to the represented array either; a dense realization requires a consumer-supplied length policy (module can't infer trailing zeros past index 5). Arithmetic incidentally threads over values (`2 as`, `as+as`, `Total` → 3. all work on VALUES), which is exactly the wrong semantics for a sparse vector (ignores implicit zeros, e.g. Mean would be wrong). Both revised admission gates fail as-intended: no faithful shape, no faithful materialization. Reject stands — rationale updated from "Dimensions fails" to "Dimensions/Normal report the entry multiset, not the represented vector." Recommend users convert to `SparseArray[Normal[as]]` (rules→SparseArray) explicitly, supplying dimensions.

---

Probe scripts (rerunnable): `/private/tmp/claude-501/-Users-swish-src-wolfram-QuantumFramework/7e8c928b-f01c-4e17-8123-ca13235727bc/scratchpad/probe1_qa.wl`, `probe2_tc.wl`, `probe3_misc.wl`.

### F.2 Tabular / data / temporal

All probes complete (WL 15.0.0, wolframscript, May 2026 kernel). Per-type mini-reports against the REVISED criterion follow.

ENVIRONMENT: "15.0.0 for Mac OS X ARM (64-bit) (May 19, 2026)". Probe scripts: /private/tmp/claude-501/-Users-swish-src-wolfram-QuantumFramework/7e8c928b-f01c-4e17-8123-ca13235727bc/scratchpad/{probe_tabular.wls,probe_rest.wls,probe_followup.wls}; raw output: .../tasks/{bn5a4rgvi,bjn44849u}.output.

================================================================
1. Tabular — REVISED VERDICT: storage-only EXPLICIT (admit)
================================================================
(1) Shape WITHOUT materializing (all verified): `Dimensions[tab]` -> {100000,4} in ~0s; `Length[tab]` -> 100000; properties `tab["RowCount"]` -> 100000, `tab["ColumnCount"]` -> 4, `tab["Dimensions"]` -> {100000,4}; `ColumnKeys[tab]` -> None for anonymous columns, {"a","b","c","d"} for named. `TabularStructure[tab]` gives per-column ElementType (Real64), NonMissingCount/MissingCount/ByteCount, and `Backend -> WolframKernel`.
(2) Materialization: construction `Tabular[mat]` from a 10^5 x 4 packed matrix costs 1.29 s (the expensive direction). `Normal[tab]` on an ANONYMOUS all-numeric tabular returns a rank-2 List that is ALREADY PACKED — 0.000995 s at 10^5 x 4; Developer`ToPackedArray after is a no-op (5e-6 s). On a NAMED tabular Normal returns list of Associations (0.113 s); the cheap named path is per-column: `Normal[tab[[All,"a"]]]` -> packed vector in 0.00009 s, so matrix = Transpose of per-column Normals, or accept the 0.13 s Values/@Normal route. Column extraction `tab[[All, 1]]` and `tab[[All, "a"]]` return head `TabularColumn` (NOT List, NOT Tabular); `Length`/`Dimensions` work on it ({100000}); `Normal` on it is packed. Row extraction `tab[[2]]` returns head `TabularRow`. `Normal[tab, "Columns"]` returns a Tabular, not an assoc of columns — don't use it.
(3) Reconstruction: `Tabular[matrix]` and `ToTabular[matrix]` both work (verified, dims preserved); with column names use `Tabular[matrix, {"a","b","c","d"}]` (verified round-trip, keys preserved). `ToTabular[matrix, names]` FAILS ($Failed) — names go through Tabular, not ToTabular.
(4) Arithmetic flag: `tab + 1` -> unevaluated Plus; `2 tab` -> unevaluated Times; `tab . m` -> unevaluated Dot; `Total[tab]` fails; `Mean[tab]` unevaluated. BUT `Map[# + 1 &, tab]` -> Tabular (native structural map works). Flag: no native arithmetic/Dot; native Map only.
(5) Numericity: per-column typed (ElementType Real64/Integer64 via TabularStructure); Missing is representable per cell (verified: `Tabular[{{1, Missing["bad"]},{2, 3.5}}, {"x","y"}]` keeps `Missing["bad"]` through Normal, structure reports MissingCount). Module-level numeric verdict = all columns numeric ElementType AND MissingCount 0 per column — decidable from TabularStructure without materializing.
(6) Recognition: `TabularQ[tab]` -> True; `Head[tab] === Tabular`. `ArrayQ[tab]` -> False (so generic ArrayQ gates would wrongly reject it — recognize by head).
(7) Lazy/out-of-core in THIS version: `Options[Tabular]` / `Options[ToTabular]` are display options only (plus MetaInformation; ToTabular adds MissingValuePattern) — no backend/memory-map option. TabularStructure exposes a `Backend` slot, observed only as `WolframKernel`. Treat as fully in-memory here.
Adapter recipe: shape = Dimensions; ingest = Normal (anonymous) or per-column Normal@tab[[All,k]] then Transpose (named); egress = Tabular[matrix, ColumnKeys[original] /. None -> Sequence[]]; numeric gate = TabularStructure column types + MissingCount.

================================================================
2. Dataset — REVISED VERDICT: storage-only EXPLICIT (admit)
================================================================
(1) Shape: `Dimensions[ds]` -> {1000,3} (and {100} for rank-1) without traversal. Type without traversal: `Dataset`GetType[ds]` -> `TypeSystem`Vector[TypeSystem`Vector[TypeSystem`Atom[Real], 3], 1000]` — gives shape AND element numericity in one call (verified; no "ColumnTypes" property needed).
(2) Materialization: `Normal[ds]` returns a List that IS packed for uniform numeric rank-2 and rank-1 (reconfirmed: Developer`PackedArrayQ True both). Cost negligible (the Dataset wraps the packed array).
(3) Reconstruction: `Dataset[matrix]` verified, dims preserved.
(4) Arithmetic flag (old survey corrected): `ds + 1` -> unevaluated Plus; `ds . m` -> unevaluated Dot — NOT compute-native. However `Total[ds]` and `Mean[ds]` DO evaluate and return a plain List (head List, correct column stats). Flag: no elementwise/Dot; a few statistical reducers work but exit the container.
(5) Numericity: read off GetType (Atom[Real] leaves) — numeric at module level for uniform numeric datasets.
(6) Recognition: `Head[ds] === Dataset` plus GetType pattern `TypeSystem`Vector[TypeSystem`Vector[TypeSystem`Atom[Real|Integer],_],_]` for the uniform-matrix case.
(7) Adapter: shape = Dimensions; ingest = Normal (already packed); egress = Dataset[matrix]; materialize-before-compute always.

================================================================
3. DataStructure — REVISED VERDICT: storage-only EXPLICIT for "DynamicArray" (admit with reference-semantics caveat); rank-1 only
================================================================
Available array-like stores in 15.0 (verified by creation attempt): "DynamicArray" (yes), "FixedArray" (yes, needs length arg: CreateDataStructure["FixedArray", n]), "LinkedList" (yes), "ExtensibleVector" (yes). "SparseVector", "BitVector", "Array": DO NOT EXIST ($Failed).
Best candidate — "DynamicArray":
(1) Shape: `da["Length"]` and also top-level `Length[da]` (both verified, 10^5). 1-D only — no rank-2 store exists, so admissible for vectors only.
(2) Materialization: `da["Elements"]` 0.0015 s at 10^5 reals, returns UNPACKED List; `Normal[da]` 0.00036 s, also unpacked; `Developer`ToPackedArray[da["Elements"]]` packs successfully, total 0.0015 s. So: cheap but always add ToPackedArray.
(3) Reconstruction: `CreateDataStructure["DynamicArray", list]` accepts a whole list directly (verified, Length 10^5) — no element-by-element loop needed. FixedArray also accepts a list: `CreateDataStructure["FixedArray", RandomReal[1,100]]` -> Length 100 (note: this arg form doubles as "from list"; bare integer means length).
(4) Arithmetic: `da + 1` -> unevaluated Plus. None native.
(5) Numericity: untyped container — numeric iff contents numeric; must inspect elements (or trust producer). No type metadata.
(6) Recognition: `DataStructureQ[da]` -> True; `MatchQ[da, DataStructure["DynamicArray", ___]]` -> True (type name is the first arg — usable as recognition pattern).
(7) Reference semantics — CONFIRMED HAZARD, but not disqualifying under the revised criterion: `db = da; da["Append", x]` mutates through the copy (both lengths 100001); explicit `da["Copy"]` detaches (da 100002 vs copy 100001). Adapter recipe: recognize by DataStructureQ + type name; ingest = ToPackedArray@ds["Elements"] IMMEDIATELY at the module boundary (snapshot defeats aliasing); egress = CreateDataStructure["DynamicArray", newList]; never hold the live handle across user code.

================================================================
4. EventSeries — REVISED VERDICT: storage-only EXPLICIT time-indexed array (admit), with a native elementwise-arithmetic capability flag
================================================================
(1) Shape: `ev["ValueDimensions"]` -> 3 (reconfirmed), `ev["PathLength"]` -> 50, `Dimensions[ev]` -> {50}. No interpolation: `ev[25.]` returns the stored row exactly; `ev[25.5]` -> `Missing[NotAvailable]` (verified — genuinely discrete, NOT lazy/interpolating).
(2) Materialization — VERSION CHANGE FOUND: in this 15.0 kernel `ev["Values"]` returns a TabularColumn-backed object (Head TabularColumn; == against the original matrix still True, Dimensions {50,3}), NOT a plain list, and it is not packed as-is. `Normal[ev["Values"]]` -> plain PACKED rank-2 List equal to the input (verified). `ev["Times"]` -> packed vector. `ev["Path"]` -> unpacked {50,2} (rank-3 values would be lost there anyway; use Values/Times). Adapter must use Normal@ev["Values"], not ev["Values"] raw.
(3) Reconstruction: `EventSeries[values, {times}]` verified — rebuilt object's Values match the originals.
(4) Arithmetic flag: `ev + 1` and `2 ev` evaluate NATIVELY, return EventSeries, values verified correct. Elementwise scalar arithmetic is native (TemporalData threading); Dot/matrix ops not native. Flag: compute-capable for elementwise scalar ops only.
(5) Numericity: values are plain numerics (packed on Normal); numeric at module level.
(6) Recognition: `Head[ev] === EventSeries`.
(7) Adapter: shape = {ev["PathLength"], ev["ValueDimensions"]}; ingest = Normal[ev["Values"]] (packed) + ev["Times"]; egress = EventSeries[newValues, {times}] preserving the time index.

================================================================
5. SymmetrizedArray — capability-flag anchors (reconfirm only)
================================================================
- Substitution opacity CONFIRMED: `sa /. aa -> 5` returns a SymmetrizedArray whose `Normal` is STILL `{{0, aa, 0}, {-aa, 0, 0}, {0, 0, 0}}` — ReplaceAll does NOT penetrate the elements (Normal[sa /. aa->5] === Normal[sa]/.nothing; equality test with the penetrated dense form was False). Any symbolic-substitution path must go Normal -> ReplaceAll -> SymmetrizedArray.
- TensorContract structure preservation CONFIRMED: `TensorContract[sb, {{3,4}}]` on a rank-4 Symmetric[{1,2}] array returns head SymmetrizedArray and its Normal equals TensorContract of the dense form (True). Flag: tensor-algebra native, substitution opaque.

### F.3 GPUArray (WL 15.0, Metal backend, Mac ARM) - EXPLICIT compute-native, v2 registry behind guards

Live-kernel probe (this machine):

- **Construction/recognition**: `GPUArray[array]` - ONE argument only; a 2-arg call
  returns an INERT expression whose head is still GPUArray, so head-matching is not a
  valid recognition pattern: **use `GPUArrayQ`** (False for the inert forms). Exact
  (`{1/2, 1/3}`) and symbolic input error loudly (`GPUArray::spec`). First use runs a
  ~2 s one-time "Checking GPU compatibility" probe.
- **Shape without transfer**: `Dimensions` works (`{512, 512}` on a device-resident
  matrix); `ArrayQ` True; `AtomQ` True.
- **Compute-native, results stay on device**: `+`, scalar `*`, `Dot`, `Transpose`,
  `Part`, `Map` all return `GPUArray`; `Total` reduces to a scalar. Measured 512^2
  complex `Dot` x10: **3 ms GPU vs 29 ms CPU fp64** (~10x).
- **Precision is the catch (Metal has no fp64 compute)**: real storage is **Real32**
  (round-trip error ~3e-8; display shows `NumericArray[<n>, Real32]`; `Real64`
  NumericArray input is REJECTED with `GPUArray::spec`), complex storage round-trips
  exactly (fp64 storage) but **compute is fp32**: GPU `Dot` vs CPU fp64 max error
  ~3e-5 at 256^2. Elementwise `Sqrt` visibly fp32 (1.41421353...).
- **Substitution-opaque** like other atoms: `ga /. rule` does not penetrate (silent
  no-op); route via `ArrayMaterialize` -> ReplaceAll.
- **Materialization**: `Normal` returns a packed array.

**Verdict**: EXPLICIT tier, compute-native flag True - but admit via the **v2
registry behind two guards**, not v1: (1) an availability guard (hardware/driver
dependent, slow first-use probe); (2) a **precision-consent guard**: fp32 compute
silently violates fp64 expectations (amplitude norms drift at the 1e-5 scale), so
ingest of fp64 data into a GPUArray container must be an explicit opt-in, mirroring
ArrayPack's exact-value fidelity posture. Recognition `_GPUArray ? GPUArrayQ`;
materialize `Normal`; rebuild `GPUArray[data]`.

**QF opportunity**: a fast approximate statevector path - fp32 is ample for sampling,
visualization, and variational sweeps; `ArrayMaterialize` returns packed CPU arrays
for exact post-processing. Not suitable as a default amplitude container.

### CORRECTION (2026-07-26, disproved during implementation)

Section B.4 and the verdict table claimed array-valued `Piecewise` has a "unique
privilege: `Dimensions` works on the UNEVALUATED form". **This is false in WL 15.0.**
`Dimensions` does not thread through the branches; it reports the argument count, which
is `{2}` for *every* `Piecewise`, because the kernel normalizes `Piecewise[pairs]` to
`Piecewise[pairs, default]`. Verified at ranks 1, 2 and 3, e.g.

    Dimensions[Piecewise[{{ConstantArray[1., {2, 3, 4}], z < 0}}, ConstantArray[0., {2, 3, 4}]]]
    (* {2}, not {2, 3, 4} *)

`Piecewise` remains admitted to the lazy tier and its shape is still obtained without
materializing, but from the branch values rather than from `Dimensions`. The paclet pins
the wrong kernel answer alongside the right container answer, so a future version change
is caught rather than silently believed.

Two further sharpenings from the same implementation pass:

- `ifn'` is already an `InterpolatingFunction` in 15.0, so `Derivative[__][_InterpolatingFunction][args]`
  never survives as an inert form and needs no separate recognition pattern.
- Storing an unapplied `Function` really is load-bearing: a plain `ReplaceAll` on the
  bound parameter produces `Function[0.5, ...]` and a `Function::flpar` error.
