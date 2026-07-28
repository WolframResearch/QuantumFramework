> **Strategy note (user decision, kernel-verified):** complex-number support in the net
> runtime is absent and is NOT a blocker for QML - **the primary route is the
> phase-space picture**, where every tensor in the network is exactly real (verified in
> QF: Bell-state Wigner vector, CNOT and RZ(0.7) phase-space kernels all have
> max|Im| = 0; QuantumWignerTransform / Picture -> "PhaseSpace", and
> TensorNetworkCompile's existing phaseSpaceQ branch routes such circuits already).
> The section-2 real/imag doubling lift below therefore becomes an optional later
> addition for Schrodinger-picture nets, not mandatory plumbing; the QML examples in
> section 4 should be built phase-space-first.

# Net-backed lazy arrays and TensorNetwork -> NetGraph compilation (QML)

> **ARCHITECTURE CORRECTION (2026-07-26, maintainer decision + live-kernel findings).**
> Three things in the body below are superseded; read this first.
>
> 1. **The compiler already exists.** `TensorNetworkToNetGraph` is implemented and
>    exported by the TensorNetworks paclet
>    (`TensorNetworks/Kernel/ToTensorNetworkGraph.wl:447`), driven by
>    `OptimalContractionPath`, lifting numeric tensors to `NetArrayLayer` and symbolic
>    ones to `FunctionLayer` with one named input port per free symbol, and lowering
>    each pairwise contraction to Transpose/Reshape/Dot/Reshape/Transpose layers.
>    Section 3's "proposed API" for building this from scratch is obsolete; the probe
>    that produced it read Contraction.wl and Paths.wl but not ToTensorNetworkGraph.wl.
>
> 2. **NetGraph is a CONTAINER TYPE, not a contraction method.** It belongs in
>    `Wolfram/Arrays` as a container the tiers already describe: shape from
>    `NetExtract[net, "Output"]` / `Information[net, "OutputPorts"]` without evaluating,
>    materialization by `net[]`, and structural ops that preserve it by appending layers
>    (`ArrayTranspose` -> `TransposeLayer`, `ReshapeArray` -> `ReshapeLayer`, the
>    pairwise contraction -> `DotLayer`). Tensor-network contraction then does not need a
>    NetGraph-specific code path at all: once contraction routes its array operations
>    through `Wolfram/Arrays`, a network whose tensors are net-backed containers
>    contracts into a net-backed container by construction. This replaces both the
>    "emit an Inactive NetGraph" backend of section 3.4 and the later idea of adding
>    `"NetGraph"` to `$TensorNetworkContractionMethods`.
>
> 3. **No parameter hacks.** A tensor the net framework cannot lift - complex values, a
>    `Conjugate` on a parameter, any head `FunctionLayer` cannot compile - must FAIL with
>    a clear message. Do not apply `ComplexExpand`, `Refine`, `Conjugate` stripping, or an
>    automatic phase-space transformation to make it compile. Section 2's real/imag
>    doubling lift stays out of scope on the same grounds.
>
> **Verified integration facts (live kernel, WL 15.0):**
> - QF's `qc["TensorNetwork"]` returns a `TensorNetwork` object for which
>   `TensorNetworkGraphQ` is **False**; `ToTensorNetworkGraph` bridges to the Graph form
>   that `TensorNetworkToNetGraph` requires.
> - `qc["TensorNetwork"]` defaults to `"Computational" -> True`, which injects **complex**
>   basis-change tensors even into a phase-space circuit (4 of 8 tensors carried
>   imaginary parts up to 4.0). With `"Computational" -> False` a Wigner circuit's
>   tensors are all real. This option, not the picture alone, is what makes the
>   phase-space route work.
> - End to end: for `qcw = QuantumWignerTransform[QuantumCircuitOperator[{"H" -> 1,
>   "CNOT" -> {1, 2}}]]`, the net built from
>   `ToTensorNetworkGraph[qcw["TensorNetwork", "Computational" -> False]]` evaluates to a
>   `{4, 4}` array matching `N @ Normal @ TensorNetworkContract[graph, Automatic]` with
>   **max|diff| = 0**.
> - `TensorNetworkContract[graph]` with **no path argument** returns an `Inactive`
>   expression instead of contracting; the path argument is mandatory for a value.
> - A parametric Wigner circuit currently fails at `FunctionLayer::compilerr: Cannot
>   interpret Conjugate[#1] & as a network`, because QF does not know the parameter is
>   real. Per decision 3 this stays a failure for now; it is a QF parameter-semantics
>   question, not something the net layer should paper over.

Status: design proposal, backed by two verified probe reports (WL `15.0.0 for Mac OS X ARM (64-bit) (May 19, 2026)`, net runtime 15.0.2; TensorNetworks paclet loaded from checkout at `/Users/swish/src/wolfram/TensorNetworks/TensorNetworks`). Probe scripts: `/private/tmp/claude-501/-Users-swish-src-wolfram-QuantumFramework/7e8c928b-f01c-4e17-8123-ca13235727bc/scratchpad/net{1,2,3,4}.wls` and `.../scratchpad/{tn2net.wl,tncompile.wl,tncompile5.wl}`. Every kernel-level claim below was run and checked; items marked "untested" or "future work" were not.

## 1. The container class: NetGraph as a GenericArray tier (Wolfram/Arrays)

A `NetGraph` whose inputs are fully bound (source-only net) is already a lazy array container in the GenericArray sense: an atomic, inert expression that knows its shape statically and materializes on demand. Verified properties:

- `NetGraph[...]` is atomic and inert: `Head === NetGraph`, `AtomQ -> True`, structurally idempotent (reconstructing the same net gives `===`-identical expression). It does nothing until applied, so it satisfies the lazy-container contract with no `Inactive` wrapper needed.
- `Inactive[NetGraph][layers, edges]` also works, and `Activate` assembles and evaluates correctly, so a symbolic tier can emit a genuinely inactive form and `Activate` at the boundary.

**Recognition protocol.** There is no System-context predicate: `System`ValidNetQ`/`NetGraphQ` do not exist (`Names["System`*NetQ*"] == {}`). Use `MatchQ[x, _NetGraph | _NetChain]` as the cheap structural check; the authoritative validity gate is `NeuralNetworks`ValidNetQ` (True on nets, False on junk), with the caveat that it is an internal context and should be wrapped defensively. The container invariant distinguishing a fully-bound array from a parametric one: `Information[x, "InputPorts"] === <||>` for source-only (array-valued) nets, non-empty for parametric ones.

**Shape protocol.** Static shape info is available with no evaluation and works on uninitialized nets. On an uninitialized net (`NetArrayLayer[{2,3}]` feeding a Transpose):

- `Information[net, "OutputPorts"]` gives `<|"Output" -> {3,2}|>`; `"InputPorts"` gives `<||>` (source-only) or e.g. `<|"T1"->{2,3},"T2"->{3,4}|>` (input-bearing).
- `NetExtract[net, "Output"]` gives `{3,2}` directly. This is the clean `ArrayDimensions` hook.
- `Information[net, "ArraysDimensions"]` gives `<|{"arr","Array"} -> {2,3}|>` for per-leaf shapes.

Shape inference propagates through the graph at construction time, so `ArrayDimensions` on this tier is free.

**Materialize protocol.** Source-only nets evaluate as both `net[]` and `net[<||>]` (identical results); a bare `NetArrayLayer["Array"->{1.,2.}][]` also works. `NetInitialize` is NOT needed when every `NetArrayLayer` has an explicit array: the net is born initialized (`NetInitialize[net] === net` is True). A shape-only `NetArrayLayer[{2,2}]` net fails on `net[]` until `NetInitialize` fills random values; that is the QML "trainable ansatz" entry point, not the array tier. Per-leaf access: `NetExtract[net, {name, "Array"}]` returns a `NumericArray` (`Normal` to unpack). `Normal[net]` gives the layer association in internal form (`NeuralNetworks`TensorT` types etc.), useful for structure, not values.

**Precision.** Everything is fp32 by default: arrays passed as machine reals are stored `Real32` at construction (`0.1` becomes `0.10000000149011612`), and default evaluation is fp32 even for Real64-stored arrays. `NetGraph[..., WorkingPrecision->"Real64"]` and `NetArrayLayer[..., WorkingPrecision->...]` are NOT accepted. The verified fp64 path: store leaves as `NumericArray[..., "Real64"]` (retained) and evaluate with `net[..., WorkingPrecision->"Real64"]`; `NetTrain[..., WorkingPrecision->"Real64"]` works and yields `Real64` stored arrays. `TargetDevice->"CPU"` and `"GPU"` both evaluate on Apple Silicon with correct values. The GenericArray materializer for this tier should therefore thread a precision option through to call time.

**Where it sits among the tiers.** Below the fully symbolic tier (Inactive expressions) and beside other lazy-numeric containers: it is shape-static, evaluation-deferred, differentiable, and device-targetable, but real-only and fp32-by-default. Recommended dispatch order: recognize with `MatchQ[x, _NetGraph|_NetChain]`, gate with `ValidNetQ`, shape via `NetExtract[net, "Output"]`, materialize via `net[]` / `net[<|inputs|>]`.

## 2. Complex-number strategy

Complex support is NOT native and this is a hard boundary, verified three ways:

- `NetArrayLayer["Array" -> {{1.+2.I, 0.}, {0., 1.-1.I}}]` fails: `NetArrayLayer::netinvtensorvals: The value specified for Array should be a real-valued tensor.`
- Feeding a `NumericArray[..., "ComplexReal64"]` also fails; `NetEncoder[{"NumericArray","ComplexReal64"}]` fails. The net runtime is real-only, period.

The verified workaround is the standard "complex as 2-channel real" lift: 2 real leaves per complex leaf, and each pairwise contraction expands to 4 real contractions plus 2 elementwise combines (ReRe - ImIm for the real part, ReIm + ImRe for the imaginary part). Verified exact against complex `Dot` (max err < 1*^-6 at fp32):

```wl
A = {{1.+2.I, 3.-1.I},{0.+1.I, 2.+0.I}}; v = {1.-1.I, 2.+3.I};
cnet = NetGraph[<|
  "Ar"->NetArrayLayer["Array"->Re[A]], "Ai"->NetArrayLayer["Array"->Im[A]],
  "vr"->NetArrayLayer["Array"->Re[v]], "vi"->NetArrayLayer["Array"->Im[v]],
  "rr"->DotLayer[], "ii"->DotLayer[], "ri"->DotLayer[], "ir"->DotLayer[],
  "re"->ThreadingLayer[Subtract], "im"->ThreadingLayer[Plus]|>,
 {{"Ar","vr"}->"rr", {"Ai","vi"}->"ii", {"Ar","vi"}->"ri", {"Ai","vr"}->"ir",
  {"rr","ii"}->"re", {"ri","ir"}->"im", "re"->NetPort["Re"], "im"->NetPort["Im"]}];
cnet[]  (* <|"Re"->{12.,5.}, "Im"->{8.,7.}|> ; Re+I*Im == A.v exactly (max err < 1*^-6) *)
```

Design decision for the container: a complex net-backed array is a net with two output ports (`"Re"`/`"Im"`) or a stacked length-2 leading axis; the materializer recombines as `Re + I Im`. The compiler (section 3) applies the doubling transform per contraction step. Since gate unitaries are complex, this is mandatory plumbing for the general TensorNetworks backend; note however that RY/CNOT/CZ-style real-gate circuits dodge it entirely (verified end to end in section 4).

## 3. TensorNetworkContract "NetGraph" backend (Wolfram`TensorNetworks`)

### 3.1 What the paclet already provides

`TensorNetworkContract[net, path]` and the Inactive variant `TensorNetworkContraction` live in `Kernel/Contraction.wl`. A canonical path is `{{i,j}..}` (pairwise merges, `CanonicalPathQ`, `Kernel/Paths.wl`); `OptimalContractionPath`/`GreedyContractionPath` produce it. `TensorNetworkContraction` folds the tree path with `contractTensorPair[{t1 -> indices1, t2 -> indices2}]`, dispatching to one of `$TensorNetworkContractionMethods = {"ArrayDotTranspose","ArrayDot","Dot","TensorContract","TableSum"}`; a final `Transpose` reorders to `"FreeIndices"`. Per-step contracted-index metadata is available standalone via `PathIndexContractions[path, indices]` (Paths.wl:109).

Crucially, the paclet already contains the exact net-shaped lowering: `einsumDot` (Contraction.wl:109-145) emits `Transpose` (A to `[free..., contracted...]`), `Transpose` (B to `[contracted..., free...]`), `ArrayReshape` fuse (only when k > 1), `Dot`, `ArrayReshape` unfuse. Observed outputs: a matrix chain with `Method -> "Dot"` gives nested `Inactive[Dot]`; a rank-3 case gives `Inactive[Dot][Inactive[ArrayReshape][T,{2,12}], Inactive[ArrayReshape][Inactive[Transpose][B, Cycles[{{1,2}}]], {12,5}]]`. The no-path form is confirmed Inactive too: `Transpose[Inactive[TensorContract][Inactive[TensorProduct][...],{{1,4}}], Cycles[{{1,2}}]]`.

### 3.2 The pairwise-step mapping (verified)

For a step contracting index set `c` (dims K = Times @@ dc) between A (free `al`) and B (free `br`):

| Inactive op | Net layer | Note (verified) |
|---|---|---|
| `Transpose[x, perm_Cycles]` | `TransposeLayer[PermutationList[perm, rank]]` | Same convention as `Transpose`: `TransposeLayer[{2,3,1}]` on {2,3,4} gives dims {4,2,3}, elementwise-equal to `Transpose[a,{2,3,1}]`. Also checked `{3,1,2}` on 2x3x4 giving dims {3,4,2}; `TransposeLayer[]` = `{2,1}` |
| `ArrayReshape[x, shape]` | `ReshapeLayer[shape]` | Row-major, total elements fixed (2x3 -> {3,2} verified); shapes static at compile time |
| `Dot[x, y]` | `DotLayer[]` (2 inputs, edge order = port order) | Exactly WL `Dot`: contracts last dim of input 1 with first dim of input 2. Verified m.v, m.m, rank3.rank2, rank3.rank3, and v.v (rank-0 scalar 14.). Arbitrary ranks, but only that one axis pair |
| `TensorProduct[x, y]` (k = 0, outer) | `ReshapeLayer[{m,1}]` + `ReshapeLayer[{1,n}]` + `DotLayer` + `ReshapeLayer[dims_x ~Join~ dims_y]` | verified on disconnected {2}x{3} |
| numeric leaf | `NetArrayLayer["Array" -> a]` | trainable parameter by default |
| `SymbolicDeltaProductArray` leaf (hyperedge binarization) | `NetArrayLayer["Array" -> Normal[delta]]` | `Normal` is the paclet's own materialization (EinsteinSummation.wl:129) |

So one step compiles to TransposeLayer(A) -> ReshapeLayer{Prod[al], K}, TransposeLayer(B) -> ReshapeLayer{K, Prod[br]}, DotLayer -> optional ReshapeLayer unfuse; identity transposes/reshapes are elidable (`einsumDot` already elides them), and k = 1 needs no reshapes at all since `DotLayer` handles higher-rank operands like `Dot`. Under the complex lift each DotLayer becomes 4. The general Transpose+Reshape+Dot+Reshape emulation was independently verified against `TensorContract`:

```wl
(* contract axis 2 of T1[2,3,4] with axis 2 of T2[5,3,6] -> [2,4,5,6] *)
enet = NetGraph[<|
  "tA"->TransposeLayer[{1,3,2}], "rA"->ReshapeLayer[{8,3}],
  "tB"->TransposeLayer[{2,1,3}], "rB"->ReshapeLayer[{3,30}],
  "dot"->DotLayer[], "rOut"->ReshapeLayer[{2,4,5,6}]|>,
 {NetPort["T1"]->"tA"->"rA", NetPort["T2"]->"tB"->"rB", {"rA","rB"}->"dot"->"rOut"}];
(* matches TensorContract[TensorProduct[T1,T2],{{2,5}}]; max abs err 3.4*^-6 (fp32) *)
```

What is NOT available and must be handled upstream: no arbitrary-axis `TensorContract` layer, no trace (repeated index on one tensor), no multi-way contraction. `FunctionLayer[TensorContract[#,{{1,2}}]&]` fails; `FunctionLayer[Dot]` and `FunctionLayer[#1 . #2 &]` fail to construct (`FunctionLayer[#Input1 . #Input2 &]` does compile a Dot, verified {3.,7.}; `FunctionLayer[Sin]` is fine). `ThreadingLayer` has no automatic rank broadcasting (`{m, v}` fails with `ThreadingLayer::nettypeinc`); `ThreadingLayer[Times, 1]` broadcasts the lower-rank input against trailing axes numpy-style, while the forms `2->1`, `{2->1}`, `"Input2"->1` all fail. `SummationLayer[]` sums all elements to a scalar; `AggregationLayer[Total, 1]` sums level 1 and `AggregationLayer[Total, All]` is the full sum. Traces must be pre-eliminated symbolically or done with a mask-and-`AggregationLayer` trick (untested).

### 3.3 Prototype results

Emission is fully mechanical: `TensorNetworkContraction[data, path, Method -> "Dot"]` already produces exactly the layer program, and the prototype compiler (`compile` in tncompile.wl, ~25 lines) is a 6-rule pattern walk of that Inactive tree. Numerical agreement of `net[]` vs `TensorNetworkContract` (fp32, so ~1e-7):

- Matrix chain {2,3}.{3,4}.{4,2}, optimal path `{{2,3},{1,2}}`: max diff 4.5e-8.
- Rank-3 with transpose+fuse (T{2,3,4} idx {a,i,j} x B{4,3,5} idx {j,i,b} -> {2,5}, two contracted indices): 1.6e-7.
- Outer product (disconnected {2}x{3}): 2.2e-8.
- Mechanical compiler on chain optimal+greedy, rank-3, outer, and a 2-qubit circuit optimal+greedy: all agree <= 1.6e-7. Random 6-vertex/9-edge hyperedge network (after `Normal`izing delta leaves): 25 layers, output dims {3,2,2,3,3}, max diff 5.2e-7.

### 3.4 Proposed API

- New function `TensorNetworkNetGraph[tn, path, opts]` in a new `Kernel/NetGraph.wl`, implemented as a post-pass over the existing `Method -> "Dot"` Inactive output (zero disturbance to Contraction.wl), plus `TensorNetworkContraction[..., Method -> "NetGraph"]` sugar delegating to it.
- Return value: the assembled `NetGraph[...]` object directly. The assembled net is already inert and is what `NetTrain`/`NetPortGradient` consume. An `"Inactive" -> True` option gives `Inactive[NetGraph][nodes, edges]` for symbolic uniformity with the paclet's Inactive pipeline (`Activate` then assembles the net) - both forms verified to work.
- Required pre-pass: `/. d: _SymbolicDeltaProductArray | _SymbolicIdentityArray :> Normal[d]` on binarized data. Side observation from the probe: the paclet's own active path-based fold and the ActivateTensors route both currently return unevaluated for hyperedge networks containing symbolic deltas (`ArrayReshape::listrp` etc.); materializing deltas first fixes both, so this pre-pass doubles as a bug workaround.
- Shape derivation inside the compiler should use the paclet's internal `symbolicTensorDimensions` (Utilities.wl:77, handles delta arrays) to get ranks for `PermutationList[perm, rank]` without evaluating; the prototype's `Dimensions@Activate` hack is the only non-production shortcut and is trivial to replace once the code lives inside the package.

## 4. QML story (QuantumFramework)

**What works today.** Real-gate circuits compile and differentiate end to end. Verified flow: 2-qubit |00> state, RY(0.7) tensor RY(-1.3), then CNOT, all real tensors. TN reference (`TensorNetworkContract`, hyperedges psi{a,b}, RY1{c,a}, RY2{d,b}, CNOT{e,f,c,d}) = `{{0.7478, -0.5685}, {-0.2075, 0.2730}}`, matching dense matrix algebra to 0. A hand-built net (NetArrayLayer gates + DotLayer/ReshapeLayer, 9 layers) and the mechanically compiled net (10 layers, both greedy and optimal paths) agree to 2.6e-8.

**Differentiability** (verified in both probes):

- `net[<|"Input"->{1.,2.}|>, NetPortGradient[{"A","Array"}]]` returns `{{1.,2.},{1.,2.}}`, exactly d(Total[A.x])/dA; `NetPortGradient["Input"]` returns `{4.,6.}` (column sums of A). Vector outputs need no SummationLayer: the gradient of a vector-output net is the gradient of the implicit total, no message.
- Source-only nets are differentiable: `sonet[<||>, NetPortGradient[{"A","Array"}]]` returns `{{1.,1.},{1.,1.}}`. The empty association is required; `sonet[NetPortGradient[...]]` alone fails.
- On the compiled circuit: with a `#^2`/SummationLayer head appended, `qscalar[<||>, NetPortGradient[{"q","ry1","Array"}]]` returned `{{1.879, 0.}, {0.686, 0.}}` - a real gradient with respect to the gate array (zero second column because |00> only engages RY's first column). NetArrayLayer arrays are trainable parameters by default, so NetTrain on gate leaves works out of the box.
- NetTrain smoke test (train A in `A.x + B.x` to fit `M.x`, B frozen):

```wl
NetTrain[tnet, <|"Input"->xs,"Output"->ys|>, LossFunction->MeanSquaredLossLayer[],
 LearningRateMultipliers->{"B"->None}, Method->"ADAM", MaxTrainingRounds->2000,
 BatchSize->128, TrainingProgressReporting->None, RandomSeeding->1]
(* A -> M - Id to 5.8*^-6; B bit-identical to IdentityMatrix[2] *)
```

`LearningRateMultipliers -> {"B" -> None}` freezes exactly. Caveat: default settings at 300 rounds left error ~1.2; explicit ADAM with enough rounds converges cleanly.

**What needs unitarity constraints.** Gradients live in unconstrained matrix space: the nonzero norm-gradient above shows updates leave the unitary manifold. Options: (a) parameterization, theta -> gate via `FunctionLayer` - the natural QML form anyway and the recommended path; (b) projection/retraction after updates; (c) a penalty term. Parameterized leaves change the leaf emission rule (NetArrayLayer holds raw tensors; theta-parameterized gates want FunctionLayer/CompiledLayer producing the gate from scalars) and simultaneously solve unitarity for standard gates.

**Concrete example flow (target).** `QuantumCircuitOperator` -> `TensorNetwork` (existing) -> `TensorNetworkNetGraph[tn, OptimalContractionPath[tn]]` -> append loss head -> `NetTrain[net, data, LearningRateMultipliers -> freeze non-ansatz leaves]` -> read trained gates back with `NetExtract[net, {gate, "Array"}] // Normal`. Every arrow except the parameterized-leaf variant was exercised in the probes.

**Batching.** TransposeLayer/ReshapeLayer/DotLayer are shape-static, so batching over input states (a data-dependent psi port instead of a NetArrayLayer leaf) works naturally in principle - but swapping a leaf for an input port changes edge emission and is untested.

## 5. Phased plan

**v1 - container class in Wolfram/Arrays (small, ~2-4 days).** Register the NetGraph tier: recognition (`MatchQ[x, _NetGraph|_NetChain]` + `ValidNetQ` gate + `InputPorts === <||>` invariant), `ArrayDimensions` via `NetExtract[net, "Output"]`, materialize via `net[]`/`net[<||>]` with a precision option threading `WorkingPrecision -> "Real64"`, per-leaf extraction. Complex convention (Re/Im ports) defined but only the real path exercised. All primitives verified; work is glue and tests.

**v2 - TN backend in Wolfram`TensorNetworks` (medium, ~1-2 weeks).** `Kernel/NetGraph.wl` with `TensorNetworkNetGraph` as the post-pass over `Method -> "Dot"` output, `Method -> "NetGraph"` sugar, delta pre-pass, `symbolicTensorDimensions`-based shapes, `"Inactive"` option, elision of identity transposes/reshapes. The prototype compiler is ~25 lines; production cost is the complex-doubling transform (mechanical but doubles the rule set: each DotLayer to 4 + Subtract/Plus threading), trace pre-elimination, and test coverage against `TensorNetworkContract` at fp32 tolerance.

**v3 - QF QML examples (medium, ~1 week + open-ended research).** Real-gate circuit training example (works today, just needs packaging); parameterized-gate leaves via FunctionLayer (new leaf emission rule, solves unitarity for standard gates); complex-circuit example once v2's doubling lands; batched-state input port variant. The unitarity-for-raw-arrays question (retraction/penalty) is research, kept out of v3 scope.

## 6. Risks and open problems

1. **Complex amplitudes.** Real/imag doubling is mandatory for generic quantum nets (4x DotLayers per contraction, 2 leaves per tensor). Verified correct for one matvec; the general graph transform is designed but not yet implemented. Real-gate circuits work today.
2. **Parameterized leaves.** NetArrayLayer holds raw tensors; theta-parameterized gates need FunctionLayer/CompiledLayer leaves, changing emission and interacting with the complex lift. Note `FunctionLayer` construction is fragile (`FunctionLayer[Dot]`, `FunctionLayer[#1 . #2 &]`, and `TensorContract` forms all fail; slot-named forms work).
3. **Delta/hyperedge leaves.** Must be materialized with `Normal`; dense copy-tensors could be large. SparseArray-aware emission or contraction-time elision (delta absorbed into DotLayer as a gather) is future work. The paclet's active hyperedge folds also fail on symbolic deltas today (pre-existing bug the pre-pass papers over; should be fixed upstream).
4. **Precision.** fp32 by default (~1e-7 agreement in all prototypes); fp64 requires Real64 NumericArray leaves plus `WorkingPrecision -> "Real64"` at call time, and neither NetGraph nor NetArrayLayer accepts a WorkingPrecision option at construction. Quantum applications sensitive to amplitude cancellation need the fp64 path threaded everywhere.
5. **Internal-context dependencies.** `NeuralNetworks`ValidNetQ` (recognition gate) and `symbolicTensorDimensions` (compiler shapes) are internal; both need defensive wrapping or upstream promotion.
6. **Traces and hyperedges beyond binarization.** No trace-type layer exists; repeated indices on one tensor must be eliminated symbolically before emission (mask-and-AggregationLayer alternative untested).
7. **Batched/state-input variant untested.** Shape-static layers should batch, but leaf-to-port conversion changes edge emission and has not been run.
8. **Unitarity under training.** Raw-array gradients leave the unitary manifold (demonstrated concretely); acceptable for the parameterized-gate path, unsolved for raw-array ansatze.