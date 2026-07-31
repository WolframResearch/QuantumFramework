(* A QuantumState holds its amplitudes in any array container that
   Wolfram/Arrays recognizes, not only in a SparseArray.  These tests pin the
   three things that has to mean: the container that arrives is the container
   that is stored, every property answers the same regardless of which one it
   is, and the tiers that carry no elements still report their shape. *)

BeginTestSection["QuantumState containers - storage"]

(* A list of amplitudes still normalizes to a SparseArray on intake; that is
   the one coercion the constructor performs, and every internal producer
   depends on it. *)
VerificationTest[
    Head @ QuantumState[{1, 0, 0, 1} / Sqrt[2], QuantumBasis[2]]["State"],
    SparseArray,
    TestID -> "Container-list-normalizes-to-SparseArray"
]

VerificationTest[
    Head @ QuantumState[Developer`ToPackedArray[{1., 0., 0., 1.} / Sqrt[2]], QuantumBasis[2]]["State"],
    SparseArray,
    TestID -> "Container-packed-list-normalizes-to-SparseArray"
]

(* Every other container is stored exactly as it arrived. *)
VerificationTest[
    Head @ QuantumState[NumericArray[{1., 0., 0., 1.} / Sqrt[2], "ComplexReal64"], QuantumBasis[2]]["State"],
    NumericArray,
    TestID -> "Container-NumericArray-preserved"
]

VerificationTest[
    Head @ QuantumState[QuantityArray[{1., 0., 0., 1.} / Sqrt[2], "Meters"], QuantumBasis[2]]["State"],
    QuantityArray,
    TestID -> "Container-QuantityArray-preserved"
]

VerificationTest[
    Head @ QuantumState[SparseArray[{1., 0., 0., 1.} / Sqrt[2]], QuantumBasis[2]]["State"],
    SparseArray,
    TestID -> "Container-SparseArray-preserved"
]

EndTestSection[]


BeginTestSection["QuantumState containers - property agreement"]

(* The same amplitudes in three containers must answer identically across the
   property surface.  Each test compares every container in one shot, so a
   container that diverges names itself.  QuantityArray is deliberately not in
   this fixture: its elements are numbers with a unit, a different
   mathematical object, and it gets its own tests below. *)

$containerStates := Block[{amps = {1., 0., 0., 1.} / Sqrt[2]},
    QuantumState[#, QuantumBasis[2]] & /@ {
        amps,
        SparseArray[amps],
        NumericArray[amps, "ComplexReal64"]
    }
]

VerificationTest[
    Equal @@ Through[$containerStates["StateType"]],
    True,
    TestID -> "Agree-StateType"
]

VerificationTest[
    Equal @@ Through[$containerStates["Dimension"]],
    True,
    TestID -> "Agree-Dimension"
]

VerificationTest[
    Equal @@ Through[$containerStates["NumericQ"]],
    True,
    TestID -> "Agree-NumericQ"
]

VerificationTest[
    Chop[# - {0.5, 0., 0., 0.5}] & /@ Through[$containerStates["ProbabilitiesList"]],
    ConstantArray[{0, 0, 0, 0}, 3],
    TestID -> "Agree-ProbabilitiesList"
]

VerificationTest[
    With[{control = Normal @ First @ Through[$containerStates["DensityMatrix"]]},
        Chop[Normal[#] - control] & /@ Through[$containerStates["DensityMatrix"]]
    ],
    ConstantArray[ConstantArray[0, {4, 4}], 3],
    TestID -> "Agree-DensityMatrix"
]

VerificationTest[
    Equal @@ Through[$containerStates["Type"]],
    True,
    TestID -> "Agree-Type"
]

VerificationTest[
    Chop[# - 1] & /@ Through[$containerStates["Purity"]],
    {0, 0, 0},
    TestID -> "Agree-Purity"
]

VerificationTest[
    Equal @@ Through[$containerStates["Amplitude"]],
    True,
    TestID -> "Agree-Amplitude"
]

(* Operations, not only readouts: applying a gate, tensoring and tracing must
   all give the same numbers whatever the input container was. *)

VerificationTest[
    Chop[(QuantumOperator["H", {1}] @ #)["ProbabilitiesList"] - {0.25, 0.25, 0.25, 0.25}] & /@ $containerStates,
    ConstantArray[{0, 0, 0, 0}, 3],
    TestID -> "Agree-apply-Hadamard"
]

VerificationTest[
    Chop[QuantumPartialTrace[#, {1}]["ProbabilitiesList"] - {0.5, 0.5}] & /@ $containerStates,
    ConstantArray[{0, 0}, 3],
    TestID -> "Agree-partial-trace"
]

VerificationTest[
    Through[(QuantumTensorProduct[#, #] & /@ $containerStates)["Dimension"]],
    {16, 16, 16},
    TestID -> "Agree-tensor-product-dimension"
]

VerificationTest[
    Chop[#["Conjugate"]["ProbabilitiesList"] - {0.5, 0., 0., 0.5}] & /@ $containerStates,
    ConstantArray[{0, 0, 0, 0}, 3],
    TestID -> "Agree-Conjugate"
]

(* A QuantityArray computes natively, so its unit survives every operation
   that keeps the amplitudes and is normalized away by every operation that
   forms a ratio.  Both halves of that are worth pinning. *)

VerificationTest[
    Chop[QuantumState[QuantityArray[{1., 0., 0., 1.} / Sqrt[2], "Meters"], QuantumBasis[2]]["ProbabilitiesList"] - {0.5, 0., 0., 0.5}],
    {0, 0, 0, 0},
    TestID -> "QuantityArray-probabilities-are-unit-free"
]

VerificationTest[
    QuantityUnit @ QuantumState[QuantityArray[{1., 0., 0., 1.} / Sqrt[2], "Meters"], QuantumBasis[2]]["Norm"],
    "Meters",
    TestID -> "QuantityArray-Norm-carries-the-unit"
]

(* Formatting reaches the amplitudes through the same route, so a container
   that formats must not emit messages doing it. *)
VerificationTest[
    Head @ ToBoxes @ QuantumState[NumericArray[{1., 0., 0., 1.} / Sqrt[2], "ComplexReal64"], QuantumBasis[2]],
    InterpretationBox,
    TestID -> "Agree-NumericArray-formats"
]

EndTestSection[]


BeginTestSection["QuantumState containers - shape-only tiers"]

(* A lazy container carries no addressable elements, so validity has to come
   from its declared shape.  An unapplied Function of the right arity is the
   cheapest lazy container to build without solving anything. *)

VerificationTest[
    stateQ[Function[{t}, {t, 1 - t, 0, 0}]],
    True,
    TestID -> "Shape-lazy-Function-passes-rank-contract"
]

VerificationTest[
    stateQ[MatrixSymbol["M", {4, 4}]],
    True,
    TestID -> "Shape-symbolic-square-passes-rank-contract"
]

(* A non-square symbolic matrix is neither a vector of amplitudes nor a density
   matrix, so the rank contract rejects it as it always did for explicit ones. *)
VerificationTest[
    stateQ[MatrixSymbol["M", {2, 3}]],
    False,
    TestID -> "Shape-symbolic-nonsquare-rejected"
]

VerificationTest[
    stateQ[{{1, 2, 3}, {4, 5, 6}}],
    False,
    TestID -> "Shape-explicit-nonsquare-rejected"
]

(* A ragged list has no shape at all and must not be mistaken for a vector of
   two elements, which is what Dimensions alone would report. *)
VerificationTest[
    stateQ[{{1, 2}, {3}}],
    False,
    TestID -> "Shape-ragged-list-rejected"
]

VerificationTest[
    stateQ["Register"[3]],
    False,
    TestID -> "Shape-constructor-name-rejected"
]


(* Classifying a state asks whether its amplitudes are all zero and whether they
   are numbers. ArrayAllZeroQ and ArrayNumericQ answer False for a lazy
   container by contract - it has no elements to inspect - so asking them
   directly classified an all-zero net as a pure state, and "Formula" then
   rendered an empty box instead of 0. That is what
   QuantumState[NetInitialize[NetArrayLayer["Output" -> 3]]]["Formula"] hit,
   NetInitialize leaving a bare array layer at zero. The zeros are set
   explicitly here so the test does not rest on that default. *)

zeroNet = NetReplacePart[NetArrayLayer["Output" -> 3], "Array" -> {0., 0., 0.}];

liveNet = NetReplacePart[NetArrayLayer["Output" -> 3], "Array" -> {0.6, 0., 0.8}];

VerificationTest[
    QuantumState[zeroNet]["Formula"],
    RawBoxes[0],
    TestID -> "Container-lazy-zero-formula-is-zero"
]

VerificationTest[
    {QuantumState[zeroNet]["Type"], QuantumState[SparseArray[{0., 0., 0.}]]["Type"]},
    {"Degenerate", "Degenerate"},
    TestID -> "Container-lazy-zero-classifies-as-degenerate"
]

VerificationTest[
    QuantumState[zeroNet]["PureStateQ"],
    False,
    TestID -> "Container-lazy-zero-is-not-pure"
]

(* A net holds Real32 numbers, so the state is numeric whichever container
   carries it. *)
VerificationTest[
    {QuantumState[liveNet]["NumericQ"], QuantumState[liveNet]["NumberQ"]},
    {True, True},
    TestID -> "Container-lazy-numeric-state-is-numeric"
]

(* A non-zero lazy container still renders its terms rather than collapsing. *)
VerificationTest[
    MatchQ[QuantumState[liveNet]["Formula"], RawBoxes[RowBox[{_, _}]]],
    True,
    TestID -> "Container-lazy-nonzero-formula-keeps-its-terms"
]

(* The contract itself is unchanged: a symbolic container cannot be shown to be
   numeric, and must not be materialized into claiming otherwise. *)
VerificationTest[
    QuantumState[{\[FormalA], \[FormalB]}]["NumericQ"],
    False,
    TestID -> "Container-symbolic-is-not-numeric"
]

EndTestSection[]


(* === tier mixing through operations ===

   Applying an operator or a circuit contracts the state tensor against the
   gate tensors, and raw TensorProduct treats a non-explicit container as a
   SCALAR: the unevaluated product then reads downstream as one amplitude of a
   state with the right dimension and wrong contents. Non-explicit operands now
   go through the ArrayContract mixing point, and a non-explicit state folds
   through direct operator composition rather than the tensor network. *)

BeginTestSection["TierMixing"]

$symState = QuantumState[VectorSymbol["mixv", 2]]

VerificationTest[
    With[{r = QuantumOperator["H"][$symState]},
        {
            ArraySymbolicQ[r["State"]],
            ArrayDimensions[r["State"]],
            Chop[N[Normal[ArrayMaterialize[ArrayReplaceAll[r["State"], VectorSymbol["mixv", 2] -> {1., 0.}]]]]] ==
                Chop[N[Normal[QuantumOperator["H"][QuantumState[{1., 0.}]]["StateVector"]]]]
        }
    ],
    {True, {2}, True},
    TestID -> "Mixing-symbolic-state-through-an-operator"
]

VerificationTest[
    With[{r = QuantumCircuitOperator[{"H", "X"}][$symState]},
        {
            ArraySymbolicQ[r["State"]],
            Chop[N[Normal[ArrayMaterialize[ArrayReplaceAll[r["State"], VectorSymbol["mixv", 2] -> {1., 0.}]]]]] ==
                Chop[N[Normal[QuantumCircuitOperator[{"H", "X"}][QuantumState[{1., 0.}]]["StateVector"]]]]
        }
    ],
    {True, True},
    TestID -> "Mixing-symbolic-state-through-a-circuit"
]

(* A deferred contraction has values, so applying a gate may compute them; what
   matters is that the values are the contraction's, not an expression taken
   for an amplitude. *)
VerificationTest[
    With[{d = QuantumState[Inactive[Dot][SparseArray[{{0., 1.}, {1., 0.}}], SparseArray[{1., 0.}]], QuantumBasis[2]]},
        Chop[N[Normal[ArrayMaterialize[QuantumOperator["H"][d]["State"]]]]] ==
            Chop[N[Normal[QuantumOperator["H"][QuantumState[{0., 1.}]]["StateVector"]]]]
    ],
    True,
    TestID -> "Mixing-deferred-state-through-an-operator"
]

(* The mixing route is only for non-explicit CONTAINERS: a scalar state tensor -
   the identity of an empty circuit - is not one, and sending it there made
   PauliStabilizer[QuantumOperator["I"]] round-trip to a projector. *)
VerificationTest[
    Normal[PauliStabilizer[QuantumOperator["I"]]["QuantumOperator"]["Matrix"]],
    {{1, 0}, {0, 1}},
    TestID -> "Mixing-scalar-operand-stays-on-the-explicit-path"
]

EndTestSection[]
