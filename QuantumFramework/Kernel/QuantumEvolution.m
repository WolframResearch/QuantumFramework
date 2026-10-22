Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Arrays`"]

PackageExport[QuantumEvolve]
PackageExport[HamiltonianTransitionRate]
PackageExport[LindbladTransitionRates]

PackageScope[$QuantumEvolveSparseThreshold]
PackageScope[ExpandInterpolatingFunction]
PackageScope[hermiteInterpolantQ]



QuantumEvolve::error = "Differential Solver failed to find a solution"
QuantumEvolve::jump = "the jump operators `1` are not QuantumOperators"
QuantumEvolve::rates = "the rates `1` do not fit the jump operators: give as many rates as there are jump operators, `2`, or a `2` by `2` rate matrix"
QuantumEvolve::args = "the arguments `1` do not fit QuantumEvolve[hamiltonian, jumps, state, time]: jumps are {L1, ...} or {L1, ...} -> {rate1, ...}, the state is a QuantumState, Automatic or None, and the time is a symbol t or {t, t0, t1} with numeric t0 and t1"
QuantumEvolve::pulseedges = "the discontinuities of `1` in [`2`, `3`] are not all given as a list of times, so restart times are read only off those that are; if it holds a pulse shorter than the steps NDSolve takes around it, bound them with MaxStepSize"
QuantumEvolve::pulsetime = "the discontinuities of `1` were not found within `2` seconds, so no restart times are read off them; if it holds a pulse shorter than the steps NDSolve takes around it, bound them with MaxStepSize"
QuantumEvolve::pulsedata = "`1` is data sampled `2` apart, closer than the steps of up to `3` that NDSolve may take, and no restart times are read off data; set MaxStepSize to `2` or less if it holds a pulse shorter than those steps"
QuantumEvolve::pulseform = "`1` is not a formula in the time, so no restart times are read off it; if it holds a pulse shorter than the steps NDSolve takes around it, bound them with MaxStepSize"

(* the number of unknowns above which a numeric solve keeps its operators sparse: NDSolve
   integrates the dense equations faster up to about a thousand unknowns and slows sharply
   beyond *)
$QuantumEvolveSparseThreshold = 1024

Options[QuantumEvolve] = DeleteDuplicatesBy[First] @ Join[
    {"AdditionalEquations" -> {}, "ReturnEquations" -> False, "ReturnSolution" -> False, "MergeInterpolatingFunctions" -> True, "Expand" -> False},
    Options[NDSolveValue], Options[DSolveValue]
]

QuantumEvolve[
    hamiltonian_,
    lindblad_,
    observable_ ? QuantumOperatorQ,
    args___
] := With[{evolution = QuantumEvolve[hamiltonian, lindblad, None, args]},
    If[FailureQ[evolution], evolution, SuperDagger[evolution][observable["MatrixQuantumState"]]["Operator"]]
]

QuantumEvolve[hamiltonian_QuantumOperator, lindblad : {___QuantumOperator}, args___] := QuantumEvolve[hamiltonian, ToList[lindblad] -> {}, args]

(* a rule whose left side is a string or a symbol is an option, not a jump operator *)
QuantumEvolve[
    hamiltonian_QuantumOperator,
    (lindblad : Except[_List | _String | _Symbol] -> gamma_) | (lindblad : Except[_String | _Symbol] -> gamma : Except[_List]),
    args___
] := QuantumEvolve[hamiltonian, ToList[lindblad] -> ToList[gamma], args]

(* a jump operator on some of the qudits acts as the identity on the rest *)
QuantumEvolve[hamiltonian_ ? QuantumOperatorQ, lindblad : {__ ? QuantumOperatorQ} -> gammas_List, args___] /;
    AnyTrue[lindblad, #["Order"] =!= hamiltonian["Order"] &] :=
    With[{ops = padJumpOperators[hamiltonian, lindblad]},
        QuantumEvolve[First[ops], Last[ops] -> gammas, args] /; AllTrue[Last[ops], #["Order"] === First[ops]["Order"] &]
    ]

(* a matrix of rates is a Kossakowski matrix: evolve with the corresponding Liouvillian superoperator *)
QuantumEvolve[hamiltonian_ ? QuantumOperatorQ, lindblad : {__ ? QuantumOperatorQ} -> gammas_List ? MatrixQ, args___] /; jumpRatesQ[lindblad, gammas] :=
    QuantumEvolve[QuantumOperator["Hamiltonian"[hamiltonian, lindblad, gammas]], args]

QuantumEvolve[
    hamiltonian_ ? QuantumOperatorQ,
    lindblad : {___ ? QuantumOperatorQ} -> gammas_List ? VectorQ,
    defaultState : _ ? QuantumStateQ | Automatic | None : Automatic,
    defaultParameter : _Symbol | {_Symbol, _ ? NumericQ, _ ? NumericQ} | Automatic : Automatic,
    opts : OptionsPattern[]
] /; jumpRatesQ[lindblad, gammas] := Enclose @ Block[{
    state = None, basis,
    matrix, jumps,
    parameter, parameterSpec,
    numericQ, solution,
    rhs, init, events = {},
    equations, additional, return, param,
    phaseSpaceQ = False,
    mergeQ = TrueQ[OptionValue["MergeInterpolatingFunctions"]]
},
    If[ defaultState =!= None,
        state = Replace[defaultState, Automatic :> QuantumState["Register"[hamiltonian["InputDimensions"]]]];
        phaseSpaceQ = state["Picture"] === "PhaseSpace" && state["VectorQ"] && AllTrue[Sqrt[state["Dimensions"]], IntegerQ]
    ];
    ConfirmAssert[AllTrue[Prepend[lindblad, hamiltonian] , #["OutputDimensions"] == #["InputDimensions"] &]];
    ConfirmAssert[state === None || phaseSpaceQ || hamiltonian["InputDimension"] == state["Dimension"]];
    If[ defaultParameter === Automatic,
        parameter = First[Join[hamiltonian["Parameters"], Catenate[Through[lindblad["Parameters"]]]], \[FormalT]];
        parameterSpec = Replace[First[Join[hamiltonian["ParameterSpec"], Catenate[Through[lindblad["ParameterSpec"]]]], {}], {} :> {parameter, 0, 1}],

        parameter = Replace[defaultParameter, {p_Symbol, _, _} :> p];
        parameterSpec = Replace[defaultParameter, _Symbol :> {parameter, 0, 1}]
    ];
    ConfirmAssert[Developer`SymbolQ[parameter]];
    numericQ = MatchQ[defaultParameter, {_Symbol, _ ? NumericQ, _ ? NumericQ}];
    basis = If[state =!= None && phaseSpaceQ, state["Basis"], hamiltonian["OutputBasis"]];
    matrix = Progress`EvaluateWithProgress[
        If[ (phaseSpaceQ || state === None && lindblad =!= {}) && ! hamiltonian["MatrixQ"],
            HamiltonianTransitionRate[hamiltonian, If[phaseSpaceQ, basis, basis["Double"]]],
            If[ hamiltonian["MatrixQ"] || phaseSpaceQ,
                QuantumOperator[hamiltonian["Double"], If[state =!= None && phaseSpaceQ, basis, basis["Double"]]],
                QuantumOperator[hamiltonian, basis]
            ]["Matrix"] / I
        ],
        <|"Text" -> "Preparing Hamiltonian"|>
    ];
    jumps = Progress`EvaluateWithProgress[
        If[ (phaseSpaceQ || state === None) && lindblad =!= {} && ! First[lindblad]["MatrixQ"],
            PadRight[gammas, Length[lindblad], 1] * LindbladTransitionRates[lindblad, If[phaseSpaceQ, basis, basis["Double"]]],
            PadRight[Sqrt[gammas], Length[lindblad], 1] * (QuantumOperator[#, basis]["Matrix"] & /@ lindblad)
        ],
        <|"Text" -> "Preparing Lindblad jumps"|>
    ];
    (* A numeric solve reads its restart times off the time dependence as written, before the
       interpolating functions in it are merged and before TrigToExp *)
    If[ numericQ,
        {matrix, jumps} = rangeConditions[{matrix, jumps}, parameter, Rest[parameterSpec]];
        events = restartEvents[
            {matrix, jumps}, parameter, Rest[parameterSpec],
            Min[Replace[OptionValue[MaxStepSize], Except[_ ? NumericQ] -> Infinity], OptionValue[MaxStepFraction] Abs[Subtract @@ Rest[parameterSpec]]],
            OptionValue[WorkingPrecision]
        ]
    ];
    (* TrigToExp helps with some examples *)
    matrix = Progress`EvaluateWithProgress[
        TrigToExp @ Confirm @ If[mergeQ, MergeInterpolatingFunctions, Identity][matrix],
        <|"Text" -> "Preparing Hamiltonian"|>
    ];
    jumps = Progress`EvaluateWithProgress[
        TrigToExp @ Confirm @ If[mergeQ, Map[MergeInterpolatingFunctions], Identity][jumps],
        <|"Text" -> "Preparing Lindblad jumps"|>
    ];
    If[ jumps =!= {},
        If[ phaseSpaceQ || state === None,
            matrix += Total[jumps],
            state = state["MatrixState"]
        ]
    ];
    rhs = With[{s = \[FormalS][parameter]}, If[
        phaseSpaceQ || state === None || jumps === {} && (state["VectorQ"] || hamiltonian["MatrixQ"]),
        matrix . s,
        matrix . s - s . matrix + Total[With[{L = #, Ldg = ConjugateTranspose[#]}, L . s . Ldg - 1 / 2 (Ldg . L . s + s . Ldg . L)] & /@ jumps]
    ]];
    If[ state =!= None && ! phaseSpaceQ,
        state = state["Split", state["Qudits"]];
        If[ hamiltonian["MatrixQ"],
            basis = basis["Double"];
            state = QuantumState[state["Double"], basis]
            ,
            state = QuantumState[state, hamiltonian["Input"]["Dual"]]
        ]
    ];
    init = If[defaultState === None, IdentityMatrix[Length[matrix], SparseArray], state["State"]];
    additional = Join[Flatten[{OptionValue["AdditionalEquations"]}], events];
    equations = Join[
        {
            \[FormalS]'[parameter] == rhs,
            \[FormalS][If[numericQ, parameterSpec[[2]], 0]] == init
        },
        additional
    ];
    return = If[ numericQ,
        \[FormalS],
        Element[
            \[FormalS][parameter],
            If[ state === None,
                Matrices[Dimensions[matrix]],
                If[state["VectorQ"], Vectors[state["Dimension"], If[phaseSpaceQ, Reals, Complexes]], Matrices[state["MatrixDimensions"]]]
            ]
        ]
    ];
    param = If[numericQ, parameterSpec, parameter];
    If[TrueQ[OptionValue["ReturnEquations"]], Return[{equations, return, param}]];

    (* hide sparse arrays over delayed definitions and inject numeric paratemer, initial condition is always dense *)
    If[ numericQ && ! mergeQ,
        Module[{frhs},
            Block[{s, t},
                With[{
                    def = rhs /. {
                        \[FormalS][_] -> s,
                        f_InterpolatingFunction[_] :> f[t],
                        parameter -> t
                    }
                },
                    SetDelayed @@ Hold[frhs[s_ ? ArrayQ, t_], Unevaluated[def] /. sa_SparseArray ? SparseArrayQ :> MapSparseArray[ReplaceAll[parameter -> t], sa]]
                ]
            ];
            equations = Join[
                {
                    \[FormalS]'[parameter] == frhs[\[FormalS][parameter], parameter],
                    \[FormalS][parameterSpec[[2]]] == Normal[init]
                },
                additional
            ]
        ],
        equations = If[ numericQ && Times @@ ArrayDimensions[init] > $QuantumEvolveSparseThreshold,
            Join[
                {
                    \[FormalS]'[parameter] == sparseNumericRHS[rhs, parameter, OptionValue[WorkingPrecision]],
                    \[FormalS][parameterSpec[[2]]] == Normal[init]
                },
                additional /. sa_SparseArray ? SparseArrayQ :> Normal[sa]
            ],
            equations /. sa_SparseArray ? SparseArrayQ :> Normal[sa]
        ]
    ];

    solution = If[numericQ,
        Module[{time, protectedQ = MemberQ[Attributes[Evaluate[parameter]], Protected]},
            WithCleanup[
                Unprotect[\[FormalS]];
                If[protectedQ, Unprotect[Evaluate[parameter]]]
                ,
                Progress`EvaluateWithProgress[
                    NDSolveValue[
                        equations,
                        return,
                        param,
                        FilterRules[{opts}, Options[NDSolveValue]],
                        StepMonitor :> (time = parameter)
                    ],
                    <|"Text" -> "Integrating", "Progress" :> time / Last[param], "Percentage" :> time / Last[param], "ElapsedTime" -> Automatic, "RemainingTime" -> Automatic|>
                ]
                ,
                Protect[\[FormalS]];
                If[protectedQ, Protect[Evaluate[parameter]]]
            ]
        ],
        If[ phaseSpaceQ && FreeQ[matrix, parameter],
            (* Phase-space evolution is the linear system s'[t] == R.s[t] with a
               time-independent rate matrix R, whose closed form is the matrix
               exponential.  DSolveValue can fail to return it under the real
               vector domain, so form it directly. *)
            MatrixExp[parameter Normal[matrix]] . Normal[init],
            DSolveValue[
                equations,
                return,
                param,
                FilterRules[{opts}, Options[DSolveValue]]
            ]
        ]
    ];
    If[ TrueQ[OptionValue["ReturnSolution"]], Return[solution]];
    (* The solver returns one array-valued InterpolatingFunction.  Applied to
       the parameter it is already a lazy array container, so it is stored
       whole: expanding it into one scalar interpolation per amplitude
       duplicates the time grid d^2 times and turns each later read into d^2
       separate interpolation calls.  "Expand" -> True restores that form. *)
    If[ MatchQ[solution, _InterpolatingFunction],
        solution = If[
            TrueQ[OptionValue["Expand"]],
            ExpandInterpolatingFunction[solution, parameter],
            solution[parameter]
        ]
    ];
    Which[
        state === None && MatchQ[ArrayDimensions[solution], {d_, d_}],
        QuantumOperator[
            rearrangedState[
                solution,
                If[ ArrayDimensions[solution] == hamiltonian["MatrixNameDimensions"] ^ 2,
                    QuantumState[ArrayVector[#], hamiltonian["Basis"]["Double"]]["Undouble"] &,
                    QuantumState[ArrayVector[#], hamiltonian["Basis"]] &
                ]
            ],
            hamiltonian["Order"],
            "ParameterSpec" -> DeleteDuplicatesBy[Append[hamiltonian["ParameterSpec"], parameterSpec], First]
        ],
        stateQ[solution],
        If[ state =!= None && phaseSpaceQ,
            QuantumState[
                QuantumState[
                    QuantumWeylTransform[QuantumState[solution, basis, "Picture" -> "PhaseSpace"]],
                    Sqrt[basis["Dimensions"]]
                ]["Double"],
                basis,
                "ParameterSpec" -> DeleteDuplicatesBy[Append[basis["ParameterSpec"], parameterSpec], First]
            ],
            If[hamiltonian["MatrixQ"], #["Undouble"], #] & @ QuantumState[solution, basis, "ParameterSpec" -> DeleteDuplicatesBy[Append[basis["ParameterSpec"], parameterSpec], First]]
        ],
        True,
        Message[QuantumEvolve::error];
        solution
    ]
]

(* a call without jump operators gets an empty list of them; the rule it adds in their
   place keeps this definition from applying a second time *)
QuantumEvolve[hamiltonian_QuantumOperator, args___] /; ! MatchQ[{args}, {_List -> _, ___}] := QuantumEvolve[hamiltonian, {} -> {}, args]

(* a call that reaches this definition matched none of the forms above: name the argument at fault *)
QuantumEvolve[_QuantumOperator, lindblad_List -> gammas_, args___] := Which[
    ! AllTrue[lindblad, QuantumOperatorQ],
    With[{bad = Select[lindblad, ! QuantumOperatorQ[#] &]},
        Message[QuantumEvolve::jump, bad];
        Failure["InvalidJumpOperators", <|"MessageTemplate" :> QuantumEvolve::jump, "MessageParameters" -> {bad}|>]
    ],
    ! jumpRatesQ[lindblad, gammas],
    Message[QuantumEvolve::rates, gammas, Length[lindblad]];
    Failure["InvalidRates", <|"MessageTemplate" :> QuantumEvolve::rates, "MessageParameters" -> {gammas, Length[lindblad]}|>],
    True,
    Message[QuantumEvolve::args, {args}];
    Failure["InvalidArguments", <|"MessageTemplate" :> QuantumEvolve::args, "MessageParameters" -> {{args}}|>]
]


MapSparseArray[f_, sa_] := SparseArray[Thread[sa["ExplicitPositions"] -> f /@ sa["ExplicitValues"]], Dimensions[sa]]

(* A sparse array whose entries depend on the parameter, as {coefficient, array} terms whose
   arrays do not: each entry expands into terms, and each term splits into the factors that
   hold the parameter and the rest. Plus leaves an exact and an inexact multiple of the same
   factor apart, so one entry can contribute several terms to an array; they are summed. *)
parameterTerms[sa_SparseArray, parameter_] := KeyValueMap[
    {#1, SparseArray[Normal @ Merge[#2, Total], Dimensions[sa]]} &,
    GroupBy[
        Catenate @ MapThread[
            {pos, value} |-> Map[
                term |-> With[{factors = If[MatchQ[term, _Times], List @@ term, {term}]},
                    Times @@ Select[factors, ! FreeQ[#, parameter] &] -> (pos -> Times @@ Select[factors, FreeQ[#, parameter] &])
                ],
                If[MatchQ[#, _Plus], List @@ #, {#}] & @ Expand[value]
            ],
            {sa["ExplicitPositions"], sa["ExplicitValues"]}
        ],
        First -> Last
    ]
]

(* The values the arrays hold, the explicit ones of a SparseArray, which ReplaceAll, FreeQ and
   Cases do not see into *)
arrayValues[arrays_] := Flatten[{arrays /. sa_SparseArray ? SparseArrayQ :> sa["ExplicitValues"]}]

(* Each condition on the time is refined on the open range of a numeric solve: one that is
   constant there becomes True or False, and one bounded at an end of the range loses that bound.
   NDSolve's discontinuity processing stalls on a Boole condition that vanishes where the
   integration stops (NDSolve::ndsz), and the values at the ends do not enter the solution. *)
rangeConditions[arrays_, t_, {t0_, t1_}] := With[{
    rules = DeleteCases[
        # -> Refine[#, Min[t0, t1] < t < Max[t0, t1]] & /@ DeleteDuplicates @ Cases[arrayValues[arrays],
            condition : (_Less | _LessEqual | _Greater | _GreaterEqual | _Inequality) /; ! FreeQ[condition, t],
            {0, Infinity}
        ],
        condition_ -> condition_
    ]
},
    If[ rules === {},
        arrays,
        ReplaceAll[arrays /. sa_SparseArray ? SparseArrayQ :> MapSparseArray[ReplaceAll[rules], sa], rules]
    ]
]

(* The factors of the time dependence: each term of each value, without its constant factors *)
timeFactors[arrays_, t_] := DeleteDuplicates @ Cases[
    Catenate[If[MatchQ[#, _Plus], List @@ #, {#}] & /@ Expand[Select[arrayValues[arrays], ! FreeQ[#, t] &]]],
    term_ /; ! FreeQ[term, t] :> Times @@ Select[If[MatchQ[term, _Times], List @@ term, {term}], ! FreeQ[#, t] &]
]

(* The parts of a factor that are no formula in the time: data, as an InterpolatingFunction, and
   functions that are not System` symbols, which FunctionDiscontinuities knows nothing about *)
opaqueParts[g_, t_] := DeleteDuplicates @ Cases[g,
    f : h_[__] /; ! MatchQ[h, _Symbol ? (Context[#] === "System`" &)] && ! FreeQ[f, t],
    {0, Infinity}
]

(* No restart times are read off data. It is named when its samples lie closer together than the
   largest step NDSolve may take, beyond rounding, since a pulse they resolve can then fit inside
   one step. *)
opaqueMessage[f : interpolation_InterpolatingFunction[_], {lower_, upper_}, maxStep_] :=
    With[{spacing = Min[Differences[Select[First[interpolation["Coordinates"]], Between[{lower, upper}]]]]},
        If[spacing (1 + 10 ^ -9) < maxStep, Message[QuantumEvolve::pulsedata, f, spacing, maxStep]]
    ]

opaqueMessage[f_, _, _] := Message[QuantumEvolve::pulseform, f]

(* A Piecewise as the sum of its values, each times the Boole of its own condition and of no
   earlier one: Piecewise holds its arguments, so ComplexExpand does not reach into it, while
   Boole is real *)
booleTerms[g_] := g //. HoldPattern[Piecewise[cases_, default_ : 0]] :> With[{conditions = cases[[All, 2]]},
    Total[MapThread[#1 Boole[#2 && ! Or @@ #3] &, {cases[[All, 1]], conditions, Most[FoldList[Append, {}, conditions]]}]] +
        default Boole[! Or @@ conditions]
]

(* The discontinuities of a factor of the time dependence in [lower, upper], exact. Over the reals
   FunctionDiscontinuities counts the times where a function is not real, so it reads the real and
   imaginary parts. The factor is rationalized first, so that no approximate zero is left in an
   imaginary part and Reduce solves an exact system. Reduce describes the whole set: a part of it
   that is not isolated times, or an analysis that does not finish, gives no restart times and a
   message. *)
edgeTimes[g_, t_, {lower_, upper_}] := With[{seconds = 5},
    With[{set = TimeConstrained[
        Reduce[
            FunctionDiscontinuities[ComplexExpand[{Re[#], Im[#]}] & @ booleTerms[Rationalize[g, 0]], {t}, Reals] &&
                Rationalize[lower, 0] <= t <= Rationalize[upper, 0],
            t,
            Reals
        ],
        seconds,
        Message[QuantumEvolve::pulsetime, g, seconds]; False
    ]},
        With[{parts = Replace[set, {False -> {}, HoldPattern[Or[p__]] :> {p}, p_ :> {p}}]},
            If[! MatchQ[parts, {(t == _) ...}], Message[QuantumEvolve::pulseedges, g, lower, upper]];
            Cases[parts, t == time_ :> time]
        ]
    ]
]

(* Times across an exponential profile, which TrigToExp makes of Gaussian, sech and tanh pulses,
   read off the real parts of the coefficients of its exponent, so that a phase or a detuning
   does not move them. For a Gaussian, its centre and two and four standard deviations either
   side, so that no step within four standard deviations of the centre is longer than two; for an
   exponential, the times where its magnitude is E^0, E^+-2 and E^+-4. *)
profileTimes[g_, t_] := Catenate @ Cases[TrigToExp[g], Power[E, q_] /; PolynomialQ[q, t] :> exponentTimes[Re @ CoefficientList[q, t]], {0, Infinity}]

exponentTimes[{_, b_, a_}] /; TrueQ[a < 0] := - b / (2 a) + {-4, -2, 0, 2, 4} / Sqrt[- 2 a]

exponentTimes[{b_, a_}] /; TrueQ[a != 0] := - b / a + {-4, -2, 0, 2, 4} / Abs[a]

exponentTimes[_] := {}

(* Restart events of a numeric solve, so that step control cannot step over a feature of the time
   dependence much shorter than the steps around it. Each discontinuity inside the range is crossed
   with "CrossDiscontinuity" at its exact time: NDSolve's own processing, which finds an edge by the
   sign change of a switching function across a step, misses one whose switching function changes
   sign twice within a step, and when a window is short against its distance from the time origin
   it carries the solution past an edge on the wrong side of it. A restart halfway between each
   pair of consecutive discontinuities keeps a grid point inside every window. The profile times
   of exponentials are restarts too, except on a discontinuity, which is crossed already. *)
restartEvents[arrays_, t_, {t0_, t1_}, maxStep_, precision_] := With[{
    lower = Min[t0, t1],
    upper = Max[t0, t1],
    tolerance = 10 ^ -9 Abs[t1 - t0],
    factors = timeFactors[arrays, t]
},
    With[{opaque = DeleteDuplicates @ Catenate[opaqueParts[#, t] & /@ factors]},
        Scan[opaqueMessage[#, {lower, upper}, maxStep] &, opaque];
        With[{edges = Union @ N[Catenate[edgeTimes[# /. Thread[opaque -> 1], t, {lower, upper}] & /@ factors], precision]},
            Join[
                WhenEvent[t == #, "CrossDiscontinuity"] & /@ Select[edges, lower < # < upper &],
                WhenEvent[t == #, "RestartIntegration"] & /@ First /@ Split[
                    Sort @ Join[
                        Mean /@ Partition[edges, 2, 1],
                        Select[N[Catenate[profileTimes[#, t] & /@ factors], precision], lower < # < upper && Min[Abs[# - edges]] > tolerance &]
                    ],
                    Abs[#2 - #1] < tolerance &
                ]
            ]
        ]
    ]
]

(* NDSolve multiplies a numeric SparseArray without densifying it but rejects one holding a
   symbol, so an array that depends on the parameter becomes a sum of numeric arrays with
   parameter-dependent coefficients *)
sparseNumericRHS[rhs_, parameter_, precision_] := ReplaceAll[
    ReplaceRepeated[rhs,
        Dot[a___, sa_SparseArray ? SparseArrayQ /; ! FreeQ[sa["ExplicitValues"], parameter], b___] :>
            Total[#1 Dot[a, #2, b] & @@@ parameterTerms[sa, parameter]]
    ],
    sa_SparseArray ? SparseArrayQ :> N[sa, precision]
]

MergeInterpolatingFunctions[array_ ? SparseArrayQ] := Enclose @ Block[{pos, values, ifs, grid, dims, dim, params, param},
    pos = array["ExplicitPositions"];
    values = array["ExplicitValues"];
	ifs = Cases[values, _InterpolatingFunction[_], {ArrayDepth[array], Infinity}, Heads -> True];
	If[Length[ifs] == 0, Return[array]];
    params = First /@ ifs;
    ConfirmAssert[SameQ @@ params];
    ifs = Head /@ ifs;
	dims = Through[ifs["OutputDimensions"]];
	ConfirmAssert[SameQ @@ dims];
    param = First[params];
	dim = First[dims];
	ConfirmAssert[dim === {}];
	grid = Catenate[Intersection @@ Through[ifs["Grid"]]];
    (* make custom linear interpolation supporting sparse arrays *)
	(* Interpolation[{#, Normal @ SparseArray[Thread[pos -> (values /. param -> #)], Dimensions[array]]} & /@ grid, InterpolationOrder -> 1][param] *)
	With[{arrays = SparseArray[Thread[pos -> (values /. param -> #)], Dimensions[array]] & /@ grid},
        Module[{f, getArray},
            Block[{t},
                With[{rhs = Piecewise @ MapIndexed[
                    With[{a = (#[[2]] - t) / (#[[2]] - #[[1]]), b = (t - #[[1]]) / (#[[2]] - #[[1]])},
                        {a getArray[#2[[1]]] + b getArray[#2[[1]] + 1], Between[t, #]}
                    ] &,
                    Partition[grid, 2, 1]
                ]},
                    SetDelayed @@ {f[t_ ? NumericQ], rhs}
                ];
                Format[f] ^:= Interpretation["SparseInterpolationFunction", f];
                getArray[i_] := arrays[[i]];
                f[param]
            ]
        ]
    ]
]

MergeInterpolatingFunctions[array_ ? ArrayQ] := MergeInterpolatingFunctions[SparseArray[array]]


(* For a first-order system, as QuantumEvolve's are, NDSolve's InterpolatingFunction at its
   default interpolation order is a cubic Hermite interpolant: between two successive grid times
   it is the cubic through the values and first derivatives at both. Each component is rebuilt
   from those values and derivatives, so it is that component of the solution. With
   InterpolationOrder -> All the rebuild is still the cubic through the same values and
   derivatives, and departs from the solution by that cubic's interpolation error. The cubic on a
   step reads only the step's two ends: the pair of grid times, apart by about 1e-13 of their
   magnitude or less, by which NDSolve brackets an edge it crosses stays a step of its own, as it
   is in the solution. A cubic through values alone reads two more grid times, and through two
   that close it swings far beyond them.

   Interpolation refuses a repeated time, which NDSolve gives where an event resets the state,
   with a different value on each side: the grid is split there, each run is interpolated on its
   own, and the runs are joined, the earlier one up to and including its last time, as the
   solution takes the value before the reset at the time itself. A repeat carrying the same value
   and derivative is no reset and is dropped. The join is built with Apply, so the bound variable
   is not renamed. *)
gridInterpolation[rows_List] := With[{
    runs = If[DuplicateFreeQ[rows[[All, 1]]], {rows}, Split[DeleteDuplicates[rows], First[#1] =!= First[#2] &]]
},
    With[{interpolants = If[ Length[#] == 1,
        With[{value = #[[1, 2]]}, value &],
        Interpolation[#, InterpolationOrder -> 3]
    ] & /@ runs},
        If[ Length[interpolants] == 1,
            First[interpolants],
            Function @@ {\[FormalT], Piecewise[
                MapThread[{#1[\[FormalT]], \[FormalT] <= #2} &, {Most[interpolants], Most[runs][[All, -1, 1, 1]]}],
                Last[interpolants][\[FormalT]]
            ]}
        ]
    ]
]

(* The values and the first derivatives on the grid, each with the grid index moved last, so
   that the components of the two pair up level by level. The rebuild reproduces what NDSolve
   returns for a first-order system, and not an interpolant built from values alone, so
   computableState asks hermiteInterpolantQ before it calls this. *)
ExpandInterpolatingFunction[f_InterpolatingFunction, parameter_] := With[{
    times = List /@ First[f["Coordinates"]],
    rank = ArrayDepth[f["ValuesOnGrid"]] - 1
},
    MapThread[
        gridInterpolation[Thread[{times, #1, #2}]][parameter] &,
        Transpose[#, InversePermutation[Cycles[{Range[rank + 1]}]]] & /@ {f["ValuesOnGrid"], f'["ValuesOnGrid"]},
        rank
    ]
]

(* On a step from a to b of width w, at s = (t - a) / w, the cubic Hermite interpolant of the
   values and first derivatives on the grid takes
   (2 s^3 - 3 s^2 + 1) f[a] + (s^3 - 2 s^2 + s) w f'[a] + (3 s^2 - 2 s^3) f[b] + (s^3 - s^2) w f'[b],
   and for a first-order system, as QuantumEvolve's are, NDSolve's InterpolatingFunction at its
   default interpolation order is that interpolant. It is checked a third and two thirds of the way
   across every step: a polynomial of degree five or less that keeps the step's end values and
   slopes departs from the cubic by s^2 (1 - s)^2 (c0 + c1 s), which vanishes at both points only
   if it vanishes everywhere. Each s is read off its point as rounded, since rounding moves the
   point off the exact third by more than the tolerance allows, and each value is compared against
   the size of the four terms that make it, so that the check holds to rounding at any scale of the
   values. One that Interpolation builds from values alone reports the same "InterpolationMethod"
   but passes a cubic through four neighbouring values instead, and generally fails the check.
   Only a function of one time, with no units on the time or the values, is checked; a step of
   zero length, at a repeated time, has no interior, and values that are not numbers fail. *)
hermiteInterpolantQ[f_InterpolatingFunction] := Length[f["Domain"]] == 1 && f["QuantityUnits"] === {None, None} && f["InterpolationMethod"] === "Hermite" && With[{
    times = First[f["Coordinates"]],
    values = f["ValuesOnGrid"],
    slopes = f'["ValuesOnGrid"]
},
    With[{steps = Flatten @ Position[Sign[Differences[times]], 1]},
        With[{a = times[[steps]], w = times[[steps + 1]] - times[[steps]]},
            AllTrue[{1 / 3, 2 / 3}, Function[fraction, With[{points = a + fraction w}, With[{s = (points - a) / w},
                With[{terms = {
                    (2 s ^ 3 - 3 s ^ 2 + 1) values[[steps]], (s ^ 3 - 2 s ^ 2 + s) w slopes[[steps]],
                    (3 s ^ 2 - 2 s ^ 3) values[[steps + 1]], (s ^ 3 - s ^ 2) w slopes[[steps + 1]]
                }},
                    TrueQ[Max[Abs[f[points] - Total[terms]] - 10 ^ -12 Total[Abs[terms]]] <= 0]
                ]
            ]]]]
        ]
    ]
]

ExpandInterpolatingFunction[f_InterpolatingFunction[parameter_]] := ExpandInterpolatingFunction[f, parameter]

ExpandInterpolatingFunction[array_ ? ArrayQ] := ExpandInterpolatingFunction[MergeInterpolatingFunctions[array]]

ExpandInterpolatingFunction[qs_QuantumState] := QuantumState[ExpandInterpolatingFunction @ qs["State"], qs["Basis"]]


(* The state that build makes of the propagator NDSolve returns, where build only moves the entries
   of the array it is given, as a reshape or a transpose does. Where it moves each entry is read off
   the array of the entries' own positions, and the values and first derivatives on the grid are
   moved alike: their cubic Hermite interpolant is the propagator so rearranged, as NDSolve's own
   interpolant is the cubic Hermite interpolant of its values and derivatives. Given the interpolant
   itself, a reshape is rebuilt by Arrays from the values alone, which swings far beyond the
   solution where NDSolve, crossing an edge, puts grid times at the edge and within 2e-11 of
   it on either side.

   Interpolation refuses a repeated time, which an event that resets the state gives: each entry is
   then expanded on its own, run by run, before it is moved. An interpolant that is no cubic Hermite
   interpolant of its values and derivatives is given to build as it is. hermiteInterpolantQ tells
   one apart by its method before it reads any derivative value, and the values of its derivative
   are read here only after that: reading those of the "Local Taylor series" interpolant that
   InterpolationOrder -> All gives NDSolve can end the kernel or abort with General::nomem. *)
rearrangedState[f_InterpolatingFunction[parameter_], build_] /; hermiteInterpolantQ[f] := With[{
    times = First[f["Coordinates"]]
},
    If[ DuplicateFreeQ[times],
        With[{positions = build[ArrayReshape[Range[Times @@ f["OutputDimensions"]], f["OutputDimensions"]]]},
            With[{order = Flatten[Normal[positions["State"]]], shape = Prepend[Dimensions[positions["State"]], Length[times]]},
                QuantumState[
                    Interpolation[
                        Thread[Prepend[
                            ArrayReshape[ArrayReshape[#, {Length[times], Length[order]}][[All, order]], shape] & /@ {f["ValuesOnGrid"], f'["ValuesOnGrid"]},
                            List /@ times
                        ]],
                        InterpolationOrder -> 3
                    ][parameter],
                    positions["Basis"]
                ]
            ]
        ],
        build[ExpandInterpolatingFunction[f, parameter]]
    ]
]

rearrangedState[solution_, build_] := build[solution]


HamiltonianTransitionRate[h_QuantumOperator, args___] := Chop[QuantumPhaseSpaceTransform[QuantumOperator["Hamiltonian"[h]], args]["Matrix"] / I]

LindbladTransitionRates[ls : {___QuantumOperator}, args___] := Chop[QuantumPhaseSpaceTransform[QuantumOperator["Hamiltonian"[None, {#}]], args]["Matrix"] / I] & /@ ls

