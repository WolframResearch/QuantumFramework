Package["Wolfram`QuantumFramework`"]

PackageScope["$QuantumMeasurementOperatorNames"]



$QuantumMeasurementOperatorNames = {"M", "RandomHermitian", "WignerMICPOVM", "GellMannMICPOVM", "TetrahedronSICPOVM", "QBismSICPOVM", "HesseSICPOVM", "HoggarSICPOVM", "RandomPOVM", "Lueders"}


QuantumMeasurementOperator["M"[args__ : 2], opts___] := QuantumMeasurementOperator[QuantumBasis[args], opts]

QuantumMeasurementOperator["RandomHermitian"[args___], target : _ ? targetQ : {1}, opts___] := With[{
    basis = QuantumBasis[args, "Label" -> "Random"]
},
    QuantumMeasurementOperator[
        QuantumOperator[
            With[{m = RandomComplex[1 + I, {basis["Dimension"], basis["Dimension"]}]}, (m + ConjugateTranspose[m]) / 2],
            target,
            basis
        ],
        opts
]
]

QuantumMeasurementOperator["WignerMICPOVM"[args___], target : _ ? targetQ : {1}, opts___] := Enclose @ Simplify @ QuantumMeasurementOperator[
    QuantumMeasurementOperator[
        QuantumMeasurementOperator[ConfirmBy[QuantumWignerMICPOVM[args], ArrayQ[#, 3] &], target],
        With[{basis = QuditBasis["WignerMIC"[args]]},
            QuantumBasis[QuantumTensorProduct[QuditBasis[basis["Names"]], QuditBasis[Sqrt[basis["Dimension"]]]], QuditBasis[Sqrt[basis["Dimension"]]], "Label" -> "WignerMIC"]
        ]
    ],
    opts
]

QuantumMeasurementOperator["GellMannMICPOVM"[d : _Integer ? Positive : 2, s_ : 0], opts___] :=
    QuantumMeasurementOperator[
        QuantumMeasurementOperator[
            GellMannMICPOVM[d, s],
            QuantumBasis[QuantumTensorProduct[QuditBasis[Subscript["\[ScriptCapitalG]", #] & /@ Range[d ^ 2]], QuditBasis[d]], QuditBasis[d], "Label" -> "GellMannMIC"]
        ],
        opts
    ]

QuantumMeasurementOperator["RandomPOVM"[d : _Integer ? Positive : 2, methodOpts : OptionsPattern[]], opts___] :=
    Enclose @ QuantumMeasurementOperator[
        With[{
            povm = ConfirmBy[Replace[OptionValue[{methodOpts, Method -> "Haar"}, Method], {name_, args___} | name_ :>
                Switch[name, "Haar", RandomHaarPOVM, "Bloch", RandomBlochMICPOVM][d, args]], TensorQ]
        },
            QuantumMeasurementOperator[
                povm,
                QuantumBasis[QuantumTensorProduct[QuditBasis[Subscript["\[ScriptCapitalR]", #] & /@ Range[Length[povm]]], QuditBasis[d]], QuditBasis[d], "Label" -> "RandomMIC"]
            ]
        ],
        opts
    ]

QuantumMeasurementOperator["TetrahedronSICPOVM"[HoldPattern[angles : PatternSequence[__] : Sequence[0, 0, 0]]], opts___] :=
    Simplify @ QuantumMeasurementOperator[
        QuantumMeasurementOperator[
            KroneckerProduct[#, Conjugate[#]] / 2 & [QuantumOperator["U"[angles]]["Matrix"] . #] & /@ {
                {1, 0},
                {1, Sqrt[2] E ^ (I 4 Pi / 3)} / Sqrt[3],
                {1, Sqrt[2] E ^ (I 2 Pi / 3)} / Sqrt[3],
                {1, Sqrt[2]} / Sqrt[3]
            },
            QuantumBasis[QuantumTensorProduct[QuditBasis[Subscript["\[ScriptCapitalT]", #] & /@ Range[4]], QuditBasis[2]], QuditBasis[2], "Label" -> "TetrahedronSIC"]
        ],
        opts
    ]

QuantumMeasurementOperator["QBismSICPOVM"[d : _Integer : 2], opts___] := Enclose @ QuantumMeasurementOperator[
    QuantumMeasurementOperator[
        Confirm @ QBismSICPOVM[d],
        QuantumBasis[QuantumTensorProduct[QuditBasis[Subscript["\[ScriptCapitalQ]", #] & /@ Range[d ^ 2]], QuditBasis[d]], QuditBasis[d], "Label" -> "QBismSIC"]
    ],
    opts
]

QuantumMeasurementOperator["HesseSICPOVM"[], opts___] := Enclose @ QuantumMeasurementOperator[
    QuantumMeasurementOperator[
        Confirm @ HesseSICPOVM[],
        QuantumBasis[QuantumTensorProduct[QuditBasis[Subscript["\[ScriptCapitalH]", #] & /@ Range[9]], QuditBasis[3]], QuditBasis[3], "Label" -> "HesseSIC"]
    ],
    opts
]

QuantumMeasurementOperator["HoggarSICPOVM"[], opts___] := Enclose @ QuantumMeasurementOperator[
    QuantumMeasurementOperator[
        Confirm @ HoggarSICPOVM[],
        QuantumBasis[QuantumTensorProduct[QuditBasis[Subscript["\[ScriptCapitalH]", #] & /@ Range[64]], QuditBasis[8]], QuditBasis[8], "Label" -> "HoggarSIC"]
    ],
    opts
]


QuantumMeasurementOperator::luedersnotsquare = "\"Lueders\" measures an observable: an operator whose output qudits are its input qudits, with the same dimensions."
QuantumMeasurementOperator::luedersnotnormal = "\"Lueders\" measures the eigenspaces of a normal operator; this operator is not normal, or is symbolic and could not be shown to be normal."
QuantumMeasurementOperator::luedersnotcomplete = "The eigenspace projectors of this operator could not be shown to sum to the identity, since its symbolic eigenvectors are not orthonormalized; give the symbols values first."
QuantumMeasurementOperator::luedersoneoutcome = "This operator has one eigenvalue on the targets, so \"Lueders\" has one outcome, and a measurement operator holds two or more."
QuantumMeasurementOperator::luederstarget = "The target `1` of \"Lueders\" should list distinct positive qudits: some of the operator's input order `2`, or as many qudits as that order to place the operator on."

(* "Lueders"[op] measures the observable op on the targets QuantumMeasurementOperator[op, args]
   gives it, with one outcome per distinct eigenvalue, whose Kraus operator is the projector
   onto that eigenvalue's eigenspace (Lueders' rule). The eigenbasis measurement
   "SuperOperator" builds is Sum_i |i> (x) K_i, one outcome i per orthonormal eigenvector v_i,
   with K_i = |v_i><v_i| on the targets and the identity on the other qudits, and the
   eigenqudit at order 0. The K_i of one eigenvalue add up to the projector onto its
   eigenspace, and a 0-1 matrix acting on the eigenqudit alone adds them for every eigenvalue
   at once: Sum_i |lambda(i)> (x) K_i = Sum_lambda |lambda> (x) P_lambda. Adding the maps
   K_i.rho.K_i^† of one eigenvalue instead would give von Neumann's rule, which dephases rho
   inside the eigenspace. *)

(* An integer target n is the target {n}. A target that does not lie inside op's input order
   names no qudits of op, so op is placed on it, as in the circuit form
   "Lueders"[op] -> order; that takes as many qudits as op's input order. *)
QuantumMeasurementOperator["Lueders"[op_], target : _Integer | _ ? targetQ, args___] := Enclose @ With[
    {qo = ConfirmBy[QuantumOperator[op], QuantumOperatorQ], t = Flatten[{target}]},
    {inside = ContainsAll[qo["InputOrder"], t]},
    ConfirmAssert[
        AllTrue[t, Positive] && DuplicateFreeQ[t] && (inside || Length[t] == Length[qo["InputOrder"]]),
        Message[QuantumMeasurementOperator::luederstarget, t, qo["InputOrder"]]
    ];
    If[ inside,
        luedersMeasurement @ ConfirmBy[QuantumMeasurementOperator[qo, t, args], QuantumMeasurementOperatorQ],
        QuantumMeasurementOperator["Lueders"[QuantumOperator[qo, t]], args]
    ]
]

QuantumMeasurementOperator["Lueders"[op_], args___] := Enclose @ luedersMeasurement @ ConfirmBy[
    QuantumMeasurementOperator[ConfirmBy[QuantumOperator[op], QuantumOperatorQ], args],
    QuantumMeasurementOperatorQ
]

(* The eigenvectors of a projective measurement come from its operator's
   "MatrixRepresentation", whose tensor factors follow the sorted qudit order. So an operator
   stored on an unsorted order, or a target listed out of order, is measured as the sorted
   operator on the sorted target, the same observable on the same qudits, and the result
   carries the sorted target. *)
luedersMeasurement[qmo_] /; Sort[qmo["OutputOrder"]] === Sort[qmo["InputOrder"]] && ! (qmo["QuantumOperator"]["SortedQ"] && OrderedQ[qmo["Target"]]) :=
    luedersMeasurement[QuantumMeasurementOperator[qmo["QuantumOperator"]["Sort"], {Sort @ qmo["Target"]}]]

luedersMeasurement[qmo_] := Enclose @ Block[{trace, observable, inexact, scale, scaled, povm, values, groups, lueders},
    ConfirmAssert[
        qmo["ProjectionQ"] && Sort[qmo["OutputOrder"]] === Sort[qmo["InputOrder"]],
        Message[QuantumMeasurementOperator::luedersnotsquare]
    ];
    (* the operator whose eigenvectors "SuperOperator" measures: qmo traced over the qudits
       outside its targets *)
    trace = DeleteCases[qmo["FullInputOrder"], Alternatives @@ qmo["Target"]];
    observable = Normal[Simplify[QuantumPartialTrace[qmo, trace]]["MatrixRepresentation"]];
    ConfirmAssert[normalQ[observable], Message[QuantumMeasurementOperator::luedersnotnormal]];
    (* "SuperOperator" chops the entries of that operator below 10^-10, which would erase an
       inexact operator whose entries are all small, so such an operator is measured
       multiplied up to largest entry 1, which has the same eigenspaces, and the eigenvalues
       are scaled back *)
    inexact = inexactArrayQ[observable];
    scale = If[inexact, Replace[Min[1, Max[Abs[observable]]], _ ? PossibleZeroQ -> 1], 1];
    scaled = If[scale === 1, qmo, QuantumMeasurementOperator[qmo["QuantumOperator"] / scale, qmo["Targets"]]];
    povm = scaled["POVM"];
    ConfirmAssert[povm["Eigenindex"] === {1}];
    values = First /@ povm["EigenvalueVectors"];
    groups = eigenspacePositions[values, If[inexact, groupingTolerance[Simplify[QuantumPartialTrace[scaled, trace]], values], 0]];
    (* a qudit of dimension 1 has no place in a QuditBasis, so one outcome leaves no eigenqudit *)
    ConfirmAssert[Length[groups] > 1, Message[QuantumMeasurementOperator::luedersoneoutcome]];
    lueders = QuantumMeasurementOperator[
        mergeEigenqudit[povm["Operator"], groups, scale (If[inexact, Mean[values[[#]]], values[[First[#]]]] & /@ groups)],
        qmo["Targets"]
    ];
    ConfirmAssert[completeQ[lueders["Operator"]], Message[QuantumMeasurementOperator::luedersnotcomplete]];
    lueders
]

(* The positions of the eigenvalues, gathered into those of one eigenvalue in order of first
   appearance: an eigenvalue joins the group whose first eigenvalue it equals, exactly when
   exact, after Simplify when symbolic, and to within tolerance when inexact. *)
eigenspacePositions[values_, tolerance_] := With[{equalQ = equalEigenvaluesTest[values, tolerance]},
    Gather[Range[Length[values]], equalQ[values[[#1]], values[[#2]]] &]
]

equalEigenvaluesTest[values_ ? (VectorQ[#, NumericQ] && ! inexactArrayQ[#] &), _] := exactZeroQ[#1 - #2] &

equalEigenvaluesTest[values_ ? (VectorQ[#, NumericQ] &), tolerance_] := Abs[#1 - #2] <= tolerance &

equalEigenvaluesTest[_, _] := Simplify[#1 - #2] === 0 &

(* The distance within which inexact eigenvalue labels count as one, where values are the
   labels "SuperOperator" gives the traced operator t: twice the Frobenius norm of what its
   Chop removes from t, which bounds how far that moves two eigenvalues apart, or the
   roundoff 100 n 10^-p of the largest label at precision p, whichever is more. The labels
   may be the eigenvalues of t divided by a number, so the first term is carried into their
   units by the largest label over the norm of t. *)
groupingTolerance[t_, values_] := With[
    {m = Normal[t["MatrixRepresentation"]], largest = Max[Abs[values]]},
    {norm = Norm[m]},
    Max[
        If[PossibleZeroQ[norm], 0, 2 Norm[Flatten[m - Normal[Chop[t]["MatrixRepresentation"]]]] largest / norm],
        100 Length[values] 10 ^ -SetPrecision[nonzeroPrecision[values], 20] largest
    ]
]

(* Sum_i |lambda(i)> (x) K_i from op = Sum_i |i> (x) K_i, whose first output qudit is the
   eigenqudit: the 0-1 matrix sending each outcome i to the group holding it acts on that
   qudit alone, and labels name the groups. *)
mergeEigenqudit[op_, groups_, labels_] := With[
    {d = First[op["OutputDimensions"]], amplitudes = op["State"]["StateVector"]},
    {group = Lookup[Association @ Catenate @ MapIndexed[Thread[#1 -> First[#2]] &, groups], Range[d]]},
    QuantumOperator[
        QuantumState[
            Flatten[SparseArray[Thread[Transpose[{group, Range[d]}] -> 1], {Length[groups], d}] . ArrayReshape[amplitudes, {d, Length[amplitudes] / d}]],
            QuantumBasis[op["Basis"], "Output" -> QuantumTensorProduct[eigenvalueBasis[labels], QuantumPartialTrace[op["Output"], {1}]]]
        ],
        op["Order"]
    ]
]

eigenvalueBasis[labels_] := QuditBasis[
    MapIndexed[Interpretation[Tooltip[Style[Subscript["\[ScriptCapitalE]", #1], Bold], StringTemplate["Eigenvalue ``"][First @ #2]], {#1, #2}] &, labels],
    IdentityMatrix[Length[labels], SparseArray]
]

(* whether m.m^† == m^†.m, with the roundoff scale of a product of m with itself *)
normalQ[m_] := vanishingQ[m . ConjugateTranspose[m] - ConjugateTranspose[m] . m, m, Norm[m, "Frobenius"] ^ 2]

(* whether the Kraus operators K_lambda of op, a measurement operator whose eigenqudit is its
   first output qudit, give Sum_lambda K_lambda^†.K_lambda = 1: in the computational basis,
   which is orthonormal on the eigenqudit too, the columns of its matrix are orthonormal *)
completeQ[op_] := With[{m = op["MatrixRepresentation"]}, {n = Last[Dimensions[m]]},
    vanishingQ[ConjugateTranspose[m] . m - IdentityMatrix[n, SparseArray], m, Sqrt[n]]
]

(* whether every entry of the array c, computed from the matrix m, is zero: after Simplify
   when m is symbolic, exactly when m is exact, and when m is inexact to within
   100 n 10^-p scale in the Frobenius norm, for the n columns of m and the precision p of
   its nonzero entries *)
vanishingQ[c_, m_, scale_] := Which[
    ! ArrayQ[m, _, NumericQ], AllTrue[Flatten[Normal[Simplify[c]]], PossibleZeroQ],
    ! inexactArrayQ[m], AllTrue[Flatten[Normal[c]], exactZeroQ],
    True, TrueQ[Norm[Flatten[c]] <= 100 Last[Dimensions[m]] 10 ^ -SetPrecision[nonzeroPrecision[m], 20] scale]
]

inexactArrayQ[a_] := ArrayQ[a, _, NumericQ] && nonzeroPrecision[a] < Infinity

(* The precision of the nonzero entries of a numeric array, or of all of them when every one
   is zero: an entry that cancels to zero carries no digits, and Precision of the whole
   array would read the array as having none. *)
nonzeroPrecision[a_] := With[{entries = Flatten[Normal[a]]}, Precision[Replace[Select[entries, # != 0 &], {} -> entries]]]


QuantumMeasurementOperator[name_String, opts___] /; MemberQ[$QuantumMeasurementOperatorNames, name] :=
    QuantumMeasurementOperator[name[], opts]


QuantumMeasurementOperator[name_String[args___], ___] /; ! MemberQ[$QuantumMeasurementOperatorNames, name] := (
    Message[QuantumMeasurementOperator::invalidName, Defer[name[args]]];
    Failure["InvalidName", <|"MessageTemplate" :> QuantumMeasurementOperator::invalidName, "MessageParameters" :> {Defer[name[args]]}|>]
)

QuantumMeasurementOperator[name_String[args___], ___] /; MemberQ[$QuantumMeasurementOperatorNames, name] := (
    Message[QuantumMeasurementOperator::invalidArgs, Defer[name[args]]];
    Failure["InvalidArguments", <|"MessageTemplate" :> QuantumMeasurementOperator::invalidArgs, "MessageParameters" :> {Defer[name[args]]}|>]
)
