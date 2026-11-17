(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageImport["Wolfram`QuantumFramework`"]
PackageImport["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport[QECBosonicCode]

PackageScope[bosonicCodeData]
PackageScope[bosonicElement]
PackageScope[bosonicNorm]
PackageScope[bosonicErrorSet]
PackageScope[bosonicBlocks]
PackageScope[ncDagger]


QECBosonicCode::usage = "QECBosonicCode[{w0, w1}] represents a single-mode bosonic code with the given logical codewords, each an association <|n -> amplitude|> in the Fock basis, a list of {coefficient, amplitude} pairs in the coherent basis, or a pure single-mode QuantumState.\nQECBosonicCode[\"Binomial\"[N, S]] builds the binomial code of Michael et al. with parameters N and S.\nQECBosonicCode[\"Cat\"[legs, alpha]] builds the rotation-symmetric cat code of amplitude alpha on an even number of coherent-state legs; 2d legs correct up to d-1 photon losses.\nQECBosonicCode[spec, b] writes the code in the ladder operators b and b\[Dagger].\ncode[prop] gives a property; code[\"Properties\"] lists them.";


annihilator[a_Association] := a["LadderVariables"][[2, 1]]
creator[a_Association] := a["LadderVariables"][[1, 1]]


ncDagger[a_, e : _Plus | _Times] := ncDagger[a, #] & /@ e
ncDagger[a_, w_NonCommutativeMultiply] := Reverse[ncDagger[a, #] & /@ w]
ncDagger[a_, GeneralizedPower[f_, x_, n_]] := GeneralizedPower[f, ncDagger[a, x], n]
ncDagger[a_, x_] := Replace[x, {annihilator[a] -> creator[a], creator[a] -> annihilator[a], _ :> Conjugate[x]}]

ncTimes[1, y_] := y
ncTimes[x_, 1] := x
ncTimes[x_, y_] := x ** y


bosonicCodeData[QECBosonicCode[a_Association]] := a


QECBosonicCode::words = "Codewords must be two associations <|n -> amplitude|> in the Fock basis, or two lists of {coefficient, amplitude} pairs in the coherent basis.";
QECBosonicCode::legs = "A cat code needs an even number of at least 2 legs; `1` was given.";
QECBosonicCode::vars = "`1` does not name the ladder operators of a single mode.";
QECBosonicCode::states = "QuantumState codewords must be pure states of a single mode.";
QECBosonicCode::channel = "`1` is not an error channel: use \"Loss\"[L], \"Gain\"[G], \"Dephasing\"[D], \"Displacement\"[L], a list of these, or a list of ladder-operator errors.";

Options[QECBosonicCode] = {Assumptions -> True};

fockWordQ[w_] := AssociationQ[w] && AllTrue[Keys[w], IntegerQ[#] && NonNegative[#] &]
coherentWordQ[w_] := MatchQ[w, {{_, _} ..}]

normalizeFock[w_Association] := w/Sqrt[Total[Abs[Values[w]]^2]]

fockWord[qs_QuantumState] :=
    DeleteCases[KeyMap[First[#] - 1 &, Association @ Most @ ArrayRules[qs["Computational"]["StateVector"]]], _ ? PossibleZeroQ]

codeObject[ws_, basis_, params_, vars_, assum_] :=
    With[{lv = Replace[vars, {Automatic :> FieldVariables[], v : _Symbol | {_Symbol} :> FieldVariables[v]}]},
        If[ MatchQ[lv, {{_}, {_}}],
            QECBosonicCode[<|
                "Codewords" -> ws, "Basis" -> basis, "Modes" -> 1, "Parameters" -> params,
                "LadderVariables" -> lv, "Assumptions" -> assum|>],
            Message[QECBosonicCode::vars, vars]; $Failed
        ]
    ]

QECBosonicCode[qs : {_QuantumState, _QuantumState}, rest___] :=
    If[ AllTrue[qs, #["PureStateQ"] && #["Qudits"] == 1 &],
        QECBosonicCode[fockWord /@ qs, rest],
        Message[QECBosonicCode::states]; $Failed
    ]

QECBosonicCode[ws : {_, _}, vars : Except[_Rule | _RuleDelayed] : Automatic, OptionsPattern[]] :=
    Which[
        AllTrue[ws, fockWordQ], codeObject[ws, "Fock", <||>, vars, OptionValue[Assumptions]],
        AllTrue[ws, coherentWordQ], codeObject[ws, "Coherent", <||>, vars, OptionValue[Assumptions]],
        True, Message[QECBosonicCode::words]; $Failed
    ]

QECBosonicCode["Binomial"[nn_Integer ? NonNegative, ss_Integer ? NonNegative],
        vars : Except[_Rule | _RuleDelayed] : Automatic, OptionsPattern[]] :=
    With[{word = normalizeFock[Association[Table[p (ss + 1) -> Sqrt[Binomial[nn + 1, p]], {p, #, nn + 1, 2}]]] &},
        codeObject[{word[0], word[1]}, "Fock", <|"N" -> nn, "S" -> ss|>, vars,
            OptionValue[Assumptions]]
    ]

QECBosonicCode["Cat"[legs_Integer, al_], vars : Except[_Rule | _RuleDelayed] : Automatic, OptionsPattern[]] :=
    If[ ! (EvenQ[legs] && legs >= 2),
        Message[QECBosonicCode::legs, legs]; $Failed,
        codeObject[
            Table[Table[{Exp[-2 Pi I r k/legs], al Exp[2 Pi I k/legs]}, {k, 0, legs - 1}],
                {r, {0, legs/2}}],
            "Coherent", <|"Legs" -> legs, "Alpha" -> al|>, vars, OptionValue[Assumptions]]
    ]


componentElement[a_Association, x_, y_, expr_] :=
    If[ a["Basis"] === "Fock",
        Conjugate[x[[2]]] y[[2]] BosonicMatrixElement[{x[[1]], y[[1]]}, expr, a["LadderVariables"]],
        Conjugate[x[[1]]] y[[1]] BosonicMatrixElement[{x[[2]], y[[2]]}, expr, a["LadderVariables"], "Basis" -> "Coherent"]
    ]

wordComponents[a_Association, w_] := If[a["Basis"] === "Fock", List @@@ Normal[w], w]

rawElement[a_Association, wi_, wj_, expr_] :=
    Total[Outer[componentElement[a, #1, #2, expr] &, wordComponents[a, wi], wordComponents[a, wj], 1], 2]

bosonicNorm[a_Association, i_Integer] := bosonicNorm[a, i] =
    Sqrt[FullSimplify[rawElement[a, a["Codewords"][[i]], a["Codewords"][[i]], 1], a["Assumptions"]]]

bosonicElement[a_Association, i_Integer, j_Integer, expr_] := bosonicElement[a, i, j, expr] =
    With[{ws = a["Codewords"]},
        Chop @ FullSimplify[
            rawElement[a, ws[[i]], ws[[j]], expr]/(bosonicNorm[a, i] bosonicNorm[a, j]),
            a["Assumptions"]]
    ]


bosonicErrorSet[a_, "Loss"[l_Integer ? NonNegative]] := NestList[ncTimes[#, annihilator[a]] &, 1, l]

bosonicErrorSet[a_, "Gain"[g_Integer ? NonNegative]] := NestList[ncTimes[#, creator[a]] &, 1, g]

bosonicErrorSet[a_, "Dephasing"[dd_Integer ? NonNegative]] :=
    NestList[ncTimes[#, creator[a] ** annihilator[a]] &, 1, dd]

bosonicErrorSet[a_, "Displacement"[l_Integer ? NonNegative]] :=
    Flatten @ Table[
        ncTimes[Nest[ncTimes[#, creator[a]] &, 1, p], Nest[ncTimes[#, annihilator[a]] &, 1, t - p]],
        {t, 0, l}, {p, 0, t}
    ]

bosonicErrorSet[a_, chs : {(_String[_]) ..}] := DeleteDuplicates @ Catenate[bosonicErrorSet[a, #] & /@ chs]

bosonicErrorSet[_, es_List] := es

bosonicErrorSet[_, ch_] := (Message[QECBosonicCode::channel, ch]; $Failed)


bosonicBlocks[a_Association, channel_] :=
    With[{es = bosonicErrorSet[a, channel]},
        Table[
            Table[bosonicElement[a, i, j, ncTimes[ncDagger[a, es[[p]]], es[[q]]]], {i, 2}, {j, 2}],
            {p, Length[es]}, {q, Length[es]}]
    ]

bosonicKLMatrix[a_Association, channel_] :=
    Map[FullSimplify[Tr[#]/2, a["Assumptions"]] &, bosonicBlocks[a, channel], {2}]

bosonicResidual[a_Association, channel_] :=
    With[{es = bosonicErrorSet[a, channel]},
        AssociationThread[Tuples[{es, es}],
            Flatten[Map[# - IdentityMatrix[2] FullSimplify[Tr[#]/2, a["Assumptions"]] &,
                bosonicBlocks[a, channel], {2}], 1]]
    ]

bosonicCorrectableQ[a_Association, channel_] :=
    AllTrue[Flatten @ Map[# - IdentityMatrix[2] Tr[#]/2 &, bosonicBlocks[a, channel], {2}],
        PossibleZeroQ[Simplify[#, a["Assumptions"]]] &]

QECBosonicCode::ordcap = "The conditions hold through order `1`, the largest order checked.";

scanOrder[a_Association, test_] :=
    With[{max = If[a["Basis"] === "Fock",
            bosonicFockSpaceSize[a],
            Ceiling[Length[DeleteDuplicates[Catenate[a["Codewords"]][[All, 2]]]]/2]]},
        With[{k = LengthWhile[Range[0, max], test] - 1},
            If[k === max, Message[QECBosonicCode::ordcap, max]];
            k
        ]
    ]

bosonicCorrectionOrder[a_Association, family_] := scanOrder[a, bosonicCorrectableQ[a, family[#]] &]

QECBosonicCode::novar = "This code has no amplitude parameter; supply the limit variable as the third argument.";
QECBosonicCode::numvar = "The amplitude `1` is not a symbol; build the code with a symbolic amplitude to take its large-amplitude limit.";

bosonicApproxOrder[a_Association, family_, var_] :=
    With[{ws = If[a["Basis"] === "Coherent", limitWord[#, var] & /@ a["Codewords"], {}]},
        scanOrder[a,
            If[ ws =!= {} && FreeQ[Keys /@ ws, var],
                limitCorrectableQ[a, ws, family[#], var] &,
                AllTrue[Simplify[Flatten @ Values @ bosonicResidual[a, family[#]], a["Assumptions"]],
                    PossibleZeroQ @ Quiet @ Limit[#, var -> Infinity] &] &
            ]
        ]
    ]

limitWord[w_, var_] := GroupBy[w, Simplify[Last[#]/var] & -> First, Total]

monomialValue[rules_, b_, bc_] := Total[#2 bc^#1[[1]] b^#1[[2]] & @@@ rules]

limitElement[a_Association, wi_Association, wj_Association, expr_, var_] :=
    With[{rules = CoefficientRules[
            BosonicNormalOrder[expr, a["LadderVariables"]] /. GeneralizedPower[_, x_, n_] :> x^n /. NonCommutativeMultiply -> Times,
            {creator[a], annihilator[a]}]},
        Total[Conjugate[wi[#]] wj[#] monomialValue[rules, var #, var Conjugate[#]] & /@ Intersection[Keys[wi], Keys[wj]]]/
            Sqrt[Total[Abs[Values[wi]]^2] Total[Abs[Values[wj]]^2]]
    ]

limitCorrectableQ[a_Association, ws_, channel_, var_] :=
    With[{es = bosonicErrorSet[a, channel]},
        AllTrue[
            Flatten @ Table[
                With[{blk = Table[limitElement[a, ws[[i]], ws[[j]], ncTimes[ncDagger[a, es[[p]]], es[[q]]], var], {i, 2}, {j, 2}]},
                    blk - IdentityMatrix[2] Tr[blk]/2],
                {p, Length[es]}, {q, Length[es]}],
            exactZeroPolynomialQ[#, var] &]
    ]

exactZeroPolynomialQ[x_, var_] := AllTrue[Flatten @ CoefficientList[Expand[x], var], PossibleZeroQ[#, Method -> "ExactAlgebraics"] &]


QECBosonicCode::numeric = "Codeword amplitudes must be numeric to choose a Fock space size; give one explicitly as the second argument.";

bosonicFockSpaceSize[a_Association] :=
    If[ a["Basis"] === "Fock",
        Max[Keys /@ a["Codewords"]] + 1,
        With[{nbar = Max[Abs[#[[2]]]^2 & /@ Catenate[a["Codewords"]]]},
            If[ ! NumericQ[nbar],
                $Failed,
                Ceiling[Quantile[PoissonDistribution[Max[N[nbar], 1]], 1 - 10^-12]] + 5]
        ]
    ]

bosonicCodewordStates[a_Association, dd_Integer] :=
    Total[
        If[ a["Basis"] === "Fock",
            #[[2]] FockState[#[[1]], dd],
            #[[1]] CoherentState[dd][#[[2]]]
        ] & /@ wordComponents[a, #]
    ]["Normalize"] & /@ a["Codewords"]


$bosonicDirectProperties = {
    "Codewords", "Basis", "Modes", "Parameters", "MeanPhotonNumber"
};

$bosonicParametrizedProperties = {
    "ErrorSet", "KnillLaflammeBlocks", "KnillLaflammeMatrix", "KnillLaflammeResidual", "KnillLaflammeRank", "Signature",
    "CorrectableQ", "CorrectionOrder", "ApproximateCorrectionOrder", "CodewordStates",
    "FockSpaceSize"
};

QECBosonicCode::noprop = "`1` is not a property of QECBosonicCode. Use code[\"Properties\"] for the list.";

QECBosonicCode[_Association]["Properties"] :=
    Join[$bosonicDirectProperties, $bosonicParametrizedProperties]

QECBosonicCode[a_Association]["Codewords"] := a["Codewords"]
QECBosonicCode[a_Association]["Basis"] := a["Basis"]
QECBosonicCode[a_Association]["Modes"] := a["Modes"]
QECBosonicCode[a_Association]["Parameters"] := a["Parameters"]

QECBosonicCode[a_Association]["MeanPhotonNumber"] :=
    Table[bosonicElement[a, i, i, creator[a] ** annihilator[a]], {i, 2}]

QECBosonicCode[a_Association]["ErrorSet", channel_] := bosonicErrorSet[a, channel]
QECBosonicCode[a_Association]["KnillLaflammeBlocks", channel_] := bosonicBlocks[a, channel]
QECBosonicCode[a_Association]["KnillLaflammeMatrix", channel_] := bosonicKLMatrix[a, channel]
QECBosonicCode[a_Association]["KnillLaflammeResidual", channel_] := bosonicResidual[a, channel]
QECBosonicCode[a_Association]["KnillLaflammeRank", channel_] := MatrixRank[bosonicKLMatrix[a, channel]]

QECBosonicCode[a_Association]["Signature", channel_] :=
    With[{h = bosonicKLMatrix[a, channel], es = bosonicErrorSet[a, channel]},
        Association[
            es[[{##}]] -> FullSimplify[h[[##]]/Sqrt[h[[#1, #1]] h[[#2, #2]]], a["Assumptions"]] & @@@
                Subsets[Select[Range[Length[es]], ! PossibleZeroQ[h[[#, #]]] &], {2}]]
    ]
QECBosonicCode[a_Association]["CorrectableQ", channel_] := bosonicCorrectableQ[a, channel]
QECBosonicCode[a_Association]["CorrectionOrder", family_ : "Loss"] :=
    bosonicCorrectionOrder[a, family]

QECBosonicCode[a_Association]["ApproximateCorrectionOrder", family_ : "Loss", var_ : Automatic] :=
    With[{v = Replace[var, Automatic :> Lookup[a["Parameters"], "Alpha", $Failed]]},
        Which[
            v === $Failed, Message[QECBosonicCode::novar]; $Failed,
            ! MatchQ[v, _Symbol], Message[QECBosonicCode::numvar, v]; $Failed,
            True, bosonicApproxOrder[a, family, v]
        ]
    ]
QECBosonicCode[a_Association]["CodewordStates", dd_Integer] := bosonicCodewordStates[a, dd]

QECBosonicCode[a_Association]["CodewordStates"] :=
    With[{dd = bosonicFockSpaceSize[a]},
        If[dd === $Failed, Message[QECBosonicCode::numeric]; $Failed, bosonicCodewordStates[a, dd]]]

QECBosonicCode[a_Association]["FockSpaceSize"] := bosonicFockSpaceSize[a]

QECBosonicCode[_Association][prop_, ___] /;
        ! MemberQ[Join[$bosonicDirectProperties, $bosonicParametrizedProperties, {"Properties"}], prop] :=
    (Message[QECBosonicCode::noprop, prop]; $Failed)


signedTerm[c_, k_] :=
    Which[
        c === 1, {" + ", k},
        c === -1, {" - ", k},
        TrueQ[Negative[c]], {" - ", Row[{-c, k}]},
        True, {" + ", Row[{c, k}]}]

codewordForm[a_Association, w_] :=
    Module[{parts},
        parts = Flatten @ If[
            a["Basis"] === "Fock",
            KeyValueMap[signedTerm[Chop[#2], Ket[{#1}]] &, w],
            signedTerm[Chop[#[[1]]], Ket[{Chop[#[[2]]]}]] & /@ w];
        parts[[1]] = If[parts[[1]] === " - ", "-", ""];
        Row[parts]
    ]

$bosonicIcon = With[{hs = {0.18, 0.45, 0.72, 0.99}},
    Graphics[
        {
            {GrayLevel[0.7], AbsoluteThickness[1.2],
             Line[Table[{x, x^2}, {x, -1.02, 1.02, 0.04}]]},
            {RGBColor[0.15, 0.5, 0.65], AbsoluteThickness[2.2],
             Line[{{-Sqrt[#], #}, {Sqrt[#], #}}] & /@ hs}
        },
        ImageSize -> {Automatic, 34},
        PlotRange -> {{-1.1, 1.1}, {-0.05, 1.12}},
        ImagePadding -> 1]
    ];

QECBosonicCode /: MakeBoxes[
        obj : QECBosonicCode[a_Association] /; KeyExistsQ[a, "Codewords"],
        form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECBosonicCode,
        obj,
        $bosonicIcon,
        {
            BoxForm`SummaryItem[{"Basis: ", a["Basis"]}],
            BoxForm`SummaryItem[{"Modes: ", a["Modes"]}],
            BoxForm`SummaryItem[{"Codewords: ", Column[codewordForm[a, #] & /@ a["Codewords"]]}]
        },
        {
            BoxForm`SummaryItem[{"Parameters: ", a["Parameters"]}],
            BoxForm`SummaryItem[{"Mean photon number: ", Dynamic[obj["MeanPhotonNumber"]]}]
        },
        form,
        "Interpretable" -> False
    ]
