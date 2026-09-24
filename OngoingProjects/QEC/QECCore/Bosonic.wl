(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageImport["Wolfram`QuantumFramework`"]
PackageImport["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport[QECBosonicCode]

PackageScope[bosonicCodeData]
PackageScope[bosonicElement]
PackageScope[bosonicErrorSet]
PackageScope[bosonicBlocks]
PackageScope[ncDagger]


QECBosonicCode::usage = "QECBosonicCode[{w0, w1}] represents a single-mode bosonic code with the given logical codewords, each an association <|n -> amplitude|> in the Fock basis or a list of {coefficient, amplitude} pairs in the coherent basis.\nQECBosonicCode[\"Binomial\", N, S] builds the binomial code of Michael et al. with parameters N and S.\nQECBosonicCode[\"Cat\", legs, alpha] builds the rotation-symmetric cat code of amplitude alpha on an even number of coherent-state legs; 2d legs correct up to d-1 photon losses.\ncode[prop] gives a property; code[\"Properties\"] lists them.";


$av := First[FieldVariables[]]
$adv := Last[FieldVariables[]]


daggerAtom[SuperDagger[v_]] := v
daggerAtom[v_] := SuperDagger[v]

ncDagger[e_Plus] := ncDagger /@ e
ncDagger[HoldPattern[Times[c_, w_NonCommutativeMultiply]]] := Conjugate[c] Reverse[daggerAtom /@ w]
ncDagger[w_NonCommutativeMultiply] := Reverse[daggerAtom /@ w]
ncDagger[c_ ? NumericQ] := Conjugate[c]
ncDagger[v_] := daggerAtom[v]

ncTimes[1, y_] := y
ncTimes[x_, 1] := x
ncTimes[x_, y_] := x ** y


bosonicCodeData[QECBosonicCode[a_Association]] := a


QECBosonicCode::words = "Codewords must be two associations <|n -> amplitude|> in the Fock basis, or two lists of {coefficient, amplitude} pairs in the coherent basis.";
QECBosonicCode::legs = "A cat code needs an even number of at least 2 legs; `1` was given.";

Options[QECBosonicCode] = {Assumptions -> True};

fockWordQ[w_] := AssociationQ[w] && AllTrue[Keys[w], IntegerQ[#] && NonNegative[#] &]
coherentWordQ[w_] := MatchQ[w, {{_, _} ..}]

normalizeFock[w_Association] := w/Sqrt[Total[Abs[Values[w]]^2]]

codeObject[ws_, basis_, params_, assum_] :=
    QECBosonicCode[<|
        "Codewords" -> ws, "Basis" -> basis, "Modes" -> 1,
        "Parameters" -> params, "Assumptions" -> assum|>]

QECBosonicCode[ws : {_, _}, OptionsPattern[]] :=
    Which[
        AllTrue[ws, fockWordQ], codeObject[ws, "Fock", <||>, OptionValue[Assumptions]],
        AllTrue[ws, coherentWordQ], codeObject[ws, "Coherent", <||>, OptionValue[Assumptions]],
        True, Message[QECBosonicCode::words]; $Failed
    ]

QECBosonicCode["Binomial", nn_Integer ? NonNegative, ss_Integer ? NonNegative,
        OptionsPattern[]] :=
    With[{word = Function[par,
            normalizeFock[Association[Table[
                p (ss + 1) -> Sqrt[Binomial[nn + 1, p]], {p, par, nn + 1, 2}]]]]},
        codeObject[{word[0], word[1]}, "Fock", <|"N" -> nn, "S" -> ss|>,
            OptionValue[Assumptions]]
    ]

QECBosonicCode["Cat", legs_Integer, al_, OptionsPattern[]] :=
    If[ ! (EvenQ[legs] && legs >= 2),
        Message[QECBosonicCode::legs, legs]; $Failed,
        codeObject[
            Table[Table[{Exp[-2 Pi I r k/legs], al Exp[2 Pi I k/legs]}, {k, 0, legs - 1}],
                {r, {0, legs/2}}],
            "Coherent", <|"Legs" -> legs, "Alpha" -> al|>, OptionValue[Assumptions]]
    ]


componentElement[a_Association, x_, y_, expr_] :=
    If[ a["Basis"] === "Fock",
        Conjugate[x[[2]]] y[[2]] BosonicMatrixElement[{x[[1]], y[[1]]}, expr],
        Conjugate[x[[1]]] y[[1]] BosonicMatrixElement[{x[[2]], y[[2]]}, expr, "Basis" -> "Coherent"]
    ]

wordComponents[a_Association, w_] := If[a["Basis"] === "Fock", List @@@ Normal[w], w]

rawElement[a_Association, wi_, wj_, expr_] :=
    Total[Flatten @ Outer[componentElement[a, #1, #2, expr] &,
        wordComponents[a, wi], wordComponents[a, wj], 1]]

wordNorm[a_Association, w_] := Sqrt[rawElement[a, w, w, 1]]

bosonicElement[a_Association, i_Integer, j_Integer, expr_] :=
    With[{ws = a["Codewords"]},
        Chop @ FullSimplify[
            rawElement[a, ws[[i]], ws[[j]], expr]/(wordNorm[a, ws[[i]]] wordNorm[a, ws[[j]]]),
            a["Assumptions"]]
    ]


(* The idealized span {I, a, ..., a^L}, not the Kraus set of the physical channel:
   its no-jump factor eta^(n/2) is not a polynomial in the ladder operators. *)
bosonicErrorSet["Loss"[l_Integer ? NonNegative]] := NestList[ncTimes[#, $av] &, 1, l]

bosonicErrorSet["Dephasing"[dd_Integer ? NonNegative]] :=
    NestList[ncTimes[#, $adv ** $av] &, 1, dd]

bosonicErrorSet[es_List] := es


bosonicBlocks[a_Association, channel_] :=
    With[{es = bosonicErrorSet[channel]},
        Table[
            Table[bosonicElement[a, i, j, ncTimes[ncDagger[es[[p]]], es[[q]]]], {i, 2}, {j, 2}],
            {p, Length[es]}, {q, Length[es]}]
    ]

bosonicKLMatrix[a_Association, channel_] :=
    Map[FullSimplify[Tr[#]/2, a["Assumptions"]] &, bosonicBlocks[a, channel], {2}]

(* The obstruction itself: what each block has left over once the codeword-independent
   part is removed. Zero exactly when the conditions hold. *)
bosonicResidual[a_Association, channel_] :=
    With[{b = bosonicBlocks[a, channel]},
        Flatten @ MapThread[#1 - IdentityMatrix[2] #2 &,
            {b, Map[FullSimplify[Tr[#]/2, a["Assumptions"]] &, b, {2}]}, 2]
    ]

bosonicCorrectableQ[a_Association, channel_] :=
    AllTrue[bosonicResidual[a, channel],
        PossibleZeroQ[FullSimplify[#, a["Assumptions"]]] &]

(* Capped: an approximate code never satisfies the conditions at finite alpha. *)
$bosonicMaxOrder = 4;

bosonicCorrectionOrder[a_Association, family_] :=
    LengthWhile[Range[0, $bosonicMaxOrder], bosonicCorrectableQ[a, family[#]] &] - 1

(* An approximate code satisfies the conditions only in the limit of large amplitude.
   Asking whether every residual vanishes there separates a code that is approximately
   correcting from one that is not correcting at all. *)
QECBosonicCode::novar = "This code has no amplitude parameter; supply the limit variable as the third argument.";

bosonicApproxOrder[a_Association, family_, var_] :=
    LengthWhile[
        Range[0, $bosonicMaxOrder],
        AllTrue[Simplify[bosonicResidual[a, family[#]], a["Assumptions"]],
            PossibleZeroQ @ Quiet @ Limit[#, var -> Infinity] &] &
    ] - 1


(* The only place a cutoff enters, so codewords can reach the phase-space tools.
   A Fock codeword is a finite sum and its size is exact; a coherent one needs enough
   levels to hold the Poisson tail of its largest amplitude. *)
QECBosonicCode::numeric = "Codeword amplitudes must be numeric to choose a Fock space size; give one explicitly as the second argument.";

bosonicFockSpaceSize[a_Association] :=
    If[ a["Basis"] === "Fock",
        Max[Union @@ (Keys /@ a["Codewords"])] + 1,
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
    "ErrorSet", "KnillLaflammeBlocks", "KnillLaflammeMatrix", "KnillLaflammeResidual",
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
    Table[bosonicElement[a, i, i, $adv ** $av], {i, 2}]

QECBosonicCode[a_Association]["ErrorSet", channel_] := bosonicErrorSet[channel]
QECBosonicCode[a_Association]["KnillLaflammeBlocks", channel_] := bosonicBlocks[a, channel]
QECBosonicCode[a_Association]["KnillLaflammeMatrix", channel_] := bosonicKLMatrix[a, channel]
QECBosonicCode[a_Association]["KnillLaflammeResidual", channel_] := bosonicResidual[a, channel]
QECBosonicCode[a_Association]["CorrectableQ", channel_] := bosonicCorrectableQ[a, channel]
QECBosonicCode[a_Association]["CorrectionOrder", family_ : "Loss"] :=
    bosonicCorrectionOrder[a, family]

QECBosonicCode[a_Association]["ApproximateCorrectionOrder", family_ : "Loss", var_ : Automatic] :=
    With[{v = Replace[var, Automatic :> Lookup[a["Parameters"], "Alpha", $Failed]]},
        If[ v === $Failed,
            Message[QECBosonicCode::novar]; $Failed,
            bosonicApproxOrder[a, family, v]
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

(* Fixed: the harmonic well and its ladder of levels, the one mode a bosonic code
   lives in. *)
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
