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


QECBosonicCode::usage = "QECBosonicCode[{w0, w1}] represents a single-mode bosonic code with the given logical codewords, each an association <|n -> amplitude|> in the Fock basis or a list of {coefficient, amplitude} pairs in the coherent basis.\nQECBosonicCode[\"Binomial\", N, S] builds the binomial code of Michael et al. with parameters N and S.\nQECBosonicCode[\"Cat\", legs, alpha] builds the 2- or 4-component cat code of amplitude alpha.\ncode[prop] gives a property; code[\"Properties\"] lists them.";


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
QECBosonicCode::legs = "Cat codes are built here for 2 or 4 legs; `1` was given.";

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
    If[ ! MemberQ[{2, 4}, legs],
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
        FullSimplify[
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

bosonicCorrectableQ[a_Association, channel_] :=
    AllTrue[
        Flatten @ MapThread[#1 - IdentityMatrix[2] #2 &,
            {bosonicBlocks[a, channel], bosonicKLMatrix[a, channel]}, 2],
        PossibleZeroQ[FullSimplify[#, a["Assumptions"]]] &
    ]

(* Capped: an approximate code never satisfies the conditions at finite alpha. *)
$bosonicMaxOrder = 4;

bosonicCorrectionOrder[a_Association, family_] :=
    LengthWhile[Range[0, $bosonicMaxOrder], bosonicCorrectableQ[a, family[#]] &] - 1


(* The only place a cutoff enters, so codewords can reach the phase-space tools. *)
fockAmplitudes[a_Association, w_, dd_Integer] :=
    If[ a["Basis"] === "Fock",
        Table[Lookup[w, n, 0], {n, 0, dd - 1}],
        Total[#[[1]] Exp[-Abs[#[[2]]]^2/2] Table[#[[2]]^n/Sqrt[n!], {n, 0, dd - 1}] & /@ w]
    ]

bosonicQuantumStates[a_Association, dd_Integer] :=
    QuantumState[Normalize[fockAmplitudes[a, #, dd]], dd] & /@ a["Codewords"]


$bosonicDirectProperties = {
    "Codewords", "Basis", "Modes", "Parameters", "MeanPhotonNumber"
};

$bosonicParametrizedProperties = {
    "ErrorSet", "KnillLaflammeBlocks", "KnillLaflammeMatrix", "CorrectableQ",
    "CorrectionOrder", "QuantumStates"
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
QECBosonicCode[a_Association]["CorrectableQ", channel_] := bosonicCorrectableQ[a, channel]
QECBosonicCode[a_Association]["CorrectionOrder", family_ : "Loss"] :=
    bosonicCorrectionOrder[a, family]
QECBosonicCode[a_Association]["QuantumStates", dd_Integer] := bosonicQuantumStates[a, dd]

QECBosonicCode[_Association][prop_, ___] /;
        ! MemberQ[Join[$bosonicDirectProperties, $bosonicParametrizedProperties, {"Properties"}], prop] :=
    (Message[QECBosonicCode::noprop, prop]; $Failed)


codewordForm[a_Association, w_] :=
    Row[Riffle[
        If[ a["Basis"] === "Fock",
            KeyValueMap[Row[{#2, Ket[{#1}]}] &, w],
            Row[{#[[1]], Ket[{#[[2]]}]}] & /@ w],
        " + "]]

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
