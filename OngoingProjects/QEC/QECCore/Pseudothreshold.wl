(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECPseudothreshold]


(* ============================================================================ *)
(* The pseudothreshold, honestly (Got26 sec. 14.7.6).                           *)
(*                                                                              *)
(* The rate machinery returns an exact expression in the noise parameter, so    *)
(* the point where one level of encoding stops helping -- the logical rate      *)
(* crossing the physical one, p_L(p) = p -- is a root of a polynomial.  It is    *)
(* found by Solve, exactly, as a Root object: no sampling, no fit, the same     *)
(* number on every run, to any precision N is asked for.                        *)
(*                                                                              *)
(* WHAT THE NUMBER IS NOT is the reason this is an object and not a number.     *)
(* Sec. 14.7.6 is emphatic, and its own example is the warning: the 7-qubit     *)
(* code's CNOT pseudothreshold, 2.88 x 10^-5, lies OUTSIDE the true threshold   *)
(* surface, because one rate goes down from level 0 to level 1 while another    *)
(* goes up.  A pseudothreshold compares one level of encoding against none, at  *)
(* one point of noise space, for one decoder, one gadget set and one experiment *)
(* length.  So the result carries every one of those -- the code, the noise     *)
(* model, the decoder, the rounds and the extraction -- together with the       *)
(* polynomial it came from and the caveats that apply to this particular case.  *)
(* There is no path that hands back a bare number.                              *)
(*                                                                              *)
(* A THRESHOLD SURFACE, in this layer's terms, is the same question asked along *)
(* rays: give the noise as a direction with one free parameter                  *)
(* (QECNoiseModel["Circuit", <|"TwoQubit" -> p, "OneQubit" -> c p|>]) and map    *)
(* over c.  Each ray is one root.  It is a pseudothreshold surface -- one level *)
(* against none -- and is reported as one.                                      *)
(*                                                                              *)
(* References: Got26 sec. 14.7.6 (threshold surfaces, eqs. 14.64-14.68).         *)
(* ============================================================================ *)

QECPseudothreshold::usage = "QECPseudothreshold[code, noise] gives the pseudothreshold of the code under a noise model with one symbolic rate p: the rate below which one level of encoding beats no encoding, found exactly as the first root of p_L(p) = p in (0, 1).\nThe result is an object carrying the value together with the code, noise, decoder, rounds, extraction, the logical-error polynomial and the caveats that apply; ps[\"Value\"] is the exact crossing and ps[\"Properties\"] lists the rest.\nThe options are those of QECLogicalErrorRate.";

QECPseudothreshold::params = "The noise must carry exactly one symbolic rate; the logical error rate depends on `1`.";
QECPseudothreshold::rate = "The logical error rate could not be computed exactly for this code and noise.";
QECPseudothreshold::noprop = "`1` is not a property of QECPseudothreshold. Use ps[\"Properties\"] for the list.";

Options[QECPseudothreshold] = Options[QECLogicalErrorRate];


(* ---- construction ---- *)

QECPseudothreshold[code_QECCode, qc_Wolfram`QuantumFramework`QuantumChannel, opts : OptionsPattern[]] :=
    Replace[QECNoiseModel[qc], {noise_QECNoiseModel :> QECPseudothreshold[code, noise, opts], _ -> $Failed}]

QECPseudothreshold[code : QECCode[a_Association], noise : QECNoiseModel[n_Association], opts : OptionsPattern[]] := Module[
    {report, rate, params, p, roots, slope, crossing, level, postSelected},

    report = QECLogicalErrorRate[code, noise, All, opts];
    If[ ! AssociationQ[report], Message[QECPseudothreshold::rate]; Return[$Failed]];
    rate = Simplify[report["Rate"]];

    params = DeleteDuplicates @ Cases[rate, s_Symbol /; ! NumericQ[s] && Context[s] =!= "System`", {0, Infinity}];
    If[ Length[params] =!= 1, Message[QECPseudothreshold::params, params]; Return[$Failed]];
    p = First[params];

    roots = p /. Solve[rate == p && 0 < p < 1, p, Reals];
    roots = If[! FreeQ[roots, p], {}, SortBy[roots, N]];

    (* A pseudothreshold is where encoding STOPS helping, so it only exists if encoding
       helps at low noise: p_L / p -> slope < 1 as p -> 0.  A code whose rate starts at
       order p with a coefficient above 1 -- the bare-ancilla hook error does this --
       may still cross p further up, but that crossing is where encoding starts to
       help, the opposite statement, and it is reported under "Crossings", not as the
       value. *)
    slope = Limit[rate / p, p -> 0, Direction -> "FromAbove"];
    crossing = Which[
        ! TrueQ[slope < 1], Missing["EncodingDoesNotHelp"],
        roots === {}, Missing["NoCrossing"],
        True, First[roots]
    ];

    level = noiseLevel[n];
    postSelected = ! TrueQ[report["Acceptance"] === 1];

    QECPseudothreshold[<|
        "Value" -> crossing,
        "Crossings" -> roots,
        "LowNoiseRatio" -> slope,
        "Parameter" -> p,
        "LogicalErrorRate" -> rate,
        "Acceptance" -> report["Acceptance"],
        "Code" -> a,
        "Noise" -> n,
        "Level" -> level,
        "Decoder" -> OptionValue[QECLogicalErrorRate, {opts}, "Decoder"],
        "Rounds" -> If[level === "CodeCapacity", None, Replace[OptionValue[QECLogicalErrorRate, {opts}, "Rounds"], Automatic :> defaultRounds[a]]],
        "Extraction" -> If[level === "CodeCapacity", None, OptionValue[QECLogicalErrorRate, {opts}, "Extraction"]],
        "PostSelectedQ" -> postSelected
    |>]
]


(* ---- what the number does and does not mean ---- *)

pseudothresholdCaveats[d_Association] := DeleteCases[{
    "A pseudothreshold is not a threshold: it compares one level of encoding against none, at one point of the noise space, and can lie outside the true threshold surface (Got26 sec. 14.7.6).",
    "It belongs to " <>
        If[d["Level"] === "CodeCapacity",
            "this code, this noise model and this decoder",
            "this code, this noise model, this decoder, this extraction and this number of rounds"] <>
        "; change any of them and the number changes.",
    If[ d["Level"] === "CodeCapacity",
        Nothing,
        "At the " <> d["Level"] <> " level the comparison sets the logical failure of a whole memory experiment of " <>
            ToString[d["Rounds"]] <> " round(s) against one physical rate p."
    ],
    If[ TrueQ[d["PostSelectedQ"]],
        "The logical rate is conditional on the heralds accepting; the acceptance is in ps[\"Acceptance\"] and is part of the cost.",
        Nothing
    ],
    Switch[d["Value"],
        Missing["EncodingDoesNotHelp"],
            "As p -> 0 the logical rate is " <> ToString[d["LowNoiseRatio"], InputForm] <>
            " p to leading order, so the encoding does not help at low noise and there is no pseudothreshold; any crossing in ps[\"Crossings\"] is where it starts to help, not where it stops.",
        Missing["NoCrossing"],
            "The encoded rate stays below p on the whole of (0, 1).",
        _, Nothing
    ]
}, Nothing]


(* ---- properties ---- *)

$pseudothresholdProperties = {
    "Value", "Crossings", "LowNoiseRatio", "Parameter", "LogicalErrorRate", "Acceptance", "Code", "Noise", "Level", "Decoder",
    "Rounds", "Extraction", "PostSelectedQ", "Caveats", "Plot", "Properties"
};

QECPseudothreshold[_Association]["Properties"] := $pseudothresholdProperties

QECPseudothreshold[d_Association][prop : "Value" | "Crossings" | "LowNoiseRatio" | "Parameter" | "LogicalErrorRate" | "Acceptance" |
    "Level" | "Decoder" | "Rounds" | "Extraction" | "PostSelectedQ"] := d[prop]
QECPseudothreshold[d_Association]["Code"] := QECCode[d["Code"]]
QECPseudothreshold[d_Association]["Noise"] := QECNoiseModel[d["Noise"]]
QECPseudothreshold[d_Association]["Caveats"] := pseudothresholdCaveats[d]

(* The two curves on log axes, with the crossing marked: below it the encoded rate is
   the lower curve, above it the higher. *)
QECPseudothreshold[d_Association]["Plot", opts___] := With[
    {p = d["Parameter"], rate = d["LogicalErrorRate"], v = d["Value"]},
    With[{hi = If[MissingQ[v], 1/2, Min[1, 10 N[v]]], lo = If[MissingQ[v], 10^-4, N[v] / 100]},
        LogLogPlot[{rate, p}, {p, lo, hi}, opts,
            PlotLegends -> {"logical", "physical"},
            AxesLabel -> {"p", "error rate"},
            Epilog -> If[MissingQ[v], {}, {PointSize[Medium], Point[Log[{N[v], N[v]}]]}],
            ImageSize -> 360
        ]
    ]
]

QECPseudothreshold[_Association][prop_String, ___] :=
    (Message[QECPseudothreshold::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECPseudothreshold /: MakeBoxes[obj : QECPseudothreshold[d_Association] /; KeyExistsQ[d, "Value"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECPseudothreshold,
        obj,
        None,
        {
            BoxForm`SummaryItem[{"Pseudothreshold: ", If[MissingQ[d["Value"]], d["Value"], N[d["Value"]]]}],
            BoxForm`SummaryItem[{"Code: ", {d["Code"]["Qubits"], codeLogicalQubits[d["Code"]]}}],
            BoxForm`SummaryItem[{"Level: ", d["Level"]}]
        },
        {
            BoxForm`SummaryItem[{"Decoder: ", d["Decoder"]}],
            BoxForm`SummaryItem[{"Rounds: ", d["Rounds"]}],
            BoxForm`SummaryItem[{"Not a threshold: ", "see \"Caveats\""}]
        },
        form,
        "Interpretable" -> False
    ]
