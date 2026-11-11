(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECThreshold]

PackageScope[concatenationResources]


(* ============================================================================ *)
(* The threshold of a code family, as the crossing of its curves                *)
(* (Got26 sec. 17.3.5; Dennis, Kitaev, Landahl, Preskill 2002).                  *)
(*                                                                              *)
(* A threshold is a property of a FAMILY, not of a code: the noise rate below   *)
(* which a larger member of the family fails less often than a smaller one,     *)
(* so that growing the code drives the logical rate to zero.  Its signature is  *)
(* that the curves p_L(p) for d = 3, 5, 7, ... cross: below the crossing the    *)
(* largest code is lowest, above it highest.  QECPseudothreshold compares one   *)
(* code against no code; this compares the codes against each other.            *)
(*                                                                              *)
(* TWO ROUTES, chosen by what can be computed.                                  *)
(*  - Exact.  At code capacity, when every member is small enough for the exact *)
(*    rate, each curve is a polynomial and each crossing a root, found by Solve *)
(*    as in QECPseudothreshold.  The repetition code under bit flips crosses at  *)
(*    exactly 1/2, for every pair.                                              *)
(*  - Sampled.  Otherwise each curve is sampled on a grid of p, with the        *)
(*    matching decoder when every member is graph-like.  Each adjacent pair     *)
(*    gives a crossing, read off where the larger code stops being better: the  *)
(*    two grid points that bracket the sign change of the difference are        *)
(*    joined by a line in log p, and the root carries the binomial error of the *)
(*    four rates it came from, propagated by AroundReplace.  With the default   *)
(*    range the grid is first laid coarsely over a range typical of the noise   *)
(*    level, then again, finely and with all the shots, around the bracket the  *)
(*    largest pair found.                                                       *)
(*                                                                              *)
(* WHICH CROSSING IS "THE" THRESHOLD.  Pairwise crossings drift with size --     *)
(* the finite-size effect -- and approach the threshold as the sizes grow.  The *)
(* value reported is the crossing of the two largest members, which is the best *)
(* estimate the data holds; the others are kept in "Crossings" so the drift can *)
(* be seen rather than averaged away.  The error bar is statistical only.       *)
(*                                                                              *)
(* As with the pseudothreshold, the result is an object: a threshold belongs to *)
(* a family, a noise model, a decoder and (below code capacity) an extraction    *)
(* and a number of rounds, and it carries all of them.                          *)
(*                                                                              *)
(* References: Got26 sec. 17.2.3, 17.3.5; Dennis, Kitaev, Landahl, Preskill,     *)
(* J. Math. Phys. 43, 4452 (2002); Wang, Harrington, Preskill, Ann. Phys. 303,   *)
(* 31 (2003).                                                                   *)
(* ============================================================================ *)

QECThreshold::usage = "QECThreshold[{code1, code2, ...}, noise] estimates the threshold of a family of codes under a noise model with one symbolic rate p: the rate at which the logical error curves of successive members cross.\nQECThreshold[{\"Surface\", {3, 5, 7}}, noise] takes the family by name and sizes; \"Repetition\" and \"Toric\" work the same way.\nAt code capacity, when every member is small enough, the crossings are exact roots; otherwise each curve is sampled and the crossings carry error bars.\nThe result is an object; th[\"Value\"] is the crossing of the two largest members and th[\"Properties\"] lists the rest.\nQECThreshold[code, \"Method\" -> \"Concatenation\"] gives instead a rigorous lower bound on the threshold of concatenating the code's fault-tolerant protocol, from its extended rectangles (Got26 sec. 14.5, 14.7.4).";

QECThreshold::params = "The noise must carry exactly one symbolic rate; it depends on `1`.";
QECThreshold::codes = "A threshold needs at least two codes of one family; `1` were given.";
QECThreshold::family = "`1` is not a family QECThreshold can size; use one of `2`, or give the codes themselves.";
QECThreshold::range = "\"Range\" must be Automatic or {pmin, pmax} with 0 < pmin < pmax < 1; got `1`.";
QECThreshold::noprop = "`1` is not a property of QECThreshold. Use th[\"Properties\"] for the list.";

QECThreshold::counting = "\"Counting\" must be \"Malignant\" or \"BadSets\"; got `1`.";
QECThreshold::method = "Without a noise model QECThreshold computes the concatenation bound; \"Method\" must be \"Concatenation\", got `1`.";

Options[QECThreshold] = {
    "Method" -> "Concatenation", "Ancilla" -> "Verified", "Counting" -> "Malignant", "GateSet" -> "Clifford",
    "Shots" -> Automatic, "Range" -> Automatic, "Points" -> 7,
    "Decoder" -> Automatic, "DecoderReach" -> Automatic, "Rounds" -> Automatic, "Extraction" -> "BareAncilla"
};

$thresholdFamilies = {"Repetition", "Surface", "Toric"};

(* The default sample size, and the range a sampled scan starts from at each level:
   the thresholds known for these families sit well inside them. *)
$thresholdShots = 20000;
thresholdDefaultRange["CodeCapacity"] = {1/50, 3/10};
thresholdDefaultRange["Phenomenological"] = {1/200, 1/5};
thresholdDefaultRange["Circuit"] = {1/2000, 1/50};


(* ---- the members ---- *)

(* The family form knows each member's size, and that size is also the memory
   experiment's default number of rounds: the topological-memory convention, and
   it spares the distance search a d = 7 code cannot afford. *)
thresholdMembers[{family_String, sizes : {__Integer}}] := If[
    MemberQ[$thresholdFamilies, family],
    With[{sorted = Sort[DeleteDuplicates[sizes]]},
        <|"Codes" -> (QECCode[family, #] & /@ sorted), "Sizes" -> sorted, "Rounds" -> sorted,
          "Family" -> family|>
    ],
    Message[QECThreshold::family, family, $thresholdFamilies]; $Failed
]

(* Codes given directly are sized by their qubit count. *)
thresholdMembers[codes : {__QECCode}] := With[{order = Ordering[#["Qubits"] & /@ codes]},
    <|"Codes" -> codes[[order]], "Sizes" -> (#["Qubits"] & /@ codes[[order]]),
      "Rounds" -> ConstantArray[Automatic, Length[codes]], "Family" -> None|>
]

thresholdMembers[other_] := (Message[QECThreshold::codes, Length[Flatten[{other}]]]; $Failed)

(* The one symbolic rate the noise carries. *)
noiseSymbols[n_Association] := DeleteDuplicates @ Cases[
    If[noiseLevel[n] === "Circuit", Values[n["Rates"]], {n["Probabilities"], Lookup[n, "MeasurementError", 0]}],
    s_Symbol /; ! NumericQ[s] && Context[s] =!= "System`",
    {0, Infinity}
]

noiseAt[n_Association, p_Symbol, value_] := QECNoiseModel[n /. p -> value]


(* ---- construction ---- *)

(* ---- concatenation: a rigorous lower bound ---- *)

(* THE OTHER MEANING OF THRESHOLD (Def 14.10): the rate below which concatenating
   a fault-tolerant protocol drives the failure to zero.  Level by level a bad
   exRec at level 1 behaves like a single fault at level 2 (level reduction,
   sec. 14.6), so if every exRec fails with probability at most eps(p) < p, each
   level squares the ratio.  The threshold of the protocol is therefore at least
   the fixed point of the worst exRec's bound, eps(p) = p -- eq. 14.63, which
   finds 3.6 x 10^-6 for the book's 7-qubit protocol from its CNOT exRec.

   This is a lower bound three ways over: the union bound counts every set of
   t + 1 locations as fatal, the adversary picks the worst error at each, and the
   protocol is one particular circuit.  It is reported as a bound, not as "the"
   threshold of the code. *)
$concatenationLocations = {"Preparation", "Measurement", "Wait", "H", "S", "CNOT"};

(* the universal set adds the pi/8 gate, done by magic-state injection *)
concatenationLocations["Clifford"] := $concatenationLocations
concatenationLocations["Universal"] := Append[$concatenationLocations, "T"]

QECThreshold::gateset = "\"GateSet\" must be \"Clifford\" or \"Universal\"; got `1`.";

QECThreshold[code_QECCode, opts : OptionsPattern[]] := If[
    OptionValue["Method"] === "Concatenation",
    concatenationThresholdObject[code, opts],
    Message[QECThreshold::method, OptionValue["Method"]]; $Failed
]

concatenationThresholdObject[code_QECCode, opts : OptionsPattern[QECThreshold]] := Module[
    {mode = OptionValue[QECThreshold, {opts}, "Ancilla"], counting = OptionValue[QECThreshold, {opts}, "Counting"],
     gates = OptionValue[QECThreshold, {opts}, "GateSet"], rects, rows, worst, countedBy},
    If[ ! MemberQ[{"Clifford", "Universal"}, gates], Message[QECThreshold::gateset, gates]; Return[$Failed]];
    If[ ! MemberQ[{"Malignant", "BadSets"}, counting], Message[QECThreshold::counting, counting]; Return[$Failed]];
    rects = DeleteCases[
        AssociationMap[Quiet[QECExtendedRectangle[code, #, "Ancilla" -> mode]] &, concatenationLocations[gates]],
        $Failed
    ];
    If[rects === <||>, Return[$Failed]];
    (* the malignant count is only available for t = 1; otherwise the crude bound *)
    If[ counting === "Malignant" && MissingQ[First[rects]["MalignantThresholdBound"]], counting = "BadSets"];
    (* an exRec the frame cannot search -- the pi/8 one -- falls back to bad sets,
       which is also how the book treats it *)
    countedBy = Function[r, If[counting === "Malignant" && ! MissingQ[r["MalignantThresholdBound"]], "Malignant", "BadSets"]] /@ rects;
    rows = AssociationThread[Keys[rects],
        MapThread[If[#2 === "Malignant", #1["MalignantThresholdBound"], #1["ThresholdBound"]] &, {Values[rects], Values[countedBy]}]];
    worst = First[Keys[TakeSmallestBy[rows, N, 1]]];
    QECThreshold[<|
        "Value" -> rows[worst],
        "Method" -> "Concatenation",
        "Counting" -> counting,
        "GateSet" -> gates,
        "Worst" -> worst,
        "ExtendedRectangles" -> (Function[r, <|
            "Locations" -> r["Locations"], "PostSelections" -> r["Verifications"],
            "PreparationLocations" -> r["PreparationLocations"],
            "MalignantPairs" -> If[countedBy[r["Location"]] === "Malignant", r["MalignantPairs"], Missing["NotCounted"]],
            "Counting" -> countedBy[r["Location"]],
            "RectangleLocations" -> r["RectangleLocations"],
            "Bound" -> If[countedBy[r["Location"]] === "Malignant", r["MalignantProbabilityBound", \[FormalP]], r["BadProbabilityBound", \[FormalP]]],
            "FixedPoint" -> rows[r["Location"]]|>] /@ rects),
        "Codes" -> {First[code]},
        "Sizes" -> None,
        "Family" -> "Concatenated",
        "Noise" -> None,
        "Parameter" -> \[FormalP],
        "Level" -> "Circuit",
        "Decoder" -> "Steane EC",
        "Ancilla" -> mode,
        "Shots" -> None, "Range" -> None, "Rounds" -> None, "Extraction" -> None,
        "Crossings" -> <||>, "Rates" -> <||>
    |>]
]

QECThreshold[members_, qc_Wolfram`QuantumFramework`QuantumChannel, opts : OptionsPattern[]] :=
    Replace[QECNoiseModel[qc], {noise_QECNoiseModel :> QECThreshold[members, noise, opts], _ -> $Failed}]

QECThreshold[members_, QECNoiseModel[n_Association], opts : OptionsPattern[]] := Module[
    {m, symbols, p, level, shots, decoder, exactQ},

    m = thresholdMembers[members];
    If[m === $Failed, Return[$Failed]];
    If[ Length[m["Codes"]] < 2, Message[QECThreshold::codes, Length[m["Codes"]]]; Return[$Failed]];

    symbols = noiseSymbols[n];
    If[ Length[symbols] =!= 1, Message[QECThreshold::params, symbols]; Return[$Failed]];
    p = First[symbols];
    level = noiseLevel[n];
    shots = OptionValue["Shots"];

    exactQ = shots === Automatic && level === "CodeCapacity" && AllTrue[m["Codes"], exactEnumerableQ[First[#]] &];

    (* The decoder: the rate's own default on the exact route, which is maximum
       likelihood; matching on the sampled one whenever every member is a graph. *)
    decoder = Replace[OptionValue["Decoder"], Automatic :> If[
        ! exactQ && AllTrue[m["Codes"], codeNonGraphlikeLetter[First[#]] === None &],
        "Matching",
        Automatic
    ]];

    If[ exactQ,
        exactThreshold[m, n, p, decoder, opts],
        sampledThreshold[m, n, p, level, Replace[shots, Automatic -> $thresholdShots], decoder, opts]
    ]
]


(* ---- the exact route ---- *)

exactThreshold[m_Association, n_Association, p_Symbol, decoder_, opts : OptionsPattern[QECThreshold]] := Module[
    {rates, pairs, crossings},
    rates = Simplify[QECLogicalErrorRate[#, QECNoiseModel[n], All, "Decoder" -> decoder,
        "DecoderReach" -> OptionValue[QECThreshold, {opts}, "DecoderReach"]]["Rate"]] & /@ m["Codes"];
    pairs = Partition[Range[Length[rates]], 2, 1];
    crossings = Association @ Table[
        m["Sizes"][[pair]] -> exactCrossing[rates[[pair[[1]]]], rates[[pair[[2]]]], p],
        {pair, pairs}
    ];
    thresholdObject[m, n, p, <|
        "Method" -> "Exact",
        "Rates" -> AssociationThread[m["Sizes"], rates],
        "Crossings" -> crossings,
        "Decoder" -> Replace[decoder, Automatic -> "MaximumLikelihood"],
        "Shots" -> None,
        "Range" -> None,
        "Rounds" -> None,
        "Extraction" -> None
    |>]
]

(* The root where the larger code stops being better: the difference goes from
   negative just below it to positive just above.  A tangency, or a root where the
   order goes the other way, is not a threshold. *)
exactCrossing[small_, large_, p_Symbol] := Module[{diff = Together[large - small], roots, eps},
    roots = p /. Solve[diff == 0 && 0 < p < 1, p, Reals];
    roots = If[! FreeQ[roots, p] || roots === {}, {}, SortBy[roots, N]];
    eps = 10^-6;
    SelectFirst[roots,
        Function[r, (diff /. p -> N[r - eps, 30]) < 0 && (diff /. p -> N[r + eps, 30]) > 0],
        If[ TrueQ[(diff /. p -> 1/1000) < 0], Missing["AboveRange"], Missing["BelowRange"]]
    ]
]


(* ---- the sampled route ---- *)

sampledThreshold[m_Association, n_Association, p_Symbol, level_String, shots_Integer, decoder_,
    opts : OptionsPattern[QECThreshold]] := Module[
    {range, points, grid, rates, crossings, coarse, last, bracket},

    range = OptionValue[QECThreshold, {opts}, "Range"];
    points = OptionValue[QECThreshold, {opts}, "Points"];
    If[ ! (range === Automatic || MatchQ[range, {lo_ ? NumericQ, hi_ ? NumericQ} /; 0 < lo < hi < 1]),
        Message[QECThreshold::range, range]; Return[$Failed]
    ];

    (* With the default range, a coarse pass with a fifth of the shots finds the
       bracket of the largest pair, and the fine pass is laid over it. *)
    If[ range === Automatic,
        range = thresholdDefaultRange[level];
        grid = thresholdGrid[range, points];
        coarse = sampleCurves[m, n, p, grid, Max[1000, Quotient[shots, 5]], decoder, opts];
        If[coarse === $Failed, Return[$Failed]];
        last = Length[m["Sizes"]];
        bracket = crossingBracket[coarse[[last - 1]], coarse[[last]]];
        If[ bracket =!= None,
            range = {grid[[Max[1, bracket - 1]]], grid[[Min[Length[grid], bracket + 2]]]}
        ]
    ];

    grid = thresholdGrid[range, points];
    rates = sampleCurves[m, n, p, grid, shots, decoder, opts];
    If[rates === $Failed, Return[$Failed]];

    crossings = Association @ Table[
        m["Sizes"][[{i, i + 1}]] -> sampledCrossing[grid, rates[[i]], rates[[i + 1]]],
        {i, Length[rates] - 1}
    ];

    thresholdObject[m, n, p, <|
        "Method" -> "Sampled",
        "Rates" -> AssociationThread[m["Sizes"], Function[curve, Transpose[{grid, curve}]] /@ rates],
        "Crossings" -> crossings,
        "Decoder" -> Replace[decoder, Automatic -> "Default"],
        "Shots" -> shots,
        "Range" -> range,
        "Rounds" -> If[level === "CodeCapacity", None,
            MapThread[Replace[Replace[OptionValue[QECThreshold, {opts}, "Rounds"], Automatic -> #1], Automatic :> defaultRounds[First[#2]]] &,
                {m["Rounds"], m["Codes"]}]],
        "Extraction" -> If[level === "CodeCapacity", None, OptionValue[QECThreshold, {opts}, "Extraction"]]
    |>]
]

(* Evenly spaced in log p: a threshold is a crossing of curves that are straight on
   log axes near it. *)
thresholdGrid[{lo_, hi_}, points_Integer] := N @ Exp[Subdivide[Log[lo], Log[hi], points - 1]]

(* One curve per member, each point an Around: the rate and its binomial standard
   error over the accepted shots.  A rate of 0 or 1 still gets the error of one
   event, not none. *)
sampleCurves[m_Association, n_Association, p_Symbol, grid_List, shots_Integer, decoder_, opts : OptionsPattern[QECThreshold]] :=
    Catch @ MapThread[
        Function[{code, rounds},
            With[{sampler = pointSampler[First[code], n, p, Replace[OptionValue[QECThreshold, {opts}, "Rounds"], Automatic -> rounds],
                    shots, decoder, opts]},
                Table[
                    With[{report = sampler[v]},
                        If[ ! AssociationQ[report] || ! NumericQ[report["Rate"]], Throw[$Failed]];
                        With[{r = report["Rate"], k = report["Accepted"]},
                            Around[r, Sqrt[Max[r (1 - r), 1 / k] / k]]
                        ]
                    ],
                    {v, grid}
                ]
            ]
        ],
        {m["Codes"], m["Rounds"]}
    ]

(* What one point of a curve costs.  At code capacity the rate is sampled directly.
   Below it, the detector model is built ONCE with p symbolic and each value is
   substituted into it: the faults and detectors do not depend on p, and building
   the model is nearly all of a point's time (13 of 15 seconds for the d = 5 surface
   code over five rounds). *)
pointSampler[a_Association, n_Association, p_Symbol, rounds_, shots_Integer, decoder_, opts : OptionsPattern[QECThreshold]] /;
        noiseLevel[n] === "CodeCapacity" :=
    Function[v, QECLogicalErrorRate[QECCode[a], noiseAt[n, p, v], shots, All,
        "Decoder" -> decoder, "DecoderReach" -> OptionValue[QECThreshold, {opts}, "DecoderReach"]]]

pointSampler[a_Association, n_Association, p_Symbol, rounds_, shots_Integer, decoder_, opts : OptionsPattern[QECThreshold]] := Module[
    {mode = OptionValue[QECThreshold, {opts}, "Extraction"], dem},
    If[ ! MemberQ[$codeExtractions, mode],
        Message[QECDetectorModel::extraction, mode, $codeExtractions]; Throw[$Failed]
    ];
    dem = codeDetectorModel[a, n, Replace[rounds, Automatic :> defaultRounds[a]], mode];
    Function[v, demSampledReport[dem /. p -> v, shots, decoder, OptionValue[QECThreshold, {opts}, "DecoderReach"]]]
]

(* The index k where the larger code is still better at grid point k and no longer
   at k + 1, or None.  The LAST such change is taken: at low p both rates are a few
   events out of the shots and their difference flips sign on noise alone, while
   above the threshold the larger code is worse at every point. *)
crossingBracket[small_List, large_List] := With[{diff = #["Value"] & /@ (large - small)},
    Replace[
        Last[Position[Partition[diff, 2, 1], {a_, b_} /; a < 0 && b >= 0, {1}, Heads -> False], None],
        {k_Integer} :> k
    ]
]

(* Linear interpolation of the difference in log p between the two bracketing
   points, with the four rates' errors propagated to the root. *)
sampledCrossing[grid_List, small_List, large_List] := Module[{k = crossingBracket[small, large], dLo, dHi, x},
    If[ k === None,
        (* No change of sign from below: where the range ends decides which side the
           threshold is on, since a run that ends with the larger code ahead never
           saw it fall behind. *)
        Return[If[Last[large - small]["Value"] < 0, Missing["AboveRange"], Missing["BelowRange"]]]
    ];
    x = AroundReplace[
        Log[grid[[k]]] + (Log[grid[[k + 1]]] - Log[grid[[k]]]) dLo / (dLo - dHi),
        {dLo -> large[[k]] - small[[k]], dHi -> large[[k + 1]] - small[[k + 1]]}
    ];
    Around[Exp[x["Value"]], Exp[x["Value"]] x["Uncertainty"]]
]


(* ---- the object ---- *)

thresholdObject[m_Association, n_Association, p_Symbol, route_Association] := With[
    {crossings = route["Crossings"]},
    QECThreshold[Join[route, <|
        "Value" -> Last[Values[crossings]],
        "Sizes" -> m["Sizes"],
        "Family" -> m["Family"],
        "Codes" -> (First /@ m["Codes"]),
        "Noise" -> n,
        "Parameter" -> p,
        "Level" -> noiseLevel[n]
    |>]]
]

thresholdCaveats[d_Association] /; d["Method"] === "Concatenation" := {
    "This is a lower bound on the threshold of concatenating this one fault-tolerant protocol (Got26 Def 14.10), not an estimate of it, under adversarial noise at one rate p for every kind of location.",
    If[ d["Counting"] === "Malignant",
        "Pairs of locations are counted as fatal only if some Pauli errors on them make the exRec incorrect in an accepted run (Def 14.7); every triple is still counted as fatal (eq. 14.32).",
        "Every set of t + 1 locations of an exRec is counted as fatal (eq. 14.29); counting only the malignant ones raises the bound (the book: from 3.6 x 10^-6 to 2.7 x 10^-5 for its 7-qubit protocol)."
    ],
    "It comes from the worst extended rectangle, the " <> ToString[d["Worst"]] <> " one, through eps(p) = p (eq. 14.63).",
    "The protocol is this layer's circuit -- Steane EC with verified ancillas, transversal Clifford gates -- laid out by an ASAP schedule; the book's hand-laid circuit has fewer waits, and a different circuit has a different bound.",
    If[ d["GateSet"] === "Universal",
        "The pi/8 gate is done by magic-state injection (figs. 13.7, 14.8b, 14.9); its exRec contains non-Clifford gates, so the Pauli frame cannot search it for malignant sets and it is bounded by bad sets, as in the book.",
        "Clifford gates only: \"GateSet\" -> \"Universal\" adds the pi/8 exRec, which makes the protocol universal and lowers the bound."
    ]
}

thresholdCaveats[d_Association] := DeleteCases[{
    If[ d["Level"] === "CodeCapacity",
        "At code capacity the checks are measured perfectly, so this is the threshold of the code family and the decoder, not a fault-tolerance threshold; a circuit-level threshold is typically an order of magnitude lower.",
        "It belongs to this extraction circuit and these numbers of rounds; at the " <> d["Level"] <> " level both enter the rate."
    ],
    If[ d["Decoder"] === "MaximumLikelihood",
        "It is the threshold of the maximum-likelihood decoder, the best any decoder can do on this noise at code capacity; a practical decoder's is lower.",
        "It belongs to this decoder (" <> ToString[d["Decoder"]] <> "); a better decoder has a higher threshold, maximum likelihood the highest."
    ],
    If[ d["Method"] === "Sampled",
        "The crossings drift with size and approach the threshold from finite sizes; the value is the crossing of the two largest members, and its error bar is statistical only, not this finite-size bias.",
        "The crossings are exact roots of exact rates; with sizes this small they can still differ from the large-size threshold."
    ],
    If[ MissingQ[d["Value"]],
        Switch[d["Value"],
            Missing["AboveRange"], "The largest code was better everywhere sampled: the threshold lies above the range.",
            Missing["BelowRange"], "The largest code was never better where sampled: the threshold lies below the range, or the family has none.",
            _, Nothing
        ],
        Nothing
    ]
}, Nothing]


(* ---- properties ---- *)

$thresholdProperties = {
    "Value", "Crossings", "Rates", "Sizes", "Family", "Codes", "Noise", "Parameter", "Level", "Method",
    "Decoder", "Shots", "Range", "Rounds", "Extraction", "Worst", "ExtendedRectangles", "Counting", "GateSet", "Resources",
    "Caveats", "Plot", "Properties"
};

QECThreshold[_Association]["Properties"] := $thresholdProperties

QECThreshold[d_Association][prop : "Value" | "Crossings" | "Rates" | "Sizes" | "Family" | "Parameter" | "Level" |
    "Method" | "Decoder" | "Shots" | "Range" | "Rounds" | "Extraction" | "Worst" | "ExtendedRectangles" | "Counting" | "GateSet"] :=
    Lookup[d, prop, Missing["NotApplicable"]]
QECThreshold[d_Association]["Codes"] := QECCode /@ d["Codes"]
QECThreshold[d_Association]["Noise"] := If[AssociationQ[d["Noise"]], QECNoiseModel[d["Noise"]], Missing["NotApplicable"]]
QECThreshold[d_Association]["Caveats"] := thresholdCaveats[d]

(* The curves on log axes, the threshold marked by a dashed line. *)
(* WHAT IT COSTS (sec. 14.7.5).  At rate p below the threshold, l levels of
   concatenation fail per location with probability at most p_T (p / p_T)^(2^l), so a
   circuit of up to eps / that many locations fails with probability at most eps;
   and every level multiplies the size by the largest rectangle B -- the
   rectangle, not the exRec, because FT(C) is cut into rectangles without overlap.
   Table 14.1 is this with the book's B = 492 and p_T = 2.7 x 10^-5. *)
concatenationResources[pT_, B_Integer, p_, eps_, levels_Integer] := Table[
    <|"Levels" -> l, "MaxLocations" -> eps / (pT (p / pT)^(2^l)), "Overhead" -> B^l|>,
    {l, 0, levels}
]

QECThreshold[d_Association]["Resources", p_ ? NumericQ, eps_ ? NumericQ, levels_Integer : 5] /; d["Method"] === "Concatenation" :=
    Dataset @ concatenationResources[N[d["Value"]],
        Max[#["RectangleLocations"] & /@ Values[d["ExtendedRectangles"]]], p, eps, levels]

(* The concatenation bound: each exRec's eps(p) against p, the fixed point of the
   worst marked. *)
QECThreshold[d_Association]["Plot", opts___] /; d["Method"] === "Concatenation" := With[
    {bounds = #["Bound"] & /@ d["ExtendedRectangles"], v = N[d["Value"]]},
    LogLogPlot[Evaluate[Append[Values[bounds], \[FormalP]]], {\[FormalP], v / 100, 3 v}, opts,
        PlotLegends -> Append[Keys[bounds], "p"],
        PlotStyle -> Append[Table[Automatic, Length[bounds]], Directive[Dashed, Gray]],
        AxesLabel -> {"p", "exRec failure bound"},
        Epilog -> {PointSize[Medium], Point[Log[{v, v}]]}, ImageSize -> 380]
]

QECThreshold[d_Association]["Plot", opts___] := With[
    {v = d["Value"], sizes = d["Sizes"], label = If[d["Family"] === None, "n = ", "d = "]},
    With[{mark = If[MissingQ[v], {}, {Dashed, Gray, InfiniteLine[{{Log[N[#]], 0}, {Log[N[#]], 1}}] & @ If[Head[v] === Around, v["Value"], v]}]},
        If[ d["Method"] === "Exact",
            LogLogPlot[Evaluate[Values[d["Rates"]] /. d["Parameter"] -> \[FormalP]], {\[FormalP], 1/100, 1/2}, opts,
                PlotLegends -> (label <> ToString[#] & /@ sizes),
                AxesLabel -> {"p", "logical error rate"}, Epilog -> mark, ImageSize -> 380],
            ListLogLogPlot[Values[d["Rates"]], opts,
                Joined -> True, PlotMarkers -> Automatic,
                PlotLegends -> (label <> ToString[#] & /@ sizes),
                AxesLabel -> {"p", "logical error rate"}, Epilog -> mark, ImageSize -> 380]
        ]
    ]
]

QECThreshold[_Association][prop_String, ___] :=
    (Message[QECThreshold::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECThreshold /: MakeBoxes[obj : QECThreshold[d_Association] /; KeyExistsQ[d, "Value"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECThreshold,
        obj,
        None,
        {
            BoxForm`SummaryItem[{"Threshold: ", Replace[d["Value"], r : Except[_Around | _Missing] :> N[r]]}],
            BoxForm`SummaryItem[{If[d["Family"] === None, "Sizes (qubits): ", d["Family"] <> ": "], d["Sizes"]}],
            BoxForm`SummaryItem[{"Level: ", d["Level"]}]
        },
        {
            BoxForm`SummaryItem[{"Method: ", d["Method"]}],
            BoxForm`SummaryItem[{"Decoder: ", d["Decoder"]}],
            BoxForm`SummaryItem[{"Shots per point: ", d["Shots"]}]
        },
        form,
        "Interpretable" -> False
    ]
