(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Matching.wlt

   The native minimum-weight perfect matching decoder (Got26 sec. 17.2),
   QECDecoder[..., "Method" -> "Matching"] and "Decoder" -> "Matching".

   What each group is checked against.
   - The matcher itself: on codes small enough to decode by every other means,
     matching must give the rate those give -- the binomial tail on the
     repetition code, the maximum-likelihood polynomial on the surface code
     under bit flips -- and every one of the 256 syndromes of the d = 3 surface
     code must come back as a correction with that syndrome.
   - The graph it builds from a detector model: on a CSS code a fault is split
     into its X part and its Z part, so no edge joins an X-check detector to a
     Z-check detector; and the decoder on that graph reproduces PyMatching on
     the same circuit (StimCrossCheck, November 2026: Z-bar fails in 0.443 % of
     shots of the d = 3 surface code at circuit level 1/500, three rounds).
   - What it refuses: a code where one single-qubit error trips three checks is
     a hypergraph, and the refusal is the same at every noise level.
   - That it scales: past d = 3, where nothing in the layer could decode before,
     the bit-flip rates of d = 3 and d = 5 cross.

   Everything here runs on bare Wolfram; the PyMatching number is pinned, not
   recomputed.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";
qecDemMatching = Symbol[qecScope <> "demMatching"];
qecDemRowKeys = Symbol[qecScope <> "demRowKeys"];
qecDemGroups = Symbol[qecScope <> "demGroups"];
qecDecoderLookup = Symbol[qecScope <> "decoderLookup"];
qecMatchingLookup = Symbol[qecScope <> "matchingLookup"];
qecPauliVector = Symbol[qecScope <> "pauliVector"];
qecPauliString = Symbol[qecScope <> "pauliString"];

rep5 = QECCode["Repetition", 5];
surface3 = QECCode["Surface", 3];
bitFlip = QECNoiseModel["BitFlip", qecP];
circuit500 = QECNoiseModel["Circuit", 1/500];

surfaceDecoder = QECDecoder[surface3, QECNoiseModel["Depolarizing", 1/100], "Method" -> "Matching"];
surfaceModel = QECDetectorModel[surface3, circuit500, 3];

(* Sample a detector model the way demSampledFailure does, decode with matching,
   and give the failure rate of each observable bit separately: bit 0 is the Zbar
   product (flipped by X errors, what a memory-Z experiment reads), bit 1 the Xbar
   one. *)
perObservableFailure[dem_QECDetectorModel, count_Integer] := Module[
    {a = First[dem], d, o, h, dKeys, oKeys, detectors, observables, predicted},
    {d, o, h} = qecDemRowKeys[a];
    dKeys = Prepend[d, 0];
    oKeys = Prepend[o, 0];
    {detectors, observables} = Transpose @ Map[
        Function[rows,
            With[{p = N[a["Probabilities"][[rows]]]},
                With[{draw = RandomChoice[Append[p, 1 - Total[p]] -> Append[rows, 0], count]},
                    {dKeys[[draw + 1]], oKeys[[draw + 1]]}
                ]
            ]
        ],
        qecDemGroups[a]
    ];
    detectors = Fold[BitXor, 0, detectors];
    observables = Fold[BitXor, 0, observables];
    predicted = Lookup[qecDecoderLookup[qecMatchingLookup[qecDemMatching[a]], detectors], detectors, -1];
    Table[
        N @ Mean[MapThread[Boole[BitAnd[#1, 2^b] =!= BitAnd[#2, 2^b]] &, {predicted, observables}]],
        {b, 0, 1}
    ]
]


(* ============================================================================
   The matcher, against decoders it must agree with
   ============================================================================ *)

(* Matching on a line is majority voting: the rate is the binomial tail. *)
VerificationTest[
    Simplify[
        QECLogicalErrorRate[rep5, bitFlip, "Decoder" -> "Matching"] -
            Sum[Binomial[5, j] qecP^j (1 - qecP)^(5 - j), {j, 3, 5}]
    ],
    0,
    TestID -> "QEC-Matching-repetition-code-is-the-binomial-tail"
]

(* Under bit flips the surface code's X errors are a graph problem exactly, and
   matching is the maximum-likelihood decoder's equal. *)
VerificationTest[
    Simplify[
        QECLogicalErrorRate[surface3, bitFlip, "Decoder" -> "Matching"] -
            QECLogicalErrorRate[surface3, bitFlip]
    ],
    0,
    TestID -> "QEC-Matching-surface-bit-flip-equals-maximum-likelihood"
]

VerificationTest[
    AllTrue[Tuples[{0, 1}, 8], surface3["Syndrome", surfaceDecoder["Decode", #]] === # &],
    True,
    TestID -> "QEC-Matching-every-syndrome-gets-a-matching-correction"
]

(* Every single-qubit X and Z is corrected up to a stabilizer. *)
VerificationTest[
    Module[{n = 9, singles},
        singles = Flatten[Table[StringJoin[ReplacePart[ConstantArray["I", n], q -> l]], {q, n}, {l, {"X", "Z"}}]];
        AllTrue[singles,
            Function[e, With[{c = surfaceDecoder["Decode", surface3["Syndrome", e]]},
                surface3["StabilizerMemberQ", qecPauliString[Append[
                    BitXor[Most[qecPauliVector[e]], Most[qecPauliVector[c]]], 0]]]
            ]]
        ]
    ],
    True,
    TestID -> "QEC-Matching-corrects-every-single-qubit-error"
]

VerificationTest[
    {surfaceDecoder["Method"], Head[surfaceDecoder["Graph"]], Length[surfaceDecoder["Table"]],
     surfaceDecoder["Decode", {0, 0, 0, 0, 1, 0, 0, 0}]},
    {"Matching", Graph, 256, "XIIIIIIII"},
    TestID -> "QEC-Matching-decoder-object"
]


(* ============================================================================
   The graph from a detector model
   ============================================================================ *)

(* Detector (round, i) is check i; surface-code checks 1-4 are X type, 5-8 Z. *)
VerificationTest[
    With[{m = qecDemMatching[First[surfaceModel]]},
        With[{zType = Function[det, Mod[det - 1, 8] + 1 > 4]},
            {m["Detectors"], m["DroppedMechanisms"],
             Count[Keys[m["Edges"]], {u_, v_} /; v <= 32 && zType[u] =!= zType[v]]}
        ]
    ],
    {32, 0, 0},
    TestID -> "QEC-Matching-detector-graph-splits-X-from-Z"
]

(* PyMatching on Stim's graph of the same circuit: 0.00443 (two million shots).
   200 000 shots puts the standard error near 0.00015, so the window is four of it. *)
VerificationTest[
    BlockRandom[SeedRandom[11];
        With[{f = perObservableFailure[surfaceModel, 200000]},
            {0.0038 < f[[1]] < 0.0050, f[[2]] < 0.006}
        ]
    ],
    {True, True},
    TestID -> "QEC-Matching-circuit-level-matches-PyMatching"
]

VerificationTest[
    With[{dec = QECDecoder[surfaceModel, "Method" -> "Matching"]},
        {dec["Method"], dec["Decode", ConstantArray[0, 32]], Head[dec["DetectorModel"]]}
    ],
    {"Matching", {0, 0}, QECDetectorModel},
    TestID -> "QEC-Matching-detector-model-decoder-object"
]

VerificationTest[
    BlockRandom[SeedRandom[12];
        0.004 < QECLogicalErrorRate[surface3, circuit500, 50000, "Rounds" -> 3, "Decoder" -> "Matching"] < 0.012
    ],
    True,
    TestID -> "QEC-Matching-rate-through-the-public-entry-point"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

VerificationTest[
    QECDecoder[QECCode["SteaneCode"], QECNoiseModel["Depolarizing", 1/100], "Method" -> "Matching"],
    $Failed,
    {QECDecoder::notgraphlike},
    TestID -> "QEC-Matching-refuses-Steane-at-code-capacity"
]

VerificationTest[
    QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Circuit", 1/500], 1000, "Rounds" -> 2, "Decoder" -> "Matching"],
    $Failed,
    {QECDecoder::notgraphlike},
    TestID -> "QEC-Matching-refuses-Steane-at-circuit-level"
]

(* Its X's and Z's trip two checks each; its Y's trip three. *)
VerificationTest[
    QECDecoder[QECCode["5QubitCode"], QECNoiseModel["Depolarizing", 1/100], "Method" -> "Matching"],
    $Failed,
    {QECDecoder::notgraphlike},
    TestID -> "QEC-Matching-refuses-the-5-qubit-code-for-its-Y"
]


(* ============================================================================
   Past d = 3
   ============================================================================ *)

(* Bit flips at code capacity: below the threshold (about 10.3 %) the larger code
   is better, above it worse. *)
VerificationTest[
    BlockRandom[SeedRandom[13];
        Table[
            Less @@ Table[QECLogicalErrorRate[QECCode["Surface", d], QECNoiseModel["BitFlip", p], 20000, "Decoder" -> "Matching"], {d, {5, 3}}],
            {p, {1/25, 4/25}}
        ]
    ],
    {True, False},
    TestID -> "QEC-Matching-surface-codes-cross-under-bit-flips"
]
