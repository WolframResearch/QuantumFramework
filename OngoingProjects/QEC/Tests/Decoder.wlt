(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Decoder.wlt

   Step 4 of the redesign: one return type for the rate, the post-selected report
   on the detector model, and the decoder as an object.

   The load-bearing tests compare DECISIONS, not rates.  QECDecoder's built-ins
   delegate to the engine functions that already decoded, and a refactor that
   changed a single entry of a table could still leave a rate within tolerance.
   The one that matters most is QEC-Decoder-detector-model-keeps-the-herald-filter:
   the lightest-set rule must not offer rows a verification check rejects, which
   is what moved the verified-cat exponent from 1.73 to 1.965.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

(* The row functions and helpers that the public QECPauli and QECCode["CSS", ...]
   stand on, in PackageScope since step 5 of the API redesign. *)
qecPauliString = Symbol["Wolfram`QuantumFramework`QEC`PackageScope`pauliString"];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecDemTable = Symbol[qecScope <> "demDecoderTable"];

bitFlip = QECCode["BitFlipCode"];
five = QECCode["5QubitCode"];
steane = QECCode["SteaneCode"];
dep = QECNoiseModel["Depolarizing", qecP];
circuit = QECNoiseModel["Circuit", 1/1000];

bitFlipDem = QECDetectorModel[bitFlip, circuit, 1];
fiveTransversal = QECDetectorModel[five, circuit, 1, "Extraction" -> "Transversal"];


(* ============================================================================
   One return type
   ============================================================================ *)

VerificationTest[
    {Head[QECLogicalErrorRate[steane, QECNoiseModel["Depolarizing", 1/10]]],
     Head[QECLogicalErrorRate[bitFlip, circuit, "Rounds" -> 1]]},
    {Rational, Rational},
    TestID -> "QEC-Decoder-rate-is-a-number-at-every-level"
]

VerificationTest[
    Keys[QECLogicalErrorRate[steane, dep, All]],
    {"Rate", "Acceptance", "Failure"},
    TestID -> "QEC-Decoder-report-has-the-same-keys-without-heralds"
]

VerificationTest[
    With[{r = QECLogicalErrorRate[steane, dep, All]},
        {Simplify[r["Rate"] - QECLogicalErrorRate[steane, dep]], r["Acceptance"], Simplify[r["Failure"] - r["Rate"]]}
    ],
    {0, 1, 0},
    TestID -> "QEC-Decoder-report-without-heralds-is-the-rate"
]

VerificationTest[
    Keys[QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", 1/10], 1000, All]],
    {"Rate", "Acceptance", "Accepted", "Shots"},
    TestID -> "QEC-Decoder-sampled-report-keys"
]


(* ============================================================================
   The report on the detector model
   ============================================================================ *)

VerificationTest[
    bitFlipDem["LogicalErrorRate"] === QECLogicalErrorRate[bitFlip, circuit, "Rounds" -> 1],
    True,
    TestID -> "QEC-Decoder-detector-model-rate-is-the-rate"
]

VerificationTest[
    {bitFlipDem["Acceptance"], bitFlipDem["Failure"] === bitFlipDem["LogicalErrorRate"]},
    {1, True},
    TestID -> "QEC-Decoder-detector-model-without-heralds-accepts-everything"
]

VerificationTest[
    MemberQ[bitFlipDem["Properties"], #] & /@ {"LogicalErrorRate", "Acceptance", "Failure"},
    {True, True, True},
    TestID -> "QEC-Decoder-detector-model-lists-its-report"
]

VerificationTest[
    With[{r = bitFlipDem["LogicalErrorRate", 2000, All]},
        {Keys[r], r["Acceptance"], 0 <= r["Rate"] < 0.05}
    ],
    {{"Rate", "Acceptance", "Accepted", "Shots"}, 1., True},
    TestID -> "QEC-Decoder-detector-model-samples"
]

VerificationTest[
    QECDetectorModel[bitFlip, QECNoiseModel["Circuit", qecP], 1]["LogicalErrorRate", 100],
    $Failed,
    {QECLogicalErrorRate::symbolicshots},
    TestID -> "QEC-Decoder-detector-model-refuses-symbolic-shots"
]


(* ============================================================================
   The three built-ins, decision for decision
   ============================================================================ *)

VerificationTest[
    With[{d = QECDecoder[steane]},
        {d["Method"], d["Reach"], d["Table"] === (qecPauliString /@ steane["Decoder"])}
    ],
    {"MinimumWeight", 1, True},
    TestID -> "QEC-Decoder-minimum-weight-is-the-code's-table"
]

VerificationTest[
    With[{noise = QECNoiseModel["Depolarizing", 1/10]},
        Values[QECDecoder[five, noise]["Table"]] === Values[five["Decoder", noise]]
    ],
    True,
    TestID -> "QEC-Decoder-maximum-likelihood-is-the-coset-map"
]

VerificationTest[
    Table[QECDecoder[steane]["Decode", s] === steane["Decode", s], {s, Tuples[{0, 1}, 6][[;; 12]]}],
    ConstantArray[True, 12],
    TestID -> "QEC-Decoder-decode-agrees-with-the-code"
]

(* THE herald test.  The object's table is the engine's filtered table, and the
   filter is doing something: the same model with its heralds ignored (every row
   offered as an explanation) decodes differently. *)
VerificationTest[
    Module[{d = QECDecoder[fiveTransversal], raw, unfiltered},
        raw = qecDemTable[First[fiveTransversal], d["Reach"]];
        unfiltered = qecDemTable[Append[First[fiveTransversal], "Heralds" -> 0], d["Reach"]];
        {fiveTransversal["Heralds"] > 0,
         KeyMap[FromDigits[#, 2] &, FromDigits[#, 2] & /@ d["Table"]] === raw,
         raw =!= unfiltered}
    ],
    {True, True, True},
    TestID -> "QEC-Decoder-detector-model-keeps-the-herald-filter"
]

VerificationTest[
    QECDecoder[five, circuit, "Rounds" -> 1, "Extraction" -> "Transversal"]["Method"],
    "DetectorModel",
    TestID -> "QEC-Decoder-circuit-noise-gives-the-detector-model-decoder"
]


(* ============================================================================
   The rate accepts the object, and the seam
   ============================================================================ *)

VerificationTest[
    {Simplify[QECLogicalErrorRate[steane, dep, "Decoder" -> QECDecoder[steane, dep]] - QECLogicalErrorRate[steane, dep]],
     Simplify[QECLogicalErrorRate[steane, dep, "Decoder" -> QECDecoder[steane]] -
         QECLogicalErrorRate[steane, dep, "Decoder" -> "MinimumWeight"]]},
    {0, 0},
    TestID -> "QEC-Decoder-rate-scores-a-built-in-object"
]

(* A decoder the layer knows nothing about -- here a function reading the table,
   standing in for PyMatching -- is scored by the same machinery. *)
VerificationTest[
    With[{f = QECDecoder[steane, Function[s, Lookup[steane["Decoder"], Key[s], Missing[]]]]},
        {f["Method"], Simplify[QECLogicalErrorRate[steane, dep, "Decoder" -> f] -
            QECLogicalErrorRate[steane, dep, "Decoder" -> "MinimumWeight"]]}
    ],
    {"Function", 0},
    TestID -> "QEC-Decoder-a-function-plugs-in"
]

(* The trivial decoder corrects nothing, so its bit-flip rate is the probability that
   anything happened on the logical class: every odd number of flips reaches it. *)
VerificationTest[
    Simplify[QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", qecP],
        "Decoder" -> QECDecoder[bitFlip, Function[s, If[s === {0, 0}, "III", Missing[]]]]]],
    Simplify[1 - (1 - qecP)^3],
    TestID -> "QEC-Decoder-undecoded-syndromes-count-as-failures"
]

VerificationTest[
    QECDecoder[steane, Function[s, "XX"]]["Decode", {1, 0, 0, 0, 0, 0}],
    $Failed,
    {QECDecoder::badcorr},
    TestID -> "QEC-Decoder-a-function-must-return-a-Pauli-of-the-right-size"
]

VerificationTest[
    QECLogicalErrorRate[steane, dep, "Decoder" -> QECDecoder[five]],
    $Failed,
    {QECDecoder::othercode},
    TestID -> "QEC-Decoder-refuses-another-code's-decoder"
]

VerificationTest[
    QECDecoder[steane]["Decode", {0, 1}],
    $Failed,
    {QECDecoder::badsyn},
    TestID -> "QEC-Decoder-refuses-a-wrong-length-syndrome"
]
