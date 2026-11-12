(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Threshold.wlt

   The threshold of a code family as the crossing of its curves (Got26 sec.
   17.3.5), QECThreshold.

   What each value is checked against.
   - Exact: the repetition code under bit flips.  Its rates are binomial tails,
     and r5 - r3 = 3 p^2 (p - 1)^2 (2 p - 1), so every pair crosses at exactly
     1/2 and nowhere else in (0, 1).
   - Sampled, same family: the estimator has to find 1/2 again from samples,
     within its own error bar.
   - Sampled, against the literature: the repetition code with measurement
     errors at the same rate is the toric code at code capacity in disguise
     (Dennis, Kitaev, Landahl, Preskill 2002), whose matching threshold is
     10.3 %; the surface code under bit flips has the same one.  Small sizes
     sit a little below it, which is the finite-size drift the object reports.
   - The detector model is built once with p symbolic and substituted at each
     point; that has to give the same samples, shot for shot, as building it at
     each point.

   Everything runs on bare Wolfram, seeded.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";
qecCodeDetectorModel = Symbol[qecScope <> "codeDetectorModel"];
qecDemSampledReport = Symbol[qecScope <> "demSampledReport"];

bitFlip = QECNoiseModel["BitFlip", qecP];
phenomenological = QECNoiseModel["BitFlip", qecP, "MeasurementError" -> qecP];

exactRepetition = QECThreshold[{"Repetition", {3, 5, 7}}, bitFlip];


(* ============================================================================
   Exact
   ============================================================================ *)

VerificationTest[
    {exactRepetition["Value"], exactRepetition["Crossings"], exactRepetition["Method"], exactRepetition["Decoder"]},
    {1/2, <|{3, 5} -> 1/2, {5, 7} -> 1/2|>, "Exact", "MaximumLikelihood"},
    TestID -> "QEC-Threshold-repetition-code-crosses-at-one-half-exactly"
]

VerificationTest[
    Simplify[exactRepetition["Rates"][5] - QECLogicalErrorRate[QECCode["Repetition", 5], bitFlip]],
    0,
    TestID -> "QEC-Threshold-exact-curves-are-the-rates"
]

(* Codes given directly are put in order of size. *)
VerificationTest[
    With[{th = QECThreshold[{QECCode["Repetition", 5], QECCode["Repetition", 3]}, bitFlip]},
        {th["Sizes"], th["Value"], th["Family"]}
    ],
    {{3, 5}, 1/2, None},
    TestID -> "QEC-Threshold-codes-given-directly"
]


(* ============================================================================
   Sampled
   ============================================================================ *)

VerificationTest[
    BlockRandom[SeedRandom[3];
        With[{th = QECThreshold[{"Repetition", {3, 5, 7}}, bitFlip, "Shots" -> 20000, "Range" -> {0.3, 0.7}]},
            {th["Method"], Head[th["Value"]], Abs[th["Value"]["Value"] - 1/2] < 3 th["Value"]["Uncertainty"],
             th["Value"]["Uncertainty"] < 0.03}
        ]
    ],
    {"Sampled", Around, True, True},
    TestID -> "QEC-Threshold-sampled-estimator-finds-the-exact-crossing"
]

(* The repetition code with measurement errors at the rate of the data errors: the
   matching threshold of the random-bond Ising model, 10.3 %. *)
VerificationTest[
    BlockRandom[SeedRandom[4];
        With[{th = QECThreshold[{"Repetition", {3, 5}}, phenomenological, "Shots" -> 4000, "Points" -> 5]},
            {th["Level"], th["Rounds"], th["Decoder"], 0.08 < th["Value"]["Value"] < 0.13}
        ]
    ],
    {"Phenomenological", {3, 5}, "Matching", True},
    TestID -> "QEC-Threshold-repetition-with-measurement-errors-near-ten-percent"
]

VerificationTest[
    BlockRandom[SeedRandom[1];
        With[{th = QECThreshold[{"Surface", {3, 5}}, bitFlip, "Shots" -> 5000]},
            {th["Decoder"], 0.06 < th["Value"]["Value"] < 0.12}
        ]
    ],
    {"Matching", True},
    TestID -> "QEC-Threshold-surface-code-under-bit-flips"
]

(* A range entirely below the threshold: the larger code wins everywhere, and the
   object says so instead of inventing a crossing. *)
VerificationTest[
    BlockRandom[SeedRandom[2];
        With[{th = QECThreshold[{"Surface", {3, 5}}, bitFlip, "Shots" -> 4000, "Range" -> {0.01, 0.04}, "Points" -> 4]},
            {th["Value"], AnyTrue[th["Caveats"], StringContainsQ[#, "above the range"] &]}
        ]
    ],
    {Missing["AboveRange"], True},
    TestID -> "QEC-Threshold-reports-a-threshold-above-the-range"
]

(* Building the model once with p symbolic and substituting gives the very samples
   that building it at the point gives. *)
VerificationTest[
    With[{r5 = QECCode["Repetition", 5], noise = QECNoiseModel["BitFlip", 1/20, "MeasurementError" -> 1/20]},
        {BlockRandom[SeedRandom[5]; QECLogicalErrorRate[r5, noise, 3000, "Rounds" -> 5, "Decoder" -> "Matching"]],
         BlockRandom[SeedRandom[5];
            qecDemSampledReport[qecCodeDetectorModel[First[r5], First[phenomenological], 5, "BareAncilla"] /. qecP -> 1/20,
                3000, "Matching", Automatic]["Rate"]]}
    ] // Apply[Equal],
    True,
    TestID -> "QEC-Threshold-symbolic-model-substituted-equals-model-built-at-the-point"
]


(* ============================================================================
   The object
   ============================================================================ *)

VerificationTest[
    {exactRepetition["Parameter"], exactRepetition["Level"], Head[exactRepetition["Noise"]],
     Head /@ exactRepetition["Codes"], Length[exactRepetition["Caveats"]] >= 3, Head[exactRepetition["Plot"]]},
    {qecP, "CodeCapacity", QECNoiseModel, {QECCode, QECCode, QECCode}, True, Legended},
    TestID -> "QEC-Threshold-carries-its-context"
]

VerificationTest[
    StringContainsQ[First[exactRepetition["Caveats"]], "not a fault-tolerance threshold"],
    True,
    TestID -> "QEC-Threshold-code-capacity-is-not-fault-tolerance"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

VerificationTest[
    QECThreshold[{"Repetition", {3}}, bitFlip],
    $Failed,
    {QECThreshold::codes},
    TestID -> "QEC-Threshold-needs-two-codes"
]

VerificationTest[
    QECThreshold[{"Repetition", {3, 5}}, QECNoiseModel["BitFlip", qecP, "MeasurementError" -> qecQ]],
    $Failed,
    {QECThreshold::params},
    TestID -> "QEC-Threshold-needs-one-symbolic-rate"
]

VerificationTest[
    QECThreshold[{"Steane", {3, 5}}, bitFlip],
    $Failed,
    {QECThreshold::family},
    TestID -> "QEC-Threshold-unknown-family"
]

VerificationTest[
    QECThreshold[{"Surface", {3, 5}}, bitFlip, "Shots" -> 1000, "Range" -> {0.2, 0.1}],
    $Failed,
    {QECThreshold::range},
    TestID -> "QEC-Threshold-refuses-a-bad-range"
]
