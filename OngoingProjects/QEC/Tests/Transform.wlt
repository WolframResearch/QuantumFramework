(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Transform.wlt

   Exact rates through the Walsh-Hadamard transform (QECCore/Transform.wl), the
   maximum-likelihood decoder it gives for a detector model, and a memory kept to
   one basis ("Basis" -> "X" | "Z").

   What each group is checked against.
   - The transform against the fold it replaced (demExactFailureFold, kept as the
     reference): the same polynomial, symbol for symbol, on circuit-level and
     phenomenological models, on a model with heralds, and on exact numeric rates.
   - surface-code.nb (November 2026), which computed the exact rate of
     the d = 3 surface code under one noisy round of CNOT depolarizing (pC/15 for
     each of the 15 Paulis) and readout errors (pM), with a decoder that knows
     every single fault and one that also knows every pair.  The layer's
     detector model is his fault list exactly, its default decoder is his
     single-fault one (same polynomial, pinned here by its leading terms and its
     value at 1/1000), and maximum likelihood has the leading terms of his pair
     decoder and a lower rate.
   - Maximum likelihood is a lower bound on every other decoder of the same model.
   - The one-basis memory: its counts, its refusals, and that it decodes.

   Everything here runs on bare Wolfram.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";
qecFold = Symbol[qecScope <> "demExactFailureFold"];
qecExact = Symbol[qecScope <> "demExactFailure"];

(* the transform and the fold on the same model and the same decoder table *)
sameAsFold[dem_QECDetectorModel] := With[{a = First[dem]},
    With[{new = qecExact[a, 1], old = qecFold[a, 1]},
        If[ AssociationQ[new],
            Simplify[{new["Failure"] - old["Failure"], new["Acceptance"] - old["Acceptance"]}] === {0, 0},
            Simplify[new - old] === 0
        ]
    ]
]

leading[expr_, order_] := Expand[Normal[Series[expr /. {pC -> t pC, pM -> t pM}, {t, 0, order}]] /. t -> 1]

surface3 = QECCode["Surface"[3]];
nbNoise = QECNoiseModel["Circuit"[<|"TwoQubit" -> pC, "Measurement" -> pM|>]];


(* ============================================================================
   The transform is the fold
   ============================================================================ *)

VerificationTest[
    sameAsFold[QECDetectorModel[QECCode["Repetition"[3]], QECNoiseModel["Circuit"[qecP]], 2]],
    True,
    TestID -> "QEC-Transform-equals-fold-repetition-circuit"
]

VerificationTest[
    sameAsFold[QECDetectorModel[QECCode["SteaneCode"], QECNoiseModel["Depolarizing"[qecP], "MeasurementError" -> qecQ], 1]],
    True,
    TestID -> "QEC-Transform-equals-fold-steane-phenomenological"
]

VerificationTest[
    sameAsFold[QECDetectorModel[surface3, QECNoiseModel["Circuit"[qecP]], 1, "Basis" -> "Z"]],
    True,
    TestID -> "QEC-Transform-equals-fold-surface-one-basis"
]

VerificationTest[
    sameAsFold[QECDetectorModel[QECCode["Repetition"[3]], QECNoiseModel["Circuit"[1/100]], 2]],
    True,
    TestID -> "QEC-Transform-equals-fold-exact-numeric-rates"
]

(* Heralds: the cat checks of the transversal extraction post-select, and both the
   joint failure and the acceptance must come out the same. *)
VerificationTest[
    sameAsFold[QECDetectorModel[QECCode["ShorCode"], QECNoiseModel["Circuit"[<|"TwoQubit" -> qecP, "Measurement" -> qecP|>]], 1,
        "Extraction" -> "Transversal", "Basis" -> "X"]],
    True,
    TestID -> "QEC-Transform-equals-fold-with-heralds"
]

VerificationTest[
    $QECExactDetectorLimit,
    2^24,
    TestID -> "QEC-Transform-exact-limit"
]


(* ============================================================================
   The surface-code.nb study
   ============================================================================ *)

(* 352 faults that fire something: 24 CNOTs x 15 Paulis and 8 readouts, less the 16
   that do nothing; 16 detectors and both logical observables. *)
VerificationTest[
    With[{a = First[QECDetectorModel[surface3, nbNoise, 1]]},
        {Length[a["Probabilities"]], a["Detectors"], a["Observables"], Union[a["Probabilities"]]}
    ],
    {352, 16, 2, Union[{pC/15, pM}]},
    TestID -> "QEC-Transform-surface-nb-detector-model"
]

(* His single-fault decoder: no linear term, and this p^2 term. *)
VerificationTest[
    leading[QECLogicalErrorRate[surface3, nbNoise, "Rounds" -> 1], 2],
    (48326 pC^2)/225 + (2213 pC pM)/15 + 28 pM^2,
    TestID -> "QEC-Transform-surface-nb-single-fault-decoder"
]

VerificationTest[
    QECLogicalErrorRate[surface3, nbNoise, "Rounds" -> 1] /. {pC -> 1/1000, pM -> 1/1000},
    0.00038328649725612645,
    SameTest -> (Abs[#1 - #2] < 10^-15 &),
    TestID -> "QEC-Transform-surface-nb-single-fault-value"
]

(* Maximum likelihood: the leading terms of his pair decoder, 2438/45 pC^2 + pC pM/5,
   and a lower rate than his 5.5497e-5 at 1/1000. *)
VerificationTest[
    With[{r = QECLogicalErrorRate[surface3, nbNoise, "Rounds" -> 1, "Decoder" -> {"MaximumLikelihood", 1/1000}]},
        {leading[r, 2], r /. {pC -> 1/1000, pM -> 1/1000}}
    ],
    {(2438 pC^2)/45 + (pC pM)/5, 0.0000541174729859932},
    SameTest -> (#1[[1]] === #2[[1]] && Abs[#1[[2]] - #2[[2]]] < 10^-15 &),
    TestID -> "QEC-Transform-surface-nb-maximum-likelihood"
]


(* ============================================================================
   Maximum likelihood
   ============================================================================ *)

(* No decoder of the same model beats it, at the rate its table is built for. *)
VerificationTest[
    With[{n = QECNoiseModel["Circuit"[1/200]]},
        Module[{ml, mw, mt},
            ml = QECLogicalErrorRate[surface3, n, "Rounds" -> 2, "Basis" -> "Z", "Decoder" -> "MaximumLikelihood"];
            mw = QECLogicalErrorRate[surface3, n, "Rounds" -> 2, "Basis" -> "Z"];
            mt = QECLogicalErrorRate[surface3, n, "Rounds" -> 2, "Basis" -> "Z", "Decoder" -> "Matching"];
            ml <= mw && ml <= mt
        ]
    ],
    True,
    TestID -> "QEC-Transform-maximum-likelihood-is-a-lower-bound"
]

(* The same table decodes sampled shots. *)
VerificationTest[
    With[{n = QECNoiseModel["Circuit"[1/50]]},
        Abs[QECLogicalErrorRate[surface3, n, 20000, "Rounds" -> 1, "Basis" -> "Z", "Decoder" -> "MaximumLikelihood"] -
            QECLogicalErrorRate[surface3, n, "Rounds" -> 1, "Basis" -> "Z", "Decoder" -> "MaximumLikelihood"]] < 0.01
    ],
    True,
    TestID -> "QEC-Transform-maximum-likelihood-sampled"
]

VerificationTest[
    FailureQ[QECLogicalErrorRate[surface3, nbNoise, "Rounds" -> 1, "Decoder" -> "Fastest"]],
    True,
    {QECLogicalErrorRate::decoder},
    TestID -> "QEC-Transform-unknown-decoder-refused"
]


(* ============================================================================
   A memory in one basis
   ============================================================================ *)

VerificationTest[
    With[{m = QECDetectorModel[surface3, QECNoiseModel["Circuit"[qecP]], 3, "Basis" -> "Z"]},
        {m["Basis"], m["Detectors"], m["Observables"]}
    ],
    {"Z", 16, 1},
    TestID -> "QEC-Transform-one-basis-counts"
]

(* Both one-basis models of the d = 3 surface code have distance 3: no linear term. *)
VerificationTest[
    Table[
        Coefficient[QECLogicalErrorRate[surface3, QECNoiseModel["Circuit"[qecP]], "Rounds" -> 2, "Basis" -> b], qecP, 1],
        {b, {"X", "Z"}}
    ],
    {0, 0},
    TestID -> "QEC-Transform-one-basis-distance-three"
]

VerificationTest[
    FailureQ[QECDetectorModel[surface3, QECNoiseModel["Circuit"[qecP]], 1, "Basis" -> "Y"]],
    True,
    {QECDetectorModel::basis},
    TestID -> "QEC-Transform-one-basis-unknown-refused"
]

VerificationTest[
    FailureQ[QECDetectorModel[QECCode[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}], QECNoiseModel["Circuit"[qecP]], 1, "Basis" -> "Z"]],
    True,
    {QECDetectorModel::basiscss},
    TestID -> "QEC-Transform-one-basis-needs-css"
]

VerificationTest[
    QECDetectorModel[surface3, QECNoiseModel["Circuit"[qecP]], 1, "Basis" -> "Z"]["Stim"],
    Missing["NotAvailable", "Z"],
    {QECDetectorModel::stimbasis},
    TestID -> "QEC-Transform-one-basis-no-stim"
]
