(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Herald.wlt

   Post-selection in the rate machinery.

   A herald ("MH") is a verification outcome that discards the shot rather than
   contributing a syndrome bit.  Once a circuit has one, a logical error rate is
   CONDITIONAL on acceptance, and a conditional rate reported as a bare number is
   meaningless -- so the two routes return an association carrying their own
   acceptance instead, and return a plain number only when there is nothing to
   condition on.

   The synthetic detector model below is deliberately hand-built: it isolates the
   conditioning algebra from the gadget that will produce heralds for real (the cat
   state measurement of sec. 12.1, phase B).  The row filter's herald clause in
   codeDetectorModel gets its first end-to-end exercise there.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecExact  = Symbol[qecScope <> "demExactFailure"];
qecSamp   = Symbol[qecScope <> "demSampledFailure"];
qecKeys   = Symbol[qecScope <> "demRowKeys"];
qecEffect = Symbol[qecScope <> "faultEffect"];
qecCodeData = Symbol[qecScope <> "codeData"];
qecCircInstr = Symbol[qecScope <> "codeCircuitInstructions"];

(* One detector, one observable, one herald, three independent mechanisms:
     A  fires the detector and nothing else   -- detectable and correctable
     B  flips the observable and nothing else -- an undetectable logical fault
     C  trips the herald and nothing else     -- a rejected preparation           *)
toy = <|
    "DetectorMatrix"   -> {{1}, {0}, {0}},
    "ObservableMatrix" -> {{0}, {1}, {0}},
    "HeraldMatrix"     -> {{0}, {0}, {1}},
    "Probabilities"    -> {pd, pl, ph},
    "Locations"        -> {{"A", 1}, {"B", 2}, {"C", 3}},
    "Detectors" -> 1, "Observables" -> 1, "Heralds" -> 1,
    "Rounds" -> 1, "Checks" -> 1
|>;

toyN = <|toy, "Probabilities" -> {0.05, 0.02, 0.30}|>;


(* ============================================================================
   Nothing changes when there are no heralds
   ============================================================================ *)

(* THE REGRESSION.  Both routes still return a bare number, and the numbers are the
   ones they always were, so the rest of the suite is untouched by this file. *)
VerificationTest[
    Expand[QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", p]]],
    3 p^2 - 2 p^3,
    TestID -> "QEC-Herald-code-capacity-unchanged"
]

VerificationTest[
    Normal @ Series[
        QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["Circuit", p], "Rounds" -> 1],
        {p, 0, 1}],
    32 p / 15,
    TestID -> "QEC-Herald-circuit-level-unchanged"
]

VerificationTest[
    Head[QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["Circuit", 0.01],
        20000, "Rounds" -> 1]],
    Real,
    TestID -> "QEC-Herald-no-heralds-gives-a-bare-number"
]

VerificationTest[
    With[{dem = QECDetectorModel[QECCode["BitFlipCode"], QECNoiseModel["Circuit", p], 1]},
        {dem["Heralds"], dem["PostSelectedQ"], dem["HeraldRates"]}
    ],
    {0, False, {}},
    TestID -> "QEC-Herald-extraction-circuit-has-none"
]

(* faultEffect grew a third component, and it is empty for a circuit with no heralds. *)
VerificationTest[
    With[{a = qecCodeData[QECCode["BitFlipCode"]]},
        With[{eff = qecEffect[a, qecCircInstr[a, 1], 5, 2, 1, {{0, 1, {1, 0}}}]},
            {Length[eff], Last[eff]}
        ]
    ],
    {3, {}},
    TestID -> "QEC-Herald-faultEffect-has-three-components"
]


(* ============================================================================
   The conditioning algebra
   ============================================================================ *)

(* The herald column reaches the packed keys as a third list. *)
VerificationTest[
    qecKeys[toy],
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
    TestID -> "QEC-Herald-keys-carry-three-columns"
]

(* THE ALGEBRA, exactly and symbolically.  Mechanism C rejects with probability ph, so
   acceptance is 1 - ph.  The joint failure is the logical fault AND acceptance.  And
   the conditional rate is just pl: post-selection here is independent of the logical
   fault, so conditioning divides it back out. *)
VerificationTest[
    With[{ex = qecExact[toy, 1]},
        Simplify[{ex["Rate"], ex["Acceptance"], ex["Failure"]} == {pl, 1 - ph, pl (1 - ph)}]
    ],
    True,
    TestID -> "QEC-Herald-exact-conditioning-is-correct"
]

(* And it comes back as an association, not a number: a conditional rate must carry
   its acceptance or it cannot be read. *)
VerificationTest[
    Sort[Keys[qecExact[toy, 1]]],
    {"Acceptance", "Failure", "Rate"},
    TestID -> "QEC-Herald-exact-returns-rate-with-acceptance"
]

VerificationTest[
    Sort[Keys[qecSamp[toyN, 2000, 1]]],
    {"Acceptance", "Accepted", "Rate", "Shots"},
    TestID -> "QEC-Herald-sampled-returns-rate-with-acceptance"
]

(* The two routes agree.  Fixed seed so the tolerance is a statement about the
   estimator and not about luck. *)
VerificationTest[
    Block[{}, SeedRandom[20260909];
        With[{s = qecSamp[toyN, 200000, 1], e = qecExact[toyN, 1]},
            {Abs[s["Rate"] - e["Rate"]] < 0.003,
             Abs[s["Acceptance"] - e["Acceptance"]] < 0.01}
        ]
    ],
    {True, True},
    TestID -> "QEC-Herald-sampled-agrees-with-exact"
]

(* Acceptance is the probability the heralds all read zero, and nothing else. *)
VerificationTest[
    Simplify[qecExact[<|toy, "Probabilities" -> {pd, pl, 0}|>, 1]["Acceptance"]],
    1,
    TestID -> "QEC-Herald-no-rejection-means-full-acceptance"
]

(* A herald that fires with certainty rejects everything, and the exact route says so
   rather than dividing by zero. *)
VerificationTest[
    qecExact[<|toy, "Probabilities" -> {0, 0, 1}|>, 1],
    <|"Rate" -> Indeterminate, "Acceptance" -> 0, "Failure" -> 0|>,
    {QECLogicalErrorRate::noacceptance},
    TestID -> "QEC-Herald-certain-rejection-gives-zero-acceptance"
]

(* The sampler refuses honestly when every shot was thrown away, instead of returning
   a rate computed from nothing. *)
VerificationTest[
    qecSamp[<|toy, "Probabilities" -> {0., 0., 1.}|>, 1000, 1],
    <|"Rate" -> Indeterminate, "Acceptance" -> 0., "Accepted" -> 0, "Shots" -> 1000|>,
    {QECLogicalErrorRate::allrejected},
    TestID -> "QEC-Herald-all-rejected-is-reported"
]

(* Heralds join the exact route's state space, since conditioning needs the joint
   distribution -- so a post-selected model hits the limit sooner, and the message
   counts them. *)
VerificationTest[
    Block[{$QECExactDetectorLimit = 4},
        qecExact[toy, 1]
    ],
    $Failed,
    {QECLogicalErrorRate::detbig},
    TestID -> "QEC-Herald-heralds-count-toward-the-exact-limit"
]

(* The decoder never sees heralds.  It reads detectors only, so its table is the same
   with and without the herald column -- which is the point: post-selection is not
   information the decoder gets to use. *)
VerificationTest[
    With[{withH = Symbol[qecScope <> "demDecoderTable"][toy, 1],
          withoutH = Symbol[qecScope <> "demDecoderTable"][
              <|toy, "HeraldMatrix" -> {{}, {}, {}}, "Heralds" -> 0|>, 1]},
        Normal[withH] === Normal[withoutH]
    ],
    True,
    TestID -> "QEC-Herald-decoder-does-not-see-heralds"
]
