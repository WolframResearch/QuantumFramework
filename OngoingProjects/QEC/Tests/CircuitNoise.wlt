(* ::Package:: *)

(* ============================================================================
   Tests/QEC/CircuitNoise.wlt

   The two noisy levels -- phenomenological and circuit -- the detector error
   model they produce, the memory experiment run on top of it, and the Stim
   export.

   Three oracles carry this file.

   The first is the code-capacity route, already tested against the binomial tail
   in Noise.wlt.  A phenomenological model with perfect readout, run for one
   round, is code-capacity noise wearing a circuit: its logical error rate must
   come out as the same polynomial.  That single identity exercises the circuit
   emitter, the frame propagator, the detector definition, the fault enumeration
   and the exact sum at once.

   The second is arithmetic done by hand on the bit-flip code, which is small
   enough that every row of its detector model can be written down.

   The third is Stim.  Its numbers are not reproduced here -- the suite must run
   without a Python install -- but the exported circuits were checked against it:
   58 detector firing rates and 5 observable flip rates across five cases,
   computed exactly here and sampled there, all agreeing within 2.7 sigma of two
   million shots, and PyMatching decodes the exported circuits unchanged.  The
   harness is StimCrossCheck/ (emit.wls, check.py, README.md).
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];


(* ============================================================================
   Levels of the noise model
   ============================================================================ *)

VerificationTest[
    QECNoiseModel["Depolarizing", qecP]["Level"],
    "CodeCapacity",
    TestID -> "QEC-CircuitNoise-level-default"
]

VerificationTest[
    QECNoiseModel["Depolarizing", qecP, "MeasurementError" -> qecQ]["Level"],
    "Phenomenological",
    TestID -> "QEC-CircuitNoise-level-phenomenological"
]

(* Readout error zero is code capacity, not a degenerate phenomenological model:
   the two are the same object and should report the same level. *)
VerificationTest[
    QECNoiseModel["BitFlip", qecP, "MeasurementError" -> 0]["Level"],
    "CodeCapacity",
    TestID -> "QEC-CircuitNoise-level-zero-readout"
]

VerificationTest[
    QECNoiseModel["Circuit", qecP]["Level"],
    "Circuit",
    TestID -> "QEC-CircuitNoise-level-circuit"
]

VerificationTest[
    QECNoiseModel["Depolarizing", qecP, "MeasurementError" -> qecQ]["MeasurementError"],
    qecQ,
    TestID -> "QEC-CircuitNoise-readout-rate"
]

(* A single number spreads over every location but idling, which has no parallel
   layer to sit in yet. *)
VerificationTest[
    QECNoiseModel["Circuit", qecP]["Rates"],
    <|"OneQubit" -> qecP, "TwoQubit" -> qecP, "Measurement" -> qecP, "Reset" -> qecP, "Idle" -> 0|>,
    TestID -> "QEC-CircuitNoise-uniform-rates"
]

VerificationTest[
    QECNoiseModel["Circuit", <|"Measurement" -> qecQ|>]["Rates"],
    <|"OneQubit" -> 0, "TwoQubit" -> 0, "Measurement" -> qecQ, "Reset" -> 0, "Idle" -> 0|>,
    TestID -> "QEC-CircuitNoise-partial-rates"
]

VerificationTest[
    QECNoiseModel["Circuit", <|"Typo" -> 1/10|>],
    $Failed,
    {QECNoiseModel::circuitrates},
    TestID -> "QEC-CircuitNoise-bad-rate-key"
]

VerificationTest[
    QECNoiseModel["Circuit", qecP]["SymbolicQ"],
    True,
    TestID -> "QEC-CircuitNoise-symbolic-circuit"
]

VerificationTest[
    QECNoiseModel["Circuit", 1/1000]["SymbolicQ"],
    False,
    TestID -> "QEC-CircuitNoise-numeric-circuit"
]


(* ============================================================================
   The detector model, checked row by row on the bit-flip code

   Generators ZZI and IZZ, so with one round there are four detectors: the two
   checks of round 1, then the two differences against the noiseless final
   readout.
   ============================================================================ *)

qecBF = QECCode["BitFlipCode"];

qecDem1 = QECDetectorModel[qecBF, QECNoiseModel["BitFlip", qecP, "MeasurementError" -> qecQ], 1];

VerificationTest[
    QECDetectorModel[qecBF, QECNoiseModel["BitFlip", qecP], 1],
    $Failed,
    {QECDetectorModel::level},
    TestID -> "QEC-CircuitNoise-dem-refuses-code-capacity"
]

VerificationTest[
    qecDem1["Detectors"],
    4,
    TestID -> "QEC-CircuitNoise-dem-detector-count"
]

(* Bit-flip noise makes no Y and no Z, so those rows carry probability exactly
   zero and are dropped: three data faults and two readout faults remain. *)
VerificationTest[
    qecDem1["Faults"],
    5,
    TestID -> "QEC-CircuitNoise-dem-drops-impossible-faults"
]

VerificationTest[
    qecDem1["LocationCounts"],
    <|"Data" -> 3, "Measurement" -> 2|>,
    TestID -> "QEC-CircuitNoise-dem-location-counts"
]

(* X on qubit 1 anticommutes with ZZI only, and a data error present in round 1 is
   still there at the final readout, so the difference detectors stay quiet. *)
VerificationTest[
    qecDem1["DetectorMatrix"][[1]],
    {1, 0, 0, 0},
    TestID -> "QEC-CircuitNoise-dem-data-x1"
]

VerificationTest[
    qecDem1["DetectorMatrix"][[2]],
    {1, 1, 0, 0},
    TestID -> "QEC-CircuitNoise-dem-data-x2"
]

(* The whole point of detectors: a readout that lies fires the check in its own
   round and again in the difference against the next one, so it looks nothing
   like a data error. *)
VerificationTest[
    qecDem1["DetectorMatrix"][[4]],
    {1, 0, 1, 0},
    TestID -> "QEC-CircuitNoise-dem-readout-fires-twice"
]

VerificationTest[
    qecDem1["DetectorMatrix"][[5]],
    {0, 1, 0, 1},
    TestID -> "QEC-CircuitNoise-dem-readout-fires-twice-second-check"
]

(* A readout error leaves nothing behind on the data. *)
VerificationTest[
    qecDem1["ObservableMatrix"][[{4, 5}]],
    {{0, 0}, {0, 0}},
    TestID -> "QEC-CircuitNoise-dem-readout-no-residual"
]

(* X on qubit 3 anticommutes with the logical Z, so it is the one data fault that
   moves the logical qubit. *)
VerificationTest[
    qecDem1["ObservableMatrix"][[3]],
    {0, 1},
    TestID -> "QEC-CircuitNoise-dem-logical-class"
]

VerificationTest[
    qecDem1["Probabilities"],
    {qecP, qecP, qecP, qecQ, qecQ},
    TestID -> "QEC-CircuitNoise-dem-probabilities"
]

(* With depolarizing noise the Z errors are possible, and the bit-flip code cannot
   see them: no detector fires, yet the logical qubit is damaged.  That is the
   circuit-level statement of "this code has no phase protection". *)
VerificationTest[
    Length @ QECDetectorModel[qecBF,
        QECNoiseModel["Depolarizing", qecP, "MeasurementError" -> qecQ], 1]["UndetectableFaults"],
    3,
    TestID -> "QEC-CircuitNoise-dem-undetectable-z"
]

(* Exact detector firing rates.  With no data noise every detector fires exactly
   when its own readout lies. *)
VerificationTest[
    Simplify @ QECDetectorModel[qecBF,
        QECNoiseModel["BitFlip", 0, "MeasurementError" -> qecQ], 1]["DetectorRates"],
    {qecQ, qecQ, qecQ, qecQ},
    TestID -> "QEC-CircuitNoise-dem-detector-rates"
]


(* ============================================================================
   The memory experiment

   The identity that ties the new machinery to the old: perfect readout for one
   round is code-capacity noise, and must give the same polynomial.
   ============================================================================ *)

VerificationTest[
    Simplify[
        QECLogicalErrorRate[qecBF, QECNoiseModel["BitFlip", qecP, "MeasurementError" -> 0], "Rounds" -> 1] ==
        QECLogicalErrorRate[qecBF, QECNoiseModel["BitFlip", qecP]]
    ],
    True,
    TestID -> "QEC-CircuitNoise-perfect-readout-is-code-capacity-bitflip"
]

VerificationTest[
    Simplify[
        QECLogicalErrorRate[qecBF,
            QECNoiseModel["BitFlip", qecP, "MeasurementError" -> 0], "Rounds" -> 1] ==
        3 qecP^2 - 2 qecP^3
    ],
    True,
    TestID -> "QEC-CircuitNoise-bitflip-polynomial"
]

(* The same identity on a code with all three Paulis and a real distance.  Both
   sides decode by minimum weight, which for the perfect 5-qubit code under
   symmetric noise is also the maximum-likelihood answer. *)
VerificationTest[
    Simplify[
        QECLogicalErrorRate[QECCode["5QubitCode"],
            QECNoiseModel["Depolarizing", qecP, "MeasurementError" -> 0], "Rounds" -> 1] ==
        QECLogicalErrorRate[QECCode["5QubitCode"],
            QECNoiseModel["Depolarizing", qecP], "Decoder" -> "MinimumWeight"]
    ],
    True,
    TestID -> "QEC-CircuitNoise-perfect-readout-is-code-capacity-5qubit"
]

(* No data noise: a single readout lie is recognised by its two-detector
   signature and corrected, so failure needs two of them. *)
VerificationTest[
    Simplify @ QECLogicalErrorRate[qecBF,
        QECNoiseModel["BitFlip", 0, "MeasurementError" -> qecQ], "Rounds" -> 1],
    qecQ^2,
    TestID -> "QEC-CircuitNoise-readout-only-rate"
]

(* Sampling and the exact sum are two routes to the same number. *)
VerificationTest[
    With[
        {exact = N @ QECLogicalErrorRate[qecBF, QECNoiseModel["Circuit", 1/200], "Rounds" -> 2],
         sampled = QECLogicalErrorRate[qecBF, QECNoiseModel["Circuit", 1/200], 200000, "Rounds" -> 2]},
        Abs[exact - sampled] < 5 Sqrt[exact (1 - exact) / 200000]
    ],
    True,
    TestID -> "QEC-CircuitNoise-exact-matches-sampled"
]

VerificationTest[
    QECLogicalErrorRate[qecBF, QECNoiseModel["Circuit", qecP], 1000, "Rounds" -> 1],
    $Failed,
    {QECLogicalErrorRate::symbolicshots},
    TestID -> "QEC-CircuitNoise-symbolic-shots-refused"
]

(* The state space of the exact route is 2^(detectors + observables), and it says
   so rather than trying. *)
VerificationTest[
    QECLogicalErrorRate[QECCode["SteaneCode"],
        QECNoiseModel["Depolarizing", 1/100, "MeasurementError" -> 1/100], "Rounds" -> 3],
    $Failed,
    {QECLogicalErrorRate::detbig},
    TestID -> "QEC-CircuitNoise-exact-limit"
]


(* ============================================================================
   The Stim export
   ============================================================================ *)

qecStim = QECStimCircuit[qecBF, QECNoiseModel["Circuit", 1/200], 2];

VerificationTest[
    StringQ[qecStim],
    True,
    TestID -> "QEC-CircuitNoise-stim-string"
]

(* Stim starts in |0...0>, which is not a codeword, so the encoder has to go first
   or the first round of detectors means nothing. *)
VerificationTest[
    StringContainsQ[qecStim, "# --- encode ---"],
    True,
    TestID -> "QEC-CircuitNoise-stim-encoder"
]

(* One detector per check per round, plus the final noiseless round. *)
VerificationTest[
    StringCount[qecStim, "DETECTOR"],
    6,
    TestID -> "QEC-CircuitNoise-stim-detector-count"
]

VerificationTest[
    StringCount[qecStim, "OBSERVABLE_INCLUDE"],
    1,
    TestID -> "QEC-CircuitNoise-stim-one-observable"
]

VerificationTest[
    {StringContainsQ[qecStim, "DEPOLARIZE2(0.005)"], StringContainsQ[qecStim, "X_ERROR(0.005)"],
     StringContainsQ[qecStim, "M(0.005)"]},
    {True, True, True},
    TestID -> "QEC-CircuitNoise-stim-noise-lines"
]

(* The final round carries no noise: it is the perfect readout that closes the
   experiment. *)
VerificationTest[
    StringFreeQ[StringDrop[qecStim, StringPosition[qecStim, "final noiseless round"][[1, 1]]], "DEPOLARIZE"],
    True,
    TestID -> "QEC-CircuitNoise-stim-final-round-clean"
]

VerificationTest[
    StringContainsQ[
        QECStimCircuit[qecBF, QECNoiseModel["Depolarizing", 1/100, "MeasurementError" -> 1/50], 1],
        "PAULI_CHANNEL_1"],
    True,
    TestID -> "QEC-CircuitNoise-stim-phenomenological-channel"
]

VerificationTest[
    QECStimCircuit[qecBF, QECNoiseModel["Depolarizing", 1/100], 1],
    $Failed,
    {QECStimCircuit::level},
    TestID -> "QEC-CircuitNoise-stim-refuses-code-capacity"
]

VerificationTest[
    QECStimCircuit[qecBF, QECNoiseModel["Circuit", 1/200], 1, "Observable" -> "Q"],
    $Failed,
    {QECStimCircuit::observable},
    TestID -> "QEC-CircuitNoise-stim-bad-observable"
]

(* The Y-generator path reaches Stim as its sqrt(X). *)
VerificationTest[
    StringContainsQ[QECStimCircuit[QECCode[{"YYI", "IYY"}], QECNoiseModel["Circuit", 1/500], 1], "SQRT_X"],
    True,
    TestID -> "QEC-CircuitNoise-stim-sqrt-x"
]
