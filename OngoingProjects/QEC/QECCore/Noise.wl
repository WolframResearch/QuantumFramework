(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECNoiseModel]

PackageScope[noiseProbabilities]
PackageScope[noiseSymbolicQ]
PackageScope[noiseErrorProbability]
PackageScope[noiseProfileProbability]
PackageScope[noiseSampleRows]
PackageScope[noiseLevel]
PackageScope[noiseMeasurementError]
PackageScope[noiseCircuitRates]
PackageScope[$QECNoiseLevels]
PackageScope[$QECCircuitLocations]


(* ============================================================================ *)
(* The noise model.                                                             *)
(*                                                                              *)
(* A code means nothing until you say what it is protecting against.  This is    *)
(* the simplest honest answer, the one the roadmap calls code-capacity noise:    *)
(* each qubit independently suffers I, X, Y or Z with fixed probabilities,       *)
(* between rounds, with perfect check measurements.  Circuit-level noise (noisy  *)
(* syndrome extraction) is the next object, not this one.                        *)
(*                                                                              *)
(* Two things make this more than a random-number generator.                     *)
(*                                                                              *)
(* First, the probabilities may be **symbolic**.  Everything downstream then      *)
(* returns an exact expression in p rather than an estimate: a logical error     *)
(* rate is a polynomial, and a threshold is where two polynomials cross.  That   *)
(* is the one move a sampling simulator structurally cannot make, and it is why  *)
(* the model is built around exact enumeration first and sampling second.        *)
(*                                                                              *)
(* Second, because the errors are independent and identically distributed, the   *)
(* probability of a specific error depends only on how many X, Y and Z it        *)
(* carries -- not on where they sit:                                             *)
(*                                                                              *)
(*     P(E) = pI^(n - nx - ny - nz) pX^nx pY^ny pZ^nz                            *)
(*                                                                              *)
(* so an exact sum over 4^n errors collapses to a sum over the far smaller set   *)
(* of weight profiles.  Both the decoder and the error rate lean on this.        *)
(* ============================================================================ *)

QECNoiseModel::usage = "QECNoiseModel[name, p] represents independent identically distributed Pauli noise on each qubit: \"Depolarizing\", \"BitFlip\", \"PhaseFlip\" or \"BitPhaseFlip\" at rate p, with perfect checks (code-capacity level).\nQECNoiseModel[<|\"X\" -> px, \"Y\" -> py, \"Z\" -> pz|>] gives a general Pauli channel.\nQECNoiseModel[name, p, \"MeasurementError\" -> q] adds check outcomes that are misreported with probability q (phenomenological level).\nQECNoiseModel[\"Circuit\", p] puts depolarizing noise on every gate, measurement and reset of the syndrome-extraction circuit at rate p; QECNoiseModel[\"Circuit\", <|\"OneQubit\" -> p1, ...|>] sets the locations separately (circuit level).\nThe rates may be symbolic, in which case downstream quantities come back as exact expressions in them.\nnoise[prop] gives a property; noise[\"Properties\"] lists them.";

QECNoiseModel::rates = "Error probabilities must be non-negative and sum to at most 1; got `1`.";
QECNoiseModel::unknown = "`1` is not a known noise model. Available: `2`.";
QECNoiseModel::noprop = "`1` is not a property of QECNoiseModel. Use noise[\"Properties\"] for the list.";

QECNoiseModel::level = "`1` is not a noise level. Available: `2`.";
QECNoiseModel::circuitrates = "Circuit-level rates must be given as a number or an Association over `1`; got `2`.";
QECNoiseModel::wronglevel = "This is a `1` noise model; `2` is only defined at the `3` level.";

$QECNoiseNames = {"Depolarizing", "BitFlip", "PhaseFlip", "BitPhaseFlip"};

(* The three levels, in the order the roadmap builds them.  Each drops an
   assumption the one before it made, but they are NOT three estimates of the same
   number -- see the warning at the end.

   CodeCapacity     errors on the data, checks read out perfectly.  Exact and
                    symbolic; this is where a logical error rate is a polynomial.
                    Not a term from Gottesman's book; it is what part I of the book
                    assumes throughout, "the errors only had one shot at the
                    quantum state".
   Phenomenological the same data errors, plus check outcomes that lie with
                    probability q.  One round is then not enough to tell a data
                    error from a misread check, which is why the experiment has to
                    be repeated -- this level is what makes rounds mean anything.
                    The book names it, and names it as the standard example of an
                    *optimistic* assumption (sec. 10.4): "an error rate per
                    physical qubit per time step (regardless of gates) and an error
                    rate on each bit of the error syndrome".
   Circuit          no more "checks": the extraction subcircuit is written out and
                    every gate, reset and readout in it can fail.  A Pauli
                    specialisation of the book's basic model (sec. 10.1.1,
                    Definitions 10.1-10.3): faults attach to circuit locations, one
                    rate per location type.  The level a hardware person would
                    accept, and the one where the number almost always has to be
                    sampled.

   The levels are not separate objects: a code-capacity model is a phenomenological
   one with q = 0, and both are consumed by the same detector machinery.

   The warning, and it is the book's, sec. 10.4: "A threshold in a phenomenological
   model should not be confused with the threshold derived from a full circuit
   model; they are not at all comparable."  The reason is that the p of a
   phenomenological model is a rate per qubit per round while the p of a circuit
   model is a rate per location, and the phenomenological model "completely ignores"
   how many gates a protocol needs to extract a syndrome bit.  Plotting the three
   against a shared p, as the tech note does, shows one code losing ground as
   assumptions are dropped; it does not compare thresholds, and must not be read
   as doing so. *)
$QECNoiseLevels = {"CodeCapacity", "Phenomenological", "Circuit"};

(* Where a fault can happen in the extraction circuit, and what it does.

   OneQubit    depolarizing after H, S, V, Vdg: X, Y or Z each with rate/3
   TwoQubit    depolarizing after CNOT: one of the 15 non-identity two-qubit
               Paulis, each with rate/15
   Measurement the readout lies.  Modelled as an X on the ancilla just before its
               measurement: the ancilla is reset immediately after, so this flips
               the record and nothing else, which is what a classical readout
               error is.
   Reset       the ancilla is prepared in the wrong state: an X just after "R"
   Idle        a qubit decohering while it waits.  Modelled: circuitSchedule groups
               the instructions into parallel time steps and a qubit no instruction
               of a step touches gets one location per step, three Paulis each.
               Still ZERO BY DEFAULT, and that default is optimistic rather than
               neutral -- serial extraction is the arrangement with the most idling.
               Set it and it is charged. *)
$QECCircuitLocations = {"OneQubit", "TwoQubit", "Measurement", "Reset", "Idle"};


(* ---- construction ---- *)

Options[QECNoiseModel] = {"MeasurementError" -> 0};

(* A numeric rate outside [0, 1] is a real error; a symbolic one cannot be checked
   and is taken on trust, which is the whole point of allowing it. *)
validRatesQ[{pI_, pX_, pY_, pZ_}] := ! AllTrue[{pI, pX, pY, pZ}, NumericQ] ||
    (AllTrue[{pI, pX, pY, pZ}, # >= 0 &] && pI + pX + pY + pZ == 1)

fromRates[rates_List, name_, params_] := fromRates[rates, name, params, 0]

fromRates[rates_List, name_, params_, measurementError_] := If[
    validRatesQ[rates],
    QECNoiseModel[<|
        "Probabilities" -> rates,
        "Name" -> name,
        "Parameters" -> params,
        "Level" -> If[TrueQ[measurementError === 0], "CodeCapacity", "Phenomenological"],
        "MeasurementError" -> measurementError
    |>],
    Message[QECNoiseModel::rates, rates]; $Failed
]

QECNoiseModel["Depolarizing", p_, opts : OptionsPattern[]] :=
    fromRates[{1 - p, p/3, p/3, p/3}, "Depolarizing", {p}, OptionValue[QECNoiseModel, {opts}, "MeasurementError"]]

QECNoiseModel["BitFlip", p_, opts : OptionsPattern[]] :=
    fromRates[{1 - p, p, 0, 0}, "BitFlip", {p}, OptionValue[QECNoiseModel, {opts}, "MeasurementError"]]

QECNoiseModel["PhaseFlip", p_, opts : OptionsPattern[]] :=
    fromRates[{1 - p, 0, 0, p}, "PhaseFlip", {p}, OptionValue[QECNoiseModel, {opts}, "MeasurementError"]]

QECNoiseModel["BitPhaseFlip", p_, opts : OptionsPattern[]] :=
    fromRates[{1 - p, 0, p, 0}, "BitPhaseFlip", {p}, OptionValue[QECNoiseModel, {opts}, "MeasurementError"]]

(* Circuit level.  A single number spreads the same rate over every location
   except idling, which stays at zero until the emitter produces parallel layers
   to idle in; an Association sets them one by one.  Missing keys are zero, so
   QECNoiseModel["Circuit", <|"Measurement" -> q|>] is a clean way to ask what
   readout noise alone does. *)
QECNoiseModel["Circuit", p_ ? (! AssociationQ[#] && (NumericQ[#] || ! FreeQ[#, _Symbol]) &)] := QECNoiseModel["Circuit",
    <|"OneQubit" -> p, "TwoQubit" -> p, "Measurement" -> p, "Reset" -> p, "Idle" -> 0|>
]

QECNoiseModel["Circuit", rates_Association] := If[
    ! ContainsOnly[Keys[rates], $QECCircuitLocations],
    Message[QECNoiseModel::circuitrates, $QECCircuitLocations, rates]; $Failed,
    QECNoiseModel[<|
        "Name" -> "Circuit",
        "Level" -> "Circuit",
        "Rates" -> AssociationMap[Lookup[rates, #, 0] &, $QECCircuitLocations],
        "Parameters" -> DeleteDuplicates[Cases[Values[rates], _Symbol, Infinity]]
    |>]
]

QECNoiseModel["Circuit", spec___] := (
    Message[QECNoiseModel::circuitrates, $QECCircuitLocations, {spec}];
    $Failed
)

QECNoiseModel[name_String, ___] := (
    Message[QECNoiseModel::unknown, name, StringRiffle[$QECNoiseNames, ", "]];
    $Failed
)

QECNoiseModel[rates_Association, opts : OptionsPattern[]] /; ContainsOnly[Keys[rates], {"I", "X", "Y", "Z"}] :=
    With[{px = Lookup[rates, "X", 0], py = Lookup[rates, "Y", 0], pz = Lookup[rates, "Z", 0]},
        fromRates[
            {Lookup[rates, "I", 1 - px - py - pz], px, py, pz},
            "Pauli",
            DeleteDuplicates[Cases[{px, py, pz}, _Symbol, Infinity]],
            OptionValue[QECNoiseModel, {opts}, "MeasurementError"]
        ]
    ]

QECNoiseModel[rates : {_, _, _, _}, opts : OptionsPattern[]] := fromRates[
    rates, "Pauli",
    DeleteDuplicates[Cases[rates, _Symbol, Infinity]],
    OptionValue[QECNoiseModel, {opts}, "MeasurementError"]
]


(* ---- accessors ---- *)

noiseProbabilities[QECNoiseModel[a_Association]] := a["Probabilities"]
noiseProbabilities[a_Association] := a["Probabilities"]

noiseSymbolicQ[a_Association] := ! AllTrue[
    If[a["Level"] === "Circuit", Values[a["Rates"]], a["Probabilities"]],
    NumericQ
]

noiseLevel[QECNoiseModel[a_Association]] := Lookup[a, "Level", "CodeCapacity"]
noiseLevel[a_Association] := Lookup[a, "Level", "CodeCapacity"]

noiseMeasurementError[a_Association] := Lookup[a, "MeasurementError", 0]

noiseCircuitRates[a_Association] := Lookup[a, "Rates", <||>]

(* P(E) from the multiplicities alone.

   The exponents are counts, so a factor with multiplicity zero is an empty product
   and contributes 1 -- including when its probability is also zero, which is exactly
   the case for the channels that switch a Pauli off (bit flip has P(Y) = P(Z) = 0).
   Writing this as a plain power would hand those models an Indeterminate from 0^0. *)
factor[_, 0] := 1
factor[q_, e_] := q^e

noiseProfileProbability[a_Association, n_Integer, {nx_, ny_, nz_}] := With[
    {p = a["Probabilities"]},
    factor[p[[1]], n - nx - ny - nz] factor[p[[2]], nx] factor[p[[3]], ny] factor[p[[4]], nz]
]

noiseErrorProbability[a_Association, pauli_] := With[{v = QECPauliVector[pauli]},
    With[{n = pauliQubits[v]},
        noiseProfileProbability[a, n, pauliProfile[v]]
    ]
]

(* {number of X, number of Y, number of Z} of a Pauli row. *)
pauliProfile[v_List] := With[{n = pauliQubits[v]},
    With[{pairs = Transpose[{v[[1 ;; n]], v[[n + 1 ;; 2 n]]}]},
        {Count[pairs, {1, 0}], Count[pairs, {1, 1}], Count[pairs, {0, 1}]}
    ]
]


(* ---- sampling ---- *)

(* count independent n-qubit errors, as Pauli rows.  Drawn per qubit in one
   RandomChoice call rather than qubit by qubit: the sampled route exists for the
   codes too big to enumerate, so it has to be cheap. *)
noiseSampleRows[a_Association, n_Integer, count_Integer] := Module[{letters, draws},
    letters = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    draws = RandomChoice[a["Probabilities"] -> letters, {count, n}];
    Join[draws[[All, All, 1]], draws[[All, All, 2]], ConstantArray[0, {count, 1}], 2]
]


(* ---- properties ---- *)

$noiseProperties = {
    "Probabilities", "Name", "Parameters", "Level", "MeasurementError", "Rates",
    "SymbolicQ", "MeanWeight", "ErrorProbability", "RandomError", "RandomErrors",
    "QuantumChannel", "Properties"
};

QECNoiseModel[_Association]["Properties"] := $noiseProperties

QECNoiseModel[a_Association]["Probabilities"] := Lookup[a, "Probabilities",
    Message[QECNoiseModel::wronglevel, "Circuit", "Probabilities", "code-capacity or phenomenological"];
    Missing["NotAvailable", "Probabilities"]
]
QECNoiseModel[a_Association]["Name"] := a["Name"]
QECNoiseModel[a_Association]["Parameters"] := a["Parameters"]
QECNoiseModel[a_Association]["SymbolicQ"] := noiseSymbolicQ[a]
QECNoiseModel[a_Association]["Level"] := noiseLevel[a]
QECNoiseModel[a_Association]["MeasurementError"] := noiseMeasurementError[a]
QECNoiseModel[a_Association]["Rates"] := noiseCircuitRates[a]

(* Expected number of qubits hit, per qubit. *)
QECNoiseModel[a_Association]["MeanWeight"] := Total[Rest[a["Probabilities"]]]

QECNoiseModel[a_Association]["ErrorProbability", pauli_] := noiseErrorProbability[a, pauli]

QECNoiseModel[a_Association]["RandomError", n_Integer] :=
    QECPauliString[First[noiseSampleRows[a, n, 1]]]

QECNoiseModel[a_Association]["RandomErrors", n_Integer, count_Integer] :=
    QECPauliString /@ noiseSampleRows[a, n, count]

(* The engine's own channel, for the qubits given: this is how the model reaches
   PauliStabilizer states, where qc[ps] returns the {probability, state} mixture.

   Mind the depolarizing convention.  The engine's QuantumChannel["Depolarizing"[q]]
   has Kraus coefficients Sqrt[1 - 3q/4] and Sqrt[q]/2, i.e. it leaves the qubit alone
   with probability 1 - 3q/4 and applies each of X, Y, Z with probability q/4.  This
   model says "something happens with probability p, and it is X, Y or Z with equal
   odds", so p = 3q/4 and the engine's parameter is q = 4p/3.  Passing p straight
   through would quietly simulate a weaker channel than the one asked for. *)
noiseChannelName["Depolarizing", p_] := "Depolarizing"[4 (1 - p[[1]]) / 3]
noiseChannelName["BitFlip", p_] := "BitFlip"[p[[2]]]
noiseChannelName["PhaseFlip", p_] := "PhaseFlip"[p[[4]]]
noiseChannelName["BitPhaseFlip", p_] := "BitPhaseFlip"[p[[3]]]
noiseChannelName[_, _] := Missing["NotAvailable", "a general Pauli channel has no named QuantumChannel form"]

QECNoiseModel[a_Association]["QuantumChannel", qubits_ : {1}] := With[
    {spec = noiseChannelName[a["Name"], a["Probabilities"]]},
    If[ MissingQ[spec],
        spec,
        Wolfram`QuantumFramework`QuantumChannel[spec, qubits]
    ]
]

QECNoiseModel[a_Association][prop_String] := (Message[QECNoiseModel::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECNoiseModel /: MakeBoxes[obj : QECNoiseModel[a_Association] /; KeyExistsQ[a, "Probabilities"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECNoiseModel,
        obj,
        BarChart[Replace[Rest[a["Probabilities"]], _Symbol -> 0.3, {1}], ChartLabels -> {"X", "Y", "Z"},
            ImageSize -> {Automatic, 34}, Axes -> False, ChartStyle -> RGBColor[0.15, 0.5, 0.65]],
        {
            BoxForm`SummaryItem[{"Model: ", a["Name"]}],
            BoxForm`SummaryItem[{"Level: ", noiseLevel[a]}],
            BoxForm`SummaryItem[{"Symbolic: ", noiseSymbolicQ[a]}]
        },
        {
            BoxForm`SummaryItem[{"P(I): ", a["Probabilities"][[1]]}],
            BoxForm`SummaryItem[{"P(X): ", a["Probabilities"][[2]]}],
            BoxForm`SummaryItem[{"P(Y): ", a["Probabilities"][[3]]}],
            BoxForm`SummaryItem[{"P(Z): ", a["Probabilities"][[4]]}],
            If[noiseLevel[a] === "CodeCapacity", Nothing,
                BoxForm`SummaryItem[{"Readout error: ", noiseMeasurementError[a]}]]
        },
        form,
        "Interpretable" -> False
    ]

(* Replace at level 1, not ReplaceAll.  A symbolic rate has to be swapped for a number
   before BarChart sees it, but ReplaceAll rewrites heads too and List is itself a
   Symbol, so {p, 0, 0} /. _Symbol -> 0.3 gives 0.3[0.3, 0, 0] and the chart fails to
   render.  The numbers listed under the chart are always the real ones. *)

QECNoiseModel /: MakeBoxes[obj : QECNoiseModel[a_Association] /; KeyExistsQ[a, "Rates"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECNoiseModel,
        obj,
        BarChart[Replace[Values[a["Rates"]], _Symbol -> 0.3, {1}], ChartLabels -> Keys[a["Rates"]],
            ImageSize -> {Automatic, 34}, Axes -> False, ChartStyle -> RGBColor[0.15, 0.5, 0.65]],
        {
            BoxForm`SummaryItem[{"Model: ", a["Name"]}],
            BoxForm`SummaryItem[{"Level: ", "Circuit"}],
            BoxForm`SummaryItem[{"Symbolic: ", noiseSymbolicQ[a]}]
        },
        KeyValueMap[BoxForm`SummaryItem[{# <> ": ", #2}] &, a["Rates"]],
        form,
        "Interpretable" -> False
    ]
