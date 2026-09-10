(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECDetectorModel]

PackageScope[detectorData]
PackageScope[codeDetectorModel]
PackageScope[circuitFaultMechanisms]
PackageScope[idleMechanisms]
PackageScope[faultEffect]
PackageScope[recordDetectors]
PackageScope[roundLength]
PackageScope[defaultRounds]
PackageScope[$oneQubitPaulis]
PackageScope[$twoQubitPaulis]


(* ============================================================================ *)
(* The detector error model.                                                    *)
(*                                                                              *)
(* Once the checks are measured by a circuit rather than by fiat, a syndrome bit *)
(* is no longer trustworthy on its own: it can be wrong because the data was hit *)
(* or because the readout lied, and one round cannot tell those apart.  Repeating *)
(* the extraction is the old fix (Gottesman, QECC book sec. 12.2.2); decoding the *)
(* *differences* rather than the syndromes is a newer one, and the two are not     *)
(* the same idea -- see WHOSE FIX THIS IS below before citing the book for it.    *)
(*                                                                              *)
(*     D_1     = S_1                                                             *)
(*     D_r     = S_r xor S_(r-1)          2 <= r <= R                            *)
(*     D_(R+1) = S_final xor S_R                                                 *)
(*                                                                              *)
(* A detector is a parity that is deterministically 0 when nothing goes wrong,   *)
(* so a fired detector is unambiguous evidence of a fault, and a readout error   *)
(* fires two detectors in a row rather than looking like a data error.  The last  *)
(* one closes the experiment against a noiseless final readout: S_final is the    *)
(* syndrome of whatever Pauli is left on the data.                               *)
(*                                                                              *)
(* "Deterministically 0" carries an assumption worth naming: that the data enters *)
(* inside the code space.  D_1 = S_1 is a detector only because a memory          *)
(* experiment prepares a codeword first, and at circuit level no input-error       *)
(* mechanism is emitted at all.  So this models a memory experiment, not an error- *)
(* correction gadget handed a block that already carries errors -- which is the     *)
(* case the book's ECRP and ECCP (sec. 10.2.4) are written for.  Feeding an        *)
(* incoming error would need a mechanism at instruction index 0 and would make     *)
(* the round-1 detectors ambiguous by construction.                               *)
(*                                                                              *)
(* The object built here is the map from individual faults to what they do.      *)
(* Two facts make it small and exact:                                            *)
(*                                                                              *)
(*   - the circuit is Clifford and the noise is Pauli, so a fault is a frame that *)
(*     propagates by conjugation (see Circuit.wl) -- no state simulation;         *)
(*   - that propagation is GF(2)-*linear*, so the effect of a set of faults is    *)
(*     the XOR of their individual effects.                                      *)
(*                                                                              *)
(* So one pass over the elementary faults produces a matrix: rows are faults,     *)
(* columns are detectors followed by logical observables, and every noisy run of  *)
(* the experiment is a GF(2) sum of rows.  That matrix is exactly Stim's detector *)
(* error model, and it is also what a matching decoder consumes, which is why it  *)
(* is the right boundary to stop at: everything downstream -- sampling, decoding,  *)
(* exporting to Stim, handing a graph to PyMatching -- reads this and nothing      *)
(* else.                                                                          *)
(*                                                                              *)
(* WHOSE FIX THIS IS.  The book's answer to an untrustworthy syndrome is repeat   *)
(* and *agree*: measure the whole syndrome, repeat, and trust a run of t+1         *)
(* consecutive identical syndromes, since at least one of them was taken with no    *)
(* faults (sec. 12.2.2, sec. 12.2.3).  Everything outside the winning run is        *)
(* discarded and a single (n-k)-bit syndrome is then looked up.  Differencing        *)
(* consecutive syndromes and decoding the whole spacetime history is a different    *)
(* procedure, from the topological-code literature -- Dennis, Kitaev, Landahl and   *)
(* Preskill, "Topological quantum memory", arXiv:quant-ph/0110143, where repeated   *)
(* noisy measurement becomes a lattice in one more dimension -- and as an explicit  *)
(* object with an error model attached it is Gidney, "Stim: a fast stabilizer       *)
(* circuit simulator", arXiv:2103.02202.                                           *)
(*                                                                              *)
(* The book knows about the strategy and declines to analyse it, sec. 12.2.2: "we  *)
(* might even hope to track changes to the syndrome as new faults occur ...         *)
(* However, it is difficult to analyze such strategies in the completely general    *)
(* case, so I'll stick with the simple rule given above."  So this file computes,   *)
(* exactly, the thing the standard reference sets aside as intractable in general.  *)
(* It also means the ECRP and ECCP proofs of sec. 12.2.3 do not transfer to it:     *)
(* what is exact here is the error rate of a stated protocol, not a fault-tolerance *)
(* guarantee.                                                                      *)
(*                                                                              *)
(* THE FAULT MODEL is a Pauli specialisation of the book's basic model, sec.       *)
(* 10.1.1: faults attach to circuit *locations* (Definition 10.1 -- preparation,    *)
(* gate, wait, measurement), and each faulty location performs the correct action   *)
(* followed by an error (Definition 10.2), with one rate per location type          *)
(* (Definition 10.3 gives p_P, p_G, p_S, p_M; sec. 10.1.1 explicitly anticipates    *)
(* splitting p_G by arity, which is the OneQubit/TwoQubit split here).  Two honest  *)
(* differences from the book.  It "makes no prescription for what happens at        *)
(* faulty locations", so depolarizing is a choice, not the definition -- proofs in  *)
(* that tradition use an adversarial model instead.  And the book counts the        *)
(* identity as an admissible outcome at a faulty location, whereas rate/3 and       *)
(* rate/15 here exclude it, which makes this marginally *less* conservative.        *)
(*                                                                              *)
(* One approximation, and it is the same one Stim makes: the three outcomes of a  *)
(* depolarizing channel are listed as three independent mechanisms rather than    *)
(* one three-way choice.  Sampling draws the true channel, so the simulated       *)
(* numbers are exact; only the probabilities attached to the rows are the         *)
(* independent-mechanism approximation, and the minimum-weight decoder does not   *)
(* read them.                                                                     *)
(* ============================================================================ *)

QECDetectorModel::usage = "QECDetectorModel[code, noise] gives the detector error model of a code under a phenomenological or circuit-level noise model, over the default number of rounds.\nQECDetectorModel[code, noise, r] uses r rounds of syndrome extraction followed by one noiseless readout.\nEach row is an elementary fault, with the detectors it fires and the logical observables it flips.\ndem[prop] gives a property; dem[\"Properties\"] lists them.";

QECDetectorModel::level = "A detector model needs a phenomenological or circuit-level noise model; this one is `1`. Use QECNoiseModel[name, p, \"MeasurementError\" -> q] or QECNoiseModel[\"Circuit\", p].";
QECDetectorModel::rounds = "The number of rounds must be a positive integer; got `1`.";
QECDetectorModel::noprop = "`1` is not a property of QECDetectorModel. Use dem[\"Properties\"] for the list.";


(* ---- elementary faults ---- *)

$oneQubitPaulis = {{1, 0}, {1, 1}, {0, 1}};

$twoQubitPaulis = DeleteCases[Tuples[{{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 2], {{0, 0}, {0, 0}}];

roundLength[a_Association, rounds_Integer] := Length[codeCircuitInstructions[a, rounds]] / rounds

(* A mechanism is <|"Faults" -> {{index, qubit, {x, z}}, ...}, "Probability" -> p,
   "Location" -> label|>.  Faults are placed *after* the instruction at the given
   index; index 0 is before the circuit starts. *)

circuitMechanisms[instr_List, rates_Association] := Catenate @ Table[
    With[{op = instr[[i, 1]]},
        Which[
            MemberQ[$circuitOneQubitOps, op],
                Table[<|
                    "Faults" -> {{i, instr[[i, 2]], p}},
                    "Probability" -> rates["OneQubit"] / 3,
                    "Location" -> {"OneQubit", i}
                |>, {p, $oneQubitPaulis}],

            MemberQ[$circuitTwoQubitOps, op],
                Table[<|
                    "Faults" -> {{i, instr[[i, 2]], pp[[1]]}, {i, instr[[i, 3]], pp[[2]]}},
                    "Probability" -> rates["TwoQubit"] / 15,
                    "Location" -> {"TwoQubit", i}
                |>, {pp, $twoQubitPaulis}],

            (* Placed before the measurement, so it flips the record and is then
               wiped by the next reset: a classical readout error. *)
            MemberQ[$circuitMeasureOps, op],
                {<|
                    "Faults" -> {{i - 1, instr[[i, 2]], {1, 0}}},
                    "Probability" -> rates["Measurement"],
                    "Location" -> {"Measurement", i}
                |>},

            op === "R",
                {<|
                    "Faults" -> {{i, instr[[i, 2]], {1, 0}}},
                    "Probability" -> rates["Reset"],
                    "Location" -> {"Reset", i}
                |>},

            True, {}
        ]
    ],
    {i, Length[instr]}
]

(* Phenomenological: the data is hit between rounds and the readout can lie, but
   the extraction itself is treated as a black box.  The data faults sit at the
   instruction index just before each round begins. *)
phenomenologicalMechanisms[instr_List, a_Association, noise_Association, rounds_Integer] := Module[
    {n = a["Qubits"], len, probs, q},
    len = Length[instr] / rounds;
    probs = noise["Probabilities"];
    q = noiseMeasurementError[noise];
    Join[
        Flatten[Table[
            <|
                "Faults" -> {{(r - 1) len, qubit, $oneQubitPaulis[[k]]}},
                "Probability" -> probs[[k + 1]],
                "Location" -> {"Data", r, qubit}
            |>,
            {r, rounds}, {qubit, n}, {k, 3}
        ], 2],
        Table[
            <|
                "Faults" -> {{i - 1, instr[[i, 2]], {1, 0}}},
                "Probability" -> q,
                "Location" -> {"Measurement", i}
            |>,
            {i, Select[Range[Length[instr]], MemberQ[$circuitMeasureOps, instr[[#, 1]]] &]}
        ]
    ]
]

(* Idle: a qubit that no instruction of a time step touches is waiting, and a
   waiting qubit decoheres.  One location per (layer, qubit), the three Paulis as its
   mutually exclusive outcomes, and the fault anchored by circuitIdleAnchor -- after
   the last instruction touching that qubit in an earlier layer, or at index 0 when
   none has run yet.  See the note there for why "the end of the layer" is the wrong
   anchor: the instruction list is not ordered by layer.

   Emitted only when the idle rate is not literally zero.  With Idle -> 0 (the
   default) this returns {} rather than a list of zero-probability rows, so the
   mechanism list, the detector model and every number downstream are bit-for-bit
   what they were before this existed.  A symbolic rate is never zero structurally,
   so it is always emitted.

   Ancilla idles are emitted too and mostly disappear in codeDetectorModel's filter,
   since a reset wipes them; the ones that survive fall between an ancilla's R and
   its M, act as readout errors, and belong in the model. *)
idleMechanisms[instr_List, nq_Integer, rates_Association] := If[
    rates["Idle"] === 0,
    {},
    Catenate[
        Function[slot,
            Table[
                <|
                    "Faults" -> {{slot[[1]], slot[[3]], pauli}},
                    "Probability" -> rates["Idle"] / 3,
                    "Location" -> {"Idle", slot[[2]], slot[[3]]}
                |>,
                {pauli, $oneQubitPaulis}
            ]
        ] /@ circuitIdleSlots[instr, nq]
    ]
]

circuitFaultMechanisms[a_Association, noise_Association, rounds_Integer] := With[
    {instr = codeCircuitInstructions[a, rounds],
     nq = a["Qubits"] + codeStabilizerCount[a]},
    Switch[noiseLevel[noise],
        "Circuit",
            With[{rates = noiseCircuitRates[noise]},
                Join[circuitMechanisms[instr, rates], idleMechanisms[instr, nq, rates]]
            ],
        _, phenomenologicalMechanisms[instr, a, noise, rounds]
    ]
]


(* ---- what a fault does ---- *)

(* Detector vector of a measurement record: first round bare, later rounds
   differenced, and a final difference against the noiseless readout of whatever
   Pauli survives on the data. *)
recordDetectors[record_List, finalSyndrome_List, m_Integer, rounds_Integer] := Join[
    record[[1 ;; m]],
    Catenate @ Table[
        BitXor[record[[(r - 1) m + 1 ;; r m]], record[[(r - 2) m + 1 ;; (r - 1) m]]],
        {r, 2, rounds}
    ],
    BitXor[finalSyndrome, record[[(rounds - 1) m + 1 ;; rounds m]]]
]

(* The residual data Pauli, read through the code's label matrix: the first m bits
   are the syndrome a perfect final round would see, the rest are the logical class,
   i.e. which logical operators the survivor anticommutes with. *)
faultEffect[a_Association, instr_List, nq_Integer, m_Integer, rounds_Integer, faults_List] := Module[
    {run, x, z, labels},
    run = framePropagate[instr, nq, faults];
    {x, z} = run["Frame"];
    labels = Mod[Join[x[[1 ;; a["Qubits"]]], z[[1 ;; a["Qubits"]]]] . Transpose[codeLabelMatrix[a]], 2];
    {
        recordDetectors[run["Record"], labels[[1 ;; m]], m, rounds],
        labels[[m + 1 ;; -1]],
        (* heralds need no differencing: a verification check reads 0 when nothing goes
           wrong, so the raw outcome already is the "something is wrong" bit *)
        run["Heralds"]
    }
]


(* ---- assembling the model ---- *)

(* Two kinds of row are dropped.

   Faults that fire no detector and flip no observable: the ones a reset wipes
   before they reach anything.  Keeping them only inflates the subset search a
   decoder has to do.  A fault that fires nothing but *does* flip an observable is
   kept and is the interesting one -- it is an undetectable logical fault, and its
   lowest weight is the circuit-level distance.

   And faults the channel cannot actually produce, i.e. probability literally zero:
   bit-flip noise never makes a Y or a Z.  These matter more than they look, because
   the minimum-weight decoder does not read probabilities -- a zero-probability row
   left in the table can claim a detector pattern that a real fault should have had,
   and quietly make the decoder worse than the noise deserves.  A symbolic rate is
   never dropped, since nothing can be concluded about it. *)
codeDetectorModel[a_Association, noise_Association, rounds_Integer] := codeDetectorModel[a, noise, rounds] = Module[
    {instr, nq, m, mechanisms, effects, keep},

    instr = codeCircuitInstructions[a, rounds];
    nq = a["Qubits"] + codeStabilizerCount[a];
    m = codeStabilizerCount[a];

    mechanisms = circuitFaultMechanisms[a, noise, rounds];
    effects = faultEffect[a, instr, nq, m, rounds, #["Faults"]] & /@ mechanisms;

    (* A row that only trips a herald is NOT effectless: it discards the shot, which is
       an outcome the rate has to account for.  So the filter looks at all three. *)
    keep = Select[
        Range[Length[mechanisms]],
        mechanisms[[#, "Probability"]] =!= 0 &&
            ! (AllTrue[effects[[#, 1]], # === 0 &] && AllTrue[effects[[#, 2]], # === 0 &] &&
               AllTrue[effects[[#, 3]], # === 0 &]) &
    ];

    <|
        "DetectorMatrix" -> effects[[keep, 1]],
        "ObservableMatrix" -> effects[[keep, 2]],
        "HeraldMatrix" -> effects[[keep, 3]],
        "Heralds" -> Count[instr, {"MH", ___}],
        "Probabilities" -> Lookup[mechanisms[[keep]], "Probability"],
        "Locations" -> Lookup[mechanisms[[keep]], "Location"],
        "Detectors" -> (rounds + 1) m,
        "Observables" -> 2 codeLogicalQubits[a],
        "Rounds" -> rounds,
        "Checks" -> m,
        "Code" -> a,
        "Noise" -> noise
    |>
]


(* ---- construction ---- *)

(* Default rounds: the code distance.  This is the topological-memory convention
   (Dennis-Kitaev-Landahl-Preskill, arXiv:quant-ph/0110143): with repeated noisy
   measurement the history is a lattice in one more dimension, and taking fewer
   rounds than the distance protects the time direction less than the space ones,
   which quietly caps the whole thing.

   It is deliberately NOT the book's fault-tolerant-error-correction criterion, and
   the two should not be conflated.  Gottesman, QECC book sec. 12.2.2 asks instead
   for a run of t+1 *consecutive agreeing* syndrome measurements, proves (t+1)^2
   repetitions sufficient to find one, and explicitly rejects the count d = 2t+1
   with a majority vote -- "since there are 2^(n-k) possible values of the syndrome,
   it is not enough to repeat 2t+1 times and take the majority ... There might not
   be a majority."  That criterion buys the ECRP and ECCP; this default buys a
   memory experiment whose time and space directions are equally protected.  Two
   different questions, and "Rounds" is there to ask either. *)
defaultRounds[a_Association] := With[{d = codeDistance[a]}, If[IntegerQ[d] && d > 0, d, 1]]

QECDetectorModel[QECCode[a_Association], noise_QECNoiseModel] :=
    QECDetectorModel[QECCode[a], noise, defaultRounds[a]]

QECDetectorModel[QECCode[a_Association], QECNoiseModel[noise_Association], rounds_] := Which[
    ! (IntegerQ[rounds] && rounds > 0),
        Message[QECDetectorModel::rounds, rounds]; $Failed,
    noiseLevel[noise] === "CodeCapacity",
        Message[QECDetectorModel::level, "code-capacity"]; $Failed,
    True,
        QECDetectorModel[codeDetectorModel[a, noise, rounds]]
]


(* ---- properties ---- *)

$detectorProperties = {
    "DetectorMatrix", "ObservableMatrix", "HeraldMatrix", "Probabilities", "Locations",
    "Detectors", "Observables", "Heralds", "Rounds", "Checks", "Faults", "Code", "Noise",
    "LocationCounts", "UndetectableFaults", "DetectorRates", "ObservableRates",
    "HeraldRates", "PostSelectedQ", "Properties"
};

(* How often each detector fires, exactly.  A detector fires when an odd number of
   the mechanisms touching it fire, and odd-parity probabilities multiply through
   the (1 - 2p) form, one factor per independent location:

       P(odd) = (1 - prod_locations (1 - 2 sum_(i in location, i touches d) p_i)) / 2

   the inner sum being over the mutually exclusive outcomes of that one location.
   Exact, symbolic if the rates are, and the natural quantity to check an external
   simulator against: it is basis-free and every detector is a separate test. *)
oddParityRate[a_Association, column_List] := With[
    {p = a["Probabilities"]},
    (1 - Product[
        1 - 2 Total[p[[Select[rows, column[[#]] === 1 &]]]],
        {rows, demGroups[a]}
    ]) / 2
]

detectorFiringRates[a_Association] := If[
    a["DetectorMatrix"] === {},
    ConstantArray[0, a["Detectors"]],
    Table[oddParityRate[a, a["DetectorMatrix"][[All, k]]], {k, a["Detectors"]}]
]

heraldFiringRates[a_Association] := If[
    a["HeraldMatrix"] === {} || a["Heralds"] === 0,
    ConstantArray[0, a["Heralds"]],
    Table[oddParityRate[a, a["HeraldMatrix"][[All, k]]], {k, a["Heralds"]}]
];

observableFlipRates[a_Association] := If[
    a["ObservableMatrix"] === {},
    ConstantArray[0, a["Observables"]],
    Table[oddParityRate[a, a["ObservableMatrix"][[All, k]]], {k, a["Observables"]}]
];

detectorData[QECDetectorModel[a_Association]] := a

QECDetectorModel[_Association]["Properties"] := $detectorProperties

QECDetectorModel[a_Association][prop : ("DetectorMatrix" | "ObservableMatrix" |
    "HeraldMatrix" | "Probabilities" | "Locations" | "Detectors" | "Observables" |
    "Heralds" | "Rounds" | "Checks")] := a[prop]

(* Whether any of the circuit's measurements is a herald, i.e. whether a rate computed
   from this model is conditional on acceptance. *)
QECDetectorModel[a_Association]["PostSelectedQ"] := a["Heralds"] > 0

QECDetectorModel[a_Association]["Faults"] := Length[a["Probabilities"]]
QECDetectorModel[a_Association]["Code"] := QECCode[a["Code"]]
QECDetectorModel[a_Association]["Noise"] := QECNoiseModel[a["Noise"]]
QECDetectorModel[a_Association]["LocationCounts"] := Counts[First /@ a["Locations"]]
QECDetectorModel[a_Association]["DetectorRates"] := detectorFiringRates[a]
QECDetectorModel[a_Association]["ObservableRates"] := observableFlipRates[a]
QECDetectorModel[a_Association]["HeraldRates"] := heraldFiringRates[a]

(* Faults that no detector can see.  Their existence is not a bug: they are the
   circuit-level analogue of a logical operator, and how many faults it takes to
   build one is the circuit-level distance. *)
QECDetectorModel[a_Association]["UndetectableFaults"] := With[
    {rows = Select[Range[Length[a["Probabilities"]]], AllTrue[a["DetectorMatrix"][[#]], # === 0 &] &]},
    a["Locations"][[rows]]
]

QECDetectorModel[a_Association][prop_String] := (Message[QECDetectorModel::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECDetectorModel /: MakeBoxes[obj : QECDetectorModel[a_Association] /; KeyExistsQ[a, "DetectorMatrix"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECDetectorModel,
        obj,
        ArrayPlot[If[a["DetectorMatrix"] === {}, {{0}}, a["DetectorMatrix"]],
            ColorRules -> {0 -> White, 1 -> RGBColor[0.15, 0.5, 0.65]},
            ImageSize -> {Automatic, 34}, Frame -> False],
        {
            BoxForm`SummaryItem[{"Faults: ", Length[a["Probabilities"]]}],
            BoxForm`SummaryItem[{"Detectors: ", a["Detectors"]}],
            BoxForm`SummaryItem[{"Rounds: ", a["Rounds"]}]
        },
        {
            BoxForm`SummaryItem[{"Level: ", noiseLevel[a["Noise"]]}],
            BoxForm`SummaryItem[{"Observables: ", a["Observables"]}],
            BoxForm`SummaryItem[{"Locations: ", Counts[First /@ a["Locations"]]}]
        },
        form,
        "Interpretable" -> False
    ]
