(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECStimCircuit]

PackageScope[stimGateName]
PackageScope[codeStimCircuit]


(* ============================================================================ *)
(* The Stim bridge.                                                             *)
(*                                                                              *)
(* Stim (Gidney, "Stim: a fast stabilizer circuit simulator", arXiv:2103.02202) is *)
(* the community's reference stabilizer simulator, and the roadmap is             *)
(* explicit that we connect to it rather than compete with it: it samples fast   *)
(* and large, we derive exactly and symbolically.  Writing the memory experiment *)
(* out in its syntax buys two things at once -- an independent check on           *)
(* everything above (its detector error model has to agree with ours mechanism   *)
(* for mechanism), and a path to the decoders built on top of it, since           *)
(* PyMatching consumes exactly that model.                                       *)
(*                                                                              *)
(* Structure of the emitted circuit:                                             *)
(*                                                                              *)
(*   1. the encoder, noiseless -- Stim starts in |0...0>, which is not a codeword *)
(*      for a general code, and a memory experiment has to start inside the code  *)
(*      space or its first round of detectors means nothing;                     *)
(*   2. r rounds of syndrome extraction with noise on every location;            *)
(*   3. one noiseless round, which is the perfect final readout;                 *)
(*   4. one noiseless measurement of the logical operator, declared as the        *)
(*      observable.                                                              *)
(*                                                                              *)
(* One honest asymmetry with the rest of the package.  Our own logical error rate *)
(* asks whether the decoder got the residual's *class* right, which covers both   *)
(* logical X and logical Z damage at once, because that is the right question for *)
(* an unknown logical state.  A Stim memory experiment is single-basis by          *)
(* construction: prepare |0_L>, keep it, measure Z-bar.  So the emitted circuit    *)
(* declares one observable and "Observable" chooses which.  The detector part is   *)
(* basis-independent and is where the two models are compared.                    *)
(* ============================================================================ *)

QECStimCircuit::usage = "QECStimCircuit[code, noise] gives, as a string, the Stim source of a memory experiment: the encoder, rounds of noisy syndrome extraction, a noiseless final round, and a logical observable.\nQECStimCircuit[code, noise, r] uses r noisy rounds.\nQECStimCircuit[code] emits the noiseless circuit.\nThe option \"Observable\" -> \"Z\" or \"X\" chooses which logical operator is kept.";

QECStimCircuit::observable = "\"Observable\" must be \"X\" or \"Z\"; got `1`.";
QECStimCircuit::herald = "This circuit contains `1` herald measurement(s) (\"MH\"). Stim has no post-selection primitive, and emitting them as ordinary measurements would shift every rec[-k] index the detectors are built from. Export a circuit without heralds, or declare them as detectors and post-select outside Stim.";
QECStimCircuit::level = "Stim has no notion of a code-capacity noise model, whose checks are perfect by assumption. Use a phenomenological or circuit-level model.";


(* ---- gate names ---- *)

(* V is the engine's sqrt(X).  Stim's SQRT_X is the same Clifford up to a Pauli,
   and a Pauli difference is invisible to detectors, which are defined relative to
   the circuit's own noiseless run. *)
$stimGates = <|
    "H" -> "H", "S" -> "S", "V" -> "SQRT_X", "Vdg" -> "SQRT_X_DAG",
    "X" -> "X", "Y" -> "Y", "Z" -> "Z",
    "CNOT" -> "CX", "CZ" -> "CZ", "SWAP" -> "SWAP"
|>;

stimGateName[op_String] := Lookup[$stimGates, op, Missing["NotAvailable", op]]

(* Stim counts qubits from zero. *)
qb[q_Integer] := ToString[q - 1]

qbs[qs_List] := StringRiffle[qb /@ qs, " "]

num[x_] := ToString[N[x], InputForm]


(* ---- pieces ---- *)

(* The encoder as Stim lines.  These are arrows name -> order, the same shape
   "ApplyCircuit" takes. *)
encoderLines[a_Association] := Replace[
    codeEncodingGates[a],
    {
        (op_String -> q_Integer) :> stimGateName[op] <> " " <> qb[q],
        (op_String -> qs_List) :> stimGateName[op] <> " " <> qbs[qs]
    },
    {1}
]

(* One instruction, with its noise attached.  A measurement carries its own flip
   probability in Stim (M(p)), which is exactly the readout error, so it needs no
   separate line. *)
stepLines[step_List, rates_Association] := Switch[First[step],
    "R",
        {"R " <> qb[step[[2]]],
         If[rates["Reset"] === 0, Nothing, "X_ERROR(" <> num[rates["Reset"]] <> ") " <> qb[step[[2]]]]},
    "M",
        {If[rates["Measurement"] === 0,
            "M " <> qb[step[[2]]],
            "M(" <> num[rates["Measurement"]] <> ") " <> qb[step[[2]]]]},
    "CNOT",
        {"CX " <> qbs[{step[[2]], step[[3]]}],
         If[rates["TwoQubit"] === 0, Nothing,
            "DEPOLARIZE2(" <> num[rates["TwoQubit"]] <> ") " <> qbs[{step[[2]], step[[3]]}]]},
    "CZ",
        {"CZ " <> qbs[{step[[2]], step[[3]]}],
         If[rates["TwoQubit"] === 0, Nothing,
            "DEPOLARIZE2(" <> num[rates["TwoQubit"]] <> ") " <> qbs[{step[[2]], step[[3]]}]]},
    _,
        {stimGateName[First[step]] <> " " <> qb[step[[2]]],
         If[rates["OneQubit"] === 0, Nothing,
            "DEPOLARIZE1(" <> num[rates["OneQubit"]] <> ") " <> qb[step[[2]]]]}
]

(* Idling.  One DEPOLARIZE1 per idle slot, at the slot's anchor -- the same anchor
   circuitIdleSlots hands the detector model, so the exported circuit's error model
   and ours have the same locations by construction.

   A qubit idle across several layers with one anchor gets several lines on the same
   qubit, and that is right: they are independent chances to decohere during one gap,
   and Stim reads consecutive DEPOLARIZE1 as independent channels exactly as the
   detector model reads them as separate locations.  One line per slot rather than
   repeated targets on one line, since a repeated target is not a documented way to
   ask Stim for a repeated channel.

   No TICK is emitted.  TICK would claim a time step, and the lines here are in
   *emission* order, not layer order -- the two disagree, and reordering to fix that
   would move measurements relative to each other and silently invalidate every
   rec[-k] index the detectors are built from. *)
idleLines[rates_Association, qs_List] := If[
    rates["Idle"] === 0 || qs === {},
    {},
    "DEPOLARIZE1(" <> num[rates["Idle"]] <> ") " <> qb[#] & /@ qs
]

instructionLines[instr_List, rates_Association] := instructionLines[instr, rates, <||>, 0]

(* offset is where this block starts in the whole circuit's instruction list, so a
   per-round block can look up anchors computed over all rounds at once. *)
instructionLines[instr_List, rates_Association, idleAt_Association, offset_Integer] := Join[
    If[offset === 0, idleLines[rates, Lookup[idleAt, 0, {}]], {}],
    Catenate @ Table[
        Join[
            stepLines[instr[[i]], rates],
            idleLines[rates, Lookup[idleAt, offset + i, {}]]
        ],
        {i, Length[instr]}
    ]
]

$noRates = <|"OneQubit" -> 0, "TwoQubit" -> 0, "Measurement" -> 0, "Reset" -> 0, "Idle" -> 0|>;

(* Detectors for the round just emitted.  The record grows by m entries per round,
   so generator j of the last round is rec[-(m - j + 1)] and its predecessor is
   another m further back. *)
detectorLines[m_Integer, first_ : False] := Table[
    If[ first,
        "DETECTOR rec[" <> ToString[-(m - j + 1)] <> "]",
        "DETECTOR rec[" <> ToString[-(m - j + 1)] <> "] rec[" <> ToString[-(2 m - j + 1)] <> "]"
    ],
    {j, m}
]


(* ---- the circuit ---- *)

codeStimCircuit[a_Association, noise_Association, rounds_Integer, basis_String] := Module[
    {n = a["Qubits"], m = codeStabilizerCount[a], instr, rates, dataNoise, idleAt, logical, lines, round},

    logical = First[codeLogicalVectors[a][basis]];

    instr = codeCircuitInstructions[a, 1];
    rates = If[noiseLevel[noise] === "Circuit", noiseCircuitRates[noise], $noRates];

    (* Idle anchors are computed over the *whole* multi-round circuit, because the
       rounds pipeline: an ancilla's reset for round 2 schedules while the other
       ancilla of round 1 is still being measured, so r rounds are not r copies of
       one round's schedule.  The per-round blocks below index into it by offset. *)
    idleAt = If[
        rates["Idle"] === 0,
        <||>,
        GroupBy[circuitIdleSlots[codeCircuitInstructions[a, rounds], n + m], First -> Last]
    ];

    (* Phenomenological data noise: one Pauli channel across the data qubits at the
       start of each round.  Circuit-level noise has none, because there the data is
       hit by the gates themselves. *)
    dataNoise = If[
        noiseLevel[noise] === "Phenomenological",
        With[{p = noise["Probabilities"]},
            {"PAULI_CHANNEL_1(" <> StringRiffle[num /@ p[[2 ;; 4]], ", "] <> ") " <> qbs[Range[n]]}
        ],
        {}
    ];
    If[ noiseLevel[noise] === "Phenomenological",
        rates = <|$noRates, "Measurement" -> noiseMeasurementError[noise]|>
    ];

    round[r_] := Join[
        {"# --- round " <> ToString[r] <> " ---"},
        dataNoise,
        instructionLines[instr, rates, idleAt, (r - 1) Length[instr]],
        detectorLines[m, r === 1]
    ];

    lines = Join[
        {"# " <> ToString[n] <> " data qubits, " <> ToString[m] <> " check ancillas, " <>
            ToString[rounds] <> " noisy rounds, " <> basis <> "-basis memory",
         "# emitted by QECStimCircuit"},
        {"# --- encode ---"},
        encoderLines[a],
        Catenate[round /@ Range[rounds]],

        (* the perfect final readout *)
        {"# --- final noiseless round ---"},
        instructionLines[instr, $noRates],
        detectorLines[m, rounds === 0],

        (* the logical operator, measured noiselessly through one more ancilla *)
        {"# --- logical " <> basis <> " ---"},
        instructionLines[generatorInstructions[symplecticPart[logical], n, n + m + 1], $noRates],
        {"OBSERVABLE_INCLUDE(0) rec[-1]"}
    ];

    StringRiffle[lines, "\n"] <> "\n"
]


(* ---- construction ---- *)

Options[QECStimCircuit] = {"Observable" -> "Z"};

QECStimCircuit[code_QECCode, opts : OptionsPattern[]] :=
    QECStimCircuit[code, QECNoiseModel["Circuit", <||>], 1, opts]

QECStimCircuit[code_QECCode, noise_QECNoiseModel, opts : OptionsPattern[]] :=
    QECStimCircuit[code, noise, defaultRounds[First[code]], opts]

QECStimCircuit[code_QECCode, noise_QECNoiseModel, rounds_, opts : OptionsPattern[]] := With[
    {basis = OptionValue[QECStimCircuit, {opts}, "Observable"]},
    Which[
        ! MemberQ[{"X", "Z"}, basis],
            Message[QECStimCircuit::observable, basis]; $Failed,
        ! (IntegerQ[rounds] && rounds > 0),
            Message[QECSyndromeCircuit::rounds, rounds]; $Failed,
        noiseLevel[First[noise]] === "CodeCapacity",
            Message[QECStimCircuit::level]; $Failed,
        Count[codeCircuitInstructions[First[code], 1], {"MH", ___}] > 0,
            Message[QECStimCircuit::herald,
                Count[codeCircuitInstructions[First[code], 1], {"MH", ___}]]; $Failed,
        True,
            codeStimCircuit[First[code], First[noise], rounds, basis]
    ]
]

QECCode[a_Association]["StimCircuit", rest___] := QECStimCircuit[QECCode[a], rest]
