(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECFaultTolerant]

PackageScope[ftGateWord]
PackageScope[ftRelocate]
PackageScope[ftPreparation]
PackageScope[ftMeasurement]
PackageScope[ftCorrection]
PackageScope[ftReadoutSupport]
PackageScope[ftAssemble]
PackageScope[ftLocations]
PackageScope[faultTolerantData]
PackageScope[$ftGadgetKinds]


(* ============================================================================ *)
(* FT(C): the ideal fault-tolerant simulation of a circuit (Got26 Def 10.6).    *)
(*                                                                              *)
(* Everything before this file is a gadget.  Def 10.5 says what a fault-tolerant *)
(* protocol is -- a code plus gadgets for preparation, measurement, a universal  *)
(* gate set, storage and error correction -- and Def 10.6 says how to spend      *)
(* them: take each location of C, replace it with the corresponding gadget,      *)
(* replace each qubit of C with a BLOCK of the code, and after every             *)
(* preparation, gate or storage gadget put an error correction gadget on each    *)
(* block involved.  Not after measurement gadgets: their output is classical,    *)
(* and classical circuits are assumed not to fail.                              *)
(*                                                                              *)
(* This file is that substitution, and almost nothing else.  The gadgets are     *)
(* already built and already checked; what is new is the bookkeeping, which is   *)
(* where a hand-assembled circuit goes wrong.                                    *)
(*                                                                              *)
(* THE GADGET IS CHOSEN BY WHAT IT DOES, NOT BY ITS NAME.  This is the one       *)
(* place where phase E pays for itself.  The gadget for a logical S on the       *)
(* 7-qubit code is the transversal S^dagger, because the transversal S performs  *)
(* the logical S^dagger (sec. 11.3, eq. 11.22): Y^7 = -Ybar.  So the assembler   *)
(* never emits "the transversal version of the gate it was asked for" -- it asks *)
(* QECTransversalGate which word has the requested LOGICAL action, and emits     *)
(* that.  A protocol assembled by name would silently compute the conjugate      *)
(* circuit.                                                                      *)
(*                                                                              *)
(* WAITS ARE LOCATIONS, so storage gadgets are real even though they emit no     *)
(* instructions: the book's storage gadget "can always be implemented by just    *)
(* putting a wait location for all physical qubits in the code", and this layer  *)
(* already treats idling that way -- circuitIdleSlots is the single source of    *)
(* truth for where a wait is, shared with the detector model and the Stim        *)
(* writer.  A live block idle in a layer therefore gets a storage gadget in the  *)
(* table and an error correction gadget after it, which is what makes the size   *)
(* overhead an honest number rather than a count of typed-out gates.             *)
(*                                                                              *)
(* THE THREE OVERHEADS are Def 10.6's own, so they are properties rather than    *)
(* something a caller recomputes: size (locations of FT(C) over locations of C), *)
(* qubits, and depth.                                                            *)
(*                                                                              *)
(* WHAT IS NOT HERE, stated rather than hidden:                                  *)
(*                                                                              *)
(*   - Clifford gates only.  Def 10.5 asks for a UNIVERSAL set; transversal      *)
(*     gates cannot give one (that is the point of ch. 13), so a circuit with a  *)
(*     T gate is refused rather than approximated.                               *)
(*   - The preparation gadget is the code's own encoder, which is not fault      *)
(*     tolerant.  Same open assumption the EC gadget already reports, inherited  *)
(*     rather than quietly dropped.                                              *)
(*   - Corrections stay in the Pauli frame (sec. 12.5.2).  No instruction here   *)
(*     is classically conditioned: for a Clifford circuit with Pauli errors,     *)
(*     tracking the correction is equivalent to applying it, and an applied      *)
(*     correction is one more location that can fail.                           *)
(*   - One logical qubit per block.  Def 10.4's own simplification, and what     *)
(*     makes "the gadget for logical g" a well-posed question.                   *)
(*                                                                              *)
(* References: Got26 Def 10.4 (gadget), Def 10.5 (protocol), Def 10.6 (FT(C)     *)
(* and the three overheads), fig. 10.2; sec. 11.3 (which word to emit);          *)
(* sec. 12.3 (the EC gadget); sec. 12.5.2 (the frame); sec. 14.5.3 (the numbers  *)
(* this is aimed at reproducing).                                                *)
(* ============================================================================ *)

QECFaultTolerant::usage = "QECFaultTolerant[circuit, code] gives the fault-tolerant simulation FT(C) of an ideal Clifford circuit, with each logical qubit encoded in a block of the code, each location replaced by its gadget, and an error correction gadget after every preparation, gate and storage gadget.\nThe circuit is a list of instructions on logical qubits, such as {{\"R\", 1}, {\"H\", 1}, {\"CNOT\", 1, 2}, {\"M\", 1}}.\nft[prop] gives a property; ft[\"Properties\"] lists them.";

QECFaultTolerant::css = "The fault-tolerant protocol assembled here uses Steane error correction, which needs a CSS code.";
QECFaultTolerant::logical = "This protocol encodes one logical qubit per block; the code given carries `1`.";
QECFaultTolerant::location = "`1` is not a location this protocol has a gadget for. Locations are {\"R\", j}, {\"M\", j}, a one-qubit Clifford {g, j}, or a two-qubit gate {g, j, l}.";
QECFaultTolerant::gate = "The code has no transversal gadget for a logical `1`. Transversal gates are limited to the Clifford group, and to the Cliffords this code's symmetry admits; QECTransversalGate[code, All] lists them.";
QECFaultTolerant::qubit = "`1` is not a logical qubit of a circuit on `2`.";
QECFaultTolerant::unprepared = "Logical qubit `1` is used before it is prepared. Start the circuit with {\"R\", `1`}.";
QECFaultTolerant::measured = "Logical qubit `1` is used after it has been measured.";
QECFaultTolerant::readout = "The code's logical Z is not a product of Z operators, so a transversal Z measurement does not read it.";
QECFaultTolerant::noprop = "`1` is not a property of QECFaultTolerant. Use ft[\"Properties\"] for the list.";


$ftGadgetKinds = {"Preparation", "Gate", "Storage", "Measurement", "ErrorCorrection"};


(* ---- locations, counted the way Def 10.1 counts them ---- *)

(* A location is an instruction OR a wait, and the waits come from the schedule
   rather than from the instruction list -- which is why this is not Length. *)
ftLocations[instr_List, nq_Integer] := gadgetLocations[instr, nq]


(* ---- moving a gadget onto a block ---- *)

(* The EC gadget is written for data on qubits 1..n and its ancilla at an offset,
   because that is all it ever needed.  Here the data is block b.  Every
   instruction of that gadget touches either a data qubit (<= n) or an ancilla
   qubit (> offset), so relocating it is a two-branch map and not a rewrite. *)
ftRelocate[instr_List, n_Integer, offset_Integer, dataBase_Integer, ancillaBase_Integer] :=
    With[{move = If[# <= n, dataBase + #, ancillaBase + # - offset] &},
        Replace[instr, {
            {op_String, q_Integer} :> {op, move[q]},
            {op_String, x_Integer, y_Integer} :> {op, move[x], move[y]}
        }, {1}]
    ]


(* ---- the gadgets, one per kind of location ---- *)

(* Preparation of an encoded |0>: reset the block and run the code's encoder.  Not
   fault tolerant, and reported as such rather than assumed away. *)
ftPreparation[a_Association, b_Integer] := With[{n = a["Qubits"], base = (b - 1) a["Qubits"]},
    Join[
        Table[{"R", base + q}, {q, n}],
        codeEncodedZeroInstructions[a, base]
    ]
]

(* Measurement of the logical Z: measure every qubit of the block.  The logical bit
   is the parity of the outcomes over the support of Zbar, and the parities of the
   Z-type generators over the same outcomes are the X syndrome -- which is why a
   transversal measurement is already a fault-tolerant one and needs no gadget of
   its own. *)
ftMeasurement[a_Association, b_Integer] := With[{n = a["Qubits"], base = (b - 1) a["Qubits"]},
    Table[{"M", base + q}, {q, n}]
]

(* Where the logical bit lives inside those n outcomes. *)
ftReadoutSupport[a_Association] := With[
    {n = a["Qubits"], z = First[codeLogicalVectors[a]["Z"]]},
    If[ Total[z[[1 ;; n]]] =!= 0,
        Message[QECFaultTolerant::readout]; $Failed,
        Flatten @ Position[z[[n + 1 ;; 2 n]], 1]
    ]
]

(* The error correction gadget, moved onto block b with that block's own ancilla
   workspace -- one workspace per block, so that the ECs of one layer can run in
   parallel rather than queueing on a shared ancilla. *)
ftCorrection[a_Association, b_Integer, blocks_Integer, order_String] := With[
    {n = a["Qubits"]},
    ftRelocate[
        steaneInstructions[a, n, order], n, n,
        (b - 1) n, (blocks + b - 1) n
    ]
]

(* The gate gadget: the transversal word whose LOGICAL action is the gate asked
   for.  On the 7-qubit code the answer for "S" is the word "Sdg", and that is the
   whole reason this is a search rather than a lookup by name. *)
(* The named gates are searched first, so that the gadget for a logical S comes
   back as the single gate "Sdg" rather than as the word "SSS" that does the same
   thing with three times the locations.  The words are the fallback, and they are
   what makes the search exhaustive: they cover all 24 one-qubit Cliffords, named
   or not. *)
ftGateWord[a_Association, g_String] := ftGateWord[a, g] = SelectFirst[
    Join[$transversalOneQubitGates, transversalCliffordWords[1]],
    transversalValidQ[a, #, 1] &&
        transversalName[transversalLogicalAction[a, #, 1], 1] === g &,
    Missing["NoGadget", g]
]

ftGateInstructions[a_Association, g_String, b_Integer] := With[
    {n = a["Qubits"], word = ftGateWord[a, g]},
    If[ MissingQ[word],
        $Failed,
        registerLift[transversalInstructions[word, n], n, b]
    ]
]

ftTwoQubitInstructions[a_Association, g_String, c_Integer, t_Integer] := With[
    {n = a["Qubits"]},
    If[ transversalValidQ[a, g, 2] && transversalName[transversalLogicalAction[a, g, 2], 2] === g,
        Table[{g, registerIndex[n, c, q], registerIndex[n, t, q]}, {q, n}],
        $Failed
    ]
]


(* ---- the assembly ---- *)

(* One pass over the SCHEDULE of C rather than over its instruction list: Def 10.6
   puts an error correction gadget between every adjacent pair of locations, and
   what "adjacent" means in a circuit with parallel gates is a layer.  A block that
   is live and untouched in a layer is waiting, which is a storage gadget, which
   also earns an EC.

   The measurement record is tracked while emitting, not recovered afterwards: the
   n outcomes of a measurement gadget sit in the middle of a record otherwise made
   of EC outcomes, and their positions are what a decoder needs. *)
ftAssemble[a_Association, circuit_List, blocks_Integer, order_String] := Catch[Module[
    {n = a["Qubits"], support, layers, status, out, rows, readouts, records, emit, touched},

    support = ftReadoutSupport[a];
    If[support === $Failed, Throw[$Failed, "QECFaultTolerant"]];

    layers = circuitSchedule[circuit, blocks];
    status = ConstantArray["Unprepared", blocks];
    out = Internal`Bag[];
    rows = Internal`Bag[];
    readouts = <||>;
    records = 0;

    emit[kind_, location_, bs_List, instr_List] := (
        Internal`StuffBag[out, instr, 1];
        Internal`StuffBag[rows, <|
            "Gadget" -> kind, "Location" -> location, "Blocks" -> bs,
            "Instructions" -> Length[instr]
        |>];
        records += Count[instr, {"M", _}];
    );

    Do[
        touched = {};
        Do[
            Module[{step = circuit[[i]], op = First[circuit[[i]]], qs, instr},
                qs = Rest[step];
                If[ ! AllTrue[qs, 1 <= # <= blocks &],
                    Message[QECFaultTolerant::qubit, First[Select[qs, ! (1 <= # <= blocks) &]], blocks];
                    Throw[$Failed, "QECFaultTolerant"]
                ];
                touched = Join[touched, qs];
                Which[
                    op === "R",
                        emit["Preparation", step, qs, ftPreparation[a, First[qs]]];
                        status[[First[qs]]] = "Live",

                    AnyTrue[qs, status[[#]] === "Unprepared" &],
                        Message[QECFaultTolerant::unprepared, First[Select[qs, status[[#]] === "Unprepared" &]]];
                        Throw[$Failed, "QECFaultTolerant"],

                    AnyTrue[qs, status[[#]] === "Measured" &],
                        Message[QECFaultTolerant::measured, First[Select[qs, status[[#]] === "Measured" &]]];
                        Throw[$Failed, "QECFaultTolerant"],

                    op === "M",
                        emit["Measurement", step, qs, ftMeasurement[a, First[qs]]];
                        readouts[First[qs]] = records - n + support;
                        status[[First[qs]]] = "Measured",

                    Length[qs] === 1,
                        instr = ftGateInstructions[a, op, First[qs]];
                        If[instr === $Failed, Message[QECFaultTolerant::gate, op]; Throw[$Failed, "QECFaultTolerant"]];
                        emit["Gate", step, qs, instr],

                    Length[qs] === 2,
                        instr = ftTwoQubitInstructions[a, op, qs[[1]], qs[[2]]];
                        If[instr === $Failed, Message[QECFaultTolerant::gate, op]; Throw[$Failed, "QECFaultTolerant"]];
                        emit["Gate", step, qs, instr],

                    True,
                        Message[QECFaultTolerant::location, step];
                        Throw[$Failed, "QECFaultTolerant"]
                ]
            ],
            {i, layers[[L]]}
        ];

        (* Storage gadgets: a live block nobody touched this layer is waiting, and a
           wait is a location.  It emits nothing -- the idle slots of the assembled
           circuit are where its faults live -- but it is a gadget, and it earns an
           error correction gadget like any other. *)
        Do[
            If[status[[b]] === "Live" && ! MemberQ[touched, b],
                emit["Storage", {"Wait", b}, {b}, {}]
            ],
            {b, blocks}
        ];

        (* And the correction, on every block that is live at the end of the layer:
           after preparation, gate and storage gadgets, never after a measurement. *)
        Do[
            If[status[[b]] === "Live",
                emit["ErrorCorrection", {"EC", b}, {b}, ftCorrection[a, b, blocks, order]]
            ],
            {b, blocks}
        ],
        {L, Length[layers]}
    ];

    <|
        "Instructions" -> Internal`BagPart[out, All],
        "Gadgets" -> Internal`BagPart[rows, All],
        "Readouts" -> readouts,
        "Measurements" -> records
    |>
], "QECFaultTolerant"]


(* ---- construction ---- *)

Options[QECFaultTolerant] = {"Order" -> "BitFlipFirst"};

QECFaultTolerant[circuit_List, code_QECCode, opts : OptionsPattern[]] := Module[
    {a = First[code], blocks, order, assembled},

    order = OptionValue["Order"];
    blocks = Max[Cases[circuit, {_String, qs__Integer} :> Max[{qs}]], 0];

    Which[
        ! codeCSSQ[a], Message[QECFaultTolerant::css]; $Failed,
        codeLogicalQubits[a] =!= 1, Message[QECFaultTolerant::logical, codeLogicalQubits[a]]; $Failed,
        blocks === 0, Message[QECFaultTolerant::location, circuit]; $Failed,
        True,
            assembled = ftAssemble[a, circuit, blocks, order];
            If[ assembled === $Failed,
                $Failed,
                QECFaultTolerant[<|
                    "Code" -> a,
                    "Circuit" -> circuit,
                    "Blocks" -> blocks,
                    "Order" -> order,
                    (* data blocks first, then one EC workspace per block *)
                    "Qubits" -> 2 blocks a["Qubits"],
                    assembled
                |>]
            ]
    ]
]


(* ---- properties ---- *)

$faultTolerantProperties = {
    "Code", "Circuit", "Blocks", "LogicalQubits", "Qubits", "Instructions",
    "InstructionCount", "Depth", "Locations", "CircuitLocations", "Measurements",
    "Readouts", "Gadgets", "GadgetCounts", "GadgetInstructions", "GateCounts",
    "QuantumCircuitOperator", "Diagram", "SizeOverhead",
    "QubitOverhead", "DepthOverhead", "Overheads", "OpenAssumptions", "Properties"
};

faultTolerantData[QECFaultTolerant[a_Association]] := a

QECFaultTolerant[_Association]["Properties"] := $faultTolerantProperties

QECFaultTolerant[a_Association][prop : ("Circuit" | "Blocks" | "Qubits" |
    "Instructions" | "Readouts" | "Measurements")] := a[prop]

QECFaultTolerant[a_Association]["Code"] := QECCode[a["Code"]]
QECFaultTolerant[a_Association]["LogicalQubits"] := a["Blocks"]
QECFaultTolerant[a_Association]["InstructionCount"] := Length[a["Instructions"]]
QECFaultTolerant[a_Association]["Depth"] := gadgetDepth[a["Instructions"], a["Qubits"]]
QECFaultTolerant[a_Association]["GateCounts"] := gadgetGateCounts[a["Instructions"]]
QECFaultTolerant[a_Association]["QuantumCircuitOperator"] :=
    gadgetCircuitOperator[a["Instructions"]]
QECFaultTolerant[a_Association]["Diagram"] := gadgetDiagram[a["Instructions"]]

QECFaultTolerant[a_Association]["Gadgets"] := Dataset[a["Gadgets"]]

(* The instructions of the gadgets of one kind, in order.  The table is in emission
   order and each row carries its own length, so the slices are cumulative sums --
   a gadget never has to be recognised by looking at its gates. *)
QECFaultTolerant[a_Association]["GadgetInstructions", kind_] := Module[
    {rows = a["Gadgets"], kinds = Flatten[{kind}], edges},
    edges = Prepend[Accumulate[#["Instructions"] & /@ rows], 0];
    Catenate @ Table[
        If[ MemberQ[kinds, rows[[j]]["Gadget"]],
            Take[a["Instructions"], {edges[[j]] + 1, edges[[j + 1]]}],
            {}
        ],
        {j, Length[rows]}
    ]
]
QECFaultTolerant[a_Association]["GadgetCounts"] := Counts[#["Gadget"] & /@ a["Gadgets"]]

(* Locations, both sides, counted the same way: instructions plus waits. *)
QECFaultTolerant[a_Association]["Locations"] := ftLocations[a["Instructions"], a["Qubits"]]
QECFaultTolerant[a_Association]["CircuitLocations"] := ftLocations[a["Circuit"], a["Blocks"]]

QECFaultTolerant[a_Association]["SizeOverhead"] :=
    QECFaultTolerant[a]["Locations"] / QECFaultTolerant[a]["CircuitLocations"]

QECFaultTolerant[a_Association]["QubitOverhead"] := a["Qubits"] / a["Blocks"]

QECFaultTolerant[a_Association]["DepthOverhead"] :=
    QECFaultTolerant[a]["Depth"] / Length[circuitSchedule[a["Circuit"], a["Blocks"]]]

QECFaultTolerant[a_Association]["Overheads"] := <|
    "Size" -> QECFaultTolerant[a]["SizeOverhead"],
    "Qubits" -> QECFaultTolerant[a]["QubitOverhead"],
    "Depth" -> QECFaultTolerant[a]["DepthOverhead"]
|>

QECFaultTolerant[a_Association]["OpenAssumptions"] := {
    "Preparation of an encoded |0> uses the code's own encoder, which is not fault \
tolerant (Got26 Procedure 6.6); so do the ancillas of every error correction \
gadget. Fault-tolerant preparation is chapter 13 work.",
    "The gate gadgets are transversal, so the protocol is Clifford only. A universal \
set (Def 10.5) needs the constructions of chapter 13.",
    "Corrections are carried in the Pauli frame (sec. 12.5.2) rather than applied, so \
the emitted circuit has no classically conditioned instruction."
}

QECFaultTolerant[a_Association][prop_String] :=
    (Message[QECFaultTolerant::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECFaultTolerant /: MakeBoxes[
    obj : QECFaultTolerant[a_Association] /; KeyExistsQ[a, "Gadgets"],
    form : (StandardForm | TraditionalForm)
] := BoxForm`ArrangeSummaryBox[
    QECFaultTolerant,
    obj,
    None,
    {
        BoxForm`SummaryItem[{"Logical qubits: ", a["Blocks"]}],
        BoxForm`SummaryItem[{"Qubits: ", a["Qubits"]}],
        BoxForm`SummaryItem[{"Gadgets: ", Length[a["Gadgets"]]}]
    },
    {
        BoxForm`SummaryItem[{"Code: ", QECCode[a["Code"]]["Parameters"]}],
        BoxForm`SummaryItem[{"Locations: ", obj["Locations"]}],
        BoxForm`SummaryItem[{"Size overhead: ", N[obj["SizeOverhead"]]}]
    },
    form,
    "Interpretable" -> False
]
