(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECExtendedRectangle]

PackageScope[exRecLiveWaits]
PackageScope[exRecLayout]
PackageScope[exRecLocations]
PackageScope[concatenationThreshold]
PackageScope[preparationAttemptLocations]


(* ============================================================================ *)
(* Extended rectangles, and the bound they give (Got26 sec. 14.2, 14.5, 14.7).  *)
(*                                                                              *)
(* The threshold theorem is proved one piece of FT(C) at a time.  An EXTENDED   *)
(* RECTANGLE is a gadget together with the error-correction gadgets right       *)
(* before and right after it (Def 14.2): a gate or wait gadget with its leading *)
(* and trailing EC, a preparation with the EC after it, a measurement with the  *)
(* EC before it.  It is GOOD if it holds at most t faults and BAD otherwise     *)
(* (Def 14.4), and a good exRec is CORRECT (Theorem 14.1): it does to the       *)
(* encoded state what the ideal location does.  So a computation can only go    *)
(* wrong where an exRec is bad, and the probability of that is bounded by a     *)
(* union over the sets of t + 1 locations it contains (eq. 14.29):              *)
(*                                                                              *)
(*     Prob(bad) <= Binomial[A, t + 1] p^(t + 1),   A = the locations it has,   *)
(*                                                                              *)
(* divided by (1 - Q)^k when k of its ancillas are post-selected, Q the chance  *)
(* of at least one fault in one preparation attempt (eq. 14.33).                *)
(*                                                                              *)
(* WHAT THIS FILE COUNTS IS THIS LAYER'S OWN CIRCUIT.  The exRec is built from  *)
(* the same gadgets FT(C) uses -- the Steane EC with verified ancillas of       *)
(* ErrorCorrection.wl, the transversal gate words of FaultTolerant.wl -- and A  *)
(* is read off it, location by location: every instruction, and every WAIT.    *)
(*                                                                              *)
(* A wait is a qubit that holds a state through a time step in which nothing     *)
(* touches it: after its first use and before its last, and not between a       *)
(* measurement and the next reset, when it holds nothing, nor between a reset   *)
(* and the first gate, since a reset can be done just in time.  The book's      *)
(* convention -- its EC gadget waits 34 times, not once per idle qubit per step  *)
(* -- and it is stricter than circuitIdleSlots, which charges every untouched   *)
(* qubit every step because the detector model wants an idling channel          *)
(* everywhere.  The two are different questions and are kept apart.            *)
(*                                                                              *)
(* To keep the data from waiting on ancillas, each EC gadget of the exRec has a *)
(* workspace of its own, so the trailing EC's ancillas are prepared while the    *)
(* leading EC runs, which is the arrangement the book's counts assume ("the     *)
(* data block doesn't have to wait").                                           *)
(*                                                                              *)
(* THE BOOK'S NUMBERS ARE FOR THE BOOK'S CIRCUIT.  Sec. 14.5.3 finds A = 371    *)
(* for a one-qubit gate exRec of the 7-qubit code and 735 for CNOT, Q = 63 p, *)
(* and p_T > 3.6 x 10^-6 from the CNOT exRec.  This layer's circuit is the same  *)
(* protocol laid out by an ASAP schedule rather than by hand, so its counts     *)
(* differ, and the comparison is made explicitly rather than by tuning.          *)
(*                                                                              *)
(* References: Got26 Def 14.2-14.4, Theorem 14.1, sec. 14.5.1-14.5.3 (eqs.       *)
(* 14.28-14.41), sec. 14.7.4 (eqs. 14.62-14.63, Def 14.10); fig. 13.3, 14.8a.    *)
(* ============================================================================ *)

QECExtendedRectangle::usage = "QECExtendedRectangle[code, location] gives the extended rectangle of a fault-tolerant protocol for a CSS code: the gadget for the location together with the error-correction gadgets around it.\nThe location is \"Preparation\", \"Measurement\", \"Wait\", a one-qubit Clifford such as \"H\", or \"CNOT\".\nThe options \"Ancilla\" (\"Verified\" or \"Encoder\") and \"Order\" are those of QECErrorCorrection.\nrect[\"BadProbabilityBound\", p] is the bound of Got26 eq. 14.29 with the post-selection correction of eq. 14.33.";

QECExtendedRectangle::css = "Extended rectangles are built here from Steane error correction, which needs a CSS code with one logical qubit.";
QECExtendedRectangle::location = "`1` is not a location with a gadget for this code; use \"Preparation\", \"Measurement\", \"Wait\", a transversal one-qubit Clifford, or \"CNOT\".";
QECExtendedRectangle::noprop = "`1` is not a property of QECExtendedRectangle. Use rect[\"Properties\"] for the list.";

Options[QECExtendedRectangle] = {"Ancilla" -> "Verified", "Order" -> "BitFlipFirst"};


(* ---- waits that hold a state ---- *)

(* Every (anchor, layer, qubit) at which a qubit holding a state is left alone for
   a layer: between two of its instructions, the first not a measurement.  The
   anchor is the instruction before, which is where a fault in that wait acts.

   Nor a reset: a fresh |0> can always be prepared just before it is first used,
   so the gap after a reset is not a wait -- the ASAP schedule puts every reset in
   the first layer only because it schedules everything as early as it can.  This
   is "designed to minimize waiting times", the book's fig. 14.8a, done by the
   counting rather than by redrawing the circuit. *)
exRecLiveWaits[instr_List, nq_Integer] := Module[{layerOf, touches},
    layerOf = circuitInstructionLayers[instr, nq];
    touches = GroupBy[
        Catenate @ Table[{q, layerOf[[i]], i, First[instr[[i]]]}, {i, Length[instr]}, {q, Rest[instr[[i]]]}],
        First -> Rest
    ];
    Catenate @ KeyValueMap[
        Function[{q, ts},
            Catenate @ Table[
                With[{now = ts[[j]], next = ts[[j + 1]]},
                    If[ MemberQ[$circuitMeasureOps, now[[3]]] || now[[3]] === "R",
                        {},
                        Table[{now[[2]], L, q}, {L, now[[1]] + 1, next[[1]] - 1}]
                    ]
                ],
                {j, Length[ts] - 1}
            ]
        ],
        SortBy[#, First] & /@ touches
    ]
]


(* ---- the layout ---- *)

(* Blocks first, then two workspaces per block, leading EC in the first and
   trailing EC in the second.  Each part of the exRec is kept with its name, so
   the leading EC, the gadget and the trailing EC can be told apart without
   reading the instructions. *)
exRecLayout[a_Association, location_, mode_String, order_String] := Module[
    {n = a["Qubits"], w = ancillaWidth[a, mode], blocks, ws, ec, ecChecks, parts, extraWaits = {}, gadget},

    blocks = If[MemberQ[$circuitTwoQubitOps, location] || location === "T", 2, 1];
    ws[b_, j_] := blocks n + (2 (b - 1) + (j - 1)) w;
    ec[b_, j_] := ftRelocate[Catenate[steanePieces[a, n, order, mode][[All, 2]]], n, n, (b - 1) n, ws[b, j]];
    ecChecks = Length[piecesHeraldChecks[steanePieces[a, n, order, mode]]];

    parts = Switch[location,
        "T",
            (* data block 1..n, magic block n+1..2n; then the leading EC's workspace,
               the magic block's EC workspace, the cat and its check, the trailing
               EC's workspace *)
            Module[{wl = blocks n, wm = blocks n + w, cat = blocks n + 2 w, chk = blocks n + 2 w + n + 1, wt = blocks n + 2 w + n + 1, ecAt},
                ecAt[base_, wsb_] := ftRelocate[Catenate[steanePieces[a, n, order, mode][[All, 2]]], n, n, base, wsb];
                {"LeadingEC" -> ecAt[0, wl],
                 "MagicPreparation" -> Join[Table[{"R", n + q}, {q, n}], magicStateInstructions[a, n, "Tdg"]],
                 "MagicEC" -> ecAt[n, wm],
                 "MagicTest" -> magicTestInstructions[a, n, cat, chk],
                 (* fig. 13.7: CNOT from the magic block into the data, the data read in
                    Z, and on outcome 1 the correction -- logical X, kept in the frame,
                    then a transversal S, which the state-vector check found to be the
                    one that completes the gate.  Counted always, as the book does. *)
                 "Injection" -> Join[Table[{"CNOT", n + q, q}, {q, n}], Table[{"M", q}, {q, n}],
                    Table[{"S", n + q}, {q, n}]],
                 "TrailingEC" -> ecAt[n, wt + 1]}
            ],
        "Preparation",
            {"Preparation" -> If[ mode === "Verified",
                verifiedZeroInstructions[a, 0, ws[1, 1]],
                Join[Table[{"R", q}, {q, n}], codeEncodedZeroInstructions[a, 0]]],
             "TrailingEC" -> ec[1, 2]},
        "Measurement",
            {"LeadingEC" -> ec[1, 1], "Measurement" -> ftMeasurement[a, 1]},
        "Wait",
            extraWaits = Range[n];
            {"LeadingEC" -> ec[1, 1], "Wait" -> {}, "TrailingEC" -> ec[1, 2]},
        _?(MemberQ[$circuitTwoQubitOps, #] &),
            gadget = ftTwoQubitInstructions[a, location, 1, 2];
            If[gadget === $Failed, Return[$Failed]];
            {"LeadingEC" -> Join[ec[1, 1], ec[2, 1]], "Gate" -> gadget, "TrailingEC" -> Join[ec[1, 2], ec[2, 2]]},
        _String,
            gadget = ftGateInstructions[a, location, 1];
            If[gadget === $Failed, Return[$Failed]];
            {"LeadingEC" -> ec[1, 1], "Gate" -> gadget, "TrailingEC" -> ec[1, 2]},
        _, Return[$Failed]
    ];

    <|
        "Parts" -> Association[parts],
        "Instructions" -> Catenate[Values[parts]],
        "Qubits" -> If[location === "T", blocks n + 3 w + n + 2, blocks n + 2 blocks w],
        "Blocks" -> blocks,
        (* the wait gadget emits nothing; its n waits are located on the data at the
           seam between the two ECs *)
        "GadgetWaits" -> extraWaits,
        "HeraldChecks" -> exRecHeraldChecks[a, Association[parts], blocks, mode, order],
        (* post-selected ancillas: those of every EC gadget, and the preparation
           gadget's own when it is verified *)
        "Verifications" -> ecChecks If[location === "T", 1, blocks] Count[Keys[parts], "LeadingEC" | "TrailingEC"] +
            If[location === "Preparation" && mode === "Verified", 1, 0]
    |>
]

(* The herald checks of the whole exRec, as {positions among its heralds, rows}:
   the verified preparation's own, then each EC gadget's, in emission order. *)
exRecHeraldChecks[a_Association, parts_Association, blocks_Integer, mode_String, order_String] := Module[
    {n = a["Qubits"], h = 0, out = {}, ecChecks, ecHeralds},
    ecChecks = piecesHeraldChecks[steanePieces[a, n, order, mode]];
    ecHeralds = Count[Catenate[steanePieces[a, n, order, mode][[All, 2]]], {"MH", _}];
    Do[
        Which[
            part === "Preparation" && mode === "Verified",
                AppendTo[out, {h + Range[n], verificationChecks[a, "Z"]}]; h += n,
            MemberQ[{"LeadingEC", "TrailingEC"}, part],
                Do[out = Join[out, {h + #[[1]], #[[2]]} & /@ ecChecks]; h += ecHeralds, {blocks}]
        ],
        {part, Keys[parts]}
    ];
    out
]

(* The locations of the exRec, each with the faults it can carry: the detector
   model's enumeration for the instructions, a single-qubit Pauli for each wait. *)
exRecLocations[layout_Association] := Module[
    {instr = layout["Instructions"], nq = layout["Qubits"], rates, ops, waits, seam},
    rates = <|"OneQubit" -> 1, "TwoQubit" -> 1, "Measurement" -> 1, "Reset" -> 1, "Idle" -> 1|>;
    ops = KeyValueMap[<|"Type" -> #1[[1]], "Index" -> #1[[2]], "Faults" -> #2|> &,
        GroupBy[circuitMechanisms[instr, rates], #["Location"] &, #[[All, "Faults"]] &]];
    waits = Function[slot, <|"Type" -> "Wait", "Index" -> slot[[1]],
        "Faults" -> Table[{{slot[[1]], slot[[3]], p}}, {p, $oneQubitPaulis}]|>] /@ exRecLiveWaits[instr, nq];
    seam = Length[Lookup[layout["Parts"], "LeadingEC", {}]];
    (* the pi/8 gates, which the detector model has no mechanism for: a location each *)
    ops = Join[ops, Function[i, <|"Type" -> "NonClifford", "Index" -> i,
        "Faults" -> Table[{{i, instr[[i, 2]], p}}, {p, $oneQubitPaulis}]|>] /@
        Flatten[Position[instr, {g_ /; MemberQ[$nonCliffordOps, g], _}, {1}, Heads -> False]]];
    Join[ops, waits,
        Function[q, <|"Type" -> "Wait", "Index" -> seam, "Gadget" -> True,
            "Faults" -> Table[{{seam, q, p}}, {p, $oneQubitPaulis}]|>] /@ layout["GadgetWaits"]]
]

$exRecTypeNames = <|"Reset" -> "Preparation", "Measurement" -> "Measurement", "OneQubit" -> "OneQubitGate",
    "TwoQubit" -> "TwoQubitGate", "Wait" -> "Wait", "NonClifford" -> "PiOver8Gate"|>;

(* The magic state's whole preparation -- encoder, its EC, the test -- is one
   post-selected attempt (eq. 14.41 divides by 1 - 282 p for it), counted on its
   own circuit the same way. *)
magicAttemptLocations[layout_Association] := With[
    {instr = Catenate[Lookup[layout["Parts"], {"MagicPreparation", "MagicEC", "MagicTest"}, {}]]},
    If[instr === {}, 0, Length[instr] + Length[exRecLiveWaits[instr, layout["Qubits"]]]]
]

(* B in Q = B p: the locations of ONE attempt at a post-selected ancilla, counted
   the same way, on that attempt's circuit alone. *)
preparationAttemptLocations[a_Association, "Verified"] := With[{n = a["Qubits"]},
    Max @ Table[
        With[{instr = Catenate[pieces[[All, 2]]]},
            Length[instr] + Length[exRecLiveWaits[instr, 2 n]]
        ],
        {pieces, {verifiedZeroPieces[a, 0, n],
            If[hadamardIsLogicalQ[a], verifiedZeroPieces[a, 0, n], verifiedPlusPieces[a, 0, n]]}}
    ]
]
preparationAttemptLocations[_Association, "Encoder"] := 0


(* ---- construction ---- *)

QECExtendedRectangle[code_QECCode, location_, opts : OptionsPattern[]] := Module[
    {a = First[code], mode = OptionValue["Ancilla"], order = OptionValue["Order"], layout, locs},
    Which[
        ! codeCSSQ[a] || codeLogicalQubits[a] =!= 1, Message[QECExtendedRectangle::css]; Return[$Failed],
        ! MemberQ[$ancillaModes, mode], Message[QECErrorCorrection::ancilla, mode, $ancillaModes]; Return[$Failed],
        ! MemberQ[{"BitFlipFirst", "PhaseFirst"}, order], Message[QECErrorCorrection::order, order]; Return[$Failed]
    ];
    layout = Quiet[exRecLayout[a, location, mode, order]];
    If[ layout === $Failed, Message[QECExtendedRectangle::location, location]; Return[$Failed]];
    locs = exRecLocations[layout];
    QECExtendedRectangle[Join[layout, <|
        "Code" -> a,
        "Location" -> location,
        "Ancilla" -> mode,
        "Order" -> order,
        "LocationList" -> locs,
        "PreparationLocations" -> preparationAttemptLocations[a, mode],
        "MagicPreparationLocations" -> magicAttemptLocations[layout],
        "Correctable" -> Quotient[codeDistance[a] - 1, 2]
    |>]]
]


(* ---- the bound ---- *)

(* Eq. 14.29 with eq. 14.33: every set of t + 1 locations, each at rate p, over
   the acceptance of each post-selected ancilla, at least 1 - Q. *)
exRecBadBound[d_Association, p_] := With[{A = Length[d["LocationList"]], t = d["Correctable"]},
    Binomial[A, t + 1] p^(t + 1) / (1 - d["PreparationLocations"] p)^d["Verifications"] / magicFactor[d, p]
]

magicFactor[d_Association, p_] := With[{qm = Lookup[d, "MagicPreparationLocations", 0]}, If[qm > 0, 1 - qm p, 1]]

(* The fixed point of a level of concatenation: a level-1 exRec fails with
   probability at most eps(p), and the threshold is where eps(p) = p.  With
   eps = C p^(t+1) / (1 - Q p)^k that is the smallest root of C p^t = (1 - Q p)^k
   in (0, 1/Q); eq. 14.63 takes the CNOT exRec's and finds 3.6 x 10^-6. *)
concatenationThreshold[{A_Integer, Q_Integer, k_Integer}, t_Integer] := concatenationThreshold[{A, Q, k, 0}, t]

concatenationThreshold[{A_Integer, Q_Integer, k_Integer, Qm_Integer}, t_Integer] := Module[{p, roots},
    roots = p /. Solve[Binomial[A, t + 1] p^t == (1 - Q p)^k If[Qm > 0, 1 - Qm p, 1] &&
        0 < p < If[Max[Q, Qm] > 0, 1 / Max[Q, Qm], 1], p, Reals];
    If[! FreeQ[roots, p] || roots === {}, Missing["NoFixedPoint"], First[SortBy[roots, N]]]
]


(* With the malignant pairs counted (t = 1): eq. 14.32, every triple still
   counted as fatal, over the same acceptance.  Only pairs are counted, so for a
   code that corrects more than one error there is nothing to replace. *)
nonCliffordQ[d_Association] := MemberQ[First /@ d["Instructions"], Alternatives @@ $nonCliffordOps]

exRecMalignantBound[d_Association, p_] /; nonCliffordQ[d] := Missing["NonClifford", "A Pauli frame cannot follow a fault through a pi/8 gate."]
malignantThreshold[d_Association] /; nonCliffordQ[d] := Missing["NonClifford", "A Pauli frame cannot follow a fault through a pi/8 gate."]

exRecMalignantBound[d_Association, p_] := If[ d["Correctable"] =!= 1,
    Missing["NotCounted", "Malignant sets are counted for t = 1 only."],
    With[{A = Length[d["LocationList"]], m = exRecMalignantCounts[d]["Pairs"]},
        (m p^2 + Binomial[A, 3] p^3) / (1 - d["PreparationLocations"] p)^d["Verifications"]
    ]
]

malignantThreshold[d_Association] := If[ d["Correctable"] =!= 1,
    Missing["NotCounted", "Malignant sets are counted for t = 1 only."],
    Module[{p, roots, A = Length[d["LocationList"]], m = exRecMalignantCounts[d]["Pairs"], Q = d["PreparationLocations"], k = d["Verifications"]},
        roots = p /. Solve[m p + Binomial[A, 3] p^2 == (1 - Q p)^k && 0 < p < If[Q > 0, 1 / Q, 1], p, Reals];
        If[! FreeQ[roots, p] || roots === {}, Missing["NoFixedPoint"], First[SortBy[roots, N]]]
    ]
]


(* ---- properties ---- *)

$exRecProperties = {
    "Location", "Code", "Ancilla", "Order", "Instructions", "Qubits", "Blocks", "Parts",
    "Locations", "LocationCounts", "Verifications", "PreparationLocations", "Correctable",
    "BadProbabilityBound", "ThresholdBound", "RectangleLocations", "MagicPreparationLocations",
    "MalignantSingles", "MalignantPairs", "MalignantPairsByType", "MalignantProbabilityBound", "MalignantThresholdBound",
    "Properties"
};

QECExtendedRectangle[_Association]["Properties"] := $exRecProperties

QECExtendedRectangle[d_Association][prop : "Location" | "Ancilla" | "Order" | "Instructions" | "Qubits" |
    "Blocks" | "Verifications" | "PreparationLocations" | "Correctable" | "MagicPreparationLocations"] := d[prop]
QECExtendedRectangle[d_Association]["Code"] := QECCode[d["Code"]]
QECExtendedRectangle[d_Association]["Parts"] := Length /@ d["Parts"]
QECExtendedRectangle[d_Association]["Locations"] := Length[d["LocationList"]]
QECExtendedRectangle[d_Association]["LocationCounts"] :=
    KeySort @ Counts[$exRecTypeNames[#["Type"]] & /@ d["LocationList"]]
QECExtendedRectangle[d_Association]["BadProbabilityBound", p_] := exRecBadBound[d, p]
QECExtendedRectangle[d_Association]["ThresholdBound"] := concatenationThreshold[
    {Length[d["LocationList"]], d["PreparationLocations"], d["Verifications"], Lookup[d, "MagicPreparationLocations", 0]}, d["Correctable"]]

(* A rectangle is the exRec without its leading EC gadgets: the unit FT(C) is cut
   into, and so the unit of overhead (sec. 14.7.5). *)
QECExtendedRectangle[d_Association]["RectangleLocations"] := With[
    {a = d["Code"], n = d["Code"]["Qubits"], leading = If[KeyExistsQ[d["Parts"], "LeadingEC"], If[MemberQ[$circuitTwoQubitOps, d["Location"]], 2, 1], 0]},
    Length[d["LocationList"]] - leading ecGadgetLocations[a, d["Ancilla"], d["Order"]]
]

(* one EC gadget on its own, counted the same way *)
ecGadgetLocations[a_Association, mode_String, order_String] := With[
    {instr = Catenate[steanePieces[a, a["Qubits"], order, mode][[All, 2]]]},
    Length[instr] + Length[exRecLiveWaits[instr, a["Qubits"] + ancillaWidth[a, mode]]]
]

(* Theorem 14.1 in this circuit: no single location is malignant.  Zero is the
   check that the protocol really is fault tolerant; with the bare encoder it is
   not, and the preparation faults show up here. *)
QECExtendedRectangle[d_Association][prop : "MalignantSingles" | "MalignantPairs" | "MalignantPairsByType"] /; nonCliffordQ[d] :=
    Missing["NonClifford", "A Pauli frame cannot follow a fault through a pi/8 gate."]
QECExtendedRectangle[d_Association]["MalignantSingles"] := exRecMalignantCounts[d]["Singles"]
QECExtendedRectangle[d_Association]["MalignantPairs"] := exRecMalignantCounts[d]["Pairs"]
QECExtendedRectangle[d_Association]["MalignantPairsByType"] :=
    KeyMap[Map[$exRecTypeNames], exRecMalignantCounts[d]["PairsByType"]]
QECExtendedRectangle[d_Association]["MalignantProbabilityBound", p_] := exRecMalignantBound[d, p]
QECExtendedRectangle[d_Association]["MalignantThresholdBound"] := malignantThreshold[d]

QECExtendedRectangle[_Association][prop_String, ___] :=
    (Message[QECExtendedRectangle::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECExtendedRectangle /: MakeBoxes[obj : QECExtendedRectangle[d_Association] /; KeyExistsQ[d, "LocationList"],
    form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECExtendedRectangle,
        obj,
        None,
        {
            BoxForm`SummaryItem[{"Location: ", d["Location"]}],
            BoxForm`SummaryItem[{"Locations: ", Length[d["LocationList"]]}],
            BoxForm`SummaryItem[{"Code: ", QECCode[d["Code"]]["Parameters"]}]
        },
        {
            BoxForm`SummaryItem[{"Ancilla: ", d["Ancilla"]}],
            BoxForm`SummaryItem[{"Post-selections: ", d["Verifications"]}],
            BoxForm`SummaryItem[{"Parts: ", Length /@ d["Parts"]}]
        },
        form,
        "Interpretable" -> False
    ]
