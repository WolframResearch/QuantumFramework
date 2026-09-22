(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECErrorCorrection]

PackageScope[codeLogicalWires]
PackageScope[codeEncodedZeroInstructions]
PackageScope[codeEncodedPlusInstructions]
PackageScope[steaneBitFlipInstructions]
PackageScope[steanePhaseInstructions]
PackageScope[steaneInstructions]
PackageScope[steaneRegions]
PackageScope[steaneResidualWeights]
PackageScope[correctionData]


(* ============================================================================ *)
(* Steane error correction.                                                     *)
(*                                                                              *)
(* Measurement.wl leaves a slot: repetition handles faults in the ancilla, but  *)
(* a data error anticommuting with P defeats every repetition alike, so Got26   *)
(* Theorem 12.1 needs an FTEC sub-gadget between them.  This is that gadget.    *)
(*                                                                              *)
(* WHY STEANE AND NOT SHOR.  Got26 sec. 12.3 opens on the cost argument: Shor   *)
(* EC "involves a lot of locations.  Lots of locations means lots of            *)
(* opportunities for errors, which eventually will translate into a lousy       *)
(* threshold."  Steane EC moves that work into preparing the ancilla, and the   *)
(* asymmetry that licenses it is worth stating because it is the whole idea of  *)
(* the chapter: an ancilla can be CHECKED and a data block cannot, "the reason  *)
(* is that we know the precise ancilla state that we are trying to create, but  *)
(* the state of the logical qubits somewhere in the middle of a long            *)
(* computation is unknown".                                                     *)
(*                                                                              *)
(* CSS ONLY, and for a reason rather than by restriction: the construction runs *)
(* on transversal CNOT being the logical CNOT, which is a property of CSS codes.*)
(*                                                                              *)
(* THE BIT-FLIP HALF (fig. 12.8a).  Ancilla block in encoded |+>, transversal   *)
(* CNOT from data to ancilla, transversal measurement of the ancilla.  On the   *)
(* encoded states CNOT|psi>|+> = |psi>|+>, so nothing happens logically -- but  *)
(* the gate still propagates errors, and bit flips in the data land on the same *)
(* positions of the ancilla, where measuring is harmless.                       *)
(*                                                                              *)
(* THE PHASE HALF (fig. 12.8b).  Ancilla block in encoded |0>, transversal CNOT *)
(* the other way, ancilla as control, so phase errors flow data to ancilla; a   *)
(* transversal Hadamard then exchanges X and Z errors, and the transversal      *)
(* measurement reads them.                                                      *)
(*                                                                              *)
(* WHAT GOES THE WRONG WAY, and it is not a defect but the price: propagation   *)
(* is symmetric.  While bit flips cross into the ancilla, phase errors in the   *)
(* ancilla cross into the data, and in the phase half the roles swap.  So the   *)
(* ancilla's own quality is the thing that matters, which is why its            *)
(* preparation is the hard part and why Got26 defers it to chapter 13.          *)
(*                                                                              *)
(* NO REPETITION IS NEEDED (sec. 12.3.3), and this is why Steane EC is the      *)
(* right gadget for Measurement.wl's slot.  The errors are ADDITIVE: the        *)
(* measured classical error is e + f + g, for e the true data error, f the      *)
(* ancilla's and g the measurement's, so correcting by e + f + g leaves f + g   *)
(* behind -- "a single-qubit error in the ancilla can only produce a            *)
(* single-qubit error in the final state.  This is in contrast to Shor EC,      *)
(* where a single ancilla error changing one bit of the error syndrome could    *)
(* totally change the error we deduce."  A wrong syndrome bit here shifts the   *)
(* answer by one qubit; there it can shift it to a different coset entirely.    *)
(*                                                                              *)
(* THE OPEN ASSUMPTION.  The ancilla states are prepared here by the code's own *)
(* encoder, which is the non-fault-tolerant Procedure 6.6 construction, so this *)
(* gadget is fault tolerant GIVEN a clean ancilla and not otherwise.  Got26 is  *)
(* explicit that this is the hard part -- "How do we make the ancilla states    *)
(* and ensure they do not have too many errors?  That is the complicated part   *)
(* of Steane EC" -- and equally explicit that it is chapter 13 work, needed for *)
(* a full FT protocol whatever EC gadget is chosen.  "OpenAssumptions" names it *)
(* rather than leaving it implied.                                              *)
(*                                                                              *)
(* References: Got26 sec. 12.3.1 (the gadget, figs. 12.8-12.9), sec. 12.3.3     *)
(* (additivity, and why no repetition), sec. 11.4 (transversal measurement of a *)
(* CSS code), Theorem 12.5 (the ECRP), ch. 13 (preparing the ancilla states).   *)
(* ============================================================================ *)

QECErrorCorrection::usage = "QECErrorCorrection[code] gives the Steane error-correction gadget for a CSS code: a bit-flip half and a phase half, each an encoded ancilla block, a transversal CNOT and a transversal measurement.\nQECErrorCorrection[code, offset] places the ancilla block after the given qubit.\nThe option \"Order\" -> \"BitFlipFirst\" or \"PhaseFirst\" chooses which half runs first.\nec[prop] gives a property; ec[\"Properties\"] lists them.";

QECErrorCorrection::css = "Steane error correction needs a CSS code; this one is not CSS. Its generators must each be purely X-type or purely Z-type.";
QECErrorCorrection::order = "\"Order\" must be \"BitFlipFirst\" or \"PhaseFirst\"; got `1`.";
QECErrorCorrection::noprop = "`1` is not a property of QECErrorCorrection. Use ec[\"Properties\"] for the list.";


(* ---- the encoded ancilla states ---- *)

(* Which wires carry the unencoded logical qubits going into the encoder.  The
   standard form puts them last, but that is derived here rather than assumed: a
   wire is logical when putting the input in |+> there makes the encoded state a
   +1 eigenstate of the corresponding logical X.  Memoised, since it costs one
   engine run per candidate wire. *)
codeLogicalWires[a_Association] := codeLogicalWires[a] = With[
    {n = a["Qubits"], gates = codeEncodingGates[a], xbars = codeLogicalVectors[a]["X"]},
    Table[
        SelectFirst[
            Range[n],
            applyGates[
                Wolfram`QuantumFramework`PauliStabilizer[n],
                Join[{"H" -> #}, gates]
            ]["Expectation", QECPauliString[xbars[[j]]]] === 1 &
        ],
        {j, Length[xbars]}
    ]
]

(* |0_L>: the encoder run on |0...0>.  |+_L>: the same with an H on each logical
   wire first, which is a logical Hadamard on the input rather than on the block,
   so it needs no transversal-gate machinery. *)
codeEncodedZeroInstructions[a_Association, offset_Integer] :=
    engineGateInstructions[codeEncodingGates[a], offset]

codeEncodedPlusInstructions[a_Association, offset_Integer] := engineGateInstructions[
    Join[("H" -> #) & /@ codeLogicalWires[a], codeEncodingGates[a]],
    offset
]


(* ---- the two halves ---- *)

(* Figure 12.8a.  Data is the control, so bit flips copy into the ancilla; the
   ancilla is read in the Z basis and the classical syndrome of the result is the
   X-error syndrome of the data. *)
steaneBitFlipInstructions[a_Association, offset_Integer] := With[{n = a["Qubits"]},
    Join[
        Table[{"R", offset + q}, {q, n}],
        codeEncodedPlusInstructions[a, offset],
        Table[{"CNOT", q, offset + q}, {q, n}],
        Table[{"M", offset + q}, {q, n}]
    ]
]

(* Figure 12.8b.  Ancilla is the control, so phase errors flow the other way; the
   transversal Hadamard exchanges X and Z so the same Z-basis readout sees them. *)
steanePhaseInstructions[a_Association, offset_Integer] := With[{n = a["Qubits"]},
    Join[
        Table[{"R", offset + q}, {q, n}],
        codeEncodedZeroInstructions[a, offset],
        Table[{"CNOT", offset + q, q}, {q, n}],
        Table[{"H", offset + q}, {q, n}],
        Table[{"M", offset + q}, {q, n}]
    ]
]

(* Figure 12.9, without the final Pauli correction -- which is not applied here at
   all, for the reason Circuit.wl already relies on: the correction is carried in
   the Pauli frame (sec. 12.5.2) rather than executed as gates that could fail.

   Either order is fault tolerant (Theorem 12.5 does not care), but the order is
   not cosmetic: the errors left at the end tend to be of the type corrected
   FIRST, since more locations follow it. *)
steaneInstructions[a_Association, offset_Integer, order_String] := With[
    {bit = steaneBitFlipInstructions[a, offset], phase = steanePhaseInstructions[a, offset]},
    If[order === "BitFlipFirst", Join[bit, phase], Join[phase, bit]]
]


(* ---- the two regions, and the proof obligation ---- *)

(* The gadget splits cleanly in two, and the split is the whole honesty of the
   file, so it is computed from the construction rather than matched out of the
   emitted list.  Each half is (reset, encoder, interaction), and the interaction
   is the only part that touches a data qubit:

     bit-flip half   n resets, k logical-wire H, the encoder | n CNOT, n M
     phase half      n resets, the encoder                   | n CNOT, n H, n M

   Index 0 goes with the interaction: a fault there is the error already on the
   data when the gadget starts, which is the input an EC gadget exists to fix,
   not a fault of its own. *)
steaneRegions[a_Association, order_String] := Module[
    {n = a["Qubits"], k = Length[codeLogicalVectors[a]["X"]], enc, bitPrep, phasePrep, first, second},
    enc = Length[codeEncodingGates[a]];
    bitPrep = n + k + enc;
    phasePrep = n + enc;
    first = If[order === "BitFlipFirst", {bitPrep, 2 n}, {phasePrep, 3 n}];
    second = If[order === "BitFlipFirst", {phasePrep, 3 n}, {bitPrep, 2 n}];
    <|
        "Preparation" -> Join[
            Range[1, first[[1]]],
            Range[Total[first] + 1, Total[first] + second[[1]]]
        ],
        "Interaction" -> Join[
            {0},
            Range[first[[1]] + 1, Total[first]],
            Range[Total[first] + second[[1]] + 1, Total[first] + Total[second]]
        ]
    |>
]

(* The residual weight left on the data by every single fault in a given region,
   read the same way Measurement.wl reads it.

   Over the INTERACTION this is the claim of sec. 12.3.3, and it is the reason
   Steane EC needs no repetition: the classical error is e + f + g, so correcting
   by it leaves f + g, and "a single-qubit error in the ancilla can only produce a
   single-qubit error in the final state".  One fault, one data error.

   Over the PREPARATION it is not, and should not be: the encoder used here is the
   non-fault-tolerant one, and a single fault in it can leave a quarter of the
   block wrong.  Reporting that separately is what turns the open assumption into
   a measured quantity instead of a caveat in a comment. *)
steaneResidualWeights[instr_List, nq_Integer, n_Integer, region_List] :=
    DeleteDuplicates @ Flatten @ Table[
        With[{fr = framePropagate[instr, nq, {{i, q, pauli}}]["Frame"]},
            Total @ Map[Max, Take[Transpose[fr], n]]
        ],
        {i, region}, {q, nq}, {pauli, $oneQubitPaulis}
    ]


(* ---- construction ---- *)

Options[QECErrorCorrection] = {"Order" -> "BitFlipFirst"};

QECErrorCorrection[code_QECCode, opts : OptionsPattern[]] :=
    QECErrorCorrection[code, code["Qubits"], opts]

QECErrorCorrection[code_QECCode, offset_Integer, opts : OptionsPattern[]] := Module[
    {a = First[code], n = code["Qubits"], order},
    order = OptionValue["Order"];
    Which[
        ! codeCSSQ[a],
            Message[QECErrorCorrection::css]; $Failed,
        ! MemberQ[{"BitFlipFirst", "PhaseFirst"}, order],
            Message[QECErrorCorrection::order, order]; $Failed,
        True,
            QECErrorCorrection[<|
                "Code" -> a,
                "DataQubits" -> n,
                "Offset" -> offset,
                "AncillaQubits" -> Range[offset + 1, offset + n],
                "Order" -> order,
                "Qubits" -> offset + n,
                "Instructions" -> steaneInstructions[a, offset, order]
            |>]
    ]
]


(* ---- properties ---- *)

$correctionProperties = {
    "Instructions", "Code", "DataQubits", "AncillaQubits", "Offset", "Order",
    "Qubits", "Measurements", "Depth", "InstructionCount", "GateCounts",
    "QuantumCircuitOperator", "Diagram",
    "BitFlipInstructions", "PhaseInstructions", "TransversalQ", "RepetitionsNeeded",
    "Regions", "DataWeights", "MaxDataWeight", "PreparationDataWeights",
    "OpenAssumptions", "Properties"
};

correctionData[QECErrorCorrection[a_Association]] := a

QECErrorCorrection[_Association]["Properties"] := $correctionProperties

QECErrorCorrection[a_Association][prop : ("Instructions" | "DataQubits" |
    "AncillaQubits" | "Offset" | "Order" | "Qubits")] := a[prop]

QECErrorCorrection[a_Association]["Code"] := QECCode[a["Code"]]
QECErrorCorrection[a_Association]["Measurements"] := gadgetMeasurements[a["Instructions"]]
QECErrorCorrection[a_Association]["InstructionCount"] := Length[a["Instructions"]]
QECErrorCorrection[a_Association]["Depth"] := gadgetDepth[a["Instructions"], a["Qubits"]]
QECErrorCorrection[a_Association]["GateCounts"] := gadgetGateCounts[a["Instructions"]]
QECErrorCorrection[a_Association]["QuantumCircuitOperator"] :=
    gadgetCircuitOperator[a["Instructions"]]
QECErrorCorrection[a_Association]["Diagram"] := gadgetDiagram[a["Instructions"]]

QECErrorCorrection[a_Association]["BitFlipInstructions"] :=
    steaneBitFlipInstructions[a["Code"], a["Offset"]]
QECErrorCorrection[a_Association]["PhaseInstructions"] :=
    steanePhaseInstructions[a["Code"], a["Offset"]]

(* The data-ancilla interaction is one CNOT per position and nothing else.  The
   ancilla PREPARATION is not transversal and is not claimed to be; it touches no
   data qubit, which is what this checks. *)
QECErrorCorrection[a_Association]["TransversalQ"] := With[
    {n = a["DataQubits"], off = a["Offset"]},
    AllTrue[
        Cases[a["Instructions"], {op_ /; MemberQ[$circuitTwoQubitOps, op], x_, y_} :> {x, y}],
        Function[pair,
            With[{d = Select[pair, # <= n &], anc = Select[pair, # > off &]},
                Length[d] === 0 || (Length[d] === 1 && Length[anc] === 1 && First[anc] - off === First[d])
            ]
        ]
    ]
]

(* Section 12.3.3: none.  Kept as a property rather than a comment because it is
   the reason this gadget and not Shor EC fills Measurement.wl's slot. *)
QECErrorCorrection[a_Association]["RepetitionsNeeded"] := 1

QECErrorCorrection[a_Association]["Regions"] := steaneRegions[a["Code"], a["Order"]]

QECErrorCorrection[a_Association]["DataWeights"] := Sort @ steaneResidualWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"],
    steaneRegions[a["Code"], a["Order"]]["Interaction"]]

QECErrorCorrection[a_Association]["MaxDataWeight"] := Max @ steaneResidualWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"],
    steaneRegions[a["Code"], a["Order"]]["Interaction"]]

(* What the missing chapter-13 preparation costs, in the same units. *)
QECErrorCorrection[a_Association]["PreparationDataWeights"] := Sort @ steaneResidualWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"],
    steaneRegions[a["Code"], a["Order"]]["Preparation"]]

QECErrorCorrection[a_Association]["OpenAssumptions"] := {
    "The ancilla blocks are prepared by the code's own encoder, which is not fault \
tolerant (Got26 Procedure 6.6). Fault-tolerant preparation of encoded |0> and |+> \
is chapter 13 work, and this gadget is fault tolerant given a clean ancilla."
}

QECErrorCorrection[a_Association][prop_String] :=
    (Message[QECErrorCorrection::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECErrorCorrection /: MakeBoxes[
    obj : QECErrorCorrection[a_Association] /; KeyExistsQ[a, "AncillaQubits"],
    form : (StandardForm | TraditionalForm)
] := BoxForm`ArrangeSummaryBox[
    QECErrorCorrection,
    obj,
    BarChart[Values[Counts[First /@ a["Instructions"]]],
        ChartLabels -> Keys[Counts[First /@ a["Instructions"]]],
        ImageSize -> {Automatic, 34}, Axes -> False,
        ChartStyle -> RGBColor[0.15, 0.5, 0.65]],
    {
        BoxForm`SummaryItem[{"Method: ", "Steane"}],
        BoxForm`SummaryItem[{"Data qubits: ", a["DataQubits"]}],
        BoxForm`SummaryItem[{"Order: ", a["Order"]}]
    },
    {
        BoxForm`SummaryItem[{"Ancilla qubits: ", a["AncillaQubits"]}],
        BoxForm`SummaryItem[{"Qubits: ", a["Qubits"]}],
        BoxForm`SummaryItem[{"Instructions: ", Length[a["Instructions"]]}]
    },
    form,
    "Interpretable" -> False
]
