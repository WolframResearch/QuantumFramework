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
PackageScope[verifiedZeroInstructions]
PackageScope[verifiedZeroPieces]
PackageScope[verifiedPlusPieces]
PackageScope[hadamardIsLogicalQ]
PackageScope[piecesHeraldChecks]
PackageScope[cssZeroInstructions]
PackageScope[cssPlusInstructions]
PackageScope[verifiedPlusInstructions]
PackageScope[verificationChecks]
PackageScope[ancillaWidth]
PackageScope[steanePieces]
PackageScope[gadgetLocationFaults]
PackageScope[acceptedQ]
PackageScope[steaneResidual]
PackageScope[cssCorrection]
PackageScope[cssCorrectionTable]
PackageScope[correctionResidualWeights]
PackageScope[$ancillaModes]


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
(* THE ANCILLA, TWO WAYS.  "Ancilla" -> "Verified" prepares each ancilla block  *)
(* fault-tolerantly (sec. 13.1.2-13.1.3, fig. 13.3): two copies are encoded by  *)
(* the non-fault-tolerant encoder, the one kept is compared against the other   *)
(* through a transversal CNOT, and the run is kept only if the comparison sees  *)
(* nothing.  "Ancilla" -> "Encoder", the default, is the bare encoder, and for  *)
(* it the paragraph below still holds.                                          *)
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
            ]["Expectation", pauliString[xbars[[j]]]] === 1 &
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


(* ---- verified ancillas (Got26 sec. 13.1.2-13.1.3, fig. 13.3) ---- *)

(* VERIFYING |0_L>.  Encode it twice, CNOT the kept copy into the checker
   transversally, and measure the checker in Z.  Bit flips on the kept copy are
   copied into the checker, so a bit flip that one fault in the kept copy's
   encoder spread over several qubits shows up in the checker's outcomes; a fault
   in the checker's encoder shows up there too.  Either way the run is discarded.
   What the comparison does NOT look at is phase errors, which flow the other
   way, from the checker into the kept copy.  For the 7-qubit code that is fine,
   and the book says why: the Hamming code is perfect, so every phase-error
   syndrome is that of a single phase error, and on |0_L> a logical phase error
   is no error at all.  For another code it may not be, which is why the claim is
   checked by enumeration ("ResidualWeights") rather than assumed.

   WHAT "SEES NOTHING" MEANS.  The checker's outcomes are a random codeword, not
   zeros, so the test is a set of parities: the checker's Z-type checks and its
   logical Z, read off the n outcomes.  The heralds are those n outcomes and
   "HeraldChecks" says which parities of them must vanish.

   VERIFYING |+_L>.  When the transversal H is the logical H -- the 7-qubit code,
   any self-dual CSS code -- the book's way is cheapest: verify a |0_L> and apply
   H to every qubit.  Otherwise the dual construction: encode |+_L> twice, CNOT
   from the checker into the kept copy so phase errors are copied out, and read
   the checker in X. *)

$ancillaModes = {"Encoder", "Verified"};

(* The verified gadget prepares all four of its ancilla blocks -- a kept copy and a
   checker for each half -- side by side, so that they are ready together and the
   data waits for none of them: four blocks of workspace, against one for the bare
   encoder, which reuses its block between the halves. *)
ancillaWidth[a_Association, "Encoder"] := a["Qubits"]
ancillaWidth[a_Association, "Verified"] := 4 a["Qubits"]

(* THE CSS ENCODER of fig. 14.8a, for any CSS code.  |0_L> is the uniform
   superposition over the X-type checks applied to |0...0>: |0...0> is already
   fixed by every Z-type check and by the logical Z's, and each X-type check is
   added by a Hadamard on one qubit of its support -- a pivot no other check uses,
   found by row-reducing the X part -- and CNOTs from that pivot onto the rest of
   the support.  For the 7-qubit code: 3 Hadamards and 9 CNOTs, the book's count,
   against the 25 gates of the general stabilizer encoder of Procedure 6.6.  |+_L>
   is the same construction for the dual: Hadamards on the non-pivots of the
   row-reduced Z part, CNOTs from them onto the pivot. *)
cssZeroInstructions[a_Association, offset_Integer] := With[{n = a["Qubits"]}, With[
    {xs = Select[a["CheckMatrix"], Total[#[[n + 1 ;; 2 n]]] === 0 &][[All, 1 ;; n]]},
    If[ xs === {}, {},
        With[{basis = gf2Basis[xs]},
            Join[
                Table[{"H", offset + piv}, {piv, basis[[2]]}],
                cnotLayers @ Catenate @ MapThread[
                    Function[{row, piv}, Table[{"CNOT", offset + piv, offset + j}, {j, DeleteCases[Flatten[Position[row, 1]], piv]}]],
                    basis
                ]
            ]
        ]
    ]
]]

cssPlusInstructions[a_Association, offset_Integer] := With[{n = a["Qubits"]}, With[
    {zs = Select[a["CheckMatrix"], Total[#[[1 ;; n]]] === 0 &][[All, n + 1 ;; 2 n]]},
    If[ zs === {}, Table[{"H", offset + q}, {q, n}],
        With[{basis = gf2Basis[zs]},
            Join[
                Table[{"H", offset + q}, {q, Complement[Range[n], basis[[2]]]}],
                cnotLayers @ Catenate @ MapThread[
                    Function[{row, piv}, Table[{"CNOT", offset + j, offset + piv}, {j, DeleteCases[Flatten[Position[row, 1]], piv]}]],
                    basis
                ]
            ]
        ]
    ]
]]

(* The CNOTs of the encoder all commute -- each copies a pivot that nothing else
   writes to -- so they can go in any order, and the order is chosen to pack them
   into few time steps: greedily, each layer takes every remaining CNOT whose two
   qubits the layer has not used yet.  For the 7-qubit code that is three layers of
   three, the book's fig. 14.8a, and fewer layers is fewer waits. *)
cnotLayers[cnots_List] := Module[{left = cnots, out = {}, used, layer},
    While[left =!= {},
        used = {}; layer = {};
        Do[
            If[ ! IntersectingQ[Rest[g], used], AppendTo[layer, g]; used = Join[used, Rest[g]]],
            {g, left}
        ];
        out = Join[out, layer];
        left = Complement[left, layer]
    ];
    out
]

hadamardIsLogicalQ[a_Association] := hadamardIsLogicalQ[a] =
    TrueQ[transversalValidQ[a, "H", 1]] && transversalName[transversalLogicalAction[a, "H", 1], 1] === "H"

(* The parities a clean comparison leaves at zero, as rows over the checker's n
   outcomes: for a Z readout, the supports of the Z-type checks and of the logical
   Z operators; for an X readout, the same with X. *)
verificationChecks[a_Association, "Z"] := With[{n = a["Qubits"], mat = a["CheckMatrix"]},
    Join[
        Select[mat, Total[#[[1 ;; n]]] === 0 &][[All, n + 1 ;; 2 n]],
        codeLogicalVectors[a]["Z"][[All, n + 1 ;; 2 n]]
    ]
]
verificationChecks[a_Association, "X"] := With[{n = a["Qubits"], mat = a["CheckMatrix"]},
    Join[
        Select[mat, Total[#[[n + 1 ;; 2 n]]] === 0 &][[All, 1 ;; n]],
        codeLogicalVectors[a]["X"][[All, 1 ;; n]]
    ]
]

(* Each piece is {region, instructions, checks}: the checks are the rows above,
   or None for a piece with no herald. *)
verifiedZeroPieces[a_Association, keep_Integer, check_Integer] := With[{n = a["Qubits"]},
    {
        {"Preparation", Join[
            Table[{"R", keep + q}, {q, n}], Table[{"R", check + q}, {q, n}],
            cssZeroInstructions[a, keep], cssZeroInstructions[a, check]], None},
        {"Verification", Join[
            Table[{"CNOT", keep + q, check + q}, {q, n}],
            Table[{"MH", check + q}, {q, n}]], verificationChecks[a, "Z"]}
    }
]

verifiedPlusPieces[a_Association, keep_Integer, check_Integer] := With[{n = a["Qubits"]},
    If[ hadamardIsLogicalQ[a],
        Append[verifiedZeroPieces[a, keep, check], {"Preparation", Table[{"H", keep + q}, {q, n}], None}],
        {
            {"Preparation", Join[
                Table[{"R", keep + q}, {q, n}], Table[{"R", check + q}, {q, n}],
                cssPlusInstructions[a, keep], cssPlusInstructions[a, check]], None},
            {"Verification", Join[
                Table[{"CNOT", check + q, keep + q}, {q, n}],
                Table[{"H", check + q}, {q, n}],
                Table[{"MH", check + q}, {q, n}]], verificationChecks[a, "X"]}
        }
    ]
]

verifiedZeroInstructions[a_Association, keep_Integer, check_Integer] := Catenate[verifiedZeroPieces[a, keep, check][[All, 2]]]
verifiedPlusInstructions[a_Association, keep_Integer, check_Integer] := Catenate[verifiedPlusPieces[a, keep, check][[All, 2]]]

(* The whole gadget as pieces, in order: the region each instruction belongs to
   and the herald checks are read off this, never matched out of the list. *)
steanePieces[a_Association, offset_Integer, order_String, "Encoder"] := With[{n = a["Qubits"]},
    With[{
        bit = {
            {"Preparation", Join[Table[{"R", offset + q}, {q, n}], codeEncodedPlusInstructions[a, offset]], None},
            {"Interaction", Join[Table[{"CNOT", q, offset + q}, {q, n}], Table[{"M", offset + q}, {q, n}]], None}},
        phase = {
            {"Preparation", Join[Table[{"R", offset + q}, {q, n}], codeEncodedZeroInstructions[a, offset]], None},
            {"Interaction", Join[Table[{"CNOT", offset + q, q}, {q, n}], Table[{"H", offset + q}, {q, n}],
                Table[{"M", offset + q}, {q, n}]], None}}},
        If[order === "BitFlipFirst", Join[bit, phase], Join[phase, bit]]
    ]
]

steanePieces[a_Association, offset_Integer, order_String, "Verified"] := With[{n = a["Qubits"], po = offset + 2 a["Qubits"]},
    With[{
        bit = Append[verifiedPlusPieces[a, offset, offset + n],
            {"Interaction", Join[Table[{"CNOT", q, offset + q}, {q, n}], Table[{"M", offset + q}, {q, n}]], None}],
        phase = Append[verifiedZeroPieces[a, po, po + n],
            {"Interaction", Join[Table[{"CNOT", po + q, q}, {q, n}], Table[{"H", po + q}, {q, n}],
                Table[{"M", po + q}, {q, n}]], None}]},
        If[order === "BitFlipFirst", Join[bit, phase], Join[phase, bit]]
    ]
]

(* Region -> instruction indices, and the herald checks as {herald positions,
   rows}, both from the pieces. *)
(* Index 0 goes with the interaction, as in steaneRegions: a fault there is the
   error already on the data when the gadget starts. *)
piecesRegions[pieces_List] := With[{edges = Prepend[Accumulate[Length /@ pieces[[All, 2]]], 0]},
    MapAt[Prepend[#, 0] &, Merge[Table[pieces[[j, 1]] -> Range[edges[[j]] + 1, edges[[j + 1]]], {j, Length[pieces]}], Catenate], "Interaction"]
]

piecesHeraldChecks[pieces_List] := Module[{h = 0, out = {}},
    Do[
        With[{k = Count[piece[[2]], {"MH", _}]},
            If[piece[[3]] =!= None, AppendTo[out, {Range[h + 1, h + k], piece[[3]]}]];
            h += k
        ],
        {piece, pieces}
    ];
    out
]

(* A run is kept when every herald check reads zero. *)
acceptedQ[heralds_List, checks_List] := AllTrue[checks, Function[c, AllTrue[Mod[c[[2]] . heralds[[c[[1]]]], 2], # === 0 &]]]


(* ---- what a single fault leaves behind, with the correction applied ---- *)

(* Every location of an instruction list and the faults each can carry: the
   detector model's own enumeration, with every rate set to one, grouped by
   location -- a reset, a measurement, a one- or two-qubit gate, or a wait. *)
gadgetLocationFaults[instr_List, nq_Integer] := With[
    {rates = <|"OneQubit" -> 1, "TwoQubit" -> 1, "Measurement" -> 1, "Reset" -> 1, "Idle" -> 1|>},
    GroupBy[Join[circuitMechanisms[instr, rates], idleMechanisms[instr, nq, rates]], #["Location"] &, #[[All, "Faults"]] &]
]

(* The residual a fault leaves on the data block once the gadget's own correction
   is applied: the frame on the data, times the X correction read from the
   bit-flip half's outcomes and the Z correction read from the phase half's, each
   the lightest error of its syndrome.  Its weight is taken modulo the stabilizer,
   since the data block carries an unknown logical state and nothing else can be
   divided out.  Missing when the heralds reject the run. *)
steaneResidual[a_Association, instr_List, nq_Integer, checks_List, order_String, faults_List] := Module[
    {n = a["Qubits"], run, rec, bitRec, phaseRec, x, z, cx, cz},
    run = framePropagate[instr, nq, faults];
    If[! acceptedQ[run["Heralds"], checks], Return[Missing["Rejected"]]];
    rec = run["Record"];
    {bitRec, phaseRec} = If[order === "BitFlipFirst", {rec[[1 ;; n]], rec[[n + 1 ;; 2 n]]}, {rec[[n + 1 ;; 2 n]], rec[[1 ;; n]]}];
    {x, z} = run["Frame"][[All, 1 ;; n]];
    cx = cssCorrection[a, "X", Mod[verificationChecks[a, "Z"][[1 ;; -1 - codeLogicalQubits[a]]] . bitRec, 2]];
    cz = cssCorrection[a, "Z", Mod[verificationChecks[a, "X"][[1 ;; -1 - codeLogicalQubits[a]]] . phaseRec, 2]];
    reducedWeight[a, BitXor[x, cx], BitXor[z, cz]]
]

(* The lightest X (or Z) error with a given syndrome under the Z-type (X-type)
   checks, found by weight, up to the weight where the table is complete. *)
cssCorrection[a_Association, letter_String, syn_List] := Lookup[cssCorrectionTable[a, letter], Key[syn], ConstantArray[0, a["Qubits"]]]

cssCorrectionTable[a_Association, letter_String] := cssCorrectionTable[a, letter] = Module[
    {n = a["Qubits"], checks = verificationChecks[a, If[letter === "X", "Z", "X"]][[1 ;; -1 - codeLogicalQubits[a]]], table = <||>},
    Do[
        With[{e = ReplacePart[ConstantArray[0, n], Thread[sub -> 1]]},
            With[{syn = Mod[checks . e, 2]}, If[! KeyExistsQ[table, syn], table[syn] = e]]
        ],
        {w, 0, n}, {sub, Subsets[Range[n], {w}]}
    ];
    table
]

reducedWeight[a_Association, x_List, z_List] := With[{n = a["Qubits"], group = codeStabilizerGroup[a]},
    Min[Total[BitOr[BitXor[x, #[[1 ;; n]]], BitXor[z, #[[n + 1 ;; 2 n]]]]] & /@ group]
]

(* ---- construction ---- *)

Options[QECErrorCorrection] = {"Order" -> "BitFlipFirst", "Ancilla" -> "Encoder"};

QECErrorCorrection::ancilla = "\"Ancilla\" must be one of `2`; got `1`.";

QECErrorCorrection[code_QECCode, opts : OptionsPattern[]] :=
    QECErrorCorrection[code, code["Qubits"], opts]

QECErrorCorrection[code_QECCode, offset_Integer, opts : OptionsPattern[]] := Module[
    {a = First[code], n = code["Qubits"], order, mode, pieces},
    order = OptionValue["Order"];
    mode = OptionValue["Ancilla"];
    Which[
        ! codeCSSQ[a],
            Message[QECErrorCorrection::css]; $Failed,
        ! MemberQ[{"BitFlipFirst", "PhaseFirst"}, order],
            Message[QECErrorCorrection::order, order]; $Failed,
        ! MemberQ[$ancillaModes, mode],
            Message[QECErrorCorrection::ancilla, mode, $ancillaModes]; $Failed,
        True,
            pieces = steanePieces[a, offset, order, mode];
            QECErrorCorrection[<|
                "Code" -> a,
                "DataQubits" -> n,
                "Offset" -> offset,
                "AncillaQubits" -> Range[offset + 1, offset + ancillaWidth[a, mode]],
                "Order" -> order,
                "Ancilla" -> mode,
                "Qubits" -> offset + ancillaWidth[a, mode],
                "Instructions" -> Catenate[pieces[[All, 2]]],
                "RegionIndices" -> piecesRegions[pieces],
                "HeraldChecks" -> piecesHeraldChecks[pieces]
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
    "OpenAssumptions", "Ancilla", "HeraldChecks", "Verifications", "ResidualWeights", "MaxResidualWeight",
    "Properties"
};

correctionData[QECErrorCorrection[a_Association]] := a

QECErrorCorrection[_Association]["Properties"] := $correctionProperties

QECErrorCorrection[a_Association][prop : ("Instructions" | "DataQubits" |
    "AncillaQubits" | "Offset" | "Order" | "Qubits")] := a[prop]

QECErrorCorrection[a_Association]["Ancilla"] := Lookup[a, "Ancilla", "Encoder"]
QECErrorCorrection[a_Association]["HeraldChecks"] := Lookup[a, "HeraldChecks", {}]
QECErrorCorrection[a_Association]["Verifications"] := Length[Lookup[a, "HeraldChecks", {}]]

(* The measured claim, region by region: the reduced weight each single fault
   leaves on the data after the gadget's own correction, over the runs the heralds
   accept.  With the bare encoder a fault in the preparation can leave a weight-2
   residual; verified, none does.  This is the property chapter 13 is for. *)
QECErrorCorrection[a_Association]["ResidualWeights"] := correctionResidualWeights[a]

QECErrorCorrection[a_Association]["MaxResidualWeight"] := Max[Values[correctionResidualWeights[a]]]

correctionResidualWeights[a_Association] := correctionResidualWeights[a] = Module[
    {instr = a["Instructions"], nq = a["Qubits"], regions, locs, checks},
    regions = correctionRegions[a];
    checks = Lookup[a, "HeraldChecks", {}];
    locs = gadgetLocationFaults[instr, nq];
    Association @ KeyValueMap[
        Function[{region, idx},
            region -> Union @ DeleteMissing @ Flatten @ KeyValueMap[
                Function[{loc, faultChoices},
                    If[ MemberQ[idx, loc[[2]]] && loc[[1]] =!= "Idle",
                        steaneResidual[a["Code"], instr, nq, checks, a["Order"], #] & /@ faultChoices,
                        {}
                    ]
                ],
                locs
            ]
        ],
        KeyDrop[regions, "Interaction"] ~Join~ KeyTake[regions, "Interaction"]
    ] ~Join~ <|
        (* waits are located by layer, not by instruction, so they get a row of their own *)
        "Waits" -> Union @ DeleteMissing @ Flatten @ KeyValueMap[
            Function[{loc, faultChoices},
                If[loc[[1]] === "Idle", steaneResidual[a["Code"], instr, nq, checks, a["Order"], #] & /@ faultChoices, {}]],
            locs]
    |>
]

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
            (* the verified gadget keeps its phase-half ancilla two blocks further up the
               workspace, so the position is read modulo the block *)
            With[{d = Select[pair, # <= n &], anc = Select[pair, # > off &]},
                Length[d] === 0 || (Length[d] === 1 && Length[anc] === 1 && Mod[First[anc] - off - 1, n] + 1 === First[d])
            ]
        ]
    ]
]

(* Section 12.3.3: none.  Kept as a property rather than a comment because it is
   the reason this gadget and not Shor EC fills Measurement.wl's slot. *)
QECErrorCorrection[a_Association]["RepetitionsNeeded"] := 1

(* the regions the gadget was built from: three for a verified one, two otherwise *)
correctionRegions[a_Association] := Lookup[a, "RegionIndices", steaneRegions[a["Code"], a["Order"]]]

QECErrorCorrection[a_Association]["Regions"] := correctionRegions[a]

QECErrorCorrection[a_Association]["DataWeights"] := Sort @ steaneResidualWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"],
    correctionRegions[a]["Interaction"]]

QECErrorCorrection[a_Association]["MaxDataWeight"] := Max @ steaneResidualWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"],
    correctionRegions[a]["Interaction"]]

(* What the missing chapter-13 preparation costs, in the same units. *)
QECErrorCorrection[a_Association]["PreparationDataWeights"] := Sort @ steaneResidualWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"],
    correctionRegions[a]["Preparation"]]

QECErrorCorrection[a_Association]["OpenAssumptions"] /; Lookup[a, "Ancilla", "Encoder"] === "Verified" := {
    "Each ancilla is verified against a second copy and the run post-selected on the \
comparison (Got26 fig. 13.3). The acceptance is part of the cost, and a protocol \
that uses the gadget has to repeat a rejected preparation."
}

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
