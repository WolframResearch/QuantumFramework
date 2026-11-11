(* ::Package:: *)

(* ============================================================================
   Tests/QEC/ExtendedRectangle.wlt

   Phase G: fault-tolerant ancilla preparation (Got26 sec. 13.1.2, fig. 13.3),
   extended rectangles (sec. 14.2), malignant sets (sec. 14.3.3) and the
   concatenation bound on the threshold (sec. 14.5, 14.7.4).

   What each value is checked against.
   - The verified ancilla: every single fault, in a run the comparison
     accepts, leaves at most one error on the data once the gadget's own
     correction is applied -- the property the bare encoder lacks (it leaves up
     to three, a logical error).
   - The exRec counts: the 7-qubit CNOT exRec has exactly the book's 263 CNOT,
     112 preparation and 112 measurement locations and its 8 post-selections
     (sec. 14.5.3); the one-qubit gates are 8 fewer because corrections stay in
     the Pauli frame, and the waits differ because the layout is an ASAP
     schedule rather than the book's hand-laid one.
   - The bound: with the book's own counts the fixed point is the book's
     3.6 x 10^-6 (eq. 14.63).
   - Malignant sets: Theorem 14.1 -- no single location of a verified exRec is
     malignant; and the fast count agrees, pair by pair, with a direct
     simulation written here independently of it.

   Everything runs on bare Wolfram.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";
qecSym[name_String] := Symbol[qecScope <> name];

steane = QECCode["SteaneCode"];
ecVerified = QECErrorCorrection[steane, "Ancilla" -> "Verified"];
ecEncoder = QECErrorCorrection[steane];
cnotRect = QECExtendedRectangle[steane, "CNOT"];
measRect = QECExtendedRectangle[steane, "Measurement"];


(* ============================================================================
   The verified ancilla
   ============================================================================ *)

(* The CSS encoder of fig. 14.8a: three Hadamards and nine CNOTs, and it prepares
   the state the code's generators and logical Z fix. *)
VerificationTest[
    With[{instr = qecSym["cssZeroInstructions"][First[steane], 0]},
        {Counts[First /@ instr],
         qecSym["applyGates"][PauliStabilizer[7],
            Replace[instr, {{"H", q_} :> "H" -> q, {"CNOT", c_, t_} :> "CNOT" -> {c, t}}, {1}]]["Expectation", #] & /@
            Join[steane["Generators"], steane["LogicalZ"]]}
    ],
    {<|"H" -> 3, "CNOT" -> 9|>, ConstantArray[1, 7]},
    TestID -> "QEC-ExRec-css-encoder-prepares-the-codeword"
]

VerificationTest[
    {ecVerified["Ancilla"], ecVerified["Verifications"], Length /@ ecVerified["HeraldChecks"][[All, 2]], ecVerified["TransversalQ"]},
    {"Verified", 2, {4, 4}, True},
    TestID -> "QEC-ExRec-verified-ec-has-two-post-selected-ancillas"
]

VerificationTest[
    {Max[ecVerified["ResidualWeights"]["Preparation"]], ecVerified["MaxResidualWeight"]},
    {1, 1},
    TestID -> "QEC-ExRec-verified-preparation-leaves-one-error"
]

VerificationTest[
    Max[ecEncoder["ResidualWeights"]["Preparation"]] > 1,
    True,
    TestID -> "QEC-ExRec-bare-encoder-leaves-more"
]


(* ============================================================================
   Extended rectangles and the bad-set bound
   ============================================================================ *)

VerificationTest[
    {cnotRect["LocationCounts"], cnotRect["Verifications"], cnotRect["PreparationLocations"], cnotRect["Parts"]},
    {<|"Measurement" -> 112, "OneQubitGate" -> 104, "Preparation" -> 112, "TwoQubitGate" -> 263, "Wait" -> 188|>,
     8, 58, <|"LeadingEC" -> 292, "Gate" -> 7, "TrailingEC" -> 292|>},
    TestID -> "QEC-ExRec-CNOT-exRec-counts"
]

VerificationTest[
    {cnotRect["BadProbabilityBound", qecP], N[cnotRect["ThresholdBound"]]},
    {Binomial[779, 2] qecP^2 / (1 - 58 qecP)^8, 3.2949505666*^-6},
    SameTest -> (Simplify[#1[[1]] - #2[[1]]] === 0 && Abs[#1[[2]] - #2[[2]]] < 10^-15 &),
    TestID -> "QEC-ExRec-CNOT-bad-bound-and-fixed-point"
]

(* The book's counts give the book's number. *)
VerificationTest[
    3.6 10^-6 < N[qecSym["concatenationThreshold"][{735, 63, 8}, 1]] < 3.71 10^-6,
    True,
    TestID -> "QEC-ExRec-book-counts-give-the-book-bound"
]


(* ============================================================================
   Malignant sets
   ============================================================================ *)

VerificationTest[
    Table[QECExtendedRectangle[steane, loc]["MalignantSingles"], {loc, {"Preparation", "Measurement", "Wait", "H", "CNOT"}}],
    {0, 0, 0, 0, 0},
    TestID -> "QEC-ExRec-no-single-location-is-malignant"
]

VerificationTest[
    QECExtendedRectangle[steane, "H", "Ancilla" -> "Encoder"]["MalignantSingles"] > 0,
    True,
    TestID -> "QEC-ExRec-bare-encoder-has-malignant-single-faults"
]

VerificationTest[
    {measRect["MalignantPairs"], cnotRect["MalignantPairs"], cnotRect["MalignantPairs"] < Binomial[779, 2] / 4},
    {507, 67051, True},
    TestID -> "QEC-ExRec-malignant-pairs-counted"
]

(* An independent oracle: simulate the exRec with both faults, apply each EC's
   correction as a Pauli at the right place, and decode ideally. *)
bruteMalignantQ[rect_, i_, j_] := Module[
    {d = First[rect], a, n, instr, nq, seam, lead, meas, hz, hx, zbar, locs, lastLead, flag, verdict},
    a = d["Code"]; n = a["Qubits"]; instr = d["Instructions"]; nq = d["Qubits"]; locs = d["LocationList"];
    seam = Length[d["Parts"]["LeadingEC"]];
    lead = Count[d["Parts"]["LeadingEC"], {"M", _}]; meas = Count[d["Parts"]["Measurement"], {"M", _}];
    hz = Most[qecSym["verificationChecks"][a, "Z"]]; hx = Most[qecSym["verificationChecks"][a, "X"]];
    zbar = First[qecSym["codeLogicalVectors"][a]["Z"]][[n + 1 ;; 2 n]];
    lastLead[q_] := Max[0, Max @ Flatten @ Position[Take[instr, seam], {_, ___, q, ___}, {1}, Heads -> False]];
    flag[loc_, f_] := loc["Type"] === "Wait" && f[[2]] <= n && f[[1]] === lastLead[f[[2]]];
    verdict[faults_, inRect_] := Module[{run, r, bit, ph, cx, cz, c1, decodedBit, emid},
        run = qecSym["framePropagate"][instr, nq, faults];
        If[! qecSym["acceptedQ"][run["Heralds"], d["HeraldChecks"]], Return[0]];
        r = run["Record"];
        bit = r[[1 ;; n]]; ph = r[[n + 1 ;; 2 n]];
        cx = qecSym["cssCorrection"][a, "X", Mod[hz . bit, 2]];
        cz = qecSym["cssCorrection"][a, "Z", Mod[hx . ph, 2]];
        c1 = DeleteCases[Table[If[cx[[q]] == 1 || cz[[q]] == 1, {seam, q, {cx[[q]], cz[[q]]}}, Nothing], {q, n}], Nothing];
        (* the decoded readout, with the leading correction applied before it *)
        decodedBit[x_] := Mod[zbar . BitXor[x, qecSym["cssCorrection"][a, "X", Mod[hz . x, 2]]], 2];
        emid = qecSym["framePropagate"][Take[instr, seam], nq, Join[Pick[faults, inRect, False], c1]]["Frame"][[1, 1 ;; n]];
        Boole[decodedBit[BitXor[r[[lead + 1 ;; lead + meas]], cx]] =!= decodedBit[emid]]
    ];
    Max @ Flatten @ Table[
        verdict[Join[ci, cj], Join[flag[locs[[i]], #] & /@ ci, flag[locs[[j]], #] & /@ cj]],
        {ci, locs[[i]]["Faults"]}, {cj, locs[[j]]["Faults"]}
    ]
]

VerificationTest[
    BlockRandom[SeedRandom[1];
        Module[{counts = qecSym["exRecMalignantCounts"][measRect], pl, nl, sample},
            pl = counts["PairList"]; nl = counts["Locations"];
            sample = Join[RandomSample[pl, 25],
                Take[Select[Table[Sort[RandomSample[Range[nl], 2]], {200}], ! MemberQ[pl, #] &], 25]];
            Count[Table[Boole[MemberQ[pl, s]] === bruteMalignantQ[measRect, s[[1]], s[[2]]], {s, sample}], True]
        ]
    ],
    50,
    TestID -> "QEC-ExRec-malignant-count-agrees-with-direct-simulation"
]


(* ============================================================================
   The concatenation threshold
   ============================================================================ *)

VerificationTest[
    With[{th = QECThreshold[steane, "Method" -> "Concatenation"]},
        {th["Method"], th["Counting"], th["Worst"], 1.45 10^-5 < N[th["Value"]] < 1.47 10^-5,
         StringContainsQ[First[th["Caveats"]], "lower bound"]}
    ],
    {"Concatenation", "Malignant", "CNOT", True, True},
    TestID -> "QEC-ExRec-concatenation-threshold-from-malignant-pairs"
]

VerificationTest[
    N[QECThreshold[steane, "Method" -> "Concatenation", "Counting" -> "BadSets"]["Value"]],
    3.2949505666*^-6,
    SameTest -> (Abs[#1 - #2] < 10^-15 &),
    TestID -> "QEC-ExRec-concatenation-threshold-from-bad-sets"
]


(* ============================================================================
   FT(C) with verified ancillas, and refusals
   ============================================================================ *)

VerificationTest[
    With[{ft = QECFaultTolerantCircuit[{{"R", 1}, {"H", 1}, {"M", 1}}, steane, "Ancilla" -> "Verified"]},
        {ft["Ancilla"], Length[ft["HeraldChecks"]], Count[ft["Instructions"], {"MH", _}], ft["Qubits"],
         StringContainsQ[First[ft["OpenAssumptions"]], "verified"]}
    ],
    {"Verified", 5, 35, 35, True},
    TestID -> "QEC-ExRec-fault-tolerant-circuit-with-verified-ancillas"
]

VerificationTest[
    QECExtendedRectangle[QECCode["5QubitCode"], "H"],
    $Failed,
    {QECExtendedRectangle::css},
    TestID -> "QEC-ExRec-refuses-a-non-CSS-code"
]

VerificationTest[
    QECExtendedRectangle[steane, "Toffoli"],
    $Failed,
    {QECExtendedRectangle::location},
    TestID -> "QEC-ExRec-refuses-a-location-without-a-gadget"
]


(* ============================================================================
   The pi/8 gate: magic state, test, injection
   ============================================================================ *)

qecVector = qecSym["instructionStateVector"];
qecMagic = qecSym["magicStateInstructions"];
qecZero = qecSym["cssZeroInstructions"];
qecPlus = qecSym["cssPlusInstructions"];
steaneRow = First[steane];
steaneXbar = Flatten @ Position[First[qecSym["codeLogicalVectors"][steaneRow]["X"]][[1 ;; 7]], 1];
zeroL = qecVector[qecZero[steaneRow, 0], 7];
oneL = qecVector[Join[qecZero[steaneRow, 0], {"X", #} & /@ steaneXbar], 7];
shiftBy[ins_, k_] := Replace[ins, {{op_String, q_Integer} :> {op, q + k}, {op_String, c_Integer, t_Integer} :> {op, c + k, t + k}}, {1}];

(* fig. 14.8b: 4 Hadamards, one T, 11 CNOTs, and the state is exactly T_L |+_L>. *)
VerificationTest[
    With[{ins = qecMagic[steaneRow, 0, "T"]},
        {Counts[First /@ ins], Chop[Abs[Conjugate[(zeroL + Exp[I Pi / 4.] oneL) / Sqrt[2.]] . qecVector[ins, 7]] - 1, 10^-10]}
    ],
    {<|"H" -> 4, "T" -> 1, "CNOT" -> 11|>, 0},
    TestID -> "QEC-ExRec-magic-encoder-makes-the-magic-state"
]

(* The transversal test operator (T X T^dagger)^(tensor 7) has the state made with T^dagger
   as its +1 eigenstate, so that is the state the test certifies -- the transversal
   map conjugates, as transversal S is logical S^dagger on this code. *)
VerificationTest[
    Table[
        With[{psi = qecVector[qecMagic[steaneRow, 0, g], 7]},
            Chop[Conjugate[psi] . qecVector[Join[qecMagic[steaneRow, 0, g], Catenate @ Table[{{"Tdg", q}, {"X", q}, {"T", q}}, {q, 7}]], 7], 10^-10]
        ],
        {g, {"T", "Tdg"}}
    ],
    {0, 1.},
    SameTest -> (Chop[#1 - #2, 10^-10] === {0, 0} &),
    TestID -> "QEC-ExRec-magic-test-certifies-the-T-dagger-state"
]

(* Fig. 13.7 on the encoded states: the outcome-0 branch is T^dagger_L applied to the
   input, and on outcome 1 the transversal X then S completes it. *)
applyOne[v_, {g_String, q_Integer}, nq_Integer] := With[
    {bit = 2^(nq - q), m = Lookup[<|"X" -> {{0, 1}, {1, 0}}, "S" -> {{1, 0}, {0, I}}|>, g]},
    With[{i0 = Select[Range[0, 2^nq - 1], BitAnd[#, bit] === 0 &] + 1},
        ReplacePart[v, Join[
            Thread[i0 -> m[[1, 1]] v[[i0]] + m[[1, 2]] v[[i0 + bit]]],
            Thread[(i0 + bit) -> m[[2, 1]] v[[i0]] + m[[2, 2]] v[[i0 + bit]]]]]
    ]
]

VerificationTest[
    Module[{dataPrep, psi, zb, proj, input, ideal, branch},
        dataPrep = Join[qecPlus[steaneRow, 0], Table[{"Sdg", q}, {q, 7}]];   (* S_L |+_L> *)
        psi = qecVector[Join[dataPrep, shiftBy[qecMagic[steaneRow, 0, "Tdg"], 7], Table[{"CNOT", 7 + q, q}, {q, 7}]], 14];
        zb = Flatten @ Position[First[qecSym["codeLogicalVectors"][steaneRow]["Z"]][[8 ;; 14]], 1];
        proj[m_] := Normalize[psi[[Select[Range[0, 2^14 - 1], BitShiftRight[#, 7] === FromDigits[m, 2] &] + 1]]];
        input = qecVector[dataPrep, 7];
        ideal = (Conjugate[zeroL] . input) zeroL + Exp[-I Pi / 4.] (Conjugate[oneL] . input) oneL;
        branch = GroupBy[Select[Tuples[{0, 1}, 7], Norm[proj[#]] > 0.5 &], Mod[Total[#[[zb]]], 2] &, First];
        {Abs[Conjugate[ideal] . proj[branch[0]]],
         Abs[Conjugate[ideal] . Fold[applyOne[#1, #2, 7] &, proj[branch[1]], Join[Table[{"X", q}, {q, 7}], Table[{"S", q}, {q, 7}]]]]}
    ],
    {1., 1.},
    SameTest -> (Max[Abs[#1 - #2]] < 10^-10 &),
    TestID -> "QEC-ExRec-magic-injection-performs-T-dagger"
]

VerificationTest[
    With[{r = QECExtendedRectangle[steane, "T"]},
        {r["LocationCounts"]["PiOver8Gate"], r["Verifications"], r["MagicPreparationLocations"], r["RectangleLocations"],
         MissingQ[r["MalignantPairs"]], N[r["ThresholdBound"]] < 3 10^-6}
    ],
    {15, 4, 357, 750, True, True},
    TestID -> "QEC-ExRec-pi-over-8-exRec"
]

VerificationTest[
    qecSym["framePropagate"][{{"H", 1}, {"T", 1}}, 1, {}],
    $Failed,
    {QECDetectorModel::nonclifford},
    TestID -> "QEC-ExRec-the-frame-refuses-T"
]

VerificationTest[
    With[{th = QECThreshold[steane, "Method" -> "Concatenation", "GateSet" -> "Universal"]},
        {th["Worst"], th["ExtendedRectangles"]["T"]["Counting"], N[th["Value"]] < N[QECThreshold[steane]["Value"]]}
    ],
    {"T", "BadSets", True},
    TestID -> "QEC-ExRec-universal-protocol-bound-is-set-by-the-pi-over-8-exRec"
]

(* Table 14.1, from its own inputs: p_T = 2.7 x 10^-5, p = p_T / 10, eps = 0.01, B = 492. *)
VerificationTest[
    With[{rows = qecSym["concatenationResources"][2.7*^-5, 492, 2.7*^-6, 0.01, 5]},
        {Round[Log10[#["MaxLocations"]], 0.01] & /@ rows, #["Overhead"] & /@ rows}
    ],
    {{3.57, 4.57, 6.57, 10.57, 18.57, 34.57}, {1, 492, 242064, 119095488, 58594980096, 28828730207232}},
    TestID -> "QEC-ExRec-resources-reproduce-table-14-1"
]
