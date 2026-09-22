(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECTransversalGate]

PackageScope[transversalMatrix]
PackageScope[transversalAction]
PackageScope[transversalConjugate]
PackageScope[transversalPairConjugate]
PackageScope[transversalPairCode]
PackageScope[transversalPairLogicals]
PackageScope[transversalStabilizerImages]
PackageScope[transversalDecompose]
PackageScope[transversalLogicalAction]
PackageScope[transversalInstructions]
PackageScope[transversalGateNames]
PackageScope[transversalGeneratorImages]
PackageScope[transversalName]
PackageScope[transversalValidQ]
PackageScope[transversalCliffordWords]
PackageScope[transversalData]
PackageScope[$transversalOneQubitGates]
PackageScope[$transversalTwoQubitGates]


(* ============================================================================ *)
(* Transversal gates: which Cliffords a code performs by acting qubit by qubit. *)
(*                                                                              *)
(* Got26 sec. 11.2: the transversal gates of a stabilizer code are the          *)
(* symmetries of its stabilizer.  U applied to every qubit separately is a      *)
(* valid gate gadget exactly when it maps the stabilizer group to itself; what  *)
(* it then does to the logical operators is the logical gate it performs.  Both *)
(* halves of that sentence are computed here rather than asserted.              *)
(*                                                                              *)
(* WHY THIS FILE NEEDS PHASES, AND THE FRAME PROPAGATOR WILL NOT DO.  Every     *)
(* conjugation up to now went through framePropagate, which drops signs -- fine *)
(* for Register.wl's question, "is this the logical CNOT", because the CNOT     *)
(* action has no signs in it.  It is not fine here.  The whole point of         *)
(* sec. 11.3 is a sign: on the 7-qubit code the transversal S sends             *)
(*                                                                              *)
(*     Xbar = X^7 -> Y^7    and    Y^7 = i^7 X^7 Z^7 = -i Xbar Zbar = -Ybar     *)
(*                                                                              *)
(* (eq. 11.22), so the transversal S performs the logical S^dagger, NOT the     *)
(* logical S.  A sign-blind conjugation reports "S" here and is wrong.  Hence   *)
(* the Z4 phase of Pauli.wl is carried all the way through, and membership is   *)
(* asked of codeStabilizerElement, which knows a generator's true phase, rather *)
(* than of gf2MemberQ, which only knows its symplectic part: U S U^dagger must  *)
(* land on +1 eigenoperators, not on their negatives.                           *)
(*                                                                              *)
(* THE ACTION TABLES ARE DERIVED, NOT TYPED IN.  A hand-written table of        *)
(* H: X -> Z, S: X -> Y, S: Y -> -X is exactly the kind of thing that acquires  *)
(* a wrong sign and then reports a plausible answer forever.  So the tables are *)
(* computed once from the gate's own matrix -- the engine's matrix, checked      *)
(* against it in the tests -- by conjugating each Pauli and reading the          *)
(* coefficient off a trace.  Adding a gate is adding a matrix.                   *)
(*                                                                              *)
(* WHAT COMES OUT, for the 7-qubit code (sec. 11.3, and the reason the chapter  *)
(* singles this code out):                                                      *)
(*                                                                              *)
(*     transversal H     -> logical H                                           *)
(*     transversal S     -> logical S^dagger      (the sign above)              *)
(*     transversal CNOT  -> logical CNOT          (eq. 11.23-11.26)             *)
(*                                                                              *)
(* and H, S and CNOT generate the Clifford group, so the whole logical Clifford *)
(* group is available transversally.  The pattern behind the S is the book's    *)
(* own remark: the transversal U gives the logical U*, the complex conjugate,   *)
(* which is checked here as a test rather than repeated as a slogan.            *)
(*                                                                              *)
(* AND THE CONTRAST, which is why the scan is worth having: on the 5-qubit code *)
(* neither H nor S nor the transversal CNOT preserves the stabilizer, but SH    *)
(* does -- the cyclic Clifford X -> Y -> Z -> X -- so the code is not without   *)
(* transversal gates, it has a different one.  Asking a code which Cliffords it *)
(* admits (QECTransversalGate[code, All]) answers that in one call.             *)
(*                                                                              *)
(* References: Got26 sec. 11.2 (transversal gates as stabilizer symmetries),    *)
(* sec. 11.3 (the 7-qubit code, eq. 11.18-11.26), sec. 11.4 (transversal CNOT   *)
(* and CSS codes), sec. 12.3.1 (Steane EC, the first consumer of the CNOT).     *)
(* ============================================================================ *)

QECTransversalGate::usage = "QECTransversalGate[code, gate] represents the gate applied qubit by qubit to a block of a code, where gate is \"H\", \"S\", \"Sdg\", \"V\", \"Vdg\", \"X\", \"Y\", \"Z\", \"I\", a list of those applied in order, or a two-block gate \"CNOT\", \"CZ\" or \"SWAP\".\ng[prop] gives a property; g[\"Properties\"] lists them.\nQECTransversalGate[code, All] gives the one-qubit Clifford gates that are transversal on the code, with the logical gate each one performs.";

QECTransversalGate::gate = "`1` is not a transversal gate specification. Use one of `2`, a list of them, or one of the two-block gates `3`.";
QECTransversalGate::noprop = "`1` is not a property of QECTransversalGate. Use g[\"Properties\"] for the list.";
QECTransversalGate::phase = "The conjugate of a Pauli by `1` is not a Pauli times a fourth root of unity; the gate is not Clifford.";


$transversalOneQubitGates = {"I", "X", "Y", "Z", "H", "S", "Sdg", "V", "Vdg"};

$transversalTwoQubitGates = {"CNOT", "CZ", "SWAP"};

(* The engine's own matrices, in the engine's convention (Kernel/QuantumOperator):
   V is Sqrt[X], so it sends Z to -Y.  Tests compare these against
   QuantumOperator[name]["MatrixRepresentation"] so that a drift in either
   convention fails loudly instead of quietly changing every sign below. *)
$transversalGateMatrix = <|
    "I"    -> IdentityMatrix[2],
    "X"    -> {{0, 1}, {1, 0}},
    "Y"    -> {{0, -I}, {I, 0}},
    "Z"    -> {{1, 0}, {0, -1}},
    "H"    -> {{1, 1}, {1, -1}} / Sqrt[2],
    "S"    -> {{1, 0}, {0, I}},
    "Sdg"  -> {{1, 0}, {0, -I}},
    "V"    -> {{1/2 + I/2, 1/2 - I/2}, {1/2 - I/2, 1/2 + I/2}},
    "Vdg"  -> {{1/2 - I/2, 1/2 + I/2}, {1/2 + I/2, 1/2 - I/2}},
    "CNOT" -> {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}},
    "CZ"   -> DiagonalMatrix[{1, 1, 1, -1}],
    "SWAP" -> {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}
|>;

$transversalRootPhase = <|1 -> 0, I -> 1, -1 -> 2, -I -> 3|>;


(* ---- the gate as a matrix ---- *)

transversalMatrix[gate_String] := Lookup[$transversalGateMatrix, gate, $Failed]

(* A list is a sequence in TIME order, so the matrices multiply the other way. *)
transversalMatrix[seq_List] := Fold[#2 . #1 &, IdentityMatrix[2], transversalMatrix /@ seq]

transversalQubits[gate_String] := If[MemberQ[$transversalTwoQubitGates, gate], 2, 1]
transversalQubits[_List] := 1


(* ---- the action on the Pauli group, computed from the matrix ---- *)

(* E(x, z) = i^(x z) X^x Z^z on each qubit, tensored: the same convention Pauli.wl
   uses for a row, so a row and a matrix mean the same operator. *)
transversalPauliMatrix[v_List] := With[{n = Length[v] / 2},
    Fold[
        KroneckerProduct,
        Table[
            I^(v[[q]] v[[n + q]]) *
                MatrixPower[$transversalGateMatrix["X"], v[[q]]] .
                MatrixPower[$transversalGateMatrix["Z"], v[[n + q]]],
            {q, n}
        ]
    ]
]

(* Every symplectic pattern on nq qubits, in a fixed order. *)
transversalPatterns[nq_Integer] := transversalPatterns[nq] = Table[
    Join[t[[All, 1]], t[[All, 2]]],
    {t, Tuples[Values[$pauliLetterXZ], nq]}
]

(* U E U^dagger = i^c E' for a Clifford U.  E' is found by the only nonzero trace
   overlap and c by the coefficient, which must be a fourth root of unity -- if it
   is not, the matrix was not Clifford and saying so beats returning nonsense. *)
transversalAction[gate_, nq_Integer] := transversalAction[gate, nq] = With[
    {u = transversalMatrix[gate], basis = transversalPatterns[nq]},
    Association @ Table[
        v -> Module[{m = u . transversalPauliMatrix[v] . ConjugateTranspose[u], hit},
            hit = SelectFirst[
                Table[
                    With[{c = Simplify[Tr[ConjugateTranspose[transversalPauliMatrix[w]] . m] / 2^nq]},
                        If[PossibleZeroQ[c], Nothing, {w, Lookup[$transversalRootPhase, Simplify[c], Missing["Phase"]]}]
                    ],
                    {w, basis}
                ],
                True &
            ];
            If[ MissingQ[hit] || MissingQ[hit[[2]]],
                Message[QECTransversalGate::phase, gate]; Return[$Failed, Module],
                Append[hit[[1]], hit[[2]]]
            ]
        ],
        {v, basis}
    ]
]


(* ---- conjugation of a Pauli row by the transversal gate ---- *)

(* Factors on different qubits commute and are each Hermitian in the E(x, z)
   convention, so the image is read off qubit by qubit and the phases simply add. *)
transversalConjugate[gate_, row_List] := With[
    {n = pauliQubits[row], t = transversalAction[gate, 1]},
    With[{img = Table[t[{row[[q]], row[[n + q]]}], {q, n}]},
        Join[img[[All, 1]], img[[All, 2]], {Mod[Last[row] + Total[img[[All, 3]]], 4]}]
    ]
]

(* Two blocks, laid out as Register.wl lays them out: qubit q of block 1 is q and
   qubit q of block 2 is n + q, so a row of 2n qubits has its four bits for the
   pair at q, n + q (the x half) and 2n + q, 3n + q (the z half). *)
transversalPairConjugate[gate_, row_List, n_Integer] := With[
    {t = transversalAction[gate, 2]},
    With[{img = Table[t[{row[[q]], row[[n + q]], row[[2 n + q]], row[[3 n + q]]}], {q, n}]},
        Join[
            img[[All, 1]], img[[All, 2]], img[[All, 3]], img[[All, 4]],
            {Mod[Last[row] + Total[img[[All, 5]]], 4]}
        ]
    ]
]


(* ---- two blocks as a code in their own right ---- *)

(* The stabilizer of two blocks is the direct sum, and building it as an ordinary
   code association is what lets codeStabilizerElement -- signs and all -- answer
   membership for the pair with no second implementation. *)
transversalPairCode[a_Association] := transversalPairCode[a] = With[
    {n = a["Qubits"], mat = a["CheckMatrix"], z = gf2Zero[a["Qubits"]]},
    <|
        "CheckMatrix" -> Join[
            Table[Join[r[[1 ;; n]], z, r[[n + 1 ;; 2 n]], z], {r, mat}],
            Table[Join[z, r[[1 ;; n]], z, r[[n + 1 ;; 2 n]]], {r, mat}]
        ],
        "Phases" -> Join[a["Phases"], a["Phases"]],
        "Qubits" -> 2 n
    |>
]

(* The logical operators of the pair, block by block rather than from the pair's
   own standard form: the question is what the gate does to THIS block's qubit and
   THAT block's qubit, so the basis has to be the block-wise one. *)
transversalPairLogicals[a_Association] := transversalPairLogicals[a] = With[
    {n = a["Qubits"], v = codeLogicalVectors[a], z = gf2Zero[a["Qubits"]]},
    Association @ Table[
        key -> Join[
            Table[Join[r[[1 ;; n]], z, r[[n + 1 ;; 2 n]], z, {Last[r]}], {r, v[key]}],
            Table[Join[z, r[[1 ;; n]], z, r[[n + 1 ;; 2 n]], {Last[r]}], {r, v[key]}]
        ],
        {key, {"X", "Z"}}
    ]
]


(* ---- is it a gate gadget at all ---- *)

(* Image of every generator, and whether it is still a stabilizer WITH ITS SIGN.
   codeStabilizerElement returns the group element with its true phase, so the
   comparison is on the phase too: an image equal to minus a generator preserves
   the group as a set of Paulis but not the code space, and is not a gate. *)
transversalStabilizerImages[a_Association, gate_, nq_Integer] := Module[{code, conj, gens, images},
    code = If[nq === 2, transversalPairCode[a], a];
    conj = If[nq === 2, transversalPairConjugate[gate, #, a["Qubits"]] &, transversalConjugate[gate, #] &];
    gens = codeVectors[code];
    images = conj /@ gens;
    Table[
        With[{elt = codeStabilizerElement[code, images[[i]]]},
            <|
                "Generator" -> QECPauliString[gens[[i]]],
                "Image" -> QECPauliString[images[[i]]],
                "StabilizerQ" -> (! MissingQ[elt] && Last[elt] === Last[images[[i]]])
            |>
        ],
        {i, Length[images]}
    ]
]

transversalValidQ[a_Association, gate_, nq_Integer] :=
    AllTrue[transversalStabilizerImages[a, gate, nq], #["StabilizerQ"] &]


(* ---- what it does to the logical qubits ---- *)

(* A Pauli in N(S) is (a product of logical operators) times (a stabilizer element)
   times a sign.  The exponents come from commutation -- anticommuting with Zbar_j
   is exactly what an Xbar_j factor does -- and the sign is what is left once the
   logical part and the stabilizer part are divided out.  That sign is the answer
   in sec. 11.3, so it is returned, not discarded. *)
transversalDecompose[a_Association, v_Association, w_List] := Module[
    {n = a["Qubits"], k = Length[v["X"]], ex, ez, factors, logical, residual, elt},

    ex = Table[symplecticProduct[w, v["Z"][[j]]], {j, k}];
    ez = Table[symplecticProduct[w, v["X"][[j]]], {j, k}];

    factors = Join[
        Table[If[ex[[j]] === 1, v["X"][[j]], Nothing], {j, k}],
        Table[If[ez[[j]] === 1, v["Z"][[j]], Nothing], {j, k}]
    ];

    (* i^(x z) per logical qubit, so that {1, 1} means Ybar and not Xbar Zbar. *)
    logical = MapAt[
        Mod[# + Total[ex ez], 4] &,
        If[factors === {}, pauliIdentity[n], Fold[QECPauliProduct, factors]],
        -1
    ];

    residual = Append[Mod[symplecticPart[w] + symplecticPart[logical], 2], 0];
    elt = codeStabilizerElement[a, residual];

    If[ MissingQ[elt],
        Missing["NotInNormalizer"],
        Join[ex, ez, {Mod[Last[w] - Last[QECPauliProduct[logical, elt]], 4]}]
    ]
]

(* The logical action as a list of Pauli strings on the logical qubits: the images
   of Xbar_1..Xbar_k first, then of Zbar_1..Zbar_k. *)
transversalLogicalAction[a_Association, gate_, nq_Integer] := Module[{code, logicals, conj},
    code = If[nq === 2, transversalPairCode[a], a];
    logicals = If[nq === 2, transversalPairLogicals[a], codeLogicalVectors[a]];
    conj = If[nq === 2, transversalPairConjugate[gate, #, a["Qubits"]] &, transversalConjugate[gate, #] &];
    Table[
        Replace[
            transversalDecompose[code, logicals, conj[row]],
            {r_List :> QECPauliString[r], m_ :> m}
        ],
        {row, Join[logicals["X"], logicals["Z"]]}
    ]
]


(* ---- naming the logical gate ---- *)

(* The same action, computed on bare qubits: the images of X_j and Z_j identify a
   Clifford up to a global phase, which is all a gate gadget determines anyway. *)
transversalGeneratorImages[gate_, nq_Integer] := With[
    {t = transversalAction[gate, nq], z = gf2Zero[nq]},
    Join[
        Table[QECPauliString[t[Join[Normal @ UnitVector[nq, j], z]]], {j, nq}],
        Table[QECPauliString[t[Join[z, Normal @ UnitVector[nq, j]]]], {j, nq}]
    ]
]

(* One representative of each of the 24 one-qubit Cliffords, as a word in H and S,
   with the named gates preferred as labels: the naming table is built by looking
   up an action, so a gate that has a name gets it and the rest keep their word. *)
transversalCliffordWords[nq_Integer] := transversalCliffordWords[nq] = DeleteDuplicatesBy[
    Join[{{}}, Catenate[Table[Tuples[{"H", "S"}, l], {l, 1, 6}]]],
    transversalGeneratorImages[#, nq] &
]

(* The naming table, images -> label.  The named gates go in LAST so that they win
   the key they share with a word: the transversal S of the 7-qubit code is to be
   reported as "Sdg", not as "SSS". *)
transversalGateNames[nq_Integer] := transversalGateNames[nq] = Association @ Which[
    nq === 1,
        Join[
            Table[transversalGeneratorImages[w, 1] -> If[w === {}, "I", StringJoin[w]], {w, transversalCliffordWords[1]}],
            Table[transversalGeneratorImages[g, 1] -> g, {g, $transversalOneQubitGates}]
        ],
    nq === 2,
        Table[transversalGeneratorImages[g, 2] -> g, {g, $transversalTwoQubitGates}],
    True, {}
]

transversalName[action_List, nq_Integer] := Lookup[transversalGateNames[nq], Key[action], Missing["NotNamed"]]


(* ---- the circuit ---- *)

(* Qubit by qubit, which is the definition.  A sequence is gate-major, so its depth
   is the length of the sequence and not the number of qubits.  "I" emits nothing:
   an identity gate is not a fault location of its own, it is a wait, and waits are
   the noise model's business (Noise.wl), not this file's. *)
transversalInstructions["I", _Integer] := {}

transversalInstructions[gate_String, n_Integer] /; MemberQ[$transversalOneQubitGates, gate] :=
    Table[{gate, q}, {q, n}]

transversalInstructions[gate_String, n_Integer] /; MemberQ[$transversalTwoQubitGates, gate] :=
    Table[{gate, q, n + q}, {q, n}]

transversalInstructions[seq_List, n_Integer] :=
    Catenate @ Table[transversalInstructions[g, n], {g, seq}]


(* ---- construction ---- *)

transversalGateQ[gate_String] :=
    MemberQ[$transversalOneQubitGates, gate] || MemberQ[$transversalTwoQubitGates, gate]

transversalGateQ[seq_List] := AllTrue[seq, StringQ[#] && MemberQ[$transversalOneQubitGates, #] &]

transversalGateQ[_] := False

QECTransversalGate[code_QECCode, gate_] := If[
    ! transversalGateQ[gate],
    Message[QECTransversalGate::gate, gate, $transversalOneQubitGates, $transversalTwoQubitGates]; $Failed,
    QECTransversalGate[<|"Code" -> First[code], "Gate" -> gate|>]
]

(* The scan: which one-qubit Cliffords are transversal on this code, and what each
   one does.  Twenty-four candidates, one representative per class. *)
QECTransversalGate[code_QECCode, All] := With[{a = First[code]},
    Association @ Table[
        With[{name = transversalName[transversalGeneratorImages[w, 1], 1]},
            If[ transversalValidQ[a, w, 1],
                name -> transversalLogicalAction[a, w, 1],
                Nothing
            ]
        ],
        {w, transversalCliffordWords[1]}
    ]
]


(* ---- properties ---- *)

$transversalProperties = {
    "Code", "Gate", "Blocks", "Qubits", "LogicalQubits", "Instructions", "Depth",
    "QuantumCircuitOperator", "Diagram",
    "TransversalQ", "StabilizerImages", "LogicalAction", "LogicalGate", "Properties"
};

transversalData[QECTransversalGate[a_Association]] := a

QECTransversalGate[_Association]["Properties"] := $transversalProperties

QECTransversalGate[a_Association]["Code"] := QECCode[a["Code"]]
QECTransversalGate[a_Association]["Gate"] := a["Gate"]
QECTransversalGate[a_Association]["Blocks"] := transversalQubits[a["Gate"]]
QECTransversalGate[a_Association]["Qubits"] :=
    transversalQubits[a["Gate"]] a["Code"]["Qubits"]
QECTransversalGate[a_Association]["LogicalQubits"] :=
    transversalQubits[a["Gate"]] codeLogicalQubits[a["Code"]]

QECTransversalGate[a_Association]["Instructions"] :=
    transversalInstructions[a["Gate"], a["Code"]["Qubits"]]

QECTransversalGate[a_Association]["Depth"] :=
    If[ListQ[a["Gate"]], Length[DeleteCases[a["Gate"], "I"]], If[a["Gate"] === "I", 0, 1]]

QECTransversalGate[a_Association]["QuantumCircuitOperator"] :=
    gadgetCircuitOperator[QECTransversalGate[a]["Instructions"]]
QECTransversalGate[a_Association]["Diagram"] :=
    gadgetDiagram[QECTransversalGate[a]["Instructions"]]

QECTransversalGate[a_Association]["TransversalQ"] :=
    transversalValidQ[a["Code"], a["Gate"], transversalQubits[a["Gate"]]]

QECTransversalGate[a_Association]["StabilizerImages"] := Dataset @ transversalStabilizerImages[
    a["Code"], a["Gate"], transversalQubits[a["Gate"]]
]

(* Reported as an association Xbar_j / Zbar_j -> image, which is how sec. 11.3
   writes it.  A gate that is not transversal has no logical action: the images are
   there to be looked at, but they are not a gate on the code space. *)
QECTransversalGate[a_Association]["LogicalAction"] := Module[{nq, k, action},
    nq = transversalQubits[a["Gate"]];
    If[! transversalValidQ[a["Code"], a["Gate"], nq], Return[Missing["NotTransversal"]]];
    k = nq codeLogicalQubits[a["Code"]];
    action = transversalLogicalAction[a["Code"], a["Gate"], nq];
    Association @ Join[
        Table[{"X", j} -> action[[j]], {j, k}],
        Table[{"Z", j} -> action[[k + j]], {j, k}]
    ]
]

QECTransversalGate[a_Association]["LogicalGate"] := Module[{nq, action},
    nq = transversalQubits[a["Gate"]];
    If[! transversalValidQ[a["Code"], a["Gate"], nq], Return[Missing["NotTransversal"]]];
    action = transversalLogicalAction[a["Code"], a["Gate"], nq];
    If[ Length[action] / 2 <= 2,
        transversalName[action, Length[action] / 2],
        Missing["NotNamed"]
    ]
]

QECTransversalGate[a_Association][prop_String] :=
    (Message[QECTransversalGate::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECTransversalGate /: MakeBoxes[
    obj : QECTransversalGate[a_Association] /; KeyExistsQ[a, "Gate"],
    form : (StandardForm | TraditionalForm)
] := BoxForm`ArrangeSummaryBox[
    QECTransversalGate,
    obj,
    None,
    {
        BoxForm`SummaryItem[{"Gate: ", If[ListQ[a["Gate"]], StringRiffle[a["Gate"], " "], a["Gate"]]}],
        BoxForm`SummaryItem[{"Transversal: ", obj["TransversalQ"]}],
        BoxForm`SummaryItem[{"Logical gate: ", obj["LogicalGate"]}]
    },
    {
        BoxForm`SummaryItem[{"Code: ", QECCode[a["Code"]]["Parameters"]}],
        BoxForm`SummaryItem[{"Blocks: ", transversalQubits[a["Gate"]]}]
    },
    form,
    "Interpretable" -> False
]
