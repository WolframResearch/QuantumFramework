(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECRegister]

PackageScope[registerIndex]
PackageScope[registerLift]
PackageScope[registerLabelMatrix]
PackageScope[registerLogicalVectors]
PackageScope[registerTransversalCNOT]
PackageScope[registerConjugate]
PackageScope[registerData]


(* ============================================================================ *)
(* Several blocks of one code, and gates between them.                          *)
(*                                                                              *)
(* Everything up to here lives on a single block: QECCode is n physical qubits  *)
(* carrying k logical ones, and a circuit numbers its qubits 1..n with ancillas *)
(* after.  That is enough for a memory experiment, which is why the p^2 of the  *)
(* previous file could be measured without this one.  It is not enough for      *)
(* anything that computes.                                                      *)
(*                                                                              *)
(* A logical qubit is a BLOCK.  A logical two-qubit gate is therefore a gate    *)
(* between two blocks, and Got26 Theorem 13.2 needs 2m blocks for a gate        *)
(* touching m of them, so blocks have to be addressable before chapter 13 can   *)
(* start.  ErrorCorrection.wl already used a second block by hand, for its      *)
(* ancilla; this file is that idea made general.                                *)
(*                                                                              *)
(* THE LAYOUT is the obvious one and is fixed here so that nothing has to guess *)
(* it: block b owns qubits (b-1)n+1 .. bn, in order.  A single-block register   *)
(* is then the identity on everything, which is the regression test that keeps  *)
(* the rest of the package honest.                                              *)
(*                                                                              *)
(* THE LABEL MATRIX IS BLOCK DIAGONAL, and its layout needs care.  A code's     *)
(* label matrix has its X and Z halves EXCHANGED, so that a plain matrix        *)
(* product computes symplectic products (ErrorRate.wl).  Lifting a block's rows *)
(* into a register therefore scatters each row into two windows, one in each    *)
(* half of the register's columns, not into one contiguous run.  Getting that   *)
(* wrong gives a matrix that still has the right shape and quietly reports the  *)
(* wrong syndromes.                                                             *)
(*                                                                              *)
(* THE TRANSVERSAL CNOT is the one two-block gate this file builds, and it is   *)
(* the one the rest of the roadmap needs: for a CSS code, a CNOT applied qubit  *)
(* by qubit between two blocks is the LOGICAL CNOT.  Steane EC already rests on *)
(* that (sec. 12.3.1) and this file is where the claim is finally checked       *)
(* rather than cited, by conjugating the logical operators through the gate and *)
(* reading the images off:                                                      *)
(*                                                                              *)
(*     Xbar_c -> Xbar_c Xbar_t      Zbar_c -> Zbar_c                            *)
(*     Xbar_t -> Xbar_t             Zbar_t -> Zbar_c Zbar_t                     *)
(*                                                                              *)
(* which is the action of a CNOT on the logical pair.  Note it is CSS-only, and *)
(* for the same reason Steane EC is: on a non-CSS code the transversal CNOT     *)
(* does not preserve the stabilizer group, so it is not a logical operation at  *)
(* all.  "LogicalAction" reports what the gate actually does, so a code for     *)
(* which it is not the CNOT says so instead of being assumed into one.          *)
(*                                                                              *)
(* References: Got26 sec. 11.2 (transversal gates and CSS codes), sec. 12.3.1   *)
(* (Steane EC, which is the first consumer), Theorem 13.2 and sec. 13.2 (gate   *)
(* teleportation, which is why blocks must be addressable at all).              *)
(* ============================================================================ *)

QECRegister::usage = "QECRegister[code, b] represents b blocks of a code laid out in order, so that block j owns qubits (j-1)n+1 through jn.\nreg[prop] gives a property; reg[\"Properties\"] lists them.\nreg[\"TransversalCNOT\", c, t] gives the instructions of a transversal CNOT from block c to block t, and reg[\"LogicalAction\", c, t] what it does to the logical operators.";

QECRegister::blocks = "The number of blocks must be a positive integer; got `1`.";
QECRegister::block = "`1` is not a block of this register; it has `2`.";
QECRegister::same = "A transversal CNOT needs two different blocks; got `1` twice.";
QECRegister::noprop = "`1` is not a property of QECRegister. Use reg[\"Properties\"] for the list.";


(* ---- the layout ---- *)

(* Block b, local qubit q, in register numbering.  One line, stated once, so that
   no consumer has to reconstruct it. *)
registerIndex[n_Integer, b_Integer, q_Integer] := (b - 1) n + q

(* A single-block instruction list, moved into block b.  Op names and arity are
   untouched; only the qubit slots move. *)
registerLift[instr_List, n_Integer, b_Integer] := Replace[instr, {
    {op_String, q_Integer} :> {op, registerIndex[n, b, q]},
    {op_String, x_Integer, y_Integer} :> {op, registerIndex[n, b, x], registerIndex[n, b, y]}
}, {1}]


(* ---- the label matrix, block by block ---- *)

(* The single-block matrix has 2n columns laid out as (z half | x half), so a row of
   block b occupies columns (b-1)n+1..bn of the first half and Bn+(b-1)n+1..Bn+bn of
   the second.  Two windows, not one. *)
registerLabelMatrix[a_Association, blocks_Integer] :=
    registerLabelMatrix[a, blocks] = With[
        {n = a["Qubits"], mat = codeLabelMatrix[a]},
        Catenate @ Table[
            Join[
                PadRight[PadLeft[mat[[All, 1 ;; n]], {Length[mat], b n}], {Length[mat], blocks n}],
                PadRight[PadLeft[mat[[All, n + 1 ;; 2 n]], {Length[mat], b n}], {Length[mat], blocks n}],
                2
            ],
            {b, blocks}
        ]
    ]

(* Each block's logical operators, lifted into register-wide Pauli rows: block b's
   X-bar has its support inside block b and identity everywhere else. *)
registerLogicalVectors[a_Association, blocks_Integer] :=
    registerLogicalVectors[a, blocks] = With[
        {n = a["Qubits"], v = codeLogicalVectors[a]},
        Association @ Table[
            key -> Catenate @ Table[
                With[{row = symplecticPart[#]},
                    Join[
                        PadRight[PadLeft[row[[1 ;; n]], b n], blocks n],
                        PadRight[PadLeft[row[[n + 1 ;; 2 n]], b n], blocks n],
                        {0}
                    ]
                ] & /@ v[key],
                {b, blocks}
            ],
            {key, {"X", "Z"}}
        ]
    ]


(* ---- the one gate between blocks ---- *)

(* Qubit by qubit, control block to target block.  Transversal by construction:
   each qubit of one block touches exactly the corresponding qubit of the other and
   nothing else, which is the property that keeps a single fault from spreading
   inside either block. *)
registerTransversalCNOT[n_Integer, c_Integer, t_Integer] :=
    Table[{"CNOT", registerIndex[n, c, q], registerIndex[n, t, q]}, {q, n}]

(* What a circuit does to a Pauli, by conjugation: put the Pauli in as the frame
   before the circuit starts and read the frame afterwards.  Phases are dropped, as
   everywhere in the frame propagator, so this is the action on the Pauli group
   modulo sign -- which is what "is this the logical CNOT" asks. *)
registerConjugate[instr_List, nq_Integer, row_List] := With[
    {run = framePropagate[
        instr, nq,
        Table[
            If[row[[q]] === 0 && row[[nq + q]] === 0, Nothing, {0, q, {row[[q]], row[[nq + q]]}}],
            {q, nq}
        ]
    ]},
    Join[run["Frame"][[1]], run["Frame"][[2]], {0}]
]


(* ---- construction ---- *)

QECRegister[code_QECCode] := QECRegister[code, 1]

QECRegister[code_QECCode, blocks_] := If[
    ! (IntegerQ[blocks] && blocks > 0),
    Message[QECRegister::blocks, blocks]; $Failed,
    QECRegister[<|"Code" -> First[code], "Blocks" -> blocks|>]
]


(* ---- properties ---- *)

$registerProperties = {
    "Code", "Blocks", "BlockQubits", "Qubits", "StabilizerCount", "LogicalQubits",
    "BlockRange", "Index", "Lift", "LabelMatrix", "LogicalVectors", "Generators",
    "TransversalCNOT", "LogicalAction", "CSSQ", "Properties"
};

registerData[QECRegister[a_Association]] := a

QECRegister[_Association]["Properties"] := $registerProperties

QECRegister[a_Association]["Blocks"] := a["Blocks"]
QECRegister[a_Association]["Code"] := QECCode[a["Code"]]
QECRegister[a_Association]["BlockQubits"] := a["Code"]["Qubits"]
QECRegister[a_Association]["Qubits"] := a["Blocks"] a["Code"]["Qubits"]
QECRegister[a_Association]["CSSQ"] := codeCSSQ[a["Code"]]
QECRegister[a_Association]["StabilizerCount"] := a["Blocks"] codeStabilizerCount[a["Code"]]
QECRegister[a_Association]["LogicalQubits"] := a["Blocks"] codeLogicalQubits[a["Code"]]

QECRegister[a_Association]["LabelMatrix"] := registerLabelMatrix[a["Code"], a["Blocks"]]
QECRegister[a_Association]["LogicalVectors"] := registerLogicalVectors[a["Code"], a["Blocks"]]

(* Every block's generators, as register-wide Pauli strings. *)
QECRegister[a_Association]["Generators"] := With[
    {n = a["Code"]["Qubits"], blocks = a["Blocks"]},
    Catenate @ Table[
        QECPauliString[
            Join[
                PadRight[PadLeft[#[[1 ;; n]], b n], blocks n],
                PadRight[PadLeft[#[[n + 1 ;; 2 n]], b n], blocks n],
                {0}
            ]
        ] & /@ a["Code"]["CheckMatrix"],
        {b, blocks}
    ]
]

QECRegister[a_Association]["BlockRange", b_Integer] := With[{n = a["Code"]["Qubits"]},
    If[ 1 <= b <= a["Blocks"],
        Range[registerIndex[n, b, 1], registerIndex[n, b, n]],
        Message[QECRegister::block, b, a["Blocks"]]; $Failed
    ]
]

QECRegister[a_Association]["Index", b_Integer, q_Integer] :=
    registerIndex[a["Code"]["Qubits"], b, q]

QECRegister[a_Association]["Lift", instr_List, b_Integer] :=
    registerLift[instr, a["Code"]["Qubits"], b]

QECRegister[a_Association]["TransversalCNOT", c_Integer, t_Integer] := Which[
    ! (1 <= c <= a["Blocks"]), Message[QECRegister::block, c, a["Blocks"]]; $Failed,
    ! (1 <= t <= a["Blocks"]), Message[QECRegister::block, t, a["Blocks"]]; $Failed,
    c === t, Message[QECRegister::same, c]; $Failed,
    True, registerTransversalCNOT[a["Code"]["Qubits"], c, t]
]

(* The claim, checked rather than cited: conjugate each block's logical operators
   through the transversal CNOT and report the images.  For a CSS code these come
   back as the CNOT action on the logical pair; for a code where they do not, this
   says so, which is the point of computing it. *)
QECRegister[a_Association]["LogicalAction", c_Integer, t_Integer] := Module[
    {instr, nq, v, k},
    instr = QECRegister[a]["TransversalCNOT", c, t];
    If[instr === $Failed, Return[$Failed]];
    nq = QECRegister[a]["Qubits"];
    v = registerLogicalVectors[a["Code"], a["Blocks"]];
    k = codeLogicalQubits[a["Code"]];
    (* Three iterators over a body that is a single Rule, so the Table is three
       levels deep and Flatten takes it to a flat list of rules.  Flatten on an
       Association is an error, not a no-op, which is how this was caught. *)
    Association @ Flatten @ Table[
        {name, b, j} -> QECPauliString[
            registerConjugate[instr, nq, v[name][[(b - 1) k + j]]]
        ],
        {name, {"X", "Z"}}, {b, {c, t}}, {j, k}
    ]
]

QECRegister[a_Association][prop_String] :=
    (Message[QECRegister::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECRegister /: MakeBoxes[
    obj : QECRegister[a_Association] /; KeyExistsQ[a, "Blocks"],
    form : (StandardForm | TraditionalForm)
] := BoxForm`ArrangeSummaryBox[
    QECRegister,
    obj,
    None,
    {
        BoxForm`SummaryItem[{"Blocks: ", a["Blocks"]}],
        BoxForm`SummaryItem[{"Qubits: ", a["Blocks"] a["Code"]["Qubits"]}],
        BoxForm`SummaryItem[{"Logical qubits: ", a["Blocks"] codeLogicalQubits[a["Code"]]}]
    },
    {
        BoxForm`SummaryItem[{"Code: ", QECCode[a["Code"]]["Parameters"]}],
        BoxForm`SummaryItem[{"CSS: ", codeCSSQ[a["Code"]]}]
    },
    form,
    "Interpretable" -> False
]
