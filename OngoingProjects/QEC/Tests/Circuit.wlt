(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Circuit.wlt

   The syndrome-extraction circuit and the Pauli frame propagator.

   The oracle here is codeSyndrome, which computes a syndrome by matrix product
   and knows nothing about circuits.  Emitting a circuit, pushing a Pauli frame
   through it, and reading the ancilla measurements has to reproduce that answer
   for every error -- including on the 5-qubit code, whose generators carry Y and
   so exercise the sqrt(X) basis change.

   The second oracle is the engine itself: the same circuit, run as an actual
   stabilizer state on n + m qubits, with the checks read off by expectation.

   Note on reading a deterministic outcome from the engine: this file uses
   (1 - ps["Expectation", p]) / 2 and never ps["M", q].  See
   OngoingProjects/QEC/EngineMeasurementBug.md -- ["M", q] returns an outcome that
   depends on which generating set the state was built from, and syndrome
   extraction builds states from generators.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecFrame = Symbol[qecScope <> "framePropagate"];
qecWeightK = Symbol[qecScope <> "weightKVectors"];
qecCompleted = Symbol[qecScope <> "codeCompletedGenerators"];
qecCodeData = Symbol[qecScope <> "codeData"];
qecEngineGates = Symbol[qecScope <> "instructionEngineGates"];

qecNamed = {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"};

(* The syndrome as the circuit sees it: put the error in as a frame before the
   circuit starts and read the measurement record. *)
qecFrameSyndrome[code_, v_List] := Module[{n = code["Qubits"], circ = code["SyndromeCircuit"]},
    qecFrame[
        circ["Instructions"], circ["Qubits"],
        Table[If[v[[q]] == 0 && v[[n + q]] == 0, Nothing, {0, q, {v[[q]], v[[n + q]]}}], {q, n}]
    ]["Record"]
];

qecFrameAgrees[code_] := With[
    {vs = Join @@ Table[qecWeightK[code["Qubits"], w], {w, 1, 2}]},
    AllTrue[vs, qecFrameSyndrome[code, #] === code["Syndrome", QECPauliString[#]] &]
];

(* The same circuit run on a real state through the engine. *)
qecPad[s_String, m_] := s <> StringRepeat["I", m];
qecAncZ[n_, m_, j_] := StringJoin[Table[If[k == n + j, "Z", "I"], {k, n + m}]];

qecEngineSyndrome[code_, err_String] := Module[
    {n = code["Qubits"], m = code["StabilizerCount"], base, ps},
    base = Join[
        qecPad[#, m] & /@ (QECPauliString /@ qecCompleted[qecCodeData[code]]),
        Table[qecAncZ[n, m, j], {j, m}]
    ];
    ps = Fold[
        Function[{s, qp}, If[qp[[2]] === "I", s, s[qp[[2]], qp[[1]]]]],
        PauliStabilizer[base],
        Table[{k, StringTake[qecPad[err, m], {k}]}, {k, n + m}]
    ];
    ps = ps["ApplyCircuit", qecEngineGates[code["SyndromeCircuit"]["Instructions"]]];
    Table[(1 - ps["Expectation", qecAncZ[n, m, j]]) / 2, {j, m}]
];

qecEngineAgrees[code_] := Module[{n = code["Qubits"], ref},
    ref = qecEngineSyndrome[code, StringRepeat["I", n]];
    AllTrue[
        Join @@ Table[qecWeightK[n, w], {w, 1, 2}],
        BitXor[ref, qecEngineSyndrome[code, QECPauliString[#]]] === code["Syndrome", QECPauliString[#]] &
    ]
];


(* ============================================================================
   Shape of the emitted circuit
   ============================================================================ *)

VerificationTest[
    Head[QECCode["BitFlipCode"]["SyndromeCircuit"]],
    QECSyndromeCircuit,
    TestID -> "QEC-Circuit-head"
]

(* ZZI needs no rotation at all: two CNOTs into the ancilla and a measurement. *)
VerificationTest[
    QECCode["BitFlipCode"]["SyndromeCircuit"]["Instructions"],
    {{"R", 4}, {"CNOT", 1, 4}, {"CNOT", 2, 4}, {"M", 4},
     {"R", 5}, {"CNOT", 2, 5}, {"CNOT", 3, 5}, {"M", 5}},
    TestID -> "QEC-Circuit-bitflip-instructions"
]

(* An X-type generator rotates its data qubits with H and back again. *)
VerificationTest[
    QECCode[{"XXI"}]["SyndromeCircuit"]["Instructions"],
    {{"R", 4}, {"H", 1}, {"H", 2}, {"CNOT", 1, 4}, {"CNOT", 2, 4}, {"H", 1}, {"H", 2}, {"M", 4}},
    TestID -> "QEC-Circuit-x-type-rotations"
]

(* A Y-type generator uses the engine's sqrt(X), which sends Y to Z, and undoes it. *)
VerificationTest[
    QECCode[{"YY"}]["SyndromeCircuit"]["Instructions"],
    {{"R", 3}, {"V", 1}, {"V", 2}, {"CNOT", 1, 3}, {"CNOT", 2, 3}, {"Vdg", 1}, {"Vdg", 2}, {"M", 3}},
    TestID -> "QEC-Circuit-y-type-rotations"
]

VerificationTest[
    QECCode["SteaneCode"]["SyndromeCircuit"]["Qubits"],
    13,
    TestID -> "QEC-Circuit-qubit-count"
]

(* The ancillas are reset and reused, so r rounds do not need r m of them. *)
VerificationTest[
    QECCode["SteaneCode"]["SyndromeCircuit", 5]["Qubits"],
    13,
    TestID -> "QEC-Circuit-ancillas-reused"
]

VerificationTest[
    QECCode["SteaneCode"]["SyndromeCircuit", 4]["MeasurementCount"],
    24,
    TestID -> "QEC-Circuit-measurement-count"
]

VerificationTest[
    Length[QECCode["ShorCode"]["SyndromeCircuit", 3]["Instructions"]],
    3 Length[QECCode["ShorCode"]["SyndromeCircuit", 1]["Instructions"]],
    TestID -> "QEC-Circuit-rounds-repeat"
]

VerificationTest[
    QECCode["BitFlipCode"]["SyndromeCircuit", 2]["MeasurementLabels"],
    {{1, 1}, {1, 2}, {2, 1}, {2, 2}},
    TestID -> "QEC-Circuit-measurement-labels"
]

VerificationTest[
    QECCode["BitFlipCode"]["SyndromeCircuit"]["GateCounts"],
    <|"R" -> 2, "CNOT" -> 4, "M" -> 2|>,
    TestID -> "QEC-Circuit-gate-counts"
]

VerificationTest[
    QECSyndromeCircuit[QECCode["BitFlipCode"], 0],
    $Failed,
    {QECSyndromeCircuit::rounds},
    TestID -> "QEC-Circuit-bad-rounds"
]

VerificationTest[
    QECCode["BitFlipCode"]["SyndromeCircuit"]["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECSyndromeCircuit::noprop},
    TestID -> "QEC-Circuit-bad-property"
]


(* ============================================================================
   The frame propagator against codeSyndrome
   ============================================================================ *)

Table[
    VerificationTest[
        qecFrameAgrees[QECCode[name]],
        True,
        TestID -> "QEC-Circuit-frame-matches-syndrome-" <> name
    ],
    {name, qecNamed}
]

(* An error and a stabilizer differ by nothing a check can see. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        qecFrameSyndrome[c, QECPauliVector["XII"]] ===
            qecFrameSyndrome[c, QECPauliProduct["XII", "ZZI"]]
    ],
    True,
    TestID -> "QEC-Circuit-frame-stabilizer-invariance"
]

(* Frames add: the record of a product is the XOR of the records. *)
VerificationTest[
    With[{c = QECCode["5QubitCode"], u = QECPauliVector["XIYIZ"], v = QECPauliVector["IZZXI"]},
        qecFrameSyndrome[c, QECPauliVector[QECPauliProduct[u, v]]] ===
            BitXor[qecFrameSyndrome[c, u], qecFrameSyndrome[c, v]]
    ],
    True,
    TestID -> "QEC-Circuit-frame-linearity"
]

(* The ancilla is reset between rounds, so a frame injected before the circuit is
   seen identically in every round. *)
VerificationTest[
    With[{c = QECCode["SteaneCode"], circ = QECCode["SteaneCode"]["SyndromeCircuit", 3]},
        With[{rec = qecFrame[circ["Instructions"], circ["Qubits"], {{0, 1, {1, 0}}}]["Record"]},
            rec === Join[#, #, #] & @ c["Syndrome", "XIIIIII"]
        ]
    ],
    True,
    TestID -> "QEC-Circuit-frame-repeats-across-rounds"
]


(* ============================================================================
   The circuit against the engine
   ============================================================================ *)

Table[
    VerificationTest[
        qecEngineAgrees[QECCode[name]],
        True,
        TestID -> "QEC-Circuit-engine-matches-syndrome-" <> name
    ],
    {name, qecNamed}
]


(* ============================================================================
   How a single ancilla fault spreads

   These lock in a mechanism that is easy to get backwards, and that the file
   headers and the tech note now assert precisely.  Gottesman's figure 12.1b has
   the ancilla as the *control*, where an X on it halfway through the ladder
   corrupts several data qubits.  Here the ancilla is the target, so the roles of
   X and Z are exchanged, and transcribing the book would be wrong.
   ============================================================================ *)

(* The Pauli left on the data, and the measurement record, when a single fault
   strikes the ancilla just after instruction i. *)
qecAncillaFault[code_, i_, ancilla_, pauli_] := Module[{n = code["Qubits"], circ, run, data},
    circ = code["SyndromeCircuit"];
    run = qecFrame[circ["Instructions"], circ["Qubits"], {{i, ancilla, pauli}}];
    data = Join[run["Frame"][[1, 1 ;; n]], run["Frame"][[2, 1 ;; n]], {0}];
    <|"Data" -> QECPauliString[data], "Record" -> run["Record"],
      "StabilizerQ" -> code["StabilizerMemberQ", data]|>
];

(* An X on the ancilla reaches no data qubit: CNOT carries X only from control to
   target, so it sits on the ancilla and flips one readout.  A pure measurement
   error, whatever point of the ladder it strikes. *)
VerificationTest[
    With[{bf = QECCode["BitFlipCode"]},
        Table[qecAncillaFault[bf, i, 4, {1, 0}], {i, {1, 2, 3}}]
    ],
    Table[<|"Data" -> "III", "Record" -> {1, 0}, "StabilizerQ" -> True|>, 3],
    TestID -> "QEC-Circuit-ancilla-X-is-readout-error"
]

(* A Z on the ancilla is the one that spreads: CNOT carries Z from target back to
   control, and the ancilla keeps its own Z, so it lands on the data qubits the
   ladder has still to touch. *)
VerificationTest[
    With[{bf = QECCode["BitFlipCode"]},
        #["Data"] & /@ Table[qecAncillaFault[bf, i, 4, {0, 1}], {i, {1, 2, 3}}]
    ],
    {"ZZI", "IZI", "III"},
    TestID -> "QEC-Circuit-ancilla-Z-spreads-to-the-tail"
]

(* Before the first CNOT it leaves the whole generator, which is a stabilizer and
   harmless; after the last it leaves nothing.  Only mid-ladder hurts. *)
VerificationTest[
    With[{bf = QECCode["BitFlipCode"]},
        #["StabilizerQ"] & /@ Table[qecAncillaFault[bf, i, 4, {0, 1}], {i, {1, 2, 3}}]
    ],
    {True, False, True},
    TestID -> "QEC-Circuit-only-mid-ladder-hurts"
]

(* And it does not fire its own generator's check in the round it happened: the
   residual is a sub-product of that generator, so it commutes with it. *)
VerificationTest[
    With[{bf = QECCode["BitFlipCode"]},
        First[qecAncillaFault[bf, 2, 4, {0, 1}]["Record"]]
    ],
    0,
    TestID -> "QEC-Circuit-hook-invisible-to-its-own-check"
]

(* On a non-CSS generator the un-rotation sends the residual Z back to whatever
   letter the generator carried, so what is left is literally the un-extracted
   tail of XZZXI. *)
VerificationTest[
    With[{c5 = QECCode["5QubitCode"]},
        #["Data"] & /@ Table[qecAncillaFault[c5, i, 6, {0, 1}], {i, {1, 4, 5, 6, 7}}]
    ],
    {"XZZXI", "IZZXI", "IIZXI", "IIIXI", "IIIII"},
    TestID -> "QEC-Circuit-hook-leaves-the-generator-tail"
]

(* The whole point of the above: one fault, a multi-qubit data error. *)
VerificationTest[
    QECPauliWeight[qecAncillaFault[QECCode["5QubitCode"], 4, 6, {0, 1}]["Data"]],
    3,
    TestID -> "QEC-Circuit-one-fault-many-errors"
]
