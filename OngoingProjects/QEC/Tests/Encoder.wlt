(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Encoder.wlt

   The encoding circuit (Gottesman, QECC book Procedure 6.6 with Table 6.1).

   The oracle is the engine: run the circuit on |0...0> and ask the resulting
   stabilizer state what it is stabilized by.  A circuit that prepares a codeword
   has expectation +1 on every check, sign included -- which is what the Pauli
   fixup at the end of the encoder exists to arrange.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

(* Transitional: the rebuilt QEC core still lives under OngoingProjects/QEC/.
   Once it moves into Kernel/QEC/ and PacletInfo.wl lists its context, this Get
   disappears and the Needs above is enough.  The repo root is found from the
   loaded paclet, which RunTests.wls points at this checkout through
   PacletDirectoryLoad, so the tests always run against the source tree. *)
Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

(* Derived properties are memoised on the code's data.  A kernel that already
   held results from an earlier load would let those outrank the definitions just
   read, so the tests would silently check the old package. *)
QECClearCache[];

qecNamed = {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"};

qecCliffordGates = {"H", "S", "CNOT", "CZ", "SWAP", "X", "Y", "Z"};


(* ============================================================================
   Shape of the output
   ============================================================================ *)

VerificationTest[Head[QECCode["5QubitCode"]["EncodingCircuit"]], QuantumCircuitOperator, TestID -> "QEC-Encoder-circuit-head"]

VerificationTest[Length[QECCode["5QubitCode"]["EncodingGates"]] > 0, True, TestID -> "QEC-Encoder-gates-nonempty"]

VerificationTest[
    AllTrue[QECCode["SteaneCode"]["EncodingGates"], MemberQ[qecCliffordGates, First[#]] &],
    True,
    TestID -> "QEC-Encoder-gates-are-clifford"
]

VerificationTest[
    AllTrue[QECCode["ShorCode"]["EncodingGates"],
        With[{targets = Flatten[{Last[#]}]},
            AllTrue[targets, IntegerQ[#] && 1 <= # <= 9 &]
        ] &
    ],
    True,
    TestID -> "QEC-Encoder-gates-address-real-qubits"
]


(* ============================================================================
   The circuit prepares a codeword
   ============================================================================ *)

VerificationTest[
    AllTrue[qecNamed, QECCode[#]["EncodingCircuitValidQ"] &],
    True,
    TestID -> "QEC-Encoder-valid-named"
]

VerificationTest[
    QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}]["EncodingCircuitValidQ"],
    True,
    TestID -> "QEC-Encoder-valid-distance-two-4"
]

VerificationTest[
    QECCode[{StringRepeat["X", 6], StringRepeat["Z", 6]}]["EncodingCircuitValidQ"],
    True,
    TestID -> "QEC-Encoder-valid-distance-two-6"
]

VerificationTest[
    AllTrue[{3, 4, 5, 6}, QECCode["Repetition", #]["EncodingCircuitValidQ"] &],
    True,
    TestID -> "QEC-Encoder-valid-repetition-family"
]

(* Directly through the engine rather than through the package's own predicate:
   prepare the state and measure every check on it. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], prepared},
        prepared = Wolfram`QuantumFramework`QEC`PackageScope`applyGates[
            PauliStabilizer[code["Qubits"]], code["EncodingGates"]
        ];
        AllTrue[code["Generators"], prepared["Expectation", #] === 1 &]
    ],
    True,
    TestID -> "QEC-Encoder-prepared-state-vs-engine"
]

(* The sign fixup is doing real work: the check expectations must be +1, not -1.
   Dropping it would leave some of these at -1 and the test would catch it. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"], prepared},
        prepared = Wolfram`QuantumFramework`QEC`PackageScope`applyGates[
            PauliStabilizer[5], code["EncodingGates"]
        ];
        DeleteDuplicates[prepared["Expectation", #] & /@ code["Generators"]]
    ],
    {1},
    TestID -> "QEC-Encoder-signs-are-fixed"
]

(* The prepared state fixes the whole stabilizer group, not merely the generators. *)
VerificationTest[
    Module[{code = QECCode["BitFlipCode"], prepared},
        prepared = Wolfram`QuantumFramework`QEC`PackageScope`applyGates[
            PauliStabilizer[3], code["EncodingGates"]
        ];
        AllTrue[Subsets[code["GeneratorVectors"], {1, 2}],
            prepared["Expectation", QECPauliString[QECPauliProduct @@ #]] === 1 &
        ]
    ],
    True,
    TestID -> "QEC-Encoder-prepared-state-fixes-the-group"
]


(* ============================================================================
   Applying Paulis to a state
   ============================================================================ *)

VerificationTest[
    Module[{state},
        state = Wolfram`QuantumFramework`QEC`PackageScope`applyPauliVector[PauliStabilizer[3], "IXI"];
        state["Expectation", "IZI"]
    ],
    -1,
    TestID -> "QEC-Encoder-apply-pauli-flips"
]

VerificationTest[
    Module[{state},
        state = Wolfram`QuantumFramework`QEC`PackageScope`applyPauliVector[PauliStabilizer[3], "III"];
        state["Stabilizers"] === PauliStabilizer[3]["Stabilizers"]
    ],
    True,
    TestID -> "QEC-Encoder-apply-identity-does-nothing"
]

(* A Pauli's overall phase is a global phase, which a stabilizer state does not
   carry: -P and P must act identically on it. *)
VerificationTest[
    Module[{plus, minus},
        plus = Wolfram`QuantumFramework`QEC`PackageScope`applyPauliVector[PauliStabilizer[3], "XYZ"];
        minus = Wolfram`QuantumFramework`QEC`PackageScope`applyPauliVector[PauliStabilizer[3], "-XYZ"];
        plus["Stabilizers"] === minus["Stabilizers"]
    ],
    True,
    TestID -> "QEC-Encoder-apply-pauli-ignores-global-phase"
]

VerificationTest[
    Module[{state},
        state = Fold[
            Wolfram`QuantumFramework`QEC`PackageScope`applyPauliVector,
            PauliStabilizer[4],
            {"XYZI", "XYZI"}
        ];
        state["Stabilizers"] === PauliStabilizer[4]["Stabilizers"]
    ],
    True,
    TestID -> "QEC-Encoder-apply-pauli-twice-is-identity"
]


(* ============================================================================
   Gate inversion
   ============================================================================ *)

VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`invertGate["S" -> 2],
    {"S" -> 2, "S" -> 2, "S" -> 2},
    TestID -> "QEC-Encoder-invert-S"
]

VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`invertGate["H" -> 1],
    {"H" -> 1},
    TestID -> "QEC-Encoder-invert-H"
]

VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`invertGate["CNOT" -> {1, 2}],
    {"CNOT" -> {1, 2}},
    TestID -> "QEC-Encoder-invert-CNOT"
]

(* Each gate followed by its inverse leaves any state alone. *)
VerificationTest[
    Module[{gates, state},
        gates = {"H" -> 1, "S" -> 2, "CNOT" -> {1, 2}, "CZ" -> {2, 3}, "SWAP" -> {1, 3}};
        AllTrue[gates,
            Function[g,
                state = Wolfram`QuantumFramework`QEC`PackageScope`applyGates[
                    PauliStabilizer[3]["H", 2]["CNOT", 2, 3],
                    Prepend[Wolfram`QuantumFramework`QEC`PackageScope`invertGate[g], g]
                ];
                state["Stabilizers"] === PauliStabilizer[3]["H", 2]["CNOT", 2, 3]["Stabilizers"]
            ]
        ]
    ],
    True,
    TestID -> "QEC-Encoder-gate-inverses"
]
