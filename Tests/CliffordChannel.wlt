(* ==========================================================================
   Tests/CliffordChannel.wlt -- Phase 8.1 CliffordChannel head per Yashin25.

   Phase 8.1 covers:
     - Head construction + predicate
     - Pure-state -> Choi-tableau constructor (PauliStabilizer -> CliffordChannel)
     - Identity-channel constructor
     - Basic accessors (UA, UB, c, InputQubits, OutputQubits, Rank, Tableau)

   Phase 8.2 (composition via vector-space intersection) is stubbed and not
   tested here; tests there will land with the composition implementation.
   ========================================================================== *)

Needs["Wolfram`QuantumFramework`"];


(* ============================================================================ *)
(* TIER A -- Predicate + identity-channel sanity                                *)
(* ============================================================================ *)

VerificationTest[
    Wolfram`QuantumFramework`PackageScope`CliffordChannelQ @ CliffordChannel["Identity", 1],
    True,
    TestID -> "Phase8.1-CliffordChannelQ-Identity-1"
]

VerificationTest[
    Wolfram`QuantumFramework`PackageScope`CliffordChannelQ @ CliffordChannel["Identity", 3],
    True,
    TestID -> "Phase8.1-CliffordChannelQ-Identity-3"
]

VerificationTest[
    CliffordChannel["Identity", 2]["InputQubits"],
    2,
    TestID -> "Phase8.1-Identity-InputQubits"
]

VerificationTest[
    CliffordChannel["Identity", 2]["OutputQubits"],
    2,
    TestID -> "Phase8.1-Identity-OutputQubits"
]

VerificationTest[
    CliffordChannel["Identity", 2]["Rank"],
    4,   (* 2*n = 4 rows for a 2-qubit identity *)
    TestID -> "Phase8.1-Identity-Rank-2qubit"
]

(* Identity-channel UA = UB = I_{2n}: the channel maps every Pauli to itself. *)
VerificationTest[
    CliffordChannel["Identity", 2]["UA"],
    IdentityMatrix[4],
    TestID -> "Phase8.1-Identity-UA-IsIdentity"
]

VerificationTest[
    CliffordChannel["Identity", 2]["UB"],
    IdentityMatrix[4],
    TestID -> "Phase8.1-Identity-UB-IsIdentity"
]

VerificationTest[
    CliffordChannel["Identity", 2]["c"],
    {0, 0, 0, 0},
    TestID -> "Phase8.1-Identity-c-IsZero"
]


(* ============================================================================ *)
(* TIER B -- PauliStabilizer -> CliffordChannel (state-preparation Choi tableau)*)
(* ============================================================================ *)

(* For a state, |A| = 0 and |B| = n; the tableau has n rows on the output side,*)
(* and UA is empty. *)
VerificationTest[
    CliffordChannel[PauliStabilizer[1]]["InputQubits"],
    0,
    TestID -> "Phase8.1-PSStateAsChannel-NoInput"
]

VerificationTest[
    CliffordChannel[PauliStabilizer[1]]["OutputQubits"],
    1,
    TestID -> "Phase8.1-PSStateAsChannel-OneOutput"
]

VerificationTest[
    CliffordChannel[PauliStabilizer[1]]["Rank"],
    1,    (* one stabilizer for |0> *)
    TestID -> "Phase8.1-PSStateAsChannel-Rank-1qubit"
]

VerificationTest[
    Wolfram`QuantumFramework`PackageScope`CliffordChannelQ @ CliffordChannel[PauliStabilizer["5QubitCode"]],
    True,
    TestID -> "Phase8.1-PSStateAsChannel-5QCode-Valid"
]

VerificationTest[
    CliffordChannel[PauliStabilizer["5QubitCode"]]["Rank"],
    5,    (* 5 stabilizers for a [[5,1,3]] code state *)
    TestID -> "Phase8.1-PSStateAsChannel-5QCode-Rank"
]


(* Bell-state Choi tableau: 2 qubits out, k = 2 rows for stabilizers {XX, ZZ}. *)
VerificationTest[
    With[{
        cc = CliffordChannel @ PauliStabilizer @ QuantumCircuitOperator[
            {"H" -> 1, "CNOT" -> {1, 2}}
        ]
    },
        {cc["InputQubits"], cc["OutputQubits"], cc["Rank"]}
    ],
    {0, 2, 2},
    TestID -> "Phase8.1-PSStateAsChannel-Bell-Shape"
]


(* ============================================================================ *)
(* TIER C -- Tableau accessor: assembled [UA | UB | c] matrix                  *)
(* ============================================================================ *)

VerificationTest[
    Dimensions @ CliffordChannel["Identity", 1]["Tableau"],
    {2, 5},   (* k=2 rows, columns = 2|A| + 2|B| + 1 = 2 + 2 + 1 = 5 *)
    TestID -> "Phase8.1-Tableau-Dimensions-Identity-1"
]

VerificationTest[
    Dimensions @ CliffordChannel["Identity", 3]["Tableau"],
    {6, 13},   (* k=6, columns = 6+6+1 *)
    TestID -> "Phase8.1-Tableau-Dimensions-Identity-3"
]

(* For a state-only channel (nA=0), the Tableau drops the empty UA columns: *)
VerificationTest[
    Dimensions @ CliffordChannel[PauliStabilizer["5QubitCode"]]["Tableau"],
    {5, 11},   (* k=5, columns = 0 + 10 + 1 = 11 *)
    TestID -> "Phase8.1-Tableau-Dimensions-5QState"
]


(* ============================================================================ *)
(* TIER D -- Composition stub: returns $Failed with explanatory message        *)
(* ============================================================================ *)

VerificationTest[
    CliffordChannel["Identity", 1] @ CliffordChannel["Identity", 1],
    $Failed,
    {CliffordChannel::compose},
    TestID -> "Phase8.1-Compose-Stub-EmitsMessage"
]


(* ============================================================================ *)
(* TIER E -- More PauliStabilizer-as-CliffordChannel cases                     *)
(* ============================================================================ *)

(* Steane code: 7 qubits, 7 stabilizers; tableau is 7 x (0 + 14 + 1) = 7 x 15. *)
VerificationTest[
    Dimensions @ CliffordChannel[PauliStabilizer["SteaneCode"]]["Tableau"],
    {7, 15},
    TestID -> "Phase8.1-Tableau-Dimensions-Steane"
]

VerificationTest[
    CliffordChannel[PauliStabilizer["SteaneCode"]]["Rank"],
    7,
    TestID -> "Phase8.1-Steane-Rank"
]

VerificationTest[
    CliffordChannel[PauliStabilizer["9QubitCode"]]["Rank"],
    9,
    TestID -> "Phase8.1-9QubitCode-Rank"
]

(* Random-Clifford state of size n -> CliffordChannel of rank n. *)
VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                With[{n = RandomInteger[{1, 5}]},
                    CliffordChannel[PauliStabilizer["Random", n]]["Rank"] == n
                ],
                {15}
            ],
            TrueQ
        ]
    ],
    True,
    TestID -> "Phase8.1-RandomClifford-RankMatchesQubits"
]


(* ============================================================================ *)
(* TIER F -- Identity-channel sanity for varied n                              *)
(* ============================================================================ *)

VerificationTest[
    Dimensions @ CliffordChannel["Identity", 5]["Tableau"],
    {10, 21},  (* k=10, columns = 10 + 10 + 1 *)
    TestID -> "Phase8.1-Identity-Tableau-Dimensions-5q"
]

VerificationTest[
    Block[{},
        AllTrue[
            Table[
                With[{cc = CliffordChannel["Identity", n]},
                    cc["Rank"] == 2 n &&
                    cc["UA"] == cc["UB"] &&
                    cc["c"] === ConstantArray[0, 2 n]
                ],
                {n, 1, 6}
            ],
            TrueQ
        ]
    ],
    True,
    TestID -> "Phase8.1-Identity-Invariants-Across-n"
]


(* ============================================================================ *)
(* TIER G -- Choi-tableau invariants on the UB side                            *)
(* ============================================================================ *)

(* For a state-preparation Choi tableau, the rows of UB should be linearly     *)
(* independent over F_2 (because they're stabilizer generators of a pure       *)
(* state). Rank check via MatrixRank[matrix, Modulus -> 2].                    *)
VerificationTest[
    Block[{},
        AllTrue[
            {
                CliffordChannel[PauliStabilizer[1]],
                CliffordChannel[PauliStabilizer["5QubitCode"]],
                CliffordChannel[PauliStabilizer["SteaneCode"]],
                CliffordChannel @ PauliStabilizer @ QuantumCircuitOperator[
                    {"H" -> 1, "CNOT" -> {1, 2}}
                ]
            },
            MatrixRank[#["UB"], Modulus -> 2] == #["Rank"] &
        ]
    ],
    True,
    TestID -> "Phase8.1-UB-Rank-MatchesGenerators"
]

(* UB shape sanity: k rows, 2*nB columns. *)
VerificationTest[
    With[{cc = CliffordChannel[PauliStabilizer["5QubitCode"]]},
        Dimensions[cc["UB"]] === {cc["Rank"], 2 cc["OutputQubits"]}
    ],
    True,
    TestID -> "Phase8.1-UB-ShapeSanity"
]


(* ============================================================================ *)
(* TIER H -- Properties contract                                               *)
(* ============================================================================ *)

VerificationTest[
    SubsetQ[
        CliffordChannel[PauliStabilizer[1]]["Properties"],
        {"UA", "UB", "c", "InputQubits", "OutputQubits", "Rank", "Tableau"}
    ],
    True,
    TestID -> "Phase8.1-Properties-Contains-Core"
]

VerificationTest[
    CliffordChannel[PauliStabilizer[1]]["Source"],
    "PauliStabilizer",
    TestID -> "Phase8.1-Source-Tag-PS"
]

VerificationTest[
    CliffordChannel["Identity", 2]["Source"],
    "Identity",
    TestID -> "Phase8.1-Source-Tag-Identity"
]


(* ============================================================================ *)
(* TIER I -- Validity guard reuses cached HoldValid                            *)
(* ============================================================================ *)

(* Predicate hits HoldValidQ on an already-validated CliffordChannel: same     *)
(* result on second call (consistent dispatch).                                 *)
VerificationTest[
    Module[{cc = CliffordChannel[PauliStabilizer[1]], q1, q2},
        q1 = Wolfram`QuantumFramework`PackageScope`CliffordChannelQ[cc];
        q2 = Wolfram`QuantumFramework`PackageScope`CliffordChannelQ[cc];
        q1 === q2 === True
    ],
    True,
    TestID -> "Phase8.1-Predicate-Idempotent"
]


(* ============================================================================ *)
(* TIER J -- Composition stub fires for mixed identity / state inputs          *)
(* ============================================================================ *)

VerificationTest[
    CliffordChannel["Identity", 2] @ CliffordChannel[PauliStabilizer[2]],
    $Failed,
    {CliffordChannel::compose},
    TestID -> "Phase8.1-Compose-Stub-Identity-on-State"
]

VerificationTest[
    CliffordChannel[PauliStabilizer["5QubitCode"]] @ CliffordChannel["Identity", 5],
    $Failed,
    {CliffordChannel::compose},
    TestID -> "Phase8.1-Compose-Stub-State-on-Identity"
]
