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
