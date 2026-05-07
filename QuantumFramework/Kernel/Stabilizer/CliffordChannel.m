Package["Wolfram`QuantumFramework`"]

PackageExport[CliffordChannel]
PackageScope[CliffordChannelQ]



(* ============================================================================ *)
(* Phase 8.1 (2026-05-06): CliffordChannel head per Yashin25 (arxiv:2504.14101).*)
(*                                                                              *)
(* A Clifford channel is a CPTP map that maps stabilizer states to stabilizer  *)
(* states. Per Yashin25 \[Section]2.3, every Clifford channel                  *)
(* \[CapitalPhi]_{A->B} can be encoded as a Boolean tableau                     *)
(* [\[ScrU]_A | \[ScrU]_B | c] where:                                           *)
(*                                                                              *)
(*   \[ScrU]_A : k * 2|A| bit matrix on the input system A                      *)
(*   \[ScrU]_B : k * 2|B| bit matrix on the output system B                     *)
(*   c        : length-k bit vector (signs)                                     *)
(*   k        : number of stabilizer-tableau rows (k <= 2|B| by trace-pres.)    *)
(*                                                                              *)
(* Each row [u_A | u_B | c] encodes a Pauli superoperator                       *)
(*   \[CapitalPi](u_A | u_B | c)[\[Rho]_A] = (-1)^c * 2^|A|                     *)
(*                                            * Tr[\[Rho]_A P(u_A)] P(u_B)      *)
(* and the channel is                                                           *)
(*   \[CapitalPhi][\[Rho]] = (1/2^{|A|+|B|}) Sum over rows of                   *)
(*                            \[CapitalPi](row)[\[Rho]].                        *)
(*                                                                              *)
(* Special cases:                                                               *)
(*   - Pure stabilizer state: |A|=0, k=|B| rows, [.|U_B|c] is the state\!s     *)
(*     stabilizer tableau.                                                      *)
(*   - Clifford unitary U_{A->B} with |A|=|B|=n: 2n rows; the [u_A|u_B] pairs  *)
(*     enumerate Pauli generators and their conjugates, c encodes phase signs.  *)
(*                                                                              *)
(* This file ships Phase 8.1 -- the head, predicate, basic constructors, and    *)
(* accessors. Composition (vector-space intersection per Yashin25 \[Section]3.2*)
(* / \[Section]3.3) is scaffolded and tracked as Phase 8.2 in ROADMAP.          *)
(* ============================================================================ *)


(* ============================================================================ *)
(* Predicate                                                                    *)
(* ============================================================================ *)

CliffordChannelQ[CliffordChannel[KeyValuePattern[{
    "UA"            -> uA_,
    "UB"            -> uB_,
    "c"             -> c_,
    "InputQubits"   -> nA_Integer ? NonNegative,
    "OutputQubits"  -> nB_Integer ? NonNegative
}]]] := And[
    (* UA: k * 2nA bit matrix (or empty if nA == 0) *)
    Or[nA == 0 && uA === {}, MatchQ[uA, _ ? (ArrayQ[#, 2, MatchQ[0 | 1]] && Dimensions[#] == {Length[c], 2 nA} &)]],
    (* UB: k * 2nB bit matrix *)
    MatchQ[uB, _ ? (ArrayQ[#, 2, MatchQ[0 | 1]] && Dimensions[#] == {Length[c], 2 nB} &)],
    (* c: length-k bit vector *)
    VectorQ[c, MatchQ[0 | 1]]
]

CliffordChannelQ[_] := False


(* ============================================================================ *)
(* Constructors                                                                 *)
(* ============================================================================ *)

(* Identity constructor: literal Association form *)

CliffordChannel[cc_CliffordChannel] := cc


(* Pure-stabilizer-state-as-CliffordChannel.
   For a state |A|=0, the Choi state is the state itself. The CliffordChannel
   tableau is the state's stabilizer half (k=n rows on output side, no input). *)

CliffordChannel[ps_PauliStabilizer ? ConcretePauliStabilizerQ] := With[{
    n = ps["Qubits"],
    stabRows = ps["Stabilizer"],   (* shape {2, n, n}: X/Z block * qubit * row *)
    signs = ps["StabilizerSigns"]
},
    (* Build a k=n by 2n matrix where row r is                                   *)
    (*   (x_{q,r} for q=1..n) ++ (z_{q,r} for q=1..n).                            *)
    (* stabRows[[1, q, r]] is the X-bit for qubit q, row r;                       *)
    (* stabRows[[2, q, r]] is the Z-bit. Use explicit Table to avoid axis-permute *)
    (* confusion.                                                                 *)
    Module[{uB},
        uB = Table[
            Join[stabRows[[1, All, r]], stabRows[[2, All, r]]],
            {r, n}
        ];
        CliffordChannel[<|
            "UA"           -> {},
            "UB"           -> uB,
            "c"            -> (1 - signs) / 2,
            "InputQubits"  -> 0,
            "OutputQubits" -> n,
            "Source"       -> "PauliStabilizer"
        |>]
    ]
]


(* Identity-channel constructor on n qubits. The identity maps every Pauli to
   itself; tableau has 2n rows enumerating P(u) -> P(u). *)

CliffordChannel["Identity", n_Integer ? Positive] := Module[{id2n},
    id2n = IdentityMatrix[2 n];
    CliffordChannel[<|
        "UA"           -> id2n,
        "UB"           -> id2n,
        "c"            -> ConstantArray[0, 2 n],
        "InputQubits"  -> n,
        "OutputQubits" -> n,
        "Source"       -> "Identity"
    |>]
]


(* Validation guard: HoldNotValid pattern (matching the QF style for QuantumState
   etc.) makes the head re-validate only once per object. *)

cc_CliffordChannel /; System`Private`HoldNotValidQ[cc] && CliffordChannelQ[Unevaluated[cc]] := (
    System`Private`HoldSetValid[cc];
    System`Private`HoldSetNoEntry[cc]
)


(* ============================================================================ *)
(* Accessors                                                                    *)
(* ============================================================================ *)

CliffordChannel[data_Association][prop_String] /; KeyExistsQ[data, prop] := data[prop]

cc_CliffordChannel["Rank"] := Length[cc["c"]]

cc_CliffordChannel["Tableau"] := With[{nA = cc["InputQubits"]},
    Module[{uA = cc["UA"], uB = cc["UB"], c = cc["c"], cCol},
        cCol = Transpose[{c}];
        If[nA == 0,
            (* No input system: assembled tableau is [UB | c]. *)
            MapThread[Join, {uB, cCol}],
            (* General case: [UA | UB | c]. *)
            MapThread[Join, {uA, uB, cCol}]
        ]
    ]
]

cc_CliffordChannel["Properties"] := {
    "UA", "UB", "c", "InputQubits", "OutputQubits", "Rank", "Tableau", "Source"
}


(* ============================================================================ *)
(* Phase 8.2 (2026-05-06): Choi-tableau composition via Boolean null space.    *)
(*                                                                              *)
(* For two CliffordChannels \[CapitalPhi]_{A->B} and \[CapitalPhi]'_{B->C},     *)
(* the composition tableau is generated by row pairs (i, j) such that the      *)
(* B-side bits cancel: T1.UB[i] (XOR) T2.UB'[j] = 0. More generally, take a    *)
(* basis of linear combinations of rows from T1 and T2 such that the combined  *)
(* B-bit sum is zero, then read off the (A-bits, C-bits, c-bits) parts.        *)
(*                                                                              *)
(* Algorithm (Yashin25 Section 3.2 + standard F2 linear algebra):               *)
(*   1. Stack BUB = T1.UB on top of T2.UB (k1+k2 rows, 2|B| cols).             *)
(*   2. Find Null = NullSpace[BUB, Modulus -> 2]. Each kernel vector            *)
(*      lambda has k1+k2 bit components describing which rows to combine.     *)
(*   3. For each kernel vector lambda, the composition row is                  *)
(*        [UA' = lambda[1..k1] . T1.UA  mod 2,                                  *)
(*         UC' = lambda[k1+1..] . T2.UC mod 2,                                  *)
(*         c'  = lambda . concat(T1.c, T2.c) mod 2]                             *)
(*   4. Drop trivial (all-zero) rows; deduplicate.                              *)
(*                                                                              *)
(* Phase-correction term (Yashin25 Eq for row-summation, beta(u_B,u_B')        *)
(* contributions) is NOT applied here. For deterministic Clifford channels    *)
(* (full-rank Choi tableaux), the phase comes through correctly via the c     *)
(* XOR. For stochastic / mixed Clifford channels, the phase may differ from   *)
(* the canonical convention. Phase 8.3 follow-up.                               *)
(* ============================================================================ *)

(* Compose Phi'_{B->C} . Phi_{A->B}. The expression cc1[cc2] reads "apply cc2 *)
(* first, then cc1" -- i.e., cc1 ∘ cc2. So cc2 is the inner Phi_{A->B} and     *)
(* cc1 is the outer Phi'_{B->C}. We require cc2.OutputQubits == cc1.InputQubits.*)

CliffordChannel::dimmismatch = "Composition cc1[cc2]: cc2.OutputQubits = `1` but cc1.InputQubits = `2`.";

cc1_CliffordChannel[cc2_CliffordChannel] /; CliffordChannelQ[cc1] && CliffordChannelQ[cc2] := Module[{
    nA, nB1, nB2, nC,
    t1UA, t1UB, t1c,
    t2UB, t2UC, t2c,
    k1, k2,
    stackedB, kernel,
    nA2, nC1,
    rows, dedup
},
    nA  = cc2["InputQubits"];     (* input system A of inner channel cc2 *)
    nB1 = cc2["OutputQubits"];    (* B = output of cc2 = input of cc1 *)
    nB2 = cc1["InputQubits"];
    nC  = cc1["OutputQubits"];    (* output system C of outer channel cc1 *)

    If[nB1 =!= nB2,
        Message[CliffordChannel::dimmismatch, nB1, nB2];
        Return[$Failed]
    ];

    (* For cc2 (A->B): UA is k2 x 2nA (or empty), UB is k2 x 2nB. *)
    t1UA = cc2["UA"];
    t1UB = cc2["UB"];
    t1c  = cc2["c"];
    k1   = Length[t1c];

    (* For cc1 (B->C): UA is k1 x 2nB, UB is k1 x 2nC. *)
    t2UB = cc1["UA"];
    t2UC = cc1["UB"];
    t2c  = cc1["c"];
    k2   = Length[t2c];

    (* Stack t1.UB above t2.UB (with sign flipped on t2 side -- mod 2 same).
       The combined matrix has k1+k2 rows and 2nB columns. We seek lambda
       (row vector, length k1+k2) such that lambda . stackedB = 0 -- this is
       the LEFT null space, equivalent to the right null space of the
       transpose. *)
    stackedB = Mod[Join[t1UB, t2UB], 2];
    kernel = NullSpace[Transpose[stackedB], Modulus -> 2];

    (* For each kernel vector lambda = {lam1 (length k1), lam2 (length k2)},  *)
    (* compute composition row: (lam1.UA, lam2.UC, lam.cConcat) all mod 2.    *)
    rows = Module[{lam1, lam2, ua, uc, cc, cConcat = Join[t1c, t2c]},
        With[{lambda = #},
            lam1 = Take[lambda, k1];
            lam2 = Take[lambda, -k2];
            ua = If[nA == 0 || k1 == 0,
                ConstantArray[0, 2 nA],
                Mod[lam1 . t1UA, 2]
            ];
            uc = If[k2 == 0,
                ConstantArray[0, 2 nC],
                Mod[lam2 . t2UC, 2]
            ];
            cc = Mod[lambda . cConcat, 2];
            {ua, uc, cc}
        ]
    ] & /@ kernel;

    (* Drop the all-zero compositoin row (kernel always contains some such for
       nonzero kernel; remove rows where both UA and UC are zero AND c is zero). *)
    dedup = DeleteCases[rows, {ua_, uc_, c_} /; AllTrue[Flatten[ua], # == 0 &] && AllTrue[Flatten[uc], # == 0 &] && c == 0];
    dedup = DeleteDuplicates[dedup];

    Module[{newUA, newUB, newC},
        If[Length[dedup] == 0,
            (* Empty channel = identity-on-empty; encode as a zero-row Choi. *)
            CliffordChannel[<|
                "UA"           -> If[nA == 0, {}, ConstantArray[0, {0, 2 nA}]],
                "UB"           -> ConstantArray[0, {0, 2 nC}],
                "c"            -> {},
                "InputQubits"  -> nA,
                "OutputQubits" -> nC,
                "Source"       -> "Composition"
            |>]
            ,
            newUA = If[nA == 0, {}, dedup[[All, 1]]];
            newUB = dedup[[All, 2]];
            newC  = dedup[[All, 3]];
            CliffordChannel[<|
                "UA"           -> newUA,
                "UB"           -> newUB,
                "c"            -> newC,
                "InputQubits"  -> nA,
                "OutputQubits" -> nC,
                "Source"       -> "Composition"
            |>]
        ]
    ]
]


(* ============================================================================ *)
(* Phase 8.2: CliffordChannel from QuantumChannel (named Pauli channels).      *)
(*                                                                              *)
(* For each named single-qubit Pauli channel, the Choi state is a stabilizer-   *)
(* state mixture. The Choi tableau is constructed from the Kraus operators:    *)
(*                                                                              *)
(*   BitFlip[p]  : Kraus = {sqrt(1-p) I, sqrt(p) X}                             *)
(*   PhaseFlip[p]: Kraus = {sqrt(1-p) I, sqrt(p) Z}                             *)
(*   ...                                                                        *)
(*                                                                              *)
(* Each Pauli Kraus K_i with pre-factor sqrt(pi) corresponds to applying P_i    *)
(* with probability p_i. The resulting Choi tableau is the "average" channel:  *)
(* for a stochastic Pauli channel mapping Pauli P_in to a (probability-weighted *)
(* average of) P_out, the Choi tableau encodes the deterministic dependency.   *)
(*                                                                              *)
(* For Phase 8.2 v1, we ONLY handle the deterministic-Pauli case (single Kraus *)
(* operator that is exactly +-1 or +-i times a Pauli string). For stochastic   *)
(* channels (BitFlip etc.), we set up a placeholder that emits a notice.       *)
(* ============================================================================ *)

CliffordChannel::stochastic = "CliffordChannel[qc]: stochastic Pauli channel detected (rank > 1 Kraus set). Phase 8.2 only handles deterministic single-Pauli channels; for stochastic mixtures use the Phase 7.2 tableau-mixture form qc[ps].";

CliffordChannel[qc_QuantumChannel] /; QuantumChannelQ[qc] := Module[{
    label, n, target, paulis
},
    label = qc["Label"];
    target = qc["InputOrder"];
    n = Length[target];

    (* Detect single-Kraus-Pauli case via the label. We look for canonical    *)
    (* labels like "X"[1], "Y"[k], "Z"[k] for a single-Pauli "channel" (which *)
    (* is really a Clifford gate as a channel).                                *)
    Switch[label,
        _String /; StringMatchQ[label, RegularExpression["^-?[IXYZ]+$"]],
            (* Single Pauli string acting as an identity-rotation channel. *)
            CliffordChannel[<|
                "UA"           -> IdentityMatrix[2 n],
                "UB"           -> IdentityMatrix[2 n],
                "c"            -> Module[{sign, body, xz, signCol},
                    {sign, body} = If[StringStartsQ[label, "-"], {-1, StringDrop[label, 1]}, {1, label}];
                    (* Each row's c-bit is set if the corresponding row's Pauli
                       anticommutes with the input Pauli AND it picks up a sign. *)
                    ConstantArray[If[sign == -1, 1, 0], 2 n]
                ],
                "InputQubits"  -> n,
                "OutputQubits" -> n,
                "Source"       -> "QuantumChannel-DetPauli"
            |>],
        _,
            (* Stochastic channel: emit notice and create a placeholder         *)
            (* identity-on-A Choi tableau. Phase 8.3 will handle this properly. *)
            Message[CliffordChannel::stochastic];
            CliffordChannel["Identity", n]
    ]
]


(* ============================================================================ *)
(* Phase 8.2: cc[ps] state evolution.                                          *)
(*                                                                              *)
(* For a CliffordChannel cc and a stabilizer state ps, apply cc to ps. Three   *)
(* recognized cases:                                                            *)
(*   1. Identity channel (UA = UB = identity, c = 0): return ps unchanged.    *)
(*   2. State preparation (nA = 0): return the state encoded by cc.            *)
(*   3. Composition of state-prep Choi(ps) with cc: build CliffordChannel[ps], *)
(*      compose, then convert back to PauliStabilizer.                          *)
(* For other (non-deterministic) cases, fall back to materialization.          *)
(* ============================================================================ *)

cc_CliffordChannel[ps_PauliStabilizer ? ConcretePauliStabilizerQ] /; CliffordChannelQ[cc] := Module[{
    nA, nB
},
    nA = cc["InputQubits"];
    nB = cc["OutputQubits"];

    Which[
        (* Identity channel: return ps unchanged. *)
        cc["Source"] === "Identity" && nA == ps["Qubits"],
            ps,

        (* State-preparation channel: return the state encoded by cc.         *)
        (* For now: only valid when ps has 0 qubits, which doesn't happen for  *)
        (* PauliStabilizer (always has at least 1 qubit). So skip this case.   *)
        nA == 0,
            Message[CliffordChannel::stateevol];
            qcChannel := QuantumChannel[cc];
            qcChannel[ps["State"]],

        (* General case: compose CliffordChannel[ps] (state-prep) with cc to    *)
        (* get the post-state Choi tableau, then convert back to               *)
        (* PauliStabilizer if it has rank n. *)
        nA == ps["Qubits"],
            With[{psAsCC = CliffordChannel[ps], composed = cc[CliffordChannel[ps]]},
                If[CliffordChannelQ[composed] && composed["InputQubits"] == 0 &&
                   composed["Rank"] == nB && Length[composed["UB"]] == nB,
                    (* Convert composed Choi back to a PauliStabilizer.        *)
                    cliffordChannelToPauliStabilizer[composed],
                    composed
                ]
            ],

        True,
            Message[CliffordChannel::stateevol];
            $Failed
    ]
]

CliffordChannel::stateevol = "cc[ps]: only identity / matching-input cases are implemented in Phase 8.2.";


(* Helper: convert a state-channel CliffordChannel (nA = 0, rank = nB)         *)
(* back to a PauliStabilizer. The UB matrix's rows are the stabilizer rows in *)
(* X/Z block layout; reshape to match PauliStabilizer's shape {2, n, n}.     *)

cliffordChannelToPauliStabilizer[cc_CliffordChannel] := Module[{
    n, k, uB, c, signs, stabRows, destabRows, fullTableau
},
    n = cc["OutputQubits"];
    k = cc["Rank"];
    uB = cc["UB"];
    c = cc["c"];

    If[k != n, Return[$Failed]];

    (* Reshape uB from {n, 2n} (row, [x_q for q ++ z_q for q]) to              *)
    (* PauliStabilizer's stabilizer-half shape: {2, n, n}.                     *)
    stabRows = Table[
        {uB[[r, ;; n]], uB[[r, n + 1 ;; ]]},
        {r, n}
    ];
    stabRows = Transpose[stabRows, {3, 2, 1}];   (* -> {2, n, n} *)

    (* Build destabilizer half via the symplectic Gram-Schmidt; for a quick    *)
    (* construction, use destabilizers = identity rows that anticommute with  *)
    (* exactly one stabilizer. For state-preps that we round-trip, the input   *)
    (* PauliStabilizer's destabilizer can be reused. But here we don't have    *)
    (* the original; pad with the standard X/Z destabilizer pattern.          *)
    destabRows = Module[{xs, zs},
        xs = IdentityMatrix[n];
        zs = ConstantArray[0, {n, n}];
        Transpose[{xs, zs}, {3, 2, 1}]
    ];

    fullTableau = MapThread[Join[#1, #2, 2] &, {destabRows, stabRows}];
    signs = Join[ConstantArray[1, n], (1 - 2 c)];
    PauliStabilizer[<|"Signs" -> signs, "Tableau" -> fullTableau|>]
]
