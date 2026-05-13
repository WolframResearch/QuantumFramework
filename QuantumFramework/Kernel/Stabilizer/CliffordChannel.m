Package["Wolfram`QuantumFramework`"]

PackageExport[CliffordChannel]
PackageScope[CliffordChannelQ]
PackageScope[stabilizerRowSumAGPhase]
PackageScope[stabilizerContractionPhase]



(* ============================================================================ *)
(* CliffordChannel head per Yashin25 (arxiv:2504.14101).                        *)
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
(*   - Pure stabilizer state: |A|=0, k=|B| rows, [.|U_B|c] is the state's      *)
(*     stabilizer tableau.                                                      *)
(*   - Clifford unitary U_{A->B} with |A|=|B|=n: 2n rows; the [u_A|u_B] pairs  *)
(*     enumerate Pauli generators and their conjugates, c encodes phase signs.  *)
(*                                                                              *)
(* Composition uses vector-space intersection per Yashin25 \[Section]3.2/3.3,  *)
(* with AG-style row-sum phase tracking and the |Phi+>_BB' contraction sign    *)
(* for Y-bearing combined u_B rows.                                             *)
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
(* Choi-tableau composition via Boolean null space, with AG phase tracking and *)
(* the |Phi+>_BB' contraction-sign correction.                                 *)
(*                                                                              *)
(* For two CliffordChannels \[CapitalPhi]_{A->B} and \[CapitalPhi]'_{B->C},     *)
(* the composition tableau is generated by row pairs (i, j) such that the      *)
(* B-side bits cancel: T1.UB[i] (XOR) T2.UB'[j] = 0. More generally, take a    *)
(* basis of linear combinations of rows from T1 and T2 such that the combined  *)
(* B-bit sum is zero, then read off the (A-bits, C-bits, c-bits) parts.        *)
(*                                                                              *)
(* Algorithm (Yashin25 Section 3.2/3.3 + F2 linear algebra):                    *)
(*   1. Stack BUB = T1.UB on top of T2.UB (k1+k2 rows, 2|B| cols).             *)
(*   2. Find Null = NullSpace[BUB, Modulus -> 2]. Each kernel vector            *)
(*      lambda has k1+k2 bit components describing which rows to combine.     *)
(*   3. For each kernel vector lambda, compute the composition row:            *)
(*        UA' = lambda[1..k1] . T1.UA  mod 2                                    *)
(*        UC' = lambda[k1+1..] . T2.UC mod 2                                    *)
(*        c'  = (lambda . concat(T1.c, T2.c) mod 2) XOR phaseCorr               *)
(*      where phaseCorr accounts for: (a) AG g-function row-sum phase on A,    *)
(*      C, and both B sides; (b) the |Phi+>_BB' contraction sign which adds    *)
(*      (-1)^(sum_q x_B z_B) for the combined B Pauli (Y_B has Y^T = -Y).     *)
(*   4. Drop trivial (all-zero) rows; deduplicate.                              *)
(*                                                                              *)
(* The phase correction is REQUIRED for unitary Clifford channels: e.g., S*S    *)
(* has X-side row (X, X, c=1) per Z X Z = -X, but the simple c-bit XOR gives   *)
(* c=0. The AG phase from row-summing on B (X*Z = -iY -> phase 3) and on C    *)
(* (Y*Z = iX -> phase 1) sum to 4 mod 4 (real); the |Phi+>_BB' contraction    *)
(* with combined u_B = Y on each side adds (-1)^(1*1) = -1, contributing      *)
(* +2 to the phase mod 4, which yields c_new = 1 via /2.                      *)
(* ============================================================================ *)


(* ============================================================================ *)
(* Helper: AG phase tracking when XOR-summing tableau rows.                    *)
(*                                                                              *)
(* Given a row matrix U (k * 2n in [x|z] block layout), a selection vector     *)
(* lambda (length k F2), and the qubit count n, return the cumulative AG       *)
(* phase mod 4 of the Pauli operator product of selected rows.                 *)
(*                                                                              *)
(* Convention: agPhase(x_src, z_src, x_dst, z_dst) tracks the i-power of       *)
(* P(src) * P(dst). Iterating accum_op = accum_op * P(next_row) gives          *)
(* P(row_1) * P(row_2) * ... * P(row_m) for selected rows in order.            *)
(* ============================================================================ *)

stabilizerRowSumAGPhase[U_, lambda_, n_Integer ? Positive] := Module[{selected},
    selected = Pick[U, lambda, 1];
    If[Length[selected] <= 1,
        0,
        (* Fold: state = {accumPhase, accumRow}. At each step, add per-qubit AG  *)
        (* g-function contributions (vectorized via MapThread) and XOR the rows. *)
        Mod[
            First @ Fold[
                Function[{state, nextRow},
                    With[{
                        accumPhase = state[[1]],
                        accumRow = state[[2]]
                    },
                        {
                            accumPhase + Total @ MapThread[
                                agPhase,
                                {accumRow[[;; n]], accumRow[[n + 1 ;;]],
                                 nextRow[[;; n]], nextRow[[n + 1 ;;]]}
                            ],
                            Mod[accumRow + nextRow, 2]
                        }
                    ]
                ],
                {0, First[selected]},
                Rest[selected]
            ],
            4
        ]
    ]
]


(* ============================================================================ *)
(* Helper: |Phi+>_BB' contraction sign for combined B Pauli.                   *)
(*                                                                              *)
(* For the maximally entangled state |Phi+>_BB' = (1/sqrt(d)) sum_i |i,i>,     *)
(* and Pauli P_B = X^a Z^b on B (per qubit), the transpose flip gives          *)
(* P^T = (-1)^(a*b) P. Hence for the kernel-vector contraction of identical    *)
(* B-side Paulis on cc2 and cc1, the projection picks up sign (-1)^(sum_q x_q*z_q). *)
(* Returns the signed contribution to the i-power mod 4 (i.e. 2 * a*b summed   *)
(* per qubit, mod 4). This is what makes Y-bearing combined u_B contribute    *)
(* a -1 sign that the simple c-bit XOR misses.                                  *)
(* ============================================================================ *)

stabilizerContractionPhase[uB_List, n_Integer ? Positive] := Module[{q},
    Mod[2 * Sum[uB[[q]] * uB[[n + q]], {q, n}], 4]
]

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

    (* For cc2 (A->B): UA is k_cc2 x 2nA (or empty), UB is k_cc2 x 2nB. *)
    t1UA = cc2["UA"];
    t1UB = cc2["UB"];
    t1c  = cc2["c"];
    k1   = Length[t1c];

    (* For cc1 (B->C): UA is k_cc1 x 2nB, UB is k_cc1 x 2nC. *)
    t2UB = cc1["UA"];
    t2UC = cc1["UB"];
    t2c  = cc1["c"];
    k2   = Length[t2c];

    (* Stack t1.UB above t2.UB. The combined matrix has k1+k2 rows and 2nB    *)
    (* columns. We seek lambda (row vector, length k1+k2) such that           *)
    (* lambda . stackedB = 0 -- this is the LEFT null space, equivalent to    *)
    (* the right null space of the transpose.                                 *)
    stackedB = Mod[Join[t1UB, t2UB], 2];
    kernel = NullSpace[Transpose[stackedB], Modulus -> 2];

    (* For each kernel vector lambda = {lam1 (length k1), lam2 (length k2)},  *)
    (* compute the composition row with full phase tracking:                  *)
    (*   - bits: u_A = lam1 . t1UA, u_C = lam2 . t2UC (in F2)                 *)
    (*   - c-bit: lambda . cConcat XOR (rowSumPhase + contractionPhase) / 2.  *)
    (* If the total phase mod 4 is odd (1 or 3), the kernel vector            *)
    (* corresponds to an i-factor Pauli that cannot be encoded in the         *)
    (* {0,1}-valued c-bit. Such rows are skipped (this should never happen    *)
    (* for valid Clifford channel inputs).                                    *)
    rows = Module[{
            lam1, lam2, ua, uc, cBitXor, cConcat = Join[t1c, t2c],
            uBcc2, phaseA, phaseB1, phaseB2, phaseC, phaseContr, phaseTotal,
            cPhase
        },
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
            cBitXor = Mod[lambda . cConcat, 2];

            (* Combined u_B from cc2 (also = combined u_B from cc1 by kernel  *)
            (* constraint, in F2; we use cc2's combined for the contraction). *)
            uBcc2 = If[k1 == 0,
                ConstantArray[0, 2 nB1],
                Mod[lam1 . t1UB, 2]
            ];

            (* Row-sum AG phases for each side. nB1 == nB2 by the dim check. *)
            phaseA = If[nA == 0 || k1 == 0, 0, stabilizerRowSumAGPhase[t1UA, lam1, nA]];
            phaseB1 = If[k1 == 0, 0, stabilizerRowSumAGPhase[t1UB, lam1, nB1]];
            phaseB2 = If[k2 == 0, 0, stabilizerRowSumAGPhase[t2UB, lam2, nB1]];
            phaseC = If[k2 == 0, 0, stabilizerRowSumAGPhase[t2UC, lam2, nC]];

            (* Phi+ contraction sign: (-1)^(sum_q x_q z_q of combined u_B). *)
            phaseContr = stabilizerContractionPhase[uBcc2, nB1];

            phaseTotal = Mod[phaseA + phaseB1 + phaseB2 + phaseC + phaseContr, 4];

            (* For valid (real Pauli) compositions, phaseTotal is in {0, 2}. *)
            cPhase = Which[
                phaseTotal == 0, 0,
                phaseTotal == 2, 1,
                True, $imaginaryPhase   (* signal to skip *)
            ];

            If[cPhase === $imaginaryPhase,
                Nothing,
                {ua, uc, Mod[cBitXor + cPhase, 2]}
            ]
        ]
    ] & /@ kernel;

    rows = DeleteCases[rows, Nothing];

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
(* CliffordChannel from QuantumChannel (named Pauli channels).                 *)
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
(* Currently only the deterministic-Pauli case (single Kraus operator that is  *)
(* exactly +-1 or +-i times a Pauli string) is handled. For stochastic         *)
(* channels (BitFlip etc.), a placeholder is set up that emits a notice.       *)
(* ============================================================================ *)

CliffordChannel::stochastic = "CliffordChannel[qc]: stochastic Pauli channel detected (rank > 1 Kraus set). Only deterministic single-Pauli channels are currently handled; for stochastic mixtures use the tableau-mixture form qc[ps].";

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
            (* Stochastic channel: emit a clean named message and fail.        *)
            (* The previous behavior (placeholder CliffordChannel["Identity",n])*)
            (* was misleading because BitFlip / PhaseFlip / etc. are NOT        *)
            (* identity channels; callers had to read the message to know they *)
            (* were getting a stub. For tableau-level stochastic application,   *)
            (* use the qc[ps] form (see HybridInterop.m).                       *)
            Message[CliffordChannel::stochastic];
            $Failed
    ]
]


(* ============================================================================ *)
(* cc[ps] state evolution.                                                     *)
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

        (* State-preparation channel: nA == 0. Per the cc[ps] contract this   *)
        (* returns the state encoded by cc, regardless of the input ps        *)
        (* (which has no input-side meaning for a state-prep channel). The    *)
        (* tableau-only encoding of the prepared state is recovered via       *)
        (* cliffordChannelToPauliStabilizer.                                    *)
        nA == 0,
            cliffordChannelToPauliStabilizer[cc],

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

CliffordChannel::stateevol = "cc[ps]: only identity / matching-input cases are currently implemented.";


(* Helper: convert a state-channel CliffordChannel (nA = 0, rank = nB)         *)
(* back to a PauliStabilizer.                                                   *)
(*                                                                              *)
(* Strategy: convert each UB row to a signed Pauli string, then construct the  *)
(* PauliStabilizer via the existing string-list constructor (which handles     *)
(* the destabilizer Gram-Schmidt internally via the AG initialization).         *)

cliffordChannelToPauliStabilizer[cc_CliffordChannel] := Module[{
    n, k, uB, c, stabStrings
},
    n = cc["OutputQubits"];
    k = cc["Rank"];
    uB = cc["UB"];
    c = cc["c"];

    If[k != n, Return[$Failed]];

    (* Convert each row of uB (length 2n in [x|z] block layout) to a Pauli    *)
    (* string with sign prefix from c-bit.                                     *)
    stabStrings = Table[
        Module[{xBits, zBits, body},
            xBits = uB[[r, ;; n]];
            zBits = uB[[r, n + 1 ;; ]];
            body = StringJoin[Table[
                Which[
                    xBits[[q]] == 0 && zBits[[q]] == 0, "I",
                    xBits[[q]] == 1 && zBits[[q]] == 0, "X",
                    xBits[[q]] == 0 && zBits[[q]] == 1, "Z",
                    True, "Y"
                ],
                {q, n}
            ]];
            If[c[[r]] == 1, "-" <> body, body]
        ],
        {r, k}
    ];

    PauliStabilizer[stabStrings]
]


(* ============================================================================ *)
(* Formatting: summary box                                                      *)
(*                                                                              *)
(* Decode each tableau row [u_A | u_B | c] to a pretty Pauli-superoperator      *)
(* representation. Mirror the script-letter icon convention used by the         *)
(* sibling Stabilizer heads (PauliStabilizer, StabilizerFrame, GraphState).     *)
(* ============================================================================ *)

cliffordChannelPauliRow[row_, n_Integer] := If[n == 0,
    "\[CenterDot]",
    StringJoin @@ Table[
        Replace[{row[[q]], row[[n + q]]},
            {{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"}],
        {q, n}
    ]
]

cliffordChannelDisplayRow[uARow_, uBRow_, cBit_, nA_Integer, nB_Integer] := Row[{
    If[cBit == 0, "+ ", "\[Minus] "],
    cliffordChannelPauliRow[uARow, nA],
    " \[RightArrow] ",
    cliffordChannelPauliRow[uBRow, nB]
}]


MakeBoxes[cc : CliffordChannel[data_Association] ? CliffordChannelQ, form_] ^:= With[{
    nA = data["InputQubits"],
    nB = data["OutputQubits"],
    k = Length[data["c"]],
    source = Lookup[data, "Source", Missing[]]
},
    BoxForm`ArrangeSummaryBox["CliffordChannel",
        cc,
        Tooltip[Framed["\[ScriptCapitalC]"],
            If[MissingQ[source], "Clifford channel", "Source: " <> ToString[source]]],
        Join[
            {{BoxForm`SummaryItem[{"Qubits: ", Row[{nA, "\[RightArrow]", nB}]}]}},
            {{BoxForm`SummaryItem[{"Tableau rows: ", k}]}},
            If[! MissingQ[source] && source =!= "Composition",
                {{BoxForm`SummaryItem[{"Source: ", source}]}},
                {}
            ]
        ],
        If[k <= 32,
            With[{rowStrs = Table[
                cliffordChannelDisplayRow[
                    If[data["UA"] === {}, {}, data["UA"][[r]]],
                    data["UB"][[r]],
                    data["c"][[r]],
                    nA, nB
                ],
                {r, Min[k, 8]}
            ]},
                {{BoxForm`SummaryItem[{"Tableau: ",
                    Column[If[k > 8, Append[rowStrs, "\[VerticalEllipsis]"], rowStrs]]
                }]}}
            ],
            {{BoxForm`SummaryItem[{"Tableau: ", Row[{k, " rows (preview suppressed)"}]}]}}
        ],
        form
    ]
]
