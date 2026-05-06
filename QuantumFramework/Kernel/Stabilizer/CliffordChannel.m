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
(* Phase 8.2 (TODO): Choi-tableau composition via vector-space intersection.    *)
(*                                                                              *)
(* For two CliffordChannels \[CapitalPhi]_{A->B} and \[CapitalPhi]'_{B->C},     *)
(* their composition is encoded by the rows                                     *)
(*   [u_A | u'_C | c \[CirclePlus] c']                                          *)
(* where u_B = u'_B (the inner system bits agree). Finding a basis for these    *)
(* is a Boolean linear-algebra intersection problem (Yashin25 \[Section]3.2).   *)
(*                                                                              *)
(* Stubbed below; Phase 8.2 will fill in the algorithm.                         *)
(* ============================================================================ *)

CliffordChannel::compose = "Composition of CliffordChannels is not yet implemented (Phase 8.2). Use CliffordChannel via PauliStabilizer / QuantumState round-trip for now.";

cc1_CliffordChannel[cc2_CliffordChannel] /; CliffordChannelQ[cc1] && CliffordChannelQ[cc2] := (
    Message[CliffordChannel::compose];
    $Failed
)
