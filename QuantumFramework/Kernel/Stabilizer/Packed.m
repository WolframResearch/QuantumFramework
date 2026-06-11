Package["Wolfram`QuantumFramework`"]

PackageScope[$chunkWidth]
PackageScope[psPackedQ]
PackageScope[psConcreteFastQ]
PackageScope[psGetPacked]
PackageScope[psFromPacked]
PackageScope[packChunks]
PackageScope[unpackChunks]
PackageScope[tableauFromPacked]
PackageScope[packedGateH]
PackageScope[packedGateS]
PackageScope[packedGateSdg]
PackageScope[packedGateCNOT]
PackageScope[packedGateSWAP]
PackageScope[packedGateX]
PackageScope[packedGateY]
PackageScope[packedGateZ]



(* ============================================================================ *)
(* Bit-packed fast path for the AG tableau.                                      *)
(*                                                                              *)
(* The canonical representation is the rank-3 integer array  ps["Tableau"]      *)
(* of dimensions {2, n_qubits, 2n_rows} together with the length-2n  ps["Signs"]*)
(* list.  For the Clifford hot path that array is mutated through per-gate       *)
(* property dispatch, which rebuilds the whole O(n^2) array and re-scans the     *)
(* sign list on every gate.                                                      *)
(*                                                                              *)
(* The packed view stores, for each of the 2n tableau rows, the qubit bits of    *)
(* the X part and of the Z part packed into machine integers, 62 qubits per      *)
(* chunk (bit b of chunk c = qubit  c*62 + b).  Storage is chunk-major: each of  *)
(* "PackedX" / "PackedZ" is a length-(#chunks) list whose c-th entry is a        *)
(* length-2n *packed array of machine integers* (one per tableau row).  Because  *)
(* a Clifford generator only touches the chunk column(s) of its qubit(s), each   *)
(* gate is a handful of Listable bitwise ops over one or two length-2n packed    *)
(* machine-int vectors -- O(n) vectorized work per gate instead of O(n^2), with  *)
(* no arbitrary-precision overhead (62-bit chunks stay inside Int64).            *)
(*                                                                              *)
(* Packed PauliStabilizer objects carry the keys                                *)
(*     "PackedX" -> {chunk vectors}, "PackedZ" -> {chunk vectors},              *)
(*     "Signs"   -> {-1|1 ...},      "Qubits"  -> n                             *)
(* and no "Tableau" key.  The rank-3 array is reconstructed lazily by the        *)
(* ps["Tableau"] handler below, so every downstream property (State, Measure,    *)
(* Formatting, ...) keeps working unchanged.                                     *)
(*                                                                              *)
(* References: AarGot04 \[Section]3 (canonical update rules); external-packages- *)
(* audit \[Section]4 (the multiplication kernel; 62-bit machine words in place   *)
(* of UInt64), \[Section]14 (correctness traps).                                 *)
(* ============================================================================ *)

$chunkWidth = 62;

(* {1-based chunk index, 0-based bit position} of qubit j (1-based). *)
qubitLoc[j_Integer] := {Quotient[j - 1, $chunkWidth] + 1, Mod[j - 1, $chunkWidth]}

chunkCount[n_Integer] := Ceiling[n / $chunkWidth]

(* qubit count carried by chunk c of an n-qubit register *)
chunkQubits[n_Integer, c_Integer] := Min[c $chunkWidth, n] - (c - 1) $chunkWidth

(* Pack a {n_qubits, 2n_rows} {0,1} matrix into chunk-major form: a list whose   *)
(* c-th entry is the length-2n machine-int vector for qubit-chunk c.             *)
packChunks[m_, n_Integer] := Table[
    Map[FromDigits[Reverse[#], 2] &, Transpose[m[[(c - 1) $chunkWidth + 1 ;; Min[c $chunkWidth, n]]]]],
    {c, chunkCount[n]}
]

(* Inverse of packChunks: chunk-major form -> {n_qubits, 2n_rows} {0,1}. *)
unpackChunks[chunks_, n_Integer] := Transpose[
    MapThread[Join, MapIndexed[Map[Reverse, IntegerDigits[#1, 2, chunkQubits[n, #2[[1]]]]] &, chunks]]
]

tableauFromPacked[n_, px_, pz_] := {unpackChunks[px, n], unpackChunks[pz, n]}


(* ---- predicates ---- *)

psPackedQ[ps_PauliStabilizer] := KeyExistsQ[First[ps], "PackedX"]
psPackedQ[_] := False

(* True when the packed bitwise path may be taken: either the object is already *)
(* packed, or it carries a concrete binary {0,1} tableau with concrete {-1,1}    *)
(* signs.  Symbolic-phase states (non-Clifford P/T residues, StabilizerFrame     *)
(* mixtures) fail this test and fall through to the canonical array path.        *)
psConcreteFastQ[ps_PauliStabilizer] := With[{a = First[ps]},
    KeyExistsQ[a, "PackedX"] ||
    (
        KeyExistsQ[a, "Tableau"] && KeyExistsQ[a, "Signs"] &&
        MatchQ[a["Signs"], {(-1 | 1) ...}] &&
        With[{t = a["Tableau"]},
            ArrayQ[t, 3, IntegerQ] && Min[t] >= 0 && Max[t] <= 1
        ]
    )
]
psConcreteFastQ[_] := False


(* ---- conversion to / from the {n, packedX, packedZ, signs} tuple ---- *)

psGetPacked[ps_PauliStabilizer] := With[{a = First[ps]},
    If[ KeyExistsQ[a, "PackedX"],
        {a["Qubits"], a["PackedX"], a["PackedZ"], a["Signs"]},
        With[{t = a["Tableau"], n = Dimensions[a["Tableau"]][[2]]},
            {n, packChunks[t[[1]], n], packChunks[t[[2]], n], a["Signs"]}
        ]
    ]
]

psFromPacked[{n_, px_, pz_, s_}] :=
    PauliStabilizer[<|"PackedX" -> px, "PackedZ" -> pz, "Signs" -> s, "Qubits" -> n|>]


(* ---- lazy rank-3 reconstruction for packed objects ---- *)

ps_PauliStabilizer["Tableau"] /; KeyExistsQ[First[ps], "PackedX"] := With[{a = First[ps]},
    tableauFromPacked[a["Qubits"], a["PackedX"], a["PackedZ"]]
]


(* ============================================================================ *)
(* Packed Clifford-generator kernels.  Each takes and returns the tuple          *)
(* {n, packedX, packedZ, signs} and mirrors the canonical rule in GateUpdates.m  *)
(* bit-for-bit (verified by the equivalence tests in PauliStabilizer.wlt).       *)
(* A single-qubit gate touches the one chunk column holding its qubit; CNOT/SWAP *)
(* touch the (up to two) columns of qubits j and k.  Sign of a masked column is  *)
(* the per-row 0/1 flip indicator.                                               *)
(* ============================================================================ *)

(* shift a masked column from bit position `from` to bit position `to`. *)
shiftCol[v_, from_, to_] := With[{d = to - from}, If[d >= 0, BitShiftLeft[v, d], BitShiftRight[v, -d]]]

(* H_j: swap X[j] <-> Z[j]; flip sign where x[j] == z[j] == 1. *)
packedGateH[{n_, px_, pz_, s_}, j_Integer] := With[{loc = qubitLoc[j]},
    With[{c = loc[[1]], m = BitShiftLeft[1, loc[[2]]]},
        With[{xc = px[[c]], zc = pz[[c]]},
            With[{bx = BitAnd[xc, m], bz = BitAnd[zc, m]},
                {n, ReplacePart[px, c -> xc - bx + bz], ReplacePart[pz, c -> zc - bz + bx],
                 s (1 - 2 Sign[BitAnd[bx, bz]])}
            ]
        ]
    ]
]

(* S_j: Z[j] := Z[j] XOR X[j]; flip sign where x[j] == z[j] == 1. *)
packedGateS[{n_, px_, pz_, s_}, j_Integer] := With[{loc = qubitLoc[j]},
    With[{c = loc[[1]], m = BitShiftLeft[1, loc[[2]]]},
        With[{xc = px[[c]], zc = pz[[c]]},
            With[{bx = BitAnd[xc, m], bz = BitAnd[zc, m]},
                {n, px, ReplacePart[pz, c -> BitXor[zc, bx]], s (1 - 2 Sign[BitAnd[bx, bz]])}
            ]
        ]
    ]
]

(* S^\[Dagger]_j: same Z-update; flip sign where x[j] == 1 and z[j] == 0. *)
packedGateSdg[{n_, px_, pz_, s_}, j_Integer] := With[{loc = qubitLoc[j]},
    With[{c = loc[[1]], m = BitShiftLeft[1, loc[[2]]]},
        With[{xc = px[[c]], zc = pz[[c]]},
            With[{bx = BitAnd[xc, m], bz = BitAnd[zc, m]},
                {n, px, ReplacePart[pz, c -> BitXor[zc, bx]], s (1 - 2 Sign[BitAnd[bx, BitXor[bz, m]]])}
            ]
        ]
    ]
]

(* CNOT j -> k: X[k] ^= X[j]; Z[j] ^= Z[k];                                      *)
(* flip sign where x[j] == z[k] == 1 and x[k] == z[j]  (AG sign rule).           *)
packedGateCNOT[{n_, px_, pz_, s_}, j_Integer, k_Integer] :=
    With[{lj = qubitLoc[j], lk = qubitLoc[k]},
        With[{cj = lj[[1]], bj = lj[[2]], ck = lk[[1]], bk = lk[[2]]},
            With[{mj = BitShiftLeft[1, lj[[2]]], mk = BitShiftLeft[1, lk[[2]]],
                  xcj = px[[lj[[1]]]], xck = px[[lk[[1]]]], zcj = pz[[lj[[1]]]], zck = pz[[lk[[1]]]]},
                With[{bxj = BitAnd[xcj, mj], bzk = BitAnd[zck, mk]},
                    {
                        n,
                        ReplacePart[px, ck -> BitXor[xck, shiftCol[bxj, bj, bk]]],
                        ReplacePart[pz, cj -> BitXor[zcj, shiftCol[bzk, bk, bj]]],
                        s (1 - 2 (Sign[bxj] Sign[bzk] (1 - BitXor[Sign[BitAnd[xck, mk]], Sign[BitAnd[zcj, mj]]])))
                    }
                ]
            ]
        ]
    ]

(* SWAP j <-> k: exchange bit j and bit k in every X row and every Z row.        *)
(* Signs unchanged (canonical SWAP is a column permutation).                     *)
swapBitCols[p_, cj_, bj_, ck_, bk_] := With[{mj = BitShiftLeft[1, bj], mk = BitShiftLeft[1, bk]},
    If[ cj === ck,
        With[{col = p[[cj]]},
            ReplacePart[p, cj -> BitXor[col, (mj + mk) BitXor[Sign[BitAnd[col, mj]], Sign[BitAnd[col, mk]]]]]
        ],
        With[{colj = p[[cj]], colk = p[[ck]]},
            With[{differ = BitXor[Sign[BitAnd[p[[cj]], mj]], Sign[BitAnd[p[[ck]], mk]]]},
                ReplacePart[p, {cj -> BitXor[colj, differ mj], ck -> BitXor[colk, differ mk]}]
            ]
        ]
    ]
]
packedGateSWAP[{n_, px_, pz_, s_}, j_Integer, k_Integer] :=
    With[{lj = qubitLoc[j], lk = qubitLoc[k]},
        {n, swapBitCols[px, lj[[1]], lj[[2]], lk[[1]], lk[[2]]],
            swapBitCols[pz, lj[[1]], lj[[2]], lk[[1]], lk[[2]]], s}
    ]

(* X_j: phase XOR with Z[j]  ->  flip sign where z[j] == 1. *)
packedGateX[{n_, px_, pz_, s_}, j_Integer] := With[{loc = qubitLoc[j]},
    {n, px, pz, s (1 - 2 Sign[BitAnd[pz[[loc[[1]]]], BitShiftLeft[1, loc[[2]]]]])}
]

(* Z_j: phase XOR with X[j]  ->  flip sign where x[j] == 1. *)
packedGateZ[{n_, px_, pz_, s_}, j_Integer] := With[{loc = qubitLoc[j]},
    {n, px, pz, s (1 - 2 Sign[BitAnd[px[[loc[[1]]]], BitShiftLeft[1, loc[[2]]]]])}
]

(* Y_j: phase XOR with X[j] XOR Z[j]  ->  flip sign where x[j] != z[j]. *)
packedGateY[{n_, px_, pz_, s_}, j_Integer] := With[{loc = qubitLoc[j]},
    With[{c = loc[[1]], m = BitShiftLeft[1, loc[[2]]]},
        {n, px, pz, s (1 - 2 BitXor[Sign[BitAnd[px[[c]], m]], Sign[BitAnd[pz[[c]], m]]])}
    ]
]
