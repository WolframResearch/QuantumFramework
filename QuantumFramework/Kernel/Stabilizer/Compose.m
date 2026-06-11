Package["Wolfram`QuantumFramework`"]

PackageScope[phaseLookup]



(* ============================================================================ *)
(* phaseLookup: precomputed rank-4 SparseArray for symplectic-mod-2             *)
(* multiplication phase corrections.                                            *)
(* Reference: AarGot04 \[Section]3 / Yashin25 Eq 4 / WinPay24b \[Section]3.     *)
(* ============================================================================ *)

phaseLookup := phaseLookup = SparseArray[
    Join[
        Thread[{{1, 2, 2, 1}, {2, 1, 2, 2}, {2, 2, 1, 2}} -> -1],
        Thread[{{1, 2, 2, 2}, {2, 1, 1, 2}, {2, 2, 2, 1}} ->  1]
    ],
    {2, 2, 2, 2}
]


(* ============================================================================ *)
(* Composition: left @ right -- apply `left` after `right`.                     *)
(* Performance: O(n^3) (matrix multiply mod 2). No bit-packed SIMD path; this  *)
(* simulator is not intended to compete with Stim on raw qubit count.          *)
(* ============================================================================ *)

left_PauliStabilizer[right_PauliStabilizer] := Enclose @ Block[{
    n = Max[left["Qudits"], right["Qudits"]],
    first, second, ifacts, x1, z1, p
},
    second = right["PadRight", n];
    first = left["PadRight", n];
    {x1, z1} = Transpose /@ first["Tableau"];
    ifacts = Total[BitAnd[second["X"], second["Z"]]] + Map[
        With[{x1s = Pick[x1, #, 1], z1s = Pick[z1, #, 1]},
            (* A row selecting 0 or 1 generators contributes no consecutive-product *)
            (* phase; guarding it also avoids Rest[{}]/Most[{}] on an all-zero row.  *)
            If[ Length[x1s] <= 1,
                0,
                Total @ Extract[phaseLookup,
                    Thread[Flatten[#] + 1 & /@ {Rest[x1s], Rest[z1s], Most @ FoldList[BitXor, x1s], Most @ FoldList[BitXor, z1s]}]
                ]
            ]
        ] &,
        second["Matrix"]
    ];
    p = Quotient[Mod[ifacts, 4], 2];
    PauliStabilizer[<|
        "Phase" -> Mod[second["Matrix"] . first["Phase"] + second["Phase"] + p, 2],
        "Matrix" -> Mod[second["Matrix"] . first["Matrix"], 2]
    |>]
]


(* ============================================================================ *)
(* Tensor product: block-diagonal merge of destabilizers and stabilizers.       *)
(* Linear in n_a + n_b.                                                         *)
(* ============================================================================ *)

QuantumTensorProduct[a_PauliStabilizer, b_PauliStabilizer] :=
    PauliStabilizer[<|
        "Signs" -> Join[a["DestabilizerSigns"], b["DestabilizerSigns"], a["StabilizerSigns"], b["StabilizerSigns"]],
        "Tableau" -> {
            Join[
                blockDiagonalMatrix[{a["DestabilizerX"], b["DestabilizerX"]}],
                blockDiagonalMatrix[{a["StabilizerX"], b["StabilizerX"]}],
                2
            ],
            Join[
                blockDiagonalMatrix[{a["DestabilizerZ"], b["DestabilizerZ"]}],
                blockDiagonalMatrix[{a["StabilizerZ"], b["StabilizerZ"]}],
                2
            ]
        }
    |>]
