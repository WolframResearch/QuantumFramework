Package["Wolfram`QuantumFramework`"]

PackageScope[agPhase]
PackageScope[rowsum]
PackageScope[packedRowSumPhase]
PackageScope[packedMeasureZ]



(* ============================================================================ *)
(* AG phase function and rowsum primitive.                                      *)
(* Reference: AarGot04 \[Section]3, Lemma 3.                                     *)
(* ============================================================================ *)

(* The AG phase increment g(x1, z1, x2, z2) for multiplying tableau rows. *)
agPhase[x1_, z1_, x2_, z2_] := Which[
    x1 == z1 == 0, 0,
    x1 == z1 == 1, z2 - x2,
    x1 == 1 && z1 == 0, z2 (2 x2 - 1),
    True, x2 (1 - 2 z2)
]

(* rowsum(h, i): replace row h with row i XOR row h, tracking phase via g.    *)
(* Per-qubit AG g-function contributions are vectorized via MapThread over     *)
(* the qubit axis (Dimensions[t][[2]] entries) instead of an explicit Table. *)
rowsum[{r_, t_}, h_Integer, i_Integer] :=
    {
        ReplacePart[r,
            h -> 1 - 2 Boole[
                Mod[
                    2 - r[[h]] - r[[i]] + Total @ MapThread[
                        agPhase,
                        {t[[1, All, i]], t[[2, All, i]], t[[1, All, h]], t[[2, All, h]]}
                    ],
                    4
                ] == 2
            ]
        ],
        (tx |-> SubsetMap[BitXor[tx[[All, i]], tx[[All, h]]] &, tx, {All, h}]) /@ t
    }

ps_PauliStabilizer["RowSum", h_Integer, i_Integer] :=
    PauliStabilizer[<|"Signs" -> #1, "Tableau" -> #2|>] & @@ rowsum[{ps["Signs"], ps["Tableau"]}, h, i]


(* ============================================================================ *)
(* Single-qubit Z-basis measurement.                                            *)
(* Returns an Association: deterministic = single key, non-deterministic =      *)
(* {0 -> ps_after_0, 1 -> ps_after_1}.                                          *)
(* Reference: AarGot04 \[Section]3 (measurement update, "destabilizer trick").  *)
(* ============================================================================ *)

(* Out-of-range index: emit a clean PauliStabilizer::partition message and    *)
(* return $Failed (parity with stabilizerEntropy in Entropy.m).               *)
ps_PauliStabilizer["Measure" | "M", a_Integer] /;
    With[{n = ps["GeneratorCount"]}, ! (1 <= a <= n)] := (
        Message[PauliStabilizer::partition, {a}, ps["GeneratorCount"]];
        $Failed
)

(* ============================================================================ *)
(* Packed AG measurement (Stabilizer/Packed.m representation).                   *)
(*                                                                              *)
(* packedRowSumPhase computes the AG i-power (mod 4) of the Pauli product        *)
(* P(row i) . P(row h) by the carry-save trick (external-packages-audit          *)
(* \[Section]4.1): two parity bits accumulated word-by-word over the packed      *)
(* chunks, then popcount via DigitCount. This is the bit-for-bit equivalent of   *)
(* Total @ MapThread[agPhase, ...] over the rank-3 tableau (verified, TIER 11).  *)
(* X, Z are {#chunks, 2n} machine-int matrices (chunk \[Times] row).             *)
(* ============================================================================ *)

packedRowSumPhase[x_, z_, i_Integer, h_Integer] := Mod[#1 + 2 #2, 4] & @@ DigitCount[
    Fold[
        Function[{cc, w},
            With[{
                anti = BitXor[BitAnd[x[[w, h]], z[[w, i]]], BitAnd[x[[w, i]], z[[w, h]]]]
            },
                {
                    BitXor[cc[[1]], anti],
                    BitXor[cc[[2]], BitAnd[BitXor[cc[[1]], x[[w, i]], x[[w, h]], z[[w, i]], z[[w, h]], BitAnd[x[[w, i]], z[[w, h]]]], anti]]
                }
            ]
        ],
        {0, 0},
        Range[Length[x]]
    ],
    2, 1
]

(* packedMeasureZ[ps, a]: Z-basis measurement of qubit a on a concrete ps,       *)
(* returning the same conditional Association as the canonical path -- a single  *)
(* key for deterministic outcomes, {0 -> ., 1 -> .} when random. Post-states are *)
(* packed. Mirrors AarGot04 \[Section]3 (destabilizer trick).                    *)
packedMeasureZ[ps_PauliStabilizer, a_Integer] := Block[{
    tup = psGetPacked[ps], n, x, z, ph, ca, ba, m, pos, p
},
    n = tup[[1]]; x = tup[[2]]; z = tup[[3]]; ph = (1 - tup[[4]]) / 2;
    ca = Quotient[a - 1, $chunkWidth] + 1; ba = Mod[a - 1, $chunkWidth]; m = BitShiftLeft[1, ba];
    pos = Lookup[PositionIndex[Sign[BitAnd[x[[ca]], m]]], 1, {}];
    With[{firstStab = SelectFirst[pos, GreaterThan[n]]},
        If[ MissingQ[firstStab],
            (* deterministic: accumulate the stabilizer product into a scratch    *)
            (* row (column 2n+1) and read its phase bit as the outcome; state      *)
            (* unchanged.                                                          *)
            Block[{xs = Map[Append[#, 0] &, x], zs = Map[Append[#, 0] &, z], phs = Append[ph, 0], sc = 2 n + 1},
                Scan[
                    Function[i,
                        phs[[sc]] = Boole[Mod[2 phs[[sc]] + 2 phs[[i + n]] + packedRowSumPhase[xs, zs, i + n, sc], 4] == 2];
                        xs[[All, sc]] = BitXor[xs[[All, sc]], xs[[All, i + n]]];
                        zs[[All, sc]] = BitXor[zs[[All, sc]], zs[[All, i + n]]]
                    ],
                    Select[pos, LessEqualThan[n]]
                ];
                <|phs[[sc]] -> ps|>
            ],
            (* non-deterministic: clear every other anticommuting row into         *)
            (* firstStab -- destabilizer rows included (AarGot04 \[Section]3:      *)
            (* "for each i != p with x_ia = 1, rowsum(i, p)"); skipping the        *)
            (* destabilizers breaks the symplectic pairing and corrupts later      *)
            (* deterministic outcomes. Then promote firstStab to a destabilizer    *)
            (* and install Z_a with random sign.                                   *)
            Block[{x2 = x, z2 = z, ph2 = ph, p = firstStab},
                Scan[
                    Function[h,
                        (* The product phase i-power is even for stabilizer rows  *)
                        (* (group closure); for a destabilizer row it can be odd, *)
                        (* and the Boole collapses it to phase bit 0, the same    *)
                        (* convention as the canonical rowsum (destabilizer signs *)
                        (* are not physical).                                     *)
                        ph2[[h]] = Boole[Mod[2 ph2[[h]] + 2 ph2[[p]] + packedRowSumPhase[x2, z2, p, h], 4] == 2];
                        x2[[All, h]] = BitXor[x2[[All, h]], x2[[All, p]]];
                        z2[[All, h]] = BitXor[z2[[All, h]], z2[[All, p]]]
                    ],
                    DeleteCases[pos, p]
                ];
                x2[[All, p - n]] = x2[[All, p]]; z2[[All, p - n]] = z2[[All, p]];
                x2[[All, p]] = 0 x2[[All, p]]; z2[[All, p]] = 0 z2[[All, p]];
                z2[[ca, p]] = m;
                Association[# -> psFromPacked[{n, x2, z2, 1 - 2 ReplacePart[ph2, p -> #]}] & /@ {0, 1}]
            ]
        ]
    ]
]


ps_PauliStabilizer["Measure" | "M", a_Integer] := If[psConcreteFastQ[ps] && 1 <= a <= ps["GeneratorCount"],
    packedMeasureZ[ps, a],
    Enclose @ Block[{
        r = ps["Signs"],
        t = ps["Tableau"],
        n = ps["GeneratorCount"],
        pos
    },
        ConfirmAssert[1 <= a <= n];
        pos = Lookup[PositionIndex[t[[1, a]]], 1, {}];
        Association @ With[{p = SelectFirst[pos, GreaterThan[n]]},
            If[ MissingQ[p],
                (* deterministic case: no stabilizer anticommutes with Z_a *)
                {(1 - Last[#1]) / 2 -> PauliStabilizer[<|"Signs" -> Most[#1], "Tableau" -> t|>]} & @@
                    Fold[
                        rowsum[#1, 2 n + 1, #2 + n] &,
                        {Append[r, 1], PadRight[t, Dimensions[t] + {0, 0, 1}]},
                        Select[pos, LessEqualThan[n]]
                    ],

                (* non-deterministic case: one stabilizer anticommutes. Every other *)
                (* row with x_a = 1 -- destabilizers included -- must be rowsummed  *)
                (* with p (AarGot04 \[Section]3), or the symplectic pairing breaks  *)
                (* and later deterministic outcomes come out wrong.                 *)
                Block[{r2, t2},
                    {r2, t2} = Fold[rowsum[#1, #2, p] &, {r, t}, DeleteCases[pos, p]];
                    t2[[All, All, p - n]] = t2[[All, All, p]];
                    t2[[All, All, p]] = 0;
                    t2[[2, a, p]] = 1;
                    # -> PauliStabilizer[<|
                        "Signs" -> ReplacePart[r2, p -> 1 - 2 #],
                        "Tableau" -> t2
                    |>] & /@ {0, 1}
                ]
            ]
        ]
    ]
]


(* Out-of-range index in a partition list: same partition message as Entropy *)
ps_PauliStabilizer["Measure" | "M", qudits : {___Integer}] /;
    With[{n = ps["GeneratorCount"]}, ! AllTrue[qudits, 1 <= # <= n &]] := (
        Message[PauliStabilizer::partition, qudits, ps["GeneratorCount"]];
        $Failed
)

(* Multi-qubit measurement: recursive single-qubit, building a flat keyed assoc *)
ps_PauliStabilizer["Measure" | "M", qudits : {___Integer}] := Enclose @
    If[qudits === {},
        <|{} -> ps|>,
        Join @@ KeyValueMap[
            {k, v} |-> KeyMap[Prepend[k]] @ Confirm @ v["M", Rest[qudits]],
            Confirm @ ps["M", First[qudits]]
        ]
    ]
