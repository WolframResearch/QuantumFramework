Package["Wolfram`QuantumFramework`"]

PackageScope[agPhase]
PackageScope[rowsum]



(* ============================================================================ *)
(* AG phase function and rowsum primitive (line :280, :282 in old kernel)       *)
(* Reference: AarGot04 \[Section]3, Lemma 3.                                     *)
(* ============================================================================ *)

(* The AG phase increment g(x1, z1, x2, z2) for multiplying tableau rows.       *)
(* (Symbol renamed from `g` to `agPhase` to avoid PackageScope name collision.) *)
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

ps_PauliStabilizer["Measure" | "M", a_Integer] := Enclose @ Block[{
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

            (* non-deterministic case: one stabilizer anticommutes *)
            Block[{r2, t2},
                {r2, t2} = Fold[rowsum[#1, #2, p] &, {r, t}, Select[pos, GreaterThan[p]]];
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


(* Multi-qubit measurement: recursive single-qubit, building a flat keyed assoc *)
ps_PauliStabilizer["Measure" | "M", qudits : {___Integer}] := Enclose @
    If[qudits === {},
        <|{} -> ps|>,
        Join @@ KeyValueMap[
            {k, v} |-> KeyMap[Prepend[k]] @ Confirm @ v["M", Rest[qudits]],
            Confirm @ ps["M", First[qudits]]
        ]
    ]
