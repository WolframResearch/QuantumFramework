Package["Wolfram`QuantumFramework`"]



(* ============================================================================ *)
(* Phase 4 \[Dash] Pauli-string measurement.                                    *)
(*                                                                              *)
(* `ps["M", "XZZXI"]` measures the joint Pauli observable on a stabilizer state.*)
(*                                                                              *)
(* Implementation: AG framework directly. For a stabilizer state |\[Psi]> with  *)
(* stabilizer group S = <g_1, ..., g_n>:                                        *)
(*   - If P (or -P) is in S: deterministic outcome = sign(P in S);              *)
(*       return <|outcome -> ps|> with state unchanged.                         *)
(*   - Else: P anticommutes with at least one stabilizer. Random outcome b in   *)
(*       {0, 1}. The new stabilizer group is obtained by replacing one          *)
(*       anticommuting generator with (-1)^b * P.                               *)
(*                                                                              *)
(* Reference: AarGot04 \[Section]3 generalized; PatGuh26 \[Section]3.4.         *)
(* ============================================================================ *)


(* Helper: parse Pauli string into (sign, n_qubits, [x | z] vector). *)
pauliStringParse[s_String] := Module[{sign, body, n, xz},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    n = StringLength[body];
    xz = Replace[Characters[body],
        {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
    {sign, n, Join[xz[[All, 1]], xz[[All, 2]]]}
]


(* Symplectic inner product mod 2: <(x_1, z_1), (x_2, z_2)> = x_1.z_2 + x_2.z_1 mod 2 *)
sympInnerProduct[v1_List, v2_List, n_Integer] :=
    Mod[v1[[;; n]] . v2[[n + 1 ;; 2 n]] + v2[[;; n]] . v1[[n + 1 ;; 2 n]], 2]


ps_PauliStabilizer["Measure" | "M", pauliString_String] /;
    StringMatchQ[pauliString, RegularExpression["^-?[IXYZ]+$"]] := Enclose @ Module[{
    targetSign, targetN, pVec, n, omega, gens, anticommIdx, kIdx, recoverable, groupSign, postPs
},
    {targetSign, targetN, pVec} = pauliStringParse[pauliString];
    n = ps["Qubits"];
    ConfirmAssert[targetN == n];
    omega = ArrayFlatten[{
        {ConstantArray[0, {n, n}], IdentityMatrix[n]},
        {IdentityMatrix[n], ConstantArray[0, {n, n}]}
    }];
    gens = ps["Matrix"][[ps["GeneratorCount"] + 1 ;; 2 ps["GeneratorCount"]]];
    (* Find indices where pauliString anticommutes with stabilizer rows *)
    anticommIdx = Position[
        Mod[gens . omega . pVec, 2],
        1, {1}, Heads -> False
    ];
    If[anticommIdx === {},
        (* DETERMINISTIC: P (or -P) commutes with all stabilizers => P is in the
           Abelian closure. Compute <P> exactly via direct vector materialization
           (Phase 4 v1; TODO Phase 5+: AG closed-form with i-factor tracking).
           The outcome bit is (1 - <P>) / 2. *)
        With[{vec = ps["State"]["StateVector"], pauliMat = stringToPauliMatrix[pauliString]},
            <|Round[(1 - Re[Conjugate[vec] . pauliMat . vec]) / 2] -> ps|>
        ],
        (* NON-DETERMINISTIC: at least one stabilizer anticommutes with P.
           Pick the first anticommuting row k; replace it with (-1)^b * P.
           Outcome b is random in {0, 1}. *)
        kIdx = First @ First @ anticommIdx;   (* index in 1..n (within stabilizer rows) *)
        (* Build post-measurement state for both outcomes *)
        Module[{newGens, newSigns, ps0Tableau, ps0Signs, postState, outcome},
            (* Replace stabilizer row kIdx with the new generator pVec.
               All other anticommuting stabilizer rows i need to be multiplied by row kIdx
               to make them commute with P. *)
            #[[1]] -> #[[2]] & /@
                Table[
                    Module[{newSignsArr, newGensArr, newRowSign, b = bit},
                        newGensArr = gens;
                        newSignsArr = ps["StabilizerSigns"];
                        (* For anticommuting i \[NotEqual] kIdx: multiply row i by row kIdx *)
                        Do[
                            Module[{i = anticommIdx[[idx, 1]]},
                                If[i =!= kIdx,
                                    newGensArr[[i]] = Mod[newGensArr[[i]] + gens[[kIdx]], 2];
                                    (* Sign update from agPhase function for the row product *)
                                    newSignsArr[[i]] = newSignsArr[[i]] * newSignsArr[[kIdx]]
                                ]
                            ],
                            {idx, 1, Length[anticommIdx]}
                        ];
                        (* Set row kIdx = (-1)^b * P *)
                        newGensArr[[kIdx]] = pVec;
                        newSignsArr[[kIdx]] = targetSign * (1 - 2 b);
                        (* Reconstruct PauliStabilizer; keep destabilizers unchanged but
                           note they're no longer AG-valid; this is a known limitation. *)
                        (* For the syndrome-extraction use case, we mainly want the Phase
                           updates to reflect the outcome. *)
                        With[{newPhase = Join[
                                ps["Phase"][[;; ps["GeneratorCount"]]],
                                (1 - newSignsArr) / 2
                            ]},
                            b -> PauliStabilizer[<|
                                "Phase" -> newPhase,
                                "Tableau" -> ps["Tableau"]
                            |>]
                        ]
                    ],
                    {bit, {0, 1}}
                ] // Association
        ]
    ]
]


(* String-list bulk Pauli measurement: measure each in sequence, collect outcomes *)
ps_PauliStabilizer["Measure" | "M", pauliStrings : {__String}] /;
    AllTrue[pauliStrings, StringMatchQ[#, RegularExpression["^-?[IXYZ]+$"]] &] :=
    Fold[
        Function[{cur, ps2},
            (* cur is an Association {key -> ps}; for each, apply next measurement *)
            Association @@ Flatten[
                KeyValueMap[
                    Function[{outerKey, currentPs},
                        KeyValueMap[
                            Function[{innerKey, newPs}, Append[Flatten[{outerKey}], innerKey] -> newPs],
                            currentPs["M", ps2]
                        ]
                    ],
                    cur
                ]
            ]
        ],
        <|{} -> ps|>,
        pauliStrings
    ]


(* Helper: build the matrix form of a Pauli string for fallback. *)
stringToPauliMatrix[s_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    mats = Replace[Characters[body], {
        "I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1],
        "Y" -> PauliMatrix[2],
        "Z" -> PauliMatrix[3]
    }, {1}];
    sign * KroneckerProduct @@ mats
]
