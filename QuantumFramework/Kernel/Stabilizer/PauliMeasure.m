Package["Wolfram`QuantumFramework`"]



(* ============================================================================ *)
(* Pauli-string measurement.                                                    *)
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
    targetSign, targetN, pVec, n, omega, gens, anticommIdx, kIdx
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
           Abelian closure. The expectation <P> = +-1 is recovered by the AG
           closed-form i-factor tracking (stabilizerExpectation, O(n^2)) rather
           than materializing the 2^n state vector. The outcome bit is (1-<P>)/2. *)
        <|(1 - stabilizerExpectation[ps, pauliString]) / 2 -> ps|>,
        (* NON-DETERMINISTIC: at least one stabilizer anticommutes with P.       *)
        (* AG §3 Case I: replace the first anticommuting generator with (-1)^b·P*)
        (* (step 1) and update the OTHER anticommuting generators by row-       *)
        (* multiplication so they commute with P (step 2). Outcome b is random  *)
        (* in {0, 1}.                                                            *)
        kIdx = First @ First @ anticommIdx;   (* index in 1..n (within stabilizer rows) *)
        (* Update both the sign AND the tableau row: stabilizer row kIdx is    *)
        (* overwritten with pVec's symplectic bits so the post-state is a       *)
        (* +-P eigenstate (AG §3 Case I, step 1).                              *)
        Module[{otherIdx, baseSigns, newTableau, gen = ps["GeneratorCount"]},
            otherIdx = DeleteCases[Flatten[anticommIdx, 1], kIdx];
            baseSigns = MapAt[# * ps["StabilizerSigns"][[kIdx]] &, ps["StabilizerSigns"], List /@ otherIdx];
            (* Replace stabilizer row kIdx of the tableau with pVec's bits.     *)
            (* Tableau shape is {2 (X/Z block), n_qubits, 2*GeneratorCount}.    *)
            (* Stabilizer rows live at indices gen+1 .. 2 gen along axis 3.    *)
            newTableau = ps["Tableau"];
            newTableau[[1, All, gen + kIdx]] = pVec[[;; n]];
            newTableau[[2, All, gen + kIdx]] = pVec[[n + 1 ;;]];
            Association @ Table[
                With[{
                    newSigns = ReplacePart[baseSigns, kIdx -> targetSign * (1 - 2 b)]
                },
                    b -> PauliStabilizer[<|
                        "Phase" -> Join[
                            ps["Phase"][[;; gen]],
                            (1 - newSigns) / 2
                        ],
                        "Tableau" -> newTableau
                    |>]
                ],
                {b, {0, 1}}
            ]
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
