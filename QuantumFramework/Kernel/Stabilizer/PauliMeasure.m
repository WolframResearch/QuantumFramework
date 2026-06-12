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
    targetSign, targetN, pVec, n, omega, anticommPos
},
    {targetSign, targetN, pVec} = pauliStringParse[pauliString];
    n = ps["Qubits"];
    ConfirmAssert[targetN == n];
    omega = ArrayFlatten[{
        {ConstantArray[0, {n, n}], IdentityMatrix[n]},
        {IdentityMatrix[n], ConstantArray[0, {n, n}]}
    }];
    (* Rows of ps["Matrix"] (destabilizers 1..n, stabilizers n+1..2n) that      *)
    (* anticommute with P, via the symplectic form.                             *)
    anticommPos = Lookup[PositionIndex[Mod[ps["Matrix"] . omega . pVec, 2]], 1, {}];
    With[{firstStab = SelectFirst[anticommPos, GreaterThan[n]]},
    If[MissingQ[firstStab],
        (* DETERMINISTIC: P (or -P) commutes with all stabilizers => P is in the
           Abelian closure. The expectation <P> = +-1 is recovered by the AG
           closed-form i-factor tracking (stabilizerExpectation, O(n^2)) rather
           than materializing the 2^n state vector. The outcome bit is (1-<P>)/2. *)
        <|(1 - stabilizerExpectation[ps, pauliString]) / 2 -> ps|>,
        (* NON-DETERMINISTIC: at least one stabilizer anticommutes with P.       *)
        (* AG §3 generalized to a Pauli string: rowsum EVERY other anticommuting *)
        (* row -- destabilizers included -- into the pivot stabilizer p, so all  *)
        (* remaining generators commute with P and the symplectic pairing is     *)
        (* preserved (skipping rows here silently corrupts later deterministic   *)
        (* outcomes). Then promote p's row to its destabilizer slot and install  *)
        (* (-1)^b * P as the new stabilizer. Outcome b is random in {0, 1}.      *)
        Block[{r2, t2},
            {r2, t2} = Fold[rowsum[#1, #2, firstStab] &, {ps["Signs"], ps["Tableau"]}, DeleteCases[anticommPos, firstStab]];
            t2[[All, All, firstStab - n]] = t2[[All, All, firstStab]];
            t2[[1, All, firstStab]] = pVec[[;; n]];
            t2[[2, All, firstStab]] = pVec[[n + 1 ;;]];
            Association[
                # -> PauliStabilizer[<|
                    "Signs" -> ReplacePart[r2, firstStab -> targetSign * (1 - 2 #)],
                    "Tableau" -> t2
                |>] & /@ {0, 1}
            ]
        ]
    ]]
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
