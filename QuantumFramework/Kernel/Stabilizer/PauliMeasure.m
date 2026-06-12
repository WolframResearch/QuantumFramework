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


(* Packed-representation Pauli-string measurement (Stabilizer/Packed.m layout,    *)
(* same physics as the canonical rule below). Row r anticommutes with P iff       *)
(* popcount(x_r AND z_P) + popcount(z_r AND x_P) is odd, computed per chunk by    *)
(* a vectorized DigitCount parity against P's packed word-vector. The AG          *)
(* row-clearing loop reuses packedRowSumPhase (Measurement.m), so each cleared    *)
(* row costs O(n/62) machine words instead of an O(n^2) tableau copy. Returns     *)
(* the same conditional Association contract as the canonical path.               *)
packedMeasurePauliString[ps_PauliStabilizer, pauliString_String] := Enclose @ Block[{
    tup = psGetPacked[ps], targetSign, targetN, pVec, n, x, z, ph, pxP, pzP, par, pos
},
    {targetSign, targetN, pVec} = pauliStringParse[pauliString];
    n = tup[[1]]; x = tup[[2]]; z = tup[[3]]; ph = (1 - tup[[4]]) / 2;
    ConfirmAssert[targetN == n];
    pxP = packChunks[Transpose[{pVec[[;; n]]}], n][[All, 1]];
    pzP = packChunks[Transpose[{pVec[[n + 1 ;;]]}], n][[All, 1]];
    par = Mod[
        Total @ MapThread[
            DigitCount[BitAnd[#1, #4], 2, 1] + DigitCount[BitAnd[#2, #3], 2, 1] &,
            {x, z, pxP, pzP}
        ],
        2
    ];
    pos = Lookup[PositionIndex[par], 1, {}];
    With[{firstStab = SelectFirst[pos, GreaterThan[n]]},
        If[ MissingQ[firstStab],
            (* DETERMINISTIC: P (or -P) commutes with all stabilizers; outcome     *)
            (* from the closed-form i-factor tracking, state unchanged.            *)
            <|(1 - stabilizerExpectation[ps, pauliString]) / 2 -> ps|>,
            (* NON-DETERMINISTIC: clear every other anticommuting row into the     *)
            (* pivot stabilizer p -- destabilizers included (AarGot04 \[Section]3),*)
            (* or the symplectic pairing breaks and later deterministic outcomes   *)
            (* come out wrong. Then promote p to its destabilizer slot and install *)
            (* (-1)^b * P as the new stabilizer.                                   *)
            Block[{x2 = x, z2 = z, ph2 = ph, p = firstStab},
                Scan[
                    Function[h,
                        (* Boole collapses an odd i-power (possible only for a    *)
                        (* destabilizer row, whose sign is not physical) to phase *)
                        (* bit 0, matching the canonical rowsum convention.       *)
                        ph2[[h]] = Boole[Mod[2 ph2[[h]] + 2 ph2[[p]] + packedRowSumPhase[x2, z2, p, h], 4] == 2];
                        x2[[All, h]] = BitXor[x2[[All, h]], x2[[All, p]]];
                        z2[[All, h]] = BitXor[z2[[All, h]], z2[[All, p]]]
                    ],
                    DeleteCases[pos, p]
                ];
                x2[[All, p - n]] = x2[[All, p]]; z2[[All, p - n]] = z2[[All, p]];
                x2[[All, p]] = pxP; z2[[All, p]] = pzP;
                Association[
                    # -> psFromPacked[{n, x2, z2, 1 - 2 ReplacePart[ph2, p -> (1 - targetSign (1 - 2 #)) / 2]}] & /@ {0, 1}
                ]
            ]
        ]
    ]
]


(* Concrete states take the packed path; the canonical rank-3-array path handles  *)
(* symbolic signs (non-Clifford P/T residues, SymPhase states), which fail the    *)
(* psConcreteFastQ gate.                                                           *)
ps_PauliStabilizer["Measure" | "M", pauliString_String] /;
    StringMatchQ[pauliString, RegularExpression["^-?[IXYZ]+$"]] := If[psConcreteFastQ[ps],
    packedMeasurePauliString[ps, pauliString],
    Enclose @ Module[{
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
]]


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
