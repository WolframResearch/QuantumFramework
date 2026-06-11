Package["Wolfram`QuantumFramework`"]

PackageScope[FromFullTableau]



(* ============================================================================ *)
(* Constructors: from QuantumState (4^n tomography).                            *)
(*                                                                              *)
(* Compute <psi|P|psi> for every Pauli P (4^n tomography), keep the 2^n Paulis *)
(* with |<P>| ~= 1 as the stabilizer group with their signs, pick n linearly-  *)
(* independent generators by greedy F_2 rank, then extend to a 2n-dim          *)
(* symplectic basis to obtain destabilizers that satisfy the AG commutation    *)
(* pattern. Cost is O(4^n) -- practical for n <= ~8.                            *)
(* ============================================================================ *)

(* Symplectic form J over F_2 used by the destabilizer search. *)
agSymplecticForm[n_Integer] := ArrayFlatten[{
    {ConstantArray[0, {n, n}], IdentityMatrix[n]},
    {IdentityMatrix[n], ConstantArray[0, {n, n}]}
}]

(* Encode a tuple of Paulis (each in 0..3) as a 2n-bit row {x_1,..,x_n,z_1,..,z_n}. *)
(* Convention matches PauliRow / FromFullTableau: I=00, X=10, Y=11, Z=01.       *)
agEncodePauli[plist_List] := Catenate[Transpose @ Replace[plist, {0 -> {0, 0}, 1 -> {1, 0}, 2 -> {1, 1}, 3 -> {0, 1}}, {1}]]

(* Find n linearly-independent rows of `bits` over F_2 by greedy rank-growth.  *)
(* Returns indices into `bits`. Catch + Internal`Bag avoids the AppendTo/Break  *)
(* anti-pattern; the algorithm itself is genuinely sequential (each accepted   *)
(* index changes the rank residual seen by later candidates).                   *)
agPickIndependent[bits_, n_Integer] := Module[{bag = Internal`Bag[]},
    Catch @ Do[
        With[{kept = Internal`BagPart[bag, All]},
            If[MatrixRank[Append[bits[[kept]], bits[[i]]], Modulus -> 2] > Length[kept],
                Internal`StuffBag[bag, i];
                If[Length[kept] + 1 >= n, Throw[Null]]
            ]
        ],
        {i, Length[bits]}
    ];
    Internal`BagPart[bag, All]
]

(* Given n stabilizer encodings (n x 2n), construct n destabilizer encodings   *)
(* by symplectic Gram-Schmidt: each D_i satisfies symp(D_i, S_j) = delta_ij    *)
(* and symp(D_i, D_k) = 0 for k < i. Each D_i depends on prior D_k via the    *)
(* augmented constraint matrix; expressed as Fold over Range[n].               *)
agExtendToSymplecticBasis[stabBits_, n_Integer] := With[{J = agSymplecticForm[n]},
    Fold[
        Function[{dests, i},
            Append[dests,
                LinearSolve[
                    Mod[Join[stabBits, dests] . J, 2],
                    Join[
                        Normal @ SparseArray[{i -> 1}, n],
                        ConstantArray[0, Length[dests]]
                    ],
                    Modulus -> 2
                ]
            ]
        ],
        {},
        Range[n]
    ]
]

(* Compute <psi|P|psi> for every Pauli P on n qubits. Returns a list of        *)
(* {encoding, sign} pairs for the 2^n Paulis with |<P>| ~= 1. The first  *)
(* element is always {identity, 1}. *)
agStabilizerGroup[v_, n_Integer] := Enclose @ Module[{paulis, expVals, group, identity, nonid},
    paulis = Tuples[Range[0, 3], n];
    expVals = Map[Chop[Conjugate[v] . (kroneckerProduct @@ Map[pauliMatrix, #]) . v] &, paulis];
    group = Cases[
        Transpose[{paulis, expVals}],
        {p_, val_} /; Re[val] >= 1 - 10^-8 :> {agEncodePauli[p],  1},
        {1}
    ] ~Join~ Cases[
        Transpose[{paulis, expVals}],
        {p_, val_} /; Re[val] <= -(1 - 10^-8) :> {agEncodePauli[p], -1},
        {1}
    ];
    ConfirmAssert[Length[group] == 2^n,
        "PauliStabilizer[QuantumState] expects a stabilizer state; got " <>
        ToString[Length[group]] <> " group elements (need 2^" <> ToString[n] <> ")"
    ];
    group
]

PauliStabilizer[qs_QuantumState] := Enclose @ Module[
    {n, v, encodings, signs, nonidIdx, genIdx, stabBits, stabSigns,
     destBits, destSigns, fullTableau, basePS, baseAssoc, baseVec, anchor, phase},
    n = qs["Qudits"];
    v = Normal @ qs["Computational"]["StateVector"];

    {encodings, signs} = Transpose @ ConfirmMatch[
        agStabilizerGroup[v, n],
        {{_, _} ..}
    ];

    (* Drop the identity row (all-zero encoding) before picking generators. *)
    nonidIdx = Position[encodings, _ ? (! AllTrue[#, # == 0 &] &), {1}, Heads -> False][[All, 1]];
    genIdx = nonidIdx[[ agPickIndependent[encodings[[nonidIdx]], n] ]];
    stabBits = encodings[[genIdx]];
    stabSigns = signs[[genIdx]];

    (* Canonicalize generator order to match the integer-constructor          *)
    (* convention PauliStabilizer[n] = {Z_1, Z_2, ..., Z_n}.                  *)
    (* Sort by FirstPosition[row, 1] ascending (row whose active bit is at   *)
    (* the lowest-index qubit comes first).                                   *)
    With[{order = Ordering[
        FirstPosition[#, 1, {2 n + 1}] & /@ stabBits
    ]},
        stabBits = stabBits[[order]];
        stabSigns = stabSigns[[order]]
    ];

    destBits = agExtendToSymplecticBasis[stabBits, n];
    destSigns = ConstantArray[1, n];

    (* Build the (2n) x (2n+1) tableau in PauliRow/FromFullTableau format:     *)
    (* rows 1..n  = destabilizers, rows n+1..2n = stabilizers.                 *)
    (* Last column = phase (0 = sign +1, 1 = sign -1).                         *)
    fullTableau = Join[
        MapThread[Append, {destBits, (1 - destSigns) / 2}],
        MapThread[Append, {stabBits, (1 - stabSigns) / 2}]
    ];

    basePS = Confirm @ FromFullTableau[fullTableau];
    baseAssoc = First[basePS];

    (* Recover the overall phase the tableau drops: ps["State"] is exact only    *)
    (* up to a complex unit. Anchor at the first nonzero amplitude of the         *)
    (* tableau-projected state and store the ratio to the original.               *)
    baseVec = Normal @ basePS["State"]["StateVector"];
    anchor = First @ FirstPosition[baseVec, x_ /; Chop[x] =!= 0, {1}, Heads -> False];
    phase = v[[anchor]] / baseVec[[anchor]];

    PauliStabilizer[<|baseAssoc, "GlobalPhase" -> phase|>]
]


(* ============================================================================ *)
(* Association-keyed constructors                                               *)
(* ============================================================================ *)

PauliStabilizer[data : KeyValuePattern[{"Phase" -> phase_}]] := PauliStabilizer[<|KeyDrop[data, "Phase"], "Signs" -> 1 - 2 phase|>]

PauliStabilizer[data : KeyValuePattern[{"Matrix" -> mat_}]] := PauliStabilizer[<|KeyDrop[data, "Matrix"], "Tableau" -> Transpose[ArrayReshape[mat, {Length[mat], 2, Length[mat] / 2}], {3, 1, 2}]|>]

PauliStabilizer[data : KeyValuePattern[{"Tableau" -> t_}]] /; ! KeyExistsQ[data, "Signs"] :=
    PauliStabilizer[<|"Signs" -> ConstantArray[1, Dimensions[t][[3]]], "Tableau" -> t|>]


(* ============================================================================ *)
(* From QuantumOperator                                                         *)
(* ============================================================================ *)

FromFullTableau[t_] := PauliStabilizer[<|
    "Signs" -> 1 - 2 t[[All, -1]],
    "Tableau" -> Transpose[ArrayReshape[t[[All, ;; -2]], {Length[t], 2, Length[t] / 2}], {3, 1, 2}]
|>]


(* Capture the global phase the AG decomposition drops. The tableau is         *)
(* invariant under U -> e^{i alpha} U, so for Y = -i X Z the AG-recovered     *)
(* circuit Z@X equals i*Y; we record -i as "GlobalPhase" so that the inverse  *)
(* path ps["QuantumOperator"] returns Y exactly.                                *)
PauliStabilizer[qo_QuantumOperator, n : _Integer : 1] /; Sort[qo["OutputOrder"]] == Sort[qo["InputOrder"]] :=
Enclose @ Module[{max, basePS, baseAssoc, m1, m2, anchor, phase},
    max = Max[n, qo["InputOrder"]];
    basePS = FromFullTableau[
        Confirm @ QuantumOperatorTableau[qo["Sort"]["Computational"]]
    ]["PadRight", max]["Permute", Join[Sort[qo["InputOrder"]], Complement[Range[max], qo["InputOrder"]]]];

    baseAssoc = First[basePS];
    m1 = Flatten[Normal @ qo["Matrix"]];
    m2 = Flatten[Normal @ basePS["QuantumOperator"]["Matrix"]];
    anchor = First @ FirstPosition[m2, x_ /; Chop[x] =!= 0, {1}, Heads -> False];
    phase = m1[[anchor]] / m2[[anchor]];

    PauliStabilizer[<|baseAssoc, "GlobalPhase" -> phase|>]
]


(* ============================================================================ *)
(* From QuantumCircuitOperator: fold over normal operators                      *)
(* ============================================================================ *)

PauliStabilizer[qco_QuantumCircuitOperator] := Fold[
    Collect[#1 /. ps_PauliStabilizer :> #2[ps], _PauliStabilizer, Simplify] &,
    PauliStabilizer[qco["Max"]],
    Replace[PauliStabilizer[#, qco["Max"]], _ ? FailureQ :> #] & /@ qco["NormalOperators"]
]


(* ============================================================================ *)
(* From rank-3 +/-1/0 array                                                     *)
(* ============================================================================ *)

PauliStabilizer[basis_ /; ArrayQ[basis, 3, MatchQ[0 | 1 | -1]] && MatchQ[Dimensions[basis], {2, n_, m_}]] :=
    Enclose @ PauliStabilizer[
        Confirm @ Replace[Union @@ #, {{0, 1} | {1} -> 1, {-1, 0} | {-1} -> -1, _ -> $Failed}] & /@ Transpose[basis, {3, 2, 1}],
        Abs[basis]
    ]


(* ============================================================================ *)
(* Single-row / paired-sign constructors (auto-pad destabilizer half)           *)
(* ============================================================================ *)

PauliStabilizer[tableau_ ? PauliTableauQ] := PauliStabilizer[1, tableau]

PauliStabilizer[sign : -1 | 1, tableau_ ? PauliTableauQ] := PauliStabilizer[{sign}, tableau]

PauliStabilizer[signs : {(-1 | 1) ...}, tableau_ ? PauliTableauQ] := With[{padSigns = PadRight[signs, Dimensions[tableau][[3]], 1]},
    PauliStabilizer[
        <|"Signs" -> Join[padSigns, padSigns], "Tableau" -> MapThread[Join[##, 2] &, {Reverse[tableau], tableau}]|>
    ]
]


(* ============================================================================ *)
(* String-list constructors (Pauli strings with optional sign prefix)           *)
(* ============================================================================ *)

$PauliString = Repeated["-" | "+", {0, 1}] ~~ ("I" | "X" | "Y" | "Z") ...

PauliStabilizer[signedPauliStrings : {__String}] /; AllTrue[signedPauliStrings, StringMatchQ[$PauliString]] :=
    PauliStabilizer @@
        MapAt[Transpose[PadRight[#, Automatic, {0, 0}], {3, 2, 1}] &, 2] @ Thread @ Replace[
            Characters[signedPauliStrings], {sign : "-" | "+" : 1, paulis___} :> {
                Replace[sign, {"-" -> -1, "+" -> 1}],
                Replace[{paulis}, {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}]
            },
            {1}
        ]

PauliStabilizer[stabString : {__String}, destabStrings : {__String}] /; AllTrue[Join[stabString, destabStrings], StringMatchQ[$PauliString]] :=
    PauliStabilizer[<|"Signs" -> #1, "Tableau" -> #2|>] & @@
        MapAt[Transpose[PadRight[#, Automatic, {0, 0}], {3, 2, 1}] &, 2] @ Thread @ Replace[
            Characters[Join[destabStrings, stabString]], {sign : "-" | "+" : 1, paulis___} :> {
                Replace[sign, {"-" -> -1, "+" -> 1}],
                Replace[{paulis}, {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}]
            },
            {1}
        ]


(* ============================================================================ *)
(* Integer constructor: |0...0> register                                        *)
(*                                                                              *)
(* The register tableau is known in closed form (destabilizers X_i in rows      *)
(* 1..q, stabilizers Z_i in rows q+1..2q, all signs +1), so it is assembled     *)
(* directly instead of routing the {2, q, 2q} array through the generic rank-3  *)
(* validating constructor, whose per-element pattern matching is O(q^2) and     *)
(* dominated circuit-simulation setup at q ~ 10^3.                              *)
(* ============================================================================ *)

PauliStabilizer[q_Integer ? Positive] := PauliStabilizer[<|
    "Signs" -> ConstantArray[1, 2 q],
    "Tableau" -> {
        PadRight[IdentityMatrix[q], {q, 2 q}],
        Join[ConstantArray[0, {q, q}], IdentityMatrix[q], 2]
    }
|>]

PauliStabilizer[q_Integer ? NonNegative] := PauliStabilizer[{ConstantArray[0, {Max[q, 1], q}], identityMatrix[q]}]


(* ============================================================================ *)
(* Empty + shortcut                                                             *)
(* ============================================================================ *)

PauliStabilizer[] := PauliStabilizer[1]

PauliStabilizer[shortcut : _String | (_String -> _ ? orderQ | _Integer) | _List] := PauliStabilizer[QuantumCircuitOperator[shortcut]]
