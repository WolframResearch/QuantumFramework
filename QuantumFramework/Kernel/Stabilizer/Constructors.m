Package["Wolfram`QuantumFramework`"]

PackageScope[FromFullTableau]



(* ============================================================================ *)
(* Local helpers used by string-list and array constructors.                    *)
(* (firstIndex / bitvector also used by PauliRow in Conversions.m)              *)
(* ============================================================================ *)

(* From-state constructor: tomography over all 4^n Pauli expectations.         *)
(* Cost: O(4^n) -- practical only for n <= 8 or so.                            *)
PauliStabilizerTableau[v_ ? ArrayQ, n_Integer ? Positive] :=
    Enclose @ Replace[#, Thread[PadRight[ConfirmBy[DeleteCases[Union[Flatten[#]], 0], Length[#] <= 2 &], 2, -1] -> {1, -1}], {3}] & @ Transpose[
        Partition[#, 2] & /@ Take[#, n] & @ ReverseSort @ ResourceFunction["RowSpace"][
            With[{pauli = Tuples[Range[0, 3], n]},
                Map[
                    (Conjugate[v] . (kroneckerProduct @@ PauliMatrix /@ #) . v) *
                        Catenate[Replace[#, {0 -> {0, 0}, 1 -> {1, 0}, 2 -> {1, 1}, 3 -> {0, 1}}, {1}]] &,
                    pauli
                ]
            ],
            "Basis"
        ],
        {3, 2, 1}
    ]



(* ============================================================================ *)
(* Constructors: from QuantumState (expensive 4^n tomography)                   *)
(* ============================================================================ *)

PauliStabilizer[qs_QuantumState] := Enclose @ With[{
    z = ConfirmBy[PauliStabilizerTableau[qs["Computational"]["StateVector"], qs["Qudits"]], PauliTableauQ],
    x = ConfirmBy[PauliStabilizerTableau[QuantumState[qs, Table["X", qs["Qudits"]]]["StateVector"], qs["Qudits"]], PauliTableauQ]
},
{
    tableau = MapThread[Join[##, 2] &, {x, z}]
},
    PauliStabilizer[<|"Signs" -> Map[First[ConfirmBy[DeleteCases[#, 0], Apply[Equal]], 1] &, Transpose[Catenate[tableau]]], "Tableau" -> Abs @ tableau|>]
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


PauliStabilizer[qo_QuantumOperator, n : _Integer : 1] /; Sort[qo["OutputOrder"]] == Sort[qo["InputOrder"]] := With[{max = Max[n, qo["InputOrder"]]},
    Enclose @ FromFullTableau[
        Confirm @ QuantumOperatorTableau[qo["Sort"]["Computational"]]
    ]["PadRight", max]["Permute", Join[Sort[qo["InputOrder"]], Complement[Range[max], qo["InputOrder"]]]]
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
(* ============================================================================ *)

PauliStabilizer[q_Integer ? NonNegative] := PauliStabilizer[{ConstantArray[0, {Max[q, 1], q}], identityMatrix[q]}]


(* ============================================================================ *)
(* Empty + shortcut                                                             *)
(* ============================================================================ *)

PauliStabilizer[] := PauliStabilizer[1]

PauliStabilizer[shortcut : _String | (_String -> _ ? orderQ | _Integer) | _List] := PauliStabilizer[QuantumCircuitOperator[shortcut]]
