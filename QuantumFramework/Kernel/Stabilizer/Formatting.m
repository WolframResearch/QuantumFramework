Package["Wolfram`QuantumFramework`"]

PackageScope[PauliForm]
PackageScope[TableauForm]


(* Pauli letter indexed by the code  x + 2 z  (I=00, X=10, Z=01, Y=11), i.e.    *)
(* {0,0}->I, {1,0}->X, {0,1}->Z, {1,1}->Y. The 1-based offset is applied below.  *)
$pauliLetters = {"I", "X", "Z", "Y"}



(* ============================================================================ *)
(* Pauli-string display form                                                    *)
(* ============================================================================ *)

PauliForm[ps_ ? PauliStabilizerQ, n : _Integer ? Positive | Infinity : Infinity, stab_ : True] := If[stab,
    PauliForm[Take[ps["StabilizerSigns"], UpTo[n]], Map[Take[#, UpTo[n]] &, ps["StabilizerTableau"], {2}]],
    PauliForm[Take[ps["DestabilizerSigns"], UpTo[n]], Map[Take[#, UpTo[n]] &, ps["DestabilizerTableau"], {2}]]
]

(* Handle symbolic signs (after SymbolicMeasure, signs may be polynomials in   *)
(* \[FormalS][k] symbols). Numeric signs format as "" / "-"; symbolic signs    *)
(* render as `ToString[s, InputForm] <> "*"` to keep the result readable.      *)
pauliFormSignString[1] := ""
pauliFormSignString[-1] := "-"
pauliFormSignString[s_] := ToString[s, InputForm] <> "*"

(* The Pauli letters are read off by arithmetic indexing  x + 2 z  into          *)
(* $pauliLetters, replacing an O(n_qubits * n_rows) element-by-element pattern    *)
(* Replace with one Part lookup per row. tableau[[1]]/[[2]] are the X/Z blocks    *)
(* {n_qubits, n_rows}; transposing the code grid gives one index vector per row.  *)
PauliForm[signs_, tableau_] := If[MatchQ[Dimensions[tableau], {2, n_, m_} /; 0 < m < n], Append["\[Ellipsis]"], Identity] @ MapThread[
    StringJoin,
    {
        pauliFormSignString /@ signs,
        StringJoin[$pauliLetters[[#]]] & /@ Transpose[tableau[[1]] + 2 tableau[[2]] + 1]
    }
]


(* ============================================================================ *)
(* Tableau grid display form                                                    *)
(* ============================================================================ *)

TableauForm[ps_ ? PauliStabilizerQ, n : _Integer ? Positive | Infinity : Infinity] :=
    TableauForm[
        Take[ps["StabilizerSigns"], UpTo[n]],
        Map[Take[#, UpTo[n]] &, ps["StabilizerTableau"], {2}]
    ]

TableauForm[signs_, tableau_, stab_ : True] := With[{
    short = stab && MatchQ[Dimensions[tableau], {2, n_, m_} /; 0 < m < n],
    r = PadRight[signs, Max[2, Length[signs]], 1],
    t = PadRight[tableau, MapAt[Max[#, If[stab, 1, 2]] &, Dimensions[tableau], -1], " "]
},
    Row[{
        Column[Replace[r, {1 -> "", -1 -> "-"}, {1}]],
        MatrixForm[{{
            Grid[Transpose[Join @@ t], Dividers -> {{Dimensions[t][[2]] + 1 -> True}, If[stab, False, {Length[r] / 2 + 1 -> True}]}]
        }, If[short, {"\[VerticalEllipsis]"}, Nothing]}]
    }]
]


(* ============================================================================ *)
(* MakeBoxes: TraditionalForm + summary box                                     *)
(* ============================================================================ *)

(* Gate TraditionalForm on Qubits <= 8 to avoid OOM.  Above 8 qubits,
   ps["State"] would materialize a 2^n vector -- prohibitive.
   Fallback: typeset the PauliForm string list. *)
MakeBoxes[ps_PauliStabilizer ? PauliStabilizerQ, TraditionalForm] ^:= If[
    ps["Qubits"] <= 8,
    With[{boxes = ToBoxes[ps["State"], TraditionalForm]},
        InterpretationBox[boxes, ps]
    ],
    With[{boxes = ToBoxes[Row[{"\[ScriptCapitalS]", Subscript["", ps["Qubits"]]}, " "], TraditionalForm]},
        InterpretationBox[boxes, ps]
    ]
]

MakeBoxes[ps_PauliStabilizer ? PauliStabilizerQ, form_] ^:= With[{n = ps["GeneratorCount"]},
    BoxForm`ArrangeSummaryBox["PauliStabilizer",
        ps,
        Framed["\[ScriptCapitalS]"],
        If[n < 32,
            {{BoxForm`SummaryItem[{PauliForm[ps, 5]}]},
                If[n <= 5, {BoxForm`SummaryItem[{TableauForm[ps, 5]}]}, Nothing]
            },
            {{BoxForm`SummaryItem[{"Qubits: ", ps["Qudits"]}]},
             {BoxForm`SummaryItem[{"Generators: ", n}]}}
        ],
        If[n < 32,
            {{BoxForm`SummaryItem[{"Destabilizers: ", PauliForm[ps, 5, False]}]},
             {BoxForm`SummaryItem[{"Tableau: ", ps["TableauForm"]}]}},
            {{}}
        ],
        form
    ]
]
