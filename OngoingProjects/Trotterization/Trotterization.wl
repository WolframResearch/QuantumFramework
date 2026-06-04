(* ::Package:: *)

(* ::Title:: *)
(* Standalone Suzuki-Trotter decomposition *)

(* ::Text:: *)
(* Pure Wolfram Language implementation of the Suzuki-Trotter product
   formula, independent from QuantumFramework.

   Given a Hamiltonian H = Sum_i A_i (each A_i a Hermitian d x d matrix),
   the time evolution operator U(t) = Exp[-I t H] is approximated by a
   product of single-term exponentials Exp[-I c_k A_{i(k)}], where the
   coefficients c_k and term ordering follow Suzuki's recursion (see
   arXiv:math-ph/0506007).

   Mirrors QF's Trotterization in NamedCircuits.m:604-629. The recursion
   for trotterCoeffs and trotterExpand is identical so that the matrix
   list returned here, multiplied right-to-left, equals QF's circuit
   matrix exactly. *)

BeginPackage["TrotterizationStandalone`"];

TrotterizationMatrices::usage =
"TrotterizationMatrices[ops, order, reps, time] returns the list of d x d \
gate matrices for the Suzuki-Trotter decomposition of Exp[-I time Total[ops]] \
at the given order with the given number of Trotter steps. The gates are \
returned in circuit-application order; the approximate unitary is \
Dot @@ Reverse[gates]. The argument ops is a list of d x d Hermitian \
matrices A_i with H = Sum_i A_i.

Options:
  \"Form\" -> \"Flat\"    -- flat list of gates (default)
  \"Form\" -> \"Steps\"   -- list of length reps; each entry the gates of one Trotter step
  \"Form\" -> \"Unitary\" -- the assembled approximate unitary as a single matrix";

TrotterizationCoefficients::usage =
"TrotterizationCoefficients[l, order, c] returns the Suzuki-Trotter coefficient \
list for l Hamiltonian terms at the given order with overall scale c. \
Equivalent to QF's internal trotterCoeffs.";

TrotterizationOrdering::usage =
"TrotterizationOrdering[ops, order] returns the term-index sequence (or the \
expanded operator list, if ops is a list of operators) used by the Suzuki-Trotter \
decomposition at the given order. Equivalent to QF's internal trotterExpand.";

Begin["`Private`"];

(* --- Suzuki recursion: coefficients ---
   Order 1: l copies of c
   Order 2: c/2 forward then c/2 reverse
   Even n >= 4: Suzuki's fractal recursion with p = 1 / (4 - 4^(1/(n-1)))
   Odd n: rounded up to nearest even *)

TrotterizationCoefficients[l_Integer ? Positive, 1, c_ : 1] :=
    ConstantArray[c, l];

TrotterizationCoefficients[l_Integer ? Positive, 2, c_ : 1] :=
    With[{s = TrotterizationCoefficients[l, 1, c / 2]},
        Join[s, Reverse[s]]
    ];

TrotterizationCoefficients[l_Integer ? Positive, n_Integer ? EvenQ, c_ : 1] /; n >= 4 :=
    With[{p = 1 / (4 - 4 ^ (1 / (n - 1)))},
        With[{
            s = TrotterizationCoefficients[l, n - 2, c p],
            sm = TrotterizationCoefficients[l, n - 2, (1 - 4 p) c]
        },
            Join[s, s, sm, s, s]
        ]
    ];

TrotterizationCoefficients[l_Integer ? Positive, n_Integer ? OddQ, c_ : 1] :=
    TrotterizationCoefficients[l, Round[n, 2], c];


(* --- Suzuki recursion: operator-list expansion --- *)

TrotterizationOrdering[l_List, 1] := l;

TrotterizationOrdering[l_List, 2] := Join[l, Reverse[l]];

TrotterizationOrdering[l_List, n_Integer ? EvenQ] /; n >= 4 :=
    Catenate @ Table[TrotterizationOrdering[l, n - 2], 5];

TrotterizationOrdering[l_List, n_Integer ? OddQ] :=
    TrotterizationOrdering[l, Round[n, 2]];


(* --- Public API --- *)

Options[TrotterizationMatrices] = {"Form" -> "Flat"};

TrotterizationMatrices[
    ops : {__ ? MatrixQ},
    order : (_Integer ? Positive) : 1,
    reps : (_Integer ? Positive) : 1,
    time_ : 1,
    OptionsPattern[]
] := Module[{coeffs, expandedOps, gates, perStep, form},
    form = OptionValue["Form"];

    (* Coefficients carry the 1/reps factor so that the gates of one
       Trotter step approximate Exp[-I (time/reps) H]. *)
    coeffs = time * TrotterizationCoefficients[Length[ops], order, 1 / reps];
    expandedOps = TrotterizationOrdering[ops, order];

    (* Gates of a single Trotter step, in circuit (application) order. *)
    gates = MapThread[MatrixExp[-I #1 #2] &, {coeffs, expandedOps}];

    (* All reps are identical; total gate count is reps * Length[expandedOps]. *)
    perStep = ConstantArray[gates, reps];

    Switch[form,
        "Flat",    Catenate[perStep],
        "Steps",   perStep,
        "Unitary", Dot @@ Reverse[Catenate[perStep]],
        _,         Message[TrotterizationMatrices::badform, form];
                   Catenate[perStep]
    ]
];

TrotterizationMatrices::badform =
"Unknown \"Form\" value `1`. Falling back to \"Flat\".";

End[];
EndPackage[];
