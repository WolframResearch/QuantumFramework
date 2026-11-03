(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[$QECStabilizerGroupLimit]

PackageScope[codeStabilizerGroup]
PackageScope[codeWeightEnumerator]
PackageScope[codeSignedWeightEnumerator]
PackageScope[codeDualWeightEnumerator]
PackageScope[codeCodespaceProbability]


(* ============================================================================ *)
(* Weight enumerators (Got26 sec. 7.4.1) and what they are good for.            *)
(*                                                                              *)
(* A_j is the number of stabilizer elements of weight j, B_j the number of      *)
(* normalizer elements of weight j (Props. 7.5, 7.6).  Two facts make them      *)
(* worth having as properties rather than as an exercise:                       *)
(*                                                                              *)
(*   - Theorem 7.7(c): the code has distance d iff A_j = B_j for j < d.  That    *)
(*     is a second, independent route to the distance, and the tests use it    *)
(*     to check the distance search against a count.                           *)
(*   - The quantum MacWilliams identity (Thm 7.8),                              *)
(*         B(x) = (K / 2^n) (1 + 3x)^n A((1 - x) / (1 + 3x)),                   *)
(*     so B never needs the normalizer enumerated: it follows from A, which     *)
(*     needs only the 2^(n-k) elements of the stabilizer.  The normalizer has   *)
(*     2^(n+k), which is the difference between 2^8 and 2^22 on the [[15,7,3]]  *)
(*     Hamming code.  The tests enumerate N(S) directly on small codes and      *)
(*     check the identity against it.                                           *)
(*                                                                              *)
(* THE SIGNED ENUMERATOR.  A(x) counts Paulis and forgets their signs, which is *)
(* right for the distance and wrong for the question sec. 13.4 asks of it: the  *)
(* probability that n copies of a one-qubit state land in the code space,       *)
(*                                                                              *)
(*     tr[rho^(x)n Pi] = 2^-(n-k) Sum_{M in S} phase(M) Prod_q tr(rho M_q),     *)
(*                                                                              *)
(* eq. 13.36.  The book can drop phase(M) because every element of the 5-qubit  *)
(* code's stabilizer has sign +1; a code built from a check matrix with signs   *)
(* cannot.  So the layer keeps both: "WeightEnumerator" (Def 7.1, unsigned) and *)
(* "SignedWeightEnumerator" (Sum phase(M) x^wt M), and "CodespaceProbability"   *)
(* computes the trace above for any one-qubit state, symbolic or not, without   *)
(* forming a 2^n matrix.  On five copies of the magic state |R> it gives the    *)
(* 1/6 of protocol 13.5, and the tests check it against the dense projector.    *)
(*                                                                              *)
(* The stabilizer group is enumerated once per code and memoised; its size is   *)
(* bounded by $QECStabilizerGroupLimit, with a message instead of a hang.       *)
(*                                                                              *)
(* References: Got26 sec. 7.4.1 (Defs 7.1-7.3, Props 7.5-7.6, Thm 7.7),         *)
(* sec. 7.4.2 (Thm 7.8), sec. 13.4 (protocol 13.5, eqs. 13.33-13.37).           *)
(* ============================================================================ *)

$QECStabilizerGroupLimit::usage = "$QECStabilizerGroupLimit is the largest stabilizer group, in elements, that a code enumerates to compute its weight enumerators and code-space probabilities.";

$QECStabilizerGroupLimit = 2^20;

QECCode::groupbig = "The stabilizer group has 2^`1` elements, above $QECStabilizerGroupLimit = `2`. Raise the limit to enumerate it.";
QECCode::blochstate = "`1` is not a one-qubit state. Give a one-qubit QuantumState or a Bloch vector {x, y, z}.";


(* ---- the group ---- *)

(* Every element of S as a row with its true phase: each generator doubles the list.
   2^(n-k) products, no search. *)
codeStabilizerGroup[a_Association] := codeStabilizerGroup[a] = Fold[
    Function[{elements, g}, Join[elements, pauliProduct[#, g] & /@ elements]],
    {pauliIdentity[a["Qubits"]]},
    codeVectors[a]
]

stabilizerGroupAllowedQ[a_Association] := With[{m = codeStabilizerCount[a]},
    If[ 2^m <= $QECStabilizerGroupLimit,
        True,
        Message[QECCode::groupbig, m, $QECStabilizerGroupLimit]; False
    ]
]

rowWeights[rows_List] := With[{n = pauliQubits[First[rows]]},
    Total[BitOr[rows[[All, 1 ;; n]], rows[[All, n + 1 ;; 2 n]]], {2}]
]

(* Hermitian elements only have phase 0 or 2; i^e is then +1 or -1. *)
rowSigns[rows_List] := I^rows[[All, -1]]


(* ---- the three enumerators, as coefficient lists {c_0, ..., c_n} ---- *)

codeWeightEnumerator[a_Association] := With[{w = rowWeights[codeStabilizerGroup[a]]},
    Table[Count[w, j], {j, 0, a["Qubits"]}]
]

codeSignedWeightEnumerator[a_Association] := With[
    {rows = codeStabilizerGroup[a]},
    With[{w = rowWeights[rows], s = rowSigns[rows]},
        Table[Total[Pick[s, w, j]], {j, 0, a["Qubits"]}]
    ]
]

(* The MacWilliams identity, applied to A.  Exact integers come out because the
   identity is exact; Expand is enough, no simplification is guessed at. *)
codeDualWeightEnumerator[a_Association] := Module[
    {n = a["Qubits"], k = codeLogicalQubits[a], x, A, B},
    A = codeWeightEnumerator[a] . x^Range[0, n];
    B = Expand[Together[2^k / 2^n (1 + 3 x)^n (A /. x -> (1 - x) / (1 + 3 x))]];
    CoefficientList[B, x, n + 1]
]

enumeratorPolynomial[coeffs_List, x_] := coeffs . x^Range[0, Length[coeffs] - 1]


(* ---- the probability that n copies of a state land in the code space ---- *)

blochVector[s_Wolfram`QuantumFramework`QuantumState] := If[
    s["Qudits"] === 1 && s["Dimensions"] === {2},
    With[{rho = Normal[s["DensityMatrix"]]},
        Simplify[Tr[rho . PauliMatrix[#]] & /@ {1, 2, 3}]
    ],
    $Failed
]
blochVector[r : {_, _, _}] := r
blochVector[_] := $Failed

(* tr[rho^(x)n Pi] = 2^-(n-k) Sum_M sign(M) Prod_q r_(M_q), with r_I = 1, r_X = x,
   r_Y = y, r_Z = z.  The row stores E(x, z) = i^(x z) X^x Z^z, which for x = z = 1 is
   Y itself, so the letter is read straight off the two bits. *)
codeCodespaceProbability[a_Association, r_List] := Module[
    {rows = codeStabilizerGroup[a], n = a["Qubits"], letter},
    letter = <|{0, 0} -> 1, {1, 0} -> r[[1]], {1, 1} -> r[[2]], {0, 1} -> r[[3]]|>;
    Simplify[
        Total[
            MapThread[
                #2 Times @@ Lookup[letter, Transpose[{#1[[1 ;; n]], #1[[n + 1 ;; 2 n]]}]] &,
                {rows, rowSigns[rows]}
            ]
        ] / Length[rows]
    ]
]


(* ---- properties ---- *)

QECCode[a_Association]["WeightEnumerator"] :=
    If[stabilizerGroupAllowedQ[a], codeWeightEnumerator[a], $Failed]
QECCode[a_Association]["WeightEnumerator", x_] :=
    If[stabilizerGroupAllowedQ[a], enumeratorPolynomial[codeWeightEnumerator[a], x], $Failed]

QECCode[a_Association]["SignedWeightEnumerator"] :=
    If[stabilizerGroupAllowedQ[a], codeSignedWeightEnumerator[a], $Failed]
QECCode[a_Association]["SignedWeightEnumerator", x_] :=
    If[stabilizerGroupAllowedQ[a], enumeratorPolynomial[codeSignedWeightEnumerator[a], x], $Failed]

QECCode[a_Association]["DualWeightEnumerator"] :=
    If[stabilizerGroupAllowedQ[a], codeDualWeightEnumerator[a], $Failed]
QECCode[a_Association]["DualWeightEnumerator", x_] :=
    If[stabilizerGroupAllowedQ[a], enumeratorPolynomial[codeDualWeightEnumerator[a], x], $Failed]

QECCode[a_Association]["CodespaceProbability", state_] := With[{r = blochVector[state]},
    Which[
        r === $Failed, Message[QECCode::blochstate, state]; $Failed,
        ! stabilizerGroupAllowedQ[a], $Failed,
        True, codeCodespaceProbability[a, r]
    ]
]
