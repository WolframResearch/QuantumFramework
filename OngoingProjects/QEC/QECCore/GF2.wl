(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[gf2Basis]
PackageScope[gf2Rank]
PackageScope[gf2NullSpace]
PackageScope[gf2Solve]
PackageScope[gf2Reduce]
PackageScope[gf2MemberQ]
PackageScope[gf2IndependentQ]
PackageScope[gf2Zero]


(* ============================================================================ *)
(* Linear algebra over GF(2), written once.                                     *)
(*                                                                              *)
(* The prototype called MatrixRank[Append[rows, v], Modulus -> 2] every time it *)
(* needed to ask "is this Pauli in the stabilizer group?", re-reducing the whole *)
(* check matrix per query.  Here a matrix is reduced to a row-echelon basis once *)
(* (gf2Basis) and every later membership / residual query reduces the vector     *)
(* against that basis in O(rank * n) with no further elimination.  QECCode caches *)
(* the basis of its check matrix, so the cost is paid once per code.             *)
(*                                                                              *)
(* A "basis" here is the pair {rows, pivots}: rows is the reduced matrix with    *)
(* zero rows dropped, pivots[[i]] is the column of the leading 1 of rows[[i]].   *)
(* ============================================================================ *)

gf2Zero[n_Integer] := ConstantArray[0, n]


(* ---- reduction to an echelon basis ---- *)

gf2Basis[{}] := {{}, {}}

gf2Basis[m_ ? MatrixQ] := With[{rows = DeleteCases[RowReduce[m, Modulus -> 2], {(0) ..}]},
    {rows, If[rows === {}, {}, First[FirstPosition[#, 1]] & /@ rows]}
]


gf2Rank[{}] := 0

gf2Rank[m_ ? MatrixQ] := Length[First[gf2Basis[m]]]


(* ---- queries against a precomputed basis ---- *)

(* Residual of v after subtracting off every basis row whose pivot v occupies. *)
gf2Reduce[{rows_, pivots_}, v_List] := Fold[
    If[#1[[pivots[[#2]]]] === 1, Mod[#1 + rows[[#2]], 2], #1] &,
    v,
    Range[Length[rows]]
]

gf2MemberQ[basis_, v_List] := gf2Reduce[basis, v] === gf2Zero[Length[v]]

(* True if v is outside the span, i.e. it would extend the basis. *)
gf2IndependentQ[basis_, v_List] := ! gf2MemberQ[basis, v]


(* ---- kernel and solving ---- *)

gf2NullSpace[{}] := {}

gf2NullSpace[m_ ? MatrixQ] := NullSpace[m, Modulus -> 2]

(* No Quiet: the Stabilizer subsystem's house rule is that a message either matters
   or should not be raised.  An inconsistent system is an ordinary outcome here (the
   encoder asks whether a sign pattern is reachable), so it is tested for rather than
   silenced: b must lie in the row space of m. *)
gf2Solve[m_ ? MatrixQ, b_List] := With[{basis = gf2Basis[Transpose[m]]},
    If[ gf2MemberQ[basis, b],
        LinearSolve[m, b, Modulus -> 2],
        $Failed
    ]
]
