(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[codeSpanData]
PackageScope[codeStabilizerElement]
PackageScope[codeLogicalVectors]
PackageScope[codeMinimumLogical]


(* ============================================================================ *)
(* Structure of the code: standard form, logical operators, distance.           *)
(*                                                                              *)
(* Gottesman's thesis sec. 4.1 for the standard form and the closed formulas for *)
(* the logical operators; sec. 3.2 for the distance as the minimum weight in     *)
(* N(S) \ S.  The mathematics is the prototype's, which was checked by hand;     *)
(* what changes is the bookkeeping around it.                                    *)
(* ============================================================================ *)


(* ---- membership, with signs ---- *)

(* [M | I] reduced over GF(2) gives [R | T] with R = T.M in echelon form.  T is
   what lets us name the actual product of generators behind a membership, and
   therefore its phase -- the prototype could only answer membership up to sign. *)
codeSpanData[a_Association] := codeSpanData[a] = Module[{m = codeStabilizerCount[a], n = a["Qubits"], red},
    red = RowReduce[Join[a["CheckMatrix"], IdentityMatrix[m], 2], Modulus -> 2];
    <|
        "Rows" -> red[[All, 1 ;; 2 n]],
        "Transform" -> red[[All, 2 n + 1 ;; 2 n + m]],
        "Pivots" -> (First[FirstPosition[#, 1]] & /@ red[[All, 1 ;; 2 n]])
    |>
]

(* The stabilizer group element whose symplectic part is v, as a full Pauli row
   with its true phase, or Missing if v is not in the span. *)
codeStabilizerElement[a_Association, v_List] := Module[{sd, rows, piv, trans, residual, coeff, gens, used},
    sd = codeSpanData[a];
    {rows, piv, trans} = Lookup[sd, {"Rows", "Pivots", "Transform"}];
    residual = symplecticPart[QECPauliVector[v]];
    coeff = gf2Zero[codeStabilizerCount[a]];
    Do[
        If[ residual[[piv[[i]]]] === 1,
            residual = Mod[residual + rows[[i]], 2];
            coeff = Mod[coeff + trans[[i]], 2]
        ],
        {i, Length[rows]}
    ];
    If[residual =!= gf2Zero[2 a["Qubits"]], Return[Missing["NotInStabilizerGroup"]]];
    gens = codeVectors[a];
    used = Pick[gens, coeff, 1];
    If[used === {}, pauliIdentity[a["Qubits"]], QECPauliProduct @@ used]
]

(* Membership in the stabilizer group proper: the symplectic part must be in the
   span AND the phase must agree.  -M for a stabilizer M is in N(S) but not in S. *)
codeStabilizerMemberQ[a_Association, p_] := With[{v = QECPauliVector[p], el = codeStabilizerElement[a, QECPauliVector[p]]},
    ! MissingQ[el] && phasePart[el] === phasePart[v]
]

(* In the normaliser but not the stabilizer group: a logical operator.  Sign is
   irrelevant here -- P and -P are logical together -- so this compares spans. *)
codeLogicalPauliQ[a_Association, p_] := With[{v = QECPauliVector[p]},
    Length[v] === 2 a["Qubits"] + 1 &&
    codeSyndromeVector[a, v] === gf2Zero[codeStabilizerCount[a]] &&
    MissingQ[codeStabilizerElement[a, v]]
]


(* ---- completion to a full stabilizer state ---- *)

QECCode::incomplete = "Could not complete the generator set to a full stabilizer state.";

(* Extend the m code generators to n commuting independent ones, pinning one
   fiducial encoded state.  Candidates come from the commutant of the code. *)
codeCompletedGenerators[a_Association] := codeCompletedGenerators[a] = Module[
    {n = a["Qubits"], mat = a["CheckMatrix"], swapped, commutant, basis, extras},
    swapped = Join[mat[[All, n + 1 ;; 2 n]], mat[[All, 1 ;; n]], 2];
    commutant = gf2NullSpace[swapped];
    basis = gf2Basis[mat];
    extras = {};
    Do[
        If[ Length[extras] + Length[mat] < n &&
            gf2IndependentQ[gf2Basis[Join[mat, extras]], c] &&
            AllTrue[extras, symplecticProduct[Append[#, 0], Append[c, 0]] === 0 &],
            extras = Append[extras, c]
        ],
        {c, commutant}
    ];
    If[ Length[extras] + Length[mat] =!= n,
        Message[QECCode::incomplete]; Return[$Failed]
    ];
    Join[codeVectors[a], Append[#, 0] & /@ extras]
]


(* ---- standard form (thesis sec. 4.1) ---- *)

(* Gaussian elimination over GF(2) with generator products as row operations and
   qubit relabelling as paired column swaps (column q and column n+q move
   together), bringing the check matrix to
       [ I  A1 A2 | B  C1 C2 ]
       [ 0  0  0  | D  I  E  ]
   and returning the qubit permutation used and r, the rank of the X block. *)
codeStandardForm[a_Association] := codeStandardForm[a] = Module[
    {mat = a["CheckMatrix"], n = a["Qubits"], m = codeStabilizerCount[a], perm, r, s, pivot, swapQubits},

    perm = Range[n];

    swapQubits[i_, j_] := (
        mat[[All, {i, j}]] = mat[[All, {j, i}]];
        mat[[All, {n + i, n + j}]] = mat[[All, {n + j, n + i}]];
        perm[[{i, j}]] = perm[[{j, i}]]
    );

    (* First block: eliminate in the X half, columns 1..n. *)
    r = 0;
    While[r < m,
        pivot = firstPivot[mat, r + 1, m, r + 1, n, 0];
        If[MissingQ[pivot], Break[]];
        If[pivot[[2]] =!= r + 1, swapQubits[r + 1, pivot[[2]]]];
        If[pivot[[1]] =!= r + 1, mat[[{pivot[[1]], r + 1}]] = mat[[{r + 1, pivot[[1]]}]]];
        Do[
            If[j =!= r + 1 && mat[[j, r + 1]] === 1, mat[[j]] = Mod[mat[[j]] + mat[[r + 1]], 2]],
            {j, m}
        ];
        r++
    ];

    (* Second block: eliminate in the Z half, columns n+r+1..2n, rows below r. *)
    s = 0;
    While[r + s < m,
        pivot = firstPivot[mat, r + s + 1, m, r + s + 1, n, n];
        If[MissingQ[pivot], Break[]];
        If[pivot[[2]] =!= r + s + 1, swapQubits[r + s + 1, pivot[[2]]]];
        If[pivot[[1]] =!= r + s + 1, mat[[{pivot[[1]], r + s + 1}]] = mat[[{r + s + 1, pivot[[1]]}]]];
        Do[
            If[j =!= r + s + 1 && j > r && mat[[j, n + r + s + 1]] === 1, mat[[j]] = Mod[mat[[j]] + mat[[r + s + 1]], 2]],
            {j, m}
        ];
        s++
    ];

    <|"Matrix" -> mat, "QubitPermutation" -> perm, "XRank" -> r|>
]

(* First {row, column} with a 1, scanning columns outermost as the thesis procedure
   requires; offset selects the X half (0) or the Z half (n). *)
firstPivot[mat_, rowFrom_, rowTo_, colFrom_, colTo_, offset_] := Catch[
    Do[
        If[mat[[i, offset + c]] === 1, Throw[{i, c}]],
        {c, colFrom, colTo}, {i, rowFrom, rowTo}
    ];
    Throw[Missing["NoPivot"]]
]


(* ---- logical operators (thesis sec. 4.1, closed form on the standard form) ---- *)

(* Xbar = (0 E^T I | E^T C1^T + C2^T 0 0),  Zbar = (0 0 0 | A2^T 0 I), then undo
   the qubit permutation the standard form introduced. *)
codeLogicalVectors[a_Association] := codeLogicalVectors[a] = Module[
    {n = a["Qubits"], m = codeStabilizerCount[a], k, sf, mat, perm, r, mid, a2, c1, c2, e, unpermute, xRows, zRows},

    k = codeLogicalQubits[a];
    If[k === 0, Return[<|"X" -> {}, "Z" -> {}|>]];

    sf = codeStandardForm[a];
    {mat, perm, r} = Lookup[sf, {"Matrix", "QubitPermutation", "XRank"}];
    mid = m - r;

    a2 = mat[[1 ;; r, r + mid + 1 ;; n]];
    c1 = mat[[1 ;; r, n + r + 1 ;; n + r + mid]];
    c2 = mat[[1 ;; r, n + r + mid + 1 ;; 2 n]];
    e = mat[[r + 1 ;; m, n + r + mid + 1 ;; 2 n]];

    unpermute[v_] := Module[{x = v[[1 ;; n]], z = v[[n + 1 ;; 2 n]], xo, zo},
        xo = zo = gf2Zero[n];
        Do[xo[[perm[[q]]]] = x[[q]]; zo[[perm[[q]]]] = z[[q]], {q, n}];
        Join[xo, zo, {0}]
    ];

    xRows = Table[
        Join[
            gf2Zero[r],
            Table[e[[l, i]], {l, mid}],
            IdentityMatrix[k][[i]],
            Table[Mod[Sum[e[[l, i]] c1[[j, l]], {l, mid}] + c2[[j, i]], 2], {j, r}],
            gf2Zero[mid],
            gf2Zero[k]
        ],
        {i, k}
    ];

    zRows = Table[
        Join[gf2Zero[n], Table[a2[[j, i]], {j, r}], gf2Zero[mid], IdentityMatrix[k][[i]]],
        {i, k}
    ];

    <|"X" -> (unpermute /@ xRows), "Z" -> (unpermute /@ zRows)|>
]

codeLogicalOperators[a_Association] := Map[QECPauliString, codeLogicalVectors[a], {2}]


(* ---- distance ---- *)

(* The minimum weight over N(S) \ S, with its witness.
   Two things the prototype did not do.  First, the logical operators from the
   standard form are themselves in N(S) \ S, so their minimum weight is an upper
   bound and the search never has to run past it -- on the Steane code that stops
   the scan at weight 3 instead of walking to 7.  Second, the syndromes of a whole
   weight class are computed as one matrix product rather than one call per
   candidate, and only the zero-syndrome survivors are tested for membership. *)
codeMinimumLogical[a_Association] := codeMinimumLogical[a] = Module[
    {n = a["Qubits"], logicals, weights, bound, witness, swapped, candidates, syndromes, kernel, found},

    If[codeLogicalQubits[a] === 0, Return[{Infinity, Missing["NoLogicalOperators"]}]];

    logicals = Join @@ Values[codeLogicalVectors[a]];
    weights = QECPauliWeight /@ logicals;
    bound = Min[weights];
    witness = logicals[[First[FirstPosition[weights, bound]]]];

    swapped = With[{mat = a["CheckMatrix"]}, Join[mat[[All, n + 1 ;; 2 n]], mat[[All, 1 ;; n]], 2]];

    Do[
        candidates = weightKVectors[n, w];
        syndromes = Mod[(symplecticPart /@ candidates) . Transpose[swapped], 2];
        kernel = Pick[candidates, Total /@ syndromes, 0];
        found = SelectFirst[kernel, MissingQ[codeStabilizerElement[a, #]] &];
        If[! MissingQ[found], Return[{w, found}, Module]],
        {w, 1, bound - 1}
    ];

    {bound, witness}
]

codeDistance[a_Association] := First[codeMinimumLogical[a]]

codeMinimumWeightLogical[a_Association] := With[{res = Last[codeMinimumLogical[a]]},
    If[MissingQ[res], res, QECPauliString[res]]
]


(* ---- CSS ---- *)

(* A generator is X-type or Z-type when one of its halves vanishes; a Y anywhere
   puts a 1 in both halves of the same column, so this also rules Y out. *)
codeCSSQ[a_Association] := With[{n = a["Qubits"]},
    AllTrue[a["CheckMatrix"], Total[#[[1 ;; n]]] === 0 || Total[#[[n + 1 ;; 2 n]]] === 0 &]
]
