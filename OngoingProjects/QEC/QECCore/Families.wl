(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[codeCatalog]
PackageScope[classicalHammingMatrix]


(* ============================================================================ *)
(* Named codes and code families.                                               *)
(*                                                                              *)
(* Everything here is a constructor for the same QECCode object, reached through *)
(* the object itself rather than through a separate exported symbol per family:  *)
(* QECCode["SteaneCode"], QECCode["Repetition", 5], QECCode["CSS", hx, hz].      *)
(* The scalable families of roadmap item (f) slot in here the same way:          *)
(* QECCode["Surface", d] and QECCode["Toric", L]; bivariate bicycle codes next.  *)
(* ============================================================================ *)


(* ---- fixed named codes ---- *)

(* PauliStabilizer's named codes give n stabilizers, which fix one codeword; the
   code is the group they generate minus the last one, which is the logical Z
   pinning that particular codeword. *)
$QECNamedCodes := $QECNamedCodes = <|
    "BitFlipCode" -> {"ZZI", "IZZ"},
    "PhaseFlipCode" -> {"XXI", "IXX"},
    "ShorCode" -> {
        "ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ",
        "XXXXXXIII", "IIIXXXXXX"
    },
    "5QubitCode" -> Most[Wolfram`QuantumFramework`PauliStabilizer["5QubitCode"]["Stabilizers"]],
    "SteaneCode" -> Most[Wolfram`QuantumFramework`PauliStabilizer["SteaneCode"]["Stabilizers"]]
|>

QECCode::unknown = "`1` is not a known code or code family. Named codes: `2`. Families: `3`.";

$QECFamilies = {"Repetition", "PhaseRepetition", "DistanceTwo", "Hamming", "CSS", "Surface", "Toric"};

QECCode[name_String] := If[
    KeyExistsQ[$QECNamedCodes, name],
    QECCode[$QECNamedCodes[name]],
    Message[QECCode::unknown, name, StringRiffle[Keys[$QECNamedCodes], ", "], StringRiffle[$QECFamilies, ", "]];
    $Failed
]


(* ---- parametric families ---- *)

QECCode::repsize = "The repetition code needs at least 2 qubits.";
QECCode::evensize = "The [[n, n-2, 2]] family needs an even number of qubits, at least 4.";
QECCode::hammingsize = "The Hamming CSS construction needs r >= 3 for the classical code to be weakly self-dual.";

(* Neighbouring-pair checks: ZZ for bit flips, XX for phase flips (thesis sec. 2.2). *)
repetitionRows[n_Integer, letter_ : "Z"] := Table[
    pauliVector[StringJoin[ReplacePart[ConstantArray["I", n], {i -> letter, i + 1 -> letter}]]],
    {i, n - 1}
]

QECCode["Repetition", n_Integer] := If[
    n >= 2, QECCode[repetitionRows[n, "Z"]], Message[QECCode::repsize]; $Failed
]

QECCode["PhaseRepetition", n_Integer] := If[
    n >= 2, QECCode[repetitionRows[n, "X"]], Message[QECCode::repsize]; $Failed
]

(* The [[n, n-2, 2]] family of thesis sec. 8.1, saturating the quantum Singleton bound. *)
QECCode["DistanceTwo", n_Integer] := If[
    EvenQ[n] && n >= 4,
    QECCode[{StringRepeat["X", n], StringRepeat["Z", n]}],
    Message[QECCode::evensize]; $Failed
]

(* Columns are the nonzero binary words of length r: the classical [2^r-1, 2^r-1-r, 3]
   Hamming code. *)
classicalHammingMatrix[r_Integer ? Positive] := Transpose[Table[IntegerDigits[j, 2, r], {j, 2^r - 1}]]

(* The CSS code built on it, [[2^r-1, 2^r-1-2r, 3]], which is the Steane code at r = 3. *)
QECCode["Hamming", r_Integer] := If[
    r >= 3, QECCode["CSS", classicalHammingMatrix[r]], Message[QECCode::hammingsize]; $Failed
]


(* ---- topological families (roadmap item (f), Got26 ch. 17) ---- *)

QECCode::surfacesize = "The surface code needs a distance d >= 2.";
QECCode::toricsize = "The toric code needs a lattice size L >= 2.";

(* A Pauli of one letter on the given (1-based) qubits, as a row. *)
letterRow[n_Integer, letter_String, qubits_List] :=
    pauliVector[StringJoin[ReplacePart[ConstantArray["I", n], Thread[qubits -> letter]]]]

(* THE ROTATED SURFACE CODE, [[d^2, 1, d]]: the layout of the memory experiments in
   the literature (and of Stim's "surface_code:rotated_memory" circuits).  Data qubits
   sit on a d x d grid, (r, c) -> r d + c + 1, row by row.  A check sits on each corner
   (a, b), 0 <= a, b <= d, of that grid and acts on the up-to-four data qubits around
   it; in the bulk the type alternates like a checkerboard, X when a + b is even.
   On the boundary only half the corners carry a weight-2 check: X checks on the top
   and bottom edges, Z checks on the left and right ones.  That choice is what makes
   the two kinds of boundary of Got26 sec. 17.1.4 -- logical X runs between the two
   Z-checked sides, logical Z between the two X-checked ones -- and the count works
   out to (d - 1)^2 bulk plus 2 (d - 1) boundary checks, d^2 - 1 in all, so one
   logical qubit.  X checks are listed first, then Z, each in reading order. *)
surfaceChecks[d_Integer] := Module[{at, corners, onBoundary, keep, xs, zs},
    at[a_, b_] := Select[{{a - 1, b - 1}, {a - 1, b}, {a, b - 1}, {a, b}}, 0 <= #[[1]] < d && 0 <= #[[2]] < d &];
    corners = Tuples[Range[0, d], 2];
    keep[{a_, b_}] := Which[
        1 <= a <= d - 1 && 1 <= b <= d - 1, True,
        (a == 0 || a == d) && 1 <= b <= d - 1, EvenQ[a + b],
        (b == 0 || b == d) && 1 <= a <= d - 1, OddQ[a + b],
        True, False
    ];
    corners = Select[corners, keep];
    xs = Select[corners, EvenQ[Total[#]] &];
    zs = Select[corners, OddQ[Total[#]] &];
    Join[
        letterRow[d^2, "X", (#[[1]] d + #[[2]] + 1) & /@ at @@ #] & /@ xs,
        letterRow[d^2, "Z", (#[[1]] d + #[[2]] + 1) & /@ at @@ #] & /@ zs
    ]
]

(* THE CHECK ORDER.  A fault on a check's ancilla halfway through its CNOTs leaves the
   check's letter on the qubits not yet touched: a pair, for a weight-four check.  If
   that pair lies along the logical operator of the same letter, two faults' worth of
   data error comes from one fault, and the circuit loses distance in that basis.
   Logical X runs down a column and logical Z along a row, so X checks go row by row
   (their hook is a horizontal pair) and Z checks column by column (a vertical pair):
   each hook is perpendicular to the logical it could have helped build.  This is the
   "N" and "Z" order of the rotated-surface-code literature. *)
surfaceCheckOrder[d_Integer] := Module[{n = d^2, rows = surfaceChecks[d]},
    Function[row, With[{support = Flatten[Position[BitOr[row[[1 ;; n]], row[[n + 1 ;; 2 n]]], 1]]},
        If[ Total[row[[1 ;; n]]] > 0,
            support,
            SortBy[support, {Mod[# - 1, d], Quotient[# - 1, d]} &]
        ]
    ]] /@ rows
]

QECCode["Surface", d_Integer] := If[
    d >= 2,
    Replace[QECCode[surfaceChecks[d]], QECCode[a_Association] :> QECCode[Append[a, "CheckOrder" -> surfaceCheckOrder[d]]]],
    Message[QECCode::surfacesize]; $Failed
]

(* THE TORIC CODE, [[2 L^2, 2, L]] (Got26 sec. 17.1.1): qubits on the edges of an
   L x L square lattice on a torus, horizontal edges first, h(r, c) -> r L + c + 1,
   then vertical, v(r, c) -> L^2 + r L + c + 1, indices taken mod L.  Each face gives
   Z on its four edges, each vertex X on its four.  The product of all faces is the
   identity and so is the product of all vertices (eq. 17.3), so one of each is
   dependent and is left out: 2 L^2 - 2 generators, two logical qubits, one per
   non-contractible direction. *)
toricChecks[l_Integer] := Module[{h, v, faces, vertices},
    h[r_, c_] := Mod[r, l] l + Mod[c, l] + 1;
    v[r_, c_] := l^2 + Mod[r, l] l + Mod[c, l] + 1;
    faces = Flatten[Table[{h[r, c], h[r + 1, c], v[r, c], v[r, c + 1]}, {r, 0, l - 1}, {c, 0, l - 1}], 1];
    vertices = Flatten[Table[{h[r, c], h[r, c - 1], v[r, c], v[r - 1, c]}, {r, 0, l - 1}, {c, 0, l - 1}], 1];
    Join[
        letterRow[2 l^2, "X", #] & /@ Most[vertices],
        letterRow[2 l^2, "Z", #] & /@ Most[faces]
    ]
]

QECCode["Toric", l_Integer] := If[
    l >= 2, QECCode[toricChecks[l]], Message[QECCode::toricsize]; $Failed
]


(* ---- the catalog ---- *)

$QECCatalogExamples = {
    {"Repetition", 5}, {"PhaseRepetition", 5}, {"DistanceTwo", 6}, {"Hamming", 4},
    {"Surface", 3}, {"Toric", 3}
};

(* The catalog is a question about the head, so it is asked of the head. *)
QECCode["Catalog"] := codeCatalog[]

codeCatalog[] := Dataset @ Association @ Join[
    Table[
        name -> AssociationThread[{"n", "k", "d"}, QECCode[name]["Parameters"]],
        {name, Keys[$QECNamedCodes]}
    ],
    Table[
        StringTemplate["QECCode[\"``\", ``]"][First[spec], Last[spec]] ->
            AssociationThread[{"n", "k", "d"}, (QECCode @@ spec)["Parameters"]],
        {spec, $QECCatalogExamples}
    ]
]
