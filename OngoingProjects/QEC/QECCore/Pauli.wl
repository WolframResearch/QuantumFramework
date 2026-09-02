(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECPauliVector]
PackageExport[QECPauliString]
PackageExport[QECPauliQ]
PackageExport[QECPauliWeight]
PackageExport[QECPauliCommuteQ]
PackageExport[QECPauliProduct]
PackageExport[QECPauliPhase]

PackageScope[pauliQubits]
PackageScope[symplecticPart]
PackageScope[phasePart]
PackageScope[symplecticProduct]
PackageScope[pauliIdentity]
PackageScope[weightKVectors]
PackageScope[weightOneVectors]
PackageScope[$pauliLetterXZ]


(* ============================================================================ *)
(* The Pauli layer: rows of the form  {x1..xn, z1..zn, e}.                      *)
(*                                                                              *)
(* The layout is the one the stabilizer engine already uses (see PauliRow in    *)
(* Kernel/Stabilizer/Conversions.m, which returns Join[xbits, zbits, {phase}]), *)
(* so codes and the engine share a row format and the bit-packed fast path of   *)
(* Kernel/Stabilizer/Packed.m stays available to us later without a rewrite.    *)
(*                                                                              *)
(* One deliberate generalisation: the engine's phase column is a single bit,    *)
(* because a stabilizer tableau only ever holds Hermitian rows.  Here e lives   *)
(* in Z4 and means an overall factor i^e, with the per-qubit convention         *)
(*                                                                              *)
(*     E(x, z) = i^(x z) X^x Z^z    so that  E(1,1) = Y  is Hermitian.          *)
(*                                                                              *)
(* Z4 is what makes the Pauli group closed under multiplication: X.Z = -i Y is  *)
(* not expressible with a sign bit.  The prototype carried signs as a "-" prefix *)
(* parsed off the string and never multiplied Paulis at all, which is why its   *)
(* CorrectionCycle could only certify residuals up to phase.  Hermitian rows    *)
(* have even e, and e/2 is exactly the engine's phase bit.                      *)
(* ============================================================================ *)

QECPauliVector::usage = "QECPauliVector[p] gives the Pauli row {x1..xn, z1..zn, e} of a Pauli, where the operator is i^e times the product of X^x_j Z^z_j taken Hermitian on each qubit. p may be a Pauli string such as \"XZZXI\" or \"-iXY\", a row (returned unchanged), or a list of strings.";

QECPauliString::usage = "QECPauliString[v] gives the Pauli string of a Pauli row, with a leading -, i or -i when the phase calls for one. Strings pass through unchanged and lists are mapped over.";

QECPauliQ::usage = "QECPauliQ[p] gives True if p is a Pauli string or a Pauli row.";

QECPauliWeight::usage = "QECPauliWeight[p] gives the number of qubits on which the Pauli acts nontrivially. The phase does not count.";

QECPauliCommuteQ::usage = "QECPauliCommuteQ[p, q] gives True if the two Paulis commute, which is the symplectic product of their rows being zero.";

QECPauliProduct::usage = "QECPauliProduct[p1, p2, ...] gives the product of the Paulis as a Pauli row, carrying the phase in Z4: QECPauliProduct[\"X\", \"Z\"] is -iY.";

QECPauliPhase::usage = "QECPauliPhase[p] gives the phase of a Pauli as an element of Z4: the operator carries an overall factor i^QECPauliPhase[p]. Hermitian Paulis give 0 or 2.";


$pauliLetterXZ = <|"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}|>;

$xzPauliLetter = <|{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"|>;

$phasePrefix = <|"" -> 0, "i" -> 1, "-" -> 2, "-i" -> 3|>;

$prefixPhase = <|0 -> "", 1 -> "i", 2 -> "-", 3 -> "-i"|>;


(* ---- row accessors ---- *)

pauliQubits[v_List] := (Length[v] - 1) / 2

symplecticPart[v_List] := Most[v]

phasePart[v_List] := Last[v]

pauliIdentity[n_Integer] := Append[gf2Zero[2 n], 0]


(* ---- strings in, strings out ---- *)

QECPauliVector::invalid = "`1` is not a Pauli string: expected an optional -, i or -i followed by letters I, X, Y, Z.";

QECPauliVector::badrow = "`1` is not a Pauli row: expected a list of an odd number of integers, {x1..xn, z1..zn, e}.";

QECPauliVector[s_String] := Module[{prefix, body, chars, pairs},
    prefix = Replace[StringCases[s, StartOfString ~~ p : ("-i" | "-" | "i") :> p], {{p_} :> p, _ -> ""}];
    body = StringDrop[s, StringLength[prefix]];
    chars = Characters[body];
    If[ chars === {} || ! SubsetQ[Keys[$pauliLetterXZ], DeleteDuplicates[chars]],
        Message[QECPauliVector::invalid, s]; Return[$Failed]
    ];
    pairs = Lookup[$pauliLetterXZ, chars];
    Join[pairs[[All, 1]], pairs[[All, 2]], {$phasePrefix[prefix]}]
]

QECPauliVector[v : {__Integer}] := If[
    OddQ[Length[v]] && SubsetQ[{0, 1}, DeleteDuplicates[Most[v]]],
    MapAt[Mod[#, 4] &, v, -1],
    Message[QECPauliVector::badrow, v]; $Failed
]

QECPauliVector[ps : {__String}] := QECPauliVector /@ ps


QECPauliString[v : {__Integer}] := Module[{n = pauliQubits[v]},
    $prefixPhase[Mod[Last[v], 4]] <> StringJoin[Lookup[$xzPauliLetter, Transpose[{v[[1 ;; n]], v[[n + 1 ;; 2 n]]}]]]
]

QECPauliString[s_String] := s

QECPauliString[l : {__List}] := QECPauliString /@ l


(* Checked structurally rather than by parsing and silencing the message. *)
QECPauliQ[s_String] := StringMatchQ[s, ("-i" | "-" | "i" | "") ~~ ("I" | "X" | "Y" | "Z") ..]

QECPauliQ[v : {__Integer}] := OddQ[Length[v]] && SubsetQ[{0, 1}, DeleteDuplicates[Most[v]]]

QECPauliQ[_] := False


(* ---- the algebra ---- *)

QECPauliWeight[p_] := With[{v = QECPauliVector[p]},
    With[{n = pauliQubits[v]}, Count[v[[1 ;; n]] + v[[n + 1 ;; 2 n]], _ ? Positive]]
]


(* The symplectic form  x.z' + z.x'  mod 2: zero exactly when the two Paulis commute. *)
symplecticProduct[u_List, v_List] := With[{n = pauliQubits[u]},
    Mod[u[[1 ;; n]] . v[[n + 1 ;; 2 n]] + u[[n + 1 ;; 2 n]] . v[[1 ;; n]], 2]
]

QECPauliCommuteQ::size = "Paulis act on different numbers of qubits.";

QECPauliCommuteQ[p_, q_] := With[{u = QECPauliVector[p], v = QECPauliVector[q]},
    If[ Length[u] =!= Length[v],
        Message[QECPauliCommuteQ::size]; $Failed,
        symplecticProduct[u, v] === 0
    ]
]

(* E(x1,z1) E(x2,z2) = i^c E(x3,z3) with x3 = x1+x2, z3 = z1+z2 mod 2 and
   c = x1 z1 + x2 z2 - x3 z3 + 2 z1 x2, summed over qubits.  Derivation: commute
   Z^z1 past X^x2 (the factor (-1)^(z1 x2)), then re-Hermitise the result. *)
pauliTimes[u_List, v_List] := Module[{n = pauliQubits[u], x1, z1, x2, z2, x3, z3, c},
    x1 = u[[1 ;; n]]; z1 = u[[n + 1 ;; 2 n]];
    x2 = v[[1 ;; n]]; z2 = v[[n + 1 ;; 2 n]];
    x3 = Mod[x1 + x2, 2];
    z3 = Mod[z1 + z2, 2];
    c = Total[x1 z1 + x2 z2 - x3 z3 + 2 z1 x2];
    Join[x3, z3, {Mod[Last[u] + Last[v] + c, 4]}]
]

QECPauliProduct::size = "Paulis act on different numbers of qubits.";

QECPauliProduct[ps__] := Module[{vs = QECPauliVector /@ {ps}},
    If[ Length[DeleteDuplicates[Length /@ vs]] =!= 1,
        Message[QECPauliProduct::size]; Return[$Failed]
    ];
    Fold[pauliTimes, vs]
]

QECPauliPhase[p_] := Mod[Last[QECPauliVector[p]], 4]


(* ---- enumeration of low-weight errors ---- *)

(* All Pauli rows of weight exactly k on n qubits, phase 0. Generated as vectors
   rather than strings: the distance search calls this on every weight. *)
weightKVectors[n_Integer ? Positive, 0] := {pauliIdentity[n]}

weightKVectors[n_Integer ? Positive, k_Integer ? Positive] /; k <= n := Flatten[
    Table[
        Join[
            Normal @ SparseArray[Thread[pos -> letters[[All, 1]]], n],
            Normal @ SparseArray[Thread[pos -> letters[[All, 2]]], n],
            {0}
        ],
        {pos, Subsets[Range[n], {k}]},
        {letters, Tuples[Values[$pauliLetterXZ][[2 ;; 4]], k]}
    ],
    1
]

weightKVectors[_Integer, _Integer] := {}

weightOneVectors[n_Integer ? Positive] := weightKVectors[n, 1]
