(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECPauli]
PackageExport[QECPauliQ]

PackageScope[pauliVector]
PackageScope[pauliString]
PackageScope[pauliWeight]
PackageScope[pauliCommuteQ]
PackageScope[pauliProduct]
PackageScope[pauliPhase]

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

QECPauliQ::usage = "QECPauliQ[p] gives True if p is a Pauli string or a Pauli row.";


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

QECPauli::invalid = "`1` is not a Pauli string: expected an optional -, i or -i followed by letters I, X, Y, Z.";

QECPauli::badrow = "`1` is not a Pauli row: expected a list of an odd number of integers, {x1..xn, z1..zn, e}.";

pauliVector[s_String] := Module[{prefix, body, chars, pairs},
    prefix = Replace[StringCases[s, StartOfString ~~ p : ("-i" | "-" | "i") :> p], {{p_} :> p, _ -> ""}];
    body = StringDrop[s, StringLength[prefix]];
    chars = Characters[body];
    If[ chars === {} || ! SubsetQ[Keys[$pauliLetterXZ], DeleteDuplicates[chars]],
        Message[QECPauli::invalid, s]; Return[$Failed]
    ];
    pairs = Lookup[$pauliLetterXZ, chars];
    Join[pairs[[All, 1]], pairs[[All, 2]], {$phasePrefix[prefix]}]
]

pauliVector[v : {__Integer}] := If[
    OddQ[Length[v]] && SubsetQ[{0, 1}, DeleteDuplicates[Most[v]]],
    MapAt[Mod[#, 4] &, v, -1],
    Message[QECPauli::badrow, v]; $Failed
]

pauliVector[ps : {__String}] := pauliVector /@ ps


pauliString[v : {__Integer}] := Module[{n = pauliQubits[v]},
    $prefixPhase[Mod[Last[v], 4]] <> StringJoin[Lookup[$xzPauliLetter, Transpose[{v[[1 ;; n]], v[[n + 1 ;; 2 n]]}]]]
]

pauliString[s_String] := s

pauliString[l : {__List}] := pauliString /@ l


(* Checked structurally rather than by parsing and silencing the message. *)
QECPauliQ[s_String] := StringMatchQ[s, ("-i" | "-" | "i" | "") ~~ ("I" | "X" | "Y" | "Z") ..]

QECPauliQ[v : {__Integer}] := OddQ[Length[v]] && SubsetQ[{0, 1}, DeleteDuplicates[Most[v]]]

QECPauliQ[_] := False


(* ---- the algebra ---- *)

pauliWeight[p_] := With[{v = pauliVector[p]},
    With[{n = pauliQubits[v]}, Count[v[[1 ;; n]] + v[[n + 1 ;; 2 n]], _ ? Positive]]
]


(* The symplectic form  x.z' + z.x'  mod 2: zero exactly when the two Paulis commute. *)
symplecticProduct[u_List, v_List] := With[{n = pauliQubits[u]},
    Mod[u[[1 ;; n]] . v[[n + 1 ;; 2 n]] + u[[n + 1 ;; 2 n]] . v[[1 ;; n]], 2]
]

QECPauli::size = "Paulis act on different numbers of qubits.";

pauliCommuteQ[p_, q_] := With[{u = pauliVector[p], v = pauliVector[q]},
    If[ Length[u] =!= Length[v],
        Message[QECPauli::size]; $Failed,
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


pauliProduct[ps__] := Module[{vs = pauliVector /@ {ps}},
    If[ Length[DeleteDuplicates[Length /@ vs]] =!= 1,
        Message[QECPauli::size]; Return[$Failed]
    ];
    Fold[pauliTimes, vs]
]

pauliPhase[p_] := Mod[Last[pauliVector[p]], 4]


(* ---- the object ---- *)

(* QECPauli is the one public face of everything above.  It WRAPS a row; it does
   not replace it (constraint C3 of the redesign plan): every function in the layer
   keeps taking and returning raw rows, because the algebra runs on them, and the
   object exists for the reader.  The seven former verbs are its properties:

       QECPauliVector[p]      ->  QECPauli[p]["Vector"]
       QECPauliString[p]      ->  QECPauli[p]["String"]
       QECPauliWeight[p]      ->  QECPauli[p]["Weight"]
       QECPauliPhase[p]       ->  QECPauli[p]["Phase"]
       QECPauliCommuteQ[p, q] ->  QECPauli[p]["CommuteQ", q]
       QECPauliProduct[p, q]  ->  QECPauli[p]["Product", q]   or   QECPauli[p] ** QECPauli[q]
       QECPauliQ[p]           ->  stays, as the guard

   The canonical form is QECPauli[row], phase reduced mod 4, which is also the
   shape the constructor settles into, so an object is a row with a head on it and
   pattern-matching on it costs nothing.  Every row function accepts the object as
   well, so a QECPauli can go anywhere a string or a row could: code["Syndrome", P]. *)

QECPauli::usage = "QECPauli[p] represents a Pauli operator, given as a string such as \"XZZXI\" or \"-iXY\" or as a row {x1..xn, z1..zn, e} meaning i^e times the product of X^x_j Z^z_j taken Hermitian on each qubit.\nP[\"String\"], P[\"Vector\"], P[\"Weight\"], P[\"Phase\"], P[\"Qubits\"], P[\"Support\"], P[\"HermitianQ\"], P[\"Matrix\"] and P[\"QuantumOperator\"] give its forms; P[\"CommuteQ\", q] and P[\"Product\", q, ...] relate it to others, and P ** Q multiplies.\nQECPauli[{p1, p2, ...}] maps over a list. P[\"Properties\"] lists the properties.";

QECPauli::noprop = "`1` is not a property of QECPauli. Use P[\"Properties\"] for the list.";

pauliRowQ[v_List] := VectorQ[v, IntegerQ] && OddQ[Length[v]] && Length[v] >= 3 && SubsetQ[{0, 1}, DeleteDuplicates[Most[v]]]
pauliRowQ[_] := False

QECPauli[p_QECPauli] := p
QECPauli[s_String] := Replace[pauliVector[s], v_List :> QECPauli[v]]
QECPauli[v : {__Integer}] /; ! pauliRowQ[v] := (Message[QECPauli::badrow, v]; $Failed)
QECPauli[v : {__Integer}] /; pauliRowQ[v] && ! 0 <= Last[v] <= 3 := QECPauli[MapAt[Mod[#, 4] &, v, -1]]
QECPauli[l : {(_String | _QECPauli | {__Integer}) ..}] /; ! VectorQ[l, IntegerQ] := QECPauli /@ l

pauliVector[QECPauli[v_List]] := v
pauliString[QECPauli[v_List]] := pauliString[v]
QECPauliQ[QECPauli[v_List]] := pauliRowQ[v]

$pauliProperties = {
    "String", "Vector", "Weight", "Phase", "Qubits", "Support", "HermitianQ",
    "Matrix", "QuantumOperator", "CommuteQ", "Product", "Properties"
};

QECPauli[_List]["Properties"] := $pauliProperties
QECPauli[v_List]["String"] := pauliString[v]
QECPauli[v_List]["Vector"] := v
QECPauli[v_List]["Weight"] := pauliWeight[v]
QECPauli[v_List]["Phase"] := pauliPhase[v]
QECPauli[v_List]["Qubits"] := pauliQubits[v]
QECPauli[v_List]["Support"] := With[{n = pauliQubits[v]}, Flatten @ Position[v[[1 ;; n]] + v[[n + 1 ;; 2 n]], _ ? Positive]]
QECPauli[v_List]["HermitianQ"] := EvenQ[Last[v]]
QECPauli[v_List]["Matrix"] := pauliRowMatrix[v]
QECPauli[v_List]["QuantumOperator"] := Wolfram`QuantumFramework`QuantumOperator[pauliRowMatrix[v], Range[pauliQubits[v]]]
QECPauli[v_List]["CommuteQ", q_] := pauliCommuteQ[v, q]
QECPauli[v_List]["Product", qs__] := Replace[pauliProduct[v, qs], w_List :> QECPauli[w]]
QECPauli[_List][prop_String, ___] := (Message[QECPauli::noprop, prop]; Missing["NotFound", prop])

QECPauli /: NonCommutativeMultiply[p_QECPauli, q_QECPauli] := p["Product", q]

QECPauli /: MakeBoxes[obj : QECPauli[v_List] /; pauliRowQ[v], form : (StandardForm | TraditionalForm)] :=
    InterpretationBox[RowBox[{"QECPauli", "[", #, "]"}], obj] & @ ToBoxes[pauliString[v], form]


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
