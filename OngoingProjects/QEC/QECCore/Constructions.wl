(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECConcatenate]
PackageExport[QECRemoveQubit]
PackageExport[QECPasteCodes]

PackageScope[pauliTensor]
PackageScope[embedPauli]


(* ============================================================================ *)
(* Constructions that build codes out of codes.                                 *)
(* Gottesman thesis sec. 3.5 (concatenation, qubit removal, pasting) and sec.    *)
(* 8.6 / QECC book ch. 5 (CSS).                                                 *)
(*                                                                              *)
(* All of these are row operations and block embeddings on Pauli rows, and here  *)
(* they carry phases through: a generator product is taken with QECPauliProduct  *)
(* instead of adding symplectic vectors and discarding the sign.  It matters     *)
(* because S and -S are different codes, and the prototype's string surgery      *)
(* silently dropped the distinction.                                            *)
(* ============================================================================ *)


(* ---- block algebra on rows ---- *)

(* Tensor product of Paulis on disjoint registers: the halves concatenate and the
   phases add, because the factors act on different qubits and so commute. *)
pauliTensor[u_List, v_List] := With[{m = pauliQubits[u], n = pauliQubits[v]},
    Join[u[[1 ;; m]], v[[1 ;; n]], u[[m + 1 ;; 2 m]], v[[n + 1 ;; 2 n]], {Mod[Last[u] + Last[v], 4]}]
]

pauliTensor[u_List, v_List, rest__List] := Fold[pauliTensor, u, {v, rest}]

(* v placed on qubits before+1 .. before+n of a register with `after` more qubits. *)
embedPauli[v_List, before_Integer, after_Integer] :=
    pauliTensor[pauliIdentity[before], v, pauliIdentity[after]]


(* ---- CSS ---- *)

QECCode::cssdims = "HX and HZ must have the same number of columns (qubits).";
QECCode::cssorth = "HX and HZ must satisfy HX.Transpose[HZ] = 0 (mod 2) so that X-type and Z-type stabilizers commute.";

QECCode["CSS", hx_ ? MatrixQ, hz_ ? MatrixQ] := Module[{n = Dimensions[hx][[2]]},
    If[Dimensions[hz][[2]] =!= n, Message[QECCode::cssdims]; Return[$Failed]];
    If[ ! AllTrue[Flatten[Mod[hx . Transpose[hz], 2]], # === 0 &],
        Message[QECCode::cssorth]; Return[$Failed]
    ];
    QECCode[Join[
        Join[#, gf2Zero[n], {0}] & /@ hx,
        Join[gf2Zero[n], #, {0}] & /@ hz
    ]]
]

QECCode["CSS", h_ ? MatrixQ] := QECCode["CSS", h, h]


(* ---- concatenation ---- *)

QECConcatenate::usage = "QECConcatenate[outer, inner] concatenates two codes, re-encoding each qubit of the outer code with the inner code. An [[n1,k,d1]] outer and an [[n2,1,d2]] inner give an [[n1 n2, k, d1 d2]] code.";

QECConcatenate::inner = "The inner code must encode exactly one logical qubit.";

QECConcatenate[QECCode[outer_Association], QECCode[inner_Association]] := Module[
    {n1, n2, logical, xbar, zbar, translate, innerGens, outerGens},

    If[codeLogicalQubits[inner] =!= 1, Message[QECConcatenate::inner]; Return[$Failed]];

    n1 = outer["Qubits"];
    n2 = inner["Qubits"];
    logical = codeLogicalVectors[inner];
    xbar = First[logical["X"]];
    zbar = First[logical["Z"]];

    (* One qubit of the outer code becomes one inner block: I is the identity on
       the block, X and Z become the inner logical operators, and Y = i X Z becomes
       i Xbar Zbar, which is Hermitian again because Xbar and Zbar anticommute. *)
    translate[{x_, z_}] := Which[
        {x, z} === {0, 0}, pauliIdentity[n2],
        {x, z} === {1, 0}, xbar,
        {x, z} === {0, 1}, zbar,
        True, MapAt[Mod[# + 1, 4] &, QECPauliProduct[xbar, zbar], -1]
    ];

    innerGens = Flatten[
        Table[embedPauli[g, (b - 1) n2, (n1 - b) n2], {b, n1}, {g, codeVectors[inner]}],
        1
    ];

    outerGens = Function[g,
        MapAt[
            Mod[# + Last[g], 4] &,
            Fold[pauliTensor, translate /@ Transpose[{g[[1 ;; n1]], g[[n1 + 1 ;; 2 n1]]}]],
            -1
        ]
    ] /@ codeVectors[outer];

    QECCode[Join[innerGens, outerGens]]
]


(* ---- qubit removal ---- *)

QECRemoveQubit::usage = "QECRemoveQubit[code] removes one qubit, turning an [[n,k,d]] code into an [[n-1,k+1,d-1]] code.\nQECRemoveQubit[code, q] removes the given qubit.";

QECRemoveQubit::nopivot = "No qubit admits the removal surgery: some qubit must carry both an X-type and a Z-type generator tail.";

(* Recombine the generators so that exactly one ends in X and one ends in Z at the
   removed qubit and all the others end in I; drop those two (they become the new
   logical pair) and truncate the rest.  Row operations are Pauli products, so the
   surviving generators keep their correct signs. *)
removeQubitAt[a_Association, q_Integer] := Module[{n = a["Qubits"], rows, xrow, zrow, keep, columns},
    rows = codeVectors[a];

    xrow = SelectFirst[Range[Length[rows]], rows[[#, q]] === 1 &];
    If[MissingQ[xrow], Return[$Failed]];
    rows = MapIndexed[
        If[First[#2] =!= xrow && #1[[q]] === 1, QECPauliProduct[#1, rows[[xrow]]], #1] &,
        rows
    ];

    zrow = SelectFirst[Range[Length[rows]], # =!= xrow && rows[[#, n + q]] === 1 &];
    If[MissingQ[zrow], Return[$Failed]];
    rows = MapIndexed[
        If[First[#2] =!= zrow && #1[[n + q]] === 1, QECPauliProduct[#1, rows[[zrow]]], #1] &,
        rows
    ];

    keep = Complement[Range[Length[rows]], {xrow, zrow}];
    columns = Join[Complement[Range[n], {q}], Complement[Range[n + 1, 2 n], {n + q}], {2 n + 1}];

    QECCode[rows[[keep, columns]]]
]

QECRemoveQubit[QECCode[a_Association], q_Integer] := With[{result = removeQubitAt[a, q]},
    If[result === $Failed, Message[QECRemoveQubit::nopivot]; $Failed, result]
]

QECRemoveQubit[QECCode[a_Association]] := Module[{q},
    q = SelectFirst[Reverse[Range[a["Qubits"]]], removeQubitAt[a, #] =!= $Failed &];
    If[MissingQ[q], Message[QECRemoveQubit::nopivot]; Return[$Failed]];
    removeQubitAt[a, q]
]


(* ---- pasting ---- *)

QECPasteCodes::usage = "QECPasteCodes[code1, code2, solo1, solo2] pastes two codes together, keeping the first solo1 generators of code1 and the first solo2 of code2 acting alone and pairing the rest across the two blocks.";

QECPasteCodes::pairing = "Both codes must leave the same number of generators to pair: `1` against `2`.";
QECPasteCodes::solo = "Solo counts must lie between 0 and the number of generators of each code.";

QECPasteCodes[QECCode[a1_Association], QECCode[a2_Association], solo1_Integer, solo2_Integer] := Module[
    {m1 = codeStabilizerCount[a1], m2 = codeStabilizerCount[a2], n1 = a1["Qubits"], n2 = a2["Qubits"], g1, g2, solo, paired},

    If[ ! (0 <= solo1 <= m1 && 0 <= solo2 <= m2), Message[QECPasteCodes::solo]; Return[$Failed]];
    If[ m1 - solo1 =!= m2 - solo2, Message[QECPasteCodes::pairing, m1 - solo1, m2 - solo2]; Return[$Failed]];

    g1 = codeVectors[a1];
    g2 = codeVectors[a2];

    solo = Join[
        pauliTensor[#, pauliIdentity[n2]] & /@ Take[g1, solo1],
        pauliTensor[pauliIdentity[n1], #] & /@ Take[g2, solo2]
    ];
    paired = MapThread[pauliTensor, {Drop[g1, solo1], Drop[g2, solo2]}];

    QECCode[Join[solo, paired]]
]
