(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECCliffordHierarchyLevel]

PackageScope[cliffordHierarchyMemberQ]
PackageScope[proportionalToPauliQ]
PackageScope[hierarchyPaulis]


(* ============================================================================ *)
(* The Clifford hierarchy (Got26 Def 13.1).                                     *)
(*                                                                              *)
(*     C_1 = the Pauli group,                                                   *)
(*     C_k = { U : U P U^dagger in C_(k-1) for every Pauli P }.                 *)
(*                                                                              *)
(* C_2 is the Clifford group, and C_3 is what chapter 13 is about: the gates a  *)
(* magic state plus Cliffords can teleport (T = R_pi/8, Toffoli, CCZ, CS).      *)
(* The level of a gate is the least k with U in C_k, and it is what decides how *)
(* a gate can be made fault tolerant -- transversally (it had better be in C_2  *)
(* on most codes), by injection (C_3), or by a chain of injections.             *)
(*                                                                              *)
(* THE TEST IS THE DEFINITION, run recursively and exactly.  Membership in C_1  *)
(* is "proportional to a Pauli", read off the matrix: a Pauli X^x Z^z has its   *)
(* one nonzero per column at row i XOR x, so the first row locates x, and       *)
(* U X^x must then be diagonal with the sign pattern (-1)^(z.i) of some z.      *)
(* Global phases are ignored throughout, as they must be: e^(i phi) U and U     *)
(* conjugate identically.                                                       *)
(*                                                                              *)
(* WHICH PAULIS ARE CHECKED.  For k <= 3 it is enough to check the 2n           *)
(* generators X_j, Z_j, because the image of a product is the product of the    *)
(* images and C_1, C_2 are groups.  C_3 is not a group, so from k = 4 on every  *)
(* non-identity Pauli is checked, as the book warns after eq. 13.3.  That makes *)
(* the cost grow like 4^n per level, which is fine for the gates the hierarchy  *)
(* is asked about (one to three qubits) and is why the search stops at          *)
(* "MaxLevel" rather than running on: a unitary in no C_k (R_theta for theta    *)
(* not a dyadic multiple of pi) would never stop.                               *)
(*                                                                              *)
(* Exactness: exact entries are compared with PossibleZeroQ, so                  *)
(* T = diag(1, e^(i pi/4)) is decided exactly; machine numbers to Chop's         *)
(* tolerance.                                                                    *)
(*                                                                              *)
(* This is a property of a unitary, not of a code, and its natural home is      *)
(* QuantumOperator in the framework.  It lives here until the layer moves into  *)
(* the kernel, with the conversion from QuantumOperator done at the entrance.   *)
(*                                                                              *)
(* References: Got26 sec. 13.2.2 (Def 13.1, eqs. 13.2-13.11).                    *)
(* ============================================================================ *)

QECCliffordHierarchyLevel::usage = "QECCliffordHierarchyLevel[u] gives the least k such that the gate u is in the k-th level of the Clifford hierarchy: 1 for a Pauli, 2 for a Clifford gate, 3 for T, Toffoli, CCZ or controlled-S, and so on. u may be a QuantumOperator, a gate name such as \"T\", or a unitary matrix on qubits.\nIf u is in no level up to \"MaxLevel\" (default 4), the result is Missing[\"AboveLevel\", max].";

QECCliffordHierarchyLevel::gate = "`1` is not a unitary on qubits. Give a QuantumOperator, a gate name, or a 2^n x 2^n unitary matrix.";

Options[QECCliffordHierarchyLevel] = {"MaxLevel" -> 4};


(* ---- the input, as an exact unitary matrix on n qubits ---- *)

gateMatrix[op_Wolfram`QuantumFramework`QuantumOperator] := Normal[op["MatrixRepresentation"]]
gateMatrix[name_String] := Quiet @ Check[gateMatrix[Wolfram`QuantumFramework`QuantumOperator[name]], $Failed]
gateMatrix[name_[args___]] /; StringQ[name] :=
    Quiet @ Check[gateMatrix[Wolfram`QuantumFramework`QuantumOperator[name[args]]], $Failed]
gateMatrix[m_ ? SquareMatrixQ] := Normal[m]
gateMatrix[_] := $Failed

qubitUnitaryQ[m_] := MatrixQ[m] && With[{d = Length[m]},
    d >= 2 && IntegerQ[Log2[d]] && AllTrue[Flatten[m . ConjugateTranspose[m] - IdentityMatrix[d]], zeroQ]
]


(* ---- level 1: proportional to a Pauli ---- *)

(* Exact entries are decided exactly; machine-precision ones to Chop's tolerance, so a
   numerically built gate is not refused for rounding in its last digit. *)
zeroQ[x_] := PossibleZeroQ[Chop[x]]

proportionalToPauliQ[m_List] := Module[{d = Length[m], n, x, c, xm, diag, z, bits},
    n = Log2[d];
    c = SelectFirst[Range[d], ! zeroQ[m[[1, #]]] &, None];
    If[c === None, Return[False]];
    x = c - 1;
    (* U X^x, with X^x the permutation |i> -> |i XOR x> *)
    xm = m[[All, BitXor[Range[0, d - 1], x] + 1]];
    If[ ! AllTrue[Flatten[xm - DiagonalMatrix[Diagonal[xm]]], zeroQ], Return[False]];
    diag = Diagonal[xm] / xm[[1, 1]];
    (* qubit 1 is the most significant bit, so qubit j sits at 2^(n - j) *)
    z = Table[If[zeroQ[diag[[2^(n - j) + 1]] - 1], 0, 1], {j, n}];
    bits = IntegerDigits[Range[0, d - 1], 2, n];
    AllTrue[diag - (-1)^(bits . z), zeroQ]
]


(* ---- level k, by the definition ---- *)

hierarchyPaulis[n_Integer, k_Integer] := hierarchyPaulis[n, k] = Module[{rows},
    rows = If[ k <= 3,
        Join[
            Table[Append[ReplacePart[ConstantArray[0, 2 n], j -> 1], 0], {j, n}],
            Table[Append[ReplacePart[ConstantArray[0, 2 n], n + j -> 1], 0], {j, n}]
        ],
        Rest[Join[#, {0}] & /@ Tuples[{0, 1}, 2 n]]
    ];
    Normal[pauliRowMatrix[#]] & /@ rows
]

cliffordHierarchyMemberQ[m_List, 1] := proportionalToPauliQ[m]

cliffordHierarchyMemberQ[m_List, k_Integer /; k >= 2] := With[{md = ConjugateTranspose[m]},
    AllTrue[hierarchyPaulis[Log2[Length[m]], k], cliffordHierarchyMemberQ[m . # . md, k - 1] &]
]


(* ---- the function ---- *)

QECCliffordHierarchyLevel[u_, opts : OptionsPattern[]] := Module[
    {m = gateMatrix[u], max = OptionValue["MaxLevel"]},
    If[ m === $Failed || ! qubitUnitaryQ[m],
        Message[QECCliffordHierarchyLevel::gate, u]; Return[$Failed]
    ];
    Replace[
        SelectFirst[Range[max], cliffordHierarchyMemberQ[m, #] &, None],
        None -> Missing["AboveLevel", max]
    ]
]
