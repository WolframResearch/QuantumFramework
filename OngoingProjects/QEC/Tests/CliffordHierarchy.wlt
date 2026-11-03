(* ::Package:: *)

(* ============================================================================
   Tests/QEC/CliffordHierarchy.wlt

   The Clifford hierarchy (Got26 Def 13.1), decided by its definition.

   The values pinned are the book's own: Paulis at level 1, the Clifford group at
   2 (all 24 one-qubit Cliffords, and CNOT, CZ, SWAP), and at 3 the gates
   chapter 13 names -- R_pi/8 (= T, eqs. 13.2-13.3), Toffoli (eqs. 13.4-13.9),
   CCZ, controlled-S.  Then the tower above them: R_(m pi / 2^k) with m odd sits
   at exactly level k (the induction after eq. 13.11, in the book's convention
   R_theta ~ diag(1, e^(2 i theta))), so sqrt(T) is at 4 and an angle that is not
   a dyadic multiple of pi is in no level at all.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

ccz = DiagonalMatrix[{1, 1, 1, 1, 1, 1, 1, -1}];
rotation[theta_] := DiagonalMatrix[{1, Exp[I theta]}];


VerificationTest[
    QECCliffordHierarchyLevel /@ {"I", "X", "Y", "Z"},
    {1, 1, 1, 1},
    TestID -> "QEC-CliffordHierarchy-Paulis-are-level-one"
]

VerificationTest[
    QECCliffordHierarchyLevel /@ {"H", "S", "CNOT", "CZ", "SWAP"},
    {2, 2, 2, 2, 2},
    TestID -> "QEC-CliffordHierarchy-Clifford-gates-are-level-two"
]

(* Every word in H and S up to length 6 -- which reaches all 24 one-qubit Cliffords up
   to phase -- is at level 1 or 2, and the number of distinct gates among them, counted
   up to a global phase, is 24. *)
VerificationTest[
    Module[{h = {{1, 1}, {1, -1}} / Sqrt[2], s = {{1, 0}, {0, I}}, words, normal},
        words = Flatten[Table[Fold[#2 . #1 &, IdentityMatrix[2], w], {len, 0, 6}, {w, Tuples[{h, s}, len]}], 1];
        normal[m_] := Simplify[m / SelectFirst[Flatten[m], # =!= 0 &]];
        {Union[QECCliffordHierarchyLevel /@ words], Length[DeleteDuplicates[normal /@ words]]}
    ],
    {{1, 2}, 24},
    TestID -> "QEC-CliffordHierarchy-the-one-qubit-Clifford-group"
]

VerificationTest[
    {QECCliffordHierarchyLevel["T"], QECCliffordHierarchyLevel["Toffoli"], QECCliffordHierarchyLevel[ccz],
     QECCliffordHierarchyLevel["CS"]},
    {3, 3, 3, 3},
    TestID -> "QEC-CliffordHierarchy-the-chapter-13-gates-are-level-three"
]

(* diag(1, e^(i pi / 2^k)) is the book's R_(pi / 2^(k+1)) up to phase, so it sits at
   level k + 1: Z at 1, S at 2, T at 3, sqrt(T) at 4, the next at 5. *)
VerificationTest[
    Table[QECCliffordHierarchyLevel[rotation[Pi / 2^k], "MaxLevel" -> 6], {k, 0, 4}],
    {1, 2, 3, 4, 5},
    TestID -> "QEC-CliffordHierarchy-dyadic-rotations-climb-one-level-each"
]

VerificationTest[
    QECCliffordHierarchyLevel[rotation[1/10]],
    Missing["AboveLevel", 4],
    TestID -> "QEC-CliffordHierarchy-a-generic-angle-is-in-no-level"
]

(* A Clifford times a level-3 gate stays at 3; the global phase never matters. *)
VerificationTest[
    {QECCliffordHierarchyLevel[QuantumOperator["H"] @ QuantumOperator["T"]],
     QECCliffordHierarchyLevel[Exp[I Pi / 7] Normal[QuantumOperator["T"]["MatrixRepresentation"]]],
     QECCliffordHierarchyLevel[Exp[0.3 I] IdentityMatrix[2]]},
    {3, 3, 1},
    TestID -> "QEC-CliffordHierarchy-Cliffords-and-phases-do-not-move-the-level"
]

VerificationTest[
    QECCliffordHierarchyLevel[{{1, 1}, {0, 1}}],
    $Failed,
    {QECCliffordHierarchyLevel::gate},
    TestID -> "QEC-CliffordHierarchy-refuses-a-non-unitary"
]
