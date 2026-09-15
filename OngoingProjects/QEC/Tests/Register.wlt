(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Register.wlt

   Several blocks of one code, and the one gate between them.

   Everything before this file lives on a single block, which is enough for a
   memory experiment and not enough for anything that computes: a logical qubit
   is a BLOCK, so a logical two-qubit gate is a gate between two of them, and
   Got26 Theorem 13.2 needs 2m blocks for a gate touching m.

   The load-bearing test is QEC-Register-transversal-CNOT-is-the-logical-CNOT.
   Steane EC already rested on that claim (sec. 12.3.1) and this is where it is
   checked rather than cited, by conjugating the logical operators through the
   gate and reading the images off.  Its companion is the negative one: on a
   non-CSS code the same gate takes eight of the eight register generators out of
   the stabilizer group, so it is not a logical operation at all.  Both halves
   are needed -- preserving the group and acting as the CNOT are different
   claims, and a gate can do the second on the Z operators while failing the
   first.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecConjugate = Symbol[qecScope <> "registerConjugate"];
qecLabelM    = Symbol[qecScope <> "codeLabelMatrix"];
qecCodeData  = Symbol[qecScope <> "codeData"];

steane = QECCode["SteaneCode"];
five = QECCode["5QubitCode"];

reg = QECRegister[steane, 2];
reg5 = QECRegister[five, 2];

(* The product of two Pauli strings, as a string. *)
qecProd[s1_, s2_] := QECPauliString[QECPauliProduct[s1, s2]]

(* Where a circuit sends each register generator. *)
qecImages[r_] := Module[{instr = r["TransversalCNOT", 1, 2], nq = r["Qubits"]},
    QECPauliString[qecConjugate[instr, nq, QECPauliVector[#]]] & /@ r["Generators"]
]


(* ============================================================================
   The layout
   ============================================================================ *)

(* Block b owns qubits (b-1)n+1 .. bn, stated once so nothing has to guess it. *)
VerificationTest[
    {reg["BlockRange", 1], reg["BlockRange", 2], reg["Index", 1, 3], reg["Index", 2, 3]},
    {Range[7], Range[8, 14], 3, 10},
    TestID -> "QEC-Register-blocks-are-laid-out-in-order"
]

VerificationTest[
    AssociationMap[reg, {"Blocks", "BlockQubits", "Qubits", "StabilizerCount", "LogicalQubits"}],
    <|"Blocks" -> 2, "BlockQubits" -> 7, "Qubits" -> 14, "StabilizerCount" -> 12,
      "LogicalQubits" -> 2|>,
    TestID -> "QEC-Register-counts-scale-with-the-blocks"
]

(* A single-block instruction list moves into a block; ops and arity untouched. *)
VerificationTest[
    reg["Lift", {{"R", 1}, {"H", 2}, {"CNOT", 1, 2}, {"CZ", 3, 7}}, 2],
    {{"R", 8}, {"H", 9}, {"CNOT", 8, 9}, {"CZ", 10, 14}},
    TestID -> "QEC-Register-lifts-a-circuit-into-a-block"
]

(* One block has to be the identity on everything, which is what keeps the rest of
   the package unaffected by this file existing. *)
VerificationTest[
    With[{one = QECRegister[steane, 1]},
        {one["Qubits"] === steane["Qubits"],
         one["Lift", {{"H", 1}, {"CNOT", 1, 2}}, 1] === {{"H", 1}, {"CNOT", 1, 2}},
         one["LabelMatrix"] === qecLabelM[qecCodeData[steane]],
         one["Generators"] === steane["Generators"]}
    ],
    {True, True, True, True},
    TestID -> "QEC-Register-of-one-block-is-the-code-itself"
]

(* The label matrix is block diagonal, and its two halves are exchanged, so each
   block's rows land in two windows rather than one contiguous run.  The shape is
   the cheap check that the scatter happened at all. *)
VerificationTest[
    Dimensions[reg["LabelMatrix"]],
    {16, 28},
    TestID -> "QEC-Register-label-matrix-is-block-diagonal"
]

(* Each block's logical operators, as register-wide Paulis with identity elsewhere. *)
VerificationTest[
    Map[QECPauliString, reg["LogicalVectors"], {2}],
    <|"X" -> {"IIXIXXIIIIIIII", "IIIIIIIIIXIXXI"},
      "Z" -> {"IZIZIZIIIIIIII", "IIIIIIIIZIZIZI"}|>,
    TestID -> "QEC-Register-logical-operators-sit-inside-their-block"
]


(* ============================================================================
   The one gate between blocks
   ============================================================================ *)

(* Transversal by construction: qubit q of one block touches qubit q of the other
   and nothing else. *)
VerificationTest[
    reg["TransversalCNOT", 1, 2],
    Table[{"CNOT", q, q + 7}, {q, 7}],
    TestID -> "QEC-Register-transversal-CNOT-is-qubit-by-qubit"
]

(* THE test.  Conjugating the logical operators through the gate gives the action
   of a CNOT on the logical pair, which is the claim Steane EC rests on. *)
VerificationTest[
    reg["LogicalAction", 1, 2],
    With[{
        x = QECPauliString /@ reg["LogicalVectors"]["X"],
        z = QECPauliString /@ reg["LogicalVectors"]["Z"]
    },
        <|{"X", 1, 1} -> qecProd[x[[1]], x[[2]]], {"X", 2, 1} -> x[[2]],
          {"Z", 1, 1} -> z[[1]],                  {"Z", 2, 1} -> qecProd[z[[1]], z[[2]]]|>
    ],
    TestID -> "QEC-Register-transversal-CNOT-is-the-logical-CNOT"
]

(* Acting correctly on the logical operators is only half of it: the gate also has
   to map the stabilizer group to itself, or it takes the state out of the code
   space and the logical action is beside the point. *)
VerificationTest[
    With[{group = QECCode[reg["Generators"]]},
        Count[qecImages[reg], g_ /; ! group["StabilizerMemberQ", g]]
    ],
    0,
    TestID -> "QEC-Register-transversal-CNOT-preserves-the-code-space"
]

(* And the negative.  On a non-CSS code the same gate is not a logical operation:
   every one of the eight register generators leaves the group.  This is the same
   restriction QECErrorCorrection enforces, seen at its source rather than at the
   gadget that inherits it. *)
VerificationTest[
    With[{group5 = QECCode[reg5["Generators"]]},
        {five["CSSQ"], Count[qecImages[reg5], g_ /; ! group5["StabilizerMemberQ", g]],
         Length[reg5["Generators"]]}
    ],
    {False, 8, 8},
    TestID -> "QEC-Register-a-non-CSS-code-has-no-transversal-CNOT"
]

(* The Z half of the action still comes out right on a non-CSS code, which is why
   the group test is the one that decides.  Checking the logical action alone
   would have passed half the time and been wrong. *)
VerificationTest[
    With[{
        act = reg5["LogicalAction", 1, 2],
        x = QECPauliString /@ reg5["LogicalVectors"]["X"],
        z = QECPauliString /@ reg5["LogicalVectors"]["Z"]
    },
        {act[{"Z", 2, 1}] === qecProd[z[[1]], z[[2]]],
         act[{"X", 1, 1}] === qecProd[x[[1]], x[[2]]]}
    ],
    {True, False},
    TestID -> "QEC-Register-the-Z-half-alone-would-have-fooled-us"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

VerificationTest[
    QECRegister[steane, 0],
    $Failed,
    {QECRegister::blocks},
    TestID -> "QEC-Register-needs-a-positive-block-count"
]

VerificationTest[
    reg["TransversalCNOT", 1, 3],
    $Failed,
    {QECRegister::block},
    TestID -> "QEC-Register-refuses-a-block-it-does-not-have"
]

VerificationTest[
    reg["TransversalCNOT", 2, 2],
    $Failed,
    {QECRegister::same},
    TestID -> "QEC-Register-a-CNOT-needs-two-different-blocks"
]

VerificationTest[
    reg["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECRegister::noprop},
    TestID -> "QEC-Register-unknown-property-is-refused"
]
