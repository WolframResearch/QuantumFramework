(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Structure.wlt

   Standard form, logical operators, distance, and membership in the stabilizer
   group.  The structural results are checked for internal consistency (the
   logical pair must commute with every check and anticommute with each other)
   rather than only against stored answers, so the tests still bite if the
   algorithm is replaced.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

(* Transitional: the rebuilt QEC core still lives under OngoingProjects/QEC/.
   Once it moves into Kernel/QEC/ and PacletInfo.wl lists its context, this Get
   disappears and the Needs above is enough.  The repo root is found from the
   loaded paclet, which RunTests.wls points at this checkout through
   PacletDirectoryLoad, so the tests always run against the source tree. *)
Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

(* Derived properties are memoised on the code's data.  A kernel that already
   held results from an earlier load would let those outrank the definitions just
   read, so the tests would silently check the old package. *)
QECClearCache[];

qecNamed = {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"};

(* The defining property of a logical pair: it commutes with the checks, it is
   not in the stabilizer group, Xbar_i anticommutes with Zbar_i and with nothing
   else, and the X's and the Z's commute among themselves. *)
qecLogicalsConsistentQ[code_] := Module[{lo = code["LogicalOperators"], xs, zs, k},
    xs = lo["X"]; zs = lo["Z"]; k = Length[xs];
    AllTrue[Flatten[Outer[QECPauliCommuteQ, Join[xs, zs], code["Generators"]]], TrueQ] &&
    And @@ Flatten[Table[QECPauliCommuteQ[xs[[i]], zs[[j]]] === (i =!= j), {i, k}, {j, k}]] &&
    AllTrue[Flatten[Outer[QECPauliCommuteQ, xs, xs]], TrueQ] &&
    AllTrue[Flatten[Outer[QECPauliCommuteQ, zs, zs]], TrueQ] &&
    AllTrue[Join[xs, zs], ! code["StabilizerMemberQ", #] &]
];


(* ============================================================================
   Standard form
   ============================================================================ *)

VerificationTest[QECCode["BitFlipCode"]["StandardForm"]["XRank"], 0, TestID -> "QEC-Structure-xrank-bitflip"]

VerificationTest[QECCode["PhaseFlipCode"]["StandardForm"]["XRank"], 2, TestID -> "QEC-Structure-xrank-phaseflip"]

VerificationTest[
    MatrixRank[QECCode["5QubitCode"]["StandardForm"]["Matrix"], Modulus -> 2],
    4,
    TestID -> "QEC-Structure-standard-form-rank"
]

VerificationTest[
    Sort[QECCode["ShorCode"]["StandardForm"]["QubitPermutation"]],
    Range[9],
    TestID -> "QEC-Structure-permutation-is-a-permutation"
]

(* Row operations and column relabelling preserve the rank of the check matrix. *)
VerificationTest[
    AllTrue[qecNamed,
        With[{code = QECCode[#]},
            MatrixRank[code["StandardForm"]["Matrix"], Modulus -> 2] === code["StabilizerCount"]
        ] &
    ],
    True,
    TestID -> "QEC-Structure-standard-form-preserves-rank"
]


(* ============================================================================
   Logical operators
   ============================================================================ *)

VerificationTest[
    QECCode["BitFlipCode"]["LogicalOperators"],
    <|"X" -> {"XXX"}, "Z" -> {"IIZ"}|>,
    TestID -> "QEC-Structure-logicals-bitflip"
]

(* Gottesman's thesis gives Xbar = ZIIZX for the five-qubit code. *)
VerificationTest[QECCode["5QubitCode"]["LogicalX"], {"ZIIZX"}, TestID -> "QEC-Structure-logical-X-five-qubit"]

VerificationTest[QECCode["5QubitCode"]["LogicalZ"], {"ZZZZZ"}, TestID -> "QEC-Structure-logical-Z-five-qubit"]

VerificationTest[
    qecLogicalsConsistentQ[QECCode[#]] & /@ qecNamed,
    {True, True, True, True, True},
    TestID -> "QEC-Structure-logicals-consistent-named"
]

VerificationTest[
    qecLogicalsConsistentQ[QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]],
    True,
    TestID -> "QEC-Structure-logicals-consistent-multi-logical"
]

VerificationTest[
    Length[QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]["LogicalX"]],
    6,
    TestID -> "QEC-Structure-logical-count"
]

(* Zbar and the all-Z operator differ by a stabilizer on both these codes. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"]},
        code["StabilizerMemberQ", QECPauliProduct[First[code["LogicalZ"]], "ZZZZZ"]]
    ],
    True,
    TestID -> "QEC-Structure-logical-Z-coset-five-qubit"
]

VerificationTest[
    Module[{code = QECCode["SteaneCode"]},
        code["StabilizerMemberQ", QECPauliProduct[First[code["LogicalZ"]], "ZZZZZZZ"]]
    ],
    True,
    TestID -> "QEC-Structure-logical-Z-coset-steane"
]

VerificationTest[
    QECCode[QECCode["5QubitCode"]["CompletedGenerators"]]["LogicalOperators"],
    <|"X" -> {}, "Z" -> {}|>,
    TestID -> "QEC-Structure-no-logicals-when-k-zero"
]


(* ============================================================================
   Membership, with and without signs

   -M is in the normaliser but not in the stabilizer group: it maps the code
   space to itself with a sign, so it is not a stabilizer of it.  The prototype
   compared spans only and answered True for both M and -M.
   ============================================================================ *)

VerificationTest[QECCode["BitFlipCode"]["StabilizerMemberQ", "ZIZ"], True, TestID -> "QEC-Structure-member-product"]

VerificationTest[QECCode["BitFlipCode"]["StabilizerMemberQ", "ZII"], False, TestID -> "QEC-Structure-member-outside"]

VerificationTest[QECCode["BitFlipCode"]["StabilizerMemberQ", "III"], True, TestID -> "QEC-Structure-member-identity"]

VerificationTest[QECCode["BitFlipCode"]["StabilizerMemberQ", "-ZIZ"], False, TestID -> "QEC-Structure-member-wrong-sign"]

VerificationTest[QECCode["BitFlipCode"]["StabilizerMemberQ", "ZZI"], True, TestID -> "QEC-Structure-member-generator"]

(* Every product of generators is in the group, with the sign the product carries. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], gens},
        gens = code["GeneratorVectors"];
        AllTrue[Subsets[gens, {1, 3}], code["StabilizerMemberQ", QECPauliProduct @@ #] &]
    ],
    True,
    TestID -> "QEC-Structure-member-all-products"
]

VerificationTest[QECCode["BitFlipCode"]["LogicalPauliQ", "IIZ"], True, TestID -> "QEC-Structure-logical-pauli-yes"]

VerificationTest[QECCode["BitFlipCode"]["LogicalPauliQ", "ZZI"], False, TestID -> "QEC-Structure-logical-pauli-stabilizer"]

VerificationTest[QECCode["BitFlipCode"]["LogicalPauliQ", "IXI"], False, TestID -> "QEC-Structure-logical-pauli-detectable"]

(* A logical operator stays logical when multiplied by a stabilizer. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"], xbar},
        xbar = First[code["LogicalX"]];
        AllTrue[code["GeneratorVectors"], code["LogicalPauliQ", QECPauliProduct[xbar, #]] &]
    ],
    True,
    TestID -> "QEC-Structure-logical-coset-stays-logical"
]


(* ============================================================================
   Completion to a full stabilizer state
   ============================================================================ *)

VerificationTest[
    AllTrue[qecNamed,
        Module[{code = QECCode[#], full},
            full = QECCode[code["CompletedGenerators"]];
            full =!= $Failed && full["Qubits"] === code["Qubits"] && full["LogicalQubits"] === 0
        ] &
    ],
    True,
    TestID -> "QEC-Structure-completion-named"
]

VerificationTest[
    Module[{code = QECCode["SteaneCode"], full},
        full = code["CompletedGenerators"];
        AllTrue[Subsets[full, {2}], QECPauliCommuteQ @@ # &]
    ],
    True,
    TestID -> "QEC-Structure-completion-commutes"
]

(* The completion extends the generators: the originals come first, unchanged. *)
VerificationTest[
    Take[QECCode["5QubitCode"]["CompletedGenerators"], 4],
    QECCode["5QubitCode"]["GeneratorVectors"],
    TestID -> "QEC-Structure-completion-extends"
]


(* ============================================================================
   Distance
   ============================================================================ *)

VerificationTest[QECCode["BitFlipCode"]["Distance"], 1, TestID -> "QEC-Structure-distance-bitflip"]

VerificationTest[QECCode["PhaseFlipCode"]["Distance"], 1, TestID -> "QEC-Structure-distance-phaseflip"]

VerificationTest[QECCode["5QubitCode"]["Distance"], 3, TestID -> "QEC-Structure-distance-five"]

VerificationTest[QECCode["SteaneCode"]["Distance"], 3, TestID -> "QEC-Structure-distance-steane"]

VerificationTest[QECCode["ShorCode"]["Distance"], 3, TestID -> "QEC-Structure-distance-shor"]

VerificationTest[
    QECCode[{StringRepeat["X", 6], StringRepeat["Z", 6]}]["Distance"],
    2,
    TestID -> "QEC-Structure-distance-two-family"
]

VerificationTest[
    QECCode[QECCode["5QubitCode"]["CompletedGenerators"]]["Distance"],
    Infinity,
    TestID -> "QEC-Structure-distance-no-logical-qubits"
]

VerificationTest[
    QECPauliWeight[QECCode["5QubitCode"]["MinimumWeightLogical"]],
    3,
    TestID -> "QEC-Structure-witness-weight"
]

VerificationTest[
    QECCode["5QubitCode"]["LogicalPauliQ", QECCode["5QubitCode"]["MinimumWeightLogical"]],
    True,
    TestID -> "QEC-Structure-witness-is-logical"
]

(* The witness realises the distance, on every code we have. *)
VerificationTest[
    AllTrue[qecNamed,
        With[{code = QECCode[#]},
            QECPauliWeight[code["MinimumWeightLogical"]] === code["Distance"]
        ] &
    ],
    True,
    TestID -> "QEC-Structure-witness-realises-distance"
]

(* No logical operator is lighter than the distance: exhaustive below d. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"], d},
        d = code["Distance"];
        AllTrue[
            Flatten[Table[Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[5, w], {w, 1, d - 1}], 1],
            ! code["LogicalPauliQ", #] &
        ]
    ],
    True,
    TestID -> "QEC-Structure-nothing-below-the-distance"
]

VerificationTest[
    Module[{code = QECCode["SteaneCode"], d},
        d = code["Distance"];
        AllTrue[
            Flatten[Table[Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[7, w], {w, 1, d - 1}], 1],
            ! code["LogicalPauliQ", #] &
        ]
    ],
    True,
    TestID -> "QEC-Structure-nothing-below-the-distance-steane"
]


(* ============================================================================
   CSS
   ============================================================================ *)

VerificationTest[QECCode["SteaneCode"]["CSSQ"], True, TestID -> "QEC-Structure-css-steane"]

VerificationTest[QECCode["5QubitCode"]["CSSQ"], False, TestID -> "QEC-Structure-css-five"]

VerificationTest[QECCode["ShorCode"]["CSSQ"], True, TestID -> "QEC-Structure-css-shor"]

VerificationTest[QECCode[{"XZ", "ZX"}]["CSSQ"], False, TestID -> "QEC-Structure-css-mixed"]
