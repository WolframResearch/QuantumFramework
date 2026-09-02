(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Pauli.wlt

   The Pauli layer of the QEC core: rows {x1..xn, z1..zn, e}, their string form,
   weight, commutation and multiplication.

   Two independent oracles are used rather than the package's own conventions.
   Dense 2^n x 2^n matrices, built here from scratch, check the algebra: a
   product is right when the matrices agree.  The engine's own PauliRow (in
   Kernel/Stabilizer/Conversions.m) checks the layout: our rows must be the rows
   the stabilizer engine already speaks, or the two halves of the framework would
   drift apart.
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

(* ---- local oracle: Pauli string -> dense matrix ---- *)

qecPauliMat["I"] = {{1, 0}, {0, 1}};
qecPauliMat["X"] = {{0, 1}, {1, 0}};
qecPauliMat["Y"] = {{0, -I}, {I, 0}};
qecPauliMat["Z"] = {{1, 0}, {0, -1}};

qecPrefixFactor["-"] = -1;
qecPrefixFactor["i"] = I;
qecPrefixFactor["-i"] = -I;
qecPrefixFactor[""] = 1;

qecDense[s_String] := Module[{prefix, body},
    prefix = Replace[StringCases[s, StartOfString ~~ p : ("-i" | "-" | "i") :> p], {{p_} :> p, _ -> ""}];
    body = StringDrop[s, StringLength[prefix]];
    qecPrefixFactor[prefix] * Fold[KroneckerProduct, qecPauliMat /@ Characters[body]]
];

qecRandomPaulis[n_, count_] := Table[StringJoin[RandomChoice[{"I", "X", "Y", "Z"}, n]], {count}];


(* ============================================================================
   Rows in, rows out
   ============================================================================ *)

VerificationTest[QECPauliVector["XI"], {1, 0, 0, 0, 0}, TestID -> "QEC-Pauli-vector-XI"]

VerificationTest[QECPauliVector["IX"], {0, 1, 0, 0, 0}, TestID -> "QEC-Pauli-vector-IX"]

VerificationTest[QECPauliVector["ZI"], {0, 0, 1, 0, 0}, TestID -> "QEC-Pauli-vector-ZI"]

VerificationTest[QECPauliVector["IZ"], {0, 0, 0, 1, 0}, TestID -> "QEC-Pauli-vector-IZ"]

VerificationTest[QECPauliVector["YZ"], {1, 0, 1, 1, 0}, TestID -> "QEC-Pauli-vector-YZ"]

VerificationTest[QECPauliVector["XZZXI"], {1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0}, TestID -> "QEC-Pauli-vector-XZZXI"]

VerificationTest[QECPauliVector["-ZI"], {0, 0, 1, 0, 2}, TestID -> "QEC-Pauli-vector-minus"]

VerificationTest[QECPauliVector["iZI"], {0, 0, 1, 0, 1}, TestID -> "QEC-Pauli-vector-i"]

VerificationTest[QECPauliVector["-iZI"], {0, 0, 1, 0, 3}, TestID -> "QEC-Pauli-vector-minus-i"]

VerificationTest[Quiet[QECPauliVector["XQZ"]], $Failed, TestID -> "QEC-Pauli-vector-bad-letter"]

VerificationTest[Quiet[QECPauliVector[""]], $Failed, TestID -> "QEC-Pauli-vector-empty"]

VerificationTest[QECPauliString[{1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0}], "XZZXI", TestID -> "QEC-Pauli-string-XZZXI"]

VerificationTest[QECPauliString[{0, 0, 1, 1, 0}], "ZZ", TestID -> "QEC-Pauli-string-ZZ"]

VerificationTest[QECPauliString[{0, 0, 1, 1, 2}], "-ZZ", TestID -> "QEC-Pauli-string-minus-ZZ"]

VerificationTest[QECPauliString[{0, 0, 1, 1, 3}], "-iZZ", TestID -> "QEC-Pauli-string-minus-i-ZZ"]

VerificationTest[
    Module[{},
        SeedRandom[1];
        AllTrue[
            Table[
                With[{s = StringJoin[RandomChoice[{"I", "X", "Y", "Z"}, RandomInteger[{1, 12}]]]},
                    QECPauliString[QECPauliVector[s]] === s
                ],
                {100}
            ],
            TrueQ
        ]
    ],
    True,
    TestID -> "QEC-Pauli-roundtrip-random"
]

VerificationTest[
    Module[{},
        SeedRandom[2];
        AllTrue[
            Table[
                With[{s = RandomChoice[{"", "-", "i", "-i"}] <> StringJoin[RandomChoice[{"I", "X", "Y", "Z"}, 6]]},
                    QECPauliString[QECPauliVector[s]] === s
                ],
                {60}
            ],
            TrueQ
        ]
    ],
    True,
    TestID -> "QEC-Pauli-roundtrip-phases"
]

VerificationTest[QECPauliQ["XYZ"], True, TestID -> "QEC-Pauli-Q-string"]

VerificationTest[QECPauliQ["XQZ"], False, TestID -> "QEC-Pauli-Q-bad"]

VerificationTest[QECPauliQ[{1, 0, 1, 0, 0}], True, TestID -> "QEC-Pauli-Q-row"]

VerificationTest[QECPauliQ[{1, 0, 1, 0}], False, TestID -> "QEC-Pauli-Q-even-row"]

VerificationTest[QECPauliQ[3], False, TestID -> "QEC-Pauli-Q-nonsense"]

VerificationTest[QECPauliPhase["-XZ"], 2, TestID -> "QEC-Pauli-phase-minus"]

VerificationTest[QECPauliPhase["XZ"], 0, TestID -> "QEC-Pauli-phase-plus"]


(* ============================================================================
   Layout agreement with the stabilizer engine

   PauliRow returns Join[xbits, zbits, {phase}] with a single phase bit.  Our
   rows must agree on the symplectic part exactly, and the phase must be that bit
   doubled, since we carry Z4 where the engine carries Z2.
   ============================================================================ *)

VerificationTest[
    Module[{samples},
        SeedRandom[3];
        samples = Join[{"XZ", "ZX", "YI", "IY", "XY", "ZZ", "XXI"}, qecRandomPaulis[4, 20]];
        AllTrue[samples, Most[PauliRow[qecDense[#], StringLength[#]]] === Most[QECPauliVector[#]] &]
    ],
    True,
    TestID -> "QEC-Pauli-layout-matches-engine"
]

VerificationTest[
    Module[{samples},
        SeedRandom[4];
        samples = ("-" <> # &) /@ qecRandomPaulis[3, 15];
        AllTrue[samples,
            With[{row = PauliRow[qecDense[#], 3], ours = QECPauliVector[#]},
                Most[row] === Most[ours] && 2 Last[row] === Last[ours]
            ] &
        ]
    ],
    True,
    TestID -> "QEC-Pauli-phase-matches-engine"
]

VerificationTest[
    With[{p = PauliStabilizer[3]["H", 1]["CNOT", 1, 2]["CNOT", 2, 3]},
        Most[QECPauliVector[#]] & /@ Join[p["Destabilizers"], p["Stabilizers"]] === Normal[p["Matrix"]]
    ],
    True,
    TestID -> "QEC-Pauli-tableau-agreement"
]


(* ============================================================================
   Weight
   ============================================================================ *)

VerificationTest[QECPauliWeight["XZZXI"], 4, TestID -> "QEC-Pauli-weight-XZZXI"]

VerificationTest[QECPauliWeight["III"], 0, TestID -> "QEC-Pauli-weight-identity"]

VerificationTest[QECPauliWeight["-ZI"], 1, TestID -> "QEC-Pauli-weight-signed"]

VerificationTest[QECPauliWeight[{1, 0, 1, 1, 0}], 2, TestID -> "QEC-Pauli-weight-row"]

VerificationTest[
    Module[{},
        SeedRandom[5];
        AllTrue[qecRandomPaulis[8, 40], QECPauliWeight[#] === StringLength[#] - StringCount[#, "I"] &]
    ],
    True,
    TestID -> "QEC-Pauli-weight-random"
]


(* ============================================================================
   Commutation, against the dense oracle
   ============================================================================ *)

VerificationTest[QECPauliCommuteQ["X", "Z"], False, TestID -> "QEC-Pauli-commute-XZ"]

VerificationTest[QECPauliCommuteQ["X", "Y"], False, TestID -> "QEC-Pauli-commute-XY"]

VerificationTest[QECPauliCommuteQ["Y", "Z"], False, TestID -> "QEC-Pauli-commute-YZ"]

VerificationTest[QECPauliCommuteQ["X", "X"], True, TestID -> "QEC-Pauli-commute-XX"]

VerificationTest[QECPauliCommuteQ["X", "I"], True, TestID -> "QEC-Pauli-commute-XI"]

VerificationTest[QECPauliCommuteQ["XX", "ZZ"], True, TestID -> "QEC-Pauli-commute-XXZZ"]

VerificationTest[QECPauliCommuteQ["XI", "ZZ"], False, TestID -> "QEC-Pauli-commute-XIZZ"]

VerificationTest[QECPauliCommuteQ["XXXX", "ZZZZ"], True, TestID -> "QEC-Pauli-commute-even"]

VerificationTest[QECPauliCommuteQ["XXXXX", "ZZZZZ"], False, TestID -> "QEC-Pauli-commute-odd"]

VerificationTest[QECPauliCommuteQ["-ZI", "XI"], False, TestID -> "QEC-Pauli-commute-signed"]

VerificationTest[
    Outer[QECPauliCommuteQ, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}],
    ConstantArray[True, {4, 4}],
    TestID -> "QEC-Pauli-commute-five-qubit-checks"
]

VerificationTest[
    QECPauliCommuteQ["IXIII", #] & /@ {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"},
    {False, True, True, True},
    TestID -> "QEC-Pauli-commute-five-qubit-error"
]

VerificationTest[Quiet[QECPauliCommuteQ["XI", "XXX"]], $Failed, TestID -> "QEC-Pauli-commute-size-mismatch"]

VerificationTest[
    Module[{samples},
        SeedRandom[6];
        samples = qecRandomPaulis[3, 25];
        AllTrue[Tuples[samples, 2],
            QECPauliCommuteQ @@ # === (qecDense[#[[1]]] . qecDense[#[[2]]] == qecDense[#[[2]]] . qecDense[#[[1]]]) &
        ]
    ],
    True,
    TestID -> "QEC-Pauli-commute-vs-dense"
]


(* ============================================================================
   Multiplication, against the dense oracle

   This is the part the prototype did not have at all: it added symplectic
   vectors and threw the phase away, so X.Z came back as Y instead of -iY.
   ============================================================================ *)

VerificationTest[QECPauliString[QECPauliProduct["X", "Z"]], "-iY", TestID -> "QEC-Pauli-product-XZ"]

VerificationTest[QECPauliString[QECPauliProduct["Z", "X"]], "iY", TestID -> "QEC-Pauli-product-ZX"]

VerificationTest[QECPauliString[QECPauliProduct["X", "X"]], "I", TestID -> "QEC-Pauli-product-XX"]

VerificationTest[QECPauliString[QECPauliProduct["Y", "Y"]], "I", TestID -> "QEC-Pauli-product-YY"]

VerificationTest[QECPauliString[QECPauliProduct["Y", "X"]], "-iZ", TestID -> "QEC-Pauli-product-YX"]

(* The prototype answered "XZY" here: it added symplectic vectors and dropped the
   phase.  X.Z on the third qubit contributes -i, and the dense oracle above agrees. *)
VerificationTest[QECPauliString[QECPauliProduct["XIX", "IZZ"]], "-iXZY", TestID -> "QEC-Pauli-product-strings"]

VerificationTest[QECPauliString[QECPauliProduct["XX", "ZZ"]], "-YY", TestID -> "QEC-Pauli-product-two-qubit"]

VerificationTest[QECPauliString[QECPauliProduct["X", "Y", "Z"]], "iI", TestID -> "QEC-Pauli-product-three-factors"]

VerificationTest[Quiet[QECPauliProduct["XI", "XXX"]], $Failed, TestID -> "QEC-Pauli-product-size-mismatch"]

VerificationTest[
    Module[{samples},
        SeedRandom[7];
        samples = qecRandomPaulis[3, 25];
        AllTrue[Tuples[samples, 2],
            qecDense[QECPauliString[QECPauliProduct @@ #]] == qecDense[#[[1]]] . qecDense[#[[2]]] &
        ]
    ],
    True,
    TestID -> "QEC-Pauli-product-vs-dense"
]

VerificationTest[
    Module[{samples},
        SeedRandom[8];
        samples = (RandomChoice[{"", "-", "i", "-i"}] <> # &) /@ qecRandomPaulis[2, 20];
        AllTrue[Tuples[samples, 2],
            qecDense[QECPauliString[QECPauliProduct @@ #]] == qecDense[#[[1]]] . qecDense[#[[2]]] &
        ]
    ],
    True,
    TestID -> "QEC-Pauli-product-phases-vs-dense"
]

VerificationTest[
    Module[{samples},
        SeedRandom[9];
        samples = qecRandomPaulis[4, 20];
        AllTrue[samples, QECPauliString[QECPauliProduct[#, #]] === StringRepeat["I", 4] &]
    ],
    True,
    TestID -> "QEC-Pauli-product-self-inverse"
]


(* ============================================================================
   Enumeration of low-weight errors
   ============================================================================ *)

VerificationTest[
    QECPauliString /@ Wolfram`QuantumFramework`QEC`PackageScope`weightOneVectors[2],
    {"XI", "YI", "ZI", "IX", "IY", "IZ"},
    TestID -> "QEC-Pauli-weight-one-enumeration"
]

VerificationTest[
    Sort[QECPauliString /@ Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[2, 1]],
    Sort[{"XI", "YI", "ZI", "IX", "IY", "IZ"}],
    TestID -> "QEC-Pauli-weight-k-one"
]

VerificationTest[
    Length[Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[5, 2]],
    9 Binomial[5, 2],
    TestID -> "QEC-Pauli-weight-k-count"
]

VerificationTest[
    QECPauliString /@ Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[3, 0],
    {"III"},
    TestID -> "QEC-Pauli-weight-zero"
]

VerificationTest[
    AllTrue[Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[4, 3], QECPauliWeight[#] === 3 &],
    True,
    TestID -> "QEC-Pauli-weight-k-all-have-weight-k"
]

VerificationTest[
    Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[3, 4],
    {},
    TestID -> "QEC-Pauli-weight-k-too-large"
]
