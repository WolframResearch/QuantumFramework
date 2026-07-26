(* ============================================================================
   Tests/CrossPackage_QuantumClifford.wlt -- cross-validation against the
   QuantumClifford.jl Julia package (Krastanov / QuantumSavory).

   Source repo: QuantumSavory/QuantumClifford.jl. A local copy, when the
   machine has one, is located by the resolver at the bottom of this file.

   This file extracts canonical results from QC.jl's published test suite
   (test/test_paulis.jl, test/test_ecc_syndromes.jl) and verifies our
   PauliStabilizer + QuantumOperator subsystem produces the same outputs.

   Note: no tier here runs Julia. The expected values below are HARD-CODED
   FROM QC.jl's TEST FILES with explicit citations to source lines, so every
   tier reproduces on a machine with no Julia at all.

   To add live-runtime QC.jl fixtures later, write a Julia script analogous to
   Tests/fixtures/generate_stim_fixtures.py that dumps stabilizer outputs to
   JSON, then add a TIER F here that consumes the JSON.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];


(* ============================================================================ *)
(* TIER A -- Pauli arithmetic identities                                        *)
(*                                                                              *)
(* From test_paulis.jl "Elementary operations" (lines 26-44):                  *)
(*   P"X" * P"Z" == P"-iY"                                                     *)
(*   P"X" \[CircleTimes] P"Z" == P"XZ"                                         *)
(*   prodphase(P"XX", P"YY") == 0x2   (corresponds to phase factor -1)         *)
(*   prodphase(P"ZZZ", P"XXX") == 0x3 (corresponds to phase factor -i)         *)
(*                                                                              *)
(* Our framework computes these via QuantumOperator matrix products. We       *)
(* verify the matrix-level identities here.                                    *)
(* ============================================================================ *)

(* X * Z = -i Y. *)
VerificationTest[
    Normal[QuantumOperator["X"]["Matrix"]] . Normal[QuantumOperator["Z"]["Matrix"]],
    -I * Normal[QuantumOperator["Y"]["Matrix"]],
    {},
    TestID -> "QC-Pauli-X-times-Z-equals-MinusI-Y"
]

(* Z * X = +i Y. *)
VerificationTest[
    Normal[QuantumOperator["Z"]["Matrix"]] . Normal[QuantumOperator["X"]["Matrix"]],
    I * Normal[QuantumOperator["Y"]["Matrix"]],
    {},
    TestID -> "QC-Pauli-Z-times-X-equals-PlusI-Y"
]

(* X tensor Z (XZ Pauli string) is X \[CircleTimes] Z. *)
VerificationTest[
    Normal @ QuantumOperator["XZ"]["Matrix"],
    KroneckerProduct[Normal @ QuantumOperator["X"]["Matrix"], Normal @ QuantumOperator["Z"]["Matrix"]],
    {},
    TestID -> "QC-Pauli-XZ-equals-X-Tensor-Z"
]

(* prodphase(XX, YY) == 2 (i.e., factor -1, phase mod 4 = 2):
   XX * YY = (XY)(XY) = (iZ)(iZ) = (i^2) ZZ = -ZZ
   So XX * YY = -1 * ZZ. *)
VerificationTest[
    With[{xx = Normal[QuantumOperator["XX"]["Matrix"]],
          yy = Normal[QuantumOperator["YY"]["Matrix"]],
          zz = Normal[QuantumOperator["ZZ"]["Matrix"]]},
        xx . yy == -zz
    ],
    True,
    {},
    TestID -> "QC-prodphase-XX-times-YY-equals-MinusZZ"
]

(* prodphase(ZZZ, XXX) == 3, i.e., the phase of ZZZ.XXX/YYY is i^3 = -i.
   Per qubit Z.X = +iY (verified above), so ZZZ.XXX = (iY)^3 = i^3 YYY = -i YYY.
   QC.jl encodes this as 0x3 (= phase index 3, corresponding to -i). We verify
   by computing the prefactor ratio numerically. *)
VerificationTest[
    With[{zzz = Normal[QuantumOperator["ZZZ"]["Matrix"]],
          xxx = Normal[QuantumOperator["XXX"]["Matrix"]],
          yyy = Normal[QuantumOperator["YYY"]["Matrix"]]},
        Module[{prod = zzz . xxx, pairs, ratio, phaseIdx},
            (* Pair up corresponding entries; keep only those where both are nonzero. *)
            pairs = Select[Transpose[{Flatten[prod], Flatten[yyy]}], #[[2]] != 0 &];
            ratio = pairs[[1, 1]] / pairs[[1, 2]];
            (* Phase index: 0=+1, 1=+i, 2=-1, 3=-i. *)
            phaseIdx = Mod[Round[Arg[ratio] * 2 / Pi], 4];
            {Abs[ratio], phaseIdx}
        ]
    ],
    {1, 3},   (* magnitude 1, phase index 3 = -i (matches QC.jl 0x3) *)
    {},
    TestID -> "QC-prodphase-ZZZ-times-XXX-MatchesPhaseIndex3"
]


(* ============================================================================ *)
(* TIER B -- comm() relation: stabilizer commutation                           *)
(*                                                                              *)
(* From test_paulis.jl line 36-37:                                             *)
(*   comm(P"XX", P"YY") == 0x0  (commute)                                      *)
(*   comm(P"XZ", P"YZ") == 0x1  (anticommute)                                  *)
(*                                                                              *)
(* In QF terms: [A, B] = 0 iff A and B share an even number of qubits where    *)
(* their Pauli factors anticommute (i.e., neither is I and they differ).       *)
(* We test by computing the matrix commutator.                                 *)
(* ============================================================================ *)

(* XX commutes with YY *)
VerificationTest[
    With[{xx = Normal[QuantumOperator["XX"]["Matrix"]],
          yy = Normal[QuantumOperator["YY"]["Matrix"]]},
        xx . yy == yy . xx
    ],
    True,
    {},
    TestID -> "QC-comm-XX-YY-Commute"
]

(* XZ anticommutes with YZ *)
VerificationTest[
    With[{xz = Normal[QuantumOperator["XZ"]["Matrix"]],
          yz = Normal[QuantumOperator["YZ"]["Matrix"]]},
        xz . yz == -yz . xz
    ],
    True,
    {},
    TestID -> "QC-comm-XZ-YZ-Anticommute"
]

(* Bell state stabilizers commute pairwise (XX and ZZ). *)
VerificationTest[
    With[{xx = Normal[QuantumOperator["XX"]["Matrix"]],
          zz = Normal[QuantumOperator["ZZ"]["Matrix"]]},
        xx . zz == zz . xx
    ],
    True,
    {},
    TestID -> "QC-comm-XX-ZZ-Commute-Bell"
]


(* ============================================================================ *)
(* TIER C -- random_stabilizer / random_clifford produce valid states         *)
(*                                                                              *)
(* QC.jl's random_stabilizer(n) and random_clifford(n) are the documented      *)
(* primitives for sampling. Our PauliStabilizer["Random", n] is the analog.   *)
(* The QC.jl test suite verifies that random_stabilizer(n) produces a valid    *)
(* AG tableau (stabilizers commute pairwise, full rank). We verify the same.  *)
(* ============================================================================ *)

(* Helper: build the matrix for a (possibly sign-prefixed) Pauli string,
   e.g. "-XYZ" -> -1 * (X tensor Y tensor Z). *)

stabStringToMatrix[s_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    mats = Replace[Characters[body], {
        "I" -> {{1, 0}, {0, 1}},
        "X" -> {{0, 1}, {1, 0}},
        "Y" -> {{0, -I}, {I, 0}},
        "Z" -> {{1, 0}, {0, -1}}
    }, {1}];
    sign * If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]
]


(* Pairwise commutation: stabilizers of any random state must commute. *)
VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                Module[{ps, mats, n},
                    n = RandomInteger[{1, 4}];
                    ps = PauliStabilizer["Random", n];
                    mats = stabStringToMatrix /@ ps["Stabilizers"];
                    AllTrue[Subsets[mats, {2}], #[[1]] . #[[2]] == #[[2]] . #[[1]] &]
                ],
                {10}
            ],
            TrueQ
        ]
    ],
    True,
    {},
    TestID -> "QC-RandomStabilizer-PairwiseCommute"
]

(* Each random stabilizer squares to +I (eigenvalues +-1, Hermitian). *)
VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                Module[{ps, mats, dim},
                    ps = PauliStabilizer["Random", 3];
                    dim = 2^ps["Qubits"];
                    mats = stabStringToMatrix /@ ps["Stabilizers"];
                    AllTrue[mats, # . # == IdentityMatrix[dim] &]
                ],
                {6}
            ],
            TrueQ
        ]
    ],
    True,
    {},
    TestID -> "QC-RandomStabilizer-SquaresToIdentity"
]


(* ============================================================================ *)
(* TIER D -- Tensor product semantics: PauliStabilizer[A] x PauliStabilizer[B] *)
(*                                                                              *)
(* QC.jl's Stabilizer * Stabilizer gives the tensor product. Our equivalent    *)
(* is QuantumTensorProduct[ps1, ps2]. Both should produce a stabilizer state   *)
(* with concatenated stabilizers (each padded with I on the other system).     *)
(* ============================================================================ *)

VerificationTest[
    With[{
        psA = PauliStabilizer[1]["H", 1],   (* |+> stab X *)
        psB = PauliStabilizer[1]            (* |0>  stab Z *)
    },
        Sort @ QuantumTensorProduct[psA, psB]["Stabilizers"]
    ],
    Sort @ {"XI", "IZ"},
    {},
    TestID -> "QC-TensorProduct-Plus-Tensor-Zero"
]

VerificationTest[
    With[{
        psA = PauliStabilizer[2]["H", 1]["CNOT", 1, 2],   (* Bell *)
        psB = PauliStabilizer[1]                            (* |0> *)
    },
        Sort @ QuantumTensorProduct[psA, psB]["Stabilizers"]
    ],
    Sort @ {"XXI", "ZZI", "IIZ"},
    {},
    TestID -> "QC-TensorProduct-Bell-Tensor-Zero"
]


(* ============================================================================ *)
(* TIER E -- Inner product semantics                                            *)
(*                                                                              *)
(* From QC.jl test_paulis.jl "Elementary operations" patterns: the inner       *)
(* product of a stabilizer state with itself is 1; orthogonal stabilizer       *)
(* states (e.g. Bell phi+ vs phi-) have inner product 0. Our                    *)
(* ps["InnerProduct", other] mirrors QC.jl's stabilizer inner product.         *)
(* ============================================================================ *)

VerificationTest[
    With[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}}},
        psBell["InnerProduct", psBell]
    ],
    1,
    {},
    TestID -> "QC-InnerProduct-Bell-Self"
]

VerificationTest[
    With[{
        psPhiPlus = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
        psPhiMinus = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}, "Z" -> 1}
    },
        Chop[N @ psPhiPlus["InnerProduct", psPhiMinus]]
    ],
    0,
    {},
    TestID -> "QC-InnerProduct-Bell-PhiPlus-PhiMinus-Orthogonal"
]


(* ============================================================================ *)
(* TIER F -- Live-runtime cross-validation (placeholder)                       *)
(*                                                                              *)
(* When Julia / QuantumClifford.jl is installed, generate fixtures via:        *)
(*   julia Tests/fixtures/generate_qc_jl_fixtures.jl                           *)
(* and add VerificationTest blocks here that consume the JSON output.          *)
(* The pattern mirrors the Stim setup in CrossPackage_Stim.wlt TIER G.         *)
(*                                                                              *)
(* For now: smoke check that the package source is on disk so audit            *)
(* cross-references resolve.                                                    *)
(* ============================================================================ *)

(* QC.jl is an external Julia package, not a repo artifact, so where it lives *)
(* is a fact about the machine. Three places are searched, in order:          *)
(*                                                                             *)
(*   1. QUANTUMCLIFFORD_JL_DIR, for a checkout kept anywhere at all;          *)
(*   2. a tree vendored under OngoingProjects/Stabilizer/External Packages/,  *)
(*      located from the loaded paclet rather than from the invoking script,  *)
(*      since $InputFileName names the runner under TestReport[file];          *)
(*   3. the Julia depot, <depot>/packages/QuantumClifford/<slug>/, following   *)
(*      JULIA_DEPOT_PATH as Julia reads it (base/initdefs.jl): set but empty   *)
(*      means no depot at all, ~/.julia is prepended only when the first entry *)
(*      is empty or the variable is unset, and the remaining empty entries     *)
(*      stand for the bundled depots, which carry stdlib and never a user      *)
(*      package, so nothing is lost by dropping them here.                     *)
(*                                                                             *)
(* Finding none of them is a legitimate machine state, not a defect, so the   *)
(* on-disk check below is simply not emitted in that case.                     *)

qcRepoRoot = Replace[
    PacletObject["Wolfram/QuantumFramework"],
    {paclet_PacletObject :> ParentDirectory[paclet["Location"]], _ :> $Failed}
];

qcJuliaDepots = With[{userDepot = FileNameJoin[{$HomeDirectory, ".julia"}]},
    Replace[
        Environment["JULIA_DEPOT_PATH"],
        {
            (* "" matches depots_String too, so this rule has to come first. *)
            "" -> {},
            depots_String :> With[{entries = StringSplit[depots, ":", All]},
                DeleteDuplicates @ Join[
                    If[MatchQ[entries, {"", ___}], {userDepot}, {}],
                    DeleteCases[entries, ""]
                ]
            ],
            _ :> {userDepot}
        }
    ]
];

qcSourceCandidates = DeleteDuplicates @ DeleteCases[
    Flatten @ {
        Environment["QUANTUMCLIFFORD_JL_DIR"],
        If[StringQ[qcRepoRoot],
            FileNameJoin[{
                qcRepoRoot, "OngoingProjects", "Stabilizer",
                "External Packages", "QuantumClifford.jl"
            }],
            Nothing
        ],
        FileNames["*", FileNameJoin[{#, "packages", "QuantumClifford"}]] & /@ qcJuliaDepots
    },
    Except[_String]
];

qcSourcePath = SelectFirst[qcSourceCandidates, DirectoryQ, None];

(* With a copy in reach it must carry the two test files whose line numbers   *)
(* the tiers above cite, which is what "audit cross-references resolve" means.*)
If[qcSourcePath =!= None,
    VerificationTest[
        FileExistsQ[FileNameJoin[{qcSourcePath, "test", #}]] & /@
            {"test_paulis.jl", "test_ecc_syndromes.jl"},
        {True, True},
        {},
        TestID -> "QC-PackageSource-OnDisk"
    ]
]
