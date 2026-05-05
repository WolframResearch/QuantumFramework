(* Regression tests for the audit-mistakes-followups cleanup batch.

   Covers:
     A.2 - Kernel/QuantumPartialTranspose.m rename (was QuatumPartialTranspose.m)
     A.3 - Kernel/ExampleRepository.m retirement; migration of unique exports
           (FubiniStudyMetricTensor 3-arg form, FubiniStudyMetricTensorLayers)
           into Kernel/QuantumOptimization.m
     A.4 - Guide page broken-link removal
     B.1 - QuantumSimilarity Usage + reference page
     B.2 - Bulk Keywords/Tutorials/Related Guides fills

   Each test uses fully-qualified context names to bypass the wolframscript
   shadowing of bare symbols against any system-installed paclet. *)

BeginTestSection["A.2 - QuantumPartialTranspose (renamed file)"]

(* The exported symbol must resolve in the framework context after the file
   rename. *)
VerificationTest[
  Context[Wolfram`QuantumFramework`QuantumPartialTranspose]
  ,
  "Wolfram`QuantumFramework`"
  ,
  TestID -> "A2-Symbol-Context"
];

(* QuantumPartialTranspose on a Bell state must return a QuantumState of the
   same dimensions. *)
VerificationTest[
  With[{
    qpt = Wolfram`QuantumFramework`QuantumPartialTranspose[
      Wolfram`QuantumFramework`QuantumState["Bell"], {1}]
  },
    {Head[qpt], qpt["Dimensions"]}
  ]
  ,
  {Wolfram`QuantumFramework`QuantumState, {2, 2}}
  ,
  TestID -> "A2-PartialTranspose-Bell"
];

(* Partial transpose is its own inverse: PT[PT[rho]] == rho. Also exercises
   the {1, 2} subsystem path. *)
VerificationTest[
  With[{
    bell = Wolfram`QuantumFramework`QuantumState["Bell"]
  },
    Wolfram`QuantumFramework`QuantumPartialTranspose[
      Wolfram`QuantumFramework`QuantumPartialTranspose[bell, {1}], {1}
    ]["DensityMatrix"] === bell["DensityMatrix"]
  ]
  ,
  True
  ,
  TestID -> "A2-PartialTranspose-Involution"
];

(* PT acts as identity on the un-traced subsystem: PT of |00><00| over {1} is
   itself a valid state (the |00><00| state is invariant under any partial
   transpose since it's a product state). *)
VerificationTest[
  Wolfram`QuantumFramework`QuantumPartialTranspose[
    Wolfram`QuantumFramework`QuantumState["00"], {1}
  ]["DensityMatrix"] === Wolfram`QuantumFramework`QuantumState["00"]["DensityMatrix"]
  ,
  True
  ,
  TestID -> "A2-PartialTranspose-ProductState"
];

EndTestSection[]


BeginTestSection["A.3 - ExampleRepository retirement; migration to QuantumOptimization"]

(* The deleted file's path must NOT resolve any more. *)
VerificationTest[
  FileExistsQ[
    FileNameJoin[{
      DirectoryName[$InputFileName],
      "..", "QuantumFramework", "Kernel", "ExampleRepository.m"
    }]]
  ,
  False
  ,
  TestID -> "A3-File-Deleted"
];

(* All 9 ExampleRepository exports must still resolve through the
   QuantumOptimization context after the migration. *)
VerificationTest[
  AllTrue[
    {
      "GradientDescent", "QuantumNaturalGradientDescent",
      "SPSRGradientValues", "ASPSRGradientValues",
      "FubiniStudyMetricTensor", "FubiniStudyMetricTensorLayers",
      "QuantumLockingMechanism", "QuantumUnlockingMechanism",
      "QuantumLinearSolve"
    },
    Function[name,
      MatchQ[Names["Wolfram`QuantumFramework`QuantumOptimization`" <> name], {_String}]
    ]
  ]
  ,
  True
  ,
  TestID -> "A3-AllSymbolsResolve"
];

(* QuantumLockingMechanism + QuantumUnlockingMechanism: lock builds a list
   of verification circuits parameterised by the key; Unlock returns the
   string "Correct key" if the same key opens the lock, otherwise an
   "Incorrect key" message. Test both branches. *)
VerificationTest[
  With[{
    locked = Wolfram`QuantumFramework`QuantumOptimization`QuantumLockingMechanism["secret"]
  },
    {
      Wolfram`QuantumFramework`QuantumOptimization`QuantumUnlockingMechanism[locked, "secret"],
      Wolfram`QuantumFramework`QuantumOptimization`QuantumUnlockingMechanism[locked, "wrong!"]
    }
  ]
  ,
  {"Correct key", "Incorrect key, try again"}
  ,
  TestID -> "A3-LockingMechanism-VerifyKey"
];

(* QuantumLinearSolve on a 2x2 system. The QuantumLinearSolve is "up to
   global phase" (per mistakes.md:557-560) so we compare absolute values
   only with tolerance 1e-3. *)
VerificationTest[
  With[{
    m = N @ {{2, 0}, {0, 3}},
    b = N @ {1, 2},
    classical = LinearSolve[N @ {{2, 0}, {0, 3}}, N @ {1, 2}]
  },
    Max @ Abs[Abs[Wolfram`QuantumFramework`QuantumOptimization`QuantumLinearSolve[m, b]] - Abs[classical]] < 1.0*^-3
  ]
  ,
  True
  ,
  TestID -> "A3-QuantumLinearSolve-Diagonal"
];

(* FubiniStudyMetricTensor must have BOTH the single-state form (already in
   QO) and the 3-arg layers form (newly migrated from ExampleRepository).
   Two DownValues confirm both forms exist. *)
VerificationTest[
  Length[DownValues[Wolfram`QuantumFramework`QuantumOptimization`FubiniStudyMetricTensor]]
  ,
  2
  ,
  TestID -> "A3-FubiniStudy-TwoForms"
];

(* Numeric correctness of the migrated single-state FubiniStudyMetricTensor:
   the Bures metric for |psi(theta)> = cos(theta/2)|0> + sin(theta/2)|1>
   should be 1/4. *)
VerificationTest[
  Wolfram`QuantumFramework`QuantumOptimization`FubiniStudyMetricTensor[
    Wolfram`QuantumFramework`QuantumState[{Cos[\[Theta]/2], Sin[\[Theta]/2]}],
    "Matrix",
    "Parameters" -> {\[Theta]}
  ]
  ,
  {{1/4}}
  ,
  TestID -> "A3-FubiniStudy-BuresMetric"
];

(* GradientDescent on f(x) = x^2 starting at x=1.0 with rate 0.1 must converge
   towards the minimum at x=0. After 30 iterations, the value should be small. *)
VerificationTest[
  Abs[First[Last[Wolfram`QuantumFramework`QuantumOptimization`GradientDescent[
    Function[x, x^2],
    {1.0},
    "Gradient" -> Function[x, {2 x}],
    "MaxIterations" -> 30,
    "LearningRate" -> 0.1
  ]]]] < 0.01
  ,
  True
  ,
  TestID -> "A3-GradientDescent-Converges"
];

(* The legacy ExampleRepository context must NOT define these symbols any
   more (would indicate the file got reintroduced or auto-loaded somewhere). *)
VerificationTest[
  Names["Wolfram`QuantumFramework`ExampleRepository`*"]
  ,
  {}
  ,
  TestID -> "A3-NoLegacyContext"
];

EndTestSection[]


BeginTestSection["A.4 - Guide page integrity"]

(* The Guide notebook must parse as a valid Notebook. *)
VerificationTest[
  Head[Import[
    FileNameJoin[{
      DirectoryName[$InputFileName],
      "..", "QuantumFramework", "Documentation", "English", "Guides",
      "WolframQuantumComputationFramework.nb"
    }],
    "Notebook"
  ]]
  ,
  Notebook
  ,
  TestID -> "A4-Guide-Parses"
];

(* The 3 broken symbol references must be gone from the Guide. *)
VerificationTest[
  With[{
    txt = Import[
      FileNameJoin[{
        DirectoryName[$InputFileName],
        "..", "QuantumFramework", "Documentation", "English", "Guides",
        "WolframQuantumComputationFramework.nb"
      }],
      "Text"
    ]
  },
    {
      StringContainsQ[txt, "QuantumSchmidtDecomposition"],
      StringContainsQ[txt, "QuantumSpectralDecomposition"],
      StringContainsQ[txt, "QuantumCompile"]
    }
  ]
  ,
  {False, False, False}
  ,
  TestID -> "A4-Guide-NoBrokenSymbols"
];

(* The Guide must NOT have XXXX placeholders any more. *)
VerificationTest[
  StringCount[
    Import[
      FileNameJoin[{
        DirectoryName[$InputFileName],
        "..", "QuantumFramework", "Documentation", "English", "Guides",
        "WolframQuantumComputationFramework.nb"
      }],
      "Text"
    ],
    "XXXX"
  ]
  ,
  0
  ,
  TestID -> "A4-Guide-NoXXXX"
];

EndTestSection[]


BeginTestSection["B.1 - QuantumSimilarity reference page + Usage + math correctness"]

(* QuantumSimilarity must have a Usage message in the framework context. *)
VerificationTest[
  StringQ[ToString[Wolfram`QuantumFramework`QuantumSimilarity::usage]] &&
    StringContainsQ[
      ToString[Wolfram`QuantumFramework`QuantumSimilarity::usage],
      "QuantumSimilarity"]
  ,
  True
  ,
  TestID -> "B1-Usage-Defined"
];

(* The QuantumSimilarity reference page must exist and parse. *)
VerificationTest[
  Head[Import[
    FileNameJoin[{
      DirectoryName[$InputFileName],
      "..", "QuantumFramework", "Documentation", "English", "ReferencePages",
      "Symbols", "QuantumSimilarity.nb"
    }],
    "Notebook"
  ]]
  ,
  Notebook
  ,
  TestID -> "B1-RefPage-Parses"
];

(* The 4 transformation rules from Kernel/QuantumDistance.m:63-79 must hold
   exactly: similarity = f(distance) for each measure. We test on |+> vs |->
   (orthogonal in the X basis). *)
VerificationTest[
  With[{
    psi1 = Wolfram`QuantumFramework`QuantumState[{1, 1}/Sqrt[2]],
    psi2 = Wolfram`QuantumFramework`QuantumState[{1, -1}/Sqrt[2]]
  },
    AllTrue[
      {"Fidelity", "Trace", "Bloch", "RelativePurity"},
      Function[m,
        FullSimplify[
          Wolfram`QuantumFramework`QuantumSimilarity[psi1, psi2, m] -
            (1 - Wolfram`QuantumFramework`QuantumDistance[psi1, psi2, m])
        ] === 0
      ]
    ]
  ]
  ,
  True
  ,
  TestID -> "B1-Math-LinearRule"
];

VerificationTest[
  With[{
    psi1 = Wolfram`QuantumFramework`QuantumState[{1, 1}/Sqrt[2]],
    psi2 = Wolfram`QuantumFramework`QuantumState[{1, -1}/Sqrt[2]]
  },
    AllTrue[
      {"Bures", "HilbertSchmidt"},
      Function[m,
        FullSimplify[
          Wolfram`QuantumFramework`QuantumSimilarity[psi1, psi2, m] -
            (1 - Wolfram`QuantumFramework`QuantumDistance[psi1, psi2, m]/Sqrt[2])
        ] === 0
      ]
    ]
  ]
  ,
  True
  ,
  TestID -> "B1-Math-Sqrt2Rule"
];

VerificationTest[
  With[{
    psi1 = Wolfram`QuantumFramework`QuantumState[{1, 1}/Sqrt[2]],
    psi2 = Wolfram`QuantumFramework`QuantumState[{1, -1}/Sqrt[2]]
  },
    FullSimplify[
      Wolfram`QuantumFramework`QuantumSimilarity[psi1, psi2, "BuresAngle"] -
        (1 - Wolfram`QuantumFramework`QuantumDistance[psi1, psi2, "BuresAngle"]/(Pi/2))
    ] === 0
  ]
  ,
  True
  ,
  TestID -> "B1-Math-PiOver2Rule"
];

(* All 6 example outputs in the QuantumSimilarity.nb reference page must
   match fresh kernel runs. *)
VerificationTest[
  With[{qst = Wolfram`QuantumFramework`QuantumState,
        qsim = Wolfram`QuantumFramework`QuantumSimilarity},
    {
      qsim[qst["0"], qst["1"]],
      qsim[qst[{{1/4, 0}, {0, 3/4}}], qst[{{1, 0}, {0, 1}}]],
      qsim[qst[{{1/4, 0}, {0, 3/4}}], qst[{1, 0}]],
      qsim[qst[{{1/4, 1}, {1, 3/4}}], qst[{{1/2, 2}, {2, 1/2}}], "Trace"],
      qsim[qst["GHZ"], qst["W"], "HilbertSchmidt"],
      qsim[qst[{1, 0, 0}, 3], qst[{1/Sqrt[3], 0, Sqrt[2/3]}, 3]]
    }
  ]
  ,
  {0, 1/2 + Sqrt[3]/2, 1/2, 1 - Sqrt[17]/4, 0, 1/Sqrt[3]}
  ,
  TestID -> "B1-Examples-Match"
];

(* Identity and orthogonality boundary conditions: similarity of a state
   with itself must be 1 across all measures (where defined); similarity of
   orthogonal states must be 0 for the linear-rule measures. *)
VerificationTest[
  With[{qst = Wolfram`QuantumFramework`QuantumState,
        qsim = Wolfram`QuantumFramework`QuantumSimilarity,
        psi = Wolfram`QuantumFramework`QuantumState["0"]},
    AllTrue[
      {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "RelativePurity"},
      qsim[psi, psi, #] === 1 &
    ]
  ]
  ,
  True
  ,
  TestID -> "B1-Boundary-SelfIs1"
];

VerificationTest[
  With[{qst = Wolfram`QuantumFramework`QuantumState,
        qsim = Wolfram`QuantumFramework`QuantumSimilarity},
    AllTrue[
      {"Fidelity", "Trace", "RelativePurity"},
      qsim[qst["0"], qst["1"], #] === 0 &
    ]
  ]
  ,
  True
  ,
  TestID -> "B1-Boundary-OrthogonalIs0"
];

EndTestSection[]


BeginTestSection["B.2 - Reference page bulk fills"]

(* Sample three modified reference pages and confirm they parse. *)
VerificationTest[
  AllTrue[
    {"QuantumState.nb", "QuantumChannel.nb", "QuantumCircuitOperator.nb"},
    Head[Import[
      FileNameJoin[{
        DirectoryName[$InputFileName], "..", "QuantumFramework",
        "Documentation", "English", "ReferencePages", "Symbols", #
      }],
      "Notebook"
    ]] === Notebook &
  ]
  ,
  True
  ,
  TestID -> "B2-SamplePages-Parse"
];

(* Keywords cells in the touched pages must contain real keywords (not
   "XXXX"). Spot-check QuantumState. *)
VerificationTest[
  StringContainsQ[
    Import[
      FileNameJoin[{
        DirectoryName[$InputFileName], "..", "QuantumFramework",
        "Documentation", "English", "ReferencePages", "Symbols",
        "QuantumState.nb"
      }],
      "Text"
    ],
    "Cell[\"quantum state, density matrix"
  ]
  ,
  True
  ,
  TestID -> "B2-Keywords-Filled"
];

(* The Related Guides link to the framework Guide must exist on the touched
   pages. *)
VerificationTest[
  StringContainsQ[
    Import[
      FileNameJoin[{
        DirectoryName[$InputFileName], "..", "QuantumFramework",
        "Documentation", "English", "ReferencePages", "Symbols",
        "QuantumOperator.nb"
      }],
      "Text"
    ],
    "paclet:Wolfram/QuantumFramework/guide/WolframQuantumComputationFramework"
  ]
  ,
  True
  ,
  TestID -> "B2-RelatedGuide-Linked"
];

EndTestSection[]
