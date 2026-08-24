(* Loud, message-sensitive assertions for every repaired shared-core contract. *)

ClearAll[regressionCheck];
issueCount = 0;
regressionCheck[label_String, held_HoldComplete] := Module[{value, messages},
  Block[{$MessageList = {}},
    value = Check[ReleaseHold[held], $Failed];
    messages = $MessageList
    ];
  If[TrueQ[value] && messages === {},
    WriteString[$Output, "PASS: ", label, "\n"],
    issueCount++;
    WriteString[$Output, "FAIL: ", label, " => ", ToString[value, InputForm],
      "; messages = ", ToString[messages, InputForm], "\n"]
    ]
  ];

regressionCheck["catalog source hash is present", HoldComplete[
  StringQ[catalogSourceHash] && StringLength[catalogSourceHash] === 64
  ]];

regressionCheck["all core definition cells were extracted", HoldComplete[
  catalogExtractedDefinitionCells === 33
  ]];

regressionCheck["density-matrix constructor makes a valid pure state", HoldComplete[
  With[{rho = densityMatrix[plus]},
    hermitianMatrixQ[rho] && Abs[Tr[rho] - 1] < 10^-12 &&
      Min[Re@Eigenvalues[N@rho]] >= -10^-12 && Abs[purity[rho] - 1] < 10^-12
    ]
  ]];

regressionCheck["expectation preserves complex values for non-Hermitian operators", HoldComplete[
  expectation[I X, densityMatrix[plus]] === I
  ]];

regressionCheck["stationary Hermitian basis helper exists", HoldComplete[
  Length[DownValues[stationaryHermitianBasis]] > 0
  ]];

regressionCheck["zero generator has a complete four-dimensional Hermitian stationary basis", HoldComplete[
  With[{basis = stationaryHermitianBasis[0 id2, {}]},
    Length[basis] === 4 && MatrixRank[Flatten /@ basis] === 4 &&
      And @@ (hermitianMatrixQ /@ basis) &&
      Max[Norm[lindbladian[0 id2, {}, #], "Frobenius"] & /@ basis] < 10^-12
    ]
  ]];

regressionCheck["pure dephasing has a complete two-dimensional Hermitian stationary basis", HoldComplete[
  With[{basis = stationaryHermitianBasis[0 id2, {Z}]},
    Length[basis] === 2 && MatrixRank[Flatten /@ basis] === 2 &&
      And @@ (hermitianMatrixQ /@ basis) &&
      Max[Norm[lindbladian[0 id2, {Z}, #], "Frobenius"] & /@ basis] < 10^-12
    ]
  ]];

regressionCheck["non-unique steady state is reported explicitly", HoldComplete[
  MatchQ[steadyState[0 id2, {}], Failure["NonUniqueSteadyState", _Association]]
  ]];

regressionCheck["unique amplitude-damping steady state is physical and stationary", HoldComplete[
  With[{rho = First@steadyState[0 id2, {lower}]},
    densityOperatorQ[rho] && Norm[rho - densityMatrix[ground], "Frobenius"] < 10^-12 &&
      Norm[lindbladian[0 id2, {lower}, rho], "Frobenius"] < 10^-12
    ]
  ]];

regressionCheck["empty monitored and unmonitored channel lists are supported", HoldComplete[
  With[{result = measurementStep[0 id2, {}, {}, {}, 0.01][id2/2, {}]},
    MatrixQ[result, NumericQ] && densityOperatorQ[result]
    ]
  ]];

regressionCheck["scalar and empty Hamiltonians are rejected without messages", HoldComplete[
  FailureQ[measurementStep[5, {}, {}, {}, 0.01]] &&
    FailureQ[measurementStep[{}, {}, {}, {}, 0.01]]
  ]];

regressionCheck["non-list channel containers are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, Z, {1.}, {}, 0.01]] &&
    FailureQ[measurementStep[0 id2, {Z}, 1., {}, 0.01]] &&
    FailureQ[measurementStep[0 id2, {Z}, {1.}, Z, 0.01]]
  ]];

regressionCheck["invalid efficiency is rejected", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {-0.1}, {}, 0.01]] &&
    FailureQ[measurementStep[0 id2, {Z}, {1.1}, {}, 0.01]]
  ]];

regressionCheck["mismatched watched and efficiency lists are rejected", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {}, {}, 0.01]] &&
    FailureQ[measurementStep[0 id2, {}, {1.}, {}, 0.01]]
  ]];

regressionCheck["mismatched record length is rejected", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01][densityMatrix[plus], {}]] &&
    FailureQ[measurementStep[0 id2, {}, {}, {}, 0.01][densityMatrix[plus], {0.}]]
  ]];

regressionCheck["non-list records are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01][densityMatrix[plus], 0.]]
  ]];

regressionCheck["non-unit-trace states are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01][2 densityMatrix[plus], {0.}]]
  ]];

regressionCheck["non-Hermitian states are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01][{{1., I}, {0., 0.}}, {0.}]]
  ]];

regressionCheck["non-positive states are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01][{{1.1, 0.}, {0., -0.1}}, {0.}]]
  ]];

regressionCheck["non-Hermitian Hamiltonians are rejected without messages", HoldComplete[
  FailureQ[measurementStep[{{0., 1.}, {0., 0.}}, {Z}, {1.}, {}, 0.01]]
  ]];

regressionCheck["complex records are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01][densityMatrix[plus], {I}]]
  ]];

regressionCheck["complex efficiencies are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {0.5 I}, {}, 0.01]]
  ]];

regressionCheck["complex time steps are rejected without messages", HoldComplete[
  FailureQ[measurementStep[0 id2, {Z}, {1.}, {}, 0.01 I]]
  ]];

regressionCheck["finite-dimensional oscillator commutator has the cutoff correction", HoldComplete[
  With[{n = 7, aa = annihilation[7]},
    aa . ConjugateTranspose[aa] - ConjugateTranspose[aa] . aa ==
      IdentityMatrix[n] - n SparseArray[{{n, n} -> 1}, {n, n}]
    ]
  ]];

regressionCheck["modal spectrum rejects a defective eigenbasis", HoldComplete[
  With[{jordan = {{0., 1., 0., 0.}, {0., 0., 0., 0.},
       {0., 0., -1., 0.}, {0., 0., 0., -2.}}},
    FailureQ[regressionSpectrum[jordan, id2/2, Z, {1.}]]
    ]
  ]];

regressionCheck["QPC spectrum is finite at the exceptional point", HoldComplete[
  With[{spectrum = qpcSpectrum[4 OmQpc, {-1., 0., 1.}]},
    VectorQ[spectrum, NumericQ] && FreeQ[spectrum, Indeterminate | ComplexInfinity]
    ]
  ]];

regressionCheck["resolvent spectrum follows the displayed exp(-i omega t) convention", HoldComplete[
  With[{big = DiagonalMatrix[{-1 + 2 I, -2., -3., -4.}],
      seed = {{1., 0.}, {0., 0.}}, readout = {{1., 0.}, {0., 0.}}},
    Max@Abs[resolventSpectrum[big, seed, readout, {-2., 2.}] - {2/17., 2.}] < 10^-12
    ]
  ]];

WriteString[$Output, "OPEN ISSUES: ", ToString[issueCount], "\n"];
Exit[If[issueCount === 0, 0, 1]];
