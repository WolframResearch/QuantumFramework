(* ::Package:: *)

QECRunTests[] := Module[{fiveQubitChecks, tests, failures, logicalsConsistentQ},
  logicalsConsistentQ[code_] := Module[{lo = LogicalOperators[code], gens = code["Generators"], xs, zs, kk},
    xs = lo["X"]; zs = lo["Z"]; kk = Length[xs];
    AllTrue[Flatten[Outer[PauliCommuteQ, Join[xs, zs], gens]], TrueQ] &&
    And @@ Flatten[Table[PauliCommuteQ[xs[[i]], zs[[j]]] === (i =!= j), {i, kk}, {j, kk}]] &&
    AllTrue[Flatten[Outer[PauliCommuteQ, xs, xs]], TrueQ] &&
    AllTrue[Flatten[Outer[PauliCommuteQ, zs, zs]], TrueQ] &&
    AllTrue[Join[xs, zs], ! InStabilizerGroupQ[code, PauliToSymplectic[#]] &]
  ];
  fiveQubitChecks = {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"};
  SeedRandom[1];
  tests = {
    VerificationTest[PauliToSymplectic["XI"], {1, 0, 0, 0}],
    VerificationTest[PauliToSymplectic["IX"], {0, 1, 0, 0}],
    VerificationTest[PauliToSymplectic["ZI"], {0, 0, 1, 0}],
    VerificationTest[PauliToSymplectic["IZ"], {0, 0, 0, 1}],
    VerificationTest[PauliToSymplectic["YZ"], {1, 0, 1, 1}],
    VerificationTest[PauliToSymplectic["XZZXI"], {1, 0, 0, 1, 0, 0, 1, 1, 0, 0}],
    VerificationTest[PauliToSymplectic["-ZI"], {0, 0, 1, 0}],
    VerificationTest[Quiet[PauliToSymplectic["XQZ"]], $Failed],
    VerificationTest[
      With[{p = PauliStabilizer[3]["H", 1]["CNOT", 1, 2]["CNOT", 2, 3]},
        PauliToSymplectic /@ Join[p["Destabilizers"], p["Stabilizers"]] === Normal[p["Matrix"]]
      ], True],
    VerificationTest[PauliStringSign["-ZI"], -1],
    VerificationTest[PauliStringSign["ZI"], 1],
    VerificationTest[PauliStringSign["+ZI"], 1],
    VerificationTest[PauliStringSign["ZZZZZ"], 1],
    VerificationTest[SymplecticToPauli[{1, 0, 0, 1, 0, 0, 1, 1, 0, 0}], "XZZXI"],
    VerificationTest[SymplecticToPauli[{0, 0, 1, 1}], "ZZ"],
    VerificationTest[Quiet[SymplecticToPauli[{1, 0, 1}]], $Failed],
    VerificationTest[
      AllTrue[
        Table[
          With[{s = StringJoin[RandomChoice[{"I", "X", "Y", "Z"}, RandomInteger[{1, 12}]]]},
            SymplecticToPauli[PauliToSymplectic[s]] === s
          ],
          {100}
        ],
        TrueQ
      ], True],
    VerificationTest[SymplecticInnerProduct[{1, 0}, {0, 1}], 1],
    VerificationTest[SymplecticInnerProduct[{1, 0}, {1, 0}], 0],
    VerificationTest[SymplecticInnerProduct[{1, 1}, {1, 1}], 0],
    VerificationTest[Quiet[SymplecticInnerProduct[{1, 0, 1}, {1, 0, 1}]], $Failed],
    VerificationTest[Quiet[SymplecticInnerProduct[{1, 0, 0, 0}, {0, 1, 0, 0, 0, 0}]], $Failed],
    VerificationTest[PauliCommuteQ["X", "Z"], False],
    VerificationTest[PauliCommuteQ["X", "Y"], False],
    VerificationTest[PauliCommuteQ["Y", "Z"], False],
    VerificationTest[PauliCommuteQ["X", "X"], True],
    VerificationTest[PauliCommuteQ["X", "I"], True],
    VerificationTest[PauliCommuteQ["XX", "ZZ"], True],
    VerificationTest[PauliCommuteQ["XI", "ZZ"], False],
    VerificationTest[PauliCommuteQ["XXXX", "ZZZZ"], True],
    VerificationTest[PauliCommuteQ["XXXXX", "ZZZZZ"], False],
    VerificationTest[Outer[PauliCommuteQ, fiveQubitChecks, fiveQubitChecks], ConstantArray[True, {4, 4}]],
    VerificationTest[PauliCommuteQ["IXIII", #] & /@ fiveQubitChecks, {False, True, True, True}],
    VerificationTest[PauliCommuteQ["-ZI", "XI"], False],
    VerificationTest[PauliCommuteQ[{1, 0, 0, 0}, {0, 0, 1, 0}], False],
    VerificationTest[QECCode[{"ZZI", "IZZ"}]["Parameters"], {3, 1}],
    VerificationTest[QECCode[fiveQubitChecks]["Parameters"], {5, 1}],
    VerificationTest[QECCode[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"}]["Parameters"], {7, 1}],
    VerificationTest[QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]["Parameters"], {8, 6}],
    VerificationTest[QECCode[{"ZZI", "IZZ"}]["StabilizerCount"], 2],
    VerificationTest[QECCode[{"-ZI", "IZ"}]["Signs"], {-1, 1}],
    VerificationTest[QECCode[{"ZZI", "IZZ"}]["Matrix"], {{0, 0, 0, 1, 1, 0}, {0, 0, 0, 0, 1, 1}}],
    VerificationTest[Quiet[QECCode[{"XXX", "ZZZ"}]], $Failed],
    VerificationTest[Quiet[QECCode[{"ZZI", "IZZ", "ZIZ"}]], $Failed],
    VerificationTest[Quiet[QECCode[{"ZQI", "IZZ"}]], $Failed],
    VerificationTest[Quiet[QECCode[{"ZZI", "IZZZ"}]], $Failed],
    VerificationTest[Head[QECCode[{"ZZI", "IZZ"}]["PauliStabilizer"]], PauliStabilizer],
    VerificationTest[Normal[QECCode[{"ZI", "IZ"}]["State"]["StateVector"]], {1, 0, 0, 0}],
    VerificationTest[QECCode[{"ZI", "IZ"}]["PauliStabilizer"]["Expectation", "ZZ"], 1],
    VerificationTest[QECCode[{"ZZI", "IZZ"}]["NoSuchKey"], Missing["NotFound", "NoSuchKey"]],
    VerificationTest[Syndrome[QECCode[{"ZZI", "IZZ"}], "XII"], {1, 0}],
    VerificationTest[Syndrome[QECCode[{"ZZI", "IZZ"}], "IXI"], {1, 1}],
    VerificationTest[Syndrome[QECCode[{"ZZI", "IZZ"}], "IIX"], {0, 1}],
    VerificationTest[Syndrome[QECCode[{"ZZI", "IZZ"}], "III"], {0, 0}],
    VerificationTest[Syndrome[QECCode[{"ZZI", "IZZ"}], "ZII"], {0, 0}],
    VerificationTest[Syndrome[QECCode[fiveQubitChecks], "IXIII"], {1, 0, 0, 0}],
    VerificationTest[
      Module[{code, errors, syndromes},
        code = QECCode[fiveQubitChecks];
        errors = Flatten[Table[StringJoin[ReplacePart[ConstantArray["I", 5], q -> p]], {p, {"X", "Y", "Z"}}, {q, 1, 5}]];
        syndromes = Syndrome[code, #] & /@ errors;
        Length[DeleteDuplicates[syndromes]] === 15 && FreeQ[syndromes, {0, 0, 0, 0}]
      ], True],
    VerificationTest[
      Module[{code, errors, mine, oracle},
        code = QECCode[fiveQubitChecks];
        errors = Flatten[Table[{p, q}, {p, {"X", "Y", "Z"}}, {q, 1, 5}], 1];
        mine = Syndrome[code, StringJoin[ReplacePart[ConstantArray["I", 5], #[[2]] -> #[[1]]]]] & /@ errors;
        oracle = Function[e, (1 - (PauliStabilizer["5QubitCode"][e[[1]], e[[2]]]["Expectation", #] & /@ fiveQubitChecks))/2] /@ errors;
        mine === oracle
      ], True],
    VerificationTest[Quiet[Syndrome[QECCode[{"ZZI", "IZZ"}], "XIII"]], $Failed],
    VerificationTest[Quiet[Syndrome[QECCode[{"ZZI", "IZZ"}], "XQI"]], $Failed],
    VerificationTest[WeightOneErrors[2], {"XI", "YI", "ZI", "IX", "IY", "IZ"}],
    VerificationTest[Length[SyndromeTable[QECCode[{"ZZI", "IZZ"}]]], 9],
    VerificationTest[SyndromeTable[QECCode[{"ZZI", "IZZ"}]]["XII"], {1, 0}],
    VerificationTest[SyndromeTable[QECCode[{"ZZI", "IZZ"}]]["ZII"], {0, 0}],
    VerificationTest[SyndromeTable[QECCode[fiveQubitChecks]]["IXIII"], {1, 0, 0, 0}],
    VerificationTest[PerfectCodeQ[QECCode[fiveQubitChecks]], True],
    VerificationTest[PerfectCodeQ[QECCode[{"ZZI", "IZZ"}]], False],
    VerificationTest[PerfectCodeQ[QECCode[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"}]], False],
    VerificationTest[PerfectCodeQ[QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]], False],
    VerificationTest[LookupDecoder[QECCode[{"ZZI", "IZZ"}]],
      <|{0, 0} -> "III", {1, 0} -> "XII", {1, 1} -> "IXI", {0, 1} -> "IIX"|>],
    VerificationTest[Length[LookupDecoder[QECCode[fiveQubitChecks]]], 16],
    VerificationTest[DecodeSyndrome[QECCode[fiveQubitChecks], {1, 0, 0, 0}], "IXIII"],
    VerificationTest[DecodeSyndrome[QECCode[fiveQubitChecks], {0, 0, 0, 0}], "IIIII"],
    VerificationTest[DecodeSyndrome[QECCode[{"ZZI", "IZZ"}], {0, 0}], "III"],
    VerificationTest[
      Module[{code},
        code = QECCode[fiveQubitChecks];
        AllTrue[WeightOneErrors[5], DecodeSyndrome[code, Syndrome[code, #]] === # &]
      ], True],
    VerificationTest[
      MissingQ[DecodeSyndrome[QECCode[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"}], {1, 1, 0, 0, 0, 1}]], True],
    VerificationTest[Quiet[DecodeSyndrome[QECCode[{"ZZI", "IZZ"}], {1, 0, 1}]], $Failed],
    VerificationTest[Quiet[DecodeSyndrome[QECCode[{"ZZI", "IZZ"}], {1, 2}]], $Failed],
    VerificationTest[InStabilizerGroupQ[QECCode[{"ZZI", "IZZ"}], PauliToSymplectic["ZIZ"]], True],
    VerificationTest[InStabilizerGroupQ[QECCode[{"ZZI", "IZZ"}], PauliToSymplectic["ZII"]], False],
    VerificationTest[InStabilizerGroupQ[QECCode[{"ZZI", "IZZ"}], ConstantArray[0, 6]], True],
    VerificationTest[CorrectionCycle[QECCode[{"ZZI", "IZZ"}], "IXI"]["Success"], True],
    VerificationTest[CorrectionCycle[QECCode[{"ZZI", "IZZ"}], "IXI"]["Correction"], "IXI"],
    VerificationTest[CorrectionCycle[QECCode[{"ZZI", "IZZ"}], "IXI"]["Residual"], "III"],
    VerificationTest[CorrectionCycle[QECCode[{"ZZI", "IZZ"}], "ZII"]["Success"], False],
    VerificationTest[CorrectionCycle[QECCode[fiveQubitChecks], "IIIII"]["Success"], True],
    VerificationTest[AllTrue[WeightOneErrors[5], CorrectionCycle[QECCode[fiveQubitChecks], #]["Success"] &], True],
    VerificationTest[CorrectionCycle[QECCode[fiveQubitChecks], "XXIII"]["Success"], False],
    VerificationTest[CorrectionCycle[QECCode[fiveQubitChecks], "XXIII"]["Correction"], "IIIZI"],
    VerificationTest[
      MissingQ[CorrectionCycle[QECCode[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"}], "IXIIZII"]["Correction"]], True],
    VerificationTest[
      Module[{},
        SeedRandom[42];
        AllTrue[Table[CorrectionCycle[QECCode[fiveQubitChecks]]["Success"], {20}], TrueQ]
      ], True],
    VerificationTest[QECCode["5QubitCode"]["Parameters"], {5, 1}],
    VerificationTest[QECCode["5QubitCode"]["Generators"], fiveQubitChecks],
    VerificationTest[QECCode["SteaneCode"]["Parameters"], {7, 1}],
    VerificationTest[QECCode["ShorCode"]["Parameters"], {9, 1}],
    VerificationTest[QECCode["BitFlipCode"]["Generators"], {"ZZI", "IZZ"}],
    VerificationTest[Quiet[QECCode["NoSuchCode"]], $Failed],
    VerificationTest[QECCode[PauliStabilizer[2]]["Parameters"], {2, 0}],
    VerificationTest[QECCode[PauliStabilizer[2]["X", 1]]["Signs"], {-1, 1}],
    VerificationTest[
      AllTrue[{"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"},
        Module[{c = QECCode[#], full},
          full = QECCode[CompleteStabilizerGenerators[c]];
          full =!= $Failed && full["Parameters"] === {c["Qubits"], 0}
        ] &], True],
    VerificationTest[PhysicalCorrectionCycle[QECCode["BitFlipCode"], "IXI"]["Success"], True],
    VerificationTest[PhysicalCorrectionCycle[QECCode["BitFlipCode"], "IXI"]["Correction"], "IXI"],
    VerificationTest[PhysicalCorrectionCycle[QECCode["BitFlipCode"], "ZII"]["Success"], True],
    VerificationTest[AllTrue[WeightOneErrors[5], PhysicalCorrectionCycle[QECCode["5QubitCode"], #]["Success"] &], True],
    VerificationTest[
      AllTrue[WeightOneErrors[3],
        PhysicalCorrectionCycle[QECCode["BitFlipCode"], #]["Syndrome"] === Syndrome[QECCode["BitFlipCode"], #] &], True],
    VerificationTest[CodeStandardForm[QECCode["BitFlipCode"]]["XRank"], 0],
    VerificationTest[MatrixRank[CodeStandardForm[QECCode["5QubitCode"]]["Matrix"], Modulus -> 2], 4],
    VerificationTest[Sort[CodeStandardForm[QECCode["ShorCode"]]["QubitPermutation"]], Range[9]],
    VerificationTest[LogicalOperators[QECCode["BitFlipCode"]], <|"X" -> {"XXX"}, "Z" -> {"IIZ"}|>],
    VerificationTest[
      logicalsConsistentQ[QECCode[#]] & /@ {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"},
      {True, True, True, True, True}],
    VerificationTest[logicalsConsistentQ[QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]], True],
    VerificationTest[
      Module[{code = QECCode["5QubitCode"], zbar},
        zbar = LogicalOperators[code]["Z"][[1]];
        InStabilizerGroupQ[code, Mod[PauliToSymplectic[zbar] + PauliToSymplectic["ZZZZZ"], 2]]
      ], True],
    VerificationTest[
      Module[{code = QECCode["SteaneCode"], zbar},
        zbar = LogicalOperators[code]["Z"][[1]];
        InStabilizerGroupQ[code, Mod[PauliToSymplectic[zbar] + PauliToSymplectic["ZZZZZZZ"], 2]]
      ], True],
    VerificationTest[LogicalOperators[QECCode[CompleteStabilizerGenerators[QECCode["5QubitCode"]]]], <|"X" -> {}, "Z" -> {}|>],
    VerificationTest[PauliWeight["XZZXI"], 4],
    VerificationTest[PauliWeight["III"], 0],
    VerificationTest[PauliWeight["-ZI"], 1],
    VerificationTest[Sort[WeightKPaulis[2, 1]], Sort[{"XI", "YI", "ZI", "IX", "IY", "IZ"}]],
    VerificationTest[Length[WeightKPaulis[5, 2]], 9 Binomial[5, 2]],
    VerificationTest[LogicalPauliQ[QECCode["BitFlipCode"], "IIZ"], True],
    VerificationTest[LogicalPauliQ[QECCode["BitFlipCode"], "ZZI"], False],
    VerificationTest[LogicalPauliQ[QECCode["BitFlipCode"], "IXI"], False],
    VerificationTest[CodeDistance[QECCode["BitFlipCode"]], 1],
    VerificationTest[CodeDistance[QECCode["PhaseFlipCode"]], 1],
    VerificationTest[CodeDistance[QECCode["5QubitCode"]], 3],
    VerificationTest[CodeDistance[QECCode["SteaneCode"]], 3],
    VerificationTest[CodeDistance[QECCode["ShorCode"]], 3],
    VerificationTest[CodeDistance[QECCode[{StringRepeat["X", 6], StringRepeat["Z", 6]}]], 2],
    VerificationTest[CodeDistance[QECCode[CompleteStabilizerGenerators[QECCode["5QubitCode"]]]], Infinity],
    VerificationTest[PauliWeight[MinimumWeightLogical[QECCode["5QubitCode"]]], 3],
    VerificationTest[LogicalPauliQ[QECCode["5QubitCode"], MinimumWeightLogical[QECCode["5QubitCode"]]], True],
    VerificationTest[QECCode["SteaneCode"]["FullParameters"], {7, 1, 3}],
    VerificationTest[
      AllTrue[{"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"},
        EncodingCircuitValidQ[QECCode[#]] &], True],
    VerificationTest[Head[EncodingCircuit[QECCode["5QubitCode"]]], QuantumCircuitOperator],
    VerificationTest[Length[EncodingGates[QECCode["5QubitCode"]]] > 0, True],
    VerificationTest[EncodingCircuitValidQ[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}]], True],
    VerificationTest[EncodingCircuitValidQ[QECCode[{StringRepeat["X", 6], StringRepeat["Z", 6]}]], True],
    VerificationTest[Wolfram`QuantumFramework`QEC`PackageScope`rowToPauli[{1, 0, 1, 1}, "X"], "XIXX"],
    VerificationTest[
      With[{hamming = {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}}},
        CSSCode[hamming]["FullParameters"]], {7, 1, 3}],
    VerificationTest[CSSCode[{{1, 1, 1, 1}}, {{1, 1, 1, 1}}]["Parameters"], {4, 2}],
    VerificationTest[CSSCode[{{1, 1, 1, 1}}]["Distance"], 2],
    VerificationTest[
      With[{hamming = {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}}},
        Sort[CSSCode[hamming]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]]], True],
    VerificationTest[Quiet[CSSCode[{{1, 1, 0}}, {{1, 0, 0}}]], $Failed, {}],
    VerificationTest[Quiet[CSSCode[{{1, 1}}, {{1, 1, 1}}]], $Failed, {}],
    VerificationTest[CSSCodeQ[QECCode["SteaneCode"]], True],
    VerificationTest[CSSCodeQ[QECCode["5QubitCode"]], False],
    VerificationTest[CSSCodeQ[QECCode["ShorCode"]], True],
    VerificationTest[Wolfram`QuantumFramework`QEC`PackageScope`pauliStringProduct["XIX", "IZZ"], "XZY"],
    VerificationTest[Wolfram`QuantumFramework`QEC`PackageScope`placeOnBlock["ZZ", 2, 2, 6], "IIZZII"],
    VerificationTest[Wolfram`QuantumFramework`QEC`PackageScope`translateOuterGenerator["XXI", "XXX", "IIZ", 3], "XXXXXXIII"],
    VerificationTest[ConcatenateCodes[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Parameters"], {9, 1}],
    VerificationTest[
      Sort[ConcatenateCodes[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Generators"]] === Sort[QECCode["ShorCode"]["Generators"]], True],
    VerificationTest[ConcatenateCodes[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Distance"], 3],
    VerificationTest[ConcatenateCodes[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], QECCode["BitFlipCode"]]["Parameters"], {12, 2}],
    VerificationTest[Quiet[ConcatenateCodes[QECCode["PhaseFlipCode"], QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}]]], $Failed, {}],
    VerificationTest[RemoveQubit[QECCode["5QubitCode"]]["Parameters"], {4, 2}],
    VerificationTest[RemoveQubit[QECCode["5QubitCode"]]["FullParameters"], {4, 2, 2}],
    VerificationTest[RemoveQubit[QECCode["SteaneCode"]]["Parameters"], {6, 2}],
    VerificationTest[RemoveQubit[QECCode["SteaneCode"]]["Distance"], 2],
    VerificationTest[Head[RemoveQubit[QECCode["5QubitCode"], 3]], QECCode],
    VerificationTest[Quiet[RemoveQubit[QECCode["BitFlipCode"]]], $Failed, {}],
    VerificationTest[
      PasteCodes[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], 1, 1]["Parameters"],
      {8, 5}],
    VerificationTest[
      PasteCodes[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], 1, 1]["Generators"],
      {"XXXXIIII", "IIIIXXXX", "ZZZZZZZZ"}],
    VerificationTest[
      PasteCodes[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], 1, 1]["Distance"],
      2],
    VerificationTest[Quiet[PasteCodes[QECCode["BitFlipCode"], QECCode["SteaneCode"], 1, 1]], $Failed, {}],
    VerificationTest[Quiet[PasteCodes[QECCode["BitFlipCode"], QECCode["PhaseFlipCode"], 5, 1]], $Failed, {}],
    VerificationTest[RepetitionCode[3]["Generators"], QECCode["BitFlipCode"]["Generators"]],
    VerificationTest[PhaseRepetitionCode[3]["Generators"], QECCode["PhaseFlipCode"]["Generators"]],
    VerificationTest[RepetitionCode[5]["Parameters"], {5, 1}],
    VerificationTest[Quiet[RepetitionCode[1]], $Failed, {}],
    VerificationTest[DistanceTwoCode[4]["FullParameters"], {4, 2, 2}],
    VerificationTest[DistanceTwoCode[6]["FullParameters"], {6, 4, 2}],
    VerificationTest[Quiet[DistanceTwoCode[5]], $Failed, {}],
    VerificationTest[Dimensions[ClassicalHammingMatrix[3]], {3, 7}],
    VerificationTest[ClassicalHammingMatrix[3], {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}}],
    VerificationTest[Sort[HammingCSSCode[3]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]], True],
    VerificationTest[HammingCSSCode[4]["Parameters"], {15, 7}],
    VerificationTest[Quiet[HammingCSSCode[2]], $Failed, {}],
    VerificationTest[Head[QECCodeCatalog[]], Dataset]
  };
  failures = Select[tests, #["Outcome"] =!= "Success" &];
  <|
    "Total" -> Length[tests],
    "Outcomes" -> Counts[#["Outcome"] & /@ tests],
    "FailedInputs" -> (HoldForm @@ {#["Input"]} & /@ failures)
  |>
]
