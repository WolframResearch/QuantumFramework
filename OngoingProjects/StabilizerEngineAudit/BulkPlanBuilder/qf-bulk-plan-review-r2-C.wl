TestReportObject[<|"Title" -> Automatic, "Aborted" -> False, 
  "TestResults" -> <|5862224438760584352 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709862249103`15.65585058241948*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "bcc2be42-a8e3-4ab3-b3e7-a7e406ed6514", "TestID" -> 
        "all exact primitives certify and bypass", "TestFileName" -> "", 
       "EvaluationID" -> "9d647428-e36f-4952-8c6a-d2eea2885b0c", 
       "Input" -> HoldForm[Module[{calls = 0, specs, ops}, 
          ops = FromCircuitOperatorShorthand /@ Join[Thread[{"I", "H", "S", 
                 SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"]} -> 
                ConstantArray[3, 9]], Thread[{"CNOT", "CZ", "SWAP"} -> 
                ConstantArray[{3, 1}, 3]]]; stabilizerClearCaches[]; 
           specs = Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`\
stabilizerOperatorClassification = Function[Null, calls++; Missing[
                  "NotClifford"]]}, stabilizerCircuitSpecs[qc[ops]]]; 
           calls === 0 && FreeQ[specs, _Missing] && denseCheck[ops, 3]]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[False], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.153031`5.3352944140431315, "CPUTimeUsed" -> 0.14691600000000005, 
       "MemoryUsed" -> 2307904, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Failure", 
       "FailureType" -> "SameTestFailure"|>], 514867276541816926 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709862402319`15.65585058243611*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "213fe1df-5129-482c-8354-ba7b4f7b5bde", "TestID" -> 
        "approximate and global phases use classifier", "TestFileName" -> "", 
       "EvaluationID" -> "641d9215-e846-47ef-b141-2c62e6a0eddc", 
       "Input" -> HoldForm[Module[{ops, calls = 0, specs}, 
          ops = (QuantumOperator[#1, {2}] & ) /@ {N[{{1, 1}, {1, -1}}/
                Sqrt[2]], I*PauliMatrix[1], Exp[I*(Pi/9)]*DiagonalMatrix[
                {1, I}]}; stabilizerClearCaches[]; 
           specs = Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`\
stabilizerOperatorClassification = Function[Null, calls++; Missing[
                  "NotUnitary"]]}, stabilizerCircuitSpecs[qc[ops]]]; 
           calls === 3 &&  !FreeQ[specs, _Missing]]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.008726`4.091330206482795, "CPUTimeUsed" -> 0.008642000000000483, 
       "MemoryUsed" -> 8072, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 4036946948771701622 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709862411217`15.655850582437076*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "041a7076-d49d-4321-93d2-e97eac964ac9", "TestID" -> 
        "independent input output permutations", "TestFileName" -> "", 
       "EvaluationID" -> "9ca2be57-e026-4c57-a110-b7d6f100c1c5", 
       "Input" -> HoldForm[And @@ 
          (denseCheck[{QuantumOperator[KroneckerProduct[{{1, 1}, {1, -1}}/
                 Sqrt[2], DiagonalMatrix[{1, I}]], #1]}, 3] & ) /@ 
           {{{1, 3}, {1, 3}}, {{3, 1}, {3, 1}}, {{1, 3}, {3, 1}}, 
            {{3, 1}, {1, 3}}}], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 1.765045`6.39727078007853, 
       "CPUTimeUsed" -> 0.6817669999999998, "MemoryUsed" -> 3212096, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    69741821599586381 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709864176552`15.655850582628661*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "e6ed74f3-7cb1-4b0a-aa7c-eacfa7a04d5f", 
       "TestID" -> "pictures remain distinct", "TestFileName" -> "", 
       "EvaluationID" -> "d7bb0fa4-bae9-4395-ae0e-a28452d5db2a", 
       "Input" -> HoldForm[With[{a = QuantumOperator["H"], 
           b = QuantumOperator["H", QuantumBasis["Heisenberg"]]}, 
          Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperat\
orData[a] =!= 
            Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOper\
atorData[b] && eq[qc[{a, b}][PauliStabilizer[1]], 
            Fold[#1["ApplyCircuit", stabilizerOperatorSpecs[#2]] & , 
             PauliStabilizer[1], {a, b}]]]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.030686`4.63745527853669, "CPUTimeUsed" -> 0.030577999999999328, 
       "MemoryUsed" -> 32544, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 7090948033447821345 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709864207438`15.655850582632015*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "32bc7cb1-db50-48f4-9119-fd2993d8c918", "TestID" -> 
        "custom four qubit Clifford independent dense", "TestFileName" -> "", 
       "EvaluationID" -> "6d502453-1d17-4bf7-80b3-7a44ac69ac76", 
       "Input" -> HoldForm[And @@ Table[With[{s = randomSpecs[4, 18]}, 
            With[{u = dense[s, 4]}, denseCheck[{QuantumOperator[
                Exp[I*(Pi/11)]*u, {4, 2, 1, 3}], "Y" -> 2}, 4]]], {4}]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        1.114978`6.19778129609253, "CPUTimeUsed" -> 1.1046060000000004, 
       "MemoryUsed" -> 1565456, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 45849542212452059 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709865322604`15.65585058275304*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "32e4ff62-fb12-4826-b3cb-a8bc89cf6097", "TestID" -> 
        "five qubit signed random dense", "TestFileName" -> "", 
       "EvaluationID" -> "961109f9-8006-4d94-ba63-fdf8b457322d", 
       "Input" -> HoldForm[And @@ Table[With[{s = randomSpecs[5, 35]}, 
            denseCheck[s, 5]], {3}]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 1.147344`6.210208646686331, 
       "CPUTimeUsed" -> 1.1384880000000006, "MemoryUsed" -> 829056, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    1544030826388916777 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 
        4.0017098664702219999999999999`15.655850582877587*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "83b8e91f-5a08-4cb9-a016-ebc4efbaaf4a", 
       "TestID" -> "arbitrary symbolic signs", "TestFileName" -> "", 
       "EvaluationID" -> "1a71eb04-a435-4ec1-8fdf-74d0aa104f59", 
       "Input" -> HoldForm[Module[{a = Array[sgn, 14], initial, specs, 
           actual, expected}, initial = PauliStabilizer[Association[
              "Tableau" -> PauliStabilizer[7]["Tableau"], "Signs" -> 1 - 
                2*a]]; specs = randomSpecs[7, 50]; actual = 
            qc[specs][initial]; expected = Fold[#1[#2] & , initial, specs]; 
           actual["Tableau"] === expected["Tableau"] && 
            And @@ Table[With[{rules = Thread[a -> RandomInteger[1, 14]]}, 
               Simplify[actual["Signs"] - expected["Signs"] /. rules] === 
                ConstantArray[0, 14]], {40}]]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.07697`5.036836484391474, "CPUTimeUsed" -> 0.07616900000000015, 
       "MemoryUsed" -> 81360, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 2624755168616745055 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709866547444`15.655850582885968*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "2cc58306-ae7c-42fa-838d-d18d3635775a", "TestID" -> 
        "symbolic phase repeated on moved local wires", "TestFileName" -> "", 
       "EvaluationID" -> "2794d98c-e469-4c46-b91e-0bcb8d2d73c2", 
       "Input" -> HoldForm[Module[{op, ops, out, u, v}, 
          op = QuantumOperator[DiagonalMatrix[{Exp[I*alpha], Exp[I*beta]}]]; 
           ops = {"H" -> 1, op, QuantumOperator[op, {2}], "CNOT" -> {1, 2}, 
             "H" -> 2}; out = qc[ops][PauliStabilizer[2]]; u = dense[ops, 2]; 
           v = Normal[out["State"]["StateVector"]]; 
           And @@ Flatten[Table[With[{rules = {alpha -> aa, beta -> bb}}, 
               Chop[N[Outer[Times, v, Conjugate[v]] - Outer[Times, u[[All,
                      1]], Conjugate[u[[All,1]]]] /. rules]] === 
                ConstantArray[0, {4, 4}]], {aa, {0, Pi/8, 2*(Pi/7)}}, 
              {bb, {0, -Pi/6, Pi}}]]]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 0.096851`5.136619107060999, 
       "CPUTimeUsed" -> 0.09540499999999952, "MemoryUsed" -> 196248, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    7740603960089025014 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709866644516`15.655850582896504*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "85c16dfc-5ae4-485c-b3d5-d037abe0c3ef", 
       "TestID" -> "occurrence symbolic token fallback", 
       "TestFileName" -> "", "EvaluationID" -> 
        "da7c5077-27a3-442d-9513-8d463b6d00e7", 
       "Input" -> HoldForm[With[{p = QuantumOperator["P"[ArcCos[t]]]}, 
          With[{ops = {p, QuantumOperator[p, {3}], QuantumOperator[p, {2}]}}, 
           stabilizerCircuitSpecs[qc[ops]] === Catenate[
             stabilizerOperatorSpecs /@ ops]]]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.26993`5.581766152527149, "CPUTimeUsed" -> 0.26220099999999924, 
       "MemoryUsed" -> 26805216, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 3087750128904625647 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709866914641`15.65585058292582*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "1a03a8b3-79ac-4793-b264-1ce688f6adbf", "TestID" -> 
        "composite misleading label refused", "TestFileName" -> "", 
       "EvaluationID" -> "201c4a08-8738-458c-8cd9-9f325e532620", 
       "Input" -> HoldForm[qc[{QuantumOperator["T"][QuantumOperator["H"]], 
            QuantumOperator[QuantumOperator["T"][QuantumOperator["H"]], 
             "Label" -> "H"]}][PauliStabilizer[1]]], "ExpectedMessages" -> 
        HoldForm[{PauliStabilizer::nonclifford}], "ActualOutput" -> 
        HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "H" -> {1}, "its \
matrix is not a Clifford, a product of one-qubit Cliffords and phase gates, \
or the product of the gates its label names"]]}, "AbsoluteTimeUsed" -> 
        0.078869`5.047421332163343, "CPUTimeUsed" -> 0.07975699999999986, 
       "MemoryUsed" -> 91544, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[$Failed], "Outcome" -> "Success", 
       "FailureType" -> None|>], 6850708628947041703 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709866993686`15.655850582934397*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "6eb34114-acc3-4751-8616-fcf5524f19ff", "TestID" -> 
        "moved wire shares data refused", "TestFileName" -> "", 
       "EvaluationID" -> "143be7f6-20f8-416c-8277-9225f0119db6", 
       "Input" -> HoldForm[qc[{QuantumOperator[PauliMatrix[1], {{1}, {1}}], 
            QuantumOperator[PauliMatrix[1], {{2}, {1}}]}][
          PauliStabilizer[2]]], "ExpectedMessages" -> 
        HoldForm[{PauliStabilizer::nonclifford}], "ActualOutput" -> 
        HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "U"[Pi, 0, Pi] -> 
            {{2}, {1}}, "it moves a qubit to another wire"]]}, 
       "AbsoluteTimeUsed" -> 0.050408`4.853014464439065, 
       "CPUTimeUsed" -> 0.0501049999999994, "MemoryUsed" -> 51104, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[$Failed], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    8694201709558214461 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709867044266`15.655850582939886*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "a00f110c-2e14-4318-afed-abd67d90d58d", 
       "TestID" -> "noncomputational basis refused", "TestFileName" -> "", 
       "EvaluationID" -> "30b007ad-09f4-4e8d-9d2c-7d14414d1dfa", 
       "Input" -> HoldForm[qc[{QuantumOperator["H"], QuantumOperator[
             QuantumOperator["H"], QuantumBasis["PauliY"]]}][
          PauliStabilizer[1]]], "ExpectedMessages" -> 
        HoldForm[{PauliStabilizer::nonclifford}], "ActualOutput" -> 
        HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "PauliY" -> {1}, 
           "it is stored in a basis that is not computational"]]}, 
       "AbsoluteTimeUsed" -> 0.026782`4.578358003537527, 
       "CPUTimeUsed" -> 0.025663000000000658, "MemoryUsed" -> 260776, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[$Failed], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    6507689714075545700 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709867071193`15.655850582942808*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "a5cc448e-7be5-4814-a947-3cf523319c29", 
       "TestID" -> "one ququart refused", "TestFileName" -> "", 
       "EvaluationID" -> "93bc639b-2a4c-4b72-911f-e571685cefd6", 
       "Input" -> HoldForm[qc[{QuantumOperator[IdentityMatrix[4], 
             QuantumBasis[4]], QuantumOperator["H"]}][PauliStabilizer[1]]], 
       "ExpectedMessages" -> HoldForm[{PauliStabilizer::nonclifford}], 
       "ActualOutput" -> HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "I"[4] -> {1}, 
           "it acts on a wire that is not a qubit"]]}, 
       "AbsoluteTimeUsed" -> 0.014331`4.306791493835011, 
       "CPUTimeUsed" -> 0.01423299999999994, "MemoryUsed" -> 31992, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[$Failed], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    7881694633785832797 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709867085698`15.655850582944382*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "9e083509-85d8-47dc-9669-37ec806d72f6", 
       "TestID" -> "nonunitary refused", "TestFileName" -> "", 
       "EvaluationID" -> "bf082de9-e05c-4ec0-8df2-63def9f43736", 
       "Input" -> HoldForm[qc[{"H" -> 1, QuantumOperator[
             2*IdentityMatrix[2]], "S" -> 1}][PauliStabilizer[1]]], 
       "ExpectedMessages" -> HoldForm[{PauliStabilizer::nonclifford}], 
       "ActualOutput" -> HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "U"[0, 0, 0] -> {1}, 
           "its matrix is not unitary"]]}, "AbsoluteTimeUsed" -> 
        0.020585`4.464065869165499, "CPUTimeUsed" -> 0.020398000000000138, 
       "MemoryUsed" -> 27888, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[$Failed], "Outcome" -> "Success", 
       "FailureType" -> None|>], 8149839959218374686 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709867106452`15.655850582946636*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "9470f578-ece7-4369-9bb5-a97cb724d4dd", "TestID" -> 
        "nested mixed empty", "TestFileName" -> "", "EvaluationID" -> 
        "8415b1b3-f4f1-4309-bbd6-d61782e69c5a", 
       "Input" -> HoldForm[With[{s = randomSpecs[3, 20]}, 
          eq[qc[{qc[Take[s, 10]], "Barrier", qc[Drop[s, 10]], qc[{}]}][
            PauliStabilizer[3]], Fold[#1[#2] & , PauliStabilizer[3], s]]]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.027183`4.584812382956502, "CPUTimeUsed" -> 0.027059000000001276, 
       "MemoryUsed" -> 19688, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 776548365169581405 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709867133786`15.6558505829496*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "14b887c2-af46-48b6-8514-920b0a5a490e", "TestID" -> 
        "empty packed identity", "TestFileName" -> "", 
       "EvaluationID" -> "e0c7c167-3a07-4074-939b-96f4cb3719e2", 
       "Input" -> HoldForm[qc[{}][PauliStabilizer[63]] === 
          PauliStabilizer[63]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 0.008596`4.0848114046514, 
       "CPUTimeUsed" -> 0.008519999999999861, "MemoryUsed" -> 12240, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    4606215408967447305 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709867142508`15.655850582950551*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "a7a4df17-4693-4dbd-8001-2637826136d4", 
       "TestID" -> "packed boundary scalar independent", 
       "TestFileName" -> "", "EvaluationID" -> 
        "27bdca4c-4fe2-4a4f-a194-954f029f9b3f", 
       "Input" -> HoldForm[And @@ Table[With[{s = randomSpecs[n, 120], 
             init = PauliStabilizer[n]["ApplyCircuit", {"H" -> n, 
                "Y" -> 1}]}, eq[qc[s][init], Fold[#1[#2] & , 
              PauliStabilizer[Association["Tableau" -> init["Tableau"], 
                "Signs" -> init["Signs"]]], s]]], {n, {31, 32, 61, 62, 63, 
             124, 125, 1024}}]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 2.821539`6.601001055403759, 
       "CPUTimeUsed" -> 3.2765189999999995, "MemoryUsed" -> 4519072, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    7992931833638538209 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709869964257`15.655850583256788*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "d3c77b5d-28e6-40bb-a88c-39f8cfc1eca1", 
       "TestID" -> "twelve thousand reused block", "TestFileName" -> "", 
       "EvaluationID" -> "804275ad-809a-42f9-81b3-b710ad5e372c", 
       "Input" -> HoldForm[With[{s = randomSpecs[128, 300]}, 
          eq[Nest[qc[s], PauliStabilizer[128], 40], PauliStabilizer[128][
            "ApplyCircuit", Catenate[ConstantArray[s, 40]]]]]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.611582`5.936969692663003, "CPUTimeUsed" -> 0.6177050000000008, 
       "MemoryUsed" -> 1445408, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 2619446515677916272 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709870576184`15.655850583323199*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "f60aef20-fbf4-4235-87bc-c33236fb2715", "TestID" -> 
        "new suite independently rerun", "TestFileName" -> "", 
       "EvaluationID" -> "6469b7ad-a867-452e-b4dc-dc8574ba0e9d", 
       "Input" -> HoldForm[Module[{report}, 
          report = TestReport["/Users/mohammadb/Documents/GitHub/QuantumFrame\
work/Tests/Stabilizer/BulkPlanBuilder.wlt"]; {report["TestsSucceededCount"], 
            report["TestsFailedCount"]}]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[{56, 0}], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        5.792375`6.913371668279162, "CPUTimeUsed" -> 5.927352000000001, 
       "MemoryUsed" -> 3377648, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[{56, 0}], "Outcome" -> "Success", 
       "FailureType" -> None|>]|>, "FailureResults" -> <||>, 
  "TestsNotEvaluatedKeys" -> {}, "TestsFailedWithMessagesKeys" -> {}, 
  "TestsFailedWrongResultsKeys" -> {5862224438760584352}, 
  "TestsSucceededKeys" -> {514867276541816926, 4036946948771701622, 
   69741821599586381, 7090948033447821345, 45849542212452059, 
   1544030826388916777, 2624755168616745055, 7740603960089025014, 
   3087750128904625647, 6850708628947041703, 8694201709558214461, 
   6507689714075545700, 7881694633785832797, 8149839959218374686, 
   776548365169581405, 4606215408967447305, 7992931833638538209, 
   2619446515677916272}|>]
