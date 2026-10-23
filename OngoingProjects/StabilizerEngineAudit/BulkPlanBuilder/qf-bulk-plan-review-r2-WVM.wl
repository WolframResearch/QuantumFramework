TestReportObject[<|"Title" -> Automatic, "Aborted" -> False, 
  "TestResults" -> <|7719989559251238061 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709940672352`15.655850590930537*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "2dedf688-c0d2-4b8a-adad-dab2779159d3", "TestID" -> 
        "all exact primitives certify and bypass", "TestFileName" -> "", 
       "EvaluationID" -> "64966946-d9dc-43e0-b863-fe26991cc31c", 
       "Input" -> HoldForm[Module[{calls = 0, specs, ops}, 
          ops = FromCircuitOperatorShorthand /@ Join[Thread[{"I", "H", "S", 
                 SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"]} -> 
                ConstantArray[3, 9]], Thread[{"CNOT", "CZ", "SWAP"} -> 
                ConstantArray[{1, 3}, 3]]]; stabilizerClearCaches[]; 
           specs = Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`\
stabilizerOperatorClassification = Function[Null, calls++; Missing[
                  "NotClifford"]]}, stabilizerCircuitSpecs[qc[ops]]]; 
           calls === 0 && FreeQ[specs, _Missing] && denseCheck[ops, 3]]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.282966`5.602249253499536, "CPUTimeUsed" -> 0.27085199999999965, 
       "MemoryUsed" -> 3154960, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 1693321360729441668 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709940955514`15.65585059096127*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "9839ba86-5c80-48e6-bd0f-ffb991a1f03f", "TestID" -> 
        "approximate and global phases use classifier", "TestFileName" -> "", 
       "EvaluationID" -> "3e14ae01-e135-4774-97b1-75b3841d504e", 
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
        0.008533`4.081616743464632, "CPUTimeUsed" -> 0.008472999999999953, 
       "MemoryUsed" -> 5944, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 875622453437243772 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709940964228`15.655850590962215*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "f133bc71-cb85-4b3c-98f9-c7f95570f306", "TestID" -> 
        "independent input output permutations", "TestFileName" -> "", 
       "EvaluationID" -> "d530f9ca-a616-4c50-87ab-eb5348e8d3f5", 
       "Input" -> HoldForm[And @@ 
          (denseCheck[{QuantumOperator[KroneckerProduct[{{1, 1}, {1, -1}}/
                 Sqrt[2], DiagonalMatrix[{1, I}]], #1]}, 3] & ) /@ 
           {{{1, 3}, {1, 3}}, {{3, 1}, {3, 1}}, {{1, 3}, {3, 1}}, 
            {{3, 1}, {1, 3}}}], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 0.444059`5.797955674415373, 
       "CPUTimeUsed" -> 0.4399560000000009, "MemoryUsed" -> 352704, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    7784459805187623493 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709941408521`15.655850591010436*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "5687353d-a999-435d-9143-2c4763742c1d", 
       "TestID" -> "pictures remain distinct", "TestFileName" -> "", 
       "EvaluationID" -> "60c575c9-eee6-4a92-b821-e1d9f3990a62", 
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
        0.027742`4.593652765227282, "CPUTimeUsed" -> 0.027511000000000507, 
       "MemoryUsed" -> 33696, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 3457800277856353931 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709941436432`15.655850591013465*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "96a59506-f5e2-401f-98fe-60c1c597079f", "TestID" -> 
        "custom four qubit Clifford independent dense", "TestFileName" -> "", 
       "EvaluationID" -> "a7f08d1d-c395-4cfc-9d07-b657ebbf8db3", 
       "Input" -> HoldForm[And @@ Table[With[{s = randomSpecs[4, 18]}, 
            With[{u = dense[s, 4]}, denseCheck[{QuantumOperator[
                Exp[I*(Pi/11)]*u, {4, 2, 1, 3}], "Y" -> 2}, 4]]], {4}]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        1.118272`6.19906244872732, "CPUTimeUsed" -> 1.0926989999999996, 
       "MemoryUsed" -> 985432, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 6388599441083085990 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709942554887`15.655850591134845*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "31991ffc-96f5-4cde-83d3-0ace172a9ced", "TestID" -> 
        "five qubit signed random dense", "TestFileName" -> "", 
       "EvaluationID" -> "84be0a33-571c-41d1-930a-c2a114d5c737", 
       "Input" -> HoldForm[And @@ Table[With[{s = randomSpecs[5, 35]}, 
            denseCheck[s, 5]], {3}]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 1.131884`6.2043169187290435, 
       "CPUTimeUsed" -> 1.1210909999999998, "MemoryUsed" -> 842432, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    8527572030742399195 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709943686974`15.655850591257707*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "8378d3d3-8d55-4b2b-9ee5-b06f963da8b2", 
       "TestID" -> "arbitrary symbolic signs", "TestFileName" -> "", 
       "EvaluationID" -> "55ff092d-bd65-4b49-b520-d55c8e5fd767", 
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
        0.074478`5.022543003605683, "CPUTimeUsed" -> 0.07351299999999927, 
       "MemoryUsed" -> 81408, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 287759855703468482 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709943761662`15.655850591265812*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "87aa5166-b4c2-4916-8093-9ba8d6790186", "TestID" -> 
        "symbolic phase repeated on moved local wires", "TestFileName" -> "", 
       "EvaluationID" -> "9b8b94c1-bf33-469c-8674-9d607ef92fcc", 
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
       "AbsoluteTimeUsed" -> 0.090264`5.106029572965999, 
       "CPUTimeUsed" -> 0.08917700000000028, "MemoryUsed" -> 197064, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    1101070687696261247 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709943852094`15.655850591275627*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "1bd30ea8-9c2f-4dbd-9e43-218e8f1fcdea", 
       "TestID" -> "occurrence symbolic token fallback", 
       "TestFileName" -> "", "EvaluationID" -> 
        "220eaeda-5217-4c81-9f29-5b0dc3c19996", 
       "Input" -> HoldForm[With[{p = QuantumOperator["P"[ArcCos[t]]]}, 
          With[{ops = {p, QuantumOperator[p, {3}], QuantumOperator[p, {2}]}}, 
           stabilizerCircuitSpecs[qc[ops]] === Catenate[
             stabilizerOperatorSpecs /@ ops]]]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.253442`5.554393584925487, "CPUTimeUsed" -> 0.24611799999999917, 
       "MemoryUsed" -> 27120288, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 7767383306710514040 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709944105717`15.655850591303155*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "2fc9726b-f034-4ebd-bcfe-da525a3c8540", "TestID" -> 
        "composite misleading label refused", "TestFileName" -> "", 
       "EvaluationID" -> "206d773c-23ef-4916-9b0a-a0c946e98089", 
       "Input" -> HoldForm[qc[{QuantumOperator["T"][QuantumOperator["H"]], 
            QuantumOperator[QuantumOperator["T"][QuantumOperator["H"]], 
             "Label" -> "H"]}][PauliStabilizer[1]]], "ExpectedMessages" -> 
        HoldForm[{PauliStabilizer::nonclifford}], "ActualOutput" -> 
        HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "H" -> {1}, "its \
matrix is not a Clifford, a product of one-qubit Cliffords and phase gates, \
or the product of the gates its label names"]]}, "AbsoluteTimeUsed" -> 
        0.076628`5.034902488398482, "CPUTimeUsed" -> 0.07729099999999978, 
       "MemoryUsed" -> 90552, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[$Failed], "Outcome" -> "Success", 
       "FailureType" -> None|>], 1990677175226400609 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709944182528`15.655850591311488*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "730017ad-011e-466f-8c3e-2b32519f6a0b", "TestID" -> 
        "moved wire shares data refused", "TestFileName" -> "", 
       "EvaluationID" -> "9c94b46b-0d59-4f06-bb16-7b7858d98faf", 
       "Input" -> HoldForm[qc[{QuantumOperator[PauliMatrix[1], {{1}, {1}}], 
            QuantumOperator[PauliMatrix[1], {{2}, {1}}]}][
          PauliStabilizer[2]]], "ExpectedMessages" -> 
        HoldForm[{PauliStabilizer::nonclifford}], "ActualOutput" -> 
        HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "U"[Pi, 0, Pi] -> 
            {{2}, {1}}, "it moves a qubit to another wire"]]}, 
       "AbsoluteTimeUsed" -> 0.049971`4.849233038291919, 
       "CPUTimeUsed" -> 0.04951199999999911, "MemoryUsed" -> 51360, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[$Failed], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    5756945234631181745 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709944232677`15.65585059131693*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "8ed28d8d-7b87-4314-bdec-05ecb9cdc2c6", 
       "TestID" -> "noncomputational basis refused", "TestFileName" -> "", 
       "EvaluationID" -> "448aaa54-da41-48b8-9a7b-63befb4a1838", 
       "Input" -> HoldForm[qc[{QuantumOperator["H"], QuantumOperator[
             QuantumOperator["H"], QuantumBasis["PauliY"]]}][
          PauliStabilizer[1]]], "ExpectedMessages" -> 
        HoldForm[{PauliStabilizer::nonclifford}], "ActualOutput" -> 
        HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "PauliY" -> {1}, 
           "it is stored in a basis that is not computational"]]}, 
       "AbsoluteTimeUsed" -> 0.027049`4.582666211728936, 
       "CPUTimeUsed" -> 0.025727999999999973, "MemoryUsed" -> 260456, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[$Failed], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    7689347875872562702 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.00170994425989`15.655850591319886*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "74b957a5-9e76-414d-b727-f852cbfeb614", 
       "TestID" -> "one ququart refused", "TestFileName" -> "", 
       "EvaluationID" -> "ef48a13f-c4b9-4787-b94b-7f12aed16e28", 
       "Input" -> HoldForm[qc[{QuantumOperator[IdentityMatrix[4], 
             QuantumBasis[4]], QuantumOperator["H"]}][PauliStabilizer[1]]], 
       "ExpectedMessages" -> HoldForm[{PauliStabilizer::nonclifford}], 
       "ActualOutput" -> HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "I"[4] -> {1}, 
           "it acts on a wire that is not a qubit"]]}, 
       "AbsoluteTimeUsed" -> 0.014185`4.302344337963865, 
       "CPUTimeUsed" -> 0.014125000000000831, "MemoryUsed" -> 31992, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[$Failed], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    6774580445475100484 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709944274243`15.655850591321443*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "e3a23520-0458-41b8-b3ee-0d0b94ac93fd", 
       "TestID" -> "nonunitary refused", "TestFileName" -> "", 
       "EvaluationID" -> "97a793d7-90db-4cab-be8a-49bcb4b5259e", 
       "Input" -> HoldForm[qc[{"H" -> 1, QuantumOperator[
             2*IdentityMatrix[2]], "S" -> 1}][PauliStabilizer[1]]], 
       "ExpectedMessages" -> HoldForm[{PauliStabilizer::nonclifford}], 
       "ActualOutput" -> HoldForm[$Failed], "ActualMessages" -> 
        {HoldForm[Message[PauliStabilizer::nonclifford, "U"[0, 0, 0] -> {1}, 
           "its matrix is not unitary"]]}, "AbsoluteTimeUsed" -> 
        0.020161`4.455027067456229, "CPUTimeUsed" -> 0.019906999999999897, 
       "MemoryUsed" -> 27952, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[$Failed], "Outcome" -> "Success", 
       "FailureType" -> None|>], 7618529951015624393 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709944294563`15.655850591323645*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "ff141dfe-7e59-4919-9f04-aaf264adf160", "TestID" -> 
        "nested mixed empty", "TestFileName" -> "", "EvaluationID" -> 
        "71f8b3ae-8da0-4efb-946f-b923d62696f0", 
       "Input" -> HoldForm[With[{s = randomSpecs[3, 20]}, 
          eq[qc[{qc[Take[s, 10]], "Barrier", qc[Drop[s, 10]], qc[{}]}][
            PauliStabilizer[3]], Fold[#1[#2] & , PauliStabilizer[3], s]]]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.02746`4.58921553073273, "CPUTimeUsed" -> 0.027091999999999672, 
       "MemoryUsed" -> 19752, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 8046671356386764080 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709944322167`15.655850591326645*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "f3cabc5f-2ecb-47b7-933a-64ca8bfcf0e8", "TestID" -> 
        "empty packed identity", "TestFileName" -> "", 
       "EvaluationID" -> "1ced8162-0216-40a3-8466-49e228a4f8e9", 
       "Input" -> HoldForm[qc[{}][PauliStabilizer[63]] === 
          PauliStabilizer[63]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 0.008581`4.084052899803698, 
       "CPUTimeUsed" -> 0.00849799999999945, "MemoryUsed" -> 12008, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    7702525012733844263 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709944330894`15.65585059132759*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "d514fe0e-034a-4e65-8c1f-815e78d427c7", 
       "TestID" -> "packed boundary scalar independent", 
       "TestFileName" -> "", "EvaluationID" -> 
        "5b6b9414-f7b7-4b02-9ecb-60261889a5b0", 
       "Input" -> HoldForm[And @@ Table[With[{s = randomSpecs[n, 120], 
             init = PauliStabilizer[n]["ApplyCircuit", {"H" -> n, 
                "Y" -> 1}]}, eq[qc[s][init], Fold[#1[#2] & , 
              PauliStabilizer[Association["Tableau" -> init["Tableau"], 
                "Signs" -> init["Signs"]]], s]]], {n, {31, 32, 61, 62, 63, 
             124, 125, 1024}}]], "ExpectedMessages" -> HoldForm[{}], 
       "ActualOutput" -> HoldForm[True], "ActualMessages" -> {}, 
       "AbsoluteTimeUsed" -> 2.808058`6.598921071667187, 
       "CPUTimeUsed" -> 3.243465999999999, "MemoryUsed" -> 4527808, 
       "IntermediateTestsTree" -> {}, "ExpectedOutput" -> HoldForm[True], 
       "Outcome" -> "Success", "FailureType" -> None|>], 
    5477293934455215111 -> TestObject[<|"MetaInformation" -> None, 
       "AbsoluteTime" -> 4.001709947139175`15.655850591632364*^9, 
       "SameTest" -> SameQ, "SameMessages" -> Testing`MessageMatchQ, 
       "MemoryConstraint" -> Infinity, "TimeConstraint" -> Infinity, 
       "CreationID" -> "635fd2c2-b179-47c9-bc03-58bfbf8be5c3", 
       "TestID" -> "twelve thousand reused block", "TestFileName" -> "", 
       "EvaluationID" -> "0f59f2a8-a583-4488-a6e2-f6b45459ef78", 
       "Input" -> HoldForm[With[{s = randomSpecs[128, 300]}, 
          eq[Nest[qc[s], PauliStabilizer[128], 40], PauliStabilizer[128][
            "ApplyCircuit", Catenate[ConstantArray[s, 40]]]]]], 
       "ExpectedMessages" -> HoldForm[{}], "ActualOutput" -> HoldForm[True], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        0.61875`5.942030209773619, "CPUTimeUsed" -> 0.6239599999999985, 
       "MemoryUsed" -> 1445608, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[True], "Outcome" -> "Success", 
       "FailureType" -> None|>], 5623373632592091034 -> 
     TestObject[<|"MetaInformation" -> None, "AbsoluteTime" -> 
        4.001709947758125`15.65585059169954*^9, "SameTest" -> SameQ, 
       "SameMessages" -> Testing`MessageMatchQ, "MemoryConstraint" -> 
        Infinity, "TimeConstraint" -> Infinity, "CreationID" -> 
        "255b4894-a873-48bb-bc16-87df2655eba7", "TestID" -> 
        "new suite independently rerun", "TestFileName" -> "", 
       "EvaluationID" -> "6f8480d0-c85b-43cf-831a-28eb7a51055c", 
       "Input" -> HoldForm[Module[{report}, 
          report = TestReport["/Users/mohammadb/Documents/GitHub/QuantumFrame\
work/Tests/Stabilizer/BulkPlanBuilder.wlt"]; {report["TestsSucceededCount"], 
            report["TestsFailedCount"]}]], "ExpectedMessages" -> 
        HoldForm[{}], "ActualOutput" -> HoldForm[{56, 0}], 
       "ActualMessages" -> {}, "AbsoluteTimeUsed" -> 
        5.750301`6.9102055762986065, "CPUTimeUsed" -> 5.899493, 
       "MemoryUsed" -> 3379208, "IntermediateTestsTree" -> {}, 
       "ExpectedOutput" -> HoldForm[{56, 0}], "Outcome" -> "Success", 
       "FailureType" -> None|>]|>, "FailureResults" -> <||>, 
  "TestsNotEvaluatedKeys" -> {}, "TestsFailedWrongResultsKeys" -> {}, 
  "TestsFailedWithMessagesKeys" -> {}, "TestsSucceededKeys" -> 
   {7719989559251238061, 1693321360729441668, 875622453437243772, 
   7784459805187623493, 3457800277856353931, 6388599441083085990, 
   8527572030742399195, 287759855703468482, 1101070687696261247, 
   7767383306710514040, 1990677175226400609, 5756945234631181745, 
   7689347875872562702, 6774580445475100484, 7618529951015624393, 
   8046671356386764080, 7702525012733844263, 5477293934455215111, 
   5623373632592091034}|>]
