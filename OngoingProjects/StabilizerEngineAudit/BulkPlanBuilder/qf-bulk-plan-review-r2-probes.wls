PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
If[Environment["BPB_REVIEW_TARGET"] === "WVM", $stabilizerCompileTarget = "WVM"];
Echo[<|"Version" -> $Version, "Target" -> $stabilizerCompileTarget, "Loaded" -> FindFile["Wolfram`QuantumFramework`"]|>, "Metadata"];
qc[ops_List] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ ops), "Label" -> None|>];
eq[a_PauliStabilizer,b_PauliStabilizer] := a["Tableau"] === b["Tableau"] && a["Signs"] === b["Signs"];
eq[_,_] := False;
paulis[n_] := Join[Table[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2],n], j -> PauliMatrix[1]],{j,n}], Table[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2],n], j -> PauliMatrix[3]],{j,n}]];
images[ps_,n_] := With[{m=Normal[ps["Matrix"]],s=ps["Signs"]}, MapThread[ #1 (KroneckerProduct@@MapThread[PauliMatrix[Switch[{#1,#2},{0,0},0,{1,0},1,{1,1},2,{0,1},3]]&,{Take[#2,n],Take[#2,-n]}]) &, {s,m}]];
dense[ops_,n_] := Normal[qc[Prepend[ops,QuantumOperator[IdentityMatrix[2^n]]]]["QuantumOperator"]["Matrix"]];
denseCheck[ops_,n_] := With[{u=dense[ops,n],out=qc[ops][PauliStabilizer[n]]}, Head[out]===PauliStabilizer && And@@MapThread[TrueQ[Simplify[#1-#2]===ConstantArray[0,{2^n,2^n}]] &, {images[out,n],u.#.ConjugateTranspose[u]& /@ paulis[n]}]];
randomSpecs[n_,count_] := Table[With[{w=RandomSample[Range[n],2],g=RandomChoice[{"H","S",SuperDagger["S"],"X","Y","Z","V",SuperDagger["V"],"CNOT","CZ","SWAP"}]},g->If[MemberQ[{"CNOT","CZ","SWAP"},g],w,First[w]]],{count}];
SeedRandom[187389];
tests = {
 VerificationTest[Module[{calls=0,specs,ops}, ops=FromCircuitOperatorShorthand /@ Join[Thread[{"I","H","S",SuperDagger["S"],"X","Y","Z","V",SuperDagger["V"]}->ConstantArray[3,9]], Thread[{"CNOT","CZ","SWAP"}->ConstantArray[{1,3},3]]]; stabilizerClearCaches[]; specs=Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification=Function[Null,calls++;Missing["NotClifford"]]},stabilizerCircuitSpecs[qc[ops]]]; calls===0 && FreeQ[specs,_Missing] && denseCheck[ops,3]],True,TestID->"all exact primitives certify and bypass"],
 VerificationTest[Module[{ops,calls=0,specs},ops=QuantumOperator[#, {2}]& /@ {N[{{1,1},{1,-1}}/Sqrt[2]],I PauliMatrix[1],Exp[I Pi/9]DiagonalMatrix[{1,I}]}; stabilizerClearCaches[]; specs=Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification=Function[Null,calls++;Missing["NotUnitary"]]},stabilizerCircuitSpecs[qc[ops]]];calls===3 && !FreeQ[specs,_Missing]],True,TestID->"approximate and global phases use classifier"],
 VerificationTest[And@@(denseCheck[{QuantumOperator[KroneckerProduct[{{1,1},{1,-1}}/Sqrt[2],DiagonalMatrix[{1,I}]],#]},3]& /@ {{{1,3},{1,3}},{{3,1},{3,1}},{{1,3},{3,1}},{{3,1},{1,3}}}),True,TestID->"independent input output permutations"],
 VerificationTest[With[{a=QuantumOperator["H"],b=QuantumOperator["H",QuantumBasis["Heisenberg"]]}, Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[a]=!=Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[b] && eq[qc[{a,b}][PauliStabilizer[1]],Fold[#1["ApplyCircuit",stabilizerOperatorSpecs[#2]]&,PauliStabilizer[1],{a,b}]]],True,TestID->"pictures remain distinct"],
 VerificationTest[And@@Table[With[{s=randomSpecs[4,18]},With[{u=dense[s,4]},denseCheck[{QuantumOperator[Exp[I Pi/11]u,{4,2,1,3}],"Y"->2},4]]],{4}],True,TestID->"custom four qubit Clifford independent dense"],
 VerificationTest[And@@Table[With[{s=randomSpecs[5,35]},denseCheck[s,5]],{3}],True,TestID->"five qubit signed random dense"],
 VerificationTest[Module[{a=Array[sgn,14],initial,specs,actual,expected},initial=PauliStabilizer[<|"Tableau"->PauliStabilizer[7]["Tableau"],"Signs"->1-2a|>];specs=randomSpecs[7,50];actual=qc[specs][initial];expected=Fold[#1[#2]&,initial,specs];actual["Tableau"]===expected["Tableau"] && And@@Table[With[{rules=Thread[a->RandomInteger[1,14]]},Simplify[(actual["Signs"]-expected["Signs"])/.rules]===ConstantArray[0,14]],{40}]],True,TestID->"arbitrary symbolic signs"],
 VerificationTest[Module[{op,ops,out,u,v},op=QuantumOperator[DiagonalMatrix[{Exp[I alpha],Exp[I beta]}]];ops={"H"->1,op,QuantumOperator[op,{2}],"CNOT"->{1,2},"H"->2};out=qc[ops][PauliStabilizer[2]];u=dense[ops,2];v=Normal[out["State"]["StateVector"]];And@@Flatten[Table[With[{rules={alpha->aa,beta->bb}},Chop[N[(Outer[Times,v,Conjugate[v]]-Outer[Times,u[[All,1]],Conjugate[u[[All,1]]]])/.rules]]===ConstantArray[0,{4,4}]],{aa,{0,Pi/8,2Pi/7}},{bb,{0,-Pi/6,Pi}}]]],True,TestID->"symbolic phase repeated on moved local wires"],
 VerificationTest[With[{p=QuantumOperator["P"[ArcCos[t]]]},With[{ops={p,QuantumOperator[p,{3}],QuantumOperator[p,{2}]}},stabilizerCircuitSpecs[qc[ops]]===Catenate[stabilizerOperatorSpecs /@ ops]]],True,TestID->"occurrence symbolic token fallback"],
 VerificationTest[qc[{QuantumOperator["T"]@QuantumOperator["H"],QuantumOperator[QuantumOperator["T"]@QuantumOperator["H"],"Label"->"H"]}][PauliStabilizer[1]],$Failed,{PauliStabilizer::nonclifford},TestID->"composite misleading label refused"],
 VerificationTest[qc[{QuantumOperator[PauliMatrix[1],{{1},{1}}],QuantumOperator[PauliMatrix[1],{{2},{1}}]}][PauliStabilizer[2]],$Failed,{PauliStabilizer::nonclifford},TestID->"moved wire shares data refused"],
 VerificationTest[qc[{QuantumOperator["H"],QuantumOperator[QuantumOperator["H"],QuantumBasis["PauliY"]]}][PauliStabilizer[1]],$Failed,{PauliStabilizer::nonclifford},TestID->"noncomputational basis refused"],
 VerificationTest[qc[{QuantumOperator[IdentityMatrix[4],QuantumBasis[4]],QuantumOperator["H"]}][PauliStabilizer[1]],$Failed,{PauliStabilizer::nonclifford},TestID->"one ququart refused"],
 VerificationTest[qc[{"H"->1,QuantumOperator[2IdentityMatrix[2]],"S"->1}][PauliStabilizer[1]],$Failed,{PauliStabilizer::nonclifford},TestID->"nonunitary refused"],
 VerificationTest[With[{s=randomSpecs[3,20]},eq[qc[{qc[Take[s,10]],"Barrier",qc[Drop[s,10]],qc[{}]}][PauliStabilizer[3]],Fold[#1[#2]&,PauliStabilizer[3],s]]],True,TestID->"nested mixed empty"],
 VerificationTest[qc[{}][PauliStabilizer[63]]===PauliStabilizer[63],True,TestID->"empty packed identity"],
 VerificationTest[And@@Table[With[{s=randomSpecs[n,120],init=PauliStabilizer[n]["ApplyCircuit",{"H"->n,"Y"->1}]},eq[qc[s][init],Fold[#1[#2]&,PauliStabilizer[<|"Tableau"->init["Tableau"],"Signs"->init["Signs"]|>],s]]],{n,{31,32,61,62,63,124,125,1024}}],True,TestID->"packed boundary scalar independent"],
 VerificationTest[With[{s=randomSpecs[128,300]},eq[Nest[qc[s],PauliStabilizer[128],40],PauliStabilizer[128]["ApplyCircuit",Catenate[ConstantArray[s,40]]]]],True,TestID->"twelve thousand reused block"],
 VerificationTest[Module[{report},report=TestReport["/Users/mohammadb/Documents/GitHub/QuantumFramework/Tests/Stabilizer/BulkPlanBuilder.wlt"];{report["TestsSucceededCount"],report["TestsFailedCount"]}],{56,0},TestID->"new suite independently rerun"]
};
report=TestReport[tests];
Echo[<|"Succeeded"->report["TestsSucceededCount"],"Failed"->report["TestsFailedCount"],"Results"->KeyValueMap[<|"ID"->#2["TestID"],"Outcome"->#2["Outcome"],"Actual"->ToString[#2["ActualOutput"],InputForm],"Expected"->ToString[#2["ExpectedOutput"],InputForm],"Messages"->ToString[#2["ActualMessages"],InputForm]|>&,report["TestResults"]]|>,"Review"];
Echo[<|"CompiledHead"->ToString[Head[compiledCliffordFold],InputForm],"CompiledTarget"->$stabilizerCompileTarget,"KernelMessages"->ToString[$MessageList,InputForm]|>,"Final"];
Put[report,"/private/tmp/qf-bulk-plan-review-r2-"<>$stabilizerCompileTarget<>".wl"];
Exit[If[report["TestsFailedCount"]===0,0,1]];
