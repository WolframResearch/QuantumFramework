repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$ProgressReporting = False;
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
Echo[<|"Version" -> $Version, "Machine" -> $SystemID, "Loaded" -> FindFile["Wolfram`QuantumFramework`"], "Paclet" -> PacletFind["Wolfram/QuantumFramework"]|>, "ENV"];
checks = Internal`Bag[];
SetAttributes[verify, HoldRest];
verify[name_, expr_] := With[{result = expr}, Internal`StuffBag[checks, <|"Name" -> name, "Pass" -> TrueQ[result], "Result" -> ToString[result, InputForm]|>]; Echo[{name, result}, "PROBE"]];
same[a_, b_] := Normal[a["Tableau"]] === Normal[b["Tableau"]] && a["Signs"] === b["Signs"];
canonical[p_] := PauliStabilizer[<|"Signs" -> p["Signs"], "Tableau" -> p["Tableau"]|>];
layouts[p_] := {canonical[p], psFromPacked[psGetPacked[p]], psToGeneratorPacked[p]};
paulis = {IdentityMatrix[2], PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};
hmat = {{1,1},{1,-1}}/Sqrt[2]; smat = DiagonalMatrix[{1,I}];
vmat = {{1+I,1-I},{1-I,1+I}}/2;
tensor[ms_] := Fold[KroneckerProduct, {{1}}, ms];
embed[m_, q_, n_] := tensor[ReplacePart[ConstantArray[IdentityMatrix[2], n], q -> m]];
gateMatrix[g_ -> wire_, n_] := With[{w = Flatten[{wire}]},
 If[MemberQ[{"CNOT","CX","CZ","SWAP"},g],
  With[{a=w[[1]],b=w[[2]], states = IntegerDigits[#,2,n]& /@ Range[0,2^n-1]},
   Which[g === "CZ", DiagonalMatrix[(-1)^(#[[a]] #[[b]])& /@ states],
    g === "SWAP", With[{out = (FromDigits[ReplacePart[#, {a -> #[[b]], b -> #[[a]]}],2]+1)& /@ states}, Normal[SparseArray[Thread[Transpose[{out, Range[2^n]}] -> 1], {2^n,2^n}]]],
    True, With[{out = (FromDigits[ReplacePart[#, b -> Mod[#[[b]]+#[[a]],2]],2]+1)& /@ states}, Normal[SparseArray[Thread[Transpose[{out, Range[2^n]}] -> 1], {2^n,2^n}]]]]],
  embed[Replace[g,{"H"->hmat,"S"->smat,SuperDagger["S"]->ConjugateTranspose[smat],"V"->vmat,SuperDagger["V"]->ConjugateTranspose[vmat],"X"->paulis[[2]],"Y"->paulis[[3]],"Z"->paulis[[4]]}], First[w], n]]];
rowMatrix[p_, j_] := With[{t=p["Tableau"]}, p["Signs"][[j]] tensor[Map[paulis[[#]] &, Replace[Transpose[{t[[1,All,j]],t[[2,All,j]]}], {{0,0}->1,{1,0}->2,{1,1}->3,{0,1}->4},{1}]]]];
wordKey[pair_] := Flatten[Simplify /@ {pair[[1]].PauliMatrix[1].ConjugateTranspose[pair[[1]]],pair[[1]].PauliMatrix[3].ConjugateTranspose[pair[[1]]]}];
cliffords = FixedPoint[DeleteDuplicatesBy[Join[#, Catenate[Function[p, {{hmat.p[[1]], Append[p[[2]], "H"->1]}, {smat.p[[1]], Append[p[[2]],"S"->1]}}] /@ #]], wordKey]&, {{IdentityMatrix[2],{}}}];
verify["DenseEnumerates24IndependentCliffords", Length[cliffords] === 24];
verify["Dense24ForwardInverseAllStorage", And @@ Map[Function[pair, With[{u=pair[[1]], p=PauliStabilizer[1]["ApplyCircuit",pair[[2]]]}, And @@ Map[Function[obj,
 With[{inv=obj["Dagger"]}, same[inv["Dagger"],obj] && inv["Dagger"] === obj && same[obj[inv],PauliStabilizer[1]] && same[inv[obj],PauliStabilizer[1]] &&
 And @@ Table[Simplify[rowMatrix[obj,j]-u.paulis[[If[j==1,2,4]]].ConjugateTranspose[u]] === ConstantArray[0,{2,2}] && Simplify[rowMatrix[inv,j]-ConjugateTranspose[u].paulis[[If[j==1,2,4]]].u] === ConstantArray[0,{2,2}],{j,2}]]], layouts[p]]]], cliffords]];
SeedRandom[48172831];
makeSpecs[n_, count_] := Table[With[{w=If[n==1,{1,1},RandomSample[Range[n],2]],g=RandomChoice[If[n==1,{"H","S",SuperDagger["S"],"V",SuperDagger["V"],"X","Y","Z"},{"H","S",SuperDagger["S"],"V",SuperDagger["V"],"X","Y","Z","CNOT","CZ","SWAP"}]]},g -> If[MemberQ[{"CNOT","CZ","SWAP"},g],w,First[w]]],{count}];
Scan[Function[n, verify["DenseEntanglingAllGeneratorAndInverse-n"<>ToString[n], And @@ Table[With[{specs=makeSpecs[n,27]}, With[{p=PauliStabilizer[n]["ApplyCircuit",specs], u=Fold[gateMatrix[#2,n].#1&,IdentityMatrix[2^n],specs]},
 Simplify[ConjugateTranspose[u].u]===IdentityMatrix[2^n] && And @@ Map[Function[obj, With[{inv=obj["Dagger"]}, And @@ Table[With[{base=embed[PauliMatrix[If[j<=n,1,3]],Mod[j-1,n]+1,n]}, Simplify[rowMatrix[obj,j]-u.base.ConjugateTranspose[u]]===ConstantArray[0,{2^n,2^n}] && Simplify[rowMatrix[inv,j]-ConjugateTranspose[u].base.u]===ConstantArray[0,{2^n,2^n}]],{j,2n}]]],layouts[p]]]],{4}]]],{2,3,4}];
verify["PhysicalVAndAdjoint", And @@ Table[With[{u=If[dag,ConjugateTranspose[vmat],vmat],g=If[dag,SuperDagger["V"],"V"]}, And @@ Table[With[{p=PauliStabilizer[1][g,1]}, Simplify[rowMatrix[p,j]-u.PauliMatrix[If[j==1,1,3]].ConjugateTranspose[u]]===ConstantArray[0,{2,2}]],{j,2}]],{dag,{False,True}}]];
(* Independent scalar AG recurrence works on ordinary unpacked arrays. *)
refApply[{x_,z_,s_}, g_ -> w_] := Block[{xx=x,zz=z,ss=s,a=First[Flatten[{w}]], b=Last[Flatten[{w}]], xa,za,xb,zb},
 xa=x[[a]]; za=z[[a]];
 Switch[g,
 "H", ss=s(1-2xa za); xx[[a]]=za;zz[[a]]=xa,
 "S",ss=s(1-2xa za);zz[[a]]=Mod[xa+za,2],
 SuperDagger["S"],ss=s(1-2xa(1-za));zz[[a]]=Mod[xa+za,2],
 "X",ss=s(1-2za),"Y",ss=s(1-2Mod[xa+za,2]),"Z",ss=s(1-2xa),
 "CNOT", xb=x[[b]];zb=z[[b]];ss=s(1-2xa zb(1-Mod[xb+za,2]));xx[[b]]=Mod[xb+xa,2];zz[[a]]=Mod[za+zb,2],
 "SWAP",xx[[{a,b}]]=xx[[{b,a}]];zz[[{a,b}]]=zz[[{b,a}]],
 "CZ", Return[Fold[refApply,{x,z,s},{"H"->b,"CNOT"->{a,b},"H"->b}]],
 "V",Return[Fold[refApply,{x,z,s},{"H"->a,"S"->a,"H"->a}]],
 SuperDagger["V"],Return[Fold[refApply,{x,z,s},{"H"->a,SuperDagger["S"]->a,"H"->a}]]];{xx,zz,ss}];
Scan[Function[n, verify["ScalarBoundaryAndSplit-n"<>ToString[n], With[{base=PauliStabilizer[n], specs=makeSpecs[n,80]}, With[{p=base["ApplyCircuit",specs], r=Fold[refApply,Join[base["Tableau"],{base["Signs"]}],specs]}, p["Tableau"]===r[[1;;2]] && p["Signs"]===r[[3]] && PauliStabilizerQ[p] && same[p,Fold[#1["ApplyCircuit",#2]&,base,Partition[specs,10]]] && same[p,psFromPacked[psGetPacked[p]]] && !KeyExistsQ[First[p],"Tableau"]]]]],{1,2,5,31,32,61,62,63,124,125,1024}];
verify["SymbolicBinaryPhasesAll16Assignments", Module[{p, specs, vars, symbolic, evolved, inv},
 p = PauliStabilizer[2]["ApplyCircuit",makeSpecs[2,19]];
 specs = makeSpecs[2,19]; vars = Array[phaseBit,4];
 symbolic = PauliStabilizer[<|"Matrix"->p["Matrix"],"Signs"->1-2vars|>];
 evolved = symbolic["ApplyCircuit",specs]; inv = symbolic["Dagger"];
 And @@ Table[With[{rules = Thread[vars -> bits]}, With[{con = canonical[symbolic /. rules]},
  same[evolved /. rules, con["ApplyCircuit", specs]] && same[inv /. rules, con["Dagger"]]
 ]], {bits, Tuples[{0,1},4]}]
]];
verify["PropertyEquivalenceAcrossAllStorage", With[{p=PauliStabilizer[4]["ApplyCircuit",makeSpecs[4,43]]},With[{ls=layouts[p]}, And @@ Map[Function[prop, SameQ @@ (#[prop]& /@ ls)], {"Qubits","GeneratorCount","Signs","Phase","X","Z","Matrix","TableauPhase","StabilizerX","StabilizerZ","DestabilizerSigns","PauliSymbols"}] && SameQ @@ (#["Entropy",{1,3}]& /@ ls) && SameQ @@ (#["Expectation","YIXZ"]& /@ ls) && SameQ @@ (Keys[#["M",2]]& /@ ls) && And @@ Map[Function[q,With[{v=Normal[q["State"]["StateVector"]]},Simplify[Conjugate[v].v]===1]],ls]]]];
verify["ImmutableValidationNoReusedBit", With[{p=PauliStabilizer[32]}, Block[{bad=p},bad[[1,"Signs"]]={};PauliStabilizerQ[p] && !PauliStabilizerQ[bad] && !psConcreteFastQ[bad]]]];
verify["MalformedAllStorageRejected", And @@ Map[Function[n,With[{p=PauliStabilizer[n],c=psFromPacked[psGetPacked[PauliStabilizer[n]]]}, And @@ (Not[PauliStabilizerQ[#]]& /@ {PauliStabilizer[KeyDrop[First[p],"GeneratorZ"]],PauliStabilizer[Append[First[p],"Signs"->{}]],PauliStabilizer[Append[First[p],"GeneratorX"->{{-1}}]],PauliStabilizer[Append[First[c],"PackedZ"->{{2^62}}]]})]],{1,32,62,63}]];
verify["EmptyIdentity", With[{p=PauliStabilizer[5]},p["ApplyCircuit",{}]===p]];
verify["ExtremeIntegerEncodingIsExact", encodeStabilizerGates[{"H"->2^63,"CNOT"->{-2^63,2^63}}] === {{0,2^63-1,0},{3,-2^63-1,2^63-1}}];
Echo[<|"OptionsApply"->Options[PauliStabilizerApply],"Contexts"->(Context/@{PauliStabilizerApply,psFromGeneratorPacked,applyCompiledFold,stabilizerExecutionPlan}),"CompileTarget"->$stabilizerCompileTarget,"Properties"->PauliStabilizer[1]["Properties"]|>,"SIGNATURES"];
reports = Map[Function[file,With[{r=TestReport[FileNameJoin[{repo,"Tests","Stabilizer",file}]]},With[{results=r["Results"]},With[{row=<|"File"->file,"Count"->Length[results],"Outcomes"->Counts[#["Outcome"]& /@ results],"Failures"->Map[Function[z,<|"ID"->ToString[z["TestID"],InputForm],"Expected"->ToString[z["ExpectedOutput"],InputForm],"Actual"->ToString[z["ActualOutput"],InputForm],"Messages"->ToString[z["ActualMessages"],InputForm]|>],Select[results,#["Outcome"]=!="Success"&]]|>},Echo[row,"FOCUSED"];row]]]],{"CompiledEvolution.wlt","InversePhase.wlt","ExecutionContracts.wlt"}];
rows = Internal`BagPart[checks,All];
Export["/private/tmp/stabilizer-style-r2-independent.json",<|"Version"->$Version,"LoadedFile"->FindFile["Wolfram`QuantumFramework`"],"Checks"->rows,"Focused"->reports|>,"RawJSON"];
Echo[<|"Checks"->Length[rows],"Failures"->Select[rows,!TrueQ[#Pass]&],"FocusedCount"->Total[Lookup[reports,"Count"]],"FocusedFailures"->Total[Length/@Lookup[reports,"Failures"]]|>,"FINAL"];
Exit[If[AllTrue[rows,TrueQ[#Pass]&] && AllTrue[reports,#Failures==={}&],0,1]];
