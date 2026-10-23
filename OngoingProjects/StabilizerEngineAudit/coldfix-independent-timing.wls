repo=Environment["VERIFY_REPO"]; phase=Environment["VERIFY_PHASE"];
PacletDirectoryLoad[FileNameJoin[{repo,"QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
SeedRandom[20261022];
rows=Internal`Bag[];
measure[label_,f_,warm_:True] := Module[{trials,row},
 If[warm,f[]];trials=Table[ClearSystemCache[];AbsoluteTiming[Head[f[]]][[1]],{7}];
 row=<|"Label"->label,"Seconds"->trials,"MedianMs"->1000 Median[trials],"MinMs"->1000 Min[trials],"MaxMs"->1000 Max[trials]|>;
 Echo[row];Internal`StuffBag[rows,row]
];
makeSpecs[n_,m_] := Table[With[{j=Mod[37 i,n]+1},Switch[Mod[i,7],0,"H"->j,1,"S"->j,2,"CNOT"->{j,Mod[j,n]+1},3,"CZ"->{j,Mod[j,n]+1},4,"X"->j,5,"SWAP"->{j,Mod[j,n]+1},6,SuperDagger["S"]->j]],{i,m}];
makeCircuit[specs_] := QuantumCircuitOperator[<|"Elements"->(FromCircuitOperatorShorthand/@specs),"Label"->None|>];
PauliStabilizer[2]["ApplyCircuit",{"H"->1}];
Scan[Function[n,With[{p=PauliStabilizer[n],specs=Catenate[ConstantArray[makeSpecs[n,1000],10]]},
 If[phase=!="original",With[{base=makeSpecs[n,1000]},With[{qc=makeCircuit[base]},
 measure["WarmDirect1000/n="<>ToString[n],Function[p["ApplyCircuit",base]]];
 measure["WarmPublic1000/n="<>ToString[n],Function[qc[p]]];
 measure["ColdPublicPlan1000/n="<>ToString[n],Function[stabilizerClearCaches[];qc[p]],False]
 ]]];
 measure["ColdRepeated/n="<>ToString[n],Function[stabilizerClearCaches[];p["ApplyCircuit",specs]],False];
 If[phase=!="original",measure["WarmDirect/n="<>ToString[n],Function[p["ApplyCircuit",specs]]]];
 If[phase=!="original",
 measure["WarmDirectTen/n="<>ToString[n],Function[Fold[#1["ApplyCircuit",#2]&,p,Partition[specs,1000]]]];
 With[{evolved=p["ApplyCircuit",specs]},measure["Materialize/n="<>ToString[n],Function[evolved["Tableau"]]]]
 ];
 If[n===128 && phase=!="original",With[{extended=Insert[specs,"V"->1,5000]},measure["WarmV/n=128",Function[p["ApplyCircuit",extended]]]]];
 If[phase=!="original",With[{qc=makeCircuit[makeSpecs[n,1000]]},measure["WarmPublicTen/n="<>ToString[n],Function[Fold[#2[#1]&,p,ConstantArray[qc,10]]]]]]
]],{128,1024}];
SeedRandom[87923];
With[{p=PauliStabilizer[128],specs=Table[With[{w=RandomSample[Range[128],2],g=RandomChoice[{"H","S","X","Y","Z","CNOT","CZ","SWAP"}]},g->If[MemberQ[{"CNOT","CZ","SWAP"},g],w,First[w]]],{10000}]},
 measure["ColdRandom/n=128",Function[stabilizerClearCaches[];p["ApplyCircuit",specs]],False]];
Export["/private/tmp/coldfix-independent-timing-"<>phase<>".json",<|"Version"->$Version,"Repo"->repo,"Loaded"->FindFile["Wolfram`QuantumFramework`"],"Rows"->Internal`BagPart[rows,All],"Messages"->ToString[$MessageList,InputForm]|>,"RawJSON"];
Exit[];
