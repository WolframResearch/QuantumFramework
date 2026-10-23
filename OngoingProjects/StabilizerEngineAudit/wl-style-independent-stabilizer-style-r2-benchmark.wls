repo = Environment["QF_REVIEW_REPO"];
phase = Environment["QF_REVIEW_LABEL"];
$ProgressReporting = False;
PacletDirectoryLoad[FileNameJoin[{repo,"QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
makeSpecs[n_,count_] := Table[With[{a=Mod[37 k,n]+1},Switch[Mod[k,7],0,"H"->a,1,"S"->a,2,"CNOT"->{a,Mod[a,n]+1},3,"CZ"->{a,Mod[a,n]+1},4,"X"->a,5,"SWAP"->{a,Mod[a,n]+1},6,SuperDagger["S"]->a]],{k,count}];
makeCircuit[s_] := QuantumCircuitOperator[<|"Elements"->(FromCircuitOperatorShorthand /@ s),"Label"->None|>];
ps128=PauliStabilizer[128]; ps1024=PauliStabilizer[1024];
s128=Flatten[ConstantArray[makeSpecs[128,1000],10],1];s1024=Flatten[ConstantArray[makeSpecs[1024,1000],10],1];
v128=Insert[s128,"V"->47,5000];block1024=makeCircuit[makeSpecs[1024,1000]];
PauliStabilizer[2]["ApplyCircuit",{"H"->1,"CNOT"->{1,2}}];
measure[name_,f_,cold_] := Module[{trials},
 stabilizerClearCaches[];f[];
 trials=Table[ClearSystemCache[];If[cold,stabilizerClearCaches[]];Block[{$MessageList={}},With[{t=AbsoluteTiming[f[]]},
 <|"Seconds"->First[t],"Head"->ToString[Head[Last[t]],InputForm],"Digest"->ToString[Hash[First[Last[t]]]],"Messages"->ToString[$MessageList,InputForm]|>]],{5}];
 With[{row=<|"Case"->name,"Seconds"->Lookup[trials,"Seconds"],"Range"->MinMax[Lookup[trials,"Seconds"]],"Median"->Median[Lookup[trials,"Seconds"]],"Trials"->trials|>},Echo[row,"TIMING"];row]
];
rows={measure["DirectCold128/10000",Function[ps128["ApplyCircuit",s128]],True],
 measure["DirectWarm1024/10000",Function[ps1024["ApplyCircuit",s1024]],False],
 measure["TenPublic1024/1000",Function[Fold[#2[#1]&,ps1024,ConstantArray[block1024,10]]],False],
 measure["V128/10001",Function[ps128["ApplyCircuit",v128]],False],
 measure["DirectCold1024/10000",Function[ps1024["ApplyCircuit",s1024]],True]};
Export["/private/tmp/stabilizer-style-r2-timing-"<>phase<>".json",<|"Version"->$Version,"Repo"->repo,"Loaded"->FindFile["Wolfram`QuantumFramework`"],"Compiler"->$stabilizerCompileTarget,"Rows"->rows|>,"RawJSON"];
Echo[<|"Complete"->phase|>];Exit[];
