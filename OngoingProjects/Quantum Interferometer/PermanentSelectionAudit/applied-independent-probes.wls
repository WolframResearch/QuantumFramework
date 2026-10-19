$AllowInternet=False;
PacletDirectoryLoad /@ {
 "/Users/mohammadb/Library/Wolfram/Paclets/Repository/Wolfram__Arrays-1.4.1",
 "/Users/mohammadb/Library/Wolfram/Paclets/Repository/Wolfram__TensorNetworks-1.2.0",
 "/Users/mohammadb/Library/Wolfram/Paclets/Repository/IBMQuantumPlatform-0.0.4"
};
PacletDirectoryLoad["/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
Echo[{$Version,$SystemID,FindFile["Wolfram`QuantumFramework`"]},"ENVIRONMENT"];
$qiLoadMessages = $MessageList;
$qiAuditDir = "/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/OngoingProjects/Quantum Interferometer/PermanentSelectionAudit";
Begin["QIOriginalVerifier`"];
Get["/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/OngoingProjects/Quantum Interferometer/PermanentSelectionAudit/OriginalHelpers.wl"];
End[];
Get[FileNameJoin[{$qiAuditDir,"Proposal.wl"}]];
$qiApplied = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`permanentAmplitude;
$qiPredicate = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`useRepeatedPermanentQ;
Echo[With[{f=$qiApplied}, {Length[DownValues[f]],FreeQ[DownValues[f],PermanentSelectionAudit`ProposedAmplitude],!FreeQ[DownValues[f],Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`useRepeatedPermanentQ]}],"ACTUAL_DISPATCH"];
$qiRows=Import[FileNameJoin[{$qiAuditDir,"benchmark.json"}],"RawJSON"];
$qiCompare=Map[Function[row,With[{w=If[row["Precision"]==="Machine",N[IdentityMatrix[row["Modes"]]],N[IdentityMatrix[row["Modes"]],ToExpression[row["Precision"]]]],s=row["Input"],t=row["Output"]},$qiPredicate[w,s,t] === (PermanentSelectionAudit`SelectionPlan[w,s,t] === "Repeated")]],$qiRows];
Echo[{Length[$qiCompare],Counts[$qiCompare]},"PREDICATE_DATASET_EQUIVALENCE"];
Echo[Trace[$qiPredicate[{{a,b},{c,d}},{1,1},{1,1}],_Count|_Total|_Times],"EXACT_SHORT_CIRCUIT"];
Echo[{$qiApplied[{}, {}, {}],$qiApplied[IdentityMatrix[2],{1,0},{0,0}],$qiApplied[{{a,b},{c,d}},{1,1},{1,1}]},"EDGE_RESULTS"];
$qiDefs=With[{f=$qiPredicate},DownValues[f]];
Echo[$qiDefs,"PRODUCTION_PREDICATE_DEFINITION"];
$qiCases={
 {"TinyFree",16,2,"Machine",{1,1}},
 {"MildCollision",16,12,"Machine",Join[{2},ConstantArray[1,10]]},
 {"HighModeBalanced",128,12,"Machine",{6,6}},
 {"StrongCollision",16,16,"Machine",{8,8}},
 {"ThirtyFiveDigitBalanced",16,8,35,{4,4}},
 {"ThirtyFiveDigitMild",24,8,35,Join[{2},ConstantArray[1,6]]},
 {"OutputCollision",16,12,"Machine",ConstantArray[1,12]},
 {"TransposedCollision",16,12,"Machine",{12}},
 {"HighModeSingle",256,16,"Machine",{16}}
};
$qiOutputs=Map[Function[case,With[{label=case[[1]],m=case[[2]],n=case[[3]],prec=case[[4]],s=PadRight[case[[5]],case[[2]]]},
 Module[{t,wref,w,old,actual,reference,abserr,relerr,newTimes,oldTimes},
 t=Switch[label,"OutputCollision",PadLeft[{n},m],"TransposedCollision",PadLeft[ConstantArray[1,n],m],"HighModeSingle",s,_,Reverse[s]];
 wref=N[FourierMatrix[m],85];w=If[prec==="Machine",N[wref],N[wref,prec]];
 actual=$qiApplied[w,s,t];old=QIOriginalVerifier`permanentAmplitude[w,s,t];
 reference=If[label==="HighModeSingle",wref[[1,1]]^n,QIOriginalVerifier`permanentAmplitude[wref,s,t]];
 abserr=N[Abs[Rationalize[actual,0]-Rationalize[reference,0]],20];
 relerr=If[Abs[reference]>10^-70,N[abserr/Abs[reference],20],Null];
 newTimes=Table[ClearSystemCache[];First[AbsoluteTiming[$qiApplied[w,s,t]]],{5}];
 oldTimes=Table[ClearSystemCache[];First[AbsoluteTiming[QIOriginalVerifier`permanentAmplitude[w,s,t]]],{5}];
 <|"Case"->label,"Modes"->m,"Photons"->n,"Precision"->prec,"Repeated"->$qiPredicate[w,s,t],"OldNewAbsoluteDifference"->N[Abs[Rationalize[actual,0]-Rationalize[old,0]],20],"ReferencePrecision"->Precision[reference],"AbsoluteError"->abserr,"RelativeError"->relerr,"Pass"->TrueQ[abserr<If[prec==="Machine",10^-12,10^-28]],"NewSeconds"->newTimes,"OldSeconds"->oldTimes,"NewMedian"->Median[newTimes],"OldMedian"->Median[oldTimes]|>
 ]]],$qiCases];
Export["/tmp/qi-applied-verifier-probes.json",$qiOutputs,"RawJSON"];
Echo[$qiOutputs,"NUMERICAL_TIMING_PROBES"];
Echo[{Counts[Lookup[$qiOutputs,"Pass"]],Drop[$MessageList,UpTo[Length[$qiLoadMessages]]]},"PROBE_SUMMARY"];
Echo[With[{f=$qiApplied}, {FreeQ[DownValues[f],PermanentSelectionAudit`ProposedAmplitude],!FreeQ[DownValues[f],Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`useRepeatedPermanentQ]}],"FINAL_ACTUAL_DISPATCH"];
