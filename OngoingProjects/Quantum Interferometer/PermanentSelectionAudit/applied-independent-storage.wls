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

SeedRandom[9911];
$qiStorageMatrices = With[{v=N[N[FourierMatrix[128],85]]},
 { {"FourierUnpacked",v}, {"FourierPacked",Developer`ToPackedArray[v]},
   {"RandomPacked",RandomVariate[CircularUnitaryMatrixDistribution[128]]} }];
$qiStorageResults=Map[Function[item,With[{label=item[[1]],w=item[[2]],s=PadRight[{6,6},128]},With[{t=Reverse[s]},
 Module[{trials,actual,old},actual=$qiApplied[w,s,t];old=QIOriginalVerifier`permanentAmplitude[w,s,t];
 trials=Table[ClearSystemCache[];With[{new=First[AbsoluteTiming[$qiApplied[w,s,t]]]},ClearSystemCache[];{new,First[AbsoluteTiming[QIOriginalVerifier`permanentAmplitude[w,s,t]]]}],{7}];
 <|"Case"->label,"Packed"->Developer`PackedArrayQ[w],"Precision"->ToString[Precision[w],InputForm],"Repeated"->$qiPredicate[w,s,t],"UnitarityError"->Max[Abs[w.ConjugateTranspose[w]-IdentityMatrix[128]]],"OldNewAbsoluteDifference"->Abs[actual-old],"NewSeconds"->trials[[All,1]],"OldSeconds"->trials[[All,2]],"NewMedian"->Median[trials[[All,1]]],"OldMedian"->Median[trials[[All,2]]]|>
 ]]]],$qiStorageMatrices];
Export["/tmp/qi-applied-verifier-storage.json",$qiStorageResults,"RawJSON"];
Echo[$qiStorageResults,"STORAGE_TIMING_PROBES"];
Echo[Drop[$MessageList,UpTo[Length[$qiLoadMessages]]],"LATE_MESSAGES"];
