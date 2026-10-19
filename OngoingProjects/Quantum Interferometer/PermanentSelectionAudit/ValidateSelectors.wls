$HistoryLength=0;
$AuditDir=DirectoryName[$InputFileName];
Get[FileNameJoin[{$AuditDir,"SourceHelpers.wl"}]];
Begin["PermanentSelectionAudit`"];
$BenchmarkLimit=1.0;
If[!NumericQ[DirectAmplitude[N[IdentityMatrix[2]],{1,0},{1,0}]] || !NumericQ[RepeatedAmplitude[N[IdentityMatrix[2]],{1,0},{1,0}]], Echo["Helper validation failed","ERROR"];Exit[2]];
$MemoryLimit=256*1024^2;
SetAttributes[Measure,HoldFirst];
Measure[expr_] := Module[{first, reps, times, value, one},
 first=MemoryConstrained[TimeConstrained[AbsoluteTiming[value=expr;],$BenchmarkLimit,{$TimedOut,Null}],$MemoryLimit,{$Aborted,Null}];
 If[!NumberQ[first[[1]]],Return[<|"Status"->If[first[[1]]===$TimedOut,"Timeout","MemoryLimit"],"Times"->{},"Median"->Null,"Value"->Null|>]];
 reps=Min[32,Max[1,Ceiling[0.005/Max[first[[1]],10^-6]]]];
 times=Table[ClearSystemCache[];
  one=MemoryConstrained[TimeConstrained[First@AbsoluteTiming[Do[value=expr,{reps}]],$BenchmarkLimit*reps,$TimedOut],$MemoryLimit,$Aborted];
  If[NumberQ[one],one/reps,Null],{3}];
 <|"Status"->If[AllTrue[times,NumberQ],"OK","Incomplete"],"Times"->times,
   "Median"->If[AllTrue[times,NumberQ],Median[times],Null],"Value"->value|>
];
PadOccupation[a_,m_] := PadRight[a,m];
Layouts[n_,m_,seed_] := BlockRandom[SeedRandom[seed];
 With[{free=ConstantArray[1,n],pair=Join[{2},ConstantArray[1,n-2]],pairs=Join[ConstantArray[2,Quotient[n,2]],ConstantArray[1,Mod[n,2]]],two={Ceiling[n/2],Floor[n/2]},
 random=Counts[RandomInteger[{1,Min[n,m]},n]]},
  DeleteDuplicatesBy[{
   {"Free",PadOccupation[free,m],Reverse@PadOccupation[free,m]},
   {"OnePair",PadOccupation[pair,m],Reverse@PadOccupation[pair,m]},
   {"Pairs",PadOccupation[pairs,m],Reverse@PadOccupation[pairs,m]},
   {"TwoBalanced",PadOccupation[two,m],Reverse@PadOccupation[two,m]},
   {"SingleToFree",PadOccupation[{n},m],Reverse@PadOccupation[free,m]},
   {"FreeToSingle",PadOccupation[free,m],Reverse@PadOccupation[{n},m]},
   {"SingleToSingle",PadOccupation[{n},m],Reverse@PadOccupation[{n},m]},
   {"Random",Lookup[random,Range[m],0],RandomSample@Lookup[random,Range[m],0]}
  },Rest]
 ]
];
MatrixFor[m_,p_,seed_] := BlockRandom[SeedRandom[seed];
 If[p==="Machine",RandomVariate[CircularUnitaryMatrixDistribution[m]],
  First@QRDecomposition[N[RandomInteger[{-10,10},{m,m}]+I RandomInteger[{-10,10},{m,m}],p]]]
];

Get[FileNameJoin[{$AuditDir,"Proposal.wl"}]];
$BenchmarkLimit=1.0;
$Cases=Join[
 Flatten[Table[{n,m,"Machine"},{n,{6,10,14,18}},{m,{24,96}}],1],
 Flatten[Table[{n,m,35},{n,{6,10,14}},{m,{24,96}}],1],
 {{2,16,"Machine"},{2,128,"Machine"},{4,16,"Machine"},{4,128,"Machine"}}
];
$NewRows={};$NewStarted=AbsoluteTime[];
Scan[Function[spec,
 With[{n=spec[[1]],m=spec[[2]],p=spec[[3]]},
  With[{w=MatrixFor[m,p,5511+m+n]},
   Scan[Function[layout,
    With[{s=layout[[2]],t=layout[[3]]},
     With[{old=Measure[permanentAmplitude[w,s,t]],proposal=Measure[ProposedAmplitude[w,s,t]],direct=Measure[DirectAmplitude[w,s,t]]},
      $NewRows=Append[$NewRows,Join[<|"Layout"->layout[[1]],"Input"->s,"Output"->t,"ProposedRoute"->SelectionPlan[w,s,t]|>,Features[w,s,t],
       <|"OriginalStatus"->old["Status"],"ProposedStatus"->proposal["Status"],"DirectStatus"->direct["Status"],
        "OriginalSeconds"->old["Median"],"ProposedSeconds"->proposal["Median"],"DirectSeconds"->direct["Median"],
        "OriginalTimes"->old["Times"],"ProposedTimes"->proposal["Times"],"DirectTimes"->direct["Times"],
        "AbsoluteDifferenceFromOriginal"->If[NumericQ[old["Value"]]&&NumericQ[proposal["Value"]],N[Abs[old["Value"]-proposal["Value"]]],Null]|>]];
      Export[FileNameJoin[{$AuditDir,"selector-validation.json"}],$NewRows,"RawJSON"];
      Echo[{Length[$NewRows],n,m,p,layout[[1]],SelectionPlan[w,s,t]},"SELECTOR_PROGRESS"]
     ]
    ]
   ],Select[Layouts[n,m,5511+n+m],MemberQ[If[n<=4,{"Free"},If[p==="Machine",{"OnePair","TwoBalanced","SingleToSingle"},{"OnePair","TwoBalanced"}]],First[#]]&]]
  ]
 ]
],$Cases];
Echo[{Length[$NewRows],AbsoluteTime[]-$NewStarted},"SELECTOR_DONE"];
End[];Quit[];
