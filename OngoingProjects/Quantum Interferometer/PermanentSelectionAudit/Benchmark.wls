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
$Rows={};
$Specs=Join[
 Flatten[Table[{seed,n,m,"Machine"},{seed,{1101,2202}},{n,{2,4,8,12,16,20}},{m,{Max[16,n],64,128}}],2],
 Flatten[Table[{3303,n,m,p},{n,{4,8,12,16}},{m,{16,64}},{p,{30,60}}],2]
];
$Started=AbsoluteTime[];
Echo[{$Version,Length[$Specs],$BenchmarkLimit},"BENCHMARK_START"];
Scan[Function[spec,
 With[{seed=spec[[1]],n=spec[[2]],m=spec[[3]],p=spec[[4]]},
  With[{w=MatrixFor[m,p,seed+m+n]},
   Scan[Function[layout,
    With[{s=layout[[2]],t=layout[[3]]},
     With[{direct=Measure[DirectAmplitude[w,s,t]],repeated=Measure[RepeatedAmplitude[w,s,t]]},
      With[{row=Join[<|"Seed"->seed,"Split"->If[seed===1101||seed===3303,"Training","Validation"],"Layout"->layout[[1]],"Input"->s,"Output"->t|>,Features[w,s,t],
       <|"DirectStatus"->direct["Status"],"RepeatedStatus"->repeated["Status"],"DirectTimes"->direct["Times"],"RepeatedTimes"->repeated["Times"],
         "DirectSeconds"->direct["Median"],"RepeatedSeconds"->repeated["Median"],
         "AbsoluteError"->If[NumericQ[direct["Value"]]&&NumericQ[repeated["Value"]],N[Abs[direct["Value"]-repeated["Value"]]],Null],
         "RelativeError"->If[NumericQ[direct["Value"]]&&NumericQ[repeated["Value"]],N[Abs[direct["Value"]-repeated["Value"]]/Max[Abs[repeated["Value"]],10^-100]],Null]
       |>]},
       $Rows=Append[$Rows,row];
       If[Mod[Length[$Rows],16]===0,
        Export[FileNameJoin[{$AuditDir,"benchmark.json"}],$Rows,"RawJSON"];
        Echo[{Length[$Rows],Round[AbsoluteTime[]-$Started,.1],n,m,p},"PROGRESS"]]
      ]
     ]
    ]
   ],Layouts[n,m,seed+n+m]]
  ]
 ]
],$Specs];
Export[FileNameJoin[{$AuditDir,"benchmark.json"}],$Rows,"RawJSON"];
$Columns={"Seed","Split","Layout","Photons","Modes","Precision","Terms","Rows","InputOccupied","OutputOccupied","Direction","DirectStatus","RepeatedStatus","DirectSeconds","RepeatedSeconds","AbsoluteError","RelativeError"};
(* CSV is generated from benchmark.json by Analyze.py, avoiding kernel formatting overhead. *)
Export[FileNameJoin[{$AuditDir,"metadata.json"}],<|"Version"->$Version,"SystemID"->$SystemID,"SourcePath"->$SourcePath,"SourceSHA256"->IntegerString[FileHash[$SourcePath,"SHA256"],16,64],"CaseCount"->Length[$Rows],"Seconds"->N[AbsoluteTime[]-$Started],"TimeoutPerEvaluation"->$BenchmarkLimit,"Trials"->3,"Measure"->"Normalized amplitudes including repeated-row/column matrix extraction, normalization, and whole-W precision promotion in repeated route"|>,"RawJSON"];
Echo[{Length[$Rows],AbsoluteTime[]-$Started},"BENCHMARK_DONE"];
End[];
Quit[];
