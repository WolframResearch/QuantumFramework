(* Isolated, unmodified helpers loaded from the pinned QI implementation. *)
Begin["PermanentSelectionAudit`"];
$SourcePath = "/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/QuantumFramework/Kernel/SecondQuantization/Interferometer.m";
Get[FileNameJoin[{DirectoryName[$InputFileName],"OriginalHelpers.wl"}]];
DirectAmplitude[w_, s_, t_] := Which[
 Total[s] != Total[t], 0,
 Total[s] == 0, 1,
 True, Permanent[w[[modeList[t], modeList[s]]]]/Sqrt[Times@@(s!) Times@@(t!)]
];
RepeatedAmplitude[w_, s_, t_] := Which[
 Total[s] != Total[t], 0,
 Total[s] == 0, 1,
 True, With[{v = If[Precision[w] === Infinity, w, SetPrecision[w, Precision[w] + Total[s]]]},
  N[If[Times@@(s+1) <= Times@@(t+1), repeatedPermanent[v,s,t], repeatedPermanent[Transpose[v],t,s]]/
   Sqrt[Times@@(s!) Times@@(t!)], Precision[w]]]
];
Features[w_, s_, t_] := With[{ks=Times@@(s+1), kt=Times@@(t+1)},
 <|"Photons"->Total[s], "Modes"->Length[s],
   "Precision"->If[Precision[w]===MachinePrecision,"Machine",ToString[Round[Precision[w]]]],
   "Terms"->Min[ks,kt], "Rows"->If[ks<=kt,Count[t,_?Positive],Count[s,_?Positive]],
   "InputOccupied"->Count[s,_?Positive], "OutputOccupied"->Count[t,_?Positive],
   "Direction"->If[ks<=kt,"Input","Output"]|>
];
End[];
