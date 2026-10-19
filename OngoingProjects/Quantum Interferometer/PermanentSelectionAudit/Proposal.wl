(* Prototype only: no QI kernel source is modified. Coefficients are calibrated
   against benchmark.json training rows on WL 15.0, MacOSX-ARM64.
   K is the smaller grouped Ryser term count; r is its distinct output-row count.
   m^2 represents promotion of the whole transfer matrix; K m (r+1) represents
   the existing interpreted summation, including its zero-occupation dimensions.
   The model compares weighted work scores, not guaranteed operation counts. *)
Get[FileNameJoin[{DirectoryName[$InputFileName],"SourceHelpers.wl"}]];
Begin["PermanentSelectionAudit`"];
SelectionPlan[w_,s_,t_] := With[{n=Total[s],m=Length[s],ks=Times@@(s+1),kt=Times@@(t+1),p=Precision[w]},
 Which[
  n!=Total[t],"Zero",
  n==0,"Vacuum",
  p===Infinity,"Direct",
  True,With[{k=Min[ks,kt],r=If[ks<=kt,Count[t,_?Positive],Count[s,_?Positive]],factor=If[p===MachinePrecision,16,1/10]},
   If[k<2^n && factor (m^2+k m (r+1))<n 2^n,"Repeated","Direct"]]
 ]
];
ProposedAmplitude[w_,s_,t_] := Switch[SelectionPlan[w,s,t],
 "Zero",0,
 "Vacuum",1,
 "Repeated",RepeatedAmplitude[w,s,t],
 "Direct",DirectAmplitude[w,s,t]
];
(* Fresh-kernel integration hook. Changes only the in-memory private helper,
   and returns its original DownValues so a caller can restore them. *)
InstallProposal[] := With[{sym=ToExpression[First@Names["Wolfram`QuantumFramework`SecondQuantization`*`permanentAmplitude"]]},
 With[{saved=DownValues[sym]},
  Clear[sym];
  sym[w_,s_,t_] := ProposedAmplitude[w,s,t];
  saved
 ]
];
End[];
