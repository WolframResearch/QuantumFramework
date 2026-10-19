Get[FileNameJoin[{DirectoryName[$InputFileName],"Proposal.wl"}]];
BeginTestSection["Permanent selection proposal"];
VerificationTest[
 With[{w=FourierMatrix[3]}, PermanentSelectionAudit`ProposedAmplitude[w,{2,1,0},{1,2,0}] === PermanentSelectionAudit`permanentAmplitude[w,{2,1,0},{1,2,0}]],
 True,TestID->"Proposal-preserves-exact-Fourier-amplitude"];
VerificationTest[
 {PermanentSelectionAudit`ProposedAmplitude[IdentityMatrix[2],{0,0},{0,0}],PermanentSelectionAudit`ProposedAmplitude[IdentityMatrix[2],{1,0},{0,0}]},
 {1,0},TestID->"Proposal-vacuum-and-number-conservation"];
VerificationTest[
 With[{w={{1,-1},{1,1}}/Sqrt[2]},
  {PermanentSelectionAudit`ProposedAmplitude[w,{1,1},{1,1}],Abs[PermanentSelectionAudit`ProposedAmplitude[w,{1,1},{2,0}]]^2,Abs[PermanentSelectionAudit`ProposedAmplitude[w,{1,1},{0,2}]]^2}],
 {0,1/2,1/2},TestID->"Proposal-exact-HOM"];
VerificationTest[
 With[{s=PadRight[Join[{2},ConstantArray[1,10]],16],t=Reverse@PadRight[Join[{2},ConstantArray[1,10]],16]},
  PermanentSelectionAudit`SelectionPlan[N[FourierMatrix[16]],s,t]],
 "Direct",TestID->"Proposal-avoids-mild-collision-slowdown"];
VerificationTest[
 With[{s=PadRight[{8,8},16],t=Reverse@PadRight[{8,8},16]},
  PermanentSelectionAudit`SelectionPlan[N[FourierMatrix[16]],s,t]],
 "Repeated",TestID->"Proposal-retains-strong-collision-route"];
VerificationTest[
 With[{s=PadRight[{6,6},128],t=Reverse@PadRight[{6,6},128]},
  PermanentSelectionAudit`SelectionPlan[N[FourierMatrix[128]],s,t]],
 "Direct",TestID->"Proposal-accounts-for-mode-conversion-cost"];
VerificationTest[
 With[{s=PadRight[{4,4},16],t=Reverse@PadRight[{4,4},16]},
  {PermanentSelectionAudit`SelectionPlan[N[FourierMatrix[16]],s,t],PermanentSelectionAudit`SelectionPlan[N[FourierMatrix[16],35],s,t]}],
 {"Direct","Repeated"},TestID->"Proposal-distinguishes-machine-and-arbitrary-precision"];
VerificationTest[
 With[{w=N[DiagonalMatrix[{Exp[37 I/100],Exp[-I/5]}]],s={24,0}},
  {PermanentSelectionAudit`SelectionPlan[w,s,s]==="Repeated",Abs[PermanentSelectionAudit`ProposedAmplitude[w,s,s]-w[[1,1]]^24]<10^-13}],
 {True,True},TestID->"Proposal-preserves-large-bunching-accuracy"];
VerificationTest[
 With[{m=16,n=12},
  With[{w=N[FourierMatrix[m]],s=PadRight[ConstantArray[1,n],m],t=PadLeft[{n},m]},
   Abs[PermanentSelectionAudit`ProposedAmplitude[w,s,t]-Sqrt[n!] Times@@w[[-1,Range[n]]]]<10^-12]],
 True,TestID->"Proposal-output-collision-factorial-normalization"];
VerificationTest[
 With[{we=DiagonalMatrix[{1,I,-1,-I}].FourierMatrix[4].DiagonalMatrix[{1,-I,-1,I}]},
  And@@Flatten@Table[
   With[{patterns=FrobeniusSolve[ConstantArray[1,4],n]},
    And@@Flatten@Table[
     With[{a=PermanentSelectionAudit`DirectAmplitude[we,s,t]},
      Abs[PermanentSelectionAudit`ProposedAmplitude[N[we],s,t]-N[a]]<10^-12 &&
      Abs[PermanentSelectionAudit`ProposedAmplitude[N[we,35],s,t]-N[a,35]]<10^-28],
     {s,patterns},{t,patterns}]],
   {n,1,4}]],
 True,TestID->"Proposal-all-1741-small-occupation-pairs-at-two-precisions"];

VerificationTest[
 With[{w=N[FourierMatrix[3]],s={2,1,1},outs=FrobeniusSolve[{1,1,1},4]},
  Abs[Total[Abs[PermanentSelectionAudit`ProposedAmplitude[w,s,#]]^2&/@outs]-1]<10^-12],
 True,TestID->"Proposal-scalar-amplitudes-normalize-complete-sector"];
VerificationTest[
 With[{w=SparseArray[N[DiagonalMatrix[{Exp[I/3],Exp[I/5],Exp[I/7],Exp[I/11]}]]],s={4,4,0,0}},
  Abs[PermanentSelectionAudit`ProposedAmplitude[w,s,s]-w[[1,1]]^4 w[[2,2]]^4]<10^-12],
 True,TestID->"Proposal-sparse-unitary-correctness"];
VerificationTest[
 With[{we=FourierMatrix[4]},
  With[{w=ReplacePart[N[we,50],{1,1}->N[we[[1,1]]]],s={2,1,1,0},t={0,1,1,2}},
   Abs[PermanentSelectionAudit`ProposedAmplitude[w,s,t]-N[PermanentSelectionAudit`DirectAmplitude[we,s,t],50]]<10^-12]],
 True,TestID->"Proposal-mixed-precision-unitary-correctness"];
EndTestSection[];
