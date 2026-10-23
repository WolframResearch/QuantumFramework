root=DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{root,"optimization","WidthOptimization.wl"}]];
source=FileNameJoin[{root,"feasibility","B-quadratic-192.json"}];
rows=Import[source,"RawJSON"]["table"];
saved=Import[FileNameJoin[{root,"optimization","width-certificates.json"}],"RawJSON"];
selected=InvSqrtWidths`CertifyWidths[rows,{30,23,16},15/(16 2^24)];
savedSelected=SelectFirst[saved["certificates"],#["fractional_bits"]=={30,23,16}&];
tests={
 VerificationTest[InvSqrtWidths`HornerError[t,e0,e1],e0-t e1,TestID->"signed-error"],
 VerificationTest[Resolve[ForAll[{t,e0,e1,b0,b1},
   t>=0 && b0>=0 && b1>=0 && 0<=e0<=b0 && 0<=e1<=b1,
   -t b1<=InvSqrtWidths`HornerError[t,e0,e1]<=b0],Reals],True,TestID->"general-error-envelope"],
 VerificationTest[Resolve[ForAll[{u,k},k<=u<k+1,0<=u-k<1],Reals],True,TestID->"floor-row-residual"],
 VerificationTest[FullSimplify[Sum[2^-f,{j,1,n}],Element[{f,n},Integers] && f>0 && n>=0],n/2^f,TestID->"arbitrary-row-count"],
 VerificationTest[InvSqrtWidths`HornerError[t,0,0],0,TestID->"exact-arithmetic-limit"],
 VerificationTest[InvSqrtWidths`HornerError[0,e0,e1],e0,TestID->"interval-origin"],
 VerificationTest[FullSimplify[Limit[n/2^f,f->Infinity],n>0],0,TestID->"precision-limit"],
 VerificationTest[FullSimplify[(y 2^(s-10)) Sqrt[x]==y Sqrt[x 4^s/2^20],
   x>0 && Element[{s,y},Reals]],True,TestID->"normalization-relative-error"],
 VerificationTest[FullSimplify[Abs[Exp[2 Pi I phi]]^2,Element[phi,Reals]],1,TestID->"unit-phase-modulus"],
 VerificationTest[selected["certified"],True,TestID->"all-intervals-and-capacities"],
 VerificationTest[selected,savedSelected,TestID->"reproduce-stored-certificate"],
 VerificationTest[selected["register_bits"],{30,22,15},TestID->"leading-zero-removal"],
 VerificationTest[Dimensions[selected["table"]],{192,3},TestID->"nontrivial-interval-count"],
 VerificationTest[Round[{1/2,3/2,5/2,7/2}],{0,2,2,4},TestID->"round-ties-to-even"],
 VerificationTest[FailureQ[InvSqrtWidths`QuantizedRows[{}, {30,23,16}]],True,TestID->"empty-table"],
 VerificationTest[FailureQ[InvSqrtWidths`CertifyWidths[rows,{0,23,16},1/2]],True,TestID->"invalid-width"],
 VerificationTest[FailureQ[InvSqrtWidths`CertifyWidths[rows,{30,23,16},0.01]],True,TestID->"reject-inexact-tolerance"],
 VerificationTest[FailureQ[InvSqrtWidths`CertifyWidths[rows,{30,23,16},1]],True,TestID->"reject-unit-tolerance"],
 VerificationTest[Count[Lookup[saved["certificates"],"certified"],False],4,TestID->"rejected-aggressive-designs"],
 VerificationTest[FileHash[source,"SHA256","HexString"],saved["source_sha256"],TestID->"frozen-source-hash"],
 VerificationTest[FileHash[FileNameJoin[{root,"optimization","WidthOptimization.wl"}],"SHA256","HexString"],
   saved["certificate_source_sha256"],TestID->"certificate-code-hash"]
};
report=TestReport[tests];
summary=<|"kernel"->$Version,"tests"->Length[tests],
  "succeeded"->report["TestsSucceededCount"],"failed"->report["TestsFailedCount"],
  "verifier_source_sha256"->FileHash[$InputFileName,"SHA256","HexString"]|>;
Export[FileNameJoin[{root,"optimization","wl-checks.json"}],summary,"RawJSON"];
Echo[summary];
If[report["TestsFailedCount"]!=0,Exit[1]];
