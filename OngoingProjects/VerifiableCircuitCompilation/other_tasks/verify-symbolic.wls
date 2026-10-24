root=DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"FieldArithmetic.wl"}]];
Get[FileNameJoin[{root,"ComponentMixing.wl"}]];
saved=Import[FileNameJoin[{root,"symbolic.json"}],"RawJSON"];
plans=OtherTasksField`TowerPlan /@ {6,4,3};
mixing=OtherTasksMixing`MixingProofs[];
field=FiniteField[2,(#^12+#^6+#^4+#+1)&];
delta=2^-128;shots=102400;
upper=1-delta^(1/shots);
repetition=N[Log[2]/-Log[1/2+2 Sqrt[upper]],40];
baseline=Import[FileNameJoin[{root,"dqi-baseline.json"}],"RawJSON"];
tests={
 VerificationTest[plans,saved["tower_plans"],TestID->"reproduce-tower-maps"],
 VerificationTest[OtherTasksField`SixProductProof[],saved["six_product"],TestID->"reproduce-six-product-proof"],
 VerificationTest[mixing,saved["mixing"],TestID->"reproduce-mixing-proofs"],
 VerificationTest[Lookup[plans,"tensor_certificate"],{True,True,True},TestID->"all-bilinear-coefficients"],
 VerificationTest[Lookup[plans,"subfield_root_certificate"],{True,True,True},TestID->"matching-subfield-embeddings"],
 VerificationTest[Lookup[plans,"subfield_products"],{3,5,7},TestID->"interpolation-product-counts"],
 VerificationTest[Dimensions /@ Lookup[plans,"evaluation_maps"],{{3,6,12},{5,4,12},{7,3,12}},TestID->"linear-map-dimensions"],
 VerificationTest[saved["six_product"]["certificate"],True,TestID->"symbolic-gf8-six-products"],
 VerificationTest[And@@(#^4094 #==field[1]& /@ (field /@ Range[1,4095])),True,TestID->"all-nonzero-inverses"],
 VerificationTest[field[0]^4094,field[0],TestID->"zero-inverse-extension"],
 VerificationTest[field[2]^4095,field[1],TestID->"native-generator-order-divides-4095"],
 VerificationTest[And@@(field[2]^(4095/#)!=field[1]& /@ {3,5,7,13}),True,TestID->"native-generator-primitive"],
 VerificationTest[FailureQ[OtherTasksField`TowerPlan[5]],True,TestID->"unsupported-extension-degree"],
 VerificationTest[FailureQ[OtherTasksField`TowerPlan[3.]],True,TestID->"inexact-field-degree"],
 VerificationTest[mixing["two_body_bilinear"],True,TestID->"general-sos-tensor-invariant"],
 VerificationTest[mixing["offset_contraction"],True,TestID->"general-h1-offset-invariant"],
 VerificationTest[mixing["offset_squared_norm"],True,TestID->"general-energy-invariant"],
 VerificationTest[mixing["exact_nontrivial_gain"],2,TestID->"nontrivial-exact-normalization-gain"],
 VerificationTest[OtherTasksMixing`PairRotation[0],IdentityMatrix[2],TestID->"identity-mixing"],
 VerificationTest[mixing["large_parameter_limit"],True,TestID->"limiting-sign-flip"],
 VerificationTest[FailureQ[OtherTasksMixing`PairRotation[I]],True,TestID->"nonreal-rotation-rejected"],
 VerificationTest[FailureQ[OtherTasksMixing`PairRotation[{}]],True,TestID->"malformed-rotation-rejected"],
 VerificationTest[N[upper,30]>0,True,TestID->"zero-observed-errors-nonzero-confidence-bound"],
 VerificationTest[Abs[N[upper,30]-baseline["effective_error_rate"]]<10^-16,True,TestID->"official-confidence-bound"],
 VerificationTest[Abs[repetition-baseline["repetition_factor"]]<10^-14,True,TestID->"official-repetition-penalty"]
};
report=TestReport[tests];
summary=<|"kernel"->$Version,"tests"->Length[tests],"succeeded"->report["TestsSucceededCount"],
 "failed"->report["TestsFailedCount"],"symbolic_sha256"->FileHash[FileNameJoin[{root,"symbolic.json"}],"SHA256","HexString"],
 "verifier_sha256"->FileHash[$InputFileName,"SHA256","HexString"]|>;
Export[FileNameJoin[{root,"wl-checks.json"}],summary,"RawJSON"];
Echo[summary];
If[report["TestsFailedCount"]!=0,Exit[1]];
