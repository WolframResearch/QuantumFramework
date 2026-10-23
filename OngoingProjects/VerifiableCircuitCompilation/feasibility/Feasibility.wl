(* Exact certificates and resource inputs for invsqrt architecture research.
   Numeric minimax is discovery only. Resolve sees dyadic coefficients, never floats. *)
BeginPackage["InvSqrtFeasibility`"];
DiscoverDesign::usage = "DiscoverDesign[name,family,r,degree,w,d] discovers dyadic piecewise polynomials and proves their continuum bounds.";
NewtonRelativeError::usage = "NewtonRelativeError[e] is the exact relative error after one inverse-square-root Newton step.";
SquareDropBound::usage = "SquareDropBound[w,b] bounds all omitted bit products in a w-bit square with b fractional output bits.";
HornerDropBound::usage = "HornerDropBound[r,degree,w,d] bounds per-row truncation in local Horner arithmetic.";
CertifyPolynomial::usage = "CertifyPolynomial[p,t,a,h,eta] proves positivity, relative error and derivative bounds with exact arithmetic.";
PhaseVariants::usage = "PhaseVariants[] exports exact parity phases for the retained final partial products.";
Begin["`Private`"];

NewtonRelativeError[e_] := Expand[(1+e) (3-(1+e)^2)/2-1];
SquareDropBound[w_Integer?Positive,b_Integer?Positive] :=
  Total[2^(-2 Select[Range[w],2 #>b&])] +
  Total[2^(-Total[#]+1)& /@ Select[Subsets[Range[w],{2}],Total[#]-1>b&]];
HornerDropBound[r_Integer?NonNegative,degree_Integer?Positive,w_Integer?Positive,d_Integer?Positive] :=
  Max[0,Min[d,20]-r] 2^-w Sum[2^(-r j),{j,0,degree-1}];

CertifyPolynomial[p_,t_Symbol,a_,h_,eta_] /;
    PolynomialQ[p,t] && VectorQ[Join[CoefficientList[p,t],{a,h,eta}],MatchQ[#,_Integer|_Rational]&] &&
    a>=1 && h>0 && 0<eta<1 := <|
  "relative" -> TimeConstrained[Resolve[ForAll[t,0<=t<=h,
      p>0 && (1-eta)^2 <= (a+t) p^2 <= (1+eta)^2],Reals],15,"TIMEOUT"],
  "derivative" -> TimeConstrained[Resolve[ForAll[t,0<=t<=h,-1<=D[p,t]<=0],Reals],15,"TIMEOUT"]|>;
CertifyPolynomial[___] := Failure["CertificateDomain",<|"MessageTemplate"->
  "Expected a rational polynomial, rational a>=1, h>0 and 0<eta<1."|>];

DiscoverDesign[name_String,family:("A"|"B"),r_Integer?Positive,degree_Integer?Positive,
    w_Integer?Positive,d_Integer?Positive] := Module[
  {h=2^-r,t,intervals,seedBound,hornerBound,deltaBound,totalSeed,nrBounds,nrBits,bound},
  Needs["FunctionApproximations`"];
  intervals=Map[Function[i,Module[{a=1+i h,approx,coefs,p,eta,cert},
    approx=TimeConstrained[FunctionApproximations`MiniMaxApproximation[1/Sqrt[a+t],
      {t,{0,h},degree,0},WorkingPrecision->50],20,$Failed];
    If[!MatchQ[approx,{_List,{_,_?NumericQ}}],Return[Failure["Discovery",<|"interval"->i|>]]];
    coefs=Round[2^w CoefficientList[approx[[2,1]],t]];
    p=coefs . t^Range[0,degree]/2^w;
    eta=Ceiling[(Abs[approx[[2,2]]] + 2^-w Sum[h^j,{j,0,degree}]) 2^40]/2^40;
    cert=CertifyPolynomial[p,t,a,h,eta];
    With[{stages=FoldList[#2-t #1&,Last[Abs[coefs]]/2^w,Reverse[Most[Abs[coefs]]]/2^w],
        hb=HornerDropBound[r,degree,w,d]},
      cert["horner_registers"]=TimeConstrained[Resolve[ForAll[t,0<=t<=h,
        And@@Thread[hb<stages] && And@@Thread[stages<1-hb]],Reals],15,"TIMEOUT"]];
    <|"index"->i,"a"->ToString[a,InputForm],"coefs"->coefs,
      "eta"->eta,"certificates"->cert,
      "alternating"->And@@Thread[Sign[coefs]==(-1)^Range[0,degree]],
      "numeric_relative_error"->N[Max[Abs[Sqrt[a+#] (p/.t->#)-1]& /@ Subdivide[0,h,128]],16]|>
    ]],Range[0,3 2^r-1]];
  If[AnyTrue[intervals,FailureQ],Return[Failure["Intervals",<|"name"->name|>]]];
  seedBound=Max[Lookup[intervals,"eta"]];
  hornerBound=HornerDropBound[r,degree,w,d];
  deltaBound=If[d>=20,0,2^-d];
  totalSeed=seedBound+2 (hornerBound+deltaBound);
  nrBounds=Association@Table[b -> ((1+totalSeed) (4 SquareDropBound[w,b]+20 2^-b)+w 2^-b),{b,28,36}];
  nrBits=SelectFirst[Range[28,36],3 totalSeed^2/2+totalSeed^3/2+nrBounds[#] < 2^-25&,36];
  bound=If[family=="A",3 totalSeed^2/2+totalSeed^3/2+nrBounds[nrBits],totalSeed];
  <|"name"->name,"family"->family,"degree"->degree,"interval_bits"->r,
    "coef_bits"->w,"delta_bits"->d,"nr_bits"->nrBits,
    "table"->Lookup[intervals,"coefs"],"intervals"->Map[Append[#, "eta"->ToString[#["eta"],InputForm]]&,intervals],
    "eta_bound"->ToString[seedBound,InputForm],
    "horner_absolute_bound"->ToString[hornerBound,InputForm],
    "delta_absolute_bound"->ToString[deltaBound,InputForm],
    "seed_relative_bound"->ToString[totalSeed,InputForm],
    "certified_relative_bound"->ToString[bound,InputForm],
    "margin_fraction"->N[1-bound 2^24,16],
    "newton_arithmetic_bound"->ToString[nrBounds[nrBits],InputForm],
    "certified"->(And@@(TrueQ /@ Flatten[Values /@ Lookup[intervals,"certificates"]]) &&
      And@@Lookup[intervals,"alternating"] && bound<2^-24)|>
];

(* A product of n bits is a signed sum of nonempty subset parities / 2^(n-1).
   Coefficients here are turns over Q, never GF(2). Equal parities are collected. *)
parityTerms[bits_List,turns_] := Map[
  {Sort[#],(-1)^(Length[#]-1) turns/2^(Length[bits]-1)}&,Subsets[bits,{1,Length[bits]}]];
collectTerms[terms_List] := Map[
  <|"bits"->First[#][[1]],"turns"->ToString[Mod[Total[#[[All,2]]],1],InputForm]|>&,
  Select[GatherBy[terms,First],Mod[Total[#[[All,2]]],1]!=0&]];
phaseBlock[s_Integer,controlled_] := With[
  {pairs=Flatten[Table[Table[{17-k,17+j,2^(j-27-k+s-11)},
      {j,Max[0,27+k-28],28}],{k,1,17}],1]},
  <|"s"->s,"retained_pairs"->Length[pairs],"terms"->collectTerms[Catenate[
    parityTerms[If[controlled,{#[[1]],#[[2]],46},#[[1;;2]]],#[[3]]]& /@ pairs]]|>];
PhaseVariants[] := <|"variants"->{
  <|"name"->"fixed-s0-optimistic","exponent_blocks"->{phaseBlock[0,False]}|>,
  <|"name"->"all-exponents","exponent_blocks"->(phaseBlock[#,True]& /@ Range[0,10])|>}|>;

End[];
EndPackage[];
