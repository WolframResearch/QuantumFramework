BeginPackage["CompilationRepresentations`"];
CachedQuotient::usage = "CachedQuotient[sigma,omega,z,h] caches omega/(z^h sigma') modulo sigma over a finite field. Requires positive degree, h>=0 and an invertible denominator; otherwise returns Failure.";
RepresentationChecks::usage = "RepresentationChecks[] returns exact symbolic and finite-field checks for three compilation representations.";
Begin["`Private`"];

CachedQuotient[sigma_, omega_, z_Symbol, h_Integer?NonNegative] := Module[
  {denominator, bezout, gcd, inverse},
  If[!PolynomialQ[sigma,z] || !PolynomialQ[omega,z] || Exponent[sigma,z]<1,
    Return[Failure["PolynomialDomain",<||>]]];
  denominator = PolynomialRemainder[z^h D[sigma,z],sigma,z];
  bezout = PolynomialExtendedGCD[sigma,denominator,z];
  gcd = First[bezout];
  If[!FreeQ[gcd,z] || TrueQ[gcd==0],
    Return[Failure["NoninvertibleDenominator",<||>]]];
  inverse = Last[Last[bezout]]/gcd;
  <|"denominator"->denominator,
    "inverse"->PolynomialRemainder[inverse,sigma,z],
    "quotient"->PolynomialRemainder[omega inverse,sigma,z]|>
];
CachedQuotient[___] := Failure["Arguments",<||>];

(* A polynomial certificate is stronger than testing roots alone. *)
quotientCase[field_, roots_List, h_Integer] := Module[
  {z, sigma, omega, cache, zero},
  zero=field[0];
  sigma=Expand[Times@@(z-#& /@ roots)];
  omega=Total[MapIndexed[field[17 First[#2]+3] z^(First[#2]-1)&,roots]];
  cache=CachedQuotient[sigma,omega,z,h];
  If[FailureQ[cache],Return[False]];
  TrueQ[PolynomialRemainder[cache["quotient"] cache["denominator"]-omega,sigma,z]==zero] &&
    AllTrue[roots,Function[root,TrueQ[
      (cache["quotient"]/.z->root)==(omega/.z->root)/(root^h (D[sigma,z]/.z->root))]]]
];

mobiusMatrix[n_Integer] := Nest[KroneckerProduct[#,{{1,0},{1,1}}]&,{{1}},n];

(* Match the released loop, not a differently indexed Forney convention. *)
nativeScanChecks[field_] := Module[
  {z,roots,sigma,evaluator,omega,cache,wrong,slotValues,trace,phaseCase},
  roots=field /@ Range[2,36];
  sigma=Expand[Times@@(z-#& /@ roots)];
  evaluator=field /@ Range[1,35];
  omega=Total[evaluator z^Range[34,0,-1]];
  cache=CachedQuotient[sigma,omega,z,1];
  slotValues[root_]:=Total[evaluator root^Range[34,0,-1]]/
    Total[(Coefficient[sigma,z,#]& /@ Range[1,35,2]) root^Range[1,35,2]];
  trace[element_]:=Total[NestList[#^2&,element,11]];
  phaseCase[degree_]:=With[{localRoots=field /@ Range[2,degree+1]},
    With[{locator=Expand[Times@@(z-#& /@ localRoots)]},
      With[{quotient=CachedQuotient[locator,omega,z,1]["quotient"]},
        AllTrue[localRoots,Function[root,
          With[{native=Total[evaluator root^Range[34,0,-1]]/
            Total[(Coefficient[locator,z,#]& /@ Range[1,35,2]) root^Range[1,35,2]],
            cached=quotient/.z->root},
            AllTrue[field /@ (2^Range[0,11]),Function[mask,
              (-1)^Boole[trace[mask native]==field[1]]==
              (-1)^Boole[trace[mask cached]==field[1]]]]]]]]]];
  wrong=CachedQuotient[sigma,field[1],z,34];
  <|"native_scan_reversed_evaluator_h1"->AllTrue[roots,
      TrueQ[(cache["quotient"]/.z->#)==slotValues[#]]&],
    "native_scan_masked_phases_degrees_1_2_17_35"->AllTrue[{1,2,17,35},phaseCase],
    "native_derivative_accumulator_is_z_sigma_prime"->TrueQ[
      Expand[Total[(Coefficient[sigma,z,#] z^#)& /@ Range[1,35,2]]-z D[sigma,z]]==field[0]],
    "incorrect_h34_adapter_counterexample_detected"->TrueQ[
      (wrong["quotient"]/.z->First[roots])!=First[roots]^33/(D[sigma,z]/.z->First[roots])]|>
];

RepresentationChecks[] := Module[
  {a,b,c,d,h,t,u,v,z,field,small,polynomialError,compensation,
   lowerBound,matrices,prefix,trace,gram,phaseSigns,linearMasks,
   om,ps,kk,ss,mm,qq},
  polynomialError=Expand[(a+d h^2/8)-(b+d h)t+(c+d)t^2-(a-b t+c t^2)];
  compensation=TrueQ[Expand[polynomialError-d(t^2-h t+h^2/8)]==0];
  lowerBound=Expand[((t^2+u t+v)/.t->0)-2 ((t^2+u t+v)/.t->1/2)+((t^2+u t+v)/.t->1)];
  field=FiniteField[2,(#^12+#^6+#^4+#+1)&];
  small=FiniteField[2,(#^4+#+1)&];
  matrices=mobiusMatrix /@ Range[1,6];
  prefix=AllTrue[matrices,Function[matrix,
    AllTrue[Range[Length[matrix]-1],Function[cut,
      Total[Flatten[matrix[[1;;cut,cut+1;;]]]]==0]]]];
  trace[element_]:=Total[NestList[#^2&,element,11]];
  gram=Outer[trace[#1 #2]&,field /@ (2^Range[0,11]),field /@ (2^Range[0,11])];
  linearMasks=Det[gram]!=field[0];
  (* Root flag times a linear output bit is quadratic Boolean phase: CZ. *)
  phaseSigns=Flatten[Outer[(-1)^(#1 #2)&,{0,1},{0,1}]];
  Join[<|
    "curvature_identity"->compensation,
    "curvature_uniform_bound"->TrueQ[Resolve[ForAll[t,h>0 && 0<=t<=h,
      -h^2/8<=t^2-h t+h^2/8<=h^2/8],Reals]],
    "curvature_equioscillation"->TrueQ[((t^2-h t+h^2/8)/.{t->#}& /@ {0,h/2,h})=={h^2/8,-h^2/8,h^2/8}],
    "affine_minimax_lower_certificate"->TrueQ[lowerBound==1/2],
    "integer_intercept_interval"->TrueQ[Resolve[ForAll[{a,b,c,u},
      Element[{a,b,c,u},Reals] && u>0,
      Equivalent[b<=(a-c) u<=h,(b/u+c)<=a<=(h/u+c)]] ,Reals]],
    "curvature_zero_change_limit"->TrueQ[Limit[polynomialError,d->0]==0],
    "GF12_degree35_general_quotient_h34"->quotientCase[field,field /@ Range[2,36],34],
    "GF16_degree1_prefactor0"->quotientCase[small,{small[2]},0],
    "GF16_degree5_prefactor71"->quotientCase[small,small /@ Range[2,6],71],
    "general_quotient_ideal_certificate"->TrueQ[Last[PolynomialReduce[
      (om ps-kk ss) mm-om-ss(om qq-kk mm),
      {ps mm-1-qq ss},{ps,mm,om,kk,ss,qq}]]==0],
    "repeated_root_rejected"->FailureQ[CachedQuotient[Expand[(z-small[2])^2],small[1],z,0]],
    "zero_root_with_prefactor_rejected"->FailureQ[CachedQuotient[small[1] z,small[1],z,1]],
    "empty_locator_rejected"->FailureQ[CachedQuotient[small[1],small[1],z,0]],
    "malformed_prefactor_rejected"->FailureQ[CachedQuotient[z+small[1],small[1],z,-1]],
    "trace_dual_basis_all_4096_masks"->TrueQ[linearMasks],
    "cached_root_phase_unitary"->TrueQ[phaseSigns=={1,1,1,-1} && phaseSigns^2=={1,1,1,1}],
    "prefix_independent_of_invalid_suffix"->TrueQ[prefix],
    "mobius_involution_sizes_2_to_64"->AllTrue[matrices,Mod[#.#,2]==IdentityMatrix[Length[#]]&],
    "invertible_output_basis_preserves_support"->TrueQ[Resolve[ForAll[{u,v},
      Element[{u,v},Integers] && 0<=u<=1 && 0<=v<=1,
      Equivalent[u==0 && v==0,Mod[u+v,2]==0 && v==0]]]],
    "unitarity_scope"->"Only phase identity checked; no compiled workspace or cleanup claim."
  |>,nativeScanChecks[field]]
];
End[];
EndPackage[];
