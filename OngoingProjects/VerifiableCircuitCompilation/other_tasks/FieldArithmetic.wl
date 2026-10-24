(* Exact arithmetic certificates. Native kickmix basis: x^12+x^6+x^4+x+1. *)
BeginPackage["OtherTasksField`"];
TowerPlan::usage="TowerPlan[d] constructs exact degree-12 multiplication via interpolation over GF(2^d), d=3,4,6.";
SixProductProof::usage="SixProductProof[] proves the six-product GF(8) multiplication identity.";
Begin["`Private`"];
bits[e_,n_] := PadRight[e["Coefficients"],n];
word[v_] := FromDigits[Reverse[v],2];
moduli=<|3->11,4->19,6->67|>;
poly[m_,x_] := Reverse[IntegerDigits[m,2]].x^Range[0,IntegerLength[m,2]-1];

SixProductProof[] := With[
 {products={a0 b0,a1 b1,a2 b2,(a0+a1)(b0+b1),(a0+a2)(b0+b2),(a1+a2)(b1+b2)},
  outputs={{1,1,1,0,0,1},{1,0,0,1,0,1},{1,1,0,0,1,0}}},
 <|"input_forms"->{{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1}},
   "output_forms"->outputs,"products"->6,
   "certificate"->TrueQ[PolynomialMod[Expand[
      (outputs.products).{1,z,z^2}-PolynomialRemainder[
       (a0+a1 z+a2 z^2)(b0+b1 z+b2 z^2),z^3+z+1,z,Modulus->2]],2]===0]|>];

TowerPlan[d:(3|4|6)] := Module[
 {field=FiniteField[2,(#^12+#^6+#^4+#+1)&],subfield,degree=12/d,
  beta,alpha,basis,basisInverse,coordinates,embed,points,vandermonde,interpolation,
  evaluations,outputs,expected,reconstructed},
 subfield=FiniteField[2,With[{m=moduli[d]},poly[m,#]&]];
 alpha=field[2];
 beta=SelectFirst[alpha^((4095/(2^d-1)) Range[1,2^d-1]),
   poly[moduli[d],#]==field[0]&];
 basis=Transpose[Catenate[Table[bits[beta^i alpha^j,12],{j,0,degree-1},{i,0,d-1}]]];
 basisInverse=Inverse[basis,Modulus->2];
 coordinates[e_] := subfield /@ (word /@ Partition[Mod[basisInverse.bits[e,12],2],d]);
 embed[e_] := Total[MapThread[Times,{field /@ bits[e,d],NestList[beta #&,field[1],d-1]}]];
 points=Append[subfield /@ Range[0,2 degree-3],"infinity"];
 vandermonde=Map[Function[p,If[p==="infinity",
   Join[ConstantArray[subfield[0],2 degree-2],{subfield[1]}],
   NestList[p #&,subfield[1],2 degree-2]]],points];
 interpolation=Inverse[vandermonde];
 evaluations=Map[Function[p,Transpose[Map[Function[e,bits[
   If[p==="infinity",Last[coordinates[e]],
     coordinates[e].NestList[p #&,subfield[1],degree-1]],d]],field /@ (2^Range[0,11])]]],points];
 outputs=Table[Transpose[Table[bits[Total[Table[
   embed[interpolation[[j,t]] subfield[2^b]] alpha^(j-1),{j,1,2 degree-1}]],12],
   {b,0,d-1}]],{t,1,Length[points]}];
 expected=Table[bits[field[2^i] field[2^j],12],{i,0,11},{j,0,11}];
 reconstructed=Table[Mod[Total[MapThread[Function[{out,ev},
   out.bits[subfield[word[ev[[All,i+1]]]] subfield[word[ev[[All,j+1]]]],d]],
   {outputs,evaluations}]],2],{i,0,11},{j,0,11}];
 <|"subfield_degree"->d,"extension_degree"->degree,"native_modulus"->4179,
   "subfield_modulus"->moduli[d],"subfield_generator_native"->beta["Index"],
   "evaluation_maps"->evaluations,"output_maps"->outputs,
   "subfield_products"->Length[points],
   "subfield_root_certificate"->TrueQ[poly[moduli[d],beta]==field[0]],
   "basis_certificate"->TrueQ[Mod[basisInverse.basis,2]===IdentityMatrix[12]],
   "interpolation_certificate"->TrueQ[
     Map[#["Index"]&,interpolation.vandermonde,{2}]===IdentityMatrix[2 degree-1]],
   "tensor_certificate"->TrueQ[reconstructed===expected],
   "coefficient_tensor_entries_checked"->12^3|>
];
TowerPlan[___] := Failure["SubfieldDegree",<|"MessageTemplate"->"Expected subfield degree 3, 4 or 6."|>];
End[];
EndPackage[];
