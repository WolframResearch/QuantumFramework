repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
PacletDirectoryLoad[FileNameJoin[{repo,"QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting=False;
Echo[<|"Version"->$Version,"Loaded"->FindFile["Wolfram`QuantumFramework`"]|>];
rows = Internal`Bag[];
record[label_, value_] := With[{r=<|"Label"->label,"Pass"->TrueQ[value]|>},Echo[r];Internal`StuffBag[rows,r]];
same[a_,b_] := a["Signs"]===b["Signs"] && a["Tableau"]===b["Tableau"];
canonical[ps_] := PauliStabilizer[<|"Signs"->ps["Signs"],"Tableau"->ps["Tableau"]|>];
replaceV[(g:("V"|SuperDagger["V"])) -> q_Integer] := {"H"->q,If[g==="V","S",SuperDagger["S"]]->q,"H"->q};
replaceV[spec_] := {spec};
SeedRandom[417522];
sg = {"H","S",SuperDagger["S"],"X","Y","Z","V",SuperDagger["V"]};
makeSpecs[n_,m_] := Table[With[{g=RandomChoice[Join[sg,If[n>1,{"CNOT","CX","CZ","SWAP"},{}]]],w=RandomSample[Range[n],Min[n,2]]},g->If[MemberQ[{"CNOT","CX","CZ","SWAP"},g],w,First[w]]],{m}];
sm["H"]={{1,1},{1,-1}}/Sqrt[2]; sm["S"]=DiagonalMatrix[{1,I}]; sm[SuperDagger["S"]]=DiagonalMatrix[{1,-I}];
sm["X"]=PauliMatrix[1];sm["Y"]=PauliMatrix[2];sm["Z"]=PauliMatrix[3];
sm["V"]={{1+I,1-I},{1-I,1+I}}/2;sm[SuperDagger["V"]]=ConjugateTranspose[sm["V"]];
embed[u_,q_,n_] := KroneckerProduct@@ReplacePart[ConstantArray[IdentityMatrix[2],n],q->u];
dense[g_->q_Integer,n_] := embed[sm[g],q,n];
dense[g_->{a_Integer,b_Integer},n_] := Module[{bits=IntegerDigits[Range[0,2^n-1],2,n],out},
 If[g==="CZ",Return[DiagonalMatrix[1-2 bits[[All,a]] bits[[All,b]]]]];
 out=Map[Function[v,Switch[g,"CNOT"|"CX",ReplacePart[v,b->Mod[v[[a]]+v[[b]],2]],"SWAP",ReplacePart[v,{a->v[[b]],b->v[[a]]}]]],bits];
 Transpose[IdentityMatrix[2^n][[FromDigits[#,2]+1&/@out]]]
];
pm[x_,z_] := KroneckerProduct@@MapThread[Function[{a,b},PauliMatrix[Switch[{a,b},{0,0},0,{1,0},1,{1,1},2,{0,1},3]]],{x,z}];
Scan[Function[trial,Module[{n=3,specs,ps,u,t,expected,actual},
 specs=makeSpecs[n,60];ps=PauliStabilizer[n]["ApplyCircuit",specs];u=Fold[dense[#2,n].#1&,IdentityMatrix[2^n],specs];t=ps["Tableau"];
 expected=Join[Table[u.embed[PauliMatrix[1],q,n].ConjugateTranspose[u],{q,n}],Table[u.embed[PauliMatrix[3],q,n].ConjugateTranspose[u],{q,n}]];
 actual=Table[ps["Signs"][[k]] pm[t[[1,All,k]],t[[2,All,k]]],{k,2 n}];
 record["DenseExactConjugation/"<>ToString[trial],Simplify[expected-actual]===ConstantArray[0,{2 n,2^n,2^n}]];
 record["DenseUnitarity/"<>ToString[trial],Simplify[u.ConjugateTranspose[u]]===IdentityMatrix[2^n]];
 record["DenseNumericReference/"<>ToString[trial],Max[Abs[N[expected-actual]]] < 10^-12]
]],Range[5]];
record["BellExact",PauliStabilizer[2]["ApplyCircuit",{"H"->1,"CNOT"->{1,2}}]["Stabilizers"]==={"XX","ZZ"}];
Scan[Function[n,Module[{specs=makeSpecs[n,60],p=PauliStabilizer[n],bulk,scalar,mat,j},
 bulk=p["ApplyCircuit",specs]; scalar=Fold[#1[#2]&,canonical[p],Catenate[replaceV /@ specs]];
 record["BoundaryExactScalar/n="<>ToString[n],same[bulk,scalar]];
 record["PackedResident/n="<>ToString[n],KeyExistsQ[First[bulk],"GeneratorX"] && Developer`PackedArrayQ[First[bulk]["GeneratorX"]] && Developer`PackedArrayQ[First[bulk]["GeneratorZ"]]];
 If[n<=125,mat=bulk["Matrix"];j=ArrayFlatten[{{ConstantArray[0,{n,n}],IdentityMatrix[n]},{IdentityMatrix[n],ConstantArray[0,{n,n}]}}];record["Symplectic/n="<>ToString[n],Mod[mat.j.Transpose[mat],2]===j]];
]],{1,2,31,32,61,62,63,124,125,1024}];
(* Exact wire arithmetic must remain exact even when packed machine subtraction overflows. *)
Scan[Function[w,record["ExactEncoding/"<>ToString[w],encodeStabilizerGates[{"H"->w,"V"->{w},"CNOT"->{w,-w}}]==={{0,w-1,0},{0,w-1,0},{1,w-1,0},{0,w-1,0},{3,w-1,-w-1}}]],{-2^100,-2^63-1,-2^63,-2^63+1,-1,0,1,2^63-1,2^63,2^63+1,2^100}];
record["LegacyControlled",encodeStabilizerGates[{{"H",2},{"V",3},{SuperDagger["V"],1},"C"["X"->{3},{2},{}],"C"["Z"->{1},{3},{}]}]===encodeStabilizerGates[{"H"->2,"V"->3,SuperDagger["V"]->1,"CNOT"->{2,3},"CZ"->{3,1}}]];
record["EmptyIdentity",With[{p=PauliStabilizer[4]},p["ApplyCircuit",{}]===p && encodeStabilizerGates[{}]==={}]];
record["EncoderFallback",And@@(encodeStabilizerGates[#]===$Failed&/@{{"T"->1},{"Unknown"->1},{"H"->1.},{"H"->s},{"H"->{1,2}},{"CNOT"->{1}},{"H"->1,"T"->2}})];
record["NonCliffordFrame",Head[PauliStabilizer[2]["ApplyCircuit",{"H"->1,"T"->1,"CNOT"->{1,2}}]]===StabilizerFrame];
Module[{p=PauliStabilizer[2]["ApplyCircuit",{"H"->1,"CNOT"->{1,2}}],sym,specs,evolved},
 sym=p["SymbolicMeasure",1];specs={"H"->2,"S"->1,"CNOT"->{1,2},"X"->1,"Y"->2};evolved=sym["ApplyCircuit",specs];
 record["SymbolicFallbackExact",!psConcreteFastQ[sym] && same[evolved,Fold[#1[#2]&,sym,specs]]];
 record["SymbolicSamplesValid",AllTrue[evolved["SampleOutcomes",20],PauliStabilizerQ]]
];
(* Enumerate the old and new gate shape guards over list, array, malformed, and exact bignum data. *)
shapes=Join[{{},{{}},{{0,0,0}},{{0,0,0},{1,1,0}},{{0,0}},{{0,0,0,0}},{{0,0,0},{1,0}},{{0.,0,0}},{{0,s,0}},{0,0,0},{{0,0,0},{}},SparseArray[{{0,0,0}}],{{2^100,0,0}},{{0,-2^100,0}}},Tuples[{{0,0,0},{1,1,0}},2]];
record["ShapeGuardEquivalence",And@@Map[Function[g,MatchQ[g,{{_Integer,_Integer,_Integer}..}]===TrueQ[ListQ[g]&&MatrixQ[g,IntegerQ]&&MatchQ[Dimensions[g],{_Integer,3}]]],shapes]];
(* Observe native entry by replacing only the compiled function in this temporary kernel. *)
Module[{entry=0,p=PauliStabilizer[3],invalids={{0,-1,0},{0,3,0},{3,1,1},{4,2,2},{0,-2^63-1,0},{0,2^63-1,0}}},
 Internal`InheritedBlock[{compiledCliffordFold},Clear[compiledCliffordFold];compiledCliffordFold[x_,z_,s_,valid_,g_]:=(entry++;Join[x,z,{s}]);
 Scan[Function[g,record["InvalidWireResult/"<>ToString[g,InputForm],applyCompiledFold[p,{g}]===$Failed]],invalids];
 record["InvalidWiresNeverNative",entry===0];
 Scan[Function[g,If[g=!={},applyCompiledFold[p,g]]],Select[shapes,Not[MatchQ[#,{{_Integer,_Integer,_Integer}..}]]&]];
 record["BadShapeNeverNative",entry===0];
 applyCompiledFold[p,{{0,0,0}}];record["ValidIntegerListEntersNative",entry===1]
 ]
];
report=TestReport[FileNameJoin[{repo,"Tests","Stabilizer","CompiledEvolution.wlt"}]];
record["CompiledEvolution65",Length[report["Results"]]===65 && AllTrue[report["Results"],#["Outcome"]==="Success"&]];
Echo[Counts[#["Outcome"]&/@report["Results"]]];
Export["/private/tmp/coldfix-independent-correctness.json",<|"Version"->$Version,"Rows"->Internal`BagPart[rows,All],"Failures"->Select[Internal`BagPart[rows,All],!#Pass&]|>,"RawJSON"];
Exit[If[AllTrue[Internal`BagPart[rows,All],#Pass&],0,1]];
