repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
Echo[<|"Version" -> $Version, "LoadedFile" -> FindFile["Wolfram`QuantumFramework`"]|>];
results = <||>;
record[name_, value_] := (AssociateTo[results, name -> value]; Echo[name -> value]);
same[a_, b_] := a["Signs"] === b["Signs"] && a["Tableau"] === b["Tableau"];
canon[ps_] := PauliStabilizer[<|"Signs" -> ps["Signs"], "Tableau" -> ps["Tableau"]|>];
specs[n_, count_] := BlockRandom[SeedRandom[1837+n+count]; Table[With[{g = RandomChoice[{"H", "S", SuperDagger["S"], "V", SuperDagger["V"], "X", "Y", "Z", "CNOT", "CX", "SWAP", "CZ"}], w = RandomSample[Range[n], Min[n, 2]]}, If[n == 1, RandomChoice[{"H", "S", SuperDagger["S"], "V", SuperDagger["V"], "X", "Y", "Z"}] -> 1, g -> If[MemberQ[{"CNOT", "CX", "SWAP", "CZ"}, g], w, First[w]]]], {count}]];
(* Independent dense reference, entirely in System matrix/bit primitives. *)
unitary[g_ -> w_, n_] := Module[{wire = Flatten[{w}], bitstrings = Tuples[{0,1}, n], images, one},
 If[MemberQ[{"CNOT", "CX", "SWAP", "CZ"}, g],
  If[g === "CZ", Return[DiagonalMatrix[(-1)^(Times@@#[[wire]]) & /@ bitstrings]]];
  images = If[g === "SWAP", ReplacePart[#, Thread[wire -> Reverse[#[[wire]]]]] &, ReplacePart[#, wire[[2]] -> BitXor[#[[wire[[1]]]], #[[wire[[2]]]]]] &] /@ bitstrings;
  Return[Transpose[UnitVector[2^n, FromDigits[#, 2]+1] & /@ images]]];
 one = Switch[g, "H", {{1,1},{1,-1}}/Sqrt[2], "S", DiagonalMatrix[{1,I}], SuperDagger["S"], DiagonalMatrix[{1,-I}], "X", PauliMatrix[1], "Y", PauliMatrix[2], "Z", PauliMatrix[3], "V", ((1+I)IdentityMatrix[2]+(1-I)PauliMatrix[1])/2, SuperDagger["V"], ((1-I)IdentityMatrix[2]+(1+I)PauliMatrix[1])/2];
 KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], n], First[wire] -> one]
];
pauliMatrices[ps_] := MapThread[Times, {ps["Signs"], (KroneckerProduct @@ (PauliMatrix[Switch[#, {0,0},0, {1,0},1, {1,1},2, {0,1},3]] & /@ #) & /@ Transpose[ps["Tableau"], {3,2,1}])}];
zeroMatrixQ[m_] := And@@Flatten[Map[TrueQ[FullSimplify[#] == 0]&, m, {2}]];
record["DenseAllGeneratorsAndState", Table[Module[{s = specs[n, 50], initial = PauliStabilizer[n], evolved, u, actual, target, v},
 evolved = initial["ApplyCircuit",s]; u = Fold[unitary[#2,n].#1 &, IdentityMatrix[2^n],s];
 actual = pauliMatrices[evolved]; target = (u.#.ConjugateTranspose[u] &) /@ pauliMatrices[initial];
 v = Normal[evolved["State"]["StateVector"]];
 <|"n" -> n,"GeneratorsExact" -> And@@MapThread[zeroMatrixQ[#1-#2]&, {actual,target}],"Unitary" -> zeroMatrixQ[u.ConjugateTranspose[u]-IdentityMatrix[2^n]],"StateNorm" -> Simplify[Conjugate[v].v],"SignsConsistent" -> And@@(zeroMatrixQ[{#.v-v}]& /@ Take[actual,-n])|>],{n,1,4}]];
record["ScalarSplitBoundaries", Table[Module[{s=(specs[n,70]/.{SuperDagger["V"]->SuperDagger["S"],"V"->"H"}), p=PauliStabilizer[n], a,b,c},a=p["ApplyCircuit",s];b=Fold[#1[#2]&,canon[p],s];c=Fold[#1["ApplyCircuit",#2]&,p,Partition[s,10]];{n,same[a,b],same[a,c],PauliStabilizerQ[a]}],{n,{1,2,31,32,61,62,63,124,125,1024}}]];
record["Bell", PauliStabilizer[2]["ApplyCircuit",{"H"->1,"CNOT"->{1,2}}]["Stabilizers"]];
record["VScalarDense", Table[With[{p=PauliStabilizer[3]["ApplyCircuit",specs[3,20]], s=g->q}, same[p["ApplyCircuit",{s}],canon[p][s]] && And@@MapThread[zeroMatrixQ[#1-#2]&,{pauliMatrices[p["ApplyCircuit",{s}]],(unitary[s,3].#.ConjugateTranspose[unitary[s,3]]&) /@ pauliMatrices[p]}]],{g,{"V",SuperDagger["V"]}},{q,3}]];
record["SymbolicSigns", Module[{base = PauliStabilizer[3]["ApplyCircuit",specs[3,20]],sym,s, out,a,b},sym=PauliStabilizer[<|"Tableau"->base["Tableau"],"Signs"->(1-2 Array[Symbol["r"<>ToString[#]]&,6])|>];s=specs[3,20];out=sym["ApplyCircuit",s];a=out/.Thread[Array[Symbol["r"<>ToString[#]]&,6]->{0,1,1,0,1,0}];b=(sym/.Thread[Array[Symbol["r"<>ToString[#]]&,6]->{0,1,1,0,1,0}])["ApplyCircuit",s];same[a,b] && !psConcreteFastQ[sym] && PauliStabilizerQ[out]]];
record["ValidityMutation", Module[{p=PauliStabilizer[32],bad},Table[bad=Switch[k,1,ReplacePart[p,{1,"Signs"}->{}],2,ReplacePart[p,{1,"Qubits"}->33],3,ReplacePart[p,{1,"GeneratorX",1,1}->-1],4,PauliStabilizer[KeyDrop[First[p],"GeneratorZ"]],5,PauliStabilizer[Append[First[p],"GeneratorX"->{{1}}]],6,MapAt[Rest,p,{1,"Signs"}],7,ReplacePart[p,{1,"GeneratorZ",32,2}->4]];{k,PauliStabilizerQ[bad],System`Private`HoldValidQ[bad]}, {k,7}]]];
record["Compatibility", Module[{p=PauliStabilizer[5]["ApplyCircuit",specs[5,20]],c,old,partial},c=canon[p];old=psFromPacked[psGetPacked[c]];partial=PauliStabilizer[{"XXI","IZZ"}];{same[p["ApplyCircuit",specs[5,12]],old["ApplyCircuit",specs[5,12]]],same[p["Dagger"],c["Dagger"]],same[p["Dagger"]["Dagger"],c["Dagger"]["Dagger"]],PauliStabilizerQ[partial],psConcreteFastQ[partial],Head[partial["ApplyCircuit",{"H"->1}]]}]];
record["MeasurementBranches", Module[{p=PauliStabilizer[4]["ApplyCircuit",specs[4,30]],c,resultsP,resultsC},c=canon[p];resultsP=p["M",{1,3,2}];resultsC=c["M",{1,3,2}];Keys[resultsP]===Keys[resultsC]&&And@@MapThread[same,{Values[resultsP],Values[resultsC]}]]];
record["SymbolicBellCorrelations",Module[{p=PauliStabilizer[2]["ApplyCircuit",{"H"->1,"CNOT"->{1,2}}],s,a,b},s=p["SymbolicMeasure",{1,2}];a=s["SubstituteOutcomes",\[FormalS][1]->0];b=s["SubstituteOutcomes",\[FormalS][1]->1];{a["Stabilizers"],b["Stabilizers"],AllTrue[s["SampleOutcomes",20],PauliStabilizerQ]}]];
record["PackTraceSuccessiveBulk",Module[{p=PauliStabilizer[63]["ApplyCircuit",specs[63,20]],trace},trace=Trace[p["ApplyCircuit",specs[63,20]],_packGenRows|_unpackGenRows,TraceInternal->True];Cases[trace,(h:packGenRows|unpackGenRows)[args___]:>{ToString[Unevaluated[h],InputForm],Dimensions[{args}[[1]]]},Infinity]]];
record["StorageBytes",Table[Module[{p=PauliStabilizer[n]}, {n,ByteCount[p],ByteCount[canon[p]],Developer`PackedArrayQ[First[p]["GeneratorX"]],Dimensions[First[p]["GeneratorX"]]}],{n,{31,62,128,1024}}]];
record["PlanSemantics",Module[{p=PauliStabilizer[3], q1,q3,qh,qx,qparam,noncomp,a,b},q1=QuantumCircuitOperator[{"H"->1}];q3=QuantumCircuitOperator[{"H"->3}];qh=QuantumCircuitOperator[{QuantumOperator[(PauliMatrix[1]+PauliMatrix[3])/Sqrt[2],"Label"->"X"]}];qx=QuantumCircuitOperator[{QuantumOperator[PauliMatrix[1],"Label"->"X"]}];qparam=QuantumCircuitOperator[{"P"[theta]->1},"Parameters"->{theta}];stabilizerClearCaches[];a=q1[p];b=q1[p];{same[a,b],q3[p]["Stabilizers"],qh[PauliStabilizer[1]]["Stabilizers"],qx[PauliStabilizer[1]]["Stabilizers"],Head[qparam[<|theta->Pi/2|>][PauliStabilizer[1]]],Head[qparam[<|theta->Pi/4|>][PauliStabilizer[1]]]}]];
record["Empty",With[{p=PauliStabilizer[3]},p["ApplyCircuit",{}]===p]];
Export["/private/tmp/stabilizer-independent-r1.json",ToString[results,InputForm],"Text"];
Echo[results];
Exit[];
