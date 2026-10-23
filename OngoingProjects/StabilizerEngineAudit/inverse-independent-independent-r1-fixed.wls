repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
PacletDirectoryLoad[repo <> "/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
Echo[<|"Version" -> $Version, "Loaded" -> FindFile["Wolfram`QuantumFramework`"], "Process" -> $ProcessID|>];
$ProgressReporting = False;
check = <||>;
record[name_, value_] := (check[name] = TrueQ[value]; Echo[{name, check[name]}]);
sameAction[a_, b_] := a["Signs"] === b["Signs"] && Normal[a["Matrix"]] === Normal[b["Matrix"]];
canonical[p_] := PauliStabilizer[<|"Signs" -> p["Signs"], "Matrix" -> Normal[p["Matrix"]]|>];
layouts[p_] := {canonical[p], psFromPacked[psGetPacked[p]], psToGeneratorPacked[p]};
shape[p_] := Which[psGeneratorPackedQ[p], "Generator", psPackedQ[p], "Chunk", True, "Canonical"];
tensor[ms_List] := Fold[KroneckerProduct, {{1}}, ms];
embed[g_, wire_, n_] := tensor[ReplacePart[ConstantArray[IdentityMatrix[2], n], wire -> g]];
rowMatrix[p_, j_] := With[{xz = Normal[p["Matrix"]][[j]], n = p["Qubits"]},
 p["Signs"][[j]] tensor[MapThread[Switch[{#1, #2}, {0,0}, PauliMatrix[0], {1,0}, PauliMatrix[1], {1,1}, PauliMatrix[2], {0,1}, PauliMatrix[3]] &, {Take[xz,n],Drop[xz,n]}]]];
generators[n_] := Join[Table[embed[PauliMatrix[1], q, n], {q,n}], Table[embed[PauliMatrix[3], q, n], {q,n}]];
(* Construct U directly from its signed X/Z conjugations; no QF circuit synthesis or composition. *)
reconstruct[p_] := Module[{n=p["Qubits"], image, projector, column, vacuum},
 image = rowMatrix[p,#]& /@ Range[2 n];
 projector = Fold[(IdentityMatrix[2^n] + #2) . #1 /2 &, IdentityMatrix[2^n], Drop[image,n]];
 column = First @ Select[Transpose[projector], !AllTrue[#,PossibleZeroQ]&];
 vacuum = Simplify[column/Sqrt[Conjugate[column].column]];
 Simplify @ Transpose[Map[Function[bits, Fold[#2 . #1 &, vacuum, Pick[Take[image,n],bits,1]]], Tuples[{0,1},n]]]
];
denseCheck[p_] := With[{u = reconstruct[p], gs = generators[p["Qubits"]], inv = p["Dagger"]},
 Simplify[ConjugateTranspose[u].u] === IdentityMatrix[Length[u]] &&
 And @@ MapThread[(Simplify[u.#1.ConjugateTranspose[u] - #2] === ConstantArray[0,Dimensions[u]])&, {gs,rowMatrix[p,#]& /@ Range[Length[gs]]}] &&
 And @@ MapThread[(Simplify[ConjugateTranspose[u].#1.u - #2] === ConstantArray[0,Dimensions[u]])&, {gs,rowMatrix[inv,#]& /@ Range[Length[gs]]}]
];
fixture = PauliStabilizer[<|"Signs"->{1,1,1,-1,1,1},"Tableau"->{{{1,1,0,1,1,1},{1,1,0,1,1,0},{1,0,0,1,0,1}},{{1,1,1,0,0,0},{0,0,1,1,0,1},{0,0,0,1,1,0}}}|>];
record["FixtureDenseAll64SignedActions", And @@ Map[Function[bits,denseCheck[PauliStabilizer[<|"Matrix"->fixture["Matrix"],"Signs"->1-2 bits|>]]],Tuples[{0,1},6]]];
one = FixedPoint[DeleteDuplicatesBy[Join[#, Catenate[Map[Function[p,{p["H",1],p["S",1]}],#]]],Function[p,{p["Matrix"],p["Signs"]}]] &, {PauliStabilizer[1]}];
record["All24DensePhysicalAdjoints", Length[one]===24 && And@@(denseCheck/@one)];
record["LayoutsFullInvolution", And @@ Map[Function[p, With[{q=p["Dagger"]},q["Dagger"]===p && p["Inverse"]===q && shape[p]===shape[q] && sameAction[p[q],PauliStabilizer[p["Qubits"]]]]],Catenate[layouts /@ Join[one,{fixture}]]]];
adjointSpec[(g_ -> wires_)] := Replace[g,{"S"->SuperDagger["S"],SuperDagger["S"]->"S","V"->SuperDagger["V"],SuperDagger["V"]->"V"}] -> wires;
SeedRandom[91732026];
specs[n_, len_] := Table[With[{g=RandomChoice[{"H","S",SuperDagger["S"],"X","Y","Z","V",SuperDagger["V"],"CNOT","CZ","SWAP"}],w=RandomSample[Range[n],2]},g->If[MemberQ[{"CNOT","CZ","SWAP"},g],w,First[w]]],{len}];
record["BoundariesExactAllSignedGenerators", And@@Map[Function[n,With[{s=specs[n,130]},With[{p=PauliStabilizer[n]["ApplyCircuit",s],q=PauliStabilizer[n]["ApplyCircuit",adjointSpec/@Reverse[s]]},And@@Map[Function[r,sameAction[r["Dagger"],q] && r["Dagger"]["Dagger"]===r && sameAction[r[r["Dagger"]],PauliStabilizer[n]] && sameAction[r["Dagger"][r],PauliStabilizer[n]]],layouts[p]]]]],{31,32,61,62,63,124,125}]];
record["RandomDenseAction2Through4",And@@Map[Function[n,And@@Table[denseCheck[PauliStabilizer[n]["ApplyCircuit",specs[n,35]]],{4}]],{2,3,4}]];
variables=Array[bit,6];
symbolic=PauliStabilizer[<|"Matrix"->fixture["Matrix"],"Signs"->1-2 variables|>];
symbolicInverse=symbolic["Dagger"];
record["BinarySymbolicAll64Assignments",And@@Map[Function[bits,With[{rules=Thread[variables->bits]},sameAction[symbolicInverse/.rules,canonical[symbolic/.rules]["Dagger"]]]],Tuples[{0,1},6]]];
record["EmptyCircuitIdentities",And@@Map[Function[n,With[{p=PauliStabilizer[n]},sameAction[p["Dagger"],p] && p["ApplyCircuit",{}]===p && p["Dagger"]["Dagger"]===p]],{1,2,31,62,63,125}]];
report=TestReport[{VerificationTest[PauliStabilizer[{"XX","XX"}]["Dagger"],$Failed,{PauliStabilizer::singular}],VerificationTest[PauliStabilizer[{"XX","XX"}]["Inverse"],$Failed,{PauliStabilizer::singular}]}];
record["SingularBothAliasesMessageContract",And@@(#["Outcome"]==="Success"&/@report["Results"])];
Echo[<|"Checks"->check,"Messages"->ToString[$MessageList,InputForm]|>];
Export["/private/tmp/qf-inverse-independent-r1.json",<|"Version"->$Version,"Checks"->check,"Messages"->ToString[$MessageList,InputForm]|>,"RawJSON"];
Exit[If[And@@Values[check],0,1]];
