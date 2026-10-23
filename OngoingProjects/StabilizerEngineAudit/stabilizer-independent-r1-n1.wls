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
 If[n == 1, one, KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], n], First[wire] -> one]]
];
kron[list_] := If[Length[list] == 1, First[list], KroneckerProduct @@ list];
pauliMatrices[ps_] := MapThread[Times, {ps["Signs"], (kron[(PauliMatrix[Switch[#, {0,0},0, {1,0},1, {1,1},2, {0,1},3]] & /@ #)] & /@ Transpose[ps["Tableau"], {3,2,1}])}];
zeroMatrixQ[m_] := And@@Flatten[Map[TrueQ[FullSimplify[#] == 0]&, m, {2}]];
record["DenseAllGeneratorsAndState", Table[Module[{s = specs[n, 50], initial = PauliStabilizer[n], evolved, u, actual, target, v},
 evolved = initial["ApplyCircuit",s]; u = Fold[unitary[#2,n].#1 &, IdentityMatrix[2^n],s];
 actual = pauliMatrices[evolved]; target = (u.#.ConjugateTranspose[u] &) /@ pauliMatrices[initial];
 v = Normal[evolved["State"]["StateVector"]];
 <|"n" -> n,"GeneratorsExact" -> And@@MapThread[zeroMatrixQ[#1-#2]&, {actual,target}],"Unitary" -> zeroMatrixQ[u.ConjugateTranspose[u]-IdentityMatrix[2^n]],"StateNorm" -> Simplify[Conjugate[v].v],"SignsConsistent" -> And@@(zeroMatrixQ[{#.v-v}]& /@ Take[actual,-n])|>],{n,1,1}]];

Exit[];
