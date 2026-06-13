PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
SetAttributes[tm, HoldFirst];
tm[e_] := Module[{b=Infinity,t}, Do[ClearSystemCache[]; t=First@AbsoluteTiming[e]; b=Min[b,t], {3}]; 1000. b];
toSpec[ops_]:=Map[Switch[#[[1]],"H","H"->#[[2]]+1,"S","S"->#[[2]]+1,"CNOT","CNOT"->{#[[2]]+1,#[[3]]+1}]&,ops];
SeedRandom[3];
mk[n_] := Fold[#1["H",RandomInteger[{1,n}]]["S",RandomInteger[{1,n}]]["CNOT",#2,Mod[#2,n]+1]&, PauliStabilizer[n], RandomInteger[{1,n},n]];
PauliStabilizer[4]["ApplyCircuit", {"H"->1}];  (* warm compile *)
Do[Module[{ops=Import["/tmp/stab_ops_"<>ToString[n]<>".json"], specs},
   specs=toSpec[ops];
   pr["gates_n"<>ToString[n]<>"_AC_ms", tm[PauliStabilizer[n]["ApplyCircuit", specs]["Signs"]]];
   pr["gates_n"<>ToString[n]<>"_pergate_ms", tm[Fold[Function[{p,op},Switch[op[[1]],"H",p["H",op[[2]]+1],"S",p["S",op[[2]]+1],"CNOT",p["CNOT",op[[2]]+1,op[[3]]+1]]],PauliStabilizer[n],ops]["Signs"]]];
], {n,{100,500,1000}}];
Do[pr["ctor_n"<>ToString[n]<>"_ms", tm[PauliStabilizer[n]]], {n,{100,500,1000}}];
Do[Module[{ps=mk[n], zstr=StringJoin[ConstantArray["Z",n]]},
   pr["measZ_n"<>ToString[n]<>"_ms", tm[ps["M",1]]];
   pr["measPauli_n"<>ToString[n]<>"_ms", tm[ps["M",zstr]]];
   pr["stabs_n"<>ToString[n]<>"_ms", tm[ps["Stabilizers"]]];
], {n,{100,500}}];
Module[{a=mk[100], b=mk[100]}, pr["compose_n100_ms", tm[a[b]]]];
Print["DONE"];
