PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 60], " ", ToString[x, InputForm]]);
kerrH[n_] := Module[{a = AnnihilationOperator[n], ad}, ad = a["Dagger"];
   {-0.5 (ad @ a) + 0.25 (ad @ ad @ a @ a) + 1.0 (a + ad), a}];
Do[Module[{H, L, eqs, ret, par, dense, r, t, sas},
   {H, L} = kerrH[n];
   {eqs, ret, par} = QuantumEvolve[H, {L} -> {0.2}, FockState[0, n], {\[FormalT], 0, 10}, "ReturnEquations" -> True];
   sas = Cases[eqs, sa_SparseArray /; SparseArrayQ[sa], Infinity];
   P["N=" <> ToString[n] <> " sparse arrays in the equations / their element kinds",
     {Length[sas], Union[Head /@ Flatten[Normal /@ sas]]}];
   {t, r} = AbsoluteTiming[NDSolveValue[eqs, ret, par]];
   P["N=" <> ToString[n] <> " A: NDSolveValue on sparse equations (s, head)", {t, Head[r]}];
   dense = eqs /. sa_SparseArray ? SparseArrayQ :> Normal[sa];
   P["N=" <> ToString[n] <> "   dense arrays packed?", Developer`PackedArrayQ /@ Cases[dense, _List ? (ArrayQ[#, 2] &), Infinity][[;; UpTo[3]]]];
   {t, r} = AbsoluteTiming[NDSolveValue[dense, ret, par]];
   P["N=" <> ToString[n] <> " B: NDSolveValue after the Normal[] QuantumEvolve applies", {t, Head[r]}];
   {t, r} = AbsoluteTiming[NDSolveValue[dense /. l_List ? (ArrayQ[#, 2, NumericQ] &) :> Developer`ToPackedArray[N[l]], ret, par]];
   P["N=" <> ToString[n] <> " C: same, arrays forced packed", {t, Head[r]}]],
  {n, {30, 35, 40}}];
