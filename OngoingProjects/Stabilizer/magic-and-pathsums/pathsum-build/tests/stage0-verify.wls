(* Stage 0 verification: native-WL refactor + exact harness + convention.       *)
(* Run: wolframscript -file stage0-verify.wls                                     *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
dir = DirectoryName[$InputFileName] // ParentDirectory;
Get[FileNameJoin[{dir, "pathsum.wl"}]];
Get[FileNameJoin[{dir, "harness.wl"}]];

(* demo circuit constructors (from the frozen POC) *)
cliffChain[k_] := Join[{{"H", 1}}, Catenate[Table[{{"Z", 1}, {"H", 1}}, k]]];
magicChain[k_] := Join[{{"H", 1}}, Catenate[Table[{{"T", 1}, {"H", 1}}, k]]];
deep2q[k_]     := Join[{{"H", 1}, {"H", 2}}, Catenate[Table[{{"CZ", 1, 2}, {"H", 1}, {"H", 2}}, k]]];
mixed          = {{"H", 1}, {"Z", 1}, {"H", 1}, {"CZ", 1, 2}, {"H", 1}, {"T", 1}, {"H", 1}, {"H", 2}};

Echo["==== convention invariant: T is full-angle diag(1, e^{i pi/4}) ====", ""];
tps = PathSum`build[1, {{"H", 1}, {"T", 1}}];
Echo[{PathSum`exactEqualQ[PathSum`amp[tps, 1, {0}], PathSum`qfAmp[1, {{"H", 1}, {"T", 1}}, {0}]],
      PathSum`exactEqualQ[PathSum`amp[tps, 1, {1}], PathSum`qfAmp[1, {{"H", 1}, {"T", 1}}, {1}]]},
  "T|+> amplitudes match QF exactly (expect {True,True}):"];

Echo["==== regression demos: {vars before, after, reduced===QF (exact)} ====", ""];
demo[label_, n_, gates_] := With[{ps = PathSum`build[n, gates], r = PathSum`reduce[PathSum`build[n, gates]]},
    Echo[{Length[ps["vars"]], Length[r["vars"]], PathSum`reducedMatchesQ[n, gates],
          IntegerQ[Log2[r["scale"]]]}, label]];
demo["Clifford H(ZH)^8 1q (expect 9,1,True,True):", 1, cliffChain[8]];
demo["Clifford 2q deep    (expect 10,2,True,True):", 2, deep2q[4]];
demo["Magic    H(TH)^8 1q (expect 9,9,True,True):", 1, magicChain[8]];
demo["Mixed Clifford+1T   (expect 5,3,True,True):", 2, mixed];

Echo["==== builder exactness (un-reduced brute force vs QF) ====", ""];
Echo[{PathSum`builderMatchesQ[1, magicChain[6]], PathSum`builderMatchesQ[2, deep2q[3]],
      PathSum`builderMatchesQ[3, {{"H",1},{"H",2},{"H",3},{"T",1},{"CCZ",1,2,3},{"H",2}}]},
  "builder === QF on 3 circuits (expect {True,True,True}):"];

Echo["==== randomized differential suite (reduced === QF, exact) ====", ""];
SeedRandom[20260622];
specs = Join[
    Table[{2, RandomInteger[{4, 10}]}, {4}],
    Table[{3, RandomInteger[{4, 8}]}, {3}],
    Table[{4, RandomInteger[{3, 6}]}, {2}]
];
results = Table[
    Module[{n = sp[[1]], d = sp[[2]], c}, c = PathSum`randomCircuit[n, d];
        {n, d, PathSum`reducedMatchesQ[n, c]}],
    {sp, specs}];
Echo[results, "per-circuit {n, depth, reduced===QF}:"];
Echo[{Count[results[[All, 3]], True], Length[results]}, "passed / total (expect all):"];

Echo["stage0 done", "STATUS:"];
