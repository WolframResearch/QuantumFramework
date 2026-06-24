(* Stage 3 verification: the path-sum rewrite engine (full rule set + brute-force     *)
(* completion), degree invariant, confluence-by-replay, structured telescoping, the    *)
(* Inactive[Sum] carrier, and the magic-residue / fallback behaviour.                   *)
(* Run: wolframscript -file stage3-verify.wls                                           *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
dir = DirectoryName[$InputFileName] // ParentDirectory;
Get[FileNameJoin[{dir, "pathsum.wl"}]];
Get[FileNameJoin[{dir, "harness.wl"}]];
Get[FileNameJoin[{dir, "elimination.wl"}]];

hl[n_] := Table[{"H", q}, {q, n}];

Echo["==== A. exact differential: engine === brute-force build === QF Schrodinger ====", ""];
(* a battery spanning Clifford, Clifford+T, interleaved magic, diagonal, mixed:      *)
(* {name, qubits, gates}. *)
battery = {
  {"Clifford H(SH)^3 1q", 1, {{"H",1},{"S",1},{"H",1},{"S",1},{"H",1},{"S",1},{"H",1}}},
  {"Clifford 2q HSCZ", 2, {{"H",1},{"H",2},{"S",1},{"CZ",1,2},{"H",1},{"H",2},{"CZ",1,2},{"H",1}}},
  {"magic H(TH)^4 1q", 1, Join[{{"H",1}}, Flatten[Table[{{"T",1},{"H",1}}, 4], 1]]},
  {"mixed 3q", 3, {{"H",1},{"H",2},{"H",3},{"S",1},{"CZ",1,2},{"T",2},{"H",2},{"CNOT",1,3},{"CZ",2,3},{"H",1}}},
  {"diagonal-on-+ 2q", 2, Join[hl[2], {{"T",1},{"S",2},{"CZ",1,2},{"CNOT",1,2}}]}
};
Echo[{#[[1]], PathSum`engineMatchesQ[#[[2]], #[[3]]]} & /@ battery, "A engineMatchesQ per circuit (expect all True):"];

Echo["==== A2. randomized {H,T,S,Z,CZ,CCZ,CNOT} (the fallback gate: reduce + brute) ====", ""];
SeedRandom[314159];
randTrials = Table[Module[{n = RandomInteger[{2, 4}], d = RandomInteger[{8, 18}], c},
    c = PathSum`randomCircuit[n, d]; {n, d, PathSum`engineMatchesQ[n, c]}], {10}];
Echo[{Count[randTrials[[All, 3]], True], Length[randTrials]}, "A2 passed / total (expect all):"];

Echo["==== B. degree-<=3 invariant holds after reduction ====", ""];
SeedRandom[2718];
degOk = Table[Module[{n = RandomInteger[{2, 4}], c}, c = PathSum`randomCircuit[n, RandomInteger[{8, 16}]];
    PathSum`pdeg[PathSum`engineReduce[PathSum`build[n, c]]["phase"]] <= 3], {15}];
Echo[{Count[degOk, True], Length[degOk]}, "B degree<=3 after reduce (expect all):"];

Echo["==== C. confluence by replay: amplitude invariant under shuffled rule order ====", ""];
SeedRandom[11];
confl = Table[Module[{n = RandomInteger[{2, 4}], c}, c = PathSum`randomCircuit[n, RandomInteger[{8, 14}]];
    PathSum`confluentQ[n, c, {1, 2, 3, 4, 5}]], {8}];
Echo[{Count[confl, True], Length[confl]}, "C confluent (shuffled === canonical) passed / total:"];

Echo["==== D. structured telescoping: Clifford H/S/CZ ladder, internal residual -> ? ====", ""];
czchain[n_] := Table[{"CZ", q, q + 1}, {q, n - 1}];
slayer[n_] := Table[{"S", q}, {q, n}];
cliffLadder[n_, d_] := Join[hl[n], Flatten[Table[Join[slayer[n], czchain[n], hl[n]], {d}], 1]];
ladderRow[n_, d_] := Module[{ps = PathSum`build[n, cliffLadder[n, d]], r},
    r = PathSum`engineReduce[ps];
    {n, d, Length[ps["vars"]], Length[r["vars"]], Length[PathSum`internalResidual[r, n]], r["steps"]}];
Echo[Prepend[Table[ladderRow[4, d], {d, 1, 5}],
    {"n", "depth", "vars(before)", "vars(after)", "internal-residual", "steps"}],
    "D Clifford ladder (n=4): internal residual should fall to 0, steps ~ linear:"];
Echo[PathSum`engineMatchesQ[3, cliffLadder[3, 2]], "D correctness on cliffLadder[3,2] (exact vs QF):"];

Echo["==== D2. magic residue: H(TH)^k leaves exactly the T-blocked variables ====", ""];
magicRow[k_] := Module[{c = Join[{{"H",1}}, Flatten[Table[{{"T",1},{"H",1}}, k], 1]], r},
    r = PathSum`engineReduce[PathSum`build[1, c]];
    {k, Length[r["vars"]], Length[PathSum`internalResidual[r, 1]]}];
Echo[Prepend[Table[magicRow[k], {k, 1, 5}], {"#T", "vars(after)", "internal-residual (the magic)"}],
    "D2 magic survives (internal residual grows with #T):"];

Echo["==== E. Inactive[Sum] formal carrier: inert, Activate === engine amplitude ====", ""];
cE = {{"H",1},{"H",2},{"S",1},{"CZ",1,2},{"H",1},{"H",2},{"CZ",1,2},{"H",1}};
psE = PathSum`engineReduce[PathSum`build[2, cE]];
fsE = PathSum`formalSum[psE, 2, {1, 0}];
Echo[{! FreeQ[fsE, Inactive], Head /@ {PathSum`PS}}, "E formal sum carries Inactive[Sum]; PS head exists:"];
Echo[{MemberQ[Attributes[PathSum`PS], HoldAll], MemberQ[Attributes[PathSum`PS], NonThreadable]},
    "E PS attributes {HoldAll, NonThreadable}:"];
Echo[PathSum`exactEqualQ[PathSum`activate[fsE], PathSum`amp[psE, 2, {1, 0}]], "E Activate[formalSum] === engine amp:"];
Echo[PathSum`exactEqualQ[PathSum`activate[fsE], PathSum`qfAmp[2, cE, {1, 0}]], "E Activate[formalSum] === QF amp:"];

Echo["stage3 done", "STATUS:"];
