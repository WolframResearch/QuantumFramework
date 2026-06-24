(* Run the path-sum benchmarks. Output is markdown-table-ready.                      *)
(* wolframscript -file run-benchmarks.wls                                             *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
dir = ParentDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{dir, "pathsum.wl"}]];
Get[FileNameJoin[{dir, "harness.wl"}]];
Get[FileNameJoin[{dir, "elimination.wl"}]];
Get[FileNameJoin[{dir, "diagonal.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "bench.wl"}]];

Echo[$MachineName <> " | " <> ToString[$ProcessorCount] <> " cores | WL " <> ToString[$VersionNumber], "MACHINE:"];

Echo["==== M1: diagonal single-amplitude scaling (diagAmp vs QF dense, seconds) ====", ""];
Echo[Prepend[Table[PathSum`benchDiag[n], {n, {4, 6, 8, 10, 12, 14}}],
    {"n", "#gates", "diagAmp s", "QF-dense s"}], "M1:"];

Echo["==== M2: engine compression, Clifford ladder (residual -> 0, reduce time) ====", ""];
Echo[Prepend[Table[PathSum`benchLadder[6, d], {d, {1, 2, 3, 4, 6, 8}}],
    {"n", "depth", "vars0", "varsR", "internalResidual", "steps", "reduce s"}], "M2:"];

Echo["==== M3: parametric capability (one symbolic traversal vs dense per angle) ====", ""];
Echo[PathSum`benchParam[Range[0, 2 Pi, 2 Pi/16]], "M3 {#angles, engine-once s, dense-per-angle s, closed form}:"];

Echo["==== M4: generic-magic HONEST LOSS (engine vs dense single amplitude) ====", ""];
Echo[Prepend[Table[PathSum`benchMagic[3, d, 100 + d], {d, {8, 10, 12, 14}}],
    {"n", "depth", "internalResidual", "engineAmp s", "QF-dense s"}], "M4:"];

Echo["==== M5: representational, eager StabilizerFrame 2^{#T} vs compressor rank ====", ""];
Echo[Prepend[Table[PathSum`benchFrame[3, k], {k, {1, 2, 3, 4, 5, 6}}],
    {"n", "#T", "frame components", "compressor rank"}], "M5:"];

Echo["benchmarks done", "STATUS:"];
