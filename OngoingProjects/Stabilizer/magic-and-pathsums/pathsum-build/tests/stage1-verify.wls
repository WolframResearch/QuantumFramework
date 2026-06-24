(* Stage 1 verification: exact span-dimension rank-cap on StabilizerFrame.        *)
(* Run: wolframscript -file stage1-verify.wls                                      *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
dir = DirectoryName[$InputFileName] // ParentDirectory;
Get[FileNameJoin[{dir, "pathsum.wl"}]];
Get[FileNameJoin[{dir, "rankcap.wl"}]];

(* exact equality of two state vectors (the cap is a linear re-expression of the   *)
(* SAME vector, so equality is exact, not just up-to-phase) *)
exactVecEqualQ[a_, b_] := AllTrue[RootReduce[Normal[a] - Normal[b]], PossibleZeroQ];
fidelity[a_, b_] := Module[{x = N[Normal[a]], y = N[Normal[b]]}, Abs[Conjugate[x] . y]/(Norm[x] Norm[y])];

check[label_, f_] := Module[{capped = PathSum`capFrameRank[f]},
    Echo[{
        "eager->capped" -> {f["Length"], capped["Length"]},
        "capped===eager (exact)" -> exactVecEqualQ[capped["StateVector"], f["StateVector"]],
        "numeric |Delta|" -> Chop[Max[Abs[N[RootReduce[Normal[capped["StateVector"]] - Normal[f["StateVector"]]]]]]]
    }, label]];

Echo["==== rank-cap: reproduces the state exactly, caps the component count ====", ""];
check["T|0> (expect 2->1, exact True, |Delta| 0):", PauliStabilizer[1]["T", 1]];
check["interleaved H,T,H,T,H,T 1q (expect 8->2):", PauliStabilizer[1][{"H", "T", "H", "T", "H", "T"}]];
check["diagonal T1.T2.CZ (expect 4->4, no over-collapse):",
    PauliStabilizer[2][{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CZ" -> {1, 2}}]];
check["4-T two-qubit (expect 16->4):",
    PauliStabilizer[2][{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CNOT" -> {1, 2}, "T" -> 1, "H" -> 2, "T" -> 2}]];

Echo["==== convention cross-check: capped frame ~ path-sum amplitude (fidelity 1) ====", ""];
(* StabilizerFrame is correct up to a global phase (kernel contract); compare via  *)
(* numeric fidelity to the path-sum state (which uses the canonical full-angle T).  *)
psT = PathSum`build[1, {{"H", 1}, {"T", 1}}];
psVec = Table[PathSum`amp[psT, 1, {b}], {b, 0, 1}];
Echo[Chop[fidelity[PathSum`capFrameRank[PauliStabilizer[1][{"H", "T"}]]["StateVector"], psVec]],
  "fidelity(capped T|+>, path-sum T|+>) (expect 1):"];

Echo["stage1 done", "STATUS:"];
