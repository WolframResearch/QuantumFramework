PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
b[e_] := TimeConstrained[MemoryConstrained[e, 2*^9, "OOM"], 120, "TIMEOUT"];
(* H = sigma_x, psi(t0) = |0>  =>  P0(t) = Cos[t - t0]^2.  Range {t, -1, 1}.
   If the state is pinned at the START (-1): P0(1) = Cos[2]^2 = 0.1732
   If it is pinned at t = 0 instead:          P0(1) = Cos[1]^2 = 0.2919 *)
p0[opts___] := b[Block[{psi = QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], -1, 1}, opts]},
    Chop @ Abs[Normal[psi["StateVector"] /. \[FormalT] -> 1][[1]]]^2]];
Print["@@ expected if pinned at start  ", N[Cos[2]^2]];
Print["@@ expected if pinned at t=0    ", N[Cos[1]^2]];
Print["@@ default                      ", p0[]];
Print["@@ MergeInterpolatingFunctions->False ", p0["MergeInterpolatingFunctions" -> False]];
Print["@@ time-dependent H (t-dependence forces the other branch?)"];
pt = b[Block[{psi = QuantumEvolve[QuantumOperator[(1 + 0 \[FormalT]) PauliMatrix[1]] , QuantumState["0"], {\[FormalT], -1, 1}]},
    Chop @ Abs[Normal[psi["StateVector"] /. \[FormalT] -> 1][[1]]]^2]];
Print["@@   (1+0t) X                   ", pt];
pt2 = b[Block[{psi = QuantumEvolve[QuantumOperator[\[FormalT] PauliMatrix[1]], QuantumState["0"], {\[FormalT], -1, 1}]},
    Chop @ Abs[Normal[psi["StateVector"] /. \[FormalT] -> 1][[1]]]^2]];
Print["@@   t X: start-pinned -> Cos[0]^2=1 (phase integral -1..1 is 0); 0-pinned -> Cos[1/2]^2=", N[Cos[1/2]^2], "  got ", pt2];
