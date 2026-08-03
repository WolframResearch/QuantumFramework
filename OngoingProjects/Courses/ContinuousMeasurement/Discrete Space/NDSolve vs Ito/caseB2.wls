SetDirectory["/private/tmp/claude-10838/-Users-mohammadb-Documents-GitHub-QuantumMachines-TransmonLab--claude-worktrees-compassionate-banach-e7090b/ae3164c5-9a0b-4ab5-bd0b-7adb86cf987e/scratchpad"];
Get["common.wl"];
ratesB = {3., 0., 1.5, 0., 1., 3.};   (* {gci,gba,gd,gphi,g1,wx}; eta=1 *)
initB = {0., 0., 1.}; tfB = 3.; ntraj = 1000; seed = 2024;
lindZtf = lindbladZfun[ratesB, initB, tfB][tfB];
lind2x2 = lindblad2x2Zfun[ratesB, initB, tfB][tfB];
Echo[{lindZtf, lind2x2, Abs[lindZtf - lind2x2]}, "Lindblad z(tf) affine, 2x2, diff: "];

fnum[x_, d_:4] := If[NumericQ[x], NumberForm[N[x], d], x];
sweep = Table[
  Module[{dt = tfB/ns, emP, emN, it},
    emP = emEnsembleVec[ratesB, initB, tfB, dt, ntraj, seed, True];
    emN = emEnsembleVec[ratesB, initB, tfB, dt, ntraj, seed, False];
    it  = itoEnsemble[ratesB, initB, tfB, dt, ntraj, seed];
    {ns, N@Sqrt[3. dt], Last[emP["meanZ"]], Last[it["meanZ"]],
     Abs[Last[emP["meanZ"]] - lindZtf], Abs[Last[it["meanZ"]] - lindZtf],
     emP["maxr"], TrueQ[emP["diverged"]], emN["maxr"], TrueQ[emN["diverged"]], it["maxr"]}],
  {ns, {60, 150, 300, 600, 1200}}];
hdr = {"nSteps","kick","EMproj_z","Ito_z","|EMp-L|","|Ito-L|","EMp_maxr","EMp_div","EMnoP_maxr","EMnoP_div","Ito_maxr"};
Echo[Grid[Prepend[Map[MapAt[fnum, #, {{3},{4},{5},{6},{7},{9},{11}}] &, sweep], hdr], Frame -> All, Alignment -> Left],
  "Case B sweep (ntraj=1000, Lindblad z(tf)=" <> ToString[lindZtf] <> ", 1/sqrt(N)=" <> ToString[N[1/Sqrt[ntraj]]] <> "): "];
Export["caseB2_result.m", <|"sweep"->sweep, "hdr"->hdr, "lindZtf"->lindZtf,
   "rates"->ratesB, "tf"->tfB, "ntraj"->ntraj, "mcFloor"->N[1/Sqrt[ntraj]]|>];
"caseB2-done"
