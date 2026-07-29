SetDirectory["/private/tmp/claude-10838/-Users-mohammadb-Documents-GitHub-QuantumMachines-TransmonLab--claude-worktrees-compassionate-banach-e7090b/ae3164c5-9a0b-4ab5-bd0b-7adb86cf987e/scratchpad"];
Get["common.wl"];
ratesDet = {0., 0., 0., 33333.33333333333, 33333.333333333336, 2 Pi*10^6};
initA = N[{Sin[0.001], 0., Cos[0.001]}];
tfA = 1.*^-5; wx = N[2 Pi*10^6];
g1 = ratesDet[[5]]; dep = depOf[ratesDet[[5]], ratesDet[[4]], ratesDet[[3]]];  (* gd=0 here => dep=1/T2 *)
nStab = tfA/(2*((dep + g1)/2)/wx^2);
lindDet = lindbladZfun[ratesDet, initA, tfA][tfA];
Echo[{nStab, lindDet}, "{stability-threshold nSteps, deterministic Lindblad z(tf)}: "];
det = Table[{ns, N[wx*(tfA/ns)], Last@First@emPath[ratesDet, initA, tfA, tfA/ns, False],
             Last@First@emPath[ratesDet, initA, tfA, tfA/ns, True]},
   {ns, {1000, 2000, 4000, 6000, 8000, 16000, 32000}}];
Echo[Grid[Prepend[det, {"nSteps", "wx*dt", "EM noProj z(tf)", "EM proj z(tf)"}], Frame -> All],
  "deterministic EM (explicit Euler on affine ODE), target " <> ToString[lindDet] <> ": "];
"det-done"
