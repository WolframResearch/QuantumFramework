SetDirectory["/private/tmp/claude-10838/-Users-mohammadb-Documents-GitHub-QuantumMachines-TransmonLab--claude-worktrees-compassionate-banach-e7090b/ae3164c5-9a0b-4ab5-bd0b-7adb86cf987e/scratchpad"];
Get["common.wl"];
ratesB = {3., 0., 1.5, 0., 1., 3.}; initB = {0., 0., 1.}; tfB = 3.; dtB = tfB/400;
lindZtf = lindbladZfun[ratesB, initB, tfB][tfB];
seOf[assoc_] := With[{zt = assoc["zTrajs"][[All, -1]]}, StandardDeviation[zt]/Sqrt[Length[zt]]];
clt = Table[
  Module[{emP = emEnsembleVec[ratesB, initB, tfB, dtB, n, 2024, True],
          it  = itoEnsemble[ratesB, initB, tfB, dtB, n, 2024]},
    {n, N[1/Sqrt[n]],
     Abs[Last[emP["meanZ"]] - lindZtf], seOf[emP], emP["maxr"],
     Abs[Last[it["meanZ"]]  - lindZtf], seOf[it],  it["maxr"]}],
  {n, {125, 250, 500, 1000, 2000}}];
(* Edge: GBA != 0 (non-optimal LO). Lindblad depends only on gd, so same target. *)
ratesBA = {1.5, 1.5, 1.5, 0., 1., 3.};
lindBA = lindbladZfun[ratesBA, initB, tfB][tfB];
emBA = emEnsembleVec[ratesBA, initB, tfB, dtB, 2000, 7, True];
itBA = itoEnsemble[ratesBA, initB, tfB, dtB, 2000, 7];
Export["caseB_clt_result.m", <|
  "clt" -> clt,
  "hdr" -> {"Ntraj","1/sqrt(N)","|EMp-L|","SE_EM","EMp_maxr","|Ito-L|","SE_Ito","Ito_maxr"},
  "lindZtf" -> lindZtf,
  "edge" -> <|"lindBA_gdonly" -> lindBA, "lindZtf_GBA0" -> lindZtf,
     "emBA_z" -> Last[emBA["meanZ"]], "itBA_z" -> Last[itBA["meanZ"]],
     "emBA_maxr" -> emBA["maxr"], "itBA_maxr" -> itBA["maxr"]|>|>];
Echo["written", "status: "];
