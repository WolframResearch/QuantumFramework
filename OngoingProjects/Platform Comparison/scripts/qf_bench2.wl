PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_,v_]:=Print["RESULT|",l,"|",v];
SetAttributes[tm, HoldFirst];
tm[expr_]:=Module[{t,best=Infinity},Do[ClearSystemCache[];t=First@AbsoluteTiming[expr];best=Min[best,t],{3}];best];

svCircuit[n_]:=QuantumCircuitOperator@Flatten@Table[
   Join[Table["H"->q,{q,n}], Table["CNOT"->{q,q+1},{q,n-1}], Table["RZ"[0.3]->q,{q,n}]],{3}];
pr["sv_circuit_n12_gatecount", svCircuit[12]["GateCount"]];
Do[Module[{qc=svCircuit[n], psi=QuantumState[Table[0,{n}],2]},
    pr["sv_n"<>ToString[n]<>"_schrodinger_ms", 1000*tm[qc[psi,Method->"Schrodinger"]["DensityMatrix"]]]
  ],{n,{12,16}}];

Do[Module[{ops=Import["/tmp/stab_ops_"<>ToString[n]<>".json"]},
    pr["stab_n"<>ToString[n]<>"_m"<>ToString[Length[ops]]<>"_ms",
      1000*tm[Fold[Function[{ps,op},
        Switch[op[[1]],"H",ps["H",op[[2]]+1],"S",ps["S",op[[2]]+1],"CNOT",ps["CNOT",op[[2]]+1,op[[3]]+1]]],
        PauliStabilizer[n], ops]["Signs"]]]
  ],{n,{100,500}}];

Module[{H,L,psi0,sol,rho,finalPop,tt},
  H=QuantumOperator[Pi*PauliMatrix[1]];
  L=QuantumOperator[Sqrt[0.3]*{{0,1},{0,0}}];
  psi0=QuantumState["1"];
  tt=tm[QuantumEvolve[H,{L},psi0,{t,0,2}]];
  sol=QuantumEvolve[H,{L},psi0,{t,0,2}];
  rho=sol["DensityMatrix"];
  finalPop=Re[(rho[[2,2]]) /. t->2.];
  pr["lindblad_ms", 1000*tt];
  pr["lindblad_final_excited_pop", finalPop]
];
