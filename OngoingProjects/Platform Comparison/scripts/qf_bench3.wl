PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_,v_]:=Print["RESULT|",l,"|",v];
SetAttributes[tm, HoldFirst];
tm[expr_]:=Module[{t,best=Infinity},Do[ClearSystemCache[];t=First@AbsoluteTiming[expr];best=Min[best,t],{2}];best];
svCircuit[n_]:=QuantumCircuitOperator@Flatten@Table[
   Join[Table["H"->q,{q,n}], Table["CNOT"->{q,q+1},{q,n-1}], Table["RZ"[0.3]->q,{q,n}]],{3}];
Do[Module[{qc=svCircuit[n], psi=QuantumState[Table[0,{n}],2]},
    pr["sv_n"<>ToString[n]<>"_default_TN_ms", 1000*tm[qc[psi]["DensityMatrix"]]]
  ],{n,{12,16}}];
