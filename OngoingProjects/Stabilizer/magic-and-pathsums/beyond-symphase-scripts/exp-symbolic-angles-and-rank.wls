PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

Echo["==== P6': symbolic-angle expectation with NON-trivial angle dependence ====", ""];
(* {H, P[theta], H} on |0> : <Z>(theta) = cos(theta/2), a real closed form     *)
(* (QF's P[theta] = diag(1, e^{i theta/2}), the half-angle convention).        *)
v1[th_] := Normal[PauliStabilizer[1][{"H", "P"[th], "H"}]["StateVector"]];
zTheta = FullSimplify[Conjugate[#] . (PauliMatrix[3] . #) &[v1[\[Theta]]], Assumptions -> \[Theta] \[Element] Reals];
Echo[zTheta, "P6' <Z>(theta) for H.P[theta].H |0> (expect Cos[theta/2]):"];
Echo[Table[{th, N[zTheta /. \[Theta] -> th]}, {th, {0, Pi/4, Pi/2, Pi}}], "P6' evaluated by substitution (no re-traversal):"];

Echo["==== P7': does the eager frame have LINEARLY DEPENDENT components? ====", ""];
compMatrix[f_] := Chop[N[Normal[#["State"]["StateVector"]] & /@ f["Stabilizers"]]];
(* (a) T|0>: two components are the SAME physical state -> rank 1 (mergeable). *)
fr0 = PauliStabilizer[1]["T", 1];
Echo[{fr0["Length"], MatrixRank[compMatrix[fr0]]},
  "P7' T|0>: {eager comps, dim of component span} -> 2 comps but rank 1 (mergeable):"];
(* (b) Diagonal T1.T2.CZ|++>: components ORTHOGONAL -> span = 4 (NOT mergeable; *)
(*     compression here is representational, via the phase polynomial).        *)
fdiag = PauliStabilizer[2][{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CZ" -> {1, 2}}];
Echo[{fdiag["Length"], MatrixRank[compMatrix[fdiag]]},
  "P7' diagonal: {eager comps, span} -> 4 comps, span 4 (component-merge does NOT help):"];
(* (c) interleaved H/T on one qubit: 1-qubit state, so its dim is <= 2:        *)
(*     the 2^k components must be massively linearly dependent (mergeable).     *)
fint = PauliStabilizer[1][{"H", "T", "H", "T", "H", "T"}];
Echo[{fint["Length"], MatrixRank[compMatrix[fint]]},
  "P7' interleaved 1-qubit: {eager comps 8, span} -> span<=2 (heavily mergeable):"];

Echo["==== P8: 2-qubit genuine-magic case -- is the eager frame compressible? ====", ""];
(* H between T's on TWO qubits: real multi-qubit magic. Compare eager count to *)
(* the dimension of the component span (an upper-rank witness).                 *)
mrow[gates_] := With[{f = PauliStabilizer[2][gates]},
   {f["Length"], MatrixRank[compMatrix[f]]}];
Echo[{
   "T1,T2 then H1: " -> mrow[{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "H" -> 1}],
   "deeper:        " -> mrow[{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CNOT" -> {1, 2}, "T" -> 1, "H" -> 2, "T" -> 2}]},
  "P8 {eager comps, span dim} (span < comps => mergeable redundancy):"];

Echo["==== P9: exact stabilizer arithmetic -- amplitudes in Z[omega,1/Sqrt2] ====", ""];
(* Every amplitude of any Clifford+T state is an exact element of the ring.    *)
amps = Normal[PauliStabilizer[2][{"H" -> 1, "H" -> 2, "T" -> 1, "CNOT" -> {1, 2}, "T" -> 2}]["StateVector"]];
Echo[RootReduce[amps], "P9 exact amplitudes:"];
Echo[AllTrue[amps, Element[#, Algebraics] &], "P9 all amplitudes exact algebraic numbers:"];

Echo["probe3 done", "STATUS:"];
