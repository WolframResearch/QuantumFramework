PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst];
safe[e_] := Check[e, "ERR"];

(* C1 projective measurement, computational basis: clean accessors *)
qmZ = QuantumMeasurementOperator[QuantumOperator["Z"]];
mPlus = qmZ[QuantumState["+"]];
c1probs = mPlus["ProbabilitiesList"];
c1eigs = QuantumOperator["Z"]["Eigenvalues"];
c1mean = mPlus["Mean"];

(* C2 measure |0> in X basis *)
mX = QuantumMeasurementOperator[QuantumOperator["X"]][QuantumState["0"]];
c2 = mX["ProbabilitiesList"];

(* C3 post-measurement (collapsed) states *)
c3vecs = (#["StateVector"] // Normal) & /@ mPlus["States"] // safe;

(* C4 SIC-POVM (tetrahedron) on |0> : 4 outcomes, sum to 1 *)
povm = QuantumMeasurementOperator["TetrahedronSICPOVM"];
mTet = povm[QuantumState["0"]];
c4probs = mTet["ProbabilitiesList"] // safe;
c4total = Total[mTet["ProbabilitiesList"]] // safe;
c4povmQ = povm["POVMQ"] // safe;

(* C5 finite-shot statistics vs Born rule *)
SeedRandom[42];
c5counts = mPlus["SimulatedCounts", 2000] // safe;
c5theory = N[mPlus["ProbabilitiesList"]];

(* C6 correlated measurement of a Bell pair *)
mBell = QuantumMeasurementOperator[{1, 2}][QuantumState["PhiPlus"]];
c6 = mBell["ProbabilitiesList"] // safe;
c6keys = Keys[mBell["Probabilities"]] /. QuditName[k_, ___] :> k // safe;

res = <|
 "C1-Zon+probs" -> c1probs, "C1-eigs" -> c1eigs, "C1-mean" -> c1mean,
 "C2-Xon0-probs" -> c2,
 "C3-postvecs" -> c3vecs,
 "C4-tet-probs" -> c4probs, "C4-total" -> c4total, "C4-povmQ" -> c4povmQ,
 "C5-counts" -> c5counts, "C5-theory" -> c5theory,
 "C6-bell-jointprobs" -> c6, "C6-keys" -> c6keys
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
