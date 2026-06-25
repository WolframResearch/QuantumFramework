PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];

(* A1 amplitudes / Born rule *)
psi = QuantumState[{Cos[Pi/6], Sin[Pi/6] Exp[I Pi/4]}];
a1a = psi["StateVector"] // Normal;
a1b = psi["Probabilities"];
a1c = psi["ProbabilitiesList"];

(* A2 normalization & global phase *)
a2a = Norm[psi["StateVector"]];
a2norm = QuantumState[{1, 1}]["Normalized"]["StateVector"] // Normal;
g1 = QuantumState[{1, 0}]; g2 = QuantumState[Exp[I Pi/3] {1, 0}];
a2b = (g1["DensityMatrix"] // Normal) == (g2["DensityMatrix"] // Normal);
a2c = QuantumDistance[g1, g2];

(* A3 Bloch vector *)
a3 = QuantumState["+"]["BlochVector"];
a3y = QuantumState["R"]["BlochVector"];

(* A4 Bloch-angle prep *)
th = Pi/3; ph = Pi/4;
bvec = {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]};
qb = QuantumState["BlochVector"[bvec]];
a4 = Chop[qb["BlochVector"] - bvec];

(* A5 overlap *)
a5a = (QuantumState["0"]["Dagger"] @ QuantumState["+"])["Scalar"];
a5b = Abs[a5a]^2;

(* A6 qudit prep + uniform superposition *)
a6a = QuantumState["0", 3]["StateVector"] // Normal;
a6b = QuantumState["UniformSuperposition", {3}]["StateVector"] // Normal;  (* qutrit uniform *)
a6c = QuantumState["UniformSuperposition"[2]]["StateVector"] // Normal;     (* 2 qubits *)
a6d = QuantumState["+"]["Dimension"];

res = <|
 "A1a-statevector" -> a1a, "A1b-probabilities" -> a1b, "A1c-problist" -> a1c,
 "A2a-norm" -> a2a, "A2norm-example" -> a2norm, "A2b-globalphase-rho-equal" -> a2b, "A2c-distance" -> a2c,
 "A3-bloch-plus" -> a3, "A3y-bloch-R" -> a3y,
 "A4-bloch-recovered-diff" -> a4,
 "A5a-overlap" -> a5a, "A5b-prob" -> a5b,
 "A6a-qutrit0" -> a6a, "A6b-unif-qutrit" -> a6b, "A6c-unif-2qubit" -> a6c, "A6d-dim-plus" -> a6d
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
