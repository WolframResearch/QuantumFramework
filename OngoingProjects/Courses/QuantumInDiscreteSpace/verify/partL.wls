PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* L1 no-cloning: a cloner would force <psi|phi> = <psi|phi>^2 *)
ov = (QuantumState["0"]["Dagger"] @ QuantumState["+"])["Scalar"];   (* 1/Sqrt2 *)
l1lhs = ov;
l1rhs = ov^2;
l1possible = (l1lhs == l1rhs);                                       (* False unless 0 or 1 *)

(* L2 CHSH Tsirelson bound 2 Sqrt2 on the singlet *)
rhoS = QuantumState["PsiMinus"]["DensityMatrix"] // Normal;
pX = {{0, 1}, {1, 0}}; pZ = {{1, 0}, {0, -1}};
Edir[A_, B_] := Re@Tr[KroneckerProduct[A, B] . rhoS];
B0 = (pZ + pX)/Sqrt[2]; B1 = (pZ - pX)/Sqrt[2];
l2zz = Edir[pZ, pZ];
S = Edir[pZ, B0] + Edir[pZ, B1] + Edir[pX, B0] - Edir[pX, B1] // FullSimplify // safe;
l2N = N[S] // safe;

(* L3 Holevo quantity of {|0> w.p. 1/2, |+> w.p. 1/2} *)
rhoEns = 1/2 QuantumState["0"]["DensityMatrix"] + 1/2 QuantumState["+"]["DensityMatrix"];
l3chi = QuantumState[rhoEns]["VonNeumannEntropy"] // safe;            (* pure parts have S=0 *)
l3N = N[QuantumState[rhoEns]["VonNeumannEntropy"]] // safe;

(* L4 GHZ / Mermin paradox: XYY=YXY=YYX=-1 but XXX=+1 *)
ghz = QuantumState["GHZ"[3]];
ev[str_] := QuantumMeasurementOperator[QuantumOperator[str]][ghz]["Mean"];
l4xxx = ev["XXX"] // safe;
l4xyy = ev["XYY"] // safe;
l4yxy = ev["YXY"] // safe;
l4yyx = ev["YYX"] // safe;
l4product = l4xyy l4yxy l4yyx // safe;     (* local-realism prediction for XXX *)

res = <|
 "L1-overlap" -> l1lhs, "L1-overlap^2" -> l1rhs, "L1-cloning-possible" -> l1possible,
 "L2-singlet-ZZ" -> l2zz, "L2-CHSH-S" -> S, "L2-N" -> l2N, "L2-2Sqrt2" -> N[2 Sqrt[2]],
 "L3-Holevo" -> l3chi, "L3-N" -> l3N,
 "L4-XXX" -> l4xxx, "L4-XYY" -> l4xyy, "L4-YXY" -> l4yxy, "L4-YYX" -> l4yyx, "L4-product(XYY,YXY,YYX)" -> l4product
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
