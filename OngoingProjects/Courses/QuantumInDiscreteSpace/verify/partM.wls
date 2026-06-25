PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
SetAttributes[rec, HoldRest];
rec[id_, e_] := Module[{v = e}, WriteString[$Output, id <> " = " <> ToString[InputForm[v]] <> "\n"]; v];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];
SetFockSpaceSize[12];
am = AnnihilationOperator[]["Matrix"] // Normal;
adm = ConjugateTranspose[am];
nm = adm . am;
fock[k_] := FockState[k]["StateVector"] // Normal;

(* M1 Fock states + ladder operators *)
rec["M1-a-dim", Dimensions[am]];
rec["M1-fock2", fock[2]];
rec["M1-a|2>=Sqrt2|1>", Chop[am . fock[2] - Sqrt[2] fock[1]]];
rec["M1-N|3>=3|3>", Chop[nm . fock[3] - 3 fock[3]]];

(* M2 coherent state |alpha> = D(alpha)|0>: eigenstate of a, mean photon |alpha|^2 *)
al = 1.2;
cohv = DisplacementOperator[al][FockState[0]]["StateVector"] // Normal;
rec["M2-eigencheck-low", Chop[(am . cohv - al cohv)][[1 ;; 6]]];
rec["M2-meanN", Chop[Re[Conjugate[cohv] . nm . cohv]]];
rec["M2-|alpha|^2", al^2];

(* M3 coherent photon distribution is Poisson *)
rec["M3-photon-probs", Chop[Abs[cohv]^2][[1 ;; 6]]];
rec["M3-poisson-pdf", Table[N@PDF[PoissonDistribution[al^2], n], {n, 0, 5}]];

(* M4 squeezing pushes a quadrature variance below the vacuum value *)
quad = QuadratureOperators[][[1]];
rec["M4-vacuum-var", N@OperatorVariance[FockState[0], quad] // safe];
rec["M4-squeezed-var", N@OperatorVariance[SqueezeOperator[1.0][FockState[0]], quad] // safe];

(* M5 g^(2)(0) photon statistics: thermal (bunching), coherent (1), Fock (antibunching) *)
g2[v_] := Re[Conjugate[v] . (adm . adm . am . am) . v]/Re[Conjugate[v] . nm . v]^2;
rec["M5-g2-thermal", N@G2Coherence[ThermalState[1.0]] // safe];
rec["M5-g2-coherent", g2[cohv]];
rec["M5-g2-fock1", N@G2Coherence[FockState[1]] // safe];
