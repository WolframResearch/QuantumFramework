(* X-009 bosons, general two-mode case from QF's exact Bogoliubov data (finite 2N x 2N matrices, no Fock space) *)
emit[k_, v_] := WriteString["stdout", k, " -> ", ToString[InputForm[v]], "\n"];
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
qi = QuantumInterferometer[{"S"[1/3 Exp[I/2]] -> 1, "S"[1/4 Exp[-I]] -> 2, "BS"[3/5, 9/10] -> {1, 2}, "S"[1/5 Exp[2 I]] -> 1}, 2];
{al, be, c} = First[qi]["Bogoliubov"];
{u, v} = N[ConjugateTranspose /@ {al, be}, 40]; c = N[c, 40];
z = Inverse[u] . v; w = Conjugate[v] . Inverse[u];
emit["Z, W symmetric", {Max[Abs[z - Transpose[z]]], Max[Abs[w - Transpose[w]]]}];
(* a^dag conjugation through the three factors must give conj(u) a^dag + conj(v) a: needs conj(u) == Transpose[Inverse[u]] + conj(v).Inverse[u].v *)
emit["conj(u) - (u^-1)^T - conj(v) u^-1 v", Max[Abs[Conjugate[u] - Transpose[Inverse[u]] - Conjugate[v] . Inverse[u] . v]]];
emit["|c| - det(u u^dag)^(-1/4)", Abs[c] - Det[u . ConjugateTranspose[u]]^(-1/4)];
(* QF conjugation of the two outer factors, symbolic symmetric Z, W *)
xz = -(z11 ad1 ** ad1 + 2 z12 ad1 ** ad2 + z22 ad2 ** ad2)/2;
yw = (w11 a1 ** a1 + 2 w12 a1 ** a2 + w22 a2 ** a2)/2;
emit["e^X a e^-X", {BosonicExpSimplify[Exp[xz] ** a1 ** Exp[-xz]], BosonicExpSimplify[Exp[xz] ** a2 ** Exp[-xz]]}];
emit["e^Y ad e^-Y", {BosonicExpSimplify[Exp[yw] ** ad1 ** Exp[-yw]], BosonicExpSimplify[Exp[yw] ** ad2 ** Exp[-yw]]}];
