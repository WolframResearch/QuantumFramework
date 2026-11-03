{a, ad} = FieldVariables[{1}];
qi = QuantumInterferometer[{"S"[1/3 Exp[I/2]] -> 1, "S"[1/4 Exp[-I]] -> 2, "BS"[3/5, 9/10] -> {1, 2}}, 2];
{al, be, c} = First[qi]["Bogoliubov"];
{u, v} = ConjugateTranspose /@ {al, be};
gk = Symbol[First[Names["Wolfram`QuantumFramework`*`gaussianKernel"]]][qi];
{FullSimplify[First[gk] - ArrayFlatten[{{-Inverse[u] . v, Inverse[u]}, {Transpose[Inverse[u]], Conjugate[v] . Inverse[u]}}]], Last[gk] === c, BosonicMatrixElement[{z, w}, Exp[h ad ** a], "Basis" -> "Coherent"]}
