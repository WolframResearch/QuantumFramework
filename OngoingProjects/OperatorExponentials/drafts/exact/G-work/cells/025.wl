{a, ad} = FieldVariables[{1}];
recon[f_] := Module[{g = Exp[al Conjugate[al]] BosonicMatrixElement[{al, al}, f, "Basis" -> "Coherent"] /. Abs[al] -> Sqrt[al Conjugate[al]] /. Conjugate[al] -> u /. al -> v},
  Table[Simplify[(D[g, {u, m}, {v, n}] /. {u -> 0, v -> 0})/Sqrt[m! n!] - BosonicMatrixElement[{m, n}, f]], {m, 0, 4}, {n, 0, 4}]];
{recon[ad ** a ** a + 3 a ** ad ** ad ** ad], recon[Exp[c ad ** a]]}
