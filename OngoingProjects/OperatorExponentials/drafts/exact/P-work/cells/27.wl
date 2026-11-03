{a, ad} = FieldVariables[{1}];
chiN = Exp[g Conjugate[be] - gb be];
{BosonicMatrixElement[{be, be}, ad ** ad ** a, "Basis" -> "Coherent"] /. Abs[be]^2 -> be Conjugate[be],
 -D[chiN, {g, 2}, {gb, 1}] /. {g -> 0, gb -> 0}}
