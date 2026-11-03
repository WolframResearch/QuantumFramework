{a, ad} = FieldVariables[{1}]; nn = ad ** a;
dg = al ad - Conjugate[al] a;
{BosonicNormalOrder[ad ** nn ** a - nn ** (nn - 1), {a, ad}], BosonicExpSimplify[Exp[-dg] ** ad ** Exp[dg]], FullSimplify[BosonicMatrixElement[{al, al}, a ** ad, "Basis" -> "Coherent"]]}
