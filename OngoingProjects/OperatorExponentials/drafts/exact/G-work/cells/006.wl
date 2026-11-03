{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[a ** nn ** nn ** nn - (nn + 1) ** (nn + 1) ** (nn + 1) ** a, {a, ad}],
 BosonicNormalOrder[nn ** nn ** a - a ** (nn - 1) ** (nn - 1), {a, ad}],
 FullSimplify[BosonicExpSimplify[Exp[-s nn] ** a ** Exp[s nn]], s > 0]}
