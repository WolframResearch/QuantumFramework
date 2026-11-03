{a, ad} = FieldVariables[{1}];
rot = BosonicExpSimplify[Exp[-I t ad ** a] ** a ** Exp[I t ad ** a]];
{rot, Simplify[rot/a - Exp[I t], t > 0]}
