{a, ad} = FieldVariables[{1}];
ada = BosonicExpSimplify[Exp[s ad ** a] ** a ** Exp[-s ad ** a]];
{ada, Simplify[ada/a - Exp[-s], Element[s, Reals]]}
