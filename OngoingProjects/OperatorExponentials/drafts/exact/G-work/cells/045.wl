{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b - bd ** a)/(2 I);
FullSimplify[{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[-g] ** b ** Exp[g]]}, th > 0]
