{a, ad, b, bd} = FieldVariables[{1, 2}];
xi = r E^(I th); g = Conjugate[xi] a ** b - xi ad ** bd;
FullSimplify[{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[-g] ** b ** Exp[g]]}, r > 0 && th \[Element] Reals]
