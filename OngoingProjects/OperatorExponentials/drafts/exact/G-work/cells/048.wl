{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b + bd ** a)/2;
FullSimplify[BosonicExpSimplify[Exp[g] ** Exp[al ad - Conjugate[al] a + be bd - Conjugate[be] b] ** Exp[-g]], th > 0]
