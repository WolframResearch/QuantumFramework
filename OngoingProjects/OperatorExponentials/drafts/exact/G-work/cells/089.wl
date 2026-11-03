{b, bd, c, cd} = FieldVariables[{1, 2}];
eta = gg E^(I ph); gen = eta bd ** c - Conjugate[eta] b ** cd;
FullSimplify[{BosonicExpSimplify[Exp[t gen] ** b ** Exp[-t gen]], BosonicExpSimplify[Exp[t gen] ** c ** Exp[-t gen]]}, gg > 0 && t > 0 && ph \[Element] Reals]
