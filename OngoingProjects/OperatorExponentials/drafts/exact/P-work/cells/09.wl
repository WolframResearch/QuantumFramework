{a, ad} = FieldVariables[{1}];
dag[x_] := BosonicExpSimplify[Exp[-(g ad - Conjugate[g] a)] ** x ** Exp[g ad - Conjugate[g] a]];
q = Sqrt[hb/2] (a + ad); p = -I Sqrt[hb/2] (a - ad);
Simplify[{dag[a] - a, dag[q] - q, dag[p] - p} /. g -> (Q + I P)/Sqrt[2 hb], Element[{P, Q}, Reals] && hb > 0]
