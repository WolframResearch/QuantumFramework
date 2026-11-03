{a, ad} = FieldVariables[{1}];
d[x_] := Exp[x ad - Conjugate[x] a];
{BosonicExpSimplify[d[al] ** d[be]], BosonicExpSimplify[d[be] ** d[al]], BosonicExpSimplify[d[al] ** d[-al]],
 ComplexExpand[(al Conjugate[be] - be Conjugate[al])/2 /. {al -> (q1 + I p1)/Sqrt[2], be -> (q2 + I p2)/Sqrt[2]}]}
