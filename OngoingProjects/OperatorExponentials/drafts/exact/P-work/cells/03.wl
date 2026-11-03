{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity] /. Conjugate[s_Plus] :> Map[Conjugate, s];
{Simplify[BosonicExpSimplify[d[al] ** d[be]]/(Exp[(al Conjugate[be] - Conjugate[al] be)/2] d[al + be])],
 BosonicExpSimplify[d[al] ** d[-al]]}
