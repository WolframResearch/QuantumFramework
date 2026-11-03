{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
ex[e_] := e /. Power[E, x_] :> Power[E, Expand[x]];
s = be Conjugate[al] + g Conjugate[al] + g Conjugate[be];
{ex[BosonicExpSimplify[d[al] ** d[be]]/BosonicExpSimplify[d[be] ** d[al]]],
 ex[BosonicExpSimplify[d[al] ** d[be] ** d[-al] ** d[-be]]],
 Simplify[ex[BosonicExpSimplify[d[g] ** d[be] ** d[al] ** d[-(al + be + g)]]]/Exp[(s - Conjugate[s])/2]]}
