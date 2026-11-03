{a, ad} = FieldVariables[{1}]; nn = ad ** a;
expo[e_] := Expand[e /. Power[E, x_] :> x];
{DisplacementOperator[-I eta, Infinity, "Ordering" -> "Normal"], BosonicExpOrder[Exp[-I eta (a + ad)]],
 FullSimplify[ExpToTrig[expo[BosonicExpSimplify[Exp[I nu t nn] ** Exp[-I eta (a + ad)] ** Exp[-I nu t nn]]] + I eta (a E^(-I nu t) + ad E^(I nu t))], nu > 0 && t > 0]}
