{a, ad} = FieldVariables[{1}];
expo[e_] := Expand[e /. Power[E, x_] :> x];
dd = BosonicExpSimplify[DisplacementOperator[al, Infinity] ** DisplacementOperator[be, Infinity]];
ddRev = BosonicExpSimplify[DisplacementOperator[be, Infinity] ** DisplacementOperator[al, Infinity]];
{Simplify[expo[dd] - ((al Conjugate[be] - Conjugate[al] be)/2 + (al + be) ad - (Conjugate[al] + Conjugate[be]) a)],
 Simplify[expo[dd] - expo[ddRev] - (al Conjugate[be] - Conjugate[al] be)],
 Simplify[ComplexExpand[(al Conjugate[be] - Conjugate[al] be)/2 - I Im[al Conjugate[be]] /. {al -> x1 + I y1, be -> x2 + I y2}]]}
