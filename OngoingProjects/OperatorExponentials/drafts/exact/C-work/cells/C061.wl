{a, ad} = FieldVariables[{1}];
d[x_, xc_] := Exp[x ad - xc a];
gen = BosonicExpSimplify[d[-al, -alc] ** d[al + e be, alc + e bec]] /. E^x_ :> x;
{gen, Simplify[D[gen, e] - (be ad - bec a + (alc be - al bec)/2)]}
