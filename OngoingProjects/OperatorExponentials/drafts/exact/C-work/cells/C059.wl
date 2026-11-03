{a, ad} = FieldVariables[{1}];
d[x_] := Exp[x ad - Conjugate[x] a];
cpx = {f -> f1 + I f2, g -> g1 + I g2};
prod = BosonicExpSimplify[d[f] ** d[g]];
ov = BosonicMatrixElement[{f, g}, 1, "Basis" -> "Coherent"];
{prod, ov, Simplify[ComplexExpand[{(f Conjugate[g] - g Conjugate[f])/2 + I Im[Conjugate[f] g],
    (ov /. E^x_ :> x) - (-Abs[f - g]^2/2 + I Im[Conjugate[f] g])} /. cpx]]}
