{a, ad} = FieldVariables[{1}];
vac = BosonicMatrixElement[{0, 0}, DisplacementOperator[al]];
ov = BosonicMatrixElement[{al, be}, 1, "Basis" -> "Coherent"];
toZ = {al -> (q1 + I p1)/Sqrt[2], be -> (q2 + I p2)/Sqrt[2]};
{vac, ov, Simplify[ComplexExpand[{vac, ov} /. E^x_ :> x /. toZ] -
   {-(q1^2 + p1^2)/4, I (q1 p2 - p1 q2)/2 - ((q1 - q2)^2 + (p1 - p2)^2)/4}]}
