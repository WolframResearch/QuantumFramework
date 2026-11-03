{a, ad} = FieldVariables[{1}];
kp = ad ** ad/2; km = a ** a/2;
gen = rr Exp[I ph] kp - rr Exp[-I ph] km;
{Simplify[BosonicExpOrder[Exp[gen]], rr > 0], Simplify[BosonicExpOrder[Exp[gen], "Ordering" -> "Antinormal"], rr > 0]}
