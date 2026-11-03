{a, ad} = FieldVariables[{1}];
qsym = BosonicMatrixElement[{al, al}, a ** ad, "Basis" -> "Coherent"];
smooth = Integrate[Exp[-((x - x0)^2 + (y - y0)^2)] (x^2 + y^2), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/Pi;
Simplify[ComplexExpand[qsym /. al -> x0 + I y0] - smooth]
