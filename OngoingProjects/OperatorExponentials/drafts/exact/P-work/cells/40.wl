{a, ad} = FieldVariables[{1}];
nf = SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"];
{nf, Simplify[BosonicExpOrder[Exp[(rr Exp[-I th] a ** a - rr Exp[I th] ad ** ad)/2]], rr > 0],
 BosonicMatrixElement[{0, 0}, nf]}
