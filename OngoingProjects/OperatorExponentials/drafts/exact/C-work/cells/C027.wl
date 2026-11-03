{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"],
 BosonicMatrixElement[{0, 0}, SqueezeOperator[xi]],
 Simplify[BosonicMatrixElement[{2, 0}, SqueezeOperator[xi]] - (1 - Abs[w]^2)^(1/4) w/Sqrt[2] /. w -> -Exp[I Arg[xi]] Tanh[Abs[xi]]]}
