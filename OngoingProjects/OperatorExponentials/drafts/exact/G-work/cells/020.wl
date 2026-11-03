{a, ad} = FieldVariables[{1}];
ov = BosonicMatrixElement[{be, al}, 1, "Basis" -> "Coherent"];
{ov, Simplify[ComplexExpand[ov Conjugate[ov] - Exp[-Abs[al - be]^2] /. {al -> x1 + I y1, be -> x2 + I y2}]],
 Integrate[2 rho^(2 m + 1) Exp[-rho^2]/m!, {rho, 0, Infinity}, Assumptions -> m >= 0],
 Simplify[Integrate[Exp[I k phi], {phi, 0, 2 Pi}], k \[Element] Integers && k != 0]}
