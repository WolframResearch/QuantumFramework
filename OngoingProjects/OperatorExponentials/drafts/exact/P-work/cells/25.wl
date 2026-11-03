{a, ad} = FieldVariables[{1}];
w1 = ComplexExpand[-2 BosonicMatrixElement[{1, 1}, DisplacementOperator[2 (x + I y)]]];
q1 = ComplexExpand[Abs[Exp[-al Conjugate[al]/2] BosonicMatrixElement[{1, 0}, Exp[al ad]]]^2 /. al -> x0 + I y0];
conv = Integrate[2 Exp[-2 ((x0 - x)^2 + (y0 - y)^2)] w1, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> Element[{x0, y0}, Reals]]/Pi;
Simplify[{q1, conv - q1}]
