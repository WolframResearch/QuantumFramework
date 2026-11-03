{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
{Simplify[BosonicExpOrder[Exp[rr Exp[I ph] ad1 ** ad2 - rr Exp[-I ph] a1 ** a2]], rr > 0],
 BosonicMatrixElement[{{3, 3}, {0, 0}}, (ad1 ** ad2) ** (ad1 ** ad2) ** (ad1 ** ad2)]/3!}
