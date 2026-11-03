{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
x = -(z11 ad1 ** ad1 + 2 z12 ad1 ** ad2 + z22 ad2 ** ad2)/2;
{BosonicExpSimplify[Exp[x] ** a1 ** Exp[-x]], BosonicExpSimplify[Exp[x] ** a2 ** Exp[-x]]}
