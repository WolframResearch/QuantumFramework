{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
{BosonicNormalOrder[a1 ** ad2 - ad2 ** a1, {a1, ad1, a2, ad2}], BosonicNormalOrder[a1 ** a2 - a2 ** a1, {a1, ad1, a2, ad2}],
 BosonicNormalOrder[a1 ** ad1 - ad1 ** a1, {a1, ad1, a2, ad2}], BosonicMatrixElement[{{2, 1}, {0, 0}}, ad1 ** ad1 ** ad2/Sqrt[2! 1!]]}
