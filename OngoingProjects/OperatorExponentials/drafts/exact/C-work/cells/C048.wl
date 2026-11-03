{a, ad} = FieldVariables[{1}];
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
{BosonicExpOrder[Exp[(r ad ** ad - r a ** a)/2], r > 0],
 BosonicExpOrder[Exp[(xi ad ** ad - Conjugate[xi] a ** a)/2]],
 BosonicExpOrder[Exp[xi a1d ** a2d - Conjugate[xi] a1 ** a2]],
 Simplify[{Tan[Sqrt[-xi Conjugate[xi]]]/Sqrt[-xi Conjugate[xi]], Sec[Sqrt[-xi Conjugate[xi]]]} /. xi -> r E^(I th), r > 0 && Element[th, Reals]]}
