{a, ad} = FieldVariables[{1}];
g = (r E^(I th) ad ** ad - r E^(-I th) a ** a)/2;
k0 = (ad ** a + 1/2)/2; kp = ad ** ad/2; km = a ** a/2;
sa = BosonicExpSimplify[Exp[g] ** a ** Exp[-g]];
sk = BosonicExpSimplify[Exp[g] ** k0 ** Exp[-g]];
{sa, sk, Simplify[BosonicNormalOrder[sk - (Cosh[2 r] k0 - Sinh[2 r]/2 (E^(I th) kp + E^(-I th) km)), {a, ad}], r > 0],
 Simplify[sa - (a - E^(I th) Tanh[r] ad)/Sqrt[1 - Tanh[r]^2], r > 0]}
