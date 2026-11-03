{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
mu = Cosh[r]; nu = E^(I th) Sinh[r]; dg = al ad - Conjugate[al] a;
step1 = BosonicExpSimplify[Exp[-dg] ** (mu a + nu ad) ** Exp[dg]];
step2 = FullSimplify[BosonicNormalOrder[step1 /. {a -> BosonicExpSimplify[sd ** a ** s], ad -> BosonicExpSimplify[sd ** ad ** s]}, {a, ad}], r > 0 && th \[Element] Reals];
y1 = (a E^(-I th/2) + ad E^(I th/2))/2; y2 = (a E^(-I th/2) - ad E^(I th/2))/(2 I);
{step2, Simplify[BosonicNormalOrder[y1 + I E^(-2 r) y2 - E^(-r - I th/2) (mu a + nu ad), {a, ad}]]}
