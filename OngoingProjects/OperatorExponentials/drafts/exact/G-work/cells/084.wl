{a, ad, b, bd, c, cd} = FieldVariables[{1, 2, 3}]; v = {a, ad, b, bd, c, cd};
g = I Pi (ad ** b + a ** bd)/4;
nbT = BosonicNormalOrder[BosonicExpSimplify[Exp[-g] ** bd ** Exp[g]] ** BosonicExpSimplify[Exp[-g] ** b ** Exp[g]], v];
gen1 = I eta (ad ** a + bd ** b)/2; gen2 = eta (ad ** b - bd ** a)/2;
{Simplify[nbT], Simplify[BosonicNormalOrder[I eta nbT ** cd ** c - (gen1 + gen2) ** cd ** c, v]],
 BosonicNormalOrder[(gen1 ** cd ** c) ** (gen2 ** cd ** c) - (gen2 ** cd ** c) ** (gen1 ** cd ** c), v]}
