{a, ad} = FieldVariables[{1}];
lhs = BosonicExpSimplify[Exp[al ad - Conjugate[al] a] ** Exp[be ad - Conjugate[be] a]];
{lhs, Simplify[lhs - Exp[(al Conjugate[be] - Conjugate[al] be)/2] Exp[(al + be) ad - (Conjugate[al] + Conjugate[be]) a]]}
