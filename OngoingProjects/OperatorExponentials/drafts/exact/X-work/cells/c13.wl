{a, ad} = FieldVariables[{1}];
q = p a ** a + qq ad ** ad + s ad ** a;
lhs = BosonicExpSimplify[Exp[t q] ** a ** Exp[-t q]];
rhs = First[MatrixExp[t {{-s, -2 qq}, {2 p, s}}] . {a, ad}];
{lhs, FullSimplify[lhs - rhs], Chop[Expand[(lhs - rhs) /. {p -> 0.3 + 0.2 I, qq -> -0.7 + 0.1 I, s -> 0.4 - 0.9 I, t -> 1.3}]]}
