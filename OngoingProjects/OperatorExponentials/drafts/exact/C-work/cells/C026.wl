{a, ad} = FieldVariables[{1}];
one = BosonicExpSimplify[Exp[(r E^(I th) ad ** ad - r E^(-I th) a ** a)/2] ** a **
    Exp[-(r E^(I th) ad ** ad - r E^(-I th) a ** a)/2]];
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
b = {{1/4 + I/10, 3/20 - I/20}, {3/20 - I/20, -1/10 + I/5}};
cr = {a1d, a2d}; an = {a1, a2}; ops = Join[an, cr];
g = (Sum[b[[i, j]] cr[[i]] ** cr[[j]], {i, 2}, {j, 2}] - Sum[Conjugate[b[[i, j]]] an[[i]] ** an[[j]], {i, 2}, {j, 2}])/2;
adg = Table[Coefficient[BosonicNormalOrder[g ** x - x ** g, ops], y], {x, ops}, {y, ops}];
mb = MatrixPower[N[ConjugateTranspose[b] . b, 40], 1/2];
w = N[b, 40] . Inverse[mb] . MatrixFunction[Tanh, mb];
blk = ArrayFlatten[{{MatrixPower[IdentityMatrix[2] - w . ConjugateTranspose[w], -1/2], -w . MatrixPower[IdentityMatrix[2] - ConjugateTranspose[w] . w, -1/2]},
    {-MatrixPower[IdentityMatrix[2] - ConjugateTranspose[w] . w, -1/2] . ConjugateTranspose[w], MatrixPower[IdentityMatrix[2] - ConjugateTranspose[w] . w, -1/2]}}];
{one, adg, Chop[Max[Abs[MatrixExp[N[adg, 40]] - blk]], 10^-30]}
