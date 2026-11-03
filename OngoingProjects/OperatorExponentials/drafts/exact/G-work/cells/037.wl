{a, ad} = FieldVariables[{1}]; nn = ad ** a;
sp = {{0, 1}, {0, 0}}; sm = Transpose[sp]; s3 = {{1, 0}, {0, -1}}; id = IdentityMatrix[2]; pe = {{1, 0}, {0, 0}};
ncm[1, y_] := y; ncm[x_, 1] := x; ncm[x_, y_] := BosonicNormalOrder[x ** y, {a, ad}];
prod[x_, y_] := Total[Flatten[Outer[#1[[1]] #2[[1]] ncm[#1[[2]], #2[[2]]] (#1[[3]] . #2[[3]]) &, x, y, 1], 1]];
H = {{w, nn, id}, {w0/2, 1, s3}, {lam, a, sp}, {lam, ad, sm}};
M = {{1, nn, id}, {1, 1, pe}};
{sp . sm - sm . sp == s3, s3 . sp - sp . s3 == 2 sp, s3 . sm - sm . s3 == -2 sm, Simplify[Expand[prod[H, M] - prod[M, H]]]}
