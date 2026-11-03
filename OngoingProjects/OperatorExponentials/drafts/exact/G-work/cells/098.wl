{a, ad} = FieldVariables[{1}]; nn = ad ** a;
me[x_, op_, y_] := BosonicMatrixElement[{x, y}, op, "Basis" -> "Coherent"];
bb = -2 I gG/ga;
lrho = -I gG (me[u, a + ad, bb] me[bb, 1, v] - me[u, 1, bb] me[bb, a + ad, v]) + ga (me[u, a, bb] me[bb, ad, v] - (me[u, nn, bb] me[bb, 1, v] + me[u, 1, bb] me[bb, nn, v])/2);
sol = alpha[t] /. First[DSolve[{alpha'[t] == -I gG - ga alpha[t]/2, alpha[0] == a0}, alpha[t], t]];
{FullSimplify[lrho, gG > 0 && ga > 0], FullSimplify[sol - (a0 E^(-ga t/2) - 2 I gG/ga (1 - E^(-ga t/2)))], Limit[sol, t -> Infinity, Assumptions -> ga > 0]}
