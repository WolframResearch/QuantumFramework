{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
an = {a1, a2}; cr = {ad1, ad2}; d = KroneckerDelta; cm[x_, y_] := x ** y - y ** x;
e[i_, j_] := cr[[i]] ** an[[j]] + d[i, j]/2; f[i_, j_] := an[[i]] ** an[[j]]; g[i_, j_] := cr[[i]] ** cr[[j]];
res = Flatten[Table[BosonicNormalOrder[#, {a1, ad1, a2, ad2}] & /@ {
     cm[e[i, j], e[k, l]] - (d[j, k] e[i, l] - d[i, l] e[k, j]),
     cm[e[i, j], f[k, l]] - (-d[i, k] f[j, l] - d[i, l] f[k, j]),
     cm[e[i, j], g[k, l]] - (d[j, k] g[i, l] + d[j, l] g[k, i]),
     cm[f[i, j], g[k, l]] - (d[j, k] e[l, i] + d[j, l] e[k, i] + d[i, k] e[l, j] + d[i, l] e[k, j]),
     cm[f[i, j], f[k, l]], cm[g[i, j], g[k, l]]}, {i, 2}, {j, 2}, {k, 2}, {l, 2}]];
{Length[res], Union[Expand[res]]}
