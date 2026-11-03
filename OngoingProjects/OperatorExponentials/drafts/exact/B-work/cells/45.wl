{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d}; an = {a1, a2}; cr = {a1d, a2d};
amat = {{A11, A12}, {A12, A22}}; abar = {{Ab11, Ab12}, {Ab12, Ab22}}; cmat = {{c11, c12}, {c21, c22}};
h = Sum[amat[[i, j]] cr[[i]] ** cr[[j]] + abar[[i, j]] an[[i]] ** an[[j]] + 2 cmat[[i, j]] cr[[i]] ** an[[j]], {i, 2}, {j, 2}]/2 + {f1, f2} . cr + {fb1, fb2} . an;
{Expand[BosonicNormalOrder[h ** # - # ** h, v] & /@ an + cmat . an + amat . cr + {f1, f2}],
 Expand[BosonicNormalOrder[h ** # - # ** h, v] & /@ cr - abar . an - Transpose[cmat] . cr - {fb1, fb2}]}
