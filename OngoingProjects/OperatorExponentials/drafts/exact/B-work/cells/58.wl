{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
cr = {a1d, a2d}; an = {a1, a2}; u = {{u11, u12}, {u12, u22}};
vq = Sum[u[[i, j]] cr[[i]] ** cr[[j]] ** an[[j]] ** an[[i]], {i, 2}, {j, 2}];
{BosonicMatrixElement[{{1, 1}, {1, 1}}, vq], BosonicMatrixElement[{{2, 0}, {2, 0}}, vq], BosonicMatrixElement[{{1, 1}, {2, 0}}, vq]}
