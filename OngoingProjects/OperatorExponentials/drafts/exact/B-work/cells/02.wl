{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
k = {{k11, k12}, {k12, k22}}; kb = {{kb11, kb12}, {kb12, kb22}};
ket = Sum[k[[i, j]] {a1d, a2d}[[i]] ** {a1d, a2d}[[j]], {i, 2}, {j, 2}]/Sqrt[2];
bra = Sum[kb[[i, j]] {a1, a2}[[j]] ** {a1, a2}[[i]], {i, 2}, {j, 2}]/Sqrt[2];
Expand[BosonicVEV[bra ** ket, {a1, a1d, a2, a2d}]]
