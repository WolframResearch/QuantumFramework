{a, ad} = FieldVariables[{1}];
qn = Refine[BosonicMatrixElement[{n, 0}, Exp[f ad]], n >= 0];
w = qn (qn /. f -> fb);
{DisplacementOperator[f, Infinity, "Ordering" -> "Normal"], qn, Sum[w, {n, 0, Infinity}], Exp[-f fb] w}
