{a, ad} = FieldVariables[{1}];
k = BosonicMatrixElement[{0, 0}, Exp[cal a] ** DisplacementOperator[g, Infinity] ** Exp[be ad]];
Simplify[k/Exp[-g Conjugate[g]/2 + cal be + cal g - be Conjugate[g]]]
