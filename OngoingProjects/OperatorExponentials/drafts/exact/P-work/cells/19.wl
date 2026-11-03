{a, ad} = FieldVariables[{1}];
q = Exp[-al Conjugate[al]] (BosonicMatrixElement[{0, 0}, Exp[cal a] ** DisplacementOperator[g, Infinity] ** Exp[al ad]] /. cal -> Conjugate[al]);
Simplify[q/Exp[-g Conjugate[g]/2 + g Conjugate[al] - Conjugate[g] al]]
