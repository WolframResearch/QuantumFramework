{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
Simplify[Exp[-I P Q/2] BosonicExpSimplify[Exp[I P q] ** Exp[-I Q p]]/DisplacementOperator[(Q + I P)/Sqrt[2], Infinity],
  Element[{P, Q}, Reals]]
