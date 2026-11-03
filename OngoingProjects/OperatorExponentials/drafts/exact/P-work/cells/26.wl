{a, ad} = FieldVariables[{1}];
dd = DisplacementOperator[g, Infinity];
Simplify[{BosonicExpSimplify[Exp[g ad] ** Exp[-Conjugate[g] a]]/dd, BosonicExpSimplify[Exp[-Conjugate[g] a] ** Exp[g ad]]/dd}]
