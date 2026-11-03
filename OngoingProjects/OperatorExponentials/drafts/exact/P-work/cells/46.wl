{a, ad} = FieldVariables[{1}];
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
gen2 = -(z11 ad1 ** ad1 + 2 z12 ad1 ** ad2 + z22 ad2 ** ad2)/2;
{BosonicExpSimplify[a ** Exp[-z/2 ad ** ad]], BosonicExpSimplify[a1 ** Exp[gen2]], BosonicExpSimplify[a2 ** Exp[gen2]]}
