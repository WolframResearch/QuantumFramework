{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
{BosonicExpSimplify[Exp[I (p0 q - q0 p)] ** q ** Exp[-I (p0 q - q0 p)]],
 BosonicExpSimplify[Exp[I (p0 q - q0 p)] ** p ** Exp[-I (p0 q - q0 p)]],
 BosonicExpSimplify[Exp[-I q0 p0/2] Exp[I p0 q] ** Exp[-I q0 p]]}
