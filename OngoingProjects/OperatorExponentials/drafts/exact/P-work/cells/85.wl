{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity] /. Conjugate[s_Plus] :> Map[Conjugate, s];
ph[x_] := Simplify[BosonicExpSimplify[d[-b] ** d[x]]/d[x - b]];
inner = BosonicMatrixElement[{al - b, alp - b}, Exp[Log[mu/(1 + mu)] ad ** a], "Basis" -> "Coherent"]/(1 + mu);
kern = Exp[(al Conjugate[al] + alp Conjugate[alp])/2] Conjugate[ph[al]] ph[alp] inner;
target = Exp[(mu Conjugate[al] alp + b Conjugate[al] + Conjugate[b] alp - b Conjugate[b])/(1 + mu)]/(1 + mu);
re = {al -> x1 + I y1, alp -> x2 + I y2, b -> u + I v};
Simplify[ComplexExpand[kern/target /. re], mu > 0]
