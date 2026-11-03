pp = Integrate[Exp[-nb (x^2 + y^2)] Exp[2 I (x q - y p)], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> nb > 0 && Element[{p, q}, Reals]]/Pi^2;
cn = Integrate[Exp[-(p^2 + q^2)/nb]/(Pi nb) Exp[2 I (y p - x q)], {p, -Infinity, Infinity}, {q, -Infinity, Infinity}, Assumptions -> nb > 0 && Element[{x, y}, Reals]];
{Simplify[pp - Exp[-(p^2 + q^2)/nb]/(Pi nb)], Simplify[cn - Exp[-nb (x^2 + y^2)]],
 ComplexExpand[((x - I y) (p + I q) - (x + I y) (p - I q)) - 2 I (x q - y p)], ComplexExpand[((x + I y) (p - I q) - (x - I y) (p + I q)) - 2 I (y p - x q)]}
