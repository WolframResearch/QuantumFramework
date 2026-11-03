sec[s_] := Integrate[Integrate[Exp[s I dl (t1 - t2)], {t2, 0, t1}], {t1, 0, t}];
nsec[s_] := Integrate[Integrate[Exp[s I dl (t1 + t2)], {t2, 0, t1}], {t1, 0, t}];
{Simplify[sec[1] - (I t/dl + (1 - E^(I dl t))/dl^2)], Simplify[sec[-1] - (-I t/dl + (1 - E^(-I dl t))/dl^2)], Simplify[nsec[1]], Simplify[nsec[-1]],
 Simplify[-g^2 (aa1 I t/dl + aa2 (-I t/dl)) - (-I g^2 t/dl) (aa1 - aa2)]}
