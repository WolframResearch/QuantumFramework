{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
c = (a + b)/Sqrt[2]; cd = (ad + bd)/Sqrt[2]; d = (a - b)/Sqrt[2]; dd = (ad - bd)/Sqrt[2];
xi = r E^(I th);
gc = (Conjugate[xi] c ** c - xi cd ** cd)/2; gd = (-Conjugate[xi] d ** d + xi dd ** dd)/2;
{BosonicNormalOrder[gc + gd - (Conjugate[xi] a ** b - xi ad ** bd), v], BosonicNormalOrder[gc ** gd - gd ** gc, v], BosonicNormalOrder[c ** dd - dd ** c, v], BosonicNormalOrder[c ** cd - cd ** c, v]}
