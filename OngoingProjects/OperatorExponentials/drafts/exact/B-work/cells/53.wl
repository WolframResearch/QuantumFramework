{a, ad} = FieldVariables[{1}];
h = (t ad ** ad + t a ** a + 2 s ad ** a)/2;
om = Sqrt[s^2 - t^2]; ch = Sqrt[(s/om + 1)/2]; sh = Sqrt[(s/om - 1)/2];
{b, bd} = {ch a + sh ad, ch ad + sh a};
Simplify[{BosonicNormalOrder[b ** bd - bd ** b, {a, ad}], BosonicNormalOrder[om bd ** b + (om - s)/2 - h, {a, ad}]}, s > t > 0]
