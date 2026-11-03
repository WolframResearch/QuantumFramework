{a, ad} = FieldVariables[{1}];
hd = eps ad + Conjugate[eps] a;
BosonicNormalOrder[-I (a ** hd - hd ** a), {a, ad}]
