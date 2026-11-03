{a, ad} = FieldVariables[{1}];
WriteString["stdout", ToString[InputForm[{BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}], BosonicNormalOrder[a ** ad, {a, ad}]}]], "\n"];
