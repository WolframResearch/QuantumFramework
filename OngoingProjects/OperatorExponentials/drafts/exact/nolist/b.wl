{a, ad} = FieldVariables[{1}];
res = <|"usage" -> BosonicNormalOrder::usage, "list" -> BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}], "single" -> BosonicNormalOrder[a ** ad, {a, ad}], "listNoVars" -> BosonicNormalOrder[{a ** ad, ad ** a}]|>;
WriteString["stdout", ToString[InputForm[res]], "\n"]
