(* Probes BosonicNormalOrder on list input. Loaded by a runner after load.wl. *)
{a, ad} = FieldVariables[{1}];
{b, bd} = FieldVariables[{2}];
results = <|
  "qfLocation" -> PacletObject["Wolfram/QuantumFramework"]["Location"],
  "qfVersion" -> PacletObject["Wolfram/QuantumFramework"]["Version"],
  "attributes" -> Attributes[BosonicNormalOrder],
  "usage" -> BosonicNormalOrder::usage,
  "downValueLHS" -> Map[First, DownValues[BosonicNormalOrder]],
  "single" -> BosonicNormalOrder[a ** ad, {a, ad}],
  "list1" -> Check[BosonicNormalOrder[{a ** ad}, {a, ad}], "MESSAGE"],
  "list2" -> Check[BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}], "MESSAGE"],
  "list2swap" -> Check[BosonicNormalOrder[{ad ** a, a ** ad}, {a, ad}], "MESSAGE"],
  "list3" -> Check[BosonicNormalOrder[{a ** ad, a ** a ** ad, ad ** a ** a ** ad}, {a, ad}], "MESSAGE"],
  "list2mixed" -> Check[BosonicNormalOrder[{a ** ad + 2 a, 3 + a ** a ** ad ** ad}, {a, ad}], "MESSAGE"],
  "list2multimode" -> Check[BosonicNormalOrder[{a ** ad ** b, b ** bd ** a}, {a, ad, b, bd}], "MESSAGE"],
  "list2blasiak" -> Check[BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}, Method -> "Blasiak"], "MESSAGE"],
  "list2noVars" -> Check[BosonicNormalOrder[{a ** ad, ad ** a}], "MESSAGE"],
  "singleNoVars" -> Check[BosonicNormalOrder[a ** ad], "MESSAGE"],
  "mapped" -> Map[BosonicNormalOrder[#, {a, ad}] &, {a ** ad, ad ** a}],
  "ncprListDirect" -> NonCommutativePolynomialReduce[{a ** ad, ad ** a},
     {a ** ad - ad ** a - 1}, NonCommutativeAlgebra[<|"Generators" -> {a, ad}|>]],
  "antinormalList2" -> Check[BosonicAntinormalOrder[{ad ** a, a ** ad}, {a, ad}], "MESSAGE"]
|>;
WriteString["stdout", ToString[InputForm[results]], "\n"]
