(* Small-fragment reference bridge, not a competition-scale simulator.
   Kickmix q0 maps to QF wire 1; registers remain little-endian.
   HMR branches retain their 1/Sqrt[2] Kraus normalization and classical fixups.
   R is a clean-wire assertion (projector), not an arbitrary physical reset.
   Unsupported instructions fail explicitly. No source text is evaluated as WL. *)
BeginPackage["QFKickmixPrototype`", {"Wolfram`QuantumFramework`"}];
KickmixBranchOperator::usage = "KickmixBranchOperator[fragment, outcomes] returns a small QF branch operator for specified HMR outcomes.";
Begin["`Private`"];

wire[s_String] /; StringMatchQ[s, RegularExpression["q[0-9]+"]] := FromDigits[StringDrop[s, 1]] + 1;
bitKey[s_String] := "b" <> IntegerString[FromDigits[StringDrop[s, 1]]];
asciiLineQ[s_String] := AllTrue[ToCharacterCode[s], # === 9 || 32 <= # <= 126 &];

argumentKinds = <|"X" -> {"q"}, "Z" -> {"q"}, "CX" -> {"q", "q"}, "CZ" -> {"q", "q"},
    "CCX" -> {"q", "q", "q"}, "CCZ" -> {"q", "q", "q"}, "SWAP" -> {"q", "q"},
    "Z_POW" -> {"q", "angle"}, "R" -> {"q"}, "HMR" -> {"q", "b"}|>;
tokenQ["q", s_String] := StringMatchQ[s, RegularExpression["q[0-9]+"]];
tokenQ["b", s_String] := StringMatchQ[s, RegularExpression["b[0-9]+"]];
tokenQ["angle", s_String] := StringMatchQ[s, RegularExpression["[0-9]+\\.[0-9]+"]];
instructionQ[tokens_List] := With[{kinds = Lookup[argumentKinds, First[tokens], Missing["Instruction"]]},
    ListQ[kinds] && Length[kinds] === Length[tokens] - 1 &&
    And @@ MapThread[tokenQ, {kinds, Rest[tokens]}] &&
    DuplicateFreeQ[wire /@ Pick[Rest[tokens], kinds, "q"]]];

decimal[s_String] /; tokenQ["angle", s] :=
    With[{parts = StringSplit[s, "."]}, FromDigits[First[parts]] + FromDigits[Last[parts]]/10^StringLength[Last[parts]]];

gate[{"X", q_}] := QuantumOperator["X", {wire[q]}];
gate[{"Z", q_}] := QuantumOperator["Z", {wire[q]}];
gate[{"CX", a_, t_}] := QuantumOperator["CNOT", wire /@ {a, t}];
gate[{"CZ", a_, t_}] := QuantumOperator["CZ", wire /@ {a, t}];
gate[{"CCX", a_, b_, t_}] := QuantumOperator["C"["X", wire /@ {a, b}, {}], {wire[t]}];
gate[{"CCZ", a_, b_, t_}] := QuantumOperator["C"["Z", wire /@ {a, b}, {}], {wire[t]}];
gate[{"SWAP", a_, b_}] := QuantumOperator["SWAP", wire /@ {a, b}];
gate[{"Z_POW", q_, exponent_}] /; tokenQ["angle", exponent] :=
    QuantumOperator["Phase"[Pi decimal[exponent]], {wire[q]}];
gate[{"R", q_}] := QuantumOperator[{{1, 0}, {0, 0}}, {wire[q]}];
gate[tokens_] := Failure["UnsupportedInstruction", <|"Tokens" -> tokens|>];

step[{bits_, outcomes_}, line_String] := Enclose @ Module[
    {tokens = StringSplit[line], condition, active, instruction, branch, nextBits = bits, remaining = outcomes, op},
    ConfirmAssert[tokens =!= {}];
    condition = MemberQ[tokens, "if"];
    If[condition,
        ConfirmAssert[Count[tokens, "if"] === 1 && Length[tokens] >= 3 && tokens[[-2]] === "if" &&
            tokenQ["b", Last[tokens]]]
    ];
    active = !condition || Lookup[bits, bitKey[Last[tokens]], 0] === 1;
    If[condition,
        tokens = Drop[tokens, -2]
    ];
    ConfirmAssert[instructionQ[tokens]];
    If[!active, Return[{bits, outcomes}]];
    instruction = First[tokens];
    If[instruction === "HMR",
        ConfirmAssert[Length[tokens] === 3 && Length[outcomes] > 0];
        branch = First[outcomes];
        remaining = Rest[outcomes];
        nextBits = Append[bits, bitKey[tokens[[3]]] -> branch];
        op = QuantumOperator[{{1, (-1)^branch}, {0, 0}}/Sqrt[2], {wire[tokens[[2]]]}],
        op = gate[tokens]
    ];
    Sow[Confirm[op], "Operators"];
    {nextBits, remaining}
];

KickmixBranchOperator[fragment_Association, outcomes_List : {}] := Enclose @ Module[
    {n = Lookup[fragment, "NumQubits"], lines = Lookup[fragment, "Lines"], folded, collected, operators, compiled},
    ConfirmAssert[IntegerQ[n] && 1 <= n <= 8];
    ConfirmAssert[VectorQ[lines, StringQ] && VectorQ[outcomes, MemberQ[{0, 1}, #] &]];
    ConfirmAssert[AllTrue[lines, asciiLineQ]];
    ConfirmAssert[AllTrue[Cases[Catenate[StringSplit /@ lines], s_String /; tokenQ["q", s]],
        1 <= wire[#] <= n &]];
    {folded, collected} = Reap[Fold[Confirm @ step[#1, #2] &, {<||>, outcomes}, lines], "Operators"];
    ConfirmAssert[folded[[2]] === {}];
    operators = If[collected === {}, {}, First[collected]];
    ConfirmAssert[AllTrue[operators, MatchQ[#, _QuantumOperator] &]];
    compiled = QuantumCircuitOperator[Join[QuantumOperator["I", {#}] & /@ Range[n], operators]][
        "QuantumOperator", Method -> "Schrodinger"];
    ConfirmAssert[compiled["InputDimension"] === 2^n && compiled["OutputDimension"] === 2^n];
    compiled
];

KickmixBranchOperator[___] := Failure["ArgumentShape", <|
    "MessageTemplate" -> "Expected an Association and a list of binary HMR outcomes."|>];

End[];
EndPackage[];
