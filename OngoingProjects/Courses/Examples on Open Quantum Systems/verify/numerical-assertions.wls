(* Seeded numerical/convergence assertions for the catalog's stochastic claims. *)

ClearAll[numericalCheck];
numericalFailures = {};
numericalCheck[name_String, held_HoldComplete] := Module[{value, messages},
  Block[{$MessageList = {}},
    value = Check[ReleaseHold[held], $Failed];
    messages = $MessageList
  ];
  If[TrueQ[value] && messages === {},
    WriteString[$Output, "PASS: ", name, "\n"],
    AppendTo[numericalFailures, <|"Name" -> name, "Value" -> ToString[value, InputForm],
      "Messages" -> ToString[messages, InputForm]|>];
    WriteString[$Output, "FAIL: ", name, " => ", ToString[value, InputForm],
      "; messages = ", ToString[messages, InputForm], "\n"]
  ]
];

biasRows = Table[{h, oneStepBias[h]}, {h, {0.1, 0.05, 0.025}}];
WriteString[$Output, "EVIDENCE: one-step biases = ", ToString[biasRows, InputForm], "\n"];
numericalCheck["finite-step ensemble bias contracts under refinement", HoldComplete[
  Greater @@ biasRows[[All, 2]] && Last[biasRows][[2]] < 10^-3
]];

rapidRows = Table[With[{n = Round[0.6/h]},
    With[{finals = Last[steeredBlur[#, n, h]] & /@ Range[80], target = 0.5 Exp[-4 kPur 0.6]},
      {h, Mean[finals], StandardDeviation[finals], Abs[Mean[finals] - target]}]],
  {h, {0.02, 0.01, 0.005}}];
WriteString[$Output, "EVIDENCE: rapid-purification refinement = ", ToString[rapidRows, InputForm], "\n"];
numericalCheck["rapid-purification path spread contracts under refinement", HoldComplete[
  Greater @@ rapidRows[[All, 3]]
]];
numericalCheck["rapid-purification mean bias contracts from coarsest to finest step", HoldComplete[
  Last[rapidRows][[4]] < First[rapidRows][[4]]
]];

fixedForPassage = fixedBlur[#, 400] & /@ Range[200];
steeredForPassage = steeredBlur[#, 400] & /@ Range[200];
fixedPassage = fptPur /@ fixedForPassage;
steeredPassage = fptPur /@ steeredForPassage;
censorCounts = {Count[fixedPassage, _Missing], Count[steeredPassage, _Missing]};
WriteString[$Output, "EVIDENCE: first-passage censor counts = ", ToString[censorCounts, InputForm], "\n"];
numericalCheck["seeded first-passage comparison is uncensored in the four-unit window", HoldComplete[censorCounts === {0, 0}]];
numericalCheck["fixed-axis strategy has the shorter seeded mean first-passage time", HoldComplete[
  Mean[fixedPassage] < Mean[steeredPassage]
]];

WriteString[$Output, "OPEN ISSUES: ", Length[numericalFailures], "\n"];
If[numericalFailures =!= {}, Exit[1], Exit[0]];
