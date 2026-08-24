(* Rebuild the generated notebook twin from the audited Markdown source. *)
Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/MarkdownToNotebook.wl"];

(* Stage counters make a long conversion auditable without changing its result. *)
progress[label_String] := (WriteString[$Output, label, "\n"]; Flush[$Output]);

DownValues[originalLitParse] = DownValues[litParse] /. HoldPattern[litParse] :> originalLitParse;
Clear[litParse];
litParse[args___] := Module[{result},
  progress["STAGE litParse begin"];
  result = originalLitParse[args];
  progress["STAGE litParse end"];
  result
  ];

$inlineTextCount = 0;
DownValues[originalInlineTextData] = DownValues[inlineTextData] /. HoldPattern[inlineTextData] :> originalInlineTextData;
Clear[inlineTextData];
inlineTextData[args___] := Module[{result},
  $inlineTextCount++;
  If[Mod[$inlineTextCount, 100] === 0,
   progress["PROGRESS inline text " <> ToString[$inlineTextCount]]
   ];
  result = originalInlineTextData[args];
  result
  ];

$texBoxCount = 0;
DownValues[originalTexBoxes] = DownValues[texBoxes] /. HoldPattern[texBoxes] :> originalTexBoxes;
Clear[texBoxes];
texBoxes[args___] := Module[{result},
  $texBoxCount++;
  If[Mod[$texBoxCount, 100] === 0,
   progress["PROGRESS TeX expressions " <> ToString[$texBoxCount]]
   ];
  result = originalTexBoxes[args];
  result
  ];

DownValues[originalBuildNotebook] = DownValues[buildNotebook] /. HoldPattern[buildNotebook] :> originalBuildNotebook;
Clear[buildNotebook];
buildNotebook[args___] := Module[{result},
  progress["STAGE buildNotebook begin"];
  result = originalBuildNotebook[args];
  progress["STAGE buildNotebook end"];
  result
  ];

DownValues[originalWithDeterministicIDs] = DownValues[withDeterministicIDs] /. HoldPattern[withDeterministicIDs] :> originalWithDeterministicIDs;
Clear[withDeterministicIDs];
withDeterministicIDs[args___] := Module[{result},
  progress["STAGE deterministic IDs begin"];
  result = originalWithDeterministicIDs[args];
  progress["STAGE deterministic IDs end"];
  result
  ];

DownValues[originalExportDeterministicNB] = DownValues[exportDeterministicNB] /. HoldPattern[exportDeterministicNB] :> originalExportDeterministicNB;
Clear[exportDeterministicNB];
exportDeterministicNB[args___] := Module[{result},
  progress["STAGE export begin"];
  result = originalExportDeterministicNB[args];
  progress["STAGE export end"];
  result
  ];

source = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/Examples on Open Quantum Systems/open-system-simulation-catalog.md";
target = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/Examples on Open Quantum Systems/open-system-simulation-catalog.nb";
Block[{$MessageList = {}},
  result = Check[MarkdownToNotebook[source, target, "Evaluate" -> False], $Failed];
  buildMessages = $MessageList
  ];
WriteString[$Output, "BUILD RESULT: ", ToString[result, InputForm], "\n"];
If[buildMessages =!= {},
  WriteString[$Output, "BUILD MESSAGES: ", ToString[buildMessages, InputForm], "\n"]
  ];
If[result === $Failed || FailureQ[result] || buildMessages =!= {} || !FileExistsQ[target],
  Exit[1],
  Exit[0]
  ];
