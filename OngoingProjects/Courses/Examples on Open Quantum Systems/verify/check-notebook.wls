source = FileNameJoin[{DirectoryName[$InputFileName], "..", "open-system-simulation-catalog.md"}];
notebook = FileNameJoin[{DirectoryName[$InputFileName], "..", "open-system-simulation-catalog.nb"}];

markdown = Import[source, "Text"];
expectedInputs = Length @ StringCases[
   markdown,
   RegularExpression["(?s)```wl\\s*(.*?)```"]
   ];

Block[{$MessageList = {}},
  heldNotebook = Check[Get[notebook], $Failed];
  notebookMessages = $MessageList
  ];
If[heldNotebook === $Failed,
  Print["FAIL notebook import"];
  Print["OPEN ISSUES: 1"];
  Exit[1]
  ];

inputCells = Cases[heldNotebook, Cell[_, "Input", ___], Infinity];
titleCells = Cases[heldNotebook, Cell[_, "Title", ___], Infinity];
sectionCells = Cases[
   heldNotebook,
   Cell[_, style : ("Chapter" | "Section" | "Subsection" | "Subsubsection"), ___] :> style,
   Infinity
   ];

checks = {
   "notebook parses" -> (Head[heldNotebook] =!= Failure),
   "notebook import emits no messages" -> (notebookMessages === {}),
   "input-cell parity" -> (Length[inputCells] === expectedInputs),
   "single title" -> (Length[titleCells] === 1),
   "section hierarchy present" -> (Length[sectionCells] > 20),
   "no conversion failures" -> FreeQ[heldNotebook, $Failed | Failure[__], Infinity]
   };

Scan[
  Function[check,
   Print[If[TrueQ[Last[check]], "PASS ", "FAIL "] <> First[check]]
   ],
  checks
  ];
Print["EXPECTED INPUT CELLS: " <> ToString[expectedInputs]];
Print["ACTUAL INPUT CELLS: " <> ToString[Length[inputCells]]];
Print["SECTION CELLS: " <> ToString[Length[sectionCells]]];
If[notebookMessages =!= {}, Print["IMPORT MESSAGES: " <> ToString[notebookMessages, InputForm]]];
issues = Count[checks, _?(Not @* TrueQ @* Last)];
Print["OPEN ISSUES: " <> ToString[issues]];
Exit[If[issues === 0, 0, 1]];
