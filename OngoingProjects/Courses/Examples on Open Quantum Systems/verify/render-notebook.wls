Needs["JLink`"];

source = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/Examples on Open Quantum Systems/open-system-simulation-catalog.nb";
target = "/tmp/open-system-simulation-catalog-preview.pdf";

Block[{$MessageList = {}},
  result = Check[
    UsingFrontEnd[
      Module[{notebookObject, exportResult},
        notebookObject = NotebookOpen[source, Visible -> False];
        If[Head[notebookObject] =!= NotebookObject, Return[$Failed, Module]];
        exportResult = Export[target, notebookObject, "PDF"];
        NotebookClose[notebookObject];
        exportResult
        ]
      ],
    $Failed
    ];
  renderMessages = $MessageList
  ];

Print["RENDER RESULT: " <> ToString[result, InputForm]];
If[renderMessages =!= {},
  Print["RENDER MESSAGES: " <> ToString[renderMessages, InputForm]]
  ];
If[result === $Failed || renderMessages =!= {} || ! FileExistsQ[target], Exit[1], Exit[0]];
