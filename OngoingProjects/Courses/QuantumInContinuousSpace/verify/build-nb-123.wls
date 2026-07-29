#!/usr/bin/env wolframscript
(* Build Solutions-Parts-1-3.nb from the merged literate-markdown twin via md2nb.
   Template: Default -> StyleDefinitions Default.nb; Title from the body H1.
   Evaluate -> False keeps Input cells only (reader-runs-cells design). *)
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/MarkdownToNotebook.wl"];
res = MarkdownToNotebook["Solutions-Parts-1-3.md", "Solutions-Parts-1-3.nb", "Evaluate" -> False];
WriteString[$Output, "RESULT: ", ToString[res], "\n"];
WriteString[$Output, "exists: ", ToString[FileExistsQ["Solutions-Parts-1-3.nb"]],
   "  bytes: ", ToString[FileByteCount["Solutions-Parts-1-3.nb"]], "\n"];
