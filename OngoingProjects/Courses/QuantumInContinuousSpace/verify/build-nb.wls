#!/usr/bin/env wolframscript
(* Build Solutions-Parts-1-2.nb from its literate-markdown twin via md2nb. *)
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/MarkdownToNotebook.wl"];
res = MarkdownToNotebook["Solutions-Parts-1-2.md", "Solutions-Parts-1-2.nb", "Evaluate" -> False];
WriteString[$Output, "RESULT: ", ToString[res], "\n"];
WriteString[$Output, "exists: ", ToString[FileExistsQ["Solutions-Parts-1-2.nb"]],
   "  bytes: ", ToString[FileByteCount["Solutions-Parts-1-2.nb"]], "\n"];
