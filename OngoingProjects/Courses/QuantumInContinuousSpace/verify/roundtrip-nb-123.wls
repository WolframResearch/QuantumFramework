#!/usr/bin/env wolframscript
(* Round-trip the built nb back to markdown via the canonical NotebookToMarkdown,
   then confirm the twin md recovers the structure and the code verbatim. *)
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/NotebookToMarkdown.wl"];
md = NotebookToMarkdown["Solutions-Parts-1-3.nb"];
Export["verify/Solutions-Parts-1-3.roundtrip.md", md, "Text"];
WriteString[$Output, "roundtrip md bytes: ", ToString[StringLength[md]], "\n"];

(* structure recovered? *)
lines = StringSplit[md, "\n"];
h1 = Select[lines, StringMatchQ[#, "# " ~~ ___] &];
h2 = Select[lines, StringMatchQ[#, "## " ~~ ___] &];
h3 = Select[lines, StringMatchQ[#, "### " ~~ ___] &];
WriteString[$Output, "H1: ", ToString[Length[h1]], "  H2: ", ToString[Length[h2]],
   "  H3: ", ToString[Length[h3]], "\n"];
Scan[WriteString[$Output, "  ", StringTake[#, UpTo[55]], "\n"] &, h3];

(* code fences recovered? *)
fences = Length[StringCases[md, "```wl"]];
WriteString[$Output, "wl code fences: ", ToString[fences], "\n"];

(* spot-check that a couple of signature code lines survived verbatim *)
probes = {"Vpt[x_, \[Lambda]_] :=", "LegendreP[\[Lambda], \[Lambda] - n, Tanh[x]]",
   "AiryAiZero[n]", "Partition[Transpose[{ks, 1 + h^2 ks/12}], 3, 1]", "tAmp[\[Lambda]_, k_]"};
Scan[WriteString[$Output, "  present? ", ToString[StringContainsQ[md, #]], "  <- ", #, "\n"] &, probes];
WriteString[$Output, "DONE\n"];
