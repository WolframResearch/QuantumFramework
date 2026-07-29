#!/usr/bin/env wolframscript
(* Verify the built Solutions-Parts-1-3.nb:
   (1) cell-style counts + which headings survived (md2nb silently drops titles
       outside its taxonomy, so confirm all 10 Subsections rendered);
   (2) dropped-LaTeX scan (empty/Null FormBox, raw backslash TeX);
   (3) hbar -> \[HBar] glyph present. *)
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
nb = Get["Solutions-Parts-1-3.nb"];

styles = Cases[nb, Cell[_, s_String, ___] :> s, Infinity];
WriteString[$Output, "=== cell-style counts ===\n"];
Scan[WriteString[$Output, "  ", #[[1]], " : ", ToString[#[[2]]], "\n"] &, SortBy[Tally[styles], First]];

(* which headings rendered? pull the plain text out of Title/Section/Subsection cells *)
headText[c_] := StringJoin[Cases[c, _String, Infinity]];
WriteString[$Output, "=== Title ===\n"];
Scan[WriteString[$Output, "  ", headText[#], "\n"] &, Cases[nb, Cell[_, "Title", ___], Infinity]];
WriteString[$Output, "=== Sections ===\n"];
Scan[WriteString[$Output, "  ", headText[#], "\n"] &, Cases[nb, Cell[_, "Section", ___], Infinity]];
WriteString[$Output, "=== Subsections (expect 1.1-1.9, 2.1-2.9, 3.1-3.10, 28 total) ===\n"];
subs = Cases[nb, Cell[_, "Subsection", ___], Infinity];
WriteString[$Output, "  count: ", ToString[Length[subs]], "\n"];
Scan[WriteString[$Output, "  ", StringTake[headText[#], UpTo[60]], "\n"] &, subs];

(* dropped LaTeX: empty / Null FormBox *)
emptyForms = Cases[nb,
   fb : FormBox[b_, ___] /; (MatchQ[b, "" | "\[Null]" | RowBox[{}] | RowBox[{""}]] ||
      (StringQ[b] && StringTrim[b] === "")) :> fb, Infinity];
WriteString[$Output, "=== math integrity ===\n"];
WriteString[$Output, "  InlineFormula cells: ", ToString[Length[Cases[nb, Cell[_, "InlineFormula", ___], Infinity]]], "\n"];
WriteString[$Output, "  FormBox total: ", ToString[Length[Cases[nb, _FormBox, Infinity]]], "\n"];
WriteString[$Output, "  EMPTY/Null FormBox (dropped math): ", ToString[Length[emptyForms]], "\n"];

(* raw un-parsed TeX: a backslash-command left in a FormBox string *)
strs = Cases[nb, FormBox[stuff_, ___] :> stuff, Infinity];
rawTeX = Select[Cases[strs, _String, Infinity],
   StringContainsQ[#, "\\"] && ! StringMatchQ[#, WhitespaceCharacter ...] &];
WriteString[$Output, "  FormBox strings with a stray backslash: ", ToString[Length[rawTeX]], "\n"];
If[Length[rawTeX] > 0, Scan[WriteString[$Output, "    RAW: ", #, "\n"] &, Take[rawTeX, UpTo[12]]]];

(* hbar glyph *)
nHBar = Length[Select[Cases[nb, _String, Infinity], StringContainsQ[#, "\[HBar]"] &]];
WriteString[$Output, "  HBar glyph occurrences: ", ToString[nHBar], "\n"];
WriteString[$Output, "  Input cells: ", ToString[Length[Cases[nb, Cell[_, "Input", ___], Infinity]]], "\n"];
WriteString[$Output, "DONE\n"];
