#!/usr/bin/env wolframscript
(* Verify the built Solutions-Parts-1-2.nb: cell-style counts + scan for dropped LaTeX. *)
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
nb = Get["Solutions-Parts-1-2.nb"];

styles = Cases[nb, Cell[_, s_String, ___] :> s, Infinity];
WriteString[$Output, "=== cell-style counts ===\n"];
Scan[WriteString[$Output, #[[1]], " : ", ToString[#[[2]]], "\n"] &, Tally[styles]];

(* Count InlineFormula (inline $...$) and any FormBox math *)
nInline = Length[Cases[nb, Cell[_, "InlineFormula", ___], Infinity]];
nForm = Length[Cases[nb, _FormBox, Infinity]];
WriteString[$Output, "InlineFormula cells: ", ToString[nInline], "\n"];
WriteString[$Output, "FormBox total: ", ToString[nForm], "\n"];

(* Dropped-LaTeX signature: a FormBox whose content is empty / just \[Null] / "" *)
emptyForms = Cases[nb,
   fb : FormBox[b_, ___] /; (
      MatchQ[b, "" | "\[Null]" | RowBox[{}] | RowBox[{""}]] ||
      (StringQ[b] && StringTrim[b] === "")
   ) :> fb, Infinity];
WriteString[$Output, "EMPTY/Null FormBox count: ", ToString[Length[emptyForms]], "\n"];

(* Scan all string leaves inside FormBoxes for a stray backslash-command (un-parsed LaTeX) *)
strs = Cases[nb, FormBox[stuff_, ___] :> stuff, Infinity];
rawTeX = Select[Cases[strs, _String, Infinity], StringContainsQ[#, "\\"] &];
WriteString[$Output, "FormBox strings containing a backslash (possible raw TeX): ", ToString[Length[rawTeX]], "\n"];
If[Length[rawTeX] > 0, Scan[WriteString[$Output, "  RAW: ", #, "\n"] &, Take[rawTeX, UpTo[10]]]];

(* HBar presence (the hbar -> \[HBar] check) *)
nHBar = Length[Cases[nb, "\[HBar]", Infinity]] + Length[Select[Cases[nb, _String, Infinity], StringContainsQ[#, "\[HBar]"] &]];
WriteString[$Output, "HBar glyph occurrences: ", ToString[nHBar], "\n"];

(* Input cell count (the wl code cells) *)
WriteString[$Output, "Input cells: ", ToString[Length[Cases[nb, Cell[_, "Input", ___], Infinity]]], "\n"];
WriteString[$Output, "DONE\n"];
