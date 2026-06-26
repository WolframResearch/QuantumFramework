#!/usr/bin/env wolframscript
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
nb = Get["Solutions-Parts-1-2.nb"];
WriteString[$Output, "styles: ", ToString[Tally[Cases[nb, Cell[_, s_String, ___] :> s, Infinity]], InputForm], "\n"];
inl = Cases[nb, Cell[BoxData[b_], "InlineFormula", ___] :> b, Infinity];
empty = Select[inl, (MatchQ[#, "" | "\[Null]" | RowBox[{}] | RowBox[{""}] | RowBox[{"\[Null]"}]] || (StringQ[#] && StringTrim[#] === "")) &];
WriteString[$Output, "InlineFormula: ", ToString[Length[inl]], "  empty: ", ToString[Length[empty]], "  Null-leaves: ", ToString[Count[inl, "\[Null]", Infinity]], "\n"];
WriteString[$Output, "has sech text: ", ToString[Length[Select[Cases[inl, _String, Infinity], StringContainsQ[#, "sech"] &]] > 0], "\n"];
WriteString[$Output, "HBar glyphs: ", ToString[Length[Select[Cases[nb, _String, Infinity], StringContainsQ[#, "\[HBar]"] &]]], "\n"];
WriteString[$Output, "Input cells: ", ToString[Length[Cases[nb, Cell[_, "Input", ___], Infinity]]], "\n"];
WriteString[$Output, "DONE\n"];
