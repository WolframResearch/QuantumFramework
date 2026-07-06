(* Per-part verify + md2nb build. Usage: wolframscript -file _check-one.wls <Part-0N.md>
   Runs every wl cell in THIS fresh kernel (flags any that error or emit a message),
   then rebuilds the sibling .nb via md2nb (Evaluate -> False) and reports cell counts. *)
part = $ScriptCommandLine[[2]];
Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/MarkdownToNotebook.wl"];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
txt = Import[part, "Text"];
cells = StringCases[txt, Shortest["```wl\n" ~~ c___ ~~ "\n```"] :> c];
bad = {};
Do[If[Check[ToExpression[cells[[i]]]; True, False] =!= True, AppendTo[bad, i]], {i, Length[cells]}];
out = StringReplace[part, ".md" -> ".nb"];
MarkdownToNotebook[part, out, "Evaluate" -> False];
nb = Get[out];
nin = Length[Cases[nb, Cell[_, "Input", ___], Infinity]];
nsub = Length[Cases[nb, Cell[_, "Subsection", ___], Infinity]];
WriteString[$Output,
  FileNameTake[part] <> ": wl-cells=" <> ToString[Length[cells]] <>
  " errored/messaged=" <> ToString[bad] <>
  "  ->  " <> FileNameTake[out] <> " (Input=" <> ToString[nin] <> ", Subsection=" <> ToString[nsub] <> ")\n"];
