(* Per-part verify + md2nb build. Usage: wolframscript -file _check-one.wls <Part-NN.md>
   Runs every wl cell in THIS fresh kernel and flags any that error or emit a message.
   A clean run then rebuilds the sibling .nb via md2nb (Evaluate -> False) and reports cell
   counts; a flagged run leaves the existing notebook untouched. Exits with status 1 when a
   cell was flagged, when no notebook was written, or when the notebook's Input-cell count
   differs from the source's wl-cell count. *)
part = $ScriptCommandLine[[2]];
Get["/Users/mohammadb/Documents/GitHub/MarkdownToNotebook/MarkdownToNotebook.wl"];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
txt = Import[part, "Text"];
cells = StringCases[txt, Shortest["```wl\n" ~~ c___ ~~ "\n```"] :> c];
bad = Flatten @ Position[Check[ToExpression[#]; True, False] & /@ cells, False];
out = FileNameJoin[{DirectoryName[part], FileBaseName[part] <> ".nb"}];
If[bad =!= {},
  WriteString[$Output,
    FileNameTake[part] <> ": wl-cells=" <> ToString[Length[cells]] <>
    " errored/messaged=" <> ToString[bad] <> "  (notebook not rebuilt)\n"];
  Exit[1]];
MarkdownToNotebook[part, out, "Evaluate" -> False];
nb = Get[out];
nin = Length[Cases[nb, Cell[_, "Input", ___], Infinity]];
nsub = Length[Cases[nb, Cell[_, "Subsection", ___], Infinity]];
WriteString[$Output,
  FileNameTake[part] <> ": wl-cells=" <> ToString[Length[cells]] <>
  " errored/messaged={}" <>
  "  ->  " <> FileNameTake[out] <> " (Input=" <> ToString[nin] <> ", Subsection=" <> ToString[nsub] <> ")\n"];
If[! FileExistsQ[out] || nin =!= Length[cells], Exit[1]];
