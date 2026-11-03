(* Runs every cell file exactly as written, each after clearing Global`, and writes
   "NN :: <InputForm of the last expression>" per cell. Optional 3rd script arg: cell list. *)
drv`cellDir = FileNameJoin[{DirectoryName[$InputFileName], "cells"}];
drv`files = Sort[FileNames["*.wl", drv`cellDir]];
drv`sel = If[Length[$ScriptCommandLine] >= 3, StringSplit[$ScriptCommandLine[[3]], ","], All];
drv`files = If[drv`sel === All, drv`files, Select[drv`files, MemberQ[drv`sel, FileBaseName[#]] &]];
drv`run[drv`f_] := Module[{drv`t0 = AbsoluteTime[], drv`v},
  ClearAll["Global`*"];
  drv`v = TimeConstrained[Get[drv`f], 900, "TIMEOUT"];
  FileBaseName[drv`f] <> " [" <> ToString[Round[AbsoluteTime[] - drv`t0, 0.1]] <> " s] :: " <> ToString[InputForm[drv`v]]];
Scan[WriteString["stdout", drv`run[#], "\n"] &, drv`files];
