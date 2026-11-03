Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
dir = FileNameJoin[{DirectoryName[$InputFileName], "cells"}];
files = Sort[FileNames["*.wl", dir]];
sel = Select[$ScriptCommandLine, StringMatchQ[#, "C" ~~ DigitCharacter .. ] &];
files = If[sel === {}, files, Select[files, MemberQ[sel, FileBaseName[#]] &]];
Scan[Function[file,
   ClearAll["Global`*"];
   CellResults`r = TimeConstrained[Get[file], 300, $TimedOut];
   WriteString["stdout", "=== ", FileBaseName[file], "\n", ToString[InputForm[CellResults`r]], "\n"]], files];
