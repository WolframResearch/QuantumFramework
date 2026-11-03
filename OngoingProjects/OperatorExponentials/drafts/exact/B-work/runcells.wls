Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
R`dir = FileNameJoin[{DirectoryName[$InputFileName], "cells"}];
R`sel = If[Length[$ScriptCommandLine] < 2, All, StringSplit[$ScriptCommandLine[[2]], ","]];
R`files = Select[FileNames["*.wl", R`dir], R`sel === All || MemberQ[R`sel, FileBaseName[#]] &];
R`out = OpenWrite[R`dir <> "/../cells-out-" <> If[R`sel === All, "all", StringRiffle[R`sel, "_"]] <> ".txt"];
Scan[
  Function[R`file,
    ClearAll["Global`*"];
    {R`t, {R`r, R`msgs}} = AbsoluteTiming[Block[{$MessageList = {}}, {TimeConstrained[Get[R`file], 300, $TimedOut], $MessageList}]];
    WriteString[R`out, "=== ", FileBaseName[R`file], " (", ToString[Round[R`t, 0.1]], " s) messages: ", ToString[R`msgs], "\n", ToString[InputForm[R`r]], "\n\n"]
  ],
  R`files];
Close[R`out];
