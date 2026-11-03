drv`dir = FileNameJoin[{DirectoryName[$InputFileName], $ScriptCommandLine[[3]]}];
drv`files = Sort[FileNames["*.wl", drv`dir]];
drv`sel = If[Length[$ScriptCommandLine] >= 4, ToExpression[$ScriptCommandLine[[4]]], All];
drv`files = If[drv`sel === All, drv`files, Select[drv`files, MemberQ[drv`sel, ToExpression[FileBaseName[#]]] &]];
Scan[Function[drv`f,
   {drv`t, {drv`res, drv`msgs}} = AbsoluteTiming[Block[{$MessageList = {}}, {TimeConstrained[Get[drv`f], 300, "TIMEOUT"], $MessageList}]];
   WriteString["stdout", "#### ", FileBaseName[drv`f], " (", ToString[Round[drv`t, 0.1]], " s) msgs=", ToString[drv`msgs, InputForm], "\n", ToString[drv`res, InputForm], "\n"];
   ClearAll["Global`*"]], drv`files];
