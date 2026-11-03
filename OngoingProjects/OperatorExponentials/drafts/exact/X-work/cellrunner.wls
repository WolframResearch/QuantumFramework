Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
(WriteString["stdout", "=== ", FileBaseName[#], "\n", ToString[InputForm[TimeConstrained[Get[#], 120]]], "\n"]; ClearAll["Global`*"]) & /@ $ScriptCommandLine[[2 ;;]];
