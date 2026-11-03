emit[name_, val_] := WriteString["stdout", name, " => ", ToString[InputForm[val]], "\n"];
SetAttributes[chk, HoldRest];
chk[name_, expr_] := emit[name, TimeConstrained[expr, 60, "TIMEOUT"]];
