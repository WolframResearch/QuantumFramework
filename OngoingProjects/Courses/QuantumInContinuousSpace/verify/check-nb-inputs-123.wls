#!/usr/bin/env wolframscript
(* Every Input cell in the built nb must reparse (no ToExpression::sntx, no the
   `_\[Name]` -> `_*\[Name]` blank-head mangling) AND run. Confirms md2nb did not
   corrupt the WL on the way into the notebook. *)
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
nb = Get["Solutions-Parts-1-3.nb"];
inputs = Cases[nb, Cell[BoxData[b_], "Input", ___] :> b, Infinity];
WriteString[$Output, "Input cells found: ", ToString[Length[inputs]], "\n"];
i = 0; nfail = 0;
Scan[(
   i++;
   exprs = Check[ToExpression[#, StandardForm, Hold], $FAILEDPARSE, {ToExpression::sntx}];
   If[exprs === $FAILEDPARSE,
     nfail++; WriteString[$Output, "cell ", ToString[i], " PARSE-FAIL\n"],
     Check[ReleaseHold[exprs], WriteString[$Output, "cell ", ToString[i], " RUN-MSG\n"]]
   ]
) &, inputs];
WriteString[$Output, "parse failures: ", ToString[nfail], "\n"];
WriteString[$Output, If[nfail == 0, "ALL-INPUT-CELLS-PARSED-AND-RAN\n", "SOME-CELLS-FAILED\n"]];
