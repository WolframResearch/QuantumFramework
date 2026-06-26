#!/usr/bin/env wolframscript
SetDirectory["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/QuantumInContinuousSpace"];
Off[Integrate::idiv];
nb = Get["Solutions-Parts-1-2.nb"];
inputs = Cases[nb, Cell[BoxData[b_], "Input", ___] :> b, Infinity];
WriteString[$Output, "Input cells found: ", ToString[Length[inputs]], "\n"];
i = 0;
Scan[(
   i++;
   exprs = Quiet[Check[ToExpression[#, StandardForm, Hold], $FAILEDPARSE], {ToExpression::sntx}];
   If[exprs === $FAILEDPARSE,
     WriteString[$Output, "cell ", ToString[i], " PARSE-FAIL\n"],
     ReleaseHold[exprs]
   ]
) &, inputs];
WriteString[$Output, "ALL-INPUT-CELLS-PARSED-AND-RAN\n"];
