$HistoryLength = 0;
Get["/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/OngoingProjects/Quantum Interferometer/PermanentSelectionAudit/Proposal.wl"];
Begin["PermanentSelectionQuality`"];
referenceMatrix[m_] := referenceMatrix[m] = BlockRandom[
 SeedRandom[99173 + m];
 First @ QRDecomposition[N[RandomInteger[{-20, 20}, {m, m}] + I RandomInteger[{-20, 20}, {m, m}], 85]]
];
coefficientAmplitude[w_, s_, t_] := With[{
 outputModes = Flatten @ Position[t, _?Positive],
 inputModes = Flatten @ Position[s, _?Positive]
},
 Module[{variables},
  variables = Array[Unique["referenceVariable"] &, Length[outputModes]];
  Fold[Coefficient[#1, First[#2], Last[#2]] &,
   Expand[Times @@ MapThread[(variables . #1)^#2 &, {Transpose[w[[outputModes, inputModes]]], s[[inputModes]]}]],
   Transpose[{variables, t[[outputModes]]}]] Sqrt[Times @@ (t!) / Times @@ (s!)]
 ]
];
rows = Map[Function[spec,
 With[{n = spec[[1]], m = spec[[2]], p = spec[[3]]},
  With[{s = PadRight[{n/2, n/2}, m], wHigh = referenceMatrix[m]},
   With[{t = Reverse[s], w = If[p === "Machine", N[wHigh], N[wHigh, p]]},
    With[{
     reference = coefficientAmplitude[wHigh, s, t],
     result = PermanentSelectionAudit`ProposedAmplitude[w, s, t]
    },
     With[{
      absoluteError = N[Abs[Rationalize[result, 0] - Rationalize[reference, 0]], 20],
      referenceMagnitude = N[Abs[reference], 20],
      tolerance = If[p === "Machine", 10^-10, 10^-25]
     },
      With[{relativeError = absoluteError / referenceMagnitude},
       <|"Photons" -> n, "Modes" -> m, "RequestedPrecision" -> p,
        "InputPrecision" -> N[Precision[w]],
        "Selection" -> PermanentSelectionAudit`SelectionPlan[w, s, t],
        "ReferencePrecision" -> N[Precision[reference]], "ResultPrecision" -> N[Precision[result]],
        "AbsReference" -> referenceMagnitude, "AbsoluteError" -> absoluteError,
        "RelativeError" -> relativeError, "RelativeTolerance" -> tolerance,
        "Pass" -> TrueQ[If[referenceMagnitude > 10^-60, relativeError < tolerance, absoluteError < 10^-60]]|>
      ]
     ]
    ]
   ]
  ]
 ]
], Tuples[{{10, 18}, {24, 96}, {"Machine", 35}}]];
Export["/tmp/permanent-rule-quality/relative-accuracy.json", rows, "RawJSON"];
Echo[rows, "RELATIVE_ACCURACY"];
Echo[{Length[rows], Count[Lookup[rows, "Pass"], True], Count[Lookup[rows, "Pass"], False]}, "RELATIVE_ACCURACY_COUNTS"];
End[];
PermanentSelectionQuality`rows
