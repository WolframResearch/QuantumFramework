(* One evidence-backed invariant gate for each of the twenty catalog examples. *)

ClearAll[exampleCheck];
exampleCheck[number_Integer, name_String, held_HoldComplete] := Module[{value, label, messages},
  label = "example " <> ToString[number] <> " — " <> name;
  Block[{$MessageList = {}},
    value = Check[ReleaseHold[held], $Failed];
    messages = $MessageList
  ];
  If[TrueQ[value] && messages === {},
    WriteString[$Output, "PASS: ", label, "\n"],
    issueCount++;
    WriteString[$Output, "FAIL: ", label, " => ", ToString[value, InputForm],
      "; messages = ", ToString[messages, InputForm], "\n"]
  ]
];

exampleCheck[1, "pure-dephasing Bloch rates and stationary family", HoldComplete[
  Simplify[rates[0 id2, {Sqrt[g] Z}] == {-2 g x, -2 g y, 0}, g > 0] &&
    lindbladian[0 id2, {Sqrt[1/2] Z}, blochState[0, 0, 2/5]] == ConstantArray[0, {2, 2}]
]];

exampleCheck[2, "driven damped atom steady state", HoldComplete[
  With[{om = 3., gg = 1., rho = First@steadyState[(3./2) X, {lower}]},
    Norm[blochVector[rho] - {0., 2 om gg/(gg^2 + 2 om^2), -gg^2/(gg^2 + 2 om^2)}] < 10^-10 &&
      Norm[lindbladian[(om/2) X, {Sqrt[gg] lower}, rho], "Frobenius"] < 10^-10
  ]
]];

exampleCheck[3, "leaky-cavity Kraus and Liouvillian solutions", HoldComplete[
  Max@Abs@Flatten[evolve[0 blankCat, {Sqrt[\[Gamma]cat] downCat}, densityMatrix[cat], 0.6] -
      photonLoss[densityMatrix[cat], 0.6]] < 10^-9
]];

exampleCheck[4, "atom-made cat conditioning and separation-dependent decay", HoldComplete[
  Abs[purity[conditionalCat[1.3, 1]] - 1] < 10^-9 &&
    Last[far/First[far]] < Last[close/First[close]]
]];

exampleCheck[5, "dispersive-readout pointer formula", HoldComplete[
  Max@Table[traceDistance[cavityCat[tt, dispRungs], cavityGiven[dispRun[tt], dispRungs, plus]],
      {tt, 0.1, 3, 0.2}] < 10^-3
]];

exampleCheck[6, "high-temperature Brownian covariance remains physical", HoldComplete[
  Min@Table[Det[spread[tt]], {tt, 0, 25, 0.1}] >= 1/4 - 10^-8 &&
    Abs[particle[[3]][25] - kT/(M \[CapitalOmega]qbm^2)] < 0.02
]];

exampleCheck[7, "photon-counting antibunching and ensemble recovery", HoldComplete[
  Abs[waitDensity[0.]] < 10^-12 &&
    Max@Abs[(excitedPop /@ Mean[Table[densityMatrix /@
          Take[countRuns[[kk]]["states"], Length[when]], {kk, nCount}]]) - averageChance] < 0.12
]];

exampleCheck[8, "ideal homodyne path purity", HoldComplete[
  Max[Abs[purity[#] - 1] & /@ smoothHistory["states"]] < 10^-8
]];

exampleCheck[9, "heterodyne channel identity and independent innovations", HoldComplete[
  Norm[dissipator[cAtom/Sqrt[2], start] + dissipator[I cAtom/Sqrt[2], start] -
      dissipator[cAtom, start], "Frobenius"] < 10^-12 &&
    Abs[Correlation[dWI, dWQ]] < 0.1
]];

exampleCheck[10, "Zeno exceptional-point eigenvalue law", HoldComplete[
  Abs[zenoFreq[1/2] - Sqrt[3]] < 10^-10 && Abs[zenoRate[1/2] - 1] < 10^-10 &&
    Abs[zenoFreq[1]] < 10^-8
]];

exampleCheck[11, "QPC spectrum through the exceptional point", HoldComplete[
  With[{atEP = qpcSpectrum[4 OmQpc, {-1., 0., 1.}], wide = Range[-40, 40, 0.04]},
    VectorQ[atEP, NumericQ] && FreeQ[atEP, Indeterminate | ComplexInfinity] &&
      Abs[0.04 Total[qpcSpectrum[1., wide]]/(2 Pi) - 1] < 0.03
  ]
]];

exampleCheck[12, "localization martingale and finite-sample readout", HoldComplete[
  Max[Abs[meanLoc - z0Loc]] < 10^-8 && Min[Abs[Last[#]] & /@ fanLoc] > 0.99 && fidLoc[tfLoc] >= 0.99
]];

exampleCheck[13, "Kalman Riccati fixed point", HoldComplete[
  With[{s = Sqrt[1 + 4 kKal^2], exact = {1/Sqrt[2 (1 + Sqrt[1 + 4 kKal^2])],
       (Sqrt[1 + 4 kKal^2] - 1)/(4 kKal),
       Sqrt[1 + 4 kKal^2]/Sqrt[2 (1 + Sqrt[1 + 4 kKal^2])]}},
    Norm[Through[riccati[12.]] - exact] < 3 10^-4 && Abs[Det[covariance[12.]] - 1/4] < 10^-5
  ]
]];

exampleCheck[14, "measurement-feedback target and ensemble", HoldComplete[
  Norm[(Z - I Y) . plus] < 10^-12 && Last[predicted] > 0.95 &&
    Max[Abs[steeredMean - predicted]] < 0.15
]];

exampleCheck[15, "thermal mean and Gibbs populations", HoldComplete[
  Max[Abs[{units[fromCold[16.]], units[fromHot[16.]]} - nT]] < 0.01 &&
    Max[Abs[pops - warmShape]] < 10^-3
]];

exampleCheck[16, "rapid-purification refinement and uncensored passage", HoldComplete[
  Greater @@ purificationConvergence[[All, 3]] &&
    {Count[fptFixed, _Missing], Count[fptSteer, _Missing]} === {0, 0} &&
    Mean[fptFixed] < Mean[fptSteer]
]];

exampleCheck[17, "feedback cooling has a finite optimal gain", HoldComplete[
  nUncCool[Gcool] < nUncCool[0.] && nUncCool[Gcool] > condFloorCool &&
    nUncCool[12.] > nUncCool[Gcool]
]];

exampleCheck[18, "linear reference-measure trajectory average", HoldComplete[
  Norm[linearAverage - evolve[(\[CapitalOmega]lin/2) X, {cLin}, densityMatrix[excited], steps step],
      "Frobenius"] < 0.15 &&
    Abs[Mean[Re[Conjugate[#] . #] & /@ linearStates] - 1] < 0.2
]];

exampleCheck[19, "general driven atom and ideal watched purity", HoldComplete[
  Max@Abs[steadyGeneral[1.3, 3., 1., 0.4, 0.25] -
      blochVector[First@steadyState[(1.3/2) Z + (3./2) X,
        {Sqrt[1.4] lower, Sqrt[0.4] ConjugateTranspose[lower], Sqrt[0.25] Z}]]] < 10^-9 &&
    Max[Abs[(1 - Norm[blochVector[#]]^2)/2] & /@ watched["states"]] < 10^-8
]];

exampleCheck[20, "Mollow triplet and homodyne squeezing", HoldComplete[
  Min[strong] > -10^-8 && Length@Select[tones[[#]] & /@ FindPeaks[strong][[All, 1]],
      Abs[Abs[#] - 6] < 0.6 &] == 2 &&
    Min[mollowHom[0.25, 1., Pi/2, tones]] < 1 &&
    Max@Table[Abs[zeroFreqSqueeze[om] - (1 + 8 om^2 (4 om^2 - 1)/(1 + 2 om^2)^3)],
      {om, 0.05, 1.2, 0.05}] < 10^-9
]];
