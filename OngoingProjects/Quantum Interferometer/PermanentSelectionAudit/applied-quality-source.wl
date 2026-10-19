modeList[occupation_] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occupation]

photonPatterns[n_, m_] := ReverseSort @ FrobeniusSolve[ConstantArray[1, m], n]

(* <T| U |S> = Perm(W_{T,S}) / Sqrt[Prod s! Prod t!]; Ryser over the multiplicities x <= s, Prod (s + 1) terms of
   Sum_x (-1)^(N - |x|) Prod Binomial[s, x] Prod_j (W_j . x)^t_j, run N digits above the input since they cancel *)
repeatedPermanent[w_, s_, t_] := With[{rows = Flatten @ Position[t, _ ? Positive]},
    Total[Times @@ Binomial[s, #] (-1) ^ (Total[s] - Total[#]) Times @@ ((w[[rows]] . #) ^ t[[rows]]) & /@ Tuples[Range[0, #] & /@ s]]
]

(* Compare calibrated work scores rather than grouped term counts alone. Promoting the whole
   transfer matrix costs m^2; the interpreted grouped sum also depends on its occupied rows.
   Weights 16 (machine) and 1/10 (arbitrary precision) were calibrated on WL 15.0 ARM;
   see OngoingProjects/Quantum Interferometer/PermanentSelectionAudit for the benchmark data. *)
useRepeatedPermanentQ[w_, s_, t_] := With[{p = Precision[w]},
    p =!= Infinity && With[{n = Total[s], m = Length[s], a = Times @@ (s + 1), b = Times @@ (t + 1)},
        With[{k = Min[a, b]},
            k < 2 ^ n && With[{r = Count[If[a <= b, t, s], _ ? Positive], weight = If[p === MachinePrecision, 16, 1/10]},
                weight (m ^ 2 + k m (r + 1)) < n 2 ^ n
            ]
        ]
    ]
]

permanentAmplitude[w_, s_, t_] := Which[
    Total[s] != Total[t], 0,
    Total[s] == 0, 1,
    useRepeatedPermanentQ[w, s, t],
        With[{v = SetPrecision[w, Precision[w] + Total[s]]},
            N[If[Times @@ (s + 1) <= Times @@ (t + 1), repeatedPermanent[v, s, t], repeatedPermanent[Transpose[v], t, s]] / Sqrt[Times @@ (s!) Times @@ (t!)], Precision[w]]
        ],
    True, Permanent[w[[modeList[t], modeList[s]]]] / Sqrt[Times @@ (s!) Times @@ (t!)]
]

