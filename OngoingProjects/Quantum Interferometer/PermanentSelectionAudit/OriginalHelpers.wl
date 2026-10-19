(* Verbatim snapshot of QI helpers; see metadata.json for source hash. *)
modeList[occupation_] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occupation]

photonPatterns[n_, m_] := ReverseSort @ FrobeniusSolve[ConstantArray[1, m], n]

(* <T| U |S> = Perm(W_{T,S}) / Sqrt[Prod s! Prod t!]; Ryser over the multiplicities x <= s, Prod (s + 1) terms of
   Sum_x (-1)^(N - |x|) Prod Binomial[s, x] Prod_j (W_j . x)^t_j, run N digits above the input since they cancel *)
repeatedPermanent[w_, s_, t_] := With[{rows = Flatten @ Position[t, _ ? Positive]},
    Total[Times @@ Binomial[s, #] (-1) ^ (Total[s] - Total[#]) Times @@ ((w[[rows]] . #) ^ t[[rows]]) & /@ Tuples[Range[0, #] & /@ s]]
]

permanentAmplitude[w_, s_, t_] := Which[
    Total[s] != Total[t], 0,
    Total[s] == 0, 1,
    Precision[w] =!= Infinity && Min[Times @@ (s + 1), Times @@ (t + 1)] < 2 ^ Total[s],
        With[{v = SetPrecision[w, Precision[w] + Total[s]]},
            N[If[Times @@ (s + 1) <= Times @@ (t + 1), repeatedPermanent[v, s, t], repeatedPermanent[Transpose[v], t, s]] / Sqrt[Times @@ (s!) Times @@ (t!)], Precision[w]]
        ],
    True, Permanent[w[[modeList[t], modeList[s]]]] / Sqrt[Times @@ (s!) Times @@ (t!)]
]

