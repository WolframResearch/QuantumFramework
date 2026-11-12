(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[matchingLookup]
PackageScope[decoderLookup]
PackageScope[codeMatching]
PackageScope[demMatching]
PackageScope[codeGraphlikeQ]
PackageScope[codeMatchingGraph]
PackageScope[demMatchingGraph]
PackageScope[codeNonGraphlikeLetter]
PackageScope[matchingDecodeIndices]
PackageScope[matchingDecodeKey]
PackageScope[matchingSpecQ]


(* ============================================================================ *)
(* Minimum-weight perfect matching (Got26 sec. 17.2.1-17.2.2).                   *)
(*                                                                              *)
(* The decoder the surface code is decoded with in practice, and the first one  *)
(* in this layer that does not enumerate.  Every decoder before it is a table   *)
(* -- over errors up to a weight, over cosets, over fault subsets -- and a      *)
(* table is what stops at d = 3.                                                *)
(*                                                                              *)
(* THE GRAPH.  When every elementary error flips at most two detectors, errors  *)
(* are edges and detectors are nodes; an error that flips one detector is an    *)
(* edge to a single virtual BOUNDARY node.  A syndrome is a set of lit nodes,   *)
(* and an error that explains it is a set of paths pairing them up (or sending  *)
(* them to the boundary).  With independent errors of probability p_e, the most  *)
(* likely such set is the one minimising Sum_e w_e, w_e = log((1 - p_e)/p_e).   *)
(* That is a minimum-weight perfect matching on the complete graph of lit       *)
(* nodes, weighted by shortest-path distance -- Edmonds' problem, solved        *)
(* exactly here by the built-in FindIndependentEdgeSet, which on a weighted     *)
(* graph returns a maximum-cardinality matching of maximum weight: give the     *)
(* edges C - w with C larger than any total and the matching is perfect and of  *)
(* minimum w.  (Checked against brute force before relying on it.)  Each lit    *)
(* node also gets its own copy of the boundary, the copies joined to each other *)
(* at weight 0, so any subset of nodes may go to the boundary independently.    *)
(*                                                                              *)
(* WHERE THE EDGES COME FROM.                                                   *)
(*  - Code capacity: the single-qubit errors X_q (probability p_X + p_Y) and     *)
(*    Z_q (p_Z + p_Y), with the syndrome bits as detectors and the logical      *)
(*    class bits as the observable.  This is the standard independent-X-and-Z  *)
(*    approximation: a Y is decoded as an X and a Z.  It needs every X_q and    *)
(*    Z_q to trip at most two checks (and every Y_q, on a code that is not     *)
(*    CSS), which holds for the repetition, surface, toric and Shor codes and   *)
(*    fails for the 5-qubit and Steane codes: those are refused, with a         *)
(*    message, at every noise level, rather than decoded wrongly.               *)
(*  - A detector model: its rows, with the rows a herald rejects left out (the  *)
(*    same filter demDecoderTable applies, for the same reason).  On a CSS code *)
(*    each row is first SPLIT into its X part and its Z part (splitByBasis),    *)
(*    since the two are seen by disjoint sets of detectors.  A part that still  *)
(*    trips more than two detectors -- a fault that spreads -- is               *)
(*    DECOMPOSED into edges the graph already has, whose observables must XOR   *)
(*    to its own; its probability is then added to each of them, which is what *)
(*    Stim's decompose_errors does for PyMatching.  A row with no consistent    *)
(*    decomposition is dropped and counted, not guessed at.                    *)
(*                                                                              *)
(* Parallel edges are merged: mechanisms with the same observable combine as    *)
(* independent flips, p + q - 2 p q; across different observables the more      *)
(* likely one is kept.  Symbolic rates are evaluated at $QECDecoderProbeRate, as *)
(* the maximum-likelihood decoder does, because a matching needs numbers.       *)
(*                                                                              *)
(* The decoder is lazy: it decodes the syndromes it is asked about and          *)
(* remembers the paths it found, so a sampled rate costs one matching per       *)
(* distinct syndrome in the sample rather than one per syndrome in 2^m.         *)
(*                                                                              *)
(* References: Got26 sec. 17.2 (decoding as matching); Dennis, Kitaev, Landahl, *)
(* Preskill, J. Math. Phys. 43, 4452 (2002); Edmonds, Can. J. Math. 17, 449     *)
(* (1965); Higgott, PyMatching, ACM TQC 3, 16 (2022).                           *)
(* ============================================================================ *)

QECDecoder::notgraphlike = "Matching needs every single-qubit error to trip at most two checks; on this code `1` trips `2`. Use \"MaximumLikelihood\" or \"MinimumWeight\" for it.";


(* ---- building the graph ---- *)

oddParity[p_, q_] := p + q - 2 p q

edgeKeyFor[dets_List, boundary_Integer] := If[Length[dets] === 1, {First[dets], boundary}, Sort[dets]]

edgeWeight[p_] := With[{q = Clip[N[p], {10^-300, 1/2 - 10^-12}]}, Log[(1 - q) / q]]

(* mechanisms: {detector indices, observable key, probability, payload}.  The payload is
   what the edge contributes to a correction: a symplectic vector at code capacity,
   0 for a detector model, where only the observable is ever reported. *)
matchingFromMechanisms[nd_Integer, mechs_List] := Module[
    {boundary = nd + 1, live, graphlike, hyper, edges, byEdge, keys, weights, graph,
     dist, cache, dropped = 0, decompose, extra = <||>},

    live = Select[mechs, Length[#[[1]]] > 0 && TrueQ[#[[3]] > 0] &];
    graphlike = Select[live, Length[#[[1]]] <= 2 &];
    hyper = Select[live, Length[#[[1]]] > 2 &];

    (* merge parallel mechanisms on each edge *)
    byEdge = GroupBy[graphlike, edgeKeyFor[#[[1]], boundary] &];
    edges = Map[
        Function[group,
            With[{byObs = GroupBy[group, #[[2]] &]},
                With[{best = First @ MaximalBy[
                        KeyValueMap[{#1, Fold[oddParity, 0, #2[[All, 3]]], First[MaximalBy[#2, #[[3]] &]][[4]]} &, byObs],
                        #[[2]] &]},
                    <|"Observable" -> best[[1]], "Probability" -> best[[2]], "Payload" -> best[[3]]|>
                ]
            ]
        ],
        byEdge
    ];

    (* split a hyperedge into existing edges whose observables XOR to its own *)
    decompose[dets_List, obs_] := Module[{go},
        go[{}, acc_] := If[BitXor @@ Append[acc[[All, 2]], 0] === obs, acc, $Failed];
        go[rest_List, acc_] := Module[{d = First[rest], tries},
            tries = Join[
                If[KeyExistsQ[edges, {d, boundary}], {{{d, boundary}, Rest[rest]}}, {}],
                Table[
                    If[KeyExistsQ[edges, Sort[{d, e}]], {Sort[{d, e}], DeleteCases[rest, d | e]}, Nothing],
                    {e, Rest[rest]}
                ]
            ];
            Catch[
                Scan[
                    Function[t, With[{r = go[t[[2]], Append[acc, {t[[1]], edges[t[[1]]]["Observable"]}]]},
                        If[r =!= $Failed, Throw[r, "decomposed"]]]],
                    tries
                ];
                $Failed,
                "decomposed"
            ]
        ];
        go[dets, {}]
    ];
    Scan[
        Function[h, With[{parts = decompose[h[[1]], h[[2]]]},
            If[ parts === $Failed,
                dropped++,
                Scan[(extra[#[[1]]] = oddParity[Lookup[extra, Key[#[[1]]], 0], h[[3]]]) &, parts]
            ]
        ]],
        hyper
    ];
    KeyValueMap[(edges[[Key[#1], "Probability"]] = oddParity[edges[[Key[#1], "Probability"]], #2]) &, extra];

    keys = Keys[edges];
    weights = edgeWeight[#["Probability"]] & /@ Values[edges];
    graph = Graph[Range[boundary], UndirectedEdge @@@ keys, EdgeWeight -> weights];
    dist = GraphDistanceMatrix[graph];

    (* shortest paths, found once per pair and remembered *)
    cache = Unique["matchingPath"];
    With[{c = cache, g = graph, e = edges},
        c[u_, v_] := c[u, v] = With[{path = FindShortestPath[g, u, v]},
            With[{steps = Sort /@ Partition[path, 2, 1]},
                {BitXor @@ Append[e[#]["Observable"] & /@ steps, 0],
                 If[steps === {}, 0, BitXor @@ (e[#]["Payload"] & /@ steps)]}
            ]
        ]
    ];

    <|
        "Detectors" -> nd, "Boundary" -> boundary, "Graph" -> graph, "Edges" -> edges,
        "Distances" -> dist, "Paths" -> cache, "DroppedMechanisms" -> dropped
    |>
]


(* ---- decoding ---- *)

(* The lit detectors in, {observable key, payload} out, or Missing when no perfect
   matching exists (a lit detector with no path to anything). *)
matchingDecodeIndices[m_Association, {}] := {0, 0}

matchingDecodeIndices[m_Association, lit_List] := Module[
    {k = Length[lit], dist = m["Distances"], b = m["Boundary"], pairs, weights, big, matching, pieces},
    pairs = Join[
        Flatten[Table[{{i, j}, dist[[lit[[i]], lit[[j]]]]}, {i, k}, {j, i + 1, k}], 1],
        Table[{{i, k + i}, dist[[lit[[i]], b]]}, {i, k}],
        Flatten[Table[{{k + i, k + j}, 0}, {i, k}, {j, i + 1, k}], 1]
    ];
    pairs = Select[pairs, NumericQ[#[[2]]] && #[[2]] < Infinity &];
    big = 1 + Total[pairs[[All, 2]]];
    matching = FindIndependentEdgeSet[
        Graph[Range[2 k], UndirectedEdge @@@ pairs[[All, 1]], EdgeWeight -> (big - pairs[[All, 2]])]
    ];
    If[Length[matching] < k, Return[Missing["UndecodableSyndrome"]]];
    pieces = Replace[List @@@ matching, {
        {i_, j_} /; i <= k && j <= k :> m["Paths"][lit[[i]], lit[[j]]],
        {i_, j_} /; i <= k :> m["Paths"][lit[[i]], b],
        {i_, j_} /; j <= k :> m["Paths"][lit[[j]], b],
        _ -> Nothing
    }, {1}];
    If[pieces === {}, {0, 0}, {BitXor @@ pieces[[All, 1]], BitXor @@ pieces[[All, 2]]}]
]

(* A syndrome or detector key -- the integer whose binary digits, most significant
   first, are the detector bits -- to the predicted observable key. *)
matchingDecodeKey[m_Association, key_Integer] := Replace[
    matchingDecodeIndices[m, Flatten @ Position[IntegerDigits[key, 2, m["Detectors"]], 1]],
    {{obs_, _} :> obs}
]

(* The rate machinery reads decoders as key -> key associations.  A matching decoder
   is materialised only on the keys asked about; a key it cannot decode is left out,
   which the rate counts as a failure, exactly as for a table that does not cover it. *)
decoderLookup[d_Association, _] := d
decoderLookup[matchingLookup[m_Association], keys_List] :=
    DeleteMissing @ AssociationMap[matchingDecodeKey[m, #] &, DeleteDuplicates[keys]]


matchingSpecQ["Matching"] := True
matchingSpecQ[QECDecoder[d_Association]] := d["Method"] === "Matching"
matchingSpecQ[_] := False


(* ---- the two sources of a graph ---- *)

numericProbabilities[noise_Association] := With[{p = noise["Probabilities"]},
    If[ noiseSymbolicQ[noise],
        N[p /. Thread[noise["Parameters"] -> $QECDecoderProbeRate]],
        N[p]
    ]
]

(* The test comes first and is not memoised with the graph, so a refused code says
   why every time it is asked, not only the first. *)
codeMatching[a_Association, noise_Association] := If[codeGraphlikeQ[a], codeMatchingGraph[a, noise], $Failed]

codeMatchingGraph[a_Association, noise_Association] := codeMatchingGraph[a, noise] = Module[
    {n = a["Qubits"], p, mechs, bad},
    p = numericProbabilities[noise];
    mechs = Flatten[Table[
        With[{row = ReplacePart[pauliIdentity[n], If[letter === "X", q, n + q] -> 1]},
            {Flatten @ Position[codeSyndromeVector[a, row], 1],
             First[Last[codeLabelKeys[a, {row}]]],
             If[letter === "X", p[[2]] + p[[3]], p[[4]] + p[[3]]],
             Most[row]}
        ],
        {q, n}, {letter, {"X", "Z"}}
    ], 1];
    matchingFromMechanisms[codeStabilizerCount[a], mechs]
]

(* WHICH CODES MATCH.  Matching pairs up defects, so it needs every single-qubit X and
   every single-qubit Z to trip at most two checks: then each is an edge, and a
   pattern of defects is a set of paths.  A code where one letter trips three -- the
   Steane code, any color code -- is a hypergraph, and pairing its defects is not
   decoding it.  On a code that is not CSS a Y is not an X and a Z seen by separate
   checks, so it has to be an edge too: that is what refuses the 5-qubit code, whose
   X's and Z's trip two checks each but whose Y's trip three.  This is a property of
   the code, not of the noise or the circuit, so both the code-capacity and the
   detector-model decoders ask it, and refuse the same codes for the same reason. *)
codeGraphlikeQ[a_Association] := With[{bad = codeNonGraphlikeLetter[a]},
    If[ bad === None,
        True,
        Message[QECDecoder::notgraphlike, pauliString[bad], Total[codeSyndromeVector[a, bad]]]; False
    ]
]

codeNonGraphlikeLetter[a_Association] := codeNonGraphlikeLetter[a] = With[
    {n = a["Qubits"]},
    SelectFirst[
        Flatten[Table[
            ReplacePart[pauliIdentity[n], Thread[positions -> 1]],
            {q, n},
            {positions, If[codeCSSQ[a], {{q}, {n + q}}, {{q}, {n + q}, {q, n + q}}]}
        ], 1],
        Total[codeSyndromeVector[a, #]] > 2 &,
        None
    ]
]

(* SPLITTING BY BASIS.  On a CSS code a fault's X part and Z part are seen by different
   checks: the X part trips Z-check detectors and flips the observables that read Zbar,
   the Z part trips X-check detectors and flips the ones that read Xbar.  A Y, or a
   two-qubit fault with both, is therefore two independent edges, not one edge between
   an X-check detector and a Z-check detector -- an edge like that is a shortcut the
   matching would take between two unrelated defects.  This is the decomposition Stim
   performs for PyMatching, done here from the detector's check type: detector
   (round, i) is check i, rounds major; the observable bits are the Xbar products first,
   then the Zbar ones.  A part that trips no detector is an undetectable logical flip,
   and no decoder can see it, so it is left out of the graph. *)
splitByBasis[mech_List, zCheckQ_List, k_Integer] := Module[
    {dets = mech[[1]], bits = IntegerDigits[mech[[2]], 2, 2 k], dz, dx, xPart, zPart},
    dz = Select[dets, zCheckQ[[#]] &];
    dx = Select[dets, ! zCheckQ[[#]] &];
    xPart = {dz, FromDigits[Join[ConstantArray[0, k], bits[[k + 1 ;; 2 k]]], 2], mech[[3]], 0};
    zPart = {dx, FromDigits[Join[bits[[1 ;; k]], ConstantArray[0, k]], 2], mech[[3]], 0};
    Select[{xPart, zPart}, #[[1]] =!= {} &]
]

demMatching[a_Association] := If[codeGraphlikeQ[a["Code"]], demMatchingGraph[a], $Failed]

demMatchingGraph[a_Association] := demMatchingGraph[a] = Module[
    {d, o, h, p, rows, mechs, code = a["Code"], m, zCheck},
    {d, o, h} = demRowKeys[a];
    p = If[AllTrue[a["Probabilities"], NumericQ], N[a["Probabilities"]],
        N[a["Probabilities"] /. Thread[DeleteDuplicates[Cases[a["Probabilities"], _Symbol, Infinity]] -> $QECDecoderProbeRate]]];
    rows = Select[Range[Length[d]], h[[#]] === 0 &];
    mechs = Table[{Flatten @ Position[a["DetectorMatrix"][[r]], 1], o[[r]], p[[r]], 0}, {r, rows}];
    If[ codeCSSQ[code],
        m = codeStabilizerCount[code];
        (* a Z check is one with no X part *)
        zCheck = With[{n = code["Qubits"]}, Total[#[[1 ;; n]]] === 0 & /@ code["CheckMatrix"]];
        zCheck = zCheck[[Mod[Range[a["Detectors"]] - 1, m] + 1]];
        mechs = Catenate[splitByBasis[#, zCheck, codeLogicalQubits[code]] & /@ mechs]
    ];
    matchingFromMechanisms[a["Detectors"], mechs]
]
