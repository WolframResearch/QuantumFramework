(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[$QECTransformChunk]

PackageScope[syndromeVector]
PackageScope[whtTransform]
PackageScope[whtSpan]
PackageScope[whtLocationTables]
PackageScope[demTransformFailure]
PackageScope[demDistribution]
PackageScope[demMaximumLikelihoodTable]
PackageScope[demTableVector]
PackageScope[demSignatures]


(* ============================================================================ *)
(* Exact rates through the Walsh-Hadamard transform.                            *)
(*                                                                              *)
(* The exact route of Memory.wl folds every fault location into a distribution  *)
(* over x = (detectors, observables, heralds), a vector of length 2^N with N the *)
(* number of those bits.  With symbolic rates every entry of that vector is a    *)
(* polynomial, and the fold expands 2^N of them once per location: it is what    *)
(* capped the exact route at N = 16.                                             *)
(*                                                                              *)
(* The fold is a convolution over GF(2)^N, and the Walsh-Hadamard transform      *)
(* turns a convolution into a product.  For a location whose mutually exclusive  *)
(* outcomes have probabilities p_r and effects s_r (bit vectors over x),          *)
(*                                                                              *)
(*     f(k) = sum_x P_location(x) (-1)^(k.x) = 1 - 2 sum_(r : k.s_r odd) p_r,     *)
(*                                                                              *)
(* and the transform of the whole distribution is the product of these factors.  *)
(* A decoder's failure probability is the overlap of its failure indicator W     *)
(* with the distribution, and Parseval moves the overlap to the other side:       *)
(*                                                                              *)
(*     P_fail = 2^-N sum_k What(k) prod_locations f_location(k).                  *)
(*                                                                              *)
(* What is a transform of a 0/1 vector, so it is an integer vector computed in    *)
(* machine arithmetic with no rounding error (|What| <= 2^N < 2^53).  The product *)
(* is where the symbols live, and it takes very few distinct values: at a two-    *)
(* qubit depolarizing location exactly 0 or 8 of the 15 Paulis are odd under any  *)
(* k, so the factor is 1 or 1 - 16 p/15, and similarly elsewhere.  So the k are   *)
(* grouped by how many locations contribute each distinct factor, the integers   *)
(* What are summed within each group, and the symbolic part is a short sum of     *)
(* monomials -- one per group, not one per k.                                   *)
(*                                                                              *)
(* The idea comes from surface-code.nb (November 2026), which used it for the   *)
(* d = 3 surface code with one noisy round; here it is general: any detector      *)
(* model, any mutually exclusive outcomes per location, heralds included.         *)
(*                                                                              *)
(* MAKING IT FAST.  Everything per k is a pass over 2^N entries, so the cost is   *)
(* the number of passes.  The factor of a location depends on k only through the *)
(* parities of k against a basis of the span of its effects (at most four bits for *)
(* a CNOT), a "code".  Splitting k into a high and a low half, a parity is the XOR *)
(* of a parity of each half, so the codes of several locations at once come from  *)
(* two short tables and one Outer, and their joint contribution from one lookup.  *)
(* Locations are batched up to ten code bits, and k is walked in blocks of         *)
(* $QECTransformChunk so that no per-k array but the transform itself is ever     *)
(* held whole.                                                                    *)
(*                                                                              *)
(* THE SAME MACHINERY GIVES THE BEST DECODER.  With numeric rates the product of  *)
(* factors is the transform of the distribution itself, and one inverse transform *)
(* gives the probability of every (syndrome, observables) outcome.  The decoder   *)
(* that picks, for each syndrome, the most probable observable class is maximum   *)
(* likelihood over the whole detector history: no decoder does better on this     *)
(* model.  Its table is a vector over all 2^(detectors) syndromes, and the exact    *)
(* route then gives its rate -- as a polynomial, for a table chosen at one rate.  *)
(* In machine arithmetic the distribution is accurate to about 1e-16 absolute, so *)
(* the choice for syndromes rarer than that is noise; they carry no weight at the *)
(* rate the table is built for.                                                  *)
(* ============================================================================ *)

$QECTransformChunk = 2^22;

$whtBatchBits = 10;


(* ---- the transform ---- *)

(* Natural (Hadamard) ordering: entry k pairs with x through (-1)^popcount(k & x).
   DiscreteHadamardTransform normalises by 2^(-N/2); undoing that leaves integers
   for an integer input, exact in machine arithmetic below 2^53. *)
whtTransform[v_List, n_Integer] :=
    Round[2.^(n/2) DiscreteHadamardTransform[Developer`ToPackedArray[N[v]], Method -> "BitComplement"]]


(* ---- one location ---- *)

(* A basis of the span of the effects, taken from the effects themselves, and each
   effect's coordinates in it as a bit mask -- Gaussian elimination on the integers. *)
whtSpan[vs_List] := Module[{red = {}, piv = {}, msk = {}, kept = {}, coords = {}, v, m, nk, o},
    Do[
        v = u; m = 0;
        Do[
            If[BitGet[v, piv[[i]]] == 1, v = BitXor[v, red[[i]]]; m = BitXor[m, msk[[i]]]],
            {i, Length[red]}
        ];
        If[ v == 0,
            AppendTo[coords, m],
            nk = Length[kept];
            AppendTo[kept, u];
            AppendTo[red, v]; AppendTo[piv, BitLength[v] - 1]; AppendTo[msk, BitXor[m, 2^nk]];
            AppendTo[coords, 2^nk];
            o = Ordering[piv, All, Greater];
            red = red[[o]]; piv = piv[[o]]; msk = msk[[o]]
        ],
        {u, vs}
    ];
    {kept, coords}
]

(* {basis, factor for every code value}: the code's bit i is the parity of k against
   basis vector i, and an outcome is odd exactly when its coordinate mask meets the
   code in an odd number of bits. *)
whtLocationTables[sigs_List, probs_List, groups_List] := Function[rows,
    Module[{basis, coords},
        {basis, coords} = whtSpan[sigs[[rows]]];
        {
            basis,
            Table[
                Expand[1 - 2 Total[Pick[probs[[rows]], BitAnd[DigitCount[BitAnd[coords, c], 2, 1], 1], 1]]],
                {c, 0, 2^Length[basis] - 1}
            ]
        }
    ]
] /@ groups


(* ---- walking k ---- *)

(* The code of one location on one half of k. *)
whtHalfCode[basis_List, shift_Integer, bits_Integer] := With[{ks = Range[0, 2^bits - 1]},
    Total[MapIndexed[
        2^(First[#2] - 1) BitAnd[DigitCount[BitAnd[ks, BitAnd[BitShiftRight[#1, shift], 2^bits - 1]], 2, 1], 1] &,
        basis
    ]]
]

(* Batches of locations, each {bits, high-half codes, low-half codes, table}; the table
   is indexed by (high code, low code) and holds the batch's combined contribution, inc
   giving one location's contribution for each of its code values. *)
whtBatches[tables_List, n_Integer, inc_, combine_] := Module[
    {lo = Floor[n/2], hi, maxBits, batches = {}, cur = {}, curBits = 0},
    hi = n - lo;
    maxBits = Min[$whtBatchBits, Max[4, Floor[(n - 4)/2]]];
    Do[
        If[ t[[1]] =!= {},
            If[curBits + Length[t[[1]]] > maxBits && cur =!= {}, AppendTo[batches, cur]; cur = {}; curBits = 0];
            AppendTo[cur, t]; curBits += Length[t[[1]]]
        ],
        {t, tables}
    ];
    If[cur =!= {}, AppendTo[batches, cur]];
    Function[batch,
        Module[{offs = Most[FoldList[Plus, 0, Length[#[[1]]] & /@ batch]], bb, xs},
            bb = Total[Length[#[[1]]] & /@ batch];
            xs = With[{rr = Range[0, 2^bb - 1]},
                BitXor[Flatten[Transpose[ConstantArray[rr, 2^bb]]], Flatten[ConstantArray[rr, 2^bb]]]
            ];
            {
                bb,
                Total[MapThread[2^#2 whtHalfCode[#1[[1]], lo, hi] &, {batch, offs}]],
                Total[MapThread[2^#2 whtHalfCode[#1[[1]], 0, lo] &, {batch, offs}]],
                combine @@ MapThread[
                    inc[#1][[BitAnd[BitShiftRight[xs, #2], 2^Length[#1[[1]]] - 1] + 1]] &,
                    {batch, offs}
                ]
            }
        ]
    ] /@ batches
]

(* f[{first, last}, values] for consecutive blocks of k, values being the combination
   over all batches. *)
whtScan[batches_List, n_Integer, init_, combine_, f_] := Module[{lo = Floor[n/2], hi, step},
    hi = n - lo;
    step = Max[1, Quotient[$QECTransformChunk, 2^lo]];
    Do[
        With[{hs = Range[s + 1, Min[s + step, 2^hi]]},
            f[
                {s 2^lo + 1, Last[hs] 2^lo},
                Fold[
                    combine[#1, #2[[4]][[Flatten[Outer[Plus, 2^#2[[1]] #2[[2]][[hs]], #2[[3]]]] + 1]]] &,
                    init,
                    batches
                ]
            ]
        ],
        {s, 0, 2^hi - 1, step}
    ]
]


(* ---- the model's side ---- *)

(* The effect of each row as one integer over x = (detectors, observables, heralds),
   detectors most significant -- the same packing as demExactFailure's state index. *)
demSignatures[a_Association] := With[{keys = demRowKeys[a], no = a["Observables"], nh = Lookup[a, "Heralds", 0]},
    BitShiftLeft[keys[[1]], no + nh] + BitShiftLeft[keys[[2]], nh] + keys[[3]]
]

demBits[a_Association] := a["Detectors"] + a["Observables"] + Lookup[a, "Heralds", 0]

(* Group the k by how many locations contribute each distinct factor, and sum each
   integer vector in ws within the groups.  Returns the factors and, for every vector,
   an Association group -> sum.  A group is a mixed-radix integer of counts. *)
whtGroupSums[a_Association, ws_List] := Module[
    {n = demBits[a], tables, factors, ids, maxCount, radix, kmax, batches, sums, old},

    tables = whtLocationTables[demSignatures[a], a["Probabilities"], demGroups[a]];
    factors = DeleteCases[DeleteDuplicates[Catenate[tables[[All, 2]]]], 1];
    ids = AssociationThread[factors, Range[Length[factors]]];
    maxCount = Table[Count[tables[[All, 2]], t_ /; MemberQ[t, f]], {f, factors}];
    radix = FoldList[Times, 1, Most[maxCount + 1]];
    kmax = If[factors === {}, 1, Last[radix] (Last[maxCount] + 1)];
    If[kmax > 2^60, Return[$Failed, Module]];

    batches = whtBatches[tables, n, Function[t, Prepend[radix, 0][[Lookup[ids, t[[2]], 0] + 1]]], Plus];
    sums = ConstantArray[<||>, Length[ws]];
    old = SystemOptions["SparseArrayOptions" -> "TreatRepeatedEntries"];
    Internal`WithLocalSettings[
        (* repeated positions of a SparseArray are summed: a bucket sum without sorting *)
        SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}],
        If[ batches === {},
            sums = Function[w, <|0 -> Total[w]|>] /@ ws,
            whtScan[batches, n, 0, Plus, Function[{span, keys},
                sums = MapThread[
                    Merge[{#1, Association[
                        Most[ArrayRules[SparseArray[Transpose[{keys + 1}] -> #2[[span[[1]] ;; span[[2]]]], {kmax}]]] /.
                            ({k_} -> v_) :> (k - 1 -> v)
                    ]}, Total] &,
                    {sums, ws}
                ]
            ]]
        ],
        SetSystemOptions[old]
    ];
    {factors, radix, maxCount, sums}
]

whtMonomial[key_Integer, factors_List, radix_List, maxCount_List] :=
    Times @@ (factors ^ Mod[Quotient[key, radix], maxCount + 1])

whtCollect[{factors_, radix_, maxCount_, _}, sums_Association, n_Integer] :=
    Expand[2^-n Total[KeyValueMap[#2 whtMonomial[#1, factors, radix, maxCount] &, sums]]]


(* The decoder as a vector over all syndromes, -1 where it gives no answer. *)
demTableVector[a_Association, syndromeVector[v_List]] := v
demTableVector[a_Association, table_] := With[{keys = Range[0, 2^a["Detectors"] - 1]},
    Lookup[decoderLookup[table, keys], keys, -1]
]

(* {failure and acceptance, acceptance} probabilities of a decoder table, exactly:
   symbolic if the rates are.  $Failed when the factors are too varied to group, which
   needs a model whose locations nearly all have different rates. *)
demTransformFailure[a_Association, table_] := Module[
    {nd = a["Detectors"], no = a["Observables"], nh = Lookup[a, "Heralds", 0], n, vec, wrong, accept, ws, g},
    n = nd + no + nh;
    vec = demTableVector[a, table];
    (* failure over (syndrome, observables), then the heralds: only x with no herald
       fired is accepted *)
    wrong = Unitize[Outer[Subtract, vec, Range[0, 2^no - 1]]];
    ws = If[ nh === 0,
        {whtTransform[Flatten[wrong], n]},
        {
            whtTransform[Flatten[Outer[Times, wrong, UnitVector[2^nh, 1]]], n],
            whtTransform[Flatten[ConstantArray[UnitVector[2^nh, 1], 2^(nd + no)]], n]
        }
    ];
    g = whtGroupSums[a, ws];
    If[FailureQ[g], Return[g, Module]];
    whtCollect[g, #, n] & /@ Last[g]
]


(* ---- numeric distribution and maximum likelihood ---- *)

(* The probability of every x at numeric rates: the product of factors per k, then
   the inverse transform. *)
demDistribution[a_Association, rules_] := Module[{n = demBits[a], tables, batches, dhat},
    tables = whtLocationTables[demSignatures[a], N[a["Probabilities"] /. rules], demGroups[a]];
    batches = whtBatches[tables, n, Function[t, Developer`ToPackedArray[N[t[[2]]]]], Times];
    dhat = ConstantArray[1., 2^n];
    whtScan[batches, n, 1., Times, Function[{span, vals}, dhat[[span[[1]] ;; span[[2]]]] = vals]];
    2.^(-n/2) DiscreteHadamardTransform[dhat, Method -> "BitComplement"]
]

(* For every syndrome the observable class with the most accepted probability, the
   first on ties.  Symbolic rates are read at $QECDecoderProbeRate, unless a reference
   point is given as rules or as one number for every rate parameter. *)
demMaximumLikelihoodTable[a_Association] := demMaximumLikelihoodTable[a, Automatic]

demMaximumLikelihoodTable[a_Association, ref_] := qecMemoKey[demMaximumLikelihoodTable[a, ref]] = Module[
    {params, rules, dist, m, best, val, better, no = a["Observables"], nh = Lookup[a, "Heralds", 0]},
    params = DeleteDuplicates[Cases[a["Probabilities"], s_Symbol /; ! NumericQ[s], Infinity]];
    rules = Switch[ref,
        Automatic, Thread[params -> $QECDecoderProbeRate],
        _?NumericQ, Thread[params -> ref],
        _, ref
    ];
    dist = demDistribution[a, rules];
    m = ArrayReshape[dist, {2^a["Detectors"], 2^no, 2^nh}][[All, All, 1]];
    best = ConstantArray[0, Length[m]];
    val = m[[All, 1]];
    Do[
        better = Ramp[Sign[m[[All, j]] - val]];
        best = best + better (j - 1 - best);
        val = val + better (m[[All, j]] - val),
        {j, 2, 2^no}
    ];
    syndromeVector[best]
]

decoderLookup[syndromeVector[v_List], keys_List] := AssociationThread[keys, v[[keys + 1]]]
