(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[$QECExactDetectorLimit]
PackageExport[$QECDecoderSubsetLimit]

PackageScope[demRowKeys]
PackageScope[demGroups]
PackageScope[demDecoderTable]
PackageScope[demSampledFailure]
PackageScope[demExactFailure]
PackageScope[demRate]


(* ============================================================================ *)
(* The memory experiment.                                                       *)
(*                                                                              *)
(* Prepare a codeword, run r rounds of noisy syndrome extraction, read the data  *)
(* out perfectly at the end, decode the whole detector history, and ask whether  *)
(* the logical qubit survived.  With the detector model already built this is     *)
(* GF(2) bookkeeping on top of it, and it comes in the same two routes as the     *)
(* code-capacity rate: exact, or sampled.                                        *)
(*                                                                              *)
(* The exact route is the interesting one, because it is the move a sampling      *)
(* simulator cannot make.  Every fault mechanism is a coin, and its effect is a    *)
(* fixed vector of detector bits and observable bits.  Rather than enumerate the   *)
(* 2^(number of faults) configurations -- hopeless -- carry the *distribution over *)
(* effects*, a vector of length 2^(detectors + observables), and fold the coins in *)
(* one at a time:                                                                 *)
(*                                                                              *)
(*     P'[s] = P[s] (1 - sum_i p_i) + sum_i p_i P[s xor sig_i]                    *)
(*                                                                              *)
(* the sum running over the mutually exclusive outcomes of one location, so a     *)
(* three-way depolarizing channel is handled exactly rather than as three          *)
(* independent coins.  Cost is (locations) x (state space), not 2^faults, and     *)
(* with symbolic rates the answer is an exact polynomial in p and q -- a           *)
(* circuit-level logical error rate in closed form.  The state space is what       *)
(* bounds it: 2^((r+1) m + 2k) grows fast, so this is for the small codes and few  *)
(* rounds, and $QECExactDetectorLimit says where it stops.                        *)
(*                                                                              *)
(* The decoder is minimum weight over fault mechanisms: the fewest faults whose    *)
(* detector pattern matches what was seen, precomputed into a lookup table the     *)
(* same way codeDecoderToWeight builds one over Pauli errors.  A detector pattern  *)
(* the table does not reach is undecodable and counts as a failure.  This is the   *)
(* honest built-in; matching decoders (PyMatching) plug in at exactly this point,  *)
(* reading the detector model and returning the same thing.                       *)
(*                                                                              *)
(* This is a different kind of decoder from the book's, not a variant of it.       *)
(* Gottesman's FTEC (QECC book sec. 12.2.3) uses the repetitions only as a         *)
(* certificate: find the last run of t+1 agreeing syndromes, conclude that one of  *)
(* them was taken with no faults, discard everything else, and look that single    *)
(* (n-k)-bit syndrome up in a table.  Correctness comes from an existence           *)
(* argument.  What happens here instead is an optimisation over the whole           *)
(* spacetime history: nothing is discarded, all rounds are used jointly, and the    *)
(* answer is the lightest fault configuration consistent with every detector.       *)
(* That is the strategy sec. 12.2.2 raises and then sets aside as hard to analyse   *)
(* in general, so it inherits none of that section's guarantees -- what is exact     *)
(* here is the error rate of a stated protocol, not a fault-tolerance proof.        *)
(* ============================================================================ *)

$QECExactDetectorLimit = 2^16;

$QECDecoderSubsetLimit = 2 * 10^6;

QECLogicalErrorRate::detbig = "The detector model has `1` detectors and `2` observables, so an exact sum needs a state space of `3`, past $QECExactDetectorLimit (`4`). Give a shot count to sample instead, or use fewer rounds.";
QECLogicalErrorRate::subsets = "Stopping the decoder table at weight `1`: weight `2` would need `3` fault subsets, past $QECDecoderSubsetLimit (`4`).";
QECLogicalErrorRate::symbolicshots = "A sampled rate needs numeric noise rates.";


(* ---- packing a detector model into integers ---- *)

(* Detector and observable patterns are bit vectors that are only ever XORed and
   compared, so they are carried as integers: one BitXor replaces a vector Mod. *)
bitsToInteger[bits_List] := If[bits === {}, 0, FromDigits[bits, 2]]

demRowKeys[a_Association] := demRowKeys[a] = {
    bitsToInteger /@ a["DetectorMatrix"],
    bitsToInteger /@ a["ObservableMatrix"]
}

(* Mechanisms that belong to the same physical location are alternatives, not
   independent coins: at most one of the three Paulis of a depolarizing channel
   happens.  The location label already carries that grouping. *)
demGroups[a_Association] := demGroups[a] = Values @ GroupBy[
    Range[Length[a["Probabilities"]]],
    a["Locations"][[#]] &
]


(* ---- the decoder ---- *)

(* Detector pattern -> the observable flip implied by the lightest fault set that
   explains it.  Built once per (model, reach). *)
demDecoderTable[a_Association, reach_Integer] := demDecoderTable[a, reach] = Module[
    {d, o, nf, table, stop},

    {d, o} = demRowKeys[a];
    nf = Length[d];
    table = <|0 -> 0|>;
    stop = reach;

    Do[
        If[ Binomial[nf, w] > $QECDecoderSubsetLimit,
            Message[QECLogicalErrorRate::subsets, w - 1, w, Binomial[nf, w], $QECDecoderSubsetLimit];
            stop = w - 1;
            Break[]
        ];
        Do[
            With[{key = BitXor @@ d[[sub]]},
                If[! KeyExistsQ[table, key], table[key] = BitXor @@ o[[sub]]]
            ],
            {sub, Subsets[Range[nf], {w}]}
        ],
        {w, 1, reach}
    ];

    table
]


(* ---- sampled ---- *)

(* One draw per location, all shots at once: the group's own probabilities plus the
   leftover mass for "nothing happened".  Mechanisms dropped from the model because
   they had no effect are simply part of that leftover, which is why the weights
   still add up. *)
demSampledFailure[a_Association, count_Integer, reach_Integer] := Module[
    {d, o, groups, table, dKeys, oKeys, detectors, observables, predicted},

    {d, o} = demRowKeys[a];
    groups = demGroups[a];
    table = demDecoderTable[a, reach];

    (* index 0 means no fault, so the key lists are padded at the front *)
    dKeys = Prepend[d, 0];
    oKeys = Prepend[o, 0];

    {detectors, observables} = Transpose @ Map[
        Function[rows,
            With[{draw = RandomChoice[
                    Append[a["Probabilities"][[rows]], 1 - Total[a["Probabilities"][[rows]]]] ->
                        Append[rows, 0],
                    count]},
                {dKeys[[draw + 1]], oKeys[[draw + 1]]}
            ]
        ],
        groups
    ];

    detectors = Fold[BitXor, ConstantArray[0, count], detectors];
    observables = Fold[BitXor, ConstantArray[0, count], observables];

    (* A pattern the table does not reach is undecodable, hence a failure; -1 is
       never a real observable key, so those shots count as one. *)
    predicted = Lookup[table, detectors, -1];

    N[Count[Transpose[{predicted, observables}], {x_, y_} /; x =!= y] / count]
]


(* ---- exact ---- *)

demExactFailure[a_Association, reach_Integer] := Module[
    {d, o, no, nd, dim, sigs, groups, dist, idx, table, mask},

    nd = a["Detectors"];
    no = a["Observables"];
    dim = 2^(nd + no);

    If[ dim > $QECExactDetectorLimit,
        Message[QECLogicalErrorRate::detbig, nd, no, dim, $QECExactDetectorLimit];
        Return[$Failed]
    ];

    {d, o} = demRowKeys[a];
    sigs = BitShiftLeft[d, no] + o;
    groups = demGroups[a];
    idx = Range[0, dim - 1];

    dist = ConstantArray[0, dim];
    dist[[1]] = 1;

    Do[
        dist = Expand[
            dist (1 - Total[a["Probabilities"][[rows]]]) +
            Sum[a["Probabilities"][[r]] dist[[BitXor[idx, sigs[[r]]] + 1]], {r, rows}]
        ],
        {rows, groups}
    ];

    table = demDecoderTable[a, reach];

    (* 1 on the states the decoder gets wrong. *)
    mask = Table[
        Boole[Lookup[table, Key[BitShiftRight[s, no]], -1] =!= BitAnd[s, 2^no - 1]],
        {s, idx}
    ];

    Simplify[mask . dist]
]


(* ---- the entry point for the two noisy levels ---- *)

demRate[a_Association, noise_Association, rounds_Integer, count_, reach_Integer] := Module[{dem},
    dem = codeDetectorModel[a, noise, rounds];
    If[ count === None,
        demExactFailure[dem, reach],
        If[ noiseSymbolicQ[noise],
            Message[QECLogicalErrorRate::symbolicshots]; $Failed,
            demSampledFailure[dem, count, reach]
        ]
    ]
]
