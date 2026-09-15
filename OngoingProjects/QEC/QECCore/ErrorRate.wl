(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECLogicalErrorRate]
PackageExport[$QECExactEnumerationLimit]

PackageScope[codeLabelMatrix]
PackageScope[$codeExtractions]        (* Measurement.wl *)
PackageScope[codeErrorLabels]
PackageScope[codeErrorTally]
PackageScope[codeCosetRepresentatives]
PackageScope[cosetProbabilities]
PackageScope[codeLabelKeys]
PackageScope[codeMaximumLikelihoodDecoder]
PackageScope[codeMinimumWeightDecoder]
PackageScope[decoderFor]
PackageScope[decoderReachFor]
PackageScope[exactLogicalErrorRate]
PackageScope[sampledLogicalErrorRate]


(* ============================================================================ *)
(* Decoding by inference, and the figure of merit.                              *)
(*                                                                              *)
(* One observation organises this whole file.  Label every n-qubit Pauli error   *)
(* by two things: which checks it anticommutes with (its syndrome, n-k bits),    *)
(* and which logical operators it anticommutes with (its logical class, 2k bits). *)
(* Both are the same kind of question, so both are one matrix product over GF(2) *)
(* against the checks stacked on top of the logical operators, halves swapped.    *)
(*                                                                              *)
(* That label decides everything:                                                *)
(*                                                                              *)
(*   - Two errors with the same label differ by a stabilizer, so a correction     *)
(*     for one corrects the other.  The label's fibre IS the coset.               *)
(*   - Decoding is: given the syndrome half, guess the logical half.  The best    *)
(*     guess is the class carrying the most probability -- maximum likelihood.    *)
(*     Minimum weight is the p -> 0 limit of that, and on a degenerate code like   *)
(*     Shor the two genuinely differ: many light errors in one class can outweigh  *)
(*     a single lighter error in another.                                         *)
(*   - Correction succeeds exactly when the guessed class was right.  So the      *)
(*     logical error rate needs no residual algebra per trial, only a comparison  *)
(*     of labels.                                                                 *)
(*                                                                                *)
(* And because the noise is i.i.d., the probability of an error depends only on   *)
(* its {number of X, of Y, of Z}.  So one pass over the errors, tallying counts   *)
(* per (syndrome, class, profile), answers every question about this code under   *)
(* every noise model and every rate, symbolic or numeric.  The pass is memoised    *)
(* on the code alone; changing p costs nothing.                                    *)
(* ============================================================================ *)

QECLogicalErrorRate::usage = "QECLogicalErrorRate[code, noise] gives the probability that a correction cycle leaves a logical error, computed exactly by enumeration. With a symbolic noise rate the result is an exact expression in it.\nQECLogicalErrorRate[code, noise, n] estimates the same quantity from n sampled errors.\nThe option \"Decoder\" chooses which decoder the number describes: \"MaximumLikelihood\" (the most probable coset, needs the full enumeration), \"MinimumWeight\" (the lookup table, cheap at any size), or Automatic, which takes maximum likelihood when the code can be enumerated and minimum weight when it cannot.";

QECLogicalErrorRate::toobig = "Exact enumeration needs 4^`1` = `2` errors, above the limit $QECExactEnumerationLimit = `3`. Give a sample count as a third argument to estimate it instead, or raise the limit.";

QECLogicalErrorRate::nological = "The code has no logical qubits, so there is no logical error to make.";

QECLogicalErrorRate::mwdecoder = "4^`1` errors is above $QECExactEnumerationLimit, so the estimate describes the minimum-weight lookup decoder rather than maximum likelihood. Pass \"Decoder\" -> \"MinimumWeight\" to choose it deliberately, or raise the limit for maximum likelihood.";

QECLogicalErrorRate::reachone = "The code has `1` qubits, too many to establish its distance cheaply, so the minimum-weight table reaches weight one only. Pass \"DecoderReach\" -> w for a deeper table.";

QECLogicalErrorRate::mlneedsall = "Maximum-likelihood decoding weighs every coset, so it needs all 4^`1` = `2` errors, above $QECExactEnumerationLimit = `3`. Use \"Decoder\" -> \"MinimumWeight\", or raise the limit.";

$QECExactEnumerationLimit = 4^10;

(* Where the decoder's argmax is resolved when the noise rate is symbolic.  A
   decoder is a fixed function; it has to be chosen before the rate is a free
   variable.  Choosing it at a small rate is both the physically relevant regime
   and the one where maximum likelihood agrees with minimum weight.  Once chosen,
   the error rate reported for it is exact for every p. *)
$QECDecoderProbeRate = 1/1000;


(* ---- the label ---- *)

(* Checks stacked on logical operators, X and Z halves exchanged, so that
   row . labelMatrix^T is the vector of symplectic products. *)
codeLabelMatrix[a_Association] := codeLabelMatrix[a] = With[
    {n = a["Qubits"], rows = Join[a["CheckMatrix"], symplecticPart /@ Join @@ Values[codeLogicalVectors[a]]]},
    Join[rows[[All, n + 1 ;; 2 n]], rows[[All, 1 ;; n]], 2]
]

(* {syndrome bits, class bits} for a list of Pauli rows, in one product. *)
codeErrorLabels[a_Association, rows_List] :=
    Mod[(symplecticPart /@ rows) . Transpose[codeLabelMatrix[a]], 2]

(* Both halves of the label as integers, from one product: {syndrome keys, class keys}.
   Computing them together matters -- the sampled route calls this on every draw. *)
codeLabelKeys[a_Association, rows_List] := With[
    {m = codeStabilizerCount[a], labels = codeErrorLabels[a, rows]},
    {FromDigits[#, 2] & /@ labels[[All, 1 ;; m]], FromDigits[#, 2] & /@ labels[[All, m + 1 ;; -1]]}
]

(* All 4^n Pauli rows on n qubits, phase 0. *)
allPauliRows[n_Integer] := With[
    {pairs = Tuples[Values[$pauliLetterXZ], n]},
    Join[pairs[[All, All, 1]], pairs[[All, All, 2]], ConstantArray[0, {4^n, 1}], 2]
]


(* ---- the one pass ---- *)

(* <|{synKey, classKey, nx, ny, nz} -> count|> over every error, plus a
   minimum-weight representative per (synKey, classKey).  Memoised on the code:
   it is a fact about the code, not about any noise model. *)
codeErrorTally[a_Association] := codeErrorTally[a] = Module[
    {n = a["Qubits"], m = codeStabilizerCount[a], rows, labels, synKeys, classKeys, x, z, nx, ny, nz, weights, order, keys},

    rows = allPauliRows[n];
    labels = codeErrorLabels[a, rows];

    synKeys = FromDigits[#, 2] & /@ labels[[All, 1 ;; m]];
    classKeys = If[m === Length[First[labels]], ConstantArray[0, Length[labels]], FromDigits[#, 2] & /@ labels[[All, m + 1 ;; -1]]];

    x = rows[[All, 1 ;; n]];
    z = rows[[All, n + 1 ;; 2 n]];
    nx = Total[x (1 - z), {2}];
    ny = Total[x z, {2}];
    nz = Total[(1 - x) z, {2}];
    weights = nx + ny + nz;

    keys = Transpose[{synKeys, classKeys, nx, ny, nz}];

    (* Heaviest first, so that building the Association leaves the lightest
       representative in place (a later key overwrites an earlier one). *)
    order = Ordering[weights, All, Greater];

    <|
        "Counts" -> Counts[keys],
        "Representatives" -> Association[Thread[Transpose[{synKeys[[order]], classKeys[[order]]}] -> rows[[order]]]],
        "Qubits" -> n,
        "Checks" -> m
    |>
]

exactEnumerableQ[a_Association] := 4^a["Qubits"] <= $QECExactEnumerationLimit


(* ---- maximum-likelihood decoding ---- *)

(* Probability of each (syndrome, class) under the noise, as an exact expression. *)
cosetProbabilities[a_Association, noise_Association] := cosetProbabilities[a, noise] = Module[
    {tally = codeErrorTally[a], n},
    n = tally["Qubits"];
    Merge[
        KeyValueMap[
            {#1[[1]], #1[[2]]} -> #2 noiseProfileProbability[noise, n, #1[[3 ;; 5]]] &,
            tally["Counts"]
        ],
        Total
    ]
]

(* syndrome -> the class carrying the most probability.  With a symbolic rate the
   comparison is made at $QECDecoderProbeRate; see the note above. *)
codeMaximumLikelihoodDecoder[a_Association, noise_Association] := codeMaximumLikelihoodDecoder[a, noise] = Module[
    {probs = cosetProbabilities[a, noise], concrete, grouped},

    concrete = If[
        noiseSymbolicQ[noise],
        With[{rules = Thread[noise["Parameters"] -> $QECDecoderProbeRate]}, N[probs /. rules]],
        N[probs]
    ];

    (* Keys of `concrete` are {syndrome, class} pairs; group them by syndrome and keep
       the class of the heaviest.  The value is the class alone, not the pair. *)
    grouped = GroupBy[Keys[concrete], First];

    Association @ KeyValueMap[
        #1 -> Last[First[MaximalBy[#2, concrete[#] &, 1]]] &,
        grouped
    ]
]

(* The same syndrome -> class map, read off the minimum-weight lookup table instead.

   This is the decoder for codes that cannot be enumerated.  The table only needs the
   errors up to the weight the code guarantees (see codeDecoderToWeight), not all 4^n,
   so it is cheap at any size; what it gives up is the weighing of whole cosets, which
   is what maximum likelihood adds and what pays off under biased noise.

   A syndrome the table does not cover has no correction, and an uncorrectable syndrome
   is a failure, so it is simply absent here and the caller counts it as one. *)
codeMinimumWeightDecoder[a_Association, reach_Integer] := codeMinimumWeightDecoder[a, reach] = With[
    {corrections = Values[codeDecoderToWeight[a, reach]]},
    With[{keys = codeLabelKeys[a, corrections]},
        AssociationThread[First[keys] -> Last[keys]]
    ]
]

(* How far up the table reaches.  Deliberately NOT the code's guaranteed weight by
   default on a large code: that weight is floor((d-1)/2), and knowing d means running
   the distance search, whose only bound is the weight of the standard-form logical
   operators.  When those come out heavy -- the [[n,n-2,2]] family is the case that
   caught this -- the bound bounds nothing and the search degenerates to nearly the
   full 4^n it was supposed to avoid.  So a code small enough to enumerate uses its
   real guaranteed weight, and a larger one defaults to weight one and says so; raise
   it explicitly with "DecoderReach" when the code deserves more, at the cost of
   enumerating every error up to that weight. *)
decoderReachFor[a_Association, Automatic] := If[
    exactEnumerableQ[a],
    codeDecoderReach[a],
    Message[QECLogicalErrorRate::reachone, a["Qubits"]]; 1
]

decoderReachFor[_Association, w_Integer ? Positive] := w

(* Which decoder a call describes.  Automatic prefers maximum likelihood and falls back
   to minimum weight when the enumeration it needs is out of reach -- saying so, because
   the two answer different questions and a silent switch would make the number mean
   something the caller did not ask for. *)
decoderFor[a_Association, noise_Association, spec_, reach_] := Switch[spec,
    "MaximumLikelihood",
        If[ exactEnumerableQ[a],
            codeMaximumLikelihoodDecoder[a, noise],
            Message[QECLogicalErrorRate::mlneedsall, a["Qubits"], 4^a["Qubits"], $QECExactEnumerationLimit];
            $Failed
        ],
    "MinimumWeight",
        codeMinimumWeightDecoder[a, reach],
    _,
        If[ exactEnumerableQ[a],
            codeMaximumLikelihoodDecoder[a, noise],
            Message[QECLogicalErrorRate::mwdecoder, a["Qubits"]];
            codeMinimumWeightDecoder[a, reach]
        ]
]

(* The correction the decoder applies for a syndrome: a minimum-weight member of
   the chosen coset. *)
codeCosetRepresentatives[a_Association, noise_Association] := codeCosetRepresentatives[a, noise] = With[
    {tally = codeErrorTally[a], decoder = codeMaximumLikelihoodDecoder[a, noise]},
    Association @ KeyValueMap[#1 -> tally["Representatives"][{#1, #2}] &, decoder]
]


(* ---- the figure of merit ---- *)

(* Exact: sum, over syndromes, of the probability the decoder guessed the class
   wrong.  Every term is exact, so a symbolic rate gives a polynomial. *)
exactLogicalErrorRate[a_Association, noise_Association, decoder_Association] :=
    Simplify @ Total @ KeyValueMap[
        If[#1[[2]] === decoder[#1[[1]]], 0, #2] &,
        cosetProbabilities[a, noise]
    ]

(* Sampled: draw errors, label them, and count the ones whose class the decoder
   would have missed.  No per-trial GF(2) work -- the label already decides it. *)
sampledLogicalErrorRate[a_Association, noise_Association, count_Integer, decoder_Association] := Module[
    {rows, synKeys, classKeys},

    rows = noiseSampleRows[noise, a["Qubits"], count];
    {synKeys, classKeys} = codeLabelKeys[a, rows];

    (* A syndrome absent from the decoder is one it cannot correct, which is a failure;
       the -1 default never equals a real class key, so those draws count as such. *)
    N @ Divide[
        Count[Transpose[{synKeys, classKeys}], {s_, c_} /; Lookup[decoder, s, -1] =!= c],
        count
    ]
]


Options[QECLogicalErrorRate] = {"Decoder" -> Automatic, "DecoderReach" -> Automatic,
    "Rounds" -> Automatic, "Extraction" -> "BareAncilla"};

QECLogicalErrorRate[code_QECCode, noise_QECNoiseModel, opts : OptionsPattern[]] :=
    logicalErrorRate[First[code], First[noise], None,
        OptionValue[QECLogicalErrorRate, {opts}, "Decoder"], OptionValue[QECLogicalErrorRate, {opts}, "DecoderReach"],
        OptionValue[QECLogicalErrorRate, {opts}, "Rounds"],
        OptionValue[QECLogicalErrorRate, {opts}, "Extraction"]]

QECLogicalErrorRate[code_QECCode, noise_QECNoiseModel, count_Integer ? Positive, opts : OptionsPattern[]] :=
    logicalErrorRate[First[code], First[noise], count,
        OptionValue[QECLogicalErrorRate, {opts}, "Decoder"], OptionValue[QECLogicalErrorRate, {opts}, "DecoderReach"],
        OptionValue[QECLogicalErrorRate, {opts}, "Rounds"],
        OptionValue[QECLogicalErrorRate, {opts}, "Extraction"]]

(* Rounds only mean something once a check can lie, so a code-capacity model ignores
   them and the other two levels default to the code distance -- the standard memory
   experiment, where the time direction is protected as well as the space ones. *)
roundsFor[a_Association, Automatic] := defaultRounds[a]
roundsFor[_Association, r_Integer ? Positive] := r

logicalErrorRate[a_Association, noise_Association, count_, spec_, reachSpec_, roundSpec_, mode_] /;
        noiseLevel[noise] =!= "CodeCapacity" := Module[{reach},
    If[ codeLogicalQubits[a] === 0,
        Message[QECLogicalErrorRate::nological]; Return[$Failed]
    ];
    If[ ! MemberQ[$codeExtractions, mode],
        Message[QECDetectorModel::extraction, mode, $codeExtractions]; Return[$Failed]
    ];
    reach = decoderReachFor[a, reachSpec];
    demRate[a, noise, roundsFor[a, roundSpec], count, reach, mode]
]

(* Code capacity has no circuit, so it has no extraction to choose; the option is
   accepted and ignored rather than refused, since a sweep across noise levels
   should not have to strip it. *)
logicalErrorRate[a_Association, noise_Association, count_, spec_, reachSpec_, _, _] := Module[{decoder, reach},
    Which[
        codeLogicalQubits[a] === 0,
            Message[QECLogicalErrorRate::nological]; Return[$Failed],
        (* the exact route sums over every coset, so it needs the enumeration whatever
           decoder it is describing *)
        count === None && ! exactEnumerableQ[a],
            Message[QECLogicalErrorRate::toobig, a["Qubits"], 4^a["Qubits"], $QECExactEnumerationLimit];
            Return[$Failed],
        IntegerQ[count] && noiseSymbolicQ[noise],
            Message[QECNoiseModel::rates, noise["Probabilities"]]; Return[$Failed]
    ];

    reach = decoderReachFor[a, reachSpec];
    decoder = decoderFor[a, noise, spec, reach];
    If[decoder === $Failed, Return[$Failed]];

    If[ count === None,
        exactLogicalErrorRate[a, noise, decoder],
        sampledLogicalErrorRate[a, noise, count, decoder]
    ]
]
