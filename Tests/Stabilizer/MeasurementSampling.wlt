(* ==========================================================================
   Tests/Stabilizer/MeasurementSampling.wlt -- Method -> {"Stabilizer",
   "Shots" -> n} on circuits with computational-basis measurements.

   "Shots" -> n returns n records (bit lists, in the default method's order:
   by measured wire, time order within a wire); "States" -> True returns n
   rules record -> PauliStabilizer. Every sampled record must be one exact mode
   gives, and its tableau must be that record's post-measurement state.
   ========================================================================== *)

Get[FileNameJoin[{ParentDirectory[PacletObject["Wolfram/QuantumFramework"]["Location"]], "Tests", "Stabilizer", "fixtures", "MeasurementHelpers.wl"}]]

stabSampM = QuantumMeasurementOperator;

(* Every sampled record is in the exact result, with the same state. The input *)
(* is a list of leading arguments: {} for the default register, {qs} or {ps}.  *)
stabSampConsistentQ[qc_QuantumCircuitOperator, input_List, shots_Integer] := With[{
    exact = qc @@ Append[input, Method -> "Stabilizer"],
    sampled = qc @@ Append[input, Method -> {"Stabilizer", "Shots" -> shots, "States" -> True}]
},
    Length[sampled] == shots &&
    AllTrue[sampled, Function[rule,
        KeyExistsQ[exact, First[rule]] &&
        stabMeasExactSameQ[stabMeasTableauDensity[Last[rule]], stabMeasTableauDensity[exact[First[rule]]]]
    ]]
]

stabSampBell = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, stabSampM[{1}], stabSampM[{2}]}];


(* ---------- records only ---------- *)

VerificationTest[
    BlockRandom[SeedRandom[11];
        With[{r = stabSampBell[Method -> {"Stabilizer", "Shots" -> 200}]},
            {Length[r], Union[r], Union[Length /@ r]}
        ]
    ],
    {200, {{0, 0}, {1, 1}}, {2}},
    TestID -> "Sampling-Bell-RecordsOnly"
]

(* impossible records never appear *)
VerificationTest[
    BlockRandom[SeedRandom[12];
        FreeQ[stabSampBell[Method -> {"Stabilizer", "Shots" -> 500}], {0, 1} | {1, 0}, {1}]
    ],
    True,
    TestID -> "Sampling-ImpossibleRecordsNeverAppear"
]

(* record order: by wire, time order within wire 3 *)
VerificationTest[
    BlockRandom[SeedRandom[13];
        Union @ QuantumCircuitOperator[{stabSampM[{3}], stabSampM[{1}], "X" -> 3, stabSampM[{3}], stabSampM[{2}], stabSampM[{4}]}][
            QuantumState["0110"], Method -> {"Stabilizer", "Shots" -> 5}]
    ],
    {{0, 1, 1, 0, 0}},
    TestID -> "Sampling-RecordOrder"
]


(* ---------- record with state ---------- *)

VerificationTest[
    BlockRandom[SeedRandom[14];
        {
            stabSampConsistentQ[stabSampBell, {QuantumState["00"]}, 40],
            stabSampConsistentQ[QuantumCircuitOperator[{"H" -> 1, stabSampM[{1}], "H" -> 1, stabSampM[{1}]}], {QuantumState["0"]}, 40],
            stabSampConsistentQ[QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CZ" -> {1, 2}, stabSampM[{1}], "S" -> 2}], {QuantumState["00"]}, 40],
            stabSampConsistentQ[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, stabSampM[{2}], "H" -> 1, stabSampM[{1}]}], {QuantumState["00"]}, 40],
            stabSampConsistentQ[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}, stabSampM[{3, 1}], "H" -> 2, stabSampM[{2}]}], {QuantumState["000"]}, 40]
        }
    ],
    {True, True, True, True, True},
    TestID -> "Sampling-StatesMatchExact-ClosedForms"
]

stabSampRandomCircuit[seed_] := BlockRandom[
    SeedRandom[seed];
    With[{n = RandomInteger[{3, 5}]},
        QuantumCircuitOperator @ Catenate @ Table[
            Append[Table[stabMeasRandomGate[n], RandomInteger[{1, 5}]],
                If[RandomReal[] < 0.3, stabSampM[RandomSample[Range[n], 2]], stabSampM[{RandomInteger[{1, n}]}]]],
            RandomInteger[{1, 4}]
        ]
    ]
]

VerificationTest[
    BlockRandom[SeedRandom[15];
        Select[Range[10], ! stabSampConsistentQ[stabSampRandomCircuit[#], {}, 30] &]
    ],
    {},
    TestID -> "Sampling-StatesMatchExact-RandomCircuits"
]


(* ---------- random stream ---------- *)

VerificationTest[
    {
        SeedRandom[21]; stabSampBell[Method -> {"Stabilizer", "Shots" -> 50}],
        SeedRandom[21]; stabSampBell[Method -> {"Stabilizer", "Shots" -> 50}]
    } // Apply[SameQ],
    True,
    TestID -> "Sampling-Reproducible"
]

(* records-only and with-states sampling draw the same records *)
VerificationTest[
    With[{
        records = (SeedRandom[22]; QuantumCircuitOperator[{"H" -> 1, stabSampM[{1}], "H" -> 1, stabSampM[{1}], "H" -> 2, stabSampM[{2}]}][Method -> {"Stabilizer", "Shots" -> 50}]),
        rules = (SeedRandom[22]; QuantumCircuitOperator[{"H" -> 1, stabSampM[{1}], "H" -> 1, stabSampM[{1}], "H" -> 2, stabSampM[{2}]}][Method -> {"Stabilizer", "Shots" -> 50, "States" -> True}])
    },
        records === Keys[rules]
    ],
    True,
    TestID -> "Sampling-OneSamplerForBothSettings"
]

(* the caller's stream advances with the draws, and zero shots leave it untouched *)
VerificationTest[
    {
        (SeedRandom[23]; stabSampBell[Method -> {"Stabilizer", "Shots" -> 5}]; RandomInteger[10^9]) =!= (SeedRandom[23]; RandomInteger[10^9]),
        (SeedRandom[23]; stabSampBell[Method -> {"Stabilizer", "Shots" -> 0}]; RandomInteger[10^9]) === (SeedRandom[23]; RandomInteger[10^9])
    },
    {True, True},
    TestID -> "Sampling-RandomStream"
]


(* ---------- frequencies (secondary) ---------- *)

(* H on three wires then measure all: 8 records at 1/8; Pearson chi-square with *)
(* 7 degrees of freedom below 30 (its 0.9999 quantile is 29.88)                 *)
VerificationTest[
    BlockRandom[SeedRandom[31];
        With[{counts = Lookup[Counts[QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "H" -> 3, stabSampM[{1}], stabSampM[{2}], stabSampM[{3}]}][
                Method -> {"Stabilizer", "Shots" -> 4000}]], Tuples[{0, 1}, 3], 0]},
            Total[(counts - 500)^2 / 500] < 30
        ]
    ],
    True,
    TestID -> "Sampling-ChiSquare"
]


(* ---------- contract ---------- *)

VerificationTest[
    {
        stabSampBell[Method -> {"Stabilizer", "Shots" -> 0}],
        stabSampBell[Method -> {"Stabilizer", "Shots" -> 0, "States" -> True}]
    },
    {{}, {}},
    TestID -> "Sampling-ZeroShots"
]

VerificationTest[
    #["Reason"] & /@ {
        QuantumCircuitOperator[{"H" -> 1}][Method -> {"Stabilizer", "Shots" -> 3}],
        stabSampBell[Method -> {"Stabilizer", "States" -> True}],
        stabSampBell[Method -> {"Stabilizer", "Shots" -> -1}],
        stabSampBell[Method -> {"Stabilizer", "Shots" -> 1.5}],
        stabSampBell[Method -> {"Stabilizer", "Shots" -> 3, "States" -> "yes"}]
    },
    {"ShotsWithoutMeasurement", "StatesWithoutShots", "InvalidOption", "InvalidOption", "InvalidOption"},
    TestID -> "Sampling-Refusals"
]

(* without a measurement, the options at their defaults change nothing, and an *)
(* invalid "MaxBranches" is refused                                             *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]},
        {
            qc[Method -> {"Stabilizer", "Shots" -> None, "States" -> False}] === qc[Method -> "Stabilizer"],
            qc[PauliStabilizer[2], Method -> {"Stabilizer", "Shots" -> None}] === qc[PauliStabilizer[2], Method -> "Stabilizer"],
            qc[Method -> {"Stabilizer", "MaxBranches" -> -3}]["Reason"]
        }
    ],
    {True, True, "InvalidOption"},
    TestID -> "Sampling-MeasurementFreeOptions"
]

(* "MaxBranches" bounds exact mode only: sampling 2^20 possible records is fine *)
VerificationTest[
    BlockRandom[SeedRandom[32];
        Length @ QuantumCircuitOperator[Join[Table["H" -> i, {i, 20}], stabSampM[{#}] & /@ Range[20]]][
            Method -> {"Stabilizer", "Shots" -> 10, "MaxBranches" -> 2}]
    ],
    10,
    TestID -> "Sampling-MaxBranchesIgnored"
]


(* ---------- large n, inside memory and time limits far below a dense state ---------- *)

stabSampLimited[expr_] := MemoryConstrained[TimeConstrained[expr, 300, "TimedOut"], 2 * 10^9, "OutOfMemory"]
SetAttributes[stabSampLimited, HoldFirst]

VerificationTest[
    BlockRandom[SeedRandom[41];
        stabSampLimited @ With[{r = QuantumCircuitOperator[Join[Table["H" -> i, {i, 50}], stabSampM[{#}] & /@ Range[50]]][
                Method -> {"Stabilizer", "Shots" -> 1000}]},
            {Length[r], Union[Length /@ r], Length[Union[r]] > 990}
        ]
    ],
    {1000, {50}, True},
    TestID -> "Sampling-Large-IndependentWires"
]

(* GHZ_50: all bits equal; for record b the state is fixed by Z_i Z_(i+1) = +1 and Z_1 = (-1)^b *)
VerificationTest[
    BlockRandom[SeedRandom[42];
        stabSampLimited @ With[{r = QuantumCircuitOperator[Join[stabMeasGHZ[50], stabSampM[{#}] & /@ Range[50]]][
                Method -> {"Stabilizer", "Shots" -> 20, "States" -> True}]},
            AllTrue[r, Function[rule,
                With[{b = First[First[rule]], ps = Last[rule]},
                    MatchQ[First[rule], {0 ..} | {1 ..}] &&
                    AllTrue[Range[49], Keys[ps["M", stabMeasZ[50, {#, # + 1}]]] === {0} &] &&
                    Keys[ps["M", stabMeasZ[50, {1}]]] === {b}
                ]
            ]]
        ]
    ],
    True,
    TestID -> "Sampling-Large-GHZ50"
]

(* graph state on 30 sites, then X measurement (H, then Z) of the 15 odd sites:  *)
(* on a ring every even site has two odd neighbours, so the product of the odd  *)
(* sites' generators is the product of their X and the parity is even; on an    *)
(* open chain site 30 has one odd neighbour and both parities occur.            *)
stabSampGraph[edges_] := Join[Table["H" -> i, {i, 30}], "CZ" -> # & /@ edges, Table["H" -> i, {i, 1, 29, 2}], stabSampM[{#}] & /@ Range[1, 29, 2]]

VerificationTest[
    BlockRandom[SeedRandom[43];
        stabSampLimited @ {
            Union[Mod[Total /@ QuantumCircuitOperator[stabSampGraph[Table[{i, Mod[i, 30] + 1}, {i, 30}]]][Method -> {"Stabilizer", "Shots" -> 200}], 2]],
            Union[Mod[Total /@ QuantumCircuitOperator[stabSampGraph[Table[{i, i + 1}, {i, 29}]]][Method -> {"Stabilizer", "Shots" -> 200}], 2]]
        }
    ],
    {{0}, {0, 1}},
    TestID -> "Sampling-Large-RingAndChainParity"
]
