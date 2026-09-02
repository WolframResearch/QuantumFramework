(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Noise.wlt

   The noise model, maximum-likelihood decoding, and the logical error rate.

   Two independent oracles carry this file.  The engine's own QuantumChannel,
   applied to a stabilizer state, returns the {probability, state} mixture, and
   our declared probabilities must equal it symbolically.  And the repetition
   code under bit-flip noise fails exactly when a majority of qubits flips, so
   its logical error rate must equal the binomial tail term for term -- a
   closed-form answer that owes nothing to this package.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

(* Transitional: the rebuilt QEC core still lives under OngoingProjects/QEC/.
   Once it moves into Kernel/QEC/ and PacletInfo.wl lists its context, this Get
   disappears and the Needs above is enough.  The repo root is found from the
   loaded paclet, which RunTests.wls points at this checkout through
   PacletDirectoryLoad, so the tests always run against the source tree. *)
Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

(* Derived properties are memoised on the code's data.  A kernel that already
   held results from an earlier load would let those outrank the definitions just
   read, so the tests would silently check the old package. *)
QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

(* Failure probability of the n-qubit repetition code under independent bit
   flips, from first principles: the majority is wrong. *)
qecBinomialTail[n_, p_] := Sum[Binomial[n, j] p^j (1 - p)^(n - j), {j, Floor[n/2] + 1, n}];


(* ============================================================================
   Building a noise model
   ============================================================================ *)

VerificationTest[Head[QECNoiseModel["Depolarizing", 1/10]], QECNoiseModel, TestID -> "QEC-Noise-head"]

VerificationTest[
    QECNoiseModel["Depolarizing", qecP]["Probabilities"],
    {1 - qecP, qecP/3, qecP/3, qecP/3},
    TestID -> "QEC-Noise-depolarizing-rates"
]

VerificationTest[
    QECNoiseModel["BitFlip", qecP]["Probabilities"],
    {1 - qecP, qecP, 0, 0},
    TestID -> "QEC-Noise-bitflip-rates"
]

VerificationTest[
    QECNoiseModel["PhaseFlip", qecP]["Probabilities"],
    {1 - qecP, 0, 0, qecP},
    TestID -> "QEC-Noise-phaseflip-rates"
]

VerificationTest[
    QECNoiseModel["BitPhaseFlip", qecP]["Probabilities"],
    {1 - qecP, 0, qecP, 0},
    TestID -> "QEC-Noise-bitphaseflip-rates"
]

VerificationTest[
    QECNoiseModel[<|"X" -> 1/100, "Z" -> 2/100|>]["Probabilities"],
    {97/100, 1/100, 0, 2/100},
    TestID -> "QEC-Noise-general-pauli"
]

VerificationTest[QECNoiseModel["Depolarizing", qecP]["SymbolicQ"], True, TestID -> "QEC-Noise-symbolic-flag"]

VerificationTest[QECNoiseModel["Depolarizing", 1/10]["SymbolicQ"], False, TestID -> "QEC-Noise-numeric-flag"]

VerificationTest[QECNoiseModel["Depolarizing", qecP]["Parameters"], {qecP}, TestID -> "QEC-Noise-parameters"]

VerificationTest[QECNoiseModel["Depolarizing", qecP]["MeanWeight"], qecP, TestID -> "QEC-Noise-mean-weight"]

VerificationTest[Quiet[QECNoiseModel["NoSuchModel", 1/10]], $Failed, TestID -> "QEC-Noise-unknown-model"]

VerificationTest[
    Quiet[Check[QECNoiseModel["NoSuchModel", 1/10], "fired", QECNoiseModel::unknown]],
    "fired",
    TestID -> "QEC-Noise-unknown-message"
]

(* Rates that cannot be probabilities are rejected; symbolic ones are taken on
   trust, since that is the point of allowing them. *)
VerificationTest[Quiet[QECNoiseModel[{1/2, 9/10, 0, 0}]], $Failed, TestID -> "QEC-Noise-reject-bad-rates"]

VerificationTest[Head[QECNoiseModel[{1 - qecP, qecP, 0, 0}]], QECNoiseModel, TestID -> "QEC-Noise-accept-symbolic-rates"]


(* ============================================================================
   Error probabilities

   P(E) depends only on how many X, Y and Z the error carries, never on where
   they sit.  Everything exact downstream rests on that.
   ============================================================================ *)

VerificationTest[
    QECNoiseModel["Depolarizing", qecP]["ErrorProbability", "IIIII"],
    (1 - qecP)^5,
    TestID -> "QEC-Noise-probability-identity"
]

VerificationTest[
    Simplify[QECNoiseModel["Depolarizing", qecP]["ErrorProbability", "XIIII"] - (qecP/3) (1 - qecP)^4],
    0,
    TestID -> "QEC-Noise-probability-weight-one"
]

(* Position cannot matter. *)
VerificationTest[
    Module[{noise = QECNoiseModel["Depolarizing", qecP], errors},
        errors = {"XYIII", "IXYII", "IIXYI", "YXIII"};
        Length[DeleteDuplicates[Simplify[noise["ErrorProbability", #]] & /@ errors]]
    ],
    1,
    TestID -> "QEC-Noise-probability-depends-only-on-profile"
]

(* A channel with a Pauli switched off still gives an ordinary product: the
   zero-probability factors have multiplicity zero, an empty product, not 0^0. *)
VerificationTest[
    QECNoiseModel["BitFlip", qecP]["ErrorProbability", "XIII"],
    qecP (1 - qecP)^3,
    TestID -> "QEC-Noise-probability-with-zero-rate"
]

VerificationTest[
    QECNoiseModel["BitFlip", qecP]["ErrorProbability", "ZIII"],
    0,
    TestID -> "QEC-Noise-impossible-error"
]

(* Over all 4^n errors the probabilities sum to one. *)
VerificationTest[
    Module[{noise = QECNoiseModel["Depolarizing", qecP], all},
        all = QECPauliString /@ Flatten[Table[Wolfram`QuantumFramework`QEC`PackageScope`weightKVectors[3, w], {w, 0, 3}], 1];
        Simplify[Total[noise["ErrorProbability", #] & /@ all]]
    ],
    1,
    TestID -> "QEC-Noise-probabilities-normalise"
]


(* ============================================================================
   Agreement with the engine's own channels

   The strongest check available: apply the engine's QuantumChannel to a
   stabilizer state and read back the mixture it produces.  Note the
   depolarizing convention differs -- the engine's parameter q gives each Pauli
   q/4, this model's p spreads p over three Paulis -- so the constructor
   converts.  If that conversion were wrong, this test fails.
   ============================================================================ *)

VerificationTest[
    Module[{noise, mixture},
        noise = QECNoiseModel["Depolarizing", qecP];
        mixture = noise["QuantumChannel", {1}][PauliStabilizer[1]];
        Simplify[Sort[mixture[[All, 1]]] == Sort[noise["Probabilities"]]]
    ],
    True,
    TestID -> "QEC-Noise-depolarizing-matches-engine"
]

VerificationTest[
    Module[{noise, mixture},
        noise = QECNoiseModel["BitFlip", qecP];
        mixture = noise["QuantumChannel", {1}][PauliStabilizer[1]];
        Simplify[Sort[mixture[[All, 1]]] == Sort[DeleteCases[noise["Probabilities"], 0]]]
    ],
    True,
    TestID -> "QEC-Noise-bitflip-matches-engine"
]

VerificationTest[
    Module[{noise, mixture},
        noise = QECNoiseModel["PhaseFlip", qecP];
        mixture = noise["QuantumChannel", {1}][PauliStabilizer[1]];
        Simplify[Sort[mixture[[All, 1]]] == Sort[DeleteCases[noise["Probabilities"], 0]]]
    ],
    True,
    TestID -> "QEC-Noise-phaseflip-matches-engine"
]

(* The bit-flip channel really does flip the sign of the Z stabilizer. *)
VerificationTest[
    Module[{mixture},
        mixture = QECNoiseModel["BitFlip", qecP]["QuantumChannel", {1}][PauliStabilizer[2]];
        Sort[mixture[[All, 2]] /. ps_PauliStabilizer :> First[ps["Stabilizers"]]]
    ],
    Sort[{"ZI", "-ZI"}],
    TestID -> "QEC-Noise-bitflip-flips-the-sign"
]


(* ============================================================================
   Sampling
   ============================================================================ *)

VerificationTest[
    Module[{},
        SeedRandom[1];
        Length[QECNoiseModel["Depolarizing", 1/2]["RandomErrors", 5, 20]]
    ],
    20,
    TestID -> "QEC-Noise-sample-count"
]

VerificationTest[
    Module[{},
        SeedRandom[2];
        AllTrue[QECNoiseModel["Depolarizing", 1/2]["RandomErrors", 6, 30], StringLength[#] === 6 &]
    ],
    True,
    TestID -> "QEC-Noise-sample-width"
]

(* A bit-flip channel can only produce X errors. *)
VerificationTest[
    Module[{},
        SeedRandom[3];
        AllTrue[QECNoiseModel["BitFlip", 1/2]["RandomErrors", 6, 50], StringFreeQ[#, "Y" | "Z"] &]
    ],
    True,
    TestID -> "QEC-Noise-sample-respects-channel"
]

(* Sampled weights track the model's mean. *)
VerificationTest[
    Module[{sample},
        SeedRandom[4];
        sample = QECNoiseModel["Depolarizing", 1/10]["RandomErrors", 20, 5000];
        Abs[Mean[N[QECPauliWeight /@ sample]] - 2] < 0.15
    ],
    True,
    TestID -> "QEC-Noise-sample-mean-weight"
]


(* ============================================================================
   The logical error rate, exactly

   The repetition code oracle: under independent bit flips the code fails
   exactly when more than half the qubits flip.  Term for term against the
   binomial tail, for three sizes, and again for the phase-flip dual.
   ============================================================================ *)

VerificationTest[
    Simplify[QECLogicalErrorRate[QECCode["Repetition", 3], QECNoiseModel["BitFlip", qecP]] == qecBinomialTail[3, qecP]],
    True,
    TestID -> "QEC-Rate-repetition-3-vs-binomial"
]

VerificationTest[
    Simplify[QECLogicalErrorRate[QECCode["Repetition", 5], QECNoiseModel["BitFlip", qecP]] == qecBinomialTail[5, qecP]],
    True,
    TestID -> "QEC-Rate-repetition-5-vs-binomial"
]

VerificationTest[
    Simplify[QECLogicalErrorRate[QECCode["Repetition", 7], QECNoiseModel["BitFlip", qecP]] == qecBinomialTail[7, qecP]],
    True,
    TestID -> "QEC-Rate-repetition-7-vs-binomial"
]

VerificationTest[
    Simplify[QECLogicalErrorRate[QECCode["PhaseRepetition", 5], QECNoiseModel["PhaseFlip", qecP]] == qecBinomialTail[5, qecP]],
    True,
    TestID -> "QEC-Rate-phase-repetition-dual"
]

(* A distance-3 code corrects one error, so it cannot fail below weight two:
   the rate must vanish to second order. *)
VerificationTest[
    Module[{rate = QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Depolarizing", qecP]]},
        {rate /. qecP -> 0, Exponent[Normal[Series[rate, {qecP, 0, 4}]], qecP, Min]}
    ],
    {0, 2},
    TestID -> "QEC-Rate-five-qubit-vanishes-to-second-order"
]

(* The [[5,1,3]] code is perfect: every weight-one error is corrected and every
   weight-two error fails, so the leading coefficient is exactly the number of
   weight-two errors times their probability, 9 Binomial[5,2] / 9 = 10. *)
VerificationTest[
    Coefficient[Normal[Series[QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Depolarizing", qecP]], {qecP, 0, 2}]], qecP, 2],
    10,
    TestID -> "QEC-Rate-five-qubit-leading-coefficient"
]

VerificationTest[
    Module[{rate = QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Depolarizing", qecP]]},
        Exponent[Normal[Series[rate, {qecP, 0, 4}]], qecP, Min]
    ],
    2,
    TestID -> "QEC-Rate-steane-vanishes-to-second-order"
]

(* Rate zero at rate zero, and no code protects what it does not hold. *)
VerificationTest[
    QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Depolarizing", 0]],
    0,
    TestID -> "QEC-Rate-no-noise-no-error"
]

VerificationTest[
    Quiet[QECLogicalErrorRate[QECCode[QECCode["5QubitCode"]["CompletedGenerators"]], QECNoiseModel["Depolarizing", qecP]]],
    $Failed,
    TestID -> "QEC-Rate-no-logical-qubits"
]


(* ============================================================================
   The logical error rate, by sampling

   The two routes are independent implementations of the same quantity -- one
   enumerates, one draws -- so they check each other.
   ============================================================================ *)

VerificationTest[
    Module[{code = QECCode["5QubitCode"], exact, sampled, n = 100000},
        exact = N[QECLogicalErrorRate[code, QECNoiseModel["Depolarizing", 1/20]]];
        SeedRandom[11];
        sampled = QECLogicalErrorRate[code, QECNoiseModel["Depolarizing", 1/20], n];
        Abs[sampled - exact] <= 4 Sqrt[exact (1 - exact)/n]
    ],
    True,
    TestID -> "QEC-Rate-sampled-agrees-with-exact-five"
]

VerificationTest[
    Module[{code = QECCode["SteaneCode"], exact, sampled, n = 100000},
        exact = N[QECLogicalErrorRate[code, QECNoiseModel["Depolarizing", 1/20]]];
        SeedRandom[12];
        sampled = QECLogicalErrorRate[code, QECNoiseModel["Depolarizing", 1/20], n];
        Abs[sampled - exact] <= 4 Sqrt[exact (1 - exact)/n]
    ],
    True,
    TestID -> "QEC-Rate-sampled-agrees-with-exact-steane"
]

VerificationTest[
    Module[{code = QECCode["Repetition", 5], exact, sampled, n = 100000},
        exact = N[QECLogicalErrorRate[code, QECNoiseModel["BitFlip", 1/10]]];
        SeedRandom[13];
        sampled = QECLogicalErrorRate[code, QECNoiseModel["BitFlip", 1/10], n];
        Abs[sampled - exact] <= 4 Sqrt[exact (1 - exact)/n]
    ],
    True,
    TestID -> "QEC-Rate-sampled-agrees-with-exact-repetition"
]

(* The sampled route needs a number to draw from. *)
VerificationTest[
    Quiet[QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Depolarizing", qecP], 100]],
    $Failed,
    TestID -> "QEC-Rate-sampling-needs-a-numeric-rate"
]


(* ============================================================================
   Maximum-likelihood decoding

   Minimum weight is the p -> 0 limit of maximum likelihood, so under symmetric
   noise the two agree and there is nothing to show.  Under biased noise -- Z
   ten times likelier than X, as on real hardware -- they part company, and the
   inference wins.  That is what the decoder is for.
   ============================================================================ *)

VerificationTest[
    Module[{decoder = QECCode["SteaneCode"]["Decoder", QECNoiseModel["Depolarizing", 1/100]]},
        {AssociationQ[decoder], Length[decoder], AllTrue[Values[decoder], StringQ]}
    ],
    {True, 64, True},
    TestID -> "QEC-Decoder-covers-every-syndrome"
]

(* Every stored correction really does produce the syndrome it is stored under. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], decoder},
        decoder = code["Decoder", QECNoiseModel["Depolarizing", 1/100]];
        AllTrue[Normal[decoder], FromDigits[code["Syndrome", Last[#]], 2] === First[#] &]
    ],
    True,
    TestID -> "QEC-Decoder-corrections-match-their-syndromes"
]

VerificationTest[
    Module[{code = QECCode["5QubitCode"], decoder},
        decoder = code["Decoder", QECNoiseModel["Depolarizing", 1/100]];
        AllTrue[Normal[decoder], FromDigits[code["Syndrome", Last[#]], 2] === First[#] &]
    ],
    True,
    TestID -> "QEC-Decoder-corrections-match-their-syndromes-five"
]

(* Under biased noise the Steane code decodes strictly better by inference than
   by weight: 651/25 against 678/25 at leading order. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], biased, ml},
        biased = QECNoiseModel[<|"X" -> qecP/10, "Y" -> qecP/10, "Z" -> qecP|>];
        ml = QECLogicalErrorRate[code, biased];
        Coefficient[Normal[Series[ml, {qecP, 0, 2}]], qecP, 2]
    ],
    651/25,
    TestID -> "QEC-Decoder-biased-noise-steane"
]

(* The property on the code object and the free function are the same thing. *)
VerificationTest[
    Module[{code = QECCode["Repetition", 5], noise = QECNoiseModel["BitFlip", qecP]},
        Simplify[code["LogicalErrorRate", noise] - QECLogicalErrorRate[code, noise]]
    ],
    0,
    TestID -> "QEC-Rate-property-matches-function"
]

(* Adding qubits helps below threshold and hurts above it -- the fact the whole
   subject is about.  For repetition codes under bit flips the crossing is at
   p = 1/2 exactly. *)
VerificationTest[
    Module[{r3, r5},
        r3 = QECLogicalErrorRate[QECCode["Repetition", 3], QECNoiseModel["BitFlip", qecP]];
        r5 = QECLogicalErrorRate[QECCode["Repetition", 5], QECNoiseModel["BitFlip", qecP]];
        {N[(r5 - r3) /. qecP -> 1/10] < 0, N[(r5 - r3) /. qecP -> 9/10] > 0,
         Simplify[(r5 - r3) /. qecP -> 1/2]}
    ],
    {True, True, 0},
    TestID -> "QEC-Rate-repetition-threshold-is-one-half"
]


(* ============================================================================
   Choosing the decoder

   The sampled route exists for codes too large to enumerate, and for a while it
   did not deliver: it built a maximum-likelihood decoder, which weighs every
   coset and so needs all 4^n errors after all. These tests pin the fix. The
   second one is the regression guard: it runs on a code above the enumeration
   limit and must finish, which the old code could not do.
   ============================================================================ *)

VerificationTest[
    Quiet[QECLogicalErrorRate[QECCode["DistanceTwo", 12], QECNoiseModel["Depolarizing", 1/100]]],
    $Failed,
    TestID -> "QEC-Decoder-exact-refuses-above-the-limit"
]

(* The code has 12 qubits, so 4^12 errors is past the limit; the estimate must
   still come back, and quickly, on the minimum-weight decoder. *)
VerificationTest[
    Module[{r},
        SeedRandom[21];
        r = TimeConstrained[
            Quiet[QECLogicalErrorRate[QECCode["DistanceTwo", 12], QECNoiseModel["Depolarizing", 1/50], 20000]],
            30, $Aborted];
        {NumericQ[r], 0 < r < 1}
    ],
    {True, True},
    TestID -> "QEC-Decoder-sampled-works-above-the-limit"
]

(* And it says which decoder the number describes rather than switching quietly. *)
VerificationTest[
    Quiet[Check[
        QECLogicalErrorRate[QECCode["DistanceTwo", 12], QECNoiseModel["Depolarizing", 1/50], 500],
        "fired", QECLogicalErrorRate::mwdecoder]],
    "fired",
    TestID -> "QEC-Decoder-announces-the-fallback"
]

(* Asking for maximum likelihood on a code that cannot supply it is refused, not
   silently downgraded. *)
VerificationTest[
    Quiet[QECLogicalErrorRate[QECCode["DistanceTwo", 12], QECNoiseModel["Depolarizing", 1/50], 500,
        "Decoder" -> "MaximumLikelihood"]],
    $Failed,
    TestID -> "QEC-Decoder-explicit-ML-refused-above-the-limit"
]

(* On a code that can be enumerated, the option reproduces both decoders exactly.
   Under biased noise they differ and maximum likelihood wins: the leading
   coefficient drops from 756/25 to 651/25. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], biased, ml, mw, lead},
        biased = QECNoiseModel[<|"X" -> qecP/10, "Y" -> qecP/10, "Z" -> qecP|>];
        ml = QECLogicalErrorRate[code, biased, "Decoder" -> "MaximumLikelihood"];
        mw = QECLogicalErrorRate[code, biased, "Decoder" -> "MinimumWeight"];
        lead[e_] := Coefficient[Normal[Series[e, {qecP, 0, 2}]], qecP, 2];
        {lead[ml], lead[mw]}
    ],
    {651/25, 756/25},
    TestID -> "QEC-Decoder-ML-beats-minimum-weight-under-bias"
]

(* On a perfect code the table is already optimal: every weight-one error has its
   own syndrome and there is nothing for the coset weighing to improve, so the two
   decoders give the same leading coefficient under symmetric noise. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"], noise, lead},
        noise = QECNoiseModel["Depolarizing", qecP];
        lead[e_] := Coefficient[Normal[Series[e, {qecP, 0, 2}]], qecP, 2];
        {lead[QECLogicalErrorRate[code, noise, "Decoder" -> "MaximumLikelihood"]],
         lead[QECLogicalErrorRate[code, noise, "Decoder" -> "MinimumWeight"]]}
    ],
    {10, 10},
    TestID -> "QEC-Decoder-agree-on-a-perfect-code"
]

(* On a degenerate code they part company even under symmetric noise, and this is
   the textbook reason degeneracy matters: many distinct light errors share a coset,
   the table sees one representative per syndrome and can pick the wrong coset,
   while weighing whole cosets adds their probabilities up correctly. On Shor that
   is 13 p^2 against 31 p^2 -- the table is more than twice as bad. *)
VerificationTest[
    Module[{code = QECCode["ShorCode"], noise, lead},
        noise = QECNoiseModel["Depolarizing", qecP];
        lead[e_] := Coefficient[Normal[Series[e, {qecP, 0, 2}]], qecP, 2];
        {lead[QECLogicalErrorRate[code, noise, "Decoder" -> "MaximumLikelihood"]],
         lead[QECLogicalErrorRate[code, noise, "Decoder" -> "MinimumWeight"]]}
    ],
    {13, 31},
    TestID -> "QEC-Decoder-degeneracy-separates-them-on-shor"
]

(* Automatic on a small code is maximum likelihood, so nothing changed there. *)
VerificationTest[
    Module[{code = QECCode["5QubitCode"], noise = QECNoiseModel["Depolarizing", qecP]},
        Simplify[
            QECLogicalErrorRate[code, noise] -
            QECLogicalErrorRate[code, noise, "Decoder" -> "MaximumLikelihood"]
        ]
    ],
    0,
    TestID -> "QEC-Decoder-automatic-is-ML-when-affordable"
]

(* A deeper table is a deliberate choice with a cost, so the reach is an option. *)
VerificationTest[
    Module[{code = QECCode["SteaneCode"], noise = QECNoiseModel["Depolarizing", 1/50], shallow, deep},
        SeedRandom[22];
        shallow = QECLogicalErrorRate[code, noise, 20000, "Decoder" -> "MinimumWeight", "DecoderReach" -> 1];
        SeedRandom[22];
        deep = QECLogicalErrorRate[code, noise, 20000, "Decoder" -> "MinimumWeight", "DecoderReach" -> 2];
        {NumericQ[shallow], NumericQ[deep], deep <= shallow}
    ],
    {True, True, True},
    TestID -> "QEC-Decoder-reach-is-an-option"
]
