(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Bosonic.wlt

   The bosonic code object and the symbolic Knill-Laflamme decision.

   Every answer here is an exact symbolic identity, not a tolerance on a truncated
   matrix.  The binomial milestone of Bosonic-QEC-PartA-Spec is pinned first:
   h = {{1, 0}, {0, 2}} for the error set {I, a}.  The cats then exercise the
   coherent-state route, where the codeword norms carry the tanh/coth factors that
   decide correctability.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

{qecAv, qecAdv} = FieldVariables[];


BeginTestSection["QECBosonicCode - construction"]

(* The spec's worked example: N = 1, S = 1 gives (|0> + |4>)/Sqrt[2] and |2>. *)
VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["Codewords"],
    {<|0 -> 1/Sqrt[2], 4 -> 1/Sqrt[2]|>, <|2 -> 1|>},
    TestID -> "QBC-Binomial-11-Codewords"
]

VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["Parameters"],
    <|"N" -> 1, "S" -> 1|>,
    TestID -> "QBC-Binomial-Parameters"
]

VerificationTest[
    {QECBosonicCode["Binomial", 1, 1]["Basis"], QECBosonicCode["Cat", 4, 2]["Basis"]},
    {"Fock", "Coherent"},
    TestID -> "QBC-Basis-Detection"
]

VerificationTest[
    QECBosonicCode["Cat", 4, \[Alpha]]["Modes"],
    1,
    TestID -> "QBC-Modes"
]

(* Any even number of legs is a rotation-symmetric cat code; odd counts and
   degenerate ones are not. *)
VerificationTest[
    QECBosonicCode["Cat", 3, \[Alpha]],
    $Failed,
    {QECBosonicCode::legs},
    TestID -> "QBC-Cat-OddLegs"
]

VerificationTest[
    QECBosonicCode["Cat", 0, \[Alpha]],
    $Failed,
    {QECBosonicCode::legs},
    TestID -> "QBC-Cat-DegenerateLegs"
]

VerificationTest[
    QECBosonicCode["Cat", 6, \[Alpha]]["Parameters"],
    <|"Legs" -> 6, "Alpha" -> \[Alpha]|>,
    TestID -> "QBC-Cat6-Parameters"
]

VerificationTest[
    QECBosonicCode[{1, 2}],
    $Failed,
    {QECBosonicCode::words},
    TestID -> "QBC-BadCodewords"
]

VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["Nonsense"],
    $Failed,
    {QECBosonicCode::noprop},
    TestID -> "QBC-UndefinedProperty"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - mean photon number"]

(* Equal across the two codewords, which is half of why loss is correctable here. *)
VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["MeanPhotonNumber"],
    {2, 2},
    TestID -> "QBC-Binomial-MeanPhotonNumber"
]

(* The two-legged cat in closed form.  These differ, so the conditions fail. *)
VerificationTest[
    QECBosonicCode["Cat", 2, \[Alpha], Assumptions -> \[Alpha] > 0]["MeanPhotonNumber"],
    {\[Alpha]^2 Tanh[\[Alpha]^2], \[Alpha]^2 Coth[\[Alpha]^2]},
    TestID -> "QBC-Cat2-MeanPhotonNumber"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - Knill-Laflamme, binomial"]

(* The milestone: h00 = <W|I|W> = 1, h11 = <W|n|W> = 2, and h01 = 0 because a flips
   photon-number parity.  Exact, with no Fock cutoff anywhere. *)
VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["KnillLaflammeMatrix", "Loss"[1]],
    {{1, 0}, {0, 2}},
    TestID -> "QBC-Binomial-KLMatrix-Loss1"
]

(* Every block is a multiple of the identity over codewords, which is the condition
   itself rather than a summary of it. *)
VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["KnillLaflammeBlocks", "Loss"[1]],
    {{{{1, 0}, {0, 1}}, {{0, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}, {{2, 0}, {0, 2}}}},
    TestID -> "QBC-Binomial-KLBlocks-Loss1"
]

VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["CorrectableQ", "Loss"[1]],
    True,
    TestID -> "QBC-Binomial-Correctable-Loss1"
]

(* Two losses are not corrected: <W|ad^2 a^2|W> is 6 on one codeword and 2 on the
   other, so the Knill-Laflamme diagonal is codeword-dependent. *)
VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["CorrectableQ", "Loss"[2]],
    False,
    TestID -> "QBC-Binomial-NotCorrectable-Loss2"
]

VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["CorrectionOrder"],
    1,
    TestID -> "QBC-Binomial-CorrectionOrder"
]

(* An explicit set of ladder words must agree with the named channel. *)
VerificationTest[
    QECBosonicCode["Binomial", 1, 1]["KnillLaflammeMatrix", {1, qecAv}],
    QECBosonicCode["Binomial", 1, 1]["KnillLaflammeMatrix", "Loss"[1]],
    TestID -> "QBC-ExplicitErrorSet-MatchesNamed"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - Knill-Laflamme, cats"]

(* a sends the mod-4 combs out of the code space, so the single-loss block vanishes
   identically and only the mean-photon gap obstructs correction. *)
VerificationTest[
    QECBosonicCode["Cat", 4, \[Alpha], Assumptions -> \[Alpha] > 0][
        "KnillLaflammeBlocks", "Loss"[1]][[1, 2]],
    {{0, 0}, {0, 0}},
    TestID -> "QBC-Cat4-LossBlockVanishes"
]

(* Neither cat corrects a single loss exactly at finite alpha.  Cat4 is an approximate
   code, exact only as alpha -> infinity; Cat2 is not approximately correcting at all,
   since its off-diagonal block grows with alpha. *)
VerificationTest[
    QECBosonicCode["Cat", 4, \[Alpha], Assumptions -> \[Alpha] > 0]["CorrectableQ", "Loss"[1]],
    False,
    TestID -> "QBC-Cat4-NotExactlyCorrectable"
]

VerificationTest[
    QECBosonicCode["Cat", 2, \[Alpha], Assumptions -> \[Alpha] > 0]["CorrectableQ", "Loss"[1]],
    False,
    TestID -> "QBC-Cat2-NotExactlyCorrectable"
]

(* Cat2 fails in the off-diagonal, not the diagonal: a maps the even cat onto the odd
   cat, so the block is nonzero exactly where the binomial code's is zero. *)
VerificationTest[
    FullSimplify[
        QECBosonicCode["Cat", 2, \[Alpha], Assumptions -> \[Alpha] > 0][
            "KnillLaflammeBlocks", "Loss"[1]][[1, 2]] == {{0, 0}, {0, 0}},
        \[Alpha] > 0],
    False,
    TestID -> "QBC-Cat2-LossBlockNonzero"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - the 2d-leg family corrects d-1 losses"]

(* The defining property of the rotation-symmetric cat family: on 2d legs, a^k takes
   the code space outside itself for every k < d, and back inside at k = d, where it
   acts as a logical error the syndrome cannot see. Derived from the codeword
   overlaps rather than assumed, so it is a check on the layer as much as on the
   codes. See Grimsmo, Combes & Baragiola on rotation-symmetric bosonic codes. *)

ClearAll[lossBlock, lossBlockZeroQ]
lossBlock[legs_, k_] :=
    QECBosonicCode["Cat", legs, \[Alpha], Assumptions -> \[Alpha] > 0][
        "KnillLaflammeBlocks", {1, Nest[# ** \[FormalA] &, \[FormalA], k - 1]}][[1, 2]]
lossBlockZeroQ[legs_, k_] :=
    TrueQ @ Simplify[lossBlock[legs, k] == {{0, 0}, {0, 0}}, \[Alpha] > 0]

(* d = 1: even/odd cat, a single loss already acts inside the code space *)
VerificationTest[
    lossBlockZeroQ[2, 1],
    False,
    TestID -> "QBC-Cat2-d1-FailsAtOne"
]

(* d = 2: one loss leaves, two land back inside *)
VerificationTest[
    {lossBlockZeroQ[4, 1], lossBlockZeroQ[4, 2]},
    {True, False},
    TestID -> "QBC-Cat4-d2-FailsAtTwo"
]

(* d = 3: the six-component cat, robust to two losses where the four-component one
   is not, which is the reason the family is generalised past four legs *)
VerificationTest[
    {lossBlockZeroQ[6, 1], lossBlockZeroQ[6, 2], lossBlockZeroQ[6, 3]},
    {True, True, False},
    TestID -> "QBC-Cat6-d3-FailsAtThree"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - approximate correction"]

(* An exact code is correctable at finite amplitude, so the two orders agree. *)
VerificationTest[
    {QECBosonicCode["Binomial", 2, 2]["CorrectionOrder"],
     QECBosonicCode["Binomial", 2, 2]["ApproximateCorrectionOrder", "Loss", \[Alpha]]},
    {2, 2},
    TestID -> "QBC-Binomial-ExactEqualsApproximate"
]

(* A cat satisfies the conditions only as alpha grows, so its exact order is 0 while its
   approximate order is d - 1.  This is the distinction the boolean cannot make: Cat2 is
   not approximately correcting at all, Cat4 is. *)
VerificationTest[
    {QECBosonicCode["Cat", 2, \[Alpha], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder"],
     QECBosonicCode["Cat", 4, \[Alpha], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder"]},
    {0, 1},
    TestID -> "QBC-Cat-ApproximateOrder"
]

(* The residual is the obstruction itself, and vanishes identically for an exact code. *)
VerificationTest[
    DeleteDuplicates @ QECBosonicCode["Binomial", 1, 1]["KnillLaflammeResidual", "Loss"[1]],
    {0},
    TestID -> "QBC-Residual-ZeroForExactCode"
]

(* A code with no amplitude parameter has no limit to take. *)
VerificationTest[
    QECBosonicCode[{<|0 -> 1|>, <|1 -> 1|>}]["ApproximateCorrectionOrder"],
    $Failed,
    {QECBosonicCode::novar},
    TestID -> "QBC-ApproximateOrder-NoAmplitude"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - truncated conversion"]

(* A Fock codeword is a finite sum, so its Fock space size is exact: highest occupied
   level plus one.  A coherent codeword is sized from the Poisson tail of its amplitude. *)
VerificationTest[
    {QECBosonicCode["Binomial", 1, 1]["FockSpaceSize"],
     QECBosonicCode["Binomial", 3, 2]["FockSpaceSize"]},
    {5, 13},
    TestID -> "QBC-FockSpaceSize-Fock-Exact"
]

VerificationTest[
    QECBosonicCode["Cat", 4, \[Alpha]]["CodewordStates"],
    $Failed,
    {QECBosonicCode::numeric},
    TestID -> "QBC-FockSpaceSize-NeedsNumeric"
]

(* Auto-sized states are normalized and the two codewords stay orthogonal. *)
VerificationTest[
    With[{qs = QECBosonicCode["Cat", 4, 2.]["CodewordStates"]},
        Chop[{#["Norm"] & /@ qs, (First[qs]["Dagger"] @ Last[qs])["Scalar"]} - {{1, 1}, 0}]],
    {{0, 0}, 0},
    TestID -> "QBC-CodewordStates-Auto-Orthonormal"
]

(* The one cutoff in the object, for handing codewords to the phase-space tools. *)
VerificationTest[
    #["Dimensions"] & /@ QECBosonicCode["Binomial", 1, 1]["CodewordStates", 12],
    {{12}, {12}},
    TestID -> "QBC-CodewordStates-Dimensions"
]

VerificationTest[
    Chop[#["Norm"] - 1] & /@ QECBosonicCode["Cat", 4, 2.]["CodewordStates", 24],
    {0, 0},
    TestID -> "QBC-CodewordStates-Normalized"
]

EndTestSection[]
