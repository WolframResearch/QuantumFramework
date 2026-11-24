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

{{qecAdv}, {qecAv}} = FieldVariables[];


BeginTestSection["QECBosonicCode - construction"]

(* The spec's worked example: N = 1, S = 1 gives (|0> + |4>)/Sqrt[2] and |2>. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["Codewords"],
    {<|0 -> 1/Sqrt[2], 4 -> 1/Sqrt[2]|>, <|2 -> 1|>},
    TestID -> "QBC-Binomial-11-Codewords"
]

VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["Parameters"],
    <|"N" -> 1, "S" -> 1|>,
    TestID -> "QBC-Binomial-Parameters"
]

VerificationTest[
    {QECBosonicCode["Binomial"[1, 1]]["Basis"], QECBosonicCode["Cat"[4, 2]]["Basis"]},
    {"Fock", "Coherent"},
    TestID -> "QBC-Basis-Detection"
]

VerificationTest[
    QECBosonicCode["Cat"[4, \[Alpha]]]["Modes"],
    1,
    TestID -> "QBC-Modes"
]

(* Any even number of legs is a rotation-symmetric cat code; odd counts and
   degenerate ones are not. *)
VerificationTest[
    QECBosonicCode["Cat"[3, \[Alpha]]],
    $Failed,
    {QECBosonicCode::legs},
    TestID -> "QBC-Cat-OddLegs"
]

VerificationTest[
    QECBosonicCode["Cat"[0, \[Alpha]]],
    $Failed,
    {QECBosonicCode::legs},
    TestID -> "QBC-Cat-DegenerateLegs"
]

VerificationTest[
    QECBosonicCode["Cat"[6, \[Alpha]]]["Parameters"],
    <|"Legs" -> 6, "Alpha" -> \[Alpha]|>,
    TestID -> "QBC-Cat6-Parameters"
]

VerificationTest[
    QECBosonicCode[{1, 2}],
    $Failed,
    {QECBosonicCode::words},
    TestID -> "QBC-BadCodewords"
]

(* An unknown property, and a channel property asked without its channel. *)
VerificationTest[
    {QECBosonicCode["Binomial"[1, 1]]["Nonsense"], QECBosonicCode["Binomial"[1, 1]]["CorrectableQ"]},
    {$Failed, $Failed},
    {QECBosonicCode::noprop, QECBosonicCode::nochannel},
    TestID -> "QBC-UndefinedProperty"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - mean photon number"]

(* Equal across the two codewords, which is half of why loss is correctable here. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["MeanPhotonNumber"],
    {2, 2},
    TestID -> "QBC-Binomial-MeanPhotonNumber"
]

(* The two-legged cat in closed form.  These differ, so the conditions fail. *)
VerificationTest[
    QECBosonicCode["Cat"[2, \[Alpha]], Assumptions -> \[Alpha] > 0]["MeanPhotonNumber"],
    {\[Alpha]^2 Tanh[\[Alpha]^2], \[Alpha]^2 Coth[\[Alpha]^2]},
    TestID -> "QBC-Cat2-MeanPhotonNumber"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - Knill-Laflamme, binomial"]

(* The milestone: h00 = <W|I|W> = 1, h11 = <W|n|W> = 2, and h01 = 0 because a flips
   photon-number parity.  Exact, with no Fock cutoff anywhere.  "Loss"[L] is the idealized
   span {I, a, ..., a^L}, not the Kraus set of the physical channel, whose no-jump factor
   eta^(n/2) is not a polynomial in the ladder operators. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeMatrix", "Loss"[1]],
    {{1, 0}, {0, 2}},
    TestID -> "QBC-Binomial-KLMatrix-Loss1"
]

(* Every block is a multiple of the identity over codewords, which is the condition
   itself rather than a summary of it. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeBlocks", "Loss"[1]],
    {{{{1, 0}, {0, 1}}, {{0, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}, {{2, 0}, {0, 2}}}},
    TestID -> "QBC-Binomial-KLBlocks-Loss1"
]

VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["CorrectableQ", "Loss"[1]],
    True,
    TestID -> "QBC-Binomial-Correctable-Loss1"
]

(* Two losses are not corrected: <W|ad^2 a^2|W> is 6 on one codeword and 2 on the
   other, so the Knill-Laflamme diagonal is codeword-dependent. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["CorrectableQ", "Loss"[2]],
    False,
    TestID -> "QBC-Binomial-NotCorrectable-Loss2"
]

VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["CorrectionOrder"],
    1,
    TestID -> "QBC-Binomial-CorrectionOrder"
]

(* An explicit set of ladder words must agree with the named channel. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeMatrix", {1, qecAv}],
    QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeMatrix", "Loss"[1]],
    TestID -> "QBC-ExplicitErrorSet-MatchesNamed"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - Knill-Laflamme, cats"]

(* a sends the mod-4 combs out of the code space, so the single-loss block vanishes
   identically and only the mean-photon gap obstructs correction. *)
VerificationTest[
    QECBosonicCode["Cat"[4, \[Alpha]], Assumptions -> \[Alpha] > 0][
        "KnillLaflammeBlocks", "Loss"[1]][[1, 2]],
    {{0, 0}, {0, 0}},
    TestID -> "QBC-Cat4-LossBlockVanishes"
]

(* Neither cat corrects a single loss exactly at finite alpha.  Cat4 is an approximate
   code, exact only as alpha -> infinity; Cat2 is not approximately correcting at all,
   since its off-diagonal block grows with alpha. *)
VerificationTest[
    QECBosonicCode["Cat"[#, \[Alpha]], Assumptions -> \[Alpha] > 0]["CorrectableQ", "Loss"[1]] & /@ {2, 4},
    {False, False},
    TestID -> "QBC-Cat-NotExactlyCorrectable"
]

(* Cat2 fails in the off-diagonal, not the diagonal: a maps the even cat onto the odd
   cat, so the block is nonzero exactly where the binomial code's is zero. *)
VerificationTest[
    FullSimplify[
        QECBosonicCode["Cat"[2, \[Alpha]], Assumptions -> \[Alpha] > 0][
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
    QECBosonicCode["Cat"[legs, \[Alpha]], Assumptions -> \[Alpha] > 0][
        "KnillLaflammeBlocks", {1, Nest[# ** \[FormalA] &, \[FormalA], k - 1]}][[1, 2]]
lossBlockZeroQ[legs_, k_] :=
    TrueQ @ Simplify[lossBlock[legs, k] == {{0, 0}, {0, 0}}, \[Alpha] > 0]

(* d = 1: even/odd cat, a single loss already acts inside the code space.
   d = 2: one loss leaves, two land back inside.
   d = 3: the six-component cat, robust to two losses where the four-component one
   is not, which is the reason the family is generalised past four legs. *)
VerificationTest[
    Table[lossBlockZeroQ[2 d, k], {d, 3}, {k, d}],
    {{False}, {True, False}, {True, True, False}},
    TestID -> "QBC-Cat-2dLegs-FailAtD"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - approximate correction"]

(* An exact code is correctable at finite amplitude, so the two orders agree. *)
VerificationTest[
    {QECBosonicCode["Binomial"[2, 2]]["CorrectionOrder"],
     QECBosonicCode["Binomial"[2, 2]]["ApproximateCorrectionOrder", "Loss", \[Alpha]]},
    {2, 2},
    TestID -> "QBC-Binomial-ExactEqualsApproximate"
]

(* A coherent code whose amplitudes are the limit variable times fixed numbers, beta_k = a z_k,
   decides its large-amplitude limit without Limit: |<a z_k|a z_l>| = Exp[-a^2 |z_k - z_l|^2/2]
   vanishes for distinct z, so only each component paired with itself survives and a residual
   tends to zero exactly when its polynomial in a vanishes identically.  The cat orders below
   were also obtained with Limit, including the 4-cat under gain, which took Limit 109 s. *)
VerificationTest[
    Quiet[QECBosonicCode["Cat"[#1, \[Alpha]], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", #2] & @@@
        {{2, "Loss"}, {2, "Gain"}, {2, "Dephasing"}, {4, "Loss"}, {4, "Gain"}, {4, "Dephasing"}}],
    {0, 0, 1, 1, 1, 2},
    TestID -> "QBC-CoherentLimit-Cats"
]

(* No symmetry needed.  Code A's four amplitudes are distinct and both codewords give
   1 + 2^d for <n^d>, so dephasing holds up to the scan's end.  Code B fails at order one:
   <n> is 1 on one codeword and 5/2 on the other, and the (a, a) entries are a^2 and 5a^2/2. *)
VerificationTest[
    Quiet[{
        QECBosonicCode[{{{1, \[Alpha]}, {1, (1 + I) \[Alpha]}}, {{1, -\[Alpha]}, {1, -(1 + I) \[Alpha]}}},
            Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", #, \[Alpha]] & /@ {"Loss", "Dephasing"},
        QECBosonicCode[{{{1, \[Alpha]}, {2, -\[Alpha]}, {1, I \[Alpha]}}, {{1, -I \[Alpha]}, {1, 2 \[Alpha]}}},
            Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", #, \[Alpha]] & /@ {"Loss", "Dephasing"}}],
    {{0, 2}, {0, 0}},
    TestID -> "QBC-CoherentLimit-Asymmetric"
]

(* A leg is a value, not a form.  The rotated 4-cat, its last leg written E^(I Pi/5) or in
   radicals, which Simplify keeps apart.  Split in halves within one codeword, the two forms
   read as orthogonal legs break the norms; written differently in the two codewords, the leg
   drops out of <W0|W1>.  Either way order 0 fails and gives -1. *)
VerificationTest[
    With[{legs = \[Alpha] E^(I Pi/5) I^Range[0, 3],
          radical = -I \[Alpha] ((1 + Sqrt[5])/4 + I Sqrt[(5 - Sqrt[5])/8]),
          signs = {1, -1, 1, -1}},
        QECBosonicCode[#, Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", "Loss", \[Alpha]] & /@ {
            {{1, #} & /@ legs,
             Join[Transpose[{Most[signs], Most[legs]}], {{-1/2, Last[legs]}, {-1/2, radical}}]},
            {{1, #} & /@ legs,
             Transpose[{signs, Append[Most[legs], radical]}]}}
    ],
    {1, 1},
    TestID -> "QBC-CoherentLimit-EqualLegsInTwoForms"
]

(* Amplitudes that are not a multiple of the limit variable still go through Limit. *)
VerificationTest[
    QECBosonicCode[{{{1, \[Alpha] + 1}, {1, -\[Alpha] - 1}}, {{1, \[Alpha] + 1}, {-1, -\[Alpha] - 1}}},
        Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", "Loss", \[Alpha]],
    0,
    TestID -> "QBC-ApproximateOrder-LimitFallback"
]

(* 2d legs correct d - 1 losses in the limit; with Limit the 6-cat alone took minutes. *)
VerificationTest[
    Table[QECBosonicCode["Cat"[l, \[Alpha]], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder"], {l, {6, 8, 10}}],
    {2, 3, 4},
    TestID -> "QBC-CatLimit-LargeCats"
]

(* The residual is the obstruction itself, and vanishes identically for an exact code. *)
VerificationTest[
    DeleteDuplicates @ Values @ QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeResidual", "Loss"[1]],
    {{{0, 0}, {0, 0}}},
    TestID -> "QBC-Residual-ZeroForExactCode"
]

(* Keyed by error pair: past its design order the binomial code fails only where a^2 enters. *)
VerificationTest[
    Keys @ Select[QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeResidual", "Loss"[2]],
        # =!= {{0, 0}, {0, 0}} &],
    {{1, qecAv ** qecAv}, {qecAv ** qecAv, 1}, {qecAv ** qecAv, qecAv ** qecAv}},
    TestID -> "QBC-Residual-KeyedByErrorPair"
]

(* Matrix elements are cached once per code, error product and codeword pair, so repeated
   properties add nothing and the next loss order adds only its new products: {I, a} gives
   4 products x 4 codeword pairs = 16 entries, and {I, a, a^2} adds 5 more products, 20 entries. *)
VerificationTest[
    With[{element = Symbol["Wolfram`QuantumFramework`QEC`PackageScope`bosonicElement"],
          code = QECBosonicCode["Binomial"[2, 1]]},
        With[{n0 = Length @ DownValues[element]},
            code["CorrectableQ", "Loss"[1]];
            code["KnillLaflammeMatrix", "Loss"[1]];
            code["KnillLaflammeResidual", "Loss"[1]];
            With[{n1 = Length @ DownValues[element]},
                code["CorrectableQ", "Loss"[2]];
                {n1 - n0, Length @ DownValues[element] - n1}]]],
    {16, 20},
    TestID -> "QBC-Blocks-Memoised"
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
    {QECBosonicCode["Binomial"[1, 1]]["FockSpaceSize"],
     QECBosonicCode["Binomial"[3, 2]]["FockSpaceSize"]},
    {5, 13},
    TestID -> "QBC-FockSpaceSize-Fock-Exact"
]

VerificationTest[
    QECBosonicCode["Cat"[4, \[Alpha]]]["CodewordStates"],
    $Failed,
    {QECBosonicCode::numeric},
    TestID -> "QBC-FockSpaceSize-NeedsNumeric"
]

(* Auto-sized states are normalized and the two codewords stay orthogonal. *)
VerificationTest[
    With[{qs = QECBosonicCode["Cat"[4, 2.]]["CodewordStates"]},
        Chop[{#["Norm"] & /@ qs, (First[qs]["Dagger"] @ Last[qs])["Scalar"]} - {{1, 1}, 0}]],
    {{0, 0}, 0},
    TestID -> "QBC-CodewordStates-Auto-Orthonormal"
]

(* The one cutoff in the object, for handing codewords to the phase-space tools. *)
VerificationTest[
    {#["Dimensions"] & /@ QECBosonicCode["Binomial"[1, 1]]["CodewordStates", 12],
     Chop[#["Norm"] - 1] & /@ QECBosonicCode["Cat"[4, 2.]]["CodewordStates", 24]},
    {{{12}, {12}}, {0, 0}},
    TestID -> "QBC-CodewordStates-ExplicitSize"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - displacement error sets"]

(* Every normally ordered monomial up to total degree L, so the set has the triangular
   count (L+1)(L+2)/2 rather than the L+1 of a loss set. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["ErrorSet", "Displacement"[1]],
    {1, qecAv, qecAdv},
    TestID -> "QBC-Displacement-ErrorSet"
]

VerificationTest[
    Table[Length[QECBosonicCode["Binomial"[1, 1]]["ErrorSet", "Displacement"[l]]], {l, 0, 3}],
    Table[(l + 1) (l + 2)/2, {l, 0, 3}],
    TestID -> "QBC-Displacement-SetSize"
]

(* A code built against photon loss does not correct displacements: a^dagger takes it out
   of the code space in a way the loss set never probes. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["CorrectableQ", "Displacement"[1]],
    False,
    TestID -> "QBC-Displacement-BinomialNotCorrectable"
]

(* The structural difference from loss.  A loss set is graded by photon number, so
   E_p^dagger E_q moves it unless p = q and the coefficient matrix comes out diagonal at
   every order.  A displacement set is not graded - pairs such as (I, a^dagger a) preserve
   photon number - so the off-diagonal coefficients survive.  They are what a signature
   in the sense of Du et al. would be read from. *)
VerificationTest[
    With[{offDiagonal = DeleteCases[Flatten[# - DiagonalMatrix[Diagonal[#]]], 0] &},
        {offDiagonal[QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeMatrix", "Loss"[2]]],
         offDiagonal[QECBosonicCode["Binomial"[1, 1]]["KnillLaflammeMatrix", "Displacement"[2]]]}],
    {{}, {2, Sqrt[3/2], 2, Sqrt[3/2]}},
    TestID -> "QBC-Displacement-OffDiagonalsSurvive"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - gain, dephasing and combined channels"]

(* Michael et al., PRX 6, 031006 (2016): the error set of Eq. (6) is the union
   {I, a..a^L, ad..ad^G, n..n^D}, and the code of Eq. (7) corrects it when S = L + G
   and N = max{L, G, 2D}.  Their Eqs. (2), (3) and (5) are the (N, S) = (1, 1), (2, 2)
   and (2, 1) codes. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["ErrorSet", "Gain"[2]],
    {1, qecAdv, qecAdv ** qecAdv},
    TestID -> "QBC-Gain-ErrorSet"
]

VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["ErrorSet", {"Loss"[1], "Dephasing"[1]}],
    {1, qecAv, qecAdv ** qecAv},
    TestID -> "QBC-Combined-ErrorSet-Union"
]

VerificationTest[
    QECBosonicCode["Binomial"[2, 1]]["Codewords"],
    {<|0 -> 1/2, 4 -> Sqrt[3]/2|>, <|2 -> Sqrt[3]/2, 6 -> 1/2|>},
    TestID -> "QBC-Binomial-21-Eq5-Codewords"
]

(* Eq. (5) corrects {I, a, n}.  Dephasing makes h non-diagonal: h(I, n) = <n> = 3 and
   h(n, n) = <n^2> = 12 on both codewords, while a flips parity and leaves h(I, a) = 0. *)
VerificationTest[
    QECBosonicCode["Binomial"[2, 1]]["KnillLaflammeMatrix", {"Loss"[1], "Dephasing"[1]}],
    {{1, 0, 3}, {0, 3, 0}, {3, 0, 12}},
    TestID -> "QBC-Binomial-21-KLMatrix-LossDephasing"
]

VerificationTest[
    QECBosonicCode["Binomial"[2, 1]]["CorrectableQ", {"Loss"[1], "Dephasing"[1]}],
    True,
    TestID -> "QBC-Binomial-21-Correctable-LossDephasing"
]

(* N = 2 >= 2D stops at one dephasing. *)
VerificationTest[
    QECBosonicCode["Binomial"[2, 1]]["CorrectableQ", "Dephasing"[2]],
    False,
    TestID -> "QBC-Binomial-21-NotCorrectable-Dephasing2"
]

(* Eq. (2) has N = 1, so no dephasing: <n^2> is 8 on one codeword and 4 on the other. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1]]["CorrectableQ", "Dephasing"[1]],
    False,
    TestID -> "QBC-Binomial-11-NotCorrectable-Dephasing1"
]

(* S = L + G = 1 buys one loss or one gain, not both: a W0 and ad W1 both land on |3>. *)
VerificationTest[
    {QECBosonicCode["Binomial"[1, 1]]["CorrectableQ", "Gain"[1]],
     QECBosonicCode["Binomial"[1, 1]]["CorrectableQ", {"Loss"[1], "Gain"[1]}]},
    {True, False},
    TestID -> "QBC-Binomial-11-GainOrLossNotBoth"
]

(* Eq. (3) corrects either {I, a, a^2, n} or {I, a, ad, n}. *)
VerificationTest[
    {QECBosonicCode["Binomial"[2, 2]]["CorrectableQ", {"Loss"[2], "Dephasing"[1]}],
     QECBosonicCode["Binomial"[2, 2]]["CorrectableQ", {"Loss"[1], "Gain"[1], "Dephasing"[1]}]},
    {True, True},
    TestID -> "QBC-Binomial-22-Eq3-Correctable"
]

(* {f1, f2, ...}[K] is every sequence of up to K jumps drawn from the families, in any order:
   the errors at order dt^K when the channels act together.  A loss then a dephasing is as
   likely as two losses, so the union {"Loss"[K], "Dephasing"[K]} is not the same set.  Words
   that add nothing to the span are dropped (a a\[Dagger] and a\[Dagger] a differ by I), so with
   all three families K jumps leave (K + 1)^2 independent errors. *)
VerificationTest[
    With[{code = QECBosonicCode["Binomial"[1, 1]]},
        {code["ErrorSet", {"Loss", "Gain"}[1]],
         code["ErrorSet", {"Loss"}[2]] === code["ErrorSet", "Loss"[2]],
         Length[code["ErrorSet", {"Loss", "Gain", "Dephasing"}[#]]] & /@ Range[0, 3]}],
    {{1, qecAv, qecAdv}, True, {1, 4, 9, 16}},
    TestID -> "QBC-Composed-ErrorSet"
]

(* Loss and gain together span the normally ordered monomials, which is "Displacement". *)
VerificationTest[
    {#["CorrectionOrder", {"Loss", "Gain"}], #["CorrectionOrder", "Displacement"]} & /@
        {QECBosonicCode["Binomial"[1, 1]], QECBosonicCode["Binomial"[2, 2]]},
    {{0, 0}, {1, 1}},
    TestID -> "QBC-Composed-LossGainIsDisplacement"
]

(* (N, S) = (2, 1) corrects each channel once, but a loss and a gain land on the same occupations,
   so any family set holding both drops to order 0. *)
VerificationTest[
    With[{code = QECBosonicCode["Binomial"[2, 1]]},
        code["CorrectionOrder", #] & /@
            {"Loss", "Gain", "Dephasing", {"Loss", "Gain"}, {"Loss", "Dephasing"}, {"Gain", "Dephasing"}, {"Loss", "Gain", "Dephasing"}}],
    {1, 1, 1, 0, 1, 1, 0},
    TestID -> "QBC-Composed-CorrectionOrder"
]

VerificationTest[
    QECBosonicCode["Cat"[#, \[Alpha]], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", #2] & @@@
        {{4, {"Loss"}}, {4, {"Loss", "Gain"}}, {4, {"Loss", "Dephasing"}}, {6, {"Loss", "Gain"}}, {6, {"Loss", "Dephasing"}}},
    {1, 0, 1, 1, 2},
    TestID -> "QBC-Composed-ApproximateOrder"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - ladder variables"]

VerificationTest[
    QECBosonicCode["Binomial"[2, 1], qecB]["ErrorSet", {"Loss"[1], "Dephasing"[1]}],
    {1, qecB, SuperDagger[qecB] ** qecB},
    TestID -> "QBC-Vars-NamedChannelInVars"
]

VerificationTest[
    QECBosonicCode["Binomial"[2, 1], {{SuperDagger[qecB]}, {qecB}}]["ErrorSet", "Loss"[1]],
    {1, qecB},
    TestID -> "QBC-Vars-FieldVariablesForm"
]

VerificationTest[
    QECBosonicCode["Binomial"[2, 1], qecB]["KnillLaflammeMatrix", {1, qecB, SuperDagger[qecB] ** qecB}],
    QECBosonicCode["Binomial"[2, 1]]["KnillLaflammeMatrix", {"Loss"[1], "Dephasing"[1]}],
    TestID -> "QBC-Vars-ExplicitErrorsMatchDefault"
]

VerificationTest[
    QECBosonicCode["Cat"[2, \[Alpha]], qecB, Assumptions -> \[Alpha] > 0]["MeanPhotonNumber"],
    {\[Alpha]^2 Tanh[\[Alpha]^2], \[Alpha]^2 Coth[\[Alpha]^2]},
    TestID -> "QBC-Vars-WithAssumptions"
]

VerificationTest[
    QECBosonicCode["Binomial"[2, 1], {qecB, qecC}],
    $Failed,
    {QECBosonicCode::vars},
    TestID -> "QBC-Vars-TwoModesRejected"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - QuantumState codewords"]

VerificationTest[
    QECBosonicCode[QECBosonicCode["Binomial"[2, 1]]["CodewordStates"]]["Codewords"],
    QECBosonicCode["Binomial"[2, 1]]["Codewords"],
    TestID -> "QBC-States-RoundTrip"
]

VerificationTest[
    QECBosonicCode[{FockState[0, 5], FockState[3, 5]}, qecB]["ErrorSet", "Loss"[1]],
    {1, qecB},
    TestID -> "QBC-States-WithVars"
]

VerificationTest[
    QECBosonicCode[{QuantumState[IdentityMatrix[3]/3, 3], FockState[1, 3]}],
    $Failed,
    {QECBosonicCode::states},
    TestID -> "QBC-States-MixedRejected"
]

VerificationTest[
    QECBosonicCode[{QuantumTensorProduct[FockState[0, 3], FockState[0, 3]],
        QuantumTensorProduct[FockState[1, 3], FockState[1, 3]]}],
    $Failed,
    {QECBosonicCode::states},
    TestID -> "QBC-States-TwoModesRejected"
]

(* A state written in another basis is read through its Fock amplitudes: (|+> + |->)/Sqrt[2] is |0>. *)
VerificationTest[
    QECBosonicCode[{QuantumState[{1, 1}/Sqrt[2], QuantumBasis["X"]],
        QuantumState[{1, -1}/Sqrt[2], QuantumBasis["X"]]}]["Codewords"],
    {<|0 -> 1|>, <|1 -> 1|>},
    TestID -> "QBC-States-NonFockBasis"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - explicit errors and edge cases"]

(* The adjoint swaps the code's two ladder operators and conjugates every scalar, so h is
   Hermitian: (2b)^dagger (i b) = 2i n averages to 4i and (i b)^dagger (2b) to -4i, while
   (b + bd)^2 averages to 2<n> + 1 = 5 because b^2 and bd^2 leave the codewords' parity class. *)
VerificationTest[
    QECBosonicCode["Binomial"[1, 1], qecB][
        "KnillLaflammeMatrix", {1, 2 qecB, I qecB, qecB + SuperDagger[qecB]}],
    {{1, 0, 0, 0}, {0, 8, 4 I, 4}, {0, -4 I, 2, -2 I}, {0, 4, 2 I, 5}},
    TestID -> "QBC-Explicit-ScalarsAndSums"
]

VerificationTest[
    QECBosonicCode["Binomial"[1, 1], qecB]["KnillLaflammeMatrix",
        {1, GeneralizedPower[NonCommutativeMultiply, qecB, 2], GeneralizedPower[NonCommutativeMultiply, SuperDagger[qecB], 2]}],
    QECBosonicCode["Binomial"[1, 1], qecB]["KnillLaflammeMatrix",
        {1, qecB ** qecB, SuperDagger[qecB] ** SuperDagger[qecB]}],
    TestID -> "QBC-Explicit-GeneralizedPower"
]

VerificationTest[
    QECBosonicCode["Binomial"[1, 1], qecB, Assumptions -> \[Epsilon] > 0][
        "KnillLaflammeMatrix", {1, \[Epsilon] qecB}],
    {{1, 0}, {0, 2 \[Epsilon]^2}},
    TestID -> "QBC-Explicit-SymbolicScalar"
]

(* The order scan runs until the conditions fail, not to a fixed cap: (N, S) = (5, 5)
   corrects min(N, S) = 5 losses. *)
VerificationTest[
    QECBosonicCode["Binomial"[5, 5]]["CorrectionOrder"],
    5,
    TestID -> "QBC-CorrectionOrder-NoFixedCap"
]

(* A bad channel is reported where the error set is built, and every property that needs the
   set fails with it: an empty set would pass the conditions vacuously.  Unknown or misspelled
   names, orders that are not non-negative integers, a bad member of a list, a channel name
   among ladder-operator errors, and a bad family in a combined channel. *)
VerificationTest[
    With[{code = QECBosonicCode["Binomial"[1, 1]]},
        Quiet[
            Outer[code[#1, #2] &,
                {"ErrorSet", "CorrectableQ", "KnillLaflammeBlocks", "KnillLaflammeMatrix",
                 "KnillLaflammeResidual", "KnillLaflammeRank", "Signature"},
                {"Foo"[1], "loss"[1], "Loss"[1.], "Loss"[-1], {"Loss"[1], "Loss"[-1]},
                 {"Loss"[1], qecAv ** qecAv}, {"Loss", "loss"}[1]},
                1],
            {QECBosonicCode::channel, General::stop}]
    ],
    ConstantArray[$Failed, {7, 7}],
    TestID -> "QBC-Channel-Rejected"
]

(* The order scans check the family before scanning, and report it. *)
VerificationTest[
    {QECBosonicCode["Binomial"[1, 1]]["CorrectionOrder", "loss"],
     QECBosonicCode["Binomial"[1, 1]]["CorrectionOrder", {"Loss", "loss"}],
     QECBosonicCode["Cat"[4, \[Alpha]], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", "loss"]},
    {$Failed, $Failed, $Failed},
    {QECBosonicCode::channel, QECBosonicCode::channel, QECBosonicCode::channel, General::stop},
    TestID -> "QBC-Channel-RejectedByOrders"
]

(* The limit variable can be any symbolic amplitude, not only a bare symbol; a number or a
   compound expression such as 2 \[Alpha] has no limit to take. *)
VerificationTest[
    {QECBosonicCode["Cat"[4, 2.]]["ApproximateCorrectionOrder"],
     QECBosonicCode["Cat"[4, Subscript[\[Alpha], 0]], Assumptions -> Subscript[\[Alpha], 0] > 0]["ApproximateCorrectionOrder"],
     QECBosonicCode["Cat"[4, \[Alpha]], Assumptions -> \[Alpha] > 0]["ApproximateCorrectionOrder", "Loss", 2 \[Alpha]]},
    {$Failed, 1, $Failed},
    {QECBosonicCode::numvar, QECBosonicCode::numvar},
    TestID -> "QBC-ApproximateOrder-AmplitudeVariable"
]

EndTestSection[]


BeginTestSection["QECBosonicCode - rank and signature"]

(* The rank of h counts the error images the code can tell apart. Eq. (5) of Michael et al.
   separates all three of {I, a, n}; an error set holding b and 2b is degenerate by construction. *)
VerificationTest[
    {QECBosonicCode["Binomial"[2, 1]]["KnillLaflammeRank", {"Loss"[1], "Dephasing"[1]}],
     QECBosonicCode["Binomial"[2, 1], qecB]["KnillLaflammeRank", {1, qecB, 2 qecB}]},
    {3, 2},
    TestID -> "QBC-Rank-FullAndDegenerate"
]

(* The signature of Du et al. (npj QI 2025) is the vector of off-diagonal coefficients of h.
   Bosonic errors carry no fixed normalization, so each is scaled to unit norm on the code:
   h(I, n)/Sqrt[h(I, I) h(n, n)] = 3/Sqrt[12] for Eq. (5), and a flips parity so the rest vanish. *)
VerificationTest[
    QECBosonicCode["Binomial"[2, 1]]["Signature", {"Loss"[1], "Dephasing"[1]}],
    <|{1, qecAv} -> 0, {1, qecAdv ** qecAv} -> Sqrt[3]/2, {qecAv, qecAdv ** qecAv} -> 0|>,
    TestID -> "QBC-Signature-Eq5"
]

VerificationTest[
    Values @ QECBosonicCode["Binomial"[2, 1], qecB]["Signature", {1, 2 qecB, 5 SuperDagger[qecB] ** qecB}],
    Values @ QECBosonicCode["Binomial"[2, 1], qecB]["Signature", {1, qecB, SuperDagger[qecB] ** qecB}],
    TestID -> "QBC-Signature-ScaleFree"
]

(* A loss set is graded by photon number, so h is diagonal and the signature vanishes. *)
VerificationTest[
    Union @ Values @ QECBosonicCode["Binomial"[2, 2]]["Signature", "Loss"[2]],
    {0},
    TestID -> "QBC-Signature-LossVanishes"
]

EndTestSection[]
