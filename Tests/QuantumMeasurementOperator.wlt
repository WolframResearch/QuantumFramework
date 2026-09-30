BeginTestSection["QuantumMeasurementOperator - constructors"]

VerificationTest[QuantumMeasurementOperator["RandomHermitian"[2]]["Dimensions"], {2, 2}, TestID -> "RandomHermitian-2"]

VerificationTest[QuantumMeasurementOperator["GellMannMICPOVM"[2]]["Dimensions"], {4, 2, 2}, TestID -> "GellMannMICPOVM-2"]

VerificationTest[QuantumMeasurementOperator["TetrahedronSICPOVM"]["Dimensions"], {4, 2, 2}, TestID -> "TetrahedronSICPOVM-bare"]

VerificationTest[QuantumMeasurementOperator["QBismSICPOVM"[2]]["Dimensions"], {4, 2, 2}, TestID -> "QBismSICPOVM-2"]

VerificationTest[QuantumMeasurementOperator["HesseSICPOVM"]["Dimensions"], {9, 3, 3}, TestID -> "HesseSICPOVM-bare"]

VerificationTest[QuantumMeasurementOperator["HoggarSICPOVM"]["Dimensions"], {64, 8, 8}, TestID -> "HoggarSICPOVM-bare"]

VerificationTest[QuantumMeasurementOperator[1]["Targets"], {{1}}, TestID -> "Integer-target"]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - Lueders measurement"]

(* "Lueders"[op] measures the observable op with one outcome per distinct eigenvalue, whose
   Kraus operator is the projector P onto that eigenvalue's eigenspace. So for an op on
   exactly the measured qudits the outcome labels are the eigenvalues, the probabilities
   ||P psi||^2 sum to 1, the branch of an outcome of nonzero probability is
   P psi / ||P psi||, and "Mean" is <psi|op|psi>. Most of the observables have repeated
   eigenvalues, where this differs from the eigenbasis measurement
   QuantumMeasurementOperator[op], and the rotated ones give Eigensystem a non-orthogonal
   basis of each eigenspace. *)

(* the outcome values of qmo[qs] minus values, then the probabilities minus ||P psi||^2,
   their sum minus 1, each branch of nonzero probability normalized minus P psi / ||P psi||,
   and "Mean" minus <psi|o|psi>, for the projectors ps listed in the order of values *)
luedersResiduals[qmo_, o_, qs_, values_, ps_] := With[
    {qm = qmo[qs], psi = Normal[qs["StateVector"]]},
    {p = qm["ProbabilitiesList"]},
    Flatten[{
        (First /@ qm["EigenvalueVectors"]) - values,
        p - (Norm[# . psi] ^ 2 & /@ ps),
        Total[p] - 1,
        MapThread[If[PossibleZeroQ[#3], {}, Normalize[Normal[#1["StateVector"]]] - Normalize[#2 . psi]] &, {qm["States"], ps, p}],
        qm["Mean"] - Conjugate[psi] . o . psi
    }]
]

(* the exact probabilities come back as nested radicals, which PossibleZeroQ decides in
   minutes unless they are simplified first *)
luedersExactQ[residuals_] := AllTrue[Simplify[residuals], PossibleZeroQ[#, Method -> "ExactAlgebraics"] &]

(* ZZ has the eigenvalues -1 and 1, each twice, and |Phi+> lies in the eigenspace of 1: it
   is found there with certainty and left unchanged *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Lueders"[QuantumOperator["ZZ"]]][QuantumState["PhiPlus"]]},
        {First /@ qm["EigenvalueVectors"], qm["ProbabilitiesList"], qm["Mean"], Normal[Last[qm["States"]]["StateVector"]]}
    ],
    {{-1, 1}, {0, 1}, 1, {1, 0, 0, 1} / Sqrt[2]},
    TestID -> "Lueders-ZZ-PhiPlus"
]

(* its Kraus operators, under the outcome labels -1 and 1, are the projectors onto the two
   eigenspaces of ZZ *)
VerificationTest[
    KeyMap[
        Replace[QuditName[Interpretation[_, {value_, _}], ___] :> value],
        Normal[#["MatrixRepresentation"]] & /@ QuantumMeasurementOperator["Lueders"[QuantumOperator["ZZ"]]]["Operators"]
    ],
    <|-1 -> DiagonalMatrix[{0, 1, 1, 0}], 1 -> DiagonalMatrix[{1, 0, 0, 1}]|>,
    TestID -> "Lueders-ZZ-Operators-are-eigenspace-projectors"
]

(* the eigenbasis measurement is still the default: one outcome per eigenvector *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator[QuantumOperator["ZZ"]][QuantumState["PhiPlus"]]},
        {First /@ qm["EigenvalueVectors"], qm["ProbabilitiesList"], Normal[#["StateVector"]] & /@ qm["States"]}
    ],
    {{-1, -1, 1, 1}, {0, 0, 1/2, 1/2}, {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 1/Sqrt[2]}, {1/Sqrt[2], 0, 0, 0}}},
    TestID -> "Lueders-default-eigenbasis-measurement-unchanged"
]

(* exact: the eigenspaces of diag(1, 1, -1, -1) rotated by the Fourier matrix, on a state
   with weight in both *)
VerificationTest[
    With[{f = FourierMatrix[4]}, {o = f . DiagonalMatrix[{1, 1, -1, -1}] . ConjugateTranspose[f]},
        luedersExactQ @ luedersResiduals[
            QuantumMeasurementOperator["Lueders"[QuantumOperator[o, {1, 2}]]],
            o,
            QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}],
            {-1, 1},
            f . DiagonalMatrix[#] . ConjugateTranspose[f] & /@ {{0, 0, 1, 1}, {1, 1, 0, 0}}
        ]
    ],
    True,
    TestID -> "Lueders-ExactRotatedDegenerate"
]

(* Machine precision: this rotation fails HermitianMatrixQ from roundoff alone, so
   Eigensystem takes its general solver, whose eigenvectors for a repeated eigenvalue are
   not orthogonal, and the two eigenvalues of each eigenspace differ in their last bits. *)
VerificationTest[
    With[{u = BlockRandom[Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {4, 4}]], RandomSeeding -> 35]},
        {o = ConjugateTranspose[u] . DiagonalMatrix[N @ {1, 1, -1, -1}] . u},
        {
            HermitianMatrixQ[o],
            Max[Abs[luedersResiduals[
                QuantumMeasurementOperator["Lueders"[QuantumOperator[o, {1, 2}]]],
                o,
                QuantumState[N @ {1, 2 I, -1, 3} / Sqrt[15], {2, 2}],
                {-1, 1},
                ConjugateTranspose[u] . DiagonalMatrix[#] . u & /@ {{0, 0, 1, 1}, {1, 1, 0, 0}}
            ]]] < 10 ^ -8
        }
    ],
    {False, True},
    TestID -> "Lueders-MachineRotatedDegenerate"
]

(* a qutrit: diag(1, 1, 0) rotated by the 3 x 3 Fourier matrix *)
VerificationTest[
    With[{f = FourierMatrix[3]}, {o = f . DiagonalMatrix[{1, 1, 0}] . ConjugateTranspose[f]},
        luedersExactQ @ luedersResiduals[
            QuantumMeasurementOperator["Lueders"[QuantumOperator[o, {1}, 3]]],
            o,
            QuantumState[{1, 2, 2 I} / 3, 3],
            {0, 1},
            f . DiagonalMatrix[#] . ConjugateTranspose[f] & /@ {{0, 0, 1}, {1, 1, 0}}
        ]
    ],
    True,
    TestID -> "Lueders-ExactRotatedDegenerateQutrit"
]

(* symbolic eigenvalues a and b on exact eigenspaces: the outcomes are labeled a and b *)
VerificationTest[
    Module[{a, b},
        With[{f = FourierMatrix[4], psi = {1, 2 I, -1, 3} / Sqrt[15]},
            {pa = f . DiagonalMatrix[{1, 1, 0, 0}] . ConjugateTranspose[f], pb = f . DiagonalMatrix[{0, 0, 1, 1}] . ConjugateTranspose[f]},
            With[{qm = QuantumMeasurementOperator["Lueders"[QuantumOperator[a pa + b pb, {1, 2}]]][QuantumState[psi, {2, 2}]]},
                {
                    KeySort[AssociationThread[First /@ qm["EigenvalueVectors"], Simplify[qm["ProbabilitiesList"]]]] ===
                        KeySort[Simplify /@ <|a -> Norm[pa . psi] ^ 2, b -> Norm[pb . psi] ^ 2|>],
                    Simplify[qm["Mean"] - Conjugate[psi] . (a pa + b pb) . psi]
                }
            ]
        ]
    ],
    {True, 0},
    TestID -> "Lueders-SymbolicEigenvalues"
]

(* a mixed state rho: the probabilities are Tr[P rho] and the branches P.rho.P *)
VerificationTest[
    With[{f = FourierMatrix[4], psi = {1, 2 I, -1, 3} / Sqrt[15]},
        {rho = 3/4 KroneckerProduct[psi, Conjugate[psi]] + IdentityMatrix[4] / 16, ps = f . DiagonalMatrix[#] . ConjugateTranspose[f] & /@ {{0, 0, 1, 1}, {1, 1, 0, 0}}},
        With[{qm = QuantumMeasurementOperator["Lueders"[QuantumOperator[f . DiagonalMatrix[{1, 1, -1, -1}] . ConjugateTranspose[f], {1, 2}]]][QuantumState[rho, {2, 2}]]},
            luedersExactQ @ Flatten[{
                qm["ProbabilitiesList"] - (Tr[# . rho] & /@ ps),
                (Normal[#["DensityMatrix"]] & /@ qm["States"]) - (# . rho . # & /@ ps)
            }]
        ]
    ],
    True,
    TestID -> "Lueders-MixedState-branches"
]

(* without a repeated eigenvalue each eigenspace is one eigenvector, and the measurement is
   the eigenbasis measurement *)
VerificationTest[
    With[{f = FourierMatrix[4], qs = QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]},
        {o = QuantumOperator[f . DiagonalMatrix[{3, 1, -1, -2}] . ConjugateTranspose[f], {1, 2}]},
        With[{lueders = QuantumMeasurementOperator["Lueders"[o]][qs], eigenbasis = QuantumMeasurementOperator[o][qs]},
            luedersExactQ @ Flatten[{
                (First /@ lueders["EigenvalueVectors"]) - (First /@ eigenbasis["EigenvalueVectors"]),
                lueders["ProbabilitiesList"] - eigenbasis["ProbabilitiesList"],
                (Normal[#["StateVector"]] & /@ lueders["States"]) - (Normal[#["StateVector"]] & /@ eigenbasis["States"])
            }]
        ]
    ],
    True,
    TestID -> "Lueders-nondegenerate-is-eigenbasis-measurement"
]

(* "SuperOperator" chops the entries of the operator it diagonalizes below 10^-10, so inexact
   eigenvalues are one outcome within the larger of twice what that removes and the
   roundoff 100 n 10^-p of the largest: 1 and 1 + 10^-14 are one, 1 and 1 + 10^-6 two *)
VerificationTest[
    With[{f = FourierMatrix[4]},
        {Head[#], Length[#["Eigenvalues"]]} & @ QuantumMeasurementOperator["Lueders"[QuantumOperator[N[f . DiagonalMatrix[{1, 1 + #, -1, -1}] . ConjugateTranspose[f]], {1, 2}]]] & /@ {10 ^ -14, 10 ^ -6}
    ],
    {{QuantumMeasurementOperator, 2}, {QuantumMeasurementOperator, 3}},
    TestID -> "Lueders-machine-eigenvalue-grouping"
]

(* where the Chop removes nothing the tolerance is the roundoff alone, so the eigenvalues
   1 - 10^-11 and 1 + 10^-11 are two outcomes *)
VerificationTest[
    Length[QuantumMeasurementOperator["Lueders"[QuantumOperator[N[IdentityMatrix[2] + 10 ^ -11 PauliMatrix[3]]]]]["Eigenvalues"]],
    2,
    TestID -> "Lueders-close-eigenvalues-nothing-chopped"
]

(* an operator with entries 10^12 and 1 is measured as it is, so the Chop leaves its
   eigenvalues -1 and 1 apart *)
VerificationTest[
    First /@ QuantumMeasurementOperator["Lueders"[QuantumOperator[N[DiagonalMatrix[{10 ^ 12, 10 ^ 12, 1, -1}]], {1, 2}]]]["EigenvalueVectors"],
    {-1., 1., 1.*^12},
    SameTest -> (Max[Abs[#1 - #2] / Abs[#2]] < 10 ^ -8 &),
    TestID -> "Lueders-wide-range-eigenvalues"
]

(* an operator whose entries are all small is measured multiplied up to largest entry 1, so
   an observable in SI units keeps its two outcomes and its eigenvalue labels *)
VerificationTest[
    With[{f = FourierMatrix[4]},
        With[{qmo = QuantumMeasurementOperator["Lueders"[QuantumOperator[N[10 ^ -12 f . DiagonalMatrix[{1, 1, -1, -1}] . ConjugateTranspose[f]], {1, 2}]]]},
            Max[Abs[(First /@ qmo["EigenvalueVectors"]) - {-10 ^ -12, 10 ^ -12}]] < 10 ^ -24
        ]
    ],
    True,
    TestID -> "Lueders-small-scale-observable"
]

(* ZX on the order {2, 1} is X on qudit 1 and Z on qudit 2, so |+0> is its eigenstate for 1 *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Lueders"[QuantumOperator["ZX", {2, 1}]]][QuantumState["+0"]]},
        {First /@ qm["EigenvalueVectors"], qm["ProbabilitiesList"]}
    ],
    {{-1, 1}, {0, 1}},
    TestID -> "Lueders-unsorted-order"
]

(* Z (x) diag(1, 1, 2) on a qubit and a qutrit, with its input legs listed in the order {2, 1}:
   the dimensions read {2, 3} out and {3, 2} in, and it is measured as its sorted form *)
VerificationTest[
    With[{sorted = QuantumOperator[KroneckerProduct[PauliMatrix[3], DiagonalMatrix[{1, 1, 2}]], {1, 2}, {2, 3}]},
        With[{permuted = QuantumOperator[sorted["PermuteInput", Cycles[{{1, 2}}]], {{1, 2}, {2, 1}}]},
            {
                permuted["InputDimensions"],
                First /@ QuantumMeasurementOperator["Lueders"[permuted]]["EigenvalueVectors"],
                First /@ QuantumMeasurementOperator["Lueders"[sorted]]["EigenvalueVectors"]
            }
        ]
    ],
    {{3, 2}, {-2, -1, 1, 2}, {-2, -1, 1, 2}},
    TestID -> "Lueders-permuted-input-order-mixed-dimensions"
]

(* an integer target inside op's order is the target {1}: Z (x) I on qudits 1 and 2, measured
   on qudit 1 of |10>, gives its negative eigenvalue with certainty (keyed by the sign of
   the label, whose size depends on how the partial trace over qudit 2 is normalized) *)
VerificationTest[
    With[{qmo = QuantumMeasurementOperator["Lueders"[QuantumOperator[KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]], {1, 2}]], 1]},
        With[{qm = qmo[QuantumState["10"]]},
            {qmo["Target"], KeySort[AssociationThread[Sign[First /@ qm["EigenvalueVectors"]], qm["ProbabilitiesList"]]]}
        ]
    ],
    {{1}, <|-1 -> 1, 1 -> 0|>},
    TestID -> "Lueders-integer-target"
]

(* in a circuit, "Lueders"[op] -> order places op on order and measures it as one observable,
   and so does a target outside op's qudits, as {2, 3} is for "ZZ", built on {1, 2} *)
VerificationTest[
    {First /@ #["EigenvalueVectors"], #["ProbabilitiesList"]} & /@ {
        QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "Lueders"["ZZ"] -> {1, 2}}][],
        QuantumCircuitOperator[{"H" -> 2, "CNOT" -> {2, 3}, "Lueders"["ZZ"] -> {2, 3}}][],
        QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1, 2} -> "Lueders"["ZZ"]}][],
        QuantumCircuitOperator[{"H" -> 2, "CNOT" -> {2, 3}, {2, 3} -> "Lueders"["ZZ"]}][]
    },
    ConstantArray[{{-1, 1}, {0, 1}}, 4],
    TestID -> "Lueders-circuit-specs"
]

(* a non-normal operator has no orthogonal eigenspaces, exact or at 30 digits *)
VerificationTest[
    FailureQ @ QuantumMeasurementOperator["Lueders"[QuantumOperator[{{1, 1}, {0, 2}}]]],
    True,
    {QuantumMeasurementOperator::luedersnotnormal},
    TestID -> "Lueders-non-normal-fails"
]

VerificationTest[
    FailureQ @ QuantumMeasurementOperator["Lueders"[QuantumOperator[N[{{1, 1}, {-1, 2}}, 30]]]],
    True,
    {QuantumMeasurementOperator::luedersnotnormal},
    TestID -> "Lueders-non-normal-arbitrary-precision-fails"
]

(* a target outside op's qudits must have as many qudits as op, and a target lists distinct
   qudits *)
VerificationTest[
    FailureQ /@ {QuantumMeasurementOperator["Lueders"["ZZ"], {3}], QuantumMeasurementOperator["Lueders"["ZZ"], {1, 1}]},
    {True, True},
    {QuantumMeasurementOperator::luederstarget, QuantumMeasurementOperator::luederstarget},
    TestID -> "Lueders-invalid-target-fails"
]

(* an operator from qudit 1 to qudit 2 is not an observable of qudit 1 *)
VerificationTest[
    FailureQ @ QuantumMeasurementOperator["Lueders"[QuantumOperator[PauliMatrix[3], {{2}, {1}}]]],
    True,
    {QuantumMeasurementOperator::luedersnotsquare},
    TestID -> "Lueders-output-qudits-not-input-qudits-fails"
]

(* a multiple of the identity has a single outcome, which leaves no eigenqudit *)
VerificationTest[
    FailureQ @ QuantumMeasurementOperator["Lueders"[QuantumOperator[3 IdentityMatrix[4], {1, 2}]]],
    True,
    {QuantumMeasurementOperator::luedersoneoutcome},
    TestID -> "Lueders-identity-one-outcome-fails"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - basic action"]

(* projector measurement: probabilities of |+> in computational basis *)
VerificationTest[
    Chop @ Values @ QuantumMeasurementOperator[QuantumOperator["I"]][QuantumState["+"]]["Probabilities"],
    {1/2, 1/2},
    TestID -> "Plus-projector-probabilities"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - Computational property"]

(* gh#3: ["Computational"] used to silently corrupt non-Z measurement statistics
   by rebuilding the QMO with the operator re-seated into the computational frame,
   stripping the eigenvalue-tagged pointer labels. It should now leave a non-
   computational measurement untouched. *)
VerificationTest[
    Values @ QuantumMeasurementOperator["X", {1}]["Computational"][QuantumState["+"]]["Probabilities"],
    {1, 0},
    TestID -> "Computational-X-plus"
]

VerificationTest[
    Values @ QuantumMeasurementOperator["X", {1}]["Computational"][QuantumState["-"]]["Probabilities"],
    {0, 1},
    TestID -> "Computational-X-minus"
]

VerificationTest[
    Values @ QuantumMeasurementOperator["Y", {1}]["Computational"][QuantumState["+"]]["Probabilities"],
    {1/2, 1/2},
    TestID -> "Computational-Y-plus"
]

(* a Z measurement is already in the computational frame and continues to be
   rebuilt through the standard path *)
VerificationTest[
    Values @ QuantumMeasurementOperator["Z", {1}]["Computational"][QuantumState["0"]]["Probabilities"],
    {1, 0},
    TestID -> "Computational-Z-zero"
]

(* circuit-level: QuantumCircuitOperator[...]["Computational"] maps the property
   over each element; the non-Z measurement inside must stay faithful *)
VerificationTest[
    Values @ QuantumCircuitOperator[{"H" -> 1, QuantumMeasurementOperator["X", {1}]}]["Computational"][QuantumState["0"]]["Probabilities"],
    {1, 0},
    TestID -> "Computational-circuit-H-then-X"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - failure"]

VerificationTest[
    QuantumMeasurementOperator["NotAMeasurement"[2]],
    Failure["InvalidName", _],
    {QuantumMeasurementOperator::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

VerificationTest[
    QuantumMeasurementOperator["HesseSICPOVM"["bad-arg"]],
    Failure["InvalidArguments", _],
    {QuantumMeasurementOperator::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-HesseSICPOVM"
]

(* A bare string here is a BASIS specification, not a measurement name: anything
   QuantumBasis accepts measures in that basis. So the unrecognized bare name is
   adjudicated by QuditBasis, which owns the basis-name registry, and the message
   is QuditBasis::invalidName rather than a second copy of the registry here. *)
VerificationTest[
    FailureQ @ QuantumMeasurementOperator["NotAMeasurement"],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-bare-form"
]

VerificationTest[
    FailureQ @ QuantumMeasurementOperator["NotAMeasurement", {2}],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-bare-form-with-target"
]

(* Which is why bare basis names must keep building measurements: a guard that
   rejected every string outside $QuantumMeasurementOperatorNames would break all
   of these. Pauli strings and "Basis" suffixes reach here through QuditBasis too. *)
VerificationTest[
    {
        QuantumMeasurementOperator["PauliZ"]["Dimension"],
        QuantumMeasurementOperator["X"]["Dimension"],
        QuantumMeasurementOperator["Bell"]["Dimension"],
        QuantumMeasurementOperator["XYZ"]["Target"],
        QuantumMeasurementOperator["PauliBasis"]["Dimension"]
    },
    {8, 8, 64, {1, 2, 3}, 64},
    {},
    TestID -> "BasisName-bare-form-still-measures"
]

(* Registered measurement names are unaffected. *)
VerificationTest[
    {
        QuantumMeasurementOperator["M"]["Dimension"],
        QuantumMeasurementOperator["TetrahedronSICPOVM"]["Dimension"]
    },
    {8, 16},
    {},
    TestID -> "ValidName-registered-measurements-unaffected"
]

EndTestSection[]


BeginTestSection["QuantumMeasurement - Entropy"]

(* Shannon entropy of the outcome distribution, H = -Sum p_i Log_b p_i, in closed form
   (mirrors QuantumState["VonNeumannEntropy", logBase]). Computing it via
   Information[CategoricalDistribution[..]] is unsafe: on symbolic / non-numeric
   probabilities (any parametric circuit) that path falls back to Monte-Carlo sampling
   and emits RandomVariate::unsdst + Extract::psl1, including on plain display. *)

(* --- known numeric values --- *)

VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT", {1, 2}}][]["Entropy"],
    Quantity[1, "Bits"],
    TestID -> "Entropy-Bell-1bit"
]

VerificationTest[
    QuantumCircuitOperator[{"H", {1}}][]["Entropy"],
    Quantity[1, "Bits"],
    TestID -> "Entropy-PlusState-1bit"
]

(* deterministic outcome carries no information *)
VerificationTest[
    QuantumMeasurementOperator["Z"][QuantumState["0"]]["Entropy"],
    Quantity[0, "Bits"],
    TestID -> "Entropy-Deterministic-0bit"
]

(* biased {1/4, 3/4}: H = 2 - (3/4) Log2[3], exact *)
VerificationTest[
    Simplify[QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["Entropy"] - Quantity[2 - 3/4 Log2[3], "Bits"]],
    Quantity[0, "Bits"],
    TestID -> "Entropy-Biased-ExactShannon"
]

(* --- logBase argument: parity with QuantumState["Entropy", logBase] (QM lacked it) --- *)

(* base e returns the entropy in nats *)
VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT", {1, 2}}][]["Entropy", E],
    Log[2],
    TestID -> "Entropy-LogBase-Nats"
]

(* {"Entropy", logBase} list-form dispatches to the bare base-b value *)
VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT", {1, 2}}][][{"Entropy", 2}],
    1,
    TestID -> "Entropy-ListForm-Base2"
]

(* base conversion is consistent on a non-trivial distribution: nats = bits * Ln[2] *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]},
        FullSimplify[qm["Entropy", E] - qm["Entropy", 2] Log[2]]
    ],
    0,
    TestID -> "Entropy-LogBase-Consistency"
]

(* --- symbolic robustness (the reported regression) --- *)

(* a bare parametric rotation measured in Z keeps symbolic probabilities (theta is not
   assumed real); qm["Entropy"] must stay a Quantity and emit NO messages. The empty
   expected-message list makes any leaked message fail the test. *)
VerificationTest[
    Module[{t}, Head @ QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["Entropy"]],
    Quantity,
    {},
    TestID -> "Entropy-Symbolic-NoMessages"
]

(* the summary box evaluates N @ qm["Entropy"] on display, so the display path must be
   quiet too *)
VerificationTest[
    Module[{t}, ToBoxes @ QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]; True],
    True,
    {},
    TestID -> "Entropy-Symbolic-Display-NoMessages"
]

(* symbolic entropy is correct, not merely quiet: at theta = Pi/2 the Z outcomes are
   equiprobable, so H = 1 bit *)
VerificationTest[
    Module[{t}, Simplify[QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["Entropy"] /. t -> Pi / 2]],
    Quantity[1, "Bits"],
    TestID -> "Entropy-Symbolic-CorrectAtPiOver2"
]

EndTestSection[]


BeginTestSection["QuantumMeasurement - Symbolic Probabilities"]

(* Outcome ordering ("TopProbabilities") and Monte-Carlo sampling ("SimulatedMeasurement",
   "SimulatedCounts", "SimulatedStateMeasurement"), together with the generic
   "DistributionInformation" passthrough, are defined only when every outcome probability
   is numeric. On a symbolic / parametric measurement the old code fed a non-numeric
   distribution into Information / RandomVariate and leaked CategoricalDistribution::elmntavsl,
   MultinomialDistribution::vprobprm, Extract::psl1, RandomVariate::unsdst, KeyMap::invak, etc.,
   including on plain display. These now return Indeterminate (mirroring "Entropy"). The empty
   expected-message list makes any leaked message fail the test. *)

(* --- the reported regression: no leaked messages, Indeterminate value --- *)

(* TopProbabilities needs to order the weights (was CategoricalDistribution::elmntavsl) *)
VerificationTest[
    Module[{t}, QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["TopProbabilities"]],
    Indeterminate,
    {},
    TestID -> "SymProb-TopProbabilities"
]

(* SimulatedCounts samples a MultinomialDistribution (was MultinomialDistribution::vprobprm) *)
VerificationTest[
    Module[{t}, QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["SimulatedCounts", 10]],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedCounts"
]

(* SimulatedMeasurement samples the CategoricalDistribution; previously returned an
   unevaluated RandomVariate[...] instead of a clean value *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedMeasurement"],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedMeasurement"
]

VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedMeasurement", 3],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedMeasurement-n"
]

(* the explicit DistributionInformation path with a sampling sub-property: the exact pair
   from the Entropy regression, Extract::psl1 + RandomVariate::unsdst *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["DistributionInformation", "Entropy"],
    Indeterminate,
    {},
    TestID -> "SymProb-DistributionInformation-Entropy"
]

(* same path under N, the route the summary box used to take on display *)
VerificationTest[
    N @ QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["DistributionInformation", "Entropy"],
    Indeterminate,
    {},
    TestID -> "SymProb-DistributionInformation-Entropy-N"
]

(* SimulatedStateMeasurement keys the state association by a simulated outcome *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedStateMeasurement"],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedStateMeasurement"
]

(* list form was RandomVariate::array + Part::pkspec1 *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedStateMeasurement", 3],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedStateMeasurement-n"
]

(* TopStateProbabilities depends on the (now guarded) TopProbabilities; was KeyMap::invak *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["TopStateProbabilities"],
    Indeterminate,
    {},
    TestID -> "SymProb-TopStateProbabilities"
]

(* display path stays quiet on a symbolic measurement *)
VerificationTest[
    Module[{t}, ToBoxes @ QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]; True],
    True,
    {},
    TestID -> "SymProb-Display"
]

(* --- properties that are well defined symbolically must NOT degrade to Indeterminate --- *)

(* Categories is just the outcome support, independent of the probability values *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["Categories"],
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["Outcomes"],
    {},
    TestID -> "SymProb-Categories-Preserved"
]

(* ProbabilityArray is the symbolic weight vector, fully resolved (no leftover Information) *)
VerificationTest[
    With[{a = QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["ProbabilityArray"]},
        Head[a] === List && FreeQ[a, _Information | _CategoricalDistribution]
    ],
    True,
    {},
    TestID -> "SymProb-ProbabilityArray-Preserved"
]

VerificationTest[
    Head @ QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["ProbabilityTable"],
    Dataset,
    {},
    TestID -> "SymProb-ProbabilityTable-Preserved"
]

(* --- numeric measurements are unchanged and stay message-free --- *)

VerificationTest[
    Length @ QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["TopProbabilities"],
    2,
    {},
    TestID -> "SymProb-Numeric-TopProbabilities"
]

VerificationTest[
    Total @ QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["SimulatedCounts", 50],
    50,
    {},
    TestID -> "SymProb-Numeric-SimulatedCounts"
]

VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]},
        MemberQ[qm["Outcomes"], qm["SimulatedMeasurement"]]
    ],
    True,
    {},
    TestID -> "SymProb-Numeric-SimulatedMeasurement"
]

VerificationTest[
    NumericQ @ QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["DistributionInformation", "Entropy"],
    True,
    {},
    TestID -> "SymProb-Numeric-DistributionInformation"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - degenerate observable eigenbasis"]

(* A rank-1 projector on a qutrit has a degenerate 0-eigenspace; at arbitrary precision
   Eigensystem returns non-orthonormal eigenvectors there, and projectors built from the
   skewed basis silently corrupted the outcome probabilities for some KCBS pentagon
   projectors (k = 2, 3) while others were fine. All five expectation values are exactly
   Cos[alpha]^2 = 1/Sqrt[5] by symmetry. *)
VerificationTest[
    Block[{sol, alpha, v, P, qs},
        sol = Solve[Cos[a]^2 + Sin[a]^2 Cos[4 Pi/5] == 0 && 0 < a < Pi/2, a];
        alpha = a /. First[sol];
        v[k_] := N[{Cos[alpha], Sin[alpha] Cos[4 Pi k/5], Sin[alpha] Sin[4 Pi k/5]}, 16];
        P[k_] := Outer[Times, v[k], v[k]];
        qs = QuantumState[{1, 0, 0}, {3}];
        Max @ Abs[N @ Table[QuantumMeasurementOperator[P[k]][qs]["Mean"], {k, 0, 4}] - N[1/Sqrt[5]]] < 1*^-10
    ],
    True,
    TestID -> "DegenerateEigenbasis-KCBS-arbitrary-precision"
]

(* same pentagon at machine precision *)
VerificationTest[
    Block[{sol, alpha, v, P, qs},
        sol = Solve[Cos[a]^2 + Sin[a]^2 Cos[4 Pi/5] == 0 && 0 < a < Pi/2, a];
        alpha = a /. First[sol];
        v[k_] := N @ {Cos[alpha], Sin[alpha] Cos[4 Pi k/5], Sin[alpha] Sin[4 Pi k/5]};
        P[k_] := Outer[Times, v[k], v[k]];
        qs = QuantumState[{1, 0, 0}, {3}];
        Max @ Abs[Table[QuantumMeasurementOperator[P[k]][qs]["Mean"], {k, 0, 4}] - N[1/Sqrt[5]]] < 1*^-8
    ],
    True,
    TestID -> "DegenerateEigenbasis-KCBS-machine-precision"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - order re-seating on direct state application"]

(* A tensor-product observable whose factors were built with off-1 wire labels (output
   wire {3} each, renumbered to {3, 4} by the tensor product) failed with
   QuantumCircuitOperator::dim when applied directly to a two-qutrit state: phantom
   wires 1-2 were padded at qubit dimension 2. A directly applied measurement has no
   circuit context, so wire labels that do not fit the state are re-seated onto the
   state's qudits. *)
VerificationTest[
    Block[{P0, opAB, qmo, w, qs},
        P0 = Outer[Times, {1., 0., 0.}, {1., 0., 0.}];
        opAB = QuantumTensorProduct[QuantumOperator[P0, {3}, {1}], QuantumOperator[P0, {3}, {2}]];
        qmo = QuantumMeasurementOperator[opAB];
        w = Normalize[{1., 2., 0.}];
        qs = QuantumState[Flatten[Outer[Times, w, w]], {3, 3}];
        {Head[qmo[qs]], Abs[N @ qmo[qs]["Mean"] - 1/25] < 1*^-8}
    ],
    {QuantumMeasurement, True},
    TestID -> "OrderReseat-tensor-product-observable"
]

(* a measurement whose order fits inside a larger state is untouched by the re-seat *)
VerificationTest[
    Round[N @ Values @ QuantumMeasurementOperator["Z", {2}][
        QuantumCircuitOperator[{"X" -> 2}][QuantumState["000"]]]["Probabilities"], 0.001],
    {0., 1.},
    TestID -> "OrderReseat-preserves-targeted-measurement"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - operator stored on an unsorted qudit order"]

(* A measurement of an operator measures that operator, whatever order its qudits are stored
   in and whatever order its target lists them in. Every reference comes from the operator's
   own "MatrixRepresentation" m, which lists the qudits in sorted order: the mean <v|m|v>,
   and the Born weight <v|P|v> of each eigenvalue lambda, with P the eigenspace projector
   Prod_{mu != lambda} (m - mu)/(lambda - mu). P needs no eigenvectors, so a degenerate
   observable has a reference too. "ZX" on {2, 1} is X on qudit 1 and Z on qudit 2: |+0> is
   its +1 eigenstate and |0+> has mean 0, where the qudit-swapped Z1 X2 gives 0 and 1. *)

unsortedOrderEigenprojectors[m_] := With[{id = IdentityMatrix[Length[m]], values = Union[Eigenvalues[m]]},
    AssociationMap[l |-> Fold[#1 . (m - #2 id) / (l - #2) &, id, DeleteCases[values, l]], values]
]

unsortedOrderEigenvalueWeights[m_, v_] := Simplify[Conjugate[v] . # . v] & /@ unsortedOrderEigenprojectors[m]

(* the measured weight of each eigenvalue: the probabilities of the outcomes it labels, summed *)
unsortedOrderMeasuredWeights[qm_] := KeySort @ Merge[
    KeyValueMap[Replace[#1, QuditName[Interpretation[_, {l_, _}], ___] :> l] -> #2 &, qm["Probabilities"]],
    Total
]

(* the eigenvalue labels of a measurement's outcomes, sorted *)
measurementOutcomeEigenvalues[qmo_] := Sort[Replace[qmo["Eigenvalues"], QuditName[Interpretation[_, {l_, _}], ___] :> l, {1}]]

(* whether the mean and the eigenvalue weights of qmo in psi are those of the operator op *)
unsortedOrderReadoutQ[qmo_, op_, psi_] := With[
    {m = Normal[op["MatrixRepresentation"]], v = Normal[psi["StateVector"]], qm = qmo[psi]},
    With[{measured = unsortedOrderMeasuredWeights[qm], expected = unsortedOrderEigenvalueWeights[m, v]},
        Simplify[qm["Mean"] - Conjugate[v] . m . v] === 0 &&
            Keys[measured] === Keys[expected] &&
            Simplify[Values[measured] - Values[expected]] === ConstantArray[0, Length[expected]]
    ]
]

VerificationTest[
    With[{op = QuantumOperator["ZX", {2, 1}]},
        With[{qmo = QuantumMeasurementOperator[op]},
            {
                Normal[op["MatrixRepresentation"]] === Normal[QuantumOperator["XZ", {1, 2}]["MatrixRepresentation"]],
                qmo[QuantumState["0+"]]["Mean"],
                qmo[QuantumState["+0"]]["Mean"],
                unsortedOrderMeasuredWeights[qmo[QuantumState["0+"]]],
                unsortedOrderMeasuredWeights[qmo[QuantumState["+0"]]]
            }
        ]
    ],
    {True, 0, 1, <|-1 -> 1/2, 1 -> 1/2|>, <|-1 -> 0, 1 -> 1|>},
    {},
    TestID -> "UnsortedOrder-ZX-on-21-MeasuresX1Z2"
]

(* the states on which each measurement departs from its operator: none, for an observable
   whose eigenvalues are all degenerate, a nondegenerate one (diag(1, 2) on qudit 2, X on
   qudit 1), and the first re-expressed in the PauliX basis *)
VerificationTest[
    Map[
        op |-> With[{qmo = QuantumMeasurementOperator[op]},
            Select[{"00", "0+", "+0", "+-", "LR", "1L", "R-"}, ! unsortedOrderReadoutQ[qmo, op, QuantumState[#]] &]
        ],
        {
            QuantumOperator["ZX", {2, 1}],
            QuantumOperator[KroneckerProduct[DiagonalMatrix[{1, 2}], PauliMatrix[1]], {2, 1}],
            QuantumOperator[QuantumOperator["ZX", {2, 1}], QuantumBasis["PauliX", 2]]
        }
    ],
    {{}, {}, {}},
    {},
    TestID -> "UnsortedOrder-TwoQubit-MatchesOwnMatrix"
]

(* X on qudit 3, Y on qudit 1, Z on qudit 2: the sorted operator is Y1 Z2 X3 *)
VerificationTest[
    With[{op = QuantumOperator["XYZ", {3, 1, 2}]},
        With[{qmo = QuantumMeasurementOperator[op]},
            Select[{"000", "L0+", "+L0", "0+L", "RR1", "-0L", "1-R"}, ! unsortedOrderReadoutQ[qmo, op, QuantumState[#]] &]
        ]
    ],
    {},
    {},
    TestID -> "UnsortedOrder-ThreeCycle-MatchesOwnMatrix"
]

(* a qutrit observable diag(1, 2, 3) on qudit 2 and X on qubit 1, stored on {2, 1}: the stored
   dimensions are {3, 2} and the sorted ones {2, 3}. |+> (|0> + |2>)/Sqrt[2] has mean 2. *)
VerificationTest[
    With[{op = QuantumOperator[KroneckerProduct[DiagonalMatrix[{1, 2, 3}], PauliMatrix[1]], {2, 1}, QuantumBasis[{3, 2}, {3, 2}]]},
        With[{qmo = QuantumMeasurementOperator[op], states = {
                QuantumState[Flatten[KroneckerProduct[{1, 1} / Sqrt[2], {1, 0, 1} / Sqrt[2]]], {2, 3}],
                QuantumState[Flatten[KroneckerProduct[{1, 0}, {1, 1, 1} / Sqrt[3]]], {2, 3}],
                QuantumState[Flatten[KroneckerProduct[{1, -I} / Sqrt[2], {0, 1, 1} / Sqrt[2]]], {2, 3}]
            }},
            {
                Normal[op["MatrixRepresentation"]] === KroneckerProduct[PauliMatrix[1], DiagonalMatrix[{1, 2, 3}]],
                qmo[First[states]]["Mean"],
                unsortedOrderReadoutQ[qmo, op, #] & /@ states
            }
        ]
    ],
    {True, 2, {True, True, True}},
    {},
    TestID -> "UnsortedOrder-MixedDimensions-MatchesOwnMatrix"
]

(* a partial target with the traced qudit between the targets: Z on 3, I on 2 and X on 1,
   measured on {3, 1}, is X1 Z3, whose mean and eigenvalue weights are those of the whole
   operator *)
VerificationTest[
    With[{op = QuantumOperator["ZIX", {3, 2, 1}]},
        With[{qmo = QuantumMeasurementOperator[op, {3, 1}]},
            Select[{"000", "+00", "00+", "-01", "L1R", "+1-"}, ! unsortedOrderReadoutQ[qmo, op, QuantumState[#]] &]
        ]
    ],
    {},
    {},
    TestID -> "UnsortedOrder-PartialTarget-TracedQuditBetween"
]

(* the target names the measured qudits: listing them as {2, 1} measures the same Z1 X2 *)
VerificationTest[
    With[{op = QuantumOperator["ZX"]},
        With[{qmo = QuantumMeasurementOperator[op, {2, 1}]},
            {
                qmo[QuantumState["0+"]]["Mean"],
                qmo[QuantumState["+0"]]["Mean"],
                qmo == QuantumMeasurementOperator[op]
            }
        ]
    ],
    {1, 0, True},
    {},
    TestID -> "UnsortedOrder-TargetListOrder-SameObservable"
]

(* the same observable stored on two orders is one measurement: equal under ==, which
   compares eigenvectors, with the same eigenvalue labels, mean and eigenvalue weights *)
VerificationTest[
    With[{
        qmo = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]],
        twin = QuantumMeasurementOperator[QuantumOperator["XZ", {1, 2}]],
        psi = QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
    },
        {
            qmo == twin,
            measurementOutcomeEigenvalues[qmo] === measurementOutcomeEigenvalues[twin],
            Simplify[qmo[psi]["Mean"] - twin[psi]["Mean"]] === 0,
            unsortedOrderMeasuredWeights[qmo[psi]] === unsortedOrderMeasuredWeights[twin[psi]]
        }
    ],
    {True, True, True, True},
    {},
    TestID -> "UnsortedOrder-SameObservableTwoOrders-Equal"
]

(* the forms built from the measurement's dilation keep that agreement: its "POVM" form is
   the same measurement, its "SuperOperator", "Dagger", "POVM", "Conjugate", "Dual" and
   "Double" forms equal those of the sorted twin, and a tensor product with another
   measurement gives the twin's outcome probabilities in the same order. The last entry is
   the qudit-swapped observable, whose "POVM" form differs. *)
VerificationTest[
    With[{
        qmo = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]],
        twin = QuantumMeasurementOperator[QuantumOperator["XZ", {1, 2}]],
        swapped = QuantumMeasurementOperator[QuantumOperator["ZX", {1, 2}]],
        psi = QuantumState[{1, 2, 0, 1, 3, 1, 2, 1} / Sqrt[21]]
    },
        {
            qmo["POVM"] == qmo,
            qmo["SuperOperator"] == twin["SuperOperator"],
            qmo["Dagger"] == twin["Dagger"],
            qmo["POVM"] == twin["POVM"],
            qmo["Conjugate"] == twin["Conjugate"],
            qmo["Dual"] == twin["Dual"],
            qmo["Double"] == twin["Double"],
            QuantumTensorProduct[QuantumMeasurementOperator[{1}], qmo][psi]["ProbabilitiesList"] ===
                QuantumTensorProduct[QuantumMeasurementOperator[{1}], twin][psi]["ProbabilitiesList"],
            qmo["POVM"] == swapped["POVM"]
        }
    ],
    {True, True, True, True, True, True, True, True, False},
    {},
    TestID -> "UnsortedOrder-DilationFormsMatchSortedTwin"
]

(* an operator stored on {2, 1} with its target listed in sorted order is the same
   measurement too, with the same "POVM" form and the same measurement channel *)
VerificationTest[
    With[{
        qmo = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}], {1, 2}],
        twin = QuantumMeasurementOperator[QuantumOperator["XZ", {1, 2}]],
        swapped = QuantumMeasurementOperator[QuantumOperator["ZX", {1, 2}]],
        psi = QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
    },
        {
            qmo == twin,
            qmo["POVM"] == qmo,
            qmo[QuantumState[#]]["Mean"] & /@ {"0+", "+0"},
            Simplify[Normal[qmo["DiscardExtraQudits"][psi]["DensityMatrix"]] - Normal[twin["DiscardExtraQudits"][psi]["DensityMatrix"]]] ===
                ConstantArray[0, {4, 4}],
            Simplify[Normal[qmo["DiscardExtraQudits"][psi]["DensityMatrix"]] - Normal[swapped["DiscardExtraQudits"][psi]["DensityMatrix"]]] ===
                ConstantArray[0, {4, 4}]
        }
    ],
    {True, True, {0, 1}, True, False},
    {},
    TestID -> "UnsortedOrder-SortedTargetOnUnsortedOperator"
]

(* "ReverseEigenQudits" puts the dilation back on the measurement's own order: it measures the
   operator, with <X1 Z2> = -2/15 on the state below where Z1 X2 gives 2/5, and equals the
   sorted twin's form *)
VerificationTest[
    With[{
        qmo = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]],
        twin = QuantumMeasurementOperator[QuantumOperator["XZ", {1, 2}]],
        swapped = QuantumMeasurementOperator[QuantumOperator["ZX", {1, 2}]],
        psi = QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
    },
        {
            Simplify[qmo["ReverseEigenQudits"][psi]["Mean"]],
            qmo["ReverseEigenQudits"] == twin["ReverseEigenQudits"],
            swapped["ReverseEigenQudits"] == twin["ReverseEigenQudits"]
        }
    ],
    {-2/15, True, False},
    {},
    TestID -> "UnsortedOrder-ReverseEigenQudits-MeasuresOperator"
]

(* composing an operator after the measurement keeps it the same measurement as the sorted
   twin's composition *)
VerificationTest[
    With[{
        qmo = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]],
        twin = QuantumMeasurementOperator[QuantumOperator["XZ", {1, 2}]],
        swapped = QuantumMeasurementOperator[QuantumOperator["ZX", {1, 2}]]
    },
        {
            QuantumOperator["H"][qmo] == QuantumOperator["H"][twin],
            QuantumOperator["H"][qmo]["QuantumOperator"] == QuantumOperator["H"][twin]["QuantumOperator"],
            QuantumOperator["H"][qmo] == QuantumOperator["H"][swapped]
        }
    ],
    {True, True, False},
    {},
    TestID -> "UnsortedOrder-OperatorAfterMeasurementMatchesTwin"
]

(* the measurement's operators are the projectors it applies: weighted by their eigenvalues
   they rebuild the operator's own matrix, and they sum to the identity *)
VerificationTest[
    With[{op = QuantumOperator["ZX", {2, 1}]},
        With[{ops = QuantumMeasurementOperator[op]["Operators"]},
            {
                Total[KeyValueMap[Replace[#1, QuditName[Interpretation[_, {l_, _}], ___] :> l] Normal[#2["MatrixRepresentation"]] &, ops]] -
                    Normal[op["MatrixRepresentation"]],
                Total[Normal[#["MatrixRepresentation"]] & /@ Values[ops]] - IdentityMatrix[4]
            }
        ]
    ],
    {ConstantArray[0, {4, 4}], ConstantArray[0, {4, 4}]},
    {},
    TestID -> "UnsortedOrder-OperatorsRebuildTheOperator"
]

(* the circuit route and the stabilizer route measure the same observable: after H on
   qudit 1 the register is |+0>, certain to give +1, outcome 0 of the tableau *)
VerificationTest[
    With[{qmo = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]]},
        {
            QuantumCircuitOperator[{"H" -> 1, qmo}][]["Mean"],
            Keys @ qmo[PauliStabilizer[2][{"H" -> 1}]]
        }
    ],
    {1, {0}},
    {},
    TestID -> "UnsortedOrder-CircuitAndStabilizerRoutesAgree"
]

(* the measurement channel of a nondegenerate observable is Sum_lambda P rho P, with the
   eigenspace projectors of the operator's own matrix: for the operator stored on {2, 1}, and
   for the sorted operator with its target listed as {2, 1} *)
VerificationTest[
    With[{op = QuantumOperator[KroneckerProduct[DiagonalMatrix[{1, 2}], PauliMatrix[1]], {2, 1}], psi = QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]},
        With[{rho = Normal[psi["DensityMatrix"]]},
            With[{channel = Total[(# . rho . #) & /@ Values[unsortedOrderEigenprojectors[Normal[op["MatrixRepresentation"]]]]]},
                Simplify[Normal[#["DiscardExtraQudits"][psi]["DensityMatrix"]] - channel] & /@ {
                    QuantumMeasurementOperator[op],
                    QuantumMeasurementOperator[op["Sort"], {2, 1}]
                }
            ]
        ]
    ],
    ConstantArray[0, {2, 4, 4}],
    {},
    TestID -> "UnsortedOrder-DiscardExtraQuditsIsTheMeasurementChannel"
]

(* a state with a free angle t: |+> on qudit 1 and Cos[t]|0> + Sin[t]|1> on qudit 2 has
   <X1 Z2> = Cos[2 t] and eigenvalue weights Sin[t]^2, Cos[t]^2, where Z1 X2 gives 0 and
   1/2, 1/2 for every t *)
VerificationTest[
    Module[{t},
        With[{qm = QuantumMeasurementOperator[QuantumOperator["ZX", {2, 1}]][
                QuantumState[Flatten[KroneckerProduct[{1, 1} / Sqrt[2], {Cos[t], Sin[t]}]]]]},
            Simplify[{qm["Mean"] - Cos[2 t], Values[unsortedOrderMeasuredWeights[qm]] - {Sin[t]^2, Cos[t]^2}}, Element[t, Reals]]
        ]
    ],
    {0, {0, 0}},
    {},
    TestID -> "UnsortedOrder-SymbolicAngle-ClosedForm"
]

(* machine precision: a Hermitian operator with a nondegenerate spectrum, stored on {2, 1} *)
VerificationTest[
    With[{op = QuantumOperator[N @ {{1, 2 I, 0, 1}, {-2 I, -1, 3, 0}, {0, 3, 2, -I}, {1, 0, I, 0}}, {2, 1}]},
        With[{qmo = QuantumMeasurementOperator[op], m = Normal[op["MatrixRepresentation"]]},
            Max @ Map[
                With[{v = Normal[#["StateVector"]]}, Abs[qmo[#]["Mean"] - Conjugate[v] . m . v]] &,
                QuantumState /@ {"00", "0+", "+0", "LR", "1-"}
            ]
        ]
    ],
    _ ? (# < 10^-8 &),
    {},
    SameTest -> MatchQ,
    TestID -> "UnsortedOrder-MachinePrecision-MeanMatchesOwnMatrix"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - partial target measures the reduced observable"]

(* A measurement whose target is part of its operator's qudits measures the partial trace
   over the other qudits divided by their dimension: A itself when the operator is A on the
   target and the identity elsewhere, and in general the part of the operator that acts as
   the identity on the traced qudits. *)

VerificationTest[
    With[{qmo = QuantumMeasurementOperator[QuantumOperator["ZI"], {1}]},
        {
            measurementOutcomeEigenvalues[qmo],
            qmo[QuantumState[#]]["Mean"] & /@ {"00", "10", "+0"}
        }
    ],
    {{-1, 1}, {1, -1, 0}},
    {},
    TestID -> "PartialTarget-ZI-on-1-MeasuresZ1"
]

(* Z (x) Z + X (x) 1 on target {1}: the Z Z part traces to zero, so the measured observable is X1 *)
VerificationTest[
    With[{qmo = QuantumMeasurementOperator[
            QuantumOperator[KroneckerProduct[PauliMatrix[3], PauliMatrix[3]] + KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]], {1, 2}],
            {1}
        ]},
        {
            measurementOutcomeEigenvalues[qmo],
            qmo[QuantumState[#]]["Mean"] & /@ {"+0", "0+", "-1"}
        }
    ],
    {{-1, 1}, {1, 0, -1}},
    {},
    TestID -> "PartialTarget-SumOperator-MeasuresIdentityPart"
]

(* X on a qubit and the identity on a traced qutrit: the division is by 3. The state
   (Cos[1/3], Sin[1/3]) on the qubit has <X> = Sin[2/3]. *)
VerificationTest[
    With[{qmo = QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[1], IdentityMatrix[3]], {1, 2}, QuantumBasis[{2, 3}, {2, 3}]], {1}]},
        {
            measurementOutcomeEigenvalues[qmo],
            Simplify[qmo[QuantumState[Flatten[KroneckerProduct[{Cos[1/3], Sin[1/3]}, {1, 1, 1} / Sqrt[3]]], {2, 3}]]["Mean"] - Sin[2/3]]
        }
    ],
    {{-1, 1}, 0},
    {},
    TestID -> "PartialTarget-TracedQutrit-DividesByThree"
]

(* an operator on the non-contiguous order {1, 3}, measured on qudit 1: the traced qudit's
   dimension is looked up by its label, not by its position *)
VerificationTest[
    With[{qmo = QuantumMeasurementOperator[QuantumOperator["XI", {1, 3}], {1}]},
        {
            measurementOutcomeEigenvalues[qmo],
            qmo[QuantumState[#]]["Mean"] & /@ {"+00", "-00", "0+1"}
        }
    ],
    {{-1, 1}, {1, -1, 0}},
    {},
    TestID -> "PartialTarget-NonContiguousOrder"
]

(* eigenvalue labels given explicitly are used as given, with no division *)
VerificationTest[
    measurementOutcomeEigenvalues[QuantumMeasurementOperator[QuantumOperator["ZI"] -> {7, 9, 11, 13}, {1}]],
    {7, 9},
    {},
    TestID -> "PartialTarget-ExplicitLabelsUnscaled"
]

(* "-> Automatic" labels the outcomes with the measured observable's own eigenvalues: those
   of the reduced observable on a partial target, and the operator's on a full target *)
VerificationTest[
    {
        measurementOutcomeEigenvalues[QuantumMeasurementOperator[QuantumOperator["ZI"] -> Automatic, {1}]],
        QuantumMeasurementOperator[QuantumOperator["ZI"] -> Automatic, {1}][QuantumState["00"]]["Mean"],
        With[{op = QuantumOperator[KroneckerProduct[DiagonalMatrix[{1, 2}], PauliMatrix[1]], {2, 1}]},
            measurementOutcomeEigenvalues[QuantumMeasurementOperator[op -> Automatic]] === Sort[Eigenvalues[Normal[op["MatrixRepresentation"]]]]
        ]
    },
    {{-1, 1}, 1, True},
    {},
    TestID -> "PartialTarget-AutomaticLabelsAreMeasuredEigenvalues"
]

EndTestSection[]


BeginTestSection["QuantumMeasurement - branch conditional states"]

(* "States" returns the conditional post-measurement branch for each outcome, the CP
   map image M_m . rho . ConjugateTranspose[M_m] (Lueders / Kraus branch), trace p_m.
   For a mixed input this was Sum_j M_m rho M_j^dag = M_m rho: the whole pointer block-row
   was summed, collapsing the bra pointer index and leaving a non-Hermitian, non-positive
   object that is not a valid state. The branch must instead be the diagonal pointer block
   M_m rho M_m^dag. These branches were undocumented and untested, so the whole contract is
   pinned here. The defect was invisible whenever rho is diagonal in the measurement
   eigenbasis (M_m rho = M_m rho M_m), which is why the common demos never caught it. *)

(* --- the reported regression: Z on a mixed state with computational-basis coherence.
   Before the fix this returned {{{1/2, 1/4}, {0, 0}}, {{0, 0}, {1/4, 1/2}}} (M_m rho). --- *)
VerificationTest[
    Normal[#["DensityMatrix"]] & /@
        QuantumMeasurementOperator["Z"][QuantumState[{{1/2, 1/4}, {1/4, 1/2}}]]["States"],
    {{{1/2, 0}, {0, 0}}, {{0, 0}, {0, 1/2}}},
    TestID -> "BranchStates-Z-mixed-Lueders-exact"
]

(* every branch is a valid density operator: Hermitian and positive semidefinite. A
   complex off-diagonal makes M_m rho manifestly non-Hermitian, so this bites hardest. *)
VerificationTest[
    With[{dms = Normal[#["DensityMatrix"]] & /@
        QuantumMeasurementOperator["Z"][QuantumState[{{1/2, I/4}, {-I/4, 1/2}}]]["States"]},
        {AllTrue[dms, HermitianMatrixQ], AllTrue[dms, PositiveSemidefiniteMatrixQ]}
    ],
    {True, True},
    TestID -> "BranchStates-Hermitian-and-PSD-complex-mixed"
]

(* branch trace equals the outcome probability, Tr(M_m rho M_m^dag) = p_m *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{{1/2, 1/4}, {1/4, 1/2}}]]},
        (Tr[Normal[#["DensityMatrix"]]] & /@ qm["States"]) == qm["ProbabilitiesList"]
    ],
    True,
    TestID -> "BranchStates-trace-equals-probability"
]

(* the branches resolve the averaged post-measurement (decohered) state:
   Sum_m M_m rho M_m^dag = PostMeasurementState. The buggy branches summed to rho itself
   (Sum_m M_m rho = rho), so this invariant fails on the old code. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{{1/2, I/4}, {-I/4, 1/2}}]]},
        Total[Normal[#["DensityMatrix"]] & /@ qm["States"]] ==
            Normal[qm["PostMeasurementState"]["DensityMatrix"]]
    ],
    True,
    TestID -> "BranchStates-sum-equals-decohered-state"
]

(* non-projective measurement: a SIC-POVM has non-orthogonal Kraus operators, so the
   correct branches M_m rho M_m^dag are genuinely non-diagonal (not just zeroed off a
   projector). They must still be Hermitian, PSD, trace p_m, and resolve the averaged
   state. Guards that the fix is CP-map correct, not merely projector-correct. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["TetrahedronSICPOVM"][
        QuantumState[{{0.6, 0.2 + 0.1 I}, {0.2 - 0.1 I, 0.4}}]]},
        Module[{dms, probs, post},
            dms = Normal[#["DensityMatrix"]] & /@ qm["States"];
            probs = qm["ProbabilitiesList"];
            post = Normal[qm["PostMeasurementState"]["DensityMatrix"]];
            {
                AllTrue[dms, Max[Abs[# - ConjugateTranspose[#]]] < 1.*^-10 &],
                AllTrue[dms, Min[Re @ Eigenvalues[#]] > -1.*^-10 &],
                Max[Abs[(Tr /@ dms) - probs]] < 1.*^-10,
                Max[Abs[Flatten[Total[dms] - post]]] < 1.*^-10
            }
        ]
    ],
    {True, True, True, True},
    TestID -> "BranchStates-POVM-tetrahedron-CP-branches"
]

(* basis independence: an X measurement of a state with coherence in the X eigenbasis
   (here rho is diagonal in Z, so it is NOT diagonal in X). Branch validity is intrinsic
   (Hermiticity, positivity, trace) and holds in any basis representation. *)
VerificationTest[
    With[{dms = Normal[#["DensityMatrix"]] & /@
        QuantumMeasurementOperator["X"][QuantumState[{{3/4, 0}, {0, 1/4}}]]["States"],
        probs = QuantumMeasurementOperator["X"][QuantumState[{{3/4, 0}, {0, 1/4}}]]["ProbabilitiesList"]},
        {AllTrue[dms, HermitianMatrixQ], AllTrue[dms, PositiveSemidefiniteMatrixQ], (Tr /@ dms) == probs}
    ],
    {True, True, True},
    TestID -> "BranchStates-X-basis-intrinsic-invariants"
]

(* two-qubit Z(x)Z on a mixed Bell state: the branches carry the Bell coherence on one
   row before the fix. All four branches must be Hermitian, PSD, and trace to the
   outcome probability. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator[{1, 2}][
        QuantumState[0.7 KroneckerProduct[{1, 0, 0, 1}/Sqrt[2], {1, 0, 0, 1}/Sqrt[2]] + 0.3 IdentityMatrix[4]/4, 2]]},
        Module[{dms, probs},
            dms = Normal[#["DensityMatrix"]] & /@ qm["States"];
            probs = qm["ProbabilitiesList"];
            {
                Length[dms],
                AllTrue[dms, Max[Abs[# - ConjugateTranspose[#]]] < 1.*^-10 &],
                AllTrue[dms, Min[Re @ Eigenvalues[#]] > -1.*^-10 &],
                Max[Abs[(Tr /@ dms) - probs]] < 1.*^-10
            }
        ]
    ],
    {4, True, True, True},
    TestID -> "BranchStates-two-qubit-ZZ-mixed-Bell"
]

(* a pure input carries no bra pointer index; its branches are M_m|psi><psi|M_m^dag
   directly and were always correct. The fix leaves this path untouched. *)
VerificationTest[
    Normal[#["DensityMatrix"]] & /@ QuantumMeasurementOperator["Z"][QuantumState["+"]]["States"],
    {{{1/2, 0}, {0, 0}}, {{0, 0}, {0, 1/2}}},
    TestID -> "BranchStates-pure-ket-unchanged"
]

(* regression anchor for the case that always worked: a state diagonal in the
   measurement eigenbasis. The fix must not disturb it. *)
VerificationTest[
    Normal[#["DensityMatrix"]] & /@
        QuantumMeasurementOperator["Z"][QuantumState[{{3/4, 0}, {0, 1/4}}]]["States"],
    {{{3/4, 0}, {0, 0}}, {{0, 0}, {0, 1/4}}},
    TestID -> "BranchStates-diagonal-input-unchanged"
]

(* the user-facing accessor "StateAssociation" (outcome -> branch) is built from
   "States" and must surface the corrected branches too *)
VerificationTest[
    Normal[#["DensityMatrix"]] & /@ Values @
        QuantumMeasurementOperator["Z"][QuantumState[{{1/2, 1/4}, {1/4, 1/2}}]]["StateAssociation"],
    {{{1/2, 0}, {0, 0}}, {{0, 0}, {0, 1/2}}},
    TestID -> "BranchStates-StateAssociation-accessor-corrected"
]

(* representation independence: the branch is fixed by storage type, not purity. A pure
   state given as a DENSITY MATRIX carries both pointer indices, so it must yield the same
   two branches as the same physical state given as a ket, not Eigendimension^2 of them.
   Before the fix this returned four branches with traces {1/4,1/4,0,0}. *)
VerificationTest[
    Normal[#["DensityMatrix"]] & /@
        QuantumMeasurementOperator["Z"][QuantumState[{{1/2, 1/2}, {1/2, 1/2}}]]["States"],
    {{{1/2, 0}, {0, 0}}, {{0, 0}, {0, 1/2}}},
    TestID -> "BranchStates-pure-state-as-density-matrix"
]

(* a deterministic input has a zero-probability outcome: its branch is the zero operator,
   trace 0, and the outcome drops out of "StateAssociation". Before the fix a pure density
   matrix produced a length mismatch (Thread::tdlen) between branches and probabilities. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{{1, 0}, {0, 0}}]]},
        {
            Normal[#["DensityMatrix"]] & /@ qm["States"],
            (Tr[Normal[#["DensityMatrix"]]] & /@ qm["States"]) == qm["ProbabilitiesList"],
            Length[qm["StateAssociation"]]
        }
    ],
    {{{{1, 0}, {0, 0}}, {{0, 0}, {0, 0}}}, True, 1},
    TestID -> "BranchStates-zero-probability-outcome"
]

(* a genuine single-qudit d=3 measurement (not the d=4 two-qubit tensor): three branches,
   each a valid conditional state, resolving the decohered state. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator[{1}, 3][
        QuantumState[{{1/3, 1/12, 0}, {1/12, 1/3, 1/12}, {0, 1/12, 1/3}}, 3]]},
        Module[{dms = Normal[#["DensityMatrix"]] & /@ qm["States"]},
            {Length[dms], AllTrue[dms, HermitianMatrixQ], AllTrue[dms, PositiveSemidefiniteMatrixQ],
             (Tr /@ dms) == qm["ProbabilitiesList"],
             Total[dms] == Normal[qm["PostMeasurementState"]["DensityMatrix"]]}
        ]
    ],
    {3, True, True, True, True},
    TestID -> "BranchStates-qutrit-single-qudit"
]

(* a rank-1 projector observable on a qutrit has a degenerate zero eigenspace; the eigenbasis
   comes back from an orthogonalized degenerate block. The branches must still be valid
   conditional states with the right traces and resolution. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator[Outer[Times, {1, 1, 0}/Sqrt[2], {1, 1, 0}/Sqrt[2]]][
        QuantumState[N @ {{1/2, 1/6, 1/12}, {1/6, 1/4, 0}, {1/12, 0, 1/4}}, 3]]},
        Module[{dms = Normal[Chop[#["DensityMatrix"]]] & /@ qm["States"], probs = qm["ProbabilitiesList"]},
            {AllTrue[dms, Max[Abs[# - ConjugateTranspose[#]]] < 1.*^-10 &],
             AllTrue[dms, Min[Re @ Eigenvalues[#]] > -1.*^-10 &],
             Max[Abs[(Tr /@ dms) - probs]] < 1.*^-10,
             Max[Abs[Flatten[Total[dms] - Normal[qm["PostMeasurementState"]["DensityMatrix"]]]]] < 1.*^-10}
        ]
    ],
    {True, True, True, True},
    TestID -> "BranchStates-degenerate-eigenspace-projective"
]

(* invariants are necessary but not sufficient: pin the branch values against the CP image
   M_m . rho . ConjugateTranspose[M_m] built independently from the measurement's own
   projectors, for a projective Z on a mixed state. *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{{1/2, 1/4}, {1/4, 1/2}}]],
          rho = {{1/2, 1/4}, {1/4, 1/2}}},
        Sort[Normal[#["DensityMatrix"]] & /@ qm["States"]] ==
            Sort[(# . rho . ConjugateTranspose[#] &) /@ (Normal /@ qm["Projectors"])]
    ],
    True,
    TestID -> "BranchStates-matches-independent-CP-image"
]

(* the four invariants are necessary but do not fix each branch of a non-orthogonal-Kraus
   measurement. Pin the SIC-POVM branches by value against the CP image
   M_m . rho . ConjugateTranspose[M_m] with M_m = Sqrt[E_m] reconstructed from the POVM
   elements: each branch must equal one such image. *)
VerificationTest[
    With[{rho = {{0.6, 0.2 + 0.1 I}, {0.2 - 0.1 I, 0.4}}, qmo = QuantumMeasurementOperator["TetrahedronSICPOVM"]},
        With[{
            states = Normal[Chop[#["DensityMatrix"]]] & /@ qmo[QuantumState[rho]]["States"],
            cp = (# . rho . ConjugateTranspose[#] &) /@ (MatrixPower[N[#], 1/2] & /@ (Normal /@ qmo["POVMElements"]))
        },
            AllTrue[states, s |-> AnyTrue[cp, Max[Abs[Flatten[s - #]]] < 1.*^-9 &]]
        ]
    ],
    True,
    TestID -> "BranchStates-POVM-matches-independent-CP-image"
]

(* partial measurement: measure qubit 1 of a mixed two-qubit state and keep the entangled
   spectator. Each branch lives in the full two-qubit space and equals
   (P_m (x) I) . rho . (P_m (x) I), and the branches resolve the FULL decohered state.
   The spectator is retained here, so the resolution is against the full state, not
   "PostMeasurementState" (which traces the spectator out and lives in a smaller space). *)
VerificationTest[
    With[{rho = 0.7 KroneckerProduct[{1, 0, 0, 1}/Sqrt[2], {1, 0, 0, 1}/Sqrt[2]] + 0.3 IdentityMatrix[4]/4},
        With[{qm = QuantumMeasurementOperator["Z", {1}][QuantumState[rho, 2]],
              p0 = KroneckerProduct[{{1, 0}, {0, 0}}, IdentityMatrix[2]],
              p1 = KroneckerProduct[{{0, 0}, {0, 1}}, IdentityMatrix[2]]},
            With[{dms = Normal[Chop[#["DensityMatrix"]]] & /@ qm["States"], cp = {p0 . rho . p0, p1 . rho . p1}},
                {
                    Dimensions /@ dms,
                    AllTrue[dms, s |-> AnyTrue[cp, Max[Abs[Flatten[s - #]]] < 1.*^-9 &]],
                    Max[Abs[Flatten[Total[dms] - Total[cp]]]] < 1.*^-9,
                    Max[Abs[(Tr /@ dms) - qm["ProbabilitiesList"]]] < 1.*^-10
                }
            ]
        ]
    ],
    {{{4, 4}, {4, 4}}, True, True, True},
    TestID -> "BranchStates-partial-measurement-keeps-spectator"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - operators of a projective measurement"]

(* "Operators" gives, under each outcome label, the operator a projective measurement
   applies for that outcome: the projector P_k onto the orthonormal eigenbasis the
   measurement itself builds. So the operators carry the outcome labels of the
   measurement in its order, sum to the identity, are idempotent, and give its
   branches: P_k.rho.P_k is the k-th post-measurement state, whose trace is the
   probability of outcome k. The observables have a repeated eigenvalue in a rotated
   basis, where Eigensystem returns a non-orthogonal basis of the eigenspace. Exact
   residuals are decided exactly, inexact ones to 10^-8. *)

(* whether the operators of qmo carry the outcome labels of qmo[psi] in order, and the
   entries of Sum_k P_k - 1, P_k.P_k - P_k and P_k.rho.P_k minus the k-th branch, all
   in the computational basis, the one "MatrixRepresentation" writes the operators in;
   a measurement returns its branches in the basis of its observable *)
eigenbasisMeasurementResiduals[qmo_, psi_] := With[
    {ops = Normal[#["MatrixRepresentation"]] & /@ qmo["Operators"], qm = qmo[psi], rho = Normal[psi["Computational"]["DensityMatrix"]]},
    {
        Keys[ops] === Keys[qm["Probabilities"]],
        Flatten[{
            Total[Values[ops]] - IdentityMatrix[Length[rho]],
            (# . # - #) & /@ Values[ops],
            MapThread[#1 . rho . #1 - Normal[#2["Computational"]["DensityMatrix"]] &, {Values[ops], qm["States"]}]
        }]
    }
]

eigenbasisMeasurementExactQ[{labelsQ_, residuals_}] := {labelsQ, AllTrue[residuals, PossibleZeroQ[#, Method -> "ExactAlgebraics"] &]}

VerificationTest[
    With[{f = FourierMatrix[4]},
        eigenbasisMeasurementExactQ @ eigenbasisMeasurementResiduals[
            QuantumMeasurementOperator[QuantumOperator[f . DiagonalMatrix[{1, 1, -1, -1}] . ConjugateTranspose[f], {1, 2}]],
            QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
        ]
    ],
    {True, True},
    TestID -> "Operators-ExactRotatedDegenerate"
]

(* Machine precision: this rotation fails HermitianMatrixQ from roundoff alone, so
   Eigensystem takes its general solver, whose eigenvectors for a repeated eigenvalue
   are not orthogonal. *)
VerificationTest[
    With[{u = BlockRandom[Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {4, 4}]], RandomSeeding -> 35]},
        With[{m = ConjugateTranspose[u] . DiagonalMatrix[N @ {1, 1, -1, -1}] . u},
            With[{r = eigenbasisMeasurementResiduals[QuantumMeasurementOperator[QuantumOperator[m, {1, 2}]], QuantumState[N @ {1, 2 I, -1, 3} / Sqrt[15], {2, 2}]]},
                {HermitianMatrixQ[m], First[r], Max[Abs[Last[r]]]}
            ]
        ]
    ],
    {False, True, _ ? (# < 10^-8 &)},
    SameTest -> MatchQ,
    TestID -> "Operators-MachineRotatedDegenerate"
]

VerificationTest[
    With[{f = FourierMatrix[3]},
        eigenbasisMeasurementExactQ @ eigenbasisMeasurementResiduals[
            QuantumMeasurementOperator[QuantumOperator[f . DiagonalMatrix[{1, 1, 0}] . ConjugateTranspose[f], {1}, 3]],
            QuantumState[{1, 2, 2 I} / 3, 3]
        ]
    ],
    {True, True},
    TestID -> "Operators-ExactRotatedDegenerateQutrit"
]

(* The outcomes come in the order of the sorted eigenvalues -1, 0, 1/2, 2, not the
   order Eigenvectors lists the eigenvectors of diag(2, -1, 1/2, 0) in; each outcome
   label must carry its own projector. *)
VerificationTest[
    eigenbasisMeasurementExactQ @ eigenbasisMeasurementResiduals[
        QuantumMeasurementOperator[QuantumOperator[DiagonalMatrix[{2, -1, 1/2, 0}], {1, 2}]],
        QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
    ],
    {True, True},
    TestID -> "Operators-OutcomeOrder"
]

(* an observable stored in the PauliX basis *)
VerificationTest[
    eigenbasisMeasurementExactQ @ eigenbasisMeasurementResiduals[
        QuantumMeasurementOperator[QuantumOperator[QuantumOperator["ZZ"], QuantumBasis["PauliX", 2]]],
        QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
    ],
    {True, True},
    TestID -> "Operators-NonComputationalBasis"
]

(* a measurement that targets qudit 1 of a two-qudit observable: one operator per
   outcome, acting on both qudits *)
VerificationTest[
    eigenbasisMeasurementExactQ @ eigenbasisMeasurementResiduals[
        QuantumMeasurementOperator[QuantumOperator["ZI"], {1}],
        QuantumState[{1, 2 I, -1, 3} / Sqrt[15], {2, 2}]
    ],
    {True, True},
    TestID -> "Operators-PartialTarget"
]

(* "POVMElements", which QuantumStateEstimate reads, are the operators squared and so,
   for a projective measurement, the projectors themselves, summing to the identity *)
VerificationTest[
    With[{f = FourierMatrix[4]},
        With[{qmo = QuantumMeasurementOperator[QuantumOperator[f . DiagonalMatrix[{1, 1, -1, -1}] . ConjugateTranspose[f], {1, 2}]]},
            AllTrue[
                Flatten[{
                    Total[Normal /@ qmo["POVMElements"]] - IdentityMatrix[4],
                    (Normal /@ qmo["POVMElements"]) - Values[Normal[#["MatrixRepresentation"]] & /@ qmo["Operators"]]
                }],
                PossibleZeroQ[#, Method -> "ExactAlgebraics"] &
            ]
        ]
    ],
    True,
    TestID -> "POVMElements-ExactRotatedDegenerate"
]

EndTestSection[]
