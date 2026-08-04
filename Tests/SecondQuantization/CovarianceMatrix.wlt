BeginTestSection["CovarianceMatrix - vacuum and Fock states"]

(* Vacuum: Serafini convention gives sigma = (1/2) I *)
VerificationTest[
    CovarianceMatrix[FockState[0]],
    {{1/2, 0}, {0, 1/2}},
    TestID -> "CM-Vacuum-Serafini"
]

(* Fock |1>: sigma = (3/2) I  [general |n>: (2n+1)/2 * I] *)
VerificationTest[
    CovarianceMatrix[FockState[1]],
    {{3/2, 0}, {0, 3/2}},
    TestID -> "CM-Fock1-Serafini"
]

(* Fock |2>: sigma = (5/2) I *)
VerificationTest[
    CovarianceMatrix[FockState[2]],
    {{5/2, 0}, {0, 5/2}},
    TestID -> "CM-Fock2-Serafini"
]

EndTestSection[]


BeginTestSection["CovarianceMatrix - mixed states"]

(* Thermal state with nbar=1: sigma = (nbar + 1/2) I = (3/2) I
   Use large Fock space so the truncation error is negligible *)
VerificationTest[
    Chop[CovarianceMatrix[ThermalState[1., 64]] - {{3/2, 0}, {0, 3/2}}],
    {{0, 0}, {0, 0}},
    TestID -> "CM-Thermal-nbar1"
]

(* Thermal state with nbar=0 recovers vacuum *)
VerificationTest[
    Chop[Quiet@CovarianceMatrix[ThermalState[10.^-12, 32]] - {{1/2, 0}, {0, 1/2}}],
    {{0, 0}, {0, 0}},
    TestID -> "CM-Thermal-nbar0-IsVacuum"
]

EndTestSection[]


BeginTestSection["CovarianceMatrix - squeezed vacuum"]

(* Squeezed vacuum with real r: sigma = diag(Exp[-2r]/2, Exp[2r]/2)
   X is squeezed, P is anti-squeezed *)
VerificationTest[
    Chop[
        CovarianceMatrix[SqueezeOperator[0.5, 40] @ FockState[0, 40]] -
        {{Exp[-1.]/2, 0}, {0, Exp[1.]/2}}
    ],
    {{0, 0}, {0, 0}},
    TestID -> "CM-SqueezedVacuum-r05"
]

(* At r=0 squeezed = vacuum *)
VerificationTest[
    Chop[Quiet@CovarianceMatrix[SqueezeOperator[10.^-12, 32] @ FockState[0, 32]] - {{1/2, 0}, {0, 1/2}}],
    {{0, 0}, {0, 0}},
    TestID -> "CM-SqueezedVacuum-r0-IsVacuum"
]

(* Uncertainty principle: det(sigma) >= 1/4 for all physical states *)
VerificationTest[
    Det[CovarianceMatrix[SqueezeOperator[1.5, 64] @ FockState[0, 64]]] >= 1/4,
    True,
    TestID -> "CM-SqueezedVacuum-UncertaintyPrinciple"
]

EndTestSection[]


BeginTestSection["CovarianceMatrix - QuadratureScaling option"]

(* hbar=1/2 convention (current QuadratureOperators default): sigma_vac = (1/4) I *)
VerificationTest[
    CovarianceMatrix[FockState[0], "QuadratureScaling" -> 1/2],
    {{1/4, 0}, {0, 1/4}},
    TestID -> "CM-Vacuum-HalfHbar"
]

(* Simon convention (s=1): sigma_vac = I *)
VerificationTest[
    CovarianceMatrix[FockState[0], "QuadratureScaling" -> 1],
    {{1, 0}, {0, 1}},
    TestID -> "CM-Vacuum-Simon"
]

(* Scaling relation: CM(s) = (s * sqrt(2))^2 * CM(1/sqrt(2)) *)
VerificationTest[
    Chop[CovarianceMatrix[FockState[1], "QuadratureScaling" -> 1] -
         2 * CovarianceMatrix[FockState[1], "QuadratureScaling" -> 1/Sqrt[2]]],
    {{0, 0}, {0, 0}},
    TestID -> "CM-ScalingRelation"
]

EndTestSection[]


BeginTestSection["CovarianceMatrix - multi-mode"]

(* Two-mode vacuum: sigma is 4x4 block-diagonal (1/2)I_4 - modes are uncorrelated *)
VerificationTest[
    CovarianceMatrix[FockState[{0, 0}], {1, 2}],
    {{1/2, 0, 0, 0}, {0, 1/2, 0, 0}, {0, 0, 1/2, 0}, {0, 0, 0, 1/2}},
    TestID -> "CM-TwoMode-Vacuum"
]

(* Two-mode Fock |1,0>: mode 1 has sigma (3/2)I, mode 2 has sigma (1/2)I, no cross terms *)
VerificationTest[
    CovarianceMatrix[FockState[{1, 0}], {1, 2}],
    {{3/2, 0, 0, 0}, {0, 3/2, 0, 0}, {0, 0, 1/2, 0}, {0, 0, 0, 1/2}},
    TestID -> "CM-TwoMode-Fock10"
]

(* Single-mode via explicit order {1} gives same result as default *)
VerificationTest[
    CovarianceMatrix[FockState[0], {1}],
    CovarianceMatrix[FockState[0]],
    TestID -> "CM-SingleMode-ConsistencyCheck"
]

EndTestSection[]
