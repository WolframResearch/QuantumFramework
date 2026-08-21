
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]



BeginTestSection["BosonicExpOrder - independence from parallel sub-kernels"]

(* Ordering a sum must not depend on the parallel subsystem: sub-kernels do
   not load the paclet, and a parallel split tied the answer to kernel
   availability (a stalled auto-launch can consume BosonicExpOrder's
   TimeConstraint and hand back the input unchanged).  The known coherent
   shift is pinned in both kernel states, and the definitions are checked
   free of parallel constructs, which is the machine- and run-order-
   independent invariant.  This section stays first in the file so that a
   single-file TestReport catches a parallel split by value; in a shared
   suite session an earlier file may absorb the first call, and the
   undeclared ParallelMap message still fails the test.  The section owns
   the parallel pool: it closes any running kernels and ends with the pool
   empty. *)
VerificationTest[
    (CloseKernels[]; BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]]]),
    (\[Alpha] + SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "PlusSplitZeroKernels"
]

VerificationTest[
    (LaunchKernels[2];
    With[{ordered = BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]]]},
    CloseKernels[]; ordered]),
    (\[Alpha] + SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "PlusSplitLaunchedKernels"
]

(* The same contract at the source level, at no runtime cost on any machine:
   the ordering heads carry no parallel constructs in their definitions. *)
VerificationTest[
    FreeQ[{DownValues[BosonicNormalOrder], DownValues[BosonicAntinormalOrder]},
    ParallelMap | Parallelize | ParallelTable],
    True,
    TestID -> "OrderingHeadsFreeOfParallelConstructs"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - ordering sums"]

(* Both ordering heads are linear over sums.  A two-mode sum with a declared
   scalar shows what one mode cannot: cross-mode factors commute into
   creators-left order and each mode's [a, a\[Dagger]] = 1 contributes its own
   c-number. *)
VerificationTest[
    BosonicNormalOrder[\[FormalA]**SuperDagger[\[FormalA]] + \[FormalB]**SuperDagger[\[FormalB]] +
    \[FormalA]**SuperDagger[\[FormalB]] + \[Lambda]*SuperDagger[\[FormalA]]**\[FormalB],
    {\[FormalA], SuperDagger[\[FormalA]], \[FormalB], SuperDagger[\[FormalB]]},
    "Scalars" -> {\[Lambda]}],
    2 + SuperDagger[\[FormalA]]**\[FormalA] + \[Lambda]*SuperDagger[\[FormalA]]**\[FormalB] +
    SuperDagger[\[FormalB]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB],
    TestID -> "TwoModeSumNormalOrders"
]

(* Reducing the sum whole collects scalar coefficients across summands. *)
VerificationTest[
    BosonicAntinormalOrder[SuperDagger[\[FormalA]]**\[FormalA] -
    4*\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA],
    {\[FormalA], SuperDagger[\[FormalA]]}, "Scalars" -> {\[Lambda]}],
    -1 + 4*\[Lambda] + (1 - 4*\[Lambda])*\[FormalA]**SuperDagger[\[FormalA]],
    TestID -> "AntinormalSumCollects"
]

(* The invariant itself: ordering a sum equals the sum of the ordered parts,
   for both heads, through degree two and a scalar-weighted term. *)
VerificationTest[
    With[{
        x = \[FormalA]**SuperDagger[\[FormalA]],
        y = GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]**\[FormalA],
        z = \[Lambda]*\[FormalA],
        vars = {\[FormalA], SuperDagger[\[FormalA]]}
    },
        {
            Simplify[BosonicNormalOrder[x + y + z, vars, "Scalars" -> {\[Lambda]}] -
                Total[BosonicNormalOrder[#, vars, "Scalars" -> {\[Lambda]}] & /@ {x, y, z}]],
            Simplify[BosonicAntinormalOrder[x + y + z, vars, "Scalars" -> {\[Lambda]}] -
                Total[BosonicAntinormalOrder[#, vars, "Scalars" -> {\[Lambda]}] & /@ {x, y, z}]]
        }
    ],
    {0, 0},
    TestID -> "OrderingLinearOverSums"
]

(* An unknown method on a sum fails once, not once per summand; Automatic
   routes to the default reduction. *)
VerificationTest[
    BosonicNormalOrder[\[FormalA] + SuperDagger[\[FormalA]],
    {\[FormalA], SuperDagger[\[FormalA]]}, Method -> "NoSuchMethod"],
    $Failed,
    {BosonicNormalOrder::unknownMethod},
    TestID -> "UnknownMethodOnSumFailsOnce"
]

VerificationTest[
    BosonicNormalOrder[\[FormalA]**SuperDagger[\[FormalA]],
    {\[FormalA], SuperDagger[\[FormalA]]}, Method -> Automatic],
    1 + SuperDagger[\[FormalA]]**\[FormalA],
    TestID -> "AutomaticMethodRoutesToDefault"
]

(* A summand the Blasiak formula rejects fails the whole sum: no $Failed may
   survive inside a Plus. *)
VerificationTest[
    BosonicNormalOrder[\[FormalA]**SuperDagger[\[FormalA]] + Sin[\[FormalA]],
    {\[FormalA], SuperDagger[\[FormalA]]}, Method -> "Blasiak"],
    $Failed,
    {Wolfram`QuantumFramework`SecondQuantization`PackageScope`ParseBlasiakMonomial::badfactor},
    TestID -> "BlasiakPartialFailureCollapses"
]

(* The general-degree coherent shift: e^{\[Alpha] a} a\[Dagger]^m normal-orders
   to the binomial shift (\[Alpha] + a\[Dagger])^m e^{\[Alpha] a}, checked
   through degree six. *)
VerificationTest[
    Table[Simplify[BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**
        GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], m]] -
        Sum[Binomial[m, j]*\[Alpha]^(m - j)*
            If[j === 0, 1, GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], j]],
            {j, 0, m}]**E^(\[FormalA]*\[Alpha])], {m, 1, 6}],
    {0, 0, 0, 0, 0, 0},
    TestID -> "GeneralDegreeShiftIsBinomial"
]

(* The Blasiak route agrees with the default reduction on a sum, not only on
   monomials. *)
VerificationTest[
    BosonicNormalOrder[\[FormalA]**SuperDagger[\[FormalA]] +
    SuperDagger[\[FormalA]]**\[FormalA]**SuperDagger[\[FormalA]],
    {\[FormalA], SuperDagger[\[FormalA]]}, Method -> "Blasiak"],
    1 + GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]**\[FormalA] +
    SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]],
    TestID -> "BlasiakSumMatchesDefaultRoute"
]

(* The reducer recovers its own canonical relation, and both orderings are
   idempotent on an already-ordered sum. *)
VerificationTest[
    With[{x = SuperDagger[\[FormalA]]**\[FormalA]**SuperDagger[\[FormalA]] + \[Lambda]*\[FormalA],
        vars = {\[FormalA], SuperDagger[\[FormalA]]}},
        {
            BosonicNormalOrder[Commutator[\[FormalA], SuperDagger[\[FormalA]]], vars],
            Simplify[BosonicNormalOrder[BosonicNormalOrder[x, vars, "Scalars" -> {\[Lambda]}],
                vars, "Scalars" -> {\[Lambda]}] -
                BosonicNormalOrder[x, vars, "Scalars" -> {\[Lambda]}]],
            Simplify[BosonicAntinormalOrder[BosonicAntinormalOrder[x, vars, "Scalars" -> {\[Lambda]}],
                vars, "Scalars" -> {\[Lambda]}] -
                BosonicAntinormalOrder[x, vars, "Scalars" -> {\[Lambda]}]]
        }
    ],
    {1, 0, 0},
    TestID -> "OrderingIdempotentAndCanonical"
]

(* Arbitrary degree: the vacuum moments <0|a^m a\[Dagger]^m|0> = m! pin the
   ladder algebra at every degree up to ten operators. *)
VerificationTest[
    Table[BosonicVEV[GeneralizedPower[NonCommutativeMultiply, \[FormalA], m]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], m],
    {\[FormalA], SuperDagger[\[FormalA]]}], {m, 1, 5}],
    {1, 2, 6, 24, 120},
    TestID -> "VacuumMomentsAreFactorials"
]

(* The beam-splitter generator conserves the total photon number and rotates
   the number difference: one commutator normal-orders to zero, the other to
   the orthogonal cross-mode pair. *)
VerificationTest[
    With[{
        gen = \[FormalA]**SuperDagger[\[FormalB]] + SuperDagger[\[FormalA]]**\[FormalB],
        vars = {\[FormalA], SuperDagger[\[FormalA]], \[FormalB], SuperDagger[\[FormalB]]}
    },
        {
            BosonicNormalOrder[Commutator[
                SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB], gen], vars],
            BosonicNormalOrder[Commutator[
                SuperDagger[\[FormalA]]**\[FormalA] - SuperDagger[\[FormalB]]**\[FormalB], gen], vars]
        }
    ],
    {0, 2*SuperDagger[\[FormalA]]**\[FormalB] - 2*SuperDagger[\[FormalB]]**\[FormalA]},
    TestID -> "NumberConservationNormalOrders"
]

(* Independent numeric check against truncated Fock matrices, as in
   Antinormal.wlt: the realization map is validated on the number operator
   first, then each ordered sum must act identically to its original on a
   sub-block safely below the truncation edge. *)
VerificationTest[
    Module[{size = 12, block = 8, am, adm, toMatrix, pN, pA, vars},
        am = AnnihilationOperator[size]["Matrix"];
        adm = ConjugateTranspose[am];
        vars = {\[FormalA], SuperDagger[\[FormalA]]};
        toMatrix[e_Plus] := Total[toMatrix /@ (List @@ e)];
        toMatrix[e_ /; FreeQ[e, \[FormalA]]] := e IdentityMatrix[size, SparseArray];
        toMatrix[e_] := e /. {
            GeneralizedPower[NonCommutativeMultiply, op_, p_] :>
                MatrixPower[op /. {SuperDagger[\[FormalA]] -> adm, \[FormalA] -> am}, p],
            SuperDagger[\[FormalA]] -> adm,
            \[FormalA] -> am,
            NonCommutativeMultiply -> Dot
        };
        pN = \[FormalA]**\[FormalA]**SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]] +
            \[FormalA]**SuperDagger[\[FormalA]];
        pA = SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]];
        {
            Normal[Diagonal[Take[toMatrix[\[FormalA]**SuperDagger[\[FormalA]]], block, block]]],
            Max[Abs[Take[toMatrix[BosonicNormalOrder[pN, vars]] - toMatrix[pN], block, block]]],
            Max[Abs[Take[toMatrix[BosonicAntinormalOrder[pA, vars]] - toMatrix[pA], block, block]]]
        }
    ],
    {Range[8], 0, 0},
    TestID -> "SumRealizationMatchesTruncatedFock"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - number-operator exponentials"]

(* Exp[lambda a^dagger a] normal-orders to the Stirling/Touchard form
   NormalOrdered[Exp[(E^lambda - 1) a^dagger a]].  The "Series" property re-expresses
   that as an inactive sum over normally-ordered monomials. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]]],
    NormalOrdered[E^((-1 + E^\[Lambda])*SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedRes1"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]]]["Series"],
    Inactive[Sum][
    ((-1 + E^\[Lambda])^k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], k])/k!, {k, 0, Infinity}],
    TestID -> "NormalOrderResolving"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]]],
    E^\[Lambda]*NormalOrdered[
    E^((-1 + E^\[Lambda])*SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedRes2"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]]]["Series"],
    E^\[Lambda]*Inactive[Sum][
    ((-1 + E^\[Lambda])^k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], k])/k!, {k, 0, Infinity}],
    TestID -> "NormalOrderedResolving2"
]

VerificationTest[
    BosonicExpOrder[Exp[(3/10)*SuperDagger[\[FormalA]]**\[FormalA]]],
    NormalOrdered[E^((-1 + E^(3/10))*SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedResExactNumber"
]

VerificationTest[
    BosonicExpOrder[Exp[0.3*SuperDagger[\[FormalB]]**\[FormalB]]],
    NormalOrdered[E^((Exp[0.3] - 1)*SuperDagger[\[FormalB]]**\[FormalB])],
    TestID -> "NormalOrderOtherVariableAndFloat"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - linear ladder combinations"]

(* A linear combination in a and a^dagger disentangles by the Weyl relation,
   picking up the Gaussian prefactor Exp[alpha beta / 2]. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA] + \[Beta]*SuperDagger[\[FormalA]]]],
    Exp[(\[Beta]*\[Lambda])/2]*Exp[\[Beta]*SuperDagger[\[FormalA]]]**
    Exp[\[Lambda]*\[FormalA]],
    TestID -> "LinearCombinationLadder"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA] + \[Beta]*SuperDagger[\[FormalA]]], 
    Assumptions -> \[Beta] == 1/\[Lambda]],
    Sqrt[E]*E^(\[Beta]*SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Lambda]),
    TestID -> "TestAssumptions"
]

VerificationTest[
    BosonicExpOrder[Exp[0.3*\[FormalA] + 0.5*SuperDagger[\[FormalA]]]],
    1.0778841508846315*E^(0.5*SuperDagger[\[FormalA]])**E^(0.3*\[FormalA]),
    TestID -> "LinearCombinationLadderFloat"
]

VerificationTest[
    BosonicExpOrder[Exp[0.8*SuperDagger[\[FormalA]] + 0.5*\[FormalA]]],
    E^((0.8*0.5)/2)*E^(0.8*SuperDagger[\[FormalA]])**E^(0.5*\[FormalA]),
    TestID -> "LinearCombinationLadderFloat2"
]

VerificationTest[
    BosonicExpOrder[DisplacementOperator[\[Alpha], Infinity]],
    DisplacementOperator[\[Alpha], Infinity, "Ordering" -> "Normal"],
    TestID -> "DisplacementOperatorNormal"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - SU(1,1) disentangling and squeezing"]

(* Quadratic (two-photon) generators close on su(1,1), so products and sums of
   Exp[alpha a a] and Exp[beta a^dagger a^dagger] disentangle into the
   Exp[a^dagger a^dagger] . (...)^(a^dagger a) . Exp[a a] normal form. *)
VerificationTest[
    BosonicExpOrder[Exp[k*\[FormalC]**\[FormalC]]**
    Exp[m*SuperDagger[\[FormalC]]**SuperDagger[\[FormalC]]]],
    E^((m*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalC]], 2])/
    (1 - 4*k*m))**(1 - 4*k*m)^(-SuperDagger[\[FormalC]]**\[FormalC])**
    E^((k*GeneralizedPower[NonCommutativeMultiply, \[FormalC], 2])/(1 - 4*k*m))/Sqrt[1 - 4*k*m],
    TestID -> "SU11RuleApply"
]

VerificationTest[
    BosonicExpOrder[Exp[Conjugate[\[Lambda]]*\[FormalC]**\[FormalC]]**
    Exp[\[Lambda]*SuperDagger[\[FormalC]]**SuperDagger[\[FormalC]]]],
    E^((\[Lambda]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalC]], 2])/
    (1 - 4*\[Lambda]*Conjugate[\[Lambda]]))**(1 - 4*\[Lambda]*Conjugate[\[Lambda]])^
    (-SuperDagger[\[FormalC]]**\[FormalC])**
    E^((Conjugate[\[Lambda]]*GeneralizedPower[NonCommutativeMultiply, \[FormalC], 2])/
    (1 - 4*\[Lambda]*Conjugate[\[Lambda]]))/Sqrt[1 - 4*\[Lambda]*Conjugate[\[Lambda]]],
    TestID -> "SU11ApplyConjugate"
]

(* The machine-float coefficients vary by an ulp across machines, so the
   comparison rounds every Real to twelve digits before matching. *)
VerificationTest[
    Evaluate[BosonicExpOrder[Exp[0.5*SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]] -
    0.2*\[FormalA]**\[FormalA]]]],
    E^(0.4425138959884551*GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 2])**(0.9103119246332295/
    E^(0.18793592776541562*SuperDagger[\[FormalA]]**\[FormalA]))**
    E^(-0.17700555839538204*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]),
    SameTest -> Function[{a, b},
        (a /. x_Real :> Round[x, 1.*^-12]) === (b /. x_Real :> Round[x, 1.*^-12])],
    TestID -> "DisentanglingWithFloats"
]

(* The same disentangling at exact rationals closes machine-independently,
   keeping the su(1,1) structure the float form collapses. *)
VerificationTest[
    BosonicExpOrder[Exp[(1/2)*SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]] -
    (1/5)*\[FormalA]**\[FormalA]]],
    E^((Sqrt[5/2]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
    Tanh[Sqrt[2/5]])/2)**Sech[Sqrt[2/5]]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(-((GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[Sqrt[2/5]])/Sqrt[10])),
    TestID -> "DisentanglingExactRationals"
]

VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
    \[Beta]*\[FormalA]**\[FormalA]]],
    E^((k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
    Tan[2*Sqrt[(-k)*\[Beta]]])/(2*Sqrt[(-k)*\[Beta]]))**
    Sec[2*Sqrt[(-k)*\[Beta]]]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(-((\[Beta]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*
    Tan[2*Sqrt[(-k)*\[Beta]]])/(2*Sqrt[(-k)*\[Beta]]))),
    TestID -> "DisentanglingSymbolicPassingGP"
]

VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
    \[Beta]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]]],
    E^((k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
    Tan[2*Sqrt[(-k)*\[Beta]]])/(2*Sqrt[(-k)*\[Beta]]))**
    Sec[2*Sqrt[(-k)*\[Beta]]]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(-((\[Beta]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*
    Tan[2*Sqrt[(-k)*\[Beta]]])/(2*Sqrt[(-k)*\[Beta]]))),
    TestID -> "DisentanglingSymbolicPassingGP2"
]

VerificationTest[
    BosonicExpOrder[Exp[(1/2)*(r/E^(I*\[Theta]))*\[FormalA]**\[FormalA] - 
    (1/2)*(r*E^(I*\[Theta]))*SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]],
    Assumptions -> r > 0],
    E^((-(1/2))*E^(I*\[Theta])*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]*Tanh[r])**Sech[r]^(SuperDagger[\[FormalA]]**\[FormalA] + 1/2)**
    E^(((1/2)*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "SqueezeOperatorConvention1"
]

VerificationTest[
    BosonicExpOrder[
    Exp[(1/2)*(r*E^(I*\[Theta]))*GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 2] - (1/2)*(r/E^(I*\[Theta]))*
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]], Assumptions -> r > 0],
    E^((1/2)*E^(I*\[Theta])*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]*Tanh[r])**Sech[r]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(((-(1/2))*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "SqueezeOperatorConvention2"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - shift rules on exponential-polynomial products"]

(* Conjugating a polynomial by Exp[alpha a] / Exp[beta a^dagger] shifts the ladder
   operator it does not commute with: a^dagger -> a^dagger + alpha and a -> a + beta. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]]],
    (\[Alpha] + SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "CoherentStateShift"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**(SuperDagger[\[FormalA]]**\[FormalA])],
    (\[FormalA]*\[Alpha] + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(\[FormalA]*\[Alpha]),
    TestID -> "NumberOperatorShift"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**(SuperDagger[\[FormalA]]**
    SuperDagger[\[FormalA]])],
    (\[Alpha]^2 + GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    2] + 2*\[Alpha]*SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "TwoPhotonNormalShift"
]

VerificationTest[
    BosonicExpOrder[\[FormalA]**Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(\[FormalA] + \[Beta]),
    TestID -> "SecondShiftAnnihilation"
]

VerificationTest[
    BosonicExpOrder[(SuperDagger[\[FormalA]]**\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(SuperDagger[\[FormalA]]**\[FormalA] + 
    \[Beta]*SuperDagger[\[FormalA]]),
    TestID -> "NumberOperatorShift2"
]

VerificationTest[
    BosonicExpOrder[(SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**
    (2*\[Beta]*SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    \[Beta]^2*SuperDagger[\[FormalA]]),
    TestID -> "MixedExpPolyCase"
]

VerificationTest[
    BosonicExpOrder[(SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA] + y*\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(y*\[FormalA] + y*\[Beta] + 
    2*\[Beta]*SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    \[Beta]^2*SuperDagger[\[FormalA]]),
    TestID -> "MixedRuleExtraScalarCase"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - generalized normal ordering"]

(* Exp[lambda (a^dagger)^r a (a^dagger)^s] falls under the Blasiak-Penson-Solomon
   theorem, giving a NormalOrdered expression with a non-trivial (a^dagger)-dependent
   prefactor and exponent. *)
VerificationTest[
    BosonicExpOrder[
    Exp[\[Kappa]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3]]],
    NormalOrdered[
    (1/(1 - 7*\[Kappa]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 7])^
    (3/7))**
    E^((-1 + 1/(1 - 7*\[Kappa]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]],
    7])^(1/7))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraForms1"
]

VerificationTest[
    BosonicExpOrder[
    Exp[2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]]],
    NormalOrdered[
    (1/Sqrt[1 - 16*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8]])**
    E^((-1 + 1/(1 - 16*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8])^
    (1/8))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraForms2"
]

VerificationTest[
    BosonicExpOrder[
    Exp[0.3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]],
    NormalOrdered[
    (1/(1 - 1.5*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5])^(1/5))**
    E^((-1 + 1/(1 - 1.5*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5])^
    (1/5))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraFormFloat"
]

VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]],
    NormalOrdered[(1/(1 - k*SuperDagger[\[FormalA]]))**
    E^((-1 + 1/(1 - k*SuperDagger[\[FormalA]]))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraForms3"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - series expansions"]

(* "Series"[lambda] returns the coefficient generator -- a DifferenceRoot when the
   recursion has no closed form.  "Series"[lambda, order] truncates at that order, and
   must agree with truncating the inactive sum. *)
VerificationTest[
    BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k],
    NormalOrdered[Inactive[Sum][k^n*SuperDagger[\[FormalA]]^n*
    DifferenceRoot[Function[{\[FormalY], \[FormalN]},
    {(-\[FormalN])*\[FormalY][\[FormalN]] + (2 + \[FormalN])*\[FormalY][3 + \[FormalN]] +
    \[FormalY][2 + \[FormalN]]*(-4 - 3*\[FormalN] - \[FormalA]*SuperDagger[\[FormalA]]) +
    \[FormalY][1 + \[FormalN]]*(2 + 3*\[FormalN] + \[FormalA]*SuperDagger[\[FormalA]]) ==
    0, \[FormalY][-1] == 0, \[FormalY][0] == 0, \[FormalY][1] == 1}]][1 + n],
    {n, 0, Infinity}]],
    TestID -> "ShowCoefficientGenerator"
]

VerificationTest[
    BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k, 3],
    1 + k^2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] + 
    k^3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3] +
    k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]**\[FormalA] +
    2*k^2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3]**\[FormalA] +
    3*k^3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]**\[FormalA] +
    (1/2)*k^2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    (3/2)*k^3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    (1/6)*k^3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 6]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 3] + k*SuperDagger[\[FormalA]],
    TestID -> "SeriesTruncation"
]

VerificationTest[
    BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k, 3],
    Activate[BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    1]**\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k]] /. Infinity -> 3,
    TestID -> "ActivateSumEqualsExplicitTruncation"
]

VerificationTest[
    Collect[BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    2]**\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k, 3], k],
    1 + k*(GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] + 
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3]**\[FormalA]) +
    k^2*((3/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4] +
    (5/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA] +
    (1/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 6]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]) +
    k^3*((5/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 6] +
    (11/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 7]**\[FormalA] +
    2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    (1/6)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 9]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 3]),
    TestID -> "SeriesTruncation2"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - vanishing left ladder block"]

(* The r = 0 case, Exp[lambda a (a^dagger)^s], and its truncated expansions. *)
VerificationTest[
    BosonicExpOrder[
    Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]]],
    NormalOrdered[(1/(-1 + k*SuperDagger[\[FormalA]])^2)**
    E^((-1 + 1/(1 - k*SuperDagger[\[FormalA]]))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "OtherCaseZeroALeft"
]

VerificationTest[
    Collect[BosonicExpOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]]]["Series", k, 2], k],
    1 + k*(GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]**
    \[FormalA] + 2*SuperDagger[\[FormalA]]) +
    k^2*(3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] +
    3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3]**\[FormalA] +
    (1/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]),
    TestID -> "ZeroALeftExpansion"
]

(* The coefficient generator's DifferenceRoot shape is not unique: equivalent
   recurrences come back depending on the kernel session's history.  Assert
   the mathematical content instead: activating and truncating the inactive
   sum reproduces the explicit order-3 truncation. *)
VerificationTest[
    Activate[BosonicExpOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 2]]]["Series", k]] /. Infinity -> 3,
    BosonicExpOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 2]]]["Series", k, 3],
    TestID -> "ZeroALeftSeriesMatchesTruncation"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - vanishing right ladder block"]

(* The s = 0 case, Exp[lambda (a^dagger)^r a]. *)
VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]]],
    NormalOrdered[
    E^((-1 + 1/(1 - 4*k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4])^
    (1/4))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "OtherCaseZeroARight"
]

VerificationTest[
    Activate[BosonicExpOrder[Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 5]**\[FormalA]]]["Series", \[Lambda]]] /. Infinity -> 3,
    1 + \[Lambda]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    5]**\[FormalA] + (5/2)*\[Lambda]^2*GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 9]**\[FormalA] + (1/2)*\[Lambda]^2*
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 10]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    (15/2)*\[Lambda]^3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 13]**
    \[FormalA] + (5/2)*\[Lambda]^3*GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 14]**GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    (1/6)*\[Lambda]^3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 15]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 3],
    TestID -> "ZeroARightExpansionFromActivate"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - two modes"]

(* Cross-mode quadratic generators: the two-mode squeezer, and the beam-splitter-like
   a b^dagger + a^dagger b generator. *)
VerificationTest[
    BosonicExpOrder[Exp[(r*\[FormalA]**\[FormalB])/E^(I*\[Theta]) - 
    r*E^(I*\[Theta])*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]], Assumptions -> r > 0],
    E^((-E^(I*\[Theta]))*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]*Tanh[r])**
    Sech[r]*Sech[r]^(SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB])**
    E^((\[FormalA]**\[FormalB]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "TwoModeSqueezeOperatorNormalOrder"
]

VerificationTest[
    Simplify[BosonicExpOrder[Exp[\[Alpha]*\[FormalA]**SuperDagger[\[FormalB]] + 
    \[Beta]*SuperDagger[\[FormalA]]**\[FormalB]]]] /. {Sqrt[\[Alpha]*\[Beta]] -> \[Gamma],
    1/Sqrt[\[Alpha]*\[Beta]] -> 1/\[Gamma]},
    E^((\[Beta]*SuperDagger[\[FormalA]]**\[FormalB]*Tanh[\[Gamma]])/\[Gamma])**
    Cosh[\[Gamma]]^(-SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB])**
    E^((\[Alpha]*\[FormalA]**SuperDagger[\[FormalB]]*Tanh[\[Gamma]])/\[Gamma]),
    TestID -> "TwoModeCrossTermNormalOrder"
]

EndTestSection[]
