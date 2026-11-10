
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

(* The summation index of a "Series" result is owned by the paclet and prints as a
   plain letter; rename it to a local symbol so the expected values below do not
   pin an internal name. *)
normalizeIndex[expr_, s_] :=
    expr /. Inactive[Sum][b_, {i_, lo_, hi_}] :> Inactive[Sum][b /. i -> s, {s, lo, hi}]

$vA = {{SuperDagger[\[FormalA]]}, {\[FormalA]}};
$vB = {{SuperDagger[\[FormalB]]}, {\[FormalB]}};
$vC = {{SuperDagger[\[FormalC]]}, {\[FormalC]}};
$vAB = {{SuperDagger[\[FormalA]], SuperDagger[\[FormalB]]}, {\[FormalA], \[FormalB]}};


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
    (CloseKernels[]; BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]], $vA]),
    (\[Alpha] + SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "PlusSplitZeroKernels"
]

VerificationTest[
    (LaunchKernels[2];
    With[{ordered = BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]], $vA]},
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
    $vAB],
    2 + SuperDagger[\[FormalA]]**\[FormalA] + \[Lambda]*SuperDagger[\[FormalA]]**\[FormalB] +
    SuperDagger[\[FormalB]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB],
    TestID -> "TwoModeSumNormalOrders"
]

(* Reducing the sum whole collects scalar coefficients across summands. *)
VerificationTest[
    BosonicAntinormalOrder[SuperDagger[\[FormalA]]**\[FormalA] -
    4*\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA],
    $vA],
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
        vars = $vA
    },
        {
            Simplify[BosonicNormalOrder[x + y + z, vars] -
                Total[BosonicNormalOrder[#, vars] & /@ {x, y, z}]],
            Simplify[BosonicAntinormalOrder[x + y + z, vars] -
                Total[BosonicAntinormalOrder[#, vars] & /@ {x, y, z}]]
        }
    ],
    {0, 0},
    TestID -> "OrderingLinearOverSums"
]

(* An unknown method on a sum fails once, not once per summand; Automatic
   routes to the default reduction. *)
VerificationTest[
    BosonicNormalOrder[\[FormalA] + SuperDagger[\[FormalA]],
    $vA, Method -> "NoSuchMethod"],
    $Failed,
    {BosonicNormalOrder::unknownMethod},
    TestID -> "UnknownMethodOnSumFailsOnce"
]

VerificationTest[
    BosonicNormalOrder[\[FormalA]**SuperDagger[\[FormalA]],
    $vA, Method -> Automatic],
    1 + SuperDagger[\[FormalA]]**\[FormalA],
    TestID -> "AutomaticMethodRoutesToDefault"
]

(* A summand the Blasiak formula rejects fails the whole sum: no $Failed may
   survive inside a Plus. *)
VerificationTest[
    BosonicNormalOrder[\[FormalA]**SuperDagger[\[FormalA]] + Sin[\[FormalA]],
    $vA, Method -> "Blasiak"],
    $Failed,
    {BosonicNormalOrder::ladpoly},
    TestID -> "BlasiakPartialFailureCollapses"
]

(* The general-degree coherent shift: e^{\[Alpha] a} a\[Dagger]^m normal-orders
   to the binomial shift (\[Alpha] + a\[Dagger])^m e^{\[Alpha] a}, checked
   through degree six. *)
VerificationTest[
    Table[Simplify[BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**
        GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], m], $vA] -
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
    $vA, Method -> "Blasiak"],
    1 + GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]**\[FormalA] +
    SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]],
    TestID -> "BlasiakSumMatchesDefaultRoute"
]

(* The reducer recovers its own canonical relation, and both orderings are
   idempotent on an already-ordered sum. *)
VerificationTest[
    With[{x = SuperDagger[\[FormalA]]**\[FormalA]**SuperDagger[\[FormalA]] + \[Lambda]*\[FormalA],
        vars = $vA},
        {
            BosonicNormalOrder[Commutator[\[FormalA], SuperDagger[\[FormalA]]], vars],
            Simplify[BosonicNormalOrder[BosonicNormalOrder[x, vars],
                vars] -
                BosonicNormalOrder[x, vars]],
            Simplify[BosonicAntinormalOrder[BosonicAntinormalOrder[x, vars],
                vars] -
                BosonicAntinormalOrder[x, vars]]
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
    $vA], {m, 1, 5}],
    {1, 2, 6, 24, 120},
    TestID -> "VacuumMomentsAreFactorials"
]

(* The beam-splitter generator conserves the total photon number and rotates
   the number difference: one commutator normal-orders to zero, the other to
   the orthogonal cross-mode pair. *)
VerificationTest[
    With[{
        gen = \[FormalA]**SuperDagger[\[FormalB]] + SuperDagger[\[FormalA]]**\[FormalB],
        vars = $vAB
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
        vars = $vA;
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
    BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]], $vA],
    NormalOrdered[E^((-1 + E^\[Lambda])*SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedRes1"
]

VerificationTest[
    normalizeIndex[
        BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]], $vA]["Series"], k],
    Inactive[Sum][
    ((-1 + E^\[Lambda])^k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], k])/k!, {k, 0, Infinity}],
    TestID -> "NormalOrderResolving"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]], $vA],
    E^\[Lambda]*NormalOrdered[
    E^((-1 + E^\[Lambda])*SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedRes2"
]

VerificationTest[
    normalizeIndex[
        BosonicExpOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]], $vA]["Series"], k],
    E^\[Lambda]*Inactive[Sum][
    ((-1 + E^\[Lambda])^k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], k])/k!, {k, 0, Infinity}],
    TestID -> "NormalOrderedResolving2"
]

VerificationTest[
    BosonicExpOrder[Exp[(3/10)*SuperDagger[\[FormalA]]**\[FormalA]], $vA],
    NormalOrdered[E^((-1 + E^(3/10))*SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedResExactNumber"
]

VerificationTest[
    BosonicExpOrder[Exp[0.3*SuperDagger[\[FormalB]]**\[FormalB]], $vB],
    NormalOrdered[E^((Exp[0.3] - 1)*SuperDagger[\[FormalB]]**\[FormalB]), $vB],
    TestID -> "NormalOrderOtherVariableAndFloat"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - linear ladder combinations"]

(* A linear combination in a and a^dagger disentangles by the Weyl relation,
   picking up the Gaussian prefactor Exp[alpha beta / 2]. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA] + \[Beta]*SuperDagger[\[FormalA]]], $vA],
    Exp[(\[Beta]*\[Lambda])/2]*Exp[\[Beta]*SuperDagger[\[FormalA]]]**
    Exp[\[Lambda]*\[FormalA]],
    TestID -> "LinearCombinationLadder"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA] + \[Beta]*SuperDagger[\[FormalA]]], $vA, 
    Assumptions -> \[Beta] == 1/\[Lambda]],
    Sqrt[E]*E^(\[Beta]*SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Lambda]),
    TestID -> "TestAssumptions"
]

VerificationTest[
    BosonicExpOrder[Exp[0.3*\[FormalA] + 0.5*SuperDagger[\[FormalA]]], $vA],
    1.0778841508846315*E^(0.5*SuperDagger[\[FormalA]])**E^(0.3*\[FormalA]),
    TestID -> "LinearCombinationLadderFloat"
]

VerificationTest[
    BosonicExpOrder[Exp[0.8*SuperDagger[\[FormalA]] + 0.5*\[FormalA]], $vA],
    E^((0.8*0.5)/2)*E^(0.8*SuperDagger[\[FormalA]])**E^(0.5*\[FormalA]),
    TestID -> "LinearCombinationLadderFloat2"
]

(* The symbolic displacement is written in the declared variables, FieldVariables[] by
   default, and normal-orders to its own "Normal" form in the same variables. *)
VerificationTest[
    {BosonicExpOrder[DisplacementOperator[\[Alpha], Infinity], FieldVariables[]] ===
        DisplacementOperator[\[Alpha], Infinity, "Ordering" -> "Normal"],
     BosonicExpOrder[DisplacementOperator[\[Alpha], Infinity, {\[FormalB]}], {\[FormalB]}] ===
        DisplacementOperator[\[Alpha], Infinity, {\[FormalB]}, "Ordering" -> "Normal"],
     And @@ Table[
        DisplacementOperator[\[Alpha], Infinity, {\[FormalB]}, "Ordering" -> o] === (DisplacementOperator[\[Alpha], Infinity, "Ordering" -> o] /. \[FormalA] -> \[FormalB]) &&
        SqueezeOperator[\[Xi], Infinity, {\[FormalB]}, "Ordering" -> o] === (SqueezeOperator[\[Xi], Infinity, "Ordering" -> o] /. \[FormalA] -> \[FormalB]),
        {o, {"Weak", "Normal", "Antinormal"}}]},
    {True, True, True},
    TestID -> "DisplacementOperatorNormal"
]

(* A displacement or squeeze acts on one mode; declaring two is refused. *)
VerificationTest[
    {DisplacementOperator[\[Alpha], Infinity, {\[FormalA], \[FormalB]}], SqueezeOperator[\[Xi], Infinity, {\[FormalA], \[FormalB]}]},
    {$Failed, $Failed},
    {DisplacementOperator::ladone, SqueezeOperator::ladone},
    TestID -> "SymbolicOperatorOneMode"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - SU(1,1) disentangling and squeezing"]

(* Quadratic (two-photon) generators close on su(1,1), so products and sums of
   Exp[alpha a a] and Exp[beta a^dagger a^dagger] disentangle into the
   Exp[a^dagger a^dagger] . (...)^(a^dagger a) . Exp[a a] normal form. *)
VerificationTest[
    BosonicExpOrder[Exp[k*\[FormalC]**\[FormalC]]**
    Exp[m*SuperDagger[\[FormalC]]**SuperDagger[\[FormalC]]], $vC],
    E^((m*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalC]], 2])/
    (1 - 4*k*m))**(1 - 4*k*m)^(-SuperDagger[\[FormalC]]**\[FormalC])**
    E^((k*GeneralizedPower[NonCommutativeMultiply, \[FormalC], 2])/(1 - 4*k*m))/Sqrt[1 - 4*k*m],
    TestID -> "SU11RuleApply"
]

VerificationTest[
    BosonicExpOrder[Exp[Conjugate[\[Lambda]]*\[FormalC]**\[FormalC]]**
    Exp[\[Lambda]*SuperDagger[\[FormalC]]**SuperDagger[\[FormalC]]], $vC],
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
    0.2*\[FormalA]**\[FormalA]], $vA]],
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
    (1/5)*\[FormalA]**\[FormalA]], $vA],
    E^((Sqrt[5/2]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
    Tanh[Sqrt[2/5]])/2)**Sech[Sqrt[2/5]]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(-((GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[Sqrt[2/5]])/Sqrt[10])),
    TestID -> "DisentanglingExactRationals"
]

VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
    \[Beta]*\[FormalA]**\[FormalA]], $vA],
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
    \[Beta]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]], $vA],
    E^((k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
    Tan[2*Sqrt[(-k)*\[Beta]]])/(2*Sqrt[(-k)*\[Beta]]))**
    Sec[2*Sqrt[(-k)*\[Beta]]]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(-((\[Beta]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*
    Tan[2*Sqrt[(-k)*\[Beta]]])/(2*Sqrt[(-k)*\[Beta]]))),
    TestID -> "DisentanglingSymbolicPassingGP2"
]

VerificationTest[
    BosonicExpOrder[Exp[(1/2)*(r/E^(I*\[Theta]))*\[FormalA]**\[FormalA] - 
    (1/2)*(r*E^(I*\[Theta]))*SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]], $vA,
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
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]], $vA, Assumptions -> r > 0],
    E^((1/2)*E^(I*\[Theta])*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]*Tanh[r])**Sech[r]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(((-(1/2))*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "SqueezeOperatorConvention2"
]

(* The polar form S(r e^{i th}) carries no Arg or Abs: its normal form is written in r and th, and
   BosonicExpOrder of its weak form gives the same product. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        {
            SqueezeOperator[{r, th}, Infinity, {a}, "Ordering" -> "Normal"],
            Simplify[BosonicExpOrder[SqueezeOperator[{r, th}, Infinity, {a}], $vA, Assumptions -> r > 0 && th \[Element] Reals] -
                SqueezeOperator[{r, th}, Infinity, {a}, "Ordering" -> "Normal"], r > 0 && th \[Element] Reals]
        }
    ],
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        {Exp[-E^(I th) Tanh[r]/2 ad ** ad] ** Sech[r]^(ad ** a + 1/2) ** Exp[E^(-I th) Tanh[r]/2 a ** a], 0}
    ],
    TestID -> "SqueezeOperatorPolar"
]

(* The polar form also builds the truncated matrix: equal to the one from r e^{i th} in every ordering, and
   with symbolic r and th written in them alone, its normal-order entries the closed forms of BosonicMatrixElement. *)
VerificationTest[
    {
        Normal[SqueezeOperator[{1/2, Pi/3}, 4, "Ordering" -> #]["Matrix"]] == Normal[SqueezeOperator[E^(I Pi/3)/2, 4, "Ordering" -> #]["Matrix"]] & /@
            {"Normal", "Weak", "Antinormal"},
        FreeQ[Normal[SqueezeOperator[{r, th}, 4, "Ordering" -> #]["Matrix"]], Abs | Arg | Conjugate] & /@ {"Normal", "Weak", "Antinormal"},
        Simplify[Normal[SqueezeOperator[{r, th}, 4]["Matrix"]][[{1, 2, 3}, {1, 2, 1}]] -
            Outer[BosonicMatrixElement[{#1, #2}, SqueezeOperator[{r, th}]] &, {0, 1, 2}, {0, 1, 0}], r > 0 && th \[Element] Reals]
    },
    {{True, True, True}, {True, True, True}, ConstantArray[0, {3, 3}]},
    TestID -> "SqueezeOperatorPolarTruncated"
]

(* Any single-mode exponent of degree two, e^{t(q a^2 + p a^dagger^2 + w N + mu a^dagger + nu a)}, disentangles
   in either order: the general one, the parabolic (a + a^dagger)^2 + mu a^dagger + nu a where 4pq = w^2 and
   the coefficients are the l -> 0 limits, and a number term beside linear ones. An ordered product's
   symbol (a^dagger -> x, a -> y) is the product of its factors' symbols, the middle one read through
   e^{g N} = :e^{(e^g - 1) N}: and e^{g a a^dagger} = antinormal e^{(1 - e^-g) a a^dagger}; to order t^3 it
   matches sum_k t^k/k! times the symbol of the ordered H^k. *)
VerificationTest[
    Module[{p, q, w, mu, nu, t, x, y, sym, fac},
        With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
            sym[e_] := e /. GeneralizedPower[NonCommutativeMultiply, o_, k_] :> o^k /. NonCommutativeMultiply -> Times /. {ad -> x, a -> y};
            fac[n_, mid_][f_] := With[{e = Expand[Replace[f, {Times[c__, Power[b_, z_ /; ! FreeQ[z, n]]] :> Log[Times[c]] + z Log[b], Power[b_, z_] :> z Log[b]}]]},
                Exp[sym[e - Coefficient[e, n] n] + mid[Coefficient[e, n]] x y]];
            Table[
                Simplify[Normal[Series[BosonicExpOrder[Exp[t h], $vA, "Ordering" -> o[[1]]] /. c_. NonCommutativeMultiply[f__] :> c Times @@ (fac[o[[2]], o[[3]]] /@ {f}), {t, 0, 3}]] -
                    Sum[t^k/k! sym[o[[4]][If[k == 0, 1, NonCommutativeMultiply @@ ConstantArray[h, k]], $vA]], {k, 0, 3}]],
                {h, {q a ** a + p ad ** ad + w ad ** a + mu ad + nu a, (a + ad) ** (a + ad) + mu ad + nu a, w ad ** a + mu ad + nu a}},
                {o, {{"Normal", ad ** a, E^# - 1 &, BosonicNormalOrder}, {"Antinormal", a ** ad, 1 - E^-# &, BosonicAntinormalOrder}}}
            ]
        ]
    ],
    ConstantArray[0, {3, 2}],
    TestID -> "DisentanglingQuadraticPlusLinear"
]

(* The root D^(-1/2) stays continuous where D crosses the negative real axis: for e^{-itH} with
   H = N + c a^dagger^2 + conj(c) a^2, c = 1/5 + i/10, D winds once per period 2 Pi/Omega, Omega = 2/Sqrt[5],
   and the vacuum element has one limit from both sides of the first crossing t = Pi Sqrt[5]/2. *)
VerificationTest[
    Module[{t},
        With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
            With[{v0 = BosonicMatrixElement[{0, 0}, BosonicExpOrder[Exp[-I t (ad ** a + (1/5 + I/10) ad ** ad + (1/5 - I/10) a ** a)], $vA], $vA]},
                Simplify[Limit[v0, t -> Pi Sqrt[5]/2, Direction -> "FromBelow"] - Limit[v0, t -> Pi Sqrt[5]/2, Direction -> "FromAbove"]]
            ]
        ]
    ],
    0,
    TestID -> "DisentanglingRootContinuousInTime"
]

(* With a symbolic frequency om the winding phase cannot be chosen, so the root stays principal and the
   vacuum element is regular at om = 0, where it is the pure squeezer's Sech[2 t/3]^(1/2). *)
VerificationTest[
    With[{b = \[FormalA], bd = SuperDagger[\[FormalA]]},
        With[{g = BosonicMatrixElement[{0, 0}, BosonicExpOrder[Exp[-I t (om bd ** b + 1/3 bd ** bd + 1/3 b ** b)], $vA], $vA]},
            FullSimplify[{g /. om -> 0, Limit[g, om -> 0]} - Sqrt[Sech[2 t/3]], t \[Element] Reals]
        ]
    ],
    {0, 0},
    TestID -> "DisentanglingRootRegularAtZeroFrequency"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - shift rules on exponential-polynomial products"]

(* Conjugating a polynomial by Exp[alpha a] / Exp[beta a^dagger] shifts the ladder
   operator it does not commute with: a^dagger -> a^dagger + alpha and a -> a + beta. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]], $vA],
    (\[Alpha] + SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "CoherentStateShift"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**(SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    (\[FormalA]*\[Alpha] + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(\[FormalA]*\[Alpha]),
    TestID -> "NumberOperatorShift"
]

VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*\[FormalA]]**(SuperDagger[\[FormalA]]**
    SuperDagger[\[FormalA]]), $vA],
    (\[Alpha]^2 + GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    2] + 2*\[Alpha]*SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "TwoPhotonNormalShift"
]

VerificationTest[
    BosonicExpOrder[\[FormalA]**Exp[\[Beta]*SuperDagger[\[FormalA]]], $vA],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(\[FormalA] + \[Beta]),
    TestID -> "SecondShiftAnnihilation"
]

VerificationTest[
    BosonicExpOrder[(SuperDagger[\[FormalA]]**\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]], $vA],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(SuperDagger[\[FormalA]]**\[FormalA] + 
    \[Beta]*SuperDagger[\[FormalA]]),
    TestID -> "NumberOperatorShift2"
]

VerificationTest[
    BosonicExpOrder[(SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]], $vA],
    E^(\[Beta]*SuperDagger[\[FormalA]])**
    (2*\[Beta]*SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    \[Beta]^2*SuperDagger[\[FormalA]]),
    TestID -> "MixedExpPolyCase"
]

VerificationTest[
    BosonicExpOrder[(SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA] + y*\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]], $vA],
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
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3]], $vA],
    NormalOrdered[
    (1/(1 - 7*\[Kappa]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 7])^
    (3/7))**
    E^((-1 + 1/(1 - 7*\[Kappa]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]],
    7])^(1/7))**SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedExtraForms1"
]

VerificationTest[
    BosonicExpOrder[
    Exp[2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]], $vA],
    NormalOrdered[
    (1/Sqrt[1 - 16*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8]])**
    E^((-1 + 1/(1 - 16*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8])^
    (1/8))**SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedExtraForms2"
]

VerificationTest[
    BosonicExpOrder[
    Exp[0.3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA],
    NormalOrdered[
    (1/(1 - 1.5*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5])^(1/5))**
    E^((-1 + 1/(1 - 1.5*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5])^
    (1/5))**SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedExtraFormFloat"
]

VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA],
    NormalOrdered[(1/(1 - k*SuperDagger[\[FormalA]]))**
    E^((-1 + 1/(1 - k*SuperDagger[\[FormalA]]))**SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "NormalOrderedExtraForms3"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - series expansions"]

(* "Series"[lambda] returns the coefficient generator -- a DifferenceRoot when the
   recursion has no closed form.  "Series"[lambda, order] truncates at that order, and
   must agree with truncating the inactive sum. *)
VerificationTest[
    normalizeIndex[
        BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
        \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA][
        "Series", k], n],
    NormalOrdered[Inactive[Sum][k^n*SuperDagger[\[FormalA]]^n*
    DifferenceRoot[Function[{\[FormalY], \[FormalN]},
    {(-\[FormalN])*\[FormalY][\[FormalN]] + (2 + \[FormalN])*\[FormalY][3 + \[FormalN]] +
    \[FormalY][2 + \[FormalN]]*(-4 - 3*\[FormalN] - \[FormalA]*SuperDagger[\[FormalA]]) +
    \[FormalY][1 + \[FormalN]]*(2 + 3*\[FormalN] + \[FormalA]*SuperDagger[\[FormalA]]) ==
    0, \[FormalY][-1] == 0, \[FormalY][0] == 0, \[FormalY][1] == 1}]][1 + n],
    {n, 0, Infinity}], $vA],
    TestID -> "ShowCoefficientGenerator"
]

VerificationTest[
    BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA][
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
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA][
    "Series", k, 3],
    Activate[BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    1]**\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA][
    "Series", k]] /. Infinity -> 3,
    TestID -> "ActivateSumEqualsExplicitTruncation"
]

VerificationTest[
    Collect[BosonicExpOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    2]**\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]], $vA][
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
    Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]], $vA],
    NormalOrdered[(1/(-1 + k*SuperDagger[\[FormalA]])^2)**
    E^((-1 + 1/(1 - k*SuperDagger[\[FormalA]]))**SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "OtherCaseZeroALeft"
]

VerificationTest[
    Collect[BosonicExpOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]], $vA]["Series", k, 2], k],
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
    SuperDagger[\[FormalA]], 2]], $vA]["Series", k]] /. Infinity -> 3,
    BosonicExpOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 2]], $vA]["Series", k, 3],
    TestID -> "ZeroALeftSeriesMatchesTruncation"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - vanishing right ladder block"]

(* The s = 0 case, Exp[lambda (a^dagger)^r a]. *)
VerificationTest[
    BosonicExpOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]], $vA],
    NormalOrdered[
    E^((-1 + 1/(1 - 4*k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4])^
    (1/4))**SuperDagger[\[FormalA]]**\[FormalA]), $vA],
    TestID -> "OtherCaseZeroARight"
]

VerificationTest[
    Activate[BosonicExpOrder[Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 5]**\[FormalA]], $vA]["Series", \[Lambda]]] /. Infinity -> 3,
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
    r*E^(I*\[Theta])*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]], $vAB, Assumptions -> r > 0],
    Sech[r]*(E^((-E^(I*\[Theta]))*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]*Tanh[r])**
    Sech[r]^(SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB])**
    E^((\[FormalA]**\[FormalB]*Tanh[r])/E^(I*\[Theta]))),
    TestID -> "TwoModeSqueezeOperatorNormalOrder"
]

(* The disentangled squeezer must stay the same operator once numbers go in, before or after
   ordering, in both orderings.  A scalar left between two factors splits the product into an
   orderless Times of operator halves, which machine-number substitution re-sorts; that was
   0.19 off the exact exponential.  Each result is realized as a truncated operator and
   compared with MatrixExp of the generator on the n_a, n_b <= 2 block.  The anti-normal
   form sums over states above the block, so each mode keeps 25 levels.  The reference is
   built in machine numbers: exact 625 x 625 products took 28 s of the test's 31. *)
VerificationTest[
    Module[{s = 25, pars = {r -> 0.3, \[Theta] -> 0.7}, am, bm, idx, gen, exact},
        am = N @ KroneckerProduct[Normal[AnnihilationOperator[s]["Matrix"]], IdentityMatrix[s]];
        bm = N @ KroneckerProduct[IdentityMatrix[s], Normal[AnnihilationOperator[s]["Matrix"]]];
        idx = Flatten[Table[s i + j + 1, {i, 0, 2}, {j, 0, 2}]];
        gen = r E^(-I \[Theta]) \[FormalA]**\[FormalB] -
            r E^(I \[Theta]) SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]];
        exact = MatrixExp[(r E^(-I \[Theta]) /. pars) am . bm -
            (r E^(I \[Theta]) /. pars) ConjugateTranspose[am] . ConjugateTranspose[bm]][[idx, idx]];
        Max[Abs[Normal[ToBosonicOperator[#, {\[FormalA], \[FormalB]}, s]["Matrix"]][[idx, idx]] - exact]] < 10^-12 & /@ {
            BosonicExpOrder[Exp[gen], $vAB, Assumptions -> r > 0] /. pars,
            BosonicExpOrder[Exp[gen], $vAB, Assumptions -> r > 0, "Ordering" -> "Antinormal"] /. pars,
            BosonicExpOrder[Exp[gen /. pars], $vAB],
            BosonicExpOrder[Exp[gen /. pars], $vAB, "Ordering" -> "Antinormal"]
        }
    ],
    {True, True, True, True},
    TestID -> "TwoModeSqueezerSurvivesNumericSubstitution"
]

VerificationTest[
    Simplify[BosonicExpOrder[Exp[\[Alpha]*\[FormalA]**SuperDagger[\[FormalB]] + 
    \[Beta]*SuperDagger[\[FormalA]]**\[FormalB]], $vAB]] /. {Sqrt[\[Alpha]*\[Beta]] -> \[Gamma],
    1/Sqrt[\[Alpha]*\[Beta]] -> 1/\[Gamma]},
    E^((\[Beta]*SuperDagger[\[FormalA]]**\[FormalB]*Tanh[\[Gamma]])/\[Gamma])**
    Cosh[\[Gamma]]^(-SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB])**
    E^((\[Alpha]*\[FormalA]**SuperDagger[\[FormalB]]*Tanh[\[Gamma]])/\[Gamma]),
    TestID -> "TwoModeCrossTermNormalOrder"
]

(* An exponential of one mode is not a polynomial in another, so the shift rules must
   leave it alone; once they treated E^(b^dagger)**E^b as a mode-a scalar and fused it
   into a single exponential with no Weyl factor.  Checked on the n_a, n_b <= 2 block
   against exact single-mode products: truncated MatrixExp of a ladder operator is exact
   there. *)
VerificationTest[
    Module[{s = 10, big = 40, am, idx, ord},
        am = Normal[AnnihilationOperator[big]["Matrix"]];
        idx = Flatten[Table[s i + j + 1, {i, 0, 2}, {j, 0, 2}]];
        ord = BosonicExpOrder[Exp[0.2 \[FormalA]]**Exp[0.3 SuperDagger[\[FormalB]]]**
            Exp[-0.25 \[FormalB]]**Exp[0.15 SuperDagger[\[FormalA]]], $vAB];
        Max[Abs[Normal[ToBosonicOperator[ord, {\[FormalA], \[FormalB]}, s]["Matrix"]][[idx, idx]] -
            KroneckerProduct[
                (MatrixExp[0.2 am] . MatrixExp[0.15 ConjugateTranspose[am]])[[;; 3, ;; 3]],
                (MatrixExp[0.3 ConjugateTranspose[am]] . MatrixExp[-0.25 am])[[;; 3, ;; 3]]
            ]]] < 10^-12
    ],
    True,
    TestID -> "TwoModeInterleavedExponentialsStayApart"
]

(* Normal order across modes: an exponential of one mode commutes with every exponential of
   another, so creators of any mode move left of annihilators of any other at no cost.  The
   ranks are 1 for creator-only factors, 3 for annihilator-only, 2 otherwise. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Alpha] \[FormalA]]**Exp[\[Beta] SuperDagger[\[FormalB]]]**Exp[\[Gamma] \[FormalB]], $vAB],
    Exp[\[Beta] SuperDagger[\[FormalB]]]**Exp[\[Alpha] \[FormalA]]**Exp[\[Gamma] \[FormalB]],
    TestID -> "TwoModeCreatorsMoveLeft"
]

VerificationTest[
    Module[{s = 10, big = 40, am, idx, in, ord, rank},
        am = Normal[AnnihilationOperator[big]["Matrix"]];
        idx = Flatten[Table[s i + j + 1, {i, 0, 2}, {j, 0, 2}]];
        rank = Which[FreeQ[#, SuperDagger], 3, FreeQ[# /. SuperDagger[_] -> 1, \[FormalA] | \[FormalB]], 1, True, 2] &;
        in = Exp[0.2 \[FormalA]]**Exp[0.3 SuperDagger[\[FormalB]]]**Exp[-0.25 \[FormalB]]**Exp[0.15 SuperDagger[\[FormalA]]];
        ord = BosonicExpOrder[in, $vAB];
        {
            rank /@ List @@ Replace[ord, c_ m_NonCommutativeMultiply /; FreeQ[c, \[FormalA] | \[FormalB]] :> m],
            Max[Abs[Normal[ToBosonicOperator[ord, {\[FormalA], \[FormalB]}, s]["Matrix"]][[idx, idx]] -
                KroneckerProduct[
                    (MatrixExp[0.2 am] . MatrixExp[0.15 ConjugateTranspose[am]])[[;; 3, ;; 3]],
                    (MatrixExp[0.3 ConjugateTranspose[am]] . MatrixExp[-0.25 am])[[;; 3, ;; 3]]
                ]]] < 10^-12
        }
    ],
    {{1, 1, 3, 3}, True},
    TestID -> "TwoModeInterleavedProductNormalOrders"
]

(* The scalar of e^(-i t H), H = Na + 2 Nb + (ad bd + a b)/2, against its su(1,1) closed form: the vacuum
   amplitude in normal order, and in antinormal order the same with -3/2 for 3/2 in front of Sin, from
   a -> ad, ad -> -a.  Det[D] circles the origin with frequency 2 Sqrt[2], so the principal root alone flips
   sign on t in [1.11, 3.33], [5.55, 7.77], [9.99, 12.2]: t = 2, 6, 10 sit inside, t = 1, 4 outside. *)
VerificationTest[
    Module[{h},
        h = SuperDagger[\[FormalA]]**\[FormalA] + 2 SuperDagger[\[FormalB]]**\[FormalB] +
            (SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]] + \[FormalA]**\[FormalB])/2;
        Max[Abs[N[#1 /. (NormalOrdered | AntinormalOrdered)[_, _] -> 1 /. NonCommutativeMultiply -> Times /.
            Thread[Flatten[$vAB] -> 0] /. t -> {1, 2, 4, 6, 10}, 30] -
            N[E^(3 I t/2)/(Cos[Sqrt[2] t] + #2 3 I/(2 Sqrt[2]) Sin[Sqrt[2] t]) /. t -> {1, 2, 4, 6, 10}, 30]]] < 10^-25 & @@@ {
            {BosonicExpOrder[Exp[-I t h], $vAB], 1},
            {BosonicExpOrder[Exp[-I t h], $vAB, "Ordering" -> "Antinormal"], -1}
        }
    ],
    {True, True},
    TestID -> "TwoModeRootContinuousInTime"
]

(* A beam splitter on top of the squeezer, ordered at machine times: each vacuum amplitude matches the
   truncated evolution, where the principal root alone gave the opposite sign at t = 2, 6.5, 11. *)
VerificationTest[
    Module[{s = 20, h, am, bm, vals, vecs},
        h = SuperDagger[\[FormalA]]**\[FormalA] + 2 SuperDagger[\[FormalB]]**\[FormalB] +
            (1/2 + I/5) SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]] + (1/2 - I/5) \[FormalA]**\[FormalB] +
            (SuperDagger[\[FormalA]]**\[FormalB] + SuperDagger[\[FormalB]]**\[FormalA])/3;
        am = N @ KroneckerProduct[Normal[AnnihilationOperator[s]["Matrix"]], IdentityMatrix[s]];
        bm = N @ KroneckerProduct[IdentityMatrix[s], Normal[AnnihilationOperator[s]["Matrix"]]];
        {vals, vecs} = Eigensystem[h /. NonCommutativeMultiply -> Dot /.
            {SuperDagger[\[FormalA]] -> ConjugateTranspose[am], SuperDagger[\[FormalB]] -> ConjugateTranspose[bm], \[FormalA] -> am, \[FormalB] -> bm}];
        Max @ Table[Abs[(BosonicExpOrder[Exp[-I t h], $vAB] /. NormalOrdered[_, _] -> 1 /. NonCommutativeMultiply -> Times /.
            Thread[Flatten[$vAB] -> 0]) - Abs[vecs[[All, 1]]]^2 . Exp[-I t vals]], {t, {1., 2., 6.5, 11.}}] < 10^-6
    ],
    True,
    TestID -> "TwoModeRootContinuousAtMachineTime"
]

(* Coupled two-mode generators with linear terms through the multimode rule, hopping in normal order and
   squeezing plus a constant in antinormal order: the symbol of the ordered form (a^dagger -> x, a -> y, the
   number middle merged before it is read) equals the series of e^(t G) to t^2. *)
VerificationTest[
    Module[{ad = SuperDagger[\[FormalA]], bd = SuperDagger[\[FormalB]], a = \[FormalA], b = \[FormalB], deq, sym, series},
        deq = # //. {ad -> x1, bd -> x2, a -> y1, b -> y2, GeneralizedPower[NonCommutativeMultiply, u_, k_] :> u^k, NonCommutativeMultiply -> Times} &;
        sym[o_, ord_] := With[{nums = If[ord === "Normal", {ad ** a, bd ** b}, {a ** ad, b ** bd}]},
            deq[(o /. {(NormalOrdered | AntinormalOrdered)[f_, _] :> deq[f],
                    Power[base_, e_] /; ! FreeQ[e, Alternatives @@ nums] && FreeQ[Expand[e] /. Thread[nums -> 0], a | b] :> Inactive[Exp][Log[base] e]}) //.
                Inactive[Exp][u_] Inactive[Exp][w_] :> Inactive[Exp][u + w] /.
                Inactive[Exp][e_] :> With[{c = Coefficient[Expand[e], #] & /@ nums},
                    Exp[If[ord === "Normal", E^c - 1, 1 - E^-c] . {x1 y1, x2 y2} + (Expand[e] /. Thread[nums -> 0])]]]];
        series[g_, ord_] := 1 + Sum[t^j/j! deq[If[ord === "Normal", BosonicNormalOrder, BosonicAntinormalOrder][NonCommutativeMultiply @@ ConstantArray[g, j], $vAB]], {j, 2}];
        Simplify[Normal[Series[sym[BosonicExpOrder[Exp[t #1], $vAB, "Ordering" -> #2], #2] - series[#1, #2], {t, 0, 2}]]] & @@@ {
            {ad ** a/3 + (ad ** b + bd ** a)/7 + bd ** b/2 + ad/5 + b/6, "Normal"},
            {ad ** bd/2 + a ** b/7 + ad/4 + b/9 + 3/7, "Antinormal"}}
    ],
    {0, 0},
    TestID -> "TwoModeCoupledQuadraticAndLinearSeries"
]

(* A displaced two-mode squeezer, e^(r (a^dagger b^dagger - a b) + g a^dagger - conj(g) a), ordered with linear
   terms on a coupled pair, realized on 16 levels per mode against MatrixExp of the generator on n_a, n_b <= 2. *)
VerificationTest[
    Module[{s = 16, am, bm, idx, gen},
        am = N @ KroneckerProduct[Normal[AnnihilationOperator[s]["Matrix"]], IdentityMatrix[s]];
        bm = N @ KroneckerProduct[IdentityMatrix[s], Normal[AnnihilationOperator[s]["Matrix"]]];
        idx = Flatten[Table[s i + j + 1, {i, 0, 2}, {j, 0, 2}]];
        gen = 0.3 (ConjugateTranspose[am] . ConjugateTranspose[bm] - am . bm) + (0.2 + 0.1 I) ConjugateTranspose[am] - (0.2 - 0.1 I) am;
        Max[Abs[Normal[ToBosonicOperator[
            BosonicExpOrder[Exp[r (SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]] - \[FormalA]**\[FormalB]) + g SuperDagger[\[FormalA]] - Conjugate[g] \[FormalA]],
                $vAB, Assumptions -> r > 0] /. {r -> 0.3, g -> 0.2 + 0.1 I}, {\[FormalA], \[FormalB]}, s]["Matrix"]][[idx, idx]] - MatrixExp[gen][[idx, idx]]]] < 10^-10
    ],
    True,
    TestID -> "TwoModeDisplacedSqueezerMatchesMatrixExp"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - products of three or more exponentials"]

(* The rules repeat until nothing changes, so a reordering that exposes a new adjacent pair is
   followed through; a single pass left E^(gamma a) to the left of E^(beta a^dagger). *)
VerificationTest[
    BosonicExpOrder[Exp[\[Gamma] \[FormalA]]**Exp[\[Alpha] \[FormalA]]**Exp[\[Beta] SuperDagger[\[FormalA]]], $vA],
    E^(\[Beta] (\[Alpha] + \[Gamma])) Exp[\[Beta] SuperDagger[\[FormalA]]]**Exp[(\[Alpha] + \[Gamma]) \[FormalA]],
    TestID -> "ThreeFactorProductNormalOrders"
]

(* Exponentials of one mode that are both creator-only or both annihilator-only commute, so they
   merge exactly; across modes, and with anything mixed, they stay apart. *)
VerificationTest[
    {
        BosonicExpOrder[Exp[\[Alpha] \[FormalA]]**Exp[\[Beta] SuperDagger[\[FormalA]]]**Exp[\[Gamma] \[FormalA]], $vA],
        BosonicExpOrder[Exp[\[Alpha] \[FormalA]]**Exp[\[Gamma] \[FormalA]], $vA],
        BosonicExpOrder[Exp[\[Beta] SuperDagger[\[FormalB]]]**Exp[\[Delta] SuperDagger[\[FormalA]]], $vAB]
    },
    {
        E^(\[Alpha] \[Beta]) Exp[\[Beta] SuperDagger[\[FormalA]]]**Exp[(\[Alpha] + \[Gamma]) \[FormalA]],
        Exp[(\[Alpha] + \[Gamma]) \[FormalA]],
        Exp[\[Beta] SuperDagger[\[FormalB]]]**Exp[\[Delta] SuperDagger[\[FormalA]]]
    },
    TestID -> "SameModeCommutingExponentialsMerge"
]

(* A product of n displacements needs of order n^2 reorderings; it must run to the end and
   collapse to one creator and one annihilator exponential. *)
VerificationTest[
    Module[{size = 60, am, in, ord},
        am = Normal[AnnihilationOperator[size]["Matrix"]];
        in = NonCommutativeMultiply @@ Table[Exp[0.1 j SuperDagger[\[FormalA]] - 0.05 j \[FormalA]], {j, 8}];
        ord = BosonicExpOrder[in, $vA];
        {
            Length[Replace[ord, c_ m_NonCommutativeMultiply /; FreeQ[c, \[FormalA]] :> m]],
            Max[Abs[Normal[ToBosonicOperator[ord, {\[FormalA]}, size]["Matrix"]][[;; 12, ;; 12]] -
                (Dot @@ Table[MatrixExp[0.1 j ConjugateTranspose[am] - 0.05 j am], {j, 8}])[[;; 12, ;; 12]]]] < 10^-10
        }
    ],
    {2, True},
    TestID -> "LongDisplacementChainCollapses"
]

(* Exact check, independent of the exponential rules: expand input and output to t^4 and
   normal-order both with BosonicNormalOrder. *)
VerificationTest[
    Module[{vars = $vAB,
            pars = {\[Alpha], \[Beta], \[Gamma], \[Delta]}, n = 4, t, nf},
        nf[e_] := Module[{x = e /. Thread[pars -> t pars]},
            x = x /. Power[E, g_] /; ! FreeQ[g, Alternatives @@ Flatten[vars]] :>
                Sum[If[k == 0, 1, NonCommutativeMultiply @@ ConstantArray[g, k]]/k!, {k, 0, n}];
            x = x /. Power[E, c_] /; FreeQ[c, Alternatives @@ Flatten[vars]] && ! FreeQ[c, t] :> Normal[Series[E^c, {t, 0, n}]];
            x = Expand[NonCommutativeExpand[x]] /. t^k_Integer /; k > n -> 0;
            Expand[BosonicNormalOrder[x, vars]] /. t^k_Integer /; k > n -> 0];
        Simplify[Expand[nf[#] - nf[BosonicExpOrder[#, vars]]]] & /@ {
            Exp[\[Alpha] \[FormalA]]**Exp[\[Beta] SuperDagger[\[FormalA]]]**Exp[\[Gamma] \[FormalA]],
            Exp[\[Alpha] \[FormalA]]**Exp[\[Beta] SuperDagger[\[FormalB]]]**Exp[\[Gamma] \[FormalB]]**Exp[\[Delta] SuperDagger[\[FormalA]]]
        }
    ],
    {0, 0},
    TestID -> "OrderingKeepsTheOperatorExactly"
]

(* the displacement group law D(alpha) D(beta) = e^((alpha beta* - alpha* beta)/2) D(alpha + beta),
   read off as <1| D(alpha) D(beta) |0> *)
VerificationTest[
    FullSimplify[
        BosonicMatrixElement[{1, 0},
            Exp[\[Alpha] SuperDagger[\[FormalA]] - Conjugate[\[Alpha]] \[FormalA]]**
                Exp[\[Beta] SuperDagger[\[FormalA]] - Conjugate[\[Beta]] \[FormalA]], $vA] -
        Exp[-(\[Alpha] Conjugate[\[Alpha]] + \[Beta] Conjugate[\[Beta]])/2 - Conjugate[\[Alpha]] \[Beta]] (\[Alpha] + \[Beta])
    ],
    0,
    TestID -> "DisplacementCompositionMatrixElement"
]

(* Repeating must not re-enter a marker: the number exponential gives one NormalOrdered, not a
   tower, and ordering an ordered product changes nothing. *)
VerificationTest[
    With[{d = Exp[\[Alpha] SuperDagger[\[FormalA]] - Conjugate[\[Alpha]] \[FormalA]]**
            Exp[\[Beta] SuperDagger[\[FormalA]] - Conjugate[\[Beta]] \[FormalA]]},
        {
            Count[BosonicExpOrder[Exp[\[Lambda] SuperDagger[\[FormalA]]**\[FormalA]], $vA], _NormalOrdered, {0, Infinity}],
            BosonicExpOrder[BosonicExpOrder[d, $vA], $vA] === BosonicExpOrder[d, $vA]
        }
    ],
    {1, True},
    TestID -> "RepeatedOrderingIsStableAndIdempotent"
]

(* The number operator rotates the ladder operators, e^{-sN} a e^{sN} = e^s a, so an
   exponential of a alone moves right past e^{sN} as e^{beta a} e^{sN} = e^{sN} e^{beta e^s a},
   and one of a^dagger moves left. Fock elements are compared with the hand-written
   <m|e^{beta a}|n> = beta^(n-m) Sqrt[n!/m!]/(n-m)! and the diagonal e^{sn}, for m, n <= 3. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        With[{low = If[#3 >= #2, #1^(#3 - #2) Sqrt[#3!/#2!]/(#3 - #2)!, 0] &},
            {Table[BosonicMatrixElement[{m, n}, BosonicExpOrder[Exp[be a] ** Exp[s ad ** a], $vA], $vA] - low[be, m, n] E^(s n), {m, 0, 3}, {n, 0, 3}],
             Table[BosonicMatrixElement[{m, n}, BosonicExpOrder[Exp[s ad ** a] ** Exp[be ad], $vA], $vA] - E^(s m) low[be, n, m], {m, 0, 3}, {n, 0, 3}]} // Simplify
        ]
    ],
    ConstantArray[0, {2, 4, 4}],
    TestID -> "LadderExponentialPastNumberExponential"
]

(* With the ladder exponentials moved past the number factor, D(-alpha) e^{-itN} is in
   normal order and its coherent element is read directly:
   <beta|D(-alpha) e^{-itN}|beta> = e^{-|alpha|^2/2 - alpha betabar + (e^{-it} - 1)|beta|^2 + alphabar e^{-it} beta}. *)
VerificationTest[
    With[{a = \[FormalA], ad = SuperDagger[\[FormalA]]},
        Simplify[(BosonicMatrixElement[{be, be}, Exp[-al ad + Conjugate[al] a] ** Exp[-I t ad ** a], $vA, "Basis" -> "Coherent"] /.
            Abs[x_]^2 :> x Conjugate[x]) -
            Exp[-al Conjugate[al]/2 - al Conjugate[be] + (E^(-I t) - 1) be Conjugate[be] + Conjugate[al] E^(-I t) be]]
    ],
    0,
    TestID -> "DisplacementTimesNumberExponentialCoherentElement"
]

EndTestSection[]
