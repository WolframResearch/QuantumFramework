
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

(* Ported from TestsExpSimplify.nb.  Ladder operators are the formal symbols
   \[FormalA], \[FormalB], \[FormalC]; SuperDagger marks a creator, ** is the
   non-commutative product, and GeneralizedPower[NonCommutativeMultiply, op, n]
   is op ** ... ** op repeated n times. *)


BeginTestSection["ExpBosonOrder - number-operator exponentials"]

(* Exp[lambda a^dagger a] normal-orders to the Stirling/Touchard form
   NormalOrdered[Exp[(E^lambda - 1) a^dagger a]].  The "Series" property re-expresses
   that as an inactive sum over normally-ordered monomials. *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]]],
    NormalOrdered[E^((-1 + E^\[Lambda])*SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedRes1"
]

VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]]]["Series"],
    Inactive[Sum][
    ((-1 + E^\[Lambda])^k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], k])/k!, {k, 0, Infinity}],
    TestID -> "NormalOrderResolving"
]

VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]]],
    E^\[Lambda]*NormalOrdered[
    E^((-1 + E^\[Lambda])*SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedRes2"
]

VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]]]["Series"],
    E^\[Lambda]*Inactive[Sum][
    ((-1 + E^\[Lambda])^k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], k])/k!, {k, 0, Infinity}],
    TestID -> "NormalOrderedResolving2"
]

VerificationTest[
    ExpBosonOrder[Exp[(3/10)*SuperDagger[\[FormalA]]**\[FormalA]]],
    NormalOrdered[E^((-1 + E^(3/10))*SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedResExactNumber"
]

VerificationTest[
    ExpBosonOrder[Exp[0.3*SuperDagger[\[FormalB]]**\[FormalB]]],
    NormalOrdered[E^((Exp[0.3] - 1)*SuperDagger[\[FormalB]]**\[FormalB])],
    TestID -> "NormalOrderOtherVariableAndFloat"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - linear ladder combinations"]

(* A linear combination in a and a^dagger disentangles by the Weyl relation,
   picking up the Gaussian prefactor Exp[alpha beta / 2]. *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*\[FormalA] + \[Beta]*SuperDagger[\[FormalA]]]],
    Exp[(\[Beta]*\[Lambda])/2]*Exp[\[Beta]*SuperDagger[\[FormalA]]]**
    Exp[\[Lambda]*\[FormalA]],
    TestID -> "LinearCombinationLadder"
]

VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*\[FormalA] + \[Beta]*SuperDagger[\[FormalA]]], 
    Assumptions -> \[Beta] == 1/\[Lambda]],
    Sqrt[E]*E^(\[Beta]*SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Lambda]),
    TestID -> "TestAssumptions"
]

VerificationTest[
    ExpBosonOrder[Exp[0.3*\[FormalA] + 0.5*SuperDagger[\[FormalA]]]],
    1.0778841508846315*E^(0.5*SuperDagger[\[FormalA]])**E^(0.3*\[FormalA]),
    TestID -> "LinearCombinationLadderFloat"
]

VerificationTest[
    ExpBosonOrder[Exp[0.8*SuperDagger[\[FormalA]] + 0.5*\[FormalA]]],
    E^((0.8*0.5)/2)*E^(0.8*SuperDagger[\[FormalA]])**E^(0.5*\[FormalA]),
    TestID -> "LinearCombinationLadderFloat2"
]

VerificationTest[
    ExpBosonOrder[DisplacementOperator[\[Alpha], Infinity]],
    DisplacementOperator[\[Alpha], Infinity, "Ordering" -> "Normal"],
    TestID -> "DisplacementOperatorNormal"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - SU(1,1) disentangling and squeezing"]

(* Quadratic (two-photon) generators close on su(1,1), so products and sums of
   Exp[alpha a a] and Exp[beta a^dagger a^dagger] disentangle into the
   Exp[a^dagger a^dagger] . (...)^(a^dagger a) . Exp[a a] normal form. *)
VerificationTest[
    ExpBosonOrder[Exp[k*\[FormalC]**\[FormalC]]**
    Exp[m*SuperDagger[\[FormalC]]**SuperDagger[\[FormalC]]]],
    E^((m*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalC]], 2])/
    (1 - 4*k*m))**(1 - 4*k*m)^(-SuperDagger[\[FormalC]]**\[FormalC])**
    E^((k*GeneralizedPower[NonCommutativeMultiply, \[FormalC], 2])/(1 - 4*k*m))/Sqrt[1 - 4*k*m],
    TestID -> "SU11RuleApply"
]

VerificationTest[
    ExpBosonOrder[Exp[Conjugate[\[Lambda]]*\[FormalC]**\[FormalC]]**
    Exp[\[Lambda]*SuperDagger[\[FormalC]]**SuperDagger[\[FormalC]]]],
    E^((\[Lambda]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalC]], 2])/
    (1 - 4*\[Lambda]*Conjugate[\[Lambda]]))**(1 - 4*\[Lambda]*Conjugate[\[Lambda]])^
    (-SuperDagger[\[FormalC]]**\[FormalC])**
    E^((Conjugate[\[Lambda]]*GeneralizedPower[NonCommutativeMultiply, \[FormalC], 2])/
    (1 - 4*\[Lambda]*Conjugate[\[Lambda]]))/Sqrt[1 - 4*\[Lambda]*Conjugate[\[Lambda]]],
    TestID -> "SU11ApplyConjugate"
]

VerificationTest[
    Evaluate[ExpBosonOrder[Exp[0.5*SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]] - 
    0.2*\[FormalA]**\[FormalA]]]],
    E^(0.4425138959884551*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2])**(0.9103119246332295/
    E^(0.18793592776541562*SuperDagger[\[FormalA]]**\[FormalA]))**
    E^(-0.17700555839538204*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]),
    TestID -> "DisentanglingWithFloats"
]

VerificationTest[
    ExpBosonOrder[
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
    ExpBosonOrder[
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
    ExpBosonOrder[Exp[(1/2)*(r/E^(I*\[Theta]))*\[FormalA]**\[FormalA] - 
    (1/2)*(r*E^(I*\[Theta]))*SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]],
    Assumptions -> r > 0],
    E^((-(1/2))*E^(I*\[Theta])*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]*Tanh[r])**Sech[r]^(SuperDagger[\[FormalA]]**\[FormalA] + 1/2)**
    E^(((1/2)*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "SqueezeOperatorConvention1"
]

VerificationTest[
    ExpBosonOrder[
    Exp[(1/2)*(r*E^(I*\[Theta]))*GeneralizedPower[NonCommutativeMultiply,
    SuperDagger[\[FormalA]], 2] - (1/2)*(r/E^(I*\[Theta]))*
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]], Assumptions -> r > 0],
    E^((1/2)*E^(I*\[Theta])*GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]*Tanh[r])**Sech[r]^(1/2 + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(((-(1/2))*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "SqueezeOperatorConvention2"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - shift rules on exponential-polynomial products"]

(* Conjugating a polynomial by Exp[alpha a] / Exp[beta a^dagger] shifts the ladder
   operator it does not commute with: a^dagger -> a^dagger + alpha and a -> a + beta. *)
VerificationTest[
    ExpBosonOrder[Exp[\[Alpha]*\[FormalA]]**SuperDagger[\[FormalA]]],
    (\[Alpha] + SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "CoherentStateShift"
]

VerificationTest[
    ExpBosonOrder[Exp[\[Alpha]*\[FormalA]]**(SuperDagger[\[FormalA]]**\[FormalA])],
    (\[FormalA]*\[Alpha] + SuperDagger[\[FormalA]]**\[FormalA])**
    E^(\[FormalA]*\[Alpha]),
    TestID -> "NumberOperatorShift"
]

VerificationTest[
    ExpBosonOrder[Exp[\[Alpha]*\[FormalA]]**(SuperDagger[\[FormalA]]**
    SuperDagger[\[FormalA]])],
    (\[Alpha]^2 + GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    2] + 2*\[Alpha]*SuperDagger[\[FormalA]])**E^(\[FormalA]*\[Alpha]),
    TestID -> "TwoPhotonNormalShift"
]

VerificationTest[
    ExpBosonOrder[\[FormalA]**Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(\[FormalA] + \[Beta]),
    TestID -> "SecondShiftAnnihilation"
]

VerificationTest[
    ExpBosonOrder[(SuperDagger[\[FormalA]]**\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(SuperDagger[\[FormalA]]**\[FormalA] + 
    \[Beta]*SuperDagger[\[FormalA]]),
    TestID -> "NumberOperatorShift2"
]

VerificationTest[
    ExpBosonOrder[(SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**
    (2*\[Beta]*SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    \[Beta]^2*SuperDagger[\[FormalA]]),
    TestID -> "MixedExpPolyCase"
]

VerificationTest[
    ExpBosonOrder[(SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA] + y*\[FormalA])**
    Exp[\[Beta]*SuperDagger[\[FormalA]]]],
    E^(\[Beta]*SuperDagger[\[FormalA]])**(y*\[FormalA] + y*\[Beta] + 
    2*\[Beta]*SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalA]]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
    \[Beta]^2*SuperDagger[\[FormalA]]),
    TestID -> "MixedRuleExtraScalarCase"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - generalized normal ordering"]

(* Exp[lambda (a^dagger)^r a (a^dagger)^s] falls under the Blasiak-Penson-Solomon
   theorem, giving a NormalOrdered expression with a non-trivial (a^dagger)-dependent
   prefactor and exponent. *)
VerificationTest[
    ExpBosonOrder[
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
    ExpBosonOrder[
    Exp[2*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]]],
    NormalOrdered[
    (1/Sqrt[1 - 16*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8]])**
    E^((-1 + 1/(1 - 16*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 8])^
    (1/8))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraForms2"
]

VerificationTest[
    ExpBosonOrder[
    Exp[0.3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]],
    NormalOrdered[
    (1/(1 - 1.5*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5])^(1/5))**
    E^((-1 + 1/(1 - 1.5*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5])^
    (1/5))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraFormFloat"
]

VerificationTest[
    ExpBosonOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**\[FormalA]**
    GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]],
    NormalOrdered[(1/(1 - k*SuperDagger[\[FormalA]]))**
    E^((-1 + 1/(1 - k*SuperDagger[\[FormalA]]))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "NormalOrderedExtraForms3"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - series expansions"]

(* "Series"[lambda] returns the coefficient generator -- a DifferenceRoot when the
   recursion has no closed form.  "Series"[lambda, order] truncates at that order, and
   must agree with truncating the inactive sum. *)
VerificationTest[
    ExpBosonOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
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
    ExpBosonOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
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
    ExpBosonOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]**
    \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k, 3],
    Activate[ExpBosonOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
    1]**\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 1]]][
    "Series", k]] /. Infinity -> 3,
    TestID -> "ActivateSumEqualsExplicitTruncation"
]

VerificationTest[
    Collect[ExpBosonOrder[Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 
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


BeginTestSection["ExpBosonOrder - vanishing left ladder block"]

(* The r = 0 case, Exp[lambda a (a^dagger)^s], and its truncated expansions. *)
VerificationTest[
    ExpBosonOrder[
    Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]]],
    NormalOrdered[(1/(-1 + k*SuperDagger[\[FormalA]])^2)**
    E^((-1 + 1/(1 - k*SuperDagger[\[FormalA]]))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "OtherCaseZeroALeft"
]

VerificationTest[
    Collect[ExpBosonOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]]]["Series", k, 2], k],
    1 + k*(GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]**
    \[FormalA] + 2*SuperDagger[\[FormalA]]) +
    k^2*(3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] +
    3*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 3]**\[FormalA] +
    (1/2)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4]**
    GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]),
    TestID -> "ZeroALeftExpansion"
]

VerificationTest[
    ExpBosonOrder[Exp[k*\[FormalA]**GeneralizedPower[NonCommutativeMultiply, 
    SuperDagger[\[FormalA]], 2]]]["Series", k],
    NormalOrdered[Inactive[Sum][k^n*SuperDagger[\[FormalA]]^n*
    ((1 + n)*DifferenceRoot[Function[{\[FormalY], \[FormalN]},
    {(-\[FormalN])*\[FormalY][\[FormalN]] + (2 + \[FormalN])*\[FormalY][3 + \[FormalN]] +
    \[FormalY][2 + \[FormalN]]*(-4 - 3*\[FormalN] - \[FormalA]*SuperDagger[
    \[FormalA]]) + \[FormalY][1 + \[FormalN]]*(2 + 3*\[FormalN] +
    \[FormalA]*SuperDagger[\[FormalA]]) == 0, \[FormalY][-1] == 0,
    \[FormalY][0] == 0, \[FormalY][1] == 1}]][1 + n] -
    DifferenceRoot[Function[{\[FormalY], \[FormalN]},
    {(-1 - \[FormalN])*\[FormalY][\[FormalN]] + (1 + \[FormalN])*
    \[FormalY][3 + \[FormalN]] + \[FormalY][2 + \[FormalN]]*(-3 - 3*\[FormalN] -
    \[FormalA]*SuperDagger[\[FormalA]]) + \[FormalY][1 + \[FormalN]]*
    (3 + 3*\[FormalN] + \[FormalA]*SuperDagger[\[FormalA]]) == 0, \[FormalY][0] == 0,
    \[FormalY][1] == 0, \[FormalY][2] == \[FormalA]*SuperDagger[\[FormalA]]}]][1 + n]),
    {n, 0, Infinity}]],
    TestID -> "ZeroALeftCoefficientForm"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - vanishing right ladder block"]

(* The s = 0 case, Exp[lambda (a^dagger)^r a]. *)
VerificationTest[
    ExpBosonOrder[
    Exp[k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 5]**\[FormalA]]],
    NormalOrdered[
    E^((-1 + 1/(1 - 4*k*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 4])^
    (1/4))**SuperDagger[\[FormalA]]**\[FormalA])],
    TestID -> "OtherCaseZeroARight"
]

VerificationTest[
    Activate[ExpBosonOrder[Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, 
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


BeginTestSection["ExpBosonOrder - two modes"]

(* Cross-mode quadratic generators: the two-mode squeezer, and the beam-splitter-like
   a b^dagger + a^dagger b generator. *)
VerificationTest[
    ExpBosonOrder[Exp[(r*\[FormalA]**\[FormalB])/E^(I*\[Theta]) - 
    r*E^(I*\[Theta])*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]], Assumptions -> r > 0],
    E^((-E^(I*\[Theta]))*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]*Tanh[r])**
    Sech[r]*Sech[r]^(SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB])**
    E^((\[FormalA]**\[FormalB]*Tanh[r])/E^(I*\[Theta])),
    TestID -> "TwoModeSqueezeOperatorNormalOrder"
]

VerificationTest[
    Simplify[ExpBosonOrder[Exp[\[Alpha]*\[FormalA]**SuperDagger[\[FormalB]] + 
    \[Beta]*SuperDagger[\[FormalA]]**\[FormalB]]]] /. {Sqrt[\[Alpha]*\[Beta]] -> \[Gamma],
    1/Sqrt[\[Alpha]*\[Beta]] -> 1/\[Gamma]},
    E^((\[Beta]*SuperDagger[\[FormalA]]**\[FormalB]*Tanh[\[Gamma]])/\[Gamma])**
    Cosh[\[Gamma]]^(-SuperDagger[\[FormalA]]**\[FormalA] + SuperDagger[\[FormalB]]**\[FormalB])**
    E^((\[Alpha]*\[FormalA]**SuperDagger[\[FormalB]]*Tanh[\[Gamma]])/\[Gamma]),
    TestID -> "TwoModeCrossTermNormalOrder"
]

EndTestSection[]
