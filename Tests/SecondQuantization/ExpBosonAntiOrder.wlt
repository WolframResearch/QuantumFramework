
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

(* The anti-normal rule set is the image of rulesNO under the automorphism
   a -> a^dagger, a^dagger -> -a, which carries creator-left words to annihilator-left
   ones and sends a^dagger a to -(a a^dagger + 1).  Each expected form below was checked
   against a truncated Fock-space evaluation of both sides before being recorded here. *)


BeginTestSection["ExpBosonOrder - Ordering option"]

(* the default must remain the normal-ordering behaviour that predates the option *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]]] ===
        ExpBosonOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]],
            "Ordering" -> "Normal"],
    True,
    TestID -> "ExpBosonAnti-DefaultIsNormal"
]

(* "Normal" and "AntiNormal" are the only accepted settings; anything else must fail
   loudly rather than fall through to one of them *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]],
        "Ordering" -> "NotAnOrdering"],
    $Failed,
    {ExpBosonOrder::ordering},
    TestID -> "ExpBosonAnti-RejectsUnknownOrdering"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - anti-normal number-operator exponentials"]

(* e^(l a a^dag) = A[e^((1 - e^-l) a a^dag)], the mirror of the Stirling/Touchard form.
   Equivalently (1-c)^-(n+1) with c = 1 - e^-l, which is where the sign flip comes from. *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]],
        "Ordering" -> "AntiNormal"],
    AntiNormalOrdered[E^((1 - E^(-\[Lambda]))*\[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "ExpBosonAnti-NumberOperatorAAdag"
]

(* the same expansion for a^dag a picks up the scalar e^-l, which the Times UpValue
   absorbs into the AntiNormalOrdered argument *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]],
        "Ordering" -> "AntiNormal"],
    AntiNormalOrdered[
        E^(-\[Lambda] + (1 - E^(-\[Lambda]))*\[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "ExpBosonAnti-NumberOperatorAdagA"
]

(* "Series" resolves the marker into annihilator-left monomials a^k (a^dag)^k *)
VerificationTest[
    ExpBosonOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]],
        "Ordering" -> "AntiNormal"]["Series"],
    Inactive[Sum][
        ((1 - E^(-\[Lambda]))^k*GeneralizedPower[NonCommutativeMultiply, \[FormalA], k]**
            GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k])/k!,
        {k, 0, Infinity}],
    TestID -> "ExpBosonAnti-NumberOperatorResolving"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - anti-normal linear and quadratic rules"]

(* Weyl splitting: the Gaussian prefactor flips sign and the factors swap order *)
VerificationTest[
    ExpBosonOrder[Exp[\[Alpha]*SuperDagger[\[FormalA]] + \[Beta]*\[FormalA]],
        "Ordering" -> "AntiNormal"],
    E^(\[FormalA]*\[Beta])**E^(\[Alpha]*SuperDagger[\[FormalA]])/E^((\[Alpha]*\[Beta])/2),
    TestID -> "ExpBosonAnti-WeylSplitting"
]

(* the reordering rule matches the opposite input shape to its normal counterpart *)
VerificationTest[
    ExpBosonOrder[Exp[\[Beta]*SuperDagger[\[FormalA]]]**Exp[\[Alpha]*\[FormalA]],
        "Ordering" -> "AntiNormal"],
    E^(\[FormalA]*\[Alpha])**E^(\[Beta]*SuperDagger[\[FormalA]])/E^(\[Alpha]*\[Beta]),
    TestID -> "ExpBosonAnti-ReorderProduct"
]

(* squeeze from a factored input; the middle factor becomes (1-4ab)^(a a^dag) rather
   than (1-4ab)^-(a^dag a) *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]]**
            Exp[\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]],
        "Ordering" -> "AntiNormal"],
    E^((\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2])/
            (1 - 4*\[Alpha]*\[Beta]))**
        (1 - 4*\[Alpha]*\[Beta])^(\[FormalA]**SuperDagger[\[FormalA]])**
        E^((\[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2])/
            (1 - 4*\[Alpha]*\[Beta]))/Sqrt[1 - 4*\[Alpha]*\[Beta]],
    TestID -> "ExpBosonAnti-SqueezeFactored"
]

(* squeeze from a summed input; 1/2 + a^dag a becomes 1/2 - a a^dag *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
            \[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]],
        "Ordering" -> "AntiNormal"],
    E^((\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*
            Tan[2*Sqrt[\[Alpha]*\[Beta]]])/(2*Sqrt[\[Alpha]*\[Beta]]))**
        Sec[2*Sqrt[\[Alpha]*\[Beta]]]^(1/2 - \[FormalA]**SuperDagger[\[FormalA]])**
        E^((\[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
            Tan[2*Sqrt[\[Alpha]*\[Beta]]])/(2*Sqrt[\[Alpha]*\[Beta]])),
    TestID -> "ExpBosonAnti-SqueezeSummed"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - anti-normal shift rules"]

(* P(a,a^dag) e^(a a) = e^(a a) A[P(a, a^dag - a)]: the exponential migrates left and the
   shift changes sign, with BosonicAntiNormalOrder reducing the shifted polynomial *)
VerificationTest[
    ExpBosonOrder[
        (SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]**\[FormalA])**
            Exp[\[Alpha]*\[FormalA]],
        "Ordering" -> "AntiNormal"],
    E^(\[FormalA]*\[Alpha])**(2*\[Alpha] + \[FormalA]*\[Alpha]^2 +
        \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
        2*\[Alpha]*\[FormalA]**SuperDagger[\[FormalA]] - 2*SuperDagger[\[FormalA]]),
    TestID -> "ExpBosonAnti-ShiftRight"
]

(* e^(b a^dag) P(a,a^dag) = A[P(a - b, a^dag)] e^(b a^dag) *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Beta]*SuperDagger[\[FormalA]]]**
            (SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA]),
        "Ordering" -> "AntiNormal"],
    (-2*\[FormalA] + 2*\[Beta] - 2*\[Beta]*\[FormalA]**SuperDagger[\[FormalA]] +
        GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]**SuperDagger[\[FormalA]] +
        \[Beta]^2*SuperDagger[\[FormalA]])**E^(\[Beta]*SuperDagger[\[FormalA]]),
    TestID -> "ExpBosonAnti-ShiftLeft"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - anti-normal Blasiak Thm 4.4"]

(* The mirrored theorem governs one creator flanked by annihilators, e = L + R - 1:
   e^(l a^L a^dag a^R) = (1 + e l a^e)^(-R/e) A[exp((1 - (1 + e l a^e)^(-1/e)) a a^dag)] *)

(* L=2, R=1  ->  e=2, so the prefactor exponent is -R/e = -1/2 *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]**
            SuperDagger[\[FormalA]]**\[FormalA]],
        "Ordering" -> "AntiNormal"],
    AntiNormalOrdered[
        (1/Sqrt[1 + 2*\[Lambda]*
            GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]])**
        E^((1 - 1/Sqrt[1 + 2*\[Lambda]*
            GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]])**
            \[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "ExpBosonAnti-Thm44-L2R1"
]

(* L=0, R=2  ->  e=1, prefactor exponent -R/e = -2 *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Lambda]*SuperDagger[\[FormalA]]**
            GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]],
        "Ordering" -> "AntiNormal"],
    AntiNormalOrdered[
        (1 + \[FormalA]*\[Lambda])^(-2)**
        E^(((\[FormalA]*\[Lambda])/(1 + \[FormalA]*\[Lambda]))**
            \[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "ExpBosonAnti-Thm44-L0R2"
]

(* L=2, R=0  ->  e=1 and R=0, so the prefactor is the identity *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]**
            SuperDagger[\[FormalA]]],
        "Ordering" -> "AntiNormal"],
    AntiNormalOrdered[
        E^(((\[FormalA]*\[Lambda])/(1 + \[FormalA]*\[Lambda]))**
            \[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "ExpBosonAnti-Thm44-L2R0"
]

EndTestSection[]


BeginTestSection["ExpBosonOrder - anti-normal two mode"]

(* Mode-mixing generators are already mixed in a and a^dagger, so these rules reverse the
   factorisation order rather than producing a strictly annihilator-left form.  Note the
   two modes exchange the sign on log Cosh relative to the normal rule. *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Alpha]*SuperDagger[\[FormalA]]**\[FormalB] +
            \[Beta]*\[FormalA]**SuperDagger[\[FormalB]]],
        "Ordering" -> "AntiNormal"],
    E^((\[Beta]*\[FormalA]**SuperDagger[\[FormalB]]*Tanh[Sqrt[\[Alpha]*\[Beta]]])/
            Sqrt[\[Alpha]*\[Beta]])**
        Cosh[Sqrt[\[Alpha]*\[Beta]]]^(\[FormalA]**SuperDagger[\[FormalA]] -
            \[FormalB]**SuperDagger[\[FormalB]])**
        E^((\[Alpha]*SuperDagger[\[FormalA]]**\[FormalB]*Tanh[Sqrt[\[Alpha]*\[Beta]]])/
            Sqrt[\[Alpha]*\[Beta]]),
    TestID -> "ExpBosonAnti-BeamSplitter"
]

(* the two-mode squeezer does come out strictly annihilator-left *)
VerificationTest[
    ExpBosonOrder[
        Exp[\[Alpha]*\[FormalA]**\[FormalB] +
            \[Beta]*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]],
        "Ordering" -> "AntiNormal"],
    E^((\[Alpha]*\[FormalA]**\[FormalB]*Tan[Sqrt[\[Alpha]*\[Beta]]])/
            Sqrt[\[Alpha]*\[Beta]])**Sec[Sqrt[\[Alpha]*\[Beta]]]*
        Sec[Sqrt[\[Alpha]*\[Beta]]]^(-\[FormalA]**SuperDagger[\[FormalA]] -
            \[FormalB]**SuperDagger[\[FormalB]])**
        E^((\[Beta]*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]*
            Tan[Sqrt[\[Alpha]*\[Beta]]])/Sqrt[\[Alpha]*\[Beta]]),
    TestID -> "ExpBosonAnti-TwoModeSqueezer"
]

EndTestSection[]
