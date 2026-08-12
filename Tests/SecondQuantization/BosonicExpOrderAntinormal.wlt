
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

(* The anti-normal rule set is the image of rulesNO under the automorphism
   a -> a^dagger, a^dagger -> -a, which carries creator-left words to annihilator-left
   ones and sends a^dagger a to -(a a^dagger + 1).  Each expected form below was checked
   against a truncated Fock-space evaluation of both sides before being recorded here. *)


BeginTestSection["BosonicExpOrder - Ordering option"]

(* the default must remain the normal-ordering behaviour that predates the option *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]]] ===
        BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]],
            "Ordering" -> "Normal"],
    True,
    TestID -> "BosonicExpAnti-DefaultIsNormal"
]

(* "Normal" and "Antinormal" are the only accepted settings; anything else must fail
   loudly rather than fall through to one of them *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]],
        "Ordering" -> "NotAnOrdering"],
    $Failed,
    {BosonicExpOrder::ordering},
    TestID -> "BosonicExpAnti-RejectsUnknownOrdering"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - anti-normal number-operator exponentials"]

(* e^(l a a^dag) = A[e^((1 - e^-l) a a^dag)], the mirror of the Stirling/Touchard form.
   Equivalently (1-c)^-(n+1) with c = 1 - e^-l, which is where the sign flip comes from. *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]],
        "Ordering" -> "Antinormal"],
    AntinormalOrdered[E^((1 - E^(-\[Lambda]))*\[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "BosonicExpAnti-NumberOperatorAAdag"
]

(* the same expansion for a^dag a picks up the scalar e^-l, which the Times UpValue
   absorbs into the AntinormalOrdered argument *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*SuperDagger[\[FormalA]]**\[FormalA]],
        "Ordering" -> "Antinormal"],
    AntinormalOrdered[
        E^(-\[Lambda] + (1 - E^(-\[Lambda]))*\[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "BosonicExpAnti-NumberOperatorAdagA"
]

(* "Series" resolves the marker into annihilator-left monomials a^k (a^dag)^k *)
VerificationTest[
    BosonicExpOrder[Exp[\[Lambda]*\[FormalA]**SuperDagger[\[FormalA]]],
        "Ordering" -> "Antinormal"]["Series"],
    Inactive[Sum][
        ((1 - E^(-\[Lambda]))^k*GeneralizedPower[NonCommutativeMultiply, \[FormalA], k]**
            GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], k])/k!,
        {k, 0, Infinity}],
    TestID -> "BosonicExpAnti-NumberOperatorResolving"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - anti-normal linear and quadratic rules"]

(* Weyl splitting: the Gaussian prefactor flips sign and the factors swap order *)
VerificationTest[
    BosonicExpOrder[Exp[\[Alpha]*SuperDagger[\[FormalA]] + \[Beta]*\[FormalA]],
        "Ordering" -> "Antinormal"],
    E^(\[FormalA]*\[Beta])**E^(\[Alpha]*SuperDagger[\[FormalA]])/E^((\[Alpha]*\[Beta])/2),
    TestID -> "BosonicExpAnti-WeylSplitting"
]

(* the reordering rule matches the opposite input shape to its normal counterpart *)
VerificationTest[
    BosonicExpOrder[Exp[\[Beta]*SuperDagger[\[FormalA]]]**Exp[\[Alpha]*\[FormalA]],
        "Ordering" -> "Antinormal"],
    E^(\[FormalA]*\[Alpha])**E^(\[Beta]*SuperDagger[\[FormalA]])/E^(\[Alpha]*\[Beta]),
    TestID -> "BosonicExpAnti-ReorderProduct"
]

(* squeeze from a factored input; the middle factor becomes (1-4ab)^(a a^dag) rather
   than (1-4ab)^-(a^dag a) *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]]**
            Exp[\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]],
        "Ordering" -> "Antinormal"],
    E^((\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2])/
            (1 - 4*\[Alpha]*\[Beta]))**
        (1 - 4*\[Alpha]*\[Beta])^(\[FormalA]**SuperDagger[\[FormalA]])**
        E^((\[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2])/
            (1 - 4*\[Alpha]*\[Beta]))/Sqrt[1 - 4*\[Alpha]*\[Beta]],
    TestID -> "BosonicExpAnti-SqueezeFactored"
]

(* squeeze from a summed input; 1/2 + a^dag a becomes 1/2 - a a^dag *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2] +
            \[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]],
        "Ordering" -> "Antinormal"],
    E^((\[Alpha]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]*
            Tan[2*Sqrt[\[Alpha]*\[Beta]]])/(2*Sqrt[\[Alpha]*\[Beta]]))**
        Sec[2*Sqrt[\[Alpha]*\[Beta]]]^(1/2 - \[FormalA]**SuperDagger[\[FormalA]])**
        E^((\[Beta]*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2]*
            Tan[2*Sqrt[\[Alpha]*\[Beta]]])/(2*Sqrt[\[Alpha]*\[Beta]])),
    TestID -> "BosonicExpAnti-SqueezeSummed"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - anti-normal shift rules"]

(* P(a,a^dag) e^(a a) = e^(a a) A[P(a, a^dag - a)]: the exponential migrates left and the
   shift changes sign, with BosonicAntinormalOrder reducing the shifted polynomial *)
VerificationTest[
    BosonicExpOrder[
        (SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]**\[FormalA])**
            Exp[\[Alpha]*\[FormalA]],
        "Ordering" -> "Antinormal"],
    E^(\[FormalA]*\[Alpha])**(2*\[Alpha] + \[FormalA]*\[Alpha]^2 +
        \[FormalA]**GeneralizedPower[NonCommutativeMultiply, SuperDagger[\[FormalA]], 2] -
        2*\[Alpha]*\[FormalA]**SuperDagger[\[FormalA]] - 2*SuperDagger[\[FormalA]]),
    TestID -> "BosonicExpAnti-ShiftRight"
]

(* e^(b a^dag) P(a,a^dag) = A[P(a - b, a^dag)] e^(b a^dag) *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Beta]*SuperDagger[\[FormalA]]]**
            (SuperDagger[\[FormalA]]**\[FormalA]**\[FormalA]),
        "Ordering" -> "Antinormal"],
    (-2*\[FormalA] + 2*\[Beta] - 2*\[Beta]*\[FormalA]**SuperDagger[\[FormalA]] +
        GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]**SuperDagger[\[FormalA]] +
        \[Beta]^2*SuperDagger[\[FormalA]])**E^(\[Beta]*SuperDagger[\[FormalA]]),
    TestID -> "BosonicExpAnti-ShiftLeft"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - anti-normal Blasiak Thm 4.4"]

(* The mirrored theorem governs one creator flanked by annihilators, e = L + R - 1:
   e^(l a^L a^dag a^R) = (1 + e l a^e)^(-R/e) A[exp((1 - (1 + e l a^e)^(-1/e)) a a^dag)] *)

(* L=2, R=1  ->  e=2, so the prefactor exponent is -R/e = -1/2 *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]**
            SuperDagger[\[FormalA]]**\[FormalA]],
        "Ordering" -> "Antinormal"],
    AntinormalOrdered[
        (1/Sqrt[1 + 2*\[Lambda]*
            GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]])**
        E^((1 - 1/Sqrt[1 + 2*\[Lambda]*
            GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]])**
            \[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "BosonicExpAnti-Thm44-L2R1"
]

(* L=0, R=2  ->  e=1, prefactor exponent -R/e = -2 *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Lambda]*SuperDagger[\[FormalA]]**
            GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]],
        "Ordering" -> "Antinormal"],
    AntinormalOrdered[
        (1 + \[FormalA]*\[Lambda])^(-2)**
        E^(((\[FormalA]*\[Lambda])/(1 + \[FormalA]*\[Lambda]))**
            \[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "BosonicExpAnti-Thm44-L0R2"
]

(* L=2, R=0  ->  e=1 and R=0, so the prefactor is the identity *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Lambda]*GeneralizedPower[NonCommutativeMultiply, \[FormalA], 2]**
            SuperDagger[\[FormalA]]],
        "Ordering" -> "Antinormal"],
    AntinormalOrdered[
        E^(((\[FormalA]*\[Lambda])/(1 + \[FormalA]*\[Lambda]))**
            \[FormalA]**SuperDagger[\[FormalA]])],
    TestID -> "BosonicExpAnti-Thm44-L2R0"
]

EndTestSection[]


BeginTestSection["BosonicExpOrder - anti-normal two mode"]

(* Mode-mixing generators are already mixed in a and a^dagger, so these rules reverse the
   factorisation order rather than producing a strictly annihilator-left form.  Note the
   two modes exchange the sign on log Cosh relative to the normal rule. *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Alpha]*SuperDagger[\[FormalA]]**\[FormalB] +
            \[Beta]*\[FormalA]**SuperDagger[\[FormalB]]],
        "Ordering" -> "Antinormal"],
    E^((\[Beta]*\[FormalA]**SuperDagger[\[FormalB]]*Tanh[Sqrt[\[Alpha]*\[Beta]]])/
            Sqrt[\[Alpha]*\[Beta]])**
        Cosh[Sqrt[\[Alpha]*\[Beta]]]^(\[FormalA]**SuperDagger[\[FormalA]] -
            \[FormalB]**SuperDagger[\[FormalB]])**
        E^((\[Alpha]*SuperDagger[\[FormalA]]**\[FormalB]*Tanh[Sqrt[\[Alpha]*\[Beta]]])/
            Sqrt[\[Alpha]*\[Beta]]),
    TestID -> "BosonicExpAnti-BeamSplitter"
]

(* the two-mode squeezer does come out strictly annihilator-left *)
VerificationTest[
    BosonicExpOrder[
        Exp[\[Alpha]*\[FormalA]**\[FormalB] +
            \[Beta]*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]],
        "Ordering" -> "Antinormal"],
    E^((\[Alpha]*\[FormalA]**\[FormalB]*Tan[Sqrt[\[Alpha]*\[Beta]]])/
            Sqrt[\[Alpha]*\[Beta]])**Sec[Sqrt[\[Alpha]*\[Beta]]]*
        Sec[Sqrt[\[Alpha]*\[Beta]]]^(-\[FormalA]**SuperDagger[\[FormalA]] -
            \[FormalB]**SuperDagger[\[FormalB]])**
        E^((\[Beta]*SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]*
            Tan[Sqrt[\[Alpha]*\[Beta]]])/Sqrt[\[Alpha]*\[Beta]]),
    TestID -> "BosonicExpAnti-TwoModeSqueezer"
]

EndTestSection[]
