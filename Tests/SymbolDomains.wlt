BeginTestSection["SymbolDomains"]

(* A basis records what is known about each symbol: "SymbolDomains" maps a symbol, plain or indexed,
   to the domains it lies in. Roles fill it, the user may add to it, and a symbol the
   user marks Complexes is never taken as real. *)

(* ========== roles ========== *)

(* A gate angle is real. An affine a theta + b with real a != 0 and b is real only for real theta, so
   theta joins; a + b says nothing about a or b; a number adds nothing. *)
VerificationTest[
    {
        QuantumOperator["RZ"[\[Theta]]]["Basis"]["SymbolDomains"],
        QuantumOperator["RZ"[2 \[Theta] + 1]]["Basis"]["SymbolDomains"],
        QuantumOperator["RZ"[\[Theta] / 2]]["Basis"]["SymbolDomains"],
        QuantumOperator["RZ"[a + b]]["Basis"]["SymbolDomains"],
        QuantumOperator["RZ"[0.3]]["Basis"]["SymbolDomains"],
        QuantumOperator["U3"[t[1], Subscript[f, 2], l]]["Basis"]["SymbolDomains"]
    },
    {
        <|\[Theta] -> {Reals}|>, <|\[Theta] -> {Reals}|>, <|\[Theta] -> {Reals}|>, <||>, <||>,
        <|t[1] -> {Reals}, Subscript[f, 2] -> {Reals}, l -> {Reals}|>
    },
    TestID -> "SymbolDomains-GateAngles"
]

(* A channel probability lies in [0, 1]; for Depolarizing in [0, 4/3], where the Kraus weight 1 - 3 p / 4
   stays nonnegative. *)
VerificationTest[
    {
        QuantumChannel["BitFlip"[p]]["Basis"]["SymbolDomains"],
        QuantumChannel["GeneralizedAmplitudeDamping"[g, p]]["ParameterAssumptions"],
        QuantumChannel["Depolarizing"[q]]["ParameterAssumptions"],
        QuantumChannel["ResetError"[p0, p1]]["ParameterAssumptions"],
        QuantumChannel["BitFlip"[0.1]]["Basis"]["SymbolDomains"]
    },
    {
        <|p -> {Interval[{0, 1}]}|>,
        0 <= g <= 1 && 0 <= p <= 1,
        0 <= q <= 4 / 3,
        0 <= p0 <= 1 && 0 <= p1 <= 1,
        <||>
    },
    TestID -> "SymbolDomains-ChannelProbabilities"
]

(* The Werner weight lies in [0, 1]; a Bloch vector is real. *)
VerificationTest[
    {
        QuantumState["Werner"[p]]["ParameterAssumptions"],
        QuantumState["BlochVector"[{x, y, z}]]["ParameterAssumptions"]
    },
    {0 <= p <= 1, Element[x | y | z, Reals]},
    TestID -> "SymbolDomains-NamedStates"
]

(* ========== merging ========== *)

(* Records join through circuits, tensor products and application; a symbol in two roles lies in both
   domains, so the angle and the probability p make p in [0, 1]. *)
VerificationTest[
    {
        QuantumCircuitOperator[{"RZ"[\[Theta]] -> 1, "BitFlip"[p] -> 1, "RX"[\[Theta]] -> 2}]["Basis"]["SymbolDomains"],
        QuantumCircuitOperator[{"RZ"[p] -> 1, "BitFlip"[p] -> 1}]["Basis"]["SymbolDomains"],
        QuantumTensorProduct[QuantumOperator["RZ"[a]], QuantumChannel["BitFlip"[p]]["QuantumOperator"]]["Basis"]["SymbolDomains"],
        QuantumChannel["BitFlip"[p]][QuantumState["0"]]["Basis"]["SymbolDomains"]
    },
    {
        <|\[Theta] -> {Reals}, p -> {Interval[{0, 1}]}|>,
        <|p -> {Interval[{0, 1}]}|>,
        <|a -> {Reals}, p -> {Interval[{0, 1}]}|>,
        <|p -> {Interval[{0, 1}]}|>
    },
    TestID -> "SymbolDomains-Merge"
]

(* ========== the user's record ========== *)

(* A user record joins the role record: rules or an association, a domain or a list of them; a domain
   containing another one listed is dropped, and domains neither of which contains the other both stay. *)
VerificationTest[
    {
        QuantumOperator["RZ"[t], "SymbolDomains" -> <|t -> NonNegativeReals|>]["Basis"]["SymbolDomains"],
        QuantumState[{a, b}, "SymbolDomains" -> {a -> {PositiveReals, Interval[{0, 1}]}, b -> {Reals, NonNegativeIntegers}}]["Basis"]["SymbolDomains"]
    },
    {
        <|t -> {NonNegativeReals}|>,
        <|a -> {PositiveReals, Interval[{0, 1}]}, b -> {NonNegativeIntegers}|>
    },
    TestID -> "SymbolDomains-UserRecord"
]

(* Complexes overrides every role and the declaration: RZ[t] on |0> then keeps the norm of a complex
   angle, E^(Im[t] / 2), which at t = I is Sqrt[E]. *)
VerificationTest[
    With[{rz = QuantumOperator["RZ"[t], "SymbolDomains" -> <|t -> Complexes|>]},
        {
            rz["Basis"]["SymbolDomains"], rz["ParameterAssumptions"], rz["Basis"]["RealSymbols"],
            rz[QuantumState[{1, 0}]]["Norm"] /. t -> I,
            QuantumState[{Cos[x], Sin[x]}, "Parameters" -> {x}, "SymbolDomains" -> <|x -> Complexes|>]["ParameterAssumptions"]
        }
    ],
    {<|t -> {Complexes}|>, True, {}, Sqrt[E], True},
    TestID -> "SymbolDomains-ComplexOverride"
]

(* Facts name symbols, plain or indexed, never patterns: a pattern key is refused. *)
VerificationTest[
    FailureQ /@ {
        QuantumState[{a, b}, "SymbolDomains" -> <|Subscript[t, _] -> Reals|>],
        QuantumState[{a, b}, "SymbolDomains" -> <|_ -> Reals|>],
        QuantumState[{a, b}, "SymbolDomains" -> {Subscript[t, 1 | 2] -> Reals}],
        QuantumState[{a, b}, "SymbolDomains" -> {Subscript[t, Except[0]] -> Reals}],
        QuantumState[{a, b}, "SymbolDomains" -> {t[_] -> Reals}]
    },
    {True, True, True, True, True},
    {QuantumBasis::domainspec, QuantumBasis::domainspec, QuantumBasis::domainspec, General::stop},
    TestID -> "SymbolDomains-PatternRefused"
]

(* Substituting a parameter spends its record. *)
VerificationTest[
    QuantumOperator["RZ"[\[Theta]], "Parameters" -> {\[Theta]}][<|\[Theta] -> 0.3|>]["Basis"]["SymbolDomains"],
    <||>,
    TestID -> "SymbolDomains-SubstitutionSpends"
]

(* An unknown domain is refused with a message. *)
VerificationTest[
    FailureQ[QuantumState[{a, b}, "SymbolDomains" -> <|a -> Foo|>]],
    True,
    {QuantumBasis::domain},
    TestID -> "SymbolDomains-UnknownDomain"
]

(* A domain named twice is kept once, and intervals meet in one interval. *)
VerificationTest[
    {
        QuantumState[{a, b}, "SymbolDomains" -> {a -> Interval[{0, 1}], a -> Interval[{0, 1}]}]["Basis"]["SymbolDomains"],
        QuantumState[{a, b}, "SymbolDomains" -> <|a -> {Reals, Interval[{0, 2}], NonNegativeReals, Interval[{1, 3}]}|>]["Basis"]["SymbolDomains"]
    },
    {<|a -> {Interval[{0, 1}]}|>, <|a -> {Interval[{1, 2}]}|>},
    TestID -> "SymbolDomains-EqualAndMeet"
]

(* Domains with no common point give one message and no fact, so nothing is simplified under false facts;
   the symbol keeps an empty list so that no later merge gives it a fact again. *)
VerificationTest[
    With[{s = QuantumState[{a, b}, "SymbolDomains" -> <|a -> {Interval[{0, 1}], Interval[{2, 3}]}|>]},
        {s["Basis"]["SymbolDomains"], s["ParameterAssumptions"]}
    ],
    {<|a -> {}|>, True},
    {QuantumBasis::domainempty},
    TestID -> "SymbolDomains-Contradiction"
]

VerificationTest[
    QuantumState[{a, b}, "SymbolDomains" -> <|a -> {Integers, Interval[{1/4, 3/4}]}|>]["Basis"]["SymbolDomains"],
    <|a -> {}|>,
    {QuantumBasis::domainempty},
    TestID -> "SymbolDomains-NoIntegerInInterval"
]

(* A malformed specification is refused with one message, wherever it is given. *)
VerificationTest[
    FailureQ[QuantumOperator["RZ"[a], "SymbolDomains" -> <|a -> Foo|>]],
    True,
    {QuantumBasis::domain},
    TestID -> "SymbolDomains-UnknownDomainOnGate"
]

VerificationTest[
    FailureQ /@ {
        QuantumState[{a, b}, "SymbolDomains" -> {"x" -> Reals}],
        QuantumState[{a, b}, "SymbolDomains" -> {a + b -> PositiveReals}],
        QuantumState[{a, b}, "SymbolDomains" -> {a :> Reals}]
    },
    {True, True, True},
    {QuantumBasis::domainspec, QuantumBasis::domainspec, QuantumBasis::domainspec, General::stop},
    TestID -> "SymbolDomains-MalformedSpec"
]

(* An empty list takes away every fact about a symbol, the role's included. *)
VerificationTest[
    With[{o = QuantumOperator["RZ"[a], "SymbolDomains" -> <|a -> {}|>]}, {o["Basis"]["SymbolDomains"], o["ParameterAssumptions"]}],
    {<|a -> {}|>, True},
    TestID -> "SymbolDomains-EmptyListNoFacts"
]

(* The record survives the operations that rebuild the basis. *)
VerificationTest[
    With[{s = QuantumState["BlochVector"[{x, y, 0}]]},
        #["Basis"]["SymbolDomains"] & /@ {
            s["Computational"],
            QuantumPartialTrace[QuantumTensorProduct[s, s], {1}],
            QuantumOperator["RZ"[a]]["Transpose"],
            QuantumEvolve[QuantumOperator["Z"], QuantumState["BlochVector"[{x, 0, 0}]], t]
        }
    ],
    {<|x -> {Reals}, y -> {Reals}|>, <|x -> {Reals}, y -> {Reals}|>, <|a -> {Reals}|>, <|x -> {Interval[{-1, 1}]}|>},
    TestID -> "SymbolDomains-Carried"
]

(* Conflicting facts from two objects: the sum has none for the symbol, in either order, and a later
   composition does not bring one back. *)
VerificationTest[
    With[{
        s1 = QuantumState[{a, b}, "SymbolDomains" -> {a -> Interval[{0, 1}]}],
        s2 = QuantumState[{a, b}, "SymbolDomains" -> {a -> Interval[{2, 3}]}],
        o1 = QuantumOperator[{{a, 0}, {0, 1}}, "SymbolDomains" -> {a -> Interval[{0, 1}]}],
        o2 = QuantumOperator[{{a, 0}, {0, 1}}, "SymbolDomains" -> {a -> Interval[{2, 3}]}]
    },
        #["Basis"]["SymbolDomains"] & /@ {s1 + s2, s2 + s1, o1[o2[o1]]}
    ],
    {<|a -> {}|>, <|a -> {}|>, <|a -> {}|>},
    {QuantumBasis::domainempty, QuantumBasis::domainempty, QuantumBasis::domainempty, General::stop},
    TestID -> "SymbolDomains-ConflictStays"
]

(* A symbol named twice in one specification gets both facts; Complexes still wins. *)
VerificationTest[
    QuantumState[{a, b}, "SymbolDomains" -> {a -> Complexes, a -> Reals, b -> Reals, b -> Interval[{0, 2}]}]["Basis"]["SymbolDomains"],
    <|a -> {Complexes}, b -> {Interval[{0, 2}]}|>,
    TestID -> "SymbolDomains-RepeatedKey"
]

(* Formal symbols are symbols too. *)
VerificationTest[
    QuantumOperator["RZ"[\[FormalT]]]["Basis"]["SymbolDomains"],
    <|\[FormalT] -> {Reals}|>,
    TestID -> "SymbolDomains-FormalSymbol"
]

(* A bad specification added to a channel fails as a whole. *)
VerificationTest[
    FailureQ @ QuantumChannel[QuantumChannel["BitFlip"[p]], "SymbolDomains" -> <|p -> Foo|>],
    True,
    {QuantumBasis::domain},
    TestID -> "SymbolDomains-ChannelBadSpec"
]

(* The direct sum of states keeps the facts of each. *)
VerificationTest[
    QuantumState[QuantumState["BlochVector"[{x, 0, 0}]], QuantumState["0"]]["Basis"]["SymbolDomains"],
    <|x -> {Interval[{-1, 1}]}|>,
    TestID -> "SymbolDomains-DirectSum"
]

(* A sum or a difference of operators whose facts conflict reports the conflict once. *)
VerificationTest[
    With[{
        o1 = QuantumOperator[{{a, 0}, {0, 1}}, "SymbolDomains" -> {a -> Interval[{0, 1}]}],
        o2 = QuantumOperator[{{a, 0}, {0, 1}}, "SymbolDomains" -> {a -> Interval[{2, 3}]}]
    },
        (o1 + o2)["Basis"]["SymbolDomains"]
    ],
    <|a -> {}|>,
    {QuantumBasis::domainempty},
    TestID -> "SymbolDomains-SumMessageOnce"
]

VerificationTest[
    With[{
        o1 = QuantumOperator[{{a, 0}, {0, 1}}, "SymbolDomains" -> {a -> Interval[{0, 1}]}],
        o2 = QuantumOperator[{{a, 0}, {0, 1}}, "SymbolDomains" -> {a -> Interval[{2, 3}]}]
    },
        (o1 - o2)["Basis"]["SymbolDomains"]
    ],
    <|a -> {}|>,
    {QuantumBasis::domainempty},
    TestID -> "SymbolDomains-DifferenceMessageOnce"
]

(* A change of basis does not give back the facts taken away from a symbol. *)
VerificationTest[
    QuantumState[
        QuantumState[{a, 1}, "SymbolDomains" -> <|a -> {}|>],
        QuantumBasis[2, "SymbolDomains" -> {a -> Reals, b -> Reals}]
    ]["Basis"]["SymbolDomains"],
    <|a -> {}, b -> {Reals}|>,
    TestID -> "SymbolDomains-ChangeOfBasisKeepsCleared"
]

(* An empty list takes away every fact in each form of a specification. *)
VerificationTest[
    QuantumState[{a, b}, "SymbolDomains" -> #]["Basis"]["SymbolDomains"] & /@ {{a -> {}}, a -> {}, {a -> Reals, a -> {}, b -> Reals}},
    {<|a -> {}|>, <|a -> {}|>, <|a -> {}, b -> {Reals}|>},
    TestID -> "SymbolDomains-EmptyListAllForms"
]

(* Numeric slots give no facts. *)
VerificationTest[
    #["Basis"]["SymbolDomains"] & /@ {QuantumState["BlochVector"[{0.1, 0, 0.2}]], QuantumState["Werner"[0.3]], QuantumOperator["RX"[0.4]], QuantumChannel["BitFlip"[0.1]]},
    {<||>, <||>, <||>, <||>},
    TestID -> "SymbolDomains-NumericNoFacts"
]

(* A value nested deeper than a list of domains is not a specification. *)
VerificationTest[
    FailureQ /@ {
        QuantumState[{a, b}, "SymbolDomains" -> {a -> {Reals, {}}}],
        QuantumState[{a, b}, "SymbolDomains" -> {a -> {{Reals, {Integers}}}}],
        QuantumState[{a, b}, "SymbolDomains" -> {a -> {{}}}]
    },
    {True, True, True},
    {QuantumBasis::domainspec, QuantumBasis::domainspec, QuantumBasis::domainspec, General::stop},
    TestID -> "SymbolDomains-NestedValueRejected"
]

(* An interval of one integer lies in the integers. *)
VerificationTest[
    QuantumState[{a, b}, "SymbolDomains" -> {a -> {Interval[{0, 0}], Integers}}]["Basis"]["SymbolDomains"],
    <|a -> {Interval[{0, 0}]}|>,
    TestID -> "SymbolDomains-IntegerPointInterval"
]

(* A change of basis joins the facts of both bases, so a conflict between them is found. *)
VerificationTest[
    QuantumState[
        QuantumState[{a, 1}, "SymbolDomains" -> {a -> Interval[{0, 1}], b -> Interval[{0, 2}]}],
        QuantumBasis[2, "SymbolDomains" -> {a -> Interval[{2, 3}], b -> Reals}]
    ]["Basis"]["SymbolDomains"],
    <|a -> {}, b -> {Interval[{0, 2}]}|>,
    {QuantumBasis::domainempty},
    TestID -> "SymbolDomains-ChangeOfBasisConflict"
]

(* t[k] and Subscript[t, k] are indexed symbols alike. *)
VerificationTest[
    {
        QuantumOperator["RZ"[t[k]]]["Basis"]["SymbolDomains"],
        QuantumOperator["RZ"[Subscript[t, k]]]["Basis"]["SymbolDomains"],
        QuantumState[{a, b}, "SymbolDomains" -> {t[k] -> Reals}]["Basis"]["SymbolDomains"]
    },
    {<|t[k] -> {Reals}|>, <|Subscript[t, k] -> {Reals}|>, <|t[k] -> {Reals}|>},
    TestID -> "SymbolDomains-IndexForms"
]

(* Spider phases and Euler angles are real; a probability slot holding a p + b puts p in the preimage. *)
VerificationTest[
    {
        QuantumOperator["ZSpider"[t], {1} -> {1}]["Basis"]["SymbolDomains"],
        QuantumOperator["WignerD"[1/2, {a, b, c}]]["Basis"]["SymbolDomains"],
        QuantumChannel["BitFlip"[1 - p]]["Basis"]["SymbolDomains"],
        QuantumChannel["BitFlip"[2 q]]["Basis"]["SymbolDomains"],
        QuantumChannel["ResetError"[r, r]]["Basis"]["SymbolDomains"]
    },
    {
        <|t -> {Reals}|>,
        <|a -> {Reals}, b -> {Reals}, c -> {Reals}|>,
        <|p -> {Interval[{0, 1}]}|>,
        <|q -> {Interval[{0, 1/2}]}|>,
        <|r -> {Interval[{0, 1/2}]}|>
    },
    TestID -> "SymbolDomains-MoreRoles"
]

(* The recorded range of a role is the exact region where the object is physical. For a channel, trace
   preservation on |0> and |1>, reduced over the reals, gives the region in closed form; for a Werner state
   and a qubit Bloch vector with one symbolic component, positivity of the eigenvalues does. Each derived
   region is equivalent to the recorded interval. The region is derived from the object with its facts taken
   away, since with them its trace and eigenvalues are computed as if the facts already held. *)
VerificationTest[
    With[{
        recorded = Function[{obj, s}, With[{i = First[obj["Basis"]["SymbolDomains"][s]]}, i[[1, 1]] <= s <= i[[1, 2]]]],
        tracePreserving = Function[{c, s}, With[{bare = QuantumChannel[c, "SymbolDomains" -> {s -> {}}]},
            Reduce[And @@ (ComplexExpand[Tr[bare[QuantumState[#]]["DensityMatrix"]]] == 1 & /@ {"0", "1"}), s, Reals]]],
        positive = Function[{st, s}, With[{bare = QuantumState[st, "SymbolDomains" -> {s -> {}}]},
            Reduce[Min[Eigenvalues[Normal[bare["DensityMatrix"]]]] >= 0, s, Reals]]]
    },
        Join[
            Reduce[Equivalent[tracePreserving[#1, #2], recorded[#1, #2]], #2, Reals] & @@@ {
                {QuantumChannel["BitFlip"[p]], p},
                {QuantumChannel["AmplitudeDamping"[g]], g},
                {QuantumChannel["PhaseDamping"[l]], l},
                {QuantumChannel["Depolarizing"[q]], q},
                {QuantumChannel["ResetError"[r, r]], r},
                {QuantumChannel["BitFlip"[1 - p]], p}
            },
            Reduce[Equivalent[positive[#1, #2], recorded[#1, #2]], #2, Reals] & @@@ {
                {QuantumState["Werner"[p]], p},
                {QuantumState["Werner"[p, 3]], p},
                {QuantumState["BlochVector"[{x, 1/2, 0}]], x}
            }
        ]
    ],
    ConstantArray[True, 9],
    TestID -> "SymbolDomains-PhysicalRegion"
]

(* An inexact coefficient gives no interval, since Interval widens inexact bounds; a real slot still gives Reals. *)
VerificationTest[
    #["Basis"]["SymbolDomains"] & /@ {QuantumChannel["BitFlip"[0.5 p]], QuantumOperator["RZ"[0.5 t + 0.1]], QuantumState["BlochVector"[{x, 0.5, 0}]]},
    {<||>, <|t -> {Reals}|>, <|x -> {Reals}|>},
    TestID -> "SymbolDomains-InexactCoefficients"
]

(* Facts from several roles on several wires meet in one record: a two-qubit Werner state, a controlled
   rotation, a channel on one wire and an evolution share symbols. *)
VerificationTest[
    With[{
        rho = QuantumState["Werner"[p]],
        cp = QuantumOperator["C"["RZ"[2 t + 1]], {1, 2}],
        noise = QuantumChannel["BitFlip"[1 - p], {2}]
    },
        With[{out = noise[cp[rho]]},
            {out["Basis"]["SymbolDomains"], out["ParameterAssumptions"],
             QuantumEvolve[QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["I"]], out, s]["Basis"]["SymbolDomains"]}
        ]
    ],
    {
        <|p -> {Interval[{0, 1}]}, t -> {Reals}|>,
        Element[t, Reals] && 0 <= p <= 1,
        <|p -> {Interval[{0, 1}]}, t -> {Reals}|>
    },
    TestID -> "SymbolDomains-ManyRoles"
]

(* Exact irrational coefficients give intervals; an impossible Bloch vector gives no warning; a phase-shift
   index is an integer. *)
VerificationTest[
    {
        QuantumChannel["BitFlip"[Pi p]]["Basis"]["SymbolDomains"],
        QuantumState["BlochVector"[{x, 1 / Sqrt[2], 0}]]["Basis"]["SymbolDomains"],
        QuantumState["BlochVector"[{x, 1, 1}]]["Basis"]["SymbolDomains"],
        QuantumOperator["PhaseShift"[k]]["Basis"]["SymbolDomains"]
    },
    {
        <|p -> {Interval[{0, 1 / Pi}]}|>,
        <|x -> {Interval[{- 1 / Sqrt[2], 1 / Sqrt[2]}]}|>,
        <|x -> {Reals}|>,
        <|k -> {Integers}|>
    },
    TestID -> "SymbolDomains-ExactConstantsAndIntegers"
]


(* ========== computing on the facts ========== *)

(* Conjugation takes the real symbols real: the dagger of RZ[theta] is RZ[-theta] entry by entry, the
   product with RZ[theta] is the identity, and the same for U3. *)
VerificationTest[
    With[{rz = QuantumOperator["RZ"[\[Theta]]], u3 = QuantumOperator["U3"[\[Theta], \[Phi], \[Lambda]]]}, {
        Normal[rz["Dagger"]["Matrix"]] === Normal[QuantumOperator["RZ"[- \[Theta]]]["Matrix"]],
        Normal[(rz["Dagger"] @ rz)["Matrix"]],
        FreeQ[Normal[u3["Dagger"]["Matrix"]], Conjugate],
        Simplify[Normal[(u3["Dagger"] @ u3)["Matrix"]]]
    }],
    {True, IdentityMatrix[2], True, IdentityMatrix[2]},
    TestID -> "SymbolDomains-DaggerOfRealAngles"
]

(* A symbol the user marks complex keeps its conjugate; so does the square root of a real symbol whose
   sign is unknown, and of a symbol with no fact at all. *)
VerificationTest[
    {
        FreeQ[Normal[QuantumOperator["RZ"[\[Theta]], "SymbolDomains" -> {\[Theta] -> Complexes}]["Dagger"]["Matrix"]], Conjugate],
        FreeQ[Normal[QuantumState[{Sqrt[x], 1}, "Parameters" -> {x}]["DensityMatrix"]], Conjugate],
        FreeQ[Normal[QuantumState[{Sqrt[q], 1}]["DensityMatrix"]], Conjugate]
    },
    {False, False, False},
    TestID -> "SymbolDomains-ConjugateKeptWithoutFacts"
]

(* A square root of a probability is real on [0, 1], so a channel's output is free of Conjugate: bit flip
   gives diag(1 - p, p), and amplitude damping the familiar populations and coherences. *)
VerificationTest[
    {
        Normal[QuantumChannel["BitFlip"[p]][QuantumState["0"]]["DensityMatrix"]],
        Expand @ Normal[QuantumChannel["AmplitudeDamping"[g]][QuantumState[{\[Alpha], \[Beta]}]]["DensityMatrix"]]
    },
    {
        {{1 - p, 0}, {0, p}},
        Expand @ {
            {\[Alpha] Conjugate[\[Alpha]] + g \[Beta] Conjugate[\[Beta]], Sqrt[1 - g] \[Alpha] Conjugate[\[Beta]]},
            {Sqrt[1 - g] \[Beta] Conjugate[\[Alpha]], (1 - g) \[Beta] Conjugate[\[Beta]]}
        }
    },
    TestID -> "SymbolDomains-ChannelOutputs"
]

(* Symbolic unitarity and Hermiticity hold for every value the facts allow, decided exactly: a three-qubit
   circuit of rotations and CNOTs is unitary for real angles; RZ is not Hermitian; a rotation with a complex
   angle is not unitary; a symbolic Hermitian matrix with real parameters is Hermitian; amplitude damping
   is trace preserving. Numeric operators answer as before. *)
VerificationTest[
    {
        QuantumOperator[QuantumCircuitOperator[{
            "RY"[Subscript[\[Theta], 1]] -> 1, "RY"[Subscript[\[Theta], 2]] -> 2, "RY"[Subscript[\[Theta], 3]] -> 3,
            "CNOT" -> {1, 2}, "CNOT" -> {2, 3},
            "RZ"[Subscript[\[Phi], 1]] -> 1, "RZ"[Subscript[\[Phi], 2]] -> 2, "RZ"[Subscript[\[Phi], 3]] -> 3
        }]]["UnitaryQ"],
        QuantumOperator["RZ"[\[Theta]]]["HermitianQ"],
        QuantumOperator["RZ"[\[Theta]], "SymbolDomains" -> {\[Theta] -> Complexes}]["UnitaryQ"],
        QuantumOperator[{{a, b - I c}, {b + I c, d}}, "Parameters" -> {a, b, c, d}]["HermitianQ"],
        QuantumChannel["AmplitudeDamping"[g]]["TracePreservingQ"],
        QuantumOperator["RZ"[0.3]]["UnitaryQ"],
        QuantumOperator["Y"]["HermitianQ"]
    },
    {True, False, False, True, True, True, True},
    TestID -> "SymbolDomains-SymbolicUnitaryHermitian"
]

(* The three answers of the zero test, each from its own stage: proven zero (TrigToExp then Together), nonzero
   at a point the facts allow, and undecided when neither applies (Log[E^x] - x vanishes for real x, but no
   rational form shows it). *)
VerificationTest[
    With[{qb = QuantumBasis["SymbolDomains" -> {x -> Reals}]}, {
        zeroOnDomains[{Sin[x] ^ 2 + Cos[x] ^ 2 - 1, Cos[x / 2] ^ 2 - (1 + Cos[x]) / 2}, qb],
        zeroOnDomains[{Sin[x] ^ 2 + Cos[x] ^ 2 - 1, x - Conjugate[y]}, qb],
        zeroOnDomains[{Log[Exp[x]] - x}, qb]
    }],
    {True, False, Indeterminate},
    TestID -> "SymbolDomains-ZeroTestStages"
]

(* A symbolic density matrix is physical when it is positive semidefinite for every value the facts allow.
   Channel outputs, Werner states for d = 2 and 3, a Bloch vector with one symbolic component, and an
   amplitude-damped state with complex amplitudes all are; diag(p, 1 - p) is physical for p in [0, 1], and for
   p in [0, 2] only where p <= 1, so there the answer is Indeterminate; a pure state written as a density
   matrix is classified pure. *)
VerificationTest[
    {
        QuantumChannel["BitFlip"[p]][QuantumState["0"]]["PhysicalQ"],
        QuantumState["Werner"[p]]["PhysicalQ"],
        QuantumState["Werner"[p, 3]]["PhysicalQ"],
        QuantumState["BlochVector"[{x, 1 / 2, 0}]]["PhysicalQ"],
        QuantumChannel["AmplitudeDamping"[g]][QuantumState[{\[Alpha], \[Beta]}]]["MixedStateQ"],
        QuantumState[{{p, 0}, {0, 1 - p}}, "SymbolDomains" -> {p -> Interval[{0, 1}]}]["Type"],
        QuantumState[{{p, 0}, {0, 1 - p}}, "SymbolDomains" -> {p -> Interval[{0, 2}]}]["PhysicalQ"],
        QuantumState[QuantumState[{Cos[t], Sin[t]}, "SymbolDomains" -> {t -> Reals}]["DensityMatrix"], "SymbolDomains" -> {t -> Reals}]["Type"]
    },
    {True, True, True, True, True, "Mixed", Indeterminate, "Pure"},
    {QuantumState::undecidedregion},
    TestID -> "SymbolDomains-SymbolicPhysicalStates"
]

(* The private symbols used during a computation never outlive it. *)
VerificationTest[
    QuantumOperator["U3"[\[Theta], \[Phi], \[Lambda]]]["Dagger"];
    QuantumChannel["BitFlip"[p]][QuantumState["0"]]["DensityMatrix"];
    Select[Names["Wolfram`QuantumFramework`*`r$*"], StringMatchQ[#, __ ~~ "r$" ~~ DigitCharacter ..] &],
    {},
    TestID -> "SymbolDomains-NoPrivateSymbolsLeft"
]

(* ========== yes/no answers that cannot be decided, and declarations ========== *)

$rabi = MatrixExp[- I t / 2 (\[CapitalDelta] PauliMatrix[3] + \[CapitalOmega] PauliMatrix[1])];
$generator = \[CapitalDelta] PauliMatrix[3] + \[CapitalOmega] PauliMatrix[1];
$split = QuantumState[DiagonalMatrix[{p, 2 - p}], "Parameters" -> {p}, "SymbolDomains" -> {p -> Interval[{0, 3}]}];

(* A density matrix with declared parameters is checked against its facts like any other: diag(p, 2 - p) is
   physical for p in [0, 2] and for no p in [3, 4]. *)
VerificationTest[
    {
        QuantumState[DiagonalMatrix[{p, 2 - p}], "Parameters" -> {p}, "SymbolDomains" -> {p -> Interval[{0, 2}]}]["PhysicalQ"],
        QuantumState[DiagonalMatrix[{p, 2 - p}], "Parameters" -> {p}, "SymbolDomains" -> {p -> Interval[{3, 4}]}]["PhysicalQ"]
    },
    {True, False},
    TestID -> "SymbolDomains-PhysicalQ-DeclaredParameters"
]

(* For p in [0, 3] it is physical only where p <= 2: the answer is Indeterminate, the message names where it
   holds, and Indeterminate combines with And and Or as an unknown does. "Physical" leaves such a state as it
   is and says nothing, since it was not asked whether the state is physical. *)
VerificationTest[
    {$split["PhysicalQ"], $split["PhysicalQ"] && False, $split["PhysicalQ"] || True, QuantumStateQ[$split["Physical"]]},
    {Indeterminate, False, True, True},
    {QuantumState::undecidedregion, QuantumState::undecidedregion, QuantumState::undecidedregion, General::stop},
    TestID -> "SymbolDomains-PhysicalQ-Split"
]

(* A pure state written by hand with Conjugate[Sin[t]], t real, is pure. *)
VerificationTest[
    QuantumState[{{Cos[t] ^ 2, Cos[t] Conjugate[Sin[t]]}, {Sin[t] Cos[t], Sin[t] ^ 2}}, "SymbolDomains" -> {t -> Reals}]["Type"],
    "Pure",
    TestID -> "SymbolDomains-Type-PureWithConjugate"
]

(* With undeclared symbols, which may be complex, a Rabi propagator is not unitary and its generator not
   Hermitian; the message names the symbols whose declaration as real would change the answer. *)
VerificationTest[
    {QuantumOperator[$rabi]["UnitaryQ"], QuantumOperator[$generator]["HermitianQ"]},
    {False, False},
    {QuantumOperator::undeclared, QuantumOperator::undeclared},
    TestID -> "SymbolDomains-Undeclared-Hint"
]

(* "SymbolDomains" -> Reals makes every symbol of the object real; it meets a role's range rather than
   replacing it. *)
VerificationTest[
    {
        QuantumOperator[$rabi, "SymbolDomains" -> Reals]["UnitaryQ"],
        QuantumOperator[$generator, "SymbolDomains" -> Reals]["HermitianQ"],
        QuantumOperator[$generator, "SymbolDomains" -> Reals]["Basis"]["SymbolDomains"],
        QuantumChannel["BitFlip"[p], "SymbolDomains" -> Reals]["Basis"]["SymbolDomains"]
    },
    {True, True, <|\[CapitalDelta] -> {Reals}, \[CapitalOmega] -> {Reals}|>, <|p -> {Interval[{0, 1}]}|>},
    TestID -> "SymbolDomains-AllReal"
]

(* A user's explicit Reals or Complexes for a symbol replaces its role's range; a user interval meets it. *)
VerificationTest[
    {
        QuantumChannel["BitFlip"[p], "SymbolDomains" -> {p -> Reals}]["Basis"]["SymbolDomains"],
        QuantumChannel["BitFlip"[p], "SymbolDomains" -> {p -> Complexes}]["Basis"]["SymbolDomains"],
        QuantumChannel["BitFlip"[p], "SymbolDomains" -> {p -> Interval[{0, 2}]}]["Basis"]["SymbolDomains"]
    },
    {<|p -> {Reals}|>, <|p -> {Complexes}|>, <|p -> {Interval[{0, 1}]}|>},
    TestID -> "SymbolDomains-ExplicitDomainReplacesRole"
]

(* A complex mark beside a range in one declaration contradicts it: one message, no fact. *)
VerificationTest[
    QuantumState[{a, b}, "SymbolDomains" -> {a -> {Complexes, Interval[{0, 1}]}}]["Basis"]["SymbolDomains"],
    <|a -> {}|>,
    {QuantumBasis::complexbounded},
    TestID -> "SymbolDomains-ComplexWithRange"
]

(* An entry that is not a number, a product of noncommuting symbols, leaves the yes/no answer undecided. *)
VerificationTest[
    QuantumOperator[{{a ** b, 0}, {0, 1}}]["UnitaryQ"],
    Indeterminate,
    {QuantumOperator::notnumber},
    TestID -> "SymbolDomains-NotANumber"
]

(* A basis with no facts lists no "SymbolDomains" among its options. *)
VerificationTest[
    {MemberQ[Keys[QuantumBasis[]["Options"]], "SymbolDomains"], MemberQ[Keys[QuantumOperator["RZ"[t]]["Basis"]["Options"]], "SymbolDomains"]},
    {False, True},
    TestID -> "SymbolDomains-OptionsHideEmpty"
]

EndTestSection[]
