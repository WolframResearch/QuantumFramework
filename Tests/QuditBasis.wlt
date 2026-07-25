BeginTestSection["QuditBasis - constructors"]

VerificationTest[QuditBasis["Bell"]["Dimension"], 4, TestID -> "Bell-bare"]

VerificationTest[QuditBasis["Fourier"]["Dimension"], 2, TestID -> "Fourier-default"]

VerificationTest[QuditBasis["Fourier"[3]]["Dimension"], 3, TestID -> "Fourier-call-3"]

VerificationTest[QuditBasis["PauliX"[3]]["Dimension"], 3, TestID -> "PauliX-3"]

VerificationTest[QuditBasis["GellMann"[3]]["Dimension"], 9, TestID -> "GellMann-3"]

VerificationTest[QuditBasis["HesseSIC"]["Dimension"], 9, TestID -> "HesseSIC-bare"]

VerificationTest[QuditBasis["HoggarSIC"]["Dimension"], 64, TestID -> "HoggarSIC-bare"]

(* tensor-product LIST stays as a list (semantic composition) *)
VerificationTest[QuditBasis["Computational"[{2, 3, 4}]]["Dimension"], 24, TestID -> "Computational-tensor-list"]

(* alternative dispatch: Pauli string of single-letter names *)
VerificationTest[QuditBasis["XYZ"]["Dimension"], 8, TestID -> "PauliString-XYZ"]

VerificationTest[QuditBasis[3]["Dimension"], 3, TestID -> "Integer-3"]

VerificationTest[QuditBasis[{"a", "b"}, {{1, 0}, {0, 1}}]["Dimension"], 2, TestID -> "Explicit-elements"]

EndTestSection[]


BeginTestSection["QuditBasis - properties"]

VerificationTest[QuditBasis["Pauli"]["Qudits"], 1, TestID -> "Pauli-Qudits"]

VerificationTest[Length[QuditBasis["Pauli"]["Names"]], 4, TestID -> "Pauli-NameCount"]

VerificationTest[QuditBasis["Computational"[2]]["ElementDimensions"], {2}, TestID -> "Comp2-ElementDims"]

VerificationTest[QuditBasis["GellMann"[2]]["Dimension"], 4, TestID -> "GellMann-2-dim"]

EndTestSection[]


BeginTestSection["QuditBasis - tensor product"]

VerificationTest[
    QuantumTensorProduct[QuditBasis["Pauli"], QuditBasis["Pauli"]]["Dimension"],
    16,
    TestID -> "Pauli-tensor-Pauli"
]

VerificationTest[
    QuditBasis["Computational"[{2, 3}]]["Dimensions"],
    {2, 3},
    TestID -> "Comp-tensor-dimensions"
]

EndTestSection[]


BeginTestSection["QuditBasis - failure"]

(* The cache wrapper turns InvalidName into a ConfirmationFailed wrapping
   the InvalidName, so check for either. *)
VerificationTest[
    QuditBasis["NotAnActualName"["bar"]],
    Failure["InvalidName" | "ConfirmationFailed", _],
    {QuditBasis::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

(* A bare unrecognized name used to fall through every rule and stay unevaluated,
   so QuditBasis["NotAnActualName"] had head QuditBasis and ["Dimension"] returned
   the unevaluated property lookup, with no message at all. It has to fail the same
   way the call form does. *)
VerificationTest[
    QuditBasis["NotAnActualName"],
    Failure["InvalidName" | "ConfirmationFailed", _],
    {QuditBasis::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-bare-form"
]

(* Nothing unevaluated may reach a numeric property. *)
VerificationTest[
    NumericQ @ QuditBasis["NotAnActualName"]["Dimension"],
    False,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-no-numeric-dimension"
]

(* The guard's open tail matches the multiplicity arity directly, ahead of the
   multiplicity rule, so that form fails too. *)
VerificationTest[
    FailureQ @ QuditBasis["NotAnActualName", 2],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-bare-form-with-multiplicity"
]

(* Near misses are what a user actually types: a truncation, a wrong case, and a
   name belonging to a sibling constructor. One test each, because three copies of
   the same message in a single test trips General::stop. *)
VerificationTest[
    FailureQ @ QuditBasis["Paul"],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-truncated"
]

VerificationTest[
    FailureQ @ QuditBasis["pauli"],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-wrong-case"
]

VerificationTest[
    FailureQ @ QuditBasis["Hadamard"],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-belongs-to-sibling"
]

(* Two letters that are not all Pauli letters are not a Pauli string, and the empty
   string names no basis: neither may slip past the guard as a shorthand. *)
VerificationTest[
    FailureQ @ QuditBasis["AB"],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-non-Pauli-letters"
]

VerificationTest[
    FailureQ @ QuditBasis[""],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-empty-string"
]

(* A "Basis" suffix is stripped before the guard runs, so an unknown name still
   fails; the message names the stripped form. *)
VerificationTest[
    FailureQ @ QuditBasis["NoSuchBasis"],
    True,
    {QuditBasis::invalidName},
    TestID -> "InvalidName-unknown-Basis-suffix"
]

(* A nested basis specification is diagnosed by name: every "Fourier"[___] call form
   is claimed by the recursive rule, so it is the inner "bad-arg" that is rejected,
   not the Fourier arity. *)
VerificationTest[
    QuditBasis["Fourier"["bad-arg"]],
    Failure["InvalidName" | "ConfirmationFailed", _],
    {QuditBasis::invalidName},
    SameTest -> MatchQ,
    TestID -> "InvalidName-nested-basis-arg"
]

(* A registered name whose call form matches no rule is a different failure: the
   arity is wrong, the name is not. *)
VerificationTest[
    QuditBasis["Bell"[3]],
    Failure["InvalidArguments" | "ConfirmationFailed", _],
    {QuditBasis::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-Bell-call-form"
]

EndTestSection[]


BeginTestSection["QuditBasis - names the invalidName guard must not swallow"]

(* Multi-character Pauli strings are split into one qudit per letter by a rule that
   never consults $QuditBasisNames, in both the bare and the dimensioned call form,
   and under multiplicity. *)
VerificationTest[
    {
        QuditBasis["XYZ"]["Dimension"],
        QuditBasis["XX"]["Dimension"],
        QuditBasis["IIZ"]["Qudits"],
        QuditBasis["XYZ"[3]]["Dimension"],
        QuditBasis["XY", 2]["Dimension"]
    },
    {8, 4, 3, 27, 16},
    {},
    TestID -> "Shorthand-Pauli-strings-not-swallowed"
]

(* Names containing "Basis" are rewritten by deleting it, which also never consults
   $QuditBasisNames. The guard holds no exemption for them, the rewrite rule
   outranking it, so "PauliBasis" and "PauliBasis", 2 are what detect a reordering
   that puts the guard first. "XYZBasis" is the composite case: the rewrite hands
   the stripped name back into dispatch, so it arrives as the Pauli string "XYZ",
   and at multiplicity arity it is the Pauli exemption that carries it. *)
VerificationTest[
    {
        QuditBasis["PauliBasis"]["Dimension"],
        QuditBasis["ComputationalBasis"]["Dimension"],
        QuditBasis["BellBasis"]["Dimension"],
        QuditBasis["FourierBasis"[3]]["Dimension"],
        QuditBasis["PauliBasis", 2]["Dimension"],
        QuditBasis["XYZBasis"]["Dimension"],
        QuditBasis["XYZBasis", 2]["Dimension"]
    },
    {4, 2, 4, 3, 16, 8, 64},
    {},
    TestID -> "Shorthand-Basis-suffix-not-swallowed"
]

(* "Properties" is a query on the symbol itself, not a basis name. *)
VerificationTest[
    MemberQ[QuditBasis["Properties"], "Elements"],
    True,
    {},
    TestID -> "Properties-query-not-a-name"
]

(* Every registered name still builds a valid basis in its bare form. Reading the
   registry rather than a copy of it keeps this honest as names are added. *)
VerificationTest[
    AllTrue[
        Wolfram`QuantumFramework`PackageScope`$QuditBasisNames,
        Wolfram`QuantumFramework`PackageScope`QuditBasisQ[QuditBasis[#]] &
    ],
    True,
    {},
    TestID -> "ValidName-every-registered-name-builds"
]

EndTestSection[]
