Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]


BeginTestSection["SOrderedFunction - reductions to the named members"]

(* s = 0 must reproduce WignerFunction and s = -1 HusimiQFunction. These are the oracle for
   the general kernel: both named functions were written independently of it, so an index or
   conjugation error in SOrderedKernel shows up here and nowhere else. *)
VerificationTest[
    Chop[SOrderedFunction[FockState[1, 6], 0.53 - 0.31 I, 0] -
        WignerFunction[FockState[1, 6], 0.53 - 0.31 I], 10^-10],
    0,
    TestID -> "SOrd-Reduces-Wigner"
]

VerificationTest[
    Chop[SOrderedFunction[FockState[1, 6], 0.53 - 0.31 I, -1] -
        HusimiQFunction[FockState[1, 6], 0.53 - 0.31 I], 10^-10],
    0,
    TestID -> "SOrd-Reduces-Husimi"
]

(* The s == -1 branch is selected by TrueQ[s == -1], so a machine real must reach it too:
   the general form is Indeterminate there, not merely inaccurate. *)
VerificationTest[
    Chop[SOrderedFunction[FockState[1, 6], 0.53 - 0.31 I, -1.] -
        HusimiQFunction[FockState[1, 6], 0.53 - 0.31 I], 10^-10],
    0,
    TestID -> "SOrd-Reduces-Husimi-MachineReal"
]

(* A state with coherences: the diagonal cases above cannot catch a swapped alpha and
   alphaStar in the kernel. *)
VerificationTest[
    Chop[SOrderedFunction[QuantumState[{1, I, 0, 0, 0, 0}/Sqrt[2], 6], 0.53 - 0.31 I, 0] -
        WignerFunction[QuantumState[{1, I, 0, 0, 0, 0}/Sqrt[2], 6], 0.53 - 0.31 I], 10^-10],
    0,
    TestID -> "SOrd-Reduces-Wigner-Coherences"
]

EndTestSection[]


BeginTestSection["SOrderedFunction - closed forms and exact identities"]

VerificationTest[
    Simplify @ SOrderedFunction[QuantumState[{0, 1, 0}, 3], {x, p}, s],
    -((E^((p^2 + x^2)/(-1 + s)) (-1 + 2 p^2 + s^2 + 2 x^2))/(Pi (-1 + s)^3)),
    TestID -> "SOrd-Fock1-ClosedForm-SymbolicS"
]

(* Normalization holds for every s < 1, which is what makes the s -> 1 blow-up a
   distributional limit rather than a defect. *)
VerificationTest[
    Simplify @ Integrate[
        SOrderedFunction[QuantumState[{0, 1, 0}, 3], {u, v}, s],
        {u, -Infinity, Infinity}, {v, -Infinity, Infinity},
        Assumptions -> {Element[{u, v}, Reals], s < 1}],
    1,
    TestID -> "SOrd-Fock1-Normalized-SymbolicS"
]

(* Interpolates between orderings: 2 at s = -1 (antinormal), 3/2 at s = 0 (symmetric), and 1
   in the s -> 1 limit, which is <a^dagger a> for |1>. *)
VerificationTest[
    Simplify @ Integrate[
        (u^2 + v^2)/2 SOrderedFunction[QuantumState[{0, 1, 0}, 3], {u, v}, s],
        {u, -Infinity, Infinity}, {v, -Infinity, Infinity},
        Assumptions -> {Element[{u, v}, Reals], s < 1}],
    (3 - s)/2,
    TestID -> "SOrd-Fock1-SecondMoment-SymbolicS"
]

(* Vacuum stays a positive Gaussian for every s, of width (1-s)/4. *)
VerificationTest[
    FullSimplify @ SOrderedFunction[FockState[0, 4], alpha, s],
    -((2 E^((2 alpha Conjugate[alpha])/(-1 + s)))/(Pi (-1 + s))),
    TestID -> "SOrd-Vacuum-ClosedForm-SymbolicS"
]

EndTestSection[]


BeginTestSection["SOrderedFunction - guards"]

(* The phase-space origin raises 0^0 in the kernel unless ZeroLimitPower intervenes. *)
VerificationTest[
    NumericQ @ SOrderedFunction[FockState[1, 6], 0, -1/2],
    True,
    TestID -> "SOrd-Origin-Finite"
]

VerificationTest[
    SOrderedFunction[FockState[1, 6], 0.5, 1],
    $Failed,
    {SOrderedFunction::sunit},
    TestID -> "SOrd-RejectsSEqualOne"
]

VerificationTest[
    SOrderedFunction[FockState[1, 6], 0.5, 1.5],
    $Failed,
    {SOrderedFunction::snum},
    TestID -> "SOrd-RejectsSAboveOne"
]

VerificationTest[
    SOrderedFunction[QuantumState["00"], 0.5, 0],
    $Failed,
    {SOrderedFunction::multimode},
    TestID -> "SOrd-RejectsMultimode"
]

EndTestSection[]


BeginTestSection["SOrderedRepresentation"]

(* -1/Pi is the closed-form origin value for |1> at s = 0 in the {x,p} normalization, half
   the -2/Pi carried by the alpha convention. *)
VerificationTest[
    Round[SOrderedRepresentation[FockState[1, 6], {-5, 5}, {-5, 5}, 0][0, 0], 0.001],
    Round[-1/Pi, 0.001],
    TestID -> "SOrdRep-Fock1-Origin-Wigner"
]

(* General origin value is -(1+s)/(Pi (1-s)^2) over dx dp. *)
VerificationTest[
    Round[SOrderedRepresentation[FockState[1, 6], {-5, 5}, {-5, 5}, 1/2][0, 0], 0.01],
    Round[-(1 + 1/2)/(Pi (1 - 1/2)^2), 0.01],
    TestID -> "SOrdRep-Fock1-Origin-SHalf"
]

(* Husimi Q of |1> vanishes at the origin. *)
VerificationTest[
    Abs[SOrderedRepresentation[FockState[1, 6], {-5, 5}, {-5, 5}, -1][0, 0]] < 10^-4,
    True,
    TestID -> "SOrdRep-Fock1-Origin-Husimi"
]

VerificationTest[
    Head @ SOrderedRepresentation[FockState[2, 6], {-4, 4}, {-4, 4}, -1/2, "GridSize" -> 40],
    InterpolatingFunction,
    TestID -> "SOrdRep-ReturnsInterpolatingFunction"
]

VerificationTest[
    SOrderedRepresentation[QuantumState["00"], {-5, 5}, {-5, 5}, 0],
    $Failed,
    {SOrderedRepresentation::multimode},
    TestID -> "SOrdRep-RejectsMultimode"
]

VerificationTest[
    SOrderedRepresentation[FockState[1, 6], {-5, 5}, {-5, 5}, 1],
    $Failed,
    {SOrderedRepresentation::sunit},
    TestID -> "SOrdRep-RejectsSEqualOne"
]

(* Structure of scale Sqrt[1-s] = 0.1 against a spacing of about 0.53. *)
VerificationTest[
    Head @ SOrderedRepresentation[FockState[1, 6], {-5, 5}, {-5, 5}, 99/100, "GridSize" -> 20],
    InterpolatingFunction,
    {SOrderedRepresentation::coarse},
    TestID -> "SOrdRep-WarnsOnCoarseGrid"
]

EndTestSection[]
