Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

$v = FieldVariables[];
{{adv}, {av}} = $v;


BeginTestSection["BosonicMatrixElement - Fock basis (default)"]

(* <m|ad a|n> = Sqrt[m] Sqrt[n] delta[m-1, n-1], the number operator with symbolic indices *)
VerificationTest[
    BosonicMatrixElement[{m, n}, adv ** av, $v],
    Sqrt[m] Sqrt[n] KroneckerDelta[m - 1, n - 1],
    TestID -> "BME-Fock-NumberOperator-Symbolic"
]

VerificationTest[
    BosonicMatrixElement[{3, 3}, adv ** av, $v],
    3,
    TestID -> "BME-Fock-NumberOperator-Integer"
]

(* a lowers: <m|a|n> is nonzero only for m = n - 1 *)
VerificationTest[
    BosonicMatrixElement[{m, n}, av, $v],
    Sqrt[n] KroneckerDelta[m, n - 1],
    TestID -> "BME-Fock-Annihilation"
]

(* A c-number is the identity times delta[m, n] *)
VerificationTest[
    BosonicMatrixElement[{m, n}, 7],
    7 KroneckerDelta[m, n],
    TestID -> "BME-Fock-Constant"
]

(* The canonical commutation relation as a matrix element: <m|[a, ad]|n> = delta[m, n].
   Commutator expands to the NonCommutativeMultiply form the reducer expects. *)
VerificationTest[
    BosonicMatrixElement[{m, n}, Commutator[av, adv], $v],
    KroneckerDelta[m, n],
    TestID -> "BME-Fock-CanonicalCommutator"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - coherent basis"]

(* Coherent states are eigenstates of a, so <alpha|a|alpha> = alpha and <alpha|alpha> = 1 *)
VerificationTest[
    Simplify[BosonicMatrixElement[{\[Alpha], \[Alpha]}, av, $v, "Basis" -> "Coherent"],
        \[Alpha] \[Element] Reals],
    \[Alpha],
    TestID -> "BME-Coherent-Annihilation-Diagonal"
]

VerificationTest[
    Simplify[BosonicMatrixElement[{\[Alpha], \[Alpha]}, 1, "Basis" -> "Coherent"],
        \[Alpha] \[Element] Reals],
    1,
    TestID -> "BME-Coherent-Normalization"
]

(* Normal ordering carries over to the coherent basis unchanged *)
VerificationTest[
    Simplify[BosonicMatrixElement[{\[Alpha], \[Alpha]}, Commutator[av, adv], $v, "Basis" -> "Coherent"],
        \[Alpha] \[Element] Reals],
    1,
    TestID -> "BME-Coherent-CanonicalCommutator"
]

(* Off-diagonal overlap <alpha|beta> = Exp[-(alpha - beta)^2/2] for real amplitudes *)
VerificationTest[
    Simplify[BosonicMatrixElement[{\[Alpha], \[Beta]}, 1, "Basis" -> "Coherent"],
        {\[Alpha], \[Beta]} \[Element] Reals],
    Exp[-(\[Alpha] - \[Beta])^2/2],
    TestID -> "BME-Coherent-Overlap"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - cat-code matrix elements"]

(* A codeword is a list of {coefficient, coherent amplitude} pairs, and comb[r, m] is the
   m-component cat comb holding photon number r mod m. These blocks are the Knill-Laflamme
   building blocks for the rotation-symmetric bosonic codes, in closed form and with no
   Fock truncation anywhere. *)

ClearAll[catRaw, catNorm, catElem, comb, catBlock]
catRaw[w1_, w2_, expr_] := Total[Flatten @ Outer[
    Conjugate[#1[[1]]] #2[[1]] BosonicMatrixElement[{#1[[2]], #2[[2]]}, expr, $v, "Basis" -> "Coherent"] &,
    w1, w2, 1]]
catNorm[w_] := Sqrt[catRaw[w, w, 1]]
catElem[w1_, w2_, expr_] := FullSimplify[catRaw[w1, w2, expr]/(catNorm[w1] catNorm[w2]),
    \[Alpha] \[Element] Reals && \[Alpha] > 0]
comb[r_, m_] := Table[{Exp[-2 Pi I r k/m], \[Alpha] Exp[2 Pi I k/m]}, {k, 0, m - 1}]
catBlock[ws_, expr_] := Table[catElem[ws[[i]], ws[[j]], expr], {i, 2}, {j, 2}]

(* Two-legged cat, logical words the even and odd cat. The mean photon numbers differ,
   so the Knill-Laflamme diagonal is not codeword-independent. *)
VerificationTest[
    FullSimplify[catBlock[{comb[0, 2], comb[1, 2]}, adv ** av],
        \[Alpha] \[Element] Reals && \[Alpha] > 0],
    {{\[Alpha]^2 Tanh[\[Alpha]^2], 0}, {0, \[Alpha]^2 Coth[\[Alpha]^2]}},
    TestID -> "BME-Cat2-NumberBlock"
]

(* Four-component cat: a sends the mod-4 combs out of the code space entirely, so the
   single-loss block vanishes identically. This is the structural reason Cat4 corrects a
   single loss in the large-alpha limit while Cat2 does not correct it at any alpha. *)
VerificationTest[
    FullSimplify[catBlock[{comb[0, 4], comb[2, 4]}, av],
        \[Alpha] \[Element] Reals && \[Alpha] > 0],
    {{0, 0}, {0, 0}},
    TestID -> "BME-Cat4-LossBlockVanishes"
]

(* Cat4 mean photon numbers in closed form. Their difference is the only residual
   obstruction, and it decays with an oscillation in alpha^2 as alpha grows. *)
VerificationTest[
    FullSimplify[Diagonal @ catBlock[{comb[0, 4], comb[2, 4]}, adv ** av],
        \[Alpha] \[Element] Reals && \[Alpha] > 0],
    FullSimplify[
        {\[Alpha]^2 (Sinh[\[Alpha]^2] - Sin[\[Alpha]^2])/(Cosh[\[Alpha]^2] + Cos[\[Alpha]^2]),
         \[Alpha]^2 (Sinh[\[Alpha]^2] + Sin[\[Alpha]^2])/(Cosh[\[Alpha]^2] - Cos[\[Alpha]^2])},
        \[Alpha] \[Element] Reals && \[Alpha] > 0],
    TestID -> "BME-Cat4-MeanPhotonNumbers"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - exponentials"]

(* A function of the number operator is diagonal.  This one is the no-jump factor of the
   photon-loss channel, which is not a polynomial and so had no route before. *)
VerificationTest[
    BosonicMatrixElement[{m, k}, \[Eta]^((adv ** av)/2), $v],
    \[Eta]^(k/2) KroneckerDelta[m, k],
    TestID -> "BME-Exp-NumberDiagonal"
]

(* A wing is an exponential of one operator alone, evaluated as a Taylor coefficient of
   Exp[P].  A quadratic wing raises two levels at a time, so odd differences vanish. *)
VerificationTest[
    {BosonicMatrixElement[{4, 0}, Exp[s adv ** adv], $v],
     BosonicMatrixElement[{3, 0}, Exp[s adv ** adv], $v]},
    {Sqrt[6] s^2, 0},
    TestID -> "BME-Exp-QuadraticWingParity"
]

(* An exponential mixing the two operators is disentangled first.  Cross-checked against
   MatrixExp on a 40-level truncation, which gives 1.8682476751142998. *)
VerificationTest[
    Chop[N[BosonicMatrixElement[{6, 2}, Exp[av + adv], $v]] - 1.8682476751142998],
    0,
    TestID -> "BME-Exp-MixedDisentangles"
]

(* The payoff: a squeeze written as its disentangled product of wings agrees with the
   closed form the SqueezeOperator clause gives for the same operator. *)
VerificationTest[
    With[{sd = BosonicExpOrder[
            Exp[(Conjugate[\[Xi]] av ** av - \[Xi] adv ** adv)/2], $v, Assumptions -> \[Xi] > 0] /. \[Xi] -> 0.4},
        Chop[BosonicMatrixElement[{2, 0}, Evaluate[sd], $v] -
             BosonicMatrixElement[{2, 0}, SqueezeOperator[0.4]]]],
    0,
    TestID -> "BME-Exp-SqueezeMatchesNamed"
]

(* The factors are assembled as one normally ordered product, so an anti-normally ordered
   product must not be taken at face value: reordering it while keeping its own scalar is
   wrong by the exponential of the commutator.  It is declined and then disentangled, so
   both spellings of the same operator agree. *)
VerificationTest[
    With[{no = BosonicExpOrder[Exp[0.7 adv - 0.4 av], $v],
          an = BosonicExpOrder[Exp[0.7 adv - 0.4 av], $v, "Ordering" -> "Antinormal"]},
        Chop[BosonicMatrixElement[{3, 2}, Evaluate[an], $v] -
             BosonicMatrixElement[{3, 2}, Evaluate[no], $v]]],
    0,
    TestID -> "BME-Exp-AntinormalNotReordered"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - vanishing powers at the vacuum"]

(* :Exp[-ad a]: = Sum[(-1)^k ad^k a^k/k!] is the vacuum projector |0><0|.  Its diagonal is
   (1 + c)^n at c = -1, and the n = 0 element is the k = 0 term alone, so it is 1, not 0^0.
   The whole n = 0 column is at stake: an off-diagonal zero there must not come out as
   0^0 times a vanishing Kronecker delta. *)
VerificationTest[
    Table[BosonicMatrixElement[{m, n}, NormalOrdered[Exp[-adv ** av], $v], $v], {m, 0, 3}, {n, 0, 3}],
    Table[KroneckerDelta[m, 0] KroneckerDelta[n, 0], {m, 0, 3}, {n, 0, 3}],
    TestID -> "BME-Vacuum-ProjectorFock"
]

(* The closed form for symbolic indices must remain right when the indices reach the vacuum. *)
VerificationTest[
    With[{e = BosonicMatrixElement[{m, n}, NormalOrdered[Exp[-adv ** av], $v], $v]},
        {e /. {m -> 0, n -> 0}, e /. {m -> 2, n -> 2}, e /. {m -> 1, n -> 0}}],
    {1, 0, 0},
    TestID -> "BME-Vacuum-ProjectorSymbolicIndex"
]

(* A machine-precision -1. leaves a base of 0. rather than an exact 0, and the elements
   stay machine numbers. *)
VerificationTest[
    Table[BosonicMatrixElement[{m, n}, NormalOrdered[Exp[-1. adv ** av], $v], $v], {m, 0, 2}, {n, 0, 2}],
    {{1., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
    TestID -> "BME-Vacuum-ProjectorMachineCoefficient"
]

(* 0^(-ad a) is unbounded off the vacuum: it must not be read as the vacuum projector. *)
VerificationTest[
    BosonicMatrixElement[{1, 1}, 0^(-adv ** av), $v],
    ComplexInfinity,
    {Power::infy},
    TestID -> "BME-Vacuum-NegativePowerUnbounded"
]

(* Wings on either side turn the projector into the unnormalized ket-bra
   Exp[z ad]|0><0|Exp[w a], whose elements are z^m w^n / Sqrt[m! n!]: every element,
   not only the vacuum one, passes through the n = 0 term of the middle sum. *)
VerificationTest[
    Simplify[
        Table[BosonicMatrixElement[{m, n}, Exp[z adv] ** NormalOrdered[Exp[-adv ** av], $v] ** Exp[w av], $v],
            {m, 0, 2}, {n, 0, 2}] -
        Table[z^m w^n/Sqrt[m! n!], {m, 0, 2}, {n, 0, 2}]],
    ConstantArray[0, {3, 3}],
    TestID -> "BME-Vacuum-ProjectorBetweenWings"
]

(* Total loss: the no-jump Kraus operator eta^(ad a/2) of the photon-loss channel at
   eta = 0 keeps only the vacuum. *)
VerificationTest[
    Table[BosonicMatrixElement[{m, n}, 0^((adv ** av)/2), $v], {m, 0, 2}, {n, 0, 2}],
    Table[KroneckerDelta[m, 0] KroneckerDelta[n, 0], {m, 0, 2}, {n, 0, 2}],
    TestID -> "BME-Vacuum-FullLossNoJump"
]

(* Between coherent states the projector gives <alpha|0><0|beta>. *)
VerificationTest[
    Simplify[BosonicMatrixElement[{\[Alpha], \[Beta]}, NormalOrdered[Exp[-adv ** av], $v], $v, "Basis" -> "Coherent"] -
        Exp[-(Conjugate[\[Alpha]] \[Alpha] + Conjugate[\[Beta]] \[Beta])/2]],
    0,
    TestID -> "BME-Vacuum-ProjectorCoherent"
]

(* With a zero amplitude the series in Conjugate[alpha] beta keeps only its first term:
   <0|0><0|beta> = <0|beta>, <0|0><0|0> = 1, and Exp[lambda ad a] fixes the vacuum. *)
VerificationTest[
    Simplify[{
        BosonicMatrixElement[{0, \[Beta]}, NormalOrdered[Exp[-adv ** av], $v], $v, "Basis" -> "Coherent"] - Exp[-Conjugate[\[Beta]] \[Beta]/2],
        BosonicMatrixElement[{0, 0}, NormalOrdered[Exp[-adv ** av], $v], $v, "Basis" -> "Coherent"] - 1,
        BosonicMatrixElement[{0, \[Beta]}, Exp[\[Lambda] adv ** av], $v, "Basis" -> "Coherent"] - Exp[-Conjugate[\[Beta]] \[Beta]/2]
    }],
    {0, 0, 0},
    TestID -> "BME-Vacuum-ProjectorCoherentZeroAmplitude"
]

(* The vacuum is the coherent state of amplitude 0, so <0|a|beta> = beta <0|beta>,
   <alpha|ad|0> = Conjugate[alpha] <alpha|0>, and <0|ad a|beta> = 0.  The zero amplitude
   carries the power 0 for every operator it does not touch. *)
VerificationTest[
    Simplify[{
        BosonicMatrixElement[{0, \[Beta]}, av, $v, "Basis" -> "Coherent"] - \[Beta] Exp[-Conjugate[\[Beta]] \[Beta]/2],
        BosonicMatrixElement[{\[Alpha], 0}, adv, $v, "Basis" -> "Coherent"] - Conjugate[\[Alpha]] Exp[-Conjugate[\[Alpha]] \[Alpha]/2],
        BosonicMatrixElement[{0, \[Beta]}, 1 + adv ** av, $v, "Basis" -> "Coherent"] - Exp[-Conjugate[\[Beta]] \[Beta]/2]
    }],
    {0, 0, 0},
    TestID -> "BME-Vacuum-CoherentZeroAmplitude"
]

(* One mode in the vacuum, the other coherent: the product of the per-mode overlaps times
   1 + 0 + Conjugate[alpha] beta. *)
VerificationTest[
    With[{b = FieldVariables[{1, 2}]},
        Simplify[
            BosonicMatrixElement[{{0, \[Alpha]}, {0, \[Beta]}}, 1 + b[[1, 1]] ** b[[2, 1]] + b[[1, 2]] ** b[[2, 2]], b,
                "Basis" -> "Coherent"] -
            (1 + Conjugate[\[Alpha]] \[Beta]) Exp[Conjugate[\[Alpha]] \[Beta] - (Conjugate[\[Alpha]] \[Alpha] + Conjugate[\[Beta]] \[Beta])/2]]],
    0,
    TestID -> "BME-Vacuum-CoherentZeroAmplitudeMultimode"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - independent modes"]

$v2 = FieldVariables[{1, 2}];
{{bd1, bd2}, {b1, b2}} = $v2;

(* Distinct modes commute and the Fock space is their tensor product, so the element of a
   normally ordered monomial is one single-mode element per mode, multiplied. *)
VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1 ** bd2 ** b2, $v2],
    Sqrt[m1] Sqrt[m2] Sqrt[n1] Sqrt[n2] KroneckerDelta[m1 - 1, n1 - 1] KroneckerDelta[m2 - 1, n2 - 1],
    TestID -> "BME-Modes-NumberProduct"
]

VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b2, $v2],
    Sqrt[m1] Sqrt[n2] KroneckerDelta[m1 - 1, n1] KroneckerDelta[m2, n2 - 1],
    TestID -> "BME-Modes-CrossMode"
]

(* Index slot k belongs to mode k of vars however the operator happens to be written. *)
VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd2 ** b1, $v2],
    Sqrt[n1] Sqrt[m2] KroneckerDelta[m1, n1 - 1] KroneckerDelta[m2 - 1, n2],
    TestID -> "BME-Modes-InferredOrderSorted"
]

(* A mode the operator leaves alone contributes the identity, and the declared order is kept. *)
VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1, $v2],
    Sqrt[m1] Sqrt[n1] KroneckerDelta[m1 - 1, n1 - 1] KroneckerDelta[m2, n2],
    TestID -> "BME-Modes-IdleModeDeclared"
]

VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1, {{bd2, bd1}, {b2, b1}}],
    Sqrt[m2] Sqrt[n2] KroneckerDelta[m1, n1] KroneckerDelta[m2 - 1, n2 - 1],
    TestID -> "BME-Modes-DeclaredOrderKept"
]

(* A product of single-mode operators factorizes into their single-mode elements. *)
VerificationTest[
    Simplify[
        BosonicMatrixElement[{{m1, m2}, {n1, n2}}, (b1 + bd1) ** (b1 + bd1) ** bd2 ** b2 ** b2, $v2] -
        BosonicMatrixElement[{m1, n1}, (b1 + bd1) ** (b1 + bd1), {{bd1}, {b1}}] *
            BosonicMatrixElement[{m2, n2}, bd2 ** b2 ** b2, {{bd2}, {b2}}]
    ],
    0,
    TestID -> "BME-Modes-ProductFactorizes"
]

(* Independent of the reduction: a mixed, not normally ordered operator against explicit
   matrices on a 6-level truncation of each mode.  The operator has degree 3, so on
   indices up to 2 no intermediate state leaves the truncation and the block is exact. *)
VerificationTest[
    With[{d = 6, op = (b2 + bd1) ** (b1 + bd2) ** b1},
        With[{low = SparseArray[{i_, j_} /; j == i + 1 :> Sqrt[i], {d, d}], one = IdentityMatrix[d]},
            With[{mat = Normal[op /. NonCommutativeMultiply -> Dot /. {
                    b1 -> KroneckerProduct[low, one], bd1 -> KroneckerProduct[Transpose[low], one],
                    b2 -> KroneckerProduct[one, low], bd2 -> KroneckerProduct[one, Transpose[low]]}]},
                Count[Tuples[Tuples[Range[0, 2], 2], 2],
                    {{i1_, i2_}, {j1_, j2_}} /;
                        Simplify[BosonicMatrixElement[{{i1, i2}, {j1, j2}}, op, $v2] - mat[[i1 d + i2 + 1, j1 d + j2 + 1]]] =!= 0]
            ]
        ]
    ],
    0,
    TestID -> "BME-Modes-MatchesTruncatedMatrices"
]

(* Coherent states factorize too: the overlap is a product of single-mode overlaps. *)
VerificationTest[
    Simplify[
        BosonicMatrixElement[{{\[Alpha]1, \[Alpha]2}, {\[Beta]1, \[Beta]2}}, b1, $v2,
            "Basis" -> "Coherent"] -
        \[Beta]1 Exp[Conjugate[\[Alpha]1] \[Beta]1 - (Conjugate[\[Alpha]1] \[Alpha]1 + Conjugate[\[Beta]1] \[Beta]1)/2] *
            Exp[Conjugate[\[Alpha]2] \[Beta]2 - (Conjugate[\[Alpha]2] \[Alpha]2 + Conjugate[\[Beta]2] \[Beta]2)/2]
    ],
    0,
    TestID -> "BME-Modes-CoherentFactorizes"
]

VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, 7],
    7 KroneckerDelta[m1, n1] KroneckerDelta[m2, n2],
    TestID -> "BME-Modes-Constant"
]

(* Operators of different modes commute, and an operator that cancels only once ordered
   is zero, not a stray term with a negative power. *)
VerificationTest[
    {BosonicMatrixElement[{{m1, m2}, {n1, n2}}, Commutator[b1, bd2], $v2],
     BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1 - b1 ** bd1 + 1, $v2],
     BosonicMatrixElement[{m, n}, adv ** av - av ** adv + 1, $v]},
    {0, 0, 0},
    TestID -> "BME-Modes-CancellationIsZero"
]

(* One mode reads the same with a scalar index, a one-component index vector, or a
   labelled field variable, on every route: polynomial, wing and closed form. *)
VerificationTest[
    {BosonicMatrixElement[{{m}, {n}}, adv ** av, $v],
     BosonicMatrixElement[{{2}, {0}}, Exp[z adv], $v],
     BosonicMatrixElement[{{1}, {0}}, DisplacementOperator[\[Alpha]]],
     BosonicMatrixElement[{m, n}, bd1 ** b1, {{bd1}, {b1}}],
     BosonicMatrixElement[{2, 0}, Exp[z bd1], {{bd1}, {b1}}]},
    {Sqrt[m] Sqrt[n] KroneckerDelta[m - 1, n - 1], z^2/Sqrt[2], \[Alpha] E^(-Abs[\[Alpha]]^2/2),
     Sqrt[m] Sqrt[n] KroneckerDelta[m - 1, n - 1], z^2/Sqrt[2]},
    TestID -> "BME-Modes-OneModeSpellings"
]

(* An exponential in one mode of several: the other mode contributes its overlap. *)
VerificationTest[
    {BosonicMatrixElement[{{2, 1}, {0, 0}}, Exp[z bd1], $v2], BosonicMatrixElement[{{2, 0}, {0, 0}}, Exp[z bd1], $v2]},
    {0, z^2/Sqrt[2]},
    TestID -> "BME-Modes-ExponentialOneMode"
]

(* A function of the number operators alone is diagonal in the multimode Fock basis,
   <m|f(N1, N2)|n> = f(n) delta_{mn}: the dispersive rotation e^{-i chi t (N1 + N3 (N1 - N2))}
   gives e^{-i chi t (n + 1)} with mode 1 excited, e^{i chi t n} with mode 2 excited, and 0
   between the two; products and powers of the N_k are read on the same indices. *)
VerificationTest[
    With[{u = Exp[-I chi t (bd1 ** b1 + bd3 ** b3 ** (bd1 ** b1 - bd2 ** b2))], v3 = {{bd1, bd2, bd3}, {b1, b2, b3}}},
        {BosonicMatrixElement[{{1, 0, n}, {1, 0, n}}, u, v3],
         BosonicMatrixElement[{{0, 1, n}, {0, 1, n}}, u, v3],
         BosonicMatrixElement[{{1, 0, n}, {0, 1, n}}, u, v3],
         BosonicMatrixElement[{{2, 3}, {2, 3}}, bd1 ** b1 ** bd2 ** b2 + (bd1 ** b1)^2 Exp[c bd2 ** b2], {{bd1, bd2}, {b1, b2}}]}
    ],
    {E^(-I chi (1 + n) t), E^(I chi n t), 0, 6 + 4 E^(3 c)},
    TestID -> "BME-Modes-NumberFunctionDiagonal"
]

(* The trace of a product of per-mode Boltzmann factors factorizes into geometric series,
   Tr x1^N1 x2^N2 = 1/((1 - x1)(1 - x2)) = det(I - diag(x1, x2))^-1. *)
VerificationTest[
    Simplify[Sum[BosonicMatrixElement[{{k1, k2}, {k1, k2}}, Exp[Log[x1] bd1 ** b1 + Log[x2] bd2 ** b2], {{bd1, bd2}, {b1, b2}}],
        {k1, 0, Infinity}, {k2, 0, Infinity}], 0 < x1 < 1 && 0 < x2 < 1],
    1/((1 - x1) (1 - x2)),
    SameTest -> (Simplify[#1 - #2] === 0 &),
    TestID -> "BME-Modes-NumberFunctionTrace"
]

(* Assumptions reach the reordering: the squeeze exponential's element comes out in hyperbolic
   functions, <2|S(r)|0> = -Tanh[r] Sqrt[Sech[r]]/Sqrt[2], as an option or through Assuming. *)
VerificationTest[
    With[{sq = Exp[r/2 (b ** b - bd ** bd)]},
        {BosonicMatrixElement[{2, 0}, sq, {{bd}, {b}}, Assumptions -> r > 0], Assuming[r > 0, BosonicMatrixElement[{2, 0}, sq, {{bd}, {b}}]]}
    ],
    {-(Sech[r]^(3/2) Sinh[r])/Sqrt[2], -(Sech[r]^(3/2) Sinh[r])/Sqrt[2]},
    TestID -> "BME-Assumptions-SqueezeElement"
]

(* A creator right of an annihilator wing moves across shifted, e^{c a} a^dagger = (a^dagger + c) e^{c a},
   and the element is linear over the resulting sum, in both bases. *)
VerificationTest[
    With[{v = {{bd}, {b}}, ov = Exp[Conjugate[al] be - (Conjugate[al] al + Conjugate[be] be)/2]},
        {
            Simplify[BosonicMatrixElement[{al, be}, Exp[c b] ** bd, v, "Basis" -> "Coherent"] - (Conjugate[al] + c) Exp[c be] ov],
            Simplify[BosonicMatrixElement[{al, be}, Exp[c b] + bd, v, "Basis" -> "Coherent"] - (Exp[c be] + Conjugate[al]) ov],
            Union[Flatten[Table[Simplify[BosonicMatrixElement[{i, k}, Exp[c b] ** bd, v] -
                BosonicMatrixElement[{i, k}, bd ** Exp[c b], v] - c BosonicMatrixElement[{i, k}, Exp[c b], v]], {i, 0, 3}, {k, 0, 3}]]]
        }
    ],
    {0, 0, {0}},
    TestID -> "BME-LadderPastWing-Linear"
]

(* <beta|D(g) a^dagger|beta> = (Conjugate[beta] - Conjugate[g]) <beta|D(g)|beta>, from D a^dagger = (a^dagger - Conjugate[g]) D;
   an unordered word after D reads as its normal-ordered form. *)
VerificationTest[
    With[{v = {{bd}, {b}}, d = Exp[g bd - Conjugate[g] b]},
        {
            Simplify[BosonicMatrixElement[{be, be}, d ** bd, v, "Basis" -> "Coherent"] -
                (Conjugate[be] - Conjugate[g]) BosonicMatrixElement[{be, be}, d, v, "Basis" -> "Coherent"]],
            Simplify[BosonicMatrixElement[{be, be}, d ** b ** bd, v, "Basis" -> "Coherent"] -
                BosonicMatrixElement[{be, be}, d ** BosonicNormalOrder[b ** bd, v], v, "Basis" -> "Coherent"]]
        }
    ],
    {0, 0},
    TestID -> "BME-DisplacementTimesCreatorCoherent"
]

(* Multimode wings: e^{z a1^dagger a2^dagger}|0,0> = Sum z^k |k,k>, and the two-mode squeezed vacuum
   <k,k|S2(r e^{i th})|0,0> = e^{i k th} Tanh[r]^k Sech[r] after BosonicExpOrder disentangles S2. *)
VerificationTest[
    With[{v2 = {{bd1, bd2}, {b1, b2}}, s2 = Exp[r E^(I th) bd1 ** bd2 - r E^(-I th) b1 ** b2], as = r > 0 && Element[th, Reals]},
        {
            Table[BosonicMatrixElement[{{m1, m2}, {0, 0}}, Exp[z bd1 ** bd2], v2], {m1, 0, 2}, {m2, 0, 2}],
            BosonicMatrixElement[{{1, 1}, {0, 0}}, Exp[z1 bd1] ** Exp[z2 bd2], v2],
            Table[Simplify[BosonicMatrixElement[{{k, k}, {0, 0}}, s2, v2, Assumptions -> as] - E^(I k th) Tanh[r]^k Sech[r], as], {k, 0, 3}],
            BosonicMatrixElement[{{1, 2}, {0, 0}}, s2, v2, Assumptions -> as]
        }
    ],
    {{{1, 0, 0}, {0, z, 0}, {0, 0, z^2}}, z1 z2, {0, 0, 0, 0}, 0},
    TestID -> "BME-Modes-WingFock"
]

(* Both wings and a middle: <1,1|e^{r(a1^dagger a2^dagger - a1 a2)}|1,1> agrees through r^6 with the
   term-by-term polynomial expansion, which uses no exponential code. *)
VerificationTest[
    With[{v2 = {{bd1, bd2}, {b1, b2}}, K = bd1 ** bd2 - b1 ** b2},
        Simplify[Normal[Series[BosonicMatrixElement[{{1, 1}, {1, 1}}, Exp[r K], v2, Assumptions -> r > 0], {r, 0, 6}]] -
            Sum[r^j/j! BosonicMatrixElement[{{1, 1}, {1, 1}}, If[j == 0, 1, NonCommutativeMultiply @@ ConstantArray[K, j]], v2], {j, 0, 6}]]
    ],
    0,
    TestID -> "BME-Modes-WingSeries"
]

(* BosonicVEV of exponentials goes through the vacuum element, in one or several modes. *)
VerificationTest[
    {
        Simplify[BosonicVEV[Exp[al bd - Conjugate[al] b], {{bd}, {b}}] - Exp[-al Conjugate[al]/2]],
        BosonicVEV[Exp[r (bd1 ** bd2 - b1 ** b2)], {{bd1, bd2}, {b1, b2}}, Assumptions -> r > 0]
    },
    {0, Sech[r]},
    TestID -> "BME-VEV-Exponential"
]

(* A power of a sum above NonCommutativeExpand's "ExpandMaxExponent" (30) is still ordered:
   <0|(b + b†)^(2j)|0> = (2j - 1)!!, the odd ones vanish. *)
VerificationTest[
    Table[BosonicVEV[GeneralizedPower[NonCommutativeMultiply, b + bd, n], {{bd}, {b}}], {n, 31, 34}],
    {0, 31!!, 0, 33!!},
    TestID -> "BME-VEV-PowerAboveExpandCap"
]

VerificationTest[
    Expand[BosonicNormalOrder[GeneralizedPower[NonCommutativeMultiply, b + bd, 31], {{bd}, {b}}] -
        BosonicNormalOrder[BosonicNormalOrder[GeneralizedPower[NonCommutativeMultiply, b + bd, 30], {{bd}, {b}}] ** (b + bd), {{bd}, {b}}]],
    0,
    TestID -> "BME-NormalOrder-PowerAboveExpandCap"
]

EndTestSection[]
