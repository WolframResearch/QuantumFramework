Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

{av, adv} = FieldVariables[];


BeginTestSection["BosonicMatrixElement - Fock basis (default)"]

(* <m|ad a|n> = Sqrt[m] Sqrt[n] delta[m-1, n-1], the number operator with symbolic indices *)
VerificationTest[
    BosonicMatrixElement[{m, n}, adv ** av],
    Sqrt[m] Sqrt[n] KroneckerDelta[m - 1, n - 1],
    TestID -> "BME-Fock-NumberOperator-Symbolic"
]

VerificationTest[
    BosonicMatrixElement[{3, 3}, adv ** av],
    3,
    TestID -> "BME-Fock-NumberOperator-Integer"
]

(* a lowers: <m|a|n> is nonzero only for m = n - 1 *)
VerificationTest[
    BosonicMatrixElement[{m, n}, av],
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
    BosonicMatrixElement[{m, n}, Commutator[av, adv]],
    KroneckerDelta[m, n],
    TestID -> "BME-Fock-CanonicalCommutator"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - coherent basis"]

(* Coherent states are eigenstates of a, so <alpha|a|alpha> = alpha and <alpha|alpha> = 1 *)
VerificationTest[
    Simplify[BosonicMatrixElement[{\[Alpha], \[Alpha]}, av, "Basis" -> "Coherent"],
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
    Simplify[BosonicMatrixElement[{\[Alpha], \[Alpha]}, Commutator[av, adv], "Basis" -> "Coherent"],
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
    Conjugate[#1[[1]]] #2[[1]] BosonicMatrixElement[{#1[[2]], #2[[2]]}, expr, "Basis" -> "Coherent"] &,
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
    BosonicMatrixElement[{m, k}, \[Eta]^((adv ** av)/2)],
    \[Eta]^(k/2) KroneckerDelta[m, k],
    TestID -> "BME-Exp-NumberDiagonal"
]

(* A wing is an exponential of one operator alone, evaluated as a Taylor coefficient of
   Exp[P].  A quadratic wing raises two levels at a time, so odd differences vanish. *)
VerificationTest[
    {BosonicMatrixElement[{4, 0}, Exp[s adv ** adv]],
     BosonicMatrixElement[{3, 0}, Exp[s adv ** adv]]},
    {Sqrt[6] s^2, 0},
    TestID -> "BME-Exp-QuadraticWingParity"
]

(* An exponential mixing the two operators is disentangled first.  Cross-checked against
   MatrixExp on a 40-level truncation, which gives 1.8682476751142998. *)
VerificationTest[
    Chop[N[BosonicMatrixElement[{6, 2}, Exp[av + adv]]] - 1.8682476751142998],
    0,
    TestID -> "BME-Exp-MixedDisentangles"
]

(* The payoff: a squeeze written as its disentangled product of wings agrees with the
   closed form the SqueezeOperator clause gives for the same operator. *)
VerificationTest[
    With[{sd = BosonicExpOrder[
            Exp[(Conjugate[\[Xi]] av ** av - \[Xi] adv ** adv)/2], \[Xi] > 0] /. \[Xi] -> 0.4},
        Chop[BosonicMatrixElement[{2, 0}, Evaluate[sd]] -
             BosonicMatrixElement[{2, 0}, SqueezeOperator[0.4]]]],
    0,
    TestID -> "BME-Exp-SqueezeMatchesNamed"
]

(* The factors are assembled as one normally ordered product, so an anti-normally ordered
   product must not be taken at face value: reordering it while keeping its own scalar is
   wrong by the exponential of the commutator.  It is declined and then disentangled, so
   both spellings of the same operator agree. *)
VerificationTest[
    With[{no = BosonicExpOrder[Exp[0.7 adv - 0.4 av]],
          an = BosonicExpOrder[Exp[0.7 adv - 0.4 av], "Ordering" -> "Antinormal"]},
        Chop[BosonicMatrixElement[{3, 2}, Evaluate[an]] -
             BosonicMatrixElement[{3, 2}, Evaluate[no]]]],
    0,
    TestID -> "BME-Exp-AntinormalNotReordered"
]

EndTestSection[]


BeginTestSection["BosonicMatrixElement - independent modes"]

{b1, bd1, b2, bd2} = FieldVariables[{1, 2}];

(* Distinct modes commute and the Fock space is their tensor product, so the element of a
   normally ordered monomial is one single-mode element per mode, multiplied. *)
VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1 ** bd2 ** b2],
    Sqrt[m1] Sqrt[m2] Sqrt[n1] Sqrt[n2] KroneckerDelta[m1 - 1, n1 - 1] KroneckerDelta[m2 - 1, n2 - 1],
    TestID -> "BME-Modes-NumberProduct"
]

VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b2],
    Sqrt[m1] Sqrt[n2] KroneckerDelta[m1 - 1, n1] KroneckerDelta[m2, n2 - 1],
    TestID -> "BME-Modes-CrossMode"
]

(* Inferred modes are sorted, so index slot 1 belongs to the first label however the
   operator happens to be written. *)
VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd2 ** b1],
    Sqrt[n1] Sqrt[m2] KroneckerDelta[m1, n1 - 1] KroneckerDelta[m2 - 1, n2],
    TestID -> "BME-Modes-InferredOrderSorted"
]

(* An operator that leaves a mode alone cannot name it, so the modes are declared; the idle
   mode then contributes the identity.  Declared generators keep the order given. *)
VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1, "Generators" -> {b1, bd1, b2, bd2}],
    Sqrt[m1] Sqrt[n1] KroneckerDelta[m1 - 1, n1 - 1] KroneckerDelta[m2, n2],
    TestID -> "BME-Modes-IdleModeDeclared"
]

VerificationTest[
    BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1, "Generators" -> {b2, bd2, b1, bd1}],
    Sqrt[m2] Sqrt[n2] KroneckerDelta[m1, n1] KroneckerDelta[m2 - 1, n2 - 1],
    TestID -> "BME-Modes-DeclaredOrderKept"
]

(* A product of single-mode operators factorizes into their single-mode elements. *)
VerificationTest[
    Simplify[
        BosonicMatrixElement[{{m1, m2}, {n1, n2}}, (b1 + bd1) ** (b1 + bd1) ** bd2 ** b2 ** b2] -
        BosonicMatrixElement[{m1, n1}, (b1 + bd1) ** (b1 + bd1), "Generators" -> {b1, bd1}] *
            BosonicMatrixElement[{m2, n2}, bd2 ** b2 ** b2, "Generators" -> {b2, bd2}]
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
                        Simplify[BosonicMatrixElement[{{i1, i2}, {j1, j2}}, op] - mat[[i1 d + i2 + 1, j1 d + j2 + 1]]] =!= 0]
            ]
        ]
    ],
    0,
    TestID -> "BME-Modes-MatchesTruncatedMatrices"
]

(* Coherent states factorize too: the overlap is a product of single-mode overlaps. *)
VerificationTest[
    Simplify[
        BosonicMatrixElement[{{\[Alpha]1, \[Alpha]2}, {\[Beta]1, \[Beta]2}}, b1,
            "Basis" -> "Coherent", "Generators" -> {b1, bd1, b2, bd2}] -
        \[Beta]1 Exp[Conjugate[\[Alpha]1] \[Beta]1 - (Abs[\[Alpha]1]^2 + Abs[\[Beta]1]^2)/2] *
            Exp[Conjugate[\[Alpha]2] \[Beta]2 - (Abs[\[Alpha]2]^2 + Abs[\[Beta]2]^2)/2]
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
    {BosonicMatrixElement[{{m1, m2}, {n1, n2}}, Commutator[b1, bd2], "Generators" -> {b1, bd1, b2, bd2}],
     BosonicMatrixElement[{{m1, m2}, {n1, n2}}, bd1 ** b1 - b1 ** bd1 + 1, "Generators" -> {b1, bd1, b2, bd2}],
     BosonicMatrixElement[{m, n}, adv ** av - av ** adv + 1]},
    {0, 0, 0},
    TestID -> "BME-Modes-CancellationIsZero"
]

(* One mode reads the same with a scalar index, a one-component index vector, or a
   labelled field variable, on every route: polynomial, wing and closed form. *)
VerificationTest[
    {BosonicMatrixElement[{{m}, {n}}, adv ** av],
     BosonicMatrixElement[{{2}, {0}}, Exp[z adv]],
     BosonicMatrixElement[{{1}, {0}}, DisplacementOperator[\[Alpha]]],
     BosonicMatrixElement[{m, n}, bd1 ** b1],
     BosonicMatrixElement[{2, 0}, Exp[z bd1]]},
    {Sqrt[m] Sqrt[n] KroneckerDelta[m - 1, n - 1], z^2/Sqrt[2], \[Alpha] E^(-Abs[\[Alpha]]^2/2),
     Sqrt[m] Sqrt[n] KroneckerDelta[m - 1, n - 1], z^2/Sqrt[2]},
    TestID -> "BME-Modes-OneModeSpellings"
]

(* Exponentials are still read one mode at a time; several modes stay unevaluated rather
   than feed an index vector to the single-mode wing formula. *)
VerificationTest[
    Head[BosonicMatrixElement[{{2, 1}, {0, 0}}, Exp[z bd1]]],
    BosonicMatrixElement,
    TestID -> "BME-Modes-ExponentialUnevaluated"
]

EndTestSection[]
