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

EndTestSection[]
