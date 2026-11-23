Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

$v = FieldVariables[];
{{adv}, {av}} = $v;
$v2 = FieldVariables[{1, 2}];
{{adv1, adv2}, {av1, av2}} = $v2;
$p[x_, k_] := GeneralizedPower[NonCommutativeMultiply, x, k];
$proj[n_] := $p[adv, n] ** NormalOrdered[Exp[-adv ** av], $v] ** $p[av, n] / n!;


BeginTestSection["SOrderedSymbol"]

(* D(g) = e^(g a^dagger - gbar a) has s-symbol e^(-s |g|^2/2 + g alphabar - gbar alpha) *)
VerificationTest[
    Simplify[SOrderedSymbol[Exp[g adv - gc av], $v, al, s] - Exp[-s g gc/2 + g Conjugate[al] - gc al]],
    0,
    TestID -> "SOrderedSymbol-Displacement-AllS"
]

(* s = 1 and s = -1 are the normal and anti-normal orderings with a^dagger -> alphabar, a -> alpha *)
VerificationTest[
    With[{p = av ** av ** adv ** adv ** av + x adv ** av ** adv, toAl = # //. {GeneralizedPower[NonCommutativeMultiply, b_, k_] :> b^k, NonCommutativeMultiply -> Times} /. {adv -> Conjugate[al], av -> al} &},
        Expand[{SOrderedSymbol[p, $v, al, 1] - toAl[BosonicNormalOrder[p, $v]], SOrderedSymbol[p, $v, al, -1] - toAl[BosonicAntinormalOrder[p, $v]]}]
    ],
    {0, 0},
    TestID -> "SOrderedSymbol-NormalAntinormalEndpoints"
]

(* Weyl symbol of e^(-tau (N + 1/2)) is sech(tau/2) e^(-2 tanh(tau/2) |alpha|^2) *)
VerificationTest[
    Simplify[SOrderedSymbol[Exp[-tau (adv ** av + 1/2)], $v, al, 0, Assumptions -> tau > 0] - Sech[tau/2] Exp[-2 Tanh[tau/2] Conjugate[al] al], tau > 0],
    0,
    TestID -> "SOrderedSymbol-Thermal-Weyl"
]

(* Cahill-Glauber: the t-symbol of |n><n| is 2/(1+t) ((t-1)/(t+1))^n e^(-2|alpha|^2/(1+t)) L_n(4|alpha|^2/(1-t^2)) *)
VerificationTest[
    Table[Simplify[SOrderedSymbol[$proj[n], $v, al, t] - 2/(1 + t) ((t - 1)/(t + 1))^n Exp[-2 Conjugate[al] al/(1 + t)] LaguerreL[n, 4 Conjugate[al] al/(1 - t^2)], -1 < t < 1], {n, 0, 3}],
    {0, 0, 0, 0},
    TestID -> "SOrderedSymbol-FockProjector-AllS"
]

(* squeeze operator: the Gaussian smoothing of the Weyl symbol is the normal symbol *)
VerificationTest[
    With[{sq = Exp[r/2 ($p[av, 2] - $p[adv, 2])]},
        Simplify[SOrderedSymbol[sq, $v, al, 0, Assumptions -> r > 0] - Sech[r/2] Exp[Tanh[r/2] (al^2 - Conjugate[al]^2)], r > 0]
    ],
    0,
    TestID -> "SOrderedSymbol-Squeeze-Weyl"
]

(* the Weyl symbol of a commutator is the Poisson bracket for a quadratic Hamiltonian *)
VerificationTest[
    With[{h = adv ** av + k (av ** av + adv ** adv), b = adv ** adv ** av + nu av ** av ** av, w = SOrderedSymbol[#, $v, al, 0] /. Conjugate[al] -> alc &},
        Expand[w[BosonicNormalOrder[h ** b - b ** h, $v]] - (D[w[h], al] D[w[b], alc] - D[w[h], alc] D[w[b], al])]
    ],
    0,
    TestID -> "SOrderedSymbol-Commutator-PoissonBracket"
]

(* two modes: a product of displacements, a two-mode squeezer, and two oscillators in real time at t = 3, past Pi/2,
   where the root of the whole determinant flips sign; sec(t/2)^2 = 4 e^(i t)/(1 + e^(i t))^2, tan(t/2) = -i (e^(i t) - 1)/(e^(i t) + 1) *)
VerificationTest[
    Simplify[{
        SOrderedSymbol[Exp[g adv1 - gc av1] ** Exp[h adv2 - hc av2], $v2, {al1, al2}, 0] - Exp[g Conjugate[al1] - gc al1 + h Conjugate[al2] - hc al2],
        SOrderedSymbol[Exp[r (av1 ** av2 - adv1 ** adv2)], $v2, {al1, al2}, 0, Assumptions -> r > 0] - Sech[r/2]^2 Exp[2 Tanh[r/2] (al1 al2 - Conjugate[al1] Conjugate[al2])],
        SOrderedSymbol[Exp[-3 I (adv1 ** av1 + adv2 ** av2 + 1)], $v2, {al1, al2}, 0] -
            4 E^(3 I)/(1 + E^(3 I))^2 Exp[-2 (E^(3 I) - 1)/(E^(3 I) + 1) (Conjugate[al1] al1 + Conjugate[al2] al2)]
    }, r > 0],
    {0, 0, 0},
    TestID -> "SOrderedSymbol-TwoModes"
]

(* the anti-normal symbol of a Fock projector is a distribution, and a non-Gaussian normal symbol is out of scope *)
VerificationTest[
    {SOrderedSymbol[$proj[2], $v, al, -1], SOrderedSymbol[Exp[k adv ** av ** adv], $v, al, 0]},
    {_SOrderedSymbol, _SOrderedSymbol},
    {SOrderedSymbol::sing, SOrderedSymbol::nosym},
    SameTest -> MatchQ,
    TestID -> "SOrderedSymbol-Messages"
]

EndTestSection[]
