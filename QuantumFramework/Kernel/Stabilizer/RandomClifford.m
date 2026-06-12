Package["Wolfram`QuantumFramework`"]

PackageScope[RandomClifford]
PackageScope[SampleMallows]



(* ============================================================================ *)
(* Mallows-distribution sampler for symmetric-group permutations.               *)
(* Used by the Bravyi-Maslov / Koenig-Smolin random Clifford construction.      *)
(* Reference: Koenig & Smolin 2014, arxiv:1406.2170, Section 3.2.               *)
(* ============================================================================ *)

SampleMallows[n_] := Block[{h = ConstantArray[0, n], perm = ConstantArray[0, n], indices = Range[n]},
    Do[
        Block[{m = n - i, eps, r, index, k},
            eps = 4 ^ (- m);
            r = RandomReal[];
            index = - Ceiling[Log2[r + (1 - r) eps]];
            h[[i + 1]] = Boole[index < m];
            k = If[index < m, index, 2 m - index - 1];
            perm[[i + 1]] = indices[[k + 1]];
            indices = Drop[indices, {Mod[k, Length[indices]] + 1}]
        ],
        {i, 0, n - 1}
    ];
    {h, perm}
]


(* Plain packed integer matrices throughout: structured-array wrappers          *)
(* (SymmetrizedArray, LowerTriangularMatrix) make the downstream Dot fall off   *)
(* the packed-matmul path onto a generic tensor route whose intermediates       *)
(* reach tens of GB at n ~ 500 and kill the kernel.                             *)
fillTril[n_, symmetric_ : False] := With[{lt = LowerTriangularize[RandomInteger[1, {n, n}], If[symmetric, 0, -1]]},
    If[symmetric,
        lt + Transpose[LowerTriangularize[lt, -1]],
        lt + IdentityMatrix[n]
    ]
]



(* ============================================================================ *)
(* RandomClifford[n] -- sample uniformly from the n-qubit Clifford group.       *)
(* Cardinality |C_n| = 2^(n^2 + 2n) Product[4^j - 1, {j, 1, n}] (KoeSmo14 Eq 2)*)
(* ============================================================================ *)

RandomClifford[n_] := Block[{h, perm, gamma1, delta1, gamma2, delta2, zero, prod1, prod2, inv1, inv2, table1, table2, table, indices},
    {h, perm} = SampleMallows[n];
    gamma1 = fillTril[n, True];
    gamma2 = fillTril[n, True];
    delta1 = fillTril[n];
    delta2 = fillTril[n];
    zero = ConstantArray[0, {n, n}];
    prod1 = Mod[gamma1 . delta1, 2];
    prod2 = Mod[gamma2 . delta2, 2];
    (* The tableau lives over F_2 (the final product is reduced mod 2), so      *)
    (* invert mod 2: the exact integer inverse of a unit lower-triangular 0/1   *)
    (* matrix has entries that grow exponentially with n, turning the final     *)
    (* 2n x 2n Dot into a big-integer matmul.                                   *)
    inv1 = Transpose[Inverse[delta1, Modulus -> 2]];
    inv2 = Transpose[Inverse[delta2, Modulus -> 2]];
    table1 = Join[Join[delta1, zero, 2], Join[prod1, inv1, 2]];
    table2 = Join[Join[delta2, zero, 2], Join[prod2, inv2, 2]];
    table = table2[[Join[perm, perm + n]]];
    indices = Pick[Range[n], h, 1];
    table[[Join[indices, indices + n]]] = table[[Join[indices + n, indices]]];
    PauliStabilizer[<|"Matrix" -> Mod[table1 . table, 2], "Phase" -> RandomInteger[1, 2 n]|>]
]


PauliStabilizer["Random", n : _Integer ? Positive : 5] := RandomClifford[n]
