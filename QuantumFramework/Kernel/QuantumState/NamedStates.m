Package["Wolfram`QuantumFramework`"]

PackageScope["$QuantumStateNames"]



$QuantumStateNames = {
    "0", "Zero", "Up", "1", "One", "Down",
    "Plus", "Minus", "Left", "Right",
    "PsiPlus", "PsiMinus", "PhiPlus", "PhiMinus",
    "BasisState", "Register",
    "UniformSuperposition",
    "UniformMixture",
    "RandomPure", "RandomMixed",
    "GHZ", "Bell", "Dicke",
    "W",
    "Werner",
    "Graph",
    "BlochVector"
}

QuantumState[] := QuantumState["0"]

QuantumState[""] := QuantumState[{}, QuantumBasis[1]]


QuantumState[s_String /; StringMatchQ[s, DigitCharacter..], args___] := With[{
    basis = QuantumBasis[args, "Label" -> s]
},
    QuantumState["BasisState"[Clip[Interpreter[DelimitedSequence["Digit", ""]] @ s, {0, basis["Dimension"] - 1}]], basis]
]

QuantumState[s_String /; StringMatchQ[s, ("0" | "1" | "+" | "-" | "L" | "R") ..], args___] :=
    QuantumState[QuantumTensorProduct[QuantumState /@ Characters[s]], args]


QuantumState[("Zero" | "Up")[args___], opts___] := QuantumState["0"[args], opts]

QuantumState[("One" | "Down")[args___], opts___] := QuantumState["1"[args], opts]

QuantumState["0"[]] := QuantumState[Normalize @ {1, 0}, "Label" -> "0"]

QuantumState["1"[]] := QuantumState[Normalize @ {0, 1}, "Label" -> "1"]


QuantumState["Plus"[]] := QuantumState[Normalize @ {1, 1}, "Label" -> "+"]
QuantumState["+"] := QuantumState["Plus"]

QuantumState["Minus"[]] := QuantumState[Normalize @ {1, -1}, "Label" -> "-"]
QuantumState["-"] := QuantumState["Minus"]

QuantumState["Left"[]] := QuantumState[Normalize @ {1, -I}, "Label" -> "L"]
QuantumState["L"] := QuantumState["Left"]
QuantumState["-i"] := QuantumState["Left"]

QuantumState["Right"[]] := QuantumState[Normalize @ {1, I}, "Label" -> "R"]
QuantumState["R"] := QuantumState["Right"]
QuantumState["+i"] := QuantumState["Right"]


QuantumState["PhiPlus"[]] := QuantumState[Normalize @ {1, 0, 0, 1}, "Label" -> "\*SubscriptBox[\[CapitalPhi], \(+\)]"]

QuantumState["PhiMinus"[]] := QuantumState[Normalize @ {1, 0, 0, -1}, "Label" -> "\*SubscriptBox[\[CapitalPhi], \(-\)]"]

QuantumState["PsiPlus"[]] := QuantumState[Normalize @ {0, 1, 1, 0}, "Label" -> "\*SubscriptBox[\[CapitalPsi], \(+\)]"]

QuantumState["PsiMinus"[]] := QuantumState[Normalize @ {0, 1, -1, 0}, "Label" -> "\*SubscriptBox[\[CapitalPsi], \(-\)]"]

QuantumState[(name : "Plus" | "Minus" | "Left" | "Right" | "PsiPlus" | "PsiMinus" | "PhiPlus" | "PhiMinus")[n_Integer ? Positive], args___] :=
    QuantumTensorProduct @ Table[QuantumState[name, args], n]


QuantumState["BasisState"[basisElement_List : {1}], args___] := Enclose @ Block[{basis, dimension, elementPosition},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[Length[basisElement] / basis["Qudits"]]];
    dimension = basis["Dimension"];
    elementPosition = FromDigits[basisElement, First[basis["Dimensions"]]] + 1;
    ConfirmAssert[1 <= elementPosition <= dimension];
    QuantumState[SparseArray[{elementPosition} -> 1, dimension], basis]
]

QuantumState[("Register" | "RandomPure" | "GHZ" | "UniformSuperposition" | "W")[0, ___], args___] := QuantumState[1, 1, args]

QuantumState[("RandomMixed" | "UniformMixture")[0, ___], args___] := QuantumState[{{1}}, 1, args]

QuantumState["Register"[subsystemCount: _Integer ? Positive : 1, state : _Integer ? NonNegative : 0], args___] := Enclose @ Block[{basis, dimension},
    basis = ConfirmBy[QuantumBasis[args, "Label" -> state], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    dimension = basis["Dimension"];
    ConfirmAssert[0 <= state < dimension];
    QuantumState[SparseArray[{{state + 1} -> 1}, dimension], basis]
]

QuantumState["Register"[basisArg_, state : _Integer ? NonNegative : 0], args___] := Enclose @ With[{basis = ConfirmBy[QuantumBasis[basisArg, args], QuantumBasisQ]},
    QuantumState["Register"[basis["Qudits"], state], basis]
]


QuantumState["UniformSuperposition"[subsystemCount : _Integer ? Positive : 1], args___] := Enclose @ Block[{basis},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    QuantumState[ConstantArray[1, basis["Dimension"]], basis]["Normalized"]
]


QuantumState["UniformMixture"[subsystemCount : _Integer ? Positive : 1], args___] := Enclose @ Block[{basis, dimension},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    dimension = basis["Dimension"];
    QuantumState[identityMatrix[dimension] / dimension, basis]
]


QuantumState["RandomPure"[subsystemCount : _Integer ? Positive : 1], args___] := Enclose @ Block[{basis},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    QuantumState[QuantumOperator["RandomUnitary"[basis["Dimensions"]], Range[basis["OutputQudits"]]]["Matrix"][[1]], basis]
]


QuantumState["RandomMixed"[subsystemCount : _Integer ? Positive : 1], args___] := Enclose @ Block[{basis, dimension, m},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    dimension = basis["Dimension"];
    m = RandomComplex[{-1 - I, 1 + I}, Table[dimension, 2]];
    QuantumState[m . ConjugateTranspose[m], basis]["Normalized"]
]


QuantumState["GHZ"[subsystemCount : _Integer ? Positive : 3], args___] := Enclose @ Block[{basis, dimension},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    dimension = basis["Dimension"];
    QuantumState[SparseArray[{{1} -> 1 / Sqrt[2], {dimension} -> 1 / Sqrt[2]}, dimension], basis]
]

QuantumState["Bell"[args___], opts___] := QuantumState["GHZ"[2, args], opts]


QuantumState["W"[subsystemCount : _Integer ? Positive : 3], args___] := Enclose @ Block[{basis, dimension},
    basis = ConfirmBy[QuantumBasis[args], QuantumBasisQ];
    basis = QuantumBasis[basis, Ceiling[subsystemCount / basis["Qudits"]]];
    ConfirmAssert[Equal @@ basis["Dimensions"]];
    dimension = First @ basis["Dimensions"];
    QuantumState[SparseArray[{element_} /; IntegerQ[Log[dimension, element - 1]] -> 1 / Sqrt[basis["Qudits"]], {basis["Dimension"]}], basis]
]



wernerState[p_, qb_QuditBasis] /; qb["Qudits"] == 2 :=
    Module[{dim = qb["Dimension"], d = First[qb["Dimensions"]], fAB, sym, as},
        fAB = QuantumOperator["Permutation"[qb["Dimensions"], Cycles[{{1, 2}}]]]["Matrix"];
        sym = IdentityMatrix[dim] + fAB;
        as = IdentityMatrix[dim] - fAB;
        QuantumState[p sym 1 / (dim + d) + 1 / (dim - d) (1 - p) as // Simplify, qb]
    ]

QuantumState["Werner"[p_ : .5, param_ : 2], args___] /; ! QuditBasisQ[param] :=
    Enclose @ QuantumState["Werner"[p, ConfirmBy[QuditBasis[param], QuditBasisQ]], args]

QuantumState["Werner"[p_, qb_ ? QuditBasisQ], args___] := With[{
    basis = QuditBasis[qb, 2]
},
    QuantumState[wernerState[p, basis], args]
]


QuantumState["Graph"[graph : _ ? GraphQ : RandomGraph[{4, 4}]], args___] := Module[{
    indexGraph, quditCount, entanglements
},
    indexGraph = IndexGraph[graph];
    quditCount = VertexCount[indexGraph];
    entanglements = OperatorApplied[Take][2] @* List @@@ EdgeList[indexGraph];

    QuantumState[
        Fold[
            QuantumOperator["CZ", #2] @ #1 &,
            QuantumState["UniformSuperposition"[quditCount]],
            entanglements
        ],
        args
    ]
]

QuantumState["BlochVector"[r_ : {0, 0, 1}], args___] /; VectorQ[r] := Enclose @ With[{d = Ceiling[Sqrt[Length[r] + 1]]},
    QuantumState[IdentityMatrix[d, SparseArray] / d + Sqrt[(d - 1) / 2 / d] PadRight[r, d ^ 2 - 1] . GellMannMatrices[d], ConfirmBy[QuantumBasis[d, args], QuantumBasisQ]]
]


QuantumState["Dicke"[n : _Integer ? Positive : 3, k : _Integer : 1], args___] /; n >= k := QuantumState["Dicke"[{n - k, k}], args]

QuantumState["Dicke"[k_List], args___] /; VectorQ[k, IntegerQ[#] && NonNegative[#] & ] := Block[{
    n = Total[k], dim = Length[k], s
},
    s = Table[QuantumState["Register"[{dim}, i]], {i, 0, dim - 1}];
    QuantumState[Total[QuantumTensorProduct /@ Permutations @ Catenate[ConstantArray[#1, #2] & @@@ Transpose[{s, k}]]]["Normalized"], args]
]


QuantumState[SuperDagger[arg_], opts___] := QuantumState[arg, opts]["Dagger"]


QuantumState[name_String, opts___] /; MemberQ[$QuantumStateNames, name] := QuantumState[name[], opts]


QuantumState[name_String[args___], ___] /; ! MemberQ[$QuantumStateNames, name] := (
    Message[QuantumState::invalidName, Defer[name[args]]];
    Failure["InvalidName", <|"MessageTemplate" :> QuantumState::invalidName, "MessageParameters" :> {Defer[name[args]]}|>]
)
