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

(* The zero-qudit state is built from an empty amplitude list and so stores a lone
   zero, which is why its norm is 0 where "Register"[0] is also zero qudits but is
   built from {1} and is normalized. Options describe it as well as any other
   state, but it carries no levels for a basis to speak about, so a basis tail is
   a bad argument to a good name rather than a bad name. *)
QuantumState[""] := QuantumState[{}, QuantumBasis[1]]
QuantumState["", opts : OptionsPattern[]] := QuantumState[{}, QuantumBasis[1], opts]
QuantumState["", args__] := namedStateTailRejected[""]


QuantumState[s_String /; StringMatchQ[s, DigitCharacter..], args___] := With[{
    basis = QuantumBasis[args, "Label" -> s]
},
    QuantumState["BasisState"[Clip[Interpreter[DelimitedSequence["Digit", ""]] @ s, {0, basis["Dimension"] - 1}]], basis]
]

(* One qudit per character, so this rule only speaks for genuine sequences. A
   single character is the business of the alias rule that names it below, which
   would otherwise be unreachable at arity 1 and, with a tail, would hand the
   state this route's basis label rather than the state's own. *)
QuantumState[s_String /; StringLength[s] > 1 && StringMatchQ[s, ("0" | "1" | "+" | "-" | "L" | "R") ..], args___] :=
    QuantumState[QuantumTensorProduct[QuantumState /@ Characters[s]], args, "Label" -> s]


QuantumState[("Zero" | "Up")[args___], opts___] := QuantumState["0"[args], opts]

QuantumState[("One" | "Down")[args___], opts___] := QuantumState["1"[args], opts]

(* A fixed-vector name takes a basis specification followed by options. The
   name's amplitudes become the coefficients in that basis, zero-padded when it
   holds more levels; the tail says what the coefficients are coefficients of
   and does not rotate them into a new frame, which is the reading the digit
   string rule above already gives QuantumState["0", basis]. A trailing integer
   is therefore a qudit dimension and not a qudit count, the count staying
   inside the name as "Plus"[n]. These names denote kets of at least two levels,
   so a tail reaching here that is not a basis, or is one carrying an input leg
   or fewer than two levels, is rejected the way an unmatched call shape is:
   otherwise it would be taken quietly, turning the ket into an operator shaped
   object or dropping amplitudes to fit and leaving the state unnormalized. The
   bare strings "0" and "1" carry their tail to the digit string rule above
   instead of here, and keep that rule's handling of a malformed one. A
   dimension that does not resolve to a number is turned down by the same arm,
   since TrueQ sends it to the reject branch. The trailing "Label" is a default
   that an explicit one in the tail overrides. *)

namedStateTailRejected[name_] := (
    Message[QuantumState::invalidArgs, Defer[name[]]];
    Failure["InvalidArguments", <|"MessageTemplate" :> QuantumState::invalidArgs, "MessageParameters" :> {Defer[name[]]}|>]
)

namedStateTail[name_, vec_, label_, opts___] := With[{
    basis = QuantumBasis[opts]
},
    If[ ! QuantumBasisQ[basis] || basis["InputQudits"] =!= 0 || ! TrueQ[basis["Dimension"] >= 2],
        namedStateTailRejected[name],

        With[{state = QuantumState[vec, opts, "Label" -> label]},
            If[QuantumStateQ[state], state, namedStateTailRejected[name]]
        ]
    ]
]

QuantumState["0"[]] := QuantumState[Normalize @ {1, 0}, "Label" -> "0"]
QuantumState["0"[], opts__] := namedStateTail["0", Normalize @ {1, 0}, "0", opts]

QuantumState["1"[]] := QuantumState[Normalize @ {0, 1}, "Label" -> "1"]
QuantumState["1"[], opts__] := namedStateTail["1", Normalize @ {0, 1}, "1", opts]


(* Every alias forwards its tail to the name it spells, so the two spellings of
   one state are the same object at any arity. This works only because the
   letter-sequence rule above declines single characters: that rule decomposes a
   string into one qudit per character and calls back with each character alone,
   so an alias that did not outrank it would re-enter it and hit $RecursionLimit. *)

QuantumState["Plus"[]] := QuantumState[Normalize @ {1, 1}, "Label" -> "+"]
QuantumState["Plus"[], opts__] := namedStateTail["Plus", Normalize @ {1, 1}, "+", opts]
QuantumState["+", args___] := QuantumState["Plus", args]

QuantumState["Minus"[]] := QuantumState[Normalize @ {1, -1}, "Label" -> "-"]
QuantumState["Minus"[], opts__] := namedStateTail["Minus", Normalize @ {1, -1}, "-", opts]
QuantumState["-", args___] := QuantumState["Minus", args]

QuantumState["Left"[]] := QuantumState[Normalize @ {1, -I}, "Label" -> "L"]
QuantumState["Left"[], opts__] := namedStateTail["Left", Normalize @ {1, -I}, "L", opts]
QuantumState["L", args___] := QuantumState["Left", args]
QuantumState["-i", args___] := QuantumState["Left", args]

QuantumState["Right"[]] := QuantumState[Normalize @ {1, I}, "Label" -> "R"]
QuantumState["Right"[], opts__] := namedStateTail["Right", Normalize @ {1, I}, "R", opts]
QuantumState["R", args___] := QuantumState["Right", args]
QuantumState["+i", args___] := QuantumState["Right", args]


QuantumState["PhiPlus"[]] := QuantumState[Normalize @ {1, 0, 0, 1}, "Label" -> "\*SubscriptBox[\[CapitalPhi], \(+\)]"]
QuantumState["PhiPlus"[], opts__] := namedStateTail["PhiPlus", Normalize @ {1, 0, 0, 1}, "\*SubscriptBox[\[CapitalPhi], \(+\)]", opts]

QuantumState["PhiMinus"[]] := QuantumState[Normalize @ {1, 0, 0, -1}, "Label" -> "\*SubscriptBox[\[CapitalPhi], \(-\)]"]
QuantumState["PhiMinus"[], opts__] := namedStateTail["PhiMinus", Normalize @ {1, 0, 0, -1}, "\*SubscriptBox[\[CapitalPhi], \(-\)]", opts]

QuantumState["PsiPlus"[]] := QuantumState[Normalize @ {0, 1, 1, 0}, "Label" -> "\*SubscriptBox[\[CapitalPsi], \(+\)]"]
QuantumState["PsiPlus"[], opts__] := namedStateTail["PsiPlus", Normalize @ {0, 1, 1, 0}, "\*SubscriptBox[\[CapitalPsi], \(+\)]", opts]

QuantumState["PsiMinus"[]] := QuantumState[Normalize @ {0, 1, -1, 0}, "Label" -> "\*SubscriptBox[\[CapitalPsi], \(-\)]"]
QuantumState["PsiMinus"[], opts__] := namedStateTail["PsiMinus", Normalize @ {0, 1, -1, 0}, "\*SubscriptBox[\[CapitalPsi], \(-\)]", opts]

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

(* The mixing parameter defaults exactly, so the default Werner state stays a
   symbolic density matrix rather than a machine-precision one. *)
QuantumState["Werner"[p_ : 1/2, param_ : 2], args___] /; ! QuditBasisQ[param] :=
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


(* Check the inner state before daggering it: a name this route cannot resolve
   would otherwise return the Missing from a property query on a Failure. Passing
   the failure straight through keeps whatever tag and message the name itself
   produced. An Enclose/Confirm pair cannot do this job here, because Confirm
   rethrows an already-ConfirmationFailed inner failure unwrapped, and the
   handler then sees the asserted condition rather than a failure object. *)
QuantumState[SuperDagger[arg_], opts___] := With[{state = QuantumState[arg, opts]},
    If[FailureQ[state], state, state["Dagger"]]
]


QuantumState[name_String, opts___] /; MemberQ[$QuantumStateNames, name] := QuantumState[name[], opts]


QuantumState[name_String[args___], ___] /; ! MemberQ[$QuantumStateNames, name] := (
    Message[QuantumState::invalidName, Defer[name[args]]];
    Failure["InvalidName", <|"MessageTemplate" :> QuantumState::invalidName, "MessageParameters" :> {Defer[name[args]]}|>]
)

(* Last-resort arm for a bare string that names nothing. Without it an
   unrecognized name reached the scalar-amplitude rule in QuantumState.m, which
   wraps its argument as {name} and pads that to the default qubit basis, so
   QuantumState["NoSuchState"]["Norm"] came back Abs["NoSuchState"].

   THIS ARM MUST STAY AFTER THE STRING RULES ABOVE IT. Only a bare arity-1
   left-hand side outranks this arm's conditioned one whatever the load order,
   which covers the "" rule here and the QuantumState["Properties"] query defined
   in Properties.m. Nothing ranks the digit-string rule, the letter-sequence rule
   or any of the six aliases against it, all of which carry an argument tail, so
   they are reachable only by sitting earlier in this file: move the arm above
   them and QuantumState["0101"], QuantumState["+0L"], QuantumState["+"] and
   QuantumState["+i"] all start reporting an invalid name. *)
QuantumState[name_String, ___] /; ! MemberQ[$QuantumStateNames, name] := (
    Message[QuantumState::invalidName, name];
    Failure["InvalidName", <|"MessageTemplate" :> QuantumState::invalidName, "MessageParameters" :> {name}|>]
)

QuantumState[name_String[args___], ___] /; MemberQ[$QuantumStateNames, name] := (
    Message[QuantumState::invalidArgs, Defer[name[args]]];
    Failure["InvalidArguments", <|"MessageTemplate" :> QuantumState::invalidArgs, "MessageParameters" :> {Defer[name[args]]}|>]
)
