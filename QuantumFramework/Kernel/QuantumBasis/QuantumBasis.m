Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Arrays`"]

PackageExport["QuantumBasis"]

PackageScope["quantumBasisQ"]
PackageScope["QuantumBasisQ"]
PackageScope["$QuantumBasisDataKeys"]
PackageScope["$QuantumBasisStructureProperties"]
PackageScope["$QuantumBasisPictures"]
PackageScope["MergeParameterSpecs"]
PackageScope["MergeSymbolDomains"]
PackageScope["symbolDomainsFacts"]
PackageScope["realDomainSymbols"]
PackageScope["roleSymbolDomains"]
PackageScope["roleDomainOption"]
PackageScope["allDomainsRecord"]
PackageScope["exactRealQ"]
PackageScope["onRealSymbols"]
PackageScope["conjugateOnDomains"]
PackageScope["mapEntries"]
PackageScope["zeroOnDomains"]
PackageScope["basisFactsQ"]
PackageScope["nonnegativeOnDomains"]
PackageScope["positiveSemidefiniteOnDomains"]
PackageScope["positiveSemidefiniteVerdict"]
PackageScope["zeroVerdict"]
PackageScope["yesNo"]
PackageScope["$verdictMessages"]
PackageScope["nonzeroEntries"]



$QuantumBasisPictures = {
    "Schrodinger", "Heisenberg", "Interaction", "PhaseSpace"
};

(* a picture string names a picture, never a basis, so it must not be swallowed by an
   output/input slot that would hand it to QuditBasis; same exclusion QuantumState makes *)
basisNameQ[name_] := nameQ[name] && ! MemberQ[$QuantumBasisPictures, name]


(* Messages *)

QuantumBasis::wrongData = "has wrong data: ``";

QuantumBasis::picture = "should have one of the following pictures " <> StringRiffle[$QuantumBasisPictures, ", "] <> ". But got: ``";

QuantumBasis::zeroDimension = "can't have zero dimension";

QuantumBasis::inconsistentNames = "element names should have the same length";

QuantumBasis::inconsistentInputs = "number of input qudits `` should be positive integer and less than or equal to number of total qudits ``";

QuantumBasis::inconsistentElements = "elements should have the same dimensions";

QuantumBasis::dependentElements = "elements should be linearly independent";

QuantumBasis::invalidSpec = "`` is not a valid basis specification";

QuantumBasis::domain = "`` is not a symbol domain; a domain is one of Reals, NonNegativeReals, PositiveReals, Integers, NonNegativeIntegers, PositiveIntegers, Complexes, or a real Interval[{a, b}]"

QuantumBasis::domainspec = "`` is not a symbol domain specification; give rules from a symbol or an indexed symbol to a domain or a list of domains; an empty list takes away every fact about the symbol"

QuantumBasis::domainempty = "the domains `2` given for `1` have no point in common; no fact about `1` is kept"

QuantumBasis::complexbounded = "`1` is declared both complex and in `2`, which bounds it as a real number; no fact about `1` is kept"

QuantumBasis::phaseSpacePicture = "a phase space basis is already a phase space representation and cannot be renamed to the picture ``";



(* "SymbolDomains": what is known about each symbol, by its role or from the user, as an association
   from a symbol, plain (θ) or indexed (θ[1], Subscript[θ, 1]), to the list of domains it lies in, such as {Reals}
   for the angle of a named gate or {Interval[{0, 1}]} for a probability; {Complexes} marks a symbol
   the user wants taken as complex, whatever its role. With the declared parameters, which are real
   unless marked complex, it makes "ParameterAssumptions". Each fact is about one symbol: a condition that
   ties two symbols together, p0 + p1 <= 1 or a Bloch vector of length at most 1, is not recorded, so a
   role records only what each of its symbols must satisfy on its own. *)
$QuantumBasisDataKeys = {"Input", "Output", "Picture", "Label", "ParameterSpec", "SymbolDomains"}

(* the basis's data and the shape it gives, which states and operators read off their basis *)
$QuantumBasisStructureProperties = Join[$QuantumBasisDataKeys, {
    "Qudits", "InputQudits", "OutputQudits",
    "Dimension", "InputDimension", "OutputDimension",
    "Dimensions", "InputDimensions", "OutputDimensions",
    "MatrixNameDimensions"
}]

$QuantumBasisDefaults = {
    "Input" -> QuditBasis[],
    "Output" -> QuditBasis[],
    "Picture" -> "Schrodinger",
    "Label" -> None,
    "ParameterSpec" -> {},
    "SymbolDomains" -> <||>
}


ParameterSpecQ[spec_] := MatchQ[spec, subSpecs : {{_, _ ? NumericQ, _ ? NumericQ}...} /; DuplicateFreeQ[subSpecs[[All, 1]]]]

quantumBasisQ[QuantumBasis[data : HoldPattern[Association[(_String -> _) ...]]]] := Enclose[
    ConfirmAssert[ContainsAll[Keys[data], $QuantumBasisDataKeys]];

    ConfirmAssert[QuditBasisQ[data["Input"]]];
    ConfirmAssert[QuditBasisQ[data["Output"]]];

    ConfirmAssert[MemberQ[$QuantumBasisPictures, data["Picture"]], Message[QuantumBasis::picture, data["Picture"]]];
    (*ConfirmAssert[Length[inputElements] + Length[outputElements] > 0, Message[QuantumBasis::zeroDimension]];*)

    ConfirmAssert[ParameterSpecQ[data["ParameterSpec"]]];
    ConfirmAssert[symbolDomainsQ[data["SymbolDomains"]]];
    True,

    False &
]

quantumBasisQ[___] := False


QuantumBasisQ[qb : QuantumBasis[_]] := System`Private`HoldValidQ[qb] || quantumBasisQ[Unevaluated @ qb]

QuantumBasisQ[___] := False


qb_QuantumBasis /; System`Private`HoldNotValidQ[qb] && quantumBasisQ[Unevaluated @ qb] := (
    System`Private`HoldSetValid[qb];
    System`Private`HoldSetNoEntry[qb]
)


(* mutation *)

QuantumBasis[qb_QuantumBasis] := qb

(* "PhaseSpace" is a tag, set below from $QuditPhaseSpaceBasisNames, and nothing ties it to
   the elements; QuantumState nonetheless routes on it. Renaming it therefore leaves the
   elements and the route disagreeing, and the state transforms an already-transformed
   basis. This closes the positional route on an already-tagged basis only. Other routes
   still put the two out of step (a "Picture" option, QuditBasis-built and list-built
   phase-space bases, QuantumTensorProduct, partial trace); closing those needs the tag
   derived from the elements rather than guarded here. The option form is deliberately
   left open: QuantumWeylTransform demotes through it (QuantumWignerTransform.m:81) once
   it has actually undone the transform *)
QuantumBasis[QuantumBasis[data_Association], picture : Alternatives @@ $QuantumBasisPictures] /;
    data["Picture"] === "PhaseSpace" && picture =!= "PhaseSpace" :=
    (Message[QuantumBasis::phaseSpacePicture, picture]; $Failed)

QuantumBasis[QuantumBasis[data_Association], picture : Alternatives @@ $QuantumBasisPictures] := QuantumBasis[<|data, "Picture" -> picture|>]

QuantumBasis[QuantumBasis[data_Association], rules : OptionsPattern[]] /; FilterRules[{rules}, "SymbolDomains"] === {} :=
    QuantumBasis @ KeyMap[Replace["Parameter" | "Parameters" -> "ParameterSpec"]][<|data, Reverse @ FilterRules[{rules}, Join[$QuantumBasisDataKeys, {"Parameter", "Parameters"}]]|>]

(* Every "SymbolDomains" given joins the basis's record. A record, as roles and operations pass it, meets the
   record already there. A declaration the user writes as rules does too, except that its plain Reals or
   Complexes for a symbol replaces whatever the symbol had, a role's range included. A single domain name,
   "SymbolDomains" -> Reals, is kept under All until the object's own symbols are known (allDomainsRecord). *)
QuantumBasis[QuantumBasis[data_Association], rules : OptionsPattern[]] := Enclose @ With[
    {given = Values @ FilterRules[Flatten[{rules}], "SymbolDomains"]},
    {bad = SelectFirst[given, ! symbolDomainsQ[#] && ! symbolDomainsSpecQ[#] && ! symbolDomainQ[#] &, None]},
    If[bad =!= None, domainSpecMessage[bad]; Confirm[$Failed]];
    QuantumBasis @ KeyMap[Replace["Parameter" | "Parameters" -> "ParameterSpec"]][<|
        data,
        Reverse @ FilterRules[{rules}, Join[DeleteCases[$QuantumBasisDataKeys, "SymbolDomains"], {"Parameter", "Parameters"}]],
        "SymbolDomains" -> If[given === {}, data["SymbolDomains"], givenDomainsRecord[data["SymbolDomains"], given]]
    |>]
]

givenDomainsRecord[record_, given_List] := With[{
    base = mergeSymbolDomains[Prepend[Select[given, symbolDomainsQ], record]],
    user = Select[given, ! symbolDomainsQ[#] && ! symbolDomainQ[#] &],
    all = Select[given, symbolDomainQ]
},
    withAllDomain[If[user === {}, base, withUserDomains[base, userDomainRecord[Flatten[Normal /@ user]]]], all]
]

withAllDomain[record_, {}] := record
withAllDomain[record_, all_List] := mergeSymbolDomains[{record, <|All -> DeleteDuplicates[all]|>}]

(* A declaration as the user writes it: for each symbol its domains meet each other; a complex mark wins over
   Reals, and beside a range or an integer domain, which bounds the symbol as a real number, contradicts it and
   leaves the symbol no fact. *)
userDomainRecord[spec_] := Association @ KeyValueMap[#1 -> userDomainList[#1, #2] &, recordDomainLists[spec]]

userDomainList[_, {}] := {}
userDomainList[s_, ds_] /; MemberQ[ds, Complexes] && ! SubsetQ[{Complexes, Reals}, ds] :=
    (Message[QuantumBasis::complexbounded, s, DeleteCases[DeleteDuplicates[ds], Complexes | Reals]]; {})
userDomainList[s_, ds_] := meetDomainLists[s, {DeleteDuplicates[ds]}]

(* the user's plain Reals or Complexes, or an empty list, replaces the symbol's facts; any other declaration
   meets them *)
withUserDomains[base_, user_] := Join[base, Association @ KeyValueMap[
    #1 -> Replace[#2, d : Except[{} | {Reals} | {Complexes}] :> If[KeyExistsQ[base, #1], meetDomainLists[#1, {base[#1], d}], d]] &,
    user
]]

(* "SymbolDomains" -> d, a single domain, made into a fact about every symbol of the object's amplitudes, met
   with what each symbol already has *)
allDomainsRecord[record_, values_] := With[{
    symbols = DeleteDuplicates @ Cases[values, $domainSymbolPattern, {0, Infinity}]
},
    mergeSymbolDomains[{KeyDrop[record, All], AssociationThread[symbols, ConstantArray[record[All], Length[symbols]]]}]
]

QuantumBasis[data_Association, args__] := Enclose @ Fold[ConfirmBy[QuantumBasis[##], QuantumBasisQ] &, ConfirmBy[QuantumBasis[data], QuantumBasisQ], Reverse @ {args}]



(* construction *)

QuantumBasis[elements : (a_ /; ArrayQ[a, d_ /; d > 1, NumericQ]) | _Association ? (Not @* KeyExistsQ["Output"]), args___] :=
    Enclose @ QuantumBasis[<|"Output" -> ConfirmBy[QuditBasis[elements], QuditBasisQ]|>, args]


QuantumBasis[names : {Except[_Integer] ...}, elements_ ? ArrayQ] /; Length[names] == Length[elements] && ArrayDepth[elements] > 1 := QuantumBasis[AssociationThread[names, Normal @ elements]]

QuantumBasis[output : _ ? basisNameQ | _QuditBasis | _List, input : _ ? basisNameQ | _QuditBasis | _List, args___] :=
    Enclose @ QuantumBasis["Output" -> ConfirmBy[QuditBasis[output], QuditBasisQ], "Input" -> ConfirmBy[QuditBasis[input], QuditBasisQ]["Dual"], args]

QuantumBasis[output_QuditBasis ? QuditBasisQ, args___] := QuantumBasis["Output" -> output, args]

QuantumBasis[arg : {(_QuditName | _Integer | _ ? propQ) ..}, args___] := Enclose @ QuantumBasis["Output" -> ConfirmBy[QuditBasis[arg], QuditBasisQ], args]

QuantumBasis[params_List, args : PatternSequence[] | PatternSequence[Except[_List], ___]] := Enclose @ QuantumTensorProduct[ConfirmBy[QuantumBasis[#, args], QuantumBasisQ] & /@ params]

QuantumBasis[{}, args___] := QuantumBasis[QuditBasis[{}], args]


(* defaults *)

QuantumBasis[data_Association ? (Keys /* Not @* ContainsAll[$QuantumBasisDataKeys]), args___] :=
    QuantumBasis[<|$QuantumBasisDefaults, data|>, args]


QuantumBasis[data_Association] /; AssociationQ[data["Input"]] := QuantumBasis[MapAt[QuditBasis, data, "Input"]]

QuantumBasis[data_Association] /; AssociationQ[data["Output"]] := QuantumBasis[MapAt[QuditBasis, data, "Output"]]


defaultParameterSpec[{p : Except[_List]}] := defaultParameterSpec[{p, 0, 1}]

defaultParameterSpec[{p_, i_ ? NumericQ}] := defaultParameterSpec[{p, i, i + 1}]

defaultParameterSpec[{p_, i_ ? NumericQ, f_ ? NumericQ}] := {{p, i, f}}

defaultParameterSpec[ps_List] := Catenate[defaultParameterSpec /@ ps]

(* defaultParameterSpec[spec : Except[{_List...}, _List]] := defaultParameterSpec[Take[spec, UpTo[3]]] *)

defaultParameterSpec[p : Except[_List]] := defaultParameterSpec[{p}]


QuantumBasis[data_Association] /; KeyExistsQ[data, "SymbolDomains"] && symbolDomainQ[data["SymbolDomains"]] :=
    QuantumBasis[<|data, "SymbolDomains" -> <|All -> {data["SymbolDomains"]}|>|>]

QuantumBasis[data_Association] /; KeyExistsQ[data, "SymbolDomains"] && ! symbolDomainsQ[data["SymbolDomains"]] && ! symbolDomainsSpecQ[data["SymbolDomains"]] :=
    (domainSpecMessage[data["SymbolDomains"]]; $Failed)

(* a valid specification merges into a valid record, so this rule fires at most once *)
QuantumBasis[data_Association] /; ! symbolDomainsQ[data["SymbolDomains"]] && symbolDomainsSpecQ[data["SymbolDomains"]] :=
    With[{record = userDomainRecord[data["SymbolDomains"]]},
        If[symbolDomainsQ[record], QuantumBasis[<|data, "SymbolDomains" -> record|>], $Failed]
    ]

(* a well-formed record with an unknown domain names that domain; anything else names the whole specification *)
domainSpecMessage[spec_] := With[{unknown = If[domainRulesQ[spec], SelectFirst[Flatten[Values @ Flatten[{Normal[spec]}]], ! symbolDomainQ[#] &, None], None]},
    If[unknown =!= None, Message[QuantumBasis::domain, unknown], Message[QuantumBasis::domainspec, spec]]
]

QuantumBasis[data_Association] /; ! MatchQ[data["ParameterSpec"], {{_, _, _}...}] || ! DuplicateFreeQ[data["ParameterSpec"]] :=
    QuantumBasis[<|data, "ParameterSpec" -> Reverse @ DeleteDuplicatesBy[Reverse @ defaultParameterSpec[data["ParameterSpec"]], First]|>]


(* multiplicity *)

QuantumBasis[qb_ ? QuantumBasisQ, 1, args___] := QuantumBasis[qb, args]

QuantumBasis[qb_ ? QuantumBasisQ, multiplicity_Integer, args___] := QuantumBasis[qb,
    "Input" -> QuditBasis[qb["Input"], multiplicity],
    "Output" -> QuditBasis[qb["Output"], multiplicity],
    args
]

QuantumBasis[arg_, multiplicity_Integer ? NonNegative, args___] := Enclose @ With[{
    bases = Table[ConfirmBy[QuantumBasis[arg, args], QuantumBasisQ], multiplicity]
},
    If[ multiplicity > 0,
        If[ Equal @@ bases,
            QuantumBasis[
                QuantumTensorProduct[bases],
                "Label" -> If[multiplicity > 1, First[bases]["Label"] ^ CircleTimes[multiplicity], First[bases]["Label"]]
            ],
            QuantumTensorProduct[bases]
        ],
        QuantumBasis[1, args]
    ]
]

QuantumBasis[output : _ ? basisNameQ | _Integer, input  : _ ? basisNameQ | _Integer, args___] :=
    Enclose @ QuantumBasis[<|"Output" -> ConfirmBy[QuditBasis[output], QuditBasisQ], "Input" -> ConfirmBy[QuditBasis[input], QuditBasisQ]["Dual"]|>, args]

QuantumBasis[output : _ ? basisNameQ | _Integer, args___] :=
    Enclose @ QuantumBasis[<|
        "Output" -> ConfirmBy[QuditBasis[output], QuditBasisQ],
        "Label" -> Replace[output, {name_String[___] | name_String :> name, 2 -> "I", i_Integer :> "I"[i]}],
        "Picture" -> If[MemberQ[$QuditPhaseSpaceBasisNames, nameString[output]], "PhaseSpace", "Schrodinger"]
        |>,
        args
    ]

(* a picture may sit anywhere among the option rules, including the rules the construction
   routes above emit: QuantumBasis["X", "Y", "Heisenberg"] arrives back here as
   QuantumBasis["Output" -> ..., "Input" -> ..., "Heisenberg"] *)
QuantumBasis[args : (_String ? (MatchQ[Alternatives @@ $QuantumBasisPictures]) | OptionsPattern[]) ...] :=
    QuantumBasis["Computational", args, "Label" -> None]

(* a picture leading an actual basis specification applies to the basis that specification
   builds, so QuantumBasis["Heisenberg", "PauliX"] is the PauliX basis in that picture. A
   tail starting with a non-negative Integer is left to the multiplicity rule above, where
   QuantumBasis["Heisenberg", 3] means 3 qudits rather than dimension 3. A tail starting
   with a rule is left to the rule above, which also claims a tail of further picture
   strings; that overlap is harmless because both paths apply the leading picture last, the
   option path by folding over Reverse @ {args} *)
QuantumBasis[
    picture : Alternatives @@ $QuantumBasisPictures,
    args : PatternSequence[Except[_Rule | _RuleDelayed | _Integer ? NonNegative], ___]
] :=
    Enclose @ With[{basis = QuantumBasis[args]},
        QuantumBasis[
            (* a spec that failed loudly has already said why; only name one that went quiet.
               A Failure among the arguments leaves QuantumBasis[args] inert rather than
               failed, so test the arguments too, not just the result *)
            ConfirmBy[basis, QuantumBasisQ,
                If[FreeQ[{basis, args}, _Failure], Message[QuantumBasis::invalidSpec, {args}]]
            ],
            picture
        ]
    ]

(* qb's picture travels through the guarded positional form, not inside qb["Options"];
   passing the whole option list would carry a picture straight into the option rule and
   around the phase-space guard above *)
QuantumBasis[qb_QuantumBasis, args__] := Enclose @ Fold[
    ConfirmBy[QuantumBasis[#1, #2], QuantumBasisQ] &,
    ConfirmBy[QuantumBasis[args], QuantumBasisQ],
    Append[FilterRules[qb["Options"], Except["Picture"]], qb["Picture"]]
]

QuantumBasis[qs_QuantumState] := qs["Basis"]

(* equality *)

QuantumBasis /: Equal[qb__QuantumBasis] := Equal @@ (#["Input"] & /@ {qb}) && Equal @@ (#["Output"] & /@ {qb})

QuantumBasis /: Unequal[qb__QuantumBasis] := ! Equal[qb]


(* N *)

N[qb_QuantumBasis, n_] /; ! AllTrue[
    Join[Values[qb["Output"]["Representations"]], Values[qb["Input"]["Representations"]]],
    InexactNumberQ[#] || ArrayQ[#, _, InexactNumberQ] &] := QuantumBasis[qb, "Output" -> N[qb["Output"], n], "Input" -> N[qb["Input"], n]]

SetAttributes[QuantumBasis, NHoldAll]


(* addition *)

QuantumBasis /: Plus[qb__QuantumBasis ? QuantumBasisQ] :=
    QuantumBasis[
        "Output" -> Plus @@ (#["Output"] & /@ {qb}),
        "Input" -> Plus @@ (#["Input"] & /@ {qb}),
        "Label" -> Plus @@ (#["Label"] & /@ {qb})
    ]

MergeParameterSpecs[objs___] := Enclose @ Block[{specs = Catenate[ConfirmBy[#["ParameterSpec"], ParameterSpecQ] & /@ {objs}]},
    KeyValueMap[Prepend[#2, #1] &, GroupBy[specs, First, {Min[#[[All, 2]]], Max[#[[All, 3]]]} & ]]
]

MergeSymbolDomains[objs___] := mergeSymbolDomains[If[Head[#] === QuantumBasis, #, #["Basis"]]["SymbolDomains"] & /@ {objs}]


(* symbol domains *)

(* A domain is a named set of numbers or a real interval. Each real one is read as its bounds, whether its
   lower bound is open, and whether it holds integers only, so that inclusion is a comparison of bounds;
   Complexes holds every other one. *)
$symbolDomainNames = {Reals, NonNegativeReals, PositiveReals, Integers, NonNegativeIntegers, PositiveIntegers, Complexes}

symbolDomainQ[d_] := MemberQ[$symbolDomainNames, d] || MatchQ[d, Interval[{a_ ? realNumberQ, b_ ? realNumberQ}] /; a <= b]

realNumberQ[x_] := NumericQ[x] && Im[x] == 0

(* exact, so Interval keeps it as written: rationals and constants such as Pi or 1/Sqrt[2] *)
exactRealQ[x_] := realNumberQ[x] && Precision[x] === Infinity

bounds[integer_, lower_, upper_, lowerOpen_] := <|"Integer" -> integer, "Lower" -> lower, "Upper" -> upper, "LowerOpen" -> lowerOpen|>

domainBounds[Reals] = bounds[False, -Infinity, Infinity, False]
domainBounds[NonNegativeReals] = bounds[False, 0, Infinity, False]
domainBounds[PositiveReals] = bounds[False, 0, Infinity, True]
domainBounds[Integers] = bounds[True, -Infinity, Infinity, False]
domainBounds[NonNegativeIntegers] = bounds[True, 0, Infinity, False]
domainBounds[PositiveIntegers] = bounds[True, 1, Infinity, False]
(* an interval holding one integer point is an integer domain *)
domainBounds[Interval[{a_, b_}]] := bounds[a == b && IntegerQ[a], a, b, False]

(* d1 is a subset of d2 *)
subdomainQ[_, Complexes] := True
subdomainQ[Complexes, _] := False
subdomainQ[d1_, d2_] := With[{b1 = domainBounds[d1], b2 = domainBounds[d2]},
    (b1["Integer"] || ! b2["Integer"]) && b1["Upper"] <= b2["Upper"] &&
        (b1["Lower"] > b2["Lower"] || b1["Lower"] == b2["Lower"] && (b1["LowerOpen"] || ! b2["LowerOpen"]))
]

(* The domains a symbol lies in, all of them holding at once. The intervals among them meet in one interval,
   and a domain that contains another one in the list says nothing more and is dropped (two different domains
   never contain each other). Complexes, the user's mark that the symbol is not to be taken as real, wins over
   any other. The result is {} when the domains have no point in common. *)
reduceDomains[{d_}] := {d}

reduceDomains[ds_List] := Which[
    MemberQ[ds, Complexes], {Complexes},
    emptyDomainsQ[ds], {},
    True, With[{u = DeleteDuplicates @ meetIntervals[ds]},
        Select[u, d |-> NoneTrue[DeleteCases[u, d], subdomainQ[#, d] &]]
    ]
]

meetIntervals[ds_List] := With[{is = Cases[ds, _Interval]},
    If[Length[is] < 2, ds, Insert[DeleteCases[ds, _Interval], IntervalIntersection @@ is, First @ FirstPosition[ds, _Interval]]]
]

(* the domains meet nowhere: the largest lower bound passes the smallest upper bound, or the integers among
   them leave no integer between the two *)
emptyDomainsQ[ds_List] := With[
    {bs = domainBounds /@ ds},
    {lo = Max[Lookup[bs, "Lower"]], hi = Min[Lookup[bs, "Upper"]], int = AnyTrue[Lookup[bs, "Integer"], TrueQ]},
    {open = AnyTrue[Select[bs, #["Lower"] == lo &], #["LowerOpen"] &]},
    lo > hi || lo == hi && open ||
        int && If[open && IntegerQ[lo], lo + 1, Ceiling[lo]] > Floor[hi]
]

(* A record maps each symbol to its reduced list of domains, or to {} for a symbol whose domains were found to
   have no common point: that symbol has no fact, and no later merge gives it one. *)
symbolDomainsQ[<||>] := True

symbolDomainsQ[a_] := AssociationQ[a] && AllTrue[Keys[a], domainKeyQ] &&
    AllTrue[Values[a], ListQ[#] && AllTrue[#, symbolDomainQ] && (Length[#] <= 1 || reduceDomains[#] === #) &]

(* a record as the user may give it: one rule, a list of rules or an association, from a symbol or an
   indexed symbol, each value one domain or a list of them; an empty list takes away every fact *)
domainRulesQ[a_] := MatchQ[a, _Association | {___Rule} | _Rule] && AllTrue[Keys @ Flatten[{Normal[a]}], domainKeyQ]

domainKeyQ[k_] := MatchQ[k, $domainSymbolPattern | All]

symbolDomainsSpecQ[a_] := domainRulesQ[a] && AllTrue[Values @ Flatten[{Normal[a]}], symbolDomainQ[#] || ListQ[#] && AllTrue[#, symbolDomainQ] &]

(* Records are merged symbol by symbol. A symbol whose domains have no point in common gets one message and
   no fact at all, so nothing is ever simplified under false facts. *)
mergeSymbolDomains[records_List] := With[{given = DeleteDuplicates @ DeleteCases[records, <||> | {}]},
    Which[
        given === {}, <||>,
        Length[given] == 1 && symbolDomainsQ[First[given]], First[given],
        True, mergeDomainRecords[given]
    ]
]

mergeDomainRecords[records_List] := Association @ KeyValueMap[
    Function[{s, dss}, s -> meetDomainLists[s, dss]],
    Merge[recordDomainLists /@ records, Identity]
]

meetDomainLists[s_, dss_] := If[MemberQ[dss, {}], {},
    With[{ds = Catenate[dss]}, {r = reduceDomains[ds]}, If[r === {}, Message[QuantumBasis::domainempty, s, ds]]; r]
]

(* one record as an association from each symbol to its list of domains; a symbol named twice gets both lists,
   and an empty list among them, given to take away every fact, stays empty *)
recordDomainLists[record_] := Merge[Flatten[{Normal[record]}],
    If[MemberQ[#, {}], {}, Catenate[Flatten @* List /@ #]] &
]

(* The record a role gives: a slot that lies in domain, such as the angle of a gate (Reals) or the
   probability of a channel (Interval[{0, 1}]), puts the symbol it holds in that domain. The symbol may
   be plain or indexed, theta, theta[1] or Subscript[theta, 1]. A slot holding a theta + b, with real
   numbers a != 0 and b, puts theta in the preimage of the domain: Reals for a real slot, and for a
   probability slot holding 1 - p, p in [0, 1] again; an interval preimage is taken only for exact a and b,
   since Interval widens inexact bounds. Any other expression, a + b say, says nothing about its symbols
   and is skipped. A symbol in two slots lies in both domains. *)
roleSymbolDomains[slots_List, domain_] := Association @ KeyValueMap[
    #1 -> meetDomainLists[#1, List /@ #2] &,
    Merge[slotDomain[#, domain] & /@ slots, Identity]
]

(* the "SymbolDomains" option for that record, or nothing when the slots hold no symbol, so a numeric gate
   or state is built exactly as without facts *)
roleDomainOption[slots_List, _] /; AllTrue[slots, NumericQ] := Sequence @@ {}

roleDomainOption[slots_List, domain_] := With[{r = roleSymbolDomains[slots, domain]},
    If[r === <||>, Sequence @@ {}, "SymbolDomains" -> r]
]

$domainSymbolPattern = _Symbol ? domainSymbolQ | (_Symbol ? domainSymbolQ)[__ ? AtomQ] | Subscript[_Symbol ? domainSymbolQ, __ ? AtomQ]

domainSymbolQ[s_Symbol] := ! NumericQ[s] &&
    (Context[s] =!= "System`" || StringStartsQ[ToString[s, InputForm, CharacterEncoding -> "PrintableASCII"], "\\[Formal"])

slotDomain[x_, domain : Reals | _Interval] := Replace[x, {
    v : $domainSymbolPattern :> v -> domain,
    (a_. (v : $domainSymbolPattern)) + b_. /; realNumberQ[a] && a != 0 && realNumberQ[b] :> preimage[v, domain, a, b],
    a_ (v : $domainSymbolPattern) /; realNumberQ[a] && a != 0 :> preimage[v, domain, a, 0],
    _ -> Nothing
}]

slotDomain[x_, domain_] := Replace[x, {v : $domainSymbolPattern :> v -> domain, _ -> Nothing}]

(* the values of theta for which a theta + b lies in the domain *)
preimage[v_, Reals, _, _] := v -> Reals
preimage[v_, Interval[{lo_, hi_}], a_ ? exactRealQ, b_ ? exactRealQ] := v -> Interval[{(lo - b) / a, (hi - b) / a}]
preimage[_, _Interval, _, _] := Nothing

(* the symbols, among the given ones and the keys of a record, taken as real; a symbol
   marked Complexes, or whose domains had no common point, is not *)
realDomainSymbols[record_, declared_List] := DeleteCases[All] @ DeleteDuplicates @ Join[
    Select[declared, ! MatchQ[Lookup[record, Key[#], None], {Complexes} | {}] &],
    Keys @ Select[record, ! MatchQ[#, {Complexes} | {}] &]
]

(* the facts a record states, for Refine or Simplify: every symbol that is only real in one Element, then the
   stronger facts *)
symbolDomainsFacts[record_, declared_List] := With[{
    reals = Select[realDomainSymbols[record, declared], MatchQ[Lookup[record, Key[#], {Reals}], {Reals}] &]
},
    And @@ Prepend[
        Catenate @ KeyValueMap[Function[{s, ds}, domainFact[s, #] & /@ ds], Select[KeyDrop[record, All], ! MatchQ[#, {Reals} | {Complexes}] &]],
        Element[reals, Reals]
    ]
]

domainFact[s_, Interval[{a_, b_}]] := a <= s <= b
domainFact[s_, d_] := Element[s, d]


(* computing on the facts *)

(* Each real symbol of a basis with the closed interval it lies in, the whole real line for a symbol known
   only to be real, and whether it is an integer. *)
symbolRanges[qb_] := With[{record = qb["SymbolDomains"]},
    AssociationMap[domainsRange[Lookup[record, Key[#], {Reals}]] &, qb["RealSymbols"]]
]

domainsRange[ds_List] := With[{bs = domainBounds /@ ds},
    {Interval[{Max[Lookup[bs, "Lower"]], Min[Lookup[bs, "Upper"]]}], AnyTrue[Lookup[bs, "Integer"], TrueQ]}
]

(* f applied entry by entry; a sparse array keeps its pattern *)
mapEntries[f_, a_SparseArray] := SparseArray[Thread[a["NonzeroPositions"] -> f /@ a["NonzeroValues"]], Dimensions[a], f[a["Background"]]]
mapEntries[f_, a_List] := Map[f, a, {ArrayDepth[a]}]
mapEntries[f_, a_] := f[a]

(* a sparse array is atomic to FreeQ, so its values are searched instead *)
entriesFreeQ[a_SparseArray, p_] := FreeQ[{a["NonzeroValues"], a["Background"]}, p]
entriesFreeQ[a_, p_] := FreeQ[a, p]

(* f applied to an array, a matrix or a scalar whose real symbols, those of qb, are first swapped for private
   symbols that are their own conjugate and real part and have imaginary part 0. Conjugate, Re and Im of any
   expression built from them with real operations (polynomials, Cos, Exp[I x]) then evaluate away while f
   runs. A Conjugate or Abs still left on an expression in these symbols alone is settled by its sign
   (signSettled), and the user's symbols are put back. The private symbols are removed on every exit, and
   none reaches the result. *)
(* whether a basis records any fact or parameter, read off its data without the property machinery, so that
   objects with neither pay nothing for the facts *)
basisFactsQ[QuantumBasis[data_Association]] := data["SymbolDomains"] =!= <||> || data["ParameterSpec"] =!= {}
basisFactsQ[_] := True

onRealSymbols[f_, a_, qb_] /; ! basisFactsQ[qb] := f[a]

onRealSymbols[f_, a_, qb_] := With[{ranges = symbolRanges[qb]},
    If[ ranges === <||> || entriesFreeQ[a, Alternatives @@ Keys[ranges]],
        f[a],
        onPrivateSymbols[f, a, Keys[ranges], Values[ranges][[All, 1]]]
    ]
]

onPrivateSymbols[f_, a_, reals_List, intervals_List] := With[{copies = Table[Module[{r}, r], Length[reals]]},
    WithCleanup[
        Scan[realSymbol, copies];
        mapEntries[
            ReplaceAll[signSettled[conjugateDistributed[#], AssociationThread[copies, intervals]], Thread[copies -> reals]] &,
            f[mapEntries[ReplaceAll[Thread[reals -> copies]], a]]
        ],
        Remove @@ copies
    ]
]

(* Conjugate taken through sums, products and integer powers, which the kernel leaves whole when a term holds
   a symbol with no fact, so the private symbols inside can evaluate it away *)
conjugateDistributed[x_] /; FreeQ[x, Conjugate[_Plus | _Times | Power[_, _Integer]]] := x
conjugateDistributed[x_] := x //. {
    Conjugate[y : _Plus | _Times] :> Conjugate /@ y,
    Conjugate[Power[y_, n_Integer]] :> Conjugate[y] ^ n
}

realSymbol[r_Symbol] := (r /: Conjugate[r] = r; r /: Re[r] = r; r /: Im[r] = 0;)

(* Conjugate[g] is g, and Abs[g] is g or -g, for g in the private symbols alone whose sign is known. The sign
   is read first off interval arithmetic, each symbol replaced by its interval: a result that is an Interval
   means g is real on the whole domain, and its ends give the sign. What the intervals leave open goes to
   Refine, under the same facts and a time limit, and its answer is kept only if it is no larger. *)
signSettled[x_, _] /; FreeQ[x, _Conjugate | _Abs] := x

signSettled[x_, intervals_] := With[{rules = Normal[intervals], symbols = Keys[intervals]},
    Block[{$Assumptions = And @@ KeyValueMap[intervalFact, intervals]},
        x /. {
            Conjugate[g_] /; onlySymbolsQ[g, symbols] :> If[MatchQ[g /. rules, _Interval], g, refined[Conjugate[g]]],
            Abs[g_] /; onlySymbolsQ[g, symbols] :> absBySign[g, g /. rules]
        }
    ]
]

absBySign[g_, Interval[{lo_, _}, ___] ? (Min[#] >= 0 &)] := g
absBySign[g_, Interval[{_, _}, ___] ? (Max[#] <= 0 &)] := - g
absBySign[g_, _] := refined[Abs[g]]

intervalFact[s_, Interval[{lo_, hi_}]] := Element[s, Reals] && If[lo === - Infinity, True, lo <= s] && If[hi === Infinity, True, s <= hi]

$signTimeConstraint = 0.2

refined[e_] := With[{r = TimeConstrained[Refine[e], $signTimeConstraint, e]}, If[LeafCount[r] <= LeafCount[e], r, e]]

(* no symbol in x but the given ones; numeric functions and constants are not symbols here *)
onlySymbolsQ[x_, symbols_List] := FreeQ[x,
    s_Symbol /; ! MemberQ[symbols, s] && ! MemberQ[{Rational, Complex, Integer, Real}, s] && ! MemberQ[Attributes[s], NumericFunction | Constant],
    {0, Infinity},
    Heads -> True
]

(* The conjugate of an array, matrix or scalar whose entries may hold the real symbols of qb. *)
conjugateOnDomains[a_, qb_] := onRealSymbols[Replace[ArrayConjugate[#], ArrayConjugate[x_] :> Conjugate[x]] &, a, qb]

(* A yes/no question about every value the facts of qb allow is answered True, False, or undecided[reason, ...]
   with the reason it is not decided: "time" (a test ran out of time), "tests" (no test applies), "fails" (it
   fails for some values, and whether it holds for any is not known), "region" with where it holds, or
   "notnumber" with an entry that is not a number. answer turns the last into Indeterminate. *)
answer[v : True | False] := v
answer[_undecided] := Indeterminate

(* Whether every entry of a list is zero for all values the facts of qb allow. Proven zero: after TrigToExp,
   Together leaves numerator 0 (TrigToExp first, since Together alone takes Cos and Sin as unknowns), within a
   time limit. Not zero: the value at a fixed point inside the domains, real symbols inside their intervals (an
   integer for an integer symbol), any other symbol at a complex point, each given to 30 digits, is larger than
   10^-10. An expression that is not identically zero vanishes only on a lower-dimensional set, so the
   property fails for generic values and the answer is False. *)
zeroVerdict[es_List, qb_] := With[{tests = {#, provenZeroTest[#]} & /@ es},
    With[{rest = Select[tests, #[[2]] =!= True &]},
        Which[
            rest === {}, True,
            AnyTrue[rest[[All, 1]], nonzeroAtPointQ[#, symbolRanges[qb]] &], False,
            MemberQ[rest[[All, 2]], $Aborted], undecided["time"],
            True, undecided["tests"]
        ]
    ]
]

zeroOnDomains[es_List, qb_] := answer[zeroVerdict[es, qb]]

$zeroTestTimeConstraint = 2

provenZeroTest[e_] := If[e === 0, True, TimeConstrained[Numerator[Together[TrigToExp[e]]] === 0, $zeroTestTimeConstraint, $Aborted]]

provenZeroQ[e_] := TrueQ[provenZeroTest[e]]

nonzeroAtPointQ[e_, ranges_] := With[{
    v = e /. N[samplePoint[DeleteDuplicates @ Cases[e, $domainSymbolPattern, {0, Infinity}], ranges], 30]
},
    NumberQ[v] && Abs[v] > 10 ^ -10
]

samplePoint[symbols_List, ranges_] := BlockRandom[
    SeedRandom[1];
    (# -> Replace[Lookup[ranges, Key[#], None], {
        {i_Interval, int_} :> pointInside[i, int],
        None :> RandomInteger[{1, 99}] / 100 + I RandomInteger[{1, 99}] / 100
    }]) & /@ symbols
]

(* a point inside an interval, away from its ends; an integer inside it for an integer symbol *)
pointInside[Interval[{lo_, hi_}], integer_] := With[{x = Which[
    lo === - Infinity && hi === Infinity, RandomInteger[{-99, 99}] / 50 + 1 / 7,
    hi === Infinity, lo + 1 / 3 + RandomInteger[{1, 99}] / 50,
    lo === - Infinity, hi - 1 / 3 - RandomInteger[{1, 99}] / 50,
    True, lo + (hi - lo) RandomInteger[{10, 90}] / 100
]}, If[integer, Clip[Round[x], {Ceiling[lo], Floor[hi]}], x]]


(* Whether every expression of a list is nonnegative for all values the facts of qb allow. A proof comes from
   interval arithmetic, or else from Simplify under the facts, with each symbol that has no fact split into its
   real and imaginary parts, within a time limit. A negative value at an allowed point means it fails there; the
   values where it holds are then sought with Reduce over the real symbols, within a time limit: none gives
   False, a region gives undecided["region", r]. *)
nonnegativeVerdict[es_List, qb_] := With[{ranges = symbolRanges[qb]},
    With[{rest = Select[es, ! provenNonnegativeQ[#, ranges, qb] &]},
        Which[
            rest === {}, True,
            AnyTrue[rest, negativeAtPointQ[#, ranges] &], whereNonnegative[es, ranges, qb],
            True, undecided["tests"]
        ]
    ]
]

nonnegativeOnDomains[es_List, qb_] := answer[nonnegativeVerdict[es, qb]]

$regionTimeConstraint = 1

whereNonnegative[es_List, ranges_, qb_] := With[{symbols = DeleteDuplicates @ Cases[es, $domainSymbolPattern, {0, Infinity}]},
    If[ ! SubsetQ[Keys[ranges], symbols],
        undecided["fails"],
        With[{r = TimeConstrained[Reduce[And @@ Thread[es >= 0] && qb["ParameterAssumptions"], symbols, Reals], $regionTimeConstraint, $Failed]},
            Which[
                r === False, False,
                r === $Failed || ! FreeQ[r, Reduce], undecided["fails"],
                True, undecided["region", r]
            ]
        ]
    ]
]

$signProofTimeConstraint = 2

provenNonnegativeQ[e_, ranges_, qb_] := With[{x = e /. KeyValueMap[#1 -> First[#2] &, ranges]},
    TrueQ[NumericQ[x] && x >= 0] || MatchQ[x, _Interval] && TrueQ[Min[x] >= 0] ||
        TrueQ @ TimeConstrained[
            Simplify[ComplexExpand[e, Complement[DeleteDuplicates @ Cases[e, $domainSymbolPattern, {0, Infinity}], Keys[ranges]]] >= 0, qb["ParameterAssumptions"]],
            $signProofTimeConstraint,
            False
        ]
]

negativeAtPointQ[e_, ranges_] := With[{
    v = e /. N[samplePoint[DeleteDuplicates @ Cases[e, $domainSymbolPattern, {0, Infinity}], ranges], 30]
},
    NumberQ[v] && (Re[v] < - 10 ^ -10 || Abs[Im[v]] > 10 ^ -10)
]

(* Whether a square matrix is positive semidefinite for all values the facts of qb allow. It must be Hermitian
   first; then its principal minors must be nonnegative (Sylvester's criterion in its form for semidefinite
   matrices, which needs every principal minor, not only the leading ones). That is 2^d - 1 determinants, so it
   is used up to d = 4; a diagonal matrix needs only its diagonal, and a larger one its eigenvalues, when they
   come within a time limit. *)
positiveSemidefiniteVerdict[m_, qb_] := With[{
    hermitian = zeroVerdict[nonzeroEntries @ onRealSymbols[# - ConjugateTranspose[#] &, m, qb], qb]
},
    If[ hermitian =!= True,
        hermitian,
        Replace[semidefiniteConditions[m], {
            cs_List :> nonnegativeVerdict[cs, qb],
            _ -> undecided["time"]
        }]
    ]
]

positiveSemidefiniteOnDomains[m_, qb_] := answer[positiveSemidefiniteVerdict[m, qb]]


(* A yes/no property of a symbolic object, answered by check[m, qb] on its matrix m: True, False, or
   Indeterminate with one message saying why. An entry that is not a number leaves it undecided. A False with
   symbols that have no fact is checked again with them real; when that gives True, the message names them. The
   messages are off while an object is formatted. *)
$verdictMessages = True

yesNo[head_, prop_, check_, m_, qb_] := With[{entries = nonzeroEntries[m]},
    With[{bad = SelectFirst[entries, ! FreeQ[#, _NonCommutativeMultiply | _MatrixSymbol | _ArraySymbol | _QuantumState | _QuantumOperator | _QuantumChannel] &, None]},
        If[ bad =!= None,
            reportVerdict[head, prop, undecided["notnumber", bad]],
            With[{v = check[m, qb]},
                If[ v === False && $verdictMessages,
                    undeclaredHint[head, prop, check, m, qb, entries]
                ];
                reportVerdict[head, prop, v]
            ]
        ]
    ]
]

undeclaredHint[head_, prop_, check_, m_, qb_, entries_] := With[{
    free = Complement[DeleteDuplicates @ Cases[entries, $domainSymbolPattern, {0, Infinity}], Keys[qb["SymbolDomains"]], qb["Parameters"]]
},
    If[ free =!= {} && check[m, QuantumBasis[qb, "SymbolDomains" -> Thread[free -> Reals]]] === True,
        Message[MessageName[head, "undeclared"], prop, free]
    ]
]

reportVerdict[_, _, v : True | False] := v
reportVerdict[head_, prop_, undecided[reason_, args___]] := (
    If[ $verdictMessages,
        With[{tag = Replace[reason, {"tests" -> "undecided", "time" -> "undecidedtime", "fails" -> "undecidedfails", "region" -> "undecidedregion"}]},
            Message[MessageName[head, tag], prop, args]
        ]
    ];
    Indeterminate
)

Scan[
    Function[head,
        MessageName[head, "undecided"] = "`1` could not be decided by the tests for this object; the answer is Indeterminate.";
        MessageName[head, "undecidedtime"] = "`1` could not be decided within the time limit; the answer is Indeterminate.";
        MessageName[head, "undecidedfails"] = "`1` fails for some values the facts allow, and whether it holds for any could not be decided; the answer is Indeterminate.";
        MessageName[head, "undecidedregion"] = "`1` holds only where `2`, not for every value the facts allow; the answer is Indeterminate.";
        MessageName[head, "notnumber"] = "`1` cannot be decided: the entry `2` is not a number; the answer is Indeterminate.";
        MessageName[head, "undeclared"] = "`1` is False because `2` may be complex; it is True when they are real, which \"SymbolDomains\" -> Reals declares."
    ],
    {QuantumState, QuantumOperator, QuantumChannel}
]

nonzeroEntries[m_SparseArray] := m["NonzeroValues"]
nonzeroEntries[m_] := Flatten[Normal[m]]

semidefiniteConditions[m_] := With[{d = Length[m]}, Which[
    AllTrue[SparseArray[m]["NonzeroPositions"], Apply[Equal]], Normal @ Diagonal[m],
    d <= 4, Det[Normal[m][[#, #]]] & /@ Subsets[Range[d], {1, d}],
    True, TimeConstrained[Eigenvalues[Normal[m]], $signProofTimeConstraint, $Failed]
]]

(* simplify *)

Scan[
    (Symbol[#][qb_QuantumBasis, args___] ^:= qb[#, args]) &,
    {"Simplify", "FullSimplify", "Chop", "ComplexExpand"}
]


(* parameterization *)

(qb_QuantumBasis ? QuantumBasisQ)[ps : PatternSequence[p : Except[_Association], ___]] /; ! MemberQ[QuantumBasis["Properties"], p] && Length[{ps}] <= qb["ParameterArity"] :=
    qb[AssociationThread[Take[qb["Parameters"], UpTo[Length[{ps}]]], {ps}]]

substituteRepresentations[qb_, rules_] :=
    KeyMap[ReplaceAll[rules]] @ Map[Map[ReplaceAll[rules], #, {ArrayDepth[#]}] &, qb["Representations"]]

(* A basis declares its parameters in ParameterSpec whether or not its elements
   depend on them: QuantumEvolve records the evolution parameter on the basis
   while the amplitudes alone carry it.  Substituting such a parameter changes
   nothing in the basis, and rebuilding it from its representations to discover
   that costs two orders of magnitude more than the substitution - the
   constructor is ~1540us of the ~1560us, flat in dimension, so it is fixed
   reconstruction overhead rather than work on the data.  When the substitution
   turns out to have been a no-op, the basis data is kept and only the spent
   parameter is dropped from the spec. *)

(qb : QuantumBasis[data_Association] ? QuantumBasisQ)[rules_ ? AssociationQ] /; ContainsOnly[Keys[rules], qb["Parameters"]] :=
    Block[{
        output = substituteRepresentations[qb["Output"], rules],
        input = substituteRepresentations[qb["Input"], rules],
        label = data["Label"] /. rules,
        spec = DeleteCases[data["ParameterSpec"], {Alternatives @@ Keys[rules], __}],
        domains = KeyDrop[data["SymbolDomains"], Keys[rules]]
    },
        If[ output === qb["Output"]["Representations"] &&
                input === qb["Input"]["Representations"] &&
                label === data["Label"],

            QuantumBasis[<|data, "ParameterSpec" -> spec, "SymbolDomains" -> domains|>],

            QuantumBasis[
                "Output" -> output,
                "Input" -> input,
                "Label" -> label,
                "Picture" -> data["Picture"],
                "ParameterSpec" -> spec,
                "SymbolDomains" -> domains
            ]
        ]
    ]

