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
PackageScope["refineOnFacts"]
PackageScope["conjugateTransposeOnDomains"]
PackageScope["conjugateOnDomains"]
PackageScope["mapEntries"]
PackageScope["zeroOnDomains"]
PackageScope["basisFactsQ"]
PackageScope["nonnegativeOnDomains"]
PackageScope["positiveSemidefiniteOnDomains"]
PackageScope["positiveSemidefiniteVerdict"]
PackageScope["zeroVerdict"]
PackageScope["nonzeroOnDomains"]
PackageScope["signOnDomains"]
PackageScope["keptStage"]
PackageScope["simplifyOnFacts"]
PackageScope["yesNo"]
PackageScope["$verdictMessages"]
PackageScope["verdictScope"]
PackageScope["symbolDomainsOf"]
PackageScope["withRecord"]
PackageScope["exponentialZeroTest"]
PackageScope["assumedSignRecord"]
PackageScope["reportVerdict"]
PackageScope["undecided"]
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

QuantumBasis::domain = "`` is not a symbol domain; a domain is one of Reals, NonNegativeReals, PositiveReals, Integers, NonNegativeIntegers, PositiveIntegers, Complexes, or a real Interval[{a, b}], whose ends may be -Infinity and Infinity for a half-line"

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
$QuantumBasisDataKeys = {"Input", "Output", "Picture", "Label", "ParameterSpec"}

(* The record is kept in the data only when it holds a fact: a basis with no facts has the same data as one
   built before the record existed, so numeric work builds and compares bases as it did. *)
withRecord[data_Association, <||>] := KeyDrop[data, "SymbolDomains"]
withRecord[data_Association, record_] := <|data, "SymbolDomains" -> record|>

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
    "ParameterSpec" -> {}
}


ParameterSpecQ[spec_] := MatchQ[spec, subSpecs : {{_, _ ? NumericQ, _ ? NumericQ}...} /; DuplicateFreeQ[subSpecs[[All, 1]]]]

quantumBasisQ[QuantumBasis[data : HoldPattern[Association[(_String -> _) ...]]]] := Enclose[
    ConfirmAssert[ContainsAll[Keys[data], $QuantumBasisDataKeys]];

    ConfirmAssert[QuditBasisQ[data["Input"]]];
    ConfirmAssert[QuditBasisQ[data["Output"]]];

    ConfirmAssert[MemberQ[$QuantumBasisPictures, data["Picture"]], Message[QuantumBasis::picture, data["Picture"]]];
    (*ConfirmAssert[Length[inputElements] + Length[outputElements] > 0, Message[QuantumBasis::zeroDimension]];*)

    ConfirmAssert[ParameterSpecQ[data["ParameterSpec"]]];
    If[KeyExistsQ[data, "SymbolDomains"], ConfirmAssert[symbolDomainsQ[data["SymbolDomains"]]]];
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

QuantumBasis[QuantumBasis[data_Association], rules : OptionsPattern[]] /; FreeQ[Flatten[{rules}], ("SymbolDomains" -> Except[<||>]) | ("SymbolDomains" :> _), {1}] :=
    QuantumBasis @ KeyMap[Replace["Parameter" | "Parameters" -> "ParameterSpec"]][<|data, Reverse @ FilterRules[{rules}, Join[$QuantumBasisDataKeys, {"Parameter", "Parameters"}]]|>]

(* Every "SymbolDomains" given joins the basis's record. A record, as roles and operations pass it, meets the
   record already there. A declaration the user writes as rules does too, except that its plain Reals or
   Complexes for a symbol replaces whatever the symbol had, a role's range included. A single domain name,
   "SymbolDomains" -> Reals, is kept under All until the object's own symbols are known (allDomainsRecord). *)
(* a record QF built itself, already valid, joins the basis's record with no further checks *)
QuantumBasis[QuantumBasis[data_Association], rules : OptionsPattern[]] /; validRecordOptionsQ[{rules}] := With[{
    given = Values @ FilterRules[Flatten[{rules}], "SymbolDomains"]
},
    QuantumBasis @ withRecord[
        KeyMap[Replace["Parameter" | "Parameters" -> "ParameterSpec"]][<|data, Reverse @ FilterRules[{rules}, Join[$QuantumBasisDataKeys, {"Parameter", "Parameters"}]]|>],
        mergeSymbolDomains[Prepend[given, Lookup[data, "SymbolDomains", <||>]]]
    ]
]

validRecordOptionsQ[rules_] := AllTrue[Values @ FilterRules[Flatten[rules], "SymbolDomains"], AssociationQ[#] && ! KeyExistsQ[#, All] && symbolDomainsQ[#] &]

QuantumBasis[QuantumBasis[data_Association], rules : OptionsPattern[]] := Enclose @ With[
    {given = Values @ FilterRules[Flatten[{rules}], "SymbolDomains"]},
    {bad = SelectFirst[given, ! symbolDomainsQ[#] && ! symbolDomainsSpecQ[#] && ! symbolDomainQ[#] &, None]},
    If[bad =!= None, domainSpecMessage[bad]; Confirm[$Failed]];
    QuantumBasis @ KeyMap[Replace["Parameter" | "Parameters" -> "ParameterSpec"]][<|
        data,
        Reverse @ FilterRules[{rules}, Join[$QuantumBasisDataKeys, {"Parameter", "Parameters"}]],
        "SymbolDomains" -> If[given === {}, Lookup[data, "SymbolDomains", <||>], givenDomainsRecord[Lookup[data, "SymbolDomains", <||>], given]]
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


(* One rule for every basis: an empty record, the common case, passes at once; a record already valid passes
   after one check; anything else is turned into a record, or refused with a message. A valid specification
   merges into a valid record, so the rule fires at most once. *)
QuantumBasis[data_Association] /; ! readyRecordQ[Lookup[data, "SymbolDomains", None]] := With[{spec = data["SymbolDomains"]},
    Which[
        spec === <||>, QuantumBasis[KeyDrop[data, "SymbolDomains"]],
        symbolDomainQ[spec], QuantumBasis[<|data, "SymbolDomains" -> <|All -> {spec}|>|>],
        ! symbolDomainsSpecQ[spec], domainSpecMessage[spec]; $Failed,
        True, With[{record = userDomainRecord[spec]},
            If[symbolDomainsQ[record], QuantumBasis[<|data, "SymbolDomains" -> record|>], $Failed]
        ]
    ]
]

readyRecordQ[None] = True
readyRecordQ[<||>] = False
readyRecordQ[r_] := symbolDomainsQ[r]

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

MergeSymbolDomains[objs___] := mergeSymbolDomains[symbolDomainsOf /@ {objs}]

(* The record of a basis, or of the object holding it, read off the data without the property machinery: a
   merge runs on every composition and tensor product, where a property call costs more than the merge. *)
symbolDomainsOf[QuantumBasis[data_Association]] := Lookup[data, "SymbolDomains", <||>]
symbolDomainsOf[QuantumState[_, qb_QuantumBasis]] := symbolDomainsOf[qb]
symbolDomainsOf[(QuantumOperator | QuantumChannel | QuantumMeasurementOperator)[x : _QuantumState | _QuantumOperator, ___]] := symbolDomainsOf[x]
symbolDomainsOf[obj_] := obj["Basis"]["SymbolDomains"]


(* symbol domains *)

(* A domain is a named set of numbers or a real interval. Each real one is read as its bounds, whether its
   lower bound is open, and whether it holds integers only, so that inclusion is a comparison of bounds;
   Complexes holds every other one. *)
$symbolDomainNames = {Reals, NonNegativeReals, PositiveReals, Integers, NonNegativeIntegers, PositiveIntegers, Complexes}

symbolDomainQ[d_] := MemberQ[$symbolDomainNames, d] || MatchQ[d, Interval[{a_ ? lowerEndQ, b_ ? upperEndQ}] /; a <= b && ! (a === - Infinity && b === Infinity)]

(* an end of a range: a real number, or -Infinity below and Infinity above for a half-line *)
lowerEndQ[a_] := realNumberQ[a] || a === - Infinity
upperEndQ[b_] := realNumberQ[b] || b === Infinity

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

(* checked once per record: records are small, few, and checked on every basis built *)
symbolDomainsQ[a_Association] := Replace[Lookup[$validRecords, Key[a], Missing[]], _Missing :> rememberedRecordQ[a]]

$validRecords = <||>

rememberedRecordQ[a_] := With[{v = recordCheckQ[a]},
    If[Length[$validRecords] > 10000, $validRecords = <||>];
    $validRecords[a] = v
]

symbolDomainsQ[_] := False

recordCheckQ[a_] := AssociationQ[a] && AllTrue[Keys[a], domainKeyQ] &&
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

domainSymbolQ[s_Symbol] := Replace[Lookup[$domainSymbols, Key[s], Missing[]], _Missing :> ($domainSymbols[s] = ! NumericQ[s] &&
    (Context[s] =!= "System`" || StringStartsQ[ToString[s, InputForm, CharacterEncoding -> "PrintableASCII"], "\\[Formal"]))]

$domainSymbols = <||>

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

domainFact[s_, Interval[{- Infinity, b_}]] := s <= b
domainFact[s_, Interval[{a_, Infinity}]] := a <= s
domainFact[s_, Interval[{a_, b_}]] := a <= s <= b
domainFact[s_, d_] := Element[s, d]


(* signs the amplitudes need *)

(* A square root, another non-integer power or a log of a real expression of one of the given real symbols has
   two values or none unless the sign of its argument is known. When the facts of a record do not give that sign, the
   argument is taken nonnegative (positive under a log or a negative power) when the facts allow it, and of the
   opposite sign otherwise; the range of the symbol this leaves joins the record, so every later step reads it
   as a fact. The arguments are taken in turn, each under the ranges the ones before left. A sign on two
   symbols at once is not a range of either, and is not assumed. The result is {record, conditions}. *)
assumedSignRecord[record_, reals_List, values_] := Fold[
    assumeSign[reals],
    {record, {}},
    signArguments[values, reals]
]

signArguments[values_, reals_] := DeleteDuplicates @ Select[
    Join[
        Cases[values, Power[a_, r_ ? NumericQ] /; ! IntegerQ[r] :> {a, If[r > 0, False, True]}, {0, Infinity}],
        Cases[values, Log[a_] :> {a, True}, {0, Infinity}]
    ],
    With[{symbols = questionSymbols[First[#]]}, Length[symbols] === 1 && MemberQ[reals, First[symbols]] && FreeQ[First[#], _Complex]] &
]

assumeSign[reals_][{record_, conditions_}, {a_, strict_}] := With[{s = First @ questionSymbols[a]},
    With[{ranges = <|s -> domainsRange[Lookup[record, Key[s], {Reals}]]|>, facts = symbolDomainsFacts[KeyTake[record, {s}], {s}]},
        If[ intervalSign[a, ranges, False] =!= None,
            {record, conditions},
            With[{chosen = SelectFirst[
                If[strict, {a > 0, a < 0}, {a >= 0, a <= 0}],
                ! MatchQ[reduceOnFacts[#, {s}, facts], False | $Aborted] &,
                None
            ]},
                With[{range = If[chosen === None, None, conditionInterval[reduceOnFacts[chosen, {s}, facts], s]]},
                    If[ range === None,
                        {record, conditions},
                        {mergeSymbolDomains[{record, <|s -> {range}|>}], Append[conditions, chosen]}
                    ]
                ]
            ]
        ]
    ]
]

(* the closed interval a one-symbol condition from Reduce states, or None *)
conditionInterval[r_, s_] := Replace[r, {
    HoldPattern[Inequality[lo_, Less | LessEqual, s, Less | LessEqual, hi_]] :> Interval[{lo, hi}],
    HoldPattern[(LessEqual | Less)[lo_, s, hi_]] :> Interval[{lo, hi}],
    HoldPattern[(GreaterEqual | Greater)[s, lo_]] :> Interval[{lo, Infinity}],
    HoldPattern[(LessEqual | Less)[s, hi_]] :> Interval[{- Infinity, hi}],
    HoldPattern[Equal[s, c_]] :> Interval[{c, c}],
    _ -> None
}]

QuantumState::assumed = "the amplitudes need the signs of `1`; `2` is taken, where the facts allow it, and kept as a fact of the object."


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

(* whether a basis records any fact or parameter, read off its data without the property machinery, so that
   objects with neither pay nothing for the facts *)
basisFactsQ[QuantumBasis[data_Association]] := KeyExistsQ[data, "SymbolDomains"] || data["ParameterSpec"] =!= {}
basisFactsQ[_] := True

(* The entries of an array, a matrix or a scalar refined under the facts of qb: one Refine over all of them,
   inside a block holding the basis's fact expression, which is the same for every operation on the object,
   so the facts are prepared once and reused. Refine settles Conjugate, Abs, Re, Im and Sqrt of the real
   symbols and of functions of them, factor by factor, and on a time limit keeps every part it settled and
   leaves the rest as it was. It runs on amplitudes and entries as they are formed, never on a finished
   result, whose size can be many times theirs. An object with no facts, or entries free of its real symbols,
   is returned as it is. *)
$refineTimeConstraint = 5

refineOnFacts[a_, qb_] /; ! basisFactsQ[qb] := a

refineOnFacts[a_, qb_] := With[{reals = qb["RealSymbols"]},
    If[ reals === {} || entriesFreeQ[a, Alternatives @@ reals],
        a,
        Block[{$Assumptions = qb["ParameterAssumptions"]}, mapValues[keptValues[#, Refine[#, TimeConstraint -> $refineTimeConstraint]] &, a]]
    ]
]

keptValues[old_List, new_List] /; Length[old] === Length[new] := MapThread[keptStage, {old, new}]
keptValues[old_, new_] := keptStage[old, new]

(* f applied once to all the values of an array; a sparse array keeps its pattern *)
mapValues[f_, a_SparseArray] := SparseArray[Thread[a["NonzeroPositions"] -> f[a["NonzeroValues"]]], Dimensions[a], a["Background"]]
mapValues[f_, a_] := f[a]

(* The conjugate, and the conjugate transpose, of an array, matrix or scalar whose entries may hold the real
   symbols of qb. Conjugate is first taken through sums, products and integer powers, which the kernel leaves
   whole when a term holds a symbol with no fact, so that Refine sees each factor on its own. *)
conjugateOnDomains[a_, qb_] := With[{direct = structuralConjugate[a, qb]},
    If[direct =!= $Failed, direct, refineOnFacts[conjugatedEntries[Replace[ArrayConjugate[a], ArrayConjugate[x_] :> Conjugate[x]], qb], qb]]
]

(* The conjugate read off the structure when every entry is built from real symbols of qb, numbers, sums,
   products, integer powers, exponentials, and trig or hyperbolic functions of such expressions: each real
   symbol is its own conjugate, I becomes -I, and nothing needs Refine. Anything else gives $Failed. *)
structuralConjugate[a_, qb_] /; ! basisFactsQ[qb] := $Failed
structuralConjugate[a_SparseArray, qb_] := With[{vs = structuralConjugate[a["NonzeroValues"], qb]},
    If[vs === $Failed, $Failed, SparseArray[Thread[a["NonzeroPositions"] -> vs], Dimensions[a], Conjugate[a["Background"]]]]
]
structuralConjugate[a_, qb_] := With[{reals = qb["RealSymbols"]},
    If[ reals =!= {} && conjugateByStructureQ[Alternatives @@ reals][a], a /. Complex[x_, y_] :> Complex[x, - y], $Failed]
]

conjugateByStructureQ[reals_][x_] := Which[
    NumericQ[x], FreeQ[x, _Complex] || MatchQ[x, Complex[_Rational | _Integer, _Rational | _Integer]],
    MatchQ[x, reals], True,
    MatchQ[x, _List | _Plus | _Times], AllTrue[List @@ x, conjugateByStructureQ[reals]],
    MatchQ[x, Power[_, _Integer]], conjugateByStructureQ[reals][First[x]],
    MatchQ[x, Power[E, _]], conjugateByStructureQ[reals][Last[x]],
    MatchQ[x, (Sin | Cos | Tan | Cot | Sec | Csc | Sinh | Cosh | Tanh)[_]], conjugateByStructureQ[reals][First[x]] && FreeQ[First[x], _Complex],
    True, False
]

conjugateTransposeOnDomains[m_, qb_] := refineOnFacts[conjugatedEntries[ConjugateTranspose[m], qb], qb]

conjugatedEntries[a_, qb_] /; ! basisFactsQ[qb] := a
conjugatedEntries[a_, _] := mapValues[conjugateDistributed, a]

conjugateDistributed[x_] /; FreeQ[x, Conjugate[_Plus | _Times | Power[_, _Integer]]] := x
conjugateDistributed[x_] := x //. {
    Conjugate[y : _Plus | _Times] :> Conjugate /@ y,
    Conjugate[Power[y_, n_Integer]] :> Conjugate[y] ^ n
}

(* A yes/no question about every value the facts of qb allow is answered True, False, or undecided[reason, ...]
   with the reason it is not decided: "time" (a test ran out of time), "tests" (no test applies), "fails" (it
   fails for some values, and whether it holds for any is not known), "region" with where it holds, or
   "notnumber" with an entry that is not a number. answer turns the last into Indeterminate. *)
answer[v : True | False] := v
answer[_undecided] := Indeterminate

(* The facts of qb about the given symbols only, for a question that holds no other symbol. *)
linkedFacts[qb_, symbols_List] := symbolDomainsFacts[KeyTake[qb["SymbolDomains"], symbols], Intersection[qb["Parameters"], symbols]]

questionSymbols[e_] := DeleteDuplicates @ Cases[e, $domainSymbolPattern, {0, Infinity}]

(* Reduce and Simplify read the facts given to them and nothing the user set outside: $Assumptions is replaced,
   and Reduce gets every option it has. *)
$reduceOptions = {Backsubstitution -> False, Cubics -> False, GeneratedParameters -> C, Method -> Automatic, Modulus -> 0, Quartics -> False, WorkingPrecision -> Infinity}

reduceOnFacts[claim_, symbols_, facts_] := Block[{$Assumptions = True},
    TimeConstrained[Reduce[claim && facts, symbols, Reals, Sequence @@ $reduceOptions], $regionTimeConstraint, $Aborted]
]

simplifyOnFacts[e_, facts_, limit_] := Block[{$Assumptions = facts}, Simplify[e, TimeConstraint -> {limit, 2 limit}]]

(* An expression Reduce decides without "cannot solve" messages: sums, products and rational powers of the
   symbols and of numbers. *)
algebraicQ[symbols_][x_] := FreeQ[x, Alternatives @@ symbols] && NumericQ[x] || MemberQ[symbols, x] ||
    MatchQ[x, _Plus | _Times] && AllTrue[List @@ x, algebraicQ[symbols]] ||
    MatchQ[x, Power[_, _Integer | _Rational]] && algebraicQ[symbols][First[x]]

(* A cleanup step's result is kept only if it brings in no Indeterminate or infinity its input did not hold. *)
keptStage[old_, new_] := If[FreeQ[old, Indeterminate | _DirectedInfinity] && ! FreeQ[new, Indeterminate | _DirectedInfinity], old, new]


(* points where a question is asked *)

(* Values at a point are computed with the point's non-integer coordinates given to 40 digits, the error tracked
   by the kernel, never exactly: an exact point can make the kernel try to decide the sign of a quantity that is
   exactly zero or one, which it cannot. A value is nonzero when it is not zero within its error. *)
$pointPrecision = 40

numericPoint[point_, precision_] := Replace[point, (s_ -> v_) /; ! IntegerQ[v] :> s -> SetPrecision[v, precision], {1}]

exactNonzeroQ[c_] := With[{n = SetPrecision[c, $pointPrecision]}, NumberQ[n] && n != 0]

(* The sign of a value: 1, 0, -1, "complex", or $Failed when it is not a finite number. *)
numberSign[n_] := Which[! NumberQ[n], $Failed, n == 0, 0, Im[n] != 0, "complex", Re[n] > 0, 1, True, -1]

(* The parts of e that must not vanish at a point for e to be defined there and smooth around it: bases of
   negative and non-integer powers, arguments of Log, Abs, Sign, Arg and the step functions. *)
criticalParts[e_] := DeleteDuplicates @ Join[
    Cases[e, Power[b_, n_] /; ! (IntegerQ[n] && n > 0) :> b, {0, Infinity}],
    Cases[e, (Log | Abs | Sign | Arg | UnitStep | Ramp | RealAbs | Unitize | Clip)[b_, ___] :> b, {0, Infinity}]
]

(* e at a point, its coordinates given to the precision named, or $Failed where it is not defined *)
valueAt[e_, point_, precision_ : $pointPrecision] := With[{p = numericPoint[point, precision]}, With[{v = e /. p},
    If[ NumberQ[v] && AllTrue[criticalParts[e] /. p, exactNonzeroQ], v, $Failed]
]]

(* Exact points inside the facts of qb for the symbols of e, drawn from a seed made of the question, so the same
   question always gets the same points: a real symbol at a point of each of five equal parts of its range (of
   [-3, 3] when it has no bound, so an angle covers a period), an integer symbol rounded inside its range, a
   symbol with no fact at a complex rational; then each finite end of a range, the others staying inside, and
   both corners. A bounded integer symbol with few values is taken at every value, and when every symbol is
   such, the points are complete: {points, complete}. *)
$enumerationLimit = 100

questionPoints[e_, qb_, count_] := With[{symbols = questionSymbols[e], ranges = symbolRanges[qb]},
    With[{counts = integerCount[Lookup[ranges, Key[#], None]] & /@ symbols},
        If[ symbols =!= {} && AllTrue[counts, IntegerQ] && Times @@ counts <= $enumerationLimit,
            {Thread[symbols -> #] & /@ Tuples[integerValues[ranges[#]] & /@ symbols], True},
            {allowedPoints[symbols, ranges, linkedFacts[qb, symbols], Hash[{e, ranges}], count], False}
        ]
    ]
]

integerCount[{Interval[{lo_, hi_}], True}] /; lo > - Infinity && hi < Infinity := Max[0, Floor[hi] - Ceiling[lo] + 1]
integerCount[_] := None

integerValues[{Interval[{lo_, hi_}], True}] := Range[Ceiling[lo], Floor[hi]]

allowedPoints[{}, __] := {{}}

allowedPoints[symbols_, ranges_, facts_, seed_, count_] := Select[
    BlockRandom[
        SeedRandom[seed];
        DeleteDuplicates @ Join[
            Thread[symbols -> #] & /@ Transpose[symbolPoints[#, Lookup[ranges, Key[#], None], count] & /@ symbols],
            If[count > 1, endPoints[symbols, ranges], {}]
        ]
    ],
    TrueQ[facts /. #] &
]

symbolPoints[_, None, count_] := Table[(RandomInteger[{-99, 99}] + I RandomInteger[{-99, 99}]) / 50 + 1 / 7, count]

symbolPoints[_, {Interval[{lo_, hi_}], integer_}, count_] := With[{
    a = If[lo === - Infinity, If[hi === Infinity, -3, hi - 6], lo],
    b = If[hi === Infinity, If[lo === - Infinity, 3, lo + 6], hi]
},
    RandomSample @ Table[
        With[{x = a + (b - a) (k + RandomInteger[{10, 90}] / 100) / count}, If[integer, Clip[Round[x], {Ceiling[lo], Floor[hi]}], x]],
        {k, 0, count - 1}
    ]
]

(* each finite end of each symbol's range with the other symbols at their first point, then all lower and all
   upper ends *)
endPoints[symbols_, ranges_] := With[{
    inner = First /@ (symbolPoints[#, Lookup[ranges, Key[#], None], 1] & /@ symbols),
    ends = Replace[Lookup[ranges, Key[#], None], {{Interval[{lo_, hi_}], _} :> {lo, hi}, None -> {- Infinity, Infinity}}] & /@ symbols
},
    Thread[symbols -> #] & /@ Select[
        Join[
            Catenate @ MapIndexed[Function[{e, i}, ReplacePart[inner, i -> #] & /@ e], ends],
            {ends[[All, 1]], ends[[All, 2]]}
        ],
        FreeQ[#, _DirectedInfinity] &
    ]
]


(* is it zero for every allowed value *)

(* Whether every entry of a list is zero for all values the facts of qb allow. Each entry is first tried by the
   proof alone; of those left, one that is nonzero at one allowed point refutes it at once, and each other is
   tested in full: refined on the facts, proven, or tried at more points, the ends of the ranges among them. *)
zeroVerdict[es_List, qb_] := With[{rest = Select[DeleteCases[es, 0], ! TrueQ[exponentialZeroTest[#]] &]},
    Which[
        rest === {}, True,
        AnyTrue[rest, nonzeroSomewhereQ[#, qb, 1] &], False,
        True, With[{tests = zeroTest[#, qb] & /@ rest},
            Which[
                AllTrue[tests, # === "zero" &], True,
                MemberQ[tests, "nonzero"], False,
                MemberQ[tests, "time"], undecided["time"],
                True, undecided["tests"]
            ]
        ]
    ]
]

zeroOnDomains[es_List, qb_] := answer[zeroVerdict[es, qb]]

nonzeroSomewhereQ[e_, qb_, count_] := AnyTrue[First @ questionPoints[e, qb, count], nonzeroAtQ[e, #] &]

nonzeroAtQ[e_, point_] := With[{v = valueAt[e, point]}, v =!= $Failed && exactNonzeroQ[v]]

(* the sign at a point; a value that reads zero is a zero only when it reads zero again at 80 digits *)
pointSign[e_, point_] := With[{sign = numberSign[valueAt[e, point]]},
    If[sign === 0 && numberSign[valueAt[e, point, 2 $pointPrecision]] =!= 0, $Failed, sign]
]

(* "zero" when proven zero for every allowed value, "nonzero" when nonzero at an allowed point, else "time" or
   "unknown". The facts settle Sqrt, Abs, Log and Conjugate
   of the real symbols first. A bounded integer symbol with few values is checked at every value. *)
zeroTest[e_, qb_] := With[{r = refineOnFacts[e, qb]},
    With[{proof = provenZeroTest[r]},
        If[ proof === True,
            "zero",
            With[{points = questionPoints[e, qb, 5]},
                Which[
                    AnyTrue[First[points], nonzeroAtQ[e, #] &], "nonzero",
                    Last[points], "zero",
                    proof === $Aborted, "time",
                    True, "unknown"
                ]
            ]
        ]
    ]
]

$zeroTestTimeConstraint = 2

(* True when e is proven zero, $Aborted on the time limit, else False. A product is zero when a factor is, a
   positive power when its base is; otherwise trig of linear angles is written in c = Cos and s = Sin of base
   angles and reduced with s^2 = 1 - c^2, and failing that the exponential form is put over one denominator. *)
provenZeroTest[0] := True
provenZeroTest[e_Times] := With[{tests = provenZeroTest /@ List @@ e}, If[MemberQ[tests, True], True, wholeZeroTest[e]]]
provenZeroTest[Power[b_, n_ ? Positive]] := provenZeroTest[b]
provenZeroTest[e_ ? NumericQ] := ! exactNonzeroQ[e] && wholeZeroTest[e]
provenZeroTest[e_] := wholeZeroTest[e]

wholeZeroTest[e_] := Replace[exponentialZeroTest[e], False :> TimeConstrained[trigZeroQ[e], $zeroTestTimeConstraint, $Aborted]]

exponentialZeroTest[0] := True
exponentialZeroTest[e_] := TimeConstrained[Numerator[Together[TrigToExp[e]]] === 0, $zeroTestTimeConstraint, $Aborted]

provenZeroQ[e_] := TrueQ[provenZeroTest[e]]

$trigHeads = Sin | Cos | Tan | Cot | Sec | Csc

trigZeroQ[e_] /; FreeQ[e, $trigHeads | Power[E, _Complex _. | _Complex]] := False
trigZeroQ[e0_] := Module[{c, s},
    With[{e = e0 /. {
        Power[E, x_] /; linearAngleQ[x / I] :> Cos[x / I] + I Sin[x / I],
        Tan[x_] :> Sin[x] / Cos[x], Cot[x_] :> Cos[x] / Sin[x], Sec[x_] :> 1 / Cos[x], Csc[x_] :> 1 / Sin[x]
    }},
        With[{angles = DeleteDuplicates @ Cases[e, (Sin | Cos)[a_] :> a, {0, Infinity}]},
            With[{symbols = questionSymbols[angles]},
                AllTrue[angles, linearAngleQ] && symbols =!= {} && With[{
                    scale = Thread[symbols -> Map[(LCM @@ Denominator[Coefficient[angles, #]]) # &, symbols]]
                },
                    With[{x = TrigExpand[e /. scale] /. {Cos[v_] /; MemberQ[symbols, v] :> c[v], Sin[v_] /; MemberQ[symbols, v] :> s[v]}},
                        Expand[Expand[Numerator[Together[x]]] /. s[v_] ^ n_Integer /; n >= 2 :> s[v] ^ Mod[n, 2] (1 - c[v] ^ 2) ^ Quotient[n, 2]] === 0
                    ]
                ]
            ]
        ]
    ]
]

(* a sum of rational multiples of symbols and a numeric constant *)
linearAngleQ[a_] := With[{symbols = questionSymbols[a]},
    If[ symbols === {},
        NumericQ[a],
        PolynomialQ[a, symbols] && AllTrue[CoefficientRules[a, symbols],
            With[{degree = Total[First[#]]}, degree === 0 && NumericQ[Last[#]] || degree === 1 && MatchQ[Last[#], _Integer | _Rational]] &
        ]
    ]
]


(* is it nonzero for every allowed value *)

(* Whether every entry of a list is nonzero for all values the facts of qb allow, the question asked of a
   denominator. Only a proof says True: a nonzero number, an exponential, a product or a power of nonzero
   parts, a strict sign, or no real root of a one-symbol algebraic expression. A zero at an allowed point says
   False. *)
nonzeroVerdict[es_List, qb_] := With[{tests = nonzeroTest[#, qb] & /@ es},
    Which[
        AllTrue[tests, TrueQ], True,
        MemberQ[tests, False], False,
        MemberQ[tests, undecided["time"]], undecided["time"],
        True, undecided["tests"]
    ]
]

nonzeroOnDomains[es_List, qb_] := answer[nonzeroVerdict[es, qb]]

nonzeroTest[e_, qb_] := With[{r = refineOnFacts[e, qb]}, Which[
    structurallyNonzeroQ[r], True,
    zeroSomewhereQ[e, qb], False,
    MatchQ[signTest[r, qb], "Positive" | "Negative"], True,
    oneSymbolAlgebraicQ[r], Replace[reduceOnFacts[r == 0, questionSymbols[r], linkedFacts[qb, questionSymbols[r]]], {False -> True, $Aborted -> undecided["time"], _ -> undecided["tests"]}],
    True, undecided["tests"]
]]

structurallyNonzeroQ[e_ ? NumericQ] := exactNonzeroQ[e]
structurallyNonzeroQ[Power[E, _]] := True
structurallyNonzeroQ[e_Times] := AllTrue[List @@ e, structurallyNonzeroQ]
structurallyNonzeroQ[Power[b_, _ ? NumericQ]] := structurallyNonzeroQ[b]
structurallyNonzeroQ[_] := False

zeroSomewhereQ[e_, qb_] := AnyTrue[First @ questionPoints[e, qb, 5], pointSign[e, #] === 0 &]

oneSymbolAlgebraicQ[e_] := With[{symbols = questionSymbols[e]}, Length[symbols] === 1 && algebraicQ[symbols][e]]


(* the sign over the allowed values *)

(* The sign of a real expression over the values the facts of qb allow: "Positive", "Nonnegative",
   "Negative", "Nonpositive", "Zero", {"Varies", point, point} with two allowed points of opposite sign,
   {"NotReal", point} with an allowed point where it is not real, or undecided. The values at the sample points
   decide which direction to prove; a proof comes from interval arithmetic on the ranges, from Reduce for an
   algebraic expression in one symbol, or from Simplify under the facts. *)
signTest[e_, qb_, strict_ : True] := With[{ranges = symbolRanges[qb]},
    Replace[intervalSign[e, ranges, strict], None :> With[{r = refineOnFacts[e, qb]},
        Replace[If[r === e, None, intervalSign[r, ranges, strict]], None :> sampledSign[e, r, qb]]
    ]]
]

(* the sign read off interval arithmetic on the ranges, when it settles it; a weak sign only when no strict one
   is asked for. A rational expression is also tried factored, so that a symbol repeated across a sum, as in
   3 p - 2 p^2 = p (3 - 2 p), is counted once per factor. *)
intervalSign[r_, ranges_, strict_] := Replace[
    intervalSignOf[r, ranges, strict],
    None :> If[ rationalQ[r, Keys[ranges]],
        With[{f = Factor[r]}, If[f === r, None, intervalSignOf[f, ranges, strict]]],
        None
    ]
]

rationalQ[r_, symbols_] := With[{x = Together[r]}, PolynomialQ[Numerator[x], symbols] && PolynomialQ[Denominator[x], symbols] && ! FreeQ[r, Alternatives @@ symbols]]

intervalSignOf[r_, ranges_, strict_] := With[{v = intervalValue[r, ranges]}, With[{lo = Min[v], hi = Max[v]},
    Which[
        ! MatchQ[v, _Interval], None,
        lo == hi == 0, "Zero",
        lo > 0, "Positive",
        hi < 0, "Negative",
        ! strict && lo >= 0, "Nonnegative",
        ! strict && hi <= 0, "Nonpositive",
        True, None
    ]
]]

sampledSign[e_, r_, qb_] := With[{pointsComplete = questionPoints[e, qb, 5]},
    With[{signs = DeleteCases[{pointSign[e, #], #} & /@ First[pointsComplete], {$Failed, _}]},
        With[{seen = DeleteDuplicates[signs[[All, 1]]]},
            Which[
                MemberQ[seen, 1] && MemberQ[seen, -1], {"Varies", SelectFirst[signs, First[#] === 1 &][[2]], SelectFirst[signs, First[#] === -1 &][[2]]},
                MemberQ[seen, "complex"], {"NotReal", SelectFirst[signs, First[#] === "complex" &][[2]]},
                Last[pointsComplete] && seen =!= {}, signClass[seen],
                seen === {0}, If[zeroTest[e, qb] === "zero", "Zero", undecided["tests"]],
                seen === {}, undecided["tests"],
                MemberQ[seen, 1], provenSign[r, qb, 1, ! MemberQ[seen, 0]],
                True, provenSign[r, qb, -1, ! MemberQ[seen, 0]]
            ]
        ]
    ]
]

signOnDomains[e_, qb_] := Replace[signTest[e, qb], _undecided -> Indeterminate]

signClass[{0}] := "Zero"
signClass[seen_] := Which[
    SubsetQ[{1}, seen], "Positive", SubsetQ[{-1}, seen], "Negative",
    SubsetQ[{0, 1}, seen], "Nonnegative", True, "Nonpositive"
]

(* the strict sign when strict is True and it is proven, else the weak sign when proven, else undecided *)
provenSign[r_, qb_, direction_, strict_] := Which[
    strict && signProof[r, qb, direction, True] === True, If[direction > 0, "Positive", "Negative"],
    With[{weak = signProof[r, qb, direction, False]}, weak === True], If[direction > 0, "Nonnegative", "Nonpositive"],
    True, undecided["tests"]
]

$signProofTimeConstraint = 2

signProof[r_, qb_, direction_, strict_] := With[{x = direction r, symbols = questionSymbols[r]},
    With[{facts = linkedFacts[qb, symbols]},
        Or[
            intervalSignQ[x, symbolRanges[qb], strict],
            oneSymbolAlgebraicQ[x] && reduceOnFacts[If[strict, x <= 0, x < 0], symbols, facts] === False,
            TrueQ @ simplifyOnFacts[
                If[strict, # > 0, # >= 0] & @ ComplexExpand[x, Complement[symbols, Keys[symbolRanges[qb]]]],
                facts,
                $signProofTimeConstraint
            ]
        ]
    ]
]

(* x on the ranges by interval arithmetic, as an interval; Abs of a part holding a symbol with no range is at
   least 0 *)
intervalValue[x_, ranges_] := Replace[
    (x /. (Abs[y_] /; ! SubsetQ[Keys[ranges], questionSymbols[y]] :> Interval[{0, Infinity}])) /. KeyValueMap[#1 -> First[#2] &, ranges],
    v_ /; NumericQ[v] && realNumberQ[v] :> Interval[{v, v}]
]

intervalSignQ[x_, ranges_, strict_] := With[{v = intervalValue[x, ranges]},
    MatchQ[v, _Interval] && TrueQ[If[strict, Min[v] > 0, Min[v] >= 0]]
]


(* is it nonnegative for every allowed value *)

(* Whether every expression of a list is nonnegative for all values the facts of qb allow. An expression
   negative everywhere gives False; one that is negative or not real somewhere sends the question to Reduce for
   the values where all hold, within a time limit and only for algebraic expressions: none gives False, a region
   gives undecided["region", r]. *)
nonnegativeVerdict[es_List, qb_] := With[{ranges = symbolRanges[qb]}, {weak = intervalSign[#, ranges, False] & /@ es}, Which[
    AllTrue[weak, MatchQ["Positive" | "Nonnegative" | "Zero"]], True,
    MemberQ[weak, "Negative"], False,
    rationalRealQ[Pick[es, weak, None | "Nonpositive"], Keys[ranges]], exactNonnegativeVerdict[es, Pick[es, weak, None | "Nonpositive"], qb],
    True, sampledNonnegativeVerdict[es, qb]
]]

(* Expressions that are ratios of polynomials in real symbols with real coefficients are decided exactly by
   Reduce: no allowed value makes one negative gives True; otherwise the values where all are nonnegative, none
   (False) or a region. *)
exactNonnegativeVerdict[es_, rest_, qb_] := With[{symbols = questionSymbols[es]}, {facts = linkedFacts[qb, symbols]},
    Replace[reduceOnFacts[Or @@ Thread[rest < 0], symbols, facts], {
        False -> True,
        $Aborted -> undecided["time"],
        _ :> With[{r = reduceOnFacts[And @@ Thread[es >= 0], symbols, facts]},
            Which[r === False, False, r === $Aborted || ! FreeQ[r, Reduce], undecided["fails"], True, undecided["region", r]]
        ]
    }]
]

rationalRealQ[es_List, reals_] := With[{symbols = questionSymbols[es]},
    symbols =!= {} && SubsetQ[reals, symbols] && AllTrue[es, rationalRealExpressionQ[symbols]]
]

rationalRealExpressionQ[symbols_][x_] := x === Infinity || MemberQ[symbols, x] || NumericQ[x] && Im[x] == 0 ||
    MatchQ[x, _Plus | _Times] && AllTrue[List @@ x, rationalRealExpressionQ[symbols]] ||
    MatchQ[x, Power[_, _Integer]] && rationalRealExpressionQ[symbols][First[x]]

sampledNonnegativeVerdict[es_List, qb_] := With[{signs = signTest[#, qb, False] & /@ es},
    Which[
        AllTrue[signs, MatchQ["Positive" | "Nonnegative" | "Zero"]], True,
        MemberQ[signs, "Negative"], False,
        MemberQ[signs, "Nonpositive" | {"Varies" | "NotReal", ___}], whereNonnegative[es, qb],
        MemberQ[signs, undecided["time"]], undecided["time"],
        True, undecided["tests"]
    ]
]

nonnegativeOnDomains[es_List, qb_] := answer[nonnegativeVerdict[es, qb]]

$regionTimeConstraint = 1

whereNonnegative[es_List, qb_] := With[{symbols = questionSymbols[es]},
    If[ ! SubsetQ[qb["RealSymbols"], symbols] || ! AllTrue[es, algebraicQ[symbols]],
        undecided["fails"],
        With[{r = reduceOnFacts[And @@ Thread[es >= 0], symbols, linkedFacts[qb, symbols]]},
            Which[
                r === False, False,
                r === $Aborted || ! FreeQ[r, Reduce], undecided["fails"],
                True, undecided["region", r]
            ]
        ]
    ]
]


(* Whether a square matrix is positive semidefinite for all values the facts of qb allow. It must be Hermitian
   first; then its principal minors must be nonnegative (Sylvester's criterion in its form for semidefinite
   matrices, which needs every principal minor, not only the leading ones). That is 2^d - 1 determinants, so it
   is used up to d = 4; a diagonal matrix needs only its diagonal, and a larger one its eigenvalues, when they
   come within a time limit. *)
positiveSemidefiniteVerdict[m0_, qb_] := With[{m = If[realSymmetricQ[m0, qb], m0, refineOnFacts[m0, qb]]}, {
    hermitian = If[realSymmetricQ[m0, qb], True, zeroVerdict[nonzeroEntries[m - conjugateTransposeOnDomains[m, qb]], qb]]
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

(* a matrix equal to its transpose whose entries are ratios of polynomials in real symbols with real
   coefficients: Hermitian with no test *)
realSymmetricQ[m_, qb_] := With[{entries = nonzeroEntries[m]},
    m === Transpose[m] && (AllTrue[entries, NumericQ[#] && Im[#] == 0 &] || rationalRealQ[Select[entries, ! NumericQ[#] &], qb["RealSymbols"]] && AllTrue[Select[entries, NumericQ], Im[#] == 0 &])
]


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
        verdictMessage[head, "undeclared", prop, free]
    ]
]

reportVerdict[_, _, v : True | False] := v
reportVerdict[head_, prop_, undecided[reason_, args___]] := (
    If[ $verdictMessages,
        verdictMessage[head, Replace[reason, {"tests" -> "undecided", "time" -> "undecidedtime", "fails" -> "undecidedfails", "region" -> "undecidedregion"}], prop, args]
    ];
    Indeterminate
)

(* The outermost yes/no property call on a QF object owns the messages of every check run for it: inside its
   scope a message is recorded, not issued, and when the call returns each distinct message is issued once,
   under the head of the object the user asked. Outside any scope a message is issued at once. Other
   properties open no scope, since every check they run is silent and a scope on every property call costs
   time. *)
$verdictScope = None
$verdictLog = {}

SetAttributes[verdictScope, HoldRest]
verdictScope[_, Except["HermitianQ" | "UnitaryQ" | "PhysicalQ" | "UnitalQ" | "TracePreservingQ"], expr_] := expr
verdictScope[_, _, expr_] /; $verdictScope =!= None := expr
verdictScope[head_, _, expr_] := Block[{$verdictScope = head, $verdictLog = {}},
    With[{result = expr},
        Scan[Apply[Message[MessageName[head, #1], ##2] &], DeleteDuplicates[$verdictLog]];
        result
    ]
]

verdictMessage[head_, tag_, args___] := If[ $verdictScope === None,
    Message[MessageName[head, tag], args],
    $verdictLog = Append[$verdictLog, {tag, args}]
]

Scan[
    Function[head,
        MessageName[head, "undecided"] = "`1` could not be decided by the tests for this object; the answer is Indeterminate.";
        MessageName[head, "undecidedtime"] = "`1` could not be decided within the time limit; the answer is Indeterminate.";
        MessageName[head, "undecidedfails"] = "`1` fails for some values the facts allow, and whether it holds for any could not be decided; the answer is Indeterminate.";
        MessageName[head, "undecidedregion"] = "`1` holds only where `2`, not for every value the facts allow; the answer is Indeterminate.";
        MessageName[head, "notnumber"] = "`1` cannot be decided: the entry `2` is not a number; the answer is Indeterminate.";
        MessageName[head, "undeclared"] = "`1` is False because `2` may be complex; it is True when they are real, which \"SymbolDomains\" -> Reals declares."
    ],
    {QuantumState, QuantumOperator, QuantumChannel, QuantumCircuitOperator, QuantumMeasurementOperator, QuantumMeasurement, QuantumHamiltonianOperator}
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
        domains = KeyDrop[Lookup[data, "SymbolDomains", <||>], Keys[rules]]
    },
        If[ output === qb["Output"]["Representations"] &&
                input === qb["Input"]["Representations"] &&
                label === data["Label"],

            QuantumBasis[withRecord[<|data, "ParameterSpec" -> spec|>, domains]],

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

