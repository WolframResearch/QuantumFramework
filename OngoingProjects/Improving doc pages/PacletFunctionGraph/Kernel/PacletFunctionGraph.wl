(* ::Package:: *)

(*
  PacletFunctionGraph.wl
  Build and render a directed function-relation graph for any Wolfram Language paclet.


  TWO-PHASE WORKFLOW
  ──────────────────────────────────────────────────────────────────────
  Phase 1 — Build (once per paclet, slow)
      g = BuildFunctionGraph["/path/to/MyPaclet"]
    Walks DownValues/SubValues of every paclet symbol, classifies argument
    patterns by head, infers flow edges. By default also writes the result
    to <pacletPath>/Resources/FunctionGraph.wl (disable with
    "WriteFile" -> False).

  Phase 2 — Render or query (cheap, repeated)
      PacletOverviewGraph[g]
      SymbolNeighborhoodGraph[g, "QuantumState"]
      PacletFlowEdges[g]
    Pure graph derivations. No file I/O. No paclet load required.

  Persistence
      g = LoadFunctionGraph["/path/to/MyPaclet/Resources/FunctionGraph.wl"]
    Equivalent to Get[file] with type assertion. Use to skip Phase 1.


  GRAPH SCHEMA
  ──────────────────────────────────────────────────────────────────────
  BuildFunctionGraph returns an Association keyed by symbol short-name
  plus a "$Meta" entry:

      <|
          "$Meta" -> <|
              "Paclet"           -> _String,
              "PacletVersion"    -> _String,
              "PrimaryContext"   -> _String,
              "GeneratedAt"      -> _DateObject,
              "GenerationIssues" -> { _String, ... }
          |>,

          symName_String -> <|
              "Kind"       -> "Constructor" | "Operator" | "Transform" |
                              "Predicate"   | "Distance" | "NamedInstance" |
                              "Utility",
              "Accepts"    -> { <|
                                  "Head"        -> _Symbol,
                                  "Role"        -> "Primary" | "Mutation" |
                                                   "Conversion" | "NamedArg" |
                                                   "NamedInstance" |
                                                   "PatternMachinery",
                                  "Description" -> _String
                                |>, ... },
              "OperatesOn" -> None |
                              <| _Symbol -> <| "Head" -> _Symbol,
                                               "Description" -> _String |>,
                                 ... |>,
              "ConsumedBy" -> { _Symbol, ... }
          |>
      |>


  GRAPH MODEL — paclet-only data flow
  ──────────────────────────────────────────────────────────────────────
  Nodes:   ONLY paclet head symbols (e.g. QuantumState, QuantumOperator).
           Built-in heads (Real, Integer, Symbol, List, String, Association,
           ...) are NEVER nodes.

  Edges:   Directed, unlabeled. `A → B` means "A flows into B" — either
           A is an input that B accepts, or A produces something that B
           consumes.

  Sinks:   Functions whose return is a non-paclet scalar (Real, Boolean)
           appear as nodes but have no outgoing edges. Vertex color
           (Kind) is the only signal that a node is a Distance / Predicate /
           scalar-returning function.

  Edge sources (per paclet symbol F):
    1. Accepts (paclet-typed heads only):   for each paclet head H,
       emit  H → F   ("H is consumed by F as a constructor input").
    2. OperatesOn keys (operator-like F):   for each key paclet head K,
       emit  K → F   ("F can be applied to a K").
    3. OperatesOn values:                   for each value paclet head V,
       emit  F → V   ("the application of F produces a V").
*)


BeginPackage["Wolfram`PacletFunctionGraph`"];


BuildFunctionGraph::usage =
    "BuildFunctionGraph[pacletPath] introspects the paclet at pacletPath and " <>
    "returns an Association keyed by symbol name plus a \"$Meta\" entry. As a " <>
    "side effect, writes the result to <pacletPath>/Resources/FunctionGraph.wl " <>
    "unless OutputFile is overridden or WriteFile -> False.";

LoadFunctionGraph::usage =
    "LoadFunctionGraph[file] reads a graph Association previously written by " <>
    "BuildFunctionGraph. Equivalent to Get[file] but with type assertion.";

PacletFlowEdges::usage =
    "PacletFlowEdges[graph] returns a deduplicated list of {fromName, toName} " <>
    "string pairs, each representing a directed flow between two paclet heads.";

SymbolNeighborhoodGraph::usage =
    "SymbolNeighborhoodGraph[graph, symName] renders a Graph[] showing every " <>
    "directed edge that touches symName as source or target. Vertex color " <>
    "encodes Kind. Accepts any Graph option.";

PacletOverviewGraph::usage =
    "PacletOverviewGraph[graph] renders the full paclet flow graph between all " <>
    "paclet head symbols. Option \"IncludeIsolated\" -> True keeps symbols with " <>
    "no edges; default False drops them. Accepts any Graph option.";

$FunctionGraphKindColors::usage =
    "$FunctionGraphKindColors is an Association mapping Kind strings " <>
    "(\"Constructor\", \"Operator\", \"Transform\", \"Predicate\", \"Distance\", " <>
    "\"NamedInstance\", \"Utility\") to RGBColor values used by the graph renderers.";


Begin["`Private`"];


(* ==================================================================
                        COLOR PALETTE
   ================================================================== *)

$FunctionGraphKindColors = <|
    "Constructor"   -> RGBColor[0.7, 0.85, 1],
    "Operator"      -> RGBColor[0.7, 1, 0.7],
    "Transform"     -> RGBColor[1, 0.9, 0.6],
    "Predicate"     -> RGBColor[1, 0.7, 0.7],
    "Distance"      -> RGBColor[0.9, 0.7, 1],
    "NamedInstance" -> RGBColor[0.9, 0.9, 0.9],
    "Utility"       -> RGBColor[0.95, 0.95, 0.95],
    ""              -> RGBColor[0.9, 0.9, 0.9]
|>;


(* ==================================================================
                        NAME UTILITIES
   ================================================================== *)

(* headDisplayName: short rendering of a head for graph node labels.    *)
headDisplayName[s_Symbol] := SymbolName[s];
headDisplayName[s_String] := s;
headDisplayName[other_]   := ToString[other];

(* headInputForm: full InputForm rendering, used in human descriptions. *)
headInputForm[s_Symbol] := SymbolName[s];
headInputForm[expr_]    := ToString[Unevaluated[expr], InputForm];


(* ==================================================================
                        PACLET-HEAD PREDICATES
   ================================================================== *)

isPacletHeadName[name_String, pacletKeys_List] := MemberQ[pacletKeys, name];

isPacletHeadSym[h_Symbol, primaryContext_String] :=
    StringStartsQ[Context[h], primaryContext];
isPacletHeadSym[_, _] := False;


(* ==================================================================
                        PATTERN-MACHINERY HEADS
   These are pattern-syntax heads we never want to emit as graph nodes.
   ================================================================== *)

$PatternMachineryHeads = {
    Pattern, Blank, BlankSequence, BlankNullSequence,
    Except, PatternTest, PatternSequence, Optional, Repeated, RepeatedNull,
    Condition, HoldPattern, Alternatives, OptionsPattern, Verbatim
};


(* ==================================================================
                        PREDICATE -> HEAD MAPPING (extendable)
   ================================================================== *)

$DefaultPredicateHeadNames = <|
    "NumericQ"      -> {"Integer", "Real", "Complex", "Rational"},
    "NumberQ"       -> {"Integer", "Real", "Complex", "Rational"},
    "IntegerQ"      -> {"Integer"},
    "AtomQ"         -> {"Integer", "Real", "Complex", "Rational", "String", "Symbol"},
    "StringQ"       -> {"String"},
    "ListQ"         -> {"List"},
    "VectorQ"       -> {"List"},
    "MatrixQ"       -> {"List"},
    "ArrayQ"        -> {"List", "SparseArray"},
    "TensorQ"       -> {"List", "SparseArray"},
    "AssociationQ"  -> {"Association"},
    "RuleQ"         -> {"Rule", "RuleDelayed"},
    "TrueQ"         -> {"Symbol"},
    "BooleanQ"      -> {"Symbol"},
    "GraphQ"        -> {"Graph"},
    "NumericArrayQ" -> {"NumericArray"}
|>;


(* ==================================================================
                        PATTERN CLASSIFIER
   classifyFirstArg: given a held first-arg pattern, return one or more
   <|"Head" -> _, "Description" -> _|> entries.
   ================================================================== *)

classifyFirstArg[HoldComplete[arg_], predTable_Association] :=
    Replace[Hold[arg],
        {
            Hold[Verbatim[Pattern][_, Verbatim[Blank][h_Symbol]]] :>
                {<|"Head" -> h, "Description" -> "_" <> SymbolName[h]|>},

            Hold[Verbatim[Blank][h_Symbol]] :>
                {<|"Head" -> h, "Description" -> "_" <> SymbolName[h]|>},

            Hold[Verbatim[PatternTest][Verbatim[Pattern][_, Verbatim[Blank][]], pred_]] :>
                classifyPredicate[pred, predTable],

            Hold[Verbatim[PatternTest][Verbatim[Pattern][_, Verbatim[Blank][h_Symbol]], pred_]] :>
                {<|"Head" -> h,
                   "Description" -> "_" <> SymbolName[h] <> "?" <> headInputForm[pred]|>},

            Hold[Verbatim[PatternTest][Verbatim[Blank][], pred_]] :>
                classifyPredicate[pred, predTable],

            Hold[Verbatim[Pattern][_, inner_]] :>
                classifyFirstArg[HoldComplete[inner], predTable],

            Hold[Verbatim[Alternatives][alts__]] :>
                Catenate[Map[classifyFirstArg[HoldComplete[#], predTable] &,
                    Unevaluated[{alts}]]],

            Hold[Verbatim[Optional][p_, ___]] :>
                classifyFirstArg[HoldComplete[p], predTable],

            Hold[Verbatim[Condition][p_, _]] :>
                classifyFirstArg[HoldComplete[p], predTable],

            Hold[Verbatim[HoldPattern][p_]] :>
                classifyFirstArg[HoldComplete[p], predTable],

            Hold[Verbatim[Rule][_String, _] | Verbatim[RuleDelayed][_String, _]] :>
                {<|"Head" -> Rule, "Description" -> "named-argument rule"|>},

            Hold[s_String] :>
                {<|"Head" -> String, "Description" -> "literal: \"" <> s <> "\""|>},

            Hold[i_Integer] :>
                {<|"Head" -> Integer, "Description" -> "literal: " <> ToString[i]|>},

            Hold[other_] :>
                {<|
                    "Head"        -> canonicalHead[Hold[other]],
                    "Description" -> heldDescription[Hold[other]]
                |>}
        }
    ];


(* ==================================================================
                        PREDICATE CLASSIFIER
   classifyPredicate: resolve `?pred` PatternTest to the head(s) it guards.

   Strategy:
     1. User-supplied + default predicate table (string-keyed by name).
     2. Cache miss → auto-detect by walking the predicate's own DownValues:
        a predicate like `somePred[_HeadName] := True` declares its accepted
        head right in the LHS pattern, so we recover {HeadName}.
     3. Still nothing → fall back to Symbol and log an Issue.

   This is what makes BuildFunctionGraph paclet-agnostic: any predicate
   defined with head-pattern DownValues is discovered automatically.
   ================================================================== *)

classifyPredicate[pred_Symbol, predTable_Association] := With[
    {name = SymbolName[pred]},
    With[
        {tableHit  = Lookup[predTable, name, Missing[]]},
        With[
            {heads = Which[
                ! MissingQ[tableHit], Map[Symbol, tableHit],
                True, Replace[autoDetectPredicateHeads[pred],
                    {} -> Missing[]
                ]
            ]},
            If[MissingQ[heads],
                Sow["unresolved predicate " <> name <> " - fallback to Symbol", "Issue"];
                {<|"Head" -> Symbol,
                   "Description" -> "predicate-guarded (?" <> name <> ")"|>},
                Map[<|"Head" -> #, "Description" -> "?" <> name|> &, heads]
            ]
        ]
    ]
];
classifyPredicate[_, _] :=
    {<|"Head" -> Symbol, "Description" -> "anonymous predicate"|>};


(*
  autoDetectPredicateHeads: walk a predicate's DownValues to extract heads.

  Two-pass:
    Pass 1 — heads from the LHS first-arg pattern. Handles `_Head`, `x_Head`,
             `_?pred`, alternatives, and call-form patterns like `MyHead[_]`
             or `x : MyHead[_]`.
    Pass 2 — if Pass 1 returns nothing, walk the RHS for:
             (a) delegated *Q predicate calls (e.g. QuantumOperatorQ[op])
             (b) MatchQ[_, pattern] — extract heads from `pattern`.

  The first assignment to `autoDetectPredicateHeads[pred]` is `{}` so a
  recursive call (predA -> predB -> predA) terminates instead of looping.
*)

autoDetectPredicateHeads[pred_Symbol] := Block[
    {dvs, fromLhs, result},
    autoDetectPredicateHeads[pred] = {};
    dvs = DownValues @@ {pred};
    result = If[dvs === {}, {},
        fromLhs = DeleteDuplicates @ Catenate @ Map[
            extractHeadFromArgPattern[pred, First[#]] &, dvs];
        If[fromLhs =!= {},
            fromLhs,
            DeleteDuplicates @ Catenate @ Map[
                autoDetectHeadsFromRhs[pred, #] &, dvs]
        ]
    ];
    autoDetectPredicateHeads[pred] = result;
    result
];

extractHeadFromArgPattern[pred_Symbol, heldLhs_] := With[
    {argHeld = firstArgOf[pred, heldLhs]},
    If[argHeld === None, {},
        Replace[argHeld,
            {
                (* x_Head / _Head — explicit head blank, possibly with PatternTest *)
                HoldComplete[Verbatim[Pattern][_, Verbatim[Blank][h_Symbol]]] :> {h},
                HoldComplete[Verbatim[Blank][h_Symbol]] :> {h},
                HoldComplete[Verbatim[PatternTest][Verbatim[Pattern][_, Verbatim[Blank][h_Symbol]], _]] :> {h},
                HoldComplete[Verbatim[PatternTest][Verbatim[Blank][h_Symbol], _]] :> {h},

                (* alternatives over head blanks *)
                HoldComplete[Verbatim[Pattern][_, Verbatim[Alternatives][alts__]]] :>
                    Cases[Hold[alts], Hold[Verbatim[Blank][hh_Symbol]] :> hh, {0, Infinity}],
                HoldComplete[Verbatim[Alternatives][alts__]] :>
                    Cases[Hold[alts], Hold[Verbatim[Blank][hh_Symbol]] :> hh, {0, Infinity}],

                (*
                  Call-form pattern: Head[...] or x : Head[...].
                  MUST exclude pattern-machinery symbols — otherwise a fallback DownValue
                  like `MyPredQ[___] := False` (first-arg = BlankNullSequence[]) would
                  cause us to emit Head -> BlankNullSequence as a fake type.
                *)
                HoldComplete[Verbatim[Pattern][_, h_Symbol[___]]] /;
                    ! MemberQ[$PatternMachineryHeads, h] :> {h},
                HoldComplete[h_Symbol[___]] /;
                    ! MemberQ[$PatternMachineryHeads, h] :> {h},

                _ :> {}
            }
        ]
    ]
];

(*
  autoDetectHeadsFromRhs: scan the held RHS of a DownValue for two patterns:
    (a) other_Symbol[___] where other ends in "Q" -- delegated predicate.
    (b) MatchQ[_, pat_] -- the heads in pat are heads this predicate accepts.
  Recurse on (a) one level via autoDetectPredicateHeads (memo breaks cycles).
*)

autoDetectHeadsFromRhs[pred_Symbol, dv_] := With[
    {p = pred, heldRhs = Replace[dv, _[_, rhs_] :> Hold[rhs]]},
    With[{
        delegated = DeleteDuplicates @ Cases[
            heldRhs,
            other_Symbol[___] /; other =!= p && StringEndsQ[SymbolName[other], "Q"] :> other,
            {0, Infinity}
        ],
        (*
          NOTE: `MatchQ[_, pat_]` as a Cases pattern would EVALUATE — MatchQ
          actually answers "does _ match pat_?" instead of matching literal
          MatchQ[...] calls. Use HoldPattern to keep it as a literal match.
        *)
        matchPats = Cases[heldRhs,
            HoldPattern[MatchQ[_, pat_]] :> Hold[pat], {0, Infinity}]
    },
        Join[
            Catenate @ Map[autoDetectPredicateHeads, delegated],
            DeleteDuplicates @ Catenate @ Map[extractHeadsFromMatchPattern, matchPats]
        ]
    ]
];

extractHeadsFromMatchPattern[Hold[pat_]] :=
    DeleteDuplicates @ Cases[
        Hold[pat],
        Verbatim[Blank][h_Symbol] :> h,
        {0, Infinity}
    ];


(* ==================================================================
                        LITERAL HEAD / DESCRIPTION
   ================================================================== *)

canonicalHead[Hold[h_Symbol[___]]]  := h;
canonicalHead[Hold[_String[___]]]   := String;
canonicalHead[Hold[_Integer[___]]]  := Integer;
canonicalHead[Hold[_Real[___]]]     := Real;
canonicalHead[Hold[s_Symbol]]       := Symbol;
canonicalHead[Hold[_String]]        := String;
canonicalHead[Hold[_Integer]]       := Integer;
canonicalHead[Hold[_Real]]          := Real;
canonicalHead[_]                     := Symbol;

heldDescription[Hold[s_String]]      := "literal: \"" <> s <> "\"";
heldDescription[Hold[i_Integer]]     := "literal: " <> ToString[i];
heldDescription[Hold[r_Real]]        := "literal: " <> ToString[r];
heldDescription[Hold[s_Symbol]]      := "symbol: " <> headInputForm[s];
heldDescription[Hold[h_String[___]]] := "named call: \"" <> h <> "\"[...]";
heldDescription[Hold[h_Symbol[___]]] := "call: " <> headInputForm[h] <> "[...]";
heldDescription[Hold[other_]]        := "call form: " <> ToString[Unevaluated[other], InputForm];
heldDescription[_]                   := "call form";


(* ==================================================================
                        LHS / RHS SLICING WITHOUT EVALUATION
   ================================================================== *)

(* firstArgOf: matches `f[arg, ___]` and returns HoldComplete[arg].          *)
firstArgOf[sym_Symbol, heldLhs_] := With[{s = sym},
    Replace[
        Cases[
            {heldLhs},
            HoldPattern[s[a_, ___]] :> HoldComplete[a],
            {0, Infinity}
        ],
        {first_, ___} :> first,
        {0}
    ] /. {} -> None
];

(* innerArgOf: matches `(f | f[_] | f[_]?p | ...)[arg, ___]` for SubValues.  *)
innerArgOf[sym_Symbol, heldLhs_] := With[{s = sym},
    Replace[
        Cases[
            {heldLhs},
            HoldPattern[
                (s | Verbatim[Pattern][_, Verbatim[Blank][s]] |
                 Verbatim[PatternTest][_, _] |
                 Verbatim[Pattern][_, Verbatim[PatternTest][Verbatim[Blank][s], _]]
                )[a_, ___]
            ] :> HoldComplete[a],
            {0, Infinity}
        ],
        {first_, ___} :> first,
        {0}
    ] /. {} -> None
];


(* ==================================================================
                        ROLE INFERENCE
   ================================================================== *)

inferRole[sym_, heldLhs_, heldRhs_, argHeld_, headEntries_, primaryContext_String] := Which[
    AllTrue[headEntries, MemberQ[$PatternMachineryHeads, #["Head"]] &],
        "PatternMachinery",
    namedArgQ[argHeld],
        "NamedArg",
    namedInstanceQ[argHeld],
        "NamedInstance",
    mutationQ[sym, argHeld],
        "Mutation",
    conversionQ[sym, argHeld, primaryContext],
        "Conversion",
    True, "Primary"
];

namedArgQ[HoldComplete[arg_]] := MatchQ[Hold[arg],
    Hold[Verbatim[Rule][_String, _] | Verbatim[RuleDelayed][_String, _]]
];

namedInstanceQ[HoldComplete[arg_]] := MatchQ[Hold[arg],
    Hold[_String] | Hold[_String[___]] | Hold[Verbatim[Alternatives][_String, ___]] |
    Hold[Verbatim[Alternatives][_String, ___][___]] |
    Hold[Verbatim[Pattern][_, Verbatim[Alternatives][_String, ___]]] |
    Hold[Verbatim[Pattern][_, Verbatim[Alternatives][_String, ___][___]]]
];

mutationQ[sym_, HoldComplete[arg_]] := With[{s = sym},
    MatchQ[Hold[arg],
        Hold[Verbatim[Blank][s]] |
        Hold[Verbatim[Pattern][_, Verbatim[Blank][s]]] |
        Hold[Verbatim[PatternTest][Verbatim[Pattern][_, Verbatim[Blank][s]], _]] |
        Hold[Verbatim[PatternTest][Verbatim[Blank][], _ ? (TrueQ[# === Symbol[SymbolName[s] <> "Q"]] &)]]
    ]
];

conversionQ[sym_, argHeld_, primaryContext_String] := With[{s = sym},
    With[{
        argHead = Replace[argHeld,
            HoldComplete[expr_] :> Replace[Hold[expr], {
                Hold[Verbatim[Blank][h_Symbol]]                            :> h,
                Hold[Verbatim[Pattern][_, Verbatim[Blank][h_Symbol]]]      :> h,
                Hold[Verbatim[Pattern][_, Verbatim[Alternatives][alts__]]] :>
                    Alternatives @@ Cases[
                        Hold[alts], Hold[Verbatim[Blank][h_Symbol]] :> h,
                        {0, Infinity}
                    ],
                Hold[Verbatim[Alternatives][alts__]] :>
                    Alternatives @@ Cases[
                        Hold[alts], Hold[Verbatim[Blank][h_Symbol]] :> h,
                        {0, Infinity}
                    ],
                _ :> None
            }]
        ]
    },
        Which[
            argHead === None || argHead === s, False,
            Head[argHead] === Symbol,
                argHead =!= s && isPacletHeadSym[argHead, primaryContext],
            Head[argHead] === Alternatives,
                Length[Cases[argHead,
                    h_Symbol /; h =!= s && isPacletHeadSym[h, primaryContext]]] > 0,
            True, False
        ]
    ]
];


(* ==================================================================
                        NamedInstance COLLAPSE
   Multiple NamedInstance rows are folded into a single catalog entry.
   ================================================================== *)

collapseNamedInstance[entries_] := With[
    {nameRows  = Select[entries, #["Role"] === "NamedInstance" &],
     otherRows = Select[entries, #["Role"] =!= "NamedInstance" &]},
    If[nameRows === {},
        entries,
        Join[
            otherRows,
            {<|
                "Head" -> String, "Role" -> "NamedInstance",
                "Description" -> "named-instance string (one of the catalog entries)",
                "Catalog" -> DeleteDuplicates @ DeleteCases[
                    Map[extractStringLiteral, nameRows], ""]
            |>}
        ]
    ]
];

extractStringLiteral[entry_Association] := With[
    {hits = StringCases[
        Lookup[entry, "Description", ""],
        "\"" ~~ s : Except["\""].. ~~ "\"" :> s, 1]},
    If[hits === {}, "", First[hits]]
];


(* ==================================================================
                        PER-SYMBOL EXTRACTORS
   ================================================================== *)

classifyDownValue[sym_, dv_, predTable_, primaryContext_String] := Block[
    {heldLhs, heldRhs, argHeld, headEntries, role},
    heldLhs = First[dv];
    heldRhs = Replace[dv, _[_, rhs_] :> Hold[rhs]];
    argHeld = firstArgOf[sym, heldLhs];
    If[argHeld === None, Return[{}, Block]];
    headEntries = classifyFirstArg[argHeld, predTable];
    role = inferRole[sym, heldLhs, heldRhs, argHeld, headEntries, primaryContext];
    Map[Append[#, "Role" -> role] &, headEntries]
];

extractAccepts[sym_Symbol, predTable_, primaryContext_String] := With[
    {dvs = DownValues @@ {sym}},
    DeleteDuplicates @ collapseNamedInstance @ DeleteCases[
        Catenate @ Map[
            classifyDownValue[sym, #, predTable, primaryContext] &,
            dvs
        ],
        Nothing
    ]
];

extractOperatesOn[sym_Symbol, predTable_] := With[
    {svs = SubValues @@ {sym}},
    If[svs === {},
        None,
        With[{
            edges = DeleteDuplicates @ Catenate @ Map[
                Function[heldLhs,
                    With[{innerHeld = innerArgOf[sym, heldLhs]},
                        If[innerHeld === None, {},
                           classifyFirstArg[innerHeld, predTable]]
                    ]
                ],
                First /@ svs
            ]
        },
            If[edges === {},
                None,
                AssociationThread[
                    #["Head"] & /@ edges,
                    Map[<|"Head" -> #["Head"], "Description" -> #["Description"]|> &, edges]
                ]
            ]
        ]
    ]
];

(* defaultKindClassifier: substring-based heuristic. Substring matching is
   brittle ("BaseOperatorMap" classifies as Operator); override with the
   "KindClassifier" option of BuildFunctionGraph for paclet-specific rules. *)

defaultKindClassifier[sym_Symbol] := With[{name = SymbolName[sym]},
    Which[
        StringEndsQ[name, "Q"], "Predicate",
        StringContainsQ[name, "Distance"] || StringContainsQ[name, "Similarity"], "Distance",
        StringContainsQ[name, "Transform"] || StringContainsQ[name, "Trace"] || StringContainsQ[name, "Product"], "Transform",
        StringContainsQ[name, "Operator"] || StringContainsQ[name, "Channel"], "Operator",
        True, "Constructor"
    ]
];


(* ==================================================================
                        ConsumedBy REVERSE PASS

   For each symbol s, ConsumedBy[s] = the list of symbols whose entry
   references s in Accepts or in OperatesOn keys.

   Rewritten as a one-pass GroupBy: emit {consumer, producer} pairs
   while walking entries, then group by producer. Avoids the original
   nested-Do / Append loop (O(n^2) growth).
   ================================================================== *)

computeConsumedBy[graph_Association, primaryContext_String] := With[
    {allSyms = DeleteCases[Keys[graph], "$Meta"]},
    With[{
        edges = Catenate @ KeyValueMap[
            Function[{sym, entry},
                With[{
                    referenced = DeleteDuplicates @ Cases[
                        Catenate[{
                            Lookup[#, "Head", Nothing] & /@
                                Lookup[entry, "Accepts", {}],
                            If[entry["OperatesOn"] === None,
                                {}, Keys[entry["OperatesOn"]]]
                        }],
                        h_Symbol /; MemberQ[allSyms, SymbolName[h]] :>
                            SymbolName[h]
                    ]
                },
                    Map[# -> sym &, referenced]
                ]
            ],
            KeyDrop[graph, "$Meta"]
        ]
    },
        With[{consumersOf = GroupBy[edges, First -> Last]},
            AssociationMap[
                Function[s,
                    Append[graph[s],
                        "ConsumedBy" -> DeleteDuplicates @ Map[
                            Symbol[primaryContext <> #] &,
                            Lookup[consumersOf, s, {}]
                        ]
                    ]
                ],
                allSyms
            ]
        ]
    ]
];


(* ==================================================================
                        BuildFunctionGraph HELPERS
   Each helper handles one phase of the build pipeline.
   ================================================================== *)

(* readPacletConfig: load PacletInfo.wl, extract <|Name, Version, Context,
   KernelRoot, SymbolNames|>. Raises via ConfirmAssert on any malformation.  *)

readPacletConfig[pacletPath_String] := Enclose @ Module[
    {pacletAbs, pacletInfoPath, pacletObj, pacletData,
     primaryContext, extensions, kernelExt, kernelOpts,
     kernelRoot, symbolFullNames},

    pacletAbs = ConfirmBy[
        AbsoluteFileName[pacletPath],
        StringQ,
        "BuildFunctionGraph::badpath: path " <> pacletPath <> " does not exist."
    ];

    pacletInfoPath = FileNameJoin[{pacletAbs, "PacletInfo.wl"}];
    ConfirmAssert[FileExistsQ[pacletInfoPath],
        "BuildFunctionGraph::noinfo: PacletInfo.wl not found at " <> pacletInfoPath
    ];

    pacletObj = Get[pacletInfoPath];
    ConfirmAssert[Head[pacletObj] === PacletObject,
        "BuildFunctionGraph::noinfo: " <> pacletInfoPath <>
        " did not return a PacletObject"
    ];
    pacletData = First[pacletObj];

    primaryContext = Lookup[pacletData, "PrimaryContext", $Failed];
    ConfirmAssert[StringQ[primaryContext],
        "BuildFunctionGraph::noctx: paclet has no PrimaryContext"
    ];

    extensions = Lookup[pacletData, "Extensions", {}];
    kernelExt  = SelectFirst[extensions, MatchQ[#, {"Kernel", ___}] &, $Failed];
    ConfirmAssert[kernelExt =!= $Failed,
        "BuildFunctionGraph::nokernel: PacletInfo.wl has no Kernel extension"
    ];
    kernelOpts = Association[Rest[kernelExt]];
    kernelRoot = FileNameJoin[{pacletAbs, Lookup[kernelOpts, "Root", "Kernel"]}];

    symbolFullNames = Lookup[kernelOpts, "Symbols", {}];

    <|
        "Path"           -> pacletAbs,
        "Name"           -> Lookup[pacletData, "Name", "Paclet"],
        "Version"        -> Lookup[pacletData, "Version", "unknown"],
        "PrimaryContext" -> primaryContext,
        "KernelRoot"     -> kernelRoot,
        "SymbolNames"    -> DeleteCases[toBareName /@ symbolFullNames, $Failed]
    |>
];

toBareName[s_Symbol] := SymbolName[s];
toBareName[s_String] := Last[StringSplit[s, "`"]];
toBareName[_] := $Failed;

(* resolveSymbolList: apply "Exclude" and "OnlySymbol" filters.              *)
resolveSymbolList[paclet_Association, excludeList_List, onlySymbol_] := With[
    {filtered = Complement[paclet["SymbolNames"], excludeList]},
    If[onlySymbol === None,
        filtered,
        Select[filtered, # === onlySymbol &]
    ]
];

(* resolveOutputFile: Automatic → <pacletRoot>/Resources/FunctionGraph.wl.    *)
resolveOutputFile[paclet_Association, Automatic] :=
    FileNameJoin[{paclet["Path"], "Resources", "FunctionGraph.wl"}];
resolveOutputFile[_, file_String] := file;

(* introspectPaclet: load the paclet, walk each requested symbol, return
   <|"Symbols" -> rawGraph, "Issues" -> issuesList|>. ConfirmAssert is lexical,
   so this helper must have its own Enclose; the caller ConfirmMatch'es the
   result.                                                                    *)
introspectPaclet[paclet_Association, symbolNames_List, predTable_, kindFn_] :=
    Enclose @ Module[{entries, issuesBag},
        PacletDirectoryLoad[paclet["Path"]];
        Needs[paclet["PrimaryContext"]];
        ConfirmAssert[
            Length[Names[paclet["PrimaryContext"] <> "*"]] > 0,
            "BuildFunctionGraph::noload: paclet " <> paclet["Name"] <>
            " loaded zero symbols"
        ];
        {entries, issuesBag} = Reap[
            DeleteCases[
                Map[
                    introspectOneSymbol[
                        #, kindFn, predTable, paclet["PrimaryContext"]] &,
                    symbolNames
                ],
                Nothing
            ],
            "Issue"
        ];
        <|
            "Symbols" -> Association[entries],
            "Issues"  -> If[issuesBag === {}, {}, First[issuesBag]]
        |>
    ];

(*
  introspectOneSymbol reads DownValues / SubValues / UpValues BEFORE binding
  sym to a value. Doing it the obvious way (`sym = Symbol[fullName]`) triggers
  the symbol's OwnValue if it has one — e.g. `$FunctionGraphKindColors = <|...|>`
  evaluates to the Association and sym is bound to the value, not the symbol.
  All downstream introspection then errors because DownValues / SubValues
  expect a Symbol head.

  We side-step that by parsing the full name through ToExpression with string
  formatting, which gives us the unbound symbol regardless of OwnValues.
  Constants (only OwnValues, no callable values) are skipped — they contribute
  no graph edges anyway.
*)
introspectOneSymbol[symName_String, kindFn_, predTable_, primaryContext_] := With[
    {fullName = primaryContext <> symName},
    With[{
        downV = ToExpression["DownValues[" <> fullName <> "]"],
        subV  = ToExpression["SubValues["  <> fullName <> "]"],
        upV   = ToExpression["UpValues["   <> fullName <> "]"]
    },
        If[downV === {} && subV === {} && upV === {},
            Sow[symName <> " has no callable values - skipping", "Issue"];
            Nothing,

            With[{sym = Symbol[fullName]},
                symName -> <|
                    "Kind"       -> kindFn[sym],
                    "Accepts"    -> extractAccepts[sym, predTable, primaryContext],
                    "OperatesOn" -> extractOperatesOn[sym, predTable],
                    "ConsumedBy" -> {}
                |>
            ]
        ]
    ]
];

(* buildMeta: assemble the "$Meta" entry.                                     *)
buildMeta[paclet_Association, issues_List] := <|
    "Paclet"           -> paclet["Name"],
    "PacletVersion"    -> paclet["Version"],
    "PrimaryContext"   -> paclet["PrimaryContext"],
    "GeneratedAt"      -> Now,
    "GenerationIssues" -> issues
|>;

(* writeIfRequested: side-effect controlled by "WriteFile" option.            *)
writeIfRequested[output_Association, outputFile_String, True] := (
    If[! DirectoryQ[DirectoryName[outputFile]],
        CreateDirectory[DirectoryName[outputFile],
            CreateIntermediateDirectories -> True]];
    Put[output, outputFile];
    output
);
writeIfRequested[output_Association, _, False] := output;
writeIfRequested[output_Association, outputFile_, other_] :=
    writeIfRequested[output, outputFile, TrueQ[other]];


(* ==================================================================
                        BuildFunctionGraph (PUBLIC)
   ================================================================== *)

Options[BuildFunctionGraph] = {
    "OutputFile"        -> Automatic,
    "Exclude"           -> {},
    "OnlySymbol"        -> None,
    "PredicateHeads"    -> <||>,
    "KindClassifier"    -> Automatic,
    "WriteFile"         -> True
};

BuildFunctionGraph[pacletPath_String, opts : OptionsPattern[]] := Enclose @ Module[
    {paclet, symbolNames, predTable, kindFn, outputFile, build, finalGraph, output},

    paclet      = ConfirmMatch[readPacletConfig[pacletPath], _Association];
    symbolNames = resolveSymbolList[
        paclet, OptionValue["Exclude"], OptionValue["OnlySymbol"]];

    predTable = Join[$DefaultPredicateHeadNames, OptionValue["PredicateHeads"]];
    kindFn    = Replace[OptionValue["KindClassifier"],
                    Automatic -> defaultKindClassifier];
    outputFile = resolveOutputFile[paclet, OptionValue["OutputFile"]];

    build      = ConfirmMatch[
        introspectPaclet[paclet, symbolNames, predTable, kindFn],
        _Association
    ];
    finalGraph = computeConsumedBy[build["Symbols"], paclet["PrimaryContext"]];

    output = Prepend[finalGraph, "$Meta" -> buildMeta[paclet, build["Issues"]]];
    writeIfRequested[output, outputFile, OptionValue["WriteFile"]]
];


(* ==================================================================
                        LoadFunctionGraph
   ================================================================== *)

LoadFunctionGraph::nofile  = "File `1` does not exist.";
LoadFunctionGraph::badfile = "File `1` did not return an Association.";

LoadFunctionGraph[file_String] /; ! FileExistsQ[file] :=
    (Message[LoadFunctionGraph::nofile, file]; $Failed);

LoadFunctionGraph[file_String] := With[
    {data = Block[
        {$Context = "Wolfram`PacletFunctionGraph`Sandbox`",
         $ContextPath = {"System`"}},
        Get[file]
    ]},
    If[AssociationQ[data],
        data,
        Message[LoadFunctionGraph::badfile, file]; $Failed
    ]
];


(* ==================================================================
                        EDGE EXTRACTION (paclet-only)
   PacletFlowEdges: dedup list of {fromName, toName} string pairs.
   ================================================================== *)

PacletFlowEdges[graph_Association] := With[
    {pacletKeys = DeleteCases[Keys[graph], "$Meta"]},
    With[{
        edges = Catenate @ KeyValueMap[
            Function[{symName, entry},
                Join[
                    (* 1. Accepts -> F (paclet-typed inputs only) *)
                    Cases[
                        Lookup[entry, "Accepts", {}],
                        a_ /; isPacletHeadName[headDisplayName[a["Head"]], pacletKeys] &&
                              headDisplayName[a["Head"]] =!= symName &&
                              MemberQ[{"Primary", "Mutation", "Conversion"},
                                      Lookup[a, "Role", "Primary"]] :>
                            {headDisplayName[a["Head"]], symName}
                    ],
                    (* 2./3. OperatesOn keys -> F, F -> OperatesOn values *)
                    If[entry["OperatesOn"] === None,
                        {},
                        Join[
                            Cases[
                                Keys[entry["OperatesOn"]],
                                k_ /; isPacletHeadName[headDisplayName[k], pacletKeys] :>
                                    {headDisplayName[k], symName}
                            ],
                            Cases[
                                Values[entry["OperatesOn"]],
                                v_ /; isPacletHeadName[headDisplayName[v["Head"]], pacletKeys] :>
                                    {symName, headDisplayName[v["Head"]]}
                            ]
                        ]
                    ]
                ]
            ],
            KeyDrop[graph, "$Meta"]
        ]
    },
        DeleteDuplicates @ Cases[edges, {_String, _String}]
    ]
];


(* ==================================================================
                        SymbolNeighborhoodGraph
   ================================================================== *)

Options[SymbolNeighborhoodGraph] = Join[
    {},   (* no SymbolNeighborhoodGraph-specific options yet *)
    Options[Graph]
];

SymbolNeighborhoodGraph[graph_Association, symName_String, opts : OptionsPattern[]] := With[
    {pacletKeys = DeleteCases[Keys[graph], "$Meta"],
     touching   = Select[PacletFlowEdges[graph],
                    #[[1]] === symName || #[[2]] === symName &]},

    If[touching === {},
        Return @ Graph[
            {symName},
            VertexLabels -> {symName -> Placed[Style[symName, Bold], Center]},
            FilterRules[{opts}, Options[Graph]]
        ]
    ];

    With[{
        vertices = Union[DeleteDuplicates @ Flatten[touching], {symName}],
        dEdges   = Map[DirectedEdge[#[[1]], #[[2]]] &, touching],
        kindOf   = Function[v,
            If[MemberQ[pacletKeys, v], Lookup[graph[v], "Kind", ""], ""]]
    },
        Graph[
            vertices,
            dEdges,
            FilterRules[{opts}, Options[Graph]],
            VertexLabels -> Map[
                # -> Placed[
                    If[# === symName,
                        Style[#, Bold, FontSize -> 13, FontFamily -> "Helvetica"],
                        Style[#, FontSize -> 11, FontFamily -> "Helvetica"]
                    ],
                    Center
                ] &,
                vertices
            ],
            VertexStyle -> Map[
                # -> Directive[
                    FaceForm @ Lookup[$FunctionGraphKindColors, kindOf[#], LightGray],
                    EdgeForm[If[# === symName,
                        Directive[Black, Thickness[0.005]],
                        Directive[Gray, Thickness[0.001]]
                    ]]
                ] &,
                vertices
            ],
            VertexSize  -> 0.45,
            GraphLayout -> "SpringElectricalEmbedding",
            ImageSize   -> 700,
            PlotLabel   -> Style[symName <> " \[LongDash] flow neighborhood", Bold, 14]
        ]
    ]
];


(* ==================================================================
                        PacletOverviewGraph
   ================================================================== *)

Options[PacletOverviewGraph] = Join[
    {"IncludeIsolated" -> False},
    Options[Graph]
];

PacletOverviewGraph[graph_Association, opts : OptionsPattern[]] := With[
    {pacletKeys      = DeleteCases[Keys[graph], "$Meta"],
     allEdges        = PacletFlowEdges[graph],
     includeIsolated = TrueQ @ OptionValue["IncludeIsolated"]},
    With[{
        vertices = If[includeIsolated,
            Union[DeleteDuplicates @ Flatten[allEdges], pacletKeys],
            DeleteDuplicates @ Flatten[allEdges]
        ],
        dEdges = Map[DirectedEdge[#[[1]], #[[2]]] &, allEdges],
        kindOf = Function[v,
            If[MemberQ[pacletKeys, v], Lookup[graph[v], "Kind", ""], ""]]
    },
        Graph[
            vertices,
            dEdges,
            FilterRules[{opts}, Options[Graph]],
            VertexLabels -> Map[
                # -> Placed[Style[#, FontSize -> 12, FontFamily -> "Helvetica"], Center] &,
                vertices
            ],
            VertexStyle -> Map[
                # -> Directive[
                    FaceForm @ Lookup[$FunctionGraphKindColors, kindOf[#], LightGray],
                    EdgeForm[Directive[Gray, Thickness[0.001]]]
                ] &,
                vertices
            ],
            VertexSize  -> 0.4,
            GraphLayout -> "DiscreteSpiralEmbedding",
            ImageSize   -> 1200,
            PlotLabel   -> Style[
                Lookup[graph["$Meta"], "Paclet", "Paclet"] <>
                    " \[LongDash] paclet flow",
                Bold, 16
            ]
        ]
    ]
];


End[];       (* `Private` *)
EndPackage[];
