(* ::Package:: *)

(*
  PacletFunctionGraph.wl
  Build and render a directed function-relation graph for any Wolfram Language paclet.

  ──────────────────────────────────────────────────────────────────────
  GRAPH MODEL — paclet-only data flow
  ──────────────────────────────────────────────────────────────────────

  Nodes:   ONLY paclet head symbols (e.g. QuantumState, QuantumOperator).
           Built-in heads (Real, Integer, Symbol, List, String, Association,
           ...) are NEVER nodes.

  Edges:   Directed, unlabeled. `A → B` means "A flows into B" — either
           A is an input that B accepts, or A produces something that B
           consumes. Both interpretations collapse to the same edge: data
           flows from A to B.

  Sinks:   Functions whose return is a non-paclet scalar (Real, Boolean)
           appear as nodes but have no outgoing edges. Visually they are
           dead-end vertices. Vertex color (Kind) is the only signal that
           a node is a Distance / Predicate / scalar-returning function.

  Edge sources (per paclet symbol F):
    1. Accepts (paclet-typed heads only):   for each paclet head H,
       emit  H → F   ("H is consumed by F as a constructor input").
    2. OperatesOn keys (operator-like F):   for each key paclet head K,
       emit  K → F   ("F can be applied to a K").
    3. OperatesOn values:                   for each value paclet head V,
       emit  F → V   ("the application of F produces a V").

  Public API:
    BuildFunctionGraph[pacletPath, opts]         introspect a paclet, return the graph Association
    LoadFunctionGraph[file]                       load a previously-built graph
    PacletFlowEdges[graph]                        list of {fromName, toName} string pairs
    SymbolNeighborhoodGraph[graph, symName]       Graph[] showing all edges touching symName
    PacletOverviewGraph[graph]                    Graph[] of the whole paclet
    $FunctionGraphKindColors                      <|"Kind" -> RGBColor|> palette

  Typical workflow:
    Needs["Wolfram`PacletFunctionGraph`"]
    g = BuildFunctionGraph["/path/to/MyPaclet"]   (* writes Resources/FunctionGraph.wl *)
    PacletOverviewGraph[g]
    SymbolNeighborhoodGraph[g, "MyConstructor"]
*)


BeginPackage["Wolfram`PacletFunctionGraph`"];


BuildFunctionGraph::usage =
    "BuildFunctionGraph[pacletPath] introspects the paclet at pacletPath and " <>
    "returns an Association keyed by symbol name plus a \"$Meta\" entry. As a " <>
    "side effect, writes the result to <pacletPath>/Resources/FunctionGraph.wl " <>
    "unless OutputFile is overridden.";

LoadFunctionGraph::usage =
    "LoadFunctionGraph[file] reads a graph Association previously written by " <>
    "BuildFunctionGraph. Equivalent to Get[file] but with type assertion.";

PacletFlowEdges::usage =
    "PacletFlowEdges[graph] returns a deduplicated list of {fromName, toName} " <>
    "string pairs, each representing a directed flow between two paclet heads.";

SymbolNeighborhoodGraph::usage =
    "SymbolNeighborhoodGraph[graph, symName] renders a Graph[] showing every " <>
    "directed edge that touches symName as source or target. Vertex color " <>
    "encodes Kind.";

PacletOverviewGraph::usage =
    "PacletOverviewGraph[graph] renders the full paclet flow graph between all " <>
    "paclet head symbols. Option \"IncludeIsolated\" -> True keeps symbols with " <>
    "no edges; default False drops them.";

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

shortName[s_Symbol] := SymbolName[s];
shortName[s_String] := s;
shortName[other_]   := ToString[other];

safeName[s_Symbol] := SymbolName[s];
safeName[expr_]    := ToString[Unevaluated[expr], InputForm];


(* ==================================================================
                        PACLET-HEAD PREDICATE
   ================================================================== *)

isPacletHeadName[name_String, pacletKeys_List] := MemberQ[pacletKeys, name];

isPacletHeadSym[h_Symbol, primaryContext_String] := With[
    {ctx = Context @@ {h}},
    StringStartsQ[ctx, primaryContext]
];
isPacletHeadSym[_, _] := False;


(* ==================================================================
                        EDGE EXTRACTION (paclet-only)
   ================================================================== *)

PacletFlowEdges[graph_Association] := Module[
    {pacletKeys, edges},
    pacletKeys = DeleteCases[Keys[graph], "$Meta"];

    edges = Catenate @ KeyValueMap[
        Function[{symName, entry},
            Module[{inEdges, opsInEdges, opsOutEdges},

                (* 1. Accepts → F (paclet-typed inputs only) *)
                inEdges = Cases[
                    Lookup[entry, "Accepts", {}],
                    a_ /; isPacletHeadName[shortName[a["Head"]], pacletKeys] &&
                          shortName[a["Head"]] =!= symName &&
                          MemberQ[{"Primary", "Mutation", "Conversion"},
                                  Lookup[a, "Role", "Primary"]] :>
                        {shortName[a["Head"]], symName}
                ];

                (* 2. OperatesOn keys → F (when F is operator-like) *)
                {opsInEdges, opsOutEdges} = If[
                    Lookup[entry, "OperatesOn"] === None,
                    {{}, {}},
                    Module[{opsAssoc = entry["OperatesOn"]},
                        {
                            Cases[
                                Keys[opsAssoc],
                                k_ /; isPacletHeadName[shortName[k], pacletKeys] :>
                                    {shortName[k], symName}
                            ],
                            Cases[
                                Values[opsAssoc],
                                v_ /; isPacletHeadName[shortName[v["Head"]], pacletKeys] :>
                                    {symName, shortName[v["Head"]]}
                            ]
                        }
                    ]
                ];

                Join[inEdges, opsInEdges, opsOutEdges]
            ]
        ],
        KeyDrop[graph, "$Meta"]
    ];

    DeleteDuplicates @ Cases[edges, {_String, _String}]
];


(* ==================================================================
                        SymbolNeighborhoodGraph
   ================================================================== *)

SymbolNeighborhoodGraph[graph_Association, symName_String, opts___] := Module[
    {pacletKeys, allEdges, touching, vertices, dEdges,
     vertexLabels, vertexStyles, kindOf},

    pacletKeys = DeleteCases[Keys[graph], "$Meta"];
    allEdges = PacletFlowEdges[graph];

    touching = Select[allEdges, #[[1]] === symName || #[[2]] === symName &];

    If[touching === {}, Return[
        Graph[{symName}, VertexLabels -> {symName -> Placed[Style[symName, Bold], Center]}]
    ]];

    vertices = DeleteDuplicates @ Flatten[touching];
    vertices = Union[vertices, {symName}];

    dEdges = Map[DirectedEdge[#[[1]], #[[2]]] &, touching];

    kindOf[v_] := If[MemberQ[pacletKeys, v], Lookup[graph[v], "Kind", ""], ""];

    vertexLabels = Map[
        # -> Placed[
            If[# === symName,
                Style[#, Bold, FontSize -> 13, FontFamily -> "Helvetica"],
                Style[#, FontSize -> 11, FontFamily -> "Helvetica"]
            ],
            Center
        ] &,
        vertices
    ];

    vertexStyles = Map[
        # -> Directive[
            FaceForm @ Lookup[$FunctionGraphKindColors, kindOf[#], LightGray],
            EdgeForm[If[# === symName,
                Directive[Black, Thickness[0.005]],
                Directive[Gray, Thickness[0.001]]
            ]]
        ] &,
        vertices
    ];

    Graph[
        vertices,
        dEdges,
        opts,
        VertexLabels -> vertexLabels,
        VertexStyle  -> vertexStyles,
        VertexSize   -> 0.45,
        GraphLayout  -> "SpringElectricalEmbedding",
        ImageSize    -> 700,
        PlotLabel    -> Style[symName <> " \[LongDash] flow neighborhood", Bold, 14]
    ]
];


(* ==================================================================
                        PacletOverviewGraph
   ================================================================== *)

Options[PacletOverviewGraph] = {"IncludeIsolated" -> False};

(*
  Match opts as a generic Rule sequence (___Rule) rather than OptionsPattern[],
  because callers can pass Graph options (GraphLayout, ImageSize, ...). Using
  OptionValue on a name not registered in Options[PacletOverviewGraph] emits a
  spurious ::nodef message; we look up IncludeIsolated directly via Lookup.
*)

PacletOverviewGraph[graph_Association, opts___Rule] := Module[
    {pacletKeys, allEdges, vertices, dEdges, vertexLabels, vertexStyles, kindOf,
     graphOpts, includeIsolated},

    pacletKeys      = DeleteCases[Keys[graph], "$Meta"];
    allEdges        = PacletFlowEdges[graph];
    includeIsolated = TrueQ @ Lookup[{opts}, "IncludeIsolated", False];

    vertices = DeleteDuplicates @ Flatten[allEdges];
    If[includeIsolated, vertices = Union[vertices, pacletKeys]];

    dEdges = Map[DirectedEdge[#[[1]], #[[2]]] &, allEdges];

    kindOf[v_] := If[MemberQ[pacletKeys, v], Lookup[graph[v], "Kind", ""], ""];

    vertexLabels = Map[
        # -> Placed[Style[#, FontSize -> 12, FontFamily -> "Helvetica"], Center] &,
        vertices
    ];

    vertexStyles = Map[
        # -> Directive[
            FaceForm @ Lookup[$FunctionGraphKindColors, kindOf[#], LightGray],
            EdgeForm[Directive[Gray, Thickness[0.001]]]
        ] &,
        vertices
    ];

    graphOpts = FilterRules[{opts}, Except["IncludeIsolated"]];

    Graph[
        vertices,
        dEdges,
        Sequence @@ graphOpts,
        VertexLabels -> vertexLabels,
        VertexStyle  -> vertexStyles,
        VertexSize   -> 0.4,
        GraphLayout  -> "DiscreteSpiralEmbedding",
        ImageSize    -> 1200,
        PlotLabel    -> Style[
            Lookup[graph["$Meta"], "Paclet", "Paclet"] <> " \[LongDash] paclet flow",
            Bold, 16
        ]
    ]
];


(* ==================================================================
                        LoadFunctionGraph
   ================================================================== *)

LoadFunctionGraph[file_String] := Module[{data},
    If[! FileExistsQ[file],
        Message[LoadFunctionGraph::nofile, file];
        Return[$Failed]
    ];
    data = Block[
        {$Context = "Wolfram`PacletFunctionGraph`Sandbox`", $ContextPath = {"System`"}},
        Get[file]
    ];
    If[! AssociationQ[data],
        Message[LoadFunctionGraph::badfile, file];
        Return[$Failed]
    ];
    data
];

LoadFunctionGraph::nofile  = "File `1` does not exist.";
LoadFunctionGraph::badfile = "File `1` did not return an Association.";


(* ==================================================================
                        BUILD-SIDE INTROSPECTION
   ================================================================== *)

(*
  Builder state is held in module-private bindings inside BuildFunctionGraph
  itself; the helpers below take those bindings as explicit arguments so they
  remain pure and re-entrant (one build at a time is fine, but state never
  leaks between calls).
*)


(* ----- predicate -> head mapping (extendable) ----- *)

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


(* ----- argument-pattern classifier ----- *)

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
                {<|"Head" -> h, "Description" -> "_" <> SymbolName[h] <> "?" <> safeName[pred]|>},

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

(*
  classifyPredicate: resolve a `?somePred` PatternTest to the head(s) it guards.

  Strategy:
    1. Check the user-supplied + default predicate table (string-keyed by name).
    2. Cache miss -> auto-detect by walking the predicate's own DownValues:
       a predicate like `somePred[_HeadName] := True` declares its accepted
       head right in the LHS pattern, so we can recover {HeadName} by
       running the same first-arg-head extraction we use for any function.
    3. Still nothing -> fall back to Symbol and log an Issue.

  This is what makes BuildFunctionGraph paclet-agnostic: any predicate
  defined with head-pattern DownValues is discovered automatically.
*)

classifyPredicate[pred_Symbol, predTable_Association] := Module[
    {name = SymbolName[pred], heads},
    heads = Lookup[predTable, name, Missing[]];
    If[! MissingQ[heads], heads = Map[Symbol, heads]];

    If[MissingQ[heads],
        heads = autoDetectPredicateHeads[pred];
        If[heads === {}, heads = Missing[]]
    ];

    If[MissingQ[heads],
        Sow["unresolved predicate " <> name <> " - fallback to Symbol", "Issue"];
        {<|"Head" -> Symbol, "Description" -> "predicate-guarded (?" <> name <> ")"|>},
        Map[<|"Head" -> #, "Description" -> "?" <> name|> &, heads]
    ]
];
classifyPredicate[_, _] := {<|"Head" -> Symbol, "Description" -> "anonymous predicate"|>};

(*
  autoDetectPredicateHeads: walk the predicate's DownValues and extract every
  paclet/built-in head that appears as the first-arg blank. Memoized on the
  symbol identity for the lifetime of the kernel.
*)

(*
  Two-pass autodetect:
    Pass 1 — extract heads from the predicate's LHS first-arg pattern.
              Handles `_Head`, `x_Head`, `_?pred`, alternatives, and call-form
              patterns like `MyHead[_]` or `x : MyHead[_]`.
    Pass 2 — if Pass 1 returns nothing, walk the RHS for:
              (a) delegated *Q predicate calls (e.g. QuantumOperatorQ[op])
              (b) MatchQ[_, pattern] — extract heads from `pattern`.

  The first assignment to `autoDetectPredicateHeads[pred]` is `{}` so a
  recursive call (predA -> predB -> predA) terminates instead of looping.
*)

autoDetectPredicateHeads[pred_Symbol] := Module[{dvs, fromLhs, result},
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

extractHeadFromArgPattern[pred_Symbol, heldLhs_] := Module[{argHeld},
    argHeld = firstArgOf[pred, heldLhs];
    If[argHeld === None, Return[{}]];
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
];

(*
  autoDetectHeadsFromRhs: scan the held RHS of a DownValue for two patterns:
    (a) other_Symbol[___] where other ends in "Q" -- delegated predicate.
    (b) MatchQ[_, pat_] -- the heads in pat are heads this predicate accepts.
  Recurse on (a) one level via autoDetectPredicateHeads (memo breaks cycles).
*)

autoDetectHeadsFromRhs[pred_Symbol, dv_] := With[{p = pred},
    Module[{heldRhs, delegated, matchPats, headsFromMatchQ},
        heldRhs = Replace[dv, _[_, rhs_] :> Hold[rhs]];

        delegated = DeleteDuplicates @ Cases[
            heldRhs,
            other_Symbol[___] /; other =!= p && StringEndsQ[SymbolName[other], "Q"] :> other,
            {0, Infinity}
        ];

        (*
          NOTE: `MatchQ[_, pat_]` as a Cases pattern would EVALUATE — MatchQ
          actually answers "does _ match pat_?" instead of matching literal
          MatchQ[...] calls. Use HoldPattern to keep it as a literal match.
        *)
        matchPats = Cases[heldRhs,
            HoldPattern[MatchQ[_, pat_]] :> Hold[pat], {0, Infinity}];
        headsFromMatchQ = DeleteDuplicates @ Catenate @ Map[
            extractHeadsFromMatchPattern, matchPats];

        Join[
            Catenate @ Map[autoDetectPredicateHeads, delegated],
            headsFromMatchQ
        ]
    ]
];

extractHeadsFromMatchPattern[Hold[pat_]] :=
    DeleteDuplicates @ Cases[
        Hold[pat],
        Verbatim[Blank][h_Symbol] :> h,
        {0, Infinity}
    ];

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
heldDescription[Hold[s_Symbol]]      := "symbol: " <> safeName[s];
heldDescription[Hold[h_String[___]]] := "named call: \"" <> h <> "\"[...]";
heldDescription[Hold[h_Symbol[___]]] := "call: " <> safeName[h] <> "[...]";
heldDescription[Hold[other_]]        := "call form: " <> ToString[Unevaluated[other], InputForm];
heldDescription[_]                   := "call form";


(* ----- LHS / RHS slicing without releasing held content ----- *)

firstArgOf[sym_Symbol, heldLhs_] := With[{s = sym},
    Module[{matches},
        matches = Cases[
            {heldLhs},
            HoldPattern[s[a_, ___]] :> HoldComplete[a],
            {0, Infinity}
        ];
        If[matches === {}, None, First[matches]]
    ]
];

innerArgOf[sym_Symbol, heldLhs_] := With[{s = sym},
    Module[{matches},
        matches = Cases[
            {heldLhs},
            HoldPattern[
                (s | Verbatim[Pattern][_, Verbatim[Blank][s]] |
                 Verbatim[PatternTest][_, _] |
                 Verbatim[Pattern][_, Verbatim[PatternTest][Verbatim[Blank][s], _]]
                )[a_, ___]
            ] :> HoldComplete[a],
            {0, Infinity}
        ];
        If[matches === {}, None, First[matches]]
    ]
];


(* ----- Role inference ----- *)

$PatternMachineryHeads = {Pattern, Blank, BlankSequence, BlankNullSequence,
    Except, PatternTest, PatternSequence, Optional, Repeated, RepeatedNull,
    Condition, HoldPattern, Alternatives, OptionsPattern, Verbatim};

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
    Module[{argHead},
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
        ];
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


(* ----- collapsing NamedInstance rows ----- *)

collapseNamedInstance[entries_] := Module[{nameRows, otherRows, catalog},
    nameRows  = Select[entries, #["Role"] === "NamedInstance" &];
    otherRows = Select[entries, #["Role"] =!= "NamedInstance" &];
    If[nameRows === {}, Return[entries]];
    catalog = DeleteDuplicates @ DeleteCases[
        Map[extractStringLiteral, nameRows],
        ""
    ];
    Join[
        otherRows,
        {<|
            "Head" -> String, "Role" -> "NamedInstance",
            "Description" -> "named-instance string (one of the catalog entries)",
            "Catalog" -> catalog
        |>}
    ]
];

extractStringLiteral[entry_Association] := Module[{desc, hits},
    desc = Lookup[entry, "Description", ""];
    hits = StringCases[desc, "\"" ~~ s : Except["\""].. ~~ "\"" :> s, 1];
    If[hits === {}, "", First[hits]]
];


(* ----- per-symbol extractors ----- *)

classifyDownValue[sym_, dv_, predTable_, primaryContext_String] := Module[
    {heldLhs, heldRhs, argHeld, headEntries, role},
    heldLhs = First[dv];
    heldRhs = Replace[dv, _[_, rhs_] :> Hold[rhs]];
    argHeld = firstArgOf[sym, heldLhs];
    If[argHeld === None, Return[{}]];
    headEntries = classifyFirstArg[argHeld, predTable];
    role = inferRole[sym, heldLhs, heldRhs, argHeld, headEntries, primaryContext];
    Map[Append[#, "Role" -> role] &, headEntries]
];

extractAccepts[sym_Symbol, predTable_, primaryContext_String] := With[{symVal = sym},
    Module[{dvs, entries},
        dvs = DownValues @@ {sym};
        entries = Catenate @ Map[
            classifyDownValue[symVal, #, predTable, primaryContext] &,
            dvs
        ];
        entries = DeleteCases[entries, Nothing];
        DeleteDuplicates @ collapseNamedInstance[entries]
    ]
];

extractOperatesOn[sym_Symbol, predTable_] := With[{symVal = sym},
    Module[{svs, heldLhsList, edges},
        svs = SubValues @@ {sym};
        If[svs === {}, Return[None]];
        heldLhsList = First /@ svs;
        edges = Catenate @ Map[
            Function[heldLhs,
                With[{innerHeld = innerArgOf[symVal, heldLhs]},
                    If[innerHeld === None, {},
                       classifyFirstArg[innerHeld, predTable]]
                ]
            ],
            heldLhsList
        ];
        edges = DeleteDuplicates[edges];
        If[edges === {}, None,
            AssociationThread[
                #["Head"] & /@ edges,
                Map[<|"Head" -> #["Head"], "Description" -> #["Description"]|> &, edges]
            ]
        ]
    ]
];

(*
  extractKind: defaults to "Constructor". Callers can supply a custom
  classifier via the "KindClassifier" option to BuildFunctionGraph if
  the paclet uses different naming conventions.
*)

defaultKindClassifier[sym_Symbol] := Module[{name = SymbolName[sym]},
    Which[
        StringEndsQ[name, "Q"], "Predicate",
        StringContainsQ[name, "Distance"] || StringContainsQ[name, "Similarity"], "Distance",
        StringContainsQ[name, "Transform"] || StringContainsQ[name, "Trace"] || StringContainsQ[name, "Product"], "Transform",
        StringContainsQ[name, "Operator"] || StringContainsQ[name, "Channel"], "Operator",
        True, "Constructor"
    ]
];


(* ----- ConsumedBy reverse pass ----- *)

computeConsumedBy[graph_Association, primaryContext_String] := Module[
    {consumers, allSyms = Keys[graph]},
    consumers = AssociationMap[Function[s, {}], allSyms];
    Do[
        Module[{entry = graph[sym], referenced},
            referenced = DeleteDuplicates @ Catenate[{
                Lookup[#, "Head", Nothing] & /@ Lookup[entry, "Accepts", {}],
                If[Lookup[entry, "OperatesOn"] === None, {},
                   Keys[Lookup[entry, "OperatesOn"]]]
            }];
            Do[
                If[Head[ref] === Symbol && MemberQ[allSyms, SymbolName[ref]],
                    consumers[SymbolName[ref]] = Append[
                        consumers[SymbolName[ref]],
                        Symbol[primaryContext <> sym]
                    ]
                ],
                {ref, referenced}
            ]
        ],
        {sym, allSyms}
    ];
    AssociationMap[
        Function[s, Append[graph[s], "ConsumedBy" -> DeleteDuplicates[consumers[s]]]],
        allSyms
    ]
];


(* ==================================================================
                        BuildFunctionGraph (public)
   ================================================================== *)

Options[BuildFunctionGraph] = {
    "OutputFile"        -> Automatic,
    "Exclude"           -> {},
    "OnlySymbol"        -> None,
    "Fixtures"          -> <||>,
    "PredicateHeads"    -> <||>,
    "KindClassifier"    -> Automatic,
    "WriteFile"         -> True
};

BuildFunctionGraph::badpath = "Paclet path `1` does not exist.";
BuildFunctionGraph::noinfo  = "PacletInfo.wl not found at `1`.";
BuildFunctionGraph::nokernel = "PacletInfo.wl has no Kernel extension.";
BuildFunctionGraph::noctx   = "Paclet has no PrimaryContext.";
BuildFunctionGraph::noload  = "Paclet `1` loaded zero symbols.";

BuildFunctionGraph[pacletPath_String, opts : OptionsPattern[]] := Module[
    {
        pacletAbs, pacletInfoPath, pacletObj, pacletData,
        pacletName, pacletVersion, primaryContext,
        extensions, kernelExt, kernelOpts, kernelRoot,
        symbolFullNames, symbolNames, excludeList, onlySymbol,
        predTable, kindFn, outputFile, writeFile,
        loadedCount, issuesBag, entries, rawGraph, finalGraph, meta, output
    },

    (* --- validate inputs ---------------------------------------------- *)
    pacletAbs = AbsoluteFileName[pacletPath];
    If[pacletAbs === $Failed,
        Message[BuildFunctionGraph::badpath, pacletPath]; Return[$Failed]];

    pacletInfoPath = FileNameJoin[{pacletAbs, "PacletInfo.wl"}];
    If[! FileExistsQ[pacletInfoPath],
        Message[BuildFunctionGraph::noinfo, pacletInfoPath]; Return[$Failed]];

    pacletObj = Get[pacletInfoPath];
    If[Head[pacletObj] =!= PacletObject,
        Message[BuildFunctionGraph::noinfo, pacletInfoPath]; Return[$Failed]];
    pacletData = First[pacletObj];

    pacletName     = Lookup[pacletData, "Name", "Paclet"];
    pacletVersion  = Lookup[pacletData, "Version", "unknown"];
    primaryContext = Lookup[pacletData, "PrimaryContext", $Failed];
    If[primaryContext === $Failed,
        Message[BuildFunctionGraph::noctx]; Return[$Failed]];

    extensions = Lookup[pacletData, "Extensions", {}];
    kernelExt  = SelectFirst[extensions, MatchQ[#, {"Kernel", ___}] &, $Failed];
    If[kernelExt === $Failed,
        Message[BuildFunctionGraph::nokernel]; Return[$Failed]];
    kernelOpts = Association[Rest[kernelExt]];
    kernelRoot = FileNameJoin[{pacletAbs, Lookup[kernelOpts, "Root", "Kernel"]}];

    symbolFullNames = Lookup[kernelOpts, "Symbols", {}];
    symbolNames     = DeleteCases[toBareName /@ symbolFullNames, $Failed];

    excludeList = OptionValue["Exclude"];
    onlySymbol  = OptionValue["OnlySymbol"];
    symbolNames = Complement[symbolNames, excludeList];
    If[onlySymbol =!= None,
        symbolNames = Select[symbolNames, # === onlySymbol &]];

    (* --- merge predicate table & kind classifier ---------------------- *)
    predTable = Join[$DefaultPredicateHeadNames, OptionValue["PredicateHeads"]];
    kindFn    = OptionValue["KindClassifier"];
    If[kindFn === Automatic, kindFn = defaultKindClassifier];

    (* --- output target ----------------------------------------------- *)
    writeFile  = TrueQ[OptionValue["WriteFile"]];
    outputFile = OptionValue["OutputFile"];
    If[outputFile === Automatic,
        outputFile = FileNameJoin[{pacletAbs, "Resources", "FunctionGraph.wl"}]];

    (* --- load the target paclet --------------------------------------- *)
    PacletDirectoryLoad[pacletAbs];
    Needs[primaryContext];
    loadedCount = Length[Names[primaryContext <> "*"]];
    If[loadedCount === 0,
        Message[BuildFunctionGraph::noload, pacletName]; Return[$Failed]];

    (* --- per-symbol extraction (sowed issues captured by Reap) -------- *)
    {entries, issuesBag} = Reap[
        DeleteCases[
            Map[
                introspectOneSymbol[#, kindFn, predTable, primaryContext] &,
                symbolNames
            ],
            Nothing
        ],
        "Issue"
    ];
    issuesBag = If[issuesBag === {}, {}, First[issuesBag]];

    rawGraph   = Association[entries];
    finalGraph = computeConsumedBy[rawGraph, primaryContext];

    meta = <|
        "Paclet"           -> pacletName,
        "PacletVersion"    -> pacletVersion,
        "PrimaryContext"   -> primaryContext,
        "GeneratedAt"      -> Now,
        "GenerationIssues" -> issuesBag
    |>;

    output = Prepend[finalGraph, "$Meta" -> meta];

    If[writeFile,
        If[! DirectoryQ[DirectoryName[outputFile]],
            CreateDirectory[DirectoryName[outputFile], CreateIntermediateDirectories -> True]];
        Put[output, outputFile]
    ];

    output
];

toBareName[s_Symbol] := SymbolName[s];
toBareName[s_String] := Last[StringSplit[s, "`"]];
toBareName[_] := $Failed;

(*
  introspectOneSymbol reads DownValues / SubValues / UpValues BEFORE binding sym
  to a value. Doing it the obvious way (`sym = Symbol[fullName]`) triggers the
  symbol's OwnValue if it has one — e.g. `$FunctionGraphKindColors = <|...|>`
  evaluates to the Association and sym is bound to the value, not the symbol.
  All downstream introspection then errors because DownValues / SubValues expect
  a Symbol head.

  We side-step that by parsing the full name through ToExpression with HoldComplete
  formatting, which gives us the unbound symbol regardless of OwnValues. Constants
  (only OwnValues, no callable values) are skipped — they contribute no graph
  edges anyway.
*)

introspectOneSymbol[symName_String, kindFn_, predTable_, primaryContext_] := Module[
    {fullName, downV, subV, upV, sym},
    fullName = primaryContext <> symName;

    (* Read value-lists via ToExpression-on-string to avoid OwnValue evaluation. *)
    downV = ToExpression["DownValues[" <> fullName <> "]"];
    subV  = ToExpression["SubValues["  <> fullName <> "]"];
    upV   = ToExpression["UpValues["   <> fullName <> "]"];

    Which[
        downV === {} && subV === {} && upV === {},
            Sow[symName <> " has no callable values - skipping", "Issue"];
            Return[Nothing],
        True, Null
    ];

    sym = Symbol[fullName];

    symName -> <|
        "Kind"       -> kindFn[sym],
        "Accepts"    -> extractAccepts[sym, predTable, primaryContext],
        "OperatesOn" -> extractOperatesOn[sym, predTable],
        "ConsumedBy" -> {}
    |>
];


End[];       (* `Private` *)
EndPackage[];
