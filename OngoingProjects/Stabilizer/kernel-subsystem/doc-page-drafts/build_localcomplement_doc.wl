#!/usr/bin/env wolframscript
(* ============================================================================
   build_localcomplement_doc.wl
   ----------------------------------------------------------------------------
   Generates QuantumFramework/Documentation/English/ReferencePages/Symbols/
       LocalComplement.nb
   Pattern B (function-like). Returns Graph or GraphState.

   Run:  wolframscript -file OngoingProjects/Stabilizer/build_localcomplement_doc.wl
   ============================================================================ *)

$repoRoot  = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$pacletDir = FileNameJoin[{$repoRoot, "QuantumFramework"}];
$docFile   = FileNameJoin[{$pacletDir, "Documentation", "English",
    "ReferencePages", "Symbols", "LocalComplement.nb"}];

Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

SeedRandom[20260511];

(* ===================== Helpers ===================== *)

cid[]  := RandomInteger[{10^6, 10^9 - 1}];
uuid[] := CreateUUID[];

TI[s_String] := StyleBox[s, "TI"];

lcBtn[] := ButtonBox["LocalComplement", BaseStyle -> "Link",
    ButtonData -> "paclet:Wolfram/QuantumFramework/ref/LocalComplement"];

paclLink[sym_String, path_String] := ButtonBox[sym, BaseStyle -> "Link",
    ButtonData -> "paclet:" <> path];

modInfo[] := Cell["   ", "ModInfo", ExpressionUUID -> uuid[]];

inlineCode[boxes_] := Cell[BoxData[boxes], "InlineFormula", ExpressionUUID -> uuid[]];
inlineMath[boxes_] := Cell[
    BoxData[FormBox[boxes, TraditionalForm]],
    "InlineFormula", ExpressionUUID -> uuid[]];

sub[base_, s_] := SubscriptBox[base, s];
sup[base_, e_] := SuperscriptBox[base, e];

SetAttributes[ioBlock, HoldFirst];
ioBlock[expr_, idx_Integer : 1] := With[{idxStr = ToString[idx]},
    Module[{srcStr, parsed, inBoxes, result, outBoxes},
        srcStr = ToString[Unevaluated[expr], InputForm, PageWidth -> Infinity];
        parsed = Quiet @ UsingFrontEnd @ MathLink`CallFrontEnd[
            FrontEnd`UndocumentedTestFEParserPacket[srcStr, False]];
        inBoxes = If[MatchQ[parsed, {_BoxData, _}], First[parsed][[1]], srcStr];
        result = expr;
        outBoxes = ToBoxes[result, StandardForm];
        {
            Cell[BoxData[inBoxes], "Input",
                 CellLabel -> "In[" <> idxStr <> "]:=",
                 CellID -> cid[],
                 ExpressionUUID -> uuid[]],
            Cell[BoxData[outBoxes], "Output",
                 CellLabel -> "Out[" <> idxStr <> "]=",
                 CellID -> cid[],
                 ExpressionUUID -> uuid[]]
        }
    ]
];

exampleSubsection[title_String, content_List] := Cell[CellGroupData[{
    Cell[title, "ExampleSubsection",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ content
}, Open], CellID -> cid[]];

exampleSection[title_String, content_List] := Cell[CellGroupData[{
    Cell[BoxData[InterpretationBox[
        Cell[title, "ExampleSection", ExpressionUUID -> uuid[]],
        $Line = 0; Null]], "ExampleSection",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ content
}, Open], CellID -> cid[]];

exampleSectionEmpty[title_String] := Cell[BoxData[InterpretationBox[
    Cell[title, "ExampleSection", ExpressionUUID -> uuid[]],
    $Line = 0; Null]], "ExampleSection",
    CellID -> cid[], ExpressionUUID -> uuid[]];

delim[] := Cell[BoxData[""], "ExampleDelimiter",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Usage Cell (2 patterns) ===================== *)

lcCall[argBoxes___] := RowBox[{lcBtn[], "[", argBoxes, "]"}];

usageContent = Flatten[{
    modInfo[],
    inlineCode[lcCall[RowBox[{TI["g"], ",", TI["v"]}]]],
    "\[LineSeparator]returns the graph obtained from ", inlineMath[TI["g"]],
    " by complementing every edge between distinct neighbours of ",
    inlineMath[TI["v"]], ".\n",

    modInfo[],
    inlineCode[lcCall[RowBox[{TI["gs"], ",", TI["v"]}]]],
    "\[LineSeparator]applies local complementation to the underlying graph of a ",
    inlineCode[paclLink["GraphState", "Wolfram/QuantumFramework/ref/GraphState"]],
    " ", inlineMath[TI["gs"]], " and returns the resulting ",
    inlineCode["GraphState"], "."
}];

usageCell = Cell[TextData[usageContent], "Usage",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Notes Cells ===================== *)

notesCells = {
    Cell[TextData[{
        inlineCode[lcBtn[]],
        " implements the local complementation operation of Anders & Briegel \
(arXiv:quant-ph/0504117, Definition 1): for the chosen vertex ",
        inlineMath[TI["v"]],
        ", every pair of distinct neighbours ",
        inlineMath[RowBox[{"(", TI["a"], ",", TI["b"], ")"}]],
        " has the edge ",
        inlineMath[RowBox[{TI["a"], "\[Tilde]", TI["b"]}]],
        " toggled (added if absent, removed if present)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Local complementation generates the equivalence of graph states ",
        "under local Clifford operations. The Anders-Briegel theorem ",
        "(AndBri05 Theorem 1) states that two graph states ",
        inlineMath[sub["|G", "1"]],
        " and ",
        inlineMath[sub["|G", "2"]],
        " are local-Clifford-equivalent if and only if their graphs lie in ",
        "the same orbit under successive local complementations."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "On a ", inlineCode[paclLink["GraphState", "Wolfram/QuantumFramework/ref/GraphState"]],
        " input, ", inlineCode[lcBtn[]], " modifies the underlying graph but ",
        "currently keeps the same vertex-operator-Pauli (VOP) assignments. ",
        "Tracked-VOP local complementation per AndBri05 Eq (8) is on the roadmap; ",
        "use the resulting ", inlineCode["GraphState"],
        "'s graph for combinatorial reasoning and the original ",
        inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
        " for unitary evolution."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Local complementation is involutive: ", inlineCode[lcBtn[]],
        " applied twice at the same vertex restores the original graph."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The local-complementation orbit of an ", inlineMath[TI["n"]],
        "-vertex graph is finite (bounded by the ", inlineMath[TI["n"]],
        "-vertex graph count) and is the equivalence class used to classify ",
        "graph states up to local Cliffords; small representatives are ",
        "tabulated in Hein et al. (arXiv:quant-ph/0602096)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "If ", inlineMath[TI["v"]], " has fewer than two neighbours, ",
        inlineCode[lcBtn[]], " returns ", inlineMath[TI["g"]], " unchanged ",
        "(no pairs to toggle)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]]
};

(* ===================== Primary Examples (5) ===================== *)

primaryExamplesContent = Flatten[{
    Cell["Local complementation at an isolated vertex is the identity:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        LocalComplement[Graph[{1, 2, 3}, {}], 1], 1],
    delim[],

    Cell["Apply local complementation to the centre of a 4-vertex star, producing a complete graph on the leaves:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        LocalComplement[StarGraph[4], 1], 2],
    delim[],

    Cell["At the middle of a path graph, local complementation closes the path's outer arms into a triangle:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        Module[{p3},
            p3 = PathGraph[{1, 2, 3}];
            LocalComplement[p3, 2]
        ], 3],
    delim[],

    Cell["Local complementation on a GraphState returns a new GraphState:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        LocalComplement[GraphState[StarGraph[4]], 1], 4],
    delim[],

    Cell["Involutive at a single vertex:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        Module[{g, twice},
            g = StarGraph[5];
            twice = LocalComplement[LocalComplement[g, 1], 1];
            IsomorphicGraphQ[g, twice]
        ], 5]
}];

(* ===================== Scope (2 subsections) ===================== *)

scopeContent = {
    exampleSubsection["Graph input", {
        Cell["Any undirected Graph can be supplied as the first argument:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            LocalComplement[CycleGraph[5], 1], 1],
        delim[],

        Cell["The vertex argument can be any value used as a vertex label, not just an integer:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            LocalComplement[
                Graph[{"a", "b", "c", "d"},
                    UndirectedEdge @@@ {{"a", "b"}, {"a", "c"}, {"a", "d"}}],
                "a"], 1]
    }],

    exampleSubsection["GraphState input", {
        Cell["A GraphState input returns a GraphState whose underlying graph has been complemented:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs, lc},
                gs = GraphState[CycleGraph[4]];
                lc = LocalComplement[gs, 1];
                {Head[lc], lc["EdgeCount"]}
            ], 1]
    }]
};

(* ===================== Generalizations & Extensions (1) ===================== *)

genExtContent = {
    exampleSubsection["Tracked VOP updates", {
        Cell[TextData[{
            "The Anders-Briegel local-complementation theorem produces an explicit ",
            "single-qubit Clifford update on the involved vertices ",
            "(AndBri05 Eq (8)). Tracked VOPs are on the roadmap; for now ",
            inlineCode[lcBtn[]],
            " preserves the input ", inlineCode["\"VOPs\""],
            " field and only modifies the underlying graph:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs, lc},
                gs = GraphState[StarGraph[4]];
                lc = LocalComplement[gs, 1];
                {gs["VOPs"], lc["VOPs"]}
            ], 1]
    }]
};

(* ===================== Options (none) ===================== *)

optionsContent = {
    exampleSubsection["No options", {
        Cell[TextData[{
            inlineCode[lcBtn[]], " takes no options."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]]
    }]
};

(* ===================== Applications (2) ===================== *)

applicationsContent = {
    exampleSubsection["Enumerating a small local-complementation orbit", {
        Cell["Apply local complementation at each vertex of a path graph and collect the resulting graphs:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, orbit},
                g = PathGraph[{1, 2, 3, 4}];
                orbit = DeleteDuplicates[
                    LocalComplement[g, #] & /@ VertexList[g],
                    IsomorphicGraphQ];
                Length[orbit]
            ], 1]
    }],

    exampleSubsection["Reducing a graph state to a canonical form", {
        Cell["Iteratively local-complement at the vertex that decreases the edge count, terminating in a sparse representative:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, reduce, history},
                g = CompleteGraph[4];
                reduce[gr_] := Module[{cands},
                    cands = LocalComplement[gr, #] & /@ VertexList[gr];
                    First[SortBy[Append[cands, gr], EdgeCount]]
                ];
                history = NestList[reduce, g, 4];
                EdgeCount /@ history
            ], 1]
    }]
};

(* ===================== Properties & Relations (4) ===================== *)

propsRelContent = {
    exampleSubsection["Involutive at a single vertex", {
        Cell[TextData[{
            inlineCode[RowBox[{lcBtn[], "[", lcBtn[], "[", TI["g"], ",", TI["v"], "],", TI["v"], "]"}]],
            " is graph-isomorphic to ", inlineMath[TI["g"]], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, twice},
                g = CycleGraph[5];
                twice = LocalComplement[LocalComplement[g, 2], 2];
                IsomorphicGraphQ[g, twice]
            ], 1]
    }],

    exampleSubsection["Preserves the vertex set", {
        Cell["The vertex list is unchanged; only edges are toggled:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, lc},
                g = StarGraph[5];
                lc = LocalComplement[g, 1];
                Sort[VertexList[g]] === Sort[VertexList[lc]]
            ], 1]
    }],

    exampleSubsection["Star -> complete-on-leaves -> star", {
        Cell["At the centre of a star, LC creates a complete graph on the leaves; applying LC at the centre of the resulting graph restores the star:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, lc1, lc2},
                g = StarGraph[5];
                lc1 = LocalComplement[g, 1];
                lc2 = LocalComplement[lc1, 1];
                {EdgeCount[g], EdgeCount[lc1], EdgeCount[lc2]}
            ], 1]
    }],

    exampleSubsection["Edge count of LC at the centre of a star", {
        Cell[TextData[{
            "On the centre of a ", inlineMath[TI["k"]],
            "-leaf star, ", inlineCode[lcBtn[]], " produces ",
            inlineMath[RowBox[{TI["k"], "+", RowBox[{"Binomial", "[", RowBox[{TI["k"], ",", "2"}], "]"}]}]],
            " edges (the original star edges plus all pairs of leaves):"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Table[
                {k, EdgeCount[LocalComplement[StarGraph[k], 1]]},
                {k, 2, 6}
            ], 1]
    }]
};

(* ===================== Possible Issues (3) ===================== *)

possIssuesContent = {
    exampleSubsection["VOPs are not yet tracked on GraphState input", {
        Cell["LocalComplement on a GraphState modifies the graph but leaves the VOP list unchanged; downstream evolution should be done on the corresponding PauliStabilizer if exact unitary tracking is required:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs, lc},
                gs = GraphState[CycleGraph[4]];
                lc = LocalComplement[gs, 1];
                gs["VOPs"] === lc["VOPs"]
            ], 1]
    }],

    exampleSubsection["Vertices with degree less than 2", {
        Cell["A vertex with 0 or 1 neighbours has no pairs to toggle; LocalComplement leaves the graph unchanged:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, lc},
                g = PathGraph[{1, 2, 3}];
                lc = LocalComplement[g, 1];
                IsomorphicGraphQ[g, lc]
            ], 1]
    }],

    exampleSubsection["Vertex must be present in the graph", {
        Cell[TextData[{
            "The vertex argument is forwarded to ",
            inlineCode["AdjacencyList"],
            ", which emits ", inlineCode["AdjacencyList::inv"],
            " if the vertex is not in the graph. Check membership with ",
            inlineCode[RowBox[{"MemberQ", "[", RowBox[{"VertexList", "[", TI["g"], "]"}], ",", TI["v"], "]"}]],
            " before calling on user-supplied input."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]]
    }]
};

(* ===================== Neat Examples (1) ===================== *)

neatContent = {
    exampleSubsection["Visual LC orbit of a 4-cycle", {
        Cell["The local-complementation orbit of a 4-cycle visited from each vertex, collapsed up to graph isomorphism:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{g, orbit},
                g = CycleGraph[4];
                orbit = DeleteDuplicates[
                    LocalComplement[g, #] & /@ VertexList[g],
                    IsomorphicGraphQ];
                GraphicsRow[
                    Graph[#, ImageSize -> 90, VertexLabels -> "Name"] & /@
                        Prepend[orbit, g],
                    ImageSize -> 400]
            ], 1]
    }]
};

(* ===================== See Also ===================== *)

seeAlsoLinks = {
    "GraphState", "PauliStabilizer", "StabilizerStateQ",
    "StabilizerFrame", "CliffordChannel",
    "QuantumState", "QuantumOperator", "QuantumCircuitOperator"};

systemSeeAlso = {"Graph", "EdgeAdd", "EdgeDelete", "AdjacencyList"};

seeAlsoCell = Cell[TextData[Riffle[
    Join[
        Cell[BoxData[ButtonBox[#, BaseStyle -> "Link",
            ButtonData -> "paclet:Wolfram/QuantumFramework/ref/" <> #]],
            "InlineSeeAlsoFunction",
            TaggingRules -> {"PageType" -> "Function"},
            ExpressionUUID -> uuid[]] & /@ seeAlsoLinks,
        Cell[BoxData[ButtonBox[#, BaseStyle -> "Link",
            ButtonData -> "paclet:ref/" <> #]],
            "InlineSeeAlsoFunction",
            TaggingRules -> {"PageType" -> "Function"},
            ExpressionUUID -> uuid[]] & /@ systemSeeAlso
    ],
    Cell[TextData[StyleBox[" \[FilledVerySmallSquare] ", "InlineSeparator"]],
        ExpressionUUID -> uuid[]]
]], "SeeAlso", CellID -> cid[], ExpressionUUID -> uuid[]];

moreAboutCell = Cell[TextData[
    Cell[BoxData[ButtonBox[
        "Wolfram Quantum Computation Framework", BaseStyle -> "Link",
        ButtonData -> "paclet:Wolfram/QuantumFramework/guide/WolframQuantumComputationFramework"
    ]], "InlineFormula", ExpressionUUID -> uuid[]]
], "MoreAbout", CellID -> cid[], ExpressionUUID -> uuid[]];

keywordsCells = Cell[#, "Keywords",
    CellID -> cid[], ExpressionUUID -> uuid[]] & /@ {
    "local complementation", "graph state", "Anders-Briegel",
    "LC orbit", "local Clifford equivalence", "Hein", "graph",
    "neighbour pair toggle"};

(* ===================== Assemble the Notebook ===================== *)

notebook = Notebook[{

    Cell[CellGroupData[{
        Cell["LocalComplement", "ObjectName",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        usageCell,
        Sequence @@ notesCells
    }, Open]],

    Cell[CellGroupData[{
        Cell[TextData[{
            "See Also",
            Cell[BoxData[TemplateBox[{"SeeAlso",
                Cell[BoxData[FrameBox[
                    Cell["Insert links to any related reference (function) pages.",
                        "MoreInfoText"], BaseStyle -> "IFrameBox"]],
                    "MoreInfoTextOuter"]},
                "MoreInfoOpenerButtonTemplate"]], ExpressionUUID -> uuid[]]
        }], "SeeAlsoSection", CellID -> cid[], ExpressionUUID -> uuid[]],
        seeAlsoCell
    }, Open]],

    Cell[CellGroupData[{
        Cell[TextData[{
            "Tech Notes",
            Cell[BoxData[TemplateBox[{"TechNotes",
                Cell[BoxData[FrameBox[
                    Cell["Insert links to related tech notes.", "MoreInfoText"],
                    BaseStyle -> "IFrameBox"]], "MoreInfoTextOuter"]},
                "MoreInfoOpenerButtonTemplate"]], ExpressionUUID -> uuid[]]
        }], "TechNotesSection", CellID -> cid[], ExpressionUUID -> uuid[]],
        Cell["XXXX", "Tutorials", CellID -> cid[], ExpressionUUID -> uuid[]]
    }, Open]],

    Cell[CellGroupData[{
        Cell["Related Guides", "MoreAboutSection",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        moreAboutCell
    }, Open]],

    Cell[CellGroupData[{
        Cell[TextData[{
            "Related Links",
            Cell[BoxData[TemplateBox[{"RelatedLinks",
                Cell[BoxData[FrameBox[
                    Cell["Insert links to any related page, including web pages.",
                        "MoreInfoText"], BaseStyle -> "IFrameBox"]],
                    "MoreInfoTextOuter"]},
                "MoreInfoOpenerButtonTemplate"]], ExpressionUUID -> uuid[]]
        }], "RelatedLinksSection", CellID -> cid[], ExpressionUUID -> uuid[]],
        Cell["XXXX", "RelatedLinks", CellID -> cid[], ExpressionUUID -> uuid[]]
    }, Open]],

    Cell[CellGroupData[{
        Cell[TextData[{
            "Examples Initialization",
            Cell[BoxData[TemplateBox[{"ExamplesInitialization",
                Cell[BoxData[FrameBox[
                    Cell["Input that is to be evaluated before any examples are run, e.g. \
Needs[\[Ellipsis]].", "MoreInfoText"], BaseStyle -> "IFrameBox"]],
                    "MoreInfoTextOuter"]},
                "MoreInfoOpenerButtonTemplate"]], ExpressionUUID -> uuid[]]
        }], "ExamplesInitializationSection",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Cell[BoxData[RowBox[{"Needs", "[", "\"Wolfram`QuantumFramework`\"", "]"}]],
            "ExampleInitialization",
            CellID -> cid[], ExpressionUUID -> uuid[]]
    }, Open]],

    Cell[CellGroupData[{
        Cell[BoxData[InterpretationBox[GridBox[{{
            StyleBox[RowBox[{"Basic", " ", "Examples"}], "PrimaryExamplesSection"],
            ButtonBox[RowBox[{RowBox[{"More", " ", "Examples"}], " ",
                "\[RightTriangle]"}],
                BaseStyle -> "ExtendedExamplesLink",
                ButtonData :> "ExtendedExamples"]
        }}], $Line = 0; Null]], "PrimaryExamplesSection",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ primaryExamplesContent
    }, Open]],

    Cell[CellGroupData[{
        Cell[TextData[{
            "More Examples",
            Cell[BoxData[TemplateBox[{"MoreExamples",
                Cell[BoxData[FrameBox[
                    Cell["Extended examples in standardized sections.",
                        "MoreInfoText"], BaseStyle -> "IFrameBox"]],
                    "MoreInfoTextOuter"]},
                "MoreInfoOpenerButtonTemplate"]], ExpressionUUID -> uuid[]]
        }], "ExtendedExamplesSection",
            CellTags -> "ExtendedExamples",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        exampleSection["Scope", scopeContent],
        exampleSection["Generalizations & Extensions", genExtContent],
        exampleSection["Options", optionsContent],
        exampleSection["Applications", applicationsContent],
        exampleSection["Properties & Relations", propsRelContent],
        exampleSection["Possible Issues", possIssuesContent],
        exampleSectionEmpty["Interactive Examples"],
        exampleSection["Neat Examples", neatContent]
    }, Open]],

    Cell[CellGroupData[{
        Cell["Metadata", "MetadataSection",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Cell[TextData[{
            "New in: ",
            Cell["14.0", "HistoryData", CellTags -> "New",
                ExpressionUUID -> uuid[]],
            " | Modified in: ",
            Cell[" ", "HistoryData", CellTags -> "Modified",
                ExpressionUUID -> uuid[]],
            " | Obsolete in: ",
            Cell[" ", "HistoryData", CellTags -> "Obsolete",
                ExpressionUUID -> uuid[]]
        }], "History", CellID -> cid[], ExpressionUUID -> uuid[]],
        Cell[CellGroupData[{
            Cell[TextData[{
                "Categorization",
                Cell[BoxData[TemplateBox[{"Metadata",
                    Cell[BoxData[FrameBox[
                        Cell["Metadata such as page URI, context, and type of \
documentation page.", "MoreInfoText"], BaseStyle -> "IFrameBox"]],
                        "MoreInfoTextOuter"]},
                    "MoreInfoOpenerButtonTemplate"]], ExpressionUUID -> uuid[]]
            }], "CategorizationSection",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Symbol", "Categorization", CellLabel -> "Entity Type",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Wolfram/QuantumFramework", "Categorization", CellLabel -> "Paclet Name",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Wolfram`QuantumFramework`", "Categorization", CellLabel -> "Context",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Wolfram/QuantumFramework/ref/LocalComplement", "Categorization",
                CellLabel -> "URI",
                CellID -> cid[], ExpressionUUID -> uuid[]]
        }, Closed]],
        Cell[CellGroupData[{
            Cell["Keywords", "KeywordsSection",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Sequence @@ keywordsCells
        }, Closed]],
        Cell[CellGroupData[{
            Cell["Syntax Templates", "TemplatesSection",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template", CellLabel -> "Additional Function Template",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template", CellLabel -> "Arguments Pattern",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template", CellLabel -> "Local Variables",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template", CellLabel -> "Color Equal Signs",
                CellID -> cid[], ExpressionUUID -> uuid[]]
        }, Closed]]
    }, Open]]
},
    WindowSize -> {700, 770},
    WindowMargins -> {{5, Automatic}, {Automatic, 0}},
    TaggingRules -> <|"Paclet" -> "Wolfram/QuantumFramework"|>,
    CellContext -> "Global`",
    FrontEndVersion -> "15.0 for Mac OS X ARM (64-bit) (February 17, 2026)",
    StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb",
        CharacterEncoding -> "UTF-8"],
    ExpressionUUID -> uuid[]
];

Print["Writing notebook to: ", $docFile];
Export[$docFile, notebook, "Notebook"];
Print["Done. File size: ", FileByteCount[$docFile], " bytes"];
