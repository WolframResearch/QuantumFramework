#!/usr/bin/env wolframscript
(* ============================================================================
   build_graphstate_doc.wl
   ----------------------------------------------------------------------------
   Generates QuantumFramework/Documentation/English/ReferencePages/Symbols/
       GraphState.nb
   Pattern A (object-head with custom display, ScriptCapitalG summary box).

   Run:  wolframscript -file OngoingProjects/Stabilizer/build_graphstate_doc.wl
   ============================================================================ *)

$repoRoot  = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$pacletDir = FileNameJoin[{$repoRoot, "QuantumFramework"}];
$docFile   = FileNameJoin[{$pacletDir, "Documentation", "English",
    "ReferencePages", "Symbols", "GraphState.nb"}];

Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

SeedRandom[20260511];

(* ===================== Helpers ===================== *)

cid[]  := RandomInteger[{10^6, 10^9 - 1}];
uuid[] := CreateUUID[];

TI[s_String] := StyleBox[s, "TI"];

gsBtn[] := ButtonBox["GraphState", BaseStyle -> "Link",
    ButtonData -> "paclet:Wolfram/QuantumFramework/ref/GraphState"];

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

gsCall[argBoxes___] := RowBox[{gsBtn[], "[", argBoxes, "]"}];

usageContent = Flatten[{
    modInfo[],
    inlineCode[gsCall[TI["g"]]],
    "\[LineSeparator]represents the graph state of the undirected graph ",
    inlineMath[TI["g"]], " with identity vertex-operator-Paulis (VOPs) on every vertex.\n",

    modInfo[],
    inlineCode[gsCall[TI["ps"]]],
    "\[LineSeparator]converts a graph-form ",
    inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
    " ", inlineMath[TI["ps"]], " into the corresponding ", inlineCode[gsBtn[]],
    "; the input must have ",
    inlineMath[sub["X", TI["i"]]],
    " at the diagonal of every stabilizer row and ",
    inlineMath[RowBox[{sub["Z", TI["j"]], " or ", sub["I", TI["j"]]}]],
    " elsewhere."
}];

usageCell = Cell[TextData[usageContent], "Usage",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Notes Cells ===================== *)

notesCells = {
    Cell[TextData[{
        "A ", inlineCode[gsBtn[]],
        " stores a stabilizer state in the graph-state representation introduced ",
        "by Anders & Briegel (arXiv:quant-ph/0504117). For a graph ",
        inlineMath[RowBox[{TI["G"], "=", RowBox[{"(", TI["V"], ",", TI["E"], ")"}]}]],
        " the associated state is ",
        inlineMath[RowBox[{
            TemplateBox[{TI["G"]}, "Ket"], "=",
            RowBox[{UnderscriptBox["\[Product]", RowBox[{TI["e"], "\[Element]", TI["E"]}]],
                sub["CZ", TI["e"]]}],
            TemplateBox[{sup["+", RowBox[{"\[CircleTimes]", TI["n"]}]]}, "Ket"]
        }]],
        ", i.e. the application of a controlled-",
        inlineMath["Z"], " on every edge to the all-plus state."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The internal representation is the Association ",
        inlineCode[RowBox[{"<|", RowBox[{
            "\"Graph\"", "\[Rule]", TI["g"], ",", " ",
            "\"VOPs\"", "\[Rule]", RowBox[{"{", RowBox[{"0", ",", "\[Ellipsis]"}], "}"}]
        }], "|>"}]],
        ". The ", inlineCode["\"VOPs\""], " field stores a list of integer indices ",
        "into the 24-element single-qubit Clifford group; currently only ",
        "the identity (index ", inlineCode["0"], ") is exercised."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The stabilizer generators of the graph state are read off the graph: ",
        "for each vertex ", inlineMath[TI["i"]], ",",
        inlineMath[RowBox[{sub["K", TI["i"]], "=",
            RowBox[{sub["X", TI["i"]], "\[CircleTimes]",
                RowBox[{UnderscriptBox["\[Product]",
                    RowBox[{TI["j"], "\[Element]", "N", "(", TI["i"], ")"}]],
                    sub["Z", TI["j"]]}]}]}]],
        ", where ", inlineMath[RowBox[{"N", "(", TI["i"], ")"}]],
        " is the neighbour set of ", inlineMath[TI["i"]], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The graph-state representation has memory cost ",
        inlineMath[RowBox[{"O", "(", TI["n"], "\[CenterDot]",
            OverscriptBox[TI["d"], "_"], ")"}]],
        " where ", inlineMath[OverscriptBox[TI["d"], "_"]],
        " is the average vertex degree, which is asymptotically smaller than the ",
        inlineMath[RowBox[{"O", "(", sup[TI["n"], "2"], ")"}]],
        " tableau for sparsely-connected codes."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The ", inlineCode[gsCall[TI["ps"]]],
        " conversion succeeds only on graph-form stabilizers (one ",
        inlineMath["X"], " on the diagonal, ",
        inlineMath["Z"], "/", inlineMath["I"], " elsewhere); other inputs emit ",
        inlineCode["GraphState::nongraph"],
        " and return ", inlineCode["$Failed"],
        ". Use a local-Clifford transformation ",
        "(AndBri05 Lemma 1) to bring an arbitrary stabilizer state into graph form ",
        "before calling."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Two graph states (with identity VOPs) are local-Clifford-equivalent if ",
        "and only if their graphs are related by a sequence of local complementations ",
        "(",
        inlineCode[paclLink["LocalComplement", "Wolfram/QuantumFramework/ref/LocalComplement"]],
        "); this is the basis of the graph-state classification of stabilizer states ",
        "(AndBri05 Theorem 1, Hein-Eisert-Briegel arXiv:quant-ph/0602096)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Property access uses the dispatch syntax ",
        inlineCode[RowBox[{TI["gs"], "[", "\"", TI["property"], "\"", "]"}]],
        "; the full list of recognised properties is given by ",
        inlineCode[RowBox[{TI["gs"], "[", "\"Properties\"", "]"}]],
        ". Conversion to a ", inlineCode["PauliStabilizer"], " is done via ",
        inlineCode[RowBox[{TI["gs"], "[", "\"PauliStabilizer\"", "]"}]], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]]
};

(* ===================== Primary Examples (5) — Pattern A ===================== *)

primaryExamplesContent = Flatten[{
    Cell["The graph state of a 4-cycle:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[GraphState[CycleGraph[4]], 1],
    delim[],

    Cell["The graph state of a 3-path:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[GraphState[PathGraph[{1, 2, 3}]], 2],
    delim[],

    Cell["The graph state of a 4-vertex star (GHZ-equivalent under local Cliffords):",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[GraphState[StarGraph[4]], 3],
    delim[],

    Cell["A graph-form PauliStabilizer is converted directly:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[GraphState[PauliStabilizer[{"XZI", "ZXZ", "IZX"}]], 4],
    delim[],

    Cell["Inspecting stabilizers and the underlying graph:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        Module[{gs},
            gs = GraphState[CycleGraph[4]];
            {gs["VertexCount"], gs["EdgeCount"], gs["Stabilizers"]}
        ], 5]
}];

(* ===================== Scope (4 subsections) ===================== *)

scopeContent = {
    exampleSubsection["Graph input", {
        Cell["Any undirected graph can be used as the underlying graph:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[GraphState[CompleteGraph[3]], 1],
        delim[],

        Cell["Disconnected graphs produce product graph states:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[Graph[{1, 2, 3, 4}, {UndirectedEdge[1, 2], UndirectedEdge[3, 4]}]], 1]
    }],

    exampleSubsection["PauliStabilizer input", {
        Cell["Stabilizer rows of a 4-cycle graph state in graph form:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[PauliStabilizer[{"XZIZ", "ZXZI", "IZXZ", "ZIZX"}]], 1]
    }],

    exampleSubsection["Property accessors", {
        Cell[TextData[{
            "Every recognised property is enumerated by ",
            inlineCode[RowBox[{TI["gs"], "[", "\"Properties\"", "]"}]], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[CycleGraph[3]]["Properties"], 1],
        delim[],

        Cell["Vertices, edges, and adjacency:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs},
                gs = GraphState[StarGraph[4]];
                <|
                    "Vertices" -> gs["Vertices"],
                    "Edges" -> gs["Edges"],
                    "VertexCount" -> gs["VertexCount"],
                    "EdgeCount" -> gs["EdgeCount"]
                |>
            ], 1]
    }],

    exampleSubsection["Stabilizer generators", {
        Cell[TextData[{
            "For a path graph, the stabilizer at vertex ", inlineMath[TI["i"]],
            " is ", inlineMath[RowBox[{sub["X", TI["i"]], "\[CircleTimes]",
                sub["Z", RowBox[{"N", "(", TI["i"], ")"}]]}]], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[PathGraph[{1, 2, 3, 4}]]["Stabilizers"], 1]
    }]
};

(* ===================== Generalizations & Extensions (1) ===================== *)

genExtContent = {
    exampleSubsection["Vertex-operator Paulis (VOPs)", {
        Cell[TextData[{
            "Every ", inlineCode[gsBtn[]],
            " also stores a ", inlineCode["\"VOPs\""],
            " list of integer indices into the 24-element single-qubit Clifford ",
            "group. Currently only the identity VOP (index ", inlineCode["0"], ") is ",
            "exercised; tracked VOPs are on the roadmap and would extend ",
            inlineCode[gsBtn[]],
            " to every local-Clifford representative of a stabilizer state."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[CycleGraph[4]]["VOPs"], 1]
    }]
};

(* ===================== Options (literal ladder of property accessors) ===================== *)

optionsContent = {
    exampleSubsection["\"Graph\"", {
        Cell[TextData[{
            inlineCode["\"Graph\""],
            " returns the underlying graph:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[StarGraph[5]]["Graph"], 1]
    }],

    exampleSubsection["\"VertexCount\" and \"EdgeCount\"", {
        Cell[TextData[{
            inlineCode["\"VertexCount\""], " and ",
            inlineCode["\"EdgeCount\""],
            " return the size of the underlying graph:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs},
                gs = GraphState[CycleGraph[5]];
                {gs["VertexCount"], gs["EdgeCount"]}
            ], 1]
    }],

    exampleSubsection["\"AdjacencyMatrix\"", {
        Cell[TextData[{
            inlineCode["\"AdjacencyMatrix\""],
            " returns the ", inlineMath[RowBox[{TI["n"], "\[Times]", TI["n"]}]],
            " adjacency matrix of the underlying graph:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[PathGraph[{1, 2, 3, 4}]]["AdjacencyMatrix"] // MatrixForm, 1]
    }],

    exampleSubsection["\"Stabilizers\"", {
        Cell[TextData[{
            inlineCode["\"Stabilizers\""],
            " returns the Pauli-string stabilizer generators of the graph state:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[CycleGraph[3]]["Stabilizers"], 1]
    }],

    exampleSubsection["\"Vertices\" and \"Edges\"", {
        Cell[TextData[{
            inlineCode["\"Vertices\""], " and ",
            inlineCode["\"Edges\""],
            " return the vertex and edge lists of the underlying graph:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs},
                gs = GraphState[CycleGraph[3]];
                {gs["Vertices"], gs["Edges"]}
            ], 1]
    }],

    exampleSubsection["\"Qubits\"", {
        Cell[TextData[{
            inlineCode["\"Qubits\""],
            " is an alias for ", inlineCode["\"VertexCount\""], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[CycleGraph[5]]["Qubits"], 1]
    }],

    exampleSubsection["\"VOPs\"", {
        Cell[TextData[{
            inlineCode["\"VOPs\""],
            " returns the (currently all-identity) vertex-operator-Pauli indices:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[CycleGraph[3]]["VOPs"], 1]
    }],

    exampleSubsection["\"PauliStabilizer\"", {
        Cell[TextData[{
            inlineCode["\"PauliStabilizer\""],
            " converts the graph state into a ",
            inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
            ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[StarGraph[3]]["PauliStabilizer"], 1]
    }]
};

(* ===================== Applications (2) ===================== *)

applicationsContent = {
    exampleSubsection["Cluster states for MBQC", {
        Cell[TextData[{
            "A 1D cluster state of length 5 — the workhorse of measurement-based ",
            "quantum computation (Briegel & Raussendorf, arXiv:quant-ph/0004051):"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs},
                gs = GraphState[PathGraph[{1, 2, 3, 4, 5}]];
                {gs["EdgeCount"], gs["Stabilizers"]}
            ], 1]
    }],

    exampleSubsection["GHZ state in graph form", {
        Cell[TextData[{
            "The ", inlineMath[TI["n"]],
            "-qubit GHZ state is local-Clifford-equivalent to the star graph state ",
            "with one centre vertex and ", inlineMath[RowBox[{TI["n"], "-", "1"}]],
            " leaves:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs},
                gs = GraphState[StarGraph[4]];
                {gs["EdgeCount"], gs["Stabilizers"]}
            ], 1]
    }]
};

(* ===================== Properties & Relations (4) ===================== *)

propsRelContent = {
    exampleSubsection["Round-trip Graph -> GraphState -> PauliStabilizer", {
        Cell["Building a GraphState then converting to PauliStabilizer recovers the canonical graph-form stabilizer rows:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs, ps},
                gs = GraphState[CycleGraph[3]];
                ps = gs["PauliStabilizer"];
                ps["Stabilizers"]
            ], 1]
    }],

    exampleSubsection["LocalComplement maps GraphState to GraphState", {
        Cell[TextData[{
            inlineCode[paclLink["LocalComplement", "Wolfram/QuantumFramework/ref/LocalComplement"]],
            " at any vertex returns another ", inlineCode[gsBtn[]],
            " in the same local-Clifford equivalence class:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs, lc},
                gs = GraphState[StarGraph[4]];
                lc = LocalComplement[gs, 1];
                {gs["EdgeCount"], lc["EdgeCount"]}
            ], 1]
    }],

    exampleSubsection["Edge count equals stabilizer Z-weight sum / 2", {
        Cell[TextData[{
            "Every edge contributes one ", inlineMath["Z"],
            " to each of its two endpoints' stabilizer rows; thus the total ",
            inlineMath["Z"], " count equals twice the edge count:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{gs, zCount},
                gs = GraphState[CycleGraph[5]];
                zCount = Total[StringCount[#, "Z"] & /@ gs["Stabilizers"]];
                {gs["EdgeCount"], zCount, zCount == 2 gs["EdgeCount"]}
            ], 1]
    }],

    exampleSubsection["Empty graph -> product state", {
        Cell[TextData[{
            "On a graph with no edges, the graph state is the all-",
            inlineMath["+"], " product state with stabilizers ",
            inlineMath[RowBox[{sub["X", "1"], ",", sub["X", "2"], ",", "\[Ellipsis]"}]],
            ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            GraphState[Graph[{1, 2, 3}, {}]]["Stabilizers"], 1]
    }]
};

(* ===================== Possible Issues (3) ===================== *)

possIssuesContent = {
    exampleSubsection["::nongraph for non-graph-form stabilizers", {
        Cell[TextData[{
            "A ", inlineCode["PauliStabilizer"],
            " whose rows do not have exactly one ", inlineMath["X"],
            " on the diagonal and ", inlineMath["Z"], "/", inlineMath["I"],
            " elsewhere triggers ", inlineCode["GraphState::nongraph"],
            " and returns ", inlineCode["$Failed"], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Quiet @ Check[GraphState[PauliStabilizer[{"YY", "ZZ"}]], "$Failed"], 1]
    }],

    exampleSubsection["VOPs not yet tracked", {
        Cell[TextData[{
            "Operations that should rotate a vertex's single-qubit Clifford (notably ",
            inlineCode[paclLink["LocalComplement", "Wolfram/QuantumFramework/ref/LocalComplement"]],
            ") currently leave the ", inlineCode["\"VOPs\""],
            " field untouched. Use the underlying ", inlineCode["PauliStabilizer"],
            " for exact unitary evolution; the ", inlineCode[gsBtn[]],
            " carries only the combinatorial graph information until VOPs ",
            "are tracked."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]]
    }],

    exampleSubsection["Stabilizer row count must equal vertex count", {
        Cell[TextData[{
            "The conversion ", inlineCode[gsCall[TI["ps"]]],
            " requires ", inlineMath[TI["ps"]],
            " to have exactly one stabilizer row per qubit; otherwise the ",
            "graph-form check fails."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]]
    }]
};

(* ===================== Neat Examples (1) ===================== *)

neatContent = {
    exampleSubsection["Cluster vs Star vs Complete on 5 vertices", {
        Cell["Compare the edge counts of three classic 5-vertex graph states:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{cluster, star, complete},
                cluster = GraphState[PathGraph[{1, 2, 3, 4, 5}]];
                star = GraphState[StarGraph[5]];
                complete = GraphState[CompleteGraph[5]];
                <|
                    "Cluster" -> cluster["EdgeCount"],
                    "Star" -> star["EdgeCount"],
                    "Complete" -> complete["EdgeCount"]
                |>
            ], 1]
    }]
};

(* ===================== See Also ===================== *)

seeAlsoLinks = {
    "LocalComplement", "PauliStabilizer", "StabilizerStateQ",
    "StabilizerFrame", "CliffordChannel",
    "QuantumState", "QuantumCircuitOperator"};

systemSeeAlso = {"Graph", "AdjacencyMatrix", "PauliMatrix"};

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
    "graph state", "Anders-Briegel", "VOP", "cluster state",
    "MBQC", "local complementation", "stabilizer formalism",
    "Hein-Eisert-Briegel", "graph form"};

(* ===================== Assemble the Notebook ===================== *)

notebook = Notebook[{

    Cell[CellGroupData[{
        Cell["GraphState", "ObjectName",
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
            Cell["Wolfram/QuantumFramework/ref/GraphState", "Categorization",
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

(* Sanity check: did Export succeed? Fail loudly if not. *)
exportResult = Export[$docFile, notebook, "Notebook"];
If[exportResult === $Failed,
    Print["ERROR: Export returned $Failed"];
    Exit[1]
];
Print["Wrote notebook to: ", $docFile];
Print["Done. File size: ", FileByteCount[$docFile], " bytes"];
