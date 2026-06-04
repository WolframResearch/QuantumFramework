#!/usr/bin/env wolframscript
(* ============================================================================
   build_stabilizerstateq_doc.wl
   ----------------------------------------------------------------------------
   Generates QuantumFramework/Documentation/English/ReferencePages/Symbols/
       StabilizerStateQ.nb
   Pattern B (function-like predicate). Outputs are Boolean values.

   Run:  wolframscript -file OngoingProjects/Stabilizer/build_stabilizerstateq_doc.wl
   ============================================================================ *)

$repoRoot  = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$pacletDir = FileNameJoin[{$repoRoot, "QuantumFramework"}];
$docFile   = FileNameJoin[{$pacletDir, "Documentation", "English",
    "ReferencePages", "Symbols", "StabilizerStateQ.nb"}];

Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

SeedRandom[20260511];

(* ===================== Helpers ===================== *)

cid[]  := RandomInteger[{10^6, 10^9 - 1}];
uuid[] := CreateUUID[];

TI[s_String] := StyleBox[s, "TI"];

ssqBtn[] := ButtonBox["StabilizerStateQ", BaseStyle -> "Link",
    ButtonData -> "paclet:Wolfram/QuantumFramework/ref/StabilizerStateQ"];

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

SetAttributes[example, HoldRest];
example[caption_String, expr_] := Flatten @ {
    Cell[caption, "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    ioBlock[expr, 1]
};

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

(* ===================== Usage Cell (1 pattern) ===================== *)

ssqCall[argBoxes___] := RowBox[{ssqBtn[], "[", argBoxes, "]"}];

usageContent = Flatten[{
    modInfo[],
    inlineCode[ssqCall[TI["expr"]]],
    "\[LineSeparator]returns ", inlineCode["True"], " if ", inlineMath[TI["expr"]],
    " is a structurally valid stabilizer-state representation \
(a ", inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
    " or a single-component ", inlineCode[paclLink["StabilizerFrame", "Wolfram/QuantumFramework/ref/StabilizerFrame"]],
    "), and ", inlineCode["False"], " otherwise."
}];

usageCell = Cell[TextData[usageContent], "Usage",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Notes Cells ===================== *)

notesCells = {
    Cell[TextData[{
        inlineCode[ssqBtn[]], " is a structural predicate: it inspects the head and ",
        "tableau shape of its argument and does not perform any state-vector ",
        "tomography or numerical comparison."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode[RowBox[{ssqBtn[], "[", TI["ps"], "]"}]], " returns ", inlineCode["True"],
        " when ", inlineMath[TI["ps"]], " is a ",
        inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
        " with the required ", inlineCode["\"Signs\""], " and ", inlineCode["\"Tableau\""],
        " keys whose dimensions match (length-",
        inlineMath[TI["m"]], " sign vector and a ",
        inlineMath[RowBox[{"2", "\[Times]", TI["n"], "\[Times]", TI["m"]}]],
        " tableau)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode[RowBox[{ssqBtn[], "[", TI["sf"], "]"}]], " returns ", inlineCode["True"],
        " when ", inlineMath[TI["sf"]], " is a ",
        inlineCode[paclLink["StabilizerFrame", "Wolfram/QuantumFramework/ref/StabilizerFrame"]],
        " whose ", inlineCode["\"Length\""], " is exactly one (a single-component frame ",
        "is equivalent to a stabilizer state up to a scalar)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode[ssqBtn[]], " returns ", inlineCode["False"], " on every other input, ",
        "including ", inlineCode[paclLink["QuantumState", "Wolfram/QuantumFramework/ref/QuantumState"]],
        " expressions. Detecting whether an arbitrary ", inlineCode["QuantumState"],
        " is a stabilizer state requires ", inlineMath[sup["4", TI["n"]]],
        " Pauli tomography and is exposed separately via ",
        inlineCode[RowBox[{paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"], "[", TI["qs"], "]"}]],
        ", which returns a ", inlineCode["PauliStabilizer"], " on success or ",
        inlineCode["$Failed"], " otherwise."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "A multi-component ", inlineCode["StabilizerFrame"],
        " (such as the result of applying a non-Clifford ",
        inlineCode["T"], " gate) is ", inlineCode["False"], " under ", inlineCode[ssqBtn[]],
        " because it represents a superposition of stabilizer states rather than ",
        "a stabilizer state itself."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Stabilizer states form a non-generic measure-zero subset of the ",
        inlineMath[sup["2", TI["n"]]], "-dimensional Hilbert space; their cardinality scales ",
        "as ", inlineMath[RowBox[{sup["2", TI["n"]], "\[CenterDot]",
            RowBox[{UnderoverscriptBox["\[Product]", RowBox[{TI["k"], "=", "1"}], TI["n"]],
                RowBox[{"(", sup["2", TI["k"]], "+", "1", ")"}]}]}]],
        " (Aaronson-Gottesman 2004). The predicate is therefore a useful guard ",
        "for fast-path dispatch in functions that accept either tableau or ",
        "state-vector inputs."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]]
};

(* ===================== Primary Examples (5) ===================== *)

primaryExamplesContent = Flatten[{
    (* Init cell renders separately; here we build I/O pairs *)
    Cell["Test the empty (1-qubit) stabilizer register:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[StabilizerStateQ[PauliStabilizer[]], 1],
    delim[],

    Cell["Test a 3-qubit |000\[RightAngleBracket] register:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[StabilizerStateQ[PauliStabilizer[3]], 2],
    delim[],

    Cell["Test a Bell-state stabilizer constructed from Pauli strings:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[StabilizerStateQ[PauliStabilizer[{"XX", "ZZ"}]], 3],
    delim[],

    Cell["A QuantumState is not implicitly tomographed; the predicate returns False:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[StabilizerStateQ[QuantumState[{1, 0}]], 4],
    delim[],

    Cell["A single-component StabilizerFrame qualifies; a non-state value does not:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        {StabilizerStateQ[StabilizerFrame[PauliStabilizer[2]]],
         StabilizerStateQ[42]}, 5]
}];

(* ===================== Scope (3 subsections) ===================== *)

scopeContent = {
    exampleSubsection["PauliStabilizer inputs", {
        Cell["Every well-formed PauliStabilizer is recognized regardless of how it was constructed:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerStateQ /@ {
                PauliStabilizer[],
                PauliStabilizer[5],
                PauliStabilizer[{"XX", "ZZ"}],
                PauliStabilizer["5QubitCode"],
                PauliStabilizer["Random", 3]
            }, 1]
    }],

    exampleSubsection["StabilizerFrame inputs", {
        Cell["A single-component frame is recognized as a stabilizer state:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerStateQ[StabilizerFrame[PauliStabilizer[3]]], 1],
        delim[],

        Cell["A multi-component frame (e.g. the result of a T gate) is rejected:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps},
                ps = PauliStabilizer[2];
                StabilizerStateQ[ps["T", 1]]
            ], 1]
    }],

    exampleSubsection["Rejection cases", {
        Cell["QuantumState, plain matrices, and non-state expressions return False:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerStateQ /@ {
                QuantumState[{1, 0}],
                {1, 0, 0, 0},
                IdentityMatrix[2],
                "PauliStabilizer",
                Null
            }, 1]
    }]
};

(* ===================== Generalizations & Extensions (1) ===================== *)

genExtContent = {
    exampleSubsection["Tomographic detection from a QuantumState", {
        Cell[TextData[{
            "Tomographic detection from a generic ",
            inlineCode[paclLink["QuantumState", "Wolfram/QuantumFramework/ref/QuantumState"]],
            " is intentionally not done inside ",
            inlineCode[ssqBtn[]], ". Use the ",
            inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
            " constructor on a QuantumState and check the result:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{qs, ps},
                qs = QuantumState["GHZ"];
                ps = Quiet @ PauliStabilizer[qs];
                {Head[ps], StabilizerStateQ[ps]}
            ], 1]
    }]
};

(* ===================== Options (empty) ===================== *)

optionsContent = {
    exampleSubsection["No options", {
        Cell[TextData[{
            inlineCode[ssqBtn[]], " takes no options."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]]
    }]
};

(* ===================== Applications (2) ===================== *)

applicationsContent = {
    exampleSubsection["Filter a mixed list of states", {
        Cell["Select stabilizer-state entries from a heterogeneous list:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Select[
                {PauliStabilizer[2], QuantumState[{1, 0}],
                 PauliStabilizer[{"XY", "ZZ"}], 42,
                 StabilizerFrame[PauliStabilizer[1]]},
                StabilizerStateQ
            ], 1]
    }],

    exampleSubsection["Verify a round-trip QuantumState -> PauliStabilizer", {
        Cell["Confirm that PauliStabilizer[qs] succeeds by checking the result:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{qs, ps},
                qs = QuantumState["Bell"];
                ps = PauliStabilizer[qs];
                StabilizerStateQ[ps]
            ], 1]
    }]
};

(* ===================== Properties & Relations (4) ===================== *)

propsRelContent = {
    exampleSubsection["Single-component frame equivalence", {
        Cell[TextData[{
            inlineCode[ssqBtn[]], " agrees on a ",
            inlineCode["PauliStabilizer"], " and its single-component-frame embedding:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps},
                ps = PauliStabilizer[{"XX", "ZZ"}];
                StabilizerStateQ[ps] === StabilizerStateQ[StabilizerFrame[ps]]
            ], 1]
    }],

    exampleSubsection["T-gate exits the predicate", {
        Cell["Applying a T gate promotes a PauliStabilizer to a multi-component frame, so the predicate flips to False:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps, after},
                ps = PauliStabilizer[2];
                after = ps["T", 1];
                {StabilizerStateQ[ps], StabilizerStateQ[after]}
            ], 1]
    }],

    exampleSubsection["Clifford gates preserve the predicate", {
        Cell["All Clifford gates on a PauliStabilizer return a PauliStabilizer, so the predicate is preserved:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps, after},
                ps = PauliStabilizer[3];
                after = ps["H", 1]["CNOT", 1, 2]["S", 3];
                StabilizerStateQ[after]
            ], 1]
    }],

    exampleSubsection["Non-state objects return False", {
        Cell["Channels, operators, and circuits are not stabilizer states even if they encode Clifford behaviour:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerStateQ /@ {
                CliffordChannel["Identity", 2],
                QuantumOperator["X"],
                QuantumCircuitOperator["GHZ"]
            }, 1]
    }]
};

(* ===================== Possible Issues (3) ===================== *)

possIssuesContent = {
    exampleSubsection["QuantumState is never True without explicit conversion", {
        Cell[TextData[{
            "Even known stabilizer states return ", inlineCode["False"], " when wrapped ",
            "in a ", inlineCode["QuantumState"], " (the predicate intentionally avoids ",
            inlineMath[sup["4", TI["n"]]], " tomography):"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerStateQ /@ {
                QuantumState["Bell"],
                QuantumState["GHZ"],
                QuantumState[{1, 0, 0, 0} / 1]
            }, 1]
    }],

    exampleSubsection["Multi-component StabilizerFrame returns False", {
        Cell["The predicate is restricted to length-1 frames; longer frames represent superpositions, not single stabilizer states:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f},
                f = StabilizerFrame[{
                    {1, PauliStabilizer[1]},
                    {1, PauliStabilizer[{"Z"}]}
                }];
                {f["Length"], StabilizerStateQ[f]}
            ], 1]
    }],

    exampleSubsection["Malformed Association arguments", {
        Cell[TextData[{
            "A ", inlineCode["PauliStabilizer"], " whose internal Association is missing ",
            "the required keys is structurally invalid and the predicate is ",
            inlineCode["False"], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerStateQ[PauliStabilizer[<|"Junk" -> 0|>]], 1]
    }]
};

(* ===================== Neat Examples (1) ===================== *)

neatContent = {
    exampleSubsection["Counting stabilizer states", {
        Cell[TextData[{
            "The number of ", inlineMath[TI["n"]],
            "-qubit pure stabilizer states is ",
            inlineMath[RowBox[{sup["2", TI["n"]], "\[CenterDot]",
                RowBox[{UnderoverscriptBox["\[Product]", RowBox[{TI["k"], "=", "1"}], TI["n"]],
                    RowBox[{"(", sup["2", TI["k"]], "+", "1", ")"}]}]}]],
            ". Compute the first few terms:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Table[2^n * Product[2^k + 1, {k, 1, n}], {n, 1, 6}], 1]
    }]
};

(* ===================== See Also ===================== *)

seeAlsoLinks = {
    "PauliStabilizer", "StabilizerFrame", "CliffordChannel",
    "GraphState", "LocalComplement",
    "QuantumState"};

systemSeeAlso = {"MatrixQ", "BooleanQ"};

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
    "stabilizer state", "predicate", "Pauli", "tableau",
    "Clifford", "Aaronson-Gottesman", "structural check",
    "type guard"};

(* ===================== Assemble the Notebook ===================== *)

notebook = Notebook[{

    Cell[CellGroupData[{
        Cell["StabilizerStateQ", "ObjectName",
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
            Cell["Wolfram/QuantumFramework/ref/StabilizerStateQ", "Categorization",
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
