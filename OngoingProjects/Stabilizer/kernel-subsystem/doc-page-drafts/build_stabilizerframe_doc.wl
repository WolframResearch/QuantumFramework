#!/usr/bin/env wolframscript
(* ============================================================================
   build_stabilizerframe_doc.wl
   ----------------------------------------------------------------------------
   Generates QuantumFramework/Documentation/English/ReferencePages/Symbols/
       StabilizerFrame.nb
   Pattern A (object-head with custom display, ScriptCapitalF summary box).

   Run:  wolframscript -file OngoingProjects/Stabilizer/build_stabilizerframe_doc.wl
   ============================================================================ *)

$repoRoot  = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$pacletDir = FileNameJoin[{$repoRoot, "QuantumFramework"}];
$docFile   = FileNameJoin[{$pacletDir, "Documentation", "English",
    "ReferencePages", "Symbols", "StabilizerFrame.nb"}];

Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

SeedRandom[20260511];

(* ===================== Helpers ===================== *)

cid[]  := RandomInteger[{10^6, 10^9 - 1}];
uuid[] := CreateUUID[];

TI[s_String] := StyleBox[s, "TI"];

sfBtn[] := ButtonBox["StabilizerFrame", BaseStyle -> "Link",
    ButtonData -> "paclet:Wolfram/QuantumFramework/ref/StabilizerFrame"];

paclLink[sym_String, path_String] := ButtonBox[sym, BaseStyle -> "Link",
    ButtonData -> "paclet:" <> path];

modInfo[] := Cell["   ", "ModInfo", ExpressionUUID -> uuid[]];

inlineCode[boxes_] := Cell[BoxData[boxes], "InlineFormula", ExpressionUUID -> uuid[]];
inlineMath[boxes_] := Cell[
    BoxData[FormBox[boxes, TraditionalForm]],
    "InlineFormula", ExpressionUUID -> uuid[]];
inlineKet[content_] := Cell[
    BoxData[FormBox[TemplateBox[{content}, "Ket"], TraditionalForm]],
    FormatType -> TraditionalForm,
    ExpressionUUID -> uuid[]];

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

(* ===================== Usage Cell (4 patterns) ===================== *)

sfCall[argBoxes___] := RowBox[{sfBtn[], "[", argBoxes, "]"}];

usageContent = Flatten[{
    modInfo[],
    inlineCode[sfCall[TI["ps"]]],
    "\[LineSeparator]wraps a ",
    inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
    " ", inlineMath[TI["ps"]], " as a single-component stabilizer frame with coefficient ",
    inlineCode["1"], ".\n",

    modInfo[],
    inlineCode[sfCall[RowBox[{"{", RowBox[{
        RowBox[{"{", RowBox[{sub["c", "1"], ",", sub["ps", "1"]}], "}"}],
        ",",
        RowBox[{"{", RowBox[{sub["c", "2"], ",", sub["ps", "2"]}], "}"}],
        ",", "\[Ellipsis]"
    }], "}"}]]],
    "\[LineSeparator]builds a frame representing the superposition ",
    inlineMath[RowBox[{
        UnderscriptBox["\[Sum]", TI["i"]], " ",
        sub["c", "i"], " ", TemplateBox[{sub["s", "i"]}, "Ket"]
    }]],
    " of stabilizer states with the given (possibly symbolic) coefficients.\n",

    modInfo[],
    inlineCode[RowBox[{TI["f"], "[", RowBox[{"\"", TI["gate"], "\"", ",", TI["args"]}], "]"}]],
    "\[LineSeparator]distributes the named gate (Clifford or non-Clifford) over each ",
    "component of the frame ", inlineMath[TI["f"]],
    ".\n",

    modInfo[],
    inlineCode[RowBox[{TI["f"], "[", RowBox[{"\"InnerProduct\"", ",", TI["other"]}], "]"}]],
    "\[LineSeparator]computes the inner product of the frame ", inlineMath[TI["f"]],
    " with another stabilizer state or frame."
}];

usageCell = Cell[TextData[usageContent], "Usage",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Notes Cells ===================== *)

notesCells = {
    Cell[TextData[{
        "A ", inlineCode[sfBtn[]],
        " represents a complex linear combination ",
        inlineMath[RowBox[{
            UnderscriptBox["\[Sum]", TI["i"]],
            " ", sub["c", "i"], " ", TemplateBox[{sub["s", "i"]}, "Ket"]}]],
        " of stabilizer states, the natural object for circuits that contain ",
        "a small number of non-Clifford gates, magic-state distillation, and ",
        "stabilizer-rank simulation. Reference: Garc\[IAcute]a-Mart\[IAcute]n & Markov ",
        "(arXiv:1712.03554)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The internal representation is the Association ",
        inlineCode[RowBox[{"<|", RowBox[{"\"Components\"", "\[Rule]",
            RowBox[{"{", RowBox[{
                RowBox[{"{", RowBox[{sub["c", "1"], ",", sub["ps", "1"]}], "}"}],
                ",", "\[Ellipsis]"}], "}"}]}], "|>"}]],
        " where each ", inlineCode["ps_i"], " is a ",
        inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
        " and the ", inlineCode["c_i"], " can be exact, numeric, or symbolic ",
        "complex coefficients."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The summary box reports the number of components, the qubit count, ",
        "and the first coefficient. Inspect any property with the dispatch syntax ",
        inlineCode[RowBox[{TI["f"], "[", "\"", TI["property"], "\"", "]"}]],
        "; the list of recognised property names is given by ",
        inlineCode[RowBox[{TI["f"], "[", "\"Properties\"", "]"}]], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Clifford gates (",
        inlineCode["\"H\""], ", ", inlineCode["\"S\""], ", ",
        inlineCode["\"X\""], ", ", inlineCode["\"CNOT\""], ", \[Ellipsis]) are ",
        "distributed over the components without changing the frame length, so a ",
        "Clifford circuit acting on a single-component frame remains a single-component ",
        "frame and is equivalent in cost to the underlying tableau update."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The non-Clifford ", inlineCode[RowBox[{"\"P\"", "[", "\[Theta]", "]"}]],
        " gate is decomposed as ",
        inlineMath[RowBox[{
            RowBox[{"P", "(", "\[Theta]", ")"}], " ", TemplateBox[{"s"}, "Ket"],
            "=", " ",
            FractionBox[RowBox[{"1", "+", sup["e", RowBox[{"i", "\[Theta]", "/", "2"}]]}], "2"],
            " ", TemplateBox[{"s"}, "Ket"], "+", " ",
            FractionBox[RowBox[{"1", "-", sup["e", RowBox[{"i", "\[Theta]", "/", "2"}]]}], "2"],
            " ", sub["Z", TI["q"]], " ", TemplateBox[{"s"}, "Ket"]
        }]],
        ", so each P-gate doubles the component count. The ",
        inlineCode["\"T\""], " gate is the alias ",
        inlineCode[RowBox[{"\"P\"", "[", RowBox[{"\[Pi]", "/", "2"}], "]"}]],
        " and ", inlineCode[SuperscriptBox["\"T\"", "\[Dagger]"]],
        " is the alias ",
        inlineCode[RowBox[{"\"P\"", "[", RowBox[{"-", "\[Pi]", "/", "2"}], "]"}]], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode[RowBox[{TI["f"], "[", "\"State\"", "]"}]], " or ",
        inlineCode[RowBox[{TI["f"], "[", "\"StateVector\"", "]"}]],
        " materialises the full quantum state by summing the dense state vectors of ",
        "the components. This is only practical for small qubit counts and frame lengths."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode[RowBox[{TI["f"], "[", "\"InnerProduct\"", ",", TI["other"], "]"}]],
        " returns the inner product of ", inlineMath[TI["f"]], " with ",
        inlineMath[TI["other"]],
        " (a frame, ", inlineCode["PauliStabilizer"],
        ", or stabilizer-compatible state) by reducing to the underlying tableau-level ",
        "inner products via Garc\[IAcute]a-Markov's recursion."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode["Plus"], ", ", inlineCode["Times"], ", and ",
        inlineCode["Equal"],
        " are defined on ", inlineCode[sfBtn[]],
        " via ", inlineCode["UpValues"],
        ": addition concatenates component lists, scalar multiplication rescales ",
        "every coefficient, and equality compares the (ordered) component lists for ",
        "structural sameness."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]]
};

(* ===================== Primary Examples (5) — Pattern A ===================== *)

primaryExamplesContent = Flatten[{
    Cell["A StabilizerFrame built from a single PauliStabilizer:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[StabilizerFrame[PauliStabilizer[2]], 1],
    delim[],

    Cell["An explicit two-component frame with equal coefficients:", "ExampleText",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        StabilizerFrame[{
            {1/Sqrt[2], PauliStabilizer[{"Z"}]},
            {1/Sqrt[2], PauliStabilizer[{"-Z"}]}
        }], 2],
    delim[],

    Cell["Applying a T gate to a single-component frame doubles the frame size:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        StabilizerFrame[PauliStabilizer[1]]["T", 1], 3],
    delim[],

    Cell["Distributing a Hadamard gate over each component preserves the frame length:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        StabilizerFrame[PauliStabilizer[2]]["H", 1], 4],
    delim[],

    Cell["Inspecting structure with the accessor syntax:",
        "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ ioBlock[
        Module[{f},
            f = StabilizerFrame[{
                {1/Sqrt[2], PauliStabilizer[1]},
                {1/Sqrt[2], PauliStabilizer[{"Z"}]}
            }];
            {f["Length"], f["Qubits"], f["Coefficients"]}
        ], 5]
}];

(* ===================== Scope (4 subsections) ===================== *)

scopeContent = {
    exampleSubsection["Constructors", {
        Cell["From a PauliStabilizer (coefficient 1):", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[StabilizerFrame[PauliStabilizer[3]], 1],
        delim[],

        Cell["From an explicit list of {coefficient, PauliStabilizer} pairs:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[{
                {Cos[\[Theta]/2], PauliStabilizer[1]},
                {I Sin[\[Theta]/2], PauliStabilizer[{"X"}]}
            }], 1]
    }],

    exampleSubsection["Properties", {
        Cell[TextData[{
            "Every recognised property is enumerated by ",
            inlineCode[RowBox[{TI["f"], "[", "\"Properties\"", "]"}]], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[2]]["Properties"], 1],
        delim[],

        Cell["Component list, length, and qubit count:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f},
                f = StabilizerFrame[{
                    {1, PauliStabilizer[2]},
                    {1, PauliStabilizer[{"XX", "ZZ"}]}
                }];
                {f["Length"], f["Qubits"], Length[f["Components"]]}
            ], 1]
    }],

    exampleSubsection["Materialization to a QuantumState", {
        Cell["Sum the dense state vectors of every component:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[{
                {1/Sqrt[2], PauliStabilizer[1]},
                {1/Sqrt[2], PauliStabilizer[{"Z"}]}
            }]["State"], 1]
    }],

    exampleSubsection["Inner product", {
        Cell["Inner product of two frames:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f1, f2},
                f1 = StabilizerFrame[PauliStabilizer[1]];
                f2 = StabilizerFrame[PauliStabilizer[{"Z"}]];
                f1["InnerProduct", f2]
            ], 1],
        delim[],

        Cell["Frames are orthogonal when their underlying stabilizer states are:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f1, f2},
                f1 = StabilizerFrame[PauliStabilizer[{"Z"}]];
                f2 = StabilizerFrame[PauliStabilizer[{"-Z"}]];
                f1["InnerProduct", f2]
            ], 1]
    }]
};

(* ===================== Generalizations & Extensions (1) ===================== *)

genExtContent = {
    exampleSubsection["Non-Clifford P-gate doubles the frame", {
        Cell[TextData[{
            "Each ", inlineCode[RowBox[{"\"P\"", "[", "\[Theta]", "]"}]],
            " gate doubles the component count by emitting both the unmodified ",
            "and the ", inlineMath[sub["Z", TI["q"]]], "-flipped branch:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f0, f1, f2, f3},
                f0 = StabilizerFrame[PauliStabilizer[1]];
                f1 = f0["T", 1];
                f2 = f1["T", 1];
                f3 = f2["T", 1];
                {f0["Length"], f1["Length"], f2["Length"], f3["Length"]}
            ], 1]
    }]
};

(* ===================== Options (Method dispatch literal ladder, 7 subsections) ===================== *)

optionsContent = {
    exampleSubsection["Clifford single-qubit gates", {
        Cell[TextData[{
            "Single-qubit Clifford gates such as ", inlineCode["\"H\""],
            ", ", inlineCode["\"S\""], ", ", inlineCode["\"X\""],
            ", ", inlineCode["\"Y\""], ", ", inlineCode["\"Z\""],
            ", and ", inlineCode["\"V\""],
            " distribute over each component:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[2]]["H", 1], 1]
    }],

    exampleSubsection["Clifford two-qubit gates", {
        Cell[TextData[{
            inlineCode["\"CNOT\""], ", ", inlineCode["\"CZ\""],
            ", and ", inlineCode["\"SWAP\""],
            " act analogously and are routed through the underlying ",
            inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
            " gate dispatch:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[2]]["CNOT", 1, 2], 1]
    }],

    exampleSubsection["SuperDagger gates", {
        Cell[TextData[{
            "Daggered gates use the ",
            inlineCode[RowBox[{SuperscriptBox["\"S\"", "\[Dagger]"]}]],
            " / ",
            inlineCode[RowBox[{SuperscriptBox["\"V\"", "\[Dagger]"]}]],
            " syntax:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[1]][SuperDagger["S"], 1], 1]
    }],

    exampleSubsection[
        "\"P\"[\[Theta]] arbitrary-angle phase gate", {
        Cell[TextData[{
            inlineCode[RowBox[{"\"P\"", "[", "\[Theta]", "]"}]],
            " applies a phase-shift gate and doubles the frame size; with symbolic ",
            inlineMath["\[Theta]"], " the coefficients become trigonometric expressions:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[1]]["P"[\[Theta]], 1], 1]
    }],

    exampleSubsection["\"T\" gate", {
        Cell[TextData[{
            inlineCode["\"T\""], " is the canonical non-Clifford gate, aliased to ",
            inlineCode[RowBox[{"\"P\"", "[", RowBox[{"\[Pi]", "/", "2"}], "]"}]],
            ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[1]]["T", 1], 1]
    }],

    exampleSubsection[
        "SuperDagger[\"T\"] (T\[Dagger]) gate", {
        Cell[TextData[{
            inlineCode[SuperscriptBox["\"T\"", "\[Dagger]"]],
            " is the inverse of ", inlineCode["\"T\""], ", aliased to ",
            inlineCode[RowBox[{"\"P\"", "[", RowBox[{"-", "\[Pi]", "/", "2"}], "]"}]],
            ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[1]][SuperDagger["T"], 1], 1]
    }],

    exampleSubsection["\"InnerProduct\" method", {
        Cell["Inner product of a frame with another frame:", "ExampleText",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f1, f2},
                f1 = StabilizerFrame[PauliStabilizer[{"Z"}]];
                f2 = StabilizerFrame[PauliStabilizer[1]]["T", 1];
                f1["InnerProduct", f2]
            ], 1]
    }],

    exampleSubsection[
        "op -> order rewrite", {
        Cell[TextData[{
            "Gate calls also accept the ",
            inlineCode[RowBox[{TI["op"], "\[Rule]", TI["order"]}]],
            " rewrite form used by ",
            inlineCode[paclLink["QuantumOperator", "Wolfram/QuantumFramework/ref/QuantumOperator"]],
            ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[2]]["CNOT" -> {1, 2}], 1]
    }]
};

(* ===================== Applications (2) ===================== *)

applicationsContent = {
    exampleSubsection["Symbolic T-gate evolution on a single qubit", {
        Cell[TextData[{
            "Track the symbolic coefficient growth after a sequence of ",
            inlineCode["\"T\""], " and ", inlineCode["\"H\""], " gates:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f},
                f = StabilizerFrame[PauliStabilizer[1]];
                f = f["T", 1]["H", 1]["T", 1];
                f["Coefficients"] // Simplify
            ], 1]
    }],

    exampleSubsection["Frame growth fingerprint", {
        Cell["Track the doubling pattern of frame length under successive T gates on disjoint qubits:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f, history},
                f = StabilizerFrame[PauliStabilizer[3]];
                history = NestList[Function[g, g["T", RandomInteger[{1, 3}]]], f, 4];
                #["Length"] & /@ history
            ], 1]
    }]
};

(* ===================== Properties & Relations (4) ===================== *)

propsRelContent = {
    exampleSubsection["Single-component frame equals its PauliStabilizer", {
        Cell["A length-1 frame and its underlying PauliStabilizer have the same state:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps, f},
                ps = PauliStabilizer[2];
                f = StabilizerFrame[ps];
                f["State"] == ps["State"]
            ], 1]
    }],

    exampleSubsection["Addition concatenates component lists", {
        Cell[TextData[{
            "The ", inlineCode["Plus"], " upvalue joins the component lists:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f1, f2, sum},
                f1 = StabilizerFrame[PauliStabilizer[1]];
                f2 = StabilizerFrame[PauliStabilizer[{"Z"}]];
                sum = f1 + f2;
                {f1["Length"], f2["Length"], sum["Length"]}
            ], 1]
    }],

    exampleSubsection["Scalar multiplication rescales every coefficient", {
        Cell[TextData[{
            "The ", inlineCode["Times"], " upvalue acts element-wise on coefficients:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f},
                f = StabilizerFrame[{
                    {1, PauliStabilizer[1]},
                    {1, PauliStabilizer[{"Z"}]}
                }];
                (1/Sqrt[2] f)["Coefficients"]
            ], 1]
    }],

    exampleSubsection["Frame length under non-Clifford circuits", {
        Cell[TextData[{
            "Frame length doubles per ", inlineCode["\"T\""],
            " gate applied to distinct components; Clifford gates leave the length unchanged:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f},
                f = StabilizerFrame[PauliStabilizer[2]];
                {
                    f["Length"],
                    f["H", 1]["Length"],
                    f["T", 1]["Length"],
                    f["T", 1]["T", 2]["Length"]
                }
            ], 1]
    }]
};

(* ===================== Possible Issues (3) ===================== *)

possIssuesContent = {
    exampleSubsection["Exponential frame growth under non-Clifford circuits", {
        Cell[TextData[{
            "Each ", inlineCode[RowBox[{"\"P\"", "[", "\[Theta]", "]"}]],
            " gate doubles the component count, so a depth-",
            inlineMath[TI["d"]], " T-rich circuit can grow the frame to ",
            inlineMath[sup["2", TI["d"]]],
            " components; downstream ",
            inlineCode["\"StateVector\""],
            " materialisation cost scales accordingly:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f, history},
                f = StabilizerFrame[PauliStabilizer[1]];
                history = NestList[#["T", 1] &, f, 5];
                #["Length"] & /@ history
            ], 1]
    }],

    exampleSubsection["Equality is component-wise structural", {
        Cell[TextData[{
            inlineCode["Equal"],
            " compares the ordered component list (coefficient and PauliStabilizer ",
            "alike); two mathematically-equivalent frames with different component ",
            "orderings are ", inlineCode["False"], " under ", inlineCode["==="],
            "/", inlineCode["Equal"], ":"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{f1, f2},
                f1 = StabilizerFrame[{
                    {1, PauliStabilizer[1]},
                    {1, PauliStabilizer[{"Z"}]}
                }];
                f2 = StabilizerFrame[{
                    {1, PauliStabilizer[{"Z"}]},
                    {1, PauliStabilizer[1]}
                }];
                {f1 == f2, f1["State"] == f2["State"]}
            ], 1]
    }],

    exampleSubsection["StateVector unavailable for symbolic components", {
        Cell[TextData[{
            inlineCode["\"StateVector\""],
            " requires every component to materialise to a numeric vector; ",
            "passing a ", inlineCode["PauliStabilizer"], " whose tableau contains ",
            "symbolic signs prevents materialisation and the result remains held."
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]]
    }]
};

(* ===================== Neat Examples (1) ===================== *)

neatContent = {
    exampleSubsection["Symbolic T-gate fingerprint on a single qubit", {
        Cell[TextData[{
            "After a single ",
            inlineCode[RowBox[{"\"P\"", "[", "\[Theta]", "]"}]],
            " gate, the frame's coefficients form a recognisable trigonometric pair:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            StabilizerFrame[PauliStabilizer[1]]["P"[\[Theta]], 1]["Coefficients"] // Simplify, 1]
    }]
};

(* ===================== See Also ===================== *)

seeAlsoLinks = {
    "PauliStabilizer", "StabilizerStateQ", "CliffordChannel",
    "GraphState", "LocalComplement",
    "QuantumState", "QuantumOperator", "QuantumCircuitOperator"};

systemSeeAlso = {"PauliMatrix", "Total"};

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
    "stabilizer frame", "superposition", "non-Clifford gate",
    "T gate", "magic state", "stabilizer rank", "Quipu",
    "Garcia-Markov", "frame growth"};

(* ===================== Assemble the Notebook ===================== *)

notebook = Notebook[{

    Cell[CellGroupData[{
        Cell["StabilizerFrame", "ObjectName",
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
            Cell["Wolfram/QuantumFramework/ref/StabilizerFrame", "Categorization",
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
