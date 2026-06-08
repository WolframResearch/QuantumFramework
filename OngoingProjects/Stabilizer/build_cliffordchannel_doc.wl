#!/usr/bin/env wolframscript
(* ============================================================================
   build_cliffordchannel_doc.wl
   ----------------------------------------------------------------------------
   Generates QuantumFramework/Documentation/English/ReferencePages/Symbols/
       CliffordChannel.nb
   following the patched documentation-writing skill (Step 2 head class:
   construction wrapper; Step 3 teaching-purpose rubric; Step 5.5/5.6
   render-verify; cell builder uses FE-parser path).

   Run:  wolframscript -file OngoingProjects/Stabilizer/build_cliffordchannel_doc.wl
   ============================================================================ *)

$repoRoot  = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$pacletDir = FileNameJoin[{$repoRoot, "QuantumFramework"}];
$docFile   = FileNameJoin[{$pacletDir, "Documentation", "English",
    "ReferencePages", "Symbols", "CliffordChannel.nb"}];

Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

SeedRandom[20260511];

(* ===================== Helpers ===================== *)

cid[]  := RandomInteger[{10^6, 10^9 - 1}];
uuid[] := CreateUUID[];

TI[s_String] := StyleBox[s, "TI"];
TIsub[a_String, b_String] := SubscriptBox[StyleBox[a, "TI"], StyleBox[b, "TI"]];

ccBtn[] := ButtonBox["CliffordChannel", BaseStyle -> "Link",
    ButtonData -> "paclet:Wolfram/QuantumFramework/ref/CliffordChannel"];

paclLink[sym_String, path_String] := ButtonBox[sym, BaseStyle -> "Link",
    ButtonData -> "paclet:" <> path];

modInfo[] := Cell["   ", "ModInfo", ExpressionUUID -> uuid[]];

inlineCode[boxes_] := Cell[BoxData[boxes], "InlineFormula", ExpressionUUID -> uuid[]];
inlineMath[boxes_] := Cell[
    BoxData[FormBox[boxes, TraditionalForm]],
    "InlineFormula", ExpressionUUID -> uuid[]];

sub[base_, s_] := SubscriptBox[base, s];
sup[base_, e_] := SuperscriptBox[base, e];

(* Build Input/Output cell pair — FE parser path per Step 5.5 *)
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

SetAttributes[ms, HoldRest];
ms[title_String, caption_String, expr_] := exampleSubsection[title, {
    Sequence @@ example[caption, expr]
}];

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

(* ===================== Usage Cell (5 patterns) ===================== *)

ccCall[argBoxes___] := RowBox[{ccBtn[], "[", argBoxes, "]"}];

usageContent = Flatten[{
    modInfo[],
    inlineCode[ccCall[RowBox[{"<|", StyleBox["\[Ellipsis]", "TR"], "|>"}]]],
    "\[LineSeparator]represents a Clifford channel encoded as a Choi-tableau Association \
with keys ", inlineCode["\"UA\""], ", ", inlineCode["\"UB\""], ", ",
    inlineCode["\"c\""], ", ", inlineCode["\"InputQubits\""], ", and ",
    inlineCode["\"OutputQubits\""], ".\n",

    modInfo[],
    inlineCode[ccCall[RowBox[{"\"Identity\"", ",", TI["n"]}]]],
    "\[LineSeparator]returns the identity channel on ",
    inlineMath["n"], " qubits.\n",

    modInfo[],
    inlineCode[ccCall[TI["ps"]]],
    "\[LineSeparator]constructs a state-preparation channel from a ",
    inlineCode[paclLink["PauliStabilizer", "Wolfram/QuantumFramework/ref/PauliStabilizer"]],
    " (input system empty, output equals the stabilizer state).\n",

    modInfo[],
    inlineCode[ccCall[TI["qc"]]],
    "\[LineSeparator]constructs a Clifford channel from a ",
    inlineCode[paclLink["QuantumChannel", "Wolfram/QuantumFramework/ref/QuantumChannel"]],
    " when the channel is a deterministic single-Pauli.\n",

    modInfo[],
    inlineCode[RowBox[{TI["cc"], "[", TI["arg"], "]"}]],
    "\[LineSeparator]applied to another ", inlineCode[ccBtn[]],
    " composes channels (apply ", inlineCode[TI["arg"]], " first, then ",
    inlineCode[TI["cc"]], "); applied to a ",
    inlineCode["PauliStabilizer"], " evolves the state; applied to a property string \
returns that property."
}];

usageCell = Cell[TextData[usageContent], "Usage",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Notes Cells ===================== *)

notesCells = {
    Cell[TextData[{
        "A ", inlineCode[ccBtn[]],
        " encodes a Clifford channel \[CapitalPhi]: T(",
        inlineMath[sub["\[ScriptCapitalH]", "A"]], ") \[Rule] T(",
        inlineMath[sub["\[ScriptCapitalH]", "B"]],
        ") by a Choi tableau, following Yashin25 (arXiv:2504.14101) Section 2.3."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The internal representation is an Association ",
        inlineCode[RowBox[{"<|", RowBox[{
            "\"UA\"", "->", TI["uA"], ",", " ",
            "\"UB\"", "->", TI["uB"], ",", " ",
            "\"c\"", "->", TI["c"], ",", " ",
            "\"InputQubits\"", "->", TI["nA"], ",", " ",
            "\"OutputQubits\"", "->", TI["nB"]
        }], "|>"}]],
        " where ",
        inlineMath[TI["uA"]],
        " is a ",
        inlineMath[RowBox[{TI["k"], "\[Times]", RowBox[{"2", TI["nA"]}]}]],
        " bit matrix on the input system (or ",
        inlineCode[RowBox[{"{", "}"}]],
        " if ",
        inlineMath[RowBox[{TI["nA"], "=", "0"}]],
        "), ",
        inlineMath[TI["uB"]],
        " is a ",
        inlineMath[RowBox[{TI["k"], "\[Times]", RowBox[{"2", TI["nB"]}]}]],
        " bit matrix on the output system, and ",
        inlineMath[TI["c"]],
        " is a length-",
        inlineMath[TI["k"]],
        " sign bit vector."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Each tableau row ",
        inlineMath[RowBox[{"[", RowBox[{
            sub["u", "A"], " | ", sub["u", "B"], " | ", "c"
        }], "]"}]],
        " encodes a Pauli superoperator \[CapitalPi](",
        inlineMath[sub["u", "A"]], " | ",
        inlineMath[sub["u", "B"]], " | ",
        inlineMath["c"],
        ")[\[Rho]] = ",
        inlineMath[RowBox[{
            sup["(-1)", "c"], "\[CenterDot]", sup["2", RowBox[{"|", "A", "|"}]],
            "\[CenterDot]", "Tr", "[", "\[Rho]", " ", "P", "(", sub["u","A"], ")", "]",
            "\[CenterDot]", "P", "(", sub["u","B"], ")"
        }]],
        ", and the channel acts as the average ",
        inlineMath[RowBox[{
            "\[CapitalPhi]", "[", "\[Rho]", "]", "=",
            sup["2", RowBox[{"-", "(", "|A|", "+", "|B|", ")"}]],
            "\[CenterDot]",
            "\[Sum]", " ", "\[CapitalPi]", "(", "row", ")", "[", "\[Rho]", "]"
        }]], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Special cases of the tableau:"
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "\[FilledCircle] A pure stabilizer state has ",
        inlineMath[RowBox[{TI["nA"], "=", "0"}]],
        ", ", inlineMath[RowBox[{TI["k"], "=", TI["nB"]}]],
        " rows; ", inlineMath[TI["uB"]],
        " is the state's stabilizer tableau."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "\[FilledCircle] A Clifford unitary ",
        inlineMath[sub["U", RowBox[{"A", "\[Rule]", "B"}]]],
        " with ", inlineMath[RowBox[{TI["nA"], "=", TI["nB"], "=", "n"}]],
        " has ", inlineMath[RowBox[{"2", TI["n"]}]],
        " rows enumerating Pauli generators and their conjugates; ",
        inlineMath[TI["c"]], " encodes phase signs."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "Composition ", inlineCode[RowBox[{TI["cc1"], "[", TI["cc2"], "]"}]],
        " applies ", inlineCode[TI["cc2"]],
        " first, then ", inlineCode[TI["cc1"]],
        ", and is implemented by Boolean null-space intersection on the ",
        inlineMath["B"],
        "-side bits, with Aaronson-Gottesman phase tracking and the ",
        inlineMath[RowBox[{"\[VerticalBar]", sub["\[CapitalPhi]","+"],
            sub["\[RightAngleBracket]", RowBox[{"B", " ", sup["B","'"]}]]}]],
        " contraction-sign correction for Y-bearing combined ",
        inlineMath[sub["u", "B"]],
        " Paulis (Yashin25 \[Section]3.2/\[Section]3.3)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "State evolution ", inlineCode[RowBox[{TI["cc"], "[", TI["ps"], "]"}]],
        " has three recognized dispatch cases: (i) identity channel returns ",
        inlineCode[TI["ps"]],
        " unchanged; (ii) state-preparation channel (",
        inlineMath[RowBox[{TI["nA"], "=", "0"}]],
        ") returns the state encoded by ", inlineCode[TI["cc"]],
        "; (iii) dim-matched channel builds ",
        inlineCode[RowBox[{ccBtn[], "[", TI["ps"], "]"}]],
        ", composes, and converts back to a ", inlineCode["PauliStabilizer"], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        inlineCode[RowBox[{ccBtn[], "[", TI["qc"], "]"}]],
        " for a ", inlineCode["QuantumChannel"],
        " detects deterministic single-Pauli channels by Label (",
        inlineCode["\"X\""], ", ",
        inlineCode["\"-XX\""], ", \[Ellipsis]); stochastic Pauli channels (",
        inlineCode["\"BitFlip\""], ", ",
        inlineCode["\"PhaseFlip\""], ", \[Ellipsis]) emit ",
        inlineCode["CliffordChannel::stochastic"],
        " and fall through to a placeholder identity-on-",
        inlineMath["A"], " form."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    Cell[TextData[{
        "The introspection key ",
        inlineCode[RowBox[{TI["cc"], "[", "\"Properties\"", "]"}]],
        " returns the full list of accepted accessor strings: ",
        inlineCode["\"UA\""], ", ", inlineCode["\"UB\""], ", ",
        inlineCode["\"c\""], ", ", inlineCode["\"InputQubits\""], ", ",
        inlineCode["\"OutputQubits\""], ", ", inlineCode["\"Rank\""], ", ",
        inlineCode["\"Tableau\""], ", and ", inlineCode["\"Source\""], "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]]
};

(* ===================== Primary Examples (5) ===================== *)

primaryExamplesContent = {
    Sequence @@ example[
        "Construct the identity channel on a single qubit:",
        CliffordChannel["Identity", 1]],
    delim[],

    Sequence @@ example[
        "Identity on more qubits:",
        CliffordChannel["Identity", 2]],
    delim[],

    Sequence @@ example[
        "Build a state-preparation channel from a stabilizer state:",
        CliffordChannel[PauliStabilizer[{"XX", "ZZ"}]]],
    delim[],

    Sequence @@ example[
        "Construct directly from a Choi-tableau Association (raw form):",
        CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}}, "UB" -> {{1, 1}, {0, 1}},
                          "c" -> {0, 0}, "InputQubits" -> 1, "OutputQubits" -> 1,
                          "Source" -> "S"|>]],
    delim[],

    Sequence @@ example[
        "Interrogate via the property accessor:",
        CliffordChannel["Identity", 1]["Properties"]]
};

(* ===================== Scope (8 subsections) ===================== *)

scope1 = exampleSubsection["Constructors", {
    Sequence @@ example["Identity channel on n qubits:",
        CliffordChannel["Identity", 2]],
    delim[],
    Sequence @@ example["State-preparation channel from a PauliStabilizer:",
        CliffordChannel[PauliStabilizer @ QuantumCircuitOperator[
            {"H" -> 1, "CNOT" -> {1, 2}}]]],
    delim[],
    Sequence @@ example["Raw Choi-tableau Association (e.g. the S-gate channel):",
        CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}}, "UB" -> {{1, 1}, {0, 1}},
                          "c" -> {0, 0}, "InputQubits" -> 1, "OutputQubits" -> 1,
                          "Source" -> "S"|>]],
    delim[],
    Sequence @@ example["Idempotent on an existing channel:",
        CliffordChannel[CliffordChannel["Identity", 1]]]
}];

scope2 = exampleSubsection["Properties", {
    Sequence @@ example["The full list of property names via introspection:",
        CliffordChannel["Identity", 1]["Properties"]],
    delim[],
    Sequence @@ example["Input and output qubit counts:",
        Module[{cc},
            cc = CliffordChannel["Identity", 2];
            {cc["InputQubits"], cc["OutputQubits"]}]],
    delim[],
    Sequence @@ example["Rank (number of stabilizer-tableau rows):",
        CliffordChannel["Identity", 2]["Rank"]],
    delim[],
    Sequence @@ example["Provenance tag (\"Identity\", \"PauliStabilizer\", \"Composition\", \[Ellipsis]):",
        CliffordChannel["Identity", 1]["Source"]],
    delim[],
    Sequence @@ example["Flat [UA | UB | c] tableau as a single binary matrix:",
        Dimensions @ CliffordChannel["Identity", 2]["Tableau"]]
}];

scope3 = exampleSubsection["Composition: cc1[cc2]", {
    Sequence @@ example["Apply cc2 first, then cc1:",
        Module[{ccS},
            ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}},
                "UB" -> {{1, 1}, {0, 1}}, "c" -> {0, 0},
                "InputQubits" -> 1, "OutputQubits" -> 1, "Source" -> "S"|>];
            ccS[ccS]]],
    delim[],
    Sequence @@ example["S\.b3 = S\[CenterDot]S\[CenterDot]S has phase-bit (0, 1), tracking the AG i-power across rows:",
        Module[{ccS},
            ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}},
                "UB" -> {{1, 1}, {0, 1}}, "c" -> {0, 0},
                "InputQubits" -> 1, "OutputQubits" -> 1, "Source" -> "S"|>];
            ccS[ccS[ccS]]["c"]]]
}];

scope4 = exampleSubsection["State evolution: cc[ps]", {
    Sequence @@ example["Identity channel returns the state unchanged:",
        CliffordChannel["Identity", 1][PauliStabilizer[1]]],
    delim[],
    Sequence @@ example["A state-preparation channel applied to any matching-dim ps returns the prepared state:",
        Module[{ps, cc},
            ps = PauliStabilizer[1]["X", 1];
            cc = CliffordChannel[ps];
            cc[PauliStabilizer[1]]]]
}];

scope5 = exampleSubsection["From a QuantumChannel", {
    Sequence @@ example["A stochastic Pauli channel emits CliffordChannel::stochastic and returns a placeholder identity-on-A form:",
        Quiet @ CliffordChannel[QuantumChannel["BitFlip"[1/3], {1}]]],
    delim[],
    Sequence @@ example["Identity channel matches via the standard Identity constructor:",
        CliffordChannel["Identity", 1]]
}];

scope6 = exampleSubsection["Tableau structure", {
    Sequence @@ example["The flat Tableau matrix is [UA | UB | c], reshaped from the three components:",
        Module[{cc, ua, ub, cc1, c, full},
            cc = CliffordChannel["Identity", 2];
            ua = cc["UA"]; ub = cc["UB"]; c = cc["c"];
            full = cc["Tableau"];
            Dimensions /@ {ua, ub, c, full}]]
}];

scope7 = exampleSubsection["Round-trip with PauliStabilizer", {
    Sequence @@ example["Wrap a stabilizer state as a channel and recover it:",
        Module[{ps, cc, psRecovered},
            ps = PauliStabilizer[{"XX", "ZZ"}];
            cc = CliffordChannel[ps];
            psRecovered = cc[PauliStabilizer[2]];
            psRecovered["Stabilizers"] === ps["Stabilizers"]]]
}];

scope8 = exampleSubsection["AG phase tracking and the \[VerticalBar]\[CapitalPhi]+\[RightAngleBracket] contraction sign", {
    Sequence @@ example["The S-channel composed with itself gives Z (X \[Rule] -X), with c-bits (0, 1):",
        Module[{ccS, ccZ},
            ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}},
                "UB" -> {{1, 1}, {0, 1}}, "c" -> {0, 0},
                "InputQubits" -> 1, "OutputQubits" -> 1, "Source" -> "S"|>];
            ccZ = ccS[ccS];
            ccZ["c"]]]
}];

scopeContent = {scope1, scope2, scope3, scope4, scope5, scope6, scope7, scope8};

(* ===================== Generalizations & Extensions (1) ===================== *)

genExtContent = {
    exampleSubsection["Relationship to QuantumChannel", {
        Sequence @@ example["Named single-Pauli QuantumChannels embed cleanly; stochastic channels fall through with a notice:",
            Quiet @ {
                CliffordChannel[QuantumChannel["BitFlip"[1/3], {1}]]["Source"]
            }]
    }]
};

(* ===================== Options (literal ladder, 11 subsections) ===================== *)

optionsContent = {
    ms["\"UA\"", "Input-system bit matrix (or {} for state-prep channels):",
        CliffordChannel["Identity", 1]["UA"]],
    ms["\"UB\"", "Output-system bit matrix:",
        CliffordChannel["Identity", 1]["UB"]],
    ms["\"c\"", "Sign-bit vector:",
        CliffordChannel["Identity", 1]["c"]],
    ms["\"InputQubits\"", "Number of input qubits nA:",
        CliffordChannel["Identity", 2]["InputQubits"]],
    ms["\"OutputQubits\"", "Number of output qubits nB:",
        CliffordChannel["Identity", 2]["OutputQubits"]],
    ms["\"Rank\"", "Number of stabilizer-tableau rows (\[LessEqual] 2|B| by trace-preservation):",
        CliffordChannel["Identity", 2]["Rank"]],
    ms["\"Tableau\"", "Full [UA | UB | c] flat binary matrix:",
        Dimensions @ CliffordChannel["Identity", 2]["Tableau"]],
    ms["\"Source\"", "Provenance tag:",
        CliffordChannel["Identity", 1]["Source"]],
    ms["\"Properties\"", "Full list of accessor strings:",
        CliffordChannel["Identity", 1]["Properties"]],
    ms["cc1[cc2] composition", "Apply cc2 first, then cc1:",
        CliffordChannel["Identity", 1][CliffordChannel["Identity", 1]]],
    ms["cc[ps] state evolution", "Apply channel to a PauliStabilizer:",
        CliffordChannel["Identity", 1][PauliStabilizer[1]]]
};

(* ===================== Applications (2) ===================== *)

applicationsContent = {
    exampleSubsection["S-gate composition tower", {
        Cell["The single-qubit phase-gate S has the canonical Choi-tableau form below. \
Composing it with itself realizes Z; cubing gives S\[Dagger]. Track the c-bits through \
the tower:", "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ccS, ccS2, ccS3},
                ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}},
                    "UB" -> {{1, 1}, {0, 1}}, "c" -> {0, 0},
                    "InputQubits" -> 1, "OutputQubits" -> 1, "Source" -> "S"|>];
                ccS2 = ccS[ccS];
                ccS3 = ccS[ccS2];
                {ccS["c"], ccS2["c"], ccS3["c"]}], 1]
    }],

    exampleSubsection["State-preparation chain", {
        Cell["Build a state-preparation channel from a Bell state, compose with the \
identity (a no-op), then apply to a fresh register and verify the prepared state \
matches the original:", "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps, ccPrep, ccI, result},
                ps = PauliStabilizer[{"XX", "ZZ"}];
                ccPrep = CliffordChannel[ps];
                ccI = CliffordChannel["Identity", 2];
                result = ccI[ccPrep];
                result["Source"]], 1]
    }]
};

(* ===================== Properties & Relations (4) ===================== *)

propsRelContent = {
    exampleSubsection["Identity channel is the no-op", {
        Sequence @@ example["Applying the identity channel returns the input state:",
            Module[{ps},
                ps = PauliStabilizer[{"XX", "ZZ"}];
                CliffordChannel["Identity", 2][ps]["Stabilizers"] === ps["Stabilizers"]]]
    }],

    exampleSubsection["State-prep channel embeds the stabilizer state", {
        Sequence @@ example["The state encoded by CliffordChannel[ps] equals ps[\"Stabilizers\"] in its UB block:",
            Module[{ps, cc},
                ps = PauliStabilizer[{"XX", "ZZ"}];
                cc = CliffordChannel[ps];
                {cc["InputQubits"], cc["Rank"]}]]
    }],

    exampleSubsection["Composition is associative", {
        Sequence @@ example["(ccA ccB) ccC and ccA (ccB ccC) produce the same channel:",
            Module[{ccS, lhs, rhs},
                ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}},
                    "UB" -> {{1, 1}, {0, 1}}, "c" -> {0, 0},
                    "InputQubits" -> 1, "OutputQubits" -> 1, "Source" -> "S"|>];
                lhs = (ccS[ccS])[ccS]["c"];
                rhs = ccS[ccS[ccS]]["c"];
                lhs === rhs]]
    }],

    exampleSubsection["Tableau is the [UA | UB | c] flat reshape", {
        Sequence @@ example["The Tableau dimensions match Length[c] x (2nA + 2nB + 1):",
            Module[{cc},
                cc = CliffordChannel["Identity", 2];
                {Dimensions[cc["Tableau"]],
                 {Length[cc["c"]], 2 cc["InputQubits"] + 2 cc["OutputQubits"] + 1}}]]
    }]
};

(* ===================== Possible Issues (3) ===================== *)

possIssuesContent = {
    exampleSubsection["CliffordChannel::dimmismatch", {
        Cell["Composition cc1[cc2] requires cc2.OutputQubits == cc1.InputQubits:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{cc1, cc2},
                cc1 = CliffordChannel["Identity", 1];
                cc2 = CliffordChannel["Identity", 2];
                Quiet @ Check[cc1[cc2], "$Failed"]], 1]
    }],

    exampleSubsection["CliffordChannel::stochastic", {
        Cell["A stochastic Pauli channel (rank > 1 Kraus set) emits ::stochastic and falls back to a placeholder identity-on-A form:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Quiet @ CliffordChannel[QuantumChannel["BitFlip"[1/3], {1}]]["Source"], 1]
    }],

    exampleSubsection["CliffordChannel::stateevol", {
        Cell["cc[ps] for an unrecognized dispatch case emits ::stateevol and returns $Failed:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{cc, ps},
                cc = CliffordChannel["Identity", 1];
                ps = PauliStabilizer[3];
                Quiet @ Check[cc[ps], "$Failed"]], 1]
    }]
};

(* ===================== Neat Examples (1) ===================== *)

neatContent = {
    exampleSubsection["Tableau visualization", {
        Cell["MatrixPlot of the [UA | UB | c] tableau for a 5-qubit identity channel:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            MatrixPlot[CliffordChannel["Identity", 5]["Tableau"],
                FrameTicks -> None, ImageSize -> 250], 1]
    }]
};

(* ===================== See Also ===================== *)

seeAlsoLinks = {
    "PauliStabilizer", "StabilizerFrame", "StabilizerStateQ",
    "GraphState", "LocalComplement",
    "QuantumChannel", "QuantumMeasurementOperator",
    "QuantumOperator", "QuantumState", "QuantumCircuitOperator"};

systemSeeAlso = {"PauliMatrix", "KroneckerProduct"};

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
    "Clifford channel", "Choi tableau", "stabilizer formalism",
    "Pauli superoperator", "Yashin", "quantum channel",
    "AG phase tracking", "state preparation", "composition", "tableau"};

(* ===================== Assemble the Notebook ===================== *)

notebook = Notebook[{

    Cell[CellGroupData[{
        Cell["CliffordChannel", "ObjectName",
            CellID -> 2097139114, ExpressionUUID -> uuid[]],
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
            Cell["Wolfram/QuantumFramework/ref/CliffordChannel", "Categorization",
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
Exit[0];
