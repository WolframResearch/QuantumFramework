#!/usr/bin/env wolframscript
(* ============================================================================
   build_paulistabilizer_doc.wl
   ----------------------------------------------------------------------------
   Generates a comprehensive PauliStabilizer reference page at
       QuantumFramework/Documentation/English/ReferencePages/Symbols/PauliStabilizer.nb

   Mirrors the QuantumEvolve.nb structural model.

   Run:  wolframscript -file OngoingProjects/Stabilizer/build_paulistabilizer_doc.wl
   ============================================================================ *)

$repoRoot  = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
$pacletDir = FileNameJoin[{$repoRoot, "QuantumFramework"}];
$docFile   = FileNameJoin[{$pacletDir, "Documentation", "English",
    "ReferencePages", "Symbols", "PauliStabilizer.nb"}];

Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

SeedRandom[20260510];   (* deterministic CellIDs *)

(* ===================== Helpers ===================== *)

cid[]  := RandomInteger[{10^6, 10^9 - 1}];
uuid[] := CreateUUID[];

(* Italic argument style *)
TI[s_String] := StyleBox[s, "TI"];

(* Subscript with italic letters *)
TIsub[a_String, b_String] := SubscriptBox[StyleBox[a, "TI"], StyleBox[b, "TI"]];

(* Function-link button for PauliStabilizer *)
psBtn[] := ButtonBox["PauliStabilizer", BaseStyle -> "Link",
    ButtonData -> "paclet:Wolfram/QuantumFramework/ref/PauliStabilizer"];

paclLink[sym_String, path_String] := ButtonBox[sym, BaseStyle -> "Link",
    ButtonData -> "paclet:" <> path];

modInfo[] := Cell["   ", "ModInfo", ExpressionUUID -> uuid[]];

(* inlineCode[boxes] — for inline code references in prose (PauliStabilizer,
   "Tableau", ps["Stabilizers"], etc.). Bare BoxData, no FormBox.  *)
inlineCode[boxes_] := Cell[BoxData[boxes], "InlineFormula", ExpressionUUID -> uuid[]];

(* inlineFormula[boxes] — legacy alias kept for the Usage cell (which uses code
   references inside InlineFormula cells per paclet convention). *)
inlineFormula[boxes_] := Cell[BoxData[boxes], "InlineFormula", ExpressionUUID -> uuid[]];

(* inlineMath[boxes] — for inline MATH in prose (O(n^2), (1-Signs)/2, 2n x 2n,
   subscripts, fractions). Wraps in FormBox[..., TraditionalForm] per nb-writer
   v1 Rule 1: math inside InlineFormula MUST have FormBox + TraditionalForm or
   it renders as plain text and subscript/superscript/fraction boxes break. *)
inlineMath[boxes_] := Cell[
    BoxData[FormBox[boxes, TraditionalForm]],
    "InlineFormula", ExpressionUUID -> uuid[]];

(* inlineKet[contentBoxes] — for inline kets |x⟩. Per nb-writer v1 Rule 16:
   kets use TemplateBox[{...}, "Ket"] wrapped in FormBox, and the Cell uses
   FormatType -> TraditionalForm (NOT the "InlineFormula" style — that
   doesn't provide the template rendering context for TemplateBox).  *)
inlineKet[contentBoxes_] := Cell[
    BoxData[FormBox[TemplateBox[{contentBoxes}, "Ket"], TraditionalForm]],
    FormatType -> TraditionalForm,
    ExpressionUUID -> uuid[]];

(* Common math shorthands *)
sub[base_, sub_] := SubscriptBox[base, sub];
sup[base_, exp_] := SuperscriptBox[base, exp];
frac[num_, den_] := FractionBox[num, den];

(* One usage row: ModInfo + InlineFormula(call) + "\\[LineSeparator]" + description *)
usageRow[callBoxes_, descParts___] := Sequence[
    modInfo[],
    inlineFormula[callBoxes],
    Sequence @@ Riffle[{"\[LineSeparator]"}, {descParts}, {1, -2, 2}],
    descParts,
    "\n"
];

(* simpler: just take description as a single TextData-compatible spec *)
uRow[callBoxes_, descSpec_] := Sequence[
    modInfo[],
    inlineFormula[callBoxes],
    "\[LineSeparator]",
    descSpec,
    "\n"
];

(* Build an Input/Output cell pair from an unevaluated expression.

   PRE-FIX BUG (root cause of "input shows summary box, output shows result$NNN"):
   - `MakeBoxes[expr, StandardForm]` was evaluating expr first (because parameter
     binding lifts the held form, then default makeboxes on PauliStabilizer[3]
     evaluates it and the custom summary-box upvalue fires).
   - `MakeBoxes[result, StandardForm]` was typesetting the SYMBOL `result$NNN`
     (Module-renamed local) because MakeBoxes has HoldAllComplete and never
     dereferences the binding.

   FIX:
   - For the input cell, extract the SOURCE STRING via
     `ToString[Unevaluated[expr], InputForm]` and parse it via the front-end
     parser (`FrontEnd\`UndocumentedTestFEParserPacket`). This guarantees the
     input cell shows the typed code, not the evaluated object.
   - For the output cell, use `ToBoxes[result, ...]` (which evaluates result
     and typesets the value) — NOT `MakeBoxes`.
*)
SetAttributes[ioBlock, HoldFirst];
ioBlock[expr_, idx_Integer : 1] := With[{idxStr = ToString[idx]},
    Module[{srcStr, parsed, inBoxes, result, outBoxes},
        srcStr = ToString[Unevaluated[expr], InputForm, PageWidth -> Infinity];
        parsed = Quiet @ UsingFrontEnd @ MathLink`CallFrontEnd[
            FrontEnd`UndocumentedTestFEParserPacket[srcStr, False]];
        inBoxes = If[MatchQ[parsed, {_BoxData, _}],
            First[parsed][[1]],
            (* Fallback if FE parser is unavailable: emit the raw string.
               Per nb-writer (v1) Rule 5 this is suboptimal but at least the
               source code is visible. *)
            srcStr
        ];
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

(* Helper for assembling an example: caption + Input/Output, with continuation idx *)
SetAttributes[example, HoldRest];
example[caption_String, expr_] := Flatten @ {
    Cell[caption, "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    ioBlock[expr, 1]
};

(* Two consecutive expressions sharing input numbering *)
SetAttributes[example2, HoldRest];
example2[caption_String, expr1_, expr2_] := Flatten @ {
    Cell[caption, "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
    ioBlock[expr1, 1],
    ioBlock[expr2, 2]
};

(* Wrap a list of cells in a CellGroupData with an ExampleSubsection header *)
exampleSubsection[title_String, content_List] := Cell[CellGroupData[{
    Cell[title, "ExampleSubsection",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ content
}, Open], CellID -> cid[]];

(* Wrap an ExampleSection (a top-level section under ExtendedExamples). *)
exampleSection[title_String, content_List, opts___] := Cell[CellGroupData[{
    Cell[BoxData[InterpretationBox[
        Cell[title, "ExampleSection", ExpressionUUID -> uuid[]],
        $Line = 0; Null]], "ExampleSection",
        CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ content
}, Open], CellID -> cid[]];

(* A "leaf" empty ExampleSection (no children) — for sections like Interactive Examples *)
exampleSectionEmpty[title_String] := Cell[BoxData[InterpretationBox[
    Cell[title, "ExampleSection", ExpressionUUID -> uuid[]],
    $Line = 0; Null]], "ExampleSection",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* 3-column option/method table cell *)
threeColRow[name_String, default_, description_] := {
    Cell["      ", "ModInfo", ExpressionUUID -> uuid[]],
    Cell[name, "TableText", ExpressionUUID -> uuid[]],
    Cell[If[StringQ[default],
        TextData[Cell[BoxData[default], "InlineFormula", ExpressionUUID -> uuid[]]],
        TextData[Cell[BoxData[ToBoxes[default]], "InlineFormula", ExpressionUUID -> uuid[]]]
    ], "TableText", ExpressionUUID -> uuid[]],
    Cell[description, "TableText", ExpressionUUID -> uuid[]]
};

threeColTable[rows : {{___}..}] := Cell[BoxData[GridBox[
    rows
]], "3ColumnTableMod",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Usage Cell (9 patterns) ===================== *)

(* Helper for "PauliStabilizer[...]" inline call *)
psCall[argBoxes___] := RowBox[{psBtn[], "[", argBoxes, "]"}];

usageContent = Flatten[{
    (* 1: PauliStabilizer[] *)
    modInfo[],
    inlineFormula[psCall[]],
    "\[LineSeparator]represents an empty stabilizer state, equivalent to the single-qubit register \
",
    inlineFormula[StyleBox["|0\[RightAngleBracket]", "TI"]],
    ".\n",

    (* 2: PauliStabilizer[n] *)
    modInfo[],
    inlineFormula[psCall[TI["n"]]],
    "\[LineSeparator]represents the ",
    inlineFormula[TI["n"]],
    "-qubit ",
    inlineFormula[StyleBox["|0\[CenterEllipsis]0\[RightAngleBracket]", "TI"]],
    " register state, with stabilizers ",
    inlineFormula[SubscriptBox[StyleBox["Z", "TI"], StyleBox["i", "TI"]]],
    ".\n",

    (* 3: PauliStabilizer[{stab1, stab2, ...}] *)
    modInfo[],
    inlineFormula[psCall[RowBox[{"{",
        RowBox[{TIsub["stab","1"], ",", TIsub["stab","2"], ",", StyleBox["\[Ellipsis]","TR"]}],
        "}"}]]],
    "\[LineSeparator]constructs a stabilizer state from a list of Pauli strings.\n",

    (* 4: PauliStabilizer[{stabs},{destabs}] *)
    modInfo[],
    inlineFormula[psCall[RowBox[{
        "{", StyleBox["stabs","TI"], "}", ",",
        "{", StyleBox["destabs","TI"], "}"
    }]]],
    "\[LineSeparator]constructs a stabilizer state with explicit stabilizers and destabilizers.\n",

    (* 5: PauliStabilizer[name] *)
    modInfo[],
    inlineFormula[psCall[TI["name"]]],
    "\[LineSeparator]returns a named stabilizer state. Names: ",
    inlineFormula["\"5QubitCode\""], ", ",
    inlineFormula["\"SteaneCode\""], ", ",
    inlineFormula["\"9QubitCode\""], ", their ",
    inlineFormula["\"\[CenterDot]1\""],
    " siblings, and ",
    inlineFormula["\"Random\""], ".\n",

    (* 6: PauliStabilizer["Random", n] *)
    modInfo[],
    inlineFormula[psCall[RowBox[{"\"Random\"", ",", TI["n"]}]]],
    "\[LineSeparator]returns a uniformly random ",
    inlineFormula[TI["n"]],
    "-qubit Clifford state via the Mallows sampler.\n",

    (* 7: PauliStabilizer[qs] *)
    modInfo[],
    inlineFormula[psCall[TI["qs"]]],
    "\[LineSeparator]converts a ",
    inlineFormula[paclLink["QuantumState", "Wolfram/QuantumFramework/ref/QuantumState"]],
    ", ",
    inlineFormula[paclLink["QuantumOperator", "Wolfram/QuantumFramework/ref/QuantumOperator"]],
    ", or ",
    inlineFormula[paclLink["QuantumCircuitOperator",
        "Wolfram/QuantumFramework/ref/QuantumCircuitOperator"]],
    " to a stabilizer state (Clifford only).\n",

    (* 8: ps[method, args] *)
    modInfo[],
    inlineFormula[RowBox[{TI["ps"], "[", TI["method"], ",", StyleBox["\[Ellipsis]", "TR"], "]"}]],
    "\[LineSeparator]applies a dispatched gate, measurement, or property accessor.\n",

    (* 9: ps[other] *)
    modInfo[],
    inlineFormula[RowBox[{TI["ps"], "[", TIsub["ps", "2"], "]"}]],
    "\[LineSeparator]composes two stabilizer states (apply ",
    inlineFormula[TIsub["ps", "2"]],
    " first, then ",
    inlineFormula[TI["ps"]],
    ")."
}];

usageCell = Cell[TextData[usageContent], "Usage",
    CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Notes Cells (10 + option table) ===================== *)

notesCells = {

    (* Note 1: definition + O(n^2) memory *)
    Cell[TextData[{
        "A ", inlineCode[psBtn[]],
        " encodes an n-qubit pure stabilizer state by a tableau of Pauli generators \
plus signs, requiring ",
        inlineMath[RowBox[{"O", "(", sup["n", "2"], ")"}]],
        " memory."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 2: internal Association form + 2n shape + ±1 signs *)
    Cell[TextData[{
        "The internal representation is an Association ",
        inlineCode[RowBox[{"<|", RowBox[{
            "\"Tableau\"", "->", TI["t"], ",",
            "\"Signs\"", "->", TI["s"]
        }], "|>"}]],
        " where ", inlineCode[TI["t"]],
        " is a rank-3 binary array of shape ",
        inlineMath[RowBox[{"{", "2", ",", " ", "n", ",", " ", RowBox[{"2", "n"}], "}"}]],
        " (X-block, Z-block; n qubits; ",
        inlineMath[RowBox[{"2", "n"}]],
        " rows: destabilizers then stabilizers) and ",
        inlineCode[TI["s"]],
        " is a length-",
        inlineMath[RowBox[{"2", "n"}]],
        " list of ",
        inlineMath["\[PlusMinus]1"],
        " signs."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 3: serialization keys with proper fraction *)
    Cell[TextData[{
        "Equivalent serialization keys are also accepted: ",
        inlineCode["\"Phase\""],
        " (= ",
        inlineMath[frac[RowBox[{"1", "-", "Signs"}], "2"]],
        ") and ",
        inlineCode["\"Matrix\""],
        " (the flat ",
        inlineMath[RowBox[{RowBox[{"2", "n"}], "\[Times]", RowBox[{"2", "n"}]}]],
        " binary tableau)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 4: stabilizer / destabilizer roles *)
    Cell[TextData[{
        "Stabilizer rows are the n parity-check generators of the codespace. Destabilizer \
rows complete the symplectic basis and are required for measurement and ",
        inlineCode["\"Dagger\""],
        ". String-list constructors auto-generate destabilizers via ",
        inlineCode["Reverse"],
        "; supply them explicitly with the two-list form to avoid singular fixtures."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 5: named codes with proper kets |0_L> and |1_L> *)
    Cell[TextData[{
        "Named codes follow the convention that ",
        inlineCode["\"5QubitCode\""],
        " encodes ",
        inlineKet[sub["0", "L"]],
        " (the +1 eigenstate of logical ",
        inlineMath["Z"],
        "), and the ",
        inlineCode["\"\[CenterDot]1\""],
        " sibling encodes ",
        inlineKet[sub["1", "L"]],
        "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 6: Random sampler *)
    Cell[TextData[{
        inlineCode[RowBox[{psBtn[], "[", "\"Random\"", ",", TI["n"], "]"}]],
        " uses the Bravyi-Maslov / Koenig-Smolin (arXiv:1406.2170) Mallows sampler for uniform \
sampling over the ",
        inlineMath["n"],
        "-qubit Clifford group modulo Paulis."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 7: ps dispatch surface *)
    Cell[TextData[{
        "A stabilizer state ",
        inlineCode[TI["ps"]],
        " supports two access patterns: properties via ",
        inlineCode[RowBox[{TI["ps"], "[", "\"Property\"", "]"}]],
        " (e.g. ", inlineCode["\"Stabilizers\""], ", ", inlineCode["\"Tableau\""],
        ", ", inlineCode["\"State\""],
        ") and methods via ",
        inlineCode[RowBox[{TI["ps"], "[", "\"Method\"", ",",
            StyleBox["\[Ellipsis]", "TR"], "]"}]],
        " (e.g. Clifford gates, measurement, ", inlineCode["\"InnerProduct\""],
        "). The full list of accepted strings is available via ",
        inlineCode[RowBox[{TI["ps"], "[", "\"Properties\"", "]"}]],
        "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 8: gate updates in O(n^2), Clifford generators *)
    Cell[TextData[{
        "Clifford gate updates run in ",
        inlineMath[RowBox[{"O", "(", sup["n", "2"], ")"}]],
        " time using the Aaronson-Gottesman tableau (arXiv:quant-ph/0406196). \
Single-qubit gates: ",
        inlineCode["\"H\""], ", ",
        inlineCode["\"S\""], ", ",
        inlineMath[sup["\"S\"", "\[Dagger]"]],
        ", ",
        inlineCode["\"X\""], ", ", inlineCode["\"Y\""], ", ",
        inlineCode["\"Z\""], ", ", inlineCode["\"V\""],
        ". Two-qubit gates: ",
        inlineCode["\"CNOT\""], " (or ", inlineCode["\"CX\""], "), ",
        inlineCode["\"CZ\""], ", ", inlineCode["\"SWAP\""],
        "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 9: Non-Clifford promotion to StabilizerFrame *)
    Cell[TextData[{
        "Non-Clifford gates ",
        inlineMath[RowBox[{"\"P\"", "[", "\[Theta]", "]"}]],
        ", ",
        inlineCode["\"T\""],
        ", and ",
        inlineMath[sup["\"T\"", "\[Dagger]"]],
        " return a ",
        inlineCode[paclLink["StabilizerFrame",
            "Wolfram/QuantumFramework/ref/StabilizerFrame"]],
        " (a coherent superposition of stabilizer states) rather than a ",
        inlineCode[psBtn[]],
        ". The frame closes under further Clifford gates."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 10: Symbolic measurement *)
    Cell[TextData[{
        "Symbolic measurement (Fang-Ying 2023, arXiv:2311.03906) via ",
        inlineCode["\"SymbolicMeasure\""],
        " allocates a fresh ",
        inlineMath[RowBox[{"\[FormalS]", "[", "k", "]"}]],
        " outcome symbol per non-deterministic measurement, yielding a single \
PauliStabilizer with symbolic phase. Resolve via ",
        inlineCode["\"SubstituteOutcomes\""], " or ",
        inlineCode["\"SampleOutcomes\""],
        "."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Note 11: Method option intro -> immediately followed by table *)
    Cell[TextData[{
        "The option ",
        inlineCode["Method"],
        " on ",
        inlineCode[RowBox[{TI["ps"], "[", "\"InnerProduct\"", ",", TI["other"], "]"}]],
        " selects between exact state-vector materialization (",
        inlineCode["\"Direct\""], ", default) and the closed-form Garcia-Markov-Cross \
algorithm (",
        inlineCode["\"ClosedForm\""],
        ", arXiv:1210.6646)."
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    threeColTable[{
        threeColRow["\"Method\"", "\"Direct\"",
            "the algorithm to use when computing the inner product"]
    }],

    (* Note 12: Method values *)
    Cell[TextData[{
        "Possible ",
        inlineCode["\"Method\""],
        " values:"
    }], "Notes", CellID -> cid[], ExpressionUUID -> uuid[]],

    threeColTable[{
        threeColRow["\"Direct\"", "default",
            "materialize state vectors and compute exactly (O(2^n), recovers complex phase)"],
        threeColRow["\"ClosedForm\"", "",
            "closed-form O(n^3) magnitude only, no complex phase (concrete signs only)"]
    }]
};

(* ===================== Primary Examples (5) ===================== *)

primaryExamplesContent = {

    (* Each example shows the actual PauliStabilizer OBJECT (its summary-box
       rendering with stabilizers / tableau panel) as the Output. Property
       access via ps["Stabilizers"] etc. is demonstrated in Scope, not here. *)

    (* Example 1: empty form *)
    Sequence @@ example[
        "The empty form constructs the single-qubit \"|0\[RightAngleBracket]\" register:",
        PauliStabilizer[]],

    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Example 2: n-qubit register *)
    Sequence @@ example[
        "For a multi-qubit register, pass the qubit count:",
        PauliStabilizer[3]],

    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Example 3: Bell state from Pauli strings *)
    Sequence @@ example[
        "Construct a Bell state from two Pauli-string stabilizers:",
        PauliStabilizer[{"XX", "ZZ"}]],

    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Example 4: named code *)
    Sequence @@ example[
        "Construct a named stabilizer code:",
        PauliStabilizer["5QubitCode"]],

    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Example 5: random Clifford *)
    Sequence @@ example[
        "Sample a random Clifford state via the Mallows sampler:",
        Module[{}, SeedRandom[42]; PauliStabilizer["Random", 3]]],

    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],

    (* Example 6: built by applying gates *)
    Sequence @@ example[
        "Build a 3-qubit GHZ state by applying Clifford gates:",
        PauliStabilizer[3]["H", 1]["CNOT", 1, 2]["CNOT", 2, 3]]
};

(* ===================== Scope (12 subsections) ===================== *)

scope1 = exampleSubsection["Constructors", {
    Sequence @@ example["The empty form returns the single-qubit register state:",
        PauliStabilizer[]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Specify a number of qubits to get an n-qubit register:",
        PauliStabilizer[5]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Build a state from a list of Pauli strings:",
        PauliStabilizer[{"XX", "ZZ"}]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Inspect the stabilizers via the property accessor:",
        PauliStabilizer[{"XX", "ZZ"}]["Stabilizers"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Provide explicit destabilizers as a second list:",
        PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Construct from a QuantumState via 4^n Pauli tomography:",
        PauliStabilizer[QuantumState[{1, 0, 0, 1}/Sqrt[2]]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Construct from a QuantumOperator that is Clifford:",
        PauliStabilizer[QuantumOperator["H", {1}]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Fold a QuantumCircuitOperator over the n-qubit register:",
        PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]]
}];

scope2 = exampleSubsection["Named stabilizer codes", {
    Sequence @@ example["The 5-qubit perfect code:",
        PauliStabilizer["5QubitCode"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["The Steane CSS code on 7 qubits:",
        PauliStabilizer["SteaneCode"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["The 9-qubit Shor code:",
        PauliStabilizer["9QubitCode"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["The \"\[CenterDot]1\" siblings encode the |1\[RightAngleBracket] logical (sign-flipped Z\[Macron]):",
        PauliStabilizer["5QubitCode1"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Read off the sign-flipped logical stabilizer:",
        Last[PauliStabilizer["5QubitCode1"]["Stabilizers"]]]
}];

scope3 = exampleSubsection["Random Clifford states", {
    Sequence @@ example["Sample a random 4-qubit Clifford state via the Mallows sampler:",
        Module[{}, SeedRandom[42]; PauliStabilizer["Random", 4]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Read off the random stabilizers:",
        Module[{}, SeedRandom[42]; PauliStabilizer["Random", 4]["Stabilizers"]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["The discrete spectrum of squared overlaps for 3-qubit pairs:",
        Module[{products},
            SeedRandom[42];
            products = Table[
                Abs[PauliStabilizer["Random", 3]["InnerProduct",
                    PauliStabilizer["Random", 3]]]^2,
                {200}];
            Tally[Round[products, 1/8]]]]
}];

scope4 = exampleSubsection["Properties and tableau views", {
    Sequence @@ example["List the available property names:",
        Length[PauliStabilizer[2]["Properties"]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["The first eight property names:",
        Take[PauliStabilizer[2]["Properties"], 8]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["The rank-3 binary tableau dimensions:",
        Dimensions[PauliStabilizer["5QubitCode"]["Tableau"]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Stabilizers and destabilizers separately:",
        Module[{ps},
            ps = PauliStabilizer[
                QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
            {ps["Stabilizers"], ps["Destabilizers"]}]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Materialize a dense state vector:",
        Normal @ PauliStabilizer[{"XX", "ZZ"}]["State"]["StateVector"]]
}];

scope5 = exampleSubsection["Single-qubit Clifford gates", {
    Sequence @@ example["Apply Hadamard:",
        PauliStabilizer[2]["H", 1]["Stabilizers"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Apply S, X, Y, Z, V to qubit 1:",
        Module[{ps},
            ps = PauliStabilizer[2];
            AssociationThread[
                {"S", "X", "Y", "Z", "V"} -> {
                    ps["S", 1]["Stabilizers"],
                    ps["X", 1]["Stabilizers"],
                    ps["Y", 1]["Stabilizers"],
                    ps["Z", 1]["Stabilizers"],
                    ps["V", 1]["Stabilizers"]}]]]
}];

scope6 = exampleSubsection["Two-qubit Clifford gates", {
    Sequence @@ example["Build a Bell state with H + CNOT:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["CZ acts symmetrically on both qubits:",
        PauliStabilizer[2]["H", 1]["H", 2]["CZ", 1, 2]["Stabilizers"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["SWAP via column permutation:",
        PauliStabilizer[3]["X", 1]["SWAP", 1, 3]["Stabilizers"]]
}];

scope7 = exampleSubsection["Permutation and padding", {
    Sequence @@ example["Permute qubit labels with a Cycles spec:",
        PauliStabilizer[3]["X", 1]["PermuteQudits", Cycles[{{1, 3}}]]["Stabilizers"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Pad to more qubits via tensor product with identity:",
        PauliStabilizer[{"X"}]["PadRight", 3]["Stabilizers"]]
}];

scope8 = exampleSubsection["Dagger / Inverse", {
    Sequence @@ example["Compute the inverse Clifford via the symplectic-matrix dagger:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Dagger"]["Stabilizers"]]
}];

scope9 = exampleSubsection["Measurement", {
    Sequence @@ example["Single-qubit Z-basis measurement returns an Association keyed by outcome:",
        Keys @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", 1]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Multi-qubit measurement via a list of qubits:",
        Keys @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", {1, 2}]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Measure an arbitrary Pauli-string observable:",
        Keys @ PauliStabilizer["5QubitCode"]["M", "XZZXI"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Symbolic measurement embeds a fresh outcome symbol in the phase:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]["Stabilizers"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Substitute the symbol back to recover concrete branches:",
        Module[{psSym},
            psSym = PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1];
            psSym["SubstituteOutcomes", \[FormalS][_] -> 0]["Stabilizers"]]]
}];

scope10 = exampleSubsection["Inner product and expectation", {
    Sequence @@ example["Self-inner-product is 1:",
        Module[{bell}, bell = PauliStabilizer[{"XX", "ZZ"}]; bell["InnerProduct", bell]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Sign-flipped Bell states are orthogonal:",
        PauliStabilizer[{"XX", "ZZ"}]["InnerProduct", PauliStabilizer[{"-XX", "ZZ"}]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Expectation of a stabilizer-group element returns +1:",
        PauliStabilizer[{"XX", "ZZ"}]["Expectation", "XX"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["An anticommuting Pauli gives expectation 0:",
        PauliStabilizer[{"XX", "ZZ"}]["Expectation", "XI"]]
}];

scope11 = exampleSubsection["Composition", {
    Sequence @@ example["Compose two stabilizer states (apply second, then first):",
        Module[{ps1, ps2},
            ps1 = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1}]];
            ps2 = PauliStabilizer[QuantumCircuitOperator[{"S" -> 1}]];
            ps1[ps2]["Stabilizers"]]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Apply a QuantumOperator via UpValue:",
        Module[{qo}, qo = QuantumOperator["H", {1}];
            qo[PauliStabilizer[1]]["Stabilizers"]]]
}];

scope12 = exampleSubsection["Round-trips and interop", {
    Sequence @@ example["Convert a stabilizer state to a QuantumState:",
        Normal @ QuantumState[PauliStabilizer[{"XX", "ZZ"}]]["StateVector"]],
    Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
    Sequence @@ example["Apply a Pauli QuantumMeasurementOperator natively (fast path):",
        Module[{bell}, bell = PauliStabilizer[{"XX", "ZZ"}];
            Keys @ QuantumMeasurementOperator["XX", {1, 2}][bell]]]
}];

scopeContent = {scope1, scope2, scope3, scope4, scope5, scope6, scope7, scope8,
    scope9, scope10, scope11, scope12};

(* ===================== Generalizations & Extensions (1) ===================== *)

genExtContent = {
    exampleSubsection["Non-Clifford promotion to StabilizerFrame", {
        Sequence @@ example["Applying T promotes a stabilizer state to a StabilizerFrame:",
            Module[{psT}, psT = PauliStabilizer[1]["T", 1]; Head @ psT]],
        Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ example["Successive non-Clifford gates double the frame size:",
            PauliStabilizer[1]["T", 1]["T", 1]["Length"]]
    }]
};

(* ===================== Options (literal: 30+ method subsections) ===================== *)

(* Helper for one-method examples *)
methodSub[title_, caption_, exprHeld_] := exampleSubsection[title, {
    Sequence @@ example[caption, Evaluate@exprHeld]
}];

SetAttributes[ms, HoldRest];
ms[title_String, caption_String, expr_] := exampleSubsection[title, {
    Sequence @@ example[caption, expr]
}];

optionsContent = {

    (* The actual Method option *)
    exampleSubsection["Method", {
        Sequence @@ example[
            "The default Method is \"Direct\" (state-vector materialization):",
            PauliStabilizer[2]["H", 1]["CNOT", 1, 2][
                "InnerProduct", PauliStabilizer[2]]],
        Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ example[
            "Explicit \"Direct\" recovers the full complex amplitude:",
            PauliStabilizer[2]["H", 1]["CNOT", 1, 2][
                "InnerProduct", PauliStabilizer[2], Method -> "Direct"]],
        Cell[BoxData[""], "ExampleDelimiter", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ example[
            "\"ClosedForm\" returns the magnitude in O(n\.b3) (Garcia-Markov-Cross 2012):",
            PauliStabilizer[2]["H", 1]["CNOT", 1, 2][
                "InnerProduct", PauliStabilizer[2], Method -> "ClosedForm"]]
    }],

    (* Single-qubit Clifford gates *)
    ms["\"H\"", "Hadamard maps X \[Rule] Z:",
        PauliStabilizer[{"X"}]["H", 1]["Stabilizers"]],
    ms["\"S\"", "S maps X \[Rule] Y:",
        PauliStabilizer[{"X"}]["S", 1]["Stabilizers"]],
    ms["SuperDagger[\"S\"]", "S\[Dagger] maps Y \[Rule] X:",
        PauliStabilizer[{"Y"}][SuperDagger["S"], 1]["Stabilizers"]],
    ms["\"X\"", "X anticommutes with Z:",
        PauliStabilizer[{"Z"}]["X", 1]["Stabilizers"]],
    ms["\"Y\"", "Y anticommutes with Z:",
        PauliStabilizer[{"Z"}]["Y", 1]["Stabilizers"]],
    ms["\"Z\"", "Z anticommutes with X:",
        PauliStabilizer[{"X"}]["Z", 1]["Stabilizers"]],
    ms["\"V\"", "V = \[Sqrt]X maps Z \[Rule] -Y:",
        PauliStabilizer[{"Z"}]["V", 1]["Stabilizers"]],
    ms["SuperDagger[\"V\"]", "V\[Dagger] maps Z \[Rule] Y:",
        PauliStabilizer[{"Z"}][SuperDagger["V"], 1]["Stabilizers"]],

    (* Two-qubit Clifford gates *)
    ms["\"CNOT\"", "Build a Bell state with H + CNOT:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]],
    ms["\"CX\"", "CX is an alias for CNOT:",
        PauliStabilizer[2]["H", 1]["CX", 1, 2]["Stabilizers"]],
    ms["\"CZ\"", "CZ propagates X-on-target into Z-on-control:",
        PauliStabilizer[{"IX", "XI"}]["CZ", 1, 2]["Stabilizers"]],
    ms["\"SWAP\"", "SWAP exchanges qubit labels:",
        PauliStabilizer[{"IX", "ZI"}]["SWAP", 1, 2]["Stabilizers"]],

    (* Permutation *)
    ms["\"Permute\"", "Permute swaps tableau rows:",
        PauliStabilizer[{"XI", "IZ"}]["Permute", Cycles[{{1, 2}}]]["Stabilizers"]],
    ms["\"PermuteQudits\"", "PermuteQudits swaps qubit columns:",
        PauliStabilizer[{"XI", "IZ"}]["PermuteQudits", Cycles[{{1, 2}}]]["Stabilizers"]],

    (* Padding and Dagger/Inverse *)
    ms["\"PadRight\"", "Tensor identity on the right:",
        PauliStabilizer[{"X"}]["PadRight", 3]["Stabilizers"]],
    ms["\"PadLeft\"", "Tensor identity on the left:",
        PauliStabilizer[{"X"}]["PadLeft", 3]["Stabilizers"]],
    ms["\"Dagger\"", "F\[NestedLessLess]2\[NestedGreaterGreater] symplectic inverse:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Dagger"]["Stabilizers"]],
    ms["\"Inverse\"", "Inverse is an alias for Dagger:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Inverse"]["Stabilizers"]],

    (* Non-Clifford boundary *)
    ms["\"P\"[\[Theta]]", "Phase rotation P[\[Theta]] returns a StabilizerFrame:",
        Head @ PauliStabilizer[1]["P"[Pi/4], 1]],
    ms["\"T\"", "T = P[\[Pi]/2] returns a StabilizerFrame:",
        Head @ PauliStabilizer[1]["T", 1]],
    ms["SuperDagger[\"T\"]", "T\[Dagger] returns a StabilizerFrame:",
        Head @ PauliStabilizer[1][SuperDagger["T"], 1]],

    (* Measurement *)
    ms["\"M\" (single qubit)", "Z-basis measurement on qubit q:",
        Keys @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", 1]],
    ms["\"M\" (qubit list)", "Joint Z measurement on a list of qubits:",
        Keys @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", {1, 2}]],
    ms["\"M\" (Pauli string)", "Measurement of an arbitrary Pauli observable:",
        Keys @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", "ZI"]],
    ms["\"SymbolicMeasure\"", "Allocates a fresh \[FormalS][k] outcome symbol:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]["Phase"]],
    ms["\"SubstituteOutcomes\"", "Resolve outcome symbols to concrete values:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1][
            "SubstituteOutcomes", \[FormalS][_] -> 1]["Phase"]],
    ms["\"SampleOutcomes\"", "Draw n random samples by independent substitution:",
        Length @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1][
            "SampleOutcomes", 4]],
    ms["\"RowSum\"", "AG row-multiplication primitive:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["RowSum", 1, 2]["Stabilizers"]],

    (* Inner product and expectation *)
    ms["\"InnerProduct\"", "Inner product against another stabilizer state:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["InnerProduct", PauliStabilizer[2]]],
    ms["\"Expectation\"", "Expectation value of a Pauli observable:",
        PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Expectation", "XX"]],

    (* Property accessor sweep *)
    exampleSubsection["Property accessors", {
        Cell[TextData[{
            "PauliStabilizer exposes ~30 string-keyed property accessors. The full list ",
            "is available via ", inlineFormula[RowBox[{TI["ps"], "[", "\"Properties\"", "]"}]],
            ". A representative subset:"
        }], "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps},
                ps = PauliStabilizer[{"XX", "ZZ"}];
                {ps["Qubits"], ps["GeneratorCount"], ps["Stabilizers"], ps["Destabilizers"],
                 Dimensions @ ps["Tableau"], Dimensions @ ps["Matrix"], Head @ ps["State"]}], 1]
    }]
};

(* ===================== Applications (3) ===================== *)

applicationsContent = {

    exampleSubsection["Stabilizer-state tomography", {
        Cell["Recover a stabilizer state from a dense state vector via 4^n tomography \
and round-trip back:", "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{qs, ps},
                qs = QuantumState[{1, 0, 0, 1}/Sqrt[2]];
                ps = PauliStabilizer[qs];
                ps["Stabilizers"]], 1]
    }],

    exampleSubsection["Random Clifford circuits and statistical fingerprint", {
        Cell["Sample 200 random 5-qubit Cliffords and tabulate their stabilizer-weight distribution:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{weights},
                SeedRandom[42];
                weights = Table[
                    Total[StringCount[#, "X" | "Y" | "Z"] & /@
                        PauliStabilizer["Random", 5]["Stabilizers"]],
                    {200}];
                Tally[weights] // Sort], 1]
    }],

    exampleSubsection["Quantum error correction with the 5-qubit code", {
        Cell["Encode |0_L\[RightAngleBracket] in the 5-qubit code, inject an X error on qubit 3, \
and observe a non-trivial syndrome:", "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps0, psErr, syndrome},
                ps0 = PauliStabilizer["5QubitCode"];
                psErr = ps0["X", 3];
                syndrome = Sign /@ Table[psErr["Expectation", g],
                    {g, ps0["Stabilizers"]}];
                syndrome], 1]
    }]
};

(* ===================== Properties & Relations (5) ===================== *)

propsRelContent = {
    exampleSubsection["State conversion preserves the stabilizer up to global phase", {
        Sequence @@ example["The QuantumState round-trip is exact for circuit-built fixtures:",
            Module[{ps, qs, ps2},
                ps = PauliStabilizer[QuantumCircuitOperator[
                    {"H" -> 1, "CNOT" -> {1, 2}}]];
                qs = ps["State"];
                ps2 = PauliStabilizer[qs];
                ps["Stabilizers"] === ps2["Stabilizers"]]]
    }],
    exampleSubsection["Inner-product squared on a discrete spectrum", {
        Sequence @@ example["For random pairs the squared overlap is in {0, 1/2^k}:",
            Module[{vals},
                SeedRandom[1];
                vals = Table[
                    Abs[PauliStabilizer["Random", 3]["InnerProduct",
                        PauliStabilizer["Random", 3]]]^2, {50}];
                Union @ Round[vals, 1/8]]]
    }],
    exampleSubsection["StabilizerStateQ on a PauliStabilizer is True", {
        Sequence @@ example["The predicate accepts PauliStabilizer values:",
            StabilizerStateQ[PauliStabilizer[{"XX", "ZZ"}]]]
    }],
    exampleSubsection["Bitflip channel as a tableau-mixture list", {
        Sequence @@ example["A BitFlip channel returns a probability-weighted mixture:",
            Length @ QuantumChannel["BitFlip"[1/3], {1}][PauliStabilizer[1]]]
    }],
    exampleSubsection["Inverse of a Clifford state", {
        Sequence @@ example["Composing ps with its Dagger yields a tableau equal to the identity \
register up to row reordering:",
            Module[{ps, dag},
                ps = PauliStabilizer[QuantumCircuitOperator[
                    {"H" -> 1, "CNOT" -> {1, 2}}]];
                dag = ps["Dagger"];
                Sort @ ps[dag]["Stabilizers"]]]
    }]
};

(* ===================== Possible Issues (6 messages) ===================== *)

possIssuesContent = {

    exampleSubsection["PauliStabilizer::nonclifford", {
        Cell["Non-Clifford gates without a tableau update emit ::nonclifford:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Quiet[Check[
                PauliStabilizerApply[QuantumCircuitOperator[{"X" -> 1, "Rx"[Pi/3] -> 1}],
                    PauliStabilizer[1]],
                "$Failed"]], 1]
    }],

    exampleSubsection["PauliStabilizer::singular", {
        Cell["Stabilizer-only fixtures may produce a singular tableau on Dagger:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Quiet[Check[
                PauliStabilizer[{"XX", "ZZ"}]["Dagger"], "$Failed"]], 1]
    }],

    exampleSubsection["PauliStabilizer::tdeprecated", {
        Cell["Modern paths return StabilizerFrame for non-Clifford T gates:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Head @ PauliStabilizer[1]["T", 1], 1]
    }],

    exampleSubsection["PauliStabilizer::expectationdim", {
        Cell["The Pauli string in Expectation must match the qubit count:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Quiet[Check[
                PauliStabilizer[2]["Expectation", "X"], "$Failed"]], 1]
    }],

    exampleSubsection["PauliStabilizer::partition", {
        Cell["Out-of-range qubit indices in measurement are rejected:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Quiet[Check[
                PauliStabilizer[2]["M", {1, 5}], "$Failed"]], 1]
    }],

    exampleSubsection["PauliStabilizer::nonpaulibasis", {
        Cell["Non-Pauli QMO basis falls back to the generic circuit path:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{qmo},
                qmo = QuantumMeasurementOperator[QuantumOperator[
                    DiagonalMatrix[{1, 2}]], {1}];
                Quiet[Check[Head @ qmo[PauliStabilizer[1]], "$Failed"]]], 1]
    }]
};

(* ===================== Neat Examples (1) ===================== *)

neatContent = {
    exampleSubsection["50-qubit random Clifford tableau as MatrixPlot", {
        Cell["Visualize the symplectic structure of a large random Clifford:",
            "ExampleText", CellID -> cid[], ExpressionUUID -> uuid[]],
        Sequence @@ ioBlock[
            Module[{ps},
                SeedRandom[1];
                ps = PauliStabilizer["Random", 50];
                MatrixPlot[ps["Matrix"], FrameTicks -> None,
                    ImageSize -> 250]], 1]
    }]
};

(* ===================== See Also ===================== *)

seeAlsoLinks = {
    "StabilizerFrame", "StabilizerStateQ", "GraphState", "LocalComplement",
    "CliffordChannel", "QuantumState", "QuantumOperator",
    "QuantumCircuitOperator", "QuantumMeasurementOperator", "QuantumChannel"
};

systemSeeAlso = {"PauliMatrix", "KroneckerProduct"};

seeAlsoCell = Cell[TextData[Riffle[
    Join[
        Cell[BoxData[
            ButtonBox[#, BaseStyle -> "Link",
                ButtonData -> "paclet:Wolfram/QuantumFramework/ref/" <> #]],
            "InlineSeeAlsoFunction",
            TaggingRules -> {"PageType" -> "Function"},
            ExpressionUUID -> uuid[]] & /@ seeAlsoLinks,
        Cell[BoxData[
            ButtonBox[#, BaseStyle -> "Link",
                ButtonData -> "paclet:ref/" <> #]],
            "InlineSeeAlsoFunction",
            TaggingRules -> {"PageType" -> "Function"},
            ExpressionUUID -> uuid[]] & /@ systemSeeAlso
    ],
    Cell[TextData[StyleBox[" \[FilledVerySmallSquare] ", "InlineSeparator"]],
        ExpressionUUID -> uuid[]]
]], "SeeAlso", CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Related Guides ===================== *)

moreAboutCell = Cell[TextData[
    Cell[BoxData[ButtonBox[
        "Wolfram Quantum Computation Framework", BaseStyle -> "Link",
        ButtonData -> "paclet:Wolfram/QuantumFramework/guide/WolframQuantumComputationFramework"
    ]], "InlineFormula", ExpressionUUID -> uuid[]]
], "MoreAbout", CellID -> cid[], ExpressionUUID -> uuid[]];

(* ===================== Keywords ===================== *)

keywordsCells = MapIndexed[
    Cell[#1, "Keywords",
        CellID -> cid[], ExpressionUUID -> uuid[]] &,
    {"stabilizer", "Clifford", "tableau", "Aaronson-Gottesman",
     "Pauli", "quantum error correction", "graph state",
     "random Clifford", "Mallows sampler", "SymPhase"}];

(* ===================== Build the full notebook ===================== *)

notebook = Notebook[{

    (* ----- Top group: ObjectName + Usage + Notes ----- *)
    Cell[CellGroupData[{
        Cell["PauliStabilizer", "ObjectName",
            CellID -> 2097139113, ExpressionUUID -> uuid[]],
        usageCell,
        Sequence @@ notesCells
    }, Open]],

    (* ----- See Also ----- *)
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

    (* ----- Tech Notes (empty) ----- *)
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

    (* ----- Related Guides ----- *)
    Cell[CellGroupData[{
        Cell["Related Guides", "MoreAboutSection",
            CellID -> cid[], ExpressionUUID -> uuid[]],
        moreAboutCell
    }, Open]],

    (* ----- Related Links (empty) ----- *)
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

    (* ----- Examples Initialization ----- *)
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

    (* ----- Primary Examples ----- *)
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

    (* ----- Extended Examples ----- *)
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

    (* ----- Metadata footer ----- *)
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
            Cell["Symbol", "Categorization",
                CellLabel -> "Entity Type",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Wolfram/QuantumFramework", "Categorization",
                CellLabel -> "Paclet Name",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Wolfram`QuantumFramework`", "Categorization",
                CellLabel -> "Context",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell["Wolfram/QuantumFramework/ref/PauliStabilizer", "Categorization",
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
            Cell[BoxData[""], "Template",
                CellLabel -> "Additional Function Template",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template",
                CellLabel -> "Arguments Pattern",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template",
                CellLabel -> "Local Variables",
                CellID -> cid[], ExpressionUUID -> uuid[]],
            Cell[BoxData[""], "Template",
                CellLabel -> "Color Equal Signs",
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

(* ===================== Write the notebook ===================== *)

Print["Writing notebook to: ", $docFile];
Export[$docFile, notebook, "Notebook"];
Print["Done. File size: ", FileByteCount[$docFile], " bytes"];

Exit[0];
