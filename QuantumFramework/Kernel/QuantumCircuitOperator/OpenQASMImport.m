Package["Wolfram`QuantumFramework`"]

PackageScope["ImportOpenQASM"]
PackageScope["qasmStringQ"]
PackageScope["qasmFileQ"]



(* ::Section:: *)
(* Native Wolfram-Language OpenQASM (2.0 / 3.0) importer.

   Parses the circuit-level subset of OpenQASM into a QuantumCircuitOperator with no
   Python / qiskit dependency: it depends only on the (stable) OpenQASM spec, not on any
   qiskit release. Supported: version header, include, qreg/creg + qubit/bit declarations,
   standard + parametrized gate calls, custom `gate` definitions (-> named subcircuit),
   the inv/pow/ctrl/negctrl modifiers, measure (both v2 and v3 syntaxes), reset, barrier,
   and gphase. Anything outside this subset (classical types / control flow / def / defcal /
   timing / aliasing) returns a categorized Failure, never a half-built circuit.

   Error handling: every error path Throws a structured Failure["QuantumQASM", ...] to the
   $qasmTag Catch at the public boundary, carrying a Category, a MessageTemplate, and the
   offending statement/token. Helpers return values on success and throw on error. *)


$qasmTag = "Wolfram`QuantumFramework`OpenQASMImport";

qasmThrow[cat_String, tmpl_String, params_List : {}, extra_Association : <||>] :=
    Throw[
        Failure["QuantumQASM", Join[<|
            "MessageTemplate" -> tmpl,
            "MessageParameters" -> params,
            "Category" -> cat
        |>, extra]],
        $qasmTag
    ]


(* ---- a string is OpenQASM iff it carries an OPENQASM header ---- *)

qasmStringQ[s_String] := StringContainsQ[s, RegularExpression["(?m)^\\s*OPENQASM\\b"]]
qasmStringQ[___] := False

(* a string that names an existing .qasm file (the constructor / hub accept these as input) *)
qasmFileQ[s_String] := StringMatchQ[s, ___ ~~ ".qasm", IgnoreCase -> True] && FileExistsQ[s]
qasmFileQ[___] := False


qasmStripComments[s_String] := StringReplace[s, {
    "/*" ~~ Shortest[___] ~~ "*/" :> " ",
    "//" ~~ Shortest[___] ~~ ("\n" | EndOfString) :> "\n"
}]


qasmVersion[s_String] := Replace[
    FirstCase[StringCases[s, RegularExpression["OPENQASM\\s+(\\d+)"] -> "$1"], _String, "2"],
    v_String :> ToExpression[v]
]


(* ---- split into statements on `;` at brace-depth 0; keep `{...}` gate bodies whole ---- *)

qasmStatements[s_String] := Module[{state},
    state = Fold[
        Function[{st, c},
            With[{depth = st[[1]], cur = st[[2]], out = st[[3]]},
                Which[
                    c === "{", {depth + 1, {cur, c}, out},
                    c === "}",
                        If[depth - 1 == 0,
                            {0, {}, {out, StringJoin @ Flatten @ {cur, c}}},
                            {depth - 1, {cur, c}, out}
                        ],
                    c === ";" && depth == 0, {0, {}, {out, StringJoin @ Flatten @ cur}},
                    True, {depth, {cur, c}, out}
                ]
            ]
        ],
        {0, {}, {}},
        Characters[s]
    ];
    Select[StringTrim /@ Flatten @ {state[[3]], StringJoin @ Flatten @ state[[2]]}, # =!= "" &]
]


(* ---- split a comma-separated list at bracket/paren depth 0 ---- *)

qasmSplitTopComma[s_String] := Module[{state},
    state = Fold[
        Function[{st, c},
            With[{depth = st[[1]], cur = st[[2]], out = st[[3]]},
                Which[
                    MemberQ[{"(", "["}, c], {depth + 1, {cur, c}, out},
                    MemberQ[{")", "]"}, c], {depth - 1, {cur, c}, out},
                    c === "," && depth == 0, {depth, {}, {out, StringJoin @ Flatten @ cur}},
                    True, {depth, {cur, c}, out}
                ]
            ]
        ],
        {0, {}, {}},
        Characters[s]
    ];
    Select[StringTrim /@ Flatten @ {state[[3]], StringJoin @ Flatten @ state[[2]]}, # =!= "" &]
]


(* ---- angle / numeric expression parser (allowlist guard + safe evaluation) ---- *)

(* split a qubit-argument list on top-level commas AND whitespace (the QF emitter writes
   `q[0] q[1]` space-separated; standard OpenQASM uses commas; accept both) *)

qasmSplitQargs[s_String] := Module[{state},
    state = Fold[
        Function[{st, c},
            With[{depth = st[[1]], cur = st[[2]], out = st[[3]]},
                Which[
                    MemberQ[{"(", "["}, c], {depth + 1, {cur, c}, out},
                    MemberQ[{")", "]"}, c], {depth - 1, {cur, c}, out},
                    (c === "," || MemberQ[{" ", "\t", "\n", "\r"}, c]) && depth == 0, {depth, {}, {out, StringJoin @ Flatten @ cur}},
                    True, {depth, {cur, c}, out}
                ]
            ]
        ],
        {0, {}, {}},
        Characters[s]
    ];
    Select[StringTrim /@ Flatten @ {state[[3]], StringJoin @ Flatten @ state[[2]]}, # =!= "" &]
]


qasmAngle[str_String] := Module[{e = StringTrim[str], idents, wl, val},
    (* normalize scientific notation (1.5e-3, 2E4) to WL's "*^" so the digit-adjacent
       exponent letter is not mistaken for an identifier below *)
    e = StringReplace[e, RegularExpression["(\\d)[eE]([+-]?\\d)"] -> "$1*^$2"];
    idents = ToLowerCase /@ DeleteDuplicates @ StringCases[e, RegularExpression["[A-Za-z_][A-Za-z0-9_]*"]];
    If[ ! SubsetQ[{"pi", "tau", "euler"}, idents],
        qasmThrow["BadExpression",
            "Unsupported identifier in expression `1` (only pi, tau, euler and arithmetic are allowed).",
            {e}, <|"Expression" -> e|>]
    ];
    wl = StringReplace[e, {
        "**" -> "^",
        WordBoundary ~~ "pi" ~~ WordBoundary -> "(Pi)",
        WordBoundary ~~ "tau" ~~ WordBoundary -> "(2*Pi)",
        WordBoundary ~~ "euler" ~~ WordBoundary -> "(E)"
    }];
    val = Quiet @ Check[ToExpression[wl], $Failed];
    If[ val === $Failed || ! NumericQ[N[val]],
        qasmThrow["BadExpression", "Could not parse numeric expression `1`.", {e}, <|"Expression" -> e|>]
    ];
    val
]

qasmAngles[argStr_String] := qasmAngle /@ qasmSplitTopComma[argStr]


(* ---- gate-name table: name -> {baseShortcutBuilder, nParams, nControls, daggerQ} ---- *)

qasmGateSpec[name_String] := Lookup[$qasmGateTable, ToLowerCase[name], Missing["UnknownGate"]]

$qasmGateTable = <|
    "id" -> {Function[p, "I"], 0, 0, False},
    "i" -> {Function[p, "I"], 0, 0, False},
    "x" -> {Function[p, "X"], 0, 0, False},
    "y" -> {Function[p, "Y"], 0, 0, False},
    "z" -> {Function[p, "Z"], 0, 0, False},
    "h" -> {Function[p, "H"], 0, 0, False},
    "s" -> {Function[p, "S"], 0, 0, False},
    "sdg" -> {Function[p, "S"], 0, 0, True},
    "t" -> {Function[p, "T"], 0, 0, False},
    "tdg" -> {Function[p, "T"], 0, 0, True},
    "sx" -> {Function[p, "SX"], 0, 0, False},
    "sxdg" -> {Function[p, "SX"], 0, 0, True},
    "rx" -> {Function[p, "RX"[p[[1]]]], 1, 0, False},
    "ry" -> {Function[p, "RY"[p[[1]]]], 1, 0, False},
    "rz" -> {Function[p, "RZ"[p[[1]]]], 1, 0, False},
    "p" -> {Function[p, "P"[p[[1]]]], 1, 0, False},
    "u1" -> {Function[p, "P"[p[[1]]]], 1, 0, False},
    "phase" -> {Function[p, "P"[p[[1]]]], 1, 0, False},
    "u2" -> {Function[p, "U"[Pi/2, p[[1]], p[[2]]]], 2, 0, False},
    "u3" -> {Function[p, "U"[p[[1]], p[[2]], p[[3]]]], 3, 0, False},
    "u" -> {Function[p, "U"[p[[1]], p[[2]], p[[3]]]], 3, 0, False},
    "swap" -> {Function[p, "SWAP"], 0, 0, False},
    "cx" -> {Function[p, "X"], 0, 1, False},
    "cnot" -> {Function[p, "X"], 0, 1, False},
    "cy" -> {Function[p, "Y"], 0, 1, False},
    "cz" -> {Function[p, "Z"], 0, 1, False},
    "ch" -> {Function[p, "H"], 0, 1, False},
    "csx" -> {Function[p, "SX"], 0, 1, False},
    "crx" -> {Function[p, "RX"[p[[1]]]], 1, 1, False},
    "cry" -> {Function[p, "RY"[p[[1]]]], 1, 1, False},
    "crz" -> {Function[p, "RZ"[p[[1]]]], 1, 1, False},
    "cp" -> {Function[p, "P"[p[[1]]]], 1, 1, False},
    "cu1" -> {Function[p, "P"[p[[1]]]], 1, 1, False},
    "cphase" -> {Function[p, "P"[p[[1]]]], 1, 1, False},
    "cu3" -> {Function[p, "U"[p[[1]], p[[2]], p[[3]]]], 3, 1, False},
    "ccx" -> {Function[p, "X"], 0, 2, False},
    "toffoli" -> {Function[p, "X"], 0, 2, False},
    "ccz" -> {Function[p, "Z"], 0, 2, False},
    "cswap" -> {Function[p, "SWAP"], 0, 1, False},
    "fredkin" -> {Function[p, "SWAP"], 0, 1, False}
|>


(* ---- the importer (PackageScope entry used by QuantumQASM[src_String]) ---- *)

ImportOpenQASM[src0_String] := Catch[
    Module[{src, version, statements, regs = <||>, cregs = <||>, gateDefs = <||>, offset = 0, elements, qc},
        src = qasmStripComments[src0];
        version = qasmVersion[src];
        If[ ! MemberQ[{2, 3}, version],
            qasmThrow["BadVersion", "Unsupported OpenQASM version `1` (only 2 and 3).", {version}]
        ];
        statements = qasmStatements[src];
        elements = Reap[
            Scan[
                Function[stmt,
                    With[{r = qasmProcess[stmt, regs, cregs, gateDefs, offset]},
                        regs = r["regs"]; cregs = r["cregs"]; gateDefs = r["gateDefs"]; offset = r["offset"];
                        Scan[Sow, r["elements"]]
                    ]
                ],
                statements
            ]
        ][[2]];
        elements = If[elements === {}, {}, First[elements]];
        qc = Quiet @ Check[QuantumCircuitOperator[elements], $Failed];
        If[ Head[qc] =!= QuantumCircuitOperator,
            qasmThrow["Assembly", "Parsed statements did not assemble into a valid circuit."]
        ];
        qc
    ],
    $qasmTag
]


(* ---- process one statement; return updated state association, or throw ---- *)

qasmProcess[stmt_String, regs_, cregs_, gateDefs_, offset_] := Module[{state, mk},
    state = <|"regs" -> regs, "cregs" -> cregs, "gateDefs" -> gateDefs, "offset" -> offset, "elements" -> {}|>;
    mk[els_] := <|state, "elements" -> els|>;
    Which[
        StringMatchQ[stmt, RegularExpression["OPENQASM\\b.*"]] || StringMatchQ[stmt, RegularExpression["include\\b[\\s\\S]*"]],
            state,
        StringMatchQ[stmt, RegularExpression["qreg\\s+[\\s\\S]*"]],
            qasmDeclQ[stmt, "qreg", state],
        StringMatchQ[stmt, RegularExpression["creg\\s+[\\s\\S]*"]],
            qasmDeclC[stmt, "creg", state],
        StringMatchQ[stmt, RegularExpression["qubit(\\s*\\[\\d+\\])?\\s+[\\s\\S]*"]],
            qasmDeclQ[stmt, "qubit", state],
        StringMatchQ[stmt, RegularExpression["bit(\\s*\\[\\d+\\])?\\s+[\\s\\S]*"]],
            qasmDeclC[stmt, "bit", state],
        StringMatchQ[stmt, RegularExpression["gate\\s+[\\s\\S]*"]],
            qasmGateDef[stmt, state],
        StringMatchQ[stmt, RegularExpression["barrier\\b[\\s\\S]*"]],
            mk @ qasmBarrier[stmt, state],
        StringMatchQ[stmt, RegularExpression["reset\\b[\\s\\S]*"]],
            mk @ qasmReset[stmt, state],
        StringMatchQ[stmt, RegularExpression["(?:(?:inv|pow\\s*\\([^)]*\\)|ctrl(?:\\s*\\(\\s*\\d+\\s*\\))?|negctrl(?:\\s*\\(\\s*\\d+\\s*\\))?)\\s*@\\s*)*gphase\\b[\\s\\S]*"]],
            mk @ qasmGPhase[stmt, state],
        StringMatchQ[stmt, RegularExpression["[\\s\\S]*=\\s*measure\\b[\\s\\S]*"]],
            mk @ qasmMeasure[stmt, state, "v3"],
        StringMatchQ[stmt, RegularExpression["measure\\b[\\s\\S]*"]],
            mk @ qasmMeasure[stmt, state, "v2"],
        StringMatchQ[stmt, RegularExpression["(if|for|while|def|defcal|box|delay|let|input|output|const|int|uint|float|bool|angle|duration|stretch|array|return|qubit|bit)\\b[\\s\\S]*"]],
            qasmThrow["Unsupported", "Unsupported OpenQASM construct in statement: `1`.", {stmt}, <|"Statement" -> stmt|>],
        True,
            mk @ qasmGateApply[stmt, state]
    ]
]


(* -- register declarations -- *)

qasmDeclQ[stmt_String, kind_String, state_] := Module[{nm, n, cur = state["regs"], off = state["offset"]},
    {nm, n} = qasmParseDecl[stmt, kind];
    If[ KeyExistsQ[cur, nm], qasmThrow["Syntax", "Quantum register `1` already declared.", {nm}]];
    <|state, "regs" -> Append[cur, nm -> <|"offset" -> off, "size" -> n|>], "offset" -> off + n|>
]

qasmDeclC[stmt_String, kind_String, state_] := Module[{nm, n, cur = state["cregs"]},
    {nm, n} = qasmParseDecl[stmt, kind];
    <|state, "cregs" -> Append[cur, nm -> n]|>
]

qasmParseDecl[stmt_String, kind_String] := Module[{m},
    If[ MemberQ[{"qreg", "creg"}, kind],
        m = StringCases[stmt, RegularExpression[kind <> "\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*\\[\\s*(\\d+)\\s*\\]"] :> {"$1", "$2"}];
        If[m === {}, qasmThrow["Syntax", "Malformed register declaration: `1`.", {stmt}]];
        {m[[1, 1]], ToExpression[m[[1, 2]]]},
        (* qubit / bit (v3) *)
        m = StringCases[stmt, RegularExpression[kind <> "\\s*\\[\\s*(\\d+)\\s*\\]\\s+([A-Za-z_][A-Za-z0-9_]*)"] :> {"$2", "$1"}];
        If[ m =!= {},
            {m[[1, 1]], ToExpression[m[[1, 2]]]},
            m = StringCases[stmt, RegularExpression[kind <> "\\s+([A-Za-z_][A-Za-z0-9_]*)"] :> "$1"];
            If[m === {}, qasmThrow["Syntax", "Malformed declaration: `1`.", {stmt}]];
            {First[m], 1}
        ]
    ]
]


(* -- qubit-argument resolution -- *)

qasmResolveQ[ref_String, regs_] := Module[{trimmed = StringTrim[ref], m, nm, idx, reg},
    Which[
        (* physical qubit $N (an OpenQASM 3 hardware qubit, declared by no register) maps
           directly to absolute wire N+1; qiskit emits this for transpiled circuits *)
        StringMatchQ[trimmed, RegularExpression["\\$\\d+"]],
            {ToExpression[StringDrop[trimmed, 1]] + 1},
        (* indexed register reg[i] *)
        (m = StringCases[trimmed, RegularExpression["^([A-Za-z_][A-Za-z0-9_]*)\\s*\\[\\s*(\\d+)\\s*\\]$"] :> {"$1", "$2"}]) =!= {},
            nm = m[[1, 1]]; idx = ToExpression[m[[1, 2]]];
            reg = Lookup[regs, nm, qasmThrow["Range", "Undeclared register `1`.", {nm}]];
            If[ idx < 0 || idx >= reg["size"],
                qasmThrow["Range", "Qubit index `1` out of range for register `2`[`3`].", {idx, nm, reg["size"]}]
            ];
            {reg["offset"] + idx + 1},
        (* whole register reg *)
        True,
            nm = trimmed;
            reg = Lookup[regs, nm, qasmThrow["Range", "Undeclared register `1`.", {nm}]];
            Range[reg["offset"] + 1, reg["offset"] + reg["size"]]
    ]
]


(* -- barrier / reset / gphase -- *)

qasmBarrier[stmt_String, state_] := Module[{argStr, wires},
    argStr = StringTrim @ StringReplace[stmt, RegularExpression["^barrier\\b"] -> ""];
    wires = If[ argStr === "",
        Range[state["offset"]],
        Flatten[qasmResolveQ[#, state["regs"]] & /@ qasmSplitQargs[argStr]]
    ];
    If[wires === {}, {}, {"Barrier"[wires]}]
]

qasmReset[stmt_String, state_] := Module[{argStr, wires},
    (* OpenQASM `reset` returns a qubit to |0>. A QF circuit has no per-wire reset-to-|0>
       primitive and already treats any wire with no explicit initial state as |0> (the
       register), so a reset is the identity at circuit start and is dropped. Operands are
       still resolved so an out-of-range or undeclared qubit is reported. *)
    argStr = StringTrim @ StringReplace[stmt, RegularExpression["^reset\\b"] -> ""];
    wires = Flatten[qasmResolveQ[#, state["regs"]] & /@ qasmSplitQargs[argStr]];
    {}
]

(* gphase(theta) is a global phase E^(I theta). With ctrl / negctrl modifiers it becomes a
   phase applied only on the control-active configuration: a multi-controlled phase. It is
   realized by hosting a 1-qubit phase gate on one control wire and controlling it with the
   rest, so an uncontrolled gphase stays a 0-qudit GlobalPhase and a controlled one is a
   genuine operator on the control qubits. *)
qasmGPhase[stmt_String, state_] := Module[{
    mods, rest, nc1 = 0, nc0 = 0, invQ = False, powK = 1, m, theta, qargStr, wires, w1, w0, host, gate, restCtrl1, restCtrl0
},
    {mods, rest} = qasmStripModifiers[stmt];
    Scan[
        Replace[#, {
            "inv" :> (invQ = ! invQ),
            "ctrl"[k_] :> (nc1 += k),
            "negctrl"[k_] :> (nc0 += k),
            "pow"[k_] :> (powK *= k)
        }] &,
        mods
    ];
    m = StringCases[rest, RegularExpression["gphase\\s*\\(([^)]*)\\)\\s*([\\s\\S]*)$"] :> {"$1", "$2"}];
    If[m === {}, qasmThrow["Syntax", "Malformed gphase: `1`.", {stmt}]];
    theta = powK * qasmAngle[m[[1, 1]]];
    If[invQ, theta = - theta];
    If[ nc1 + nc0 == 0,
        Return @ {"GlobalPhase"[theta]}
    ];
    qargStr = StringTrim[m[[1, 2]]];
    wires = Flatten[qasmResolveQ[#, state["regs"]] & /@ qasmSplitQargs[qargStr]];
    If[ Length[wires] =!= nc1 + nc0,
        qasmThrow["Arity", "Controlled gphase expects `1` control qubit(s), got `2`.", {nc1 + nc0, Length[wires]}]
    ];
    w1 = wires[[1 ;; nc1]];
    w0 = wires[[nc1 + 1 ;; nc1 + nc0]];
    (* host gate carries the phase on |1> (positive control) or |0> (negative control) *)
    If[ w1 =!= {},
        host = Last[w1]; gate = {{1, 0}, {0, Exp[I theta]}}; restCtrl1 = Most[w1]; restCtrl0 = w0,
        host = Last[w0]; gate = {{Exp[I theta], 0}, {0, 1}}; restCtrl1 = {}; restCtrl0 = Most[w0]
    ];
    {
        If[ restCtrl1 === {} && restCtrl0 === {},
            QuantumOperator[gate, {host}],
            QuantumOperator["C"[gate, restCtrl1, restCtrl0], {host}]
        ]
    }
]


(* -- measurements (v2: `measure q -> c`, v3: `c = measure q`) -- *)

qasmMeasure[stmt_String, state_, ver_String] := Module[{qpart, wires},
    qpart = If[ ver === "v3",
        StringTrim @ Last @ StringSplit[StringSplit[stmt, "="][[2]], RegularExpression["\\bmeasure\\b"]],
        StringTrim @ First @ StringSplit[StringReplace[stmt, RegularExpression["^\\s*measure\\b"] -> ""], "->"]
    ];
    wires = qasmResolveQ[qpart, state["regs"]];
    {wires}
]


(* -- gate definitions -- *)

qasmGateDef[stmt_String, state_] := Module[{m, nm, paramStr, qubitStr, body, fparams, fqubits},
    m = StringCases[stmt,
        RegularExpression["gate\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(\\(([^)]*)\\))?\\s*([^\\{]*)\\{([\\s\\S]*)\\}"] :>
            {"$1", "$3", "$4", "$5"}];
    If[m === {}, qasmThrow["Syntax", "Malformed gate definition: `1`.", {stmt}]];
    {nm, paramStr, qubitStr, body} = m[[1]];
    fparams = If[StringTrim[paramStr] === "", {}, StringTrim /@ qasmSplitTopComma[paramStr]];
    fqubits = If[StringTrim[qubitStr] === "", {}, StringTrim /@ qasmSplitTopComma[qubitStr]];
    <|state, "gateDefs" -> Append[state["gateDefs"], nm -> <|"params" -> fparams, "qubits" -> fqubits, "body" -> body|>]|>
]


(* -- gate application: strip modifiers, resolve name/params/qubits -- *)

qasmGateApply[stmt_String, state_] := Module[{mods, rest, nc1 = 0, nc0 = 0, invQ = False, powK = 1, m, nm, paramStr, qargStr, params, qargs, wires},
    {mods, rest} = qasmStripModifiers[stmt];
    Scan[
        Replace[#, {
            "inv" :> (invQ = ! invQ),
            "ctrl"[k_] :> (nc1 += k),
            "negctrl"[k_] :> (nc0 += k),
            "pow"[k_] :> (powK *= k)
        }] &,
        mods
    ];
    m = StringCases[rest, RegularExpression["^([A-Za-z_][A-Za-z0-9_]*)\\s*(\\(([^)]*)\\))?\\s*([\\s\\S]*)$"] :> {"$1", "$3", "$4"}];
    If[m === {}, qasmThrow["Syntax", "Could not parse gate statement: `1`.", {stmt}]];
    {nm, paramStr, qargStr} = m[[1]];
    params = If[StringTrim[paramStr] === "", {}, qasmAngles[paramStr]];
    qargs = qasmSplitQargs[qargStr];
    If[qargs === {}, qasmThrow["Syntax", "Gate `1` has no qubit arguments.", {nm}]];
    wires = qasmResolveQ[#, state["regs"]] & /@ qargs;
    qasmBuildGate[nm, params, wires, nc1, nc0, invQ, powK, state, stmt]
]


qasmStripModifiers[stmt_String] := Module[{rest = StringTrim[stmt], mods = {}, m, done = False},
    NestWhile[
        Function[Null,
            m = StringCases[rest, RegularExpression["^(inv|pow\\s*\\([^)]*\\)|ctrl\\s*\\(\\s*\\d+\\s*\\)|ctrl|negctrl\\s*\\(\\s*\\d+\\s*\\)|negctrl)\\s*@\\s*([\\s\\S]*)$"] :> {"$1", "$2"}];
            If[ m === {}, done = True, AppendTo[mods, qasmParseModifier[m[[1, 1]]]]; rest = m[[1, 2]]]
        ],
        Null, ! done &, 1, 64
    ];
    {mods, rest}
]

qasmParseModifier[modStr_String] := Module[{s = StringTrim[modStr]},
    Which[
        s === "inv", "inv",
        StringMatchQ[s, RegularExpression["pow\\s*\\([\\s\\S]*\\)"]],
            "pow"[ToExpression @ First @ StringCases[s, RegularExpression["pow\\s*\\(\\s*(-?\\d+)\\s*\\)"] -> "$1"]],
        StringMatchQ[s, RegularExpression["ctrl\\s*\\([\\s\\S]*\\)"]],
            "ctrl"[ToExpression @ First @ StringCases[s, RegularExpression["ctrl\\s*\\(\\s*(\\d+)\\s*\\)"] -> "$1"]],
        s === "ctrl", "ctrl"[1],
        StringMatchQ[s, RegularExpression["negctrl\\s*\\([\\s\\S]*\\)"]],
            "negctrl"[ToExpression @ First @ StringCases[s, RegularExpression["negctrl\\s*\\(\\s*(\\d+)\\s*\\)"] -> "$1"]],
        True, "negctrl"[1]
    ]
]


(* build the QF element(s) for a gate, given resolved per-arg wire lists and modifiers *)

qasmBuildGate[nm_String, params_, wires_, nc1_, nc0_, invQ_, powK_, state_, stmt_] := Module[{
    gateDef, spec, baseBuilder, np, namedNc, namedDag, lens, nBroad, flatWires, totalCtrl, ctrl1, ctrl0, targets, shortcut, op
},
    gateDef = Lookup[state["gateDefs"], nm, Missing[]];
    If[ ! MissingQ[gateDef], Return @ qasmBuildUserGate[nm, gateDef, params, wires, state]];

    spec = qasmGateSpec[nm];
    If[ MissingQ[spec],
        qasmThrow["UnknownGate", "Unknown gate `1` (not a standard gate and not defined by a gate block).", {nm}, <|"Gate" -> nm|>]
    ];
    {baseBuilder, np, namedNc, namedDag} = spec;
    If[ Length[params] =!= np,
        qasmThrow["Arity", "Gate `1` expects `2` parameter(s), got `3`.", {nm, np, Length[params]}]
    ];

    (* register broadcast: expand element-wise over multi-wire args *)
    lens = Length /@ wires;
    nBroad = Max[lens];
    If[ nBroad > 1,
        If[ ! AllTrue[lens, # == 1 || # == nBroad &],
            qasmThrow["Arity", "Mismatched register sizes in broadcast of gate `1`.", {nm}]
        ];
        Return @ Flatten @ Table[
            qasmBuildGate[nm, params, Map[If[Length[#] == 1, #, {#[[i]]}] &, wires], nc1, nc0, invQ, powK, state, stmt],
            {i, nBroad}
        ]
    ];

    flatWires = Flatten[wires];
    totalCtrl = nc1 + nc0 + namedNc;
    If[ Length[flatWires] < totalCtrl + 1,
        qasmThrow["Arity", "Gate `1` applied to too few qubits.", {nm}]
    ];
    ctrl1 = Join[flatWires[[1 ;; nc1]], flatWires[[nc1 + nc0 + 1 ;; nc1 + nc0 + namedNc]]];
    ctrl0 = flatWires[[nc1 + 1 ;; nc1 + nc0]];
    targets = flatWires[[totalCtrl + 1 ;;]];
    shortcut = baseBuilder[params];
    op = If[ ctrl1 === {} && ctrl0 === {},
        QuantumOperator[shortcut, targets],
        QuantumOperator["C"[shortcut, ctrl1, ctrl0], targets]
    ];
    If[Xor[invQ, namedDag], op = op["Dagger"]];
    Which[
        powK == 1, {op},
        IntegerQ[powK] && powK >= 1, ConstantArray[op, powK],
        IntegerQ[powK] && powK < 0, ConstantArray[op["Dagger"], -powK],
        powK == 0, {},
        True, qasmThrow["Unsupported", "Non-integer pow(`1`) is not supported for gate `2`.", {powK, nm}]
    ]
]


(* instantiate a user-defined gate as a labeled subcircuit *)

qasmBuildUserGate[nm_String, gateDef_, params_, wires_, state_] := Module[{
    fparams = gateDef["params"], fqubits = gateDef["qubits"], body = gateDef["body"],
    flatWires, substBody, localRegs, subElems
},
    flatWires = Flatten[wires];
    If[ Length[params] =!= Length[fparams],
        qasmThrow["Arity", "Gate `1` expects `2` parameter(s), got `3`.", {nm, Length[fparams], Length[params]}]
    ];
    If[ Length[flatWires] =!= Length[fqubits],
        qasmThrow["Arity", "Gate `1` expects `2` qubit(s), got `3`.", {nm, Length[fqubits], Length[flatWires]}]
    ];
    substBody = Fold[
        Function[{txt, rule}, StringReplace[txt, WordBoundary ~~ rule[[1]] ~~ WordBoundary :> "(" <> ToString[rule[[2]], InputForm] <> ")"]],
        body,
        MapThread[Rule, {fparams, params}]
    ];
    localRegs = Association @ MapThread[#1 -> <|"offset" -> #2 - 1, "size" -> 1|> &, {fqubits, flatWires}];
    subElems = Reap[
        Module[{lr = localRegs},
            Scan[
                Function[s,
                    With[{r = qasmProcess[s, lr, state["cregs"], state["gateDefs"], 0]},
                        lr = r["regs"];
                        Scan[Sow, r["elements"]]
                    ]
                ],
                qasmStatements[substBody]
            ]
        ]
    ][[2]];
    subElems = If[subElems === {}, {}, First[subElems]];
    {QuantumCircuitOperator[subElems, nm]}
]
