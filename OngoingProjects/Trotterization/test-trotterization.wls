#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Test the standalone Trotterization implementation against QF.

   T1. Per-gate equality:    standalone gate k == QF gate k  (max |diff|)
   T2. Total unitary match:  Dot @@ Reverse[standalone gates] == QF circuit matrix
   T3. Coefficient parity:   our trotterCoeffs == QF's
   T4. Convergence sanity:   Trotter error -> 0 as reps grows / order rises
   T5. Form options:         "Flat", "Steps", "Unitary" return consistent shapes
   T6. Notebook replication: trace-distance table for the user's 3q Ising
                             Hamiltonian over (order, steps) = {1,2,4} x 1..10
   T7. Notebook replication: spectral 2-norm table for same grid
   T8. 4-qubit transverse Ising: agreement on a larger system
   T9. Symbolic time:        gates with symbolic t match QF symbolically
   T10. Notebook gate count: c["Gates"] = reps * Length[trotterExpand[ops, order]]

   Tolerance for floating-point comparisons: 1*^-12.
*)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Trotterization/Trotterization.wl"];

tol = 1.*^-12;

(* Render a list-of-rows as a plain text table for wolframscript stdout. *)
ClearAll[textTable, fmt];
textTable[rows_List] := Module[{strRows, widths},
    strRows = Map[fmt, rows, {2}];
    widths = Max /@ Transpose[Map[StringLength, strRows, {2}]];
    StringJoin @ Riffle[
        StringJoin @ Riffle[MapThread[StringPadRight[#1, #2 + 2] &, {#, widths}], ""] & /@ strRows,
        "\n"
    ]
];

(* Format any value as a single-line string. *)
fmt[x_Real] /; x === 0. := "0";
fmt[x_Real] := Module[{e, m, mantissa},
    e = Floor[Log10[Abs[x]]];
    If[-3 <= e <= 3,
        ToString[N[Round[x, 10.^(e - 2)]]],
        m = x / 10.^e;
        mantissa = ToString[N[Round[m, 0.01]]];
        StringJoin[mantissa, "e", If[e >= 0, "+", ""], ToString[e]]
    ]
];
fmt[x_String] := x;
fmt[x_Rational] := ToString[x, InputForm];
fmt[x_] := ToString[x, InputForm];

(* ==================================================================== *)
(* Helper: build (QF ops, plain-matrix ops) for a Hamiltonian spec.     *)
(* ==================================================================== *)

ClearAll[buildOps];
buildOps[hSpec_] := Module[{h, qfOps, matOps},
    h = N @ hSpec;
    qfOps = KeyValueMap[QuantumOperator[#2 #1] &, h["PauliDecompose"]];
    matOps = #["Matrix"] & /@ qfOps;
    matOps = Normal /@ matOps;  (* dense *)
    {h, qfOps, matOps}
];

(* QF's gate matrices in circuit order. For reps > 1, QF wraps each Trotter
   step as a sub-QuantumCircuitOperator; ["Flatten"]["Operators"] unwraps. *)
(* NOTE on syntax: QF's dispatch rule is
       QuantumCircuitOperator["Trotterization"[opArgs_, args___], opts___]
   so the named circuit must be passed as the curried head form
       QuantumCircuitOperator["Trotterization"[ops, order, reps, t]]
   Wrapping the arguments in a List instead -- e.g.
       QuantumCircuitOperator[{"Trotterization", ops, order, reps, t}]  (* WRONG *)
   does NOT call the named circuit; QuantumCircuitOperator interprets the list
   as a sequence of gate elements ("Trotterization" alone triggers the default
   3q X/Y/Z Trotterization, the operator list becomes a sub-circuit, and the
   numeric arguments are silently turned into PhaseShift/OverHat gates). The
   user-facing notebook at 411282707.nb has the same bug and produces garbage
   unitaries for that reason. *)

qfGateMatrices[qfOps_, order_, reps_, t_] := Module[{c},
    c = QuantumCircuitOperator["Trotterization"[qfOps, order, reps, t]];
    Normal[#["Matrix"]] & /@ c["Flatten"]["Operators"]
];

qfCircuitMatrix[qfOps_, order_, reps_, t_] := Module[{c},
    c = QuantumCircuitOperator["Trotterization"[qfOps, order, reps, t]];
    Normal @ c["Matrix"]
];


(* ==================================================================== *)
(* Test runner                                                          *)
(* ==================================================================== *)

results = {};
failures = {};

runCase[name_, hSpec_, order_, reps_, t_] := Module[{
    h, qfOps, matOps,
    qfGates, sGates,
    qfU, sU, exactU,
    gateDiff, unitaryDiff, exactDiff
},
    {h, qfOps, matOps} = buildOps[hSpec];

    qfGates = qfGateMatrices[qfOps, order, reps, t];
    sGates = TrotterizationMatrices[matOps, order, reps, t];

    If[Length[qfGates] =!= Length[sGates],
        AppendTo[failures, {name, "gate count mismatch", Length[qfGates], Length[sGates]}];
        Return[$Failed]
    ];

    gateDiff = Max[Abs[Flatten[qfGates - sGates]]];

    qfU = qfCircuitMatrix[qfOps, order, reps, t];
    sU  = TrotterizationMatrices[matOps, order, reps, t, "Form" -> "Unitary"];
    unitaryDiff = Max[Abs[Flatten[qfU - sU]]];

    exactU = MatrixExp[-I t Normal[h["Matrix"]]];
    exactDiff = Max[Abs[Flatten[sU - exactU]]];

    AppendTo[results, <|
        "name" -> name, "order" -> order, "reps" -> reps, "t" -> t,
        "nGates" -> Length[sGates],
        "gateDiff" -> gateDiff,
        "unitaryDiff" -> unitaryDiff,
        "trotterError" -> exactDiff
    |>];

    If[gateDiff > tol,
        AppendTo[failures, {name, "T1 gate match", order, reps, t, gateDiff}]];
    If[unitaryDiff > tol,
        AppendTo[failures, {name, "T2 unitary match", order, reps, t, unitaryDiff}]];
];


(* ==================================================================== *)
(* Test cases                                                           *)
(* ==================================================================== *)

Print["=== Tests T1, T2: standalone vs QF, gate-by-gate and full unitary ==="];

(* Case A: 1-qubit, H = X + Z *)
hA = QuantumOperator["X"] + QuantumOperator["Z"];
runCase["1q X+Z", hA, 1, 1, 0.3];
runCase["1q X+Z", hA, 1, 4, 0.3];
runCase["1q X+Z", hA, 2, 1, 0.3];
runCase["1q X+Z", hA, 2, 4, 0.3];
runCase["1q X+Z", hA, 4, 1, 0.3];
runCase["1q X+Z", hA, 4, 3, 0.3];

(* Case B: 2-qubit, H = X1 + X2 + Z1Z2 *)
hB = QuantumOperator["X" -> {1}] + QuantumOperator["X" -> {2}] + QuantumOperator["ZZ"];
runCase["2q TFIM", hB, 1, 1, 0.5];
runCase["2q TFIM", hB, 2, 2, 0.5];
runCase["2q TFIM", hB, 4, 1, 0.5];

(* Case C: 3-qubit Ising-style (the user's notebook example).
   The string-sum-inside-QuantumOperator form preserves all five Pauli terms. *)
hC = QuantumOperator["X" + ("X" -> {2}) + ("X" -> {3}) - "ZZ" - ("ZZ" -> {2, 3})];
runCase["3q Ising-like", hC, 1, 1, 1.0];
runCase["3q Ising-like", hC, 1, 5, 1.0];
runCase["3q Ising-like", hC, 2, 1, 1.0];
runCase["3q Ising-like", hC, 2, 4, 1.0];
runCase["3q Ising-like", hC, 4, 1, 1.0];
runCase["3q Ising-like", hC, 4, 2, 1.0];

(* Print results table *)
Print[];
Print[textTable[
    Prepend[
        {#name, #order, #reps, #t, #nGates, #gateDiff, #unitaryDiff, #trotterError} & /@ results,
        {"case", "order", "reps", "t", "nGates", "gateDiff", "unitaryDiff", "trotterError"}
    ]
]];


(* ==================================================================== *)
(* T3: coefficient list parity with QF's internal helper                *)
(* ==================================================================== *)

Print[];
Print["=== T3: coefficient parity with QF's trotterCoeffs ==="];

(* Reach into QF's per-file PackagePrivate context to compare directly. *)
qfCoeffs = Wolfram`QuantumFramework`NamedCircuits`PackagePrivate`trotterCoeffs;

coeffCases = {{3, 1, 1}, {3, 2, 1}, {3, 4, 1}, {5, 2, 1/2}, {5, 4, 1/3}, {2, 6, 1}};
coeffParity = Table[
    With[{ours = TrotterizationCoefficients @@ c, theirs = qfCoeffs @@ c},
        <|"args" -> c, "match" -> ours === theirs, "ours" -> ours, "theirs" -> theirs|>
    ],
    {c, coeffCases}
];
Print[textTable[
    Prepend[
        {ToString[#args, InputForm], #match, Length[#ours], Length[#theirs]} & /@ coeffParity,
        {"args (l,n,c)", "match?", "len(ours)", "len(theirs)"}
    ]
]];
If[!AllTrue[coeffParity, #match &],
    AppendTo[failures, {"T3 coeff parity", Cases[coeffParity, x_ /; !x["match"] :> x["args"]]}]];


(* ==================================================================== *)
(* T4: convergence to Exp[-I t H] as reps -> infinity, order rises      *)
(* ==================================================================== *)

Print[];
Print["=== T4: convergence sanity (Trotter error vs reps for fixed t = 1) ==="];

{hConv, qfOpsConv, matOpsConv} = buildOps[hC];
exactConv = MatrixExp[-I 1.0 Normal[hConv["Matrix"]]];

convTable = Table[
    Module[{u},
        u = TrotterizationMatrices[matOpsConv, order, reps, 1.0, "Form" -> "Unitary"];
        {order, reps, Norm[u - exactConv, 2]}
    ],
    {order, {1, 2, 4}}, {reps, {1, 2, 4, 8, 16}}
];
Print[textTable[
    Prepend[
        Map[{#[[1]], #[[2]], #[[3]]} &, Catenate[convTable]],
        {"order", "reps", "spectral error"}
    ]
]];

(* Convergence check: at each order, error at reps=16 must be substantially
   smaller than at reps=1. We require >= 5x reduction (slope expected: order=1
   gives 16x, order=2 gives 256x; high-order with small step can saturate near
   machine epsilon and grow due to floating-point noise, which is expected). *)
convergence = Table[
    Module[{eFirst = row[[1, 3]], eLast = row[[-1, 3]]},
        {row[[1, 1]], eFirst, eLast, eFirst / Max[eLast, $MachineEpsilon]}
    ],
    {row, convTable}
];
Print[];
Print["Reduction factors (reps=1 -> reps=16):"];
Print[textTable[
    Prepend[
        {#[[1]], #[[2]], #[[3]], #[[4]]} & /@ convergence,
        {"order", "err(reps=1)", "err(reps=16)", "ratio"}
    ]
]];
convergedOK = AllTrue[convergence, (#[[2]] / Max[#[[3]], $MachineEpsilon]) >= 5 || #[[2]] < 1.*^-12 &];
Print["Converges (>=5x reduction or already at noise floor): ", convergedOK];
If[! convergedOK, AppendTo[failures, {"T4 convergence", convergence}]];


(* ==================================================================== *)
(* T5: form options shape sanity                                        *)
(* ==================================================================== *)

Print[];
Print["=== T5: form-option shape checks ==="];

matOps5 = buildOps[hC][[3]];
flat   = TrotterizationMatrices[matOps5, 2, 3, 1.0, "Form" -> "Flat"];
steps  = TrotterizationMatrices[matOps5, 2, 3, 1.0, "Form" -> "Steps"];
uMat   = TrotterizationMatrices[matOps5, 2, 3, 1.0, "Form" -> "Unitary"];

t5Checks = {
    {"flat is list of matrices", AllTrue[flat, MatrixQ]},
    {"steps has length reps=3", Length[steps] === 3},
    {"flat == Catenate[steps]", flat === Catenate[steps]},
    {"unitary equals reverse-Dot of flat",
        Max[Abs[Flatten[uMat - Dot @@ Reverse[flat]]]] < tol}
};
Print[textTable[Prepend[t5Checks, {"check", "pass?"}]]];
Scan[If[! #[[2]], AppendTo[failures, {"T5", #[[1]]}]] &, t5Checks];


(* ==================================================================== *)
(* T6, T7: replicate the user's notebook                                *)
(*    h = QuantumOperator["X"+("X"->{2})+("X"->{3})-"ZZ"-("ZZ"->{2,3})] *)
(*    ops = KeyValueMap[QuantumOperator[#2 #1]&, h["PauliDecompose"]]   *)
(*    For order in {1,2,4} and steps in 1..10, build                    *)
(*       cQF       = QuantumCircuitOperator["Trotterization"[N@ops, order, steps, 1]] *)
(*       cStandalone = standalone unitary from N matrices               *)
(*    (Note: the notebook itself uses the LIST form which is buggy;     *)
(*    we use the curried form here to actually exercise the algorithm.) *)
(*    Compare unitaries (T6a), trace distances to MatrixExp[-I H] (T6), *)
(*    and spectral 2-norm distances (T7).                               *)
(* ==================================================================== *)

Print[];
Print["=== T6, T7: notebook replication (order in {1,2,4}, steps 1..10, t=1) ==="];

{hNb, qfOpsNb, matOpsNb} = buildOps[hC];
exactNb = N @ MatrixExp[-I 1.0 Normal[hNb["Matrix"]]];
exactState = QuantumState[exactNb];  (* for QuantumDistance comparison *)

orderList = {1, 2, 4};
stepList  = Range[10];

nbResults = Reap[
    Do[
        Module[{cQF, qfMat, sMat, traceDist, twoNorm, gateCnt, expectedCnt},
            cQF = N @ QuantumCircuitOperator["Trotterization"[N @ qfOpsNb, order, steps, 1]];
            qfMat = Normal @ cQF["Matrix"];
            sMat  = TrotterizationMatrices[matOpsNb, order, steps, 1.0, "Form" -> "Unitary"];

            (* T6a: standalone vs QF unitary equality (machine eps) *)
            Sow[<|
                "kind" -> "match", "order" -> order, "steps" -> steps,
                "diff" -> Max[Abs[Flatten[qfMat - sMat]]]
            |>];

            (* T6: trace distance from notebook (uses standalone matrix) *)
            traceDist = QuantumDistance[QuantumState[sMat], exactState, "Trace"];
            twoNorm   = Norm[sMat - exactNb, 2];
            gateCnt   = cQF["Gates"];
            expectedCnt = steps * Length[TrotterizationOrdering[matOpsNb, order]];
            Sow[<|
                "kind" -> "scaling", "order" -> order, "steps" -> steps,
                "gates" -> gateCnt, "expectedGates" -> expectedCnt,
                "trace" -> traceDist, "2norm" -> twoNorm
            |>];
        ],
        {order, orderList}, {steps, stepList}
    ]
][[2, 1]];

matchRows  = Cases[nbResults, KeyValuePattern["kind" -> "match"]];
scaleRows  = Cases[nbResults, KeyValuePattern["kind" -> "scaling"]];

(* T6a: print maximum standalone-vs-QF gap per order *)
Print["Max |standalone - QF| per order, across all steps:"];
Print[textTable[
    Prepend[
        Function[order,
            {order, Max[#diff & /@ Cases[matchRows, KeyValuePattern["order" -> order]]]}
        ] /@ orderList,
        {"order", "max gap"}
    ]
]];
nbMatchFailures = Cases[matchRows, x_ /; x["diff"] > tol];
If[nbMatchFailures =!= {},
    AppendTo[failures, {"T6a notebook unitary mismatch", Length[nbMatchFailures]}]];

(* T6/T7: trace distance and 2-norm tables.
   Pivot to (order x steps) grid layout matching the notebook. *)
buildGrid[key_] := Module[{header, body},
    header = Prepend[ToString /@ stepList, "order \\ steps"];
    body = Function[order,
        Prepend[
            #[key] & /@ Cases[scaleRows, KeyValuePattern["order" -> order]],
            order
        ]
    ] /@ orderList;
    Prepend[body, header]
];

Print[];
Print["T6: Trace distance  D(QuantumState[trotter], QuantumState[exact])"];
Print[textTable[buildGrid["trace"]]];

Print[];
Print["T7: Spectral 2-norm  ||trotter - exact||_2"];
Print[textTable[buildGrid["2norm"]]];

(* T7 sanity: at each order, the 2-norm error should monotonically decrease
   with steps (the notebook's plot is on log-log axes and shows a clean line). *)
nbMonotone = AllTrue[
    orderList,
    Function[order,
        With[{errs = #["2norm"] & /@ Cases[scaleRows, KeyValuePattern["order" -> order]]},
            OrderedQ[Reverse[errs]] || Last[errs] < 1.*^-12
        ]
    ]
];
Print["Monotone 2-norm decrease with steps (per order): ", nbMonotone];
If[! nbMonotone, AppendTo[failures, {"T7 monotonicity violated"}]];

(* T10: gate count formula  c["Gates"] = steps * Length[trotterExpand[ops, order]] *)
Print[];
Print["T10: gate count = steps * Length[trotterExpand[ops, order]]"];
Print[textTable[buildGrid["gates"]]];
gateCntCheck = AllTrue[scaleRows, #["gates"] === #["expectedGates"] &];
Print["All match expected formula: ", gateCntCheck];
If[! gateCntCheck, AppendTo[failures, {"T10 gate count mismatch"}]];


(* ==================================================================== *)
(* T8: 4-qubit transverse Ising (a more demanding system)               *)
(*    H = -J Sum Z_i Z_{i+1} - g Sum X_i, J=1, g=0.7, t=0.4             *)
(* ==================================================================== *)

Print[];
Print["=== T8: 4-qubit transverse Ising H ==="];

g = 0.7;
(* Use the user's notebook string-sum form so PauliDecompose returns all 7 terms.
   The arithmetic-on-QuantumOperator form silently drops the X qubit-permuted
   pieces in some versions of QF. *)
hT = QuantumOperator[
        -("ZZ" -> {1, 2}) - ("ZZ" -> {2, 3}) - ("ZZ" -> {3, 4})
        - g ("X" -> {1}) - g ("X" -> {2}) - g ("X" -> {3}) - g ("X" -> {4})
    ];

{hTI, qfOpsTI, matOpsTI} = buildOps[hT];

Print["nTerms = ", Length[matOpsTI], ", dim = ", Length[matOpsTI[[1]]]];

t8Cases = {
    {"4q TI", 1, 1, 0.4}, {"4q TI", 1, 8, 0.4},
    {"4q TI", 2, 1, 0.4}, {"4q TI", 2, 4, 0.4},
    {"4q TI", 4, 1, 0.4}, {"4q TI", 4, 2, 0.4}
};
Do[
    runCase[c[[1]] <> "(" <> ToString[Length[matOpsTI]] <> "t)", hT, c[[2]], c[[3]], c[[4]]],
    {c, t8Cases}
];

(* Print only the new T8 rows (the last 6 added to results) *)
t8Slice = Take[results, -Length[t8Cases]];
Print[textTable[
    Prepend[
        {#name, #order, #reps, #t, #nGates, #gateDiff, #unitaryDiff, #trotterError} & /@ t8Slice,
        {"case", "order", "reps", "t", "nGates", "gateDiff", "unitaryDiff", "trotterError"}
    ]
]];


(* ==================================================================== *)
(* T9: symbolic time -- gates contain unbound parameter t, must agree   *)
(* ==================================================================== *)

Print[];
Print["=== T9: symbolic time parameter ==="];

(* Use the 1-qubit X+Z system; a small system keeps the symbolic
   MatrixExp tractable without an explicit closed form. *)
qfOpsSym = buildOps[hA][[2]];
matOpsSym = buildOps[hA][[3]];
sym = Symbol["tSym"];  (* avoid colliding with global t *)

(* QF circuit at symbolic time *)
cSymQF = QuantumCircuitOperator["Trotterization"[qfOpsSym, 2, 2, sym]];
gatesQFSym = Normal[#["Matrix"]] & /@ cSymQF["Flatten"]["Operators"];
gatesStdSym = TrotterizationMatrices[matOpsSym, 2, 2, sym];

(* Substitute several real values and check numeric equality. *)
testTimes = {0.1, 0.3, 0.7, 1.5};
symResults = Table[
    Module[{qfNum, sNum},
        qfNum = gatesQFSym /. sym -> tval;
        sNum  = gatesStdSym /. sym -> tval;
        {tval, Length[qfNum], Length[sNum],
         Max[Abs[Flatten[N[qfNum] - N[sNum]]]]}
    ],
    {tval, testTimes}
];
Print[textTable[
    Prepend[symResults, {"t", "len(QF gates)", "len(std gates)", "max gap"}]
]];
symFails = Cases[symResults, {_, _, _, d_} /; d > tol];
If[symFails =!= {}, AppendTo[failures, {"T9 symbolic time gap", symFails}]];


(* ==================================================================== *)
(* Final report                                                         *)
(* ==================================================================== *)

Print[];
Print[StringRepeat["=", 60]];
If[failures === {},
    Print["ALL TESTS PASSED (", Length[results], " gate/unitary cases + T3..T10)"],
    Print["FAILURES: ", Length[failures]];
    Scan[Print["  ", #] &, failures]
];
Print[StringRepeat["=", 60]];

Quit[];
