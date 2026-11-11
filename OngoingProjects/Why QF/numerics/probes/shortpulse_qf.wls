(* Short-pulse head-to-head, QuantumFramework side. Every QuantumEvolve option is left at its
   default (no MaxStepSize, no Method). The cases and their exact values are those of
   shortpulse_qutip.py; each result is scored against its exact value, and the wall time is the
   median of three evaluations. Writes shortpulse_qf.tsv next to this file.

   Besides the default row, each case gets a control row: NDSolveValue at its defaults on the
   equations QuantumEvolve builds ("ReturnEquations" -> True), with QuantumEvolve's restart events
   removed. QuantumEvolve passes NDSolveValue no options of its own, so the control is QF without
   the pulse landmarks, and the two rows together say which of them captures a pulse. Rows marked
   "diagnostic" set one option away from its default to test what captures a pulse that has no
   landmarks.

   Run: wolframscript -file shortpulse_qf.wls *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];

repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
here = DirectoryName[$InputFileName];
out = OpenWrite[FileNameJoin[{here, "shortpulse_qf.tsv"}]];

versions = StringJoin[
    "# Wolfram Language ", $Version,
    ", QuantumFramework ", First[PacletFind["Wolfram/QuantumFramework"]]["Version"],
    " at ", StringTrim @ RunProcess[{"git", "-C", repo, "rev-parse", "--short", "HEAD"}, "StandardOutput"],
    ", uncommitted kernel files: {", StringRiffle[StringSplit @ StringTrim @ RunProcess[{"git", "-C", repo, "status", "--porcelain", "--", "QuantumFramework/Kernel"}, "StandardOutput"], " "], "}",
    ", Arrays ", First[PacletFind["Wolfram/Arrays"]]["Version"],
    ", TensorNetworks ", First[PacletFind["Wolfram/TensorNetworks"]]["Version"],
    ", ", DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute"}]
];

num[x_ ? NumericQ] := ToString[CForm[N[x]]]
num[x : {__ ? NumericQ}] := StringRiffle[num /@ x, ", "]
num[x_String] := x
num[x_] := ToString[x, InputForm]

emit[line_String] := (WriteLine[out, line]; WriteString["stdout", line, "\n"])

row[case_, variant_, setting_, observable_, value_, exact_, seconds_, note_] := emit @ StringRiffle[
    {case, variant, "QF", setting, observable, num[value], num[exact],
        If[NumericQ[value] && NumericQ[exact], num[Abs[value - exact]], ""], num[seconds], note},
    "\t"
]

(* the last value of three evaluations, their median wall time, and the names of the messages they issued *)
SetAttributes[run, HoldFirst];
run[expr_] := With[{data = Table[EvaluationData[expr], 3]},
    <|
        "Value" -> Last[data]["Result"],
        "Seconds" -> Median[#["AbsoluteTiming"] & /@ data],
        "Messages" -> StringRiffle[DeleteDuplicates[ToString /@ Catenate[#["Messages"] & /@ data]], "; "]
    |>
]

notes[items___] := StringRiffle[DeleteCases[{items}, ""], "; "]

(* probabilities of a state vector or of a density matrix's diagonal *)
probabilities[v_ ? VectorQ] := Abs[v] ^ 2
probabilities[m_ ? MatrixQ] := Re[Diagonal[m]]

(* NDSolveValue at its defaults on QuantumEvolve's equations without the restart events *)
withoutRestarts[h_, psi0_, {t0_, t1_}] := With[{eqs = QuantumEvolve[h, psi0, {\[FormalT], t0, t1}, "ReturnEquations" -> True]},
    NDSolveValue[
        DeleteCases[First[eqs], WhenEvent[_, "RestartIntegration"]] /. {sa_SparseArray ? SparseArrayQ :> Normal[sa], \[FormalS] -> ss},
        ss,
        Last[eqs]
    ]
]

(* One case: QF at its defaults, the control without restart events, and any diagnostic option
   sets. score maps a function from a time to the probabilities at that time onto the result. *)
caseRows[case_, variant_, observable_, exact_, {h_, psi0_, range_}, score_, diagnostics_ : {}] := (
    With[{r = run @ score[With[{psi = QuantumEvolve[h, psi0, Prepend[range, \[FormalT]]]}, psi[#]["ProbabilitiesList"] &]]},
        row[case, variant, "default", observable, r["Value"], exact, r["Seconds"],
            notes["restart times " <> ToString[Length[Cases[First @ QuantumEvolve[h, psi0, Prepend[range, \[FormalT]], "ReturnEquations" -> True], _WhenEvent]]], r["Messages"]]]
    ];
    With[{r = run @ score[With[{sol = withoutRestarts[h, psi0, range]}, probabilities[sol[#]] &]]},
        row[case, variant, "control: NDSolve defaults, restart events removed", observable, r["Value"], exact, r["Seconds"], r["Messages"]]
    ];
    Scan[
        With[{r = run @ score[With[{psi = QuantumEvolve[h, psi0, Prepend[range, \[FormalT]], #]}, psi[#]["ProbabilitiesList"] &]]},
            row[case, variant, "diagnostic: " <> ToString[#, InputForm], observable, r["Value"], exact, r["Seconds"], r["Messages"]]
        ] &,
        diagnostics
    ]
)

emit[versions];
emit["case\tvariant\tpackage\tsetting\tobservable\tresult\texact\tabs_error\tseconds\tnote"];

(* the first QuantumEvolve of a kernel loads code that later calls do not *)
warmup = QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 1}][1]["ProbabilitiesList"];


(* 1. QuTiP docs, "Working with pulses": Z + 10 Pi Boole[0.7 < t < 0.75] X from |1>. The pulse is
   a rotation about Z + a X, a = 10 Pi, for 0.05, so <n>(1) = Cos[w/20]^2 + Sin[w/20]^2 / w^2, w = Sqrt[1 + a^2]. *)
With[{w = Sqrt[1 + (10 Pi) ^ 2]},
    caseRows["1 docs", "pulse 0.7<t<0.75", "<n>(1)", Cos[w / 20] ^ 2 + Sin[w / 20] ^ 2 / w ^ 2,
        {QuantumOperator["Z"] + 10 Pi Boole[0.7 < \[FormalT] < 0.75] QuantumOperator["X"], QuantumState["1"], {0, 1}},
        #[1][[2]] &
    ]
];


(* 2a. Issue #1265: two cosine pi/2 pulses of width 20, the second starting `wait` after the first
   ends, evolved as a density matrix (QuTiP's mesolve with zero collapse operators). Pe(99.5) = 1. *)
piLen = 20;
piAmp = Pi 1.0 / 20;
cosinePulse[t_, t0_, amp_, w_ : piLen, phase_ : 0, df_ : 0] := With[{x = t - t0},
    amp 0.5 (1 + Cos[2 Pi x / w]) Boole[x + w / 2. > 0] Boole[- x + w / 2. > 0] Exp[I phase - 2 I Pi df x]
]
seq1265[wait_] := cosinePulse[\[FormalT], piLen / 2., 0.5 piAmp] + cosinePulse[\[FormalT], wait + 1.0 piLen + piLen / 2., 0.5 piAmp]

Scan[
    caseRows["2a #1265", "wait " <> ToString[#] <> Lookup[<|0 -> " (c_x)", 40 -> " (c_x_bug)", 23.0018 -> " (BoxiLi threshold)"|>, #, ""], "Pe(99.5)", 1,
        {seq1265[#] QuantumOperator["X"], QuantumState["0"]["MatrixState"], {0, 99.5}},
        #[99.5][[2]] &,
        If[# === 40, {MaxStepFraction -> 1}, {}]
    ] &,
    {0, 10, 20, 23, 23.0017, 23.0018, 24, 30, 40, 50}
];


(* 2b. Issue #1945: 2 Pi omegaR X on [100, 600], <Z> on Subdivide[0, 700, n - 1] (numpy's
   linspace(0, 700, n)), scored by the largest error over the grid against Cos[2 Omega (t - 100)]. *)
func1945 = Piecewise[{{0, \[FormalT] < 100}, {0, \[FormalT] > 600}}, 2 Pi 2*^-3];
err1945[n_] := With[{times = Subdivide[0., 700., n - 1]},
    Function[p, Max @ Abs[(p[#] . {1, -1} & /@ times) - Cos[2 (2 Pi 2*^-3) Clip[times - 100, {0, 500}]]]]
]
Scan[
    caseRows["2b #1945", "tlist linspace(0,700," <> ToString[#] <> ")", "max_t |<sz> - exact|", 0,
        {func1945 QuantumOperator["X"], QuantumState["0"], {0, 700}},
        err1945[#]
    ] &,
    {267, 268, 800}
];

(* the same solve, read through the state vector before its time is bound *)
With[{
    r = run @ With[{
        vec = Normal @ QuantumEvolve[func1945 QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 700}]["StateVector"]
    },
        err1945[268][Abs[vec /. \[FormalT] -> #] ^ 2 &]
    ]
},
    row["2b #1945", "tlist linspace(0,700,268)", "default, read via psi[\"StateVector\"] /. t -> time", "max_t |<sz> - exact|", r["Value"], 0, r["Seconds"], r["Messages"]]
];


(* 2d. Issue #2814 has no code; its sketch is a pulse train with long waits. Five Gaussian Pi/5
   pulses, sigma = 0.01, centred at 10, 30, 50, 70, 90: with k captured, P1 = Sin[k Pi/10]^2. *)
caseRows["2d #2814 train", "5 Gaussian pi/5, w 0.04", "P1(100)", 1,
    {(Total[(Pi / 5) / (Sqrt[2 Pi] 0.01) Exp[- (\[FormalT] - #) ^ 2 / (2 0.01 ^ 2)] & /@ {10, 30, 50, 70, 90}] / 2) QuantumOperator["X"], QuantumState["0"], {0, 100}},
    #[100][[2]] &
];


(* 3. A pi pulse at t = 50 on [0, 100] under H = Omega(t) X/2, four shapes of width w. P1(100) = 1. *)
shapes[w_] := With[{s = w / 4}, <|
    "square Piecewise" -> Piecewise[{{Pi / w, Abs[\[FormalT] - 50] < w / 2}}, 0],
    "square UnitStep" -> (Pi / w) UnitStep[\[FormalT] - (50 - w / 2)] UnitStep[(50 + w / 2) - \[FormalT]],
    "Gaussian sigma=w/4" -> (Pi / (Sqrt[2 Pi] s)) Exp[- (\[FormalT] - 50) ^ 2 / (2 s ^ 2)],
    "sech tau=w/4" -> (1 / s) Sech[(\[FormalT] - 50) / s]
|>]

Scan[
    Function[w, KeyValueMap[
        caseRows["3 pi pulse", #1 <> ", w " <> StringReplace[ToString[w], RegularExpression["\\.$"] -> ".0"], "P1(100)", 1, {(#2 / 2) QuantumOperator["X"], QuantumState["0"], {0, 100}}, #[100][[2]] &] &,
        shapes[w]
    ]],
    {1., 0.1, 0.01, 0.001}
];


(* 4. Three-level transmon, alpha/2pi = -200 MHz (time in ns), rotating frame,
   H = alpha/2 a+a+aa + Omega(t)/2 (a + a+), Gaussian pi pulse after an idle gap, no DRAG.
   "windowed": sigma = T/4, T = 16, truncated to [gap, gap + T] with area Pi over the window, t_end = gap + T.
   "untruncated": sigma = 4, centre gap + 48, t_end = gap + 96. The result is {P1, P2} at t_end; the
   exact values are in shortpulse_reference.tsv. *)
alpha = -2 Pi 0.2;
tT = 16.;
sig = 4.;
aWin = Pi / (sig Sqrt[2 Pi] Erf[tT / (2 Sqrt[2] sig)]);
bb = AnnihilationOperator[3];
bd = bb["Dagger"];
transmonH[om_] := (alpha / 2) (bd @ bd @ bb @ bb) + (om / 2) (bb + bd)
pulse4["windowed", gap_] := {Piecewise[{{aWin Exp[- (\[FormalT] - (gap + tT / 2)) ^ 2 / (2 sig ^ 2)], gap <= \[FormalT] <= gap + tT}}, 0], gap + tT}
pulse4["untruncated", gap_] := {(Pi / (Sqrt[2 Pi] sig)) Exp[- (\[FormalT] - (gap + 48)) ^ 2 / (2 sig ^ 2)], gap + 96}

Scan[
    Function[{kind}, Scan[
        Function[{gap}, With[{p = pulse4[kind, gap]},
            caseRows["4 transmon", kind <> ", gap " <> ToString[gap], "P1, P2", "",
                {transmonH[First[p]], FockState[0, 3], {0, Last[p]}},
                With[{t1 = Last[p]}, #[t1][[2 ;; 3]] &],
                If[kind === "windowed" && gap === 0, {MaxStepFraction -> 1}, {}]
            ]
        ]],
        {0, 100, 1000, 10000}
    ]],
    {"windowed", "untruncated"}
];


(* 5. A Pi/2 square pulse on [0, 1], then a Gaussian Pi/2 pulse (sigma = 0.025) at t = 50, under
   H = Omega(t) X/2: P1(100) = 1, and a skipped Gaussian leaves 1/2. Written as a product of two step
   functions, as #1265 writes its window, the factor whose only edge is the start of the range is
   where the landmark search fails; the single condition Boole[0 <= t <= 1] is the same function. *)
With[{gauss = (Pi / 2) / (Sqrt[2 Pi] 0.025) Exp[- (\[FormalT] - 50) ^ 2 / (2 0.025 ^ 2)]},
    KeyValueMap[
        caseRows["5 square from t0, then Gaussian", #1, "P1(100)", 1, {((#2 + gauss) / 2) QuantumOperator["X"], QuantumState["0"], {0, 100}}, #[100][[2]] &] &,
        <|
            "Boole[t>=0] Boole[t<=1]" -> (Pi / 2) Boole[\[FormalT] >= 0] Boole[\[FormalT] <= 1],
            "UnitStep[t] UnitStep[1-t]" -> (Pi / 2) UnitStep[\[FormalT]] UnitStep[1 - \[FormalT]],
            "Boole[0<=t<=1]" -> (Pi / 2) Boole[0 <= \[FormalT] <= 1]
        |>
    ]
];


(* 2c. Issue #2552: two Pi/2 pulses given as samples, here as Interpolation of the real and imaginary
   parts. QuantumEvolve's default "MergeInterpolatingFunctions" -> True joins them into a linear
   interpolation on the common grid, so the exact value for the problem QF solves is
   Cos[Theta]^2 with Theta the trapezoidal integral of the samples (the drive's imaginary part is
   below 1e-14 of its real part). Run last: the sampled drive is the slowest case. *)
omega0 = 10*^6 2 Pi;
readComplex[file_] := N[ToExpression[StringReplace[#, {"e+" -> "*^", "e-" -> "*^-", "j" -> " I", "(" -> "", ")" -> ""}]] & /@ ReadList[file, String]]

Scan[
    Function[{name}, With[{
        c = readComplex[FileNameJoin[{here, "shortpulse_data2552", name <> "_coef.txt"}]],
        tl = N @ ReadList[FileNameJoin[{here, "shortpulse_data2552", name <> "_tlist.txt"}], Number]
    },
        With[{
            fr = Interpolation[Transpose[{tl, Re[c]}]],
            fi = Interpolation[Transpose[{tl, Im[c]}]],
            theta = omega0 / 2 Total[Differences[tl] (Most[Re[c]] + Rest[Re[c]]) / 2]
        },
            caseRows["2c #2552", name, "|<0|psi(T)>|^2", Cos[theta] ^ 2,
                {(omega0 / 2) (fr[\[FormalT]] QuantumOperator["X"] + fi[\[FormalT]] QuantumOperator["Y"]), QuantumState["0"], {First[tl], Last[tl]}},
                With[{t1 = Last[tl]}, #[t1][[1]] &],
                If[name === "index93", {MaxStepFraction -> 1}, {}]
            ]
        ]
    ]],
    {"index92", "index93", "index94"}
];

Close[out];
