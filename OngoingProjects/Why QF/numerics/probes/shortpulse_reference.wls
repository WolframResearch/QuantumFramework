(* Exact values for the transmon case of the short-pulse head-to-head (case 4 of
   shortpulse_qutip.py and shortpulse_qf.wls). NDSolveValue integrates the three-level
   Schroedinger equation over the pulse alone, with exact parameters, at WorkingPrecision 40 and a
   step bound of 1/20 ns that no feature of a sigma = 4 ns pulse can hide between. |0> is
   stationary while the drive is off, so the values hold for every idle gap. A second solve at
   WorkingPrecision 50 with half the step bound gives the converged digits.
   Writes shortpulse_reference.tsv next to this file. Run: wolframscript -file shortpulse_reference.wls *)

here = DirectoryName[$InputFileName];

alpha = -2 Pi 2 / 10;
tT = 16;
sig = 4;
a = {{0, 1, 0}, {0, 0, Sqrt[2]}, {0, 0, 0}};
ad = ConjugateTranspose[a];
h0 = alpha / 2 ad . ad . a . a;
hd = (a + ad) / 2;

drive["windowed"] = Pi / (sig Sqrt[2 Pi] Erf[tT / (2 Sqrt[2] sig)]) Exp[- (t - tT / 2) ^ 2 / (2 sig ^ 2)];
drive["untruncated"] = Pi / (Sqrt[2 Pi] sig) Exp[- (t - 48) ^ 2 / (2 sig ^ 2)];
tEnd["windowed"] = tT;
tEnd["untruncated"] = 96;

probabilities[kind_, wp_, step_] := Abs[NDSolveValue[
    {psi'[t] == - I (h0 + drive[kind] hd) . psi[t], psi[0] == {1, 0, 0}},
    psi[tEnd[kind]],
    {t, 0, tEnd[kind]},
    WorkingPrecision -> wp, PrecisionGoal -> wp / 2, AccuracyGoal -> wp / 2, MaxStepSize -> step, MaxSteps -> Infinity
]] ^ 2

digits[x_] := StringReplace[ToString[N[x, 20], InputForm, NumberMarks -> False], "*^" -> "e"]

lines = Map[
    Function[{kind}, With[{p40 = probabilities[kind, 40, 1 / 20], p50 = probabilities[kind, 50, 1 / 40]},
        StringRiffle[{
            "4 transmon", kind, "reference", "NDSolveValue, WorkingPrecision 40, MaxStepSize 1/20, pulse only", "P0, P1, P2",
            StringRiffle[digits /@ p40, ", "], "", "", "",
            "largest change at WorkingPrecision 50, MaxStepSize 1/40: " <> digits[Max[Abs[p50 - p40]]] <> "; total probability - 1: " <> digits[Total[p40] - 1]
        }, "\t"]
    ]],
    {"windowed", "untruncated"}
];

Export[FileNameJoin[{here, "shortpulse_reference.tsv"}],
    StringRiffle[Join[
        {"# Wolfram Language " <> $Version <> ", " <> DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute"}],
         "case\tvariant\tpackage\tsetting\tobservable\tresult\texact\tabs_error\tseconds\tnote"},
        lines
    ], "\n"] <> "\n",
    "Text"
];
WriteString["stdout", StringRiffle[lines, "\n"], "\n"];
