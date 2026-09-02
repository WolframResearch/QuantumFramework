#!/usr/bin/env wolframscript

(* Evaluate the example cells of a documentation source, in order, in one kernel,
   and print each result next to the expected value recorded in its `<!-- => ... -->`
   comment.  The point is to catch a stale hint before the build bakes it into a
   notebook: MarkdownToNotebook caches outputs by content hash, so a page whose
   cells were never run against the current package can build clean and be wrong.

   Cells are evaluated dynamically rather than written out as a script, because the
   first cell is what puts the QEC context on $ContextPath -- a static script would
   resolve every later symbol into the wrong context at parse time.

   Usage: wolframscript -file docs/validate.wls [pathOrPattern]
   Default is every docs/**/*.md. *)

$here = AbsoluteFileName[DirectoryName[$InputFileName]];
$repo = ParentDirectory[$here, 3];

PacletDirectoryLoad[FileNameJoin[{$repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

$arg = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], ""];

$sources = Select[
    FileNames["*.md", $here, Infinity],
    $arg === "" || StringContainsQ[#, $arg, IgnoreCase -> True] &
];

(* A cell is a ```wl fence; the hint, if any, is the first `<!-- => ... -->` after it.
   Directive lines (#| ...) are stripped, and #| eval: false marks a cell not to run. *)
cells[file_] := Module[{text, blocks},
    text = Import[file, "Text"];
    blocks = StringCases[
        text,
        "```wl" ~~ Shortest[body__] ~~ "```" ~~ Shortest[after___] ~~ ("```wl" | EndOfString) :>
            {body, StringCases[after, "<!-- =>" ~~ Shortest[h__] ~~ "-->" :> StringTrim[h], 1]},
        Overlaps -> True
    ];
    (* the pattern above double-counts; take the fences directly instead *)
    blocks = StringCases[text, "```wl\n" ~~ Shortest[body__] ~~ "\n```" :> body];
    MapThread[
        <|"Code" -> #1, "Hint" -> #2|> &,
        {blocks, hintsFor[text, blocks]}
    ]
];

(* The hint belonging to a cell is the first `<!-- => ... -->` between the end of that
   cell and the start of the next one. *)
hintsFor[text_, blocks_] := Module[{positions, spans},
    positions = First[StringPosition[text, "```wl\n" <> # <> "\n```", 1]] & /@ blocks;
    spans = Table[
        {positions[[i, 2]] + 1, If[i < Length[positions], positions[[i + 1, 1]] - 1, StringLength[text]]},
        {i, Length[positions]}
    ];
    Function[span,
        Replace[
            StringCases[StringTake[text, span], "<!-- =>" ~~ Shortest[h__] ~~ "-->" :> StringTrim[h], 1],
            {{h_} :> h, _ -> None}
        ]
    ] /@ spans
];

report[file_] := Module[{cs, skipped = 0},
    cs = cells[file];
    Print["=== " <> FileNameTake[file] <> ": " <> ToString[Length[cs]] <> " cells ==="];
    Do[
        Module[{code = cs[[i, "Code"]], hint = cs[[i, "Hint"]], body, result},
            If[ StringContainsQ[code, "#| eval: false"],
                skipped++;
                Print["[", i, "] SKIPPED (eval: false)"];
                Continue[]
            ];
            body = StringRiffle[Select[StringSplit[code, "\n"], ! StringStartsQ[StringTrim[#], "#|"] &], "\n"];
            result = Quiet @ Check[ToExpression[body], "$EVALUATION FAILED"];
            Print["[", i, "] ", StringReplace[StringTake[body, UpTo[70]], "\n" -> " / "]];
            If[hint =!= None, Print["    expected: ", hint]];
            Print["    got:      ", StringTake[ToString[Short[result, 6], InputForm], UpTo[300]]];
            Print[];
        ],
        {i, Length[cs]}
    ];
    Print["--- " <> ToString[Length[cs]] <> " cells, " <> ToString[skipped] <> " skipped ---"];
    Print[];
];

report /@ $sources;
