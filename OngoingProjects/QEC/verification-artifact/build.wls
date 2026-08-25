#!/usr/bin/env wolframscript
(* Build the computational essay .md into a .nb with MarkdownToNotebook (md2nb).
   Portable: finds the repo root from this script's own location, sets the working
   directory there so the essay's first cell resolves the package by walking up,
   then evaluates and caches every cell. *)

artifactDir = DirectoryName[$InputFileName];
repoRoot = SelectFirst[
    NestList[ParentDirectory, artifactDir, 10],
    FileExistsQ[FileNameJoin[{#, "QuantumFramework", "PacletInfo.wl"}]] &,
    $Failed
];
If[ repoRoot === $Failed,
    WriteString[$Output, "ERROR: could not locate the QuantumFramework repo root above ", artifactDir, "\n"];
    Exit[1]
];
SetDirectory[repoRoot];

md2nb = FileNameJoin[{$HomeDirectory, "Documents", "GitHub", "MarkdownToNotebook", "MarkdownToNotebook.wl"}];
If[ ! FileExistsQ[md2nb],
    WriteString[$Output, "ERROR: MarkdownToNotebook.wl not found at ", md2nb, "\n"];
    WriteString[$Output, "Clone https://github.com/WolframResearch/MarkdownToNotebook or set md2nb.\n"];
    Exit[1]
];
Get[md2nb];

src = FileNameJoin[{artifactDir, "exact-stabilizer-codes-cross-checked.md"}];
out = FileNameJoin[{artifactDir, "exact-stabilizer-codes-cross-checked.nb"}];
result = MarkdownToNotebook[src, out];

WriteString[$Output, "Built: ", ToString[result], "\n"];
