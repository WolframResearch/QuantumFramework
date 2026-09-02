#!/usr/bin/env wolframscript

(* Build the QEC documentation pages from their markdown sources.

   Usage: wolframscript -file OngoingProjects/QEC/docs/build.wls [pattern]

   The markdown is the source of truth; the .nb are generated and should not be
   hand-edited. Every docs/**/*.md is converted to the paclet's documentation
   layout, keyed by the frontmatter Template:

     Template: TechNote  ->  Tutorials/<Name>.nb
     Template: Symbol    ->  ReferencePages/Symbols/<Name>.nb
     Template: Guide     ->  Guides/<Name>.nb

   While the QEC layer lives outside the paclet the output stays beside the
   source, under docs/. After the migration into QuantumFramework/Kernel/QEC/,
   point $outRoot at the paclet's Documentation/English/ instead and run
   DocumentationBuild on the result -- MarkdownToNotebook emits an *authoring*
   notebook, and a page deployed without that build step renders with a doubled
   section rule and its examples folded under "Examples Initialization".

   RUN THIS IN A FRESH KERNEL. The pages evaluate their own example cells, and
   a kernel that has previously loaded the superseded OngoingProjects/QEC/QEC
   package will render QECCode summary boxes from *its* definition: the old
   guard KeyExistsQ[a, "Qubits"] matches the new data association too, and the
   old box then asks for a "Matrix" key the new representation does not have.
   The collision disappears once that package is deleted.
*)

$here = AbsoluteFileName[DirectoryName[$InputFileName]];
$qec = ParentDirectory[$here];
$repo = ParentDirectory[$qec, 2];

PacletDirectoryLoad[FileNameJoin[{$repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

(* Prefer a local checkout of the converter: it is where fixes land first. *)
$mtnLocal = FileNameJoin[{ParentDirectory[$repo], "MarkdownToNotebook", "MarkdownToNotebook.wl"}];

If[ FileExistsQ[$mtnLocal],
    Get[$mtnLocal];
    $convert = MarkdownToNotebook;
    Print["  Converter:  " <> $mtnLocal],
    $convert = ResourceFunction[ResourceObject[
        "https://www.wolframcloud.com/obj/nikm/DeployedResources/Function/MarkdownToNotebook"]];
    Print["  Converter:  cloud deployment (no local checkout found)"]
];

$templateDir = <|
    "TechNote" -> {"Tutorials"},
    "Symbol" -> {"ReferencePages", "Symbols"},
    "Guide" -> {"Guides"},
    "Overview" -> {"Guides"}
|>;

(* Read the Template: and Name: keys out of the frontmatter. *)
frontmatter[file_] := Module[{lines, grab},
    lines = TakeWhile[Rest[ReadList[file, "String"]], # =!= "---" &];
    grab[key_] := Replace[
        SelectFirst[lines, StringStartsQ[#, key <> ":"] &],
        {s_String :> StringTrim[StringDrop[s, StringLength[key] + 1]], _ -> Missing[key]}
    ];
    <|"Template" -> grab["Template"], "Name" -> grab["Name"]|>
];

$pattern = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], ""];

$sources = Select[
    FileNames["*.md", $here, Infinity],
    $pattern === "" || StringContainsQ[FileNameTake[#], $pattern, IgnoreCase -> True] &
];

Print["  Sources:    " <> ToString[Length[$sources]]];
Print[""];

build[src_] := Module[{fm, dir, out},
    fm = frontmatter[src];
    If[ MissingQ[fm["Template"]] || MissingQ[fm["Name"]],
        Print["  SKIP  " <> FileNameTake[src] <> "  (no Template:/Name: frontmatter)"];
        Return[False]
    ];
    dir = Lookup[$templateDir, fm["Template"], Missing["UnknownTemplate"]];
    If[ MissingQ[dir],
        Print["  SKIP  " <> FileNameTake[src] <> "  (unknown Template: " <> fm["Template"] <> ")"];
        Return[False]
    ];
    out = FileNameJoin[Join[{$here}, dir, {fm["Name"] <> ".nb"}]];
    Quiet[CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
    $convert[src, out];

    (* The DocumentationTools templates set CellContext -> CellGroup, which gives every
       cell group its own private context. The front end saves and restores $Context and
       $ContextPath per group, so a Needs or Get in one cell is discarded when that cell
       finishes: the page builds correctly (the converter evaluates it in one pass) but a
       reader running it cell by cell finds every symbol after the setup cell unresolved,
       and the cells come out unevaluated. Since these pages are meant to be run, not only
       read, the generated notebook is switched to one shared Global context. *)
    With[{text = Import[out, "Text"]},
        If[ StringContainsQ[text, "CellContext->CellGroup"],
            Export[out, StringReplace[text, "CellContext->CellGroup" -> "CellContext->\"Global`\""], "Text"]
        ]
    ];

    Print["  BUILT " <> fm["Template"] <> "  " <> FileNameTake[src] <> "  ->  " <>
        StringDrop[out, StringLength[$qec] + 1] <> "  (" <>
        ToString[Round[FileByteCount[out]/1024]] <> " KB)"];
    True
];

$ok = build /@ $sources;

Print[""];
Print["  Built " <> ToString[Count[$ok, True]] <> " of " <> ToString[Length[$sources]] <> "."];
Exit[If[AllTrue[$ok, TrueQ], 0, 1]]
