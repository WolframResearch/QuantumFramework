#!/usr/bin/env wolframscript
(* Build the Wolfram/QuantumFramework documentation notebooks from the literate
   markdown that sits beside them under
   QuantumFramework/Documentation/English/. The markdown is the source of
   truth; each page is written back to a .nb of the same base name.

   Run with:

       wolframscript -f build_docs.wls

   Optionally naming pages to rebuild, by base name:

       wolframscript -f build_docs.wls QuantumMachineLearning
*)

(* Pin the front end appearance to Light: a headless session resolves the
   Automatic appearance to Dark, which would invert colors in rasterized
   example output. *)

Quiet @ UsingFrontEnd[CurrentValue[$FrontEnd, LightDark] = "Light"]

repoDir = DirectoryName[$InputFileName]

pacletDir = FileNameJoin[{repoDir, "QuantumFramework"}]

$mtnSource = "/Users/swish/src/wolfram/AISkills/MarkdownToNotebook.wl"

PacletDirectoryLoad[pacletDir]

Needs["Wolfram`QuantumFramework`"]

(* The example cells of a page are evaluated as it is converted, and pages
   reach into the dependencies directly, so those contexts have to be on the
   path exactly as a reader's would be. *)

Needs["Wolfram`TensorNetworks`"]

Needs["Wolfram`Arrays`"]

mtn = If[ FileExistsQ[$mtnSource],
    Get[$mtnSource];
    MarkdownToNotebook,
    Quiet @ Check[If[! TrueQ[$CloudConnected], CloudConnect[]], Null];
    ResourceFunction["https://www.wolframcloud.com/obj/nikm/DeployedResources/Function/MarkdownToNotebook"]
]

(* NO DocumentationBuild stage.  The repository builds documentation itself and
   wants AUTHORING notebooks; a page already built in place does not survive
   PacletBuild, which drops it from the paclet.  Every other page in this
   paclet is in authoring form and these have to match. *)

pages = Select[
    FileNames["*.md", FileNameJoin[{pacletDir, "Documentation"}], Infinity],
    StringStartsQ[StringTrim @ Import[#, "Text"], "---"] &
]

If[ Length[$ScriptCommandLine] > 1,
    With[{wanted = Rest[$ScriptCommandLine]},
        pages = Select[pages, MemberQ[wanted, FileBaseName[#]] &]
    ]
]

Print["building ", Length[pages], " documentation page(s)"]

Scan[
    Function[src,
        mtn[src, FileNameJoin[{DirectoryName[src], FileBaseName[src] <> ".nb"}]];
        Print["  ", FileNameTake[src], " -> ", FileBaseName[src], ".nb"]
    ],
    Sort[pages]
]

Print["done"]
