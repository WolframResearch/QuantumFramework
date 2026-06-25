#!/usr/bin/env wolframscript
(* Render each stabilizer-formalism summary box (expanded) to a white-background,
   tightly cropped PNG. *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

outDir = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/summary-box-collage/boxes";

(* One representative per stabilizer-formalism head. *)
objs = <|
    "01-PauliStabilizer-Bell"   -> PauliStabilizer[{"XX", "ZZ"}],
    "02-PauliStabilizer-GHZ"    -> PauliStabilizer[{"XXX", "ZZI", "IZZ"}],
    "03-PauliStabilizer-Steane" -> PauliStabilizer["SteaneCode"],
    "04-GraphState-Ring"        -> GraphState[CycleGraph[5]],
    "05-CliffordChannel-Bell"   -> CliffordChannel[PauliStabilizer[{"XX", "ZZ"}]],
    "06-StabilizerFrame-Magic"  -> StabilizerFrame[PauliStabilizer[1]["H", 1]]["T", 1]
|>;

(* Toggle the summary-box opener (Typeset`open$$ = False) to True so the expanded
   grid (tableau / destabilizers / graph thumbnail / coefficient) is shown. *)
openBoxes[expr_] := ToBoxes[expr] /.
    HoldPattern[Typeset`open$$ = False] :> (Typeset`open$$ = True);

ras[expr_] := UsingFrontEnd @ Rasterize[
    RawBoxes @ openBoxes[expr], "Image", Background -> White, ImageResolution -> 300];

(* Trim the white border, then re-pad with a small uniform white margin. *)
tidy[img_] := ImagePad[ImageCrop[img, Padding -> White], 24, White];

paths = Association @ KeyValueMap[
    Function[{k, o},
        Block[{img = tidy @ ras[o], file = FileNameJoin[{outDir, k <> ".png"}]},
            Export[file, img];
            k -> ImageDimensions[img]
        ]
    ],
    objs
];

(* Script result: per-box exported dimensions (the PNGs are the real product). *)
paths
