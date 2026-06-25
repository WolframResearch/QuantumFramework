#!/usr/bin/env wolframscript
(* Assemble the six stabilizer-formalism summary boxes into a LinkedIn collage. *)

dir = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Stabilizer/summary-box-collage";
boxDir = FileNameJoin[{dir, "boxes"}];

img[name_] := Import[FileNameJoin[{boxDir, name <> ".png"}]];

red  = RGBColor[0.78, 0.13, 0.09];
dark = GrayLevel[0.13];
gray = GrayLevel[0.42];

cap[txt_] := Style[txt, FontFamily -> "Helvetica Neue", FontSize -> 23,
    FontWeight -> Bold, dark];

(* A captioned box: caption on top, the rendered summary box below. *)
cell[name_, caption_, w_] := Column[
    {cap[caption], Image[img[name], ImageSize -> w]},
    Alignment -> Center, Spacings -> 0.7
];

title = Column[{
    Style["The Stabilizer Formalism, at a Glance",
        FontFamily -> "Helvetica Neue", FontSize -> 44, FontWeight -> Bold, dark],
    Style["Stabilizer states, error-correcting codes, graph states, Clifford channels,\n"
        <> "and magic-state frames as first-class objects in the Wolfram Quantum Framework.",
        FontFamily -> "Helvetica Neue", FontSize -> 23, gray, LineSpacing -> {1, 6}]
    }, Alignment -> Center, Spacings -> 0.9];

footer = Style[
    "Wolfram  QuantumFramework   \[Bullet]   PauliStabilizer   \[Bullet]   GraphState   \[Bullet]   CliffordChannel   \[Bullet]   StabilizerFrame",
    FontFamily -> "Helvetica Neue", FontSize -> 19, gray];

rule = Image[ConstantArray[List @@ red, {8, 380}], ImageSize -> 380];

collage = Column[{
    title,
    rule,
    Row[{cell["01-PauliStabilizer-Bell", "Bell state  |\[CapitalPhi]\:207a\[RightAngleBracket]", 460],
         Spacer[44],
         cell["02-PauliStabilizer-GHZ", "GHZ state  |000\[RightAngleBracket]+|111\[RightAngleBracket]", 460]},
        Alignment -> Top],
    cell["03-PauliStabilizer-Steane", "Steane [[7,1,3]] quantum error-correcting code", 1000],
    Row[{cell["04-GraphState-Ring", "Graph (cluster) state", 316],
         Spacer[34],
         cell["05-CliffordChannel-Bell", "Clifford channel", 316],
         Spacer[34],
         cell["06-StabilizerFrame-Magic", "Magic-state frame (T gate)", 316]},
        Alignment -> Top],
    footer
    }, Alignment -> Center, Spacings -> 1.9];

panel = Framed[collage, FrameMargins -> 64, FrameStyle -> None, Background -> White,
    RoundingRadius -> 0];

out = UsingFrontEnd @ Rasterize[panel, "Image", Background -> White, ImageResolution -> 96];
Export[FileNameJoin[{dir, "stabilizer-formalism-collage.png"}], out];

(* A width-1080 copy, convenient for direct LinkedIn upload. *)
web = ImageResize[out, 1080];
Export[FileNameJoin[{dir, "stabilizer-formalism-collage-1080.png"}], web];

(* Script result: exported dimensions (the collage PNGs are the real product). *)
<|"full" -> ImageDimensions[out], "web1080" -> ImageDimensions[web]|>
