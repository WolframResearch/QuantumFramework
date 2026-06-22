#!/usr/bin/env wolframscript
(* Offline comparison of vertex (T1) color schemes against fixed temperature-map
   edges. Uses the cached real ibm_fez calibration in /tmp/ibm/ -- NO network,
   NO credentials. Renders a 2x2 montage. *)

Block[{$ibmErrorMapLibraryOnly = True},
  Get[FileNameJoin[{DirectoryName[$InputFileName], "ibm-error-map.wls"}]]];

cfg   = Import["/tmp/ibm/config.json", "RawJSON"];
props = Import["/tmp/ibm/props.json", "RawJSON"];
m     = ibmDeviceModel[cfg, props];

(* vertex palettes deliberately chosen to share no hue with the blue->red edge
   map; T1 reversed so low (bad) T1 -> each scheme's bright/salient end *)
schemes = {"AvocadoColors", "SunsetColors", "BlueGreenYellow", "DeepSeaColors"};

panel[s_] := Labeled[
  ImageResize[
    Rasterize[ibmErrorMap[m, "VertexColorScheme" -> s, Background -> White],
      Background -> White],
    760],
  Style["vertex T1: " <> s, 18, Bold], Top];

grid = Grid[Partition[panel /@ schemes, 2], Spacings -> {1, 1}, Frame -> All,
  FrameStyle -> GrayLevel[0.7]];

out = FileNameJoin[{DirectoryName[$InputFileName], "vertex-color-comparison.png"}];
Export[out, grid, Background -> White];
Print["Wrote ", out];
