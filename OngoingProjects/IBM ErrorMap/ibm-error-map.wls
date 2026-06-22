#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ibm-error-map.wls

   Standalone, kernel-free proof of concept for an IBM analogue of QUALink's
   QUAConnect[...]["ErrorMap"]. NOTHING here touches QuantumFramework kernel
   functions: it is pure built-in HTTP + Graph + graphics, so it can be tested
   in isolation. Once the visual is confirmed, the parse + accessor + ErrorMap
   layers port into the IBMQuantumPlatform connection / a QF backend handle.

   It uses only the FREE, read-only metadata endpoints (no circuit is run, no
   quantum-seconds are spent):
     GET /backends
     GET /backends/<name>/configuration   (coupling map + qubit grid coords)
     GET /backends/<name>/properties      (live calibration: T1/T2, readout,
                                           per-gate error/duration)

   Credentials. By design they are NOT hardcoded. The script reads, in order:
     - environment variables IBM_APIKEY / IBM_CRN, else
     - SystemCredential["IBMQuantumPlatform_APIKEY"] and LocalSymbol["ibm_crn"]
       (the same slots the QF IBMQuantumPlatform connection uses).

   Usage:
     wolframscript -file ibm-error-map.wls                 (default: ibm_fez)
     wolframscript -file ibm-error-map.wls ibm_marrakesh   (named backend)
*)

(* ============================================================================ *)
(*  Credentials                                                                  *)
(* ============================================================================ *)

ibmAPIKey[] := With[{env = Environment["IBM_APIKEY"]},
  If[StringQ[env], env, SystemCredential["IBMQuantumPlatform_APIKEY"]]];

ibmCRN[] := With[{env = Environment["IBM_CRN"]},
  If[StringQ[env], env, LocalSymbol["ibm_crn"]]];

(* ============================================================================ *)
(*  Transport: IAM token + read-only GET                                         *)
(* ============================================================================ *)

(* IAM bearer token from the API key. Mirrors IBMQuantumPlatform.m IBMToken. *)
ibmToken[apikey_String] := Enclose @ Module[{raw, resp},
  raw = URLRead @ HTTPRequest["https://iam.cloud.ibm.com/identity/token", <|
    Method -> "POST",
    "Content-Type" -> "application/x-www-form-urlencoded",
    "Body" -> <|
      "grant_type" -> "urn:ibm:params:oauth:grant-type:apikey",
      "apikey" -> apikey
    |>
  |>];
  ConfirmAssert[raw["StatusCode"] === 200, "IAM token HTTP " <> ToString[raw["StatusCode"]]];
  resp = ConfirmBy[ImportByteArray[raw["BodyByteArray"], "RawJSON"], AssociationQ];
  "Bearer " <> ConfirmBy[Lookup[resp, "access_token"], StringQ]
];

$ibmEndpoint = "https://quantum.cloud.ibm.com/api/v1";

(* read-only GET returning associations (RawJSON), never rules *)
ibmGet[token_String, crn_String, pathElems__] := Enclose @ Module[{resp},
  resp = URLRead @ HTTPRequest[
    URLBuild[{$ibmEndpoint, pathElems}],
    <|Method -> "GET", "Headers" -> {
      "Authorization" -> token, "Service-CRN" -> crn, "Content-Type" -> "application/json"}|>];
  ConfirmAssert[resp["StatusCode"] === 200, "HTTP " <> ToString[resp["StatusCode"]]];
  ImportByteArray[resp["BodyByteArray"], "RawJSON"]
];

ibmBackends[token_, crn_] := Lookup[ibmGet[token, crn, "backends"], "devices", {}];

(* ============================================================================ *)
(*  Parse: properties + configuration -> normalized, provider-neutral model      *)
(*  (the WL analogue of QUALink's iParseDeviceModel / iDevProps)                  *)
(* ============================================================================ *)

(* one named value out of IBM's [{name,value,unit,date}, ...] parameter lists *)
ibmParam[plist_List, name_String] := Lookup[
  SelectFirst[plist, Lookup[#, "name", ""] === name &, <||>], "value", Missing["NotAvailable"]];

(* IBM's general-block pair metrics are named by a CONCATENATED digit string,
   "zz_7273" / "zz_717" / "zz_100101", whose split into two physical indices is
   ambiguous on its own (7273 -> 72,73 but 717 -> 7,17). Disambiguate against the
   coupling map: the unique split whose sorted pair is an actual coupling edge.
   Reject a part with a leading zero (except "0"). Returns the sorted pair or
   Missing when no/ambiguous split matches. *)
ibmValidPart[s_String] := s === "0" || StringTake[s, 1] =!= "0";
ibmPairFromName[name_String, coupSet_] := Module[{digits, cands},
  digits = Last[StringSplit[name, "_"], ""];
  cands = DeleteDuplicates @ Table[
    With[{a = StringTake[digits, k], b = StringTake[digits, {k + 1, -1}]},
      If[ibmValidPart[a] && ibmValidPart[b],
        With[{p = Sort[ToExpression /@ {a, b}]}, If[MemberQ[coupSet, p], p, Nothing]], Nothing]],
    {k, 1, StringLength[digits] - 1}];
  If[Length[cands] === 1, First[cands], Missing["UnresolvedPair", name]]
];

ibmDeviceModel[config_Association, props_Association] := Module[{
  nq, coupling, coupSet, coords, qubits, gates, czPairs, sxErr, czModel, perQubit,
  general, zzModel
},
  nq       = Lookup[config, "n_qubits", 0];
  (* undirected coupling map: IBM lists both [i,j] and [j,i] *)
  coupling = DeleteDuplicates[Sort /@ Lookup[config, "coupling_map", {}]];
  coupSet  = coupling;
  (* per-qubit plot coordinate {col,row}; flip row so q index grows downward like IBM's plot *)
  coords   = AssociationThread[Range[0, nq - 1],
    MapIndexed[{#1[[1]], -#1[[2]]} &, Lookup[config, "coords", ConstantArray[{0, 0}, nq]]]];

  qubits = Lookup[props, "qubits", {}];
  gates  = Lookup[props, "gates", {}];

  (* per-qubit calibration *)
  perQubit = Association @ MapIndexed[
    With[{q = #2[[1]] - 1, pl = #1}, q -> <|
      "T1"            -> ibmParam[pl, "T1"],            (* us *)
      "T2"            -> ibmParam[pl, "T2"],            (* us *)
      "ReadoutError"  -> ibmParam[pl, "readout_error"],
      "ProbMeas0Prep1"-> ibmParam[pl, "prob_meas0_prep1"],
      "ProbMeas1Prep0"-> ibmParam[pl, "prob_meas1_prep0"],
      "ReadoutLength" -> ibmParam[pl, "readout_length"] (* ns *)
    |>] &,
    qubits];

  (* per-qubit single-qubit (sx) gate error *)
  sxErr = Association @ Cases[gates,
    g_ /; Lookup[g, "gate", ""] === "sx" && MatchQ[Lookup[g, "qubits", {}], {_Integer}] :>
      (Lookup[g, "qubits"][[1]] -> ibmParam[Lookup[g, "parameters", {}], "gate_error"])];

  (* per-pair cz error + duration, keyed by the sorted physical pair.
     error == 1 is IBM's uncalibrated sentinel (analogue of QUA's error == 0). *)
  czPairs = Cases[gates,
    g_ /; Lookup[g, "gate", ""] === "cz" && MatchQ[Lookup[g, "qubits", {}], {_Integer, _Integer}] :>
      (Sort[Lookup[g, "qubits"]] -> <|
        "Error"    -> ibmParam[Lookup[g, "parameters", {}], "gate_error"],
        "Duration" -> ibmParam[Lookup[g, "parameters", {}], "gate_length"]
      |>)];
  czModel = Association[czPairs];

  (* per-pair static ZZ crosstalk (|residual coupling|, GHz) from the general block,
     keyed by sorted pair. The always-on coupling a pair carries even while idle:
     an error channel independent of the cz gate_error. names disambiguated against
     the coupling map (ibmPairFromName). *)
  general = Lookup[props, "general", {}];
  zzModel = Association @ Cases[general,
    g_ /; StringStartsQ[Lookup[g, "name", ""], "zz_"] :>
      With[{p = ibmPairFromName[Lookup[g, "name"], coupSet]},
        If[MissingQ[p], Nothing, p -> Abs[Lookup[g, "value", Missing[]]]]]];

  <|
    "Backend"        -> Lookup[config, "backend_name", Lookup[props, "backend_name", "?"]],
    "NumQubits"      -> nq,
    "BasisGates"     -> Lookup[config, "basis_gates", {}],
    "ProcessorType"  -> Lookup[config, "processor_type", <||>],
    "LastUpdate"     -> Lookup[props, "last_update_date", Missing[]],
    "CouplingMap"    -> coupling,
    "Coords"         -> coords,
    "Qubits"         -> perQubit,
    "SXErrors"       -> sxErr,
    "CZ"             -> czModel,
    "ZZ"             -> zzModel
  |>
];

(* convenience projections (the QUAConnect-style pure accessors) *)
ibmCZErrors[m_]      := Map[#["Error"] &, m["CZ"]];
ibmCZDurations[m_]   := Map[#["Duration"] &, m["CZ"]];
ibmZZ[m_]            := Lookup[m, "ZZ", <||>];          (* sorted-pair -> |ZZ| crosstalk, GHz *)
ibmReadoutErrors[m_] := Map[#["ReadoutError"] &, m["Qubits"]];
ibmT1[m_]            := Map[#["T1"] &, m["Qubits"]];

(* ============================================================================ *)
(*  ErrorMap: the device lattice annotated with the live error model             *)
(*  Channels (mirroring QUAConnect[...]["ErrorMap"]):                            *)
(*    edge color     = log10 cz gate error/gate                                  *)
(*    edge thickness = static ZZ crosstalk (default; "EdgeThicknessMetric").     *)
(*                     IBM publishes no Bell fidelity, so the second 2q channel   *)
(*                     is the always-on residual ZZ coupling (an idle-error       *)
(*                     source independent of the gate error). "Duration" (cz gate *)
(*                     length) and None (uniform) are the other choices.          *)
(*    vertex size    = readout error ("VertexSizeReversed" flips it so big=good)  *)
(*    vertex color   = T1 (cooler = longer-lived; a bonus channel IBM affords)   *)
(*    "EdgeArrows" adds double-headed arrows (cz is symmetric -> bidirectional).  *)
(*  Uncalibrated cz (error == 1) and qubits with no readout number are drawn as  *)
(*  dashed defect sites, so the physical lattice is shown complete.              *)
(* ============================================================================ *)

(* Color schemes (and the other channels) are options so vertex and edge can be
   given DISTINCT built-in gradients (sharing a gradient makes a red edge and a
   red vertex ambiguous). The defaults are the curated look: avocado vertices
   (T1) over a black field, the cz-error temperature map on the bonds, size
   inverted so a big disk is a GOOD (low-readout) qubit, double-headed arrows on
   the symmetric cz bonds, and the physical qubit index printed in each disk. *)
Options[ibmErrorMap] = {
  "EdgeColorScheme"     -> "TemperatureMap",  (* cz error/gate *)
  "VertexColorScheme"   -> "AvocadoColors",   (* T1 *)
  "VertexColorReversed" -> True,              (* True: low T1 -> scheme's high (bright) end *)
  "VertexSizeReversed"  -> True,              (* True: big disk = LOW readout error (good) *)
  "EdgeThicknessMetric" -> "ZZ",              (* "ZZ" crosstalk | "Duration" cz length | None (uniform) *)
  "EdgeThicknessRange"  -> {4., 12.},         (* AbsoluteThickness pts, mapped from the metric *)
  "EdgeArrows"          -> True,              (* double-headed arrows (cz is symmetric, so bidirectional) *)
  "ShowQubitLabels"     -> True,              (* print the physical qubit index inside each disk *)
  Background            -> Black
};

ibmErrorMap[m_Association, opts : OptionsPattern[]] := Module[{
  eScheme = OptionValue["EdgeColorScheme"],
  vScheme = OptionValue["VertexColorScheme"],
  vRev = TrueQ @ OptionValue["VertexColorReversed"],
  vSizeRev = TrueQ @ OptionValue["VertexSizeReversed"],
  eMetric = OptionValue["EdgeThicknessMetric"],
  eThick = OptionValue["EdgeThicknessRange"],
  arrowsQ = TrueQ @ OptionValue["EdgeArrows"],
  showLabels = TrueQ @ OptionValue["ShowQubitLabels"],
  bg = OptionValue[Background],
  coords, cmap, czErr, czDur, zz, roErr, t1, sxErr, g,
  validErr, lo, hi, thickData, tvVals, tvlo, tvhi, roVals, rlo, rhi, t1Vals, tlo, thi,
  edgeColor, edgeThick, vertRadius, vertColor, pct, ns, us, khz, thickFmt, thickLabel,
  bgLum, dark, fg, vEdge, esf, vsf, qTip, qLabel, qLabelColor, defect, vlist
},
  coords = m["Coords"];   cmap = m["CouplingMap"];
  czErr  = ibmCZErrors[m]; czDur = ibmCZDurations[m]; zz = ibmZZ[m];
  roErr  = ibmReadoutErrors[m]; t1 = ibmT1[m]; sxErr = m["SXErrors"];

  (* light/dark theme: pick foreground (title/legend text) + vertex outline from
     the background luminance, so a dark Background keeps text and outlines visible *)
  bgLum = Replace[ColorConvert[Replace[bg, {Automatic | None -> White}], "GrayLevel"],
    {GrayLevel[y_, ___] :> y, _ -> 1.}];
  dark  = bgLum < 0.5;
  fg    = If[dark, GrayLevel[0.92], Black];
  vEdge = If[dark, GrayLevel[0.85], GrayLevel[0.25]];

  (* which per-pair quantity drives edge thickness *)
  khz[x_] := If[NumericQ[x], ToString[Round[x*10.^6, 0.1]] <> " kHz", "n/a"]; (* GHz -> kHz *)
  ns[x_]  := If[NumericQ[x], ToString[Round[x]] <> " ns", "n/a"];
  {thickData, thickFmt, thickLabel} = Switch[eMetric,
    "Duration", {czDur, ns,  "cz duration (edge thickness)"},
    None,       {<||>,  ns,  None},
    _,          {zz,    khz, "ZZ crosstalk (edge thickness)"}];  (* default: ZZ *)

  (* per-channel ranges, guarded against degenerate spreads + sentinels *)
  validErr = Select[Values[czErr], NumericQ[#] && 0 < # < 1 &]; If[validErr === {}, validErr = {0.001, 0.05}];
  {lo, hi} = MinMax[validErr]; If[lo == hi, hi = lo (1 + 1.*^-3)];
  tvVals = Select[Values[thickData], NumericQ]; If[tvVals === {}, tvVals = {0., 1.}];
  {tvlo, tvhi} = MinMax[tvVals]; If[tvlo == tvhi, tvhi = tvlo + Abs[tvlo] 0.01 + 1.*^-9];
  roVals = Select[Values[roErr], NumericQ]; If[roVals === {}, roVals = {0., 1.}];
  {rlo, rhi} = MinMax[roVals]; If[rlo == rhi, rhi = rlo + 1.*^-3];
  t1Vals = Select[Values[t1], NumericQ]; If[t1Vals === {}, t1Vals = {1., 2.}];
  {tlo, thi} = MinMax[t1Vals]; If[tlo == thi, thi = tlo + 1.];

  (* encoders *)
  edgeColor[e_] := If[NumericQ[e] && 0 < e < 1,
    ColorData[eScheme][Rescale[Log10[e], Log10[{lo, hi}]]], GrayLevel[0.6]];
  (* thickness from the chosen metric; None or a missing value -> midpoint width *)
  edgeThick[v_] := AbsoluteThickness[
    If[eMetric === None || ! NumericQ[v], Mean[eThick], Rescale[v, {tvlo, tvhi}, eThick]]];
  (* big disk = high readout error by default; reversed -> big disk = LOW error (good).
     Radius span kept well under the 0.5 half-spacing so adjacent disks never touch,
     leaving the edges (and arrowheads) a clear gap. *)
  vertRadius[ro_] := If[NumericQ[ro],
    Rescale[ro, {rlo, rhi}, If[vSizeRev, {0.30, 0.11}, {0.11, 0.30}]], 0.11];
  vertColor[t_] := If[NumericQ[t],
    ColorData[vScheme][With[{r = Rescale[t, {tlo, thi}]}, If[vRev, 1 - r, r]]], GrayLevel[0.75]];
  pct[x_] := If[NumericQ[x], ToString[Round[100 x, 0.01]] <> "%", "n/a"];
  us[x_]  := If[NumericQ[x], ToString[Round[x, 0.1]] <> " \[Mu]s", "n/a"];

  qTip[q_] := Grid[{
    {Style["q" <> ToString[q], Bold], SpanFromLeft},
    {"readout err", pct[Lookup[roErr, q, Missing[]]]},
    {"sx-gate err", pct[Lookup[sxErr, q, Missing[]]]},
    {"T1", us[Lookup[t1, q, Missing[]]]},
    {"T2", us[Lookup[Map[#["T2"] &, m["Qubits"]], q, Missing[]]]},
    {"meas0|prep1", pct[Lookup[Map[#["ProbMeas0Prep1"] &, m["Qubits"]], q, Missing[]]]},
    {"meas1|prep0", pct[Lookup[Map[#["ProbMeas1Prep0"] &, m["Qubits"]], q, Missing[]]]}
  }, Alignment -> Left, Frame -> All, FrameStyle -> GrayLevel[0.8]];

  defect[c_] := {FaceForm[None], EdgeForm[Directive[GrayLevel[0.55], Dashing[{0.01}]]], Disk[c, 0.16],
    GrayLevel[0.55], AbsoluteThickness[1.],
    Line[{c + {-.08, -.08}, c + {.08, .08}}], Line[{c + {-.08, .08}, c + {.08, -.08}}]};

  (* edge shape function: color by cz error, thickness by the chosen metric; sentinel -> dashed gray *)
  esf = Function[{pts, e}, With[{ij = Sort[List @@ e]},
    With[{err = Lookup[czErr, Key[ij], Missing[]], dur = Lookup[czDur, Key[ij], Missing[]],
          zzv = Lookup[zz, Key[ij], Missing[]], tv = Lookup[thickData, Key[ij], Missing[]]},
      Tooltip[
        If[NumericQ[err] && err >= 1 || MissingQ[err],
          {GrayLevel[0.6], AbsoluteThickness[0.7 First[eThick]], Dashing[{0.008}], Line[pts]},
          If[arrowsQ,
            (* double-headed (cz is symmetric): heads point outward near mid-edge,
               away from the disks so they stay visible *)
            {edgeColor[err], edgeThick[tv], Arrowheads[{{-0.014, 0.42}, {0.014, 0.58}}], Arrow[pts]},
            {edgeColor[err], edgeThick[tv], Line[pts]}]],
        Grid[{{Style["cz  q" <> ToString[ij[[1]]] <> " - q" <> ToString[ij[[2]]], Bold], SpanFromLeft},
          {"gate error", If[NumericQ[err] && err >= 1, "uncalibrated", pct[err]]},
          {"ZZ crosstalk", khz[zzv]},
          {"duration", ns[dur]}}, Alignment -> Left, Frame -> All, FrameStyle -> GrayLevel[0.8]]]]]];

  (* qubit-index label: black on a light fill, white on a dark fill (defect uses
     the theme foreground); font scales with the disk so it fits small disks *)
  qLabelColor[fill_] := If[fill === None, fg,
    If[Replace[ColorConvert[fill, "GrayLevel"], {GrayLevel[y_, ___] :> y, _ -> 0.5}] < 0.5,
      GrayLevel[0.95], GrayLevel[0.12]]];
  qLabel[c_, q_, r_, fill_] := If[showLabels,
    Text[Style[ToString[q], FontFamily -> "Helvetica", Bold,
      FontSize -> Clip[Rescale[r, {0.11, 0.30}, {6.5, 11.}], {6., 11.}],
      FontColor -> qLabelColor[fill]], c],
    Nothing];

  (* vertex shape function: disk sized by readout error, colored by T1; defect when no readout *)
  vsf = Function[{c, q}, With[{ro = Lookup[roErr, q, Missing[]], t = Lookup[t1, q, Missing[]]},
    If[NumericQ[ro],
      With[{r = vertRadius[ro], fill = vertColor[t]},
        Tooltip[{EdgeForm[vEdge], FaceForm[fill], Disk[c, r], qLabel[c, q, r, fill]}, qTip[q]]],
      Tooltip[{defect[c], qLabel[c, q, 0.16, None]}, qTip[q]]]]];

  g = Graph[Sort @ Keys[coords], UndirectedEdge @@@ cmap,
    VertexCoordinates -> KeyValueMap[#1 -> #2 &, KeySort @ coords]];
  vlist = VertexList[g];

  Legended[
    Graph[g,
      EdgeShapeFunction -> esf, VertexShapeFunction -> vsf,
      VertexCoordinates -> (coords[#] & /@ vlist),
      ImageSize -> 1100, ImagePadding -> 30, Background -> bg,
      PlotLabel -> Style[Row[{m["Backend"], "  \[Bullet]  ", m["NumQubits"], " qubits  \[Bullet]  calibrated ",
        m["LastUpdate"]}], 13, fg]],
    {
      Placed[BarLegend[{eScheme, {lo, hi}}, ScalingFunctions -> "Log10",
        LegendLabel -> Style["cz error/gate (edge color)", fg], LabelStyle -> fg, LegendLayout -> "Row"], Below],
      If[thickLabel === None, Nothing,
        Placed[LineLegend[{edgeThick[tvlo], edgeThick[(tvlo + tvhi)/2], edgeThick[tvhi]},
          {thickFmt[tvlo], thickFmt[(tvlo + tvhi)/2], thickFmt[tvhi]},
          LegendLabel -> Style[thickLabel, fg], LabelStyle -> fg, LegendLayout -> "Row"], Below]],
      Placed[BarLegend[{vertColor, {tlo, thi}},
        LegendLabel -> Style["T1 / \[Mu]s (vertex color)", fg], LabelStyle -> fg, LegendLayout -> "Row"], Below],
      Placed[PointLegend[{GrayLevel[0.75], GrayLevel[0.75]},
        If[vSizeRev, {pct[rhi], pct[rlo]}, {pct[rlo], pct[rhi]}],
        LegendMarkers -> {Graphics[{EdgeForm[vEdge], Disk[]}], Graphics[{EdgeForm[vEdge], Disk[]}]},
        LegendMarkerSize -> {8, 20}, LegendLabel -> Style["readout error (vertex size)", fg], LabelStyle -> fg, LegendLayout -> "Row"], Below]
    }]
];

(* ============================================================================ *)
(*  Driver                                                                       *)
(* ============================================================================ *)

ibmFetchModel[backend_String] := Enclose @ Module[{key, crn, tok, cfg, props},
  key = ibmAPIKey[]; crn = ibmCRN[];
  ConfirmAssert[StringQ[key], "No IBM API key (set $IBM_APIKEY or SystemCredential[\"IBMQuantumPlatform_APIKEY\"])."];
  ConfirmAssert[StringQ[crn], "No IBM CRN (set $IBM_CRN or LocalSymbol[\"ibm_crn\"])."];
  tok   = Confirm @ ibmToken[key];
  cfg   = Confirm @ ibmGet[tok, crn, "backends", backend, "configuration"];
  props = Confirm @ ibmGet[tok, crn, "backends", backend, "properties"];
  ibmDeviceModel[cfg, props]
];

If[!TrueQ[$ibmErrorMapLibraryOnly],
  Module[{backend, model, map, out},
    backend = Replace[$ScriptCommandLine, {{_, b_String, ___} :> b, _ :> "ibm_fez"}];
    Print["Fetching device model for ", backend, " ..."];
    model = ibmFetchModel[backend];
    If[FailureQ[model], Print[model]; Exit[1]];
    Print["  qubits: ", model["NumQubits"], " | basis: ", model["BasisGates"],
      " | cz pairs: ", Length[model["CZ"]], " | last update: ", model["LastUpdate"]];
    map = ibmErrorMap[model];
    out = FileNameJoin[{DirectoryName[$InputFileName], backend <> "-error-map.png"}];
    Export[out, map, ImageResolution -> 144];
    Print["Wrote ", out];
  ]
];
