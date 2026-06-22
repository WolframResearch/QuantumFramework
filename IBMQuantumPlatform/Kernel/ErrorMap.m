(* ErrorMap.m -- IBM device calibration model + error-map graphic.

   Get-loaded from IBMQuantumPlatform.m inside Begin["`Private`"], so every
   symbol here lives in IBMQuantumPlatform`Private`. Pure built-in graphics
   (Graph / Legended / BarLegend / ColorData): no QuantumFramework dependency.

   The IBM analogue of QUALink's QUAConnect[backend]["ErrorMap"]. The parse +
   draw layers are lifted verbatim from the verified proof of concept
   OngoingProjects/IBM ErrorMap/ibm-error-map.wls (QuantumFramework repo); only
   the transport is different: iIBMFetchModel reads the two free, read-only
   metadata endpoints through the service connection's own raw requests
   (RawBackendConfiguration / RawBackendProperties), so no circuit is run and no
   quantum-seconds are spent. *)

(* ============================================================================ *)
(*  Parse: configuration + properties -> normalized, provider-neutral model      *)
(* ============================================================================ *)

(* one named value out of IBM's [{name,value,unit,date}, ...] parameter lists *)
iIBMParam[plist_List, name_String] := Lookup[
  SelectFirst[plist, Lookup[#, "name", ""] === name &, <||>], "value", Missing["NotAvailable"]];

(* IBM's general-block pair metrics are named by a CONCATENATED digit string,
   "zz_7273" / "zz_717" / "zz_100101", whose split into two physical indices is
   ambiguous on its own (7273 -> 72,73 but 717 -> 7,17). Disambiguate against the
   coupling map: the unique split whose sorted pair is an actual coupling edge.
   Reject a part with a leading zero (except "0"). Returns the sorted pair or
   Missing when no/ambiguous split matches. *)
iIBMValidPart[s_String] := s === "0" || StringTake[s, 1] =!= "0";
iIBMPairFromName[name_String, coupSet_] := Module[{digits, cands},
  digits = Last[StringSplit[name, "_"], ""];
  cands = DeleteDuplicates @ Table[
    With[{a = StringTake[digits, k], b = StringTake[digits, {k + 1, -1}]},
      If[iIBMValidPart[a] && iIBMValidPart[b],
        With[{p = Sort[ToExpression /@ {a, b}]}, If[MemberQ[coupSet, p], p, Nothing]], Nothing]],
    {k, 1, StringLength[digits] - 1}];
  If[Length[cands] === 1, First[cands], Missing["UnresolvedPair", name]]
];

iIBMDeviceModel[config_Association, props_Association] := Module[{
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
      "T1"            -> iIBMParam[pl, "T1"],            (* us *)
      "T2"            -> iIBMParam[pl, "T2"],            (* us *)
      "ReadoutError"  -> iIBMParam[pl, "readout_error"],
      "ProbMeas0Prep1"-> iIBMParam[pl, "prob_meas0_prep1"],
      "ProbMeas1Prep0"-> iIBMParam[pl, "prob_meas1_prep0"],
      "ReadoutLength" -> iIBMParam[pl, "readout_length"] (* ns *)
    |>] &,
    qubits];

  (* per-qubit single-qubit (sx) gate error *)
  sxErr = Association @ Cases[gates,
    g_ /; Lookup[g, "gate", ""] === "sx" && MatchQ[Lookup[g, "qubits", {}], {_Integer}] :>
      (Lookup[g, "qubits"][[1]] -> iIBMParam[Lookup[g, "parameters", {}], "gate_error"])];

  (* per-pair cz error + duration, keyed by the sorted physical pair.
     error == 1 is IBM's uncalibrated sentinel (analogue of QUA's error == 0). *)
  czPairs = Cases[gates,
    g_ /; Lookup[g, "gate", ""] === "cz" && MatchQ[Lookup[g, "qubits", {}], {_Integer, _Integer}] :>
      (Sort[Lookup[g, "qubits"]] -> <|
        "Error"    -> iIBMParam[Lookup[g, "parameters", {}], "gate_error"],
        "Duration" -> iIBMParam[Lookup[g, "parameters", {}], "gate_length"]
      |>)];
  czModel = Association[czPairs];

  (* per-pair static ZZ crosstalk (|residual coupling|, GHz) from the general block,
     keyed by sorted pair. The always-on coupling a pair carries even while idle:
     an error channel independent of the cz gate_error. names disambiguated against
     the coupling map (iIBMPairFromName). *)
  general = Lookup[props, "general", {}];
  zzModel = Association @ Cases[general,
    g_ /; StringStartsQ[Lookup[g, "name", ""], "zz_"] :>
      With[{p = iIBMPairFromName[Lookup[g, "name"], coupSet]},
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
iIBMCZErrors[m_]      := Map[#["Error"] &, m["CZ"]];
iIBMCZDurations[m_]   := Map[#["Duration"] &, m["CZ"]];
iIBMZZ[m_]            := Lookup[m, "ZZ", <||>];          (* sorted-pair -> |ZZ| crosstalk, GHz *)
iIBMReadoutErrors[m_] := Map[#["ReadoutError"] &, m["Qubits"]];
iIBMT1[m_]            := Map[#["T1"] &, m["Qubits"]];
iIBMCoherence[m_]     := Map[KeyTake[#, {"T1", "T2"}] &, m["Qubits"]];
(* per-gate, per-locus error map for an arbitrary basis gate (cz / sx / x / ...) *)
iIBMGateErrors[m_, "cz"] := iIBMCZErrors[m];
iIBMGateErrors[m_, "sx"] := m["SXErrors"];
iIBMGateErrors[m_, _]    := <||>;

(* ============================================================================ *)
(*  ErrorMap: the device lattice annotated with the live error model             *)
(*    edge color     = log10 cz gate error/gate                                  *)
(*    edge thickness = static ZZ crosstalk (default; "EdgeThicknessMetric"):     *)
(*                     IBM publishes no Bell fidelity, so the second 2q channel   *)
(*                     is the always-on residual ZZ coupling, independent of the  *)
(*                     gate error. "Duration" (cz length) and None are the others.*)
(*    vertex size    = readout error ("VertexSizeReversed" flips it so big=good)  *)
(*    vertex color   = T1 (a bonus channel IBM affords)                          *)
(*    "EdgeArrows" adds double-headed arrows (cz is symmetric -> bidirectional).  *)
(*  Uncalibrated cz (error == 1) and qubits with no readout number are drawn as  *)
(*  dashed defect sites, so the physical lattice is shown complete.              *)
(* ============================================================================ *)

Options[iIBMErrorMap] = {
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

iIBMErrorMap[m_Association, opts : OptionsPattern[]] := Module[{
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
  czErr  = iIBMCZErrors[m]; czDur = iIBMCZDurations[m]; zz = iIBMZZ[m];
  roErr  = iIBMReadoutErrors[m]; t1 = iIBMT1[m]; sxErr = m["SXErrors"];

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
(*  Transport: fetch the model through the service connection's raw requests      *)
(*  (no new HTTP; reuses RawBackendConfiguration / RawBackendProperties, whose    *)
(*   Bearer-token + Service-CRN auth the connection already attaches)             *)
(* ============================================================================ *)

(* default backend: the JobRun PreprocessingFunction idiom, proven in this paclet *)
iIBMDefaultBackend[] := Enclose @ First[
  ConfirmMatch[SF`GetDefaultServiceObject["IBMQuantumPlatform"]["Backends"], {__String}],
  Confirm @ Failure["IBMQuantumPlatform", <|"MessageTemplate" -> "No backends available on this connection."|>]];

iIBMResolveBackend[b_] := If[StringQ[b], b, iIBMDefaultBackend[]];

(* fetch both free metadata endpoints via the connection's own raw requests (the
   IBMJob idiom: so["RawX", "BackendID" -> b], which applies the request's
   HTTPResponseProcessing and the connection's auth) *)
iIBMFetchModel[backend_String] := Enclose @ Module[{so, cfg, props},
  so    = ConfirmMatch[SF`GetDefaultServiceObject["IBMQuantumPlatform"], _ServiceObject];
  cfg   = ConfirmBy[so["RawBackendConfiguration", "BackendID" -> backend], AssociationQ];
  props = ConfirmBy[so["RawBackendProperties", "BackendID" -> backend], AssociationQ];
  iIBMDeviceModel[cfg, props]
];

(* resolve the "Backend" parameter (a list of rules or an association, whichever the
   service framework hands us) to a concrete device name *)
iIBMBackendFromParams[p_] := iIBMResolveBackend[Lookup[Association[p], "Backend", Automatic]];

(* the 9 styling keys are all STRING service-parameter names; map them to the option
   rules iIBMErrorMap consumes (OptionsPattern matches a bare list). Two names need a
   rename: iIBMErrorMap takes the canvas as the SYMBOL option Background, and the
   service parameter is "BackgroundColor" (the bare name "Background" collides with the
   ServiceExecute option namespace). *)
$iIBMErrorMapStyleKeys = {
  "EdgeColorScheme", "VertexColorScheme", "VertexColorReversed", "VertexSizeReversed",
  "EdgeThicknessMetric", "EdgeThicknessRange", "EdgeArrows", "ShowQubitLabels", "BackgroundColor"};

iIBMStyleRules[p_] := Normal @ KeyMap[
  Replace["BackgroundColor" -> Background],
  KeyTake[Association[p], $iIBMErrorMapStyleKeys]];
