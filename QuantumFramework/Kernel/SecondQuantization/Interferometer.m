(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["QuantumInterferometer"]

PackageScope["QuantumInterferometerQ"]



(* ============================================================================ *)
(* QuantumInterferometer: a passive linear-optical network on m modes.          *)
(*                                                                              *)
(* The primary data is the m x m transfer matrix W, with the convention         *)
(*       W[[j, k]] = <1_j| U |1_k>,   i.e.   U a_k^dag U^dag = Sum_j W[[j,k]] a_j^dag *)
(* so a gate g1 followed by g2 has W = W2 . W1. A Fock-space circuit is built    *)
(* from W only on request ("CircuitOperator"), and photon statistics come from    *)
(* permanents of submatrices of W, with no Fock-space truncation.                *)
(*                                                                              *)
(* Internal representation:                                                     *)
(*   QuantumInterferometer[<|"Unitary" -> W, "Modes" -> m, "Method" -> method,   *)
(*       "Elements" -> {"PS"[a] -> {p}, "BS"[t, f] -> {p, q}, ...},              *)
(*       "Layers" -> {x1, x2, ...}|>]                                           *)
(* Elements are time ordered, in QuantumCircuitOperator shorthand; Layers hold   *)
(* each element's horizontal position in the mesh: integer k for the beam        *)
(* splitter of column k, k -/+ 1/4 for the phase of a mesh cell, and half        *)
(* integers for phases outside cells.                                           *)
(*                                                                              *)
(* Decompositions: Reck et al., PRL 73, 58 (1994); Clements et al., Optica 3,    *)
(* 1460 (2016). Both meshes use Clements' cell T(t, f) of Eq. 1, a phase f on    *)
(* the upper mode followed by a real beam splitter BS(t, 0).                    *)
(* ============================================================================ *)


QuantumInterferometer::usage = "QuantumInterferometer[u] decomposes the m\[Times]m unitary u into a mesh of beam splitters and phase shifters.
QuantumInterferometer[u, Method -> method] uses the \"Clements\" (default), \"ClementsPhaseEnd\" or \"Reck\" mesh.
QuantumInterferometer[{\"BS\"[\[Theta], \[Phi]] -> {p, q}, \"PS\"[\[Alpha]] -> p, \[Ellipsis]}, m] builds an m-mode interferometer from beam splitters and phase shifters.
QuantumInterferometer[\"Random\"[m]] and QuantumInterferometer[\"Fourier\"[m]] give a Haar-random and a discrete Fourier interferometer.";

QuantumInterferometer::symbolic = "The matrix has non-numeric entries; only numeric (exact or approximate) unitaries can be decomposed. Build the interferometer from a list of beam splitters and phase shifters instead.";
QuantumInterferometer::nonunitary = "The matrix is not unitary within the tolerance `1`.";
QuantumInterferometer::method = "Method `1` is not one of \"Clements\", \"ClementsPhaseEnd\" or \"Reck\".";
QuantumInterferometer::elem = "`1` is not a valid element; use \"BS\"[\[Theta], \[Phi]] -> {p, q} or \"PS\"[\[Alpha]] -> p.";
QuantumInterferometer::modes = "The elements act on mode `1`, outside the `2` modes of the interferometer.";
QuantumInterferometer::occ = "`1` is not a list of `2` nonnegative photon numbers.";
QuantumInterferometer::levels = "\"CircuitOperator\" needs the number of Fock levels per mode, as in interferometer[\"CircuitOperator\", d]; n photons need d >= n + 1.";
QuantumInterferometer::state = "The state must have `1` modes with the same number of levels each.";
QuantumInterferometer::cutoff = "The state holds up to `1` photons but each mode keeps only `2` levels; at least `3` levels are needed for a faithful result.";
QuantumInterferometer::compose ="Cannot compose interferometers on `1` and `2` modes.";
QuantumInterferometer::undefprop = "property `` is undefined for this interferometer";

Options[QuantumInterferometer] = {Method -> Automatic, Tolerance -> 10.^-10, "DropIdentities" -> True};



(* ============================================================================ *)
(* Blocks: single-photon actions of the optical elements                        *)
(* ============================================================================ *)

(* BeamSplitterOperator[{t, f}] on modes p, q *)
bsBlock[t_, f_] := {{Cos[t], - Exp[- I f] Sin[t]}, {Exp[I f] Sin[t], Cos[t]}}

(* Clements' cell T(t, f) = BS(t, 0) . PS(f on the upper mode), and its inverse for real t *)
cellBlock[t_, f_] := {{Exp[I f] Cos[t], - Sin[t]}, {Exp[I f] Sin[t], Cos[t]}}

cellInverseBlock[t_, f_] := {{Exp[- I f] Cos[t], Exp[- I f] Sin[t]}, {- Sin[t], Cos[t]}}

gateMatrix["BS"[t_, f_] -> {p_, q_}, m_] := ReplacePart[IdentityMatrix[m], Thread[Tuples[{p, q}, 2] -> Flatten[bsBlock[t, f]]]]

gateMatrix["PS"[a_] -> {p_}, m_] := ReplacePart[IdentityMatrix[m], {p, p} -> Exp[I a]]

elementsUnitary[elements_List, m_Integer] := Fold[gateMatrix[#2, m] . #1 &, IdentityMatrix[m], elements]


zeroQ[x_, tol_] := If[Precision[x] === Infinity, PossibleZeroQ[x], Abs[x] <= tol]

identityGateQ["BS"[t_, _] -> _, tol_] := zeroQ[t, tol]

identityGateQ["PS"[a_] -> _, tol_] := zeroQ[Exp[I a] - 1, tol]

unitaryQ[w_, tol_] := If[Precision[w] === Infinity,
    AllTrue[Flatten[w . ConjugateTranspose[w] - IdentityMatrix[Length[w]]], PossibleZeroQ],
    UnitaryMatrixQ[w, Tolerance -> tol]
]

normalAngle[x_] := If[NumericQ[x], Mod[x, 2 Pi, - Pi], x]



(* ============================================================================ *)
(* Nulling steps                                                                *)
(* ============================================================================ *)

(* Angles of the cell whose inverse, multiplied from the right on columns c and  *)
(* c + 1, zeroes a = u[[r, c]] against its right neighbour b = u[[r, c + 1]].   *)
rightNullAngles[a_, b_, tol_] := Which[
    zeroQ[a, tol], {0, 0},
    zeroQ[b, tol], {Pi / 2, 0},
    True, {ArcTan[Abs[b], Abs[a]], Arg[a Conjugate[b]]}
]

(* Angles of the cell that, multiplied from the left on rows r - 1 and r, zeroes *)
(* b = u[[r, c]] against the entry a = u[[r - 1, c]] above it.                 *)
leftNullAngles[a_, b_, tol_] := Which[
    zeroQ[b, tol], {0, 0},
    zeroQ[a, tol], {Pi / 2, 0},
    True, {ArcTan[Abs[a], Abs[b]], Arg[- b Conjugate[a]]}
]



(* ============================================================================ *)
(* Decompositions                                                               *)
(* ============================================================================ *)

(* One nulling step: the matrix with entry {r, c} zeroed, sowing the cell used.  *)
(* "Right" multiplies columns c, c + 1 from the right by a cell inverse, against *)
(* the right neighbour; "Left" multiplies rows r - 1, r from the left by a cell, *)
(* against the entry above.                                                     *)
nullStep[u_, {"Right", {r_, c_}}, tol_, simp_] := With[{angles = simp /@ rightNullAngles[u[[r, c]], u[[r, c + 1]], tol]},
    Sow[Append[angles, c], "Right"];
    Module[{v = u},
        v[[All, {c, c + 1}]] = simp[v[[All, {c, c + 1}]] . (cellInverseBlock @@ angles)];
        v[[r, c]] = 0;
        v
    ]
]

nullStep[u_, {"Left", {r_, c_}}, tol_, simp_] := With[{angles = simp /@ leftNullAngles[u[[r - 1, c]], u[[r, c]], tol]},
    Sow[Append[angles, r - 1], "Left"];
    Module[{v = u},
        v[[{r - 1, r}]] = simp[(cellBlock @@ angles) . v[[{r - 1, r}]]];
        v[[r, c]] = 0;
        v
    ]
]

(* Fold the steps over W: {right cells, diagonal, left cells} with              *)
(* L_k ... L_1 . W . R_1^-1 ... R_n^-1 = D.                                    *)
nullingSweep[w_, steps_, tol_, simp_] := With[{reaped = Reap[Fold[nullStep[#1, #2, tol, simp] &, w, steps], {"Right", "Left"}]},
    {Catenate[reaped[[2, 1]]], simp /@ Diagonal[reaped[[1]]], Catenate[reaped[[2, 2]]]}
]

(* Reck: the entries below the diagonal, row by row from the bottom, each against *)
(* its right neighbour, so W = D . R_K ... R_1: the cells, then the phase column. *)
reckSteps[m_] := Catenate @ Table[{"Right", {r, c}}, {r, m, 2, -1}, {c, r - 1}]

(* Clements, supplementary material: along the i-th antidiagonal, from the right *)
(* for odd i and from the left for even i, so W = L_1^-1 ... L_k^-1 . D . R_n ... R_1. *)
clementsSteps[m_] := Catenate @ Table[
    If[ OddQ[i],
        Table[{"Right", {m - j, i - j}}, {j, 0, i - 1}],
        Table[{"Left", {m + j - i, j}}, {j, i}]
    ],
    {i, m - 1}
]

(* Move the phase column through one inverse cell, sowing the forward cell left  *)
(* behind: T^-1(t, f) . diag(e^ia, e^ib) = diag(e^i(b - f + Pi), e^ib) . T(t, a - b + Pi). *)
pushPhases[angles_, {t_, f_, p_}, simp_] := With[{a = angles[[p]], b = angles[[p + 1]]},
    Sow[{t, simp @ normalAngle[a - b + Pi], p}, "Moved"];
    ReplacePart[angles, p -> simp @ normalAngle[b - f + Pi]]
]

(* Pushing through the inverse cells, last one first, leaves                     *)
(* W = D' . T'_1 ... T'_k . R_n ... R_1 (Clements Eq. 5).                        *)
clementsPhaseEnd[{rights_, diagonal_, lefts_}, simp_] := With[{
    reaped = Reap[Fold[pushPhases[#1, #2, simp] &, simp /@ Arg[diagonal], Reverse[lefts]], "Moved"]
},
    Join[
        {"Cell", ##} & @@@ rights,
        {"Cell", ##} & @@@ Catenate[reaped[[2]]],
        MapIndexed[{"PS", #1, First[#2]} &, reaped[[1]]]
    ]
]

(* A decomposition gives a time-ordered list of items: {"Cell", t, f, p},        *)
(* {"InverseCell", t, f, p}, {"BS", t, f, {p, q}} and {"PS", a, p}.              *)
decompositionItems["Reck", w_, tol_, simp_] := With[{dec = nullingSweep[w, reckSteps[Length[w]], tol, simp]},
    Join[{"Cell", ##} & @@@ dec[[1]], MapIndexed[{"PS", simp[Arg[#1]], First[#2]} &, dec[[2]]]]
]

decompositionItems["Clements", w_, tol_, simp_] := With[{dec = nullingSweep[w, clementsSteps[Length[w]], tol, simp]},
    Join[
        {"Cell", ##} & @@@ dec[[1]],
        MapIndexed[{"PS", simp[Arg[#1]], First[#2]} &, dec[[2]]],
        {"InverseCell", ##} & @@@ Reverse[dec[[3]]]
    ]
]

decompositionItems["ClementsPhaseEnd", w_, tol_, simp_] :=
    clementsPhaseEnd[nullingSweep[w, clementsSteps[Length[w]], tol, simp], simp]



(* ============================================================================ *)
(* Mesh layout: items to elements and positions                                 *)
(* ============================================================================ *)

(* Place one item given the last occupied column of every mode: each beam       *)
(* splitter goes to the first column free on all the modes it spans, and a phase *)
(* outside a cell half a column after the last beam splitter on its mode.        *)
(* Returns {placed elements, updated last columns}, the FoldPairList contract.   *)
placeItem[last_, {"Cell", t_, f_, p_}] := With[{k = 1 + Max[last[[{p, p + 1}]]]},
    {{{"PS"[f] -> {p}, k - 1/4}, {"BS"[t, 0] -> {p, p + 1}, k}}, ReplacePart[last, {p -> k, p + 1 -> k}]}
]

placeItem[last_, {"InverseCell", t_, f_, p_}] := With[{k = 1 + Max[last[[{p, p + 1}]]]},
    {{{"BS"[- t, 0] -> {p, p + 1}, k}, {"PS"[- f] -> {p}, k + 1/4}}, ReplacePart[last, {p -> k, p + 1 -> k}]}
]

placeItem[last_, {"BS", t_, f_, {p_, q_}}] := With[{span = Range[Min[p, q], Max[p, q]]},
    With[{k = 1 + Max[last[[span]]]},
        {{{"BS"[t, f] -> {p, q}, k}}, ReplacePart[last, Thread[span -> k]]}
    ]
]

placeItem[last_, {"PS", a_, p_}] := {{{"PS"[a] -> {p}, last[[p]] + 1/2}}, last}

(* A phase outside a cell with no beam splitter after it on its mode moves to    *)
(* the end of the mesh, so a trailing phase column lines up.                     *)
alignTrailingPhases[placed_, m_] := With[{
    lastSplitter = Fold[
        ReplacePart[#1, Thread[Range @@ MinMax[#2[[2]]] -> #2[[1]]]] &,
        ConstantArray[0, m],
        Cases[MapIndexed[{First[#2], First[#1]} &, placed], {i_, "BS"[__] -> modes_} :> {i, modes}]
    ],
    end = Max[0, Cases[placed, {"BS"[__] -> _, x_} :> x]] + 1/2
},
    MapIndexed[
        Replace[#1, {gate : ("PS"[_] -> {p_}), x_} /; IntegerQ[x - 1/2] && First[#2] > lastSplitter[[p]] :> {gate, end}] &,
        placed
    ]
]

layoutItems[items_List, m_Integer, drop_, tol_] := With[{
    placed = alignTrailingPhases[Catenate @ FoldPairList[placeItem, ConstantArray[0, m], items], m]
},
    With[{kept = If[TrueQ[drop], DeleteCases[placed, {gate_, _} /; identityGateQ[gate, tol]], placed]},
        If[kept === {}, {{}, {}}, Transpose[kept]]
    ]
]

makeInterferometer[w_, m_, method_, items_, drop_, tol_] := With[{layout = layoutItems[items, m, drop, tol]},
    QuantumInterferometer[<|
        "Unitary" -> w,
        "Modes" -> m,
        "Method" -> method,
        "Elements" -> layout[[1]],
        "Layers" -> layout[[2]]
    |>]
]



(* ============================================================================ *)
(* Constructors                                                                 *)
(* ============================================================================ *)

QuantumInterferometerQ[QuantumInterferometer[KeyValuePattern[{
    "Unitary" -> _ ? SquareMatrixQ, "Modes" -> _Integer, "Method" -> _, "Elements" -> _List, "Layers" -> _List
}]]] := True

QuantumInterferometerQ[_] := False


QuantumInterferometer[w_ ? SquareMatrixQ, opts : OptionsPattern[]] := Module[{
    m = Length[w],
    method = Replace[OptionValue[Method], Automatic -> "Clements"],
    tol = OptionValue[Tolerance]
},
    If[! MatrixQ[w, NumericQ], Message[QuantumInterferometer::symbolic]; Return[$Failed]];
    If[! MemberQ[{"Clements", "ClementsPhaseEnd", "Reck"}, method], Message[QuantumInterferometer::method, method]; Return[$Failed]];
    If[! unitaryQ[w, tol], Message[QuantumInterferometer::nonunitary, tol]; Return[$Failed]];
    makeInterferometer[
        w, m, method,
        decompositionItems[method, w, tol, If[Precision[w] === Infinity, Simplify, Identity]],
        OptionValue["DropIdentities"], tol
    ]
]


normalizeElement[("BS" | "BeamSplitter")[{t_, f_}] -> modes_] := normalizeElement["BS"[t, f] -> modes]

normalizeElement[("BS" | "BeamSplitter")[t_, f_ : 0] -> {p_Integer ? Positive, q_Integer ? Positive}] /; p != q := {"BS", t, f, {p, q}}

normalizeElement[("PS" | "PhaseShift")[a_] -> (p_Integer ? Positive | {p_Integer ? Positive})] := {"PS", a, p}

normalizeElement[element_] := (Message[QuantumInterferometer::elem, element]; $Failed)

QuantumInterferometer[elements : {___Rule}, m_Integer ? Positive, opts : OptionsPattern[]] := Module[{items = normalizeElement /@ elements, layout},
    If[MemberQ[items, $Failed], Return[$Failed]];
    With[{top = Max[0, Flatten[items[[All, -1]]]]},
        If[top > m, Message[QuantumInterferometer::modes, top, m]; Return[$Failed]]
    ];
    layout = layoutItems[items, m, OptionValue["DropIdentities"], OptionValue[Tolerance]];
    QuantumInterferometer[<|
        "Unitary" -> elementsUnitary[layout[[1]], m],
        "Modes" -> m,
        "Method" -> None,
        "Elements" -> layout[[1]],
        "Layers" -> layout[[2]]
    |>]
]

QuantumInterferometer[elements : {__Rule}, opts : OptionsPattern[]] :=
    QuantumInterferometer[elements, Max[1, Cases[Flatten[{elements[[All, 2]]}], _Integer]], opts]


QuantumInterferometer["Random"[m_Integer ? Positive], opts : OptionsPattern[]] :=
    QuantumInterferometer[RandomVariate[CircularUnitaryMatrixDistribution[m]], opts]

QuantumInterferometer["Fourier"[m_Integer ? Positive], opts : OptionsPattern[]] :=
    QuantumInterferometer[FourierMatrix[m], opts]

QuantumInterferometer[name : "Random" | "Fourier", m_Integer ? Positive, opts : OptionsPattern[]] :=
    QuantumInterferometer[name[m], opts]

(* re-mesh with another method *)
QuantumInterferometer[qi_QuantumInterferometer ? QuantumInterferometerQ, opts : OptionsPattern[]] :=
    QuantumInterferometer[qi["Unitary"], opts]



(* ============================================================================ *)
(* Properties                                                                   *)
(* ============================================================================ *)

$QuantumInterferometerProperties = {
    "Unitary", "TransferMatrix", "Modes", "Method", "Elements", "Layers",
    "Depth", "BeamSplitterCount", "PhaseShifterCount", "Parameters",
    "Diagram", "MatrixPlot", "CircuitOperator",
    "Amplitude", "Probability", "Probabilities", "Dagger"
};

(qi_QuantumInterferometer ? QuantumInterferometerQ)[prop_String, args___] := With[{result = QuantumInterferometerProp[qi, prop, args]},
    result /; ! MatchQ[Unevaluated @ result, _QuantumInterferometerProp] || Message[QuantumInterferometer::undefprop, prop]
]

QuantumInterferometerProp[_, "Properties"] := $QuantumInterferometerProperties

QuantumInterferometerProp[QuantumInterferometer[data_], key : "Unitary" | "Modes" | "Method" | "Elements" | "Layers"] := data[key]

QuantumInterferometerProp[qi_, "TransferMatrix"] := qi["Unitary"]

QuantumInterferometerProp[qi_, "Depth"] := Max[0, Pick[qi["Layers"], MatchQ["BS"[__] -> _] /@ qi["Elements"]]]

QuantumInterferometerProp[qi_, "BeamSplitterCount"] := Count[qi["Elements"], "BS"[__] -> _]

QuantumInterferometerProp[qi_, "PhaseShifterCount"] := Count[qi["Elements"], "PS"[_] -> _]

QuantumInterferometerProp[qi_, "Parameters"] := Dataset @ MapThread[
    <|"Element" -> Head[#1[[1]]], "Parameters" -> List @@ #1[[1]], "Modes" -> #1[[2]], "Position" -> #2|> &,
    {qi["Elements"], qi["Layers"]}
]

QuantumInterferometerProp[qi_, "Diagram", opts : OptionsPattern[]] := interferometerDiagram[qi, opts]

QuantumInterferometerProp[qi_, "MatrixPlot", opts : OptionsPattern[]] := interferometerMatrixPlot[qi, opts]

QuantumInterferometerProp[qi_, "CircuitOperator", d_Integer ? Positive, opts : OptionsPattern[]] := interferometerCircuit[qi, d, opts]

QuantumInterferometerProp[qi_, "CircuitOperator"] := (Message[QuantumInterferometer::levels]; $Failed)

(* Photon patterns are occupation lists {n1, n2, ...} or Ket[{n1, n2, ...}];    *)
(* results are keyed by Ket so they read in Dirac notation.                     *)

(* single entries cost one permanent each, for when the full distribution is too large *)
QuantumInterferometerProp[qi_, "Amplitude", s_ -> t_] := With[{in = occupation[s], out = occupation[t]},
    If[occupationQ[qi, in] && occupationQ[qi, out], permanentAmplitude[qi["Unitary"], in, out], $Failed]
]

QuantumInterferometerProp[qi_, "Probability", s_ -> t_] := With[{in = occupation[s], out = occupation[t]},
    If[occupationQ[qi, in] && occupationQ[qi, out], Abs[permanentAmplitude[qi["Unitary"], in, out]] ^ 2, $Failed]
]

QuantumInterferometerProp[qi_, "Probabilities", s_] := With[{in = occupation[s]},
    If[ occupationQ[qi, in],
        With[{amplitudes = photonAmplitudes[qi["Unitary"], in]},
            Association @ Map[Ket[#] -> Abs[Lookup[amplitudes, Key[#], 0]] ^ 2 &, photonPatterns[Total[in], qi["Modes"]]]
        ],
        $Failed
    ]
]

inverseElement["BS"[t_, f_] -> modes_] := "BS"[- t, f] -> modes

inverseElement["PS"[a_] -> modes_] := "PS"[- a] -> modes

QuantumInterferometerProp[qi_, "Dagger"] := With[{depth = qi["Depth"]},
    QuantumInterferometer[<|
        "Unitary" -> ConjugateTranspose[qi["Unitary"]],
        "Modes" -> qi["Modes"],
        "Method" -> None,
        "Elements" -> Reverse[inverseElement /@ qi["Elements"]],
        "Layers" -> Reverse[depth + 1 - qi["Layers"]]
    |>]
]



(* ============================================================================ *)
(* Photon statistics in the photon-number sector                                *)
(* ============================================================================ *)

occupation[Ket[s_List]] := s

occupation[s_] := s

occupationQ[qi_, s_] := If[
    VectorQ[s, IntegerQ[#] && # >= 0 &] && Length[s] == qi["Modes"],
    True,
    Message[QuantumInterferometer::occ, s, qi["Modes"]]; False
]

(* the row (or column) list of an occupation pattern: mode k repeated s_k times *)
modeList[occupation_] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occupation]

photonPatterns[n_, m_] := ReverseSort @ FrobeniusSolve[ConstantArray[1, m], n]

(* <T| U |S> = Perm(W_{T,S}) / Sqrt[Prod s! Prod t!] *)
permanentAmplitude[w_, s_, t_] := Which[
    Total[s] != Total[t], 0,
    Total[s] == 0, 1,
    True, Permanent[w[[modeList[t], modeList[s]]]] / Sqrt[Times @@ (s!) Times @@ (t!)]
]

(* All output amplitudes at once: U |S> = Prod_k (Sum_j W[[j, k]] a_j^dag)^s_k / Sqrt[s_k!] |0>, *)
(* so the amplitude of T is its monomial's coefficient times Sqrt[Prod t! / Prod s!]. *)
photonAmplitudes[w_, s_] := Module[{x},
    With[{vars = Array[x, Length[w]]},
        Association @ ReverseSortBy[First] @ Map[
            #[[1]] -> #[[2]] Sqrt[Times @@ (#[[1]]!) / Times @@ (s!)] &,
            CoefficientRules[Expand[Times @@ (vars . w[[All, #]] & /@ modeList[s])], vars]
        ]
    ]
]

(qi_QuantumInterferometer ? QuantumInterferometerQ)[s : _List | Ket[_List]] /; VectorQ[occupation[s], IntegerQ] := With[{in = occupation[s]},
    If[occupationQ[qi, in], KeyMap[Ket, photonAmplitudes[qi["Unitary"], in]], $Failed]
]



(* ============================================================================ *)
(* Fock-space circuit                                                           *)
(* ============================================================================ *)

elementOperator["BS"[t_, f_] -> {p_, q_}, d_, opts___] :=
    BeamSplitterOperator[{t, f}, d, {p, q}, Sequence @@ FilterRules[{opts}, Options[BeamSplitterOperator]]]

elementOperator["PS"[a_] -> {p_}, d_, ___] := PhaseShiftOperator[a, d, {p}]

interferometerCircuit[qi_, d_, opts___] := With[{
    idle = Complement[Range[qi["Modes"]], Flatten[qi["Elements"][[All, 2]]]]
},
    QuantumCircuitOperator @ Join[
        QuantumOperator[IdentityMatrix[d, SparseArray], {#}, d] & /@ idle,
        elementOperator[#, d, opts] & /@ qi["Elements"]
    ]
]

(* the Fock components {occupation, amplitude} of a state vector with d levels per mode *)
fockComponents[qs_, d_, m_] := {IntegerDigits[#[[1, 1]] - 1, d, m], #[[2]]} & /@
    Most @ ArrayRules[SparseArray[qs["Computational"]["StateVector"]]]

(* U Sum_S c_S |S> from the permanents, keeping the patterns that fit in d levels *)
photonSectorState[w_, components_, d_, m_] := With[{
    amplitudes = KeySelect[Merge[#2 photonAmplitudes[w, #1] & @@@ components, Total], Max[#] < d &]
},
    QuantumState[SparseArray[KeyValueMap[FromDigits[#1, d] + 1 -> #2 &, amplitudes], d ^ m], ConstantArray[d, m]]
]

(* A state vector goes through the permanents; a density matrix through the     *)
(* Fock-space circuit. Either way the result keeps the input's d levels per mode. *)
(qi_QuantumInterferometer ? QuantumInterferometerQ)[qs_QuantumState, opts : OptionsPattern[]] := Module[{
    dims = qs["Dimensions"], m = qi["Modes"], d, components, photons
},
    If[Length[dims] != m || ! Equal @@ dims, Message[QuantumInterferometer::state, m]; Return[$Failed]];
    d = First[dims];
    components = If[qs["StateType"] === "Vector", fockComponents[qs, d, m], None];
    photons = If[ components === None,
        Max[0, Total[IntegerDigits[# - 1, d, m]] & /@
            Flatten @ Position[Normal[qs["ProbabilitiesList"]], p_ /; ! PossibleZeroQ[p], {1}, Heads -> False]],
        Max[0, Total /@ components[[All, 1]]]
    ];
    If[photons >= d, Message[QuantumInterferometer::cutoff, photons, d, photons + 1]];
    If[ components === None,
        qi["CircuitOperator", d, opts][qs],
        photonSectorState[qi["Unitary"], components, d, m]
    ]
]



(* ============================================================================ *)
(* Composition                                                                  *)
(* ============================================================================ *)

(* qi2[qi1] is qi1 followed by qi2 *)
(qi2_QuantumInterferometer ? QuantumInterferometerQ)[qi1_QuantumInterferometer ? QuantumInterferometerQ] := If[
    qi1["Modes"] != qi2["Modes"],
    Message[QuantumInterferometer::compose, qi1["Modes"], qi2["Modes"]]; $Failed,
    QuantumInterferometer[<|
        "Unitary" -> qi2["Unitary"] . qi1["Unitary"],
        "Modes" -> qi1["Modes"],
        "Method" -> None,
        "Elements" -> Join[qi1["Elements"], qi2["Elements"]],
        "Layers" -> Join[qi1["Layers"], qi2["Layers"] + Ceiling[Max[0, qi1["Layers"]]]]
    |>]
]



(* ============================================================================ *)
(* Mesh diagram                                                                 *)
(* ============================================================================ *)

$couplerHalfWidth = 0.3;

reflectivityColor[t_] := If[NumericQ[t],
    Blend[{RGBColor[0.62, 0.74, 0.88], RGBColor[0.05, 0.25, 0.6]}, N[Sin[t] ^ 2]],
    GrayLevel[0.3]
]

phaseColor[a_] := If[NumericQ[a], Hue[Mod[N[a] / (2 Pi), 1], 0.7, 0.95], White]

parameterForm[v_] := If[NumericQ[v], NumberForm[Chop[N[v]], {3, 2}], v]

(* a waveguide that leaves the height y0 at x - w, runs at y1 through the coupling *)
(* region in the middle half, and returns to y0 at x + w                         *)
waveguide[x_, w_, y0_, y1_] := Line @ Table[
    {x + u w, y0 + (y1 - y0) (1 - Cos[Pi Min[1, 2 (1 - Abs[u])]]) / 2},
    {u, -1, 1, 1/20}
]

(* the pieces of the wire on mode p in [x0, x1] that no coupler covers *)
wireSegments[p_, couplers_, {x0_, x1_}] := With[{
    blocked = SortBy[Cases[couplers, {x_, p, _, _, _} | {x_, _, p, _, _} :> {x - $couplerHalfWidth, x + $couplerHalfWidth}], First]
},
    Partition[Flatten[{x0, blocked, x1}], 2]
]

(* A coupler is labelled by its pair (t, f): for a mesh cell f is the cell's     *)
(* phase on the upper mode, which is drawn as part of the coupler; for a plain   *)
(* beam splitter it is the beam splitter's own phase.                           *)
couplerPrimitive[{x_, p_, q_, t_, f_}, labelQ_] := With[{w = $couplerHalfWidth, mid = - (p + q) / 2, gap = 0.08},
    Tooltip[
        {
            reflectivityColor[t], AbsoluteThickness[2],
            waveguide[x, w, - p, mid + gap],
            waveguide[x, w, - q, mid - gap],
            If[ labelQ,
                {Black, Text[Style[Row[{"(", parameterForm[t], ", ", parameterForm[f], ")"}], 8], {x, - p + 0.32}]},
                Nothing
            ]
        },
        Row[{"(\[Theta], \[Phi]) = ", {t, f}, ", modes ", {p, q}}]
    ]
]

phasePrimitive[{x_, p_, a_}, labelQ_] := Tooltip[
    {
        EdgeForm[{GrayLevel[0.25], AbsoluteThickness[1]}], FaceForm[GrayLevel[0.95]],
        Rectangle[{x - 0.08, - p - 0.08}, {x + 0.08, - p + 0.08}],
        If[labelQ, {Black, Text[Style[parameterForm[a], 8], {x, - p - 0.3}]}, Nothing]
    },
    Row[{"\[Phi] = ", a, ", mode ", p}]
]

Options[interferometerDiagram] = {"ShowParameters" -> False, "ModeLabels" -> Automatic};

interferometerDiagram[qi_, opts : OptionsPattern[{interferometerDiagram, Graphics}]] := Module[{
    m = qi["Modes"], placed = Transpose[{qi["Elements"], qi["Layers"]}], cellPhases, couplers, phases,
    labelQ = TrueQ[OptionValue["ShowParameters"]], x0 = 0.25, x1
},
    x1 = Max[0, qi["Layers"]] + 0.75;
    (* a phase at k -/+ 1/4 belongs to the cell whose beam splitter is at column k *)
    cellPhases = Association @ Cases[placed, {"PS"[a_] -> {p_}, x_} /; ! IntegerQ[2 x] :> {Round[x], p} -> a];
    couplers = Cases[placed, {"BS"[t_, f_] -> {p_, q_}, x_} :>
        {x, Min[p, q], Max[p, q], t, Lookup[cellPhases, Key[{x, Min[p, q]}], f]}
    ];
    phases = Cases[placed, {"PS"[a_] -> {p_}, x_} /; IntegerQ[2 x] :> {x, p, a}];
    Graphics[
        {
            {GrayLevel[0.35], AbsoluteThickness[1.5], Table[Line[{{#1, - p}, {#2, - p}} & @@@ wireSegments[p, couplers, {x0, x1}]], {p, m}]},
            couplerPrimitive[#, labelQ] & /@ couplers,
            phasePrimitive[#, labelQ] & /@ phases,
            If[ Replace[OptionValue["ModeLabels"], Automatic -> m <= 30],
                Table[Text[Style[p, 10, GrayLevel[0.4]], {x0 - 0.1, - p}, {1, 0}], {p, m}],
                Nothing
            ]
        },
        FilterRules[{opts}, Options[Graphics]],
        PlotRange -> {{x0 - 0.45, x1 + 0.1}, {- m - 0.5, - 0.4}},
        ImageSize -> Clip[If[labelQ, 75, 40] (x1 - x0 + 0.55), {120, If[labelQ, 1200, 800]}]
    ]
]

interferometerMatrixPlot[qi_, opts : OptionsPattern[]] := With[{w = N[qi["Unitary"]]},
    Row[{
        ArrayPlot[Abs[w] ^ 2, opts,
            ColorFunction -> (Blend[{White, RGBColor[0.05, 0.25, 0.6]}, #] &), ColorFunctionScaling -> False,
            PlotLabel -> "|\!\(\*SubscriptBox[\(W\), \(jk\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)", Frame -> True, FrameTicks -> Automatic, Mesh -> All, MeshStyle -> GrayLevel[0.85], ImageSize -> 220
        ],
        ArrayPlot[Map[If[Abs[#] < 10^-12, White, phaseColor[Arg[#]]] &, w, {2}], opts,
            PlotLabel -> "Arg \!\(\*SubscriptBox[\(W\), \(jk\)]\)", Frame -> True, FrameTicks -> Automatic, Mesh -> All, MeshStyle -> GrayLevel[0.85], ImageSize -> 220
        ]
    }, "  "] /; MatrixQ[w, NumericQ]
]



(* ============================================================================ *)
(* Formatting                                                                   *)
(* ============================================================================ *)

interferometerIcon[qi_] := If[
    qi["Modes"] <= 16 && Length[qi["Elements"]] <= 400,
    interferometerDiagram[qi, "ModeLabels" -> False, ImageSize -> {36, 36}, AspectRatio -> 1],
    Framed[Style["\[ScriptCapitalI]", 16], FrameStyle -> GrayLevel[0.6]]
]

MakeBoxes[qi_QuantumInterferometer ? QuantumInterferometerQ, form_] ^:= BoxForm`ArrangeSummaryBox["QuantumInterferometer",
    qi,
    interferometerIcon[qi],
    {
        {BoxForm`SummaryItem[{"Modes: ", qi["Modes"]}], BoxForm`SummaryItem[{"Method: ", qi["Method"]}]},
        {BoxForm`SummaryItem[{"Beam splitters: ", qi["BeamSplitterCount"]}], BoxForm`SummaryItem[{"Depth: ", qi["Depth"]}]}
    },
    {
        {BoxForm`SummaryItem[{"Phase shifters: ", qi["PhaseShifterCount"]}]},
        If[qi["Modes"] <= 6, {BoxForm`SummaryItem[{"Unitary: ", MatrixForm[qi["Unitary"]]}]}, Nothing]
    },
    form
]
