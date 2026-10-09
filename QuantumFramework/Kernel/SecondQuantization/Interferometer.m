(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["QuantumInterferometer"]

PackageScope["QuantumInterferometerQ"]



(* QuantumInterferometer: a passive linear-optical network on m modes.

   The primary data is the m x m transfer matrix W, with the convention
       W[[j, k]] = <1_j| U |1_k>,   i.e.   U a_k^dag U^dag = Sum_j W[[j,k]] a_j^dag
   so a gate g1 followed by g2 has W = W2 . W1. A Fock-space circuit is built
   from W only on request ("CircuitOperator"), and photon statistics come from
   permanents of submatrices of W, with no Fock-space truncation.

   Internal representation:
     QuantumInterferometer[<|"Unitary" -> W, "Modes" -> m, "Method" -> method,
         "Elements" -> {"PS"[a] -> {p}, "BS"[t, f] -> {p, q}, ...},
         "Layers" -> {x1, x2, ...}|>]
   Elements are time ordered, in QuantumCircuitOperator shorthand; Layers hold
   each element's horizontal position in the mesh: integer k for the beam
   splitter of column k and half integers for phases.

   Decompositions: Reck et al., PRL 73, 58 (1994); Clements et al., Optica 3,
   1460 (2016); Dhand and Goyal, PRA 92, 043813 (2015). Reck and Clements null W
   with Clements' cell T(t, f) of Eq. 1, a phase f on the upper mode followed by
   a real beam splitter BS(t, 0); the cosine-sine mesh recurses on the CS
   decomposition. Every mesh then moves its phases into the beam splitters' own
   phi, which leaves a single phase column at the output. *)


QuantumInterferometer::usage = "QuantumInterferometer[u] decomposes the m\[Times]m unitary u into a mesh of beam splitters and phase shifters.
QuantumInterferometer[u, Method -> method] uses the \"Clements\" (default), \"Reck\" or \"CosineSine\" mesh.
QuantumInterferometer[{\"BS\"[\[Theta], \[Phi]] -> {p, q}, \"PS\"[\[Alpha]] -> p, \[Ellipsis]}, m] builds an m-mode interferometer from beam splitters and phase shifters.
QuantumInterferometer[\"Random\"[m]] and QuantumInterferometer[\"Fourier\"[m]] give a Haar-random and a discrete Fourier interferometer.";

QuantumInterferometer::symbolic = "The matrix has non-numeric entries; only numeric (exact or approximate) unitaries can be decomposed. Build the interferometer from a list of beam splitters and phase shifters instead.";
QuantumInterferometer::nonunitary = "The matrix is not unitary within the tolerance `1`.";
QuantumInterferometer::method = "Method `1` is not one of \"Clements\", \"Reck\" or \"CosineSine\".";
QuantumInterferometer::elem = "`1` is not a valid element; use \"BS\"[\[Theta], \[Phi]] -> {p, q} or \"PS\"[\[Alpha]] -> p.";
QuantumInterferometer::modes = "The elements act on mode `1`, outside the `2` modes of the interferometer.";
QuantumInterferometer::occ = "`1` is not a list of `2` nonnegative photon numbers.";
QuantumInterferometer::levels = "\"CircuitOperator\" needs the number of Fock levels per mode, as in interferometer[\"CircuitOperator\", d]; n photons need d >= n + 1.";
QuantumInterferometer::state = "The state must have `1` modes with the same number of levels each.";
QuantumInterferometer::cutoff = "The state holds up to `1` photons but each mode keeps only `2` levels; at least `3` levels are needed for a faithful result.";
QuantumInterferometer::compose ="Cannot compose interferometers on `1` and `2` modes.";
QuantumInterferometer::plot = "\"MatrixPlot\" needs a numeric transfer matrix.";
QuantumInterferometer::undefprop = "property `` is undefined for this interferometer";

Options[QuantumInterferometer] = {Method -> Automatic, Tolerance -> 10.^-10, "DropIdentities" -> True};



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



(* Angles of the cell whose inverse, multiplied from the right on columns c and
   c + 1, zeroes a = u[[r, c]] against its right neighbour b = u[[r, c + 1]]. *)
rightNullAngles[a_, b_, tol_] := Which[
    zeroQ[a, tol], {0, 0},
    zeroQ[b, tol], {Pi / 2, 0},
    True, {ArcTan[Abs[b], Abs[a]], Arg[a Conjugate[b]]}
]

(* Angles of the cell that, multiplied from the left on rows r - 1 and r, zeroes
   b = u[[r, c]] against the entry a = u[[r - 1, c]] above it. *)
leftNullAngles[a_, b_, tol_] := Which[
    zeroQ[b, tol], {0, 0},
    zeroQ[a, tol], {Pi / 2, 0},
    True, {ArcTan[Abs[a], Abs[b]], Arg[- b Conjugate[a]]}
]


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

(* {right cells, diagonal, left cells} with L_k ... L_1 . W . R_1^-1 ... R_n^-1 = D *)
nullingSweep[w_, steps_, tol_, simp_] := With[{reaped = Reap[Fold[nullStep[#1, #2, tol, simp] &, w, steps], {"Right", "Left"}]},
    {Catenate[reaped[[2, 1]]], simp /@ Diagonal[reaped[[1]]], Catenate[reaped[[2, 2]]]}
]

(* Reck: the entries below the diagonal, row by row from the bottom, each against
   its right neighbour, so W = D . R_K ... R_1. *)
reckSteps[m_] := Catenate @ Table[{"Right", {r, c}}, {r, m, 2, -1}, {c, r - 1}]

(* Clements, supplementary material: along the i-th antidiagonal, from the right
   for odd i and from the left for even i, so W = L_1^-1 ... L_k^-1 . D . R_n ... R_1. *)
clementsSteps[m_] := Catenate @ Table[
    If[ OddQ[i],
        Table[{"Right", {m - j, i - j}}, {j, 0, i - 1}],
        Table[{"Left", {m + j - i, j}}, {j, i}]
    ],
    {i, m - 1}
]

(* A decomposition gives a time-ordered list of items: {"Cell", t, f, p},
   {"InverseCell", t, f, p}, {"BS", t, f, {p, q}} and {"PS", a, p}. *)
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

(* U = (L1 + L2) . CS . (R1^dag + R2^dag) on the first p = Ceiling[m/2] and the last q
   modes, CS coupling mode p - q + i with p + i by BS(t_i, 0), cos t_i the singular
   values of U11. L2 is the polar factor of U21 . R1, which is L2 . S, and
   R2^dag = C . L2^dag . U22 - S . (L1^dag . U12) needs no division by C or S. *)
csItems[{{z_}}, k_] := {{"PS", Arg[z], k + 1}}

csItems[u_, k_] := Module[{p = Ceiling[Length[u] / 2], q = Floor[Length[u] / 2], l1, sigma, r1, y, l2, c, s},
    {l1, sigma, r1} = SingularValueDecomposition[u[[;; p, ;; p]]];
    y = (u[[p + 1 ;;, ;; p]] . r1)[[All, p - q + 1 ;;]];
    l2 = #1 . ConjugateTranspose[#3] & @@ SingularValueDecomposition[y];
    c = Diagonal[sigma][[p - q + 1 ;;]];
    s = Norm /@ Transpose[y];
    Join[
        csItems[ConjugateTranspose[r1], k],
        csItems[c (ConjugateTranspose[l2] . u[[p + 1 ;;, p + 1 ;;]]) - s (ConjugateTranspose[l1] . u[[;; p, p + 1 ;;]])[[p - q + 1 ;;]], k + p],
        MapThread[{"BS", ArcTan[#1, #2], 0, {k + p - q + #3, k + p + #3}} &, {c, s, Range[q]}],
        csItems[l1, k],
        csItems[l2, k + p]
    ]
]

decompositionItems["CosineSine", w_, _, _] := csItems[N[w], 0]

(* BS(t, g) . diag(e^ia, e^ib) = diag(e^ia, e^ib) . BS(t, g + a - b): the phases met
   so far ride along to the output, and each beam splitter takes up their difference.
   BS(-t, g) = BS(t, g + Pi) keeps every t in [0, Pi/2]. *)
absorbItem[phases_, {"Cell", t_, f_, p_}, simp_] := With[{new = ReplacePart[phases, p -> phases[[p]] + f]},
    Sow[{"BS", t, simp @ normalAngle[new[[p]] - new[[p + 1]]], {p, p + 1}}];
    new
]

absorbItem[phases_, {"InverseCell", t_, f_, p_}, simp_] := (
    Sow[{"BS", t, simp @ normalAngle[phases[[p]] - phases[[p + 1]] + Pi], {p, p + 1}}];
    ReplacePart[phases, p -> phases[[p]] - f]
)

absorbItem[phases_, {"BS", t_, f_, {p_, q_}}, simp_] := (
    Sow[{"BS", t, simp @ normalAngle[f + phases[[p]] - phases[[q]]], {p, q}}];
    phases
)

absorbItem[phases_, {"PS", a_, p_}, _] := ReplacePart[phases, p -> phases[[p]] + a]

meshItems[items_, m_, simp_] := With[{reaped = Reap[Fold[absorbItem[#1, #2, simp] &, ConstantArray[0, m], items]]},
    Join[Catenate[reaped[[2]]], MapIndexed[{"PS", simp @ normalAngle[#1], First[#2]} &, reaped[[1]]]]
]



(* last holds the last occupied position of every mode: a beam splitter takes the
   first column after it on all the modes it spans, a phase the next free
   half-integer position. *)
placeItem[last_, {"BS", t_, f_, {p_, q_}}] := With[{span = Range[Min[p, q], Max[p, q]]},
    With[{k = 1 + Floor[Max[last[[span]]]]},
        {{{"BS"[t, f] -> {p, q}, k}}, ReplacePart[last, Thread[span -> k]]}
    ]
]

placeItem[last_, {"PS", a_, p_}] := With[{x = Floor[last[[p]] + 1/2] + 1/2}, {{{"PS"[a] -> {p}, x}}, ReplacePart[last, p -> x]}]

modeSpan[_[__] -> modes_] := Range @@ MinMax[modes]

(* Phases outside cells with no beam splitter after them on their mode shift
   together, so that the first one lands on a common column after the mesh. *)
alignTrailingPhases[placed_, m_] := With[{
    lastSplitter = Fold[
        ReplacePart[#1, Thread[modeSpan[#2[[2]]] -> #2[[1]]]] &,
        ConstantArray[0, m],
        Cases[MapIndexed[{First[#2], First[#1]} &, placed], {i_, gate : ("BS"[__] -> _)} :> {i, gate}]
    ],
    end = Max[0, Cases[placed, {"BS"[__] -> _, x_} :> x]] + 1/2
},
    With[{trailing = Cases[MapIndexed[{First[#2], #1} &, placed], {i_, {"PS"[_] -> {p_}, x_}} /; IntegerQ[x - 1/2] && i > lastSplitter[[p]] :> {i, p, x}]},
        With[{first = GroupBy[trailing, #[[2]] & -> Last, Min]},
            ReplacePart[placed, {#1, 2} -> end + #3 - first[#2] & @@@ trailing]
        ]
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
    If[! MemberQ[{"Clements", "Reck", "CosineSine"}, method], Message[QuantumInterferometer::method, method]; Return[$Failed]];
    If[! unitaryQ[w, tol], Message[QuantumInterferometer::nonunitary, tol]; Return[$Failed]];
    makeInterferometer[
        w, m, method,
        With[{simp = If[Precision[w] === Infinity, Simplify, Identity]}, meshItems[decompositionItems[method, w, tol, simp], m, simp]],
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

(* the exact decomposition of FourierMatrix[m] stalls from m = 5 on *)
QuantumInterferometer["Fourier"[m_Integer ? Positive], opts : OptionsPattern[]] :=
    QuantumInterferometer[N @ FourierMatrix[m], opts]

QuantumInterferometer[name : "Random" | "Fourier", m_Integer ? Positive, opts : OptionsPattern[]] :=
    QuantumInterferometer[name[m], opts]

QuantumInterferometer[qi_QuantumInterferometer ? QuantumInterferometerQ, opts : OptionsPattern[]] :=
    QuantumInterferometer[qi["Unitary"], opts]



$QuantumInterferometerProperties = {
    "Unitary", "TransferMatrix", "Modes", "Method", "Elements", "Layers",
    "Depth", "BeamSplitterCount", "PhaseShifterCount", "ElementTable",
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

QuantumInterferometerProp[qi_, "ElementTable"] := Dataset @ MapThread[
    <|"Element" -> Head[#1[[1]]], "Parameters" -> List @@ #1[[1]], "Modes" -> #1[[2]], "Position" -> #2|> &,
    {qi["Elements"], qi["Layers"]}
]

QuantumInterferometerProp[qi_, "Diagram", opts : OptionsPattern[]] := interferometerDiagram[qi, opts]

QuantumInterferometerProp[qi_, "MatrixPlot", opts : OptionsPattern[]] :=
    If[MatrixQ[qi["Unitary"], NumericQ], interferometerMatrixPlot[qi, opts], Message[QuantumInterferometer::plot]; $Failed]

QuantumInterferometerProp[qi_, "CircuitOperator", d_Integer ? Positive, opts : OptionsPattern[]] := interferometerCircuit[qi, d, opts]

QuantumInterferometerProp[qi_, "CircuitOperator"] := (Message[QuantumInterferometer::levels]; $Failed)

(* Photon patterns are occupation lists {n1, n2, ...} or Ket[{n1, n2, ...}];
   results are keyed by Ket so they read in Dirac notation. *)
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

inverseElement["BS"[t_, f_] -> modes_] := "BS"[t, normalAngle[f + Pi]] -> modes

inverseElement["PS"[a_] -> modes_] := "PS"[- a] -> modes

QuantumInterferometerProp[qi_, "Dagger"] := With[{last = Max[qi["Depth"], Ceiling[Max[0, qi["Layers"]] - 1/2]]},
    QuantumInterferometer[<|
        "Unitary" -> ConjugateTranspose[qi["Unitary"]],
        "Modes" -> qi["Modes"],
        "Method" -> None,
        "Elements" -> Reverse[inverseElement /@ qi["Elements"]],
        "Layers" -> Reverse[last + 1 - qi["Layers"]]
    |>]
]



occupation[Ket[s_List]] := s

occupation[s_] := s

occupationQ[qi_, s_] := If[
    VectorQ[s, IntegerQ[#] && # >= 0 &] && Length[s] == qi["Modes"],
    True,
    Message[QuantumInterferometer::occ, s, qi["Modes"]]; False
]

modeList[occupation_] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, occupation]

photonPatterns[n_, m_] := ReverseSort @ FrobeniusSolve[ConstantArray[1, m], n]

(* <T| U |S> = Perm(W_{T,S}) / Sqrt[Prod s! Prod t!] *)
permanentAmplitude[w_, s_, t_] := Which[
    Total[s] != Total[t], 0,
    Total[s] == 0, 1,
    True, Permanent[w[[modeList[t], modeList[s]]]] / Sqrt[Times @@ (s!) Times @@ (t!)]
]

(* U |S> = Prod_k (Sum_j W[[j, k]] a_j^dag)^s_k / Sqrt[s_k!] |0>, so the amplitude
   of T is its monomial's coefficient times Sqrt[Prod t! / Prod s!]. *)
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

fockComponents[qs_, d_, m_] := {IntegerDigits[#[[1, 1]] - 1, d, m], #[[2]]} & /@
    Most @ ArrayRules[SparseArray[qs["Computational"]["StateVector"]]]

photonSectorState[w_, components_, d_, m_] := With[{
    amplitudes = KeySelect[Merge[#2 photonAmplitudes[w, #1] & @@@ components, Total], Max[#] < d &]
},
    QuantumState[SparseArray[KeyValueMap[FromDigits[#1, d] + 1 -> #2 &, amplitudes], d ^ m], ConstantArray[d, m]]
]

(* A state vector goes through the permanents; a density matrix through the
   Fock-space circuit. Either way the result keeps the input's d levels per mode. *)
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



modePositions[qi_] := Merge[MapThread[Thread[modeSpan[#1] -> #2] &, {qi["Elements"], qi["Layers"]}], Identity]

(* qi2[qi1] is qi1 followed by qi2, shifted by whole columns until it clears qi1 on every shared mode *)
(qi2_QuantumInterferometer ? QuantumInterferometerQ)[qi1_QuantumInterferometer ? QuantumInterferometerQ] := If[
    qi1["Modes"] != qi2["Modes"],
    Message[QuantumInterferometer::compose, qi1["Modes"], qi2["Modes"]]; $Failed,
    QuantumInterferometer[<|
        "Unitary" -> qi2["Unitary"] . qi1["Unitary"],
        "Modes" -> qi1["Modes"],
        "Method" -> None,
        "Elements" -> Join[qi1["Elements"], qi2["Elements"]],
        "Layers" -> Join[
            qi1["Layers"],
            qi2["Layers"] + Max[0, Values @ Merge[KeyIntersection[{Max /@ modePositions[qi1], Min /@ modePositions[qi2]}], Floor[Subtract @@ #] + 1 &]]
        ]
    |>]
]



$couplerHalfWidth = 0.3;

reflectivityColor[t_] := If[NumericQ[t],
    Blend[{RGBColor[0.62, 0.74, 0.88], RGBColor[0.05, 0.25, 0.6]}, N[Sin[t] ^ 2]],
    GrayLevel[0.3]
]

phaseColor[a_] := If[NumericQ[a], Hue[Mod[N[a] / (2 Pi), 1], 0.7, 0.95], White]

parameterForm[v_] := If[NumericQ[v], NumberForm[Chop[N[v]], {3, 2}], v]

waveguide[x_, w_, y0_, y1_] := Line @ Table[
    {x + u w, y0 + (y1 - y0) (1 - Cos[Pi Min[1, 2 (1 - Abs[u])]]) / 2},
    {u, -1, 1, 1/20}
]

wireSegments[p_, couplers_, {x0_, x1_}] := With[{
    blocked = SortBy[Cases[couplers, {x_, p, _, _, _} | {x_, _, p, _, _} :> {x - $couplerHalfWidth, x + $couplerHalfWidth}], First]
},
    Partition[Flatten[{x0, blocked, x1}], 2]
]

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
    m = qi["Modes"], placed = Transpose[{qi["Elements"], qi["Layers"]}], couplers,
    labelQ = TrueQ[OptionValue["ShowParameters"]], x0 = 0.25, x1 = Max[0, qi["Layers"]] + 0.75
},
    couplers = Cases[placed, {"BS"[t_, f_] -> {p_, q_}, x_} :> {x, Min[p, q], Max[p, q], t, f}];
    Graphics[
        {
            {GrayLevel[0.35], AbsoluteThickness[1.5], Table[Line[{{#1, - p}, {#2, - p}} & @@@ wireSegments[p, couplers, {x0, x1}]], {p, m}]},
            couplerPrimitive[#, labelQ] & /@ couplers,
            phasePrimitive[#, labelQ] & /@ Cases[placed, {"PS"[a_] -> {p_}, x_} :> {x, p, a}],
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
    }, "  "]
]



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
