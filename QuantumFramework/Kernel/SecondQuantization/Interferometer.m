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
   phi, which leaves a single phase column at the output.

   Squeezers "S"[xi] -> {p}, S(xi) = Exp[(xi^* a^2 - xi a^dag^2) / 2] as SqueezeOperator,
   make the network active: "Unitary" is None and "Bogoliubov" holds {alpha, beta, c}
   with U a_k^dag U^dag = Sum_j (alpha_jk a_j^dag + beta_jk a_j) and c = <0|U|0>.
   Passive networks carry no "Bogoliubov" key and keep the permanent paths. *)


QuantumInterferometer::usage = "QuantumInterferometer[u] decomposes the m\[Times]m unitary u into a mesh of beam splitters and phase shifters.
QuantumInterferometer[u, Method -> method] uses the \"Clements\" (default), \"Reck\" or \"CosineSine\" mesh.
QuantumInterferometer[{\"BS\"[\[Theta], \[Phi]] -> {p, q}, \"PS\"[\[Alpha]] -> p, \[Ellipsis]}, m] builds an m-mode interferometer from beam splitters and phase shifters.
QuantumInterferometer[\"Random\"[m]] and QuantumInterferometer[\"Fourier\"[m]] give a Haar-random and a discrete Fourier interferometer.
Elements \"S\"[\[Xi]] -> p add squeezers; their photon distributions take \"Cutoff\" -> n, the largest total photon number kept.
\"Amplitude\" and \"Probability\" of a passive interferometer accept symbolic photon numbers, as in {n1, n2} -> {t1, t2}, and give nested sums of Wigner d-functions.";

QuantumInterferometer::symbolic = "The matrix has non-numeric entries; only numeric (exact or approximate) unitaries can be decomposed. Build the interferometer from a list of beam splitters and phase shifters instead.";
QuantumInterferometer::nonunitary = "The matrix is not unitary within the tolerance `1`.";
QuantumInterferometer::method = "Method `1` is not one of \"Clements\", \"Reck\" or \"CosineSine\".";
QuantumInterferometer::elem = "`1` is not a valid element; use \"BS\"[\[Theta], \[Phi]] -> {p, q}, \"PS\"[\[Alpha]] -> p or \"S\"[\[Xi]] -> p.";
QuantumInterferometer::modes = "The elements act on mode `1`, outside the `2` modes of the interferometer.";
QuantumInterferometer::occ = "`1` is not a list of `2` nonnegative photon numbers.";
QuantumInterferometer::levels = "\"CircuitOperator\" needs the number of Fock levels per mode, as in interferometer[\"CircuitOperator\", d]; n photons need d >= n + 1.";
QuantumInterferometer::state = "The state must have `1` modes with the same number of levels each.";
QuantumInterferometer::truncation = "The state holds up to `1` photons but each mode keeps only `2` levels; at least `3` levels are needed for a faithful result.";
QuantumInterferometer::compose ="Cannot compose interferometers on `1` and `2` modes.";
QuantumInterferometer::plot = "\"MatrixPlot\" needs a numeric transfer matrix.";
QuantumInterferometer::sample = "Sampling needs a numeric transfer matrix.";
QuantumInterferometer::squeezed = "`1` is not defined for an interferometer with squeezers.";
QuantumInterferometer::autocutoff = "Patterns with up to `1` photons hold probability `2`; set \"Cutoff\" -> n to keep more.";
QuantumInterferometer::symcutoff = "\"Cutoff\" -> Automatic needs numeric squeezing; set \"Cutoff\" -> n.";
QuantumInterferometer::undefprop = "property `` is undefined for this interferometer";

Options[QuantumInterferometer] = {Method -> Automatic, Tolerance -> 10.^-10, "DropIdentities" -> True};



(* BeamSplitterOperator[{t, f}] on modes p, q *)
bsBlock[t_, f_] := {{Cos[t], - Exp[- I f] Sin[t]}, {Exp[I f] Sin[t], Cos[t]}}

(* Clements' cell T(t, f) = BS(t, 0) . PS(f on the upper mode), and its inverse for real t *)
cellBlock[t_, f_] := {{Exp[I f] Cos[t], - Sin[t]}, {Exp[I f] Sin[t], Cos[t]}}

cellInverseBlock[t_, f_] := {{Exp[- I f] Cos[t], Exp[- I f] Sin[t]}, {- Sin[t], Cos[t]}}

(* The nulling cells in closed rational form in the entries themselves, equal to
   cellInverseBlock @@ rightNullAngles[a, b] and cellBlock @@ leftNullAngles[a, b] but
   without the Cos[ArcTan[..]] / Exp[I Arg[..]] round-trip, which on exact input swells
   into nested radicals. The zero branches match the angle functions' zeroQ branches. *)
rightInverseCell[a_, b_, tol_] := Which[
    zeroQ[a, tol], IdentityMatrix[2],
    zeroQ[b, tol], {{0, 1}, {- 1, 0}},
    True, With[{A = Abs[a], B = Abs[b], r = Sqrt[Abs[a] ^ 2 + Abs[b] ^ 2]},
        {{Conjugate[a] b / (A r), Conjugate[a] b / (B r)}, {- A / r, B / r}}
    ]
]

leftCell[a_, b_, tol_] := Which[
    zeroQ[b, tol], IdentityMatrix[2],
    zeroQ[a, tol], {{0, - 1}, {1, 0}},
    True, With[{A = Abs[a], B = Abs[b], r = Sqrt[Abs[a] ^ 2 + Abs[b] ^ 2]},
        {{- b Conjugate[a] / (B r), - B / r}, {- b Conjugate[a] / (A r), A / r}}
    ]
]

gateMatrix["BS"[t_, f_] -> {p_, q_}, m_] := ReplacePart[IdentityMatrix[m], Thread[Tuples[{p, q}, 2] -> Flatten[bsBlock[t, f]]]]

gateMatrix["PS"[a_] -> {p_}, m_] := ReplacePart[IdentityMatrix[m], {p, p} -> Exp[I a]]

(* Apply one gate as a two-row update of the running matrix: gate . w changes only rows
   p, q, so this is O(m) per gate instead of the O(m^2) of a full m x m product. The result
   is repacked once at the end, since the per-step ReplacePart unpacks a packed matrix. *)
applyGate[w_, "BS"[t_, f_] -> {p_, q_}] := With[{r = bsBlock[t, f] . w[[{p, q}]]},
    ReplacePart[w, {p -> r[[1]], q -> r[[2]]}]
]

applyGate[w_, "PS"[a_] -> {p_}] := ReplacePart[w, p -> Exp[I a] w[[p]]]

elementsUnitary[elements_List, m_Integer] := Developer`ToPackedArray[Fold[applyGate, IdentityMatrix[m], elements]]


zeroQ[x_, tol_] := If[Precision[x] === Infinity, PossibleZeroQ[x], Abs[x] <= tol]

(* Dropping no-op gates is cosmetic, so its identity test is a numeric tolerance check:
   message-free on exact high-degree angles that no symbolic zero-test can rigorously
   decide. The exact zero detection that steers the decomposition stays in zeroQ. *)
identityGateQ["BS"[t_, _] -> _, tol_] := TrueQ[Abs[N[t]] <= tol]

identityGateQ["PS"[a_] -> _, tol_] := TrueQ[Abs[N[Exp[I a] - 1]] <= tol]

identityGateQ["S"[xi_] -> _, tol_] := TrueQ[Abs[N[xi]] <= tol]

unitaryQ[w_, tol_] := If[Precision[w] === Infinity,
    AllTrue[Flatten[w . ConjugateTranspose[w] - IdentityMatrix[Length[w]]], PossibleZeroQ],
    UnitaryMatrixQ[w, Tolerance -> tol]
]

normalAngle[x_] := If[NumericQ[x] && Precision[x] =!= Infinity, Mod[x, 2 Pi, - Pi], x]


passiveQ[QuantumInterferometer[data_]] := data["Unitary"] =!= None

bogoliubov[element : (("BS"[__] | "PS"[_]) -> _), m_] := {gateMatrix[element, m], ConstantArray[0, {m, m}]}

bogoliubov["S"[xi_] -> {p_}, m_] := {
    ReplacePart[IdentityMatrix[m], {p, p} -> Cosh[Abs[xi]]],
    ReplacePart[ConstantArray[0, {m, m}], {p, p} -> Exp[- I Arg[xi]] Sinh[Abs[xi]]]
}

vacuumMatrix[a_, b_] := - Inverse[ConjugateTranspose[a]] . ConjugateTranspose[b]

(* {alpha, beta} compose as (a2, b2) o (a1, b1) = (a2 a1 + b2^* b1, b2 a1 + a2^* b1). A squeezer
   on p acting on c Exp[a^dag . B . a^dag / 2] |0> leaves the vacuum amplitude
   c Sqrt[Sech[r]] (1 - tau^* B_pp)^(-1/2), tau = e^(i phi) Tanh[r]; |tau|, |B_pp| < 1 keep the root on its principal branch. *)
gaussianData[elements_, m_] := Fold[
    With[{a = #1[[1]], b = #1[[2]], c = #1[[3]], e = bogoliubov[#2, m]},
        {
            e[[1]] . a + Conjugate[e[[2]]] . b,
            e[[2]] . a + Conjugate[e[[1]]] . b,
            Replace[#2, {("S"[xi_] -> {p_}) :> c Sqrt[Sech[Abs[xi]]] (1 - Exp[- I Arg[xi]] Tanh[Abs[xi]] vacuumMatrix[a, b][[p, p]]) ^ (-1/2), _ :> c}]
        }
    ] &,
    {IdentityMatrix[m], ConstantArray[0, {m, m}], 1},
    elements
]

(* <0| e^(x.a) U e^(y.a^dag) |0> = c Exp[z.K.z / 2], z = (x, y), so
   <t| U |s> = c Haf(K_{t+s}) / Sqrt[t! s!] with output rows first *)
kMatrix[{a_, b_, _}] := With[{bv = vacuumMatrix[a, b]},
    ArrayFlatten[{{bv, a + bv . b}, {Transpose[a + bv . b], (Transpose[a] . b + Transpose[b] . a) / 2 + Transpose[b] . bv . b}}]
]

gaussianAmplitude[{k_, c_}, s_, t_] := c hafnian[k[[#, #]] &[Join[modeList[t], Length[s] + modeList[s]]]] / Sqrt[Times @@ (s!) Times @@ (t!)]

gaussianKernel[QuantumInterferometer[data_]] := {kMatrix[data["Bogoliubov"]], data["Bogoliubov"][[3]]}


hafnian[{}] := 1

hafnian[a_ /; OddQ[Length[a]]] := 0

hafnian[a_] := Which[
    Length[a] <= 8 || Precision[a] === Infinity, hafnianPairing[a],
    Precision[a] === MachinePrecision, With[{m = Length[a] / 2},
        With[{e = Round[Subdivide[0, 2 ^ m - 1, Min[2 ^ m - 1, 64]]]},
            Total[hafnianTraceTerms[N[a[[All, Join[Range[m + 1, 2 m], Range[m]]]]], Most[e] + 1, Rest[e]]]
        ]
    ],
    True, hafnianTrace[a]
]

hafnianPairing[{}] := 1

hafnianPairing[a_] := Sum[a[[1, j]] hafnianPairing[a[[#, #]] &[Complement[Range[2, Length[a]], {j}]]], {j, 2, Length[a]}]

(* Bjorklund, Gupt and Quesada, arXiv:1805.12498: Haf(A) = Sum_Z (-1)^(n/2 - |Z|) f((A X)_Z), f(C) the
   eta^(n/2) coefficient of Exp[Sum_k Tr[C^k] eta^k / (2 k)] *)
hafnianTrace[a_] := With[{m = Length[a] / 2},
    With[{b = a[[All, Join[Range[m + 1, 2 m], Range[m]]]]},
        Total[(-1) ^ (m - Length[#]) traceCoefficient[b[[Join[#, # + m], Join[#, # + m]]], m] & /@ Subsets[Range[m]]]
    ]
]

traceCoefficient[{}, _] := 0

traceCoefficient[c_, m_] := With[{t = Tr /@ NestList[c . # &, c, m - 1]},
    Last @ Fold[Append[#1, Total[t[[;; #2]] Reverse[#1]] / (2 #2)] &, {1}, Range[m]]
]


hafnianTraceTerms := hafnianTraceTerms = Compile[{{b, _Complex, 2}, {lo, _Integer}, {hi, _Integer}},
    Module[{m = Quotient[Length[b], 2], tot = 0. + 0. I, idx, c, p, t, f, s},
        Do[
            idx = Flatten[Position[IntegerDigits[mask, 2, m], 1]];
            idx = Join[idx, idx + m];
            c = b[[idx, idx]];
            p = c;
            t = Table[0. + 0. I, {m}];
            Do[s = 0. + 0. I; Do[s += p[[i, i]], {i, Length[p]}]; t[[k]] = s; If[k < m, p = c . p], {k, m}];
            f = Table[0. + 0. I, {m + 1}];
            f[[1]] = 1. + 0. I;
            Do[s = 0. + 0. I; Do[s += t[[i]] f[[j - i + 1]], {i, j}]; f[[j + 1]] = s / (2 j), {j, m}];
            tot += If[EvenQ[m - Quotient[Length[idx], 2]], 1., -1.] f[[m + 1]],
            {mask, lo, hi}
        ];
        tot
    ],
    RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> "Speed"
]



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


(* nullStep eliminates one entry of the working matrix v, mutating it in place. The
   sweep owns v, so after the first assignment v has a single reference and the
   two-column (Right) or two-row (Left) update writes in place instead of copying the
   whole m x m matrix each step; HoldFirst keeps v the sweep's own symbol rather than
   a copy of its value. The eliminated entry is set to 0. v[[r, c]], a zero of the
   entry's own type, which leaves a machine matrix packed (an exact 0 would unpack a
   packed array); on the exact path (exact = True) it is set to exact 0. *)
SetAttributes[nullStep, HoldFirst]

nullStep[v_, {"Right", {r_, c_}}, tol_, simp_, asimp_, exact_] := With[{a = v[[r, c]], b = v[[r, c + 1]]},
    Sow[Append[asimp /@ rightNullAngles[a, b, tol], c], "Right"];
    v[[All, {c, c + 1}]] = simp[v[[All, {c, c + 1}]] . rightInverseCell[a, b, tol]];
    v[[r, c]] = If[exact, 0, 0. v[[r, c]]]
]

nullStep[v_, {"Left", {r_, c_}}, tol_, simp_, asimp_, exact_] := With[{a = v[[r - 1, c]], b = v[[r, c]]},
    Sow[Append[asimp /@ leftNullAngles[a, b, tol], r - 1], "Left"];
    v[[{r - 1, r}]] = simp[leftCell[a, b, tol] . v[[{r - 1, r}]]];
    v[[r, c]] = If[exact, 0, 0. v[[r, c]]]
]

(* {right cells, diagonal, left cells} with L_k ... L_1 . W . R_1^-1 ... R_n^-1 = D.
   v is owned here and mutated in place by nullStep across the whole sweep. *)
nullingSweep[w_, steps_, tol_, simp_, asimp_] := Module[{v = w, exact = Precision[w] === Infinity, reaped},
    reaped = Reap[Scan[nullStep[v, #, tol, simp, asimp, exact] &, steps], {"Right", "Left"}];
    {Catenate[reaped[[2, 1]]], simp /@ Diagonal[v], Catenate[reaped[[2, 2]]]}
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
decompositionItems["Reck", w_, tol_, simp_, asimp_] := With[{dec = nullingSweep[w, reckSteps[Length[w]], tol, simp, asimp]},
    Join[{"Cell", ##} & @@@ dec[[1]], MapIndexed[{"PS", asimp[Arg[#1]], First[#2]} &, dec[[2]]]]
]

decompositionItems["Clements", w_, tol_, simp_, asimp_] := With[{dec = nullingSweep[w, clementsSteps[Length[w]], tol, simp, asimp]},
    Join[
        {"Cell", ##} & @@@ dec[[1]],
        MapIndexed[{"PS", asimp[Arg[#1]], First[#2]} &, dec[[2]]],
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

decompositionItems["CosineSine", w_, _, _, _] := csItems[N[w], 0]

(* BS(t, g) . diag(e^ia, e^ib) = diag(e^ia, e^ib) . BS(t, g + a - b): the phases met
   so far ride along to the output, and each beam splitter takes up their difference.
   BS(-t, g) = BS(t, g + Pi) keeps every t in [0, Pi/2]. meshItems owns the phase
   vector and absorbStep mutates it in place (HoldFirst), so no length-m vector is
   rebuilt per item. *)
SetAttributes[absorbStep, HoldFirst]

absorbStep[phases_, {"Cell", t_, f_, p_}, simp_] := (
    phases[[p]] += f;
    Sow[{"BS", t, simp @ normalAngle[phases[[p]] - phases[[p + 1]]], {p, p + 1}}]
)

absorbStep[phases_, {"InverseCell", t_, f_, p_}, simp_] := (
    Sow[{"BS", t, simp @ normalAngle[phases[[p]] - phases[[p + 1]] + Pi], {p, p + 1}}];
    phases[[p]] -= f
)

absorbStep[phases_, {"BS", t_, f_, {p_, q_}}, simp_] :=
    Sow[{"BS", t, simp @ normalAngle[f + phases[[p]] - phases[[q]]], {p, q}}]

absorbStep[phases_, {"PS", a_, p_}, _] := phases[[p]] += a

meshItems[items_, m_, simp_] := Module[{phases = ConstantArray[0, m], reaped},
    reaped = Reap[Scan[absorbStep[phases, #, simp] &, items]];
    Join[Catenate[reaped[[2]]], MapIndexed[{"PS", simp @ normalAngle[#1], First[#2]} &, phases]]
]



(* last holds the last occupied position of every mode: a beam splitter takes the
   first column after it on all the modes it spans, a phase the next free
   half-integer position. layoutItems owns last and placeStep mutates it in place
   (HoldFirst), Sowing each placed element, so no length-m vector is rebuilt per item. *)
SetAttributes[placeStep, HoldFirst]

placeStep[last_, {"BS", t_, f_, {p_, q_}}] := With[{span = Range[Min[p, q], Max[p, q]]},
    With[{k = 1 + Floor[Max[last[[span]]]]},
        last[[span]] = k;
        Sow[{"BS"[t, f] -> {p, q}, k}]
    ]
]

placeStep[last_, {kind : "PS" | "S", a_, p_}] := With[{x = Floor[last[[p]] + 1/2] + 1/2},
    last[[p]] = x;
    Sow[{kind[a] -> {p}, x}]
]

modeSpan[_[__] -> modes_] := Range @@ MinMax[modes]

(* Phases outside cells with no beam splitter after them on their mode shift
   together, so that the first one lands on a common column after the mesh. *)
alignTrailingPhases[placed_, m_] := With[{
    lastSplitter = Fold[
        ReplacePart[#1, Thread[modeSpan[#2[[2]]] -> #2[[1]]]] &,
        ConstantArray[0, m],
        Cases[MapIndexed[{First[#2], First[#1]} &, placed], {i_, gate : (("BS"[__] | "S"[_]) -> _)} :> {i, gate}]
    ],
    end = Max[0, Cases[placed, {"BS"[__] -> _, x_} :> x]] + 1/2
},
    With[{trailing = Cases[MapIndexed[{First[#2], #1} &, placed], {i_, {"PS"[_] -> {p_}, x_}} /; IntegerQ[x - 1/2] && i > lastSplitter[[p]] :> {i, p, x}]},
        With[{first = GroupBy[trailing, #[[2]] & -> Last, Min]},
            ReplacePart[placed, {#1, 2} -> end + #3 - first[#2] & @@@ trailing]
        ]
    ]
]

layoutItems[items_List, m_Integer, drop_, tol_] := Module[{last = ConstantArray[0, m], reaped},
    reaped = Reap[Scan[placeStep[last, #] &, items]];
    With[{placed = alignTrailingPhases[If[reaped[[2]] === {}, {}, Catenate[reaped[[2]]]], m]},
        With[{kept = If[TrueQ[drop], DeleteCases[placed, {gate_, _} /; identityGateQ[gate, tol]], placed]},
            If[kept === {}, {{}, {}}, Transpose[kept]]
        ]
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
    "Unitary" -> _ ? SquareMatrixQ | None, "Modes" -> _Integer, "Method" -> _, "Elements" -> _List, "Layers" -> _List
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
        With[{simp = If[Precision[w] === Infinity, RootReduce, Identity], asimp = If[Precision[w] === Infinity, Simplify, Identity]}, meshItems[decompositionItems[method, w, tol, simp, asimp], m, asimp]],
        OptionValue["DropIdentities"], tol
    ]
]


normalizeElement[("BS" | "BeamSplitter")[{t_, f_}] -> modes_] := normalizeElement["BS"[t, f] -> modes]

normalizeElement[("BS" | "BeamSplitter")[t_, f_ : 0] -> {p_Integer ? Positive, q_Integer ? Positive}] /; p != q := {"BS", t, f, {p, q}}

normalizeElement[("PS" | "PhaseShift")[a_] -> (p_Integer ? Positive | {p_Integer ? Positive})] := {"PS", a, p}

normalizeElement[("S" | "Squeeze")[xi_] -> (p_Integer ? Positive | {p_Integer ? Positive})] := {"S", xi, p}

normalizeElement[element_] := (Message[QuantumInterferometer::elem, element]; $Failed)

QuantumInterferometer[elements : {___Rule}, m_Integer ? Positive, opts : OptionsPattern[]] := Module[{items = normalizeElement /@ elements, layout},
    If[MemberQ[items, $Failed], Return[$Failed]];
    With[{top = Max[0, Flatten[items[[All, -1]]]]},
        If[top > m, Message[QuantumInterferometer::modes, top, m]; Return[$Failed]]
    ];
    layout = layoutItems[items, m, OptionValue["DropIdentities"], OptionValue[Tolerance]];
    elementsInterferometer[layout[[1]], layout[[2]], m, None]
]

elementsInterferometer[elements_, layers_, m_, unitary_] := With[{passive = FreeQ[elements, "S"[_] -> _]},
    QuantumInterferometer @ Join[
        <|
            "Unitary" -> If[passive, Replace[unitary, None :> elementsUnitary[elements, m]], None],
            "Modes" -> m,
            "Method" -> None,
            "Elements" -> elements,
            "Layers" -> layers
        |>,
        If[passive, <||>, <|"Bogoliubov" -> gaussianData[elements, m]|>]
    ]
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
    If[passiveQ[qi], QuantumInterferometer[qi["Unitary"], opts], Message[QuantumInterferometer::squeezed, "Re-meshing"]; $Failed]



$QuantumInterferometerProperties = {
    "Unitary", "TransferMatrix", "Modes", "Method", "Elements", "Layers",
    "Depth", "BeamSplitterCount", "PhaseShifterCount", "SqueezerCount", "ElementTable",
    "Diagram", "MatrixPlot", "CircuitOperator",
    "Amplitude", "Probability", "Probabilities", "Samples", "SampleCounts", "Dagger"
};

(qi_QuantumInterferometer ? QuantumInterferometerQ)[prop_String, args___] := With[{result = QuantumInterferometerProp[qi, prop, args]},
    result /; ! MatchQ[Unevaluated @ result, _QuantumInterferometerProp] || Message[QuantumInterferometer::undefprop, prop]
]

QuantumInterferometerProp[_, "Properties"] := $QuantumInterferometerProperties

QuantumInterferometerProp[QuantumInterferometer[data_], key : "Modes" | "Method" | "Elements" | "Layers"] := data[key]

QuantumInterferometerProp[qi : QuantumInterferometer[data_], "Unitary"] :=
    If[passiveQ[qi], data["Unitary"], Message[QuantumInterferometer::squeezed, "\"Unitary\""]; $Failed]

QuantumInterferometerProp[qi_, "TransferMatrix"] := qi["Unitary"]

QuantumInterferometerProp[qi_, "Depth"] := Max[0, Pick[qi["Layers"], MatchQ["BS"[__] -> _] /@ qi["Elements"]]]

QuantumInterferometerProp[qi_, "BeamSplitterCount"] := Count[qi["Elements"], "BS"[__] -> _]

QuantumInterferometerProp[qi_, "PhaseShifterCount"] := Count[qi["Elements"], "PS"[_] -> _]

QuantumInterferometerProp[qi_, "SqueezerCount"] := Count[qi["Elements"], "S"[_] -> _]

QuantumInterferometerProp[qi_, "ElementTable"] := Dataset @ Map[
    <|"Element" -> Head[#[[1]]], "Parameters" -> List @@ #[[1]], "Modes" -> #[[2]]|> &,
    qi["Elements"]
]

QuantumInterferometerProp[qi_, "Diagram", opts : OptionsPattern[]] := interferometerDiagram[qi, opts]

QuantumInterferometerProp[qi_, "MatrixPlot", opts : OptionsPattern[]] := Which[
    ! passiveQ[qi], Message[QuantumInterferometer::squeezed, "\"MatrixPlot\""]; $Failed,
    MatrixQ[qi["Unitary"], NumericQ], interferometerMatrixPlot[qi, opts],
    True, Message[QuantumInterferometer::plot]; $Failed
]

QuantumInterferometerProp[qi_, "CircuitOperator", d_Integer ? Positive] := interferometerCircuit[qi, d]

QuantumInterferometerProp[qi_, "CircuitOperator"] := (Message[QuantumInterferometer::levels]; $Failed)

(* Photon patterns are occupation lists {n1, n2, ...}, Ket[{n1, n2, ...}] or QuditName[{n1, n2, ...}];
   results are keyed by QuditName, as QuantumState["Probabilities"] is. *)
fockAmplitude[qi_, s_, t_] := If[passiveQ[qi], permanentAmplitude[qi["Unitary"], s, t], gaussianAmplitude[gaussianKernel[qi], s, t]]

QuantumInterferometerProp[qi_, "Amplitude", s_ -> t_] := With[{in = occupation[s], out = occupation[t]},
    Which[
        symbolicPatternsQ[qi, in, out], If[passiveQ[qi], symbolicAmplitude[qi["Elements"], in, out], Message[QuantumInterferometer::squeezed, "Symbolic photon numbers"]; $Failed],
        occupationQ[qi, in] && occupationQ[qi, out], fockAmplitude[qi, in, out],
        True, $Failed
    ]
]

QuantumInterferometerProp[qi_, "Probability", s_ -> t_] := With[{in = occupation[s], out = occupation[t]},
    Which[
        symbolicPatternsQ[qi, in, out], Replace[qi["Amplitude", in -> out], a : Except[$Failed] :> Abs[a] ^ 2],
        occupationQ[qi, in] && occupationQ[qi, out], Abs[fockAmplitude[qi, in, out]] ^ 2,
        True, $Failed
    ]
]

QuantumInterferometerProp[qi_, "Probabilities", s_, opts : OptionsPattern[{"Cutoff" -> Automatic}]] := With[{in = occupation[s]},
    Which[
        ! occupationQ[qi, in], $Failed,
        passiveQ[qi], With[{amplitudes = photonAmplitudes[qi["Unitary"], in]},
            Association @ Map[QuditName[#] -> Abs[Lookup[amplitudes, Key[#], 0]] ^ 2 &, photonPatterns[Total[in], qi["Modes"]]]
        ],
        True, Replace[gaussianAmplitudes[qi, in, OptionValue["Cutoff"]], a_Association :> Abs[a] ^ 2]
    ]
]

(* Clifford and Clifford, arXiv:1706.01260, Algorithm B: with the input photons in random order, photon k
   lands on mode i with weight |Perm(A[r ++ {i}, 1..k])|^2, expanded along row i into k minors of A[r, 1..k] *)
QuantumInterferometerProp[qi_, "Samples", s_, k_Integer ? NonNegative] := With[{in = occupation[s]},
    Which[
        ! occupationQ[qi, in], $Failed,
        ! passiveQ[qi], Message[QuantumInterferometer::squeezed, "Sampling"]; $Failed,
        ! MatrixQ[N[qi["Unitary"]], NumericQ], Message[QuantumInterferometer::sample]; $Failed,
        True, With[{w = N[qi["Unitary"]][[All, modeList[in]]], m = qi["Modes"], n = Total[in]},
            Table[
                With[{a = w[[All, RandomSample[Range[n]]]]},
                    QuditName @ Lookup[
                        Counts @ Fold[
                            Append[#1, RandomChoice[Abs[a[[All, ;; #2]] . If[#2 == 1, {1}, Table[Permanent[a[[#1, Delete[Range[#2], l]]]], {l, #2}]]] ^ 2 -> Range[m]]] &,
                            {},
                            Range[n]
                        ],
                        Range[m],
                        0
                    ]
                ],
                k
            ]
        ]
    ]
]

QuantumInterferometerProp[qi_, "SampleCounts", s_, k_Integer ? NonNegative] := Replace[qi["Samples", s, k],
    samples_List :> KeyMap[QuditName, Reverse @ KeySort @ Counts[#["Name"] & /@ samples]]
]

inverseElement["BS"[t_, f_] -> modes_] := "BS"[t, normalAngle[f + Pi]] -> modes

inverseElement[(kind : "PS" | "S")[a_] -> modes_] := kind[- a] -> modes

QuantumInterferometerProp[qi_, "Dagger"] := With[{last = Max[qi["Depth"], Ceiling[Max[0, qi["Layers"]] - 1/2]]},
    elementsInterferometer[
        Reverse[inverseElement /@ qi["Elements"]], Reverse[last + 1 - qi["Layers"]], qi["Modes"],
        If[passiveQ[qi], ConjugateTranspose[qi["Unitary"]], None]
    ]
]



occupation[Ket[s_List]] := s

occupation[s_QuditName] := s["Name"]

occupation[s_] := s

occupationQ[qi_, s_] := If[
    VectorQ[s, IntegerQ[#] && # >= 0 &] && Length[s] == qi["Modes"],
    True,
    Message[QuantumInterferometer::occ, s, qi["Modes"]]; False
]

symbolicPatternsQ[qi_, s_, t_] := VectorQ[s] && VectorQ[t] && Length[s] == Length[t] == qi["Modes"] &&
    ! VectorQ[Join[s, t], IntegerQ] && VectorQ[Join[s, t], ! NumericQ[#] || IntegerQ[#] && # >= 0 &]

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


wignerd[b_, j_ ? NumericQ, m1_ ? NumericQ, m2_ ? NumericQ] := If[j < 0 || Abs[m1] > j || Abs[m2] > j, 0, WignerD[{j, m1, m2}, b]]

wignerd /: MakeBoxes[wignerd[b_, j_, m1_, m2_], form_] := RowBox[{
    SubsuperscriptBox["d", RowBox[{MakeBoxes[m1, form], ",", MakeBoxes[m2, form]}], MakeBoxes[j, form]], "(", MakeBoxes[b, form], ")"
}]

symbolicCounts[c_, {"BS"[__] -> {p_, q_}, k_}] := ReplacePart[c, {p -> k, q -> c[[p]] + c[[q]] - k}]

symbolicCounts[c_, _] := c

symbolicFactor[c_, {"BS"[th_, ph_] -> {p_, q_}, k_}] :=
    (-1) ^ (c[[p]] - k) E^(I ph (c[[p]] - k)) wignerd[2 th, (c[[p]] + c[[q]]) / 2, k - (c[[p]] + c[[q]]) / 2, (c[[p]] - c[[q]]) / 2]

symbolicFactor[c_, {"PS"[al_] -> {p_}, _}] := E^(I al c[[p]])

symbolicAmplitude[elements_, n_, t_] := Module[{steps, counts, ks, eqs, rules, free, cs},
    steps = MapIndexed[{#1, Subscript[\[FormalK], First[#2]]} &, elements];
    counts = FoldList[symbolicCounts, n, steps];
    ks = Cases[steps, {"BS"[__] -> _, k_} :> k];
    eqs = Thread[Most[Last[counts]] == Most[t]];
    rules = If[ks === {}, {},
        With[{order = Reverse[ks]},
            First[Solve[eqs, Cases[RowReduce[Outer[Coefficient, Subtract @@@ eqs, order]], row_ /; AnyTrue[row, # != 0 &] :> order[[FirstPosition[row, x_ /; x != 0][[1]]]]], MaxExtraConditions -> All], {}] /.
                ConditionalExpression[e_, _] :> e
        ]
    ];
    free = DeleteCases[ks, Alternatives @@ Keys[rules]];
    cs = Select[DeleteDuplicates[Expand[Flatten[Rest[counts]] /. rules]], ! NumericQ[#] &];
    With[{angles = Select[Cases[elements, (("BS"[_, a_] | "PS"[a_]) -> _) :> a], ! NumericQ[#] &]},
        KroneckerDelta[Total[n], Total[t]] (Eliminate[eqs, ks] /. {Equal -> KroneckerDelta, And -> Times, True -> 1, False -> 0}) *
        Inactive[Sum][
            Times @@ MapThread[symbolicFactor, {Most[counts], steps}] /. rules,
            Sequence @@ Table[
                With[{k = free[[j]], c = Select[cs, ! FreeQ[#, free[[j]]] && FreeQ[#, Alternatives @@ Drop[free, j]] &]},
                    {k, Max[0, Sequence @@ (k - Select[c, Coefficient[#, k] == 1 &])], Min @@ (k + Select[c, Coefficient[#, k] == -1 &])}
                ],
                {j, Length[free]}
            ]
        ] //. {
            Inactive[Sum][b_, pre___, {k_, lo_, hi_} /; PossibleZeroQ[lo - hi], post___] :> Inactive[Sum][b /. k -> lo, pre, Sequence @@ ({post} /. k -> lo)],
            Inactive[Sum][b_] :> b
        } /. {
            Power[-1, x_ /; ! NumericQ[x]] :> (-1) ^ With[{v = Variables[x]}, Mod[x /. Thread[v -> 0], 2] + Mod[Coefficient[x, v], 2] . v],
            Power[E, x_ /; ! NumericQ[x]] :> E ^ (I Chop[Collect[Expand[x / I], angles]]),
            wignerd[b_, j_, m1_, m2_] :> If[PossibleZeroQ[j], 1, wignerd[b, Together[Expand[j]], Together[Expand[m1]], Together[Expand[m2]]]]
        }
    ]
]

shellAmplitudes[g_, s_, n_] := AssociationMap[gaussianAmplitude[g, s, #] &, photonPatterns[n, Length[s]]]

(* squeezers change the photon number in pairs, so the output shells share the input's parity *)
gaussianAmplitudes[qi_, s_, cutoff_] := With[{g = gaussianKernel[qi], start = Mod[Total[s], 2]},
    Which[
        IntegerQ[cutoff], KeyMap[QuditName, Association[shellAmplitudes[g, s, #] & /@ Range[start, cutoff, 2]]],
        MatrixQ[N[g[[1]]], NumericQ], KeyMap[QuditName, autoShells[g, s, start]],
        True, Message[QuantumInterferometer::symcutoff]; $Failed
    ]
]

autoShells[g_, s_, start_] := Module[{n = start, amplitudes = shellAmplitudes[g, s, start], held},
    held := N[Total[Abs[Values[amplitudes]] ^ 2]];
    While[held < 1 - 10^-4 && n < 16, n += 2; amplitudes = Join[amplitudes, shellAmplitudes[g, s, n]]];
    If[held < 1 - 10^-4, Message[QuantumInterferometer::autocutoff, n, held]];
    amplitudes
]

(qi_QuantumInterferometer ? QuantumInterferometerQ)[s : _List | Ket[_List] | _QuditName, opts : OptionsPattern[{"Cutoff" -> Automatic}]] /;
    VectorQ[occupation[s], IntegerQ] := With[{in = occupation[s]},
    Which[
        ! occupationQ[qi, in], $Failed,
        passiveQ[qi], KeyMap[QuditName, photonAmplitudes[qi["Unitary"], in]],
        True, gaussianAmplitudes[qi, in, OptionValue["Cutoff"]]
    ]
]



elementOperator["BS"[t_, f_] -> {p_, q_}, d_] := BeamSplitterOperator[{t, f}, d, {p, q}, Method -> "Recurrence"]

elementOperator["PS"[a_] -> {p_}, d_] := PhaseShiftOperator[a, d, {p}]

elementOperator["S"[xi_] -> {p_}, d_] := SqueezeOperator[xi, d, {p}]

interferometerCircuit[qi_, d_] := With[{
    idle = Complement[Range[qi["Modes"]], Flatten[qi["Elements"][[All, 2]]]]
},
    QuantumCircuitOperator @ Join[
        QuantumOperator[IdentityMatrix[d, SparseArray], {#}, d] & /@ idle,
        elementOperator[#, d] & /@ qi["Elements"]
    ]
]

fockComponents[qs_, d_, m_] := {IntegerDigits[#[[1, 1]] - 1, d, m], #[[2]]} & /@
    Most @ ArrayRules[SparseArray[qs["Computational"]["StateVector"]]]

photonSectorState[w_, components_, d_, m_] := With[{
    amplitudes = KeySelect[Merge[#2 photonAmplitudes[w, #1] & @@@ components, Total], Max[#] < d &]
},
    QuantumState[SparseArray[KeyValueMap[FromDigits[#1, d] + 1 -> #2 &, amplitudes], d ^ m], ConstantArray[d, m]]
]

(* A state vector goes through the permanents; a density matrix, or any state when
   there are squeezers, through the Fock-space circuit. Either way the result keeps the
   input's d levels per mode. *)
(qi_QuantumInterferometer ? QuantumInterferometerQ)[qs_QuantumState] := Module[{
    dims = qs["Dimensions"], m = qi["Modes"], d, components, photons
},
    If[Length[dims] != m || ! Equal @@ dims, Message[QuantumInterferometer::state, m]; Return[$Failed]];
    d = First[dims];
    components = If[passiveQ[qi] && qs["StateType"] === "Vector", fockComponents[qs, d, m], None];
    photons = If[ components === None,
        Max[0, Total[IntegerDigits[# - 1, d, m]] & /@
            Flatten @ Position[Normal[qs["ProbabilitiesList"]], p_ /; ! PossibleZeroQ[p], {1}, Heads -> False]],
        Max[0, Total /@ components[[All, 1]]]
    ];
    If[photons >= d, Message[QuantumInterferometer::truncation, photons, d, photons + 1]];
    If[ components === None,
        qi["CircuitOperator", d][qs],
        photonSectorState[qi["Unitary"], components, d, m]
    ]
]



modePositions[qi_] := Merge[MapThread[Thread[modeSpan[#1] -> #2] &, {qi["Elements"], qi["Layers"]}], Identity]

(* qi2[qi1] is qi1 followed by qi2, shifted by whole columns until it clears qi1 on every shared mode *)
(qi2_QuantumInterferometer ? QuantumInterferometerQ)[qi1_QuantumInterferometer ? QuantumInterferometerQ] := If[
    qi1["Modes"] != qi2["Modes"],
    Message[QuantumInterferometer::compose, qi1["Modes"], qi2["Modes"]]; $Failed,
    elementsInterferometer[
        Join[qi1["Elements"], qi2["Elements"]],
        Join[
            qi1["Layers"],
            qi2["Layers"] + Max[0, Values @ Merge[KeyIntersection[{Max /@ modePositions[qi1], Min /@ modePositions[qi2]}], Floor[Subtract @@ #] + 1 &]]
        ],
        qi1["Modes"],
        If[passiveQ[qi1] && passiveQ[qi2], qi2["Unitary"] . qi1["Unitary"], None]
    ]
]



$beamSplitterHalfWidth = 0.25;

reflectivityColor[t_] := If[NumericQ[t],
    Blend[{RGBColor[0.62, 0.74, 0.88], RGBColor[0.05, 0.25, 0.6]}, N[Sin[t] ^ 2]],
    GrayLevel[0.3]
]

phaseColor[a_] := If[NumericQ[a], Hue[Mod[N[a] / (2 Pi), 1], 0.7, 0.95], White]

parameterForm[v_] := If[NumericQ[v], NumberForm[Chop[N[v]], {3, 2}], v]

wireSegments[p_, splitters_, {x0_, x1_}] := With[{
    blocked = SortBy[Cases[splitters, {x_, p, _, _, _} | {x_, _, p, _, _} :> {x - $beamSplitterHalfWidth, x + $beamSplitterHalfWidth}], First]
},
    Partition[Flatten[{x0, blocked, x1}], 2]
]

beamSplitterPrimitive[{x_, p_, q_, t_, f_}, labelQ_] := With[{w = $beamSplitterHalfWidth, mid = - (p + q) / 2},
    Tooltip[
        {
            reflectivityColor[t], AbsoluteThickness[2], CapForm["Round"],
            Line[{{{x - w, - p}, {x + w, - q}}, {{x - w, - q}, {x + w, - p}}}],
            {GrayLevel[0.15], AbsoluteThickness[3], CapForm["Butt"], Line[{{x - 0.12, mid}, {x + 0.12, mid}}]},
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
        Rectangle[{x - 0.11, - p - 0.11}, {x + 0.11, - p + 0.11}],
        {Black, Text[Style["\[Phi]", 8], {x, - p}]},
        If[labelQ, {Black, Text[Style[parameterForm[a], 8], {x, - p - 0.3}]}, Nothing]
    },
    Row[{"\[Phi] = ", a, ", mode ", p}]
]

squeezerPrimitive[{x_, p_, xi_}, labelQ_] := Tooltip[
    {
        EdgeForm[{GrayLevel[0.25], AbsoluteThickness[1]}], FaceForm[RGBColor[0.96, 0.72, 0.38]],
        Rectangle[{x - 0.11, - p - 0.11}, {x + 0.11, - p + 0.11}],
        {Black, Text[Style["S", 7], {x, - p}]},
        If[labelQ, {Black, Text[Style[parameterForm[xi], 8], {x, - p - 0.3}]}, Nothing]
    },
    Row[{"\[Xi] = ", xi, ", mode ", p}]
]

Options[interferometerDiagram] = {"ShowParameters" -> False, "ModeLabels" -> Automatic};

interferometerDiagram[qi_, opts : OptionsPattern[{interferometerDiagram, Graphics}]] := Module[{
    m = qi["Modes"], placed = Transpose[{qi["Elements"], qi["Layers"]}], splitters,
    labelQ = TrueQ[OptionValue["ShowParameters"]], x0 = 0.25, x1 = Max[0, qi["Layers"]] + 0.75
},
    splitters = Cases[placed, {"BS"[t_, f_] -> {p_, q_}, x_} :> {x, Min[p, q], Max[p, q], t, f}];
    Graphics[
        {
            {GrayLevel[0.35], AbsoluteThickness[1.5], Table[Line[{{#1, - p}, {#2, - p}} & @@@ wireSegments[p, splitters, {x0, x1}]], {p, m}]},
            beamSplitterPrimitive[#, labelQ] & /@ splitters,
            phasePrimitive[#, labelQ] & /@ Cases[placed, {"PS"[a_] -> {p_}, x_} :> {x, p, a}],
            squeezerPrimitive[#, labelQ] & /@ Cases[placed, {"S"[xi_] -> {p_}, x_} :> {x, p, xi}],
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
        {BoxForm`SummaryItem[{"Phase shifters: ", qi["PhaseShifterCount"]}], BoxForm`SummaryItem[{"Squeezers: ", qi["SqueezerCount"]}]},
        If[passiveQ[qi] && qi["Modes"] <= 6, {BoxForm`SummaryItem[{"Unitary: ", MatrixForm[qi["Unitary"]]}]}, Nothing]
    },
    form
]
