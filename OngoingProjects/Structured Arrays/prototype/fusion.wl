(* Gate-class fusion prototype for QuantumFramework circuits on qubits.

   A circuit becomes a list of steps applied to a state vector of n qubits
   (qubit 1 is the most significant bit of the amplitude index):
     {"Dense", wires, matrix}     a gate applied along its qubit axes
     {"Monomial", data}           a fused run of diagonal / permutation / monomial gates
     {"Fourier", wires, sign}     a QFT (sign 1) or inverse QFT (sign -1) block

   A monomial map sends v to f * v[[p]]: p is a gather list (0-based source
   index of each amplitude) and f the factor on each amplitude. Two
   representations of the same fused run are kept, so they can be compared:
     vector form:     <|"p" -> p, "f" -> f|>
     structured form: {P, D} with P a PermutationMatrix, D a DiagonalMatrix, map = D . P *)

ts = TargetStructure -> "Structured";

(* ---------- gates from a QF circuit ---------- *)

gateRecord[op_] := Enclose[
    ConfirmAssert[op["Input"]["ComputationalQ"] && op["Output"]["ComputationalQ"]];
    ConfirmAssert[op["InputOrder"] === op["OutputOrder"]];
    ConfirmAssert[MatchQ[op["Dimensions"], {2 ..}]];
    (* QF indexes a gate's "Matrix" by its qubits in increasing order, whatever order "InputOrder" lists them in *)
    (* machine precision: the prototype targets numeric state vectors *)
    {Sort[op["InputOrder"]], N[Normal[op["Matrix"]]]}
];

(* the same circuit without its measurements, whose result qudits would widen QF's output state;
   a circuit with no measurements is returned as it is, keeping its label *)
withoutMeasurements[qc_] := If[FreeQ[qc["Elements"], _QuantumMeasurementOperator, {1}], qc,
    QuantumCircuitOperator[DeleteCases[qc["Elements"], _QuantumMeasurementOperator]]];

isQFTBlock[e_QuantumCircuitOperator] := MatchQ[e["Label"], "QFT" | SuperDagger["QFT"]];
isQFTBlock[_] := False;

(* top-level elements -> {"Gate", wires, matrix} | {"Fourier", wires, sign}; measurements are dropped *)
circuitItems[qc_, routeQFT_] := Catenate[circuitItem[#, routeQFT] & /@ qc["Elements"]];
circuitItem[e_QuantumCircuitOperator, True] /; isQFTBlock[e] := {{"Fourier", Sort[e["InputOrder"]], If[e["Label"] === "QFT", 1, -1]}};
circuitItem[e_QuantumCircuitOperator, routeQFT_] := circuitItems[e, routeQFT];
circuitItem[e_QuantumOperator, _] := {Prepend[gateRecord[e], "Gate"]};
circuitItem[_QuantumMeasurementOperator, _] := {};

(* ---------- gate classes ---------- *)

gateClass[m_] := With[{nz = Unitize[Chop[m]]},
    Which[
        nz == DiagonalMatrix[Diagonal[nz]], "Diagonal",
        Total[nz] == ConstantArray[1, Length[m]] && Total[nz, {2}] == ConstantArray[1, Length[m]],
            If[AllTrue[Pick[Flatten[m], Flatten[nz], 1], # == 1 &], "Permutation", "Monomial"],
        True, "Dense"
    ]
];

(* ---------- embedding a monomial gate in the full register ---------- *)

(* bits of every basis state, and the local index of the gate's wires in each *)
basisBits[n_] := basisBits[n] = Tuples[{0, 1}, n];
localIndex[w_, n_] := basisBits[n][[All, w]] . (2^Range[Length[w] - 1, 0, -1]);

(* v'[y] = f[y] v[src(y)]: row r of the local matrix has one nonzero, in column c(r) *)
monomialFull[w_, m_, n_] := Module[{k = Length[w], cols, vals, weights, delta, loc},
    cols = (First @ FirstPosition[#, x_ /; x != 0, {1}, Heads -> False]) - 1 & /@ m;
    vals = MapThread[#1[[#2 + 1]] &, {m, cols}];
    weights = 2^(n - w);
    delta = MapThread[(IntegerDigits[#2, 2, k] - IntegerDigits[#1, 2, k]) . weights &, {Range[0, 2^k - 1], cols}];
    loc = localIndex[w, n] + 1;
    <|"p" -> Developer`ToPackedArray[Range[0, 2^n - 1] + delta[[loc]]], "f" -> Developer`ToPackedArray[N[vals[[loc]]]]|>
];

(* fusing: first a, then b *)
composeVec[a_, b_] := <|"p" -> a["p"][[b["p"] + 1]], "f" -> b["f"] a["f"][[b["p"] + 1]]|>;
applyVec[mono_, v_] := mono["f"] v[[mono["p"] + 1]];

(* structured form: map = D . P with P . v == v[[p + 1]] *)
toStructured[mono_] := {PermutationMatrix[mono["p"] + 1, ts], DiagonalMatrix[mono["f"], ts]};
(* first {P1, D1}, then {P2, D2}: D2 . P2 . D1 . P1 = D2 . (P2 D1 P2^-1) . P2 . P1, and P2 D1 P2^-1 is D1's diagonal read through P2 *)
composeStructured[{p1_, d1_}, {p2_, d2_}] := {p2 . p1, d2 . DiagonalMatrix[Diagonal[d1][[p2 . Range[Length[d1]]]], ts]};
applyStructured[{pm_, dm_}, v_] := dm . (pm . v);

(* ---------- dense gates and Fourier blocks along qubit axes ---------- *)

axisOrder[w_, n_] := Join[w, Complement[Range[n], w]];
applyAlong[v_, w_, n_, f_] := With[{sig = axisOrder[w, n], k = Length[w]},
    Flatten[Transpose[
        ArrayReshape[f[ArrayReshape[Transpose[ArrayReshape[v, ConstantArray[2, n]], Ordering[sig]], {2^k, 2^(n - k)}]], ConstantArray[2, n]],
        sig]]
];
applyDense[v_, w_, m_, n_] := applyAlong[v, w, n, m . # &];
applyFourier[v_, w_, sign_, n_] := applyAlong[v, w, n, Transpose[If[sign == 1, Fourier, InverseFourier] /@ Transpose[#]] &];

(* ---------- compiling a circuit into steps ---------- *)

(* mode "Vector" or "Structured": the form in which monomial runs are fused *)
compileSteps[qc_, n_, mode_, routeQFT_] := Module[{items, runs},
    items = If[routeQFT && isQFTBlock[qc], {{"Fourier", Sort[qc["InputOrder"]], If[qc["Label"] === "QFT", 1, -1]}}, circuitItems[qc, routeQFT]];
    runs = SplitBy[items, #[[1]] === "Gate" && gateClass[#[[3]]] =!= "Dense" &];
    Catenate[compileRun[#, n, mode] & /@ runs]
];
compileRun[run_, n_, mode_] /; run[[1, 1]] === "Gate" && gateClass[run[[1, 3]]] =!= "Dense" :=
    With[{monos = monomialFull[#[[2]], #[[3]], n] & /@ run},
        {{"Monomial", If[mode === "Structured", Fold[composeStructured, toStructured /@ monos], Fold[composeVec, monos]], Length[run]}}];
compileRun[run_, n_, mode_] := Map[If[#[[1]] === "Fourier", #, {"Dense", #[[2]], #[[3]]}] &, run];

applySteps[steps_, v0_, n_, mode_] := Fold[
    Function[{v, s}, Switch[s[[1]],
        "Dense", applyDense[v, s[[2]], s[[3]], n],
        "Fourier", applyFourier[v, s[[2]], s[[3]], n],
        "Monomial", If[mode === "Structured", applyStructured[s[[2]], v], applyVec[s[[2]], v]]]],
    v0, steps];

(* gate-by-gate baseline: every gate applied along its axes, no classes, no fusion *)
applyGateByGate[qc_, v0_, n_] := Fold[applyDense[#1, #2[[2]], #2[[3]], n] &, v0, circuitItems[qc, False]];

(* share of gates absorbed into monomial runs *)
fusedShare[qc_] := With[{items = circuitItems[qc, False]}, N[Count[items, g_ /; gateClass[g[[3]]] =!= "Dense"] / Length[items]]];
