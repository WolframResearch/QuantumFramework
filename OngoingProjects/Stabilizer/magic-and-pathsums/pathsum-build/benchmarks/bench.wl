(* ::Package:: *)

(* bench.wl : benchmark harness for the path-sum simulator. Honest, like-for-like.  *)
(* Assumes Wolfram`QuantumFramework`, pathsum.wl, harness.wl, elimination.wl,        *)
(* diagonal.wl loaded.                                                               *)
(*                                                                                  *)
(* WHAT IS AND IS NOT MEASURED. We compare only commensurable quantities. The        *)
(* engine's general single-amplitude path completes its residual by a brute-force    *)
(* Gauss sum of size 2^{residual}, so it is NOT a single-amplitude speed win by       *)
(* itself; the genuine measured speed win is the DIAGONAL (CNOT-dihedral) compressor  *)
(* diagAmp, which solves the affine output map and evaluates one closed-form term.    *)
(* The engine's wins are (a) COMPRESSION (residual-variable count) and (b) CAPABILITY *)
(* (a closed-form amplitude as a function of gate angles; exact Z[omega] arithmetic). *)
(* We also report where it LOSES: generic magic, where the residual stays large.      *)

SetAttributes[PathSum`bench, HoldFirst];
PathSum`bench[e_, k_ : 5] := Median[Table[ClearSystemCache[]; First @ AbsoluteTiming[e], {k}]];

(* ---- circuit families ---- *)
PathSum`hl[n_] := Table[{"H", q}, {q, n}];
(* diagonal (CNOT-dihedral) fragment on |+...+>: S, CZ, CCZ, CNOT, T, depth ~ c n *)
PathSum`diagFamily[n_, c_ : 3] := Join[
    Table[{"S", q}, {q, n}], Table[{"CZ", q, Mod[q, n] + 1}, {q, n}],
    Table[{"CCZ", q, Mod[q, n] + 1, Mod[q + 1, n] + 1}, {q, n}],
    Table[{"T", q}, {q, n}], Table[{"CNOT", q, Mod[q, n] + 1}, {q, n}]][[;; UpTo[c n]]];
(* Clifford-heavy ladder (H/S/CZ), telescopes to 0 internal residual *)
PathSum`ladder[n_, d_] := Join[PathSum`hl[n],
    Flatten[Table[Join[Table[{"S", q}, {q, n}], Table[{"CZ", q, q + 1}, {q, n - 1}], PathSum`hl[n]], {d}], 1]];
(* generic random magic (the honest-loss family) *)
PathSum`magicFamily[n_, d_] := PathSum`randomCircuit[n, d];

(* ---- metric 1: diagonal single-amplitude, diagAmp vs QF dense single entry ---- *)
PathSum`benchDiag[n_] := Module[{g = PathSum`diagFamily[n], y = ConstantArray[0, n], full},
    full = Join[PathSum`hl[n], g];
    {n, Length[g],
     PathSum`bench[PathSum`diagAmp[n, g, y], 5],                 (* engine (diagonal) *)
     PathSum`bench[PathSum`qfAmp[n, full, y], 5]}];              (* QF dense, one amplitude *)

(* ---- metric 2: engine compression on the Clifford ladder ---- *)
PathSum`benchLadder[n_, d_] := Module[{ps = PathSum`build[n, PathSum`ladder[n, d]], r, t},
    t = PathSum`bench[PathSum`engineReduce[ps], 3];
    r = PathSum`engineReduce[ps];
    {n, d, Length[ps["vars"]], Length[r["vars"]], Length[PathSum`internalResidual[r, n]], r["steps"], t}];

(* ---- metric 3: parametric capability: one symbolic traversal vs dense per-angle -- *)
(* engine: closed-form <Z>(theta) once; dense: rebuild per theta value.               *)
PathSum`benchParam[mvals_] := Module[{tEngine, tDense, zexpr},
    (* a 1-qubit H P[theta] H expectation, carried symbolically through the frame *)
    tEngine = PathSum`bench[
        zexpr = FullSimplify[Conjugate[#] . (PauliMatrix[3] . #) &[
            Normal[PauliStabilizer[1][{"H", "P"[\[Theta]], "H"}]["StateVector"]]], Assumptions -> \[Theta] \[Element] Reals];, 3];
    tDense = PathSum`bench[Table[Re[Conjugate[#] . (PauliMatrix[3] . #)] &[
        Normal[QuantumCircuitOperator[{"H" -> 1, "P"[t] -> 1, "H" -> 1}][QuantumState["0"], Method -> "Schrodinger"]["StateVector"]]],
        {t, mvals}], 3];
    {Length[mvals], tEngine, tDense, zexpr}];

(* ---- metric 4: generic-magic honest loss: engine vs dense single amplitude ---- *)
PathSum`benchMagic[n_, d_, seed_] := Module[{g, r},
    SeedRandom[seed]; g = PathSum`magicFamily[n, d];
    r = PathSum`engineReduce[PathSum`build[n, g]];
    {n, d, Length[PathSum`internalResidual[r, n]],
     PathSum`bench[PathSum`engineAmp[n, g, ConstantArray[0, n]], 3],
     PathSum`bench[PathSum`qfAmp[n, g, ConstantArray[0, n]], 3]}];

(* ---- metric 4': the honest loss. Narrow-but-deep magic H(TH)^k on ONE qubit:    *)
(* the residual is k magic variables, so the path sum brute-forces 2^k PATHS where  *)
(* dense is 2^1. When #paths > #basis-states, the path sum is the wrong tool.        *)
PathSum`narrowMagic[k_] := Join[{{"H", 1}}, Flatten[Table[{{"T", 1}, {"H", 1}}, k], 1]];
PathSum`benchNarrow[k_] := Module[{g = PathSum`narrowMagic[k], r},
    r = PathSum`engineReduce[PathSum`build[1, g]];
    {k, Length[PathSum`internalResidual[r, 1]],
     PathSum`bench[PathSum`engineAmp[1, g, {0}], 3],
     PathSum`bench[PathSum`qfAmp[1, g, {0}], 3]}];

(* ---- metric 5: representational, eager StabilizerFrame components vs compressor -- *)
PathSum`benchFrame[n_, k_] := Module[{arrows, frame, diag},
    (* k T-gates after an H-layer, no CCZ so QF builds a frame *)
    arrows = Join[Table["H" -> q, {q, n}], Table["T" -> Mod[i, n] + 1, {i, 0, k - 1}],
                  Table["CZ" -> {q, Mod[q, n] + 1}, {q, n}]];
    frame = PauliStabilizer[n][arrows];
    diag = PathSum`diagonalCompress[n,
        Join[Table[{"T", Mod[i, n] + 1}, {i, 0, k - 1}], Table[{"CZ", q, Mod[q, n] + 1}, {q, n}]]];
    {n, k, If[MatchQ[frame, _StabilizerFrame], frame["Length"], 1], diag["rank"]}];
