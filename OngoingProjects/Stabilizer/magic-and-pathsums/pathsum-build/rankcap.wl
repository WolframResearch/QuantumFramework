(* ::Package:: *)

(* rankcap.wl : exact span-dimension cap on a gate-built StabilizerFrame.         *)
(*                                                                              *)
(* A gate-built StabilizerFrame stores 2^(#T) components but they can be          *)
(* linearly dependent; the stabilizer rank is at most 2^n for any T-count. This   *)
(* caps the eager component list at the EXACT span dimension over Z[omega]:        *)
(*   |s> = Sum_i c_i u_i,  u_i = (relating Pauli_i) . |reference>,                 *)
(* keep a maximal linearly-independent subset {u_j} (always including the          *)
(* reference, component 1) and re-express |s> = Sum_j d_j u_j by exact LinearSolve.*)
(*                                                                              *)
(* This is a SPAN-DIMENSION cap (an upper-rank witness, <= 2^n), NOT the           *)
(* NP-hard optimal stabilizer rank. Standalone: uses only the public               *)
(* StabilizerFrame API; no kernel edit. Assumes QF + pathsum.wl loaded.            *)

(* Apply a relating Pauli {xbits, zbits, coeff} to a dense vector (qubit 1 = MSB); *)
(* mirrors the kernel's sfApplyPauli so the reconstruction is phase-coherent.       *)
PathSum`applyPauli[{xb_, zb_, coeff_}, v_, n_] := Module[
    {dim = 2^n, xmask = FromDigits[xb, 2], iPow = Total[xb zb], targets, signs},
    targets = 1 + BitXor[Range[0, dim - 1], xmask];
    signs = Table[(-1)^Mod[zb . IntegerDigits[b, 2, n], 2], {b, 0, dim - 1}];
    Normal @ SparseArray[Thread[targets -> (coeff I^iPow) signs Normal[v]], dim]
];

(* coherent component vectors u_i and coefficients c_i; |s> = Sum c_i u_i. *)
PathSum`coherentComponents[f_StabilizerFrame] := Module[{assoc = First[f], ref, n, v1},
    If[! KeyExistsQ[assoc, "Paulis"], Return[$Failed]];
    ref = assoc["Components"][[1, 2]]; n = ref["Qubits"];
    v1 = Normal[ref["State"]["StateVector"]];
    <|"u" -> (PathSum`applyPauli[#, v1, n] & /@ assoc["Paulis"]),
      "c" -> assoc["Components"][[All, 1]],
      "ps" -> assoc["Components"][[All, 2]], "paulis" -> assoc["Paulis"]|>
];

(* exact span-dimension cap. Returns a reduced StabilizerFrame reproducing |s>. *)
PathSum`capFrameRank[f_StabilizerFrame] := Module[
    {cc = PathSum`coherentComponents[f], us, full, basis, d},
    If[cc === $Failed, Return[f]];
    us = RootReduce[cc["u"]];
    full = RootReduce[Total[cc["c"] us]];                       (* = f["StateVector"] *)
    (* greedy maximal independent subset, reference (component 1) added first *)
    basis = {}; Module[{acc = {}},
        Do[If[MatrixRank[N[Append[acc, us[[i]]]]] > Length[acc],
              AppendTo[acc, us[[i]]]; AppendTo[basis, i]], {i, Length[us]}]];
    d = LinearSolve[Transpose[us[[basis]]], full];              (* |s> = Sum_j d_j u_{basis_j} *)
    StabilizerFrame[<|
        "Components" -> Table[{d[[j]], cc["ps"][[basis[[j]]]]}, {j, Length[basis]}],
        "Paulis" -> cc["paulis"][[basis]]
    |>]
];

(* convenience: {eager count, capped count} *)
PathSum`rankCapCounts[f_StabilizerFrame] := {f["Length"], PathSum`capFrameRank[f]["Length"]};
