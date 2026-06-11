Package["Wolfram`QuantumFramework`"]

PackageScope[$stabilizerGateCodes]
PackageScope[$stabilizerCompileTarget]
PackageScope[$bitsPerWord]
PackageScope[compiledCliffordFold]
PackageScope[encodeStabilizerGates]
PackageScope[applyCompiledFold]
PackageScope[packGenRows]
PackageScope[unpackGenRows]
PackageScope[stabilizerValidMasks]



(* ============================================================================ *)
(* Compiled bulk Clifford gate-fold kernel.                                      *)
(*                                                                              *)
(* The interpreted packed path (Stabilizer/Packed.m) removed the O(n^2)-per-gate *)
(* algorithmic cost, but each gate still pays ~30-50 us of Wolfram-kernel        *)
(* dispatch. This kernel folds an entire integer-encoded gate list over the      *)
(* tableau in one compiled (C) call.                                             *)
(*                                                                              *)
(* It uses a GENERATOR-MAJOR packing (distinct from Packed.m's qubit-major one): *)
(* for each qubit a, the X (resp. Z) bits across all 2n generators are packed    *)
(* into ceil(2n/62) machine words, and the 2n sign bits are one packed word      *)
(* vector. A single-qubit gate then touches the O(2n/62) words of its qubit's    *)
(* column -- the layout that makes Clifford updates word-parallel (Stim / QC.jl  *)
(* use the same idea). The sign update is one vectorized XOR over the sign words. *)
(* Validated bit-for-bit against the interpreted packed path (PauliStabilizer.   *)
(* wlt TIER 10).                                                                  *)
(* Reference: external-packages-audit \[Section]4; "Need for Speed" hints 21-24. *)
(* ============================================================================ *)


(* Single source of truth for the integer gate encoding. The compiled body and  *)
(* the encoder both read these codes (injected into the Compile body via With),  *)
(* so the mapping is never duplicated as literals. CZ is expanded into H/CNOT/H  *)
(* by the encoder; CX is an alias of CNOT.                                       *)
$stabilizerGateCodes = <|
    "H" -> 0, "S" -> 1, SuperDagger["S"] -> 2, "CNOT" -> 3,
    "SWAP" -> 4, "X" -> 5, "Y" -> 6, "Z" -> 7
|>

$bitsPerWord = 62
$packPowers := $packPowers = Developer`ToPackedArray[2^Range[0, $bitsPerWord - 1]]
$bitSelectors := $bitSelectors = Developer`ToPackedArray[2^Range[0, $bitsPerWord - 1]]


(* Compilation target. "C" when a C compiler is installed (closest to Stim),     *)
(* "WVM" bytecode otherwise -- still free of per-gate interpreter dispatch.       *)
$stabilizerCompileTarget := $stabilizerCompileTarget = (
    Needs["CCompilerDriver`"];
    If[Length[CCompilerDriver`CCompilers[]] > 0, "C", "WVM"]
)


(* ---- generator-major pack / unpack (vectorized, packed-array fast) ---- *)

(* {rows, len}-bit matrix -> {rows, #words} machine-int words.                   *)
(* The word values reach 2^62, so the Dot result comes back UNPACKED; force it    *)
(* back to a packed Int64 array, otherwise the compiled kernel silently falls     *)
(* back to interpreted evaluation on an unpacked argument.                        *)
packGenRows[bm_, nw_] := Developer`ToPackedArray @ Dot[
    ArrayReshape[
        PadRight[Developer`ToPackedArray[bm], {Length[bm], nw $bitsPerWord}],
        {Length[bm], nw, $bitsPerWord}
    ],
    $packPowers
]

(* {rows, #words} words -> {rows, len}-bit matrix. *)
unpackGenRows[wm_, len_] := Take[
    ArrayReshape[
        Transpose[Table[Unitize[BitAnd[wm, $bitSelectors[[j]]]], {j, $bitsPerWord}], {3, 1, 2}],
        {Length[wm], Length[First[wm]] $bitsPerWord}
    ],
    All, len
]

(* per-word "all valid bits" masks for a 2n-generator register. *)
stabilizerValidMasks[n_Integer, nw_Integer] :=
    Developer`ToPackedArray[(BitShiftLeft[1, #] - 1 &) /@ Append[ConstantArray[$bitsPerWord, nw - 1], Mod[2 n - 1, $bitsPerWord] + 1]]


(* ============================================================================ *)
(* The compiled kernel (lazily built once, then memoized).                      *)
(*                                                                              *)
(* Inputs:  x0, z0  -- {n, #words} generator-major X / Z word matrices          *)
(*          s0      -- length-#words sign word vector (phase bits over gens)     *)
(*          valid   -- length-#words per-word valid-bit masks                    *)
(*          gates   -- {m, 3} integer matrix, rows {code, a, b} (0-based qubits) *)
(* Output:  Join[x, z, {s}] -- a {2n + 1, #words} integer matrix.               *)
(* ============================================================================ *)

compiledCliffordFold := compiledCliffordFold = With[{
    cH = $stabilizerGateCodes["H"], cS = $stabilizerGateCodes["S"],
    cSdg = $stabilizerGateCodes[SuperDagger["S"]], cCNOT = $stabilizerGateCodes["CNOT"],
    cSWAP = $stabilizerGateCodes["SWAP"], cX = $stabilizerGateCodes["X"],
    cY = $stabilizerGateCodes["Y"], cZ = $stabilizerGateCodes["Z"]
},
    Compile[{{x0, _Integer, 2}, {z0, _Integer, 2}, {s0, _Integer, 1}, {valid, _Integer, 1}, {gates, _Integer, 2}},
        Module[{x = x0, z = z0, s = s0, code = 0, a = 1, b = 1, xa = s0, za = s0, xb = s0, zb = s0},
            Do[
                code = gates[[g, 1]]; a = gates[[g, 2]] + 1; b = gates[[g, 3]] + 1;
                Which[
                    code == cH,
                        xa = x[[a]]; za = z[[a]];
                        s = BitXor[s, BitAnd[xa, za]];
                        x[[a]] = za; z[[a]] = xa,
                    code == cS,
                        xa = x[[a]]; za = z[[a]];
                        s = BitXor[s, BitAnd[xa, za]];
                        z[[a]] = BitXor[za, xa],
                    code == cSdg,
                        xa = x[[a]]; za = z[[a]];
                        s = BitXor[s, BitAnd[xa, BitXor[za, valid]]];
                        z[[a]] = BitXor[za, xa],
                    code == cCNOT,
                        xa = x[[a]]; xb = x[[b]]; za = z[[a]]; zb = z[[b]];
                        s = BitXor[s, BitAnd[BitAnd[xa, zb], BitXor[BitXor[xb, za], valid]]];
                        x[[b]] = BitXor[xb, xa]; z[[a]] = BitXor[za, zb],
                    code == cSWAP,
                        xa = x[[a]]; x[[a]] = x[[b]]; x[[b]] = xa;
                        za = z[[a]]; z[[a]] = z[[b]]; z[[b]] = za,
                    code == cX, s = BitXor[s, z[[a]]],
                    code == cY, s = BitXor[s, BitXor[x[[a]], z[[a]]]],
                    code == cZ, s = BitXor[s, x[[a]]]
                ],
                {g, 1, Length[gates]}
            ];
            Join[x, z, {s}]
        ],
        CompilationTarget :> $stabilizerCompileTarget, RuntimeOptions -> "Speed"
    ]
]


(* ============================================================================ *)
(* Gate encoder: gate spec(s) -> {code, j, k} integer rows (0-based qubits).     *)
(* Returns $Failed for any gate outside the compiled set, so the caller can fall *)
(* back to the interpreted fold. CZ expands to H_k . CNOT . H_k; CX aliases      *)
(* CNOT. Accepts the `op -> order` spec shape used by QuantumCircuitOperator.    *)
(* ============================================================================ *)

(* One flat Dispatch table, applied in a single Replace pass over the spec list  *)
(* (~5x faster than layered per-spec function dispatch). All codes are injected   *)
(* from $stabilizerGateCodes via With, so the encoding stays single-sourced.      *)
(* Recognized spec shapes:                                                        *)
(*   op -> q  /  op -> {q}                  single-qubit (QCO emits the {q} form) *)
(*   op -> {j, k}                           two-qubit; CZ expands to H.CNOT.H     *)
(*   "C"["NOT"|"X"|"Z" -> {t}, {c}, ___]    QuantumShortcut controlled form       *)
(*   {op, q...}                             legacy list form                      *)
(* Anything else maps to the $Failed sentinel and the whole encode returns        *)
(* $Failed, sending the caller to the per-gate fold (non-Clifford etc.).          *)
$stabilizerEncodeDispatch := $stabilizerEncodeDispatch = With[{
    cH = $stabilizerGateCodes["H"], cS = $stabilizerGateCodes["S"],
    cSdg = $stabilizerGateCodes[SuperDagger["S"]], cCNOT = $stabilizerGateCodes["CNOT"],
    cSWAP = $stabilizerGateCodes["SWAP"], cX = $stabilizerGateCodes["X"],
    cY = $stabilizerGateCodes["Y"], cZ = $stabilizerGateCodes["Z"]
},
    Dispatch @ {
        ("H" -> {j_Integer} | j_Integer) :> {{cH, j - 1, 0}},
        ("S" -> {j_Integer} | j_Integer) :> {{cS, j - 1, 0}},
        ("X" -> {j_Integer} | j_Integer) :> {{cX, j - 1, 0}},
        ("Y" -> {j_Integer} | j_Integer) :> {{cY, j - 1, 0}},
        ("Z" -> {j_Integer} | j_Integer) :> {{cZ, j - 1, 0}},
        (SuperDagger["S"] -> {j_Integer} | j_Integer) :> {{cSdg, j - 1, 0}},
        (("CNOT" | "CX") -> {j_Integer, k_Integer}) :> {{cCNOT, j - 1, k - 1}},
        ("SWAP" -> {j_Integer, k_Integer}) :> {{cSWAP, j - 1, k - 1}},
        ("CZ" -> {j_Integer, k_Integer}) :> {{cH, k - 1, 0}, {cCNOT, j - 1, k - 1}, {cH, k - 1, 0}},
        "C"[("NOT" | "X") -> {k_Integer}, {j_Integer}, ___] :> {{cCNOT, j - 1, k - 1}},
        "C"["Z" -> {k_Integer}, {j_Integer}, ___] :> {{cH, k - 1, 0}, {cCNOT, j - 1, k - 1}, {cH, k - 1, 0}},
        {op_, qudits__Integer} :> With[{f = Flatten[{qudits}]},
            Replace[op -> If[Length[f] == 1, First[f], f], $stabilizerEncodeDispatch]
        ],
        _ :> $Failed
    }
]

encodeStabilizerGates[gateSpecs_List] := With[{rows = Replace[gateSpecs, $stabilizerEncodeDispatch, {1}]},
    If[FreeQ[rows, $Failed, {1}], Developer`ToPackedArray[Catenate[rows]], $Failed]
]


(* ============================================================================ *)
(* applyCompiledFold[ps, gates]: pack ps generator-major, run the compiled       *)
(* kernel once on the encoded gate matrix, unpack to a canonical PauliStabilizer. *)
(* Developer`ToPackedArray at every boundary keeps the compiled kernel on its     *)
(* fast path (an unpacked argument silently falls back to interpreted).           *)
(* ============================================================================ *)

applyCompiledFold[ps_PauliStabilizer, gates_ /; MatchQ[gates, {{_Integer, _Integer, _Integer} ..}]] :=
    Module[{n = ps["Qubits"], t = ps["Tableau"], phase, nw, res},
        nw = Ceiling[2 n / $bitsPerWord];
        phase = Developer`ToPackedArray[(1 - ps["Signs"]) / 2];
        res = compiledCliffordFold[
            packGenRows[t[[1]], nw],
            packGenRows[t[[2]], nw],
            First @ packGenRows[{phase}, nw],
            stabilizerValidMasks[n, nw],
            Developer`ToPackedArray[gates]
        ];
        PauliStabilizer[<|
            "Signs" -> 1 - 2 First @ unpackGenRows[{res[[-1]]}, 2 n],
            "Tableau" -> {unpackGenRows[res[[1 ;; n]], 2 n], unpackGenRows[res[[n + 1 ;; 2 n]], 2 n]}
        |>]
    ]

(* No gates: identity. *)
applyCompiledFold[ps_PauliStabilizer, {}] := ps


(* ============================================================================ *)
(* Public method: apply a whole Clifford circuit in one compiled call.          *)
(*                                                                              *)
(*   ps["ApplyCircuit", {"H" -> 1, "CNOT" -> {1, 2}, "S" -> 2, ...}]            *)
(*                                                                              *)
(* Each spec is a QuantumCircuitOperator-style gate (`op -> order`). When every  *)
(* gate is in the compiled set and ps is concrete, the whole list folds in one   *)
(* compiled kernel call; otherwise it falls back to the per-gate dispatch (which *)
(* also handles non-Clifford gates / StabilizerFrame).                           *)
(* ============================================================================ *)

ps_PauliStabilizer["ApplyCircuit", gateSpecs_List] :=
    With[{gates = If[psConcreteFastQ[ps], encodeStabilizerGates[gateSpecs], $Failed]},
        If[ ListQ[gates],
            applyCompiledFold[ps, gates],
            Fold[#1[#2] &, ps, gateSpecs]
        ]
    ]
