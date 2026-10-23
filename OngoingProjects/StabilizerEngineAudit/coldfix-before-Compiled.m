Package["Wolfram`QuantumFramework`"]

PackageScope[$stabilizerGateCodes]
PackageScope[$stabilizerCliffordNames]
PackageScope[$stabilizerTwoQubitNames]
PackageScope[$stabilizerCompileTarget]
PackageScope[$bitsPerWord]
PackageScope[compiledCliffordFold]
PackageScope[encodeStabilizerGates]
PackageScope[applyCompiledFold]
PackageScope[packGenRows]
PackageScope[unpackGenRows]
PackageScope[stabilizerValidMasks]
PackageScope[psGeneratorPackedQ]
PackageScope[psGeneratorPackedDataQ]
PackageScope[psGetGeneratorPacked]
PackageScope[psToGeneratorPacked]
PackageScope[psZeroGeneratorPacked]
PackageScope[$stabilizerEncodeMemo]



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

(* Gate-name catalogs for the multi-gate dispatch forms ps[{specs..}] / ps[spec, ..]    *)
(* (Stabilizer/PauliStabilizer.m). These name what the top-level dispatch RECOGNIZES as  *)
(* a gate token; they are deliberately distinct from the two sibling gate sets:          *)
(*   - $stabilizerGateCodes (above): the gates the compiled C kernel encodes. "V" is    *)
(*     expanded to H/S/H by the encoder; "T" is non-Clifford and routes to a frame.       *)
(*     The integer kernel itself only needs the primitive gate set.                                                  *)
(*   - sfGate1 / sfGate2 (Stabilizer/StabilizerFrame.m): the gates whose Pauli-          *)
(*     conjugation tables the frame uses, keyed by normalized names ("Sdg", not          *)
(*     SuperDagger["S"]).                                                                  *)
(* A new gate must be added to whichever of the three its dispatch / engine / frame path *)
(* actually supports; they are not interchangeable.                                       *)
$stabilizerTwoQubitNames = {"CNOT", "CX", "CZ", "SWAP"};
$stabilizerCliffordNames = Join[{"H", "S", "X", "Y", "Z", "V", "T"}, $stabilizerTwoQubitNames];

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
(* Normal first: Developer`ToPackedArray is a NO-OP on a SparseArray (it stays    *)
(* sparse), and a sparse word matrix makes the compiled kernel reject its         *)
(* argument (CompiledFunction::cfta) and run interpreted; Normal on an already-   *)
(* dense list returns it unchanged, so the dense path pays nothing. The word      *)
(* values reach 2^62, so the Dot result comes back UNPACKED; force it back to a   *)
(* packed Int64 array, otherwise the compiled kernel silently falls back to       *)
(* interpreted evaluation on an unpacked argument.                                *)
packGenRows[bm_, nw_] := Developer`ToPackedArray @ Dot[
    ArrayReshape[
        PadRight[Developer`ToPackedArray[Normal @ bm], {Length[bm], nw $bitsPerWord}],
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
        ("V" -> {j_Integer} | j_Integer) :> {{cH, j - 1, 0}, {cS, j - 1, 0}, {cH, j - 1, 0}},
        (SuperDagger["V"] -> {j_Integer} | j_Integer) :> {{cH, j - 1, 0}, {cSdg, j - 1, 0}, {cH, j - 1, 0}},
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

$stabilizerEncodeMemo = <||>;

(* A hash selects the bucket; equality with the retained specification list
   checks the key, so a hash collision cannot reuse another circuit's gates. *)
encodeStabilizerGates[gateSpecs_List] := With[{key = Hash[gateSpecs]},
    With[{cached = Lookup[$stabilizerEncodeMemo, key, None]},
        If[ListQ[cached] && First[cached] === gateSpecs,
            Last[cached],
            With[{rows = Replace[gateSpecs, $stabilizerEncodeDispatch, {1}]},
                With[{gates = If[FreeQ[rows, $Failed, {1}], Developer`ToPackedArray[Catenate[rows]], $Failed]},
                    If[ListQ[gates],
                        If[Length[$stabilizerEncodeMemo] >= 16, $stabilizerEncodeMemo = <||>];
                        $stabilizerEncodeMemo[key] = {gateSpecs, gates}];
                    gates
                ]
            ]
        ]
    ]
]

(* Bulk results retain the generator-major layout. Signs remain a small ordinary
   vector for compatibility; the quadratic X/Z bit arrays stay packed. *)
psGeneratorPackedQ[ps_PauliStabilizer] := KeyExistsQ[First[ps], "GeneratorX"]
psGeneratorPackedQ[_] := False

psGeneratorPackedDataQ[a_Association] := With[{n = Lookup[a, "Qubits", 0],
    x = Lookup[a, "GeneratorX", None], z = Lookup[a, "GeneratorZ", None], s = Lookup[a, "Signs", None]},
    IntegerQ[n] && n > 0 && ListQ[s] && Length[s] === 2 n &&
    MatrixQ[x, IntegerQ] && MatrixQ[z, IntegerQ] &&
    Dimensions[x] === {n, Ceiling[2 n / $bitsPerWord]} && Dimensions[z] === Dimensions[x] &&
    Min[x, z] >= 0 && Max[x, z] < 2 ^ $bitsPerWord &&
    Max[x[[All, -1]], z[[All, -1]]] < 2 ^ (Mod[2 n - 1, $bitsPerWord] + 1)
]

ps_PauliStabilizer["Tableau"] /; psGeneratorPackedQ[ps] := With[{a = First[ps]},
    {unpackGenRows[a["GeneratorX"], 2 a["Qubits"]], unpackGenRows[a["GeneratorZ"], 2 a["Qubits"]]}
]

psGetGeneratorPacked[ps_PauliStabilizer] := With[{n = ps["Qubits"], a = First[ps]},
    If[psGeneratorPackedQ[ps],
        Developer`ToPackedArray /@ {a["GeneratorX"], a["GeneratorZ"]},
        With[{t = ps["Tableau"], nw = Ceiling[2 n / $bitsPerWord]},
            {packGenRows[t[[1]], nw], packGenRows[t[[2]], nw]}]
    ]
]


(* ============================================================================ *)
(* applyCompiledFold[ps, gates]: pack ps generator-major, run the compiled       *)
(* kernel once on the encoded gate matrix, retain its generator-packed result. *)
(* Developer`ToPackedArray at every boundary keeps the compiled kernel on its     *)
(* fast path (an unpacked argument silently falls back to interpreted).           *)
(* ============================================================================ *)

applyCompiledFold[ps_PauliStabilizer, gates_ /; MatchQ[gates, {{_Integer, _Integer, _Integer} ..}]] :=
    Module[{n = ps["Qubits"], phase, nw, xz, res, codes, controls, targets},
        (* RuntimeOptions -> "Speed" needs a checked boundary before native indexing. *)
        codes = gates[[All, 1]]; controls = gates[[All, 2]]; targets = gates[[All, 3]];
        If[Min[controls, targets] < 0 || Max[controls, targets] >= n ||
            MemberQ[Pick[controls - targets, Unitize[(codes - $stabilizerGateCodes["CNOT"])
                (codes - $stabilizerGateCodes["SWAP"])] , 0], 0],
            Message[PauliStabilizer::invalidorder, n]; Return[$Failed]];
        nw = Ceiling[2 n / $bitsPerWord];
        phase = Developer`ToPackedArray[(1 - ps["Signs"]) / 2];
        xz = psGetGeneratorPacked[ps];
        res = compiledCliffordFold[
            xz[[1]], xz[[2]], First @ packGenRows[{phase}, nw],
            stabilizerValidMasks[n, nw], Developer`ToPackedArray[gates]
        ];
        System`Private`HoldSetValid[PauliStabilizer[<|
            "Signs" -> 1 - 2 First @ unpackGenRows[{res[[-1]]}, 2 n],
            "GeneratorX" -> res[[1 ;; n]], "GeneratorZ" -> res[[n + 1 ;; 2 n]], "Qubits" -> n
        |>]]
    ]

(* Preserve storage without applying a gate, for operations such as inverse. *)
psToGeneratorPacked[ps_PauliStabilizer] := With[{xz = psGetGeneratorPacked[ps]},
    System`Private`HoldSetValid[PauliStabilizer[<|"Signs" -> ps["Signs"],
        "GeneratorX" -> xz[[1]], "GeneratorZ" -> xz[[2]], "Qubits" -> ps["Qubits"]|>]]
]

(* Identity tableau directly in generator-major words: no expanded n x 2n array. *)
psZeroGeneratorPacked[n_Integer ? Positive] := With[{
    rows = Range[n], dest = Range[0, n - 1], stab = Range[n, 2 n - 1],
    nw = Ceiling[2 n / $bitsPerWord]
},
    System`Private`HoldSetValid[PauliStabilizer[<|
        "Signs" -> ConstantArray[1, 2 n],
        "GeneratorX" -> Developer`ToPackedArray[Normal[SparseArray[
            Thread[Transpose[{rows, Quotient[dest, $bitsPerWord] + 1}] -> BitShiftLeft[1, Mod[dest, $bitsPerWord]]], {n, nw}]]],
        "GeneratorZ" -> Developer`ToPackedArray[Normal[SparseArray[
            Thread[Transpose[{rows, Quotient[stab, $bitsPerWord] + 1}] -> BitShiftLeft[1, Mod[stab, $bitsPerWord]]], {n, nw}]]],
        "Qubits" -> n
    |>]]
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
(* also handles non-Clifford gates / StabilizerFrame and returns a measurement's *)
(* outcomes), and a spec the engine cannot apply is refused with                *)
(* PauliStabilizer::nonclifford.                                                *)
(* ============================================================================ *)

ps_PauliStabilizer["ApplyCircuit", gateSpecs_List] :=
    With[{gates = If[psConcreteFastQ[ps], encodeStabilizerGates[gateSpecs], $Failed]},
        If[ ListQ[gates],
            applyCompiledFold[ps, gates],
            stabilizerGateFold[ps, gateSpecs, True]
        ]
    ]
