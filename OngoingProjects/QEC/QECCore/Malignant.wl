(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[exRecMalignantCounts]
PackageScope[exRecAtomEffects]
PackageScope[blockTables]


(* ============================================================================ *)
(* Malignant sets of locations (Got26 sec. 14.3.3, 14.5.1).                     *)
(*                                                                              *)
(* A set S of locations of an exRec is MALIGNANT if some choice of Pauli errors *)
(* at those locations, in a run the heralds accept, makes the exRec incorrect;  *)
(* BENIGN otherwise (Def 14.7).  Theorem 14.1 says every set of at most t       *)
(* locations is benign; the bad-exRec bound assumes every set of t + 1 is       *)
(* malignant, and most are not.  Counting the ones that are replaces            *)
(* Binomial[A, t + 1] by M_(t+1) in the bound (eq. 14.31-14.32).  Pauli errors  *)
(* suffice (end of sec. 14.3.3), and the identity is allowed at a location --   *)
(* but a set with an identity in it is a smaller set, already counted.          *)
(*                                                                              *)
(* CORRECT, computed.  The rectangle is correct when decoding after it gives    *)
(* the ideal gate applied to the decoding of what entered it (Def 14.5).  In    *)
(* the Pauli frame: E_mid is the residual on the data at the seam, the leading  *)
(* EC's correction included; E_ref is E_mid pushed through the gate gadget     *)
(* alone; E_out is the residual at the end, both corrections included.  The      *)
(* exRec is correct iff the ideal decoder puts E_out and E_ref in the same      *)
(* logical class.  For a preparation only the logical X part counts (a logical  *)
(* Z does nothing to |0_L>); for a measurement, the decoded readout is compared *)
(* with the ideal decoding of E_mid.                                            *)
(*                                                                              *)
(* WHY THE INPUT CAN BE A PERFECT CODEWORD, though Def 14.7 quantifies over     *)
(* faults outside the exRec, i.e. over the input.  The leading EC measures      *)
(* s = syn(E_in) + f and leaves E_in + F + c(s), whose syndrome is              *)
(* syn(F) + f: free of E_in.  Two residuals with one syndrome differ by a       *)
(* stabilizer times a logical; the stabilizer does nothing and the logical      *)
(* passes through the transversal gadget and both decoders alike, so both      *)
(* sides of the comparison shift together.  The verdict depends on the input    *)
(* only through syn(E_in), and not at all.  (This is the ECRP at work.)         *)
(*                                                                              *)
(* HOW IT IS FAST ENOUGH.  Everything but the decoders is linear over GF(2).    *)
(* Each single-qubit X or Z at each place is propagated once, and its effect is *)
(* kept as a handful of integers -- herald parities, leading and trailing       *)
(* syndromes, the seam frame, the final frame, the readout.  A fault is the XOR *)
(* of its parts, a pair the XOR of two faults, and the decoders are table        *)
(* lookups, so all pairs against one location are one vectorised pass.         *)
(* ============================================================================ *)


(* ---- bits and tables ---- *)

bitsToInt[bits_List] := FromDigits[Reverse[bits], 2]
intParity[x_Integer] := Mod[DigitCount[x, 2, 1], 2]

(* The per-block tables, over all 2^n frames: the syndrome a frame gives, the
   lightest correction of a syndrome, and the logical bit of the decoded frame. *)
blockTables[a_Association] := blockTables[a] = Module[
    {n = a["Qubits"], hz, hx, zbar, xbar, corrX, corrZ, synZ, synX},
    hz = verificationChecks[a, "Z"][[1 ;; -1 - codeLogicalQubits[a]]];
    hx = verificationChecks[a, "X"][[1 ;; -1 - codeLogicalQubits[a]]];
    zbar = First[codeLogicalVectors[a]["Z"]][[n + 1 ;; 2 n]];
    xbar = First[codeLogicalVectors[a]["X"]][[1 ;; n]];
    synZ = Table[bitsToInt[Mod[hz . IntegerDigits[x, 2, n][[-1 ;; 1 ;; -1]], 2]], {x, 0, 2^n - 1}];
    synX = Table[bitsToInt[Mod[hx . IntegerDigits[z, 2, n][[-1 ;; 1 ;; -1]], 2]], {z, 0, 2^n - 1}];
    corrX = Table[bitsToInt[cssCorrection[a, "X", IntegerDigits[s, 2, Length[hz]][[-1 ;; 1 ;; -1]]]], {s, 0, 2^Length[hz] - 1}];
    corrZ = Table[bitsToInt[cssCorrection[a, "Z", IntegerDigits[s, 2, Length[hx]][[-1 ;; 1 ;; -1]]]], {s, 0, 2^Length[hx] - 1}];
    <|
        "n" -> n, "mz" -> Length[hz], "mx" -> Length[hx],
        "SynZ" -> synZ, "SynX" -> synX, "CorrX" -> corrX, "CorrZ" -> corrZ,
        (* the logical X-type bit (anticommutes with Zbar) of the decoded X frame,
           and the logical Z-type bit of the decoded Z frame *)
        "ClassX" -> Table[intParity[BitAnd[BitXor[x, corrX[[synZ[[x + 1]] + 1]]], bitsToInt[zbar]]], {x, 0, 2^n - 1}],
        "ClassZ" -> Table[intParity[BitAnd[BitXor[z, corrZ[[synX[[z + 1]] + 1]]], bitsToInt[xbar]]], {z, 0, 2^n - 1}]
    |>
]


(* ---- the effect of one single-qubit Pauli ---- *)

(* For a list of elementary faults {index, qubit, {x, z}}, the effect of each as an
   association of integers.  Blocks are packed side by side, n bits (or m bits of
   syndrome) per block. *)
exRecAtomEffects[d_Association, elementary_List] := Module[
    {a = d["Code"], instr = d["Instructions"], nq = d["Qubits"], n = d["Code"]["Qubits"], B = d["Blocks"],
     seam = Length[Lookup[d["Parts"], "LeadingEC", {}]], tb, rec, lead, trail, meas, order = d["Order"], hc = d["HeraldChecks"],
     prefix, packFrame, packSyn},
    tb = blockTables[a];
    lead = Count[Lookup[d["Parts"], "LeadingEC", {}], {"M", _}];
    meas = Count[Lookup[d["Parts"], "Measurement", {}], {"M", _}];
    prefix = Take[instr, seam];
    packFrame[fr_] := Total[Table[bitsToInt[fr[[(b - 1) n + 1 ;; b n]]] 2^((b - 1) n), {b, B}]];
    (* the 2n records of each EC on each block -> the packed X and Z syndromes *)
    packSyn[r_List] := If[ r === {}, {0, 0},
        With[{perBlock = Partition[r, 2 n]},
            With[{halves = If[order === "BitFlipFirst", {#[[1 ;; n]], #[[n + 1 ;; 2 n]]}, {#[[n + 1 ;; 2 n]], #[[1 ;; n]]}] & /@ perBlock},
                {Total[Table[tb["SynZ"][[bitsToInt[halves[[b, 1]]] + 1]] 2^((b - 1) tb["mz"]), {b, Length[halves]}]],
                 Total[Table[tb["SynX"][[bitsToInt[halves[[b, 2]]] + 1]] 2^((b - 1) tb["mx"]), {b, Length[halves]}]]}
            ]
        ]
    ];
    Table[
        Module[{run = framePropagate[instr, nq, {e}], r, x, z, s1, s2, w},
            r = run["Record"];
            {x, z} = run["Frame"];
            s1 = packSyn[r[[1 ;; lead]]];
            s2 = packSyn[r[[lead + meas + 1 ;;]]];
            w = If[ e[[1]] <= seam && lead > 0 && ! TrueQ[e[[4]]],
                With[{pr = framePropagate[prefix, nq, {Most[e]}]["Frame"]}, {packFrame[pr[[1]]], packFrame[pr[[2]]]}],
                {0, 0}
            ];
            <|
                "H" -> bitsToInt[Catenate[Mod[#[[2]] . run["Heralds"][[#[[1]]]], 2] & /@ hc]],
                "S1x" -> s1[[1]], "S1z" -> s1[[2]], "S2x" -> s2[[1]], "S2z" -> s2[[2]],
                "Wx" -> w[[1]], "Wz" -> w[[2]],
                "Fx" -> packFrame[x], "Fz" -> packFrame[z],
                "R" -> bitsToInt[r[[lead + 1 ;; lead + meas]]]
            |>
        ],
        {e, elementary}
    ]
]


(* ---- counting ---- *)

(* The fields of a fault, as integers, from its elementary parts.  An elementary
   part is {index, qubit, {x, z}} and is split into its X and its Z part, each of
   which was propagated once. *)
$atomFields = {"H", "S1x", "S1z", "S2x", "S2z", "Wx", "Wz", "Fx", "Fz", "R"};

exRecMalignantCounts[rect_QECExtendedRectangle] := exRecMalignantCounts[First[rect]]

exRecMalignantCounts[d_Association] := exRecMalignantCounts[d] = Module[
    {a = d["Code"], n = d["Code"]["Qubits"], B = d["Blocks"], tb, locs, seam, layerOf, lastLead, elem, elemIndex,
     effects, faultVec, faults, corrXP, corrZP, classXP, classZP, synZP, synXP, c1x, c1z, gadget, refTables,
     verdict, singles, pairs, kind = d["Location"], size, nloc},

    tb = blockTables[a];
    locs = d["LocationList"];
    nloc = Length[locs];
    seam = Length[Lookup[d["Parts"], "LeadingEC", {}]];
    size = 2^(n B);

    (* A data-qubit wait after that qubit's last leading-EC instruction belongs to
       the rectangle, not to the leading EC, even though it is anchored there. *)
    lastLead = Association @ Table[
        q -> Max[0, Max @ Flatten @ Position[Take[d["Instructions"], seam], {_, ___, q, ___}, {1}, Heads -> False]],
        {q, n B}];

    (* every elementary X and Z, with a flag for parts that do not touch the seam *)
    elem = DeleteDuplicates @ Flatten[Table[
            With[{rectQ = loc["Type"] === "Wait" && f[[2]] <= n B && f[[1]] >= Lookup[lastLead, f[[2]], 0] && TrueQ[loc["Gadget"] || f[[1]] === Lookup[lastLead, f[[2]], -1]]},
                DeleteCases[{
                    If[f[[3, 1]] === 1, {f[[1]], f[[2]], {1, 0}, rectQ}, Nothing],
                    If[f[[3, 2]] === 1, {f[[1]], f[[2]], {0, 1}, rectQ}, Nothing]}, Nothing]
            ],
            {loc, locs}, {choice, loc["Faults"]}, {f, choice}
        ], 3];
    elemIndex = AssociationThread[elem -> Range[Length[elem]]];
    effects = exRecAtomEffects[d, elem];
    effects = Transpose[Lookup[effects, $atomFields]];  (* field -> list over elem *)

    (* each fault choice as a vector of field integers *)
    faultVec[loc_, choice_] := With[{ids = Catenate @ Table[
            With[{rectQ = loc["Type"] === "Wait" && f[[2]] <= n B && f[[1]] >= Lookup[lastLead, f[[2]], 0] && TrueQ[loc["Gadget"] || f[[1]] === Lookup[lastLead, f[[2]], -1]]},
                DeleteCases[{
                    If[f[[3, 1]] === 1, elemIndex[{f[[1]], f[[2]], {1, 0}, rectQ}], Nothing],
                    If[f[[3, 2]] === 1, elemIndex[{f[[1]], f[[2]], {0, 1}, rectQ}], Nothing]}, Nothing]],
            {f, choice}]},
        BitXor @@@ (effects[[All, ids]] /. {} -> {0})
    ];
    faults = Table[faultVec[loc, ch], {loc, locs}, {ch, loc["Faults"]}];   (* loc -> choice -> 10 ints *)

    (* packed tables over the blocks *)
    corrXP = Table[Total[Table[tb["CorrX"][[BitAnd[BitShiftRight[s, (b - 1) tb["mz"]], 2^tb["mz"] - 1] + 1]] 2^((b - 1) n), {b, B}]], {s, 0, 2^(tb["mz"] B) - 1}];
    corrZP = Table[Total[Table[tb["CorrZ"][[BitAnd[BitShiftRight[s, (b - 1) tb["mx"]], 2^tb["mx"] - 1] + 1]] 2^((b - 1) n), {b, B}]], {s, 0, 2^(tb["mx"] B) - 1}];
    classXP = Table[Total[Table[tb["ClassX"][[BitAnd[BitShiftRight[x, (b - 1) n], 2^n - 1] + 1]] 2^(b - 1), {b, B}]], {x, 0, size - 1}];
    classZP = Table[Total[Table[tb["ClassZ"][[BitAnd[BitShiftRight[z, (b - 1) n], 2^n - 1] + 1]] 2^(b - 1), {b, B}]], {z, 0, size - 1}];

    (* the leading correction's own effect downstream: a data Pauli placed at the
       seam, propagated through the rest -- linear, so tabulated over all values *)
    {c1x, c1z} = If[ seam === 0, {ConstantArray[{0, 0, 0, 0}, size], ConstantArray[{0, 0, 0, 0}, size]},
        Module[{basis = exRecAtomEffects[d, Catenate @ Table[{seam, q, p, True}, {p, {{1, 0}, {0, 1}}}, {q, n B}]], bx, bz, span},
            bx = Lookup[basis[[1 ;; n B]], {"S2x", "S2z", "Fx", "Fz"}];
            bz = Lookup[basis[[n B + 1 ;;]], {"S2x", "S2z", "Fx", "Fz"}];
            span[bv_] := Table[With[{on = Flatten @ Position[Reverse[IntegerDigits[v, 2, n B]], 1]}, If[on === {}, {0, 0, 0, 0}, BitXor @@ bv[[on]]]], {v, 0, size - 1}];
            {span[bx], span[bz]}
        ]
    ];

    (* the reference: E_mid through the gate gadget alone *)
    gadget = Lookup[d["Parts"], "Gate", {}];
    refTables = If[ gadget === {}, None,
        Module[{gq = d["Qubits"], basis, bx, bz, span},
            basis = Table[With[{fr = framePropagate[gadget, gq, {{0, q, p}}]["Frame"]},
                    {Total[Table[bitsToInt[fr[[1, (b - 1) n + 1 ;; b n]]] 2^((b - 1) n), {b, B}]],
                     Total[Table[bitsToInt[fr[[2, (b - 1) n + 1 ;; b n]]] 2^((b - 1) n), {b, B}]]}],
                {p, {{1, 0}, {0, 1}}}, {q, n B}];
            span[bv_] := Table[With[{on = Flatten @ Position[Reverse[IntegerDigits[v, 2, n B]], 1]}, If[on === {}, {0, 0}, BitXor @@ bv[[on]]]], {v, 0, size - 1}];
            {span[basis[[1]]], span[basis[[2]]]}
        ]
    ];

    (* the verdict on a batch of combined faults, given as field arrays: True where
       the run is accepted and the exRec incorrect *)
    verdict[v_List] := Module[{h, s1x, s1z, s2x, s2z, wx, wz, fx, fz, r, cx1, cz1, ex, ez, dx, dz, out, ref, mx, mz},
        {h, s1x, s1z, s2x, s2z, wx, wz, fx, fz, r} = v;
        cx1 = corrXP[[s1x + 1]]; cz1 = corrZP[[s1z + 1]];
        (* downstream of the leading correction *)
        dx = c1x[[cx1 + 1]]; dz = c1z[[cz1 + 1]];
        s2x = BitXor[s2x, dx[[All, 1]], dz[[All, 1]]];
        s2z = BitXor[s2z, dx[[All, 2]], dz[[All, 2]]];
        fx = BitXor[fx, dx[[All, 3]], dz[[All, 3]], corrXP[[s2x + 1]]];
        fz = BitXor[fz, dx[[All, 4]], dz[[All, 4]], corrZP[[s2z + 1]]];
        mx = BitXor[wx, cx1]; mz = BitXor[wz, cz1];
        Switch[kind,
            "Preparation",
                out = classXP[[fx + 1]]; ref = 0 out,
            "Measurement",
                out = classXP[[BitXor[r, cx1] + 1]]; ref = classXP[[mx + 1]],
            _,
                If[ refTables === None,
                    {ex, ez} = {mx, mz},
                    ex = BitXor[refTables[[1, mx + 1, 1]], refTables[[2, mz + 1, 1]]];
                    ez = BitXor[refTables[[1, mx + 1, 2]], refTables[[2, mz + 1, 2]]]
                ];
                out = classXP[[fx + 1]] + 2^B classZP[[fz + 1]];
                ref = classXP[[ex + 1]] + 2^B classZP[[ez + 1]]
        ];
        (1 - Unitize[h]) Unitize[BitXor[out, ref]]
    ];

    (* single locations: Theorem 14.1 says none is malignant, and this checks it *)
    singles = Table[
        Max[verdict[Transpose[faults[[i]]]]],
        {i, nloc}
    ];

    (* pairs: every choice at i against every choice at every j > i, in one pass per i *)
    pairs = Catenate @ Table[
        Module[{js = Range[i + 1, nloc], lens, flat, combo, flags, hit},
            If[js === {}, {},
                lens = Length /@ faults[[js]];
                flat = Catenate[faults[[js]]];   (* list of 10-int vectors *)
                hit = Total @ Table[
                    combo = BitXor[ConstantArray[fi, Length[flat]], flat];
                    flags = verdict[Transpose[combo]];
                    Max /@ TakeList[flags, lens],
                    {fi, faults[[i]]}
                ];
                {i, #} & /@ Pick[js, Unitize[hit], 1]
            ]
        ],
        {i, nloc}
    ];

    <|"Singles" -> Total[singles], "Pairs" -> Length[pairs], "Locations" -> nloc,
      "SingleList" -> Flatten @ Position[singles, 1], "PairList" -> pairs,
      "PairsByType" -> KeySort @ Counts[Sort[{locs[[#[[1]]]]["Type"], locs[[#[[2]]]]["Type"]}] & /@ pairs]|>
]
