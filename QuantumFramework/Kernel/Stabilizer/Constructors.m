Package["Wolfram`QuantumFramework`"]

PackageScope[FromFullTableau]



(* ============================================================================ *)
(* Constructors: from QuantumState (4^n tomography).                            *)
(*                                                                              *)
(* Compute <psi|P|psi> for every Pauli P (4^n tomography), keep the 2^n Paulis *)
(* with |<P>| ~= 1 as the stabilizer group with their signs, pick n linearly-  *)
(* independent generators by greedy F_2 rank, then extend to a 2n-dim          *)
(* symplectic basis to obtain destabilizers that satisfy the AG commutation    *)
(* pattern. Cost is O(4^n) -- practical for n <= ~8.                            *)
(* ============================================================================ *)

(* Encode a tuple of Paulis (each in 0..3) as a 2n-bit row {x_1,..,x_n,z_1,..,z_n}. *)
(* Convention matches PauliRow / FromFullTableau: I=00, X=10, Y=11, Z=01.       *)
agEncodePauli[plist_List] := Catenate[Transpose @ Replace[plist, {0 -> {0, 0}, 1 -> {1, 0}, 2 -> {1, 1}, 3 -> {0, 1}}, {1}]]

(* Find n linearly-independent rows of `bits` over F_2 by greedy rank-growth.  *)
(* Returns indices into `bits`. Catch + Internal`Bag avoids the AppendTo/Break  *)
(* anti-pattern; the algorithm itself is genuinely sequential (each accepted   *)
(* index changes the rank residual seen by later candidates).                   *)
agPickIndependent[bits_, n_Integer] := Module[{bag = Internal`Bag[]},
    Catch @ Do[
        With[{kept = Internal`BagPart[bag, All]},
            If[MatrixRank[Append[bits[[kept]], bits[[i]]], Modulus -> 2] > Length[kept],
                Internal`StuffBag[bag, i];
                If[Length[kept] + 1 >= n, Throw[Null]]
            ]
        ],
        {i, Length[bits]}
    ];
    Internal`BagPart[bag, All]
]

(* Symplectic Gram matrix of a set of 2n-bit Pauli encodings over F_2. The      *)
(* (i, j) entry is 1 exactly when encodings i and j anticommute; the diagonal   *)
(* is always 0, so a mutually commuting set has an all-zero Gram matrix.        *)
agSymplecticGram[rows_, n_Integer] := agModTwoDot[Join[rows[[All, n + 1 ;;]], rows[[All, ;; n]], 2], Transpose[rows]]

(* Product of two 0/1 matrices mod 2. Entries of the product are at most the    *)
(* inner dimension, exact in machine reals, so the product runs as a real       *)
(* matrix product and is rounded back.                                          *)
agModTwoDot[a_, b_] := Mod[Round[N[a] . N[b]], 2]

(* Given n stabilizer encodings S (n x 2n) that form a square, independent,     *)
(* isotropic set, construct n destabilizer encodings D with                     *)
(* symp(D_i, S_j) = delta_ij and symp(D_i, D_k) = 0. One mod-2 solve with n     *)
(* right-hand sides gives D0 with S J D0^T = I. Its Gram matrix G = D0 J D0^T   *)
(* is symmetric with zero diagonal, and D = D0 + L S with L the strictly lower  *)
(* triangle of G gives D J D^T = G + L + L^T = 0 while S J D^T stays I, since   *)
(* S J S^T = 0. Cost: one O(n^3) solve. Any input outside the precondition       *)
(* returns $Failed (no message), so the solve only runs on a solvable system.   *)
agExtendToSymplecticBasis[stabBits_ ? MatrixQ, n_Integer] /;
    Dimensions[stabBits] === {n, 2 n} && MatrixRank[stabBits, Modulus -> 2] === n && Total[agSymplecticGram[stabBits, n], 2] === 0 :=
    With[{d0 = Transpose[LinearSolve[Join[stabBits[[All, n + 1 ;;]], stabBits[[All, ;; n]], 2], IdentityMatrix[n], Modulus -> 2]]},
        Mod[d0 + agModTwoDot[LowerTriangularize[agModTwoDot[Join[d0[[All, n + 1 ;;]], d0[[All, ;; n]], 2], Transpose[d0]], -1], stabBits], 2]
    ]

agExtendToSymplecticBasis[_, _] := $Failed

(* Genuine destabilizer half {Xdest, Zdest} for a full (square) stabilizer      *)
(* tableau {Xstab, Zstab}. Destabilizer D_i is the symplectic partner of        *)
(* stabilizer S_i: symp(D_i, S_j) = KroneckerDelta[i, j] and symp(D_i, D_k) = 0.*)
(* A per-qubit X<->Z swap of the stabilizers satisfies these relations only for *)
(* the all-Z register, so a general code (5-qubit, Steane, Shor) needs the       *)
(* symplectic Gram-Schmidt to give ps["Circuit"] the valid full tableau its      *)
(* Aaronson-Gottesman decomposition assumes. agExtendToSymplecticBasis returns   *)
(* the destabilizer encodings for a square/independent/isotropic set and         *)
(* $Failed otherwise, so a non-square, rank-deficient, or non-commuting          *)
(* generating set leaves the caller on the swap padding. Normal[] first, so a    *)
(* SparseArray-backed tableau (also PauliTableauQ) is handled like a dense one.  *)
agDestabilizerTableau[tableau_ ? PauliTableauQ] := With[{stab = Normal[tableau]},
    With[{
        stabBits = MapThread[Join, {Transpose[stab[[1]]], Transpose[stab[[2]]]}],
        nQubits = Length[stab[[1]]]
    },
        Replace[
            agExtendToSymplecticBasis[stabBits, nQubits],
            {
                destBits_ ? MatrixQ :> {Transpose[destBits[[All, ;; nQubits]]], Transpose[destBits[[All, nQubits + 1 ;;]]]},
                _ :> $Failed
            }
        ]
    ]
]

(* Compute <psi|P|psi> for every Pauli P on n qubits. Returns a list of        *)
(* {encoding, sign} pairs for the 2^n Paulis with |<P>| ~= 1. The first  *)
(* element is always {identity, 1}. *)
agStabilizerGroup[v_, n_Integer] := Enclose @ Module[{paulis, expVals, group, identity, nonid},
    paulis = Tuples[Range[0, 3], n];
    expVals = Map[Chop[Conjugate[v] . (kroneckerProduct @@ Map[pauliMatrix, #]) . v] &, paulis];
    group = Cases[
        Transpose[{paulis, expVals}],
        {p_, val_} /; Re[val] >= 1 - 10^-8 :> {agEncodePauli[p],  1},
        {1}
    ] ~Join~ Cases[
        Transpose[{paulis, expVals}],
        {p_, val_} /; Re[val] <= -(1 - 10^-8) :> {agEncodePauli[p], -1},
        {1}
    ];
    ConfirmAssert[Length[group] == 2^n,
        "PauliStabilizer[QuantumState] expects a stabilizer state; got " <>
        ToString[Length[group]] <> " group elements (need 2^" <> ToString[n] <> ")"
    ];
    group
]

PauliStabilizer[qs_QuantumState] := Enclose @ Module[
    {n, v, encodings, signs, nonidIdx, genIdx, stabBits, stabSigns,
     destBits, destSigns, fullTableau, basePS, baseAssoc, baseVec, anchor, phase},
    n = qs["Qudits"];
    v = Normal @ qs["Computational"]["StateVector"];

    {encodings, signs} = Transpose @ ConfirmMatch[
        agStabilizerGroup[v, n],
        {{_, _} ..}
    ];

    (* Drop the identity row (all-zero encoding) before picking generators. *)
    nonidIdx = Position[encodings, _ ? (! AllTrue[#, # == 0 &] &), {1}, Heads -> False][[All, 1]];
    genIdx = nonidIdx[[ agPickIndependent[encodings[[nonidIdx]], n] ]];
    stabBits = encodings[[genIdx]];
    stabSigns = signs[[genIdx]];

    (* Canonicalize generator order to match the integer-constructor          *)
    (* convention PauliStabilizer[n] = {Z_1, Z_2, ..., Z_n}.                  *)
    (* Sort by FirstPosition[row, 1] ascending (row whose active bit is at   *)
    (* the lowest-index qubit comes first).                                   *)
    With[{order = Ordering[
        FirstPosition[#, 1, {2 n + 1}] & /@ stabBits
    ]},
        stabBits = stabBits[[order]];
        stabSigns = stabSigns[[order]]
    ];

    destBits = agExtendToSymplecticBasis[stabBits, n];
    destSigns = ConstantArray[1, n];

    (* Build the (2n) x (2n+1) tableau in PauliRow/FromFullTableau format:     *)
    (* rows 1..n  = destabilizers, rows n+1..2n = stabilizers.                 *)
    (* Last column = phase (0 = sign +1, 1 = sign -1).                         *)
    fullTableau = Join[
        MapThread[Append, {destBits, (1 - destSigns) / 2}],
        MapThread[Append, {stabBits, (1 - stabSigns) / 2}]
    ];

    basePS = Confirm @ FromFullTableau[fullTableau];
    baseAssoc = First[basePS];

    (* Recover the overall phase the tableau drops: ps["State"] is exact only    *)
    (* up to a complex unit. Anchor at the first nonzero amplitude of the         *)
    (* tableau-projected state and store the ratio to the original.               *)
    baseVec = Normal @ basePS["State"]["StateVector"];
    anchor = First @ FirstPosition[baseVec, x_ /; Chop[x] =!= 0, {1}, Heads -> False];
    phase = v[[anchor]] / baseVec[[anchor]];

    PauliStabilizer[<|baseAssoc, "GlobalPhase" -> phase|>]
]


(* ============================================================================ *)
(* Association-keyed constructors                                               *)
(* ============================================================================ *)

PauliStabilizer[data : KeyValuePattern[{"Phase" -> phase_}]] := PauliStabilizer[<|KeyDrop[data, "Phase"], "Signs" -> 1 - 2 phase|>]

PauliStabilizer[data : KeyValuePattern[{"Matrix" -> mat_}]] := PauliStabilizer[<|KeyDrop[data, "Matrix"], "Tableau" -> Transpose[ArrayReshape[mat, {Length[mat], 2, Length[mat] / 2}], {3, 1, 2}]|>]

PauliStabilizer[data : KeyValuePattern[{"Tableau" -> t_}]] /; ! KeyExistsQ[data, "Signs"] :=
    PauliStabilizer[<|"Signs" -> ConstantArray[1, Dimensions[t][[3]]], "Tableau" -> t|>]


(* ============================================================================ *)
(* From QuantumOperator                                                         *)
(* ============================================================================ *)

FromFullTableau[t_] := PauliStabilizer[<|
    "Signs" -> 1 - 2 t[[All, -1]],
    "Tableau" -> Transpose[ArrayReshape[t[[All, ;; -2]], {Length[t], 2, Length[t] / 2}], {3, 1, 2}]
|>]


(* Capture the global phase the AG decomposition drops. The tableau is         *)
(* invariant under U -> e^{i alpha} U, so for Y = -i X Z the AG-recovered     *)
(* circuit Z@X equals i*Y; we record -i as "GlobalPhase" so that the inverse  *)
(* path ps["QuantumOperator"] returns Y exactly.                                *)
PauliStabilizer[qo_QuantumOperator, n : _Integer : 1] /; Sort[qo["OutputOrder"]] == Sort[qo["InputOrder"]] :=
Enclose @ Module[{max, basePS, baseAssoc, m1, m2, anchor, phase},
    max = Max[n, qo["InputOrder"]];
    basePS = FromFullTableau[
        Confirm @ QuantumOperatorTableau[qo["Sort"]["Computational"]]
    ]["PadRight", max]["Permute", Join[Sort[qo["InputOrder"]], Complement[Range[max], qo["InputOrder"]]]];

    baseAssoc = First[basePS];
    m1 = Flatten[Normal @ qo["Matrix"]];
    m2 = Flatten[Normal @ basePS["QuantumOperator"]["Matrix"]];
    anchor = First @ FirstPosition[m2, x_ /; Chop[x] =!= 0, {1}, Heads -> False];
    phase = m1[[anchor]] / m2[[anchor]];

    PauliStabilizer[<|baseAssoc, "GlobalPhase" -> phase|>]
]


(* ============================================================================ *)
(* From QuantumCircuitOperator                                                  *)
(* ============================================================================ *)

(* A circuit whose every gate is in the compiled Clifford set, on wires 1..n,   *)
(* folds on the n-qubit tableau in one compiled call: the cost is polynomial in *)
(* n and the result is defined up to a global phase, as for every gate applied *)
(* to a tableau.                                                                *)
PauliStabilizer[qco_QuantumCircuitOperator] /; Min[qco["FullInputOrder"], 1] >= 1 :=
    With[{gates = stabilizerExecutionPlan[qco]["Gates"]},
        applyCompiledFold[PauliStabilizer[qco["Max"]], gates] /; ListQ[gates]
    ]

(* Any other circuit folds over its normal operators. Each operator becomes a   *)
(* stabilizer on all qubits that records its global phase, which builds a dense *)
(* matrix on every qubit of the circuit.                                        *)
PauliStabilizer[qco_QuantumCircuitOperator] := Fold[
    Collect[#1 /. ps_PauliStabilizer :> #2[ps], _PauliStabilizer, Simplify] &,
    PauliStabilizer[qco["Max"]],
    Replace[PauliStabilizer[#, qco["Max"]], _ ? FailureQ :> #] & /@ qco["NormalOperators"]
]


(* ============================================================================ *)
(* From rank-3 +/-1/0 array                                                     *)
(* ============================================================================ *)

PauliStabilizer[basis_ /; ArrayQ[basis, 3, MatchQ[0 | 1 | -1]] && MatchQ[Dimensions[basis], {2, n_, m_}]] :=
    Enclose @ PauliStabilizer[
        Confirm @ Replace[Union @@ #, {{0, 1} | {1} -> 1, {-1, 0} | {-1} -> -1, _ -> $Failed}] & /@ Transpose[basis, {3, 2, 1}],
        Abs[basis]
    ]


(* ============================================================================ *)
(* Single-row / paired-sign constructors (auto-pad destabilizer half)           *)
(* ============================================================================ *)

PauliStabilizer[tableau_ ? PauliTableauQ] := PauliStabilizer[1, tableau]

PauliStabilizer[sign : -1 | 1, tableau_ ? PauliTableauQ] := PauliStabilizer[{sign}, tableau]

PauliStabilizer[signs : {(-1 | 1) ...}, tableau_ ? PauliTableauQ] := With[
    {padSigns = PadRight[signs, Dimensions[tableau][[3]], 1], destHalf = agDestabilizerTableau[tableau]},
    If[ destHalf === $Failed,
        (* Under-determined or degenerate generating set: pad the destabilizer   *)
        (* half by the per-qubit X<->Z swap. This is a valid symplectic partner   *)
        (* only for the all-Z register, but it is the only closed form available  *)
        (* when no full symplectic completion exists.                            *)
        PauliStabilizer[<|
            "Signs" -> Join[padSigns, padSigns],
            "Tableau" -> MapThread[Join[##, 2] &, {Reverse[tableau], tableau}]
        |>],
        (* Full stabilizer state: pair each generator with a genuine destabilizer *)
        (* (signs canonically +1) so ps["Circuit"] reproduces ps["State"] from    *)
        (* |0...0> for named codes as it already does for random Cliffords.        *)
        PauliStabilizer[<|
            "Signs" -> Join[ConstantArray[1, Length[padSigns]], padSigns],
            "Tableau" -> MapThread[Join[##, 2] &, {destHalf, tableau}]
        |>]
    ]
]


(* ============================================================================ *)
(* String-list constructors (Pauli strings with optional sign prefix)           *)
(* ============================================================================ *)

$PauliString = Repeated["-" | "+", {0, 1}] ~~ ("I" | "X" | "Y" | "Z") ...

(* Each Pauli string as a 2n-bit row {x_1..x_n, z_1..z_n}, signs dropped; a      *)
(* shorter string is padded with I on the right, as the tableau builder does.     *)
psStringBits[signedPauliStrings_List] := With[{bodies = StringDelete[signedPauliStrings, StartOfString ~~ ("-" | "+")]},
    With[{pairs = Replace[Characters[StringPadRight[bodies, Max[StringLength[bodies]], "I"]], {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {2}]},
        Join[pairs[[All, All, 1]], pairs[[All, All, 2]], 2]
    ]
]

(* A list of signed Pauli strings is a stabilizer state only when the strings    *)
(* commute and have F2 rank n: the rank, not the count, fixes the state. Fewer   *)
(* independent strings fix a 2^(n - rank)-dimensional space of states, which no  *)
(* tableau represents, so the constructor refuses them. More strings than n with *)
(* rank n are reduced to n independent ones, each dropped string checked against *)
(* the sign the kept ones give it.                                               *)
PauliStabilizer[signedPauliStrings : {__String}] /; AllTrue[signedPauliStrings, StringMatchQ[$PauliString]] :=
    With[{bits = psStringBits[signedPauliStrings]},
    With[{n = Length[First[bits]] / 2, m = Length[bits]},
    With[{bad = FirstPosition[Mod[bits[[All, ;; n]] . Transpose[bits[[All, n + 1 ;;]]] + bits[[All, n + 1 ;;]] . Transpose[bits[[All, ;; n]]], 2], 1, Missing["None"], {2}],
          rank = MatrixRank[bits, Modulus -> 2]},
        Which[
            ! MissingQ[bad],
                Message[PauliStabilizer::noncomm, signedPauliStrings[[bad[[1]]]], signedPauliStrings[[bad[[2]]]]]; $Failed,
            rank < n,
                Message[PauliStabilizer::rank, m, n, rank, 2 ^ (n - rank)]; $Failed,
            m == n,
                psFromSignedStrings[signedPauliStrings],
            True,
                psFromIndependentStrings[signedPauliStrings, agPickIndependent[bits, n], n]
        ]
    ]]]

(* A signed string as "-" or "" followed by its letters padded with I to n.       *)
psNormalizeSignedString[s_String, n_Integer] :=
    If[StringStartsQ[s, "-"], "-", ""] <> StringPadRight[StringDelete[s, StartOfString ~~ ("-" | "+")], n, "I"]

(* The kept strings are independent, so every dropped string is +-1 times their   *)
(* product; it is consistent exactly when the kept tableau gives it expectation 1. *)
psFromIndependentStrings[signedPauliStrings_, keep_, n_] := With[{
    ps = psFromSignedStrings[signedPauliStrings[[keep]]],
    dropped = psNormalizeSignedString[#, n] & /@ Delete[signedPauliStrings, List /@ keep]
},
    With[{wrong = SelectFirst[dropped, stabilizerExpectation[ps, #] === -1 &]},
        If[ MissingQ[wrong],
            ps,
            Message[PauliStabilizer::inconsistent, wrong]; $Failed
        ]
    ]
]

psFromSignedStrings[signedPauliStrings_] :=
    PauliStabilizer @@
        MapAt[Transpose[PadRight[#, Automatic, {0, 0}], {3, 2, 1}] &, 2] @ Thread @ Replace[
            Characters[signedPauliStrings], {sign : "-" | "+" : 1, paulis___} :> {
                Replace[sign, {"-" -> -1, "+" -> 1}],
                Replace[{paulis}, {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}]
            },
            {1}
        ]

PauliStabilizer[stabString : {__String}, destabStrings : {__String}] /; AllTrue[Join[stabString, destabStrings], StringMatchQ[$PauliString]] :=
    PauliStabilizer[<|"Signs" -> #1, "Tableau" -> #2|>] & @@
        MapAt[Transpose[PadRight[#, Automatic, {0, 0}], {3, 2, 1}] &, 2] @ Thread @ Replace[
            Characters[Join[destabStrings, stabString]], {sign : "-" | "+" : 1, paulis___} :> {
                Replace[sign, {"-" -> -1, "+" -> 1}],
                Replace[{paulis}, {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}]
            },
            {1}
        ]


(* ============================================================================ *)
(* Integer constructor: |0...0> register in generator-major packed storage. *)
(* The public "Tableau" property still reconstructs the canonical packed-dense *)
(* integer arrays. Clifford evolution does not need to allocate those arrays. *)
(* ============================================================================ *)

PauliStabilizer[q_Integer ? Positive] := psZeroGeneratorPacked[q]

PauliStabilizer[q_Integer ? NonNegative] := PauliStabilizer[{ConstantArray[0, {Max[q, 1], q}], identityMatrix[q]}]


(* ============================================================================ *)
(* Empty + shortcut                                                             *)
(* ============================================================================ *)

PauliStabilizer[] := PauliStabilizer[1]

(* Registered PauliStabilizer names (NamedCodes.m, RandomClifford.m) have      *)
(* their own dispatch rules; the shortcut must not capture them, in either     *)
(* bare-string or "Name"[args] call form.                                      *)
pauliStabilizerShortcutQ[(name_String)[___]] := ! MemberQ[$PauliStabilizerNames, name]
pauliStabilizerShortcutQ[name_String] := ! MemberQ[$PauliStabilizerNames, name]
pauliStabilizerShortcutQ[_] := True

(* A bare rule ("CNOT" -> {1, 2}) is wrapped in a list because the circuit     *)
(* constructor only accepts shorthand rules inside a list. ConfirmBy turns an  *)
(* unrecognized name into a propagated Failure instead of an unevaluated       *)
(* PauliStabilizer[Failure[...]].                                              *)
PauliStabilizer[shortcut : (_String | _String[___] | (_String -> _ ? orderQ | _Integer) | _List) ? pauliStabilizerShortcutQ] :=
    Enclose @ PauliStabilizer @ ConfirmBy[QuantumCircuitOperator[Replace[shortcut, r_Rule :> {r}]], QuantumCircuitOperatorQ]
