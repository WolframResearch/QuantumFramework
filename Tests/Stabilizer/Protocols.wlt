(* ::Package:: *)

(* ============================================================================
   Tests/Stabilizer/Protocols.wlt -- textbook protocols run entirely through
   the two bracket entrances of the tableau:

       qco_QuantumCircuitOperator[ps]        a circuit applied to a tableau
       qmo_QuantumMeasurementOperator[ps]    a Pauli measurement operator on
                                             a few wires of a wider tableau

   Every protocol below is a chain of these two calls with feed-forward in
   between: state preparation by a circuit, a partial Pauli measurement on a
   few wires of a larger register, an outcome-dependent Pauli correction by
   another circuit, and a final partial readout. The expected values are
   derived from theory and computed in the expected slot, never read off the
   kernel: the Heisenberg conjugation rules of H, S, X, CNOT and CZ, Bell-state
   algebra, the parity-check matrix of the Steane code, the Mermin argument,
   and the cluster-state stabilizers K_k = Z_{k-1} X_k Z_{k+1}. Registers run
   up to 50 qubits, beyond any dense reference. Wherever a protocol ends in a
   definite state, that state is written down from its generators as a
   reference tableau and pinned by the closed-form overlap (Method ->
   "ClosedForm", O(n^3), no state vector): 1 for the same state up to a
   phase, 0 for an orthogonal one.

   Outcome convention: a measurement of a Pauli P returns <|bit -> tableau|>
   with bit b <-> eigenvalue 1 - 2 b, so a +1 eigenstate reads b = 0 and a -1
   eigenstate reads b = 1. A deterministic outcome is a single key; a uniformly
   random one is the key set {0, 1}. Every tableau that comes out of a chain
   is also checked to keep the Aaronson-Gottesman pairing, M Omega M^T = Omega
   over F_2, the invariant a measurement update can silently break.
   ========================================================================== *)

Needs["Wolfram`QuantumFramework`"];


(* ---------------------------------------------------------------------------- *)
(* Fixtures and classical bookkeeping                                            *)
(* ---------------------------------------------------------------------------- *)

(* The six single-qubit stabilizer states: the Clifford preparation on wire 1    *)
(* from |0>, the Pauli letter that stabilizes the state, and the bit that letter *)
(* reads on it (0 on the +1 eigenstate, 1 on the -1 eigenstate).                 *)
$stabilizerInputs = {
    {"Z+", {},                             "Z", 0},
    {"Z-", {"X" -> 1},                     "Z", 1},
    {"X+", {"H" -> 1},                     "X", 0},
    {"X-", {"X" -> 1, "H" -> 1},           "X", 1},
    {"Y+", {"H" -> 1, "S" -> 1},           "Y", 0},
    {"Y-", {"X" -> 1, "H" -> 1, "S" -> 1}, "Y", 1}
};
stabilizerInput[name_String] := SelectFirst[$stabilizerInputs, First[#] === name &]

(* Sorted outcome keys of a Pauli measurement operator on the named wires. *)
outcomeBits[ps_, label_String, wires_List] := Sort @ Keys @ QuantumMeasurementOperator[label, wires][ps]

(* Signed Pauli strings for reference tableaux: an n-letter string from a       *)
(* wire -> letter association (identity elsewhere), and a sign prefix.           *)
letterString[n_Integer, letters_Association] := StringJoin @ Lookup[letters, Range[n], "I"]
signPrefix[s_] := If[s === -1, "-", ""]

(* Overlap magnitude |<psi|phi>| of two tableaux by the closed form. *)
closedFormOverlap[ps_, other_] := ps["InnerProduct", other, Method -> "ClosedForm"]

(* The Aaronson-Gottesman pairing: the 2n generator rows form a symplectic       *)
(* matrix over F_2, M Omega M^T = Omega with Omega = [[0, I], [I, 0]].           *)
tableauSymplecticQ[ps_] := With[{n = ps["Qubits"]},
    With[{omega = ArrayFlatten[{{ConstantArray[0, {n, n}], IdentityMatrix[n]}, {IdentityMatrix[n], ConstantArray[0, {n, n}]}}]},
        Mod[ps["Matrix"] . omega . Transpose[ps["Matrix"]], 2] === omega
    ]
]

(* Classical Pauli-frame bookkeeping: the Heisenberg action of H and of X on a   *)
(* signed single-qubit Pauli {sign, letter}. H swaps X and Z and negates Y; a   *)
(* byproduct X fixes X and negates Y and Z. A measured wire with outcome m      *)
(* pushes the frame through H and then through X^m.                             *)
frameH[{s_, "X"}] := {s, "Z"}
frameH[{s_, "Z"}] := {s, "X"}
frameH[{s_, "Y"}] := {-s, "Y"}
frameX[{s_, "X"}] := {s, "X"}
frameX[{s_, p_}] := {-s, p}
framePush[frame_, outcomes_List] := Fold[If[#2 === 1, frameX[frameH[#1]], frameH[#1]] &, frame, outcomes]

(* Sign to outcome bit, and letter to Bloch axis. *)
signBit[s_] := (1 - s) / 2
blochAxis[letter_] := UnitVector[3, Replace[letter, {"X" -> 1, "Y" -> 2, "Z" -> 3}]]


(* ============================================================================ *)
(* TELEPORTATION -- |psi> on wire 1, Bell pair on wires 2,3 (H2, CNOT23), Bell   *)
(* measurement of wires 1,2 as CNOT12, H1, then Z on wire 1 (bit a) and Z on    *)
(* wire 2 (bit b). Wire 3 then holds X^b Z^a |psi>; the correction X^b Z^a       *)
(* restores |psi>. Before the correction the byproduct flips the readout of     *)
(* exactly the letter it anticommutes with: Z^a flips X and Y, X^b flips Z and  *)
(* Y. All four (a, b) appear for every input: the Bell measurement reveals      *)
(* nothing about |psi>. The corrected register is |a>|b>|psi>, written down as  *)
(* a reference tableau.                                                          *)
(* ============================================================================ *)

byproductFlip["X", a_, b_] := a
byproductFlip["Z", a_, b_] := b
byproductFlip["Y", a_, b_] := BitXor[a, b]

teleportationReference[a_, b_, letter_, sign_] :=
    PauliStabilizer[{signPrefix[1 - 2 a] <> "ZII", signPrefix[1 - 2 b] <> "IZI", signPrefix[sign] <> "II" <> letter}]

(* One run per Bell-measurement branch: {a, b, bit read on wire 3 before the    *)
(* correction, bit read after it, overlap with |a>|b>|psi>, pairing invariant}. *)
teleportationBranches[{name_, prep_, letter_, bit_}] := With[{
    m1 = QuantumMeasurementOperator["Z", {1}][
        QuantumCircuitOperator[Join[prep, {"H" -> 2, "CNOT" -> {2, 3}, "CNOT" -> {1, 2}, "H" -> 1}]][PauliStabilizer[3]]
    ]
},
    {name, Flatten[Table[
        With[{m2 = QuantumMeasurementOperator["Z", {2}][m1[a]]},
            Table[
                With[{corrected = QuantumCircuitOperator[{If[b === 1, "X" -> 3, Nothing], If[a === 1, "Z" -> 3, Nothing]}][m2[b]]},
                    {a, b, outcomeBits[m2[b], letter, {3}], outcomeBits[corrected, letter, {3}],
                        closedFormOverlap[corrected, teleportationReference[a, b, letter, 1 - 2 bit]], tableauSymplecticQ[corrected]}
                ],
                {b, Sort @ Keys @ m2}
            ]
        ],
        {a, Sort @ Keys @ m1}
    ], 1]}
]

VerificationTest[
    teleportationBranches /@ $stabilizerInputs,
    Replace[$stabilizerInputs, {name_, _, letter_, bit_} :>
        {name, Flatten[Table[{a, b, {BitXor[bit, byproductFlip[letter, a, b]]}, {bit}, 1, True}, {a, 0, 1}, {b, 0, 1}], 1]},
        {1}
    ],
    {},
    TestID -> "Protocol-Teleportation-SixInputs-AllBranches-ByproductRule-Correction"
]

(* The byproduct rule as one object: measuring wires 1 and 2 symbolically       *)
(* leaves the outcomes as free signs s1, s2 on Z1 and Z2 and one sign-free      *)
(* generator that correlates the wire-3 Pauli with them: a Z on wire 1 iff the  *)
(* byproduct Z^a flips the letter, a Z on wire 2 iff X^b does, so Z1 Z2 Y3 for  *)
(* a Y input, Z1 X3 for X, Z2 Z3 for Z, with the input's sign. Substituting     *)
(* the outcomes recovers each branch as the uncorrected register                *)
(* |a>|b> X^b Z^a |psi>.                                                         *)
teleportationSymbolic[{name_, prep_, letter_, bit_}] := With[{
    sym = QuantumCircuitOperator[Join[prep, {"H" -> 2, "CNOT" -> {2, 3}, "CNOT" -> {1, 2}, "H" -> 1}]][PauliStabilizer[3]]["SymbolicMeasure", {1, 2}]
},
    With[{
        outcomeSymbol = Function[letters,
            First @ Cases[
                Cases[Transpose[{sym["Stabilizers"], sym["StabilizerSigns"]}], {s_String, sign_} /; StringEndsQ[s, letters] :> sign],
                _\[FormalS], Infinity
            ]
        ]
    },
        {name,
            Select[sym["Stabilizers"], StringMatchQ[#, RegularExpression["^-?[IXYZ]+$"]] &],
            Flatten[Table[
                {a, b, closedFormOverlap[
                    sym["SubstituteOutcomes", {outcomeSymbol["ZII"] -> a, outcomeSymbol["IZI"] -> b}],
                    teleportationReference[a, b, letter, (1 - 2 bit) (-1)^byproductFlip[letter, a, b]]
                ]},
                {a, 0, 1}, {b, 0, 1}
            ], 1]}
    ]
]

VerificationTest[
    teleportationSymbolic /@ $stabilizerInputs,
    Replace[$stabilizerInputs, {name_, _, letter_, bit_} :>
        {name,
            {signPrefix[1 - 2 bit] <> If[byproductFlip[letter, 1, 0] === 1, "Z", "I"] <> If[byproductFlip[letter, 0, 1] === 1, "Z", "I"] <> letter},
            Flatten[Table[{a, b, 1}, {a, 0, 1}, {b, 0, 1}], 1]},
        {1}
    ],
    {},
    TestID -> "Protocol-Teleportation-SymbolicOutcomes-ByproductRuleIsOneGenerator"
]


(* ============================================================================ *)
(* ENTANGLEMENT SWAPPING -- Bell pairs on (1,2) and (3,4); Bell measurement of   *)
(* wires 2,3 as CNOT23, H2, Z2 (bit a), Z3 (bit b). Projecting wires 2,3 onto   *)
(* the Bell state with XX = (-1)^a and ZZ = (-1)^b leaves wires 1,4 in that     *)
(* same Bell state, read here as two-letter measurements on the non-adjacent    *)
(* wires {1, 4} and pinned as the reference |a>|b> on 2,3 beside that Bell     *)
(* state on 1,4. One ebit crosses the cut {1,2}|{3,4} only after the swap, and  *)
(* each measured wire ends in a definite Z eigenstate, entangled with nothing.  *)
(* ============================================================================ *)

swappingReference[a_, b_] := PauliStabilizer[{
    signPrefix[1 - 2 a] <> "XIIX", signPrefix[1 - 2 b] <> "ZIIZ",
    signPrefix[1 - 2 a] <> "IZII", signPrefix[1 - 2 b] <> "IIZI"
}]

VerificationTest[
    With[{pairs = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 3, "CNOT" -> {3, 4}}][PauliStabilizer[4]]},
        With[{m2 = QuantumMeasurementOperator["Z", {2}][QuantumCircuitOperator[{"CNOT" -> {2, 3}, "H" -> 2}][pairs]]},
            {
                {pairs["Entropy", {1}], pairs["Entropy", {1, 2}]},
                Flatten[Table[
                    With[{m3 = QuantumMeasurementOperator["Z", {3}][m2[a]]},
                        Table[
                            {a, b, outcomeBits[m3[b], "XX", {1, 4}], outcomeBits[m3[b], "ZZ", {1, 4}],
                                {m3[b]["Entropy", {1, 2}], m3[b]["Entropy", {2}], m3[b]["Entropy", {3}]},
                                closedFormOverlap[m3[b], swappingReference[a, b]], tableauSymplecticQ[m3[b]]},
                            {b, Sort @ Keys @ m3}
                        ]
                    ],
                    {a, Sort @ Keys @ m2}
                ], 1]
            }
        ]
    ],
    {
        {1, 0},
        Flatten[Table[{a, b, {a}, {b}, {1, 0, 0}, 1, True}, {a, 0, 1}, {b, 0, 1}], 1]
    },
    {},
    TestID -> "Protocol-EntanglementSwapping-OutcomesLabelTheSwappedPair"
]


(* ============================================================================ *)
(* MERMIN-GHZ -- local X and Y measurements on GHZ-n, one wire at a time. A     *)
(* string of X and Y letters with an even number k of Y's equals i^k times     *)
(* X^{(x)n} times a Z string of even weight, both stabilizers of GHZ-n, so its  *)
(* eigenvalue is (-1)^(k/2): the product of the n local outcomes is fixed in    *)
(* every branch, while the first n-1 outcomes are free (2^(n-1) branches) and   *)
(* the last is determined. No assignment of +-1 to the 2n local observables    *)
(* satisfies all these products at once (brute force over 4^n assignments).    *)
(* ============================================================================ *)

(* Every branch of n sequential single-wire measurements, as <|"Bits", "Tableau"|>. *)
sequentialBranches[ps_, letters_String] := Fold[
    Function[{branches, measurement}, Flatten[Table[
        With[{r = measurement[branch["Tableau"]]},
            Table[<|"Bits" -> Append[branch["Bits"], key], "Tableau" -> r[key]|>, {key, Sort @ Keys @ r}]
        ],
        {branch, branches}
    ], 1]],
    {<|"Bits" -> {}, "Tableau" -> ps|>},
    MapIndexed[QuantumMeasurementOperator, Characters[letters]]
]

(* The X/Y strings of length n with an even number of Y letters. *)
evenYStrings[n_Integer] := Select[StringJoin /@ Tuples[{"X", "Y"}, n], EvenQ[StringCount[#, "Y"]] &]

(* Number of classical +-1 assignments to (x_1, .., x_n, y_1, .., y_n) that     *)
(* reproduce the quantum product (-1)^(k/2) for every even-Y string.            *)
classicalAssignmentCount[n_Integer] := Count[Tuples[{1, -1}, 2 n],
    v_ /; AllTrue[evenYStrings[n],
        Function[s, Product[If[StringTake[s, {i}] === "X", v[[i]], v[[n + i]]], {i, n}] === (-1)^(StringCount[s, "Y"] / 2)]
    ]
]

VerificationTest[
    Table[
        With[{ghz = QuantumCircuitOperator["GHZ"[n]][PauliStabilizer[n]]},
            {n,
                Table[
                    With[{branches = sequentialBranches[ghz, s]},
                        {s, DeleteDuplicates[Times @@ (1 - 2 #["Bits"]) & /@ branches], Length[branches],
                            AllTrue[branches, tableauSymplecticQ[#["Tableau"]] &]}
                    ],
                    {s, evenYStrings[n]}
                ],
                classicalAssignmentCount[n]}
        ],
        {n, {3, 4, 5}}
    ],
    Table[{n, Table[{s, {(-1)^(StringCount[s, "Y"] / 2)}, 2^(n - 1), True}, {s, evenYStrings[n]}], 0}, {n, {3, 4, 5}}],
    {},
    TestID -> "Protocol-MerminGHZ-LocalOutcomeProducts-NoClassicalAssignment-n3-n4-n5"
]


(* ============================================================================ *)
(* ONE-BIT TELEPORTATION -- the measurement-based primitive: |psi> on wire 1,   *)
(* |+> on wire 2, CZ(1,2), then X on wire 1 with outcome m leaves wire 2 in    *)
(* X^m H |psi>. For a stabilizer input the output Pauli is the input Pauli      *)
(* pushed through H and then through X^m, read as a deterministic outcome on   *)
(* wire 2 alone. Both m appear for every input.                                 *)
(* ============================================================================ *)

oneBitTeleportation[{name_, prep_, letter_, bit_}] := With[{
    m = QuantumMeasurementOperator["X", {1}][QuantumCircuitOperator[Join[prep, {"H" -> 2, "CZ" -> {1, 2}}]][PauliStabilizer[2]]]
},
    {name, Table[
        Replace[framePush[{1 - 2 bit, letter}, {k}], {_, pauli_} :> {k, pauli, outcomeBits[m[k], pauli, {2}], tableauSymplecticQ[m[k]]}],
        {k, Sort @ Keys @ m}
    ]}
]

VerificationTest[
    oneBitTeleportation /@ $stabilizerInputs,
    Replace[$stabilizerInputs, {name_, _, letter_, bit_} :>
        {name, Table[Replace[framePush[{1 - 2 bit, letter}, {k}], {sign_, pauli_} :> {k, pauli, {signBit[sign]}, True}], {k, 0, 1}]},
        {1}
    ],
    {},
    TestID -> "Protocol-OneBitTeleportation-XMeasurementOnClusterPair-PauliFrame"
]


(* ============================================================================ *)
(* MEASUREMENT-BASED WIRE -- a linear cluster of n wires carries |psi> from      *)
(* wire 1 to wire n under X measurements of wires 1 .. n-1 with outcomes m_k:   *)
(* the output is the ordered product of (X^{m_k} H) applied to |psi>, so its    *)
(* Pauli is the input Pauli pushed through the frame rule once per measured     *)
(* wire. In closed form: a Y input keeps its letter and picks up the sign       *)
(* (-1)^(n - 1 + Sum m_k); an X input ends as X or Z by the parity of n-1 with  *)
(* the sign (-1)^(Sum of m_k over odd k), a Z input the same with even k. The   *)
(* whole register ends as |x_{m_1}> ... |x_{m_{n-1}}> beside the output Pauli,  *)
(* written down as a reference tableau, and the Bloch vector of the last wire  *)
(* is read exactly at lengths where the dense state would have 2^n amplitudes. *)
(* ============================================================================ *)

outcomeClosedForm[{s_, "Y"}, outcomes_List] := {s (-1)^(Length[outcomes] + Total[outcomes]), "Y"}
outcomeClosedForm[{s_, "X"}, outcomes_List] := {s (-1)^Total[outcomes[[1 ;; ;; 2]]], If[EvenQ[Length[outcomes]], "X", "Z"]}
outcomeClosedForm[{s_, "Z"}, outcomes_List] := {s (-1)^Total[outcomes[[2 ;; ;; 2]]], If[EvenQ[Length[outcomes]], "Z", "X"]}

clusterWireReference[n_, outcomes_List, {sign_, pauli_}] := PauliStabilizer[Join[
    Table[signPrefix[1 - 2 outcomes[[k]]] <> letterString[n, <|k -> "X"|>], {k, n - 1}],
    {signPrefix[sign] <> letterString[n, <|n -> pauli|>]}
]]

(* {name, outcome pattern, bit read on the last wire for the predicted letter,   *)
(* Bloch vector of the last wire, overlap with the reference, pairing invariant}. *)
clusterWire[{name_, prep_, letter_, bit_}, outcomes_List] := With[{
    n = Length[outcomes] + 1,
    predicted = framePush[{1 - 2 bit, letter}, outcomes]
},
    With[{
        cluster = QuantumCircuitOperator[Join[prep, Table["H" -> k, {k, 2, n}], Table["CZ" -> {k, k + 1}, {k, n - 1}]]][PauliStabilizer[n]],
        xMeasurements = Table[QuantumMeasurementOperator["X", {k}], {k, n - 1}]
    },
        With[{
            output = Fold[
                Function[{state, wireOutcome}, Replace[wireOutcome, {k_, m_} :> xMeasurements[[k]][state][m]]],
                cluster,
                Transpose[{Range[n - 1], outcomes}]
            ]
        },
            Replace[predicted, {sign_, pauli_} :>
                {name, outcomes, outcomeBits[output, pauli, {n}],
                    Table[output["Expectation", StringRepeat["I", n - 1] <> p], {p, {"X", "Y", "Z"}}],
                    closedFormOverlap[output, clusterWireReference[n, outcomes, predicted]],
                    tableauSymplecticQ[output]}
            ]
        ]
    ]
]

clusterWirePrediction[{name_, _, letter_, bit_}, outcomes_List] :=
    Replace[framePush[{1 - 2 bit, letter}, outcomes], {sign_, pauli_} :> {name, outcomes, {signBit[sign]}, sign blochAxis[pauli], 1, True}]

(* Every outcome pattern on a 5-wire cluster for every input, against the       *)
(* closed form in the outcome string.                                           *)
VerificationTest[
    Table[Rest @ clusterWire[input, outcomes], {input, $stabilizerInputs}, {outcomes, Tuples[{0, 1}, 4]}],
    Table[
        Replace[outcomeClosedForm[{1 - 2 input[[4]], input[[3]]}, outcomes], {sign_, pauli_} :> {outcomes, {signBit[sign]}, sign blochAxis[pauli], 1, True}],
        {input, $stabilizerInputs}, {outcomes, Tuples[{0, 1}, 4]}
    ],
    {},
    TestID -> "Protocol-ClusterWire5-AllOutcomePatterns-ClosedFormInOutcomes"
]

(* Closed form in the length: all outcomes 0 down chains of n = 2 .. 13 wires,  *)
(* where the wire applies H^(n-1).                                              *)
VerificationTest[
    Table[Rest @ clusterWire[stabilizerInput[name], ConstantArray[0, n - 1]], {name, {"Y+", "X+"}}, {n, 2, 13}],
    {
        Table[{ConstantArray[0, n - 1], {signBit[(-1)^(n - 1)]}, (-1)^(n - 1) blochAxis["Y"], 1, True}, {n, 2, 13}],
        Table[{ConstantArray[0, n - 1], {0}, blochAxis[If[EvenQ[n - 1], "X", "Z"]], 1, True}, {n, 2, 13}]
    },
    {},
    TestID -> "Protocol-ClusterWire-AllZeroOutcomes-ClosedFormInLength"
]

(* All six inputs down a 12-wire cluster along their own random outcome         *)
(* patterns, and the Y+ input along the all-zero and all-one patterns at both   *)
(* parities of the length (12 and 13 wires).                                    *)
$wirePatterns = BlockRandom[SeedRandom[20260921];
    Join[
        Table[{input, RandomInteger[1, 11]}, {input, $stabilizerInputs}],
        Catenate @ Table[{stabilizerInput["Y+"], ConstantArray[m, len - 1]}, {m, 0, 1}, {len, {12, 13}}]
    ]
];

VerificationTest[
    clusterWire @@@ $wirePatterns,
    clusterWirePrediction @@@ $wirePatterns,
    {},
    TestID -> "Protocol-ClusterWire12-13-XChain-OutputPauliFollowsFrame"
]

(* A 41-wire cluster carrying |+i> along one random outcome pattern. *)
VerificationTest[
    With[{outcomes = BlockRandom[SeedRandom[20260922]; RandomInteger[1, 40]]},
        clusterWire[stabilizerInput["Y+"], outcomes] === clusterWirePrediction[stabilizerInput["Y+"], outcomes]
    ],
    True,
    {},
    TestID -> "Protocol-ClusterWire41-XChain-BeyondDense"
]

(* The closed form for every outcome string at once: an X measurement is a Z    *)
(* measurement after H on the measured wire, so one symbolic measurement of     *)
(* wires 1 .. n-1 leaves n-1 free outcome signs s_k on Z_k and gives the wire-n *)
(* Y the expectation (-1)^(n-1) Product_k (1 - 2 s_k), the closed form in all   *)
(* 2^(n-1) strings; substituting one string reproduces outcomeClosedForm.       *)
clusterWireSymbolic[n_Integer, outcomes_List] := With[{
    sym = QuantumCircuitOperator[Table["H" -> k, {k, n - 1}]][
        QuantumCircuitOperator[Join[{"H" -> 1, "S" -> 1}, Table["H" -> k, {k, 2, n}], Table["CZ" -> {k, k + 1}, {k, n - 1}]]][PauliStabilizer[n]]
    ]["SymbolicMeasure", Range[n - 1]]
},
    With[{
        symbols = DeleteDuplicates @ Cases[sym["StabilizerSigns"], _\[FormalS], Infinity],
        wireSymbols = Association @ Cases[Transpose[{sym["Stabilizers"], sym["StabilizerSigns"]}],
            {s_String, sign_} /; StringMatchQ[StringTake[s, -n], RegularExpression["^I*ZI*$"]] :>
                (First @ First @ StringPosition[StringTake[s, -n], "Z"] -> First @ Cases[sign, _\[FormalS], Infinity])]
    },
        {n, Length[symbols],
            Expand[sym["Expectation", StringRepeat["I", n - 1] <> "Y"] - (-1)^(n - 1) Times @@ (1 - 2 symbols)],
            sym["SubstituteOutcomes", Table[wireSymbols[k] -> outcomes[[k]], {k, n - 1}]]["Expectation", StringRepeat["I", n - 1] <> "Y"]}
    ]
]

$symbolicWirePatterns = BlockRandom[SeedRandom[20260923]; Table[{n, RandomInteger[1, n - 1]}, {n, {5, 12, 13, 41}}]];

VerificationTest[
    clusterWireSymbolic @@@ $symbolicWirePatterns,
    Replace[$symbolicWirePatterns, {n_, outcomes_} :> {n, n - 1, 0, First @ outcomeClosedForm[{1, "Y"}, outcomes]}, {1}],
    {},
    TestID -> "Protocol-ClusterWire-SymbolicOutcomes-ClosedFormForAllStrings"
]


(* ============================================================================ *)
(* CLUSTER STABILIZERS ON UNSORTED WIRES -- K_k = Z_{k-1} X_k Z_{k+1} measured  *)
(* as the label "XZZ" on the wires {k, k-1, k+1}: the letters follow the wires  *)
(* named, so this is K_k and reads 0 on every interior wire of a 30-wire chain. *)
(* The same letters on the sorted wires {k-1, k, k+1} name X_{k-1} Z_k Z_{k+1}, *)
(* which anticommutes with K_{k+1}: a uniformly random outcome on every         *)
(* interior wire, both branches keeping the pairing. The wire order is physics, *)
(* not bookkeeping. The same chain carries the graph-state area law: the       *)
(* entanglement entropy of a cut A is the F_2 rank of the adjacency block       *)
(* between A and its complement (Hein, Eisert, Briegel 2004), which for         *)
(* contiguous blocks counts the block boundaries (1, 2, 3 ebits here) and for   *)
(* the alternating cut reaches 15.                                              *)
(* ============================================================================ *)

$clusterCuts = {Range[10], Range[8, 17], Join[Range[5], Range[15, 20]], {1, 3, 5}, Range[1, 29, 2], {2, 3, 4, 10, 20, 21}, Range[30]};

pathCutRank[n_Integer, cut_List] := With[{rest = Complement[Range[n], cut]},
    If[rest === {}, 0, MatrixRank[Normal[AdjacencyMatrix[PathGraph[Range[n]]]][[cut, rest]], Modulus -> 2]]
]

VerificationTest[
    With[{n = 30},
        With[{cluster = QuantumCircuitOperator[Join[Table["H" -> k, {k, n}], Table["CZ" -> {k, k + 1}, {k, n - 1}]]][PauliStabilizer[n]]},
            {
                tableauSymplecticQ[cluster],
                DeleteDuplicates @ Table[outcomeBits[cluster, "XZZ", {k, k - 1, k + 1}], {k, 2, n - 1}],
                DeleteDuplicates @ Table[outcomeBits[cluster, "XZZ", {k - 1, k, k + 1}], {k, 2, n - 1}],
                AllTrue[Flatten @ Table[Values @ QuantumMeasurementOperator["XZZ", {k - 1, k, k + 1}][cluster], {k, 2, n - 1}], tableauSymplecticQ],
                Table[cluster["Entropy", cut], {cut, $clusterCuts}]
            }
        ]
    ],
    {True, {{0}}, {{0, 1}}, True, Table[pathCutRank[30, cut], {cut, $clusterCuts}]},
    {},
    TestID -> "Protocol-Cluster30-StabilizerOnUnsortedWires-AreaLaw"
]


(* ============================================================================ *)
(* STEANE CODE -- the [[7,1,3]] CSS code on the Hamming supports {1,3,5,7},     *)
(* {2,3,6,7}, {4,5,6,7}: column j of the parity-check matrix is the binary      *)
(* expansion of j. Encoding by the synthesized circuit through the bracket      *)
(* form; syndrome extraction by six four-letter measurements, each on its own   *)
(* four wires of the seven, carrying the post-measurement tableau forward;      *)
(* single-qubit Pauli errors injected as one-gate circuits on each of the       *)
(* logical states |0_L>, |1_L> and |+_L> (transversal H on |0_L>); correction   *)
(* by the decoded Pauli. The syndrome never depends on the logical state.       *)
(* ============================================================================ *)

$hammingSupports = {{1, 3, 5, 7}, {2, 3, 6, 7}, {4, 5, 6, 7}};

(* The three Z-type checks, then the three X-type checks. *)
$steaneChecks = Join[
    QuantumMeasurementOperator["ZZZZ", #] & /@ $hammingSupports,
    QuantumMeasurementOperator["XXXX", #] & /@ $hammingSupports
];

(* Sequential syndrome extraction: <|"Syndrome" -> six sorted key lists,        *)
(* "Tableau" -> the post-measurement tableau|>.                                 *)
steaneExtract[ps_] := Fold[
    With[{r = #2[#1["Tableau"]]}, <|"Syndrome" -> Append[#1["Syndrome"], Sort @ Keys @ r], "Tableau" -> First @ Values @ r|>] &,
    <|"Syndrome" -> {}, "Tableau" -> ps|>,
    $steaneChecks
]

(* Column j of the parity-check matrix: which supports hold wire j. *)
hammingColumn[j_] := Boole[MemberQ[#, j]] & /@ $hammingSupports

(* The syndrome of a Pauli error: an X component trips the Z-type checks whose   *)
(* support holds the wire, a Z component trips the X-type checks, Y trips both. *)
errorSyndrome[p_String, j_Integer] :=
    List /@ Join[If[p === "Z", {0, 0, 0}, hammingColumn[j]], If[p === "X", {0, 0, 0}, hammingColumn[j]]]

(* The decoder: the 21 single-qubit syndromes are distinct, so each names its     *)
(* error; recovery applies the named Pauli to the post-extraction tableau.       *)
$steaneDecoder = Association @ Flatten @ Table[errorSyndrome[p, j] -> {p -> j}, {p, {"X", "Y", "Z"}}, {j, 7}];
steaneRecover[extracted_Association] := QuantumCircuitOperator[$steaneDecoder[extracted["Syndrome"]]][extracted["Tableau"]]

$steaneLogicalStates = With[{code = PauliStabilizer["SteaneCode"]},
    <|"0L" -> code, "1L" -> PauliStabilizer["SteaneCode1"], "+L" -> QuantumCircuitOperator[Table["H" -> k, {k, 7}]][code]|>
];

VerificationTest[
    With[{code = $steaneLogicalStates["0L"]},
        With[{encoded = code["Circuit"][PauliStabilizer[7]]},
            {closedFormOverlap[encoded, code], steaneExtract[encoded]["Syndrome"], outcomeBits[encoded, "ZZZZZZZ", Range[7]], tableauSymplecticQ[encoded]}
        ]
    ],
    {1, Table[{0}, {6}], {0}, True},
    {},
    TestID -> "Protocol-Steane-EncodingCircuit-PreparesLogicalZero"
]

(* Three distinct code states with clean syndromes: |0_L> and |1_L> orthogonal,  *)
(* |+_L> at overlap 1/Sqrt[2] with each, the logical X reading 0 on |+_L> and    *)
(* the logical Z reading 1 on |1_L>.                                             *)
VerificationTest[
    With[{states = $steaneLogicalStates},
        {
            closedFormOverlap[states["0L"], states["1L"]],
            closedFormOverlap[states["+L"], states["0L"]], closedFormOverlap[states["+L"], states["1L"]],
            outcomeBits[states["+L"], "XXXXXXX", Range[7]], outcomeBits[states["1L"], "ZZZZZZZ", Range[7]],
            Table[steaneExtract[states[name]]["Syndrome"], {name, {"0L", "1L", "+L"}}]
        }
    ],
    {0, 1/Sqrt[2], 1/Sqrt[2], {0}, {1}, Table[Table[{0}, {6}], {3}]},
    {},
    TestID -> "Protocol-Steane-LogicalStates-CleanSyndromes"
]

(* On each logical state and for every single-qubit Pauli error the syndrome is *)
(* deterministic and equals the parity-check prediction, the decoder names the  *)
(* injected error from the syndrome alone, the pairing survives the             *)
(* six-measurement chain, and the recovery returns the same logical state.      *)
VerificationTest[
    KeyValueMap[
        Function[{name, state},
            {name, Flatten[Table[
                With[{extracted = steaneExtract[QuantumCircuitOperator[{p -> j}][state]]},
                    {p, j, extracted["Syndrome"], $steaneDecoder[extracted["Syndrome"]], tableauSymplecticQ[extracted["Tableau"]],
                        closedFormOverlap[steaneRecover[extracted], state]}
                ],
                {p, {"X", "Y", "Z"}}, {j, 7}
            ], 1]}
        ],
        $steaneLogicalStates
    ],
    Table[{name, Flatten[Table[{p, j, errorSyndrome[p, j], {p -> j}, True, 1}, {p, {"X", "Y", "Z"}}, {j, 7}], 1]}, {name, {"0L", "1L", "+L"}}],
    {},
    TestID -> "Protocol-Steane-SingleQubitErrors-SyndromeIsParityCheckColumn-OnEveryLogicalState"
]

(* Distance 3: a weight-2 error X_i X_j has the syndrome of the single error     *)
(* X_{i xor j} (the columns are binary numbers and add over F_2), so the decoder *)
(* names that wire, clears the syndrome, and leaves X_i X_j X_{i xor j}, a       *)
(* weight-3 logical X, on the code: the logical Z reads 1 and the state is       *)
(* |1_L>, a flip no later syndrome can see. Every one of the 21 pairs does this. *)
VerificationTest[
    With[{code = $steaneLogicalStates["0L"], pairs = Subsets[Range[7], {2}]},
        {
            AllTrue[pairs, Mod[Total[hammingColumn /@ #], 2] === hammingColumn[BitXor @@ #] &],
            Table[
                With[{hidden = steaneExtract[QuantumCircuitOperator[Thread["X" -> pair]][code]]},
                    With[{decoded = steaneRecover[hidden]},
                        {pair, hidden["Syndrome"], $steaneDecoder[hidden["Syndrome"]], steaneExtract[decoded]["Syndrome"],
                            outcomeBits[decoded, "ZZZZZZZ", Range[7]],
                            closedFormOverlap[decoded, $steaneLogicalStates["1L"]], closedFormOverlap[decoded, code]}
                    ]
                ],
                {pair, pairs}
            ]
        }
    ],
    {
        True,
        Table[{pair, errorSyndrome["X", BitXor @@ pair], {"X" -> BitXor @@ pair}, Table[{0}, {6}], {1}, 1, 0}, {pair, Subsets[Range[7], {2}]}]
    },
    {},
    TestID -> "Protocol-Steane-WeightTwoErrors-AllDecodeToLogicalFlip"
]


(* ============================================================================ *)
(* GHZ-n -- at n = 9 and at n = 50, beyond any dense reference (2^50            *)
(* amplitudes), measuring an end wire, an interior wire and the other end wire. *)
(* A Z measurement of wire w collapses every wire to the same bit: the state is *)
(* the product state |m m ... m>, the X-parity coherence is gone and no cut     *)
(* carries entropy. An X measurement of wire w detaches it in the X eigenstate  *)
(* (-1)^m and leaves the other n-1 wires in the GHZ state whose X-parity sign   *)
(* is (-1)^m: every adjacent ZZ among them still reads 0, and a cut A carries   *)
(* one ebit exactly when A minus {w} is a nonempty proper subset of the other   *)
(* wires. Both post-measurement states are pinned by the overlap with the       *)
(* tableau written down from these generators, and distinguished from the       *)
(* tableau with the opposite X-parity sign.                                     *)
(* ============================================================================ *)

(* Reference tableaux: the product state |m ... m>, and wire w in the X          *)
(* eigenstate (-1)^m beside the GHZ state of the other wires with X-parity       *)
(* sign paritySign.                                                              *)
ghzProductReference[n_, m_] := PauliStabilizer[Table[signPrefix[1 - 2 m] <> letterString[n, <|k -> "Z"|>], {k, n}]]
ghzDetachedReference[n_, w_, m_, paritySign_] := With[{others = Delete[Range[n], w]},
    PauliStabilizer[Join[
        Table[letterString[n, <|others[[i]] -> "Z", others[[i + 1]] -> "Z"|>], {i, n - 2}],
        {signPrefix[paritySign] <> letterString[n, AssociationThread[others -> "X"]], signPrefix[1 - 2 m] <> letterString[n, <|w -> "X"|>]}
    ]]
]

(* The cuts examined: one other wire, the measured wire, both, the first half,   *)
(* all other wires, the whole register.                                          *)
ghzCuts[n_, w_] := With[{others = Delete[Range[n], w]},
    {{others[[1]]}, {w}, {w, others[[1]]}, Range[Floor[n / 2]], others, Range[n]}
]

ghzMeasurementReport[n_, w_] := With[{
    ghz = QuantumCircuitOperator["GHZ"[n]][PauliStabilizer[n]],
    others = Delete[Range[n], w],
    cuts = ghzCuts[n, w]
},
    With[{mz = QuantumMeasurementOperator["Z", {w}][ghz], mx = QuantumMeasurementOperator["X", {w}][ghz]},
        {n, w,
            Sort @ Keys @ mz,
            Table[
                {DeleteDuplicates @ Table[mz[m]["Expectation", letterString[n, <|k -> "Z"|>]], {k, n}],
                    outcomeBits[mz[m], "Z", {others[[1]]}], outcomeBits[mz[m], "Z", {others[[-1]]}],
                    closedFormOverlap[mz[m], ghzProductReference[n, m]],
                    mz[m]["Expectation", StringRepeat["X", n]],
                    Table[mz[m]["Entropy", cut], {cut, cuts}],
                    tableauSymplecticQ[mz[m]]},
                {m, Sort @ Keys @ mz}
            ],
            Sort @ Keys @ mx,
            Table[
                {closedFormOverlap[mx[m], ghzDetachedReference[n, w, m, (-1)^m]],
                    closedFormOverlap[mx[m], ghzDetachedReference[n, w, m, -(-1)^m]],
                    DeleteDuplicates @ Table[mx[m]["Expectation", letterString[n, <|others[[i]] -> "Z", others[[i + 1]] -> "Z"|>]], {i, n - 2}],
                    outcomeBits[mx[m], "ZZ", others[[1 ;; 2]]],
                    mx[m]["Expectation", letterString[n, AssociationThread[others -> "X"]]],
                    Table[mx[m]["Entropy", cut], {cut, cuts}],
                    tableauSymplecticQ[mx[m]]},
                {m, Sort @ Keys @ mx}
            ]}
    ]
]

ghzMeasurementPrediction[n_, w_] := With[{cuts = ghzCuts[n, w]},
    {n, w,
        {0, 1},
        Table[{{1 - 2 m}, {m}, {m}, 1, 0, Table[0, {Length[cuts]}], True}, {m, 0, 1}],
        {0, 1},
        Table[{1, 0, {1}, {0}, (-1)^m, Table[Boole[0 < Length[DeleteCases[cut, w]] < n - 1], {cut, cuts}], True}, {m, 0, 1}]}
]

$ghzMeasuredWires = Flatten[Table[{n, w}, {n, {9, 50}}, {w, {1, Ceiling[n / 3], n}}], 1];

VerificationTest[
    ghzMeasurementReport @@@ $ghzMeasuredWires,
    ghzMeasurementPrediction @@@ $ghzMeasuredWires,
    {},
    TestID -> "Protocol-GHZ9-GHZ50-ZMeasurementCollapsesAll-XMeasurementDetachesOne-ThreeWires"
]
