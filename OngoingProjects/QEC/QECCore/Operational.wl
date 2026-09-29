(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[$QECDenseQubitLimit]

PackageScope[codeEncoderMatrix]
PackageScope[codeEncoderOperator]
PackageScope[codeCodewords]
PackageScope[codeCodespace]
PackageScope[codeSyndromeProjectors]
PackageScope[codeSyndromeMeasurement]
PackageScope[codeRecovery]
PackageScope[codeLogicalClassPaulis]
PackageScope[codeLogicalPauliProbabilities]
PackageScope[codeLogicalChannel]
PackageScope[codeKnillLaflammeMatrix]
PackageScope[pauliRowMatrix]
PackageScope[denseAllowedQ]


(* ============================================================================ *)
(* The code as the objects quantum theory uses for error correction.            *)
(*                                                                              *)
(* The operational formulation (Knill-Laflamme-Viola; Rahn-Doherty-Mabuchi) is  *)
(* built from five things QuantumFramework already ships: an encoding ISOMETRY  *)
(* V from k logical qubits into n physical ones, the code-space PROJECTOR       *)
(* P = V V^dagger, noise as a CHANNEL, syndrome extraction as an INSTRUMENT     *)
(* {M_s}, and recovery as a syndrome-indexed family of channels {R_s}.  The     *)
(* whole cycle composes to one more channel, the effective logical channel      *)
(*                                                                              *)
(*     L = D o ( Sum_s R_s o M_s ) o N o E .                                    *)
(*                                                                              *)
(* This file makes a QECCode present those objects, so the cycle composes with  *)
(* the rest of the framework.  It is a change of FACE, not of engine: every     *)
(* object here is either small, or computed by machinery that already exists.   *)
(*                                                                              *)
(* THE LOGICAL CHANNEL IS NOT BUILT BY COMPOSING PHYSICAL CHANNELS, and this is *)
(* the one decision in the file that matters.  Taken literally, the formula     *)
(* above puts a 4^n x 4^n superoperator in the middle -- about 2.7 x 10^8       *)
(* entries on the 7-qubit code, each a polynomial in p -- which neither fits    *)
(* nor keeps the rate symbolic.  But for Pauli noise on a stabilizer code the   *)
(* logical channel is ALWAYS a Pauli channel on the k logical qubits: an error  *)
(* E has syndrome s, the recovery R_s is a Pauli, so the residue R_s E lies in  *)
(* the normalizer -- a logical Pauli times a stabilizer -- and its whole effect *)
(* on the encoded state is that logical Pauli.  The channel is therefore 4^k    *)
(* numbers, the probability of each logical residue, and those are exactly the  *)
(* coset probabilities the maximum-likelihood decoder already computes, exactly *)
(* and symbolically.  So the logical channel comes from the coset engine, and   *)
(* the rate is 1 - q_I, a functional of it, with nothing lost: on the Steane     *)
(* code the two agree identically.                                             *)
(*                                                                              *)
(* THE ENCODER IS STATED IN THE CODE'S OWN LOGICAL BASIS.  The encoding circuit *)
(* (Encoder.wl) is a valid encoder, but its labelling of the logical basis is   *)
(* its own: on the 5-qubit code it sends |0> to the -1 eigenstate of the Zbar   *)
(* that code["LogicalZ"] reports.  An isometry whose codewords disagree with    *)
(* the code's declared logical operators would make every property built on it *)
(* quietly relabelled, so the codewords are rebuilt here from the circuit's     *)
(* |0...0> output the textbook way: |0_L> is that output moved into the +1      *)
(* eigenspace of every Zbar_j, and |c_L> = Xbar^c |0_L>.  Then V^dagger Zbar V  *)
(* = Z and V^dagger Xbar V = X exactly, on every code, and a test says so.      *)
(*                                                                              *)
(* THE KNILL-LAFLAMME MATRIX IS COMPUTED SYMPLECTICALLY, NOT DENSELY.  For      *)
(* Pauli errors E_a, E_b the overlap <W_i| E_a^dagger E_b |W_j> is decided by   *)
(* where the Pauli E_a^dagger E_b sits: anticommuting with some check gives 0;  *)
(* in the stabilizer with phase i^c gives i^c delta_ij; a nontrivial logical    *)
(* breaks the condition.  No 2^n vector is formed, so it runs at any size the   *)
(* rest of the layer does, and the dense definition is checked against it in    *)
(* the tests.  It is also the property the bosonic sibling is specified around  *)
(* (Bosonic-QEC-PartA-Spec.md: "KnillLaflammeMatrix", "CorrectableQ"), so both  *)
(* branches answer "does this code correct this error set" with one spelling.   *)
(*                                                                              *)
(* WHAT IS DENSE, and guarded: the encoder, the codewords, the projector and    *)
(* the syndrome instrument are 2^n-sized.  $QECDenseQubitLimit bounds them, and *)
(* the refusal names the cheap route instead of just refusing.                  *)
(*                                                                              *)
(* References: Knill, Laflamme, Viola, PRL 84, 2525 (2000); Rahn, Doherty,      *)
(* Mabuchi, PRA 66, 032304 (2002); Got26 ch. 2-3 (the stabilizer formalism and  *)
(* the normalizer); QEC-API-Audit-and-Redesign.md sec. 3.                       *)
(* ============================================================================ *)

$QECDenseQubitLimit::usage = "$QECDenseQubitLimit is the largest number of physical qubits for which a code materializes its dense operational objects (\"Encoder\", \"Codewords\", \"Codespace\", \"SyndromeMeasurement\", \"Recovery\").";

$QECDenseQubitLimit = 10;

QECCode::dense = "`1` materializes a dense object on `2` qubits, above $QECDenseQubitLimit = `3`. The logical channel, the Knill-Laflamme matrix and the detector model carry the same information without it; raise $QECDenseQubitLimit to force it.";
QECCode::nolevel = "The logical channel is built from code-capacity noise, where it is exact. At the `1` level the rate, its acceptance and its failure are read from the detector model: QECDetectorModel[code, noise].";
QECCode::klerrors = "`1` is not an error set. Give a list of Pauli strings or rows on the code's qubits, or \"Weight\"[t] for every Pauli of weight at most t.";
QECCode::zbarsign = "The encoding circuit's output is not an eigenstate of Zbar_`1`; the codewords cannot be put in the code's logical basis.";


denseAllowedQ[a_Association, prop_String] := If[
    a["Qubits"] <= $QECDenseQubitLimit,
    True,
    Message[QECCode::dense, prop, a["Qubits"], $QECDenseQubitLimit]; False
]


(* ---- a Pauli row as a matrix ---- *)

(* i^e times the tensor product of E(x, z) = i^(x z) X^x Z^z, qubit 1 most significant,
   which is the framework's own ordering.  Sparse, since every Pauli is. *)
pauliRowMatrix[row_List] := With[{n = pauliQubits[row]},
    I^Last[row] Fold[
        KroneckerProduct,
        Table[
            SparseArray[I^(row[[q]] row[[n + q]]) *
                MatrixPower[PauliMatrix[1], row[[q]]] . MatrixPower[PauliMatrix[3], row[[n + q]]]],
            {q, n}
        ]
    ]
]


(* ---- the isometry and its codewords ---- *)

(* The circuit's |0...0> output, moved into the +1 eigenspace of every Zbar_j, and the
   other codewords generated from it by the Xbar's.  Column c (0-based) is the codeword
   for logical input c, logical qubit 1 most significant. *)
codeEncoderMatrix[a_Association] := codeEncoderMatrix[a] = Module[
    {n = a["Qubits"], k = codeLogicalQubits[a], u, v0, xs, zs, sign},

    u = Wolfram`QuantumFramework`QuantumCircuitOperator[
        Join[Table["I" -> q, {q, n}], codeEncodingGates[a]]
    ]["QuantumOperator"]["MatrixRepresentation"];
    v0 = Normal[u[[All, 1]]];

    If[k === 0, Return[Transpose[{v0}]]];

    xs = pauliRowMatrix /@ codeLogicalVectors[a]["X"];
    zs = pauliRowMatrix /@ codeLogicalVectors[a]["Z"];

    Do[
        sign = Simplify[Conjugate[v0] . (zs[[j]] . v0)];
        Which[
            sign === 1, Null,
            sign === -1, v0 = xs[[j]] . v0,
            True, Message[QECCode::zbarsign, j]; Return[$Failed, Module]
        ],
        {j, k}
    ];

    Transpose @ Table[
        Fold[#2 . #1 &, v0, Pick[xs, IntegerDigits[c, 2, k], 1]],
        {c, 0, 2^k - 1}
    ]
]

codeEncoderOperator[a_Association] := With[{v = codeEncoderMatrix[a]},
    If[ v === $Failed,
        $Failed,
        Wolfram`QuantumFramework`QuantumOperator[
            v, {Range[a["Qubits"]], Range[Max[codeLogicalQubits[a], 1]]}
        ]
    ]
]

codeCodewords[a_Association] := With[{v = codeEncoderMatrix[a]},
    If[v === $Failed, $Failed, Wolfram`QuantumFramework`QuantumState /@ Transpose[v]]
]

(* P = V V^dagger.  The tests check it against the product of (1 + g)/2 over the
   generators, which is the same projector reached without the encoder. *)
codeCodespace[a_Association] := With[{v = codeEncoderMatrix[a]},
    If[ v === $Failed,
        $Failed,
        Wolfram`QuantumFramework`QuantumOperator[v . ConjugateTranspose[v], Range[a["Qubits"]]]
    ]
]


(* ---- the syndrome instrument ---- *)

(* One projector per syndrome, P_s = Prod_i (1 + (-1)^(s_i) g_i) / 2.  Outcome i of the
   measurement is the syndrome IntegerDigits[i - 1, 2, m], the order code["Syndrome", e]
   reports its bits in. *)
codeSyndromeProjectors[a_Association] := codeSyndromeProjectors[a] = Module[
    {n = a["Qubits"], gens = pauliRowMatrix /@ codeVectors[a], one},
    one = IdentityMatrix[2^n, SparseArray];
    Table[
        Fold[#1 . ((one + (-1)^#2[[1]] #2[[2]]) / 2) &, one, Transpose[{s, gens}]],
        {s, Tuples[{0, 1}, Length[gens]]}
    ]
]

codeSyndromeMeasurement[a_Association] := With[{n = a["Qubits"]},
    Wolfram`QuantumFramework`QuantumMeasurementOperator[
        Wolfram`QuantumFramework`QuantumOperator[#, Range[n]] & /@ codeSyndromeProjectors[a],
        Range[n]
    ]
]


(* ---- recovery ---- *)

(* The unitary channel of the correction the code's decoder applies for syndrome s:
   one Kraus operator, a Pauli, so it goes in as that operator (see krausChannel). *)
codeRecovery[a_Association, syn_List] := With[{corr = codeDecode[a, syn]},
    If[ corr === $Failed || MissingQ[corr],
        corr,
        krausChannel[{pauliRowMatrix[pauliVector[corr]]}, Range[a["Qubits"]]]
    ]
]


(* ---- the logical channel, from the coset engine ---- *)

(* Which logical Pauli each class key names.  Read off by labelling a representative of
   every logical Pauli with the same function that labels errors, so the mapping cannot
   drift from the one the decoder and the rate use. *)
codeLogicalClassPaulis[a_Association] := codeLogicalClassPaulis[a] = Module[
    {n = a["Qubits"], k = codeLogicalQubits[a], xs, zs, reps, keys},
    xs = codeLogicalVectors[a]["X"];
    zs = codeLogicalVectors[a]["Z"];
    reps = Table[
        With[{ex = IntegerDigits[c, 2, 2 k][[1 ;; k]], ez = IntegerDigits[c, 2, 2 k][[k + 1 ;; 2 k]]},
            {Join[ex, ez, {0}],
             Fold[pauliProduct, pauliIdentity[n], Join[Pick[xs, ex, 1], Pick[zs, ez, 1]]]}
        ],
        {c, 0, 4^k - 1}
    ];
    keys = Last @ codeLabelKeys[a, reps[[All, 2]]];
    AssociationThread[keys, pauliString[#[[1]]] & /@ reps]
]

(* The 4^k probabilities of the logical residue, keyed by the logical Pauli (phase
   dropped: a Pauli channel does not see it).  The residue of an error in class c that
   the decoder assigns to class d is the logical Pauli c XOR d -- labels are linear.

   A syndrome the decoder does not cover (a minimum-weight table of limited reach) has
   no recovery, so its weight is not a logical Pauli; it is kept apart under
   "Undecoded", which is exactly the part QECLogicalErrorRate counts as failure on top
   of the wrong-class part.  So 1 - q_I equals the rate for every decoder. *)
codeLogicalPauliProbabilities[a_Association, noise_Association, decoder_Association] := Module[
    {names = logicalLetters /@ codeLogicalClassPaulis[a], probs, undecoded},
    probs = Merge[
        KeyValueMap[
            With[{d = Lookup[decoder, #1[[1]], Missing[]]},
                If[MissingQ[d], "Undecoded", BitXor[#1[[2]], d]] -> #2
            ] &,
            cosetProbabilities[a, noise]
        ],
        Total
    ];
    undecoded = Lookup[probs, "Undecoded", 0];
    Join[
        KeySort @ Association @ KeyValueMap[#2 -> Simplify[Lookup[probs, #1, 0]] &, names],
        If[PossibleZeroQ[undecoded], <||>, <|"Undecoded" -> Simplify[undecoded]|>]
    ]
]

logicalLetters[s_String] := StringDelete[s, "+" | "-" | "i"]

QECCode::undecoded = "The decoder leaves syndromes carrying probability `1` without a recovery, so the cycle does not return to the code space and has no logical channel. Use a decoder that covers every syndrome (the default on a code this size), or read that weight from code[\"LogicalPauliProbabilities\", noise].";

(* sqrt(q_L) L for every logical Pauli with nonzero weight, on the k logical qubits. *)
codeLogicalChannel[a_Association, noise_Association, decoder_Association] := Module[
    {probs = codeLogicalPauliProbabilities[a, noise, decoder], k = codeLogicalQubits[a], kraus},
    If[ KeyExistsQ[probs, "Undecoded"],
        Message[QECCode::undecoded, probs["Undecoded"]]; Return[$Failed]
    ];
    kraus = KeyValueMap[
        If[ PossibleZeroQ[#2], Nothing,
            Sqrt[#2] Normal[pauliRowMatrix[pauliVector[#1]]]
        ] &,
        probs
    ];
    krausChannel[kraus, Range[k]]
]


(* ---- Knill-Laflamme, symplectically ---- *)

knillLaflammeErrors[a_Association, "Weight"[t_Integer ? NonNegative]] :=
    Catenate @ Table[weightKVectors[a["Qubits"], w], {w, 0, Min[t, a["Qubits"]]}]

knillLaflammeErrors[a_Association, errs_List] := With[{rows = pauliVector /@ errs},
    If[ MemberQ[rows, $Failed] || ! AllTrue[rows, Length[#] === 2 a["Qubits"] + 1 &],
        $Failed,
        rows
    ]
]

knillLaflammeErrors[_, _] := $Failed

(* E^dagger for E = i^e E(x, z): the Hermitian part is its own dagger, the phase
   conjugates. *)
pauliDagger[row_List] := MapAt[Mod[-#, 4] &, row, -1]

(* h_ab, or Missing naming the first pair whose product is a nontrivial logical. *)
codeKnillLaflammeMatrix[a_Association, spec_] := Module[{rows, entry},
    rows = knillLaflammeErrors[a, spec];
    If[rows === $Failed, Message[QECCode::klerrors, spec]; Return[$Failed]];
    entry[ea_, eb_] := With[{prod = pauliProduct[pauliDagger[ea], eb]},
        If[ Total[codeSyndromeVector[a, prod]] > 0,
            0,
            With[{elt = codeStabilizerElement[a, prod]},
                If[ MissingQ[elt],
                    Throw[Missing["NotCorrectable", <|
                        "Errors" -> pauliString /@ {ea, eb},
                        "Product" -> pauliString[prod]|>], "kl"],
                    I^Mod[Last[prod] - Last[elt], 4]
                ]
            ]
        ]
    ];
    Catch[Table[entry[ea, eb], {ea, rows}, {eb, rows}], "kl"]
]


(* ---- the properties ---- *)

(* The noise a property is given: a QECNoiseModel, or a single-qubit Pauli
   QuantumChannel, which QECNoiseModel reads back into one (and refuses otherwise). *)
operationalNoise[QECNoiseModel[noise_Association]] := noise
operationalNoise[qc_Wolfram`QuantumFramework`QuantumChannel] := Replace[QECNoiseModel[qc], {
    QECNoiseModel[noise_Association] :> noise,
    _ -> $Failed
}]
operationalNoise[_] := $Failed

Options[codeOperationalChannel] = {"Decoder" -> Automatic, "DecoderReach" -> Automatic};

(* The checks the exact rate makes, in the same order and with the same messages, so
   the channel and the rate refuse the same inputs. *)
operationalChannelData[a_Association, spec_, opts : OptionsPattern[codeOperationalChannel]] := Module[
    {noise = operationalNoise[spec], reach, decoder},
    Which[
        noise === $Failed, Return[$Failed],
        codeLogicalQubits[a] === 0, Message[QECLogicalErrorRate::nological]; Return[$Failed],
        noiseLevel[noise] =!= "CodeCapacity",
            Message[QECCode::nolevel, noiseLevel[noise]]; Return[Missing["NotAvailable", noiseLevel[noise]]],
        ! exactEnumerableQ[a],
            Message[QECLogicalErrorRate::toobig, a["Qubits"], 4^a["Qubits"], $QECExactEnumerationLimit];
            Return[$Failed]
    ];
    reach = decoderReachFor[a, OptionValue[codeOperationalChannel, {opts}, "DecoderReach"]];
    decoder = decoderFor[a, noise, OptionValue[codeOperationalChannel, {opts}, "Decoder"], reach];
    If[decoder === $Failed, $Failed, {noise, decoder}]
]

QECCode[a_Association]["Encoder"] := If[denseAllowedQ[a, "Encoder"], codeEncoderOperator[a], $Failed]
QECCode[a_Association]["Codewords"] := If[denseAllowedQ[a, "Codewords"], codeCodewords[a], $Failed]
QECCode[a_Association]["Codespace"] := If[denseAllowedQ[a, "Codespace"], codeCodespace[a], $Failed]
QECCode[a_Association]["SyndromeMeasurement"] :=
    If[denseAllowedQ[a, "SyndromeMeasurement"], codeSyndromeMeasurement[a], $Failed]

QECCode[a_Association]["Recovery", syn_List] := If[denseAllowedQ[a, "Recovery"], codeRecovery[a, syn], $Failed]

QECCode[a_Association]["LogicalPauliProbabilities", noise_, opts : OptionsPattern[codeOperationalChannel]] :=
    Replace[operationalChannelData[a, noise, opts], {
        {n_Association, d_Association} :> codeLogicalPauliProbabilities[a, n, d]
    }]

QECCode[a_Association]["LogicalChannel", noise_, opts : OptionsPattern[codeOperationalChannel]] :=
    Replace[operationalChannelData[a, noise, opts], {
        {n_Association, d_Association} :> codeLogicalChannel[a, n, d]
    }]

QECCode[a_Association]["KnillLaflammeMatrix", errors_] := codeKnillLaflammeMatrix[a, errors]

QECCode[a_Association]["CorrectableQ", errors_] := With[{h = codeKnillLaflammeMatrix[a, errors]},
    If[h === $Failed, $Failed, ! MissingQ[h]]
]
