(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECPauliMeasurement]

PackageScope[catMeasureInstructions]
PackageScope[pauliMeasureInstructions]
PackageScope[pauliMeasureRepetitions]
PackageScope[pauliMeasureOutcome]
PackageScope[pauliMeasureDataWeights]
PackageScope[pauliMeasureSkipIndices]
PackageScope[codeExtraction]
PackageScope[recordSyndromes]
PackageScope[measurementData]


(* ============================================================================ *)
(* Fault-tolerant measurement of a Pauli through cat states.                    *)
(*                                                                              *)
(* Circuit.wl measures a check with one bare ancilla shared by every data qubit *)
(* of the check, and that is the construction Got26 titles "Non-Fault-Tolerant  *)
(* Measurement of Paulis": one fault on the ancilla reaches several data qubits *)
(* at once, which is the hook error, and it is what turns a distance-three      *)
(* code's rate from p^2 into O(p).                                             *)
(*                                                                              *)
(* The fix is to make the controlled-P transversal.  Cat.wl builds and verifies *)
(* the m = wt(P) qubit ancilla; this file spends it.                           *)
(*                                                                              *)
(* ONE ATTEMPT (Got26 sec. 12.1.2).  Cat qubit j interacts with exactly one     *)
(* data qubit, so a fault on it becomes at most one data error.  Reading out    *)
(* the eigenvalue takes a Hadamard transform rather than disentangling the cat: *)
(*                                                                              *)
(*     H^(tensor m) (|0...0> + (-1)^b |1...1>) = sum over wt(x) = b (mod 2)     *)
(*                                                                              *)
(* so measuring all m qubits and taking the PARITY of the weight gives b        *)
(* (eqs. 12.1-12.3).  m bits of record, one bit of answer.                      *)
(*                                                                              *)
(* THE SANDWICH IS THE SAME ONE Circuit.wl USES, in the other direction.  To    *)
(* apply a controlled-L, conjugate by the U that sends L to Z:                  *)
(*                                                                              *)
(*     controlled-L = (I tensor U^-1) CZ (I tensor U),      U L U^-1 = Z        *)
(*                                                                              *)
(* which in time order is U, then CZ, then U^-1 -- exactly rotation[basisGate], *)
(* the ladder, rotation[unbasisGate] from generatorInstructions, with CZ where  *)
(* that had CNOT.  For L = X this is H CZ H = CNOT, so nothing is wasted; the   *)
(* rotations are shared rather than re-derived, and a change to the letter      *)
(* vocabulary cannot leave the two constructions disagreeing.                   *)
(*                                                                              *)
(* REPETITION AND MAJORITY (Got26 sec. 12.1.4, Theorem 12.1).  One attempt      *)
(* still fails the measurement correctness property: a bit flip on a cat qubit  *)
(* just before its readout flips the reported b, and a phase error earlier does *)
(* the same.  The answer is 2t+1 attempts with a FRESH cat each time and a      *)
(* majority vote, since one fault can corrupt only one cat.                     *)
(*                                                                              *)
(* WHAT THIS GADGET DOES NOT YET DO, and it is the honest boundary of this      *)
(* file.  Repetition handles faults in the ancilla.  It does NOT handle a fault *)
(* in the DATA: if the codeword carries an error E with E P = -P E, then every  *)
(* repetition reads the flipped eigenvalue and the majority is confidently      *)
(* wrong -- "no matter how many times we repeat it" (sec. 12.1.4).  Got26's fix *)
(* is an FTEC sub-gadget after each repetition, and figure 12.5 is the gadget   *)
(* with those interspersed.  That sub-gadget is Steane EC, which is the next    *)
(* file, so "ErrorCorrection" is an option here that currently defaults to      *)
(* None and splices nothing.  QECPauliMeasurement therefore satisfies the MCP   *)
(* against ancilla faults only, and says so: "MeasurementCorrectQ" is False     *)
(* while the slot is empty, rather than quietly claiming a property it lacks.   *)
(*                                                                              *)
(* Note also that the choice of REPRESENTATIVE matters, and Got26 flags it as   *)
(* an example of the general phenomenon: Z-bar for the five-qubit code can be   *)
(* taken as ZZZZZ, with which a single X_i anticommutes, or as another element  *)
(* of the same logical coset with which it commutes.  Equivalent operators,     *)
(* inequivalent protocols.  The Pauli handed in here is used as given.          *)
(*                                                                              *)
(* References: Got26 sec. 12.1.2 (cat measurement, eqs. 12.1-12.3, fig. 12.2),  *)
(* sec. 12.1.4 (repetition, fig. 12.5), sec. 12.1.5 and Theorem 12.1 (the MCP), *)
(* sec. 12.2 (Shor EC, the matching EC gadget).                                 *)
(* ============================================================================ *)

QECPauliMeasurement::usage = "QECPauliMeasurement[code, P] gives the fault-tolerant gadget that measures the Pauli P on a block of the code, through verified cat states repeated with a majority vote.\nQECPauliMeasurement[code, P, r] uses r repetitions instead of the 2t+1 the code's distance calls for.\nThe options \"Pairs\" and \"CatRepetitions\" are passed to the cat state, and \"ErrorCorrection\" names the sub-gadget spliced between repetitions.\nm[prop] gives a property; m[\"Properties\"] lists them.";

QECPauliMeasurement::pauli = "`1` is not a Pauli on `2` qubits.";
QECPauliMeasurement::identity = "The identity has no eigenvalue to measure.";
QECPauliMeasurement::reps = "The number of repetitions must be a positive integer; got `1`.";
QECPauliMeasurement::ec = "\"ErrorCorrection\" -> `1` is not a known sub-gadget; use None or \"Steane\".";
QECPauliMeasurement::eccss = "Steane error correction needs a CSS code, and this one is not. Use \"ErrorCorrection\" -> None, or a CSS code; Shor error correction, which has no such restriction, is not built.";
QECPauliMeasurement::noprop = "`1` is not a property of QECPauliMeasurement. Use m[\"Properties\"] for the list.";


(* ---- one attempt ---- *)

(* The transversal controlled-P, then the Hadamard transform and the readout.  The
   cat qubits are handed in already prepared and verified; this is only the part
   that spends them.  cat[[j]] pairs with the j-th qubit of P's support, which is
   the transversality the whole construction exists for. *)
catMeasureInstructions[v_List, n_Integer, cat_List] := With[
    {support = Select[Range[n], letterAt[v, n, #] =!= {0, 0} &]},
    Join[
        rotation[basisGate, v, n, support],
        Table[{"CZ", cat[[j]], support[[j]]}, {j, Length[support]}],
        rotation[unbasisGate, v, n, support],
        Table[{"H", q}, {q, cat}],
        Table[{"M", q}, {q, cat}]
    ]
]


(* ---- the whole gadget ---- *)

(* A code correcting t errors needs 2t+1 repetitions (Theorem 12.1).  t is taken
   from the distance as usual, and a distance-one code gets one repetition, which
   is the honest answer: there is nothing to be fault tolerant about. *)
pauliMeasureRepetitions[d_Integer] := 2 Floor[(d - 1) / 2] + 1

(* Figure 12.5: repetitions of (fresh verified cat, transversal measurement,
   readout), with an error-correction sub-gadget between them.  The cat qubits are
   reset by each preparation, so the repetitions reuse them rather than needing
   2t+1 separate cats worth of hardware.

   Catenate, not Flatten at level 2: the Table has ONE iterator over a body that is
   already a flat list of instructions, so exactly one level has to come off.  The
   sibling functions in Cat.wl and DetectorModel.wl take Flatten[..., 2] because
   their Tables carry two iterators.  Counting the iterators is the rule; copying
   the neighbouring call is how this goes wrong, and it did here first time. *)
pauliMeasureInstructions[
    v_List, n_Integer, cat_List, check_Integer, pairs_List,
    catReps_Integer, reps_Integer, ec_
] := Catenate @ Table[
    Join[
        catInstructions[cat, check, pairs, catReps],
        catMeasureInstructions[v, n, cat],
        Replace[ec, None -> {}]
    ],
    {reps}
]


(* ---- reading the answer ---- *)

(* One repetition contributes m record bits whose parity is its reported eigenvalue
   (eqs. 12.1-12.3); the gadget's answer is the majority over repetitions.  With an
   odd number of repetitions there is always a majority -- which is exactly why
   Got26 can use 2t+1 here while sec. 12.2.2 rejects majority voting for a whole
   syndrome, where the 2^(n-k) possible values may have no majority at all. *)
pauliMeasureOutcome[record_List, m_Integer, reps_Integer] := With[
    {bits = Mod[Total /@ Partition[Take[record, m reps], m], 2]},
    Boole[Total[bits] > reps / 2]
]


(* ---- the proof obligation ---- *)

(* The residual weight left on the data by every single fault the gadget admits.
   This is the property the whole construction exists for, and it needs two
   corrections before the number means anything.

   FIRST, the frame is read modulo P.  X^(tensor m) stabilises the cat state, so a
   cat X-pattern and its complement are the same physical error -- the degeneracy
   catCanonicalPattern already exploits on the cat side.  Through the transversal
   CZ those two patterns differ on the data by a Z on every qubit of the support,
   which after the un-rotation is exactly P.  The frame propagator tracks literal
   Paulis and cannot know the ancilla's own stabilizer, so without this quotient it
   reports a weight-m residual for a fault that does nothing at all, and a
   transversal gadget comes out looking worse than the bare-ancilla one it replaces.

   SECOND, only ACCEPTED runs count.  A fault early in the cat preparation does
   propagate along the chain and would reach several data qubits -- that is what
   checking is for, and the check fires, so the shot is discarded.  Rejected runs
   are therefore filtered out rather than counted against the gadget; what is
   claimed is the conditional statement, which is the one Got26 Theorem 12.1 needs.

   THIRD, when an error-correction sub-gadget is spliced in, its own ancilla
   PREPARATION is skipped, for exactly the reason ErrorCorrection.wl skips it
   there: that preparation is the code's non-fault-tolerant encoder, and making it
   fault tolerant is chapter 13 work.  Its cost is reported by the EC gadget's
   "PreparationDataWeights" rather than charged to the measurement.

   With all three in place a single fault leaves at most one data error, and
   dropping the herald filter is enough to see weight two reappear, so the checks
   are visibly load-bearing rather than decorative. *)
pauliMeasureDataWeights[instr_List, nq_Integer, n_Integer, v_List] :=
    pauliMeasureDataWeights[instr, nq, n, v, {}]

pauliMeasureDataWeights[instr_List, nq_Integer, n_Integer, v_List, skip_List] := Module[
    {rows, weight},
    weight[r_] := Total @ Map[Max, Transpose[{Take[r, n], Take[r, -n]}]];
    rows = Table[
        With[{run = framePropagate[instr, nq, {{i, q, pauli}}]},
            If[ Total[run["Heralds"]] > 0,
                Nothing,
                Join[Take[run["Frame"][[1]], n], Take[run["Frame"][[2]], n]]
            ]
        ],
        {i, Complement[Range[0, Length[instr]], skip]}, {q, nq}, {pauli, $oneQubitPaulis}
    ];
    DeleteDuplicates[Min[weight[#], weight[BitXor[#, Take[v, 2 n]]]] & /@ Flatten[rows, 2]]
]


(* ---- construction ---- *)

Options[QECPauliMeasurement] = {
    "Pairs" -> Automatic,
    "CatRepetitions" -> 1,
    "ErrorCorrection" -> None
};

QECPauliMeasurement[code_QECCode, p_, opts : OptionsPattern[]] :=
    QECPauliMeasurement[code, p, pauliMeasureRepetitions[code["Distance"]], opts]

QECPauliMeasurement[code_QECCode, p_, reps_, opts : OptionsPattern[]] := Module[
    {n = code["Qubits"], v, support, m, cat, check, pairs, catReps, ec, ecGadget, ecInstr},
    ec = OptionValue["ErrorCorrection"];
    catReps = OptionValue["CatRepetitions"];
    (* QECPauliQ is the message-free guard the Pauli layer provides for exactly this
       question, so asking it is what the house rule wants instead of silencing
       QECPauliVector's complaint.  The size and phase conditions are separate because
       they are different refusals: a Pauli on the wrong number of qubits, and a Pauli
       carrying a phase, whose Hermitian representative is what this gadget measures. *)
    v = If[QECPauliQ[p], QECPauliVector[p], $Failed];
    Which[
        v === $Failed || Length[v] =!= 2 n + 1 || Last[v] =!= 0,
            Message[QECPauliMeasurement::pauli, p, n]; $Failed,
        ! (IntegerQ[reps] && reps > 0),
            Message[QECPauliMeasurement::reps, reps]; $Failed,
        ! MemberQ[{None, "Steane"}, ec],
            Message[QECPauliMeasurement::ec, ec]; $Failed,
        ec === "Steane" && ! codeCSSQ[First[code]],
            Message[QECPauliMeasurement::eccss]; $Failed,
        True,
            support = Select[Range[n], letterAt[v, n, #] =!= {0, 0} &];
            m = Length[support];
            If[ m === 0,
                Message[QECPauliMeasurement::identity]; Return[$Failed]
            ];
            cat = Range[n + 1, n + m];
            check = n + m + 1;
            pairs = Replace[OptionValue["Pairs"], Automatic :> catCheckPairs[m]];
            (* The EC sub-gadget brings its own block of n ancillas, placed past the
               cat's check qubit.  It is spliced AFTER each repetition and not before
               the first: Got26 sec. 12.1.4 notes that one at the start would be
               redundant, since the standard FT simulation (Def 10.6) already puts an
               EC between every adjacent pair of gadgets. *)
            ecGadget = If[ec === "Steane", QECErrorCorrection[code, check], None];
            ecInstr = If[ecGadget === None, None, ecGadget["Instructions"]];
            QECPauliMeasurement[<|
                "Code" -> First[code],
                "Pauli" -> Take[v, 2 n],
                "Support" -> support,
                "DataQubits" -> n,
                "CatQubits" -> cat,
                "CheckQubit" -> check,
                "Pairs" -> pairs,
                "CatRepetitions" -> catReps,
                "Repetitions" -> reps,
                "ErrorCorrection" -> ec,
                "Qubits" -> If[ecGadget === None, check, ecGadget["Qubits"]],
                "Instructions" -> pauliMeasureInstructions[
                    v, n, cat, check, pairs, catReps, reps, ecInstr
                ]
            |>]
    ]
]


(* ---- properties ---- *)

$measurementProperties = {
    "Instructions", "Pauli", "Support", "Weight", "Repetitions", "CatRepetitions",
    "DataQubits", "CatQubits", "CheckQubit", "Qubits", "Pairs", "Code", "Cat",
    "ErrorCorrection", "Heralds", "Measurements", "Depth", "InstructionCount",
    "GateCounts", "QuantumCircuitOperator", "Diagram",
    "TransversalQ", "DataWeights", "MaxDataWeight", "MeasurementCorrectQ",
    "OpenAssumptions", "Properties"
};

measurementData[QECPauliMeasurement[a_Association]] := a

QECPauliMeasurement[_Association]["Properties"] := $measurementProperties

QECPauliMeasurement[a_Association][prop : ("Instructions" | "Support" | "Repetitions" |
    "CatRepetitions" | "DataQubits" | "CatQubits" | "CheckQubit" | "Qubits" | "Pairs" |
    "ErrorCorrection")] := a[prop]

QECPauliMeasurement[a_Association]["Pauli"] := QECPauliString[Append[a["Pauli"], 0]]
QECPauliMeasurement[a_Association]["Weight"] := Length[a["Support"]]
QECPauliMeasurement[a_Association]["Code"] := QECCode[a["Code"]]
QECPauliMeasurement[a_Association]["Heralds"] := gadgetHeralds[a["Instructions"]]
QECPauliMeasurement[a_Association]["Measurements"] := gadgetMeasurements[a["Instructions"]]
QECPauliMeasurement[a_Association]["InstructionCount"] := Length[a["Instructions"]]
QECPauliMeasurement[a_Association]["Depth"] := gadgetDepth[a["Instructions"], a["Qubits"]]
QECPauliMeasurement[a_Association]["QuantumCircuitOperator"] :=
    gadgetCircuitOperator[a["Instructions"]]
QECPauliMeasurement[a_Association]["Diagram"] := gadgetDiagram[a["Instructions"]]
QECPauliMeasurement[a_Association]["GateCounts"] := Counts[First /@ a["Instructions"]]

QECPauliMeasurement[a_Association]["Cat"] :=
    QECCatState[a["CatQubits"], a["CheckQubit"],
        "Pairs" -> a["Pairs"], "Repetitions" -> a["CatRepetitions"]]

(* Transversality is a structural claim and is checked structurally: no cat qubit
   shares a two-qubit gate with more than one data qubit. *)
QECPauliMeasurement[a_Association]["TransversalQ"] := With[
    {n = a["DataQubits"]},
    AllTrue[
        a["CatQubits"],
        Function[c,
            Length @ DeleteDuplicates @ Cases[
                a["Instructions"],
                {op_ /; MemberQ[$circuitTwoQubitOps, op], x_, y_} /;
                    (x === c && y <= n) || (y === c && x <= n) :> If[x === c, y, x]
            ] <= 1
        ]
    ]
]

(* The number this gadget exists to make small: one over every single fault, against
   three for the bare-ancilla extraction of the same weight. *)
(* The instruction indices belonging to a spliced EC gadget's ancilla preparation,
   one copy per repetition.  Each repetition is (cat, cat measurement, EC), so the
   EC block of repetition r starts at a fixed stride, and the EC gadget names which
   of its own indices are preparation. *)
pauliMeasureSkipIndices[a_Association] := If[
    a["ErrorCorrection"] === None,
    {},
    Module[{ec, prep, block, stride},
        ec = QECErrorCorrection[QECCode[a["Code"]], a["CheckQubit"]];
        prep = steaneRegions[a["Code"], ec["Order"]]["Preparation"];
        block = ec["InstructionCount"];
        stride = Length[a["Instructions"]] / a["Repetitions"];
        Catenate @ Table[
            (r - 1) stride + (stride - block) + prep,
            {r, a["Repetitions"]}
        ]
    ]
]

QECPauliMeasurement[a_Association]["DataWeights"] := Sort @ pauliMeasureDataWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"], a["Pauli"],
    pauliMeasureSkipIndices[a]]

QECPauliMeasurement[a_Association]["MaxDataWeight"] := Max @ pauliMeasureDataWeights[
    a["Instructions"], a["Qubits"], a["DataQubits"], a["Pauli"],
    pauliMeasureSkipIndices[a]]

(* Theorem 12.1 needs both halves: enough repetitions for the code's distance, and
   an FTEC sub-gadget between them.  The second is not built yet, so this is False
   and will stay False until Steane EC can be spliced in. *)
QECPauliMeasurement[a_Association]["MeasurementCorrectQ"] :=
    a["ErrorCorrection"] =!= None &&
    a["Repetitions"] >= pauliMeasureRepetitions[QECCode[a["Code"]]["Distance"]]

(* What the structure still rests on.  The EC sub-gadget carries its own, and the
   measurement gadget inherits them rather than hiding them behind a True above. *)
QECPauliMeasurement[a_Association]["OpenAssumptions"] := If[
    a["ErrorCorrection"] === None,
    {"No error correction is spliced between repetitions, so a data error anticommuting with the measured Pauli defeats every repetition alike (Got26 sec. 12.1.4). The gadget handles ancilla faults only."},
    QECErrorCorrection[QECCode[a["Code"]], a["CheckQubit"]]["OpenAssumptions"]
]

QECPauliMeasurement[a_Association][prop_String] :=
    (Message[QECPauliMeasurement::noprop, prop]; Missing["NotFound", prop])


(* ---- the extraction a detector model can actually consume ---- *)

(* THE CONSTRAINT THAT SHAPES THIS, and it is worth stating plainly because it
   separates two things that look like they should compose and do not.

   A detector error model is a MATRIX: the effect of a set of faults is the XOR of
   their rows, and everything downstream -- the exact fold, the decoder table, the
   Stim export -- rests on that linearity over GF(2).

   Theorem 12.1's gadget takes a MAJORITY over 2t+1 repetitions.  Majority is not
   linear, so it cannot sit inside a detector model at all: two faults whose
   individual effects are known do not have an effect given by their XOR once a
   vote is in the path.  This is not an implementation gap; it is the same fork the
   package already stands on one side of.  Got26 sec. 12.2.2 answers a noisy
   syndrome by repeat-and-agree and sets differencing aside; DetectorModel.wl takes
   differencing.  QECPauliMeasurement is the book's gadget and keeps the vote.

   What the rate machinery takes instead is the half of the construction that fixes
   the thing majority voting never fixed: TRANSVERSALITY.  The hook error is cured
   by one cat qubit per letter, not by repeating; and a cat readout's syndrome bit
   is the PARITY of its m bits (eqs. 12.1-12.3), which is linear.  Repetition is
   then handled where it already was -- across rounds, by the detectors.  That is
   also what the surface-code practice of sec. 12.5.1 does, and what Stim can
   express.

   So "Transversal" is not a weaker "FaultTolerant"; it is the composition that the
   detector formalism admits, and the one whose rate can be computed. *)

$codeExtractions = {"BareAncilla", "Transversal"};

(* Instructions, qubit count, and how many record bits each generator contributes
   per round.  The block sizes are what lets recordSyndromes reduce a record to a
   syndrome without knowing which extraction produced it. *)
codeExtraction[a_Association, rounds_Integer, "BareAncilla"] := <|
    "Instructions" -> codeCircuitInstructions[a, rounds],
    "Qubits" -> a["Qubits"] + codeStabilizerCount[a],
    "BlockSizes" -> ConstantArray[1, codeStabilizerCount[a]]
|>

codeExtraction[a_Association, rounds_Integer, "Transversal"] :=
    codeExtraction[a, rounds, "Transversal"] = Module[
        {n = a["Qubits"], mat = a["CheckMatrix"], supports, widest, cat, check},
        supports = Table[
            Select[Range[n], letterAt[v, n, #] =!= {0, 0} &],
            {v, mat}
        ];
        widest = Max[Length /@ supports];
        (* One cat block and one check qubit, reset and reused by every generator of
           every round, exactly as the bare ancillas are.  A generator of weight w
           takes the first w of them. *)
        cat = Range[n + 1, n + widest];
        check = n + widest + 1;
        <|
            "Instructions" -> Catenate @ Table[
                Catenate @ Table[
                    With[{w = Length[supports[[j]]], sub = Take[cat, Length[supports[[j]]]]},
                        Join[
                            catInstructions[sub, check, catCheckPairs[w], 1],
                            catMeasureInstructions[mat[[j]], n, sub]
                        ]
                    ],
                    {j, Length[mat]}
                ],
                {rounds}
            ],
            "Qubits" -> check,
            "BlockSizes" -> Length /@ supports
        |>
    ]

(* A raw measurement record, reduced to one syndrome bit per generator per round.
   With one bit per block this is the identity, so the bare-ancilla path is
   untouched; with w bits it is their parity, which is the eigenvalue the Hadamard
   transform encodes and is linear, which is the whole point. *)
recordSyndromes[record_List, blocks_List, rounds_Integer] := With[
    {edges = Prepend[Accumulate[blocks], 0], per = Total[blocks]},
    Catenate @ Table[
        Table[
            Mod[Total[record[[(r - 1) per + edges[[j]] + 1 ;; (r - 1) per + edges[[j + 1]]]]], 2],
            {j, Length[blocks]}
        ],
        {r, rounds}
    ]
]


(* ---- formatting ---- *)

QECPauliMeasurement /: MakeBoxes[
    obj : QECPauliMeasurement[a_Association] /; KeyExistsQ[a, "Support"],
    form : (StandardForm | TraditionalForm)
] := BoxForm`ArrangeSummaryBox[
    QECPauliMeasurement,
    obj,
    BarChart[Values[Counts[First /@ a["Instructions"]]],
        ChartLabels -> Keys[Counts[First /@ a["Instructions"]]],
        ImageSize -> {Automatic, 34}, Axes -> False,
        ChartStyle -> RGBColor[0.15, 0.5, 0.65]],
    {
        BoxForm`SummaryItem[{"Pauli: ", QECPauliString[Append[a["Pauli"], 0]]}],
        BoxForm`SummaryItem[{"Weight: ", Length[a["Support"]]}],
        BoxForm`SummaryItem[{"Repetitions: ", a["Repetitions"]}]
    },
    {
        BoxForm`SummaryItem[{"Qubits: ", a["Qubits"]}],
        BoxForm`SummaryItem[{"Cat qubits: ", a["CatQubits"]}],
        BoxForm`SummaryItem[{"Error correction: ", a["ErrorCorrection"]}],
        BoxForm`SummaryItem[{"Instructions: ", Length[a["Instructions"]]}]
    },
    form,
    "Interpretable" -> False
]
