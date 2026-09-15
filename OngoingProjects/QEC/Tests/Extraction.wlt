(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Extraction.wlt

   Handing the rate machinery a different syndrome extraction.

   The constraint that shapes this file: a detector error model is a MATRIX, so
   the effect of a set of faults must be the XOR of their rows.  Theorem 12.1's
   gadget takes a majority over repetitions, and majority is not linear, so it
   cannot sit inside a detector model at all.  What can is the other half of the
   construction -- transversality -- whose readout is a PARITY, which is linear.
   QEC-Extraction-is-linear-over-GF2 is the test that says so, and it is what
   licenses everything else here.

   The load-bearing result is QEC-Extraction-ambiguous-signatures-collapse: with
   the bare ancilla, 288 circuit faults collapse onto 71 detector signatures of
   which 29 carry conflicting logical effects; with a transversal extraction,
   624 faults -- more than twice as many locations -- collapse onto 72 signatures
   of which only 4 do.  That is the hook error priced at the level of the rate
   rather than the gadget.

   And the four that remain are all cat-PREPARATION faults, which is the same
   open assumption Steane EC reports: verified ancillas are chapter 13 work.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecExtraction = Symbol[qecScope <> "codeExtraction"];
qecSyndromes  = Symbol[qecScope <> "recordSyndromes"];
qecFaultEffect = Symbol[qecScope <> "faultEffect"];
qecCodeData   = Symbol[qecScope <> "codeData"];
qecPaulis     = Symbol[qecScope <> "$oneQubitPaulis"];
qecModes      = Symbol[qecScope <> "$codeExtractions"];

bitFlip = QECCode["BitFlipCode"];
five = QECCode["5QubitCode"];

qecAmbiguity[dem_] := Module[{d, o, g},
    d = dem["DetectorMatrix"]; o = dem["ObservableMatrix"];
    g = GroupBy[Range[Length[d]], d[[#]] &, DeleteDuplicates[o[[#]] & /@ #] &];
    {Length[d], Length[g], Count[g, alt_ /; Length[alt] > 1]}
]


(* ============================================================================
   Reducing a record to a syndrome
   ============================================================================ *)

(* One record bit per check is the identity, which is what keeps the bare-ancilla
   path byte for byte what it was. *)
VerificationTest[
    qecSyndromes[{1, 0, 1, 1}, {1, 1}, 2],
    {1, 0, 1, 1},
    TestID -> "QEC-Extraction-one-bit-per-check-is-the-identity"
]

(* Several bits per check reduce by parity, which is the eigenvalue the Hadamard
   transform of eqs. 12.1-12.3 encodes. *)
VerificationTest[
    qecSyndromes[{1, 1, 0, 1, 0, 0, 1, 0}, {2, 2}, 2],
    {0, 1, 0, 1},
    TestID -> "QEC-Extraction-many-bits-per-check-reduce-by-parity"
]

(* Checks of different weights, which is the general case: the five-qubit code's
   generators all have weight four, Steane's mix. *)
VerificationTest[
    qecSyndromes[{1, 1, 1, 0, 0, 1}, {3, 1, 2}, 1],
    {1, 0, 1},
    TestID -> "QEC-Extraction-blocks-may-differ-in-width"
]


(* ============================================================================
   Shape
   ============================================================================ *)

VerificationTest[
    Sort[qecModes],
    {"BareAncilla", "Transversal"},
    TestID -> "QEC-Extraction-two-extractions-are-offered"
]

(* The bare-ancilla extraction is the circuit that was already there, described in
   the new vocabulary: one ancilla per check, one record bit per check. *)
VerificationTest[
    With[{ex = qecExtraction[qecCodeData[bitFlip], 1, "BareAncilla"]},
        {ex["Qubits"], ex["BlockSizes"], ex["Instructions"] === bitFlip["SyndromeCircuit"]["Instructions"]}
    ],
    {5, {1, 1}, True},
    TestID -> "QEC-Extraction-bare-ancilla-is-the-old-circuit"
]

(* The transversal one measures each check through a cat, so its block sizes are
   the generator weights and it needs a cat block plus one check qubit. *)
VerificationTest[
    With[{ex = qecExtraction[qecCodeData[bitFlip], 1, "Transversal"]},
        {ex["Qubits"], ex["BlockSizes"], Count[ex["Instructions"], {"CZ", _, _}]}
    ],
    {6, {2, 2}, 4},
    TestID -> "QEC-Extraction-transversal-measures-each-check-through-a-cat"
]

(* One cat block reused by every check of every round, as the bare ancillas are:
   the five-qubit code's four weight-four checks share four cat qubits, not
   sixteen. *)
VerificationTest[
    With[{ex = qecExtraction[qecCodeData[five], 1, "Transversal"]},
        {ex["Qubits"], ex["BlockSizes"]}
    ],
    {10, {4, 4, 4, 4}},
    TestID -> "QEC-Extraction-the-cat-block-is-reused"
]


(* ============================================================================
   The property that makes it usable at all
   ============================================================================ *)

(* A detector model is a matrix, so the effect of two faults has to be the XOR of
   their rows.  Parity readout keeps that; a majority vote would not, which is why
   Theorem 12.1's gadget cannot be the thing the rate machinery consumes. *)
VerificationTest[
    Module[{a, ex, instr, nq, blocks, eff, faults, pairs},
        a = qecCodeData[bitFlip];
        ex = qecExtraction[a, 1, "Transversal"];
        instr = ex["Instructions"]; nq = ex["Qubits"]; blocks = ex["BlockSizes"];
        eff[f_] := qecFaultEffect[a, instr, nq, 2, 1, f, blocks];
        faults = Flatten[
            Table[{i, q, pauli}, {i, 0, Length[instr]}, {q, nq}, {pauli, qecPaulis}], 2];
        pairs = BlockRandom[SeedRandom[20260914]; RandomSample[Subsets[faults, {2}], 200]];
        AllTrue[pairs,
            With[{c = eff[#], u = eff[{#[[1]]}], v = eff[{#[[2]]}]},
                c[[1]] === BitXor[u[[1]], v[[1]]] &&
                c[[2]] === BitXor[u[[2]], v[[2]]] &&
                c[[3]] === BitXor[u[[3]], v[[3]]]] &]
    ],
    True,
    TestID -> "QEC-Extraction-is-linear-over-GF2"
]

(* And it measures the right thing: with no faults every detector is quiet, and a
   data error reproduces the syndrome codeSyndrome computes by matrix product. *)
VerificationTest[
    Module[{a, ex, instr, nq, blocks},
        a = qecCodeData[bitFlip];
        ex = qecExtraction[a, 1, "Transversal"];
        instr = ex["Instructions"]; nq = ex["Qubits"]; blocks = ex["BlockSizes"];
        {
            qecFaultEffect[a, instr, nq, 2, 1, {}, blocks][[1]],
            Take[qecFaultEffect[a, instr, nq, 2, 1, {{0, 1, {1, 0}}}, blocks][[1]], 2],
            bitFlip["Syndrome", "XII"]
        }
    ],
    {{0, 0, 0, 0}, {1, 0}, {1, 0}},
    TestID -> "QEC-Extraction-transversal-reproduces-the-syndrome"
]


(* ============================================================================
   No regression
   ============================================================================ *)

(* The default is the bare ancilla, so every number the package produced before
   this option existed is unchanged. *)
VerificationTest[
    {QECDetectorModel[bitFlip, QECNoiseModel["Circuit", 1/1000], 1]["Extraction"],
     QECDetectorModel[bitFlip, QECNoiseModel["Circuit", 1/1000], 1]["Faults"],
     Normal @ Series[
         QECLogicalErrorRate[bitFlip, QECNoiseModel["Circuit", p], "Rounds" -> 1],
         {p, 0, 1}]},
    {"BareAncilla", 60, 32 p / 15},
    TestID -> "QEC-Extraction-defaults-to-the-bare-ancilla"
]

(* Code capacity has no circuit and so no extraction to choose; the option is
   accepted and ignored rather than refused, so a sweep across levels need not
   strip it. *)
VerificationTest[
    Expand @ QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", p], "Extraction" -> "Transversal"],
    3 p^2 - 2 p^3,
    TestID -> "QEC-Extraction-code-capacity-ignores-it"
]


(* ============================================================================
   What it buys
   ============================================================================ *)

(* THE result.  Bare: 288 faults onto 71 signatures, 29 of them carrying
   conflicting logical effects -- which is the hook error, and why the rate starts
   at 29p/5 rather than p^2.  Transversal: 624 faults, more than twice the
   locations, onto 72 signatures of which only 4 conflict. *)
VerificationTest[
    {qecAmbiguity[QECDetectorModel[five, QECNoiseModel["Circuit", 1/1000], 1]],
     qecAmbiguity[QECDetectorModel[five, QECNoiseModel["Circuit", 1/1000], 1,
         "Extraction" -> "Transversal"]]},
    {{288, 71, 29}, {624, 72, 4}},
    TestID -> "QEC-Extraction-ambiguous-signatures-collapse"
]

(* And the four that survive are all rejected shots.  Every one of them is a fault
   in the cat's own preparation -- the reset of a cat qubit, a chain CNOT -- that
   leaves a weight-two X pattern, which is exactly what the verification check is
   there to catch.  Condition on the checks ACCEPTING, which is the only case the
   experiment keeps, and no single fault is ambiguous at all.

   This is the same conditional statement QECPauliMeasurement makes about its
   residual data weight, and it has to be made on both sides or not at all. *)
VerificationTest[
    Module[{dem, d, o, h, keep, g},
        dem = QECDetectorModel[five, QECNoiseModel["Circuit", 1/1000], 1,
            "Extraction" -> "Transversal"];
        d = dem["DetectorMatrix"]; o = dem["ObservableMatrix"]; h = dem["HeraldMatrix"];
        keep = Select[Range[Length[d]], Total[h[[#]]] === 0 &];
        g = GroupBy[keep, d[[#]] &, DeleteDuplicates[o[[#]] & /@ #] &];
        {Length[keep], Length[g], Count[g, alt_ /; Length[alt] > 1]}
    ],
    {484, 67, 0},
    TestID -> "QEC-Extraction-accepted-faults-are-never-ambiguous"
]

(* Which is why the decoder must not offer a rejected fault as an explanation.  It
   only ever runs on an accepted shot, so a fault that trips a verification check
   cannot have happened; leaving it in the hypothesis space lets the lightest-set
   rule claim a detector pattern on behalf of a fault the experiment threw away, and
   that claim displaces the real explanation.  The table shrinks to the accepted
   rows, and the bare-ancilla path -- which has no heralds -- is untouched. *)
VerificationTest[
    Module[{demT, demB},
        demT = First[QECDetectorModel[five, QECNoiseModel["Circuit", 1/1000], 1,
            "Extraction" -> "Transversal"]];
        demB = First[QECDetectorModel[five, QECNoiseModel["Circuit", 1/1000], 1]];
        {demT["Heralds"] > 0, demB["Heralds"] === 0}
    ],
    {True, True},
    TestID -> "QEC-Extraction-only-the-transversal-path-post-selects"
]


(* ============================================================================
   The milestone
   ============================================================================ *)

(* The whole point of chapter 12, measured in this layer's own units.  The
   five-qubit code has distance three, so at code capacity and phenomenologically
   its rate starts at p^2: any one error is corrected.  With the bare-ancilla
   extraction it starts at p instead -- that is the hook error, and it is a defect
   of the gadget rather than of circuit-level noise.  With a transversal extraction
   through verified cat states, the exponent comes back.

   Measured rather than asserted: the rate is evaluated at two rates a factor of two
   apart and the exponent read off.  Two is what a distance-three code owes. *)
VerificationTest[
    Module[{r1, r2, exponent},
        r1 = QECLogicalErrorRate[five, QECNoiseModel["Circuit", 1/1000], "Rounds" -> 1,
            "Extraction" -> "Transversal"]["Rate"];
        r2 = QECLogicalErrorRate[five, QECNoiseModel["Circuit", 1/2000], "Rounds" -> 1,
            "Extraction" -> "Transversal"]["Rate"];
        exponent = N[Log[r1 / r2] / Log[2]];
        1.9 < exponent < 2.1
    ],
    True,
    TestID -> "QEC-Extraction-restores-the-p-squared"
]

(* Against the same code with the bare ancilla, whose exponent is one. *)
VerificationTest[
    Module[{exponent},
        exponent = N @ Log[
            QECLogicalErrorRate[five, QECNoiseModel["Circuit", 1/1000], "Rounds" -> 1] /
            QECLogicalErrorRate[five, QECNoiseModel["Circuit", 1/2000], "Rounds" -> 1]
        ] / Log[2];
        0.9 < exponent < 1.1
    ],
    True,
    TestID -> "QEC-Extraction-the-bare-ancilla-exponent-is-one"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

VerificationTest[
    QECDetectorModel[bitFlip, QECNoiseModel["Circuit", 1/1000], 1, "Extraction" -> "Shor"],
    $Failed,
    {QECDetectorModel::extraction},
    TestID -> "QEC-Extraction-an-unknown-extraction-is-refused"
]

VerificationTest[
    QECLogicalErrorRate[bitFlip, QECNoiseModel["Circuit", 1/1000], "Rounds" -> 1,
        "Extraction" -> "Knill"],
    $Failed,
    {QECDetectorModel::extraction},
    TestID -> "QEC-Extraction-the-rate-refuses-it-too"
]
