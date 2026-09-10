(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Cat.wlt

   Cat states and their verification (book sec. 12.1.2 - 12.1.3).

   The point of this file is that the check set is DERIVED, not assumed.  The book
   says "if we do this on enough pairs of qubits, any single fault ... will be picked
   up by the checks" and leaves the set open; catXPatterns propagates every single
   fault of the preparation through the frame propagator and catChecksCoverQ is the
   proof obligation.  The tests below pin the derivation against the book's own
   worked error (figure 12.3b) and then assert coverage.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecCatPrep   = Symbol[qecScope <> "catPrepInstructions"];
qecCatCheck  = Symbol[qecScope <> "catCheckInstructions"];
qecCatX      = Symbol[qecScope <> "catXPatterns"];
qecCatPairs  = Symbol[qecScope <> "catCheckPairs"];
qecCatCover  = Symbol[qecScope <> "catChecksCoverQ"];
qecCatCanon  = Symbol[qecScope <> "catCanonicalPattern"];
qecFrame     = Symbol[qecScope <> "framePropagate"];


(* ============================================================================
   Preparation
   ============================================================================ *)

(* Book figure 12.3a: one qubit in |+>, the rest in |0>, then CNOTs.  The chain, so
   that the book's own error example belongs to this circuit. *)
VerificationTest[
    qecCatPrep[{1, 2, 3, 4}],
    {{"R", 1}, {"R", 2}, {"R", 3}, {"R", 4}, {"H", 1},
     {"CNOT", 1, 2}, {"CNOT", 2, 3}, {"CNOT", 3, 4}},
    TestID -> "QEC-Cat-prep-is-the-chain"
]

(* It can be placed anywhere, since a cat state lives inside a bigger circuit. *)
VerificationTest[
    qecCatPrep[{5, 6, 7}],
    {{"R", 5}, {"R", 6}, {"R", 7}, {"H", 5}, {"CNOT", 5, 6}, {"CNOT", 6, 7}},
    TestID -> "QEC-Cat-prep-on-given-qubits"
]


(* ============================================================================
   The canonical form
   ============================================================================ *)

(* X^(tensor m) stabilises the cat state, so a pattern and its complement are the
   same error and must reduce to the same representative. *)
VerificationTest[
    qecCatCanon[{0, 0, 1, 1}] === qecCatCanon[{1, 1, 0, 0}],
    True,
    TestID -> "QEC-Cat-pattern-and-complement-agree"
]

VerificationTest[
    {qecCatCanon[{1, 1, 1, 1}], qecCatCanon[{1, 0, 1, 1}]},
    {{0, 0, 0, 0}, {0, 1, 0, 0}},
    TestID -> "QEC-Cat-canonical-takes-the-lighter-half"
]


(* ============================================================================
   The derived patterns
   ============================================================================ *)

(* THE BOOK'S OWN EXAMPLE.  Figure 12.3b shows a single fault leaving |0011> + |1100>
   instead of |0000> + |1111>, which is the X pattern {0,0,1,1}.  The derivation finds
   exactly that, and nothing else, for a four-qubit cat. *)
VerificationTest[
    qecCatX[{1, 2, 3, 4}, 5],
    {{0, 0, 1, 1}},
    TestID -> "QEC-Cat-reproduces-figure-12-3b"
]

(* Two and three qubit cats have no dangerous pattern at all: every single fault
   leaves something of canonical weight at most one, which is a single data error and
   the code corrects it.  So they need no checking, and the default set is empty. *)
VerificationTest[
    {qecCatX[{1, 2}, 3], qecCatX[{1, 2, 3}, 4], qecCatPairs[2], qecCatPairs[3]},
    {{}, {}, {}, {}},
    TestID -> "QEC-Cat-small-cats-need-no-checks"
]

(* A chain fault reaching qubit j leaves a suffix of length L = m - j + 1, canonical
   weight Min[L, m - L], so the dangerous L run from 2 to m - 2: there are m - 3. *)
VerificationTest[
    Table[Length[qecCatX[Range[m], m + 1]], {m, 4, 8}],
    {1, 2, 3, 4, 5},
    TestID -> "QEC-Cat-dangerous-count-is-m-minus-three"
]

(* Every dangerous pattern has canonical weight at least two, by construction: weight
   one is a single data error, which is the whole thing the code is for. *)
VerificationTest[
    AllTrue[Catenate[Table[qecCatX[Range[m], m + 1], {m, 4, 8}]], Total[#] >= 2 &],
    True,
    TestID -> "QEC-Cat-dangerous-means-weight-two-or-more"
]

(* Lowering the threshold does find the weight-one patterns, so their absence above is
   the filter doing its job and not the derivation missing them. *)
VerificationTest[
    Length[qecCatX[{1, 2, 3, 4}, 5, "Weight" -> 1]] > Length[qecCatX[{1, 2, 3, 4}, 5]],
    True,
    TestID -> "QEC-Cat-weight-one-patterns-exist-but-are-excluded"
]


(* ============================================================================
   The checks cover them
   ============================================================================ *)

(* THE PROOF OBLIGATION.  This is the test the whole file exists for: the default
   pairs detect every dangerous pattern, for every size.  A pair detects a pattern
   when the pattern has odd weight on it. *)
VerificationTest[
    AllTrue[Range[2, 9], qecCatCover[qecCatPairs[#], qecCatX[Range[#], # + 1]] &],
    True,
    TestID -> "QEC-Cat-default-pairs-cover-every-size"
]

(* And they are not more than needed: each dangerous suffix is caught by exactly one
   chain pair, so dropping any one pair loses coverage. *)
VerificationTest[
    With[{m = 6, pairs = qecCatPairs[6], pat = qecCatX[Range[6], 7]},
        AllTrue[pairs, ! qecCatCover[DeleteCases[pairs, #], pat] &]
    ],
    True,
    TestID -> "QEC-Cat-no-pair-is-redundant"
]

(* A wrong pair set is reported as not covering, rather than passing quietly. *)
VerificationTest[
    qecCatCover[{{1, 2}}, qecCatX[{1, 2, 3, 4}, 5]],
    False,
    TestID -> "QEC-Cat-wrong-pairs-do-not-cover"
]


(* ============================================================================
   The check circuit is itself safe
   ============================================================================ *)

(* The cat qubits are the CONTROLS and the check qubit the target, which is what makes
   the non-transversal check circuit harmless: X travels control to target only. *)
VerificationTest[
    qecCatCheck[{1, 2, 3, 4}, 5, {2, 3}],
    {{"R", 5}, {"CNOT", 2, 5}, {"CNOT", 3, 5}, {"MH", 5}},
    TestID -> "QEC-Cat-check-has-cat-as-control"
]

(* So an X on the check qubit reaches no cat qubit at all -- it flips the herald and
   nothing else, exactly as a readout error does in the bare-ancilla circuit. *)
VerificationTest[
    With[{instr = qecCatCheck[{1, 2, 3, 4}, 5, {2, 3}]},
        With[{run = qecFrame[instr, 5, {{1, 5, {1, 0}}}]},
            {Take[run["Frame"][[1]], 4], run["Heralds"]}
        ]
    ],
    {{0, 0, 0, 0}, {1}},
    TestID -> "QEC-Cat-check-qubit-X-cannot-reach-the-cat"
]

(* The outcome is a herald, not a syndrome bit: it post-selects the attempt and never
   enters the record the decoder reads. *)
VerificationTest[
    With[{instr = qecCatCheck[{1, 2, 3, 4}, 5, {2, 3}]},
        qecFrame[instr, 5, {{1, 5, {1, 0}}}]["Record"]
    ],
    {},
    TestID -> "QEC-Cat-checks-are-heralds-not-syndrome"
]

(* And the check really measures Z_i Z_j: an X on one of the pair flips it, an X on
   both does not, an X outside the pair does not. *)
VerificationTest[
    With[{instr = qecCatCheck[{1, 2, 3, 4}, 5, {2, 3}]},
        First /@ {
            qecFrame[instr, 5, {{0, 2, {1, 0}}}]["Heralds"],
            qecFrame[instr, 5, {{0, 3, {1, 0}}}]["Heralds"],
            qecFrame[instr, 5, {{0, 2, {1, 0}}, {0, 3, {1, 0}}}]["Heralds"],
            qecFrame[instr, 5, {{0, 1, {1, 0}}}]["Heralds"]
        }
    ],
    {1, 1, 0, 0},
    TestID -> "QEC-Cat-check-measures-the-pair-parity"
]


(* ============================================================================
   The object
   ============================================================================ *)

VerificationTest[
    QECCatState[4]["ChecksCoverQ"],
    True,
    TestID -> "QEC-Cat-object-checks-cover"
]

VerificationTest[
    With[{cat = QECCatState[4]},
        {cat["Size"], cat["CatQubits"], cat["CheckQubit"], cat["Pairs"], cat["Heralds"]}
    ],
    {4, {1, 2, 3, 4}, 5, {{2, 3}}, 1},
    TestID -> "QEC-Cat-object-shape"
]

(* Repetitions guard against a faulty check qubit tricking us into keeping a bad cat,
   so they multiply the heralds and nothing else. *)
VerificationTest[
    With[{a = QECCatState[6], b = QECCatState[6, "Repetitions" -> 3]},
        {a["Heralds"], b["Heralds"], b["Pairs"] === a["Pairs"]}
    ],
    {3, 9, True},
    TestID -> "QEC-Cat-repetitions-multiply-the-heralds"
]

(* A hand-given pair set is honoured, and reported as not covering if it does not. *)
VerificationTest[
    With[{cat = QECCatState[4, "Pairs" -> {{1, 2}}]},
        {cat["Pairs"], cat["ChecksCoverQ"]}
    ],
    {{{1, 2}}, False},
    TestID -> "QEC-Cat-hand-given-pairs-are-checked-honestly"
]

VerificationTest[
    QECCatState[1],
    $Failed,
    {QECCatState::size},
    TestID -> "QEC-Cat-too-small-refused"
]

VerificationTest[
    QECCatState[4, "Repetitions" -> 0],
    $Failed,
    {QECCatState::reps},
    TestID -> "QEC-Cat-bad-repetitions-refused"
]

VerificationTest[
    QECCatState[4, "Pairs" -> {{1, 9}}],
    $Failed,
    {QECCatState::pairs},
    TestID -> "QEC-Cat-out-of-range-pairs-refused"
]

VerificationTest[
    QECCatState[4]["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECCatState::noprop},
    TestID -> "QEC-Cat-bad-property"
]
