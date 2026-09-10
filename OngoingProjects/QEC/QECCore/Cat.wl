(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECCatState]

PackageScope[catPrepInstructions]
PackageScope[catCheckInstructions]
PackageScope[catInstructions]
PackageScope[catXPatterns]
PackageScope[catCheckPairs]
PackageScope[catChecksCoverQ]
PackageScope[catCanonicalPattern]
PackageScope[catData]


(* ============================================================================ *)
(* Cat states, and checking them.                                               *)
(*                                                                              *)
(* The bare-ancilla extraction of Circuit.wl is not fault tolerant, and the      *)
(* diagnosis in book sec. 12.1.2 is precise: one ancilla shared by every data    *)
(* qubit of a check lets a single fault reach several of them.  The fix is to    *)
(* make the controlled-P transversal -- "one ancilla qubit for each single-qubit *)
(* Pauli gate making up P, and each ancilla qubit interacts with only one qubit  *)
(* in the data block.  That limits error propagation, just as for transversal    *)
(* gate gadgets."                                                                *)
(*                                                                              *)
(* For that the m = wt(P) ancillas have to be in a coherent superposition of     *)
(* "do nothing" and "apply P to all of them", which is                           *)
(*                                                                              *)
(*     |0...0> + |1...1>                                                          *)
(*                                                                              *)
(* the GHZ state, called a cat state in this literature.  This file builds one   *)
(* and checks it.  The measurement gadget that consumes it is the next file.     *)
(*                                                                              *)
(* WHAT GOES WRONG, and it is worth being exact because the two failure modes     *)
(* have two different fixes.                                                      *)
(*                                                                              *)
(* Z errors are not the problem.  On a cat state a single Z gives |0...0> -       *)
(* |1...1>, and two Z errors on different qubits give the correct state back, so   *)
(* "an odd number of Z errors gives us the same state as just one Z error and an   *)
(* even number of Z errors is the same as no Z errors.  That is, for the cat       *)
(* state, it's only possible to have a single Z error!"  One Z does flip the       *)
(* measured eigenvalue, but that is fixed by REPEATING the measurement (sec.       *)
(* 12.1.4), not by checking the state.                                             *)
(*                                                                              *)
(* X errors are the problem, and only when there are several of them: multiple X  *)
(* errors on the cat become multiple errors in the data block, from one fault.    *)
(* That is what checking is for.  Y is X times Z and needs no separate treatment. *)
(*                                                                              *)
(* WHICH PAIRS TO CHECK.  The book says "if we do this on enough pairs of qubits, *)
(* any single fault in the circuit originally constructing the cat state will be  *)
(* picked up by the checks" and leaves the set open.  Here it is derived instead: *)
(* catXPatterns propagates every single fault of the preparation circuit through  *)
(* the frame propagator and collects the X patterns it can leave on the cat, and  *)
(* catChecksCoverQ is the proof obligation that the chosen pairs detect all the    *)
(* dangerous ones.  Nothing is assumed by analogy.                                *)
(*                                                                              *)
(* Two reductions make that set small.                                            *)
(*                                                                              *)
(*   - X^(tensor m) is a stabilizer of the cat state, so an X pattern and its     *)
(*     complement are the same error.  Patterns are canonicalised accordingly.    *)
(*   - a pattern of canonical weight one is a single data error, which the code   *)
(*     corrects.  Only weight two and above is dangerous.                          *)
(*                                                                              *)
(* For the chain preparation below, a fault reaching qubit j propagates along the *)
(* rest of the chain and leaves a suffix of 1s of length L = m - j + 1, which     *)
(* canonicalises to weight Min[L, m - L].  So the dangerous L are 2 <= L <= m-2,  *)
(* there are m-3 of them, and each is caught by exactly one chain pair -- which   *)
(* is why the default check set is the chain pairs from 2 to m-2 and not all m-1. *)
(* An m = 4 cat needs one check, not three.                                        *)
(*                                                                              *)
(* WHY THE CHECK CIRCUIT IS ITSELF SAFE, since it is not transversal either.  The *)
(* cat qubits are the CONTROLS and the check qubit is the target, so an X can go   *)
(* from a cat qubit into the check qubit but not back out into a second one -- X   *)
(* only travels control to target.  A fault in the check qubit's preparation or in *)
(* the first CNOT can only put a phase error into the cat.  The worst case is a    *)
(* fault in the first CNOT leaving a bit flip on one cat qubit and a phase on the  *)
(* check qubit that becomes a phase on a second cat qubit -- two errors from one   *)
(* fault -- except that Z errors on a cat state are degenerate, so X (tensor) Z    *)
(* rewrites as -i Y (tensor) I and it is one error after all.                      *)
(*                                                                              *)
(* The check outcomes are heralds, not syndrome bits: a nonzero one means discard  *)
(* the attempt.  That is why the instruction language has "MH".  A wrong check     *)
(* outcome does not damage the cat but can make us keep a bad one, so the checks   *)
(* are repeated for codes correcting more than one error -- the "Repetitions"      *)
(* option.                                                                         *)
(* ============================================================================ *)

QECCatState::usage = "QECCatState[m] gives the m-qubit cat state gadget: a preparation circuit followed by the parity checks that verify it.\nQECCatState[qubits, check] places it on given qubits, with check the qubit reused for every parity check.\nThe option \"Pairs\" sets which pairs are checked and \"Repetitions\" how many times.\ncat[prop] gives a property; cat[\"Properties\"] lists them.";

QECCatState::size = "A cat state needs at least two qubits; got `1`.";
QECCatState::pairs = "\"Pairs\" must be a list of pairs of positions in 1..`1`; got `2`.";
QECCatState::reps = "\"Repetitions\" must be a positive integer; got `1`.";
QECCatState::noprop = "`1` is not a property of QECCatState. Use cat[\"Properties\"] for the list.";


(* ---- preparation ---- *)

(* Book figure 12.3a: one qubit in |+>, the rest in |0>, then CNOTs to entangle.
   The chain is the arrangement the book's own error example belongs to -- figure
   12.3b's |0011> + |1100> is a suffix of two, which is what a fault mid-chain
   leaves.  A star (all CNOTs from the first qubit) is equally valid and gives the
   same set of canonical patterns; the chain is used so the figure matches. *)
catPrepInstructions[qubits_List] := Join[
    {{"R", #}} & /@ qubits // Catenate,
    {{"H", First[qubits]}},
    Table[{"CNOT", qubits[[i]], qubits[[i + 1]]}, {i, Length[qubits] - 1}]
]


(* ---- the dangerous X patterns, derived ---- *)

(* X^(tensor m) stabilises the cat state, so v and its complement are one error. *)
catCanonicalPattern[v_List] := With[{w = 1 - v},
    Which[
        Total[v] < Total[w], v,
        Total[w] < Total[v], w,
        OrderedQ[{v, w}], v,
        True, w
    ]
]

(* Every X pattern a single fault in the preparation can leave on the cat, reduced
   and filtered to the ones that matter.  "Weight" -> k keeps patterns of canonical
   weight at least k; the default 2 is the dangerous set, since weight one is a
   single data error and the code corrects that. *)
Options[catXPatterns] = {"Weight" -> 2};

catXPatterns[qubits_List, nq_Integer, OptionsPattern[]] := Module[
    {instr = catPrepInstructions[qubits], w = OptionValue["Weight"], raw},
    (* three iterators, so Flatten at level 2 -- Catenate would leave it two deep and
       the mapping below would silently run on lists of patterns instead of patterns *)
    raw = Flatten[
        Table[
            Part[framePropagate[instr, nq, {{i, q, pauli}}]["Frame"], 1][[qubits]],
            {i, 0, Length[instr]}, {q, qubits}, {pauli, $oneQubitPaulis}
        ],
        2
    ];
    Select[DeleteDuplicates[catCanonicalPattern /@ raw], Total[#] >= w &]
]


(* ---- the checks ---- *)

(* A check on the pair {i, j} of cat positions measures Z_i Z_j: the cat qubits are
   the controls, the check qubit the target, and the outcome is a herald. *)
catCheckInstructions[qubits_List, check_Integer, {i_Integer, j_Integer}] := {
    {"R", check},
    {"CNOT", qubits[[i]], check},
    {"CNOT", qubits[[j]], check},
    {"MH", check}
}

(* The chain pairs that suffice for the chain preparation, derived above: a suffix of
   length L is caught only by the pair (m-L, m-L+1), and the dangerous L run from 2
   to m-2.  For m = 2 and m = 3 no single fault leaves a dangerous pattern at all, so
   the set is empty and nothing needs checking. *)
catCheckPairs[m_Integer] := Table[{k, k + 1}, {k, 2, m - 2}]

(* A pair detects a pattern when the pattern has odd weight on it.  This is the proof
   obligation: the chosen pairs must detect every dangerous pattern. *)
catChecksCoverQ[pairs_List, patterns_List] := AllTrue[
    patterns,
    Function[v, AnyTrue[pairs, OddQ[v[[#[[1]]]] + v[[#[[2]]]]] &]]
]


(* ---- the whole gadget ---- *)

(* Flatten at level 2, not Catenate: two iterators over a body that is itself a list
   of instructions, so Catenate would leave the checks nested one level down and the
   instruction list silently malformed.  Third time this shape has bitten this
   package -- see phenomenologicalMechanisms and idleMechanisms. *)
catInstructions[qubits_List, check_Integer, pairs_List, reps_Integer] := Join[
    catPrepInstructions[qubits],
    Flatten[
        Table[
            catCheckInstructions[qubits, check, pair],
            {r, reps}, {pair, pairs}
        ],
        2
    ]
]


(* ---- construction ---- *)

Options[QECCatState] = {"Pairs" -> Automatic, "Repetitions" -> 1};

QECCatState[m_Integer, opts : OptionsPattern[]] :=
    QECCatState[Range[m], m + 1, opts]

QECCatState[qubits_List, check_Integer, opts : OptionsPattern[]] := Module[
    {m = Length[qubits], pairs, reps, nq},
    reps = OptionValue["Repetitions"];
    pairs = Replace[OptionValue["Pairs"], Automatic :> catCheckPairs[m]];
    nq = Max[Append[qubits, check]];
    Which[
        m < 2,
            Message[QECCatState::size, m]; $Failed,
        ! (IntegerQ[reps] && reps > 0),
            Message[QECCatState::reps, reps]; $Failed,
        ! MatchQ[pairs, {{_Integer, _Integer} ...}] ||
            ! AllTrue[Catenate[pairs], 1 <= # <= m &],
            Message[QECCatState::pairs, m, pairs]; $Failed,
        True,
            QECCatState[<|
                "CatQubits" -> qubits,
                "CheckQubit" -> check,
                "Pairs" -> pairs,
                "Repetitions" -> reps,
                "Qubits" -> nq,
                "Instructions" -> catInstructions[qubits, check, pairs, reps]
            |>]
    ]
]


(* ---- properties ---- *)

$catProperties = {
    "Instructions", "CatQubits", "CheckQubit", "Pairs", "Repetitions", "Qubits",
    "Size", "PrepInstructions", "Heralds", "DangerousPatterns", "ChecksCoverQ",
    "GateCounts", "Properties"
};

catData[QECCatState[a_Association]] := a

QECCatState[_Association]["Properties"] := $catProperties

QECCatState[a_Association][prop : ("Instructions" | "CatQubits" | "CheckQubit" |
    "Pairs" | "Repetitions" | "Qubits")] := a[prop]

QECCatState[a_Association]["Size"] := Length[a["CatQubits"]]
QECCatState[a_Association]["PrepInstructions"] := catPrepInstructions[a["CatQubits"]]
QECCatState[a_Association]["Heralds"] := Count[a["Instructions"], {"MH", _}]
QECCatState[a_Association]["GateCounts"] := Counts[First /@ a["Instructions"]]

(* The X patterns a single fault in the preparation can leave, and whether the checks
   catch them.  A False here means the gadget is not fault tolerant and says so. *)
QECCatState[a_Association]["DangerousPatterns"] := catXPatterns[a["CatQubits"], a["Qubits"]]
QECCatState[a_Association]["ChecksCoverQ"] :=
    catChecksCoverQ[a["Pairs"], catXPatterns[a["CatQubits"], a["Qubits"]]]

QECCatState[a_Association][prop_String] := (Message[QECCatState::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECCatState /: MakeBoxes[obj : QECCatState[a_Association] /; KeyExistsQ[a, "CatQubits"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECCatState,
        obj,
        BarChart[Values[Counts[First /@ a["Instructions"]]],
            ChartLabels -> Keys[Counts[First /@ a["Instructions"]]],
            ImageSize -> {Automatic, 34}, Axes -> False, ChartStyle -> RGBColor[0.15, 0.5, 0.65]],
        {
            BoxForm`SummaryItem[{"Size: ", Length[a["CatQubits"]]}],
            BoxForm`SummaryItem[{"Checks: ", Length[a["Pairs"]] a["Repetitions"]}],
            BoxForm`SummaryItem[{"Qubits: ", a["Qubits"]}]
        },
        {
            BoxForm`SummaryItem[{"Pairs: ", a["Pairs"]}],
            BoxForm`SummaryItem[{"Repetitions: ", a["Repetitions"]}],
            BoxForm`SummaryItem[{"Instructions: ", Length[a["Instructions"]]}]
        },
        form,
        "Interpretable" -> False
    ]
