(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECSyndromeCircuit]

PackageScope[circuitData]
PackageScope[circuitInstructions]
PackageScope[circuitQubitCount]
PackageScope[codeCircuitInstructions]
PackageScope[generatorInstructions]
PackageScope[framePropagate]
PackageScope[frameZero]
PackageScope[instructionEngineGates]
PackageScope[$circuitOneQubitOps]


(* ============================================================================ *)
(* The syndrome-extraction circuit.                                             *)
(*                                                                              *)
(* Up to here a check was measured by fiat: codeSyndrome multiplies the error by *)
(* the check matrix and reads off anticommutations.  That is the right object at *)
(* code-capacity level, where measurement is assumed perfect.  It is also the    *)
(* thing that has to go if the noise is going to be honest, because in a real    *)
(* device a check is not a matrix product -- it is a subcircuit, and every gate  *)
(* and every readout in it can fail.                                            *)
(*                                                                              *)
(* So: one ancilla per generator.  This is the construction of Gottesman, QECC     *)
(* book sec. 12.1.1, figure 12.1a -- and that section is titled "Non-Fault-        *)
(* Tolerant Measurement of Paulis", which is the honest name for what this is.  The *)
(* book's version puts the ancilla in |+> and applies controlled-P from ancilla to  *)
(* data; conjugating each data qubit by the rotation that sends its letter to Z      *)
(* turns that into a controlled-Z ladder, and since H.CZ.H = CNOT the ancilla's two *)
(* Hadamards are absorbed into the ladder.  Gate for gate it is the same circuit,    *)
(* with the basis change written out instead of folded into "controlled-P".  To      *)
(* measure a                                                                        *)
(* Pauli generator g, rotate each data qubit in g's support so that g's letter    *)
(* there becomes Z, CNOT each of those data qubits into the ancilla, rotate back, *)
(* and measure the ancilla in Z.  The ancilla ends up holding the parity of the   *)
(* rotated qubits, which is exactly the eigenvalue of g.                          *)
(*                                                                              *)
(*   X on a data qubit -> H       (H X H = Z)                                    *)
(*   Y on a data qubit -> V       (V Y V^-1 = Z, V = sqrt(X), the engine's "V")   *)
(*   Z on a data qubit -> nothing                                                 *)
(*                                                                              *)
(* Layout: data qubits are 1..n and the ancilla of generator j is n + j.  The     *)
(* ancillas are reset and reused every round, so an r-round circuit still runs    *)
(* on n + m qubits, not n + r m.                                                  *)
(*                                                                              *)
(* WHAT THIS CIRCUIT IS NOT, and it is not a rough edge.                          *)
(*                                                                              *)
(* Bare ancillas are the difference between fault tolerant and not, not a         *)
(* limitation alongside gate scheduling.  Book sec. 12.1.1: "a single faulty gate *)
(* can cause multiple data errors, even if the fault doesn't directly affect the  *)
(* data block", and "if the code only corrects 1 error, it isn't strong enough to *)
(* function properly".  The fix is a verified ancilla: cat states (sec. 12.1.2 -   *)
(* 12.1.3) feeding Shor EC (sec. 12.2), or Steane EC (sec. 12.3), or Knill EC     *)
(* (sec. 12.4).  Roadmap item 5.                                                  *)
(*                                                                              *)
(* Using it anyway is a deliberate, sourced choice rather than a shortcut.  Book  *)
(* sec. 12.5.1: "Frequently people using surface codes don't even bother with     *)
(* Shor EC, instead using the non-FT measurement technique of section 12.1.1      *)
(* instead.  They're willing to put up with the fact that one fault can cause     *)
(* multiple errors."  That concession is specific to the surface/LDPC setting; it  *)
(* is not a general licence, and nothing here should be read as claiming the       *)
(* extraction is fault tolerant for an arbitrary stabilizer code.                  *)
(*                                                                              *)
(* THE MECHANISM, worked out for *this* CNOT direction, because the book's is      *)
(* mirrored and the naive transcription is wrong.  Figure 12.1b has an X on the    *)
(* ancilla halfway through the ladder corrupting several data qubits, which holds  *)
(* when the ancilla is the control.  Here the ancilla is the *target*, so:         *)
(*                                                                              *)
(*   - an X on the ancilla reaches no data qubit at all.  X lives on the target    *)
(*     and CNOT only carries X from control to target, so it sits there until the  *)
(*     measurement and flips that one outcome.  It is a pure readout error.        *)
(*   - a Z (hence also a Y) on the ancilla is the one that spreads.  CNOT carries  *)
(*     Z from target back to control, and the ancilla keeps its Z, so it lands on  *)
(*     every data qubit the ladder has *still to touch*.  After the un-rotation    *)
(*     the residual is exactly the un-extracted tail of the generator: a Z mid-    *)
(*     ladder on the 5-qubit code's XZZXI leaves IZZXI, IIZXI or IIIXI depending   *)
(*     on where it struck -- and, being a sub-product of the generator, it does    *)
(*     not fire that generator's own check in that round.  A Z before the first    *)
(*     CNOT leaves the whole generator, which is a stabilizer and harmless; a Z    *)
(*     after the last leaves nothing.  Only mid-ladder hurts.                      *)
(*                                                                              *)
(* This is what turns a distance-three code's logical error rate from order p^2    *)
(* into order p, and it is a defect of the *gadget*, not a property of circuit-     *)
(* level noise: a gadget satisfying the book's gate and error-correction           *)
(* propagation properties (sec. 10.2) keeps the p^2.                               *)
(*                                                                              *)
(* SCHEDULING, and why "Idle" is zero.  The generators are extracted one after     *)
(* another rather than interleaved into parallel layers.  Do not read the zero      *)
(* idle rate as "there is no time step to idle in" -- that inference is exactly     *)
(* the one book sec. 15.5.2 forbids, since it resolves an ill-defined time step by  *)
(* *defining* one (the longest gate) and charging the padding as storage error.     *)
(* Worse, sequential extraction means *more* waiting, not less: sec. 15.5.1 shows   *)
(* partial parallelism multiplying the storage rate (p_S -> 3 p_S and              *)
(* p_G -> p_G + 2 p_S at one-third parallelism) and states that "to have a          *)
(* threshold, it is essential to do parallel gates, at least when the storage       *)
(* error rate p_S is non-zero".  So: Idle = 0 is an optimistic simplification that  *)
(* is not yet modelled, and it is optimistic in the direction this circuit is       *)
(* already weakest.  A fixed-size code extracted serially is a constant number of  *)
(* steps, which is what keeps it inside sec. 15.5.1's constant-fraction clause;     *)
(* that clause, not the absence of idle locations, is the excuse.                   *)
(*                                                                              *)
(* Reusing the ancillas across rounds is fine and is sourced.  Book sec. 15.4      *)
(* re-examines the assumption that fresh qubits can be prepared mid-computation,    *)
(* and names measure-and-flip reset as the remedy rather than the problem: reuse    *)
(* is legitimate given non-destructive mid-circuit measurement and classical        *)
(* conditioning, both of which are assumed here.  Sec. 15.4.3 adds the one caveat   *)
(* that matters -- a reset ancilla is not truly fresh, and the residual should be   *)
(* carried as a preparation error rate -- which is exactly what the noise model's   *)
(* "Reset" rate is.                                                                *)
(* ============================================================================ *)

QECSyndromeCircuit::usage = "QECSyndromeCircuit[code] gives the syndrome-extraction circuit of a stabilizer code: one ancilla per generator, measured once.\nQECSyndromeCircuit[code, r] repeats the extraction for r rounds, resetting the ancillas between them.\ncirc[prop] gives a property; circ[\"Properties\"] lists them.";

QECSyndromeCircuit::rounds = "The number of rounds must be a positive integer; got `1`.";
QECSyndromeCircuit::noprop = "`1` is not a property of QECSyndromeCircuit. Use circ[\"Properties\"] for the list.";


(* ---- instructions ---- *)

(* An instruction is a plain list {op, qubits...}, deliberately not an expression
   with a head: the frame propagator walks millions of these and a Switch on a
   string part is the cheapest dispatch available.

     {"R", q}          reset q to |0>
     {"H"|"S"|"V"|"Vdg", q}
     {"CNOT", c, t}
     {"M", q}          measure q in Z and append the outcome to the record
*)

$circuitOneQubitOps = {"H", "S", "V", "Vdg"};

(* The rotation that sends a generator's letter on a data qubit to Z, and the one
   that undoes it.  V has order four, so V^-1 = V^3; the two are separate ops here
   because the noise model counts instructions, and a rotation is one location
   whatever the hardware needs to realise it. *)
basisGate[{1, 0}] := "H"
basisGate[{1, 1}] := "V"
basisGate[_] := None

unbasisGate[{1, 0}] := "H"
unbasisGate[{1, 1}] := "Vdg"
unbasisGate[_] := None

letterAt[v_List, n_Integer, q_Integer] := {v[[q]], v[[n + q]]}

rotation[f_, v_List, n_Integer, support_List] :=
    Table[With[{g = f[letterAt[v, n, q]]}, If[g === None, Nothing, {g, q}]], {q, support}]

(* One generator, one ancilla, one round. *)
generatorInstructions[v_List, n_Integer, ancilla_Integer] := With[
    {support = Select[Range[n], letterAt[v, n, #] =!= {0, 0} &]},
    Join[
        {{"R", ancilla}},
        rotation[basisGate, v, n, support],
        Table[{"CNOT", q, ancilla}, {q, support}],
        rotation[unbasisGate, v, n, support],
        {{"M", ancilla}}
    ]
]

codeCircuitInstructions[a_Association, rounds_Integer] := codeCircuitInstructions[a, rounds] = With[
    {n = a["Qubits"], mat = a["CheckMatrix"]},
    Catenate @ Table[
        Catenate @ Table[generatorInstructions[mat[[j]], n, n + j], {j, Length[mat]}],
        {rounds}
    ]
]


(* ---- construction ---- *)

QECSyndromeCircuit[QECCode[a_Association]] := QECSyndromeCircuit[QECCode[a], 1]

QECSyndromeCircuit[QECCode[a_Association], rounds_] := If[
    ! (IntegerQ[rounds] && rounds > 0),
    Message[QECSyndromeCircuit::rounds, rounds]; $Failed,
    QECSyndromeCircuit[<|
        "Instructions" -> codeCircuitInstructions[a, rounds],
        "DataQubits" -> a["Qubits"],
        "Ancillas" -> codeStabilizerCount[a],
        "Rounds" -> rounds,
        "Code" -> a
    |>]
]


(* ---- accessors ---- *)

circuitData[QECSyndromeCircuit[a_Association]] := a
circuitInstructions[a_Association] := a["Instructions"]
circuitQubitCount[a_Association] := a["DataQubits"] + a["Ancillas"]

(* The record is filled in emission order: round 1 generator 1, round 1 generator 2,
   ... so measurement index (r-1) m + j is generator j of round r. *)
circuitMeasurementLabels[a_Association] := Catenate @ Table[
    {r, j}, {r, a["Rounds"]}, {j, a["Ancillas"]}
]


(* ---- the frame propagator ---- *)

(* The whole point of the circuit being Clifford and the noise being Pauli: a
   fault does not have to be simulated as a state, it can be carried as a Pauli
   *frame* and pushed through the gates by conjugation, which over GF(2) is four
   XORs.  The noiseless run produces a fixed outcome record; a frame's only effect
   is to flip a subset of those outcomes and to leave a residual Pauli on the data.
   That is exactly Stim's model, it is linear in the frame, and it means the whole
   noisy experiment is GF(2) bookkeeping rather than state simulation.
   Phases never enter: the frame's contribution to an outcome is whether it
   anticommutes with the measured Z, which is the ancilla's X bit.

     H     x <-> z
     S     z ^= x           (Z-rotation: X -> Y)
     V     x ^= z           (X-rotation: Y -> Z);  V^-1 acts the same on the
                            symplectic part, since V^2 is the Pauli X
     CNOT  x_t ^= x_c, z_c ^= z_t
     R     x, z := 0        (the ancilla is re-prepared, so its frame is discarded)
     M     record x         (an X on the ancilla flips the readout)

   Faults are given as {instructionIndex, qubit, {x, z}} and are applied *after*
   the instruction at that index; index 0 means before the circuit starts, which
   is where a code-capacity error on the data lands.  Readout noise is a fault on the
   ancilla placed just before its "M": the ancilla is reset immediately after, so
   flipping it there flips the record and nothing else, which is precisely what a
   classical readout error is. *)

frameZero[nq_Integer] := {ConstantArray[0, nq], ConstantArray[0, nq]}

framePropagate[instr_List, nq_Integer, faults_List] := Module[
    {x, z, out, byIndex, step, q, c, t, f},

    x = ConstantArray[0, nq];
    z = ConstantArray[0, nq];
    out = Internal`Bag[];
    byIndex = GroupBy[faults, First];

    Do[
        q = f[[2]];
        x[[q]] = BitXor[x[[q]], f[[3, 1]]];
        z[[q]] = BitXor[z[[q]], f[[3, 2]]],
        {f, Lookup[byIndex, 0, {}]}
    ];

    Do[
        step = instr[[i]];
        Switch[First[step],
            "R",    q = step[[2]]; x[[q]] = 0; z[[q]] = 0,
            "H",    q = step[[2]]; {x[[q]], z[[q]]} = {z[[q]], x[[q]]},
            "S",    q = step[[2]]; z[[q]] = BitXor[z[[q]], x[[q]]],
            "V",    q = step[[2]]; x[[q]] = BitXor[x[[q]], z[[q]]],
            "Vdg",  q = step[[2]]; x[[q]] = BitXor[x[[q]], z[[q]]],
            "CNOT", c = step[[2]]; t = step[[3]];
                    x[[t]] = BitXor[x[[t]], x[[c]]];
                    z[[c]] = BitXor[z[[c]], z[[t]]],
            "M",    Internal`StuffBag[out, x[[step[[2]]]]]
        ];
        Do[
            q = f[[2]];
            x[[q]] = BitXor[x[[q]], f[[3, 1]]];
            z[[q]] = BitXor[z[[q]], f[[3, 2]]],
            {f, Lookup[byIndex, i, {}]}
        ],
        {i, Length[instr]}
    ];

    <|"Record" -> Internal`BagPart[out, All], "Frame" -> {x, z}|>
]


(* ---- handing the circuit to the engine ---- *)

(* Only for validation: the engine simulates states, so this is how we check that
   the emitted circuit really measures the generators we think it does.  "R" has no
   unitary form and is dropped -- validation runs a single round from ancillas
   already in |0>.  "M" is not a gate either and is handled by the caller. *)
instructionEngineGates[instr_List] := Catenate[
    Replace[instr, {
        {"R", _} -> {},
        {"M", _} -> {},
        {"Vdg", q_} :> {"V" -> q, "V" -> q, "V" -> q},
        {op_, q_} :> {op -> q},
        {"CNOT", c_, t_} :> {"CNOT" -> {c, t}}
    }, {1}]
]


(* ---- properties ---- *)

$circuitProperties = {
    "Instructions", "DataQubits", "Ancillas", "Qubits", "Rounds", "Code",
    "MeasurementCount", "MeasurementLabels", "Depth", "GateCounts", "Properties"
};

QECSyndromeCircuit[_Association]["Properties"] := $circuitProperties

QECSyndromeCircuit[a_Association]["Instructions"] := a["Instructions"]
QECSyndromeCircuit[a_Association]["DataQubits"] := a["DataQubits"]
QECSyndromeCircuit[a_Association]["Ancillas"] := a["Ancillas"]
QECSyndromeCircuit[a_Association]["Qubits"] := circuitQubitCount[a]
QECSyndromeCircuit[a_Association]["Rounds"] := a["Rounds"]
QECSyndromeCircuit[a_Association]["Code"] := QECCode[a["Code"]]
QECSyndromeCircuit[a_Association]["MeasurementCount"] := a["Rounds"] a["Ancillas"]
QECSyndromeCircuit[a_Association]["MeasurementLabels"] := circuitMeasurementLabels[a]
QECSyndromeCircuit[a_Association]["Depth"] := Length[a["Instructions"]]
QECSyndromeCircuit[a_Association]["GateCounts"] := Counts[First /@ a["Instructions"]]

QECSyndromeCircuit[a_Association][prop_String] := (Message[QECSyndromeCircuit::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECSyndromeCircuit /: MakeBoxes[obj : QECSyndromeCircuit[a_Association] /; KeyExistsQ[a, "Instructions"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECSyndromeCircuit,
        obj,
        BarChart[Values[Counts[First /@ a["Instructions"]]],
            ChartLabels -> Keys[Counts[First /@ a["Instructions"]]],
            ImageSize -> {Automatic, 34}, Axes -> False, ChartStyle -> RGBColor[0.15, 0.5, 0.65]],
        {
            BoxForm`SummaryItem[{"Data qubits: ", a["DataQubits"]}],
            BoxForm`SummaryItem[{"Ancillas: ", a["Ancillas"]}],
            BoxForm`SummaryItem[{"Rounds: ", a["Rounds"]}]
        },
        {
            BoxForm`SummaryItem[{"Instructions: ", Length[a["Instructions"]]}],
            BoxForm`SummaryItem[{"Measurements: ", a["Rounds"] a["Ancillas"]}],
            BoxForm`SummaryItem[{"Gate counts: ", Counts[First /@ a["Instructions"]]}]
        },
        form,
        "Interpretable" -> False
    ]
