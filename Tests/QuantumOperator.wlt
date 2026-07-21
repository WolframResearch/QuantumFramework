BeginTestSection["QuantumOperator - constructors"]

VerificationTest[QuantumOperator[]["Dimensions"], {2, 2}, TestID -> "Empty"]

VerificationTest[QuantumOperator["X"]["Dimensions"], {2, 2}, TestID -> "X-bare"]

VerificationTest[QuantumOperator["X"[3]]["Dimensions"], {3, 3}, TestID -> "X-3"]

VerificationTest[QuantumOperator["Hadamard"]["Dimensions"], {2, 2}, TestID -> "Hadamard"]

VerificationTest[QuantumOperator["CNOT"]["Dimensions"], {2, 2, 2, 2}, TestID -> "CNOT"]

(* CNOT[3] gives a 2-dim control with a 3-dim target (the qutrit "shift" gate),
   matching pre-refactor semantics. *)
VerificationTest[QuantumOperator["CNOT"[3]]["Dimensions"], {2, 3, 2, 3}, TestID -> "CNOT-3"]

VerificationTest[QuantumOperator["Toffoli"]["Dimensions"], {2, 2, 2, 2, 2, 2}, TestID -> "Toffoli"]

VerificationTest[QuantumOperator["Fourier"[3]]["Dimensions"], {3, 3}, TestID -> "Fourier-3"]

VerificationTest[QuantumOperator["XRotation"[Pi/3]]["Dimensions"], {2, 2}, TestID -> "XRotation"]

VerificationTest[QuantumOperator["Phase"[Pi/4]]["Dimensions"], {2, 2}, TestID -> "Phase"]

VerificationTest[QuantumOperator["Permutation"[{2, 2}, Cycles[{{1, 2}}]]]["Dimensions"], {2, 2, 2, 2}, TestID -> "Permutation"]

VerificationTest[QuantumOperator["Spider"]["Dimensions"], {2, 2}, TestID -> "Spider-bare"]

VerificationTest[QuantumOperator["Curry"]["Dimensions"], {4, 2, 2}, TestID -> "Curry-default"]

VerificationTest[QuantumOperator["XSpider"[Pi/2], {{1, 2}, {3}}]["Dimensions"], {2, 2, 2}, TestID -> "XSpider-with-order"]

VerificationTest[
    QuantumOperator["Liouvillian"[QuantumOperator["X"]]]["Dimensions"],
    {2, 2},
    TestID -> "Liouvillian"
]

(* explicit matrix *)
VerificationTest[
    QuantumOperator[{{0, 1}, {1, 0}}]["Dimensions"],
    {2, 2},
    TestID -> "Explicit-matrix"
]

VerificationTest[QuantumOperator["XYZ"]["Dimensions"], {2, 2, 2, 2, 2, 2}, TestID -> "PauliString-XYZ"]

EndTestSection[]


BeginTestSection["QuantumOperator - properties"]

VerificationTest[QuantumOperator["X"]["UnitaryQ"], True, TestID -> "X-Unitary"]

VerificationTest[QuantumOperator["Hadamard"]["UnitaryQ"], True, TestID -> "Hadamard-Unitary"]

VerificationTest[QuantumOperator["CNOT"]["Arity"], 2, TestID -> "CNOT-Arity"]

VerificationTest[
    Total @ Total[QuantumOperator["X"]["MatrixRepresentation"] - {{0, 1}, {1, 0}}],
    0,
    TestID -> "X-Matrix"
]

EndTestSection[]


BeginTestSection["QuantumOperator - composition"]

VerificationTest[
    Total @ Total[(QuantumOperator["X"] @ QuantumOperator["X"])["MatrixRepresentation"] - IdentityMatrix[2]],
    0,
    TestID -> "XX=I"
]

VerificationTest[
    Total @ Chop @ Flatten[QuantumOperator["X"][QuantumState["0"]]["StateVector"] - {0, 1}],
    0,
    TestID -> "X-on-Zero"
]

EndTestSection[]


BeginTestSection["QuantumOperator - shortcut roundtrip"]

(* For each common named gate, QuantumShortcut should emit the new call-form
   shorthand, and feeding it through the head should reconstruct an operator
   with the same matrix. *)

roundtripQ[op_] := Module[{recovered = Head[op][First[QuantumShortcut[op]]]},
    Chop[Norm[Flatten[op["Matrix"] - recovered["Matrix"]]]] === 0
]

VerificationTest[roundtripQ[QuantumOperator["X"]], True, TestID -> "Roundtrip-X"]
VerificationTest[roundtripQ[QuantumOperator["X"[3]]], True, TestID -> "Roundtrip-X-3"]
VerificationTest[roundtripQ[QuantumOperator["Y"]], True, TestID -> "Roundtrip-Y"]
VerificationTest[roundtripQ[QuantumOperator["Z"]], True, TestID -> "Roundtrip-Z"]
VerificationTest[roundtripQ[QuantumOperator["Hadamard"]], True, TestID -> "Roundtrip-Hadamard"]
VerificationTest[roundtripQ[QuantumOperator["NOT"]], True, TestID -> "Roundtrip-NOT"]
VerificationTest[roundtripQ[QuantumOperator["S"]], True, TestID -> "Roundtrip-S"]
VerificationTest[roundtripQ[QuantumOperator["T"]], True, TestID -> "Roundtrip-T"]
VerificationTest[roundtripQ[QuantumOperator["Phase"[Pi/4]]], True, TestID -> "Roundtrip-Phase"]
VerificationTest[roundtripQ[QuantumOperator["PhaseShift"[2]]], True, TestID -> "Roundtrip-PhaseShift"]
VerificationTest[roundtripQ[QuantumOperator["U"[Pi/3, Pi/4, Pi/5]]], True, TestID -> "Roundtrip-U"]
VerificationTest[roundtripQ[QuantumOperator["U2"[Pi/4, Pi/3]]], True, TestID -> "Roundtrip-U2"]

(* QuantumShortcut emits the new call form, not the legacy list form, for any
   parameterized gate. *)
VerificationTest[
    QuantumShortcut[QuantumOperator["U"[Pi/3, Pi/4, Pi/5]]],
    {"U"[Pi/3, Pi/4, Pi/5] -> {1}},
    TestID -> "Shortcut-U-callform"
]

VerificationTest[
    QuantumShortcut[QuantumOperator["X"[3]]],
    {"X"[3] -> {1}},
    TestID -> "Shortcut-X3-callform"
]

EndTestSection[]


BeginTestSection["QuantumOperator - composition fast path"]

(* Fast path for qo1[qo2] avoids QuantumCircuitOperator/TensorNetwork by direct
   tensor contraction. Each test compares fast.Sort.Matrix to the slow route
   through QuantumCircuitOperator. *)

slowCompose[qo1_, qo2_] := QuantumOperator @ QuantumCircuitOperator[{qo2, qo1}]
cmpCompose[qo1_, qo2_] := Chop @ Norm @ Flatten[N @ qo1[qo2]["Sort"]["Matrix"] - N @ slowCompose[qo1, qo2]["Sort"]["Matrix"]]

(* Aligned 1-qubit *)
VerificationTest[cmpCompose[QuantumOperator["H"], QuantumOperator["H"]], 0, TestID -> "FastCompose-H-H"]
VerificationTest[cmpCompose[QuantumOperator["X"], QuantumOperator["H"]], 0, TestID -> "FastCompose-X-H"]
VerificationTest[cmpCompose[QuantumOperator["Y"], QuantumOperator["Z"]], 0, TestID -> "FastCompose-Y-Z"]
VerificationTest[cmpCompose[QuantumOperator["S"], QuantumOperator["T"]], 0, TestID -> "FastCompose-S-T"]
VerificationTest[cmpCompose[QuantumOperator["RX"[Pi/3]], QuantumOperator["RY"[Pi/4]]], 0, TestID -> "FastCompose-RX-RY"]

(* Disjoint qudit positions - tensor product *)
VerificationTest[cmpCompose[QuantumOperator["X"], QuantumOperator["H", {2}]], 0, TestID -> "FastCompose-X1-H2"]
VerificationTest[cmpCompose[QuantumOperator["H", {2}], QuantumOperator["X"]], 0, TestID -> "FastCompose-H2-X1"]
VerificationTest[cmpCompose[QuantumOperator["X"], QuantumOperator["Y", {2}]], 0, TestID -> "FastCompose-X1-Y2"]

(* 2-qubit on 2-qubit aligned *)
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["CNOT"]], 0, TestID -> "FastCompose-CNOT-CNOT"]
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["SWAP"]], 0, TestID -> "FastCompose-CNOT-SWAP"]
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["CZ"]], 0, TestID -> "FastCompose-CNOT-CZ"]

(* 1-qubit on 2-qubit, partial overlap *)
VerificationTest[cmpCompose[QuantumOperator["X"], QuantumOperator["CNOT"]], 0, TestID -> "FastCompose-X1-CNOT"]
VerificationTest[cmpCompose[QuantumOperator["X", {2}], QuantumOperator["CNOT"]], 0, TestID -> "FastCompose-X2-CNOT"]
VerificationTest[cmpCompose[QuantumOperator["H"], QuantumOperator["SWAP"]], 0, TestID -> "FastCompose-H1-SWAP"]
VerificationTest[cmpCompose[QuantumOperator["H", {2}], QuantumOperator["SWAP"]], 0, TestID -> "FastCompose-H2-SWAP"]

(* 2-qubit on 1-qubit, symmetric direction *)
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["X"]], 0, TestID -> "FastCompose-CNOT-X1"]
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["H"]], 0, TestID -> "FastCompose-CNOT-H1"]
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["H", {2}]], 0, TestID -> "FastCompose-CNOT-H2"]
VerificationTest[cmpCompose[QuantumOperator["SWAP"], QuantumOperator["X"]], 0, TestID -> "FastCompose-SWAP-X1"]

(* Mixed-position 2-qubit operators *)
VerificationTest[cmpCompose[QuantumOperator["CNOT", {1, 3}], QuantumOperator["H"]], 0, TestID -> "FastCompose-CNOT13-H1"]
VerificationTest[cmpCompose[QuantumOperator["CNOT", {2, 3}], QuantumOperator["H", {2}]], 0, TestID -> "FastCompose-CNOT23-H2"]
VerificationTest[cmpCompose[QuantumOperator["CNOT", {1, 3}], QuantumOperator["CNOT", {2, 3}]], 0, TestID -> "FastCompose-CNOT13-CNOT23"]
VerificationTest[cmpCompose[QuantumOperator["CNOT", {2, 3}], QuantumOperator["CNOT", {1, 2}]], 0, TestID -> "FastCompose-CNOT23-CNOT12"]

(* 3-qubit operators *)
VerificationTest[cmpCompose[QuantumOperator["Toffoli"], QuantumOperator["CNOT"]], 0, TestID -> "FastCompose-Toffoli-CNOT"]
VerificationTest[cmpCompose[QuantumOperator["Toffoli"], QuantumOperator["H", {3}]], 0, TestID -> "FastCompose-Toffoli-H3"]
VerificationTest[cmpCompose[QuantumOperator["Toffoli"], QuantumOperator["Toffoli"]], 0, TestID -> "FastCompose-Toffoli-Toffoli"]

(* Disjoint composition with multi-qubit gates *)
VerificationTest[cmpCompose[QuantumOperator["CNOT"], QuantumOperator["H", {3}]], 0, TestID -> "FastCompose-CNOT-H3"]
VerificationTest[cmpCompose[QuantumOperator["H", {3}], QuantumOperator["CNOT"]], 0, TestID -> "FastCompose-H3-CNOT"]

(* Bra/ket as a QuantumOperator *)
VerificationTest[cmpCompose[QuantumOperator[QuantumState[{1, 0}]["Dagger"]], QuantumOperator["X"]], 0, TestID -> "FastCompose-bra-X"]
VerificationTest[cmpCompose[QuantumOperator["X"], QuantumOperator[QuantumState[{1, 0}]]], 0, TestID -> "FastCompose-X-ket"]

(* PhaseSpace pictures - both operators in PhaseSpace, the Picture-equality
   guard fires and the result picks up the matching picture. *)
With[{hp = QuantumOperator[QuantumOperator["H"], "Picture" -> "PhaseSpace"]},
    VerificationTest[cmpCompose[hp, hp], 0, TestID -> "FastCompose-PhaseSpace-Hp-Hp"];
]

With[{hp = QuantumOperator[QuantumOperator["H"], "Picture" -> "PhaseSpace"], xp = QuantumOperator[QuantumOperator["X"], "Picture" -> "PhaseSpace"]},
    VerificationTest[cmpCompose[xp, hp], 0, TestID -> "FastCompose-PhaseSpace-Xp-Hp"];
]

With[{cnotp = QuantumOperator[QuantumOperator["CNOT"], "Picture" -> "PhaseSpace"]},
    VerificationTest[cmpCompose[cnotp, cnotp], 0, TestID -> "FastCompose-PhaseSpace-CNOTp-CNOTp"];
]

(* Matrix-form operators (channels) - dispatched via Bend/Unbend. *)
mkMatrixOp[op_] := QuantumOperator[QuantumState[op["State"]]["MatrixState"], op["Order"]]

With[{hMix = mkMatrixOp @ QuantumOperator["H"]},
    VerificationTest[cmpCompose[hMix, hMix], 0, TestID -> "FastCompose-Matrix-H-H"];
]

With[{hMix = mkMatrixOp @ QuantumOperator["H"], xMix = mkMatrixOp @ QuantumOperator["X"]},
    VerificationTest[cmpCompose[xMix, hMix], 0, TestID -> "FastCompose-Matrix-X-H"];
]

With[{cnotMix = mkMatrixOp @ QuantumOperator["CNOT"]},
    VerificationTest[cmpCompose[cnotMix, cnotMix], 0, TestID -> "FastCompose-Matrix-CNOT-CNOT"];
]

With[{x2Mix = mkMatrixOp @ QuantumOperator["X", {2}], cnotMix = mkMatrixOp @ QuantumOperator["CNOT"]},
    VerificationTest[cmpCompose[x2Mix, cnotMix], 0, TestID -> "FastCompose-Matrix-X2-CNOT"];
]

With[{cnotMix = mkMatrixOp @ QuantumOperator["CNOT"], h2Mix = mkMatrixOp @ QuantumOperator["H", {2}]},
    VerificationTest[cmpCompose[cnotMix, h2Mix], 0, TestID -> "FastCompose-Matrix-CNOT-H2-symmetric"];
]

(* Vector + matrix mixed *)
With[{hVec = QuantumOperator["H"], hMix = mkMatrixOp @ QuantumOperator["H"]},
    VerificationTest[cmpCompose[hVec, hMix], 0, TestID -> "FastCompose-Mixed-vector-matrix"];
    VerificationTest[cmpCompose[hMix, hVec], 0, TestID -> "FastCompose-Mixed-matrix-vector"];
]

(* Disjoint matrix tensor product *)
With[{xMix = mkMatrixOp @ QuantumOperator["X"], h2Mix = mkMatrixOp @ QuantumOperator["H", {2}]},
    VerificationTest[cmpCompose[xMix, h2Mix], 0, TestID -> "FastCompose-Matrix-disjoint"];
]

EndTestSection[]


BeginTestSection["QuantumOperator - broadcast"]

(* Small multiplicity broadcasts: an operator given an order longer than its
   qudit count tensors copies of itself across the order. *)
VerificationTest[
    QuantumOperator["H", {1, 2, 3}]["Dimensions"],
    {2, 2, 2, 2, 2, 2},
    TestID -> "Broadcast-H-3qubits"
]

VerificationTest[
    QuantumOperator["H", {1, 2, 3}]["MatrixRepresentation"],
    QuantumTensorProduct[QuantumOperator["H", {1}], QuantumOperator["H", {2}], QuantumOperator["H", {3}]]["MatrixRepresentation"],
    TestID -> "Broadcast-H-3qubits-matrix"
]

VerificationTest[
    QuantumOperator["CNOT", {1, 2, 3, 4}]["Dimensions"],
    {2, 2, 2, 2, 2, 2, 2, 2},
    TestID -> "Broadcast-CNOT-2x"
]

(* A state-shaped operator (ket, no input qudits) broadcasts too, as long as
   the implied tensor power stays below $QuantumOperatorBroadcastLimit. *)
VerificationTest[
    QuantumOperator[QuantumOperator[QuantumState["Register"[2]]], Range[4]]["Dimensions"],
    {2, 2, 2, 2, 2, 2, 2, 2},
    TestID -> "Broadcast-ket-small"
]

(* A 4-qubit ket over a length-12 order implies a 16^12-amplitude tensor power;
   the guard turns the kernel-killing materialization into a message. *)
VerificationTest[
    QuantumOperator[QuantumOperator[QuantumState["Register"[4]]], Range[12]],
    $Failed,
    {QuantumOperator::broadcast},
    TestID -> "Broadcast-limit-ket-12"
]

(* Same guard on a plain gate object: H over 14 wires implies dimension 4^14 > 2^24.
   (The named form QuantumOperator["H", order] takes the Fourier construction
   path instead and never reaches the broadcast rule.) *)
VerificationTest[
    QuantumOperator[QuantumOperator["H"], Range[14]],
    $Failed,
    {QuantumOperator::broadcast},
    TestID -> "Broadcast-limit-H-14"
]

EndTestSection[]


BeginTestSection["QuantumOperator - reorder"]

(* A single flat order is a *placement*: it relabels the operator's qudit
   footprint consistently, so an operator whose output and input orders differ
   (a permutation / swap) keeps its action. A two-element {out, in} order is the
   low-level *positional* form and may intentionally rewire the legs. *)

ordPair[o_] := {o["OutputOrder"], o["InputOrder"]}
orderedMat[o_] := Normal[o["OrderedMatrixRepresentation"]]
$swap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}

(* "Permutation" -> {2, 1} is a SWAP encoded as output order {2, 1} vs input
   order {1, 2}: the swap lives in the order mismatch, not the raw tensor. *)
VerificationTest[
    ordPair[QuantumOperator["Permutation" -> {2, 1}]],
    {{2, 1}, {1, 2}},
    TestID -> "Reorder-permutation-native-orders"
]

(* Placement onto {2, 3} must preserve the swap (regression: it used to collapse
   to the identity because both orders were clobbered to the same list). *)
With[{perm = QuantumOperator["Permutation" -> {2, 1}]},
    VerificationTest[
        ordPair[QuantumOperator[perm, {2, 3}]],
        {{3, 2}, {2, 3}},
        TestID -> "Reorder-placement-permutation-orders"
    ];
    VerificationTest[
        orderedMat[QuantumOperator[perm, {2, 3}]],
        $swap,
        TestID -> "Reorder-placement-permutation-keeps-SWAP"
    ];
    (* non-adjacent target qudits *)
    VerificationTest[
        orderedMat[QuantumOperator[perm, {5, 9}]],
        $swap,
        TestID -> "Reorder-placement-permutation-nonadjacent"
    ];
    (* placement onto its own footprint is a no-op *)
    VerificationTest[
        ordPair[QuantumOperator[perm, {1, 2}]],
        {{2, 1}, {1, 2}},
        TestID -> "Reorder-placement-identity-noop"
    ];
    (* round-trip: place away then back preserves the matrix *)
    VerificationTest[
        orderedMat[QuantumOperator[QuantumOperator[perm, {5, 9}], {1, 2}]],
        $swap,
        TestID -> "Reorder-placement-roundtrip-SWAP"
    ]
]

(* A 3-qudit cyclic permutation relabels its whole footprint and keeps its
   matrix. *)
With[{p3 = QuantumOperator["Permutation" -> {2, 3, 1}]},
    VerificationTest[
        ordPair[QuantumOperator[p3, {5, 6, 7}]],
        {{6, 7, 5}, {5, 6, 7}},
        TestID -> "Reorder-placement-3cycle-orders"
    ];
    VerificationTest[
        orderedMat[QuantumOperator[p3, {5, 6, 7}]] === orderedMat[p3],
        True,
        TestID -> "Reorder-placement-3cycle-keeps-matrix"
    ]
]

(* Ordinary operators (output order == input order) are unaffected: placement is
   an ordinary relabel and matches the historical {order, order} behavior. *)
VerificationTest[
    ordPair[QuantumOperator[QuantumOperator["X"], {3}]],
    {{3}, {3}},
    TestID -> "Reorder-placement-single-qubit"
]

With[{cnot = QuantumOperator["CNOT"]},
    VerificationTest[
        ordPair[QuantumOperator[cnot, {2, 3}]],
        {{2, 3}, {2, 3}},
        TestID -> "Reorder-placement-CNOT-orders"
    ];
    VerificationTest[
        orderedMat[QuantumOperator[cnot, {2, 3}]] === orderedMat[cnot],
        True,
        TestID -> "Reorder-placement-CNOT-keeps-matrix"
    ]
]

(* Explicit two-element {out, in} order is positional/low-level: literally
   setting output order == input order == {2, 3} rewires the swap into the
   identity. This is intentional and distinct from single-order placement. *)
With[{perm = QuantumOperator["Permutation" -> {2, 1}]},
    VerificationTest[
        ordPair[QuantumOperator[perm, {{2, 3}, {2, 3}}]],
        {{2, 3}, {2, 3}},
        TestID -> "Reorder-explicit-positional-orders"
    ];
    VerificationTest[
        orderedMat[QuantumOperator[perm, {{2, 3}, {2, 3}}]],
        IdentityMatrix[4],
        TestID -> "Reorder-explicit-positional-identity"
    ]
]

(* Output and input orders set independently. *)
VerificationTest[
    ordPair[QuantumOperator[QuantumOperator["X"], {{2}, {3}}]],
    {{2}, {3}},
    TestID -> "Reorder-explicit-out-in"
]

(* Arrow form order1 -> order2 sets {output -> order2, input -> order1}. *)
VerificationTest[
    ordPair[QuantumOperator[QuantumOperator["X"], {2} -> {3}]],
    {{3}, {2}},
    TestID -> "Reorder-arrow-swaps-out-in"
]

(* "Reorder" property: a bare order changes only the output order. *)
With[{cnot = QuantumOperator["CNOT"]},
    VerificationTest[
        ordPair[cnot["Reorder", {2, 3}]],
        {{2, 3}, {1, 2}},
        TestID -> "Reorder-property-output-only"
    ];
    VerificationTest[
        ordPair[cnot["Reorder", {{2, 3}, Automatic}]],
        {{2, 3}, {1, 2}},
        TestID -> "Reorder-property-Automatic-input"
    ]
]

(* "Shift" adds a constant to every qudit index and preserves the matrix. *)
VerificationTest[
    ordPair[QuantumOperator["CNOT"]["Shift", 5]],
    {{6, 7}, {6, 7}},
    TestID -> "Reorder-shift-CNOT"
]

With[{perm = QuantumOperator["Permutation" -> {2, 1}]},
    VerificationTest[
        ordPair[perm["Shift", 1]],
        {{3, 2}, {2, 3}},
        TestID -> "Reorder-shift-permutation-orders"
    ];
    VerificationTest[
        orderedMat[perm["Shift", 1]],
        $swap,
        TestID -> "Reorder-shift-permutation-keeps-SWAP"
    ]
]

(* An order longer than the footprint is a broadcast, not a placement: the guard
   on the single-order rule lets it fall through to the multiplicity path. *)
VerificationTest[
    QuantumOperator[QuantumOperator["X"], {1, 2, 3}]["Dimensions"],
    {2, 2, 2, 2, 2, 2},
    TestID -> "Reorder-longer-order-broadcasts"
]

EndTestSection[]


BeginTestSection["QuantumOperator - failure"]

VerificationTest[
    Quiet @ QuantumOperator["NotAnActualGate"[3]],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

(* Known name with unmatched call shape returns Failure["InvalidArguments"]
   rather than recursing or falling through unevaluated. *)
VerificationTest[
    Quiet @ QuantumOperator["X"["bad-arg"]],
    Failure["InvalidArguments", _],
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-X"
]

EndTestSection[]


BeginTestSection["QuantumOperator - cross-basis composition"]

(* Composition is composition of physical operators: whatever basis frame each
   operand carries, rep(A @ B) = rep(A) . rep(B) in the computational frame.
   The direct contraction fast path glues A's input wires to B's output wires,
   which is faithful only when both sides carry the same basis on every shared
   wire; a mismatched frame must rebase (via the circuit route) instead of
   contracting raw tensors. The tagged diagonals below are sigma_x and sigma_z,
   so their product is sigma_x . sigma_z = -I sigma_y. *)

With[{
    xp = QuantumOperator[DiagonalMatrix[{1, -1}], QuantumBasis["PauliX"]],
    zp = QuantumOperator[DiagonalMatrix[{1, -1}], QuantumBasis["PauliZ"]]
},
    VerificationTest[
        Normal @ Simplify @ xp[zp]["MatrixRepresentation"],
        {{0, -1}, {1, 0}},
        TestID -> "CrossBasis-PauliXZ-product"
    ];
    VerificationTest[
        Normal @ Simplify @ zp[xp]["MatrixRepresentation"],
        Simplify[Normal[zp["MatrixRepresentation"]] . Normal[xp["MatrixRepresentation"]]],
        TestID -> "CrossBasis-PauliZX-product"
    ];
    (* symbolic angle: the contract is basis-independent for parametric operators *)
    With[{rxp = QuantumOperator[QuantumOperator["RX"[\[Theta]]]["Matrix"], QuantumBasis["PauliX"]]},
        VerificationTest[
            Simplify[Normal[rxp[zp]["MatrixRepresentation"]] - Normal[rxp["MatrixRepresentation"]] . Normal[zp["MatrixRepresentation"]]],
            {{0, 0}, {0, 0}},
            TestID -> "CrossBasis-symbolic-RX"
        ]
    ];
    (* same-basis composition stays exact and keeps its frame (the fast path) *)
    With[{xflip = QuantumOperator[{{0, 1}, {1, 0}}, QuantumBasis["PauliX"]]},
        VerificationTest[
            Normal @ Simplify @ xp[xflip]["MatrixRepresentation"],
            Simplify[Normal[xp["MatrixRepresentation"]] . Normal[xflip["MatrixRepresentation"]]],
            TestID -> "CrossBasis-same-frame-faithful"
        ];
        VerificationTest[
            xp[xflip]["Output"]["ComputationalQ"],
            False,
            TestID -> "CrossBasis-same-frame-preserved"
        ]
    ];
    (* partial wire overlap: a 2-qubit X-frame operator glued to a 1-qubit Z-frame
       operator on wire 2 only *)
    With[{
        a = QuantumOperator[KroneckerProduct[DiagonalMatrix[{1, -1}], DiagonalMatrix[{1, 1}]], {1, 2}, QuantumBasis[{"PauliX", "PauliX"}]],
        b = QuantumOperator[DiagonalMatrix[{1, -1}], {2}, QuantumBasis["PauliZ"]]
    },
        VerificationTest[
            Normal @ Simplify @ a[b]["Sort"]["MatrixRepresentation"],
            Simplify[Normal[a["MatrixRepresentation"]] . KroneckerProduct[IdentityMatrix[2], Normal[b["MatrixRepresentation"]]]],
            TestID -> "CrossBasis-partial-overlap"
        ]
    ];
    (* matrix-form operand: composing then applying equals sequential application *)
    With[{
        m = QuantumOperator[QuantumState[{{2, 0, 0, 1 + I}, {0, 1, 0, 0}, {0, 0, 1, 0}, {1 - I, 0, 0, 2}}, QuantumBasis[QuditBasis["PauliX"], QuditBasis["PauliX"]]]],
        psi = QuantumState[{3, 4 I}/5]
    },
        VerificationTest[
            Simplify[Normal[m[zp][psi]["DensityMatrix"]] - Normal[m[zp[psi]]["DensityMatrix"]]],
            {{0, 0}, {0, 0}},
            TestID -> "CrossBasis-matrix-form-sequential"
        ]
    ];
    (* disjoint wires: no glue, plain tensor product of the two frames *)
    With[{
        c1 = QuantumOperator[DiagonalMatrix[{1, -1}], {1}, QuantumBasis["PauliX"]],
        c2 = QuantumOperator[DiagonalMatrix[{1, -1}], {2}, QuantumBasis["PauliZ"]]
    },
        VerificationTest[
            Normal @ Simplify @ c1[c2]["Sort"]["MatrixRepresentation"],
            Simplify[KroneckerProduct[Normal[c1["MatrixRepresentation"]], Normal[c2["MatrixRepresentation"]]]],
            TestID -> "CrossBasis-disjoint-tensor"
        ]
    ]
]

EndTestSection[]


BeginTestSection["QuantumOperator - transpose with complex-element bases"]

(* comp(op^T) must equal Transpose[comp(op)] in any frame; complex eigenbases
   (PauliY, Fourier) are the discriminating cases, a real frame cannot see a
   dropped conjugation. *)

VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["Y"], QuantumBasis["PauliY"]]},
        Simplify[Normal[op["Transpose"]["MatrixRepresentation"]] - Transpose[Normal[op["MatrixRepresentation"]]]]
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Transpose-PauliY-frame-faithful"
]

(* Closed form: Y^T = -Y in every faithful frame. *)
VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["Y"], QuantumBasis["PauliY"]]},
        Simplify[Normal[op["Transpose"]["MatrixRepresentation"]]]
    ],
    {{0, I}, {-I, 0}},
    TestID -> "Transpose-PauliY-YT-is-minus-Y"
]

VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["Z"[3]], QuantumBasis["Fourier"[3]]]},
        Simplify[Normal[op["Transpose"]["MatrixRepresentation"]] - Transpose[Normal[op["MatrixRepresentation"]]]]
    ],
    ConstantArray[0, {3, 3}],
    TestID -> "Transpose-Fourier3-frame-faithful"
]

(* Transpose is an involution. *)
VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["Y"], QuantumBasis["PauliY"]]},
        Simplify[Normal[op["Transpose"]["Transpose"]["MatrixRepresentation"]] - Normal[op["MatrixRepresentation"]]]
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Transpose-involution-PauliY"
]

(* Conjugate after Transpose is the dagger, elementwise in the computational
   frame. *)
VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["Y"], QuantumBasis["PauliY"]]},
        Simplify[Normal[op["Transpose"]["Conjugate"]["MatrixRepresentation"]] - ConjugateTranspose[Normal[op["MatrixRepresentation"]]]]
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Transpose-Conjugate-is-dagger-PauliY"
]

(* Real frame: transposition needs no conjugation there, a regression guard. *)
VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["X"], QuantumBasis["PauliX"]]},
        Simplify[Normal[op["Transpose"]["MatrixRepresentation"]] - Transpose[Normal[op["MatrixRepresentation"]]]]
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Transpose-real-frame-regression"
]

(* Pair-form partial transpose of a two-qubit operator on the second pair. *)
VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["CNOT"], QuantumBasis[{"PauliY", "PauliY"}]]},
        With[{
            m = Normal[op["MatrixRepresentation"]],
            mt = Normal[op["Transpose", {{2, 2}}]["MatrixRepresentation"]]
        },
            Simplify[mt - ArrayReshape[Transpose[ArrayReshape[m, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]]
        ]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "PartialTranspose-operator-PauliY-pair2"
]

(* CHARACTERIZATION (contract boundary, not an endorsement): the conjugation
   convention makes transposition frame-faithful for orthonormal frames, where
   Inverse[Conjugate[E]] equals Transpose[E]. A non-orthonormal frame is
   outside the contract; this test flips if that boundary moves. *)
VerificationTest[
    With[{op = QuantumOperator[QuantumOperator["Y"], QuantumBasis[QuditBasis[{{1, I}, {1, 2}}]]]},
        Simplify[Normal[op["Transpose"]["MatrixRepresentation"]] - Transpose[Normal[op["MatrixRepresentation"]]]] === {{0, 0}, {0, 0}}
    ],
    False,
    TestID -> "Transpose-nonorthonormal-out-of-contract"
]

(* Bending a ket leg into an operator input transposes that leg: the
   computational matrix of the bent operator is the reshape of the
   computational vector. *)
VerificationTest[
    With[{qs = QuantumState[Normalize[{1, 2 I, 3, 4 I}], QuantumBasis[{"PauliY", "PauliY"}]]},
        Simplify[
            Normal[QuantumOperator[qs, {{1}, {1}}]["MatrixRepresentation"]] -
            ArrayReshape[Normal[QuantumState[qs, QuantumBasis[4]]["StateVector"]], {2, 2}]
        ]
    ],
    ConstantArray[0, {2, 2}],
    TestID -> "SplitDual-bend-reshape-contract-PauliY"
]

EndTestSection[]
