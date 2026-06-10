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
