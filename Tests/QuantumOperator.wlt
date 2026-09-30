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


BeginTestSection["QuantumOperator - normalized eigenvectors of exact operators"]

(* "Eigenvalues", "Eigenvectors" and "Eigensystem" scale each eigenvector to unit norm
   and fix its phase by first dividing it by one of its entries, so for an exact
   operator that entry has to be told from zero exactly. The first entry of an
   eigenvector of the 5 x 5 Fourier transform for the eigenvalue I or -I is zero
   without being written 0; a test that leaves x != 0 undecided puts an unevaluated If
   inside the vector, on which "Projectors" fails with SparseArray::rect. The
   normalized vectors must also stay near the size of the vectors Eigensystem returns,
   and their machine values must keep the digits of their exact values. Identities
   between the vectors are decided exactly; sums of products of them are zero without
   being written 0 and take minutes to decide exactly, so those are evaluated to 30
   digits of accuracy. *)

VerificationTest[
    With[{qo = QuantumOperator["Fourier"[5]]},
        {m = Normal[qo["MatrixRepresentation"]], es = qo["Eigensystem"]},
        {
            FreeQ[es, _If],
            Sort[Pick[First[es], (First[#] === 0) & /@ Last[es]]],
            AllTrue[Last[es], PossibleZeroQ[Conjugate[#] . # - 1, Method -> "ExactAlgebraics"] &],
            AllTrue[Flatten[MapThread[m . #2 - #1 #2 &, es]], PossibleZeroQ[#, Method -> "ExactAlgebraics"] &],
            AllTrue[Last[es], With[{x = FirstCase[#, Except[0]]},
                ! PossibleZeroQ[x, Method -> "ExactAlgebraics"] && PossibleZeroQ[x - Abs[x], Method -> "ExactAlgebraics"]] &]
        }
    ],
    {True, {-I, I}, True, True, True},
    {},
    TestID -> "Eigensystem-Fourier5-UnitEigenvectorsWithPositiveFirstNonzeroEntry"
]

(* The eigenvalue 1 of the Fourier transform is repeated, and "Projectors" takes an
   orthonormal basis of its eigenspace. Each projector P is checked against its own
   eigenvalue Tr[m . P], since two eigensystem calls on one exact matrix need not list
   the eigenvalues in the same order. A zero whose terms cancel can need more working
   precision than the default $MaxExtraPrecision allows to reach 30 digits of accuracy. *)
VerificationTest[
    With[{qo = QuantumOperator["Fourier"[5]]},
        {m = Normal[qo["MatrixRepresentation"]], p = Normal /@ qo["Projectors"]},
        Block[{$MaxExtraPrecision = 1000},
            Max[Abs[N[Flatten[{Total[p] - IdentityMatrix[5], (# . # - #) & /@ p, (m . # - Tr[m . #] #) & /@ p}], {Infinity, 30}]]]
        ]
    ],
    _ ? (# < 10^-25 &),
    {},
    SameTest -> MatchQ,
    TestID -> "Projectors-Fourier5-ResolveIdentity"
]

VerificationTest[
    With[{vectors = Last @ QuantumOperator["Fourier"[5]]["Eigensystem"]},
        Max[Abs[N[vectors] - N[vectors, 40]]]
    ],
    _ ? (# < 10^-12 &),
    {},
    SameTest -> MatchQ,
    TestID -> "Eigensystem-Fourier5-MachineValuesKeepTheirDigits"
]

(* diag(1, 1, 1, -1, -1, 2, 2, 0) in the Fourier basis of three qubits, whose entries lie
   in Q(i, Sqrt[2]) *)
VerificationTest[
    With[{f = FourierMatrix[8]},
        {m = f . DiagonalMatrix[{1, 1, 1, -1, -1, 2, 2, 0}] . ConjugateTranspose[f]},
        {qo = QuantumOperator[m, {1, 2, 3}]},
        {es = qo["Eigensystem"], unnormalized = Last @ qo["Eigensystem", "Normalize" -> False]},
        {
            LeafCount[Last[es]] < 4 LeafCount[unnormalized],
            AllTrue[Last[es], PossibleZeroQ[Conjugate[#] . # - 1, Method -> "ExactAlgebraics"] &],
            AllTrue[Flatten[MapThread[m . #2 - #1 #2 &, es]], PossibleZeroQ[#, Method -> "ExactAlgebraics"] &]
        }
    ],
    {True, True, True},
    {},
    TestID -> "Eigensystem-FourierRotated8-VectorsStaySmall"
]

(* An eigenvector of Gaussian rationals, a symbolic one and an inexact one are divided by
   their first entry when it is a nonzero number and are otherwise only scaled; so is an
   exact one whose first entry is not zero and whose norm Normalize already writes
   compactly, as for the eigenbases of the qudit Pauli X and of the spin matrices. The eigenvalue and eigenvector pairs are
   compared as sets, since two eigensystem calls on one exact matrix need not list them in
   the same order. *)
VerificationTest[
    Function[{m, sort}, With[{es = eigensystem[m, "Sort" -> sort]},
        Sort[Transpose[eigensystem[m, "Normalize" -> True, "Sort" -> sort]]] ===
            Sort[Transpose[{First[es], Normalize[If[NumericQ[First[#]] && First[#] != 0, # / First[#], #]] & /@ Last[es]}]]
    ]] @@@ {
        {{{0, -I}, {I, 0}}, False},
        {{{1, 0, 0}, {0, 0, -I}, {0, I, 0}}, False},
        {FourierMatrix[4], False},
        {{{\[FormalA], 1}, {1, -\[FormalA]}}, False},
        {N[FourierMatrix[5]], False},
        {N[{{1, 2}, {2, 1/3}}, 60], False},
        {pauliMatrix[1, 3], True},
        {pauliMatrix[1, 6], True},
        {spinMatrix[1, 4], Identity},
        {spinMatrix[2, 5], Identity}
    },
    ConstantArray[True, 10],
    {},
    TestID -> "Eigensystem-FirstEntryRuleWhereItAlreadyGivesCompactVectors"
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


BeginTestSection["QuantumOperator - irreducible shorthands"]

(* An expression no shorthand pass can reduce fails loudly instead of bouncing
   between FromOperatorShorthand's catch-all and the NumericFunction constructor
   rule until $RecursionLimit; the message shows operators as their labels. *)
VerificationTest[
    QuantumOperator[QuantumOperator["I"]^QuantumOperator["X"]],
    Failure["InvalidArguments", _],
    {QuantumOperator::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "Irreducible-operator-power-fails"
]

VerificationTest[
    QuantumOperator[QuantumOperator["Y"] QuantumOperator["Z"]],
    Failure["InvalidArguments", _],
    {QuantumOperator::invalidArgs},
    SameTest -> MatchQ,
    TestID -> "Irreducible-operator-times-fails"
]

(* A fully numeric scalar, atomic or compound, builds the scalar-eigenvalue
   diagonal operator. *)
VerificationTest[
    Normal[QuantumOperator[#]["Matrix"]] & /@ {Sqrt[2], Pi/2, 1 + Pi, E, E^2},
    # IdentityMatrix[2] & /@ {Sqrt[2], Pi/2, 1 + Pi, E, E^2},
    TestID -> "Numeric-scalar-diagonal-family"
]

(* A Failure produced inside a constructor chain propagates instead of being
   absorbed as a scalar eigenvalue by the diagonal catch-all. *)
VerificationTest[
    QuantumOperator[Failure["InvalidName", <||>]],
    Failure["InvalidName", <||>],
    TestID -> "Failure-passthrough"
]

(* Numeric constants are scalar coefficients in shorthands, never lifted to
   operators; a bare non-numeric symbol is the scalar-eigenvalue shorthand,
   directly and inside a product. *)
VerificationTest[
    Normal[QuantumOperator[Pi "PauliX"]["Matrix"]],
    {{0, Pi}, {Pi, 0}},
    TestID -> "Shorthand-constant-coefficient"
]

VerificationTest[
    Normal[QuantumOperator[2 \[FormalX]]["Matrix"]],
    {{2 \[FormalX], 0}, {0, 2 \[FormalX]}},
    TestID -> "Shorthand-symbol-scalar-diagonal"
]

(* A nonzero scalar base with an endomorphism exponent is the one-parameter
   group element base^A = Exp[Log[base] A]; an operator base with a scalar
   exponent stays MatrixPower; the unit base is Exp[0 A], the identity. *)
VerificationTest[
    {Simplify[Normal[(E^QuantumOperator["Z"])["Matrix"]]],
     Simplify[Normal[(2^QuantumOperator["Z"])["Matrix"]]]},
    {{{E, 0}, {0, 1/E}}, {{2, 0}, {0, 1/2}}},
    TestID -> "Power-scalar-base-matrix-exponential"
]

VerificationTest[
    FullSimplify[
        Normal[(\[FormalB]^QuantumOperator["PauliX"])["Matrix"]] -
        (Cosh[Log[\[FormalB]]] IdentityMatrix[2] + Sinh[Log[\[FormalB]]] PauliMatrix[1])
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Power-symbolic-base-closed-form"
]

(* The same symbolic base across Hilbert-space dimensions: b^diag(0..d-1) is
   diag(b^0, ..., b^(d-1)) for every d. *)
VerificationTest[
    Table[Simplify[Normal[(\[FormalB]^QuantumOperator[DiagonalMatrix[Range[0, d - 1]], d])["Matrix"]] - DiagonalMatrix[\[FormalB]^Range[0, d - 1]]], {d, 2, 6}],
    Table[ConstantArray[0, {d, d}], {d, 2, 6}],
    TestID -> "Power-symbolic-base-dimension-family"
]

VerificationTest[
    Normal[(QuantumOperator["X"]^2)["Matrix"]],
    {{1, 0}, {0, 1}},
    TestID -> "Power-operator-base-matrix-power"
]

VerificationTest[
    Normal[(1^QuantumOperator["Z"])["Matrix"]],
    {{1, 0}, {0, 1}},
    TestID -> "Power-unit-base-identity"
]

(* A zero base is the limit of b^op as b -> 0: Z has the eigenvalue -1, where b^-1
   diverges, so 0^Z fails rather than returning Z^0 = 1. A non-square exponent
   keeps the generic reading. A square exponent takes the scalar-base reading
   whatever frames it carries, so 2^Z is the exponential of the stored matrix,
   diag(2, 1/2), under any frame pair. *)
VerificationTest[
    MatchQ[0^QuantumOperator["Z"], Failure["ZeroBasePowerNoLimit", _]],
    True,
    TestID -> "Power-zero-base-limit"
]

VerificationTest[
    With[{mm = QuantumOperator[PauliMatrix[3], QuantumBasis[QuditBasis["PauliX"], QuditBasis[2]]]},
        Normal[(2^mm)["Matrix"]]
    ],
    {{2, 0}, {0, 1/2}},
    TestID -> "Power-square-cross-frame-stored-exponential"
]

VerificationTest[
    FailureQ[2^QuantumOperator["Cup"]],
    True,
    {MatrixPower::matsq},
    TestID -> "Power-nonsquare-exponent-generic-failure"
]

(* The exponential of I times a general Hermitian generator is unitary; the
   commuting family satisfies the group law; a non-commuting pair breaks it at
   exactly the commutator term of second order. *)
VerificationTest[
    With[{u = Normal[(E^(I QuantumOperator[\[FormalA] "PauliX" + \[FormalB] "PauliY" + \[FormalC] "PauliZ"]))["Matrix"]]},
        FullSimplify[ConjugateTranspose[u] . u, Element[\[FormalA] | \[FormalB] | \[FormalC], Reals]]
    ],
    IdentityMatrix[2],
    TestID -> "Power-exponential-unitary-general-generator"
]

VerificationTest[
    With[{
        ua = Normal[(E^(I \[FormalA] QuantumOperator["PauliX"]))["Matrix"]],
        ub = Normal[(E^(I \[FormalB] QuantumOperator["PauliX"]))["Matrix"]],
        uab = Normal[(E^(I (\[FormalA] + \[FormalB]) QuantumOperator["PauliX"]))["Matrix"]]
    },
        FullSimplify[ua . ub - uab]
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Power-exponential-group-law"
]

VerificationTest[
    With[{
        ea = Normal[(E^(I \[FormalE] QuantumOperator["PauliX"]))["Matrix"]],
        eb = Normal[(E^(I \[FormalE] QuantumOperator["PauliZ"]))["Matrix"]],
        eab = Normal[(E^(I \[FormalE] QuantumOperator["PauliX" + "PauliZ"]))["Matrix"]],
        px = Normal[PauliMatrix[1]], pz = Normal[PauliMatrix[3]]
    },
        Normal[Series[ea . eb - eab - (1/2) (I \[FormalE])^2 (px . pz - pz . px), {\[FormalE], 0, 2}]]
    ],
    ConstantArray[0, {2, 2}],
    TestID -> "Power-exponential-noncommuting-BCH"
]

(* A nilpotent generator is the case an eigendecomposition route gets wrong;
   MatrixExp handles the defective matrix exactly. *)
VerificationTest[
    Normal[(E^QuantumOperator[{{0, 1}, {0, 0}}])["Matrix"]],
    {{1, 1}, {0, 1}},
    TestID -> "Power-nilpotent-generator"
]

(* The SU(2) sign: a 2 Pi rotation is -1, only the 4 Pi rotation is the
   identity. An elementwise misreading of the exponential cannot produce it. *)
VerificationTest[
    {Simplify[Normal[(E^(I Pi QuantumOperator["PauliX"]))["Matrix"]]],
     Simplify[Normal[(E^(2 I Pi QuantumOperator["PauliX"]))["Matrix"]]]},
    {-IdentityMatrix[2], IdentityMatrix[2]},
    TestID -> "Power-exponential-su2-sign"
]

(* A non-diagonal qutrit generator with a degenerate spectrum containing zero:
   E^(-I theta Jy) at spin 1 is the Wigner small-d matrix. *)
VerificationTest[
    With[{jy = QuantumOperator[{{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2], 3]},
        FullSimplify[
            Normal[(E^(-I \[FormalT] jy))["Matrix"]] -
            {{(1 + Cos[\[FormalT]])/2, -Sin[\[FormalT]]/Sqrt[2], (1 - Cos[\[FormalT]])/2},
             {Sin[\[FormalT]]/Sqrt[2], Cos[\[FormalT]], -Sin[\[FormalT]]/Sqrt[2]},
             {(1 - Cos[\[FormalT]])/2, Sin[\[FormalT]]/Sqrt[2], (1 + Cos[\[FormalT]])/2}}
        ]
    ],
    ConstantArray[0, {3, 3}],
    TestID -> "Power-spin1-wigner-d"
]

(* Symbolic base and symbolic angle together: b^(t A) = Exp[t Log[b] A]. *)
VerificationTest[
    FullSimplify[
        Normal[(\[FormalB]^(\[FormalT] QuantumOperator["PauliX"]))["Matrix"]] -
        (Cosh[\[FormalT] Log[\[FormalB]]] IdentityMatrix[2] + Sinh[\[FormalT] Log[\[FormalB]]] PauliMatrix[1])
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Power-symbolic-base-and-angle"
]

(* A two-qubit operator on genuinely non-contiguous wires, checked against the
   involutive closed form b^A = ((b + 1/b)/2) Id + ((b - 1/b)/2) A for A^2 = Id,
   independent of MatrixExp. *)
VerificationTest[
    With[{c = QuantumOperator["CNOT", {3, 1}]},
        {Simplify[Normal[(2^c)["Matrix"]] - ((5/4) IdentityMatrix[4] + (3/4) Normal[c["Matrix"]])], (2^c)["Order"]}
    ],
    {ConstantArray[0, {4, 4}], {{1, 3}, {1, 3}}},
    TestID -> "Power-two-qubit-nonadjacent-order"
]

EndTestSection[]


BeginTestSection["QuantumOperator - scalar sum"]

(* x + qo adds x on the diagonal of the represented operator: x I + A. *)
VerificationTest[
    Normal[(2 + QuantumOperator["X"])["Matrix"]],
    2 IdentityMatrix[2] + PauliMatrix[1],
    TestID -> "ScalarSum-vector-type"
]

(* Plus is orderless, so qo + x agrees with x + qo. *)
VerificationTest[
    Normal[(QuantumOperator["X"] + 2)["Matrix"]] === Normal[(2 + QuantumOperator["X"])["Matrix"]],
    True,
    TestID -> "ScalarSum-orderless"
]

(* A density-matrix-type operator stores the d^2 x d^2 superoperator, so the added
   identity must be sized d^2.  Sizing it from the d-sized name dimensions raised a
   Thread::tdlen dimension error and returned a Failure. *)
VerificationTest[
    With[{xm = QuantumOperator["X"]["ToMatrix"]},
        {(2 + xm)["StateType"], Normal[(2 + xm)["Matrix"]] === 2 IdentityMatrix[4] + Normal[xm["Matrix"]]}
    ],
    {"Matrix", True},
    TestID -> "ScalarSum-matrix-type"
]

(* The same on two qubits: a 16 x 16 superoperator. *)
VerificationTest[
    With[{cm = QuantumOperator["CNOT"]["ToMatrix"]},
        Normal[(2 + cm)["Matrix"]] === 2 IdentityMatrix[16] + Normal[cm["Matrix"]]
    ],
    True,
    TestID -> "ScalarSum-matrix-type-two-qubit"
]

(* The scalar lands on the represented map, not the stored coefficients: on the
   non-orthonormal PauliX basis the running operator R represents X, and 2 + R
   represents 2 I + X. *)
VerificationTest[
    Simplify[Normal[(2 + QuantumOperator[{{1, 0}, {0, -1}}, "PauliX"])["MatrixRepresentation"]] - (2 IdentityMatrix[2] + PauliMatrix[1])],
    {{0, 0}, {0, 0}},
    TestID -> "ScalarSum-nonorthonormal-basis"
]

(* A symbolic scalar stays symbolic. *)
VerificationTest[
    Simplify[Normal[(a + QuantumOperator["X"])["Matrix"]] - (a IdentityMatrix[2] + PauliMatrix[1])],
    {{0, 0}, {0, 0}},
    TestID -> "ScalarSum-symbolic-scalar"
]

EndTestSection[]


BeginTestSection["QuantumOperator - orthonormal eigenbasis of a normal operator"]

(* Eigensystem returns an arbitrary basis inside the eigenspace of a repeated
   eigenvalue. That basis need not be orthogonal for exact input, above machine
   precision, or for a machine matrix that fails HermitianMatrixQ, which roundoff alone
   can make it fail. A normal operator has an orthonormal eigenbasis, and its projectors
   P_k onto it satisfy m.P_k = lambda_k P_k, P_k.P_k = P_k, Sum_k P_k = 1 and
   Sum_k lambda_k P_k = m, with lambda_k = Tr[m.P_k] the eigenvalue of P_k, and the
   basis that diagonalizes it is orthonormal. Most normal operators below are diagonal
   matrices of known spectrum written in a rotated basis, so these identities hold
   whatever basis the eigensolver picks. Exact residuals are decided exactly; inexact
   residuals of normal operators are compared to 10^-8, above the 10^-10 below which the
   eigensystem sets entries to zero. *)

(* diag(spectrum) in the Fourier basis *)
eigenbasisFourierRotated[spectrum_] := With[{f = FourierMatrix[Length[spectrum]]}, f . DiagonalMatrix[spectrum] . ConjugateTranspose[f]]

(* diag(spectrum) in the basis of H (x) ... (x) H, a rational matrix for qubits *)
eigenbasisHadamardRotated[spectrum_] := With[
    {h = KroneckerProduct @@ ConstantArray[{{1, 1}, {1, -1}} / Sqrt[2], Log2[Length[spectrum]]]},
    h . DiagonalMatrix[spectrum] . h
]

(* diag(spectrum) in a random unitary basis of the given seed, the Gram-Schmidt
   orthonormalization of random complex vectors, at the precision of the spectrum *)
eigenbasisRandomRotated[spectrum_, seed_] := With[
    {u = BlockRandom[
        Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {Length[spectrum], Length[spectrum]}, WorkingPrecision -> Precision[spectrum]]],
        RandomSeeding -> seed
    ]},
    ConjugateTranspose[u] . DiagonalMatrix[spectrum] . u
]

(* the entries of m.P_k - lambda_k P_k, P_k.P_k - P_k, Sum_k P_k - 1 and
   Sum_k lambda_k P_k - m for the projectors P_k of qo, each with its own eigenvalue
   lambda_k = Tr[m.P_k]: a separate call of qo["Eigenvalues"] can list the eigenvalues
   of an exact matrix in another order *)
eigenbasisProjectorResiduals[qo_] := With[
    {m = Normal[qo["MatrixRepresentation"]], p = Normal /@ qo["Projectors"]},
    {values = Tr[m . #] & /@ p},
    Flatten[{
        MapThread[m . #2 - #1 #2 &, {values, p}],
        (# . # - #) & /@ p,
        Total[p] - IdentityMatrix[Length[m]],
        Total[values p] - m
    }]
]

(* the entries of v^*.v^T - 1 for the basis vectors v of qo["Diagonalize"], and of the
   operator it represents minus qo *)
eigenbasisDiagonalizeResiduals[qo_] := With[
    {d = qo["Diagonalize"]},
    {v = Normal /@ d["Basis"]["Output"]["Elements"]},
    Flatten[{
        Conjugate[v] . Transpose[v] - IdentityMatrix[Length[v]],
        Normal[d["MatrixRepresentation"]] - Normal[qo["MatrixRepresentation"]]
    }]
]

eigenbasisExactZeroQ[residuals_] := AllTrue[residuals, PossibleZeroQ[#, Method -> "ExactAlgebraics"] &]

VerificationTest[
    eigenbasisExactZeroQ @ eigenbasisProjectorResiduals[QuantumOperator[eigenbasisFourierRotated[{1, 1, -1, -1}], {1, 2}]],
    True,
    TestID -> "Projectors-ExactRotatedDegenerate"
]

(* three qubits with three repeated eigenvalues: 1 three times, -1 and 2 twice each *)
VerificationTest[
    eigenbasisExactZeroQ @ eigenbasisProjectorResiduals[QuantumOperator[eigenbasisHadamardRotated[{1, 1, 1, -1, -1, 2, 2, 0}], {1, 2, 3}]],
    True,
    TestID -> "Projectors-ExactRotatedThreeRepeatedEigenvalues"
]

(* The same spectrum in the Fourier basis of three qubits, whose entries lie in
   Q(i, Sqrt[2]), where exact expressions for the eigenvectors can grow large, so the
   time constraint is part of the test. The projectors of each repeated eigenvalue,
   identified by their eigenvalue Tr[m.P_k], add up to its spectral projector
   F.diag(1 on that eigenvalue).F^†. *)
VerificationTest[
    With[{f = FourierMatrix[8], spectrum = {1, 1, 1, -1, -1, 2, 2, 0}},
        {m = f . DiagonalMatrix[spectrum] . ConjugateTranspose[f]},
        {p = Normal /@ QuantumOperator[m, {1, 2, 3}]["Projectors"]},
        {values = RootReduce[Tr[m . #]] & /@ p},
        eigenbasisExactZeroQ @ Flatten[Map[
            Total[Pick[p, values, #]] - f . DiagonalMatrix[Boole[Thread[spectrum == #]]] . ConjugateTranspose[f] &,
            {1, -1, 2}
        ]]
    ],
    True,
    TimeConstraint -> 120,
    TestID -> "Projectors-ExactAlgebraicEntries"
]

(* The spectrum 1 + e, 1 - e, -1 + e, -1 - e has no repeated eigenvalue for 0 < e < 1,
   where the eigenvectors of a normal matrix are orthogonal whatever the solver; the
   projectors must stay right as e reaches 0 and the eigenvalues pair up. *)
VerificationTest[
    Map[
        eigenbasisExactZeroQ @ eigenbasisProjectorResiduals[QuantumOperator[eigenbasisFourierRotated[{1 + #, 1 - #, -1 + #, -1 - #}], {1, 2}]] &,
        {1/10, 10^-6, 0}
    ],
    {True, True, True},
    TestID -> "Projectors-ExactSplitPairLimit"
]

(* The reflection 1 - 2 u.u^†/(u^†.u) with u = (1, r, r^2, r^3), for the negative root r
   of x^4 = x + 1: the eigenvalue -1 on u and 1 on its three-dimensional complement,
   with entries in the field of r. Two calls of Eigensystem on this matrix can list its
   eigenvalues in different orders. *)
VerificationTest[
    With[{r = Root[#^4 - # - 1 &, 1]},
        {u = {1, r, r^2, r^3}},
        {qo = QuantumOperator[IdentityMatrix[4] - 2 Outer[Times, u, u] / (u . u), {1, 2}]},
        eigenbasisExactZeroQ @ Join[eigenbasisProjectorResiduals[qo], eigenbasisDiagonalizeResiduals[qo]]
    ],
    True,
    TestID -> "Eigenbasis-ExactQuarticFieldReflection"
]

(* Degenerate perturbation theory: e V, with V = |f1><f3| + |f3><f1| in the Fourier basis
   f1, ..., f4, couples the eigenvalue-1 vector f1 to the eigenvalue -1 vector f3 of
   F.diag(1, 1, -1, -1).F^†. For e != 0 the eigenvalues are 1, -1 and +-Sqrt[1 + e^2], so
   the projectors depend on e. As e -> 0 each tends to the projector onto one Fourier
   vector: the perturbation picks that basis of the eigenspaces, whatever basis the
   eigensolver takes at e = 0. The projectors whose eigenvalue tends to 1 add up to the
   spectral projector F.diag(1, 1, 0, 0).F^† of the eigenvalue 1 at e = 0. *)
VerificationTest[
    With[{f = FourierMatrix[4], v = Outer[Times, UnitVector[4, 1], UnitVector[4, 3]] + Outer[Times, UnitVector[4, 3], UnitVector[4, 1]]},
        {m = f . (DiagonalMatrix[{1, 1, -1, -1}] + e v) . ConjugateTranspose[f]},
        {pe = Normal /@ QuantumOperator[m, {1, 2}]["Projectors"]},
        {
            FreeQ[pe, e],
            Sort @ Map[
                Function[p, SelectFirst[Range[4], Simplify[Limit[p, e -> 0] - f . DiagonalMatrix[UnitVector[4, #]] . ConjugateTranspose[f]] === ConstantArray[0, {4, 4}] &]],
                pe
            ],
            Simplify[Limit[Total[Select[pe, Limit[Tr[m . #], e -> 0] === 1 &]], e -> 0] - f . DiagonalMatrix[{1, 1, 0, 0}] . ConjugateTranspose[f]]
        }
    ],
    {False, {1, 2, 3, 4}, ConstantArray[0, {4, 4}]},
    TestID -> "Projectors-ExactDegeneratePerturbationLimit"
]

(* The quantum Fourier transform on three qubits, a unitary whose eigenvalues 1, i, -1
   and -i repeat 3, 2, 2 and 1 times *)
VerificationTest[
    With[{qo = QuantumOperator[FourierMatrix[8], {1, 2, 3}]},
        eigenbasisExactZeroQ @ Join[eigenbasisProjectorResiduals[qo], eigenbasisDiagonalizeResiduals[qo]]
    ],
    True,
    TestID -> "Projectors-ExactQuantumFourierTransform"
]

(* The same on four qubits in machine numbers. A unitary is not Hermitian, so Eigensystem
   takes its general solver, whose eigenvectors for a repeated eigenvalue are far from
   orthogonal. *)
VerificationTest[
    With[{f = N[FourierMatrix[16]]},
        {
            Max[Abs[Conjugate[#] . Transpose[#] - IdentityMatrix[16]]] & @ Last[Eigensystem[f]],
            Max[Abs[Join[eigenbasisProjectorResiduals[#], eigenbasisDiagonalizeResiduals[#]]]] & @ QuantumOperator[f, Range[4]]
        }
    ],
    {_ ? (# > 10^-1 &), _ ? (# < 10^-8 &)},
    SameTest -> MatchQ,
    TestID -> "Projectors-MachineQuantumFourierTransform"
]

(* The Heisenberg ring H = Sum_i sigma_i . sigma_(i+1) of n spins, whose SU(2) and
   translation symmetries leave its levels degenerate: for n = 4 they repeat 1, 3, 7 and
   5 times, and for n = 5 the eigenvectors of some lie in Q(Sqrt[5]). *)
eigenbasisHeisenbergRing[n_] := Sum[
    KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], n], {i -> s, Mod[i, n] + 1 -> s}],
    {i, n}, {s, PauliMatrix /@ {1, 2, 3}}
]

VerificationTest[
    eigenbasisExactZeroQ @ eigenbasisProjectorResiduals[QuantumOperator[eigenbasisHeisenbergRing[#], Range[#]]] & /@ {4, 5},
    {True, True},
    TestID -> "Projectors-ExactHeisenbergRing"
]

(* Machine precision: this rotation fails HermitianMatrixQ from roundoff alone, so
   Eigensystem takes its general solver, whose eigenvectors for a repeated eigenvalue
   are not orthogonal. *)
VerificationTest[
    With[{m = eigenbasisRandomRotated[N @ {1, 1, -1, -1}, 35]},
        {HermitianMatrixQ[m], Max[Abs[eigenbasisProjectorResiduals[QuantumOperator[m, {1, 2}]]]]}
    ],
    {False, _ ? (# < 10^-8 &)},
    SameTest -> MatchQ,
    TestID -> "Projectors-MachineRotatedDegenerate"
]

(* 60 digits: above machine precision Eigensystem normalizes the eigenvectors of a
   repeated eigenvalue but does not orthogonalize them, Hermitian input or not. *)
VerificationTest[
    Max[Abs[eigenbasisProjectorResiduals[QuantumOperator[eigenbasisRandomRotated[N[{1/2, 1/2, 0, 0}, 60], 3], {1, 2}]]]],
    _ ? (# < 10^-8 &),
    SameTest -> MatchQ,
    TestID -> "Projectors-60DigitRotatedDegenerate"
]

VerificationTest[
    eigenbasisExactZeroQ @ eigenbasisDiagonalizeResiduals[QuantumOperator[eigenbasisFourierRotated[{1, 1, -1, -1}], {1, 2}]],
    True,
    TestID -> "Diagonalize-ExactRotatedDegenerate"
]

VerificationTest[
    With[{m = eigenbasisRandomRotated[N @ {1, 1, -1, -1}, 35]},
        {HermitianMatrixQ[m], Max[Abs[eigenbasisDiagonalizeResiduals[QuantumOperator[m, {1, 2}]]]]}
    ],
    {False, _ ? (# < 10^-8 &)},
    SameTest -> MatchQ,
    TestID -> "Diagonalize-MachineRotatedDegenerate"
]

(* A non-normal operator has no orthonormal eigenbasis, and orthogonalizing its
   eigenvectors would give vectors that are not eigenvectors. When it can be
   diagonalized, its projectors stay onto eigenvectors, m.P_k = lambda_k P_k, and its
   diagonalization still represents it, also with a repeated eigenvalue:
   s.diag(1, 1, 2).s^-1 for a non-unitary s. *)
(* the entries of m.P_k - Tr[m.P_k] P_k for the projectors P_k of qo, and of the operator
   that qo["Diagonalize"] represents minus qo: what still holds when the eigenvectors of
   qo are kept *)
eigenbasisKeptResiduals[qo_] := With[{m = Normal[qo["MatrixRepresentation"]]},
    Flatten[{(m . # - Tr[m . #] #) & /@ (Normal /@ qo["Projectors"]), Normal[qo["Diagonalize"]["MatrixRepresentation"]] - m}]
]

VerificationTest[
    eigenbasisExactZeroQ @ eigenbasisKeptResiduals[QuantumOperator[#]] & /@ {
        {{1, 1}, {0, 2}},
        With[{s = {{1, 1, 0}, {0, 1, 1}, {1, 0, 2}}}, s . DiagonalMatrix[{1, 1, 2}] . Inverse[s]]
    },
    {True, True},
    TestID -> "Eigenbasis-NonNormalKeepsEigenvectors"
]

(* Non-normal by 10^-7 with eigenvalues 10^-8 apart: m.m^† - m^†.m is only of order
   10^-14 there, but the eigenvectors (1, 0) and nearly (0.995, 0.0995) are far from
   orthogonal, and orthogonalizing them would leave a residual of order 10^-7. *)
VerificationTest[
    Max[Abs[eigenbasisKeptResiduals[QuantumOperator[{{1., 1.*^-7}, {0., 1.00000001}}]]]],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TestID -> "Eigenbasis-NearlyNonNormalKeepsEigenvectors"
]

(* Adding a multiple of the identity changes neither the eigenvectors nor the departure
   from normality, so it must not change the projectors either: m, non-normal by 10^-6
   with eigenvalues 10^-3 apart, keeps its eigenvectors with 100, 1000 or 10^4 added. *)
VerificationTest[
    With[{m = {{0, 1.*^-6}, {0, 1.*^-3}}},
        {sorted = SortBy[Normal /@ QuantumOperator[#]["Projectors"], Re[Tr[m . #]] &] &},
        Max[Abs[Flatten[{
            eigenbasisKeptResiduals[QuantumOperator[m]],
            sorted[m + # IdentityMatrix[2]] - sorted[m] & /@ {100, 1000, 10^4}
        }]]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TestID -> "Eigenbasis-ShiftKeepsEigenvectors"
]

(* The propagator exp(-i h dt) over a short time dt shares the eigenvectors of its
   generator h = X - (i/2)|0><0|, which is not normal: its eigenvectors overlap. However
   close to the identity the propagator is, its departure from normality, of order dt,
   is far above roundoff, so its projectors must be those of h, to the accuracy with
   which roundoff lets the eigensolver resolve them, which falls as dt does. *)
VerificationTest[
    With[{h = {{-I/2, 1}, {1, 0}}},
        {byEigenvalue = SortBy[#, Re[Tr[h . #]] &] &},
        {hp = byEigenvalue[Outer[Times, #, Conjugate[#]] / (Conjugate[#] . #) & /@ Eigenvectors[h]]},
        Max[Abs[Flatten[byEigenvalue[Normal /@ QuantumOperator[MatrixExp[-I h #]]["Projectors"]] - hp & /@ {1.*^-8, 1.*^-9, 1.*^-10}]]]
    ],
    _ ? (# < 10^-4 &),
    SameTest -> MatchQ,
    TestID -> "Eigenbasis-ShortTimePropagatorKeepsEigenvectors"
]

(* The same s.diag(1, 1, 2).s^-1 computed in 60-digit arithmetic. The entries that cancel
   become zeros known only to some accuracy, of precision 0, so the matrix has precision
   0; that must not make it pass for normal. *)
VerificationTest[
    With[{s = N[{{1, 1, 0}, {0, 1, 1}, {1, 0, 2}}, 60]},
        {m = s . DiagonalMatrix[N[{1, 1, 2}, 60]] . Inverse[s]},
        {qo = QuantumOperator[m]},
        {Precision[m], Max[Abs[eigenbasisKeptResiduals[qo]]]}
    ],
    {0., _ ? (# < 10^-40 &)},
    SameTest -> MatchQ,
    TestID -> "Eigenbasis-AccuracyOnlyZeroKeepsEigenvectors"
]

(* The evolution operator exp(-i h t) of the machine matrix h of the test
   Projectors-MachineRotatedDegenerate, at a long time t, is unitary, but roundoff in the
   matrix exponential leaves it further from normal than one decomposition leaves; it
   must still get a spectral decomposition. *)
VerificationTest[
    Max[Abs[eigenbasisProjectorResiduals[QuantumOperator[MatrixExp[-I 1000.7 eigenbasisRandomRotated[N @ {1, 1, -1, -1}, 35]], {1, 2}]]]],
    _ ? (# < 10^-8 &),
    SameTest -> MatchQ,
    TestID -> "Projectors-LongTimeUnitary"
]

(* Eigenvectors that are orthonormal already come back unchanged, bit for bit: a machine
   operator whose eigenvalues are all distinct, and a Hermitian one whose eigensolver
   returns an orthonormal basis. *)
VerificationTest[
    Map[
        With[{qo = #},
            {
                Normal /@ qo["Projectors"] === Normal /@ qo["Projectors", "Orthogonalize" -> False],
                qo["Diagonalize"] === qo["Diagonalize", "Orthogonalize" -> False]
            }
        ] &,
        {
            QuantumOperator["RX"[0.3]],
            QuantumOperator[eigenbasisRandomRotated[N @ {1, 2, 3, 4}, 7], {1, 2}],
            QuantumOperator[BlockRandom[Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {4, 4}]], RandomSeeding -> 11], {1, 2}]
        }
    ],
    ConstantArray[{True, True}, 3],
    TestID -> "Eigenbasis-OrthonormalEigenvectorsUnchanged"
]

(* Symbolic eigenvalues a, a, b, b with numeric eigenvectors: the eigenspaces are those of
   the distinct symbols, and their eigenvectors are orthonormalized as for numbers. *)
VerificationTest[
    With[{qo = QuantumOperator[eigenbasisFourierRotated[{a, a, b, b}], {1, 2}]},
        {v = Normal /@ qo["Diagonalize"]["Basis"]["Output"]["Elements"]},
        Union @ Flatten @ Simplify[{eigenbasisProjectorResiduals[qo], Conjugate[v] . Transpose[v] - IdentityMatrix[4]}]
    ],
    {0},
    TestID -> "Eigenbasis-SymbolicEigenvaluesNumericEigenvectors"
]

EndTestSection[]
