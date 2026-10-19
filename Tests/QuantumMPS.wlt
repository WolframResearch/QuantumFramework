BeginTestSection["QuantumMPSApply"]

(* QuantumMPSApply acts with a gate on the MPS circuit QuantumMPS builds for a state, and
   returns an MPS circuit again. Every test compares the contracted MPS, amplitude by amplitude
   (phase included), with the gate applied to the dense state. The gate's wires may be adjacent
   or not, in either order, and may number more than two; the sites are read from the flattened
   circuit, so the nesting QuantumMPS adds does not matter. *)

Needs["Wolfram`QuantumFramework`Experimental`"];

mpsAmplitudes[mps_] := Normal[mps[]["StateVector"]];
mpsMatchesDenseQ[qs_, op_] := With[{mps = QuantumMPSApply[QuantumMPS[qs], op]},
    MatchQ[mps, _QuantumCircuitOperator] && TrueQ[Chop[Norm[N[mpsAmplitudes[mps] - Normal[op[qs]["StateVector"]]]]] == 0]
];

mpsRandom4 = BlockRandom[QuantumState["RandomPure"[4]], RandomSeeding -> 7];
mpsRandom5 = BlockRandom[QuantumState["RandomPure"[5]], RandomSeeding -> 7];
mpsRandomGate = BlockRandom[QuantumOperator["RandomUnitary"[{2, 2}], {2, 5}], RandomSeeding -> 11];
mpsRandomGateReversed = BlockRandom[QuantumOperator["RandomUnitary"[{2, 2}], {3, 1}], RandomSeeding -> 13];

(* ---- a CNOT on product, GHZ and random states ---- *)

VerificationTest[mpsMatchesDenseQ[QuantumState["1000"], QuantumOperator["CNOT", {1, 2}]], True, TestID -> "mps-cnot-product-state"]

VerificationTest[mpsMatchesDenseQ[QuantumState["GHZ"[4]], QuantumOperator["CNOT", {1, 2}]], True, TestID -> "mps-cnot-ghz-sites-12"]

VerificationTest[mpsMatchesDenseQ[QuantumState["GHZ"[4]], QuantumOperator["CNOT", {2, 3}]], True, TestID -> "mps-cnot-ghz-sites-23"]

VerificationTest[mpsMatchesDenseQ[mpsRandom4, QuantumOperator["CNOT", {2, 3}]], True, TestID -> "mps-cnot-random-sites-23"]

(* ---- qubit order: the control stays on the operator's first wire ---- *)

VerificationTest[mpsMatchesDenseQ[mpsRandom4, QuantumOperator["CNOT", {3, 2}]], True, TestID -> "mps-cnot-reversed-order"]

(* ---- non-adjacent sites: the sites in between are carried through the update ---- *)

VerificationTest[mpsMatchesDenseQ[mpsRandom4, QuantumOperator["CNOT", {1, 4}]], True, TestID -> "mps-cnot-nonadjacent"]

VerificationTest[mpsMatchesDenseQ[mpsRandom4, QuantumOperator["CNOT", {4, 2}]], True, TestID -> "mps-cnot-nonadjacent-reversed"]

VerificationTest[mpsMatchesDenseQ[mpsRandom5, mpsRandomGate], True, TestID -> "mps-random-two-qubit-unitary-nonadjacent"]

VerificationTest[mpsMatchesDenseQ[mpsRandom5, QuantumOperator["Toffoli", {5, 1, 3}]], True, TestID -> "mps-three-qubit-gate"]

(* ---- a one-qubit gate keeps matching the dense state ---- *)

VerificationTest[mpsMatchesDenseQ[mpsRandom4, QuantumOperator["H", {2}]], True, TestID -> "mps-one-qubit-gate"]

(* ---- the result is an MPS again, so gates chain ---- *)

VerificationTest[
    With[{gates = {
        QuantumOperator["H", {1}], QuantumOperator["CNOT", {1, 2}], QuantumOperator["CNOT", {2, 4}],
        mpsRandomGateReversed, QuantumOperator["CZ", {4, 3}]
    }},
        TrueQ[Chop[Norm[N[
            mpsAmplitudes[Fold[QuantumMPSApply, QuantumMPS[QuantumState["0000"]], gates]] -
                Normal[Fold[#2[#1] &, QuantumState["0000"], gates]["StateVector"]]
        ]]] == 0]
    ],
    True,
    TestID -> "mps-gates-chain"
]

(* ---- the split back into sites: a CNOT on |+000> leaves one bond, of dimension 2, between
   sites 1 and 2 (the Bell pair) and no bond elsewhere (the |00> factor) ---- *)

VerificationTest[
    With[{ops = QuantumMPSApply[QuantumMPS[QuantumState["+000"]], QuantumOperator["CNOT", {1, 2}]]["Operators"]},
        {#["Order"], #["Dimensions"]} & /@ ops
    ],
    {
        {{{0, 1}, {}}, {2, 2}},
        {{{2}, {0}}, {2, 2}},
        {{{3}, {}}, {2}},
        {{{4}, {}}, {2}}
    },
    TestID -> "mps-bell-pair-bond"
]

(* exact input stays exact through the split, with amplitudes about as simple as the dense ones:
   the exact split is a rank factorization by row reduction; an exact singular value
   decomposition builds algebraic singular values that swell the amplitudes *)
VerificationTest[
    With[{amps = mpsAmplitudes[QuantumMPSApply[QuantumMPS[QuantumState["+000"]], QuantumOperator["CNOT", {1, 3}]]]},
        {FreeQ[amps, _Real], Simplify[amps - Normal[QuantumOperator["CNOT", {1, 3}][QuantumState["+000"]]["StateVector"]]]}
    ],
    {True, ConstantArray[0, 16]},
    TestID -> "mps-exact-input-stays-exact"
]

VerificationTest[
    With[{qs = QuantumState["W"[3]], op = QuantumOperator["CNOT", {1, 3}]},
        With[{amps = mpsAmplitudes[QuantumMPSApply[QuantumMPS[qs], op]], dense = Normal[op[qs]["StateVector"]]},
            {FreeQ[amps, _Real], Simplify[amps - dense] === ConstantArray[0, 8], LeafCount[amps] <= 2 LeafCount[dense]}
        ]
    ],
    {True, True, True},
    TestID -> "mps-exact-w-state-stays-simple",
    TimeConstraint -> 60
]

(* QuantumMPS of a generic exact state carries nested radicals, which exact row reduction swells
   without bound; those sites are split numerically, with a message, and still match *)
VerificationTest[
    mpsMatchesDenseQ[QuantumState[Normalize[Range[8]]], QuantumOperator["CNOT", {1, 3}]],
    True,
    {QuantumMPSApply::numsplit},
    TestID -> "mps-nested-radicals-split-numerically",
    TimeConstraint -> 60
]

(* ---- symbolic gates: a one-site gate needs no split and stays symbolic; a gate across sites
   needs a singular value decomposition, so a symbolic angle there is refused ---- *)

VerificationTest[
    Simplify[
        mpsAmplitudes[QuantumMPSApply[QuantumMPS[QuantumState["GHZ"[3]]], QuantumOperator["RX"[\[Theta]], {2}]]] -
            Normal[QuantumOperator["RX"[\[Theta]], {2}][QuantumState["GHZ"[3]]]["StateVector"]]
    ],
    ConstantArray[0, 8],
    TestID -> "mps-symbolic-one-site-gate"
]

VerificationTest[
    QuantumMPSApply[QuantumMPS[QuantumState["+000"]], QuantumOperator["R"[\[Theta], "ZZ"], {1, 2}]]["Message"],
    "Splitting the updated sites needs numeric tensors",
    TestID -> "mps-symbolic-two-site-gate-refused"
]

(* ---- a wire of the gate that carries no MPS site starts in |0>, as a circuit wire does, and
   the result is still an MPS that takes further gates ---- *)

VerificationTest[
    With[{mps = QuantumMPSApply[QuantumMPS[QuantumState["+0"]], QuantumOperator["CNOT", {1, 3}]]},
        {
            Chop[Norm[N[mpsAmplitudes[mps] - Normal[QuantumOperator["CNOT", {1, 3}][QuantumState["+00"]]["StateVector"]]]]],
            Chop[Norm[N[
                mpsAmplitudes[QuantumMPSApply[mps, QuantumOperator["X", {4}]]] -
                    Normal[QuantumOperator["X", {4}][QuantumOperator["CNOT", {1, 3}][QuantumState["+000"]]]["StateVector"]]
            ]]]
        }
    ],
    {0, 0},
    TestID -> "mps-new-wire-starts-in-zero"
]

EndTestSection[]
