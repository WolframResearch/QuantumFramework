(* A controlled gate placed on new wires moves its output and input wires with one map, so
   it stays an operator on one set of wires. Given every wire, QuantumOperator[op, order]
   sends the sorted wires of op to order, as for any operator; given an {output, input}
   pair or fewer wires, a controlled gate reads the order with its controls first, sorted,
   then its targets.

   References are built by hand: P1 (x) X (x) P0 + (1 - P1 (x) 1 (x) P0) for X on the middle
   wire, controlled on 1 by the first wire and on 0 by the last, and the gate built directly
   on its new wires.

   Regimes covered:
     exactly solvable     C_X controlled on 1 and on 0, moved to {5, 7, 6}, against the
                          matrix built by hand
     general              controlled gates whose stored order does not list the controls
                          first (C_X[{1}, {2}], C_RX(theta)[{1}, {3}], Toffoli on {3, 1, 2}),
                          moved to every permutation of three new wires, against the gate
                          built directly there
     same wires           output and input orders agree after every move, for full lists,
                          {output, input} pairs and single wires
     unchanged            a controlled gate stored controls first moves as before *)

sortedMatrix[qo_] := Normal[qo["Sort"]["Matrix"]]

p0 = {{1, 0}, {0, 0}};
p1 = {{0, 0}, {0, 1}};

VerificationTest[
    With[{qo = QuantumOperator[QuantumOperator["C"["X", {1}, {2}]], {5, 7, 6}]},
        {
            qo["Order"],
            sortedMatrix[qo] ===
                KroneckerProduct[p1, PauliMatrix[1], p0] + IdentityMatrix[8] - KroneckerProduct[p1, IdentityMatrix[2], p0]
        }
    ],
    {{{7, 5, 6}, {7, 5, 6}}, True},
    TestID -> "ControlledReorder-CX-matrix"
]

(* each gate is {op, target, control-on-1 wires, control-on-0 wires}; moving it with the map
   from its sorted wires to new wires gives the gate built on the mapped wires *)
$gates = {
    {"X", 3, {1}, {2}},
    {"RX"[\[Theta]], 2, {1}, {3}},
    {"NOT", 2, {3, 1}, {}}
};

build[{op_, t_, c1_, c0_}] := QuantumOperator["C"[op -> {t}, c1, c0]]

Scan[
    Function[{g},
        Scan[
            Function[{wires},
                VerificationTest[
                    With[{
                        moved = QuantumOperator[build[g], wires],
                        repl = Thread[Sort[Flatten[Rest[g]]] -> wires]
                    },
                        With[{direct = build[MapAt[# /. repl &, g, {{2}, {3}, {4}}]]},
                            {Sort[moved["OutputOrder"]] === Sort[moved["InputOrder"]], moved["OutputOrder"] === moved["InputOrder"], sortedMatrix[moved] === sortedMatrix[direct]}
                        ]
                    ],
                    {True, True, True},
                    TestID -> "ControlledReorder-" <> ToString[g[[1]], InputForm] <> "-" <> ToString[wires]
                ]
            ],
            Permutations[{4, 6, 9}]
        ]
    ],
    $gates
]

VerificationTest[
    Map[
        Function[qo, Map[#["OutputOrder"] === #["InputOrder"] &, {QuantumOperator[qo, qo["Order"]], QuantumOperator[qo, {7}], QuantumOperator[qo, {{8, 9, 7}, {8, 9, 7}}]}]],
        build /@ $gates
    ],
    ConstantArray[{True, True, True}, 3],
    TestID -> "ControlledReorder-pairs-and-single-wires-same-order"
]

(* read with the controls first, an {output, input} pair moves the sorted controls to its
   first wires and the target to the last *)
VerificationTest[
    With[{qo = QuantumOperator[QuantumOperator["Toffoli", {3, 1, 2}], {{7, 8, 9}, {7, 8, 9}}]},
        {qo["OutputOrder"] === qo["InputOrder"], sortedMatrix[qo] === sortedMatrix[QuantumOperator["Toffoli", {7, 8, 9}]]}
    ],
    {True, True},
    TestID -> "ControlledReorder-pair-controls-first"
]

VerificationTest[
    With[{qo = QuantumOperator[QuantumOperator["CNOT", {2, 1}], {5, 4}]}, {qo["Order"], qo["Label"]}],
    {{{4, 5}, {4, 5}}, Subscript["C", "NOT"][{4}, {}]},
    TestID -> "ControlledReorder-controls-first-unchanged"
]
