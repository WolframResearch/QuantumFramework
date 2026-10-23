Begin["Wolfram`QuantumFramework`GateMatrix`PackagePrivate`"]
stabilizerOperatorGroupSpecs[operators_List, data_, positions_List] := With[{
    orders = operators[[All, 2]], first = First[operators]
},
    With[{local = If[AllTrue[orders, Sort[First[#]] === Sort[Last[#]] &],
        Replace[Replace[data, $stabilizerPrimitiveDispatch], $Failed :>
            With[{key = stabilizerGateKey[first]},
                Lookup[$stabilizerGateMemo, key,
                    stabilizerRemember[$stabilizerGateMemo, key, stabilizerOperatorClassification[first]]]
            ]],
        $Failed
    ]},
        MapThread[List, {positions,
            If[ListQ[local],
                stabilizerRemapSpecs[local, Sort /@ orders[[All, 2]]],
                (* Symbolic/composite token recovery and refusal details may
                   depend on each occurrence's label and absolute wires. *)
                stabilizerOperatorSpecs /@ operators
            ]
        }]
    ]
]

End[]
