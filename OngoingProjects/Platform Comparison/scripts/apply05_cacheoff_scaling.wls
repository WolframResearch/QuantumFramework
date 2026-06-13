(* Separate cache pathology from algorithmic superlinearity:
   depth sweep with property cache OFF, contraction-only scaling, sparse-vs-packed micro. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`TensorNetworks`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
ms[t_] := Round[1000. t, 0.1];

resolve[name_] := With[{hits = Names["*`" <> name]}, If[hits === {}, $Failed, Symbol[First[hits]]]];
TNCompile = resolve["TensorNetworkCompile"];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_, L_ : 3] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {L}], 1]];
reg[n_] := QuantumState["Register"[n]];

(* warm up the standard circuit *)
n = 12; psi = reg[n]; qc = circ[n, 0.3];
qc[psi];

(* ---------- depth sweep, cache OFF ---------- *)
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
Do[
    Module[{qcL = circ[12, 0.3 + L 10.^-6, L], t1, t2},
        t1 = First @ AbsoluteTiming[qcL[psi]];
        t2 = First @ AbsoluteTiming[qcL[psi]];
        pr["depthOFF_L" <> ToString[L] <> "_gates" <> ToString[qcL["GateCount"]] <> "_run1_run2_ms", {ms[t1], ms[t2]}]
    ],
    {L, {1, 2, 4, 8, 16, 24}}
]
(* width sweep, cache OFF *)
Do[
    Module[{qcW = circ[w, 0.3 + w 10.^-6], psiW = reg[w], t1, t2},
        t1 = First @ AbsoluteTiming[qcW[psiW]];
        t2 = First @ AbsoluteTiming[qcW[psiW]];
        pr["widthOFF_n" <> ToString[w] <> "_gates" <> ToString[qcW["GateCount"]] <> "_run1_run2_ms", {ms[t1], ms[t2]}]
    ],
    {w, {8, 12, 16}}
]
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = True;

(* ---------- contraction-only scaling on the real net (state included) ---------- *)
Do[
    Module[{w = ww, qcW, psiW, circuit, net, t},
        qcW = circ[w, 0.3]; psiW = reg[w];
        circuit = QuantumCircuitOperator[{psiW -> qcW["Flatten"]["Sort"]["FullInputOrder"]} /* qcW["Flatten"]["Sort"], "Label" -> None];
        net = TNCompile[circuit, "ReturnTensorNetwork" -> True];
        Do[
            ClearSystemCache[];
            t = First @ AbsoluteTiming[TensorNetworkContract[net, Automatic]];
            pr["contract_n" <> ToString[w] <> "_run" <> ToString[k] <> "_ms", ms[t]],
            {k, 2}
        ]
    ],
    {ww, {8, 12, 16}}
]

(* ---------- SparseArray-as-dense-vector overhead micro ---------- *)
v = RandomComplex[{-1 - I, 1 + I}, 2^12];
sv = SparseArray[v];
m = RandomComplex[{-1 - I, 1 + I}, {16, 16}];
pv = Developer`ToPackedArray[v];
pr["dense_vec_bytes", ByteCount[pv]];
pr["sparse_vec_bytes", ByteCount[sv]];
t = First[AbsoluteTiming[Do[pv . pv, {100}]]]/100; pr["packed_dot_us", Round[10.^6 t, 0.1]];
t = First[AbsoluteTiming[Do[sv . sv, {100}]]]/100; pr["sparse_dot_us", Round[10.^6 t, 0.1]];

(* ---------- per-dispatch arithmetic: how many dispatches fit in the fast-path budget ---------- *)
h = QuantumOperator["H", {1}]; h["Order"];
t = First[AbsoluteTiming[Do[h["Order"], {1000}]]]/1000;
pr["dispatch_warm_us", Round[10.^6 t, 0.1]];
