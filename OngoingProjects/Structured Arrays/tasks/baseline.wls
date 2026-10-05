(* Baseline for the structured-array task briefs, on the QF checkout:
     wolframscript -file baseline.wls group
   group: basis | spectrum | power | liouvillian | permutation | qft *)
group = Last[$ScriptCommandLine];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, group, "\t", label, "\t", ToString[value, InputForm, PageWidth -> Infinity], "\n"];
SetAttributes[{timed, msgs}, HoldFirst];
(* one run after ClearSystemCache[], within a time limit: {seconds or "TIMEOUT>limit", value} *)
timed[e_, limit_ : 120] := (ClearSystemCache[]; Replace[
    AbsoluteTiming[TimeConstrained[e, limit, $timedOut]],
    {{_, $timedOut} :> {"TIMEOUT>" <> ToString[limit], Missing["TimedOut"]}, {t_, v_} :> {Round[t, 0.0001], v}}
]);
(* the messages an evaluation issued, by name *)
msgs[e_] := Block[{$MessageList = {}}, e; Union[ToString /@ $MessageList]];

isingEnergies[n_, seed_] := BlockRandom[
    With[{spins = 1 - 2 Tuples[{0, 1}, n]}, Total[(spins . UpperTriangularize[RandomReal[{-1, 1}, {n, n}], 1]) spins, {2}]],
    RandomSeeding -> seed
];
diagonalOperator[d_, n_] := QuantumOperator["Diagonal"[d], Range[n], "Label" -> None];
pauliOn[n_, k_, p_] := QuantumOperator[StringReplacePart[StringRepeat["I", n], p, {k, k}]];

run["basis"] := (
    Scan[
        Function[d, With[{qb = QuantumBasis["Fourier"[d]], v = BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, d]], RandomSeeding -> 1]},
            With[{m = qb["Output"]["ReducedMatrix"]},
                With[{qf = timed[QuantumState[QuantumState[v], qb]]},
                    report["Fourier d=" <> ToString[d] <> ": QuantumState[qs, qb] seconds", First[qf]];
                    report["Fourier d=" <> ToString[d] <> ": Inverse[m] seconds", First[timed[Inverse[m]]]];
                    report["Fourier d=" <> ToString[d] <> ": ConjugateTranspose[m] seconds", First[timed[ConjugateTranspose[m]]]];
                    report["Fourier d=" <> ToString[d] <> ": UnitaryMatrixQ[m] exact {seconds, value}", timed[UnitaryMatrixQ[m], 60]];
                    report["Fourier d=" <> ToString[d] <> ": UnitaryMatrixQ[N[m]] {seconds, value}", timed[UnitaryMatrixQ[N[m]]]];
                    If[NumberQ[First[qf]],
                        report["Fourier d=" <> ToString[d] <> ": |QF result - ConjugateTranspose[m].v|", Norm[N[Normal[Last[qf]["StateVector"]] - ConjugateTranspose[m] . v]]]
                    ]
                ]
            ]
        ]],
        {4, 8, 16}
    ];
    report["Fourier d=16: QuditBasis Inverse property {seconds}", First[timed[QuditBasis["Fourier"[16]]["Inverse"]]]]
    (* which named bases are unitary: supplement.wls *)
);

run["spectrum"] := (
    Scan[
        Function[n, With[{qo = diagonalOperator[isingEnergies[n, 7], n]},
            report["diagonal n=" <> ToString[n] <> ": Eigenvalues {seconds}, messages", {First[timed[qo["Eigenvalues"], 300]], msgs[diagonalOperator[isingEnergies[n, 7], n]["Eigenvalues"]]}];
            report["diagonal n=" <> ToString[n] <> ": Eigensystem seconds", First[timed[diagonalOperator[isingEnergies[n, 7], n]["Eigensystem"], 300]]];
            With[{m = qo["MatrixRepresentation"]},
                report["diagonal n=" <> ToString[n] <> ": bare Eigensystem[sparse] seconds", First[timed[Eigensystem[m], 300]]];
                report["diagonal n=" <> ToString[n] <> ": Simplify of that Eigensystem seconds", With[{es = Eigensystem[m]}, First[timed[Simplify[es], 300]]]]
            ]
        ]],
        {8, 10, 12}
    ];
    report["exact diag(2, 1, 2): Eigenvalues, Eigenvectors", With[{qo = QuantumOperator[DiagonalMatrix[{2, 1, 2}]]}, {qo["Eigenvalues"], Normal /@ qo["Eigenvectors"]}]];
    report["machine diag(2., 1., 2.): Eigenvalues", QuantumOperator[DiagonalMatrix[{2., 1., 2.}]]["Eigenvalues"]]
);

run["power"] := (
    With[{qo = diagonalOperator[isingEnergies[12, 7], 12]},
        report["machine diagonal n=12: qo^0.5 seconds", First[timed[qo ^ 0.5]]];
        report["machine diagonal n=12: qo^(1/2) seconds", First[timed[qo ^ (1/2)]]];
        report["machine diagonal n=12: Sqrt[qo] seconds", First[timed[Sqrt[qo]]]];
        report["machine diagonal n=12: |diag(qo^0.5) - entries^0.5|", Max[Abs[Normal[Diagonal[(qo ^ 0.5)["Matrix"]]] - Normal[Diagonal[qo["Matrix"]]] ^ 0.5]]]
    ];
    With[{qo = QuantumOperator[DiagonalMatrix[Range[256]], Range[8]]},
        report["exact diagonal 256: qo^(1/2) {seconds}", First[timed[qo ^ (1/2)]]];
        report["exact diagonal 256: Sqrt[qo] {seconds}", First[timed[Sqrt[qo]]]]
    ]
);

run["liouvillian"] := Scan[
    Function[n, With[{h = diagonalOperator[isingEnergies[n, 7], n], zs = pauliOn[n, #, "Z"] & /@ Range[n]},
        report["n=" <> ToString[n] <> ": full dephasing Liouvillian seconds", First[timed[QuantumOperator["Liouvillian"[h, zs, ConstantArray[0.2, n]]], 600]]];
        report["n=" <> ToString[n] <> ": Hamiltonian part only seconds", First[timed[QuantumOperator["Liouvillian"[h]], 600]]];
        report["n=" <> ToString[n] <> ": one jump operator only seconds", First[timed[QuantumOperator["Liouvillian"[None, {First[zs]}, {0.2}]], 600]]];
        report["n=" <> ToString[n] <> ": the n Pauli-Z operators seconds", First[timed[pauliOn[n, #, "Z"] & /@ Range[n], 600]]]
    ]],
    {5, 6, 7}
];

run["permutation"] := (
    Scan[
        Function[d, With[{qo = QuantumOperator["X"[d]]},
            report["X[" <> ToString[d] <> "]: Eigenvalues seconds", First[timed[qo["Eigenvalues"], 300]]];
            report["X[" <> ToString[d] <> "]: Eigenvalues[PermutationMatrix structured] seconds", First[timed[Eigenvalues[PermutationMatrix[RotateRight[Range[d]], TargetStructure -> "Structured"]]]]];
            report["X[" <> ToString[d] <> "]: matrix equals PermutationMatrix[RotateRight[Range[d]]]", Normal[qo["Matrix"]] === Normal[PermutationMatrix[RotateRight[Range[d]]]]]
        ]],
        {16, 64, 128}
    ];
    Scan[
        Function[n, With[{qo = QuantumCircuitOperator["Fourier"[n]]["QuantumOperator"]},
            report["QFT n=" <> ToString[n] <> ": Eigenvalues seconds", First[timed[qo["Eigenvalues"], 300]]];
            report["QFT n=" <> ToString[n] <> ": Eigenvalues[FourierMatrix structured] seconds", First[timed[Eigenvalues[FourierMatrix[2 ^ n, TargetStructure -> "Structured"]]]]]
        ]],
        {3, 4}
    ]
);

run["qft"] := Scan[
    Function[n, With[{qc = QuantumCircuitOperator["Fourier"[n]], v = BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, 2 ^ n]], RandomSeeding -> 1]},
        With[{qs = QuantumState[v], f = FourierMatrix[2 ^ n, TargetStructure -> "Structured", WorkingPrecision -> MachinePrecision]},
            With[{qf = timed[qc[qs], 300]},
                report["QFT n=" <> ToString[n] <> ": qc[qs] default seconds", First[qf]];
                report["QFT n=" <> ToString[n] <> ": qc[qs] greedy path seconds", First[timed[qc[qs, Method -> {"TensorNetwork", "Path" -> "Greedy"}], 300]]];
                report["QFT n=" <> ToString[n] <> ": FourierMatrix . v seconds", First[timed[f . v]]];
                report["QFT n=" <> ToString[n] <> ": Fourier[v] seconds", First[timed[Fourier[v]]]];
                If[NumberQ[First[qf]], report["QFT n=" <> ToString[n] <> ": |QF - FourierMatrix . v|", Norm[Normal[Last[qf]["StateVector"]] - f . v]]]
            ]
        ]
    ]],
    {14, 16}
];

run[group]
