PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
report[label_, value_] := WriteString[$Output, "supplement\t", label, "\t", ToString[value, InputForm, PageWidth -> Infinity], "\n"];
SetAttributes[timed, HoldFirst];
timed[e_, limit_ : 120] := (ClearSystemCache[]; Replace[AbsoluteTiming[TimeConstrained[e, limit, $timedOut]],
    {{_, $timedOut} :> {"TIMEOUT>" <> ToString[limit], Missing["TimedOut"]}, {t_, v_} :> {Round[t, 0.0001], v}}]);
basisQ = Wolfram`QuantumFramework`PackageScope`QuditBasisQ;
report["named bases at default arguments: {name, square, unitary numerically}",
    Map[Function[name, With[{qb = Check[QuditBasis[name], $Failed]},
        If[TrueQ[basisQ[qb]],
            With[{m = Normal[qb["ReducedMatrix"]]}, {name, SquareMatrixQ[m], SquareMatrixQ[m] && UnitaryMatrixQ[N[m], Tolerance -> 10^-10]}],
            {name, "no default"}]]],
        Wolfram`QuantumFramework`PackageScope`$QuditBasisNames]];
Scan[Function[d, With[{m = QuantumBasis["Fourier"[d]]["Output"]["ReducedMatrix"]},
    report["Fourier d=" <> ToString[d] <> ": exact UnitaryMatrixQ {seconds, value}", timed[UnitaryMatrixQ[m], 300]];
    report["Fourier d=" <> ToString[d] <> ": Inverse[N[m]] seconds", First[timed[Inverse[N[m]]]]];
    report["Fourier d=" <> ToString[d] <> ": ConjugateTranspose[m] equals FourierMatrix-based conjugate transpose", Simplify[Normal[ConjugateTranspose[m]] - ConjugateTranspose[FourierMatrix[d]]] === ConstantArray[0, {d, d}]]
]], {16, 32, 64}];
With[{m = QuantumBasis["Fourier"[QuditBasis["PauliX"[3]]]]["Output"]["ReducedMatrix"]},
    report["Fourier over PauliX[3]: exact UnitaryMatrixQ", UnitaryMatrixQ[m]]];
With[{m = {{a, b}, {-Conjugate[b], Conjugate[a]}}}, report["symbolic 2x2: UnitaryMatrixQ {seconds, value}", timed[UnitaryMatrixQ[m], 60]]];
report["MatrixInverse definition", Definition[Wolfram`QuantumFramework`PackageScope`MatrixInverse] // InputForm // ToString // StringLength];
