PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
report[label_, value_] := WriteString[$Output, "supplement2\t", label, "\t", ToString[value, InputForm, PageWidth -> Infinity], "\n"];
SetAttributes[timed, HoldFirst];
timed[e_, limit_ : 120] := (ClearSystemCache[]; Replace[AbsoluteTiming[TimeConstrained[e, limit, $timedOut]],
    {{_, $timedOut} :> {"TIMEOUT>" <> ToString[limit], Missing["TimedOut"]}, {t_, v_} :> {Round[t, 0.0001], v}}]);
v = BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, 2^10]], RandomSeeding -> 1];
report["|Fourier[v] - FourierMatrix[2^10] . v|", Norm[Fourier[v] - FourierMatrix[2^10, WorkingPrecision -> MachinePrecision] . v]];
With[{m = QuditBasis["HoggarSIC"]["ReducedMatrix"]}, report["HoggarSIC: dimensions, exact UnitaryMatrixQ {seconds, value}", {Dimensions[m], timed[UnitaryMatrixQ[m], 120]}]];
With[{qo = QuantumOperator[DiagonalMatrix[{2., 1., 2.}]]}, report["machine diag(2., 1., 2.): Eigenvectors", Normal /@ qo["Eigenvectors"]]];
With[{qo = QuantumOperator[DiagonalMatrix[{1, -1, I, -I, 1}]]}, report["exact diag(1, -1, i, -i, 1): Eigenvalues, Eigenvectors", {qo["Eigenvalues"], Normal /@ qo["Eigenvectors"]}]];
With[{qo = QuantumOperator["Z"[3]]}, report["Z[3] clock: Eigenvalues, diagonal", {qo["Eigenvalues"], Normal[Diagonal[qo["MatrixRepresentation"]]]}]];
report["Eigenvalues of Z measurement operator outcomes order", QuantumMeasurementOperator[QuantumOperator["Z"]]["Eigenvalues"]];
