PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 60], " ", ToString[x, InputForm]]);
(* seconds of one call after ClearSystemCache[], or "timeout" / "memory" *)
SetAttributes[sec, HoldFirst];
sec[e_, limit_ : 60] := (ClearSystemCache[]; Replace[MemoryConstrained[TimeConstrained[First[AbsoluteTiming[e;]], limit, "timeout"], 3 10^9, "memory"], t_Real :> N[Round[t, 1 / 10000]]]);
SeedRandom[7];

(* The structure of every named operator, basis and circuit at its default arguments, read off its matrix in the computational basis. *)
class[m0_] := With[{m = SparseArray[m0]}, With[{n = Length[m], pos = m["NonzeroPositions"], vals = m["NonzeroValues"]}, Which[
    ! SquareMatrixQ[m], "non-square",
    AllTrue[pos, #[[1]] == #[[2]] &], If[Length[pos] == n && AllTrue[vals, # === 1 &], "identity", "diagonal"],
    Length[pos] == n && DuplicateFreeQ[pos[[All, 1]]] && DuplicateFreeQ[pos[[All, 2]]], If[AllTrue[vals, # === 1 &], "permutation", "phase permutation"],
    n <= 64 && AllTrue[Flatten[Normal[m] - FourierMatrix[n]], PossibleZeroQ], "Fourier",
    n <= 64 && AllTrue[Flatten[Normal[m] - ConjugateTranspose[FourierMatrix[n]]], PossibleZeroQ], "inverse Fourier",
    True, "other"]]];
classes[names_, matrix_] := Counts @ Map[Replace[TimeConstrained[Check[class[matrix[#]], "error"], 20, "timeout"], Except[_String] -> "error"] &, names];
P["named operators by structure", classes[Wolfram`QuantumFramework`PackageScope`$QuantumOperatorNames, QuantumOperator[#]["MatrixRepresentation"] &]];
P["named circuits by structure", classes[Wolfram`QuantumFramework`PackageScope`$QuantumCircuitOperatorNames, QuantumCircuitOperator[#]["QuantumOperator"]["MatrixRepresentation"] &]];

(* Spectra: QF's general eigensolver against what the structure gives. *)
With[{h = RandomReal[{-1, 1}, 2^10]}, With[{qo = QuantumOperator["Diagonal"[h], Range[10], "Label" -> None]},
    P["Diagonal[h], 10 qubits: Eigenvalues s", sec[qo["Eigenvalues"]]];
    P["  equal to the diagonal", Max[Abs[Sort[Re[qo["Eigenvalues"]]] - Sort[h]]]];
    P["  the diagonal itself s", sec[Diagonal[qo["Matrix"]]]]]];
P["\"Z\"[128]: Eigenvalues s", sec[QuantumOperator["Z"[128]]["Eigenvalues"]]];
P["\"X\"[64]: Eigenvalues s", sec[QuantumOperator["X"[64]]["Eigenvalues"]]];
P["\"X\"[64] @ \"Z\"[64]: Eigenvalues s", sec[(QuantumOperator["X"[64]] @ QuantumOperator["Z"[64]])["Eigenvalues"]]];
P["QFT operator, 4 qubits: Eigenvalues s", sec[QuantumCircuitOperator["Fourier"[4]]["QuantumOperator"]["Eigenvalues"]]];

(* Change of basis into the Fourier basis: exact inverse of the basis matrix, against its conjugate transpose. *)
Do[With[{s = QuantumState[Normalize[RandomComplex[{-1 - I, 1 + I}, d]]], b = QuantumBasis["Fourier"[d]]},
    P["Fourier[" <> ToString[d] <> "] change of basis s", sec[QuantumState[s, b], 30]];
    P["  ConjugateTranspose of the basis matrix s", sec[ConjugateTranspose[b["Output"]["ReducedMatrix"]]]]], {d, {8, 16}}];

(* The QFT applied to a state: the circuit against Fourier on the amplitudes. *)
With[{v = Normalize[RandomComplex[{-1 - I, 1 + I}, 2^14]]}, Module[{r},
    P["Fourier[14] circuit on a state s", sec[r = QuantumCircuitOperator["Fourier"[14]][QuantumState[v]]]];
    P["  Fourier[v] s", sec[Fourier[v]]];
    P["  agreement", Norm[Normal[r["StateVector"]] - Fourier[v]]]]];

(* A phase oracle: the circuit of multi-controlled Z gates against the sign vector of the truth table. *)
With[{n = 10}, Module[{bf = BooleanFunction[RandomInteger[2^(2^n) - 1], n], qc, r, v = ConstantArray[1 / Sqrt[2.^n], 2^n], signs},
    P["PhaseOracle, random function of 10 variables: build s", sec[qc = QuantumCircuitOperator["PhaseOracle"[bf]]]];
    P["  gates", Length[qc["Elements"]]];
    P["  apply to the uniform state s", sec[r = qc[QuantumState[v, ConstantArray[2, n]]]]];
    P["  sign vector of the truth table s", sec[signs = (-1) ^ Boole[bf @@@ Tuples[{False, True}, n]]]];
    P["  agreement", With[{w = Normal[r["StateVector"]]}, Min[Norm[w - signs v], Norm[w + signs v]]]]]];

(* Measuring a diagonal observable: the eigensystem is a small share; applying a measurement with 64 outcomes is the cost. *)
With[{n = 6}, With[{h = RandomReal[{-1, 1}, 2^n], v = Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]]},
    With[{q = QuantumOperator["Diagonal"[h], Range[n], "Label" -> None], s = QuantumState[v, ConstantArray[2, n]]},
        P["Diagonal[h], 6 qubits: Eigensystem s", sec[q["Eigensystem"]]];
        P["  measurement of it applied to a state s", sec[QuantumMeasurementOperator[q][s]]];
        P["  computational measurement, 64 outcomes s", sec[QuantumMeasurementOperator[Range[n]][s]]]]]];
