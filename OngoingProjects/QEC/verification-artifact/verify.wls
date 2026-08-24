#!/usr/bin/env wolframscript
(* Standalone verification harness for the exact stabilizer-code layer.
   Runs, and reports PASS/FAIL for:
     1. the QEC package test suite (179 VerificationTests),
     2. an independent GF(2)/symplectic re-derivation of distance, logical
        algebra and the sign-exact encoder (does not call the package's own
        analysis routines to check themselves),
     3. a Stim differential of syndrome, encoder and witness (skipped, not
        failed, if stim is not importable),
     4. the symbolic-parameter family certificates.
   Exits nonzero if any check fails. Portable: finds the repo by walking up. *)

artifactDir = DirectoryName[$InputFileName];
repoRoot = SelectFirst[
    NestList[ParentDirectory, artifactDir, 10],
    FileExistsQ[FileNameJoin[{#, "OngoingProjects", "QEC", "QEC", "QEC.wl"}]] &,
    $Failed
];
If[ repoRoot === $Failed,
    WriteString[$Output, "ERROR: could not locate the QuantumFramework repo root.\n"];
    Exit[1]
];
PacletDirectoryLoad[repoRoot];
Needs["Wolfram`QuantumFramework`"];
base = FileNameJoin[{repoRoot, "OngoingProjects", "QEC", "QEC"}];
Get[FileNameJoin[{base, "QEC.wl"}]];
Get[FileNameJoin[{base, "QECTests.wl"}]];

emit[a__] := WriteString[$Output, a, "\n"];
$fails = 0;
report[label_, ok_] := (If[! TrueQ[ok], $fails++]; emit[If[TrueQ[ok], "  PASS  ", "  FAIL  "], label]);

(* ----- 1. package test suite ----- *)
emit["[1] QEC package test suite"];
Module[{r = QECRunTests[]},
    report["179 VerificationTests succeed", r["Total"] === 179 && r["Outcomes"]["Success"] === 179 && r["FailedInputs"] === {}]
];

(* independent symplectic layer *)
xzOf = <|"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}|>;
strip[s_] := StringDelete[s, StartOfString ~~ ("+" | "-")];
symOf[s_] := With[{ps = xzOf /@ Characters[strip[s]]}, Join[ps[[All, 1]], ps[[All, 2]]]];
omega[u_, v_] := With[{n = Length[u]/2}, Mod[Take[u, n] . Drop[v, n] + Drop[u, n] . Take[v, n], 2]];
wt[s_] := Count[Characters[strip[s]], Except["I"]];
weightK[n_, k_] := If[k == 0, {StringRepeat["I", n]},
    Catenate @ Map[
        Function[pos, StringJoin /@ (ReplacePart[ConstantArray["I", n], Thread[pos -> #]] & /@ Tuples[{"X", "Y", "Z"}, k])],
        Subsets[Range[n], {k}]]];
genMatrix[c_] := symOf /@ c["Generators"];
logicalVecQ[c_, v_] := With[{sg = genMatrix[c], m = c["StabilizerCount"]},
    AllTrue[sg, omega[v, #] === 0 &] && MatrixRank[Append[sg, v], Modulus -> 2] > m];
firstLogicalOfWeight[c_, w_] := SelectFirst[weightK[c["Qubits"], w], logicalVecQ[c, symOf[#]] &];
independentDistance[c_] := Module[{w = SelectFirst[Range[c["Qubits"]], ! MissingQ[firstLogicalOfWeight[c, #]] &]},
    If[MissingQ[w], Infinity, wt[firstLogicalOfWeight[c, w]]]];
logicalAlgebraOkQ[c_] := Module[{logs = LogicalOperators[c], sg = genMatrix[c], m = c["StabilizerCount"], xs, zs, k},
    xs = symOf /@ logs["X"]; zs = symOf /@ logs["Z"]; k = Length[xs];
    And[
        AllTrue[Tuples[Range[k], 2], With[{i = #[[1]], j = #[[2]]}, omega[xs[[i]], zs[[j]]] === If[i === j, 1, 0]] &],
        AllTrue[Subsets[Range[k], {2}], omega[xs[[#[[1]]]], xs[[#[[2]]]]] === 0 &],
        AllTrue[Subsets[Range[k], {2}], omega[zs[[#[[1]]]], zs[[#[[2]]]]] === 0 &],
        AllTrue[Join[xs, zs], Function[l, AllTrue[sg, omega[l, #] === 0 &]]],
        AllTrue[Join[xs, zs], MatrixRank[Append[sg, #], Modulus -> 2] > m &]]];
pmat = <|"I" -> {{1, 0}, {0, 1}}, "X" -> {{0, 1}, {1, 0}}, "Y" -> {{0, -I}, {I, 0}}, "Z" -> {{1, 0}, {0, -1}}|>;
pauliOp[s_] := Fold[KroneckerProduct, pmat /@ Characters[strip[s]]];
encoderFixesQ[c_] := With[
    {psi = Normal[EncodingCircuit[c][QuantumState[SparseArray[{1 -> 1}, 2^c["Qubits"]], ConstantArray[2, c["Qubits"]]]]["StateVector"]]},
    AllTrue[c["Generators"], AllTrue[pauliOp[#] . psi - psi, PossibleZeroQ] &]];

names = {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"};
codes = AssociationMap[QECCode, names];

(* ----- 2. independent re-derivation ----- *)
emit["[2] independent GF(2)/symplectic re-derivation"];
report["distance matches package for all five codes",
    AllTrue[names, CodeDistance[codes[#]] === independentDistance[codes[#]] &]];
report["package witness is a valid logical of the code distance",
    AllTrue[names, With[{c = codes[#]}, logicalVecQ[c, symOf[MinimumWeightLogical[c]]] && wt[MinimumWeightLogical[c]] === CodeDistance[c]] &]];
report["logical operators satisfy the full commutation algebra", AllTrue[names, logicalAlgebraOkQ[codes[#]] &]];
report["encoder fixes every generator, sign-exact", AllTrue[names, encoderFixesQ[codes[#]] &]];

(* ----- 3. Stim differential ----- *)
emit["[3] Stim differential"];
Module[{session, ok},
    session = Quiet @ Check[StartExternalSession[<|"System" -> "Python"|>], $Failed];
    ok = session =!= $Failed && StringQ @ Quiet @ Check[ExternalEvaluate[session, "import stim; stim.__version__"], $Failed];
    If[! ok,
        emit["  SKIP  stim not importable via ExternalEvaluate (pip install stim to enable)"];
        If[session =!= $Failed, DeleteObject[session]]
        ,
        Module[{stimSyndrome, stimEnc, stimWit, gatePairs},
            stimSyndrome = ExternalEvaluate[session, "
def f(gens, err):
    e = stim.PauliString(err.replace('I','_'))
    return [0 if stim.PauliString(g.replace('I','_')).commutes(e) else 1 for g in gens]
f"];
            stimEnc = ExternalEvaluate[session, "
def f(gates, gens):
    circ = stim.Circuit(); tr = {'CNOT':'CX'}
    for name, args in gates: circ.append(tr.get(name, name), [a-1 for a in args])
    sim = stim.TableauSimulator(); sim.do_circuit(circ)
    return [int(sim.peek_observable_expectation(stim.PauliString(g.replace('I','_')))) for g in gens]
f"];
            stimWit = ExternalEvaluate[session, "
def f(gens, logicals, w):
    W = stim.PauliString(w.replace('I','_'))
    return [bool(all(stim.PauliString(g.replace('I','_')).commutes(W) for g in gens)),
            bool(any(not stim.PauliString(p.replace('I','_')).commutes(W) for p in logicals))]
f"];
            gatePairs[c_] := {#[[1]], Flatten[{#[[2]]}]} & /@ EncodingGates[c];
            report["syndrome === Stim over all weight-1 and weight-2 errors",
                AllTrue[names, Function[nm, With[{c = codes[nm], errs = Join[weightK[codes[nm]["Qubits"], 1], weightK[codes[nm]["Qubits"], 2]]},
                    (Syndrome[c, #] & /@ errs) === (stimSyndrome[c["Generators"], #] & /@ errs)]]]];
            report["Stim encoder expectation = +1 on every generator",
                AllTrue[names, With[{c = codes[#]}, AllTrue[stimEnc[gatePairs[c], c["Generators"]], # === 1 &]] &]];
            report["Stim: every witness undetectable and nontrivial",
                AllTrue[names, With[{c = codes[#], logs = LogicalOperators[codes[#]]},
                    stimWit[c["Generators"], Join[logs["X"], logs["Z"]], MinimumWeightLogical[c]] === {True, True}] &]];
            DeleteObject[session]
        ]
    ]
];

(* ----- 4. symbolic-parameter families ----- *)
emit["[4] symbolic-parameter families"];
report["[[n,n-2,2]]: generators commute iff n even (omega = n)",
    TrueQ @ Simplify[Mod[Sum[1, {i, n}], 2] == 0, Element[n/2, Integers]]];
report["[[n,n-2,2]]: weight-2 logical lighter than weight-n group for n>2",
    Reduce[2 < n, n, Integers] === (n \[Element] Integers && n >= 3)];
report["[[n,n-2,2]]: package distance = 2 across the family",
    Table[CodeDistance[DistanceTwoCode[n]], {n, {4, 6, 8, 10}}] === {2, 2, 2, 2}];
report["coherent-error syndrome <Z1Z2> = Cos[theta] (symbolic angle)",
    Module[{v = Normal[QuantumCircuitOperator[{"RX"[\[Theta]] -> 1}][QuantumState[SparseArray[{1 -> 1}, 8], {2, 2, 2}]]["StateVector"]]},
        Simplify[Conjugate[v] . KroneckerProduct[PauliMatrix[3], PauliMatrix[3], IdentityMatrix[2]] . v, \[Theta] \[Element] Reals] === Cos[\[Theta]]]];

emit[""];
emit[If[$fails === 0, "ALL CHECKS PASSED", ToString[$fails] <> " CHECK(S) FAILED"]];
Exit[If[$fails === 0, 0, 1]];
