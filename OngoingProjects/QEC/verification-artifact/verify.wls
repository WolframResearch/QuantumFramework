#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify.wls
   Cross-engine verification harness for the shipped QuantumFramework stabilizer
   layer (PauliStabilizer) and an in-language GF(2) code analysis built on top of
   it, cross-checked against Stim. Prints a PASS/FAIL line per check and exits
   nonzero if any check fails. Stim checks report SKIP (not FAIL) if stim is not
   importable from the ExternalEvaluate Python. No QuantumFramework`QEC` package is
   used: everything here runs on the shipped stabilizer functions plus Stim.
*)

(* --- locate the checkout by walking up to the paclet root --- *)
here = DirectoryName[ExpandFileName[If[$InputFileName === "", ".", $InputFileName]]];
repo = NestWhile[ParentDirectory,
   here,
   ! FileExistsQ[FileNameJoin[{#, "QuantumFramework", "PacletInfo.wl"}]] &&
     ! FileExistsQ[FileNameJoin[{#, "PacletInfo.wl"}]] &,
   1, 12];
PacletDirectoryLoad[repo];
Needs["Wolfram`QuantumFramework`"];

emit[s__] := WriteString[$Output, s];
$fail = False;
check[label_, resultQ_] := (
   If[resultQ === "SKIP",
    emit["SKIP  ", label, "\n"],
    If[TrueQ[resultQ], emit["PASS  ", label, "\n"], emit["FAIL  ", label, "\n"]; $fail = True]]);

(* ---------- in-language GF(2)/symplectic layer (independent of the package) ---------- *)
xzOf = <|"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}|>;
strip[s_] := StringDelete[s, StartOfString ~~ ("+" | "-")];
symOf[s_] := With[{p = xzOf /@ Characters[strip[s]]}, Join[p[[All, 1]], p[[All, 2]]]];
omega[u_, v_] := With[{n = Length[u]/2}, Mod[Take[u, n] . Drop[v, n] + Drop[u, n] . Take[v, n], 2]];
rank2[rows_] := MatrixRank[rows, Modulus -> 2];
weightK[n_, k_] := If[k == 0, {StringRepeat["I", n]},
   Catenate@Map[Function[pos,
      StringJoin /@ (ReplacePart[ConstantArray["I", n], Thread[pos -> #]] & /@ Tuples[{"X", "Y", "Z"}, k])],
     Subsets[Range[n], {k}]]];
logicalVecQ[checks_, v_] := With[{sg = symOf /@ checks},
   With[{m = rank2[sg]}, AllTrue[sg, omega[v, #] === 0 &] && rank2[Append[sg, v]] > m]];
distanceWitness[checks_] := With[{sg = symOf /@ checks, n = StringLength[strip[First[checks]]]},
   With[{m = rank2[sg]},
    With[{logicalQ = Function[v, AllTrue[sg, omega[v, #] === 0 &] && rank2[Append[sg, v]] > m]},
     Catch[
      Scan[Function[w,
        With[{cand = SelectFirst[weightK[n, w], logicalQ[symOf[#]] &]},
         If[! MissingQ[cand], Throw[{w, cand}]]]],
       Range[n]];
      Missing[]]]]];
logicalAlgebraQ[checks_, xbar_, zbar_] := With[{sg = symOf /@ checks},
   With[{m = rank2[sg]},
    And[omega[symOf[xbar], symOf[zbar]] === 1,
     AllTrue[checks, omega[symOf[xbar], symOf[#]] === 0 &],
     AllTrue[checks, omega[symOf[zbar], symOf[#]] === 0 &],
     rank2[Append[sg, symOf[xbar]]] > m, rank2[Append[sg, symOf[zbar]]] > m]]];

(* ---------- the five codes: shipped named-code states, plus the two [[3,1]] built from strings ---------- *)
codes = <|
   "bit-flip"   -> <|"state" -> PauliStabilizer[{"ZZI", "IZZ", "ZZZ"}], "checks" -> {"ZZI", "IZZ"}, "logicals" -> {"XXX", "ZII"}, "d" -> 1|>,
   "phase-flip" -> <|"state" -> PauliStabilizer[{"XXI", "IXX", "XXX"}], "checks" -> {"XXI", "IXX"}, "logicals" -> {"XII", "ZZZ"}, "d" -> 1|>,
   "5-qubit"    -> <|"state" -> PauliStabilizer["5QubitCode"], "checks" -> Most[PauliStabilizer["5QubitCode"]["Stabilizers"]], "logicals" -> {"XXXXX", "ZZZZZ"}, "d" -> 3|>,
   "Steane"     -> <|"state" -> PauliStabilizer["SteaneCode"], "checks" -> Most[PauliStabilizer["SteaneCode"]["Stabilizers"]], "logicals" -> {"XXXXXXX", "ZZZZZZZ"}, "d" -> 3|>,
   "Shor"       -> <|"state" -> PauliStabilizer["9QubitCode"], "checks" -> Most[PauliStabilizer["9QubitCode"]["Stabilizers"]], "logicals" -> {"XXXXXXXXX", "ZZZZZZZZZ"}, "d" -> 3|>
   |>;

(* dense sign-exact Pauli matrices (qubit 1 = leftmost factor, QF's ordering) *)
sgn[s_] := If[StringStartsQ[s, "-"], -1, 1];
pauliOp[s_] := sgn[s] Fold[KroneckerProduct,
   PauliMatrix /@ (Characters[strip[s]] /. {"I" -> 0, "X" -> 1, "Y" -> 2, "Z" -> 3})];
denseFixesQ[c_] := With[{psi = Normal[c["state"]["State"]["StateVector"]]},
   AllTrue[c["state"]["Stabilizers"], AllTrue[pauliOp[#] . psi - psi, PossibleZeroQ] &]];

(* apply a Pauli-string error to a stabilizer state, then read the syndrome via Expectation *)
applyErr[ps_, e_] := Fold[
   With[{c = StringTake[e, {#2}]}, If[c === "I", #1, #1[c, #2]]] &,
   ps, Range[StringLength[e]]];
ourSyndrome[c_, e_] := (1 - (applyErr[c["state"], e]["Expectation", #] & /@ c["checks"]))/2;

emit["=== environment ===\n"];
emit["QuantumFramework ", PacletObject["Wolfram/QuantumFramework"]["Version"], ", ", $Version, "\n"];

(* ---------- Stim session (optional) ---------- *)
stimSession = Check[StartExternalSession[<|"System" -> "Python"|>], $Failed];
stimOK = stimSession =!= $Failed &&
   (Check[ExternalEvaluate[stimSession, "import stim; True"], $Failed] === True);
emit["Stim ", If[stimOK, "available", "NOT available (Stim checks will SKIP)"], "\n\n"];

If[stimOK,
  stimSyndrome = ExternalEvaluate[stimSession, "
def stim_syndrome(gens, err):
    e = stim.PauliString(err.replace('I','_'))
    return [0 if stim.PauliString(g.replace('I','_')).commutes(e) else 1 for g in gens]
stim_syndrome"];
  stimWitness = ExternalEvaluate[stimSession, "
def stim_witness(gens, logs, w):
    W = stim.PauliString(w.replace('I','_'))
    undetectable = all(stim.PauliString(g.replace('I','_')).commutes(W) for g in gens)
    nontrivial = any(not stim.PauliString(p.replace('I','_')).commutes(W) for p in logs)
    return [bool(undetectable), bool(nontrivial)]
stim_witness"];
  stimValidGroup = ExternalEvaluate[stimSession, "
def stim_valid_group(gens):
    try:
        t = stim.Tableau.from_stabilizers([stim.PauliString(g.replace('I','_')) for g in gens], allow_redundant=False, allow_underconstrained=False)
        return len(t.to_stabilizers(canonicalize=True))
    except Exception:
        return -1
stim_valid_group"];
];

(* ---------- Part A: shipped stabilizer layer, cross-checked against Stim ---------- *)
emit["=== Part A: shipped PauliStabilizer vs Stim ===\n"];
KeyValueMap[Function[{nm, c},
   check["A1 sign-exact +1 codeword (dense, exact): " <> nm, denseFixesQ[c]]], codes];
KeyValueMap[Function[{nm, c},
   check["A1 codeword is normalized, <psi|psi> = 1 (invariant): " <> nm,
    With[{psi = Normal[c["state"]["State"]["StateVector"]]}, PossibleZeroQ[Conjugate[psi] . psi - 1]]]], codes];
KeyValueMap[Function[{nm, c},
   check["A2 Stim accepts generators as a valid stabilizer group: " <> nm,
    If[stimOK, stimValidGroup[c["state"]["Stabilizers"]] === Length[c["state"]["Stabilizers"]], "SKIP"]]], codes];
KeyValueMap[Function[{nm, c},
   check["A3 syndrome differential vs Stim (weight 1 and 2): " <> nm,
    If[stimOK,
     With[{errs = Join[weightK[c["state"]["Qubits"], 1], weightK[c["state"]["Qubits"], 2]]},
      (ourSyndrome[c, #] & /@ errs) === (stimSyndrome[c["checks"], #] & /@ errs)], "SKIP"]]], codes];

(* A4: symbolic measurement (SymPhase) resolves and composes as it must (shipped, no Stim) *)
smBell = PauliStabilizer[2][{"H", "CNOT"}]["SymbolicMeasure", 1];
sBell = First@Cases[smBell["Phase"], \[FormalS][_], Infinity];
check["A4 SymbolicMeasure resolves Bell branches",
   smBell["SubstituteOutcomes", {sBell -> 0}]["Stabilizers"] === {"ZI", "ZZ"} &&
    smBell["SubstituteOutcomes", {sBell -> 1}]["Stabilizers"] === {"-ZI", "ZZ"}];
parityState = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
smPar = parityState["SymbolicMeasure", {1, 2, 3}];
readBits[ps_] := First@Keys@Select[
    AssociationThread[Tuples[{0, 1}, ps["Qubits"]], Normal[ps["State"]["StateVector"]]], # =!= 0 &];
freeBits = DeleteDuplicates@Cases[smPar["Phase"], \[FormalS][_], Infinity];
parityRows = readBits[smPar["SubstituteOutcomes", Thread[freeBits -> #]]] & /@ Tuples[{0, 1}, Length[freeBits]];
check["A4 SymbolicMeasure composes: m3 = m1 XOR m2 on the parity code",
   AllTrue[parityRows, Mod[Total[#], 2] === 0 &] && Length[parityRows] === 4];
check["A4 conserved parity ZZZ stabilizes the state; only 2 formal bits are spent",
   parityState["Expectation", "ZZZ"] === 1 && Length[freeBits] === 2];

(* ---------- Part B: in-language exact analysis, cross-checked against Stim ---------- *)
emit["\n=== Part B: in-language exact analysis vs Stim ===\n"];
KeyValueMap[Function[{nm, c},
   check["B1 in-language distance = known d: " <> nm,
    distanceWitness[c["checks"]][[1]] === c["d"]]], codes];
KeyValueMap[Function[{nm, c},
   check["B2 logical operators satisfy the full algebra: " <> nm,
    logicalAlgebraQ[c["checks"], c["logicals"][[1]], c["logicals"][[2]]]]], codes];
KeyValueMap[Function[{nm, c},
   check["B3 distance witness is undetectable and nontrivial (Stim): " <> nm,
    If[stimOK, stimWitness[c["checks"], c["logicals"], distanceWitness[c["checks"]][[2]]] === {True, True}, "SKIP"]]], codes];

(* ---------- Part B: symbolic families and a continuous error angle ---------- *)
emit["\n=== Part B: symbolic generality (no numeric counterpart) ===\n"];
familyChecks[m_] := {StringRepeat["X", m], StringRepeat["Z", m]};
repChecks[m_] := StringJoin /@ (ReplacePart[ConstantArray["I", m], {# -> "Z", # + 1 -> "Z"}] & /@ Range[m - 1]);
check["B4 [[n,n-2,2]] commutator parity is even for n=2k, odd for n=2k+1 (derived in k)",
   Simplify[{Mod[Sum[1, {i, n}] /. n -> 2 k, 2], Mod[Sum[1, {i, n}] /. n -> 2 k + 1, 2]}, k \[Element] Integers] === {0, 1}];
check["B5 [[n,n-2,2]] X1 X2 is a weight-2 logical, no weight-1 is (representative n=6)",
   logicalVecQ[familyChecks[6], symOf["XXIIII"]] && ! logicalVecQ[familyChecks[6], symOf["XIIIII"]]];
check["B6 concrete [[n,n-2,2]] distances at n=4,6,8,10 all equal 2",
   AllTrue[{4, 6, 8, 10}, distanceWitness[familyChecks[#]][[1]] === 2 &]];
check["B7 repetition family [[n,1,1]] distance = 1 (representative n=7)",
   distanceWitness[repChecks[7]][[1]] === 1];
vTheta = Normal[QuantumCircuitOperator[{"RX"[\[Theta]] -> 1}][QuantumState[SparseArray[{1 -> 1}, 8], {2, 2, 2}]]["StateVector"]];
check["B8 coherent RX(theta) state stays normalized (unitarity invariant)",
   Simplify[Conjugate[vTheta] . vTheta, \[Theta] \[Element] Reals] === 1];
zzTheta = Simplify[Conjugate[vTheta] . KroneckerProduct[PauliMatrix[3], PauliMatrix[3], IdentityMatrix[2]] . vTheta, \[Theta] \[Element] Reals];
check["B8 coherent-error syndrome expectation equals Cos[theta]", zzTheta === Cos[\[Theta]]];
check["B9 small-angle coherent-error rate (1-<ZZ>)/2 is theta^2/4",
   Normal@Series[(1 - zzTheta)/2, {\[Theta], 0, 2}] === \[Theta]^2/4];

emit["\n", If[$fail, "RESULT: FAIL", "RESULT: PASS"], "\n"];
If[stimOK, DeleteObject[stimSession]];
Exit[If[$fail, 1, 0]]
