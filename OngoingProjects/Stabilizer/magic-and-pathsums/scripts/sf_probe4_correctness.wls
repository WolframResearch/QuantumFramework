#!/usr/bin/env wolframscript
(* Probe 4: correctness anchor.
   (4a) the per-gate decomposition rule (1+c)/2 I + (1-c)/2 Z is an EXACT matrix identity
        (= T for c=e^{iPi/4}; = diag(1,e^{i th/2}) = P[th/2] for "P"[th] -> the half-angle
        convention mismatch vs the gate-level QuantumOperator["P"[th]] = diag(1,e^{i th})).
   (4b) BUT the dense materialization fr["StateVector"]/["State"] is WRONG for >=2-component
        frames: each component's PauliStabilizer["State"] global phase is an arbitrary
        NullSpace representative (Conversions.m:69-80) with no recovered "GlobalPhase" key
        (set only by PauliStabilizer[qs]/[qo]), so the coherent sum mixes incoherent phases.
        Fidelity vs the dense Schrodinger state is 1 only for single-component (rank-1) frames.
   (4c) the one guarding test (Tests/Stabilizer/Formulas.wlt:811, S11-T-on-Plus-equals-T-State)
        masks this: it checks Total[diff]==0 (sum of a difference vector), which is blind to
        the swapped-amplitude error; Max[Abs[diff]] (or Norm) reveals it.
   No Quiet, no Print; final expression is the result Association.
   View with:  wolframscript -code 'Print[Get["sf_probe4_correctness.wls"]]' *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

id2 = IdentityMatrix[2]; Zmat = {{1, 0}, {0, -1}};
cmpFid[a0_, b0_] := Module[{a = N @ Normal[a0], b = N @ Normal[b0]},
    Round[Abs[Conjugate[a] . b]^2 / (Norm[a]^2 Norm[b]^2), 10.^-12]];

cT = Exp[I Pi/4];                 (* "T" routes to "P"[Pi/2] -> c = e^{iPi/4} *)
theta = 0.7; cP = Exp[I theta/2];

backbone[n_] := Join[Table["H" -> i, {i, n}], Table["CNOT" -> {i, i + 1}, {i, n - 1}]];
fidRow[n_, k_] := Module[{specs, fr, svF, svD},
    specs = Join[backbone[n], Table["T" -> (Mod[j - 1, n] + 1), {j, k}]];
    fr = QuantumCircuitOperator[specs][Method -> "Stabilizer"];
    svF = If[Head[fr] === StabilizerFrame, fr["StateVector"], fr["State"]["StateVector"]];
    svD = QuantumCircuitOperator[specs][Method -> "Schrodinger"]["StateVector"];
    <|"n" -> n, "k" -> k, "components" -> If[Head[fr] === StabilizerFrame, fr["Length"], 1],
      "fidelity(frame,dense)" -> cmpFid[svF, svD]|>];

svTplus = N @ Normal @ PauliStabilizer[1]["H", 1]["T", 1]["StateVector"];
trueTplus = N @ ({1, Exp[I Pi/4]}/Sqrt[2]);

<|
    "(4a) (1+c)/2 I + (1-c)/2 Z == T" ->
        (Chop[((1 + cT)/2 id2 + (1 - cT)/2 Zmat) - DiagonalMatrix[{1, Exp[I Pi/4]}]] == {{0, 0}, {0, 0}}),
    "(4a) frame P[theta] yields [2,2]=" -> Chop[((1 + cP)/2 id2 + (1 - cP)/2 Zmat)[[2, 2]]],
    "(4a) gate-level P[theta] [2,2]="   -> Chop @ QuantumOperator["P"[theta]]["Matrix"][[2, 2]],

    "(4b) fidelity frame-vs-dense table" -> Flatten[Table[fidRow[n, k], {n, {1, 2, 3}}, {k, {0, 1, 2, 3}}], 1],

    "(4c) frame T|+> StateVector"        -> Chop[svTplus],
    "(4c) true  T|+>"                    -> trueTplus,
    "(4c) Total[diff]==0 (test's check)" -> (Chop[Total[svTplus - trueTplus]] == 0),
    "(4c) Max|diff|==0 (correct check)"  -> (Chop[Max[Abs[svTplus - trueTplus]]] == 0),
    "(4c) fidelity frame T|+> vs true"   -> cmpFid[svTplus, trueTplus],

    "rank-1 T|0> fidelity"               -> cmpFid[PauliStabilizer[1]["T", 1]["StateVector"], {1, 0}]
|>
