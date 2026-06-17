#!/usr/bin/env wolframscript
(* Probe 6: verification of the coherent-materialization fix (Section 4b/4c).
   A gate-built StabilizerFrame now carries a relating Pauli per component, so
   fr["StateVector"] materializes the reference component once and obtains the
   others by applying the actual Pauli operators that relate them. This recovers
   the relative phases that the old per-component +1-eigenvector readout dropped.

   Checks, all against the dense Method->"Schrodinger" ground truth:
     - T|+> = {1, e^{iPi/4}}/Sqrt2 exactly (the case the old code swapped).
     - 216 random Clifford+T circuits (n=1..3): frame StateVector == dense, up
       to one global phase (max amplitude error over all of them).
     - <Y> on T|+> = +0.7071 (was sign-flipped to -0.7071).
     - generic "P"[theta] still realizes the half-angle gate diag(1,e^{i th/2})
       (a separate, documented convention), now materialized coherently: frame
       "P"[theta] matches dense "P"[theta/2].
     - manual frames (no relating Paulis, possibly distinct tableaux) are
       unchanged: independent per-component materialization.
   No Quiet, no Print; final expression is the result Association.
   View with:  wolframscript -code 'Print[Get["sf_probe6_fix_verification.wls"]]' *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

(* amplitude error after aligning the single global phase that rotates b onto a *)
phaseAlignedErr[a0_, b0_] := Module[{a = N @ Normalize[Normal[a0]], b = N @ Normalize[Normal[b0]], ph},
    If[Length[a] != Length[b], Return["dim mismatch"]];
    ph = Conjugate[b] . a; ph = If[Abs[ph] < 10.^-12, 1, ph / Abs[ph]];
    Chop[Max[Abs[a - ph b]], 10.^-9]];

backbone[n_] := Join[Table["H" -> i, {i, n}], Table["CNOT" -> {i, i + 1}, {i, n - 1}]];
frameSV[specs_] := Module[{fr = QuantumCircuitOperator[specs][Method -> "Stabilizer"]},
    If[Head[fr] === StabilizerFrame, fr["StateVector"], fr["State"]["StateVector"]]];
denseSV[specs_] := QuantumCircuitOperator[specs][Method -> "Schrodinger"]["StateVector"];

randGate[n_] := With[{g = RandomChoice[If[n >= 2,
        {"H", "S", "X", "Y", "Z", SuperDagger["S"], "CNOT", "CZ", "V", SuperDagger["V"]},
        {"H", "S", "X", "Y", "Z", SuperDagger["S"], "V", SuperDagger["V"]}]]},
    If[MemberQ[{"CNOT", "CZ"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]];
randSpec[n_, gates_, tcount_, seed_] := (SeedRandom[seed];
    Join[backbone[n], Table[randGate[n], {gates}],
        Table[RandomChoice[{"T", SuperDagger["T"]}] -> RandomInteger[{1, n}], {tcount}]]);
randErr[n_, gates_, tcount_, seed_] := Module[{specs = randSpec[n, gates, tcount, seed], fr},
    fr = QuantumCircuitOperator[specs][Method -> "Stabilizer"];
    If[Head[fr] =!= StabilizerFrame, Nothing, phaseAlignedErr[fr["StateVector"], denseSV[specs]]]];
randErrs = DeleteCases[Flatten @ Table[randErr[n, gates, tcount, seed],
    {n, {1, 2, 3}}, {gates, {2, 4, 6}}, {tcount, {1, 2, 3}}, {seed, 1, 8}], Nothing];

tPlus = N @ Normal @ PauliStabilizer[1]["H", 1]["T", 1]["StateVector"];
manual = StabilizerFrame[{{1, PauliStabilizer[1]}, {I, PauliStabilizer[1]["H", 1]}}];

<|
    "T|+> StateVector (== {1, e^{iPi/4}}/Sqrt2)" -> Chop[tPlus],
    "T|+> exact amplitude error" -> Chop[Max[Abs[tPlus - {1, Exp[I Pi/4]} / Sqrt[2]]], 10.^-9],
    "random circuits tested" -> Length[randErrs],
    "random worst amplitude error" -> Max[randErrs],
    "<Y> on T|+> (was -0.7071, now +0.7071)" -> Re[Conjugate[tPlus] . PauliMatrix[2] . tPlus],
    "frame P[0.7] vs dense P[0.35] (half-angle, coherent)" -> phaseAlignedErr[frameSV[{"H" -> 1, "P"[0.7] -> 1}], denseSV[{"H" -> 1, "P"[0.35] -> 1}]],
    "manual frame keeps independent readout" -> Chop[Max[Abs[N @ Normal[manual["StateVector"]] - ({1, 0} + I {1, 1} / Sqrt[2])]], 10.^-9]
|>
