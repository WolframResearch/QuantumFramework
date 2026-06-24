#!/usr/bin/env wolframscript
(* Probe 2: headline scaling. Fix an n-qubit Clifford backbone, insert k T gates,
   apply via Method -> "Stabilizer", measure StabilizerFrame component count, wall time,
   and memory vs k. ClearSystemCache[] before each timing. No Quiet, no Print; final
   expression is the result (list of Associations) -> a Dataset.
   View with:  wolframscript -code 'Print[Get["sf_probe2_scaling.wls"]]'
   Default ks = 0..16 (149 MB / ~2 s peak). The report also ran k=17,18 (298/596 MB). *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

n = 4;
backbone = Join[Table["H" -> i, {i, n}], Table["CNOT" -> {i, i + 1}, {i, n - 1}]];
tGates[k_] := Table["T" -> (Mod[j - 1, n] + 1), {j, k}];

run[k_] := Module[{qc, res, t},
    qc = QuantumCircuitOperator[Join[backbone, tGates[k]]];
    ClearSystemCache[];
    {t, res} = AbsoluteTiming[qc[Method -> "Stabilizer"]];
    <|
        "k" -> k,
        "Head" -> Head[res],
        "components" -> Which[Head[res] === StabilizerFrame, res["Length"], Head[res] === PauliStabilizer, 1, True, Missing[]],
        "2^k" -> 2^k,
        "==2^k" -> (Which[Head[res] === StabilizerFrame, res["Length"], True, 1] === 2^k),
        "t(ms)" -> Round[1000. t, 0.1],
        "MB" -> Round[ByteCount[res]/1024.^2, 0.01]
    |>
];

Dataset[run /@ Range[0, 16]]
