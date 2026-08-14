Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumDistance"]
PackageExport["QuantumSimilarity"]

PackageScope["$QuantumDistances"]



$QuantumDistances = {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}


QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ] := QuantumDistance[qs1, qs2, "Fidelity"]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "Fidelity"] /; qs1["Dimension"] == qs2["Dimension"] :=
    1 - Re[Tr[MatrixPower[qs1["Computational"]["DensityMatrix"] . qs2["Computational"]["DensityMatrix"], 1 / 2]]]

(* Umegaki relative entropy, S(s||t) = Tr[s log s] - Tr[s log t], where the
   cross term is taken as Sum_j <t_j|s|t_j> log q_j over the eigenbasis of t
   alone - a density matrix is routinely singular and 0 log 0 = 0 is what
   decides the answer there. A zero eigenvalue of s drops out of the entropy,
   which keeps a pure s against a full-rank t finite; a zero q_j carrying weight
   sends its term to -Infinity, so a violated support condition
   supp(s) subset of supp(t) surfaces as Quantity[Infinity, "Bits"], which
   QuantumSimilarity carries through as 2^-Infinity == 0. Chop supplies the
   tolerance an inexact spectrum needs and leaves exact input untouched. *)

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "RelativeEntropy"] /; qs1["Dimension"] == qs2["Dimension"] :=
    Block[{
        s = qs1["Computational"]["DensityMatrix"],
        positive = ! TrueQ[Re[#] <= 0] &,
        vals, vecs
    },
        {vals, vecs} = Eigensystem[Normal@qs2["Computational"]["DensityMatrix"]];
        Quantity[
            Chop[(
                Total[# Log[#] & /@ Select[Chop @ Eigenvalues[Normal@s], positive]] -
                Total @ MapThread[If[positive[#2], #2 Log[#1], 0] &, Chop @ {
                    If[positive[#], #, 0] & /@ vals,
                    Conjugate[#] . s . # & /@ Normalize /@ vecs
                }]
            ) / Log[2]],
            "Bits"
        ]
    ]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "RelativePurity"] /; qs1["Dimension"] == qs2["Dimension"] := With[{
    s = qs1["Computational"]["DensityMatrix"],
    t = qs2["Computational"]["DensityMatrix"]
},
    1 - Chop[Tr[s . t]]
]


QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "Trace"] /; qs1["Dimension"] == qs2["Dimension"] := With[{
    m = qs1["Computational"]["DensityMatrix"] - qs2["Computational"]["DensityMatrix"]
},
    Re @ Tr[MatrixPower[ConjugateTranspose[m] . m, 1 / 2]] / 2
]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "Bures"] /; qs1["Dimension"] == qs2["Dimension"]  :=
    Sqrt[2 QuantumDistance[qs1, qs2, "Fidelity"]]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "BuresAngle"] /; qs1["Dimension"] == qs2["Dimension"]  :=
    Re @ ArcCos[1 - QuantumDistance[qs1, qs2, "Fidelity"]]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "HilbertSchmidt"] /; qs1["Dimension"] == qs2["Dimension"] :=
    Norm[qs1["Computational"]["DensityMatrix"] - qs2["Computational"]["DensityMatrix"], "Frobenius"]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "Bloch"] /; qs1["Dimension"] == qs2["Dimension"] :=
    Re @ EuclideanDistance[qs1["BlochCartesianCoordinates"], qs2["BlochCartesianCoordinates"]] / 2


QuantumSimilarity[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, distance_String : "Fidelity"] :=
    Block[{d = QuantumDistance[qs1, qs2, distance]},
        Switch[distance,
            "Fidelity" | "Trace" | "Bloch" | "RelativePurity",
                1 - d,
            "Bures",
                1 - d / Sqrt[2],
            "BuresAngle",
                1 - d / (Pi / 2),
            "HilbertSchmidt",
                1 - d / Sqrt[2],
            "RelativeEntropy",
                2 ^ (-Replace[d, q_Quantity :> QuantityMagnitude[q]]), (* exponential decay for unbounded metric *)
            _,
                1 - d
        ]
    ]

