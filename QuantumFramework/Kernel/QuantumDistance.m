Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumDistance"]
PackageExport["QuantumSimilarity"]

PackageScope["$QuantumDistances"]



$QuantumDistances = {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}


QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ] := QuantumDistance[qs1, qs2, "Fidelity"]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "Fidelity"] /; qs1["Dimension"] == qs2["Dimension"] :=
    1 - Re[Tr[MatrixPower[qs1["Computational"]["DensityMatrix"] . qs2["Computational"]["DensityMatrix"], 1 / 2]]]

(* Umegaki relative entropy, S(s||t) = Tr[s log s] - Tr[s log t], taken from the
   two eigendecompositions because a density matrix is routinely singular and
   0 log 0 = 0 is what decides the answer there: a zero eigenvalue of s
   contributes nothing, which keeps a pure s against a full-rank t finite, while
   a zero eigenvalue of t diverges only where s has weight - the support
   condition supp(s) subset of supp(t), reported as Quantity[Infinity, "Bits"],
   which QuantumSimilarity carries through as 2^-Infinity == 0. *)

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "RelativeEntropy"] /; qs1["Dimension"] == qs2["Dimension"] :=
    Block[{
        s = Normal @ qs1["Computational"]["DensityMatrix"],
        t = Normal @ qs2["Computational"]["DensityMatrix"],
        exact, zeroQ, sv, svec, tv, tvec, entropy, crossTerm, supported = True
    },
        (* an inexact eigenvalue carries rounding noise, so it is tested against
           a tolerance and its stray imaginary part dropped; an exact one needs
           neither, and either would spend the exactness the input still has *)
        exact = Precision[{s, t}] === Infinity;
        zeroQ = If[exact, TrueQ[# <= 0] &, TrueQ[Re[#] <= 1*^-10] &];

        {sv, svec} = Eigensystem[s];
        {tv, tvec} = Eigensystem[t];
        If[! exact, sv = Re[sv]; tv = Re[tv]];
        svec = Normalize /@ svec; tvec = Normalize /@ tvec;

        entropy = Total[If[zeroQ[#], 0, # Log[#]] & /@ sv];

        crossTerm = Total @ Flatten @ Table[
            With[{p = sv[[i]], q = tv[[j]], overlap = Abs[Conjugate[svec[[i]]] . tvec[[j]]] ^ 2},
                Which[
                    zeroQ[p] || zeroQ[overlap], 0,
                    zeroQ[q], (supported = False; 0),
                    True, p overlap Log[q]
                ]
            ],
            {i, Length[sv]}, {j, Length[tv]}
        ];

        (* Chop only touches approximate numbers, so it leaves exact input alone *)
        Quantity[If[supported, Chop[(entropy - crossTerm) / Log[2]], Infinity], "Bits"]
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

