Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumEntangledQ"]
PackageExport["QuantumEntanglementMonotone"]
PackageScope["$QuantumEntanglementMonotones"]



$QuantumEntanglementMonotones = {
    "Concurrence", "Negativity", "LogNegativity", "EntanglementEntropy", "RenyiEntropy", "Realignment",
    "MutualInformationI", "MutualInformationJ", "Discord"
}

QuantumEntangledQ[qs_ ? QuantumStateQ, biPartition_ : Automatic, method_String : "Realignment"] /; MemberQ[$QuantumEntanglementMonotones, method] :=
    Enclose[Chop[ConfirmMatch[QuantumEntanglementMonotone[qs, biPartition, method], _ ? NumericQ]] > 0, Indeterminate &]



QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition : Except[_String] : Automatic] :=
    QuantumEntanglementMonotone[qs, biPartition, "Concurrence"]


y[{j_Integer, k_Integer}, n_Integer] /; 1 <= j < k <= n := SparseArray[{{j, k} -> -I, {k, j} -> I}, {n, n}]

Y[n_] := Y[n] = Catenate @ Table[y[{j, k}, n], {k, 2, n}, {j, k - 1}]

(* Rungta-Buzek-Caves-Hillery-Milburn I-concurrence (PRA 64, 042315 (2001)): the concurrence-vector
   components are the per-generator-pair Wootters combinations lambda1 - Sum[rest] of the singular
   values of Sqrt[rho].Sqrt[rho-tilde], with the spin-flipped rho-tilde = o.Conjugate[rho].o built
   from the su(d) generator pair o = yn (x) ym. Those singular values are the square roots of the
   eigenvalues of rho.rho-tilde. Norm of the vector is exact for a pure state of any dimension
   (= Sqrt[2 (1 - Tr[rhoA^2])]) and exact for two qubits (the textbook Wootters concurrence); for a
   d1 d2 > 4 mixed state it is a lower bound on the convex-roof concurrence, which has no closed form,
   and (unlike the two-qubit Wootters concurrence) it is frame-dependent: not invariant under local
   unitaries, though it stays below the LU-invariant convex-roof value in every frame.

   Two qubits (one generator pair, d1 d2 == 4) and any symbolic reduction keep the exact operator
   square-root route. For d1 d2 > 4 on a numeric reduction, rho.rho-tilde is similar to the positive-
   semidefinite Sqrt[rho].rho-tilde.Sqrt[rho], so its spectrum is real and non-negative: reading it off
   Eigenvalues sidesteps the eigenvector inverse the operator square root Sqrt[rho] (MatrixPower) takes,
   which leaks RowReduce::luc on an ill-conditioned machine reduction and blows up into huge unsimplified
   Root objects on exact input. The eigenvalue route keeps native precision, so machine stays machine and
   exact stays exact (and tractable, unlike the operator square root). *)

wootterCombination[lambda_] := Max[0, 2 Max[lambda] - Total[lambda]] (* lambda1 - Sum[rest]: largest minus the rest, order-independent (the singular values are not guaranteed sorted for symbolic input) *)

ConcurrenceVector[qs_ ? QuantumStateQ, biPartition_ : Automatic] := Block[{
	rho = qs["Bipartition", biPartition]["Operator"], component, d1, d2, y1, y2
},
	{d1, d2} = rho["OutputDimensions"];
	y1 = Y[d1];
	y2 = Y[d2];
	component = If[ d1 d2 > 4 && MatrixQ[rho["Matrix"], NumericQ],
		With[{mat = Normal @ rho["Matrix"]}, With[{cmat = Conjugate[mat]},
			Function[{ya, yb}, With[{omat = Normal @ KroneckerProduct[ya, yb]},
				wootterCombination @ Sqrt @ Clip[Re @ Eigenvalues[mat . omat . cmat . omat], {0, Infinity}]
			]]
		]],
		With[{sr = Sqrt[rho], rc = rho["Conjugate"]},
			Function[{ya, yb}, With[{o = QuantumTensorProduct[QuantumOperator[ya, d1], QuantumOperator[yb, d2]]},
				wootterCombination @ SingularValueList[(sr @ Sqrt[o @ rc @ o])["Matrix"]]
			]]
		]
	];
	Catenate @ Table[component[yn, ym], {yn, y1}, {ym, y2}]
]

Concurrence[qs_ ? QuantumStateQ, biPartition_ : Automatic] := Norm @ ConcurrenceVector[qs, biPartition]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "ConcurrenceVector"] := ConcurrenceVector[qs, biPartition]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Concurrence"] :=
    If[ qs["VectorQ"],
        With[{val = 2 (1 - (QuantumPartialTrace[qs["Bipartition", biPartition], {1}] ^ 2)["Norm"])},
            Sqrt[If[NumericQ[val], Max[0, Re[val]], val]]
        ],
        Concurrence[qs, biPartition]
    ]


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Negativity"] :=
    Enclose[(ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"] - 1) / 2]


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "LogNegativity"] :=
    Enclose @ Log2 @ ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"]


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "EntanglementEntropy"] := Enclose @ With[{
    bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
},
    If[ bp["VectorQ"],
        Quantity[Total[-# Log2[#] & @ Select[Confirm @ bp["SchmidtBasis"]["Probability"], # > 0 &]], "Bits"],
        QuantumPartialTrace[bp, {1}]["VonNeumannEntropy"]
    ]
]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "RenyiEntanglementEntropy" | "RenyiEntropy"] :=
    QuantumEntanglementMonotone[qs, biPartition, {"RenyiEntanglementEntropy", 1 / 2}]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, {"RenyiEntanglementEntropy" | "RenyiEntropy", alpha_}] :=
    Enclose[
        With[{val = (1 / (1 - alpha)) Log[2, Tr @ MatrixPower[
            QuantumPartialTrace[
                ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &],
                {1}
            ]["DensityMatrix"], alpha]]},
            If[NumericQ[val], Re[val], val]
        ]
    ]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Realignment"] :=
    With[{bqs = qs["Bipartition", biPartition]["Normalized"]},
        Total @ SingularValueList @ ArrayReshape[Transpose[bqs["Bend"]["Tensor"], 2 <-> 3], bqs["Dimensions"] ^ 2] - 1
    ]


(* Quantum Discord *)

MutualInformationI[rho_QuantumState, biPartition_ : Automatic] := With[
    {s = rho["Bipartition", biPartition]},

	QuantumPartialTrace[s, {1}]["Entropy"] + QuantumPartialTrace[s, {2}]["Entropy"] - s["Entropy"]
]

MutualInformationJ[rho_QuantumState, qm : _QuantumMeasurementOperator | Automatic : Automatic,biPartition_ : Automatic] := With[
    {s = rho["Bipartition", biPartition]},
    {m = QuantumMeasurementOperator[Replace[qm, Automatic :> QuantumMeasurementOperator[]], {2}][s]}, 

    QuantumPartialTrace[s, {2}]["Entropy"] - m["ProbabilitiesList"] . (Simplify[QuantumPartialTrace[#, {2}]]["Entropy"] & /@ m["States"])
]

QuantumDiscord[rho_QuantumState, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic] :=
	MutualInformationI[rho, biPartition] - MutualInformationJ[rho, qm, biPartition]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "MutualInformationI"] :=
    MutualInformationI[qs, biPartition]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic, "MutualInformationJ"] :=
    MutualInformationJ[qs, qm, biPartition]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic, "Discord"] :=
    QuantumDiscord[qs, qm, biPartition]

