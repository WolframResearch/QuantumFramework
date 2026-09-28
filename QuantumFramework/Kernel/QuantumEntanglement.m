Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumEntangledQ"]
PackageExport["QuantumEntanglementMonotone"]
PackageScope["$QuantumEntanglementMonotones"]
PackageScope["numericStateNotPSDQ"]



$QuantumEntanglementMonotones = {
    "Concurrence", "Negativity", "LogNegativity", "EntanglementEntropy", "RenyiEntropy", "Realignment",
    "MutualInformationI", "MutualInformationJ", "Discord"
}

QuantumEntanglementMonotone::mixedentropy =
    "The entanglement entropy of a subsystem measures entanglement only for a pure state; the reduced von Neumann or Renyi entropy of a mixed state is not an entanglement measure. Use Negativity, LogNegativity, or Concurrence for a mixed state."

QuantumEntanglementMonotone::notphysical =
    "The input is not a positive-semidefinite density matrix; an entanglement monotone of a non-physical state may be meaningless."

warnUnphysicalMonotone[qs_] := If[numericStateNotPSDQ[qs], Message[QuantumEntanglementMonotone::notphysical]]

(* The Chop-ed monotone as a number, or $Failed when it is not numeric (symbolic or undefined). The
   self-contained Enclose keeps the monotone's own Confirm from escaping, so several criteria can be read
   and combined below without a stray Confirm. *)
monotoneValue[qs_, biPartition_, method_] :=
    Enclose[Chop @ ConfirmMatch[QuantumEntanglementMonotone[qs, biPartition, method], _ ? NumericQ], $Failed &]

(* A positive value certifies entanglement; a non-positive one does not certify separability in general,
   and a value that could not be computed is Indeterminate. *)
certifyEntangled[$Failed] := Indeterminate
certifyEntangled[val_] := Positive[val]

(* The default criterion is dimension-aware. For 2 (x) 2 and 2 (x) 3 the negativity (positive partial
   transpose) criterion is necessary and sufficient, so there a False certifies separability. In every
   larger bipartition no single efficiently computable criterion is complete: negativity and the
   realignment (computable cross norm) criterion are complementary, each certifying entangled states the
   other misses (realignment catches positive-partial-transpose bound entangled states negativity cannot),
   so a positive value from either certifies entanglement. An explicitly named method bypasses this and
   uses that monotone alone. *)
QuantumEntangledQ[qs_ ? QuantumStateQ, biPartition_ : Automatic, method : _String | Automatic : Automatic] /;
        method === Automatic || MemberQ[$QuantumEntanglementMonotones, method] :=
    Enclose[
        If[ method =!= Automatic,
            certifyEntangled @ monotoneValue[qs, biPartition, method],
            If[ MatchQ[Sort @ ConfirmMatch[qs["Bipartition", biPartition]["Dimensions"], {_Integer, _Integer}], {2, 2} | {2, 3}],
                certifyEntangled @ monotoneValue[qs, biPartition, "Negativity"],
                With[{neg = monotoneValue[qs, biPartition, "Negativity"], re = monotoneValue[qs, biPartition, "Realignment"]},
                    Which[
                        NumericQ[neg] && Positive[neg], True,
                        NumericQ[re]  && Positive[re],  True,
                        neg === $Failed || re === $Failed, Indeterminate,
                        True, False
                    ]
                ]
            ]
        ],
        Indeterminate &
    ]



(* Default monotone is the concurrence. A method reaches this head either as a bare string in the second
   slot (the optional biPartition then defaults) or as the {"RenyiEntropy", alpha} list; both are excluded
   here so a method spec is never mistaken for a bipartition, which is always Automatic, an integer, or a
   list of qudit indices, never a method name. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ,
    biPartition : Except[_String | {"RenyiEntanglementEntropy" | "RenyiEntropy", _}] : Automatic] :=
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
   and (unlike the two-qubit Wootters concurrence) it is basis-dependent: not invariant under local
   unitaries, though it stays below the LU-invariant convex-roof value in every basis.

   Two qubits (one generator pair, d1 d2 == 4) and any symbolic reduction keep the exact operator
   square-root route. For d1 d2 > 4 on a numeric reduction, rho.rho-tilde is similar to the positive-
   semidefinite Sqrt[rho].rho-tilde.Sqrt[rho], so its spectrum is real and non-negative: reading it off
   Eigenvalues sidesteps the eigenvector inverse the operator square root Sqrt[rho] (MatrixPower) takes,
   which leaks RowReduce::luc on an ill-conditioned machine reduction and blows up into huge unsimplified
   Root objects on exact input. The eigenvalue route keeps native precision, so machine stays machine and
   exact stays exact (and tractable, unlike the operator square root). *)

wootterCombination[lambda_] := Max[0, 2 Max[lambda] - Total[lambda]] (* lambda1 - Sum[rest]: largest minus the rest, order-independent (the singular values are not guaranteed sorted for symbolic input) *)

ConcurrenceVector[qs_ ? QuantumStateQ, biPartition_ : Automatic] := Block[{
	rho = qs["Bipartition", biPartition]["Normalized"]["Operator"], component, d1, d2, y1, y2
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

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "ConcurrenceVector"] :=
    (warnUnphysicalMonotone[qs]; ConcurrenceVector[qs, biPartition])

(* The concurrence depends only on the state's direction, not its overall scale: the bipartition is
   normalized before the reduced purity is read, so an input with trace or vector-norm != 1 gives the
   same value as its normalized form. Every monotone in this file shares that contract (the mixed route
   normalizes inside ConcurrenceVector). Without it the reduced purity picks up the scale, and the
   Max[0, ...] clamp on the now-shifted 2 (1 - Purity) silently returns 0 for a scaled pure state. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Concurrence"] := (
    warnUnphysicalMonotone[qs];
    If[ qs["VectorQ"],
        With[{val = 2 (1 - (QuantumPartialTrace[qs["Bipartition", biPartition]["Normalized"], {1}] ^ 2)["Norm"])},
            Sqrt[If[NumericQ[val], Max[0, Re[val]], val]]
        ],
        Concurrence[qs, biPartition]
    ]
)


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Negativity"] := (
    warnUnphysicalMonotone[qs];
    Enclose[(ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"] - 1) / 2]
)


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "LogNegativity"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ Log2 @ ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"]
)


(* Entanglement entropy is the von Neumann entropy of a reduced state, an entanglement measure only for
   a pure global state. A pure state reads it through its Schmidt weights (vector) or its reduced state
   (density-matrix form); a genuinely mixed state's reduced von Neumann entropy mixes classical ignorance
   into the count, so it is not returned as entanglement. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "EntanglementEntropy"] := Enclose @ With[{
    bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
},
    warnUnphysicalMonotone[qs];
    Which[
        bp["VectorQ"],
            Quantity[Total[-# Log2[#] & @ Select[Confirm @ bp["SchmidtBasis"]["Probability"], If[NumericQ[#], # > 0, True] &]], "Bits"],
        TrueQ[bp["PureStateQ"]],
            QuantumPartialTrace[bp, {1}]["VonNeumannEntropy"],
        True,
            Message[QuantumEntanglementMonotone::mixedentropy]; Indeterminate
    ]
]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "RenyiEntanglementEntropy" | "RenyiEntropy"] :=
    QuantumEntanglementMonotone[qs, biPartition, {"RenyiEntanglementEntropy", 1 / 2}]

(* The Renyi entanglement entropy is the Renyi entropy of the reduced state, an entanglement measure only
   for a pure global state, so a genuinely mixed input is guarded exactly as EntanglementEntropy is. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, {"RenyiEntanglementEntropy" | "RenyiEntropy", alpha_}] :=
    Enclose @ With[{
        bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
    },
        warnUnphysicalMonotone[qs];
        If[ bp["VectorQ"] || TrueQ[bp["PureStateQ"]],
            With[{val = (1 / (1 - alpha)) Log[2, Tr @ MatrixPower[QuantumPartialTrace[bp, {1}]["DensityMatrix"], alpha]]},
                If[NumericQ[val], Re[val], val]
            ],
            Message[QuantumEntanglementMonotone::mixedentropy]; Indeterminate
        ]
    ]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Realignment"] :=
    With[{bqs = qs["Bipartition", biPartition]["Normalized"]},
        warnUnphysicalMonotone[qs];
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
    (warnUnphysicalMonotone[qs]; MutualInformationI[qs, biPartition])

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic, "MutualInformationJ"] :=
    (warnUnphysicalMonotone[qs]; MutualInformationJ[qs, qm, biPartition])

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic, "Discord"] :=
    (warnUnphysicalMonotone[qs]; QuantumDiscord[qs, qm, biPartition])

