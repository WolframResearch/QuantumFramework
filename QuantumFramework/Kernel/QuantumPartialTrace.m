Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Arrays`"]

PackageExport["QuantumPartialTrace"]



QuantumPartialTrace[arg_, {}] := arg

QuantumPartialTrace[qb_QuditBasis, qudits : {_Integer ..}] := qb["Delete", qudits]


QuantumPartialTrace[qb_QuantumBasis, qudits : {{_Integer, _Integer} ..}] := Enclose @ Module[{},
    ConfirmAssert[Length[qudits] <= Min[qb["OutputRank"], qb["InputRank"]]];
    ConfirmAssert[qb["OutputDimensions"][[qudits[[All, 1]]]] == qb["InputDimensions"][[qudits[[All, -1]]]]];
    QuantumBasis[qb,
        "Output" -> QuantumPartialTrace[qb["Output"], qudits[[All, 1]]],
        "Input" -> QuantumPartialTrace[qb["Input"], qudits[[All, -1]]]
    ]
]

QuantumPartialTrace[qb_QuantumBasis, qudits : {_Integer ..}] := QuantumBasis[qb, "Output" -> QuantumPartialTrace[qb["Output"], qudits]]

QuantumPartialTrace[qs_QuantumState, outQudits : {_Integer ...}, inQudits : {_Integer ...}] :=
    QuantumPartialTrace[qs, Join[outQudits, qs["OutputQudits"] + inQudits]]

QuantumPartialTrace[qs_QuantumState, {}] := qs

QuantumPartialTrace[qs_QuantumState, qudits : {_Integer ..}] :=
    QuantumState[
        MatrixPartialTrace[qs["DensityMatrix"], qudits, qs["Dimensions"]],
        QuantumBasis["Output" -> #1, "Input" -> #2] & @@ QuantumPartialTrace[qs["QuditBasis"], qudits][
            "Split", qs["OutputQudits"] - Count[qudits, q_ /; q <= qs["OutputQudits"]]
        ]
    ]

(* A trace pair contracts a stored output leg against a stored input leg.  A plain index
   contraction sums the diagonal of the stored coefficient matrix C, Tr[C], but the physical trace
   is that of the represented map A = F_out . C . F_in^-1, i.e. Tr[A] = Tr[C . F_in^-1 . F_out].
   Fold the local change of basis M = F_in^-1 . F_out into each traced input leg before contracting;
   when the wire's two sides carry the same basis M is the identity and the leg is left untouched. *)

tracedLegBasisMatrix[qb_, q_] := Normal @ qb["Delete", Complement[Range[qb["Qudits"]], {q}]]["Matrix"]

tracedLegModeMultiply[T_, k_, M_] := With[{perm = Append[DeleteCases[Range[ArrayDepth[T]], k], k]},
    Transpose[Transpose[T, perm] . SparseArray[M], Ordering[perm]]
]

reconcileTracedLegs[qs_, qudits_] := With[{outB = qs["Basis"]["Output"], inB = qs["Basis"]["Input"], oq = qs["OutputQudits"]},
    Fold[
        With[{fOut = tracedLegBasisMatrix[outB, #2[[1]]], fIn = tracedLegBasisMatrix[inB, #2[[2]]]},
            If[fOut === fIn, #1, tracedLegModeMultiply[#1, oq + #2[[2]], Inverse[fIn] . fOut]]
        ] &,
        qs["StateTensor"],
        qudits
    ]
]

QuantumPartialTrace[qs_QuantumState, qudits : {{_Integer, _Integer} ..}] := Enclose[
    ConfirmAssert[DuplicateFreeQ[qudits[[All, 1]]] && DuplicateFreeQ[qudits[[All, 2]]], "a wire can be traced at most once"];
    With[{basis = QuantumPartialTrace[qs["Basis"], qudits]},
        QuantumState[
            If[ qs["VectorQ"],
                ArrayVector @ TensorContract[reconcileTracedLegs[qs, qudits], MapAt[qs["OutputQudits"] + # &, qudits, {All, 2}]],
                ReshapeArray[{TensorContract[qs["DensityTensor"], Join[#, # + qs["Qudits"]] & @ MapAt[qs["OutputQudits"] + # &, qudits, {All, 2}]]}, {#, #} & @ basis["Dimension"]]
            ],
            basis
        ]
    ]
]

QuantumPartialTrace[qs_QuantumState] := QuantumPartialTrace[qs, Range @ qs["Qudits"]]


QuantumPartialTrace[qo_QuantumOperator, qudits : {{_Integer, _Integer} ..}] := Enclose @ With[{
    outputIdx = qo["OutputOrderQuditMapping"],
    inputIdx  = qo["InputOrderQuditMapping"]
},
    QuantumOperator[
        ConfirmBy[QuantumPartialTrace[qo["State"], {Lookup[outputIdx, #[[1]]], Lookup[inputIdx, #[[2]]]} & /@ qudits], QuantumStateQ],
        {
            DeleteElements[qo["OutputOrder"], qudits[[All, 1]]],
            DeleteElements[qo["InputOrder"], qudits[[All, 2]]]
        }
    ]
]

QuantumPartialTrace[qo_QuantumOperator] := QuantumPartialTrace[qo, Intersection @@ qo["Order"]]


QuantumPartialTrace[qm_ ? QuantumMeasurementQ, qudits_] := QuantumMeasurement[QuantumPartialTrace[qm["State"], qudits]]

QuantumPartialTrace[qm_ ? QuantumMeasurementQ] := QuantumPartialTrace[qm, Range[qm["Eigenqudits"]]]


QuantumPartialTrace[qc_ ? QuantumCircuitOperatorQ, qudits : {{_Integer, _Integer} ..}] := Enclose @ Block[{
    outputIdx = qc["OutputOrderQuditMapping"],
    inputIdx  = qc["InputOrderQuditMapping"],
    out = qudits[[All, 1]], in = qudits[[All, 2]],
    min,
    outDims, inDims
},
    ConfirmAssert[ContainsAll[qc["FullOutputOrder"], out] && ContainsAll[qc["FullInputOrder"], in]];
    ConfirmAssert[DuplicateFreeQ[out] && DuplicateFreeQ[in], "a wire can be traced at most once"];
    outDims = qc["OutputDimensions"][[Lookup[outputIdx, out]]];
    inDims = qc["InputDimensions"][[Lookup[inputIdx, in]]];
    ConfirmAssert[outDims == inDims];
    min = Min[Keys[outputIdx], Keys[inputIdx]];
    QuantumCircuitOperator[Reverse @ MapIndexed["Cup"[#1[[2]]] -> {min - #2[[1]], #1[[1]]} &, Thread[{in, inDims}]]] /*
        qc /*
    QuantumCircuitOperator[MapIndexed["Cap"[#1[[2]]] -> {min - #2[[1]], #1[[1]]} &, Thread[{out, outDims}]]]
]

QuantumPartialTrace[qc_QuantumCircuitOperator] := QuantumPartialTrace[qc, Intersection @@ qc["Order"]]


QuantumPartialTrace[op_ ? QuantumFrameworkOperatorQ, qudits : {{_Integer, _Integer} ..}] :=
    Head[op][
        QuantumPartialTrace[op["QuantumOperator"], qudits]
    ]

QuantumPartialTrace[op_ ? QuantumFrameworkOperatorQ, qudits : {_Integer ..}] := QuantumPartialTrace[op, {#, #} & /@ qudits]

