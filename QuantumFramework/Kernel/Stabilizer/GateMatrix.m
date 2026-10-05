Package["Wolfram`QuantumFramework`"]

PackageScope[stabilizerCliffordRows]
PackageScope[stabilizerCircuitSpecs]
PackageScope[stabilizerOperatorSpecs]
PackageScope[stabilizerShortcutPairs]
PackageScope[stabilizerLabelsHonestQ]
PackageScope[stabilizerArrayZeroQ]
PackageScope[stabilizerRefusalReason]
PackageScope[stabilizerMatrixRows]
PackageScope[stabilizerClearCaches]
PackageScope[$stabilizerMatrixMaxQubits]



(* ============================================================================ *)
(* Which operation a gate performs is read from its matrix, never from its     *)
(* label: a label is display metadata a user can set to anything.              *)
(*                                                                              *)
(* A gate U on k qubits is a Clifford exactly when conjugation maps each       *)
(* generator P = X_j, Z_j of its wires to a signed Pauli string G = U P U^dag.  *)
(* G has one nonzero entry per row: G e_b = c (-1)^(z.b) e_(b xor x), with     *)
(* c = s i^(x.z) and s = +1 or -1. So G applied to |0> gives x and c, and G    *)
(* applied to the k one-bit states gives z: k + 1 matrix-vector products per    *)
(* generator. The candidate (x, z, s) is then certified entry by entry as      *)
(* U P = G U, which never forms G. When every generator passes, U^dag U         *)
(* commutes with every X_j and Z_j, so it is a multiple of the identity, and    *)
(* the norm of one column decides unitarity; that norm is checked first. The   *)
(* images are the gate's tableau rows, in the row format of FromFullTableau.   *)
(* The cost grows as 4^k.                                                      *)
(*                                                                              *)
(* Candidates are read with machine numbers and screened loosely, the bound    *)
(* growing with the size and tolerance of the matrix; certification uses the  *)
(* matrix's own arithmetic. Exact entries are decided exactly: a machine screen,*)
(* confirmed by RootReduce (PossibleZeroQ and high-precision N warn on          *)
(* unsimplified algebraic zeros), and symbolic entries by Simplify. Approximate *)
(* entries are decided to 10 units of their last digit, never finer than        *)
(* 10^-9. That digit is read from Accuracy, because a computed zero such as     *)
(* 0``29.5 has Precision 0 and a matrix has the precision of its worst entry.  *)
(* ============================================================================ *)

(* The Clifford test stops at 12 qubits: a 12-qubit matrix holds 4^12 entries.  *)
(* A diagonal matrix on more wires is still read when it factors into one-qubit *)
(* phases, which costs only its 2^k diagonal entries.                           *)
$stabilizerMatrixMaxQubits = 12;

(* From these sizes the Clifford test of one gate takes noticeable time, four   *)
(* times more per added qubit, and PauliStabilizer::largegate says so: 8 qubits *)
(* for exact entries, 10 for approximate ones.                                 *)
$stabilizerLargeGateQubits = <|True -> 8, False -> 10|>;

(* Accuracy is capped at 10, where the 10^-9 floor takes over: an array of machine zeros has *)
(* Accuracy near 323, and 10^-323 underflows (General::munfl).                              *)
stabilizerNumericSpec[u_] := If[Precision[u] === Infinity, {True, 0}, {False, Max[10^-9, 10 10^-Min[Accuracy[u], 10]]}]

stabilizerZeroQ[{True, _}][x_] := TrueQ[Abs[N[x]] < 10^-8] && (RootReduce[x] === 0 || Simplify[x] === 0)
stabilizerZeroQ[{False, tol_}][x_] := TrueQ[Abs[x] < tol]

(* Every entry of an array of differences is zero. Exact zeros are skipped, so a  *)
(* difference of two equal exact arrays costs no zero tests, and approximate      *)
(* entries are compared with the tolerance in one step. Symbolic entries are       *)
(* simplified, with Chop for the 0. a machine-precision symbolic entry leaves. The *)
(* background of a SparseArray counts as an entry.                                *)
stabilizerArrayZeroQ[array_] := With[{sparse = stabilizerSparse[array]},
    With[{values = DeleteCases[Append[sparse["NonzeroValues"], sparse["Background"]], 0]},
        Which[
            ! AllTrue[values, NumericQ], AllTrue[values, TrueQ[Chop[Simplify[#]] == 0] &],
            Precision[values] === Infinity, AllTrue[values, stabilizerZeroQ[{True, 0}]],
            True, TrueQ[Max[Abs[values], 0] < Last[stabilizerNumericSpec[values]]]
        ]
    ]
]

(* SparseArray[m] rescans every entry of a SparseArray m, seconds for an exact *)
(* matrix on 11 qubits, so an existing SparseArray is used as it is.          *)
stabilizerSparse[m_SparseArray] := m
stabilizerSparse[m_] := SparseArray[m]

(* a === b decides equal exact arrays at once, where subtracting them costs    *)
(* seconds on 11 qubits; the zero test of the difference decides the rest.     *)
stabilizerSameArrayQ[a_, b_] := a === b || stabilizerArrayZeroQ[a - b]

stabilizerNumericArrayQ[m_] := With[{sparse = stabilizerSparse[m]}, AllTrue[sparse["NonzeroValues"], NumericQ] && NumericQ[sparse["Background"]]]

stabilizerExactQ[m_] := Precision[m] === Infinity

(* |x| = 1 shown for symbolic x, its symbols taken as real (ComplexExpand). *)
stabilizerUnitModulusQ[x_] := TrueQ[Simplify[ComplexExpand[Abs[x] ^ 2] == 1]]

(* An exact matrix divided by the phase of its first entry that is not zero: the *)
(* Clifford test does not see a global phase, and an algebraic phase such as     *)
(* (-1)^(1/4) left in every entry makes its exact zero tests slow.               *)
stabilizerStripPhase[u_] := If[ stabilizerExactQ[u],
    Replace[SelectFirst[stabilizerSparse[u]["NonzeroValues"], TrueQ[Abs[N[#]] > 10^-8] &], {p_ ? NumericQ :> u Abs[p] / p, _ :> u}],
    u
]

(* Qubit 1 is the most significant bit of a basis index. *)
stabilizerBit[k_, j_] := 2 ^ (k - j)
stabilizerBitVector[k_, j_] := BitAnd[Quotient[Range[0, 2 ^ k - 1], stabilizerBit[k, j]], 1]

(* P v for a vector v, and u P for a matrix u, with P = X_j or Z_j. *)
stabilizerLeftGenerator[{"X", j_}, k_, v_] := v[[BitXor[Range[0, 2 ^ k - 1], stabilizerBit[k, j]] + 1]]
stabilizerLeftGenerator[{"Z", j_}, k_, v_] := (1 - 2 stabilizerBitVector[k, j]) v
stabilizerRightGenerator[{"X", j_}, k_, u_] := u[[All, BitXor[Range[0, 2 ^ k - 1], stabilizerBit[k, j]] + 1]]
stabilizerRightGenerator[{"Z", j_}, k_, u_] := Transpose[(1 - 2 stabilizerBitVector[k, j]) Transpose[u]]

(* The bound a candidate's coefficient c is screened against. Certification decides, *)
(* so the bound is loose: c sums 2^k products of entries, and its error grows as     *)
(* 2^(k/2) times the tolerance of the entries.                                       *)
stabilizerScreen[u_, k_] := Min[1 / 2, Max[10^-6, 2 ^ (k / 2 + 1) Last[stabilizerNumericSpec[u]]]]

(* The candidate {x, z, s} for G = un P un^dag, from G applied to |0> and the one-bit *)
(* states (probes holds un^dag e_b as rows), or Missing when none fits.              *)
stabilizerImageCandidate[un_, generator_, k_, probes_, screen_] := With[{
    images = Normal[un . Transpose[stabilizerLeftGenerator[generator, k, #] & /@ probes]]
},
    With[{x0 = First[Ordering[Abs[images[[All, 1]]], -1]] - 1},
        With[{c = images[[x0 + 1, 1]], x = IntegerDigits[x0, 2, k]},
            If[ Abs[Abs[c] - 1] >= screen,
                Missing["NotClifford"],
                With[{z = Table[Boole[Re[images[[BitXor[stabilizerBit[k, i], x0] + 1, i + 1]] / c] < 0], {i, k}]},
                    With[{s = c / I ^ (x . z)},
                        If[Abs[s - Round[Re[s]]] < screen, {x, z, Round[Re[s]]}, Missing["NotClifford"]]
                    ]
                ]
            ]
        ]
    ]
]

(* u P = G u entry by entry: row r of G u is s i^(x.z) (-1)^(z.(r xor x)) times row r xor x of u. *)
stabilizerImageCertifiedQ[u_, generator_, k_, {x_, z_, s_}] := With[{rx = BitXor[Range[0, 2 ^ k - 1], FromDigits[x, 2]]},
    stabilizerArrayZeroQ[stabilizerRightGenerator[generator, k, u] - s I ^ (x . z) (1 - 2 Mod[IntegerDigits[rx, 2, k] . z, 2]) u[[rx + 1]]]
]

(* For a refusal's reason only: un keeps, to the screen, the norms of its columns and  *)
(* of two fixed pseudo-random vectors. A matrix that is not unitary changes the norm  *)
(* of all vectors but a set of measure zero.                                          *)
stabilizerNormPreservingQ[un_, k_, screen_] := With[{vs = BlockRandom[SeedRandom[1]; RandomComplex[{-1 - I, 1 + I}, {2, 2 ^ k}]]},
    TrueQ[Max[Abs[Total[Abs[un] ^ 2] - 1]] <= screen] && AllTrue[vs, TrueQ[Abs[Norm[un . #] / Norm[#] - 1] <= screen] &]
]

stabilizerGeneratorRow[u_, un_, generator_, k_, probes_, screen_] := Replace[stabilizerImageCandidate[un, generator, k, probes, screen], {
    candidate : {x_, z_, s_} /; stabilizerImageCertifiedQ[u, generator, k, candidate] :> Join[x, z, {Boole[s == -1]}],
    _ :> Missing["NotClifford"]
}]

(* Tableau rows {x bits, z bits, sign bit} of the images of X_1..X_k, Z_1..Z_k under the *)
(* 2^k x 2^k matrix m, or Missing with a reason. m may be a SparseArray.                *)
stabilizerCliffordRows[m_] := With[{k = Log2[Length[m]]},
    Which[
        ! IntegerQ[k] || Dimensions[m] =!= {2 ^ k, 2 ^ k}, Missing["NotQubit"],
        ! stabilizerNumericArrayQ[m], Missing["Symbolic"],
        k > $stabilizerMatrixMaxQubits, Missing["TooLarge"],
        True, With[{u = stabilizerStripPhase[m]}, With[{un = N[u], screen = stabilizerScreen[u, k]},
            If[ ! stabilizerArrayZeroQ[{Total[Abs[u[[All, 1]]] ^ 2] - 1}],
                Missing["NotUnitary"],
                With[{
                    probes = Normal[Conjugate[un[[Prepend[stabilizerBit[k, #] & /@ Range[k], 0] + 1]]]],
                    generators = Join[Table[{"X", j}, {j, k}], Table[{"Z", j}, {j, k}]]
                },
                    Replace[
                        Enclose[ConfirmBy[stabilizerGeneratorRow[u, un, #, k, probes, screen], ListQ] & /@ generators, Missing["NotClifford"] &],
                        _Missing :> If[stabilizerNormPreservingQ[un, k, screen], Missing["NotClifford"], Missing["NotUnitary"]]
                    ]
                ]
            ]
        ]]
    ]
]


(* ============================================================================ *)
(* From the matrix to specs the engine applies, on the gate's sorted wires     *)
(* 1..k (the matrix of an operator is indexed by its sorted wires).            *)
(*                                                                              *)
(*   diagonal product of one-qubit phases -> each factor on its own wire, on   *)
(*                                   any number of wires, numeric or symbolic  *)
(*   Clifford on up to 12 qubits  -> the named gates of its own tableau's      *)
(*                                   decomposition (ps["Circuit"])             *)
(*   one-qubit diagonal           -> "P"[phase], the phase read from the matrix *)
(*                                   (symbolic: when both diagonal entries     *)
(*                                   have modulus 1)                           *)
(*   product of one-qubit gates   -> each factor read as a one-qubit gate on    *)
(*                                   its own wire, up to 12 qubits             *)
(*   anything else                -> Missing[reason]                           *)
(*                                                                              *)
(* The decomposition is rewritten into the forms the per-gate fold has rules   *)
(* for: daggers of self-inverse gates dropped. The controlled form             *)
(* "C"["NOT" -> {t}, {c}, {}] that QuantumShortcut writes is also written as   *)
(* "CNOT" -> {c, t}, the form of a named circuit; the fold and the compiled    *)
(* encoder read both.                                                          *)
(* ============================================================================ *)

$stabilizerEngineForm = {
    (SuperDagger[g : "H" | "X" | "Y" | "Z" | "SWAP"] -> wires_) :> (g -> wires),
    "C"[g : "NOT" | "X" | "Z" -> target_List, control_List, {}] :> (("C" <> g) -> Join[control, target])
};

stabilizerCliffordSpecs[{}] := {}
stabilizerCliffordSpecs[rows_List] := With[{
    specs = Replace[QuantumShortcut[FromFullTableau[rows]["Circuit"]], $stabilizerEngineForm, {1}]
},
    If[MatchQ[specs, {(_ -> {___Integer}) ...}], specs, Missing["UnexpectedToken"]]
]

stabilizerDiagonalPhase[u_, spec_] := Which[
    Length[u] =!= 2 || ! AllTrue[{u[[1, 2]], u[[2, 1]]}, stabilizerZeroQ[spec]], Missing["NotClifford"],
    ! AllTrue[{Abs[u[[1, 1]]] - 1, Abs[u[[2, 2]]] - 1}, stabilizerZeroQ[spec]], Missing["NotUnitary"],
    (* plain Arg emits N::meprec on an exact phase written as Cos[a] + I Sin[a] *)
    True, Arg[TrigToExp[u[[2, 2]] / u[[1, 1]]]]
]

(* A symbolic one-qubit gate is accepted only when its off-diagonal entries are zero *)
(* and both diagonal entries have modulus 1. "P"[phase] applies Exp[I phase] =        *)
(* u22 / u11 exactly: PowerExpand keeps that identity for any branch of Log, and Chop *)
(* clears the 0. parts a machine-precision matrix leaves in the phase.               *)
stabilizerSymbolicPhase[u_] := If[
    Length[u] == 2 && TrueQ[u[[1, 2]] == 0] && TrueQ[u[[2, 1]] == 0] && stabilizerUnitModulusQ[u[[1, 1]]] && stabilizerUnitModulusQ[u[[2, 2]]],
    Chop[- I PowerExpand[Log[u[[2, 2]] / u[[1, 1]]]]],
    Missing["Symbolic"]
]

stabilizerOneQubitSpecs[u_] := With[{rows = stabilizerCliffordRows[u]},
    Which[
        ListQ[rows], stabilizerCliffordSpecs[rows],
        rows === Missing["NotClifford"], Replace[stabilizerDiagonalPhase[u, stabilizerNumericSpec[u]], phase : Except[_Missing] :> {"P"[phase] -> {1}}],
        rows === Missing["Symbolic"], Replace[stabilizerSymbolicPhase[u], phase : Except[_Missing] :> {"P"[phase] -> {1}}],
        True, rows
    ]
]

(* m may be a SparseArray; label only names the gate in PauliStabilizer::largegate. *)
stabilizerMatrixSpecs[m_, label_ : None] := Which[
    Length[m] == 1, With[{x = Normal[m][[1, 1]]},
        If[If[NumericQ[x], stabilizerArrayZeroQ[{Abs[x] - 1}], stabilizerUnitModulusQ[x]], {}, Missing["NotUnitary"]]
    ],
    Length[m] == 2, stabilizerOneQubitSpecs[Normal[m]],
    True, stabilizerManyQubitSpecs[m, Log2[Length[m]], label]
]

(* The one-qubit factors of a product, or Missing: a diagonal matrix on any number of *)
(* wires, a numeric one on up to $stabilizerMatrixMaxQubits.                          *)
stabilizerFactors[m_, k_Integer, diagonalQ_] := Which[
    diagonalQ, stabilizerDiagonalFactors[Normal[Diagonal[stabilizerSparse[m]]], k],
    stabilizerNumericArrayQ[m] && k <= $stabilizerMatrixMaxQubits, stabilizerTensorFactors[m, k],
    True, Missing["NotProduct"]
]

(* PauliStabilizer::largegate, before the full Clifford test of a matrix that is not diagonal. *)
stabilizerLargeGateWarning[m_, k_, diagonalQ_, label_] :=
    If[! diagonalQ && k >= $stabilizerLargeGateQubits[stabilizerExactQ[m]], Message[PauliStabilizer::largegate, k, Replace[label, None -> "with no label"]]]

(* The factors of a matrix on several wires and the rows of its Clifford test are kept *)
(* per distinct matrix (its Hash; see the caches below): the PauliStabilizer           *)
(* constructor reads a gate through QuantumOperatorTableau and, when that refuses, again *)
(* through the gate fold. PauliStabilizer::largegate comes with the Clifford test, so   *)
(* it fires once per matrix.                                                           *)
$stabilizerMatrixMemo = <||>;

SetAttributes[stabilizerMatrixRemember, HoldRest]
stabilizerMatrixRemember[tag_, key_, value_] := With[{full = Hash[{tag, key}]},
    Lookup[$stabilizerMatrixMemo, full, stabilizerRemember[$stabilizerMatrixMemo, full, value]]
]

stabilizerTestedRows[m_, k_, diagonalQ_, label_, key_] :=
    stabilizerMatrixRemember["Rows", key, (stabilizerLargeGateWarning[m, k, diagonalQ, label]; stabilizerCliffordRows[m])]

stabilizerManyQubitSpecs[m_, k_Integer, label_] := With[{diagonalQ = stabilizerStructurallyDiagonalQ[m], key = Hash[m]},
    With[{factors = stabilizerMatrixRemember["Factors", key, stabilizerFactors[m, k, diagonalQ]]},
        Which[
            (* a product is read factor by factor, and refused when a factor is refused *)
            ListQ[factors], stabilizerProductSpecs[factors],
            ! stabilizerNumericArrayQ[m], Replace[factors, Missing["NotProduct" | "NotClifford"] -> Missing["Symbolic"]],
            (* a diagonal entry, or a product's overall scale, off modulus 1 *)
            factors === Missing["NotUnitary"], factors,
            k > $stabilizerMatrixMaxQubits, Missing["TooLarge"],
            True, Replace[stabilizerTestedRows[m, k, diagonalQ, label, key], rows_List :> stabilizerCliffordSpecs[rows]]
        ]
    ]
]
stabilizerManyQubitSpecs[_, _, _] := Missing["NotQubit"]

(* Tableau rows of a qubit matrix, for QuantumOperatorTableau: a product of one-qubit *)
(* Cliffords row by row from its factors, anything else by the full Clifford test.     *)
stabilizerMatrixRows[m_, label_ : None] := With[{k = Log2[Length[m]]},
    Which[
        ! IntegerQ[k] || k < 1, Missing["NotQubit"],
        k == 1, stabilizerCliffordRows[Normal[m]],
        ! stabilizerNumericArrayQ[m], Missing["Symbolic"],
        True, With[{diagonalQ = stabilizerStructurallyDiagonalQ[m], key = Hash[m]},
            With[{factors = stabilizerMatrixRemember["Factors", key, stabilizerFactors[m, k, diagonalQ]]},
                Which[
                    ListQ[factors], stabilizerProductRows[factors],
                    factors === Missing["NotUnitary"], factors,
                    k > $stabilizerMatrixMaxQubits, Missing["TooLarge"],
                    True, stabilizerTestedRows[m, k, diagonalQ, label, key]
                ]
            ]
        ]
    ]
]

(* The rows of a product: the image of X_j or Z_j is its factor's image on wire j. *)
stabilizerProductRows[factors_List] := With[{k = Length[factors], rows = stabilizerCliffordRows /@ factors},
    If[ FreeQ[rows, _Missing],
        Join @@ Table[
            Table[With[{r = rows[[j, generator]]}, Join[ReplacePart[ConstantArray[0, k], j -> r[[1]]], ReplacePart[ConstantArray[0, k], j -> r[[2]]], {r[[3]]}]], {j, k}],
            {generator, 2}
        ],
        Missing["NotClifford"]
    ]
]


(* ============================================================================ *)
(* An operator on several wires that is a product of one-qubit gates (QF builds *)
(* "H" -> Range[n] and "T" -> {1, 2} as one operator) is read factor by factor, *)
(* each factor as a one-qubit gate on its own wire.                             *)
(*                                                                              *)
(* A diagonal matrix factors when every diagonal entry d(x) is d(0) times the  *)
(* phases r_j = d(e_j) / d(0) of the wires j set in x; this is decided for     *)
(* symbolic entries too. The factors are diag(1, r_j), and d(0) must have      *)
(* modulus 1. A numeric matrix that is not diagonal factors when its entries,  *)
(* regrouped wire by wire into a k-index array of 2 x 2 blocks, form an outer  *)
(* product. The blocks are read off through one entry that is not zero, the    *)
(* pivot: the block of wire j holds the entries whose row and column differ    *)
(* from the pivot's only in bit j. Their Kronecker product is compared with    *)
(* the matrix, and each block is scaled by the norm of a column; the product   *)
(* of the scales must have modulus 1.                                          *)
(* ============================================================================ *)

(* Every entry off the diagonal is zero (machine zeros included). *)
stabilizerStructurallyDiagonalQ[m_] := With[{sparse = stabilizerSparse[m]},
    TrueQ[sparse["Background"] == 0] && With[{positions = sparse["NonzeroPositions"]},
        positions === {} || AllTrue[Pick[sparse["NonzeroValues"], Unitize[positions[[All, 1]] - positions[[All, 2]]], 1], TrueQ[# == 0] &]
    ]
]

stabilizerDiagonalFactors[d_List, k_Integer] := With[{d0 = First[d]},
    If[ If[NumericQ[d0], stabilizerArrayZeroQ[{Abs[d0] - 1}], stabilizerUnitModulusQ[d0]],
        With[{ratios = d[[2 ^ (k - Range[k]) + 1]] / d0},
            If[ stabilizerArrayZeroQ[d - d0 Flatten[Outer[Times, Sequence @@ ({1, #} & /@ ratios)]]],
                DiagonalMatrix[{1, #}] & /@ ratios,
                Missing["NotProduct"]
            ]
        ],
        Missing["NotUnitary"]
    ]
]

(* The pivot is the largest entry for approximate entries, and for exact ones the *)
(* first entry whose numeric value is not zero (an exact entry can be a zero that  *)
(* is not simplified). Exact blocks are compared exactly.                          *)
stabilizerTensorFactors[m_, k_Integer] := With[{u = Normal[m]},
    Replace[
        If[ stabilizerExactQ[u],
            FirstPosition[u, x_ /; TrueQ[Abs[N[x]] > 10^-8], Missing["NotUnitary"], {2}, Heads -> False],
            With[{p = QuotientRemainder[First[Ordering[Abs[Flatten[u]], -1]] - 1, Length[u]] + 1},
                If[TrueQ[Abs[Extract[u, p]] > 10^-8], p, Missing["NotUnitary"]]
            ]
        ],
        p : {_Integer, _Integer} :> stabilizerTensorBlocks[u, p, k]
    ]
]

(* The blocks through the pivot at {r, c}, the first unscaled and the others divided *)
(* by the pivot, checked against u and scaled to unitaries.                          *)
stabilizerTensorBlocks[u_, {r_, c_}, k_] := With[{
    blocks = Table[
        With[{bit = stabilizerBit[k, j]},
            u[[BitAnd[r - 1, BitNot[bit]] + {1, bit + 1}, BitAnd[c - 1, BitNot[bit]] + {1, bit + 1}]] / If[j == 1, 1, u[[r, c]]]
        ],
        {j, k}
    ]
},
    With[{scales = Norm[First[SortBy[Transpose[#], - Abs[N[Norm[#]]] &]]] & /@ blocks},
        Which[
            ! stabilizerSameArrayQ[u, KroneckerProduct @@ blocks], Missing["NotProduct"],
            ! stabilizerArrayZeroQ[{Abs[Times @@ scales] - 1}], Missing["NotUnitary"],
            True, blocks / scales
        ]
    ]
]

stabilizerProductSpecs[factors_List] := With[{specs = stabilizerOneQubitSpecs /@ factors},
    If[ FreeQ[specs, _Missing],
        Catenate @ MapIndexed[Function[{one, j}, Replace[one, (g_ -> w_List) :> (g -> w + First[j] - 1), {1}]], specs],
        First[Cases[specs, _Missing]]
    ]
]
stabilizerProductSpecs[missing_Missing] := missing


(* ============================================================================ *)
(* Per operator, with the classification cached per distinct gate.             *)
(*                                                                              *)
(* An operator is QuantumOperator[QuantumState[amplitudes, basis], order]. Its *)
(* classification depends only on the stored amplitudes, the basis parts and   *)
(* the order patterns of its input and output wires, so it is computed once    *)
(* per distinct gate and kept across calls: a circuit applied many times reads *)
(* each gate once, and PauliStabilizer::largegate fires once per gate. The     *)
(* label is left out of the key: a controlled gate's label names its control   *)
(* wire, and a label decides nothing here. A QuantumBasis is atomic once       *)
(* validated, so its parts are read by property. A cache is emptied when it    *)
(* reaches 4096 entries.                                                       *)
(*                                                                              *)
(* Refused before the matrix is read: a wire of dimension other than 2 (the    *)
(* tableau is binary, and a 4-level wire would otherwise read as two qubits),  *)
(* an operator whose output wires are not its input wires (it moves a qubit    *)
(* rather than acting on fixed wires), and an operator stored in a basis that  *)
(* is not computational. For H re-expressed in the "PauliX" basis the           *)
(* "Computational" matrix is H, yet QF applies it to |0> as |0>, so the matrix *)
(* does not say what QF does with such an operator.                            *)
(*                                                                              *)
(* A composite operator (QuantumOperator["T"] @ QuantumOperator["H"]) whose     *)
(* matrix is none of the above is read through the gates its label names. They *)
(* are used only when QF's own composition of them reproduces the operator's   *)
(* matrix up to a global phase, and each is then read from its own matrix.      *)
(* A symbolic matrix the steps above cannot read ("P"[ArcCos[x]], whose phase  *)
(* is not provably real) is applied as the one gate its label names, under the *)
(* same condition. This step depends on the label, so its cache key includes   *)
(* the label.                                                                  *)
(* ============================================================================ *)

$stabilizerGateMemo = <||>;
$stabilizerTokenMemo = <||>;
$stabilizerHonestMemo = <||>;

stabilizerClearCaches[] := ($stabilizerGateMemo = <||>; $stabilizerTokenMemo = <||>; $stabilizerHonestMemo = <||>; $stabilizerMatrixMemo = <||>;)

SetAttributes[stabilizerRemember, HoldAll]
stabilizerRemember[memo_Symbol, key_, value_] := (If[Length[memo] >= 4096, memo = <||>]; memo[key] = value)

stabilizerGateKey[op_QuantumOperator] := With[{basis = op[[1, 2]]},
    Hash[{op[[1, 1]], basis["Input"], basis["Output"], basis["Picture"], Ordering[op["InputOrder"]], Ordering[op["OutputOrder"]]}]
]

stabilizerOperatorClassification[op_QuantumOperator] := Which[
    ! AllTrue[op["Dimensions"], # == 2 &], Missing["NotQubit"],
    ! TrueQ[op[[1, 2]]["ComputationalQ"]], Missing["NonComputationalBasis"],
    True, stabilizerMatrixSpecs[op["Sort"]["Matrix"], op["Label"]]
]

(* The same matrix up to a global phase, read at the first entry of b that is not zero. *)
stabilizerSameUpToPhaseQ[a_, b_] := Dimensions[a] === Dimensions[b] && (a === b || With[{sa = stabilizerSparse[a], sb = stabilizerSparse[b]},
    With[{i = SelectFirst[Range[Length[sb["NonzeroValues"]]], TrueQ[Abs[N[sb["NonzeroValues"][[#]]]] > 10^-8] &]},
        IntegerQ[i] && With[{phase = Extract[sa, sb["NonzeroPositions"][[i]]] / sb["NonzeroValues"][[i]]},
            If[NumericQ[phase], stabilizerArrayZeroQ[{Abs[phase] - 1}], stabilizerUnitModulusQ[phase]] && stabilizerArrayZeroQ[sa - phase sb]
        ]
    ]
])

(* Several tokens are each read from their own matrix; the one token of a symbolic *)
(* gate is applied as it is, since its matrix is the operator's. That token must   *)
(* name a gate: an operator with no label has a token holding its own matrix.      *)
stabilizerTokenSpecs[op_QuantumOperator, wires_List, reason_] := With[{tokens = QuantumShortcut[op]},
    If[ ListQ[tokens] && (Length[tokens] >= 2 || reason === "Symbolic" && MatchQ[tokens, {(_String | _String[___]) -> {___Integer}} | {_String[___]}]) &&
            Length[wires] <= $stabilizerMatrixMaxQubits,
        With[{composed = QuantumCircuitOperator[Prepend[tokens, "I" -> wires]]["QuantumOperator"]},
            If[ Head[composed] === QuantumOperator && Sort[composed["InputOrder"]] === wires && Sort[composed["OutputOrder"]] === wires &&
                    stabilizerSameUpToPhaseQ[op["Sort"]["Matrix"], composed["Sort"]["Matrix"]],
                If[Length[tokens] == 1, Replace[tokens, $stabilizerEngineForm, {1}], stabilizerElementSpecs[QuantumCircuitOperator[tokens]]],
                Missing["NotClifford"]
            ]
        ],
        Missing["NotClifford"]
    ]
]

stabilizerOperatorSpecs[op_QuantumOperator] := With[{in = op["InputOrder"], out = op["OutputOrder"]},
    If[ Sort[in] =!= Sort[out],
        {Missing["MovesWires", op]},
        With[{key = stabilizerGateKey[op], wires = Sort[in]},
            With[{local = Lookup[$stabilizerGateMemo, key, stabilizerRemember[$stabilizerGateMemo, key, stabilizerOperatorClassification[op]]]},
                Which[
                    ListQ[local], Replace[local, (g_ -> w_List) :> (g -> wires[[w]]), {1}],
                    MatchQ[local, Missing["NotClifford" | "NotProduct" | "Symbolic"]], With[{tokenKey = Hash[{key, op["Label"], in}]},
                        Replace[
                            Lookup[$stabilizerTokenMemo, tokenKey, stabilizerRemember[$stabilizerTokenMemo, tokenKey, stabilizerTokenSpecs[op, wires, First[local]]]],
                            _Missing :> {Missing[First[local], op]}
                        ]
                    ],
                    True, {Missing[First[local], op]}
                ]
            ]
        ]
    ]
]

(* A circuit's specs for the engine: unitary operators by their matrices, nested *)
(* circuits element by element, and every other element (a measurement, a channel)*)
(* by its QuantumShortcut token, which the fold refuses.                         *)
stabilizerCircuitSpecs[qco_QuantumCircuitOperator] := stabilizerElementSpecs[qco]

stabilizerElementSpecs[qco_QuantumCircuitOperator] := Catenate[stabilizerElementSpecs /@ qco["Operators"]]
stabilizerElementSpecs[op_QuantumOperator] := stabilizerOperatorSpecs[op]
stabilizerElementSpecs[other_] := QuantumShortcut[other]

(* Why a gate was refused, for PauliStabilizer::nonclifford. *)
stabilizerRefusalReason["NotUnitary"] := "its matrix is not unitary"
stabilizerRefusalReason["NotClifford" | "NotProduct"] := "its matrix is not a Clifford, a product of one-qubit Cliffords and phase gates, or the product of the gates its label names"
stabilizerRefusalReason["Symbolic"] := "its matrix is symbolic, and is neither a phase gate whose entries have modulus 1 nor the matrix of the gates its label names"
stabilizerRefusalReason["TooLarge"] := "its matrix acts on more than " <> ToString[$stabilizerMatrixMaxQubits] <> " qubits and is not a diagonal product of one-qubit phases"
stabilizerRefusalReason["NotQubit"] := "it acts on a wire that is not a qubit"
stabilizerRefusalReason["MovesWires"] := "it moves a qubit to another wire"
stabilizerRefusalReason["NonComputationalBasis"] := "it is stored in a basis that is not computational"
stabilizerRefusalReason["UnexpectedToken"] := "its matrix decomposes into a gate the engine has no rule for"
stabilizerRefusalReason[reason_] := ToString[reason]


(* ============================================================================ *)
(* Label checks for the phase-polynomial build, which reads gates from their   *)
(* QuantumShortcut tokens: T, S, Z, their daggers, controlled-Z and CNOT.      *)
(*                                                                              *)
(* stabilizerShortcutPairs lists QuantumShortcut[qco] element by element, as   *)
(* {element, its tokens}, nested circuits flattened in order, so the tokens    *)
(* joined are QuantumShortcut[qco]. stabilizerLabelsHonestQ is True when every *)
(* operator's matrix is the product of its tokens (exactly for an exact matrix, *)
(* within the zero test's tolerance for an approximate one), each a gate of    *)
(* that set on the operator's own wires ("S" -> Range[3] is three tokens, and  *)
(* QuantumOperator["T"] @ QuantumOperator["T"] is two on one wire). The tokens' *)
(* matrices are written out on the local wires 1..k, qubit 1 the most          *)
(* significant factor, as sparse arrays. The verdict depends on the tokens     *)
(* only through their wires' relative order, so it is cached per distinct gate *)
(* and local tokens, and each distinct gate's matrix is read once.             *)
(* ============================================================================ *)

$stabilizerFragmentPhase = <|"T" -> Exp[I Pi / 4], SuperDagger["T"] -> Exp[- I Pi / 4], "S" -> I, SuperDagger["S"] -> - I, "Z" -> -1|>;

stabilizerFragmentOperator[g_ -> {j_Integer}, k_Integer] /; KeyExistsQ[$stabilizerFragmentPhase, g] && 1 <= j <= k :=
    DiagonalMatrix[SparseArray[1 + ($stabilizerFragmentPhase[g] - 1) stabilizerBitVector[k, j]]]
stabilizerFragmentOperator["C"["Z" -> {t_Integer}, controls : {___Integer}, {}], k_Integer] /;
    DuplicateFreeQ[Append[controls, t]] && SubsetQ[Range[k], Append[controls, t]] :=
    DiagonalMatrix[SparseArray[1 - 2 Times @@ (stabilizerBitVector[k, #] & /@ Append[controls, t])]]
stabilizerFragmentOperator["C"["NOT" -> {t_Integer}, {c_Integer}, {}], k_Integer] /; c != t && SubsetQ[Range[k], {c, t}] :=
    With[{r = Range[0, 2 ^ k - 1]},
        SparseArray[Thread[Transpose[{BitXor[r, stabilizerBit[k, t] stabilizerBitVector[k, c]] + 1, r + 1}] -> 1], {2 ^ k, 2 ^ k}]
    ]
stabilizerFragmentOperator[_, _] := Missing["NotInFragment"]

(* The tokens applied in order: the first token's matrix acts first. *)
stabilizerFragmentMatrix[tokens_List, k_Integer] := With[{ops = stabilizerFragmentOperator[#, k] & /@ tokens},
    If[tokens =!= {} && FreeQ[ops, _Missing], Dot @@ Reverse[ops], Missing["NotInFragment"]]
]

(* A fragment token with its wires replaced by their positions among the gate's sorted wires. *)
stabilizerLocalToken[token_, wires_List] := With[{position = AssociationThread[wires -> Range[Length[wires]]]},
    Replace[token, {
        (g_ -> w : {___Integer}) :> (g -> Lookup[position, w]),
        "C"[g_ -> t : {___Integer}, c : {___Integer}, a : {___Integer}] :> "C"[g -> Lookup[position, t], Lookup[position, c], Lookup[position, a]]
    }]
]

stabilizerShortcutPairs[qco_QuantumCircuitOperator] := Catenate[stabilizerShortcutPairs /@ qco["Operators"]]
stabilizerShortcutPairs[element_] := {{element, QuantumShortcut[element]}}

stabilizerLabelsHonestQ[pairs_List] := AllTrue[pairs, stabilizerTokenHonestQ]

stabilizerTokenHonestQ[{op_QuantumOperator, tokens_}] := With[{wires = Sort[op["InputOrder"]], basis = op[[1, 2]]},
    wires === Sort[op["OutputOrder"]] && TrueQ[basis["ComputationalQ"]] && ListQ[tokens] &&
        With[{local = stabilizerLocalToken[#, wires] & /@ tokens},
            With[{key = Hash[{stabilizerGateKey[op], local}]},
                Lookup[$stabilizerHonestMemo, key, stabilizerRemember[$stabilizerHonestMemo, key, stabilizerMatrixIsTokenQ[op, local, Length[wires]]]]
            ]
        ]
]
stabilizerTokenHonestQ[_] := True

stabilizerMatrixIsTokenQ[op_, local_, k_] := With[{expected = stabilizerFragmentMatrix[local, k], matrix = op["Sort"]["Matrix"]},
    MatrixQ[expected] && Dimensions[matrix] === Dimensions[expected] && stabilizerArrayZeroQ[stabilizerSparse[matrix] - expected]
]
