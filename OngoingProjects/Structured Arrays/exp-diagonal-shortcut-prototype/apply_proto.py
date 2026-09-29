# Apply the prototype of the diagonal exponential to a paclet tree exported from main a1799b61:
#   git archive a1799b61 QuantumFramework Tests | tar -x -C <dir>
#   python3 apply_proto.py <dir> [sparse-state-matrix]
# Every edit must match exactly once.
import sys
root = sys.argv[1]
K = root + '/QuantumFramework/Kernel/'

def rep(s, a, b):
    assert s.count(a) == 1, a[:90]
    return s.replace(a, b)

u = K + 'Utilities.m'; s = open(u).read()
s = rep(s, 'PackageScope["zeroBasePower"]\n', 'PackageScope["zeroBasePower"]\nPackageScope["diagonalMatrixQ"]\nPackageScope["finiteArrayQ"]\nPackageScope["matrixExponential"]\n')
block = '''(* For the exponential, a matrix counts as diagonal when it is square and every
   off-diagonal entry is zero, with no tolerance: a coupling at roundoff relative to
   the largest entry still mixes two degenerate levels over a long time, and
   MatrixExp resolves that mixing. A symbolic off-diagonal entry counts as zero when
   DiagonalMatrixQ's zero test proves it zero; that test is best effort, and an
   identically zero entry it cannot prove sends the matrix to MatrixExp. *)
diagonalMatrixQ[mat_] := SquareMatrixQ[mat] && DiagonalMatrixQ[mat, Tolerance -> 0]

(* No entry is infinite or indeterminate. A SparseArray is atomic, so its stored
   values and its background are read. *)
finiteArrayQ[a_SparseArray] := finiteArrayQ[{a["NonzeroValues"], a["Background"]}]
finiteArrayQ[a_] := FreeQ[a, Indeterminate | _DirectedInfinity]

(* The exponential of a diagonal matrix is the diagonal of the exponentials of its
   entries, and its action on a vector multiplies the vector by them. Any other
   matrix, and a diagonal with an infinite or indeterminate entry, goes to
   MatrixExp, which fails on the latter. *)
matrixExponential[mat_ ? diagonalMatrixQ, v___] := With[{d = Normal[Diagonal[mat]]},
    diagonalAction[diagonalExp[d], v] /; finiteArrayQ[d]
]
matrixExponential[mat_, v___] := MatrixExp[mat, v]

diagonalAction[values_] := DiagonalMatrix[values, TargetStructure -> "Sparse"]
diagonalAction[values_, v_] := values v

(* A numeric diagonal at machine precision, exact entries included, is exponentiated
   in machine numbers, as MatrixExp does, and as one packed array. On a packed array
   Exp returns an exponential below the smallest normalized machine number as a
   subnormal number or zero, and one above the largest as an arbitrary-precision
   number, without a message; on an unpacked list it raises General::munfl. A numeric
   diagonal at a higher finite precision is exponentiated at that precision. *)
diagonalExp[d_ ? machineVectorQ] := Exp[packedVector[N[d]]]
diagonalExp[d_ ? (VectorQ[#, NumericQ] && Precision[#] < Infinity &)] := Exp[N[d, Precision[d]]]
diagonalExp[d_] := Exp[d]

machineVectorQ[d_] := VectorQ[d, NumericQ] && Precision[d] === MachinePrecision

packedVector[y_ ? (FreeQ[#, _Complex] &)] := Developer`ToPackedArray[y, Real]
packedVector[y_] := Developer`ToPackedArray[y, Complex]


'''
anchor = 'SetPrecisionNumeric[x_ /; NumericQ[x] || ArrayQ[x, _, NumericQ]]'
s = rep(s, anchor, block + anchor)
open(u, 'w').write(s)

q = K + 'QuantumOperator/QuantumOperator.m'; s = open(q).read()
s = rep(s, 'scalarBasePower[base_, mat_] := MatrixExp[Log[base] mat]', 'scalarBasePower[base_, mat_] := matrixExponential[Log[base] mat]')
s = rep(s, 'QuantumOperator /: MatrixExp[qo_QuantumOperator] := matrixMapOperator[MatrixExp, qo, Exp]', 'QuantumOperator /: MatrixExp[qo_QuantumOperator] := matrixMapOperator[matrixExponential, qo, Exp]')
s = rep(s, '    ConfirmBy[Confirm[g[ConfirmBy[mat, FreeQ[#, Indeterminate | _DirectedInfinity] &]]], MatrixQ],', '    ConfirmBy[Confirm[g[ConfirmBy[mat, finiteArrayQ]]], MatrixQ],')
s = rep(s, '''   form; the sorting happens inside the body instead. Otherwise it is the sorted
   matrix as a List: a SparseArray is atomic, so neither Function application nor
   ReplaceAll would reach the parameters inside it. The basis written into the body
   carries no parameter specification, which a substitution would otherwise
   overwrite. *)''', '''   form; the sorting happens inside the body instead. Otherwise it is the sorted
   matrix: for a diagonal matrix the List of its diagonal entries, made into a sparse
   diagonal matrix only after the substitution, and for any other the matrix as a
   List. A SparseArray is atomic, so neither Function application nor ReplaceAll
   would reach the parameters inside one. The basis written into the body carries no
   parameter specification, which a substitution would otherwise overwrite. *)''')
s = rep(s, 'heldOperatorMatrix[_, op_] := With[{mat = Normal[op["Matrix"]]}, Hold[mat]]', '''heldOperatorMatrix[_, op_] := heldMatrix[op["Matrix"]]

heldMatrix[mat_ ? diagonalMatrixQ] := With[{d = Normal[Diagonal[mat]]}, Hold[DiagonalMatrix[d, TargetStructure -> "Sparse"]]]
heldMatrix[mat_] := With[{m = Normal[mat]}, Hold[m]]''')
s = rep(s, '''QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] := Enclose @ With[{op = qo["Sort"]},
    QuantumState[
        If[ op["VectorQ"] && qs["VectorQ"],
            MatrixExp[op["Matrix"], QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]["StateVector"]],
            ArrayReshape[
                MatrixExp[op["ToMatrix"]["Matrix"], QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]["DensityVector"]],
                {#, #} & @ op["OutputDimension"]
            ]
        ],
        QuantumBasis[
            op["Output"],
            "Label" -> If[op["Label"] === None || qs["Label"] === None, None, Exp[op["Label"]][qs["Label"]]]
        ]
    ]
]
''', '''(* The exponential e^M of the operator acting on the state, as Exp[qo][qs] gives: a
   pure state is multiplied by e^M, a mixed state rho becomes e^M rho e^(M^dagger),
   and a superoperator acts on the density vector. An operator on any other qudits
   than exactly the state's, with their dimensions, goes through Exp[qo][qs], which
   extends the operator by the identity on qudits it does not act on, extends the
   state to qudits it lacks, and fails on a dimension mismatch. *)
QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] /; ! wholeRegisterQ[qo, qs] := Exp[qo][qs]

QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] := Enclose @ With[{op = qo["Sort"]},
    QuantumState[
        ConfirmBy[exponentialAction[op, QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]], ArrayQ],
        QuantumBasis[
            op["Output"],
            "Label" -> If[op["Label"] === None || qs["Label"] === None, None, Exp[op["Label"]][qs["Label"]]]
        ]
    ]
]

wholeRegisterQ[qo_, qs_] := Sort[Transpose[{qo["InputOrder"], qo["InputDimensions"]}]] === Transpose[{Range[qs["OutputQudits"]], qs["OutputDimensions"]}]

exponentialAction[op_ ? (#["VectorQ"] &), qs_ ? (#["VectorQ"] &)] := matrixExponential[op["Matrix"], qs["StateVector"]]
exponentialAction[op_ ? (#["VectorQ"] &), qs_] := With[{u = matrixExponential[op["Matrix"]]}, u . qs["DensityMatrix"] . ConjugateTranspose[u]]
exponentialAction[op_, qs_] := ArrayReshape[matrixExponential[op["ToMatrix"]["Matrix"], qs["DensityVector"]], {#, #} & @ op["OutputDimension"]]
''')
open(q, 'w').write(s)

if len(sys.argv) > 2 and sys.argv[2] == 'sparse-state-matrix':
    p = K + 'QuantumState/Properties.m'; s = open(p).read()
    s = rep(s, '''        Transpose[ReshapeArray[{qs["StateTensor"]}, Join[#, #] & @ qs["MatrixNameDimensions"]], 2 <-> 3],
        qs["MatrixNameDimensions"] ^ 2
    ]
]
''', '''        Transpose[ReshapeArray[stateTensorArray[qs["StateTensor"]], Join[#, #] & @ qs["MatrixNameDimensions"]], 2 <-> 3],
        qs["MatrixNameDimensions"] ^ 2
    ]
]

(* The state tensor as an array for ReshapeArray: the tensor of a state with no
   qudits, such as a full trace, is a scalar. *)
stateTensorArray[t_ /; ArrayDepth[t] == 0] := {t}
stateTensorArray[t_] := t
''')
    open(p, 'w').write(s)
print('applied to', root)
