(* Physics tests for the named angular-momentum operators:
   QuantumOperator["JX"/"JY"/"JZ"[j]], their AngularMomentumX/Y/Z aliases, and
   the ladder operators "J+"/"J-".

   These assert matrix content and Lie-algebra structure, not just operator
   dimensions. A dimensions-only check or a serialize-and-reconstruct roundtrip
   cannot see a component that returns the wrong matrix of the right shape (the
   buggy and the correct JX are both 2x2, and a wrong-but-self-consistent
   matrix roundtrips fine). That blind spot is how JX and JY once silently
   returned the diagonal JZ matrix. The su(2) algebra, Casimir, and distinctness
   tests below fail on that defect and pass only for a genuine spin-j
   representation. *)

jMat[name_, j_] := Normal[QuantumOperator[name[j]]["Matrix"]]
comm[a_, b_] := a . b - b . a
zeroMat[j_] := ConstantArray[0, {2 j + 1, 2 j + 1}]

$spins = {1/2, 1, 3/2, 2};


BeginTestSection["AngularMomentum - dimensions"]

Do[
    VerificationTest[QuantumOperator["JX"[j]]["Dimensions"], {2 j + 1, 2 j + 1}, TestID -> "Dim-JX-" <> ToString[j]];
    VerificationTest[QuantumOperator["JZ"[j]]["Dimensions"], {2 j + 1, 2 j + 1}, TestID -> "Dim-JZ-" <> ToString[j]],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - Hermiticity"]

(* Each Cartesian component is a physical observable, hence a Hermitian matrix. *)
Do[
    VerificationTest[HermitianMatrixQ[jMat["JX", j]], True, TestID -> "Herm-JX-" <> ToString[j]];
    VerificationTest[HermitianMatrixQ[jMat["JY", j]], True, TestID -> "Herm-JY-" <> ToString[j]];
    VerificationTest[HermitianMatrixQ[jMat["JZ", j]], True, TestID -> "Herm-JZ-" <> ToString[j]],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - su(2) algebra"]

(* The defining relations [J_a, J_b] = I eps_abc J_c. This is the core check: it
   fails outright if any component returns the wrong matrix (e.g. a diagonal JZ
   handed to JX makes the left side vanish while the right side does not). *)
Do[
    VerificationTest[Simplify[comm[jMat["JX", j], jMat["JY", j]] - I jMat["JZ", j]], zeroMat[j], TestID -> "Comm-XY-" <> ToString[j]];
    VerificationTest[Simplify[comm[jMat["JY", j], jMat["JZ", j]] - I jMat["JX", j]], zeroMat[j], TestID -> "Comm-YZ-" <> ToString[j]];
    VerificationTest[Simplify[comm[jMat["JZ", j], jMat["JX", j]] - I jMat["JY", j]], zeroMat[j], TestID -> "Comm-ZX-" <> ToString[j]],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - Casimir"]

(* J^2 = JX^2 + JY^2 + JZ^2 = j(j+1) 1 on the spin-j irreducible representation. *)
Do[
    VerificationTest[
        Simplify[jMat["JX", j] . jMat["JX", j] + jMat["JY", j] . jMat["JY", j] + jMat["JZ", j] . jMat["JZ", j] - j (j + 1) IdentityMatrix[2 j + 1]],
        zeroMat[j],
        TestID -> "Casimir-" <> ToString[j]
    ],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - JZ form and spectra"]

(* JZ is diagonal in the weight basis, ordered from the highest weight m = +j
   down to -j, so state |0> is the highest weight. This makes JZ[1/2] =
   diag(1/2, -1/2), sign-consistent with QF's PauliZ = diag(1, -1). Every
   component carries the full spin-j spectrum {-j, ..., j}. *)
Do[
    VerificationTest[jMat["JZ", j], DiagonalMatrix[Range[j, -j, -1]], TestID -> "JZ-diag-" <> ToString[j]];
    VerificationTest[Sort[Eigenvalues[jMat["JX", j]]], Range[-j, j], TestID -> "Spec-JX-" <> ToString[j]];
    VerificationTest[Sort[Eigenvalues[jMat["JY", j]]], Range[-j, j], TestID -> "Spec-JY-" <> ToString[j]];
    VerificationTest[Sort[Eigenvalues[jMat["JZ", j]]], Range[-j, j], TestID -> "Spec-JZ-" <> ToString[j]],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - spin-1/2 reduces to Pauli/2"]

VerificationTest[jMat["JX", 1/2], Normal[QuantumOperator["PauliX"]["Matrix"]] / 2, TestID -> "Pauli-JX"]

VerificationTest[jMat["JY", 1/2], Normal[QuantumOperator["PauliY"]["Matrix"]] / 2, TestID -> "Pauli-JY"]

VerificationTest[jMat["JZ", 1/2], Normal[QuantumOperator["PauliZ"]["Matrix"]] / 2, TestID -> "Pauli-JZ"]

EndTestSection[]


BeginTestSection["AngularMomentum - aliases"]

Do[
    VerificationTest[jMat["AngularMomentumX", j], jMat["JX", j], TestID -> "Alias-X-" <> ToString[j]];
    VerificationTest[jMat["AngularMomentumY", j], jMat["JY", j], TestID -> "Alias-Y-" <> ToString[j]];
    VerificationTest[jMat["AngularMomentumZ", j], jMat["JZ", j], TestID -> "Alias-Z-" <> ToString[j]],
    {j, {1/2, 1}}
]

EndTestSection[]


BeginTestSection["AngularMomentum - components are distinct (regression)"]

(* Direct guard against the specific defect where JX, JY, and JZ all returned
   the same diagonal matrix. JX and JY must differ from JZ and must carry the
   off-diagonal (raising/lowering) structure that a diagonal matrix lacks. *)
Do[
    VerificationTest[jMat["JX", j] =!= jMat["JZ", j], True, TestID -> "Distinct-XZ-" <> ToString[j]];
    VerificationTest[jMat["JY", j] =!= jMat["JZ", j], True, TestID -> "Distinct-YZ-" <> ToString[j]];
    VerificationTest[jMat["JX", j] =!= jMat["JY", j], True, TestID -> "Distinct-XY-" <> ToString[j]];
    VerificationTest[DiagonalMatrixQ[jMat["JX", j]], False, TestID -> "Offdiag-JX-" <> ToString[j]];
    VerificationTest[DiagonalMatrixQ[jMat["JY", j]], False, TestID -> "Offdiag-JY-" <> ToString[j]],
    {j, {1/2, 1, 3/2}}
]

EndTestSection[]


BeginTestSection["AngularMomentum - ladder operators"]

(* "J+" and "J-" are a Hermitian-conjugate pair and rebuild the Cartesian
   components: JX = (J+ + J-)/2 and JY = (J- - J+)/(2 I). They are genuine
   ladder operators, eigenoperators of ad_JZ with eigenvalue +/-1.

   Convention note: QF orders the weight basis with |0> = highest weight, and in
   that basis "J+" == JX - I JY, so QF's "J+" lowers the JZ eigenvalue and "J-"
   raises it. This is the reverse of the textbook labels J_pm = JX +/- I JY. The
   last two tests pin QF's actual convention and will flip if the labels are
   ever swapped. *)
Do[
    VerificationTest[ConjugateTranspose[jMat["J+", j]], jMat["J-", j], TestID -> "Ladder-dagger-" <> ToString[j]];
    VerificationTest[jMat["JX", j], (jMat["J+", j] + jMat["J-", j]) / 2, TestID -> "Ladder-JX-" <> ToString[j]];
    VerificationTest[Simplify[jMat["JY", j] - (jMat["J-", j] - jMat["J+", j]) / (2 I)], zeroMat[j], TestID -> "Ladder-JY-" <> ToString[j]];
    VerificationTest[Simplify[comm[jMat["JZ", j], jMat["J+", j]] + jMat["J+", j]], zeroMat[j], TestID -> "Ladder-JplusLowers-" <> ToString[j]];
    VerificationTest[Simplify[comm[jMat["JZ", j], jMat["J-", j]] - jMat["J-", j]], zeroMat[j], TestID -> "Ladder-JminusRaises-" <> ToString[j]],
    {j, {1/2, 1, 3/2}}
]

(* The axis-labeled ladder aliases all resolve to the same raising/lowering
   matrices as "J+"/"J-". *)
VerificationTest[jMat["JX+", 1], jMat["J+", 1], TestID -> "Ladder-alias-JXplus"]

VerificationTest[jMat["JZ+", 1], jMat["J+", 1], TestID -> "Ladder-alias-JZplus"]

VerificationTest[jMat["JX-", 1], jMat["J-", 1], TestID -> "Ladder-alias-JXminus"]

EndTestSection[]
