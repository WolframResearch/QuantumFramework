(* Physics tests for the named angular-momentum operators:
   QuantumOperator["JX"/"JY"/"JZ"[j]], their AngularMomentumX/Y/Z aliases, and
   the ladder operators "J+"/"J-".

   Two layers are asserted, and the distinction is load-bearing:

   1. The FRAME-FAITHFUL layer, "MatrixRepresentation": the operator expressed
      in the computational basis. This is what states, expectation values, and
      compositions actually use, so every physics assertion (su(2), Casimir,
      spectra, ladder relations, expectation values) runs here. A content check
      on the raw stored matrix cannot distinguish a correct operator from a
      wrong-frame operator whose stored numbers merely look right; the
      frame-faithful layer can.

   2. The STORAGE-CONTRACT layer, raw "Matrix": each Cartesian component is
      stored in its own eigenbasis (QuditBasis["JX"[j]] etc., ordered by
      ascending m), where its matrix is DiagonalMatrix[Range[-j, j]]. The raw
      matrices of JX, JY, JZ therefore coincide by design; only the basis tags
      differ. The contract tier pins that storage so the display behavior is a
      documented invariant, not a surprise. *)

compMat[name_, j_] := Simplify[Normal[QuantumOperator[name[j]]["MatrixRepresentation"]]]
rawMat[name_, j_] := Normal[QuantumOperator[name[j]]["Matrix"]]
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


BeginTestSection["AngularMomentum - own-eigenbasis storage contract"]

(* Raw storage: diagonal ascending in the operator's own eigenbasis. The
   ascending order matches the named bases, whose k-th element carries
   m = -j + k - 1 (QuditBasis["JZ"[j]] is the computational basis reversed). *)
Do[
    VerificationTest[rawMat["JX", j], Normal[DiagonalMatrix[Range[-j, j]]], TestID -> "Raw-diag-JX-" <> ToString[j]];
    VerificationTest[rawMat["JY", j], Normal[DiagonalMatrix[Range[-j, j]]], TestID -> "Raw-diag-JY-" <> ToString[j]];
    VerificationTest[rawMat["JZ", j], Normal[DiagonalMatrix[Range[-j, j]]], TestID -> "Raw-diag-JZ-" <> ToString[j]],
    {j, $spins}
]

VerificationTest[
    Normal /@ Values[QuditBasis["JZ"[1]]["Association"]],
    Reverse[IdentityMatrix[3]],
    TestID -> "Basis-JZ-reversed"
]

EndTestSection[]


BeginTestSection["AngularMomentum - Hermiticity"]

(* Each Cartesian component is a physical observable, hence a Hermitian matrix
   in the computational frame. *)
Do[
    VerificationTest[HermitianMatrixQ[compMat["JX", j]], True, TestID -> "Herm-JX-" <> ToString[j]];
    VerificationTest[HermitianMatrixQ[compMat["JY", j]], True, TestID -> "Herm-JY-" <> ToString[j]];
    VerificationTest[HermitianMatrixQ[compMat["JZ", j]], True, TestID -> "Herm-JZ-" <> ToString[j]],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - su(2) algebra"]

(* The defining relations [J_a, J_b] = I eps_abc J_c, asserted on the
   frame-faithful layer. A wrong-frame component (right numbers, wrong basis
   tag) breaks these while leaving its raw stored matrix looking plausible. *)
Do[
    VerificationTest[Simplify[comm[compMat["JX", j], compMat["JY", j]] - I compMat["JZ", j]], zeroMat[j], TestID -> "Comm-XY-" <> ToString[j]];
    VerificationTest[Simplify[comm[compMat["JY", j], compMat["JZ", j]] - I compMat["JX", j]], zeroMat[j], TestID -> "Comm-YZ-" <> ToString[j]];
    VerificationTest[Simplify[comm[compMat["JZ", j], compMat["JX", j]] - I compMat["JY", j]], zeroMat[j], TestID -> "Comm-ZX-" <> ToString[j]],
    {j, $spins}
]

(* Same-basis products through QF's own operator composition, not hand Dot. *)
VerificationTest[
    Simplify[
        Normal[(QuantumOperator["JX"[1]] @ QuantumOperator["JX"[1]])["MatrixRepresentation"]] -
        compMat["JX", 1] . compMat["JX", 1]
    ],
    zeroMat[1],
    TestID -> "Square-JX-object-1"
]

(* CHARACTERIZATION (present behavior, not an endorsement): cross-basis
   composition qo1[qo2] with two different tagged bases does not reconcile
   frames, on any operator family, so the object product of two different
   components is NOT the matrix product of their computational representations.
   Algebra between components is therefore asserted on the matrix layer above.
   This test flips and announces the change if composition ever learns frames. *)
VerificationTest[
    Simplify[Normal[(QuantumOperator["JX"[1]] @ QuantumOperator["JY"[1]])["MatrixRepresentation"]]] =!=
        Simplify[compMat["JX", 1] . compMat["JY", 1]],
    True,
    TestID -> "CrossBasis-product-characterization"
]

EndTestSection[]


BeginTestSection["AngularMomentum - Casimir"]

(* J^2 = JX^2 + JY^2 + JZ^2 = j(j+1) 1 on the spin-j irreducible representation. *)
Do[
    VerificationTest[
        Simplify[compMat["JX", j] . compMat["JX", j] + compMat["JY", j] . compMat["JY", j] + compMat["JZ", j] . compMat["JZ", j] - j (j + 1) IdentityMatrix[2 j + 1]],
        zeroMat[j],
        TestID -> "Casimir-" <> ToString[j]
    ],
    {j, $spins}
]

(* Symbolic-j statements live on the weight spectrum: an irrep's matrix has
   fixed dimension 2j+1, so no operator matrix exists at symbolic j, but its
   moments do. Second and fourth moments, Tr[JZ^2] and Tr[JZ^4], as closed
   forms in j, each tied to a concrete trace at j = 2. *)
VerificationTest[
    Simplify[Sum[m^2, {m, -\[FormalJ], \[FormalJ]}] - \[FormalJ] (\[FormalJ] + 1) (2 \[FormalJ] + 1) / 3],
    0,
    TestID -> "Weight-sum-symbolic-j"
]

VerificationTest[Tr[compMat["JZ", 2] . compMat["JZ", 2]], 10, TestID -> "Weight-sum-trace-2"]

VerificationTest[
    Simplify[Sum[m^4, {m, -\[FormalJ], \[FormalJ]}] - \[FormalJ] (\[FormalJ] + 1) (2 \[FormalJ] + 1) (3 \[FormalJ]^2 + 3 \[FormalJ] - 1) / 15],
    0,
    TestID -> "Weight-sum4-symbolic-j"
]

VerificationTest[Tr[MatrixPower[compMat["JZ", 2], 4]], 34, TestID -> "Weight-sum4-trace-2"]

EndTestSection[]


BeginTestSection["AngularMomentum - rotation generators"]

(* The components are the generators of rotations: the named WignerD operator at
   a single (symbolic) Euler angle equals exp(-I b JY) on the same
   representation. This ties the static algebra to the dynamics it generates. *)
Do[
    VerificationTest[
        FullSimplify[
            Normal[QuantumOperator["WignerD"[j, \[FormalB]]]["MatrixRepresentation"]] -
            MatrixExp[-I \[FormalB] compMat["JY", j]]
        ],
        zeroMat[j],
        TestID -> "Generator-WignerD-Ry-" <> ToString[j]
    ],
    {j, {1, 3/2}}
]

EndTestSection[]


BeginTestSection["AngularMomentum - JZ form and spectra"]

(* Computationally, JZ is diagonal in the weight basis ordered from the highest
   weight m = +j down to -j, so state |0> is the highest weight. This makes
   JZ[1/2] = diag(1/2, -1/2), sign-consistent with QF's PauliZ = diag(1, -1).
   Every component carries the full spin-j spectrum {-j, ..., j}. *)
Do[
    VerificationTest[compMat["JZ", j], Normal[DiagonalMatrix[Range[j, -j, -1]]], TestID -> "JZ-diag-" <> ToString[j]];
    VerificationTest[Sort[Eigenvalues[compMat["JX", j]]], Range[-j, j], TestID -> "Spec-JX-" <> ToString[j]];
    VerificationTest[Sort[Eigenvalues[compMat["JY", j]]], Range[-j, j], TestID -> "Spec-JY-" <> ToString[j]];
    VerificationTest[Sort[Eigenvalues[compMat["JZ", j]]], Range[-j, j], TestID -> "Spec-JZ-" <> ToString[j]],
    {j, $spins}
]

EndTestSection[]


BeginTestSection["AngularMomentum - spin-1/2 reduces to Pauli/2"]

VerificationTest[compMat["JX", 1/2], Normal[QuantumOperator["PauliX"]["MatrixRepresentation"]] / 2, TestID -> "Pauli-JX"]

VerificationTest[compMat["JY", 1/2], Normal[QuantumOperator["PauliY"]["MatrixRepresentation"]] / 2, TestID -> "Pauli-JY"]

VerificationTest[compMat["JZ", 1/2], Normal[QuantumOperator["PauliZ"]["MatrixRepresentation"]] / 2, TestID -> "Pauli-JZ"]

EndTestSection[]


BeginTestSection["AngularMomentum - aliases"]

Do[
    VerificationTest[compMat["AngularMomentumX", j], compMat["JX", j], TestID -> "Alias-X-" <> ToString[j]];
    VerificationTest[compMat["AngularMomentumY", j], compMat["JY", j], TestID -> "Alias-Y-" <> ToString[j]];
    VerificationTest[compMat["AngularMomentumZ", j], compMat["JZ", j], TestID -> "Alias-Z-" <> ToString[j]],
    {j, {1/2, 1}}
]

EndTestSection[]


BeginTestSection["AngularMomentum - components are distinct (regression)"]

(* The three components must be distinct PHYSICAL operators. Their raw stored
   matrices coincide by the storage contract, so distinctness is a
   frame-faithful statement: computationally JX and JY carry off-diagonal
   (raising/lowering) structure that the diagonal JZ lacks. *)
Do[
    VerificationTest[compMat["JX", j] =!= compMat["JZ", j], True, TestID -> "Distinct-XZ-" <> ToString[j]];
    VerificationTest[compMat["JY", j] =!= compMat["JZ", j], True, TestID -> "Distinct-YZ-" <> ToString[j]];
    VerificationTest[compMat["JX", j] =!= compMat["JY", j], True, TestID -> "Distinct-XY-" <> ToString[j]];
    VerificationTest[DiagonalMatrixQ[compMat["JX", j]], False, TestID -> "Offdiag-JX-" <> ToString[j]];
    VerificationTest[DiagonalMatrixQ[compMat["JY", j]], False, TestID -> "Offdiag-JY-" <> ToString[j]],
    {j, {1/2, 1, 3/2}}
]

EndTestSection[]


BeginTestSection["AngularMomentum - eigenstate expectations"]

(* The basis-aware round trip that a content check cannot see: the k-th element
   of the named eigenbasis carries m = -j + k - 1, so the expectation of the
   matching operator on that eigenstate must return exactly that eigenvalue.
   Includes the degenerate m = 0 eigenstate. The trace form runs for all three
   components; the bra-op-ket form is additionally pinned for JX and JZ, whose
   eigenbases are real (the Dagger/Scalar contraction mishandles a complex dual
   basis, a defect in the bra route itself, reproduced on non-J operators too,
   so JY is asserted through the trace form only). *)
expectJ[name_, j_, k_] := With[{op = QuantumOperator[name[j]], st = QuantumState[UnitVector[2 j + 1, k], name[j]]},
    Simplify[Tr[Normal[op["MatrixRepresentation"]] . Normal[QuantumState[st, QuantumBasis[2 j + 1]]["DensityMatrix"]]]]
]
braketJ[name_, j_, k_] := With[{op = QuantumOperator[name[j]], st = QuantumState[UnitVector[2 j + 1, k], name[j]]},
    Simplify[st["Dagger"][op[st]]["Scalar"]]
]

Do[
    VerificationTest[expectJ[name, j, k], -j + k - 1, TestID -> "Expect-" <> name <> "-" <> ToString[j] <> "-m" <> ToString[-j + k - 1]],
    {name, {"JX", "JY", "JZ"}}, {j, {1/2, 1, 3/2}}, {k, 2 j + 1}
]

Do[
    VerificationTest[braketJ[name, j, k], -j + k - 1, TestID -> "Braket-" <> name <> "-" <> ToString[j] <> "-m" <> ToString[-j + k - 1]],
    {name, {"JX", "JZ"}}, {j, {1/2, 1, 3/2}}, {k, 2 j + 1}
]

(* CHARACTERIZATION (present behavior, not an endorsement): the Dagger/Scalar
   contraction drops the bra conjugation when the state's basis has complex
   elements, so the bra route returns the wrong value on the JY eigenbasis
   (the trace-form tier above carries the actual physics). This test flips and
   announces the change if the bra route learns complex dual bases. *)
VerificationTest[braketJ["JY", 1, 3] =!= 1, True, TestID -> "Braket-JY-complex-dual-characterization"]

EndTestSection[]


BeginTestSection["AngularMomentum - numerical reference (independent construction)"]

(* An independent numeric route: the spin-3/2 matrices rebuilt from the ladder
   construction in machine arithmetic, no QF calls, compared entrywise against
   the operators' numericized computational representations. *)
ladderN[j_] := With[{ms = Range[j, -j, -1]},
    N @ Table[Sqrt[j (j + 1) - m2 (m2 + 1)] KroneckerDelta[m1, m2 + 1], {m1, ms}, {m2, ms}]
]
numRef = With[{jp = ladderN[3/2]},
    <|"JX" -> (jp + ConjugateTranspose[jp]) / 2,
      "JY" -> (jp - ConjugateTranspose[jp]) / (2 I),
      "JZ" -> DiagonalMatrix[N @ Range[3/2, -3/2, -1]]|>
];
Do[
    VerificationTest[Max[Abs[N[compMat[name, 3/2]] - numRef[name]]] < 10^-12, True, TestID -> "NumRef-" <> name <> "-3/2"],
    {name, {"JX", "JY", "JZ"}}
]

EndTestSection[]


BeginTestSection["AngularMomentum - addition of angular momenta"]

(* Two coupled spins, Jtot_k = J_k x 1 + 1 x J_k, built from the frame-faithful
   component matrices. Jtot^2 block-decomposes the product space into irreps
   with eigenvalues jtot(jtot+1): 1/2 x 1/2 -> triplet + singlet {2,2,2,0};
   1/2 x 1 -> jtot = 3/2, 1/2 with {15/4 (x4), 3/4 (x2)}. A wrong component
   shifts these multiplets, and no single irrep can exhibit this structure. *)
totalJ2[j1_, j2_] := With[{d1 = 2 j1 + 1, d2 = 2 j2 + 1},
    Total[Map[
        With[{a = compMat[#, j1], b = compMat[#, j2]},
            MatrixPower[KroneckerProduct[a, IdentityMatrix[d2]] + KroneckerProduct[IdentityMatrix[d1], b], 2]
        ] &,
        {"JX", "JY", "JZ"}
    ]]
]

VerificationTest[Sort[Eigenvalues[totalJ2[1/2, 1/2]]], {0, 2, 2, 2}, TestID -> "Coupling-half-half-spectrum"]

VerificationTest[Sort[Eigenvalues[totalJ2[1/2, 1]]], {3/4, 3/4, 15/4, 15/4, 15/4, 15/4}, TestID -> "Coupling-half-one-spectrum"]

EndTestSection[]


BeginTestSection["AngularMomentum - zero-field-splitting spectrum"]

(* The spin-1 ZFS Hamiltonian H = 2 d JZ^2 + e (JX^2 - JY^2), couplings kept
   symbolic, must carry the textbook spectrum {0, 2d - e, 2d + e}. This is the
   canonical composite consumer of all three components at once: any wrong
   component shifts the spectrum. *)
VerificationTest[
    Simplify[
        Sort[Eigenvalues[2 \[FormalD] compMat["JZ", 1] . compMat["JZ", 1] + \[FormalE] (compMat["JX", 1] . compMat["JX", 1] - compMat["JY", 1] . compMat["JY", 1])]] -
        Sort[{0, 2 \[FormalD] - \[FormalE], 2 \[FormalD] + \[FormalE]}]
    ],
    {0, 0, 0},
    TestID -> "ZFS-spectrum-symbolic"
]

EndTestSection[]


BeginTestSection["AngularMomentum - classical large-j limit"]

(* Normalized components S_k = J_k/j commute in the large-j (classical) limit at
   the exact rate ||[S_x, S_y]|| = ||J_z||/j^2 = 1/j. *)
Do[
    VerificationTest[
        Simplify[Norm[comm[compMat["JX", j], compMat["JY", j]] / j^2] - 1/j],
        0,
        TestID -> "Classical-rate-" <> ToString[j]
    ],
    {j, {1/2, 3/2, 7/2}}
]

(* The highest-weight state |j, j> (computational e1) is a spin coherent state:
   minimum-uncertainty at every j, Delta JX Delta JY = |<JZ>|/2 = j/2 exactly,
   and its RELATIVE spread Delta JX / j = 1/Sqrt[2 j] vanishes as j -> Infinity:
   the multiplet's angular momentum becomes a classical vector. *)
varTop[name_, j_] := With[{m = compMat[name, j], v = UnitVector[2 j + 1, 1]},
    Simplify[v . m . m . v - (v . m . v)^2]
]

Do[
    VerificationTest[
        {varTop["JX", j], varTop["JY", j], Simplify[UnitVector[2 j + 1, 1] . compMat["JZ", j] . UnitVector[2 j + 1, 1]]},
        {j/2, j/2, j},
        TestID -> "Coherent-minimum-uncertainty-" <> ToString[j]
    ],
    {j, {1/2, 3/2, 7/2}}
]

VerificationTest[Limit[Sqrt[\[FormalJ]/2]/\[FormalJ], \[FormalJ] -> Infinity], 0, TestID -> "Coherent-relative-spread-limit"]

EndTestSection[]


BeginTestSection["AngularMomentum - ladder operators"]

(* "J+" and "J-" are a Hermitian-conjugate pair and rebuild the Cartesian
   components: JX = (J+ + J-)/2 and JY = (J+ - J-)/(2 I). On the frame-faithful
   layer they follow the textbook convention: J+ raises the JZ eigenvalue
   ([JZ, J+] = +J+) and J- lowers it. Only their RAW stored matrices look
   swapped, because the storage frame orders m ascending while the
   computational frame orders m descending. The JplusRaises/JminusLowers tests
   pin the physical convention and will flip if the labels are ever swapped. *)
Do[
    VerificationTest[ConjugateTranspose[compMat["J+", j]], compMat["J-", j], TestID -> "Ladder-dagger-" <> ToString[j]];
    VerificationTest[compMat["JX", j], (compMat["J+", j] + compMat["J-", j]) / 2, TestID -> "Ladder-JX-" <> ToString[j]];
    VerificationTest[Simplify[compMat["JY", j] - (compMat["J+", j] - compMat["J-", j]) / (2 I)], zeroMat[j], TestID -> "Ladder-JY-" <> ToString[j]];
    VerificationTest[Simplify[comm[compMat["JZ", j], compMat["J+", j]] - compMat["J+", j]], zeroMat[j], TestID -> "Ladder-JplusRaises-" <> ToString[j]];
    VerificationTest[Simplify[comm[compMat["JZ", j], compMat["J-", j]] + compMat["J-", j]], zeroMat[j], TestID -> "Ladder-JminusLowers-" <> ToString[j]],
    {j, {1/2, 1, 3/2}}
]

(* The axis-labeled ladders share the "J+"/"J-" raw content but live in their
   axis's eigenbasis, so only the z-axis pair coincides with "J+"/"J-" as a
   physical operator; the x-axis pair are the ad_JX eigenoperators instead. *)
VerificationTest[compMat["JZ+", 1], compMat["J+", 1], TestID -> "Ladder-alias-JZplus"]

VerificationTest[rawMat["JX+", 1], rawMat["J+", 1], TestID -> "Ladder-raw-JXplus"]

VerificationTest[Simplify[comm[compMat["JX", 1], compMat["JX+", 1]] - compMat["JX+", 1]], zeroMat[1], TestID -> "Ladder-JXplus-raises-JX"]

VerificationTest[Simplify[comm[compMat["JX", 1], compMat["JX-", 1]] + compMat["JX-", 1]], zeroMat[1], TestID -> "Ladder-JXminus-lowers-JX"]

EndTestSection[]
