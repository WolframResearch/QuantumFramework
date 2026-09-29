(* The exponential of a diagonal operator is the diagonal of the exponentials of its
   entries, e^(c diag(h)) = diag(e^(c h_k)). Every reference is built without the
   route under test: a closed form, MatrixExp of an explicit matrix, the 40-digit
   exponential of the stored entries, or another QF route (Exp[qo][qs], a circuit of
   rotation gates, "MatrixRepresentation").

   Physics brief, written after the prototype and before these tests:
     symbolic object  U(theta) = e^(-i theta H) for H = sum_(i<j) J_ij Z_i Z_j +
                      sum_k b_k Z_k on n qubits with symbolic J, b and theta, acting on
                      operators and on states; the Gibbs operator e^(beta J sum_k
                      Z_k Z_(k+1)) of the open Ising chain with symbolic beta J; the
                      dephasing map e^(t L), L rho = -i [H, rho] + gamma sum_k
                      (Z_k rho Z_k - rho), with symbolic gamma, t and rho; e^M rho
                      e^(M^dagger) with symbolic theta, g and rho; the depth-one QAOA
                      state on the 4-ring with symbolic gamma and beta
     invariants       unitarity, the norm of a state and the group law U(theta1)
                      U(theta2) = U(theta1 + theta2); [H, Z_k] = 0; trace, populations
                      and complete positivity under dephasing; the trace a + (1 - a)
                      e^(-g) of the no-jump evolution; Gibbs weights unchanged by a
                      constant energy shift
     refused          theta = 0, a single qubit, and random couplings compared with
                      Exp of the same list

   Regimes covered:
     general symbolic     the Ising closed form for n = 2..5 and, for n = 2..4, its
                          factorization into the QAOA cost-layer gates R_ZZ(2 theta
                          J_ij) and R_Z(2 theta b_k); unitarity and the group law at
                          n = 3; e^(-i theta H) on a symbolic 3-qubit state; 2-qubit
                          dephasing with symbolic gamma, t and rho; e^M rho e^(M^dagger)
                          for M = -i theta Z, X, Y, ZZ and the no-jump generator, on a
                          whole register, on part of one and beyond it; the QAOA edge
                          term (2 + sin 4 beta sin 2 gamma)/4
     exactly solvable     dephasing coherences e^(-i (E_x - E_y) t - 2 gamma d(x, y) t)
                          and their completely positive multiplier; collective
                          dephasing, where |01><10| is decoherence-free; one half of a
                          Bell pair dephased; the partition function 2 (2 cosh(beta
                          J))^(n - 1) and the correlations tanh(beta J)^(r - 1) of the
                          open chain; the concurrence |sin 2g| of e^(-i g ZZ)|++>; the
                          spin-1 e^(-i theta J_x), stored diagonal in its own eigenbasis
     limiting             t -> Infinity gives the pinching map rho -> diag(rho);
                          gamma -> 0 gives U rho U^dagger; beta J -> Infinity leaves
                          the chain in its two ground states; theta up to 10^6, where
                          theta ||H|| reaches 10^7, keeps U unitary to the last place
     numerical reference  QuantumEvolve's numerical integration of the Lindblad
                          equation against Exp[t L]; the QAOA edge term from the cost
                          layer with gamma declared, at 25 values of gamma; and, as a
                          rounding test of the stored entries, their 40-digit
                          exponential, within one rounding of theta h_k plus one unit
     failure / edge       e^-2000 underflow, e^800 and e^750 overflow, subnormal
                          Boltzmann factors, a complex entry with a subnormal real part,
                          a 30-digit diagonal, infinite and indeterminate entries,
                          zeros that simplify away (proved and not proved), a coupling
                          that vanishes at some parameter values, a resonant coupling
                          at roundoff, colliding parameters, a matrix-type object with
                          no qudits, an operator on a qudit the state lacks, and a
                          qubit operator on a qutrit

   The symbolic, exactly solvable and limiting rows pass on main too, where MatrixExp
   reaches the same closed forms: they guard against regressions, and those meant to
   cover the entry-by-entry route assert that the generator is diagonal with no
   tolerance. The rows main fails are the rounding test, the scaling tests, most edge
   cases, and the exponential acting on a mixed state, on part of a register or
   beyond it. *)

BeginTestSection["QuantumOperator - exponential of a diagonal operator"]

deSpins[n_] := 1 - 2 Tuples[{0, 1}, n]

(* The Ising energy E_x = sum_(i<j) J_ij z_i z_j + sum_k b_k z_k of every basis state. *)
deIsingEnergies[couplings_, fields_] := With[{s = deSpins[Length[fields]]},
    Total[(s . UpperTriangularize[couplings, 1]) s, {2}] + s . fields
]

(* sum_(i<j) J_ij Z_i Z_j, built from QF's two-qubit ZZ operators, and with the fields
   sum_k b_k Z_k added. *)
deZZOperator[couplings_] := Total[Apply[couplings[[#1, #2]] QuantumOperator["ZZ", {#1, #2}] &, Subsets[Range[Length[couplings]], {2}], {1}]]
deIsingOperator[couplings_, fields_] := deZZOperator[couplings] + Total[MapIndexed[#1 QuantumOperator["Z", #2] &, fields]]

(* sum_k Z_k Z_(k+1), the open Ising chain on n qubits. *)
deChain[n_] := Total[QuantumOperator["ZZ", {#, # + 1}] & /@ Range[n - 1]]

(* The error of each entry against the 40-digit exponential of the exact product of
   theta and the stored entry, in units of one rounding of that product plus one unit
   of the result. The difference is taken at 40 digits, so the reference is not
   rounded to a machine number first. *)
deReferenceError[values_, theta_, stored_] := With[{phases = SetPrecision[theta, 40] SetPrecision[stored, 40]},
    Max[Abs[SetPrecision[values, 40] - Exp[-I phases]] / ((Abs[phases] / 2 + 1) $MachineEpsilon)]
]

deHamiltonian2 = deJ12 QuantumOperator["ZZ"] + deB1 QuantumOperator["ZI"] + deB2 QuantumOperator["IZ"];
deEnergies2 = deIsingEnergies[{{0, deJ12}, {0, 0}}, {deB1, deB2}];
deDephasers = {QuantumOperator["ZI"], QuantumOperator["IZ"]};
deRho0 = Array[deR, {4, 4}];
deRho0State = QuantumState[deRho0, QuantumBasis[{2, 2}]];
deReals = Element[{deJ12, deB1, deB2, deGam, deT}, Reals];

(* The two-qubit Liouvillian of that Hamiltonian, dephased by Z_1 and Z_2 with the
   given rates, and rho(t) = e^(t L) rho(0). *)
deLiouvillian[rates_] := QuantumOperator["Liouvillian"[deHamiltonian2, deDephasers, rates]]
deDephase[rates_] := Normal[Exp[deT deLiouvillian[rates]][deRho0State]["DensityMatrix"]]

deHamiltonian8 = deZZOperator[BlockRandom[RandomReal[{-1, 1}, {8, 8}], RandomSeeding -> 7]];
deStored8 = Normal[Diagonal[deHamiltonian8["Sort"]["Matrix"]]];

deEnergies12 = BlockRandom[deIsingEnergies[RandomReal[{-1, 1}, {12, 12}], ConstantArray[0, 12]], RandomSeeding -> 5];
deEnergies16 = BlockRandom[deIsingEnergies[RandomReal[{-1, 1}, {16, 16}], ConstantArray[0, 16]], RandomSeeding -> 7];


(* For n = 2..5 qubits with symbolic couplings, fields and time, the generator
   -i theta H is diagonal with no tolerance, and e^(-i theta H) is the diagonal of
   e^(-i theta E_x). *)
VerificationTest[
    Table[
        With[{couplings = Array[deJ, {n, n}], fields = Array[deB, n]},
            With[{hamiltonian = deIsingOperator[couplings, fields]},
                {
                    DiagonalMatrixQ[(-I deTheta hamiltonian)["Sort"]["Matrix"], Tolerance -> 0],
                    Simplify[Normal[Exp[-I deTheta hamiltonian]["Matrix"]] - DiagonalMatrix[Exp[-I deTheta deIsingEnergies[couplings, fields]]]]
                }
            ]
        ],
        {n, 2, 5}
    ],
    Table[{True, ConstantArray[0, {2 ^ n, 2 ^ n}]}, {n, 2, 5}],
    TestID -> "DiagonalExp-Ising-symbolic-closed-form"
]

(* For real couplings, fields and times, e^(-i theta H) is unitary and obeys the
   group law U(theta1) U(theta2) = U(theta1 + theta2). *)
VerificationTest[
    With[{couplings = Array[deJ, {3, 3}], fields = Array[deB, 3]},
        With[{
            hamiltonian = deIsingOperator[couplings, fields],
            reals = Element[Join[Flatten[UpperTriangularize[couplings, 1]], fields, {deTheta, deTheta1, deTheta2}], Reals]
        },
            With[{u = Normal[Exp[-I deTheta hamiltonian]["Matrix"]]},
                Simplify[{
                    u . ConjugateTranspose[u] - IdentityMatrix[8],
                    Normal[Exp[-I deTheta1 hamiltonian]["Matrix"] . Exp[-I deTheta2 hamiltonian]["Matrix"]] -
                        Normal[Exp[-I (deTheta1 + deTheta2) hamiltonian]["Matrix"]]
                }, reals]
            ]
        ]
    ],
    ConstantArray[0, {2, 8, 8}],
    TestID -> "DiagonalExp-Ising-symbolic-unitary-group-law"
]

(* e^(-i theta H) acting on a pure 3-qubit state, all symbolic: MatrixExp[qo, qs]
   multiplies each amplitude by its phase e^(-i theta E_x), and the norm stays. *)
VerificationTest[
    With[{couplings = Array[deJ, {3, 3}], fields = Array[deB, 3], psi = Array[deV, 8]},
        With[{phi = Normal[MatrixExp[-I deTheta deIsingOperator[couplings, fields], QuantumState[psi, QuantumBasis[{2, 2, 2}]]]["StateVector"]]},
            Simplify[
                {phi - Exp[-I deTheta deIsingEnergies[couplings, fields]] psi, ComplexExpand[Total[Abs[phi] ^ 2] - Total[Abs[psi] ^ 2], psi]},
                Element[Join[Flatten[UpperTriangularize[couplings, 1]], fields, {deTheta}], Reals]
            ]
        ]
    ],
    {ConstantArray[0, 8], 0},
    TestID -> "DiagonalExp-exponential-on-a-pure-register"
]

(* The QAOA cost layer: e^(-i theta H) is the circuit of one R_ZZ(2 theta J_ij) per
   pair of qubits and one R_Z(2 theta b_k) per qubit, which all commute. *)
VerificationTest[
    Table[
        With[{couplings = Array[deJ, {n, n}], fields = Array[deB, n]},
            Simplify[
                Normal[Exp[-I deTheta deIsingOperator[couplings, fields]]["Matrix"]] -
                    Normal[QuantumCircuitOperator[Join[
                        Apply[QuantumOperator["R"[2 deTheta couplings[[#1, #2]], "ZZ"], {#1, #2}] &, Subsets[Range[n], {2}], {1}],
                        MapIndexed[QuantumOperator["RZ"[2 deTheta #1], #2] &, fields]
                    ]]["QuantumOperator"]["Matrix"]]
            ]
        ],
        {n, 2, 4}
    ],
    Table[ConstantArray[0, {2 ^ n, 2 ^ n}], {n, 2, 4}],
    TestID -> "DiagonalExp-Ising-QAOA-gate-factorization"
]

(* For theta up to 10^6, where theta ||H|| reaches 10^7, the exponential of an 8-qubit
   Ising Hamiltonian stays unitary to the last place: every diagonal entry is a phase
   within one machine epsilon of modulus one, and nothing appears off the diagonal. *)
VerificationTest[
    Map[
        With[{u = Exp[-I # deHamiltonian8]["Matrix"]},
            {Max[Abs[Abs[Diagonal[u]] - 1]] <= $MachineEpsilon, AllTrue[u["NonzeroPositions"], Apply[Equal]]}
        ] &,
        {0.3, 1000., 1.*^6}
    ],
    ConstantArray[{True, True}, 3],
    TestID -> "DiagonalExp-long-time-unitarity"
]

(* A rounding test of the stored entries: against their 40-digit exponential, each
   entry is within one rounding of the product theta h_k plus one unit of the result,
   from theta = 0.3 to 10^6. *)
VerificationTest[
    Map[deReferenceError[Normal[Diagonal[Exp[-I # deHamiltonian8]["Matrix"]]], #, deStored8] <= 1 &, {0.3, 7.1, 1000., 1.*^6}],
    ConstantArray[True, 4],
    TestID -> "DiagonalExp-40-digit-reference"
]

(* The same for a 16-qubit Ising Hamiltonian built with QF's diagonal constructor:
   65536 exponentials, within the reference bound, in well under the time limit. *)
VerificationTest[
    deReferenceError[Normal[Diagonal[Exp[-I 0.3 QuantumOperator["Diagonal"[deEnergies16], Range[16]]]["Matrix"]]], 0.3, deEnergies16] <= 1,
    True,
    TimeConstraint -> 30,
    TestID -> "DiagonalExp-16-qubit-Ising"
]

(* The QAOA cost layer e^(-i gamma H) with gamma declared: its closed form is the
   diagonal of e^(-i gamma E_x), gamma = 0 gives the identity, and 25 values on 12
   qubits each fall within the reference bound. *)
VerificationTest[
    With[{
        u4 = Exp[QuantumOperator["Diagonal"[-I deGamma deEnergies12[[;; 16]]], Range[4], "Parameters" -> {deGamma}]],
        u12 = Exp[QuantumOperator["Diagonal"[-I deGamma deEnergies12], Range[12], "Parameters" -> {deGamma}]]
    },
        {
            Normal[u4["Matrix"]] == DiagonalMatrix[Exp[-I deGamma deEnergies12[[;; 16]]]],
            u12[0]["Matrix"] == IdentityMatrix[4096, SparseArray],
            Map[deReferenceError[Normal[Diagonal[u12[#]["Matrix"]]], #, deEnergies12] <= 1 &, Subdivide[0.05, 1.25, 24]]
        }
    ],
    {True, True, ConstantArray[True, 25]},
    TimeConstraint -> 60,
    TestID -> "DiagonalExp-parameter-scan-12-qubit"
]

(* Pure dephasing of two qubits under a diagonal H, all symbolic: the superoperator
   is diagonal with no tolerance, each coherence rho_xy picks up e^(-i (E_x - E_y) t)
   and decays at 2 gamma per differing qubit, the populations and the trace stay, and
   H commutes with every Z_k. *)
VerificationTest[
    With[{rhoT = deDephase[{deGam, deGam}], bits = Tuples[{0, 1}, 2]},
        {
            DiagonalMatrixQ[deLiouvillian[{deGam, deGam}]["Matrix"], Tolerance -> 0],
            Simplify[{
                rhoT - deRho0 Exp[-I deT Outer[Subtract, deEnergies2, deEnergies2] - 2 deGam deT Outer[HammingDistance, bits, bits, 1]],
                Diagonal[rhoT] - Diagonal[deRho0],
                Tr[rhoT] - Tr[deRho0],
                Normal[Commutator[deHamiltonian2, #]["Matrix"]] & /@ deDephasers
            }, deReals]
        }
    ],
    {True, {ConstantArray[0, {4, 4}], ConstantArray[0, 4], 0, ConstantArray[0, {2, 4, 4}]}},
    TestID -> "DiagonalExp-dephasing-symbolic"
]

(* Its limits: at long times the coherences are gone and rho(t) is the pinching
   diag(rho(0)); without dephasing it is the unitary U rho(0) U^dagger. *)
VerificationTest[
    With[{rhoT = deDephase[{deGam, deGam}], u = DiagonalMatrix[Exp[-I deT deEnergies2]]},
        Simplify[{
            Limit[rhoT, deT -> Infinity, Assumptions -> deGam > 0 && Element[{deJ12, deB1, deB2}, Reals]] - DiagonalMatrix[Diagonal[deRho0]],
            (rhoT /. deGam -> 0) - u . deRho0 . ConjugateTranspose[u]
        }, deReals]
    ],
    ConstantArray[0, {2, 4, 4}],
    TestID -> "DiagonalExp-dephasing-limits"
]

(* Pure dephasing multiplies rho entry by entry by c_xy = q^d(x, y), q = e^(-2 gamma t),
   so c is the tensor product of the one-qubit multipliers {{1, q}, {q, 1}}, whose
   eigenvalues 1 - q and 1 + q are nonnegative for gamma t >= 0: c, and with it the
   Choi matrix of the map, is positive semidefinite, and the map completely positive. *)
VerificationTest[
    With[{
        c = Simplify[Normal[Exp[deT QuantumOperator["Liouvillian"[None, deDephasers, {deGam, deGam}]]][deRho0State]["DensityMatrix"]] / deRho0],
        k = {{1, Exp[-2 deGam deT]}, {Exp[-2 deGam deT], 1}}
    },
        {
            Simplify[c - KroneckerProduct[k, k]],
            Simplify[Times @@ (\[FormalX] - Eigenvalues[k]) - (\[FormalX] - (1 - Exp[-2 deGam deT])) (\[FormalX] - (1 + Exp[-2 deGam deT]))]
        }
    ],
    {ConstantArray[0, {4, 4}], 0},
    TestID -> "DiagonalExp-dephasing-complete-positivity"
]

(* Dephasing by the collective Z_1 + Z_2 (a Kossakowski matrix of ones): the
   superoperator is still diagonal, a coherence decays at gamma/2 (s_x - s_y)^2 with
   s the total spin, and |01><10| is decoherence-free. *)
VerificationTest[
    With[{rhoT = deDephase[deGam ConstantArray[1, {2, 2}]], collective = Total /@ deSpins[2]},
        {
            DiagonalMatrixQ[deLiouvillian[deGam ConstantArray[1, {2, 2}]]["Matrix"], Tolerance -> 0],
            Simplify[rhoT - deRho0 Exp[-I deT Outer[Subtract, deEnergies2, deEnergies2] - deGam / 2 deT Outer[(#1 - #2)^2 &, collective, collective]], deReals],
            Simplify[rhoT[[2, 3]] - deRho0[[2, 3]] Exp[-I deT (deEnergies2[[2]] - deEnergies2[[3]])], deReals]
        }
    ],
    {True, ConstantArray[0, {4, 4}], 0},
    TestID -> "DiagonalExp-collective-dephasing"
]

(* A superoperator acting on a state: e^(t L) on rho(0) through MatrixExp[qo, qs] is
   the same rho(t) as Exp[qo][qs]. *)
VerificationTest[
    Simplify[
        Normal[MatrixExp[deT deLiouvillian[{deGam, deGam}], deRho0State]["DensityMatrix"]] - deDephase[{deGam, deGam}],
        deReals
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "DiagonalExp-superoperator-acting-on-state"
]

(* An independent numerical reference: integrating the Lindblad equation with
   QuantumEvolve (NDSolve) at numeric couplings, rate and time gives the same rho(t)
   as Exp[t L], to NDSolve's accuracy. *)
VerificationTest[
    With[{
        h = 0.7 QuantumOperator["ZZ"] - 0.3 QuantumOperator["ZI"] + 0.45 QuantumOperator["IZ"],
        rho0 = QuantumState[BlockRandom[With[{m = RandomComplex[{-1 - I, 1 + I}, {4, 4}]}, m . ConjugateTranspose[m] / Tr[m . ConjugateTranspose[m]]], RandomSeeding -> 2], QuantumBasis[{2, 2}]],
        t1 = 1.3
    },
        Max[Abs[
            Normal[QuantumEvolve[h, deDephasers -> {0.2, 0.2}, rho0, {\[FormalT], 0, t1}][t1]["DensityMatrix"]] -
                Normal[Exp[t1 QuantumOperator["Liouvillian"[h, deDephasers, {0.2, 0.2}]]][rho0]["DensityMatrix"]]
        ]] < 10^-6
    ],
    True,
    TestID -> "DiagonalExp-dephasing-against-QuantumEvolve"
]

(* e^M acting on a mixed state rho gives e^M rho e^(M^dagger): for M = -i theta Z, X
   and Y, and for the non-Hermitian no-jump generator diag(0, -g/2) of amplitude
   damping, whose unnormalized state keeps the trace a + (1 - a) e^(-g), and for the
   entangling M = -i theta ZZ on a generic two-qubit rho; on a pure state it gives
   e^M psi. The reference is MatrixExp of the explicit matrix. *)
VerificationTest[
    With[{
        rho = QuantumState[{{deA, deC}, {Conjugate[deC], 1 - deA}}],
        psi = QuantumState[{deC0, deC1}],
        generators = {-I deTheta QuantumOperator["Z"], -I deTheta QuantumOperator["X"], -I deTheta QuantumOperator["Y"], QuantumOperator[{{0, 0}, {0, -deG / 2}}]}
    },
        FullSimplify[{
            With[{u = MatrixExp[Normal[#["Matrix"]]]}, Normal[MatrixExp[#, rho]["DensityMatrix"]] - u . Normal[rho["DensityMatrix"]] . ConjugateTranspose[u]] & /@ generators,
            Tr[Normal[MatrixExp[Last[generators], rho]["DensityMatrix"]]] - (deA + (1 - deA) Exp[-deG]),
            Normal[MatrixExp[-I deTheta QuantumOperator["X"], psi]["StateVector"]] - MatrixExp[-I deTheta {{0, 1}, {1, 0}}] . {deC0, deC1},
            With[{u = DiagonalMatrix[Exp[-I deTheta {1, -1, -1, 1}]]}, Normal[MatrixExp[-I deTheta QuantumOperator["ZZ"], deRho0State]["DensityMatrix"]] - u . deRho0 . ConjugateTranspose[u]]
        }, Element[{deTheta, deG, deA}, Reals]]
    ],
    {ConstantArray[0, {4, 2, 2}], 0, {0, 0}, ConstantArray[0, {4, 4}]},
    TestID -> "DiagonalExp-exponential-acting-on-states"
]

(* An operator on qubit 2 of a two-qubit state acts there as 1 (x) e^M: a mixed state
   becomes (1 (x) e^M) rho (1 (x) e^M)^dagger and a pure state (1 (x) e^M) psi, for
   M = -i theta X and for the diagonal M = -i theta Z. Dephasing qubit 1 of a Bell
   pair, a superoperator on part of the register, leaves the coherence e^(-2 gamma t)/2
   between 00 and 11, from the mixed and from the pure Bell state. *)
VerificationTest[
    {
        Map[
            With[{m = -I deTheta QuantumOperator[#, {2}], u = KroneckerProduct[IdentityMatrix[2], MatrixExp[-I deTheta Normal[QuantumOperator[#]["Matrix"]]]]},
                FullSimplify[{
                    Normal[MatrixExp[m, deRho0State]["DensityMatrix"]] - u . deRho0 . ConjugateTranspose[u],
                    Normal[MatrixExp[m, QuantumState[Array[deV, 4], QuantumBasis[{2, 2}]]]["StateVector"]] - u . Array[deV, 4]
                }, Element[deTheta, Reals]]
            ] &,
            {"X", "Z"}
        ],
        With[{
            dephaseOne = deT QuantumOperator["Liouvillian"[None, {QuantumOperator["Z", {1}]}, {deGam}]],
            expected = {{1/2, 0, 0, Exp[-2 deGam deT] / 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {Exp[-2 deGam deT] / 2, 0, 0, 1/2}}
        },
            Simplify[Normal[MatrixExp[dephaseOne, #]["DensityMatrix"]] - expected] & /@
                {QuantumState[QuantumState["PhiPlus"]["DensityMatrix"], QuantumBasis[{2, 2}]], QuantumState["PhiPlus"]}
        ]
    },
    {ConstantArray[{ConstantArray[0, {4, 4}], ConstantArray[0, 4]}, 2], ConstantArray[0, {2, 4, 4}]},
    TestID -> "DiagonalExp-exponential-on-part-of-the-register"
]

(* An operator on a qudit the state lacks extends the state, as Exp[qo][qs] does: X on
   qubit 3 of |0> gives cos(theta) |000> - i sin(theta) |001>. A qubit operator on a
   qutrit fails, as Exp[qo][qs] does. *)
VerificationTest[
    {
        FullSimplify[Normal[MatrixExp[-I deTheta QuantumOperator["X", {3}], QuantumState["0"]]["StateVector"]] - (Cos[deTheta] UnitVector[8, 1] - I Sin[deTheta] UnitVector[8, 2]), Element[deTheta, Reals]],
        MatrixExp[-I deTheta QuantumOperator["Z"], QuantumState["0", 3]]
    },
    {ConstantArray[0, 8], $Failed},
    {QuantumCircuitOperator::dim},
    TestID -> "DiagonalExp-exponential-beyond-the-register"
]

(* e^(-i g ZZ) entangles |++>: the concurrence of the result is |sin 2g|. *)
VerificationTest[
    FullSimplify[QuantumEntanglementMonotone[Exp[-I deG QuantumOperator["ZZ"]][QuantumState["++"]], {1}, "Concurrence"] - Abs[Sin[2 deG]], Element[deG, Reals]],
    0,
    TestID -> "DiagonalExp-ZZ-entangler-concurrence"
]

(* The spin-1 J_x is stored diagonal in its own eigenbasis, so its exponential goes
   entry by entry there; read in the computational basis it is the exponential of the
   3 x 3 matrix J_x. *)
VerificationTest[
    {
        DiagonalMatrixQ[QuantumOperator["JX"[1]]["Matrix"], Tolerance -> 0],
        FullSimplify[
            Normal[Exp[-I deTheta QuantumOperator["JX"[1]]]["MatrixRepresentation"]] - MatrixExp[-I deTheta {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}} / Sqrt[2]],
            Element[deTheta, Reals]
        ]
    },
    {True, ConstantArray[0, {3, 3}]},
    TestID -> "DiagonalExp-spin-one-JX-eigenbasis"
]

(* The open Ising chain H = -J sum_k Z_k Z_(k+1) on n qubits has the partition function
   Tr e^(-beta H) = 2 (2 cosh(beta J))^(n - 1), for symbolic beta J and n = 2..5. *)
VerificationTest[
    Table[Simplify[TrigToExp[Tr[Normal[Exp[deBetaJ deChain[n]]["Matrix"]]] - 2 (2 Cosh[deBetaJ]) ^ (n - 1)]], {n, 2, 5}],
    ConstantArray[0, 4],
    TestID -> "DiagonalExp-Ising-chain-partition-function"
]

(* The Gibbs correlations of the open 5-qubit chain: <Z_1 Z_r> = tanh(beta J)^(r - 1),
   for symbolic beta J. *)
VerificationTest[
    With[{w = Normal[Diagonal[Exp[deBetaJ deChain[5]]["Matrix"]]], s = deSpins[5]},
        Simplify[Table[Total[w s[[All, 1]] s[[All, r]]] / Total[w] - Tanh[deBetaJ] ^ (r - 1), {r, 2, 5}]]
    ],
    ConstantArray[0, 4],
    TestID -> "DiagonalExp-Ising-chain-correlations"
]

(* As beta J -> Infinity the ferromagnetic 4-qubit chain is left in its two ground
   states |0000> and |1111>, with weight 1/2 each. At beta J = 250 the same weights
   come through an overflow: the ground states' Boltzmann factor e^750 is beyond the
   largest machine number and the factor e^-750 of the most excited states is below
   the smallest. Those weights are normalized at 40 digits, since machine arithmetic
   on factors this far apart underflows. *)
VerificationTest[
    {
        With[{w = Normal[Diagonal[Exp[deBetaJ deChain[4]]["Matrix"]]]}, Limit[w / Total[w], deBetaJ -> Infinity]],
        With[{w = SetPrecision[Normal[Diagonal[Exp[250. deChain[4]]["Matrix"]]], 40]},
            With[{p = w / Total[w]}, {Max[Abs[p[[{1, 16}]] - 1/2]] < 10^-14, Total[p[[2 ;; 15]]] < 10^-200}]
        ]
    },
    {ReplacePart[ConstantArray[0, 16], {1 -> 1/2, 16 -> 1/2}], {True, True}},
    TestID -> "DiagonalExp-Gibbs-ground-states"
]

(* The Gibbs weights e^(-E_x) / Z do not change when a constant is added to the
   energies. At energies 720 and 721 the Boltzmann factors are subnormal numbers,
   spaced 6.6 10^-11 of their size near e^-721, and their weights still match those of
   energies 0 and 1 within 10^-9. *)
VerificationTest[
    With[{weights = (# / Total[#] &) @ Normal[Diagonal[Exp[-QuantumOperator[DiagonalMatrix[#]]]["Matrix"]]] &},
        Max[Abs[weights[{720., 721.}] - weights[{0., 1.}]]] < 10^-9
    ],
    True,
    TestID -> "DiagonalExp-Gibbs-shift-invariance"
]

(* QAOA at depth one on the 4-ring: the cost layer e^(-i gamma C), C = sum over the
   edges of (1 - Z_i Z_j)/2, is diagonal; after the mixer e^(-i beta sum_k X_k) the
   edge term <(1 - Z_1 Z_2)/2> is (2 + sin 4 beta sin 2 gamma)/4, for symbolic gamma
   and beta. *)
VerificationTest[
    With[{ring = Total[QuantumOperator["ZZ", #] & /@ {{1, 2}, {2, 3}, {3, 4}, {1, 4}}], s = deSpins[4]},
        With[{psi = QuantumCircuitOperator[Table[QuantumOperator["RX"[2 deBeta], {k}], {k, 4}]][Exp[I deGamma / 2 ring][QuantumState["++++"]]]},
            FullSimplify[
                TrigReduce[Total[ComplexExpand[Abs[Normal[psi["StateVector"]]] ^ 2] (1 - s[[All, 1]] s[[All, 2]]) / 2]] - (2 + Sin[4 deBeta] Sin[2 deGamma]) / 4,
                Element[{deGamma, deBeta}, Reals]
            ]
        ]
    ],
    0,
    TestID -> "DiagonalExp-QAOA-depth-one-ring"
]

(* The same edge term from the cost layer with gamma declared, read at 25 values of
   gamma at beta = pi/8, where it is (2 + sin 2 gamma)/4 and reaches 3/4 at
   gamma = pi/4. *)
VerificationTest[
    With[{
        s = deSpins[4],
        cost = Exp[QuantumOperator["Diagonal"[I deGamma / 2 Total[deSpins[4][[All, #1]] deSpins[4][[All, #2]] & @@@ {{1, 2}, {2, 3}, {3, 4}, {1, 4}}]], Range[4], "Parameters" -> {deGamma}, "Label" -> None]],
        mixer = QuantumCircuitOperator[Table[QuantumOperator["RX"[Pi / 4], {k}], {k, 4}]],
        gammas = Subdivide[0.05, 1.25, 24]
    },
        With[{edge = Total[Abs[Normal[mixer[cost[#][QuantumState["++++"]]]["StateVector"]]] ^ 2 (1 - s[[All, 1]] s[[All, 2]]) / 2] &},
            {Max[Abs[edge /@ gammas - (2 + Sin[2 gammas]) / 4]] < 10^-14, Abs[edge[N[Pi / 4]] - 3/4] < 10^-14}
        ]
    ],
    {True, True},
    TestID -> "DiagonalExp-QAOA-scan-with-gamma-declared"
]

(* Symbolic and exact diagonals keep their closed form. *)
VerificationTest[
    {
        Normal[Exp[QuantumOperator[DiagonalMatrix[{deA, deB, deC, deD}], {1, 2}]]["Matrix"]],
        Normal[Exp[QuantumOperator[DiagonalMatrix[{I Pi, Log[2], 0, -I Pi / 2}], {1, 2}]]["Matrix"]]
    },
    {DiagonalMatrix[Exp[{deA, deB, deC, deD}]], DiagonalMatrix[{-1, 2, 1, -I}]},
    TestID -> "DiagonalExp-symbolic-and-exact"
]

(* An off-diagonal entry that the zero test proves zero leaves a diagonal
   exponential, with no pole where the two diagonal entries coincide. *)
VerificationTest[
    Normal[Exp[QuantumOperator[{{deA, Sin[deY]^2 + Cos[deY]^2 - 1}, {0, deB}}]]["Matrix"]] /. deB -> deA,
    {{E^deA, 0}, {0, E^deA}},
    TestID -> "DiagonalExp-hidden-zero-no-pole"
]

(* An off-diagonal entry that is zero only through Gamma[g + 1] = g Gamma[g] is not
   proved zero by DiagonalMatrixQ, so the matrix goes to MatrixExp; the result is
   still the diagonal of exponentials. *)
VerificationTest[
    With[{m = {{deA, Gamma[deG + 1] - deG Gamma[deG]}, {0, deB}}},
        {DiagonalMatrixQ[m, Tolerance -> 0], FullSimplify[Normal[Exp[QuantumOperator[m]]["Matrix"]] - DiagonalMatrix[{E^deA, E^deB}]]}
    ],
    {False, ConstantArray[0, {2, 2}]},
    TestID -> "DiagonalExp-hidden-zero-unproved-falls-back"
]

(* e^(-2000) underflows to zero, and e^0 is exactly the machine number 1. *)
VerificationTest[
    With[{m = Normal[Exp[-2000. QuantumOperator[DiagonalMatrix[{1., 0.}]]]["Matrix"]]},
        {m[[1, 1]] == 0, m[[2, 2]] - 1 == 0, Precision[m[[2, 2]]]}
    ],
    {True, True, MachinePrecision},
    TestID -> "DiagonalExp-underflow"
]

(* e^800 is beyond the largest machine number, so that entry comes back as an
   arbitrary-precision number, compared with e^800 at 40 digits, while e^(-800)
   underflows to zero, both silently. *)
VerificationTest[
    With[{m = Normal[Exp[-800. QuantumOperator[DiagonalMatrix[{-1., 1.}]]]["Matrix"]]},
        {Abs[SetPrecision[m[[1, 1]], 40] / Exp[800] - 1] < 10^-14, m[[2, 2]] == 0}
    ],
    {True, True},
    TestID -> "DiagonalExp-overflow"
]

(* e^(-708 + 1.5 i) has a modulus just above the smallest normalized machine number
   and a subnormal real part: it comes back silently, correct to within one unit of
   its modulus, and e^0 beside it is exactly 1. The error relative to the modulus is
   computed at 30 digits, since machine arithmetic on numbers this small underflows. *)
VerificationTest[
    With[{m = Normal[Exp[QuantumOperator[DiagonalMatrix[{-708. + 1.5 I, 0.}]]]["Matrix"]]},
        {Abs[SetPrecision[m[[1, 1]], 30] - N[Exp[-708 + 3/2 I], 30]] / N[Exp[-708], 30] <= $MachineEpsilon, m[[2, 2]] - 1 == 0}
    ],
    {True, True},
    TestID -> "DiagonalExp-complex-subnormal-part"
]

(* A diagonal at 30 digits, exact entries among it, is exponentiated at that
   precision: e^(-i x) for x = 1, 2, 1/3 and Pi agrees with the exact value within the
   precision the result carries. *)
VerificationTest[
    With[{d = Normal[Diagonal[Exp[-I QuantumOperator[DiagonalMatrix[{1`30, 2`30, 1/3, Pi}], {1, 2}]]["Matrix"]]]},
        {Precision[d] > 29, Max[Abs[d - Exp[-I {1, 2, 1/3, Pi}]]] == 0}
    ],
    {True, True},
    TestID -> "DiagonalExp-thirty-digit-diagonal"
]

(* A diagonal with an infinite or indeterminate entry has no exponential: it goes
   to MatrixExp and fails, as it did before. *)
VerificationTest[
    Head /@ {Exp[QuantumOperator[DiagonalMatrix[{Indeterminate, 1.}]]], Exp[QuantumOperator[DiagonalMatrix[{-Infinity, 1.}]]]},
    {Failure, Failure},
    TestID -> "DiagonalExp-non-finite-diagonal-fails"
]

(* A machine diagonal with exact entries among it, integers or constants such as Pi,
   is exponentiated in machine numbers, each entry within one unit of the 30-digit
   exponential of its input, and an underflowing entry gives zero without a
   message. *)
VerificationTest[
    Map[
        With[{d = Normal[Diagonal[Exp[QuantumOperator[#]]["Matrix"]]], reference = Exp[SetPrecision[Diagonal[#], 30]]},
            {Precision[d], Max[Abs[SetPrecision[d, 30] - reference] / Clip[Abs[reference], {1, Infinity}]] <= $MachineEpsilon}
        ] &,
        {{{1., 0}, {0, 2}}, {{-800., 0}, {0, 0}}, {{1., 0}, {0, Pi}}, {{-800., 0}, {0, Pi}}}
    ],
    ConstantArray[{MachinePrecision, True}, 4],
    TestID -> "DiagonalExp-mixed-exact-and-machine"
]

(* The off-diagonal a b vanishes only at b = 0, so the operator keeps the full
   matrix; a substitution that makes it diagonal is exponentiated entry by entry, one
   that does not by MatrixExp. *)
VerificationTest[
    With[{u = Exp[QuantumOperator[{{deA, deA deB}, {0, 1}}, "Parameters" -> {deA, deB}]]},
        {
            Max[Abs[Normal[u[0.5, 0.]["Matrix"]] - DiagonalMatrix[Exp[{0.5, 1.}]]]] == 0,
            Max[Abs[Normal[u[0.5, 0.2]["Matrix"]] - MatrixExp[{{0.5, 0.1}, {0, 1}}]]] == 0
        }
    ],
    {True, True},
    TestID -> "DiagonalExp-substitution-decides-route"
]

(* A coupling of 10^-9 between two degenerate levels, beside an entry 10^6, is at
   roundoff relative to the largest entry, yet over t = 10^6 it rotates the pair by
   10^-3. The exponential must not treat the matrix as diagonal. *)
VerificationTest[
    With[{m = {{1.*^6, 0., 0.}, {0., 0., 1.*^-9}, {0., 1.*^-9, 0.}}, t = 1.*^6},
        Abs[Abs[Normal[Exp[-I t QuantumOperator[m, 3]]["Matrix"]][[2, 3]]] - Sin[10^-3]] < 10^-4
    ],
    True,
    TestID -> "DiagonalExp-resonant-roundoff-coupling-kept"
]

(* Two parameters that coincide: the exponential of diag(a, b) at a = b = 1 is e I. *)
VerificationTest[
    Normal[Exp[QuantumOperator[DiagonalMatrix[{deA, deB}], "Parameters" -> {deA, deB}]][1, 1]["Matrix"]],
    E IdentityMatrix[2],
    TestID -> "DiagonalExp-parameter-collision"
]

(* A matrix-type object with no qudits keeps its 1 x 1 matrix. *)
VerificationTest[
    Normal[QuantumPartialTrace[QuantumState[{{3/10, 0}, {0, 7/10}}]]["Matrix"]],
    {{1}},
    TestID -> "DiagonalExp-rank-zero-matrix"
]

EndTestSection[]
