# Exponential of a diagonal operator: entry by entry instead of MatrixExp

Investigation and prototype, 2026-09-27 and 28. WL 15.0.1, macOS ARM, 12 cores, with other kernels running throughout. Nothing was committed and no tracked file was edited. The prototype lives in scratch copies of the paclet exported from main. Its diffs, a script that applies them to a fresh export, the new test file, the scripts that produced every number below and their raw outputs are in `exp-diagonal-shortcut-prototype/` beside this report; section 12 lists them.

Main moved during the work, from `2ff7d0ce` through `81df7633` and `b13a9077` to `a1799b61`. The final prototype and every measurement below are on `a1799b61`, whose `0 ^ qo` rule the prototype leaves as it is. Earlier versions and their logs are in `exp-diagonal-shortcut-prototype/earlier/`. Main has since moved to `bdbe7165`, where `1f95b354` made the substitution guard read the stored values of a `SparseArray` (`valuelessEntriesQ`, `Utilities.m:361-363`). The prototype's guard edit (section 5.2) therefore no longer applies and is no longer needed. With it dropped, every other edit applies to `bdbe7165` and the 38 tests pass there (`file-bdbe7165.out`); on that main, `matrixExponential` would use `! valuelessEntriesQ` in place of `finiteArrayQ`.

Verification status (section 11): three rounds of `/wl-verify` ran on the `81df7633` version and did not converge. Their findings are fixed, but the code changed after the last round, so its correctness rests on the full suite, the 38 new tests and the correctness battery, not on a fresh verifier. `/wl-quality` ran its three rounds and did not converge: the third reported 13 open issues. Eleven were fixed after it, among them one change to the kernel code (`MatrixExp[qo, qs]` for an operator on other qudits than exactly the state's); one, the `f[qo]` diagonal builder, is left with decision 1; and one is the fact that the physics brief was written after the code. Nothing changed after the third round has been reviewed by a fresh critic.

**Decisions needed from Mads:**
1. The detection rule. For the exponential the prototype treats a matrix as diagonal only when every off-diagonal entry is exactly zero (`Tolerance -> 0`), because `DiagonalMatrixQ`'s default tolerance drops a coupling that `MatrixExp` resolves (section 3). The existing `f[qo]` branch keeps its default-tolerance test, since sharing the strict one would slow `f[qo]` without changing its results. So there are two tests of "diagonal", one per route, and they disagree where it matters: on section 3's resonant matrix A, e^(iA) and cos A + i sin A differ by 10⁻³, on main as on the prototype, because `Sin` and `Cos` drop the coupling that `Exp` keeps. One rule could serve every function: drop an off-diagonal entry m_ij only when its first-order effect on f(M), |m_ij| times the divided difference |f[m_ii, m_jj]| (f'(m_ii) for degenerate entries), is below the roundoff of f of the diagonal. That rule keeps the resonant coupling, whose first-order effect is 10⁻³, and drops a roundoff coupling between well-separated levels, whose divided difference is at most 2/|m_ii - m_jj| for a unitary e^(-iθM). It would not restore Euler's formula by itself, because `f[qo]`'s Schur route also sets the eigenvalues ±10⁻³ of that matrix to zero (open question 5), so it belongs with a change to `f[qo]` and is not in this prototype.
2. `MatrixExp[qo, qs]`. Today it disagrees with `Exp[qo][qs]` on a mixed state, where it does not return e^M ρ e^(M†), and for any operator not on exactly the state's qudits (section 2); the prototype makes the two agree. That changes outputs, so it may belong in its own commit.

## 1. Answer

Yes. The exponential of a diagonal operator is the diagonal of the exponentials of its entries, and QF already stores these operators with off-diagonal entries that are exactly zero, so recognizing one costs about 0.1 ms for a 2^16 `SparseArray`. For the 16-qubit Ising exponential the prototype replaces a 124 s `MatrixExp` with a 0.034 s map over 65536 entries (section 4.1).

The result is also closer to the exact exponential. Against the 40-digit exponential of the stored entries, each entry of the prototype's result lies within one rounding of the product θh_k plus one unit of the result for θ up to 10⁶, where θ‖H‖ reaches 10⁷, and the modulus of every diagonal entry of e^(-iθH) stays within one machine number of 1. `MatrixExp` on the same diagonal misses that bound by up to a factor 47, and its moduli drift from 1 by up to 5e-10 at θ = 10⁶ (section 4.4). Non-diagonal operators take today's route and give results identical to the last bit (section 6). The prototype passes the full suite on `a1799b61` (2892 of 2892: main's 2854 tests plus the 38 new ones), and 13 of the 38 new tests fail on current main.

Two further changes carry the gain to the cases where QF was slowest:

- **Declared parameters.** For an operator with declared `"Parameters"`, QF holds the whole matrix as a dense `List` inside a function of the parameters: at 12 qubits the exponential occupies 412 MB and each substitution takes 23.5 s. Holding only the diagonal gives 0.56 MB and 0.014 s (section 4.3).
- **Matrix-type operators** such as the Liouvillian are stored as a `SparseArray`, but `"StateMatrix"` rebuilds their matrix as a dense `List`, because it wraps the stored tensor in braces before reshaping it. Keeping it sparse makes the n = 6 dephasing exponential 51 times faster on top of the entrywise exponential, and building that Liouvillian 2.6 times faster (section 4.2). This is a separate small change (section 5.4); the full suite passes with it too.

Three defects turned up on the way:

- `MatrixExp[qo, qs]` on a mixed state is not e^M ρ e^(M†), and for an operator not on exactly the state's qudits it returns an invalid state or ignores the operator's order (section 2). The prototype fixes both, because its diff rewires that code.
- `qo["Eigenvalues"]` and `qo["Eigensystem"]` of a diagonal operator convert the sparse matrix to a dense one, print `Eigensystem::arh` doing so at 8 and 10 qubits, and take seconds at 12 qubits (section 8.2).
- QF's `"Diagonal"` constructor labels the operator with its whole diagonal. At 16 qubits the exponential takes 0.20 s with that label and 0.015 s without it (section 8.5).

## 2. Where QF exponentiates an operator matrix

Line numbers are for main at `a1799b61`.

| Entry point | Location | Matrix exponentiated |
|---|---|---|
| `base ^ qo` (and `E ^ qo`) | `QuantumOperator.m:587-588`, through `scalarBasePower` (`:583-585`) and `matrixMapOperator` (`:604-611`) | `Log[base] op["Matrix"]` with `op = qo["Sort"]`; a zero base takes the limit route `zeroBasePower` (`Utilities.m:422-491`), which the prototype leaves as it is |
| `Exp[qo]` | `QuantumOperator.m:595` | delegates to `E ^ qo`, so the same matrix |
| `MatrixExp[qo]` | `QuantumOperator.m:600` | `op["Matrix"]`, through `matrixMapOperator` |
| `MatrixExp[qo, qs]` | `QuantumOperator.m:689-703` | for a vector-type operator on a pure state, `op["Matrix"]` acting on the state vector; otherwise `op["ToMatrix"]["Matrix"]` acting on the density vector |
| The spellings above with declared parameters | `QuantumOperator.m:629-687` (`lazyMatrixMapQ`, `heldOperatorMatrix` at `:657` and `:664`, `parametricMatrixMapOperator` at `:672`) | the same map, applied at each substitution to `heldOperatorMatrix`, which is `Normal[op["Matrix"]]`, or for a nested map the inner result read back as a matrix |
| `"RX"`, `"RY"`, `"RZ"` | `NamedOperators.m:163-185` | `Exp[-I angle/2 QuantumOperator["Pauli*"[d]]]`; `"RZ"`, on a qubit or a qudit, is diagonal |
| `"R"[angle, ops...]` | `NamedOperators.m:187-195` | `Exp[-I angle/2 op]`, diagonal for any string of `Z` and `I` |
| `"Trotterization"` | `NamedCircuits.m:597` | `Exp[-I c op]` for each term; the `ZZ` terms of an Ising Trotter circuit are diagonal |
| Adiabatic and parameter-shift helpers | `QuantumOptimization.m:500-558` | `Exp[I s QuantumOperator[...]]` |
| `ToBosonicOperator` | `SecondQuantization/Utils.m:144` | `MatrixExp[Log[base] x]` with `x` a `QuantumOperator`, so the `MatrixExp[qo]` rule |
| `QuantumEvolve` | `QuantumEvolution.m:150-186` | no exponential: `NDSolveValue` for a numeric time range, `DSolveValue` for symbolic time; only a time-independent phase-space rate matrix is exponentiated directly (`:179`), and those rate matrices are not diagonal for a diagonal H |
| `"EvolutionOperator"`, `"NEvolutionOperator"`, `QuantumChannel` `"EvolutionChannel"` | `QuantumOperator/Properties.m:643-645`, `QuantumChannel/Properties.m` | through `QuantumEvolve` |

`op["Matrix"]` is `op["Sort"]["StateMatrix"]`. For a vector-type operator it is the stored amplitudes reshaped to output by input. For a matrix-type operator (`"StateType"` `"Matrix"`, such as the Liouvillian) it is the d² by d² superoperator rebuilt from the stored tensor (`QuantumState/Properties.m:1051-1058`). The prototype tests and exponentiates exactly this matrix at every call site, so the superoperator is what gets exponentiated, as today. The dephasing Liouvillian's superoperator is diagonal whenever H and every jump operator are diagonal.

On base, `MatrixExp` is where the time goes. For the 12-qubit Ising exponential (0.16 s), `MatrixExp` itself took 0.15 s when the stages were profiled on `2ff7d0ce`, and scaling, sorting, reading the matrix and rebuilding the operator took the rest. For the n = 5 Liouvillian, reading the matrix takes 0.012 s of the 0.31 s exponential on `a1799b61`.

### `MatrixExp[qo, qs]` on a mixed state and outside the whole register

For a vector-type operator M and a mixed state ρ, `MatrixExp[qo, qs]` exponentiates `op["ToMatrix"]["Matrix"]`, and `"ToMatrix"` turns the operator into the projector onto its own vectorization, so the result is exp(vec(M) vec(M)†) applied to vec ρ. For M = -iθZ and ρ = ((a, b), (b*, 1 - a)) the kernel returns ((a e^(θ²), b e^(-θ²)), (b* e^(-θ²), (1 - a) e^(θ²))): the trace grows as e^(θ²). `Exp[qo][qs]` gives the physical e^M ρ e^(M†). For an operator on only some of the state's qudits, such as X on qubit 2 of a two-qubit state, `MatrixExp[qo, qs]` prints `MatrixExp::lslc` and `ArrayReshape::listrp` and returns a `QuantumState` whose stored state is an unevaluated `ArrayReshape[MatrixExp[...], ...]` for a mixed state and `MatrixExp[...]` for a pure one, while `Exp[qo][qs]` gives (1 ⊗ e^M) ρ (1 ⊗ e^M)†. For an operator on a qudit the state lacks, X on qubit 3 of the one-qubit |0⟩, it ignores the order and returns the one-qubit state cos θ |0⟩ - i sin θ |1⟩, where `Exp[qo][qs]` extends the state and returns cos θ |000⟩ - i sin θ |001⟩. For a qubit operator on a qutrit it returns another invalid state, after `MatrixExp::lslc`, where `Exp[qo][qs]` fails with `QuantumCircuitOperator::dim`. No test in the suite calls the two-argument form. The prototype makes `MatrixExp[qo, qs]` give what `Exp[qo][qs]` gives (section 5.2): u ρ u† with u = e^M for a mixed state, and `Exp[qo][qs]` itself whenever the operator's qudits and dimensions are not exactly the state's.

## 3. The detection rule and why it has no tolerance

### What `DiagonalMatrixQ`'s default tolerance does

Measured in the installed kernel (`tolerance.wls`), and checked against the `DiagonalMatrixQ` doc page:

- The default tolerance is relative to the **largest entry of the whole matrix**, not to the entries a coupling connects: an off-diagonal 5.6e-14 counts as zero next to a largest entry 2, and a 1e-9 coupling between two entries of size 1 counts as zero when an entry 1e6 sits elsewhere in the matrix.
- An explicit `Tolerance -> t` is an absolute threshold ("all entries `Abs[m_ij] <= t` are taken to be zero"), so `Tolerance -> 0` accepts only entries equal to zero. The machine zeros `0.`, `-0.` and `0. + 0. I` count as zero. A symbolic entry counts as zero when the zero test proves it, as it does for `Sin[y]^2 + Cos[y]^2 - 1`; that test is best effort, and it does not prove `Gamma[g + 1] - g Gamma[g]`.
- The doc page says `DiagonalMatrixQ` works for a non-square matrix, and it does, so the test must also require `SquareMatrixQ`.

### A coupling the default tolerance drops but MatrixExp keeps

Take levels 2 and 3 degenerate at energy 0 with a coupling 1e-9 between them, beside a level at 1e6: `M = {{1e6, 0, 0}, {0, 0, 1e-9}, {0, 1e-9, 0}}`. The coupling is 1e-15 of the largest entry, so the default tolerance calls `M` diagonal. Over `t = 1e6` the coupling rotates the degenerate pair by `1e-9 t = 1e-3`, and that rotation is the whole of the dynamics in that block. The exact 2-3 block of e^(-iMt) is cos(10⁻³) 1 - i sin(10⁻³) X (`resonant.wls`):

| Route | largest error in the 2-3 block |
|---|---|
| `MatrixExp` of the full matrix | 2.5e-5 |
| QF's `Exp[-I t qo]`, base and prototype | 4.8e-5 |
| entry by entry (drops the coupling) | 1.0e-3 |
| `MatrixExp` of the 2 x 2 block alone | 1.1e-16 |

`MatrixExp` resolves the rotation, to within the error that the large entry (10¹² in Mt) costs it; the entrywise route loses the rotation altogether. Dropping a coupling ε changes e^(cM) by about |c|ε, which stays below `MatrixExp`'s own error only when the coupled levels are far apart compared with ε. Between degenerate levels the error grows with t while `MatrixExp` keeps resolving the mixing. So a tolerance relative to the largest entry is not safe for the exponential.

### The rule the prototype uses

`diagonalMatrixQ[mat] := SquareMatrixQ[mat] && DiagonalMatrixQ[mat, Tolerance -> 0]`, in `Utilities.m`, used by the exponential and by the decision, for an operator with declared parameters, to hold only the diagonal (section 5.3). It never changes the physics: a matrix it accepts has exactly zero off-diagonal entries, so its exponential is exactly diagonal. The diagonal operators tried here all pass it: a diagonal `SparseArray`, QF's `"Diagonal"` constructor, a sum of `"ZZ"` and `"Z"` operators, the dense Liouvillian superoperator (its zeros are `0. + 0. I`), `"PauliZ"[3]` (the clock matrix), and the spin-1 `"JX"`, which QF stores diagonal in its own eigenbasis. A matrix that is diagonal only up to roundoff, such as `U.D.U†` computed numerically, goes to `MatrixExp` as today: no gain and no harm. An identically zero symbolic entry that the zero test cannot prove also goes to `MatrixExp`, whose result is correct but keeps the removable pole (e^a - e^b)/(a - b).

The `f[qo]` diagonal branch (`Utilities.m:331`) keeps `DiagonalMatrixQ`'s default tolerance. A version of the prototype that shared the strict test with it sent roundoff-contaminated diagonals to the Schur route: `Sin` of a 10-qubit diagonal with a 1e-17 coupling went from 2.9 ms to between 0.1 and 1.2 s (two measurements on earlier versions of main, one of them the round-3 verifier's; that version of the prototype no longer exists), and to seconds at 12 qubits. It gained nothing. On the resonant matrix above, `Cos[qo]` loses the 10⁻³ rotation through either route (error 5.0e-7 in the cosine on base and prototype): the diagonal branch drops the coupling, and `roundoffEigenvalues` (`Utilities.m:359`), which both that branch (`:335`) and the Schur route (`:373`) apply, sets any eigenvalue within 10 machine epsilons of zero, relative to the largest, to zero. It does so by design, so that `Sqrt` and `Log` stay finite. That is a separate, pre-existing question about `f[qo]` (open question 5).

### What detection costs

`DiagonalMatrixQ[m, Tolerance -> 0]` stops at the first nonzero off-diagonal entry (`bench.wls`, group `detect`, and `liouvillian-matrix.wls`; best of three):

| Matrix | strict test |
|---|---|
| dense random 1024 x 1024 (not diagonal) | below 0.1 ms |
| sparse transverse field, 12 qubits (not diagonal) | below 0.1 ms |
| symbolic dense 64 x 64 (not diagonal) | below 0.1 ms |
| `SparseArray` diagonal, 2^16 | 0.1 ms |
| symbolic 64 x 64 diagonal whose 4032 off-diagonal entries are all `Sin[y]^2 + Cos[y]^2 - 1` | 0.7 ms |
| unpacked complex `List` diagonal, 4096 x 4096 | 0.36 s |
| the n = 6 dephasing Liouvillian's `"Matrix"` on the prototype: an unpacked `List`, 4096 x 4096 | 0.26 s |
| the same with section 5.4: a `SparseArray` | 7 µs |

The test is negligible when a nonzero off-diagonal entry comes early, and it reads every entry when the matrix is diagonal. The only case that costs noticeably is a diagonal stored as an unpacked dense `List`, which is what `"StateMatrix"` returns for a matrix-type operator; there the entrywise exponential saves far more than the scan costs, and section 5.4 removes the scan. The non-diagonal controls in section 4.5 take the same time on base and prototype.

## 4. Measured gain

A time marked "one run" is a single run after `ClearSystemCache[]`, kept together with its value; every other time is the best of three runs after `ClearSystemCache[]` when one run takes under 5 s, else that one run. Times are in seconds. "Base" is main at `a1799b61`; "prototype" is main with sections 5.1-5.3; "prototype + sparse StateMatrix" adds section 5.4. The scripts are `bench.wls` and `accuracy.wls`; they build every operator with QF's constructors, and they drop the label of `"Diagonal"` operators, whose cost is section 8.5's subject, not the exponential's.

### 4.1 Diagonal Ising exponential, `Exp[-I 0.3 qo]`

`qo = QuantumOperator["Diagonal"[h], Range[n], "Label" -> None]`, with h the energies Σ J_ij z_i z_j of the basis states and J_ij uniform in [-1, 1].

| n | base, `Exp` (one run) | prototype, `Exp` (one run) | base, `MatrixExp[qo, qs]` | prototype, `MatrixExp[qo, qs]` |
|---|---|---|---|---|
| 12 | 0.163 | 0.024 | 0.010 | 0.0071 |
| 13 | 0.582 | 0.024 | 0.037 | 0.0076 |
| 14 | 2.21 | 0.025 | 0.020 | 0.0094 |
| 15 | 22.6 | 0.028 | 0.028 | 0.013 |
| 16 | 124 | 0.034 | 0.049 | 0.019 |

`MatrixExp` does not use the structure of a sparse diagonal, so the base cost grows much faster than the number of basis states; the prototype's cost tracks it. Every base result lies within 1.1e-15 of `Exp` of the energies, and the prototype's is exactly that; section 4.4 measures both against an independent reference. `MatrixExp[qo, qs]`, which applies the exponential to a state, was already fast on base and gains less.

### 4.2 Pure-dephasing Liouvillian, `Exp[0.7 L]`

`L = QuantumOperator["Liouvillian"[H, {Z_1, ..., Z_n}, {0.2, ...}]]` with H the diagonal Ising operator of section 4.1 and each Z_k the Pauli string with `Z` on qubit k; the superoperator acts on 4^n entries.

| n | 4^n | base (one run) | prototype (one run) | prototype + sparse StateMatrix (one run) |
|---|---|---|---|---|
| 4 | 256 | 0.018 | 0.011 | 0.0083 |
| 5 | 1024 | 0.307 | 0.043 | 0.0085 |
| 6 | 4096 | 16.3 | 0.52 | 0.010 |
| 7 | 16384 | not run | not run | 0.013 |

At n = 7, base and prototype would first build `"Matrix"` as a dense `List` of 4^14 complex entries, at least 4 GB, so they were not run. Every result lies within 5.6e-16 of `Exp` of the superoperator's diagonal on base, and is exactly that on the prototype. Reading `"Matrix"` alone at n = 6 takes 0.14 s as a dense `List` and 0.0027 s as a `SparseArray`.

Building `L` itself, one run each:

| n | base | prototype | prototype + sparse StateMatrix |
|---|---|---|---|
| 4 | 1.33 | 1.33 | 1.30 |
| 5 | 2.36 | 2.32 | 1.97 |
| 6 | 8.92 | 9.13 | 3.54 |
| 7 | not run | not run | 9.59 |

At n = 6 the prototype's exponential spends most of its 0.52 s building the dense `List` that `"StateMatrix"` returns and scanning it; with the sparse `"StateMatrix"` both disappear, and building `L` becomes the whole cost of an open-system evolution of this kind.

### 4.3 Parameter scan: `Exp` of a diagonal with a declared parameter

The QAOA cost layer e^(-iγH), with γ declared (`Exp[QuantumOperator["Diagonal"[-I γ h], Range[n], "Parameters" -> {γ}, "Label" -> None]]`), read at 25 or 100 values of γ in [0.05, 1.25].

| n | copy | build (one run) | `ByteCount` | one value | 25 values | 100 values |
|---|---|---|---|---|---|---|
| 8 | base | 0.41 | 1.7 MB | 0.051 | 1.44 | 6.10 |
| 8 | prototype | 0.11 | 74 kB | 0.0044 | 0.11 | 0.46 |
| 10 | base | 4.89 | 25.9 MB | 0.91 | 24.2 | 95.6 |
| 10 | prototype | 0.062 | 176 kB | 0.0062 | 0.16 | 0.65 |
| 12 | base | 76.0 | 412 MB | 23.5 | not run | not run |
| 12 | prototype | 0.100 | 557 kB | 0.014 | 0.35 | 1.39 |
| 14 | prototype | 0.25 | 2.0 MB | 0.044 | 1.09 | 4.39 |
| 16 | prototype | 0.83 | 7.9 MB | 0.16 | 4.04 | 16.9 |

Base holds `Normal[op["Matrix"]]`, a dense list of 4^n entries, inside the function of γ, and every substitution copies it and runs `MatrixExp` on it. The scans were skipped where one value takes over 2 s, and base was not run at n = 14, where the held list would take about 6.6 GB (16 times the n = 12 figure). The prototype holds the 2^n diagonal entries, and each substitution exponentiates them. At γ = 0.3 the base result lies within 6.2e-16 of `Exp` of the energies, and the prototype's is exactly that.

### 4.4 Accuracy against an independent reference

8-qubit Ising H = Σ J_ij Z_i Z_j built from QF's `"ZZ"` operators (`accuracy.wls`). Unitarity is measured by the largest deviation of the modulus of a diagonal entry of e^(-iθH) from 1, in units of ε = `$MachineEpsilon`:

| θ | 0.3 | 1000 | 10⁶ |
|---|---|---|---|
| base | 1.5 | 2601 | 2.4 × 10⁶ |
| prototype | 0.5 | 0.5 | 0.5 |

The spacing of machine numbers just below 1 is ε/2, so the prototype's moduli are 1 or its nearest neighbor. The reference error is the largest distance to the 40-digit exponential of the exact product of θ and each stored entry, in units of (|θh_k|/2 + 1)ε, which allows one rounding of that product, made by both routes, plus one unit of the result:

| θ | 0.3 | 7.1 | 1000 | 10⁶ |
|---|---|---|---|---|
| base | 1.04 | 8.2 | 47 | 15.5 |
| prototype | 0.39 | 0.85 | 0.91 | 0.92 |

Both routes round the product θh_k once, so both lose |θh_k|ε/2 at large θ, and the reference allows for that. The prototype stays within the allowance at every θ; `MatrixExp` does not, and it also loses unitarity, which the entrywise route keeps because each entry is a single complex exponential of modulus one.

### 4.5 Non-diagonal controls (must not get slower, must not change)

| Control | base | prototype |
|---|---|---|
| transverse field Σ X_k, 10 qubits, `Exp` | 3.69 | 4.10 |
| random Hermitian, 8 qubits, `Exp` | 0.313 | 0.306 |
| random Hermitian, 10 qubits, `Exp` | 19.2 | 19.3 |
| diagonal plus a 1e-12 corner coupling, 10 qubits, `Exp` | 0.0292 | 0.0293 |
| symbolic 4 x 4 blocks, `Exp` | 0.0039 | 0.0041 |
| build `"RX"[θ]` | 0.0097 | 0.0100 |
| build `"R"[θ, "XX"]` | 0.0218 | 0.0215 |
| build `"RZ"[θ]` | 0.0079 | 0.0077 |
| build `"R"[θ, "ZZ"]` | 0.0179 | 0.0178 |
| build `"R"[0.3, "ZZZZ"]` | 0.0307 | 0.0301 |
| build a 4-qubit order-2 Trotter circuit, 5 steps | 0.157 | 0.154 |

The differences are the spread of runs on a loaded machine, in both directions. The transverse-field control came out slower on the prototype in this run and in the two before it, so it was measured again alternately on the two copies (`ab.wls`): the smallest of seven runs is 3.71 and 3.69 s on base and 3.72 and 3.70 s on the prototype, and the strict test on its matrix takes about a microsecond. The non-diagonal results are identical to the last bit (section 6).

## 5. The change

### 5.1 `Utilities.m`: the exponential and its diagonal test

```diff
 PackageScope["zeroBasePower"]
+PackageScope["diagonalMatrixQ"]
+PackageScope["finiteArrayQ"]
+PackageScope["matrixExponential"]
@@ before SetPrecisionNumeric, Utilities.m:546
+(* For the exponential, a matrix counts as diagonal when it is square and every
+   off-diagonal entry is zero, with no tolerance: a coupling at roundoff relative to
+   the largest entry still mixes two degenerate levels over a long time, and
+   MatrixExp resolves that mixing. A symbolic off-diagonal entry counts as zero when
+   DiagonalMatrixQ's zero test proves it zero; that test is best effort, and an
+   identically zero entry it cannot prove sends the matrix to MatrixExp. *)
+diagonalMatrixQ[mat_] := SquareMatrixQ[mat] && DiagonalMatrixQ[mat, Tolerance -> 0]
+
+(* No entry is infinite or indeterminate. A SparseArray is atomic, so its stored
+   values and its background are read. *)
+finiteArrayQ[a_SparseArray] := finiteArrayQ[{a["NonzeroValues"], a["Background"]}]
+finiteArrayQ[a_] := FreeQ[a, Indeterminate | _DirectedInfinity]
+
+(* The exponential of a diagonal matrix is the diagonal of the exponentials of its
+   entries, and its action on a vector multiplies the vector by them. Any other
+   matrix, and a diagonal with an infinite or indeterminate entry, goes to
+   MatrixExp, which fails on the latter. *)
+matrixExponential[mat_ ? diagonalMatrixQ, v___] := With[{d = Normal[Diagonal[mat]]},
+    diagonalAction[diagonalExp[d], v] /; finiteArrayQ[d]
+]
+matrixExponential[mat_, v___] := MatrixExp[mat, v]
+
+diagonalAction[values_] := DiagonalMatrix[values, TargetStructure -> "Sparse"]
+diagonalAction[values_, v_] := values v
+
+(* A numeric diagonal at machine precision, exact entries included, is exponentiated
+   in machine numbers, as MatrixExp does, and as one packed array. On a packed array
+   Exp returns an exponential below the smallest normalized machine number as a
+   subnormal number or zero, and one above the largest as an arbitrary-precision
+   number, without a message; on an unpacked list it raises General::munfl. A numeric
+   diagonal at a higher finite precision is exponentiated at that precision. *)
+diagonalExp[d_ ? machineVectorQ] := Exp[packedVector[N[d]]]
+diagonalExp[d_ ? (VectorQ[#, NumericQ] && Precision[#] < Infinity &)] := Exp[N[d, Precision[d]]]
+diagonalExp[d_] := Exp[d]
+
+machineVectorQ[d_] := VectorQ[d, NumericQ] && Precision[d] === MachinePrecision
+
+packedVector[y_ ? (FreeQ[#, _Complex] &)] := Developer`ToPackedArray[y, Real]
+packedVector[y_] := Developer`ToPackedArray[y, Complex]
```

`DiagonalMatrix[d, TargetStructure -> "Sparse"]` is the documented sparse diagonal. At 2^16 entries it builds in under 0.1 ms, against about 0.7 ms for index pairs and 30 ms for `Band` (`bench.wls`, group `detect`), and its zero is a machine zero for a machine diagonal, so `Normal` of the result stays packed. The condition inside `With` lets the finiteness test and the exponential share one reading of the diagonal; a diagonal with an infinite or indeterminate entry fails the condition and goes to `MatrixExp`.

`diagonalExp` treats the diagonal as a whole. A diagonal such as `{-800., 0}` or `{1., Pi}` mixes machine and exact entries; exponentiating it entry by entry would return exact `E^Pi` and raise `General::munfl` on `-800.`, where `MatrixExp` returns machine numbers silently. Packing the diagonal first matters too. On an unpacked list `Exp` raises `General::munfl` for an exponential below the smallest normalized machine number, including a complex entry such as `-708. + 1.5 I` whose real part is subnormal. On a packed array it returns a subnormal number or zero silently, as `MatrixExp` does (e^-720 is 2.0e-313 either way, `edge.wls`), so Boltzmann factors as small as e^-720 keep their ratios (`DiagonalExp-Gibbs-shift-invariance`). An exponential beyond the largest machine number comes back as an arbitrary-precision number: e^800 with precision 13.05, which is what significance arithmetic assigns to e^x for the machine number x = 800.

### 5.2 `QuantumOperator.m`: the call sites

On `a1799b61` the scalar-base power and `MatrixExp[qo]` each hand a map to `matrixMapOperator`, which applies it to `op["Matrix"]` at once or, for declared parameters, at each substitution. Changing the two maps covers both cases:

```diff
@@ QuantumOperator.m:585  (base ^ qo for a nonzero base, so also E ^ qo and Exp[qo])
-scalarBasePower[base_, mat_] := MatrixExp[Log[base] mat]
+scalarBasePower[base_, mat_] := matrixExponential[Log[base] mat]
@@ QuantumOperator.m:600  (MatrixExp[qo])
-QuantumOperator /: MatrixExp[qo_QuantumOperator] := matrixMapOperator[MatrixExp, qo, Exp]
+QuantumOperator /: MatrixExp[qo_QuantumOperator] := matrixMapOperator[matrixExponential, qo, Exp]
@@ QuantumOperator.m:639  (the guard on a substituted matrix)
-    ConfirmBy[Confirm[g[ConfirmBy[mat, FreeQ[#, Indeterminate | _DirectedInfinity] &]]], MatrixQ],
+    ConfirmBy[Confirm[g[ConfirmBy[mat, finiteArrayQ]]], MatrixQ],
```

`Log[base]` times a diagonal matrix is diagonal, so `base ^ qo` needs nothing else, and the zero-base branch is untouched. Labels are untouched, and `matrixMapOperator`'s check of the result accepts a `SparseArray`. The guard reads the stored values of a `SparseArray` because `FreeQ` does not look inside one; without that, `Exp` of `diag(1/a, b)` with declared parameters at `a = 0` exponentiated `ComplexInfinity` instead of failing at the guard as on base (`edge.wls` shows the `Failure` on both copies now). Main at `bdbe7165` makes the same guard read a `SparseArray`'s values through `valuelessEntriesQ`, so there this hunk is dropped.

`MatrixExp[qo, qs]` (`:689-703`) is rewritten so that it gives what `Exp[qo][qs]` gives:

```diff
+(* The exponential e^M of the operator acting on the state, as Exp[qo][qs] gives: a
+   pure state is multiplied by e^M, a mixed state rho becomes e^M rho e^(M^dagger),
+   and a superoperator acts on the density vector. An operator on any other qudits
+   than exactly the state's, with their dimensions, goes through Exp[qo][qs], which
+   extends the operator by the identity on qudits it does not act on, extends the
+   state to qudits it lacks, and fails on a dimension mismatch. *)
+QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] /; ! wholeRegisterQ[qo, qs] := Exp[qo][qs]
+
 QuantumOperator /: MatrixExp[qo_QuantumOperator, qs_QuantumState] := Enclose @ With[{op = qo["Sort"]},
     QuantumState[
-        If[ op["VectorQ"] && qs["VectorQ"],
-            MatrixExp[op["Matrix"], QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]["StateVector"]],
-            ArrayReshape[
-                MatrixExp[op["ToMatrix"]["Matrix"], QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]["DensityVector"]],
-                {#, #} & @ op["OutputDimension"]
-            ]
-        ],
+        ConfirmBy[exponentialAction[op, QuantumState[qs, QuantumBasis[op["Input"], qs["Input"]]]], ArrayQ],
@@ after it
+wholeRegisterQ[qo_, qs_] := Sort[Transpose[{qo["InputOrder"], qo["InputDimensions"]}]] === Transpose[{Range[qs["OutputQudits"]], qs["OutputDimensions"]}]
+
+exponentialAction[op_ ? (#["VectorQ"] &), qs_ ? (#["VectorQ"] &)] := matrixExponential[op["Matrix"], qs["StateVector"]]
+exponentialAction[op_ ? (#["VectorQ"] &), qs_] := With[{u = matrixExponential[op["Matrix"]]}, u . qs["DensityMatrix"] . ConjugateTranspose[u]]
+exponentialAction[op_, qs_] := ArrayReshape[matrixExponential[op["ToMatrix"]["Matrix"], qs["DensityVector"]], {#, #} & @ op["OutputDimension"]]
```

A pure state keeps WL's `MatrixExp[m, v]`, which acts on the vector. A mixed state under a vector-type operator is conjugated by u = e^M, which for a diagonal M is a sparse diagonal, so u ρ u† costs a pass over ρ; a superoperator keeps its route. The rule for an operator outside the whole register is defined first, and `DiagonalExp-exponential-on-part-of-the-register` and `DiagonalExp-exponential-beyond-the-register` show that it is the one that applies there. A superoperator on the whole register (`"InputOrder"` {1, 2} and `"InputDimensions"` {2, 2} for the two-qubit Liouvillian) keeps its route. The `ConfirmBy` turns anything but an array into a `Failure` inside the `Enclose` the rule already had.

### 5.3 `QuantumOperator.m:664`: an operator with declared parameters holds only its diagonal

```diff
-heldOperatorMatrix[_, op_] := With[{mat = Normal[op["Matrix"]]}, Hold[mat]]
+heldOperatorMatrix[_, op_] := heldMatrix[op["Matrix"]]
+
+heldMatrix[mat_ ? diagonalMatrixQ] := With[{d = Normal[Diagonal[mat]]}, Hold[DiagonalMatrix[d, TargetStructure -> "Sparse"]]]
+heldMatrix[mat_] := With[{m = Normal[mat]}, Hold[m]]
```

The list of diagonal entries stays inside the held body, where substituting the parameters reaches it (a built `SparseArray` is atomic and would not be reached, which is why base holds a `List`); the sparse diagonal is built after the substitution, and the map then sees a diagonal. The comment above `heldOperatorMatrix` says so. The test runs on the symbolic matrix when the operator is built and needs exact zeros, so an off-diagonal entry that vanishes only at some parameter values (say `a b` at `b = 0`) keeps the full matrix; a substitution that makes that matrix diagonal is still exponentiated entry by entry, since the map tests the substituted matrix. The same holding serves `f[qo]` with declared parameters. A nested map (`Exp[Sin[H]]`, the rule at `:657`) still reads the inner result back as a dense 4^n matrix, but the outer exponential then goes entry by entry: one substitution at 10 qubits takes 0.14 s against 1.25 s on base (`nested.wls`; open question 4).

### 5.4 Separate change: `QuantumState/Properties.m:1055`, keep the superoperator sparse

```diff
-        Transpose[ReshapeArray[{qs["StateTensor"]}, Join[#, #] & @ qs["MatrixNameDimensions"]], 2 <-> 3],
+        Transpose[ReshapeArray[stateTensorArray[qs["StateTensor"]], Join[#, #] & @ qs["MatrixNameDimensions"]], 2 <-> 3],
@@ after the "StateMatrix" property
+(* The state tensor as an array for ReshapeArray: the tensor of a state with no
+   qudits, such as a full trace, is a scalar. *)
+stateTensorArray[t_ /; ArrayDepth[t] == 0] := {t}
+stateTensorArray[t_] := t
```

The braces make `ReshapeArray` act on a `List` holding the stored `SparseArray`, and the result comes back as a dense `List`; without them it stays a `SparseArray` with the same entries. They exist for one case: for a state with no qudits (a full trace `QuantumPartialTrace[rho]`, or `QuantumState[{{3/10}}, QuantumBasis[1]]`) the tensor is a scalar, and without the braces `"Matrix"` came back as an unevaluated `ReshapeArray[...]`. The braces came with `0cf50215` ("fix equality and trace of mixed operators"), and the change keeps them exactly for a tensor of rank 0 (`DiagonalExp-rank-zero-matrix`). It reaches every matrix-type operator's `"Matrix"`, and also `"Matrix"` and `"StateMatrix"` of sparse mixed states (a density matrix, a partial trace), which come back as a `SparseArray` with the same entries instead of a `List`, and non-diagonal exponentials move in the last bit: the amplitude-damping exponential by 5.6e-17 (`reach.wls`). That reach is why it belongs in its own commit.

## 6. Correctness

`correct.wls` runs the same cases on base and prototype (`correct-base.out`, `correct-proto.out`). A SAME line hashes the binary serialization of the whole result, so equal hashes mean results identical to the last bit; a plain `Hash` gives machine numbers that differ only in their last bits the same value, and an earlier version of this battery, which used it, reported `"R"[0.7, "ZIZ"]` as unchanged. An ERR line is the largest entry error against the named reference.

**Identical to the last bit on base and prototype:** the exact diagonal diag(iπ, log 2, 0, -iπ/2), whose exponential is diag(-1, 2, 1, -i); the symbolic diag(a, b, c, d) under `Exp`, `x ^ qo` and `MatrixExp`, and the label of its exponential; `"RZ"[θ]`, `"RZ"[θ, 3]`, `"R"[θ, "ZZ"]`, `"R"[θ, "ZZZ"]`, `"RX"[θ]`, `"R"[θ, "XX"]` and the `"RZ"` label; every non-diagonal control (a 4-qubit transverse field, a random 16 x 16 Hermitian matrix, exact and symbolic 2 x 2 matrices, and 2 x 2 matrices with couplings 1e-12 and 1e-15); and `E ^ qo` of a non-square operator, which fails with the same `MatrixPower::matsq` message on both.

**Changed:**

| Case | base | prototype |
|---|---|---|
| `"RZ"[0.4]`, `"R"[0.7, "ZIZ"]` | `MatrixExp`'s phases | `Exp`'s phases, one or two units in the last place away |
| closed form of `Exp` with γ declared, 8 qubits | exact `0` off the diagonal | machine `0.` off the diagonal; the same diagonal |
| 8-qubit Ising e^(-0.3iH), error against `Exp` of the entries | 2.5e-16 | 0 |
| the same, against `MatrixExp` of the dense matrix | 0 | 2.5e-16 |
| the same, largest entry of U U† - 1 | 5.6e-16 | 2.2e-16 |
| off-diagonal a b at b = 0, error against `Exp` of the entries | 4.4e-16 | 0 |
| `MatrixExp[-0.3 i XI, ρ]`, error against `Exp[-0.3 i XI][ρ]` | 0.070 | 0 |
| `MatrixExp[-0.3 i ZZ, ρ]`, error against `Exp[-0.3 i ZZ][ρ]` | 0.13 | 2.8e-17 |
| trace of `MatrixExp[-0.3 i XI, ρ]` | 1.094 | 1 to 2.2e-16 |
| 8-qubit `MatrixExp[qo, ψ]`, error against `MatrixExp` of the matrix on the vector | 0 | 1.0e-16 |
| 3-qubit dephasing e^(0.7L)ρ, error against the closed form | 5.6e-17 | 3.6e-18 |
| the same, change of the trace and of the populations | 2.2e-16 and 5.6e-17 | 0 and 0 |
| `Exp[-2000 qo]` for qo = diag(1., 0.) | e^0 comes back as 0.99999999999995791 with precision 12.95 | the machine number 1. |
| `Exp` of diag(-2000. + 3i, i) | arbitrary-precision entries, precision 15.1 and 15.4 | machine numbers |
| e^0 beside e^-800 in diag(-800., 0), distance from 1 | 1.4e-13 | 0 |
| modulus of the phase e^(i 10²⁰) in diag(10²⁰ i, 1.) | 1.6 × 10⁶⁵¹ | 0.9999999999999999 |
| 30-digit diagonal (1, 2) | precision 30. | precision 29.70, which is 30 - log10 2, what significance arithmetic assigns to e^2 |
| diag(a, b) with the provable zero `Sin[y]^2 + Cos[y]^2 - 1` above the diagonal | e^a, e^b and the entry (e^a - e^b)(Sin[y]^2 + Cos[y]^2 - 1)/(a - b), which is 0/0 at a = b | diag(e^a, e^b) |

**Unchanged, and checked on both:** with γ declared, the value at γ = 0.3 equals `Exp[-0.3 i qo]` exactly and γ = 0 gives the identity; colliding parameters give e times the identity for `Exp` at a = b = 1, zero for `Sin` at 0, the identity for `Exp[Sin[...]]` at 0 and e² times the identity for `MatrixExp` at 2; the non-diagonal operator with declared parameters gives the identity at (0, 0) and exactly `MatrixExp` at (1/3, 1/5); the superoperator's `MatrixExp[qo, qs]` equals `Exp[qo][qs]` to 5e-18; the resonant matrix of section 3 goes to `MatrixExp` on both, with the same result; `Cos` of it loses the rotation on both (error 5.0e-7); the underflow cases and the 8-qubit exponential show no message on either; the zero operator's exponential is the identity; and `Normal` of a machine exponential is a packed array.

## 7. Risks

1. **Numeric results move in the last bits, towards the exact value, and exact ones can change form.** Numeric diagonal exponentials change by one or two units in the last place: `"RZ"[0.4]` and `"R"[0.7, "ZIZ"]` both do, and so does every diagonal Trotter step with numeric coefficients. They move by more where `MatrixExp` was worse: e^0 beside e^-800 moves by 1.4e-13, and the phase e^(i 10²⁰), which `MatrixExp` returns with modulus 1.6 × 10⁶⁵¹, gets modulus one (section 6). An exact result can be written differently with the same value: `Exp[I t QuantumOperator["PauliZ"[3]]]` has the entry `E^(-((-1)^(5/6) t))` on base and `E^((I t)/E^((2 I Pi)/3))` on the prototype. The suite does not notice; user code that compares such results with stored values would.
2. **`MatrixExp[qo, qs]` changes value** wherever it differed from `Exp[qo][qs]` (section 2): on a mixed state, from the non-physical exp(vec(M) vec(M)†) vec ρ to e^M ρ e^(M†); on part of a register, from an invalid state to (1 ⊗ e^M) ρ (1 ⊗ e^M)†; for an operator on a qudit the state lacks, from a state on the state's own qudits, with the operator's order ignored, to a state extended to the operator's qudits; and for a dimension mismatch, from an invalid state to `$Failed` with `QuantumCircuitOperator::dim`. Code that relied on the old values was relying on wrong ones, but they change.
3. **A tolerance would drop physics.** With the default tolerance, section 3's resonant coupling is lost. The strict test avoids this, at the price that a matrix diagonal only up to roundoff goes to `MatrixExp`, as today.
4. **Symbolic zero tests.** `DiagonalMatrixQ` applies its zero test to symbolic off-diagonal entries. It stops at the first entry it can show is nonzero, so a symbolic non-diagonal matrix costs next to nothing; the slow case is a symbolic diagonal whose off-diagonal entries are all zeros that must be proved, 0.7 ms for 4032 of them at 64 x 64 (section 3), where `MatrixExp` would be slower still.
5. **Dense matrix-type operators.** Without section 5.4, a diagonal superoperator is scanned in full as a dense `List` (0.26 s at n = 6; at n = 7 the list takes gigabytes before any test runs).
6. **Precision labels move beyond the machine range.** An exponential above the largest machine number now carries the precision significance arithmetic gives it (13.05 digits for e^800); base's `MatrixExp` labels its value with 15.95 digits although it is 2.9e-13 off (`edge.out`). Code that tests `Precision` of such results sees the change.
7. **Structured arrays are deliberately not used for the result.** In WL 15.0.1, `MatrixExp` of a structured `DiagonalMatrix` with a complex entry that underflows (`-1000. + 3. I`) raises `General::munfl` and returns unevaluated, and `Eigensystem` of a structured `DiagonalMatrix` returns a malformed `IdentityMatrix[WorkingPrecision -> MachinePrecision, List, {3, 2, 1}]` (`structured.wls`). The prototype returns a `SparseArray`, which every downstream QF path already handles.
8. **`MatrixExp[qo, qs]` on a pure state** changes route too, so its results move in the last bit as in risk 1; its gain is smaller than `Exp`'s (section 4.1).

## 8. Follow-ons (measured on base, not prototyped)

`followon.wls` on base, and `accuracy.wls` on both copies.

### 8.1 Non-integer powers

`qo ^ x` for a scalar x goes through the generic numeric-function rule (`QuantumOperator.m:597`) to `MatrixPower` on the full matrix (`Utilities.m:316`). For a diagonal Ising operator, `qo ^ 0.5` takes 0.08 s at 10 qubits and 3.8 s at 12, while the principal branch entry by entry, built as a sparse diagonal, takes 0.1 ms and 0.5 ms and agrees with it to 4.4e-16, negative entries included. Integer powers are already fast (`qo ^ 3`: 3 ms). A diagonal branch in `matrixFunction[Power, ...]` would give the same gain as the exponential.

### 8.2 The spectrum

`qo["Eigenvalues"]` and `qo["Eigensystem"]` (`QuantumOperator/Properties.m:506-510`) call QF's `eigensystem` (`Utilities.m:163-205`) on the sparse `"MatrixRepresentation"`, and `Eigensystem` converts it to a dense matrix, printing `Eigensystem::arh` at 8 and 10 qubits (at 12 the message had already been suppressed by `General::stop` in the same run). The first call on a fresh operator:

| n | `"Eigenvalues"` | `"Eigensystem"` | `Eigenvalues` of a structured `DiagonalMatrix` |
|---|---|---|---|
| 8 | 0.067 | 0.016 | 0.0002 |
| 10 | 0.19 | 0.17 | 0.0005 |
| 12 | 4.8 | 5.1 | 0.002 |

QF returns `Eigensystem`'s order, decreasing absolute value with the later index first among equal eigenvalues, at all three sizes, where the Ising energies come in equal pairs, since flipping every spin leaves E_x unchanged. WL's `Eigensystem` of the exact diagonal (2, 1, 2) likewise gives the eigenvalues (2, 2, 1) with eigenvectors e_3, e_1, e_2. A diagonal branch must reproduce that order or document a change (open question 3), and it cannot use the structured head, whose `Eigensystem` is malformed (risk 7).

### 8.3 `QuantumEvolve`

`QuantumEvolve` does not exponentiate. For a time-independent diagonal H and a numeric time range it integrates with `NDSolveValue`:

| n | `QuantumEvolve` over t in [0, 1] | error at t = 1 against e^(-iht)ψ(0) |
|---|---|---|
| 6 | 0.018 | 5.5e-7 |
| 8 | 0.053 | 1.3e-6 |
| 10 | 0.61 | 3.2e-6 |

The error is `NDSolve`'s default accuracy; e^(-iht)ψ(0) is the entrywise exponential acting on the state, which `MatrixExp[qo, qs]` on the prototype gives to machine precision in milliseconds (section 4.1). With symbolic t, `QuantumEvolve` already returns the closed form: for h = (1, -1, 2, 0) on the product state `++` it gives (e^(-it), e^(it), e^(-2it), 1)/2. A closed-form route for a constant diagonal H would change what `QuantumEvolve` returns for a numeric range, an interpolating function today, so it is a design question rather than a speed-up.

### 8.4 The QAOA mixer

The mixer e^(-iβΣX_k) is not diagonal and keeps `MatrixExp`: 0.58 s at 9 qubits on base and on the prototype, against 3 ms for the Kronecker product of nine 2 x 2 exponentials, which agrees with it to 6.9e-15. A sum of commuting terms on different qudits could be exponentiated factor by factor; that is a separate change.

### 8.5 The `"Diagonal"` constructor's label

`QuantumOperator["Diagonal"[d], order]` (`NamedOperators.m:243-252`) sets the label `OverHat[d]` with the whole list d, and the label arithmetic that `Exp` performs dominates at 16 qubits: on the prototype, `Exp[-I 0.3 qo]` takes 0.20 s with that label and 0.015 s with `"Label" -> None`. A short label would remove the cost (open question 7).

### 8.6 The `f[qo]` diagonal branch

`scalarMatrixFunction`'s diagonal branch (`Utilities.m:331-339`) builds its result with `SparseArray[Band[{1, 1}] -> values, dims]`, whose background is an exact 0, so `Normal` of a machine result is not packed; that is the defect `/wl-verify` round 3 found in the exponential's first form. `Sin` of the 16-qubit diagonal Ising operator: 39 ms, against 15 ms for the prototype's `Exp` of the same operator (`accuracy.wls`), and its matrix has an exact zero background and an unpacked `Normal` (`followon.wls`). Building with `DiagonalMatrix[values, TargetStructure -> "Sparse"]`, as `diagonalAction` does, would give a machine zero and a packed `Normal` and build in under 0.1 ms instead of 30 ms at 2^16 entries (section 5.1). It changes `f[qo]`'s output form, so it belongs with decision 1.

## 9. Tests

`QuantumOperatorDiagonalExp.wlt`, 38 tests, in the style of `QuantumOperatorMatrixFunction.wlt`. Its header states the physics brief and the regimes; the brief was written after the prototype and before these tests, not before the code, and it says which rows pass on main too. Every fixture is built with QF's constructors (`"ZZ"`, `"Z"`, Pauli strings, `"Diagonal"`, `"Liouvillian"`, `"JX"`, `"R"`, `"RZ"`, `"RX"`), and the references are closed forms, `MatrixExp` of an explicit matrix, `QuantumEvolve`'s numerical integration, the 40-digit exponential of the stored entries, or another QF route. Every numeric comparison is exact (`Max[Abs[...]] == 0`, or a difference compared with 0) or a stated bound, since `==` and `===` between machine numbers absorb differences in the last bits. On current main 13 of the 38 fail; with the prototype all 38 pass.

| TestID, after `DiagonalExp-` | Regime | What it asserts | On main at `a1799b61` |
|---|---|---|---|
| `Ising-symbolic-closed-form` | general symbolic | n = 2..5 with symbolic couplings, fields and θ: the generator is diagonal with no tolerance, and e^(-iθH) is the diagonal of e^(-iθE_x) | passes |
| `Ising-symbolic-unitary-group-law` | general symbolic | 3 qubits, symbolic: unitarity and U(θ1)U(θ2) = U(θ1 + θ2) | passes |
| `exponential-on-a-pure-register` | general symbolic | e^(-iθH) on a symbolic 3-qubit state through `MatrixExp[qo, qs]`: each amplitude gets e^(-iθE_x), and the norm stays | passes |
| `Ising-QAOA-gate-factorization` | general symbolic | n = 2..4, symbolic: e^(-iθH) equals the circuit of R_ZZ(2θJ_ij) and R_Z(2θb_k) gates | passes |
| `long-time-unitarity` | limiting | 8-qubit Ising at θ = 0.3, 1000, 10⁶: every diagonal modulus within ε of 1, nothing off the diagonal | fails: moduli off by 1.5ε to 2.4 × 10⁶ ε |
| `40-digit-reference` | rounding test | within the section 4.4 allowance at θ = 0.3, 7.1, 1000, 10⁶ | fails at all four |
| `16-qubit-Ising` | rounding test, scaling | 16-qubit Ising within the allowance, in 30 s | fails: time limit |
| `parameter-scan-12-qubit` | general symbolic, scaling | γ declared: the closed form, the identity at γ = 0, 25 values within the allowance, in 60 s | fails: time limit |
| `dephasing-symbolic` | exactly solvable | 2-qubit dephasing with symbolic γ, t and ρ: the superoperator is diagonal with no tolerance; each coherence gets e^(-i(E_x - E_y)t - 2γ d(x, y) t) with d the Hamming distance; populations, trace, and `Commutator[H, Z_k] = 0` on the operators | passes |
| `dephasing-limits` | limiting | t -> ∞ gives diag(ρ(0)); γ = 0 gives U ρ(0) U† | passes |
| `dephasing-complete-positivity` | exactly solvable | the pure-dephasing multiplier is {{1, q}, {q, 1}} ⊗ {{1, q}, {q, 1}}, q = e^(-2γt), with eigenvalues 1 ± q ≥ 0, so the map is completely positive | passes |
| `collective-dephasing` | exactly solvable | dephasing by Z_1 + Z_2: the superoperator is diagonal, coherences decay at γ(s_x - s_y)²/2, and the coherence between 01 and 10 is decoherence-free | passes |
| `superoperator-acting-on-state` | exactly solvable | `MatrixExp[t L, ρ]` equals `Exp[t L][ρ]` | passes |
| `dephasing-against-QuantumEvolve` | numerical reference | `Exp[t L]` against `QuantumEvolve`'s integration of the Lindblad equation (`NDSolve`), within 10⁻⁶ | passes |
| `exponential-acting-on-states` | invariant | `MatrixExp[qo, ρ]` equals e^M ρ e^(M†) for M = -iθZ, -iθX, -iθY, the no-jump generator diag(0, -g/2), whose trace is a + (1 - a) e^(-g), and the entangling -iθZZ on a generic two-qubit ρ; e^M ψ on a pure state | fails: the trace grows as e^(θ²) |
| `exponential-on-part-of-the-register` | invariant | -iθX or -iθZ on qubit 2 of a two-qubit state acts as 1 ⊗ e^M, on a mixed and on a pure state; dephasing qubit 1 of a Bell pair leaves the coherence e^(-2γt)/2, from the mixed and the pure state | fails: unevaluated `ArrayReshape[MatrixExp[...], ...]`, with `MatrixExp::lslc` and `ArrayReshape::listrp` |
| `exponential-beyond-the-register` | edge | X on qubit 3 of the one-qubit state `0` gives the three-qubit state (cos θ, -i sin θ, 0, ..., 0), and a qubit operator on a qutrit gives `$Failed` with `QuantumCircuitOperator::dim`, as `Exp[qo][qs]` does | fails: a one-qubit state, and an invalid state after `MatrixExp::lslc` |
| `ZZ-entangler-concurrence` | exactly solvable | the concurrence of e^(-igZZ) applied to `++` is the absolute value of sin 2g | passes |
| `spin-one-JX-eigenbasis` | exactly solvable | the spin-1 J_x is stored diagonal, and its exponential read in the computational basis is `MatrixExp` of the 3 x 3 J_x | passes |
| `Ising-chain-partition-function` | exactly solvable | Tr e^(βJ Σ Z_k Z_(k+1)) = 2 (2 cosh βJ)^(n - 1) for symbolic βJ, n = 2..5 | passes |
| `Ising-chain-correlations` | exactly solvable | ⟨Z_1 Z_r⟩ = tanh(βJ)^(r - 1) in the Gibbs state of the open 5-qubit chain, symbolic βJ | passes |
| `Gibbs-ground-states` | limiting, edge | βJ -> ∞ leaves weight 1/2 on each of 0000 and 1111 (a symbolic `Limit`); at βJ = 250 the same weights come through e^750, beyond the machine range | passes |
| `Gibbs-shift-invariance` | edge | the Gibbs weights at energies 720 and 721, from subnormal Boltzmann factors, match those at 0 and 1 within 10⁻⁹ | passes |
| `QAOA-depth-one-ring` | general symbolic | depth-one QAOA on the 4-ring: the edge term ⟨(1 - Z_1 Z_2)/2⟩ is (2 + sin 4β sin 2γ)/4 for symbolic γ and β | passes |
| `QAOA-scan-with-gamma-declared` | numerical reference | the same edge term from the cost layer with γ declared, at 25 values of γ and β = π/8: (2 + sin 2γ)/4, and 3/4 at γ = π/4 | passes |
| `symbolic-and-exact` | general symbolic | diag(a, b, c, d) and diag(iπ, log 2, 0, -iπ/2) keep their closed forms | passes |
| `hidden-zero-no-pole` | edge | a provable zero off the diagonal leaves no pole at a = b | fails: `Indeterminate`, with `Power::infy` and `Infinity::indet` |
| `hidden-zero-unproved-falls-back` | edge | `Gamma[g + 1] - g Gamma[g]` is not proved zero, goes to `MatrixExp`, and the result is still correct | passes |
| `underflow` | edge | e^-2000 is 0 and e^0 is exactly the machine number 1 | fails: e^0 has precision 12.95 |
| `overflow` | edge | e^800 within 10⁻¹⁴ of its value at 40 digits and e^-800 zero, silently | fails: e^800 off by 2.9e-13 relative |
| `complex-subnormal-part` | edge | e^(-708 + 1.5i) within ε of its modulus, silently, and e^0 exactly 1 | fails both |
| `thirty-digit-diagonal` | edge | a 30-digit diagonal with exact entries is exponentiated at that precision, equal to the exact e^(-ix) within the precision it carries | passes |
| `non-finite-diagonal-fails` | edge | `Indeterminate` and `-Infinity` on the diagonal give a `Failure` | passes |
| `mixed-exact-and-machine` | edge | machine diagonals with exact entries give machine results within ε | fails on the two with -800. |
| `substitution-decides-route` | edge | an off-diagonal a b that vanishes at b = 0 is exponentiated entry by entry there and by `MatrixExp` elsewhere | fails: `MatrixExp` is 4.4e-16 off at b = 0 |
| `resonant-roundoff-coupling-kept` | edge | section 3's matrix keeps its 10⁻³ rotation | passes |
| `parameter-collision` | edge | diag(a, b) at a = b = 1 gives e times the identity | passes |
| `rank-zero-matrix` | edge | a full trace keeps its 1 x 1 `"Matrix"` (guards section 5.4) | passes |

Full suite, `wolframscript -file Tests/RunTests.wls` on the scratch copies of `a1799b61`: the prototype passes 2892 of 2892, main's 2854 tests plus the 38 new ones, and so does the prototype with section 5.4.

Plan for landing it, if approved: commit 1 is sections 5.1-5.3 without the `MatrixExp[qo, qs]` rewiring, plus the new test file minus the three tests that main's `MatrixExp[qo, qs]` fails (`DiagonalExp-exponential-acting-on-states`, `DiagonalExp-exponential-on-part-of-the-register` and `DiagonalExp-exponential-beyond-the-register`); commit 2 is the `MatrixExp[qo, qs]` rewrite and those tests; commit 3 is section 5.4 plus a test that the Liouvillian's `"Matrix"` is a `SparseArray` and that building and exponentiating the n = 6 dephasing Liouvillian finishes within a limit base misses. Run the full suite serially after each. Diagonal branches for `qo ^ x` and `"Eigenvalues"` go in separate commits, the latter with a test of the eigenvalue order.

## 10. Open questions

1. Two tests of "diagonal": strict for the exponential, the default tolerance for `f[qo]`. Accept the two, or change `f[qo]`'s roundoff rule first (question 5)?
2. Land section 5.4 with this change or on its own? It is the larger win for open systems, and it reaches every matrix-type operator and every sparse mixed state.
3. Should a diagonal `"Eigenvalues"` keep `Eigensystem`'s order (decreasing absolute value, later index first among equal values) or return the diagonal order and document it?
4. A nested map of a diagonal operator with declared parameters (`Exp[Sin[H]]`) still reads the inner result back as a dense matrix. Worth keeping a diagonal inner result diagonal?
5. The `f[qo]` routes set eigenvalues within 10 machine epsilons of zero, relative to the largest, to zero. That protects `Sqrt` and `Log` but drops a resolved small eigenvalue for an analytic f (`Cos` of section 3's matrix loses the 10⁻³ rotation). Should the setting to zero apply only to functions singular at zero?
6. Building the dephasing Liouvillian is the remaining cost after section 5.4 (9.6 s at n = 7, against 0.013 s for its exponential).
7. QF's `"Diagonal"` constructor labels an operator with its whole diagonal (section 8.5). A short label would remove a cost that, at 16 qubits, exceeds the exponential's.

## 11. Verification history

**`/wl-verify`, three rounds on the `81df7633` version; it did not converge.** Round 1 found the mixed exact and machine diagonal defect (`{{1., 0}, {0, 2}}` gave exact `E^2`; `{-800., 0}` raised `General::munfl`). Round 2 found that exact constants such as `Pi` beside machine numbers still gave exact results; that the first form of section 5.4 broke `"Matrix"` for rank-0 matrix-type objects; that the diagonal body bypassed the `Indeterminate` guard; that sharing the strict test slowed `f[qo]`; that two tests used `===`, which absorbs one unit in the last place, and one failed on base for a different reason than stated; and several wrong statements in the report. Round 3, the last the loop allows, reported 11 open issues. Three were changes in behavior, each with a fix the verifier tested in a kernel: new `General::munfl` messages for a complex entry with a subnormal real or imaginary part (hence `packedVector`), `Indeterminate` and infinite diagonals returning an operator where base fails (hence `finiteArrayQ` in `diagonalExponential`), and an exact zero background that left `Normal` of the result unpacked (hence `DiagonalMatrix[..., TargetStructure -> "Sparse"]`). The rest concerned the report and the tests. All are addressed here, after the last verification round, so no fresh verifier has checked these changes.

**`/wl-quality` round 1 on the `81df7633` version: 17 open issues.** What changed for each:

- A1, a hand-built sparse diagonal: now `DiagonalMatrix[d, TargetStructure -> "Sparse"]` (the critic's `Band` form builds some 40 times slower than index pairs at 2^16, section 5.1).
- A2, fixtures that bypassed QF's constructors: the tests and scripts build with `"ZZ"`, `"Z"`, Pauli strings, `"Diagonal"`, `"Liouvillian"` and `"JX"`.
- A3, the hidden-zero behavior rests on a best-effort zero test: stated in the kernel comment and section 3, with `hidden-zero-unproved-falls-back` added.
- B1, `If` dispatch at four kernel sites: pattern-dispatched rules (`diagonalExp`, `packedVector`, `heldMatrix`, `finiteArrayQ`); the last inline `If`, in section 5.4, went in round 2.
- B2, a closed form built by walking indices: `Outer[Subtract, ...]` and `Outer[HammingDistance, ...]`.
- B3, imperative patterns in the scripts: `Block[{$MessageList = {}}, ...]`, `run[group]` definitions, and one Ising energy function everywhere.
- C1, a stale comment and duplicated checks: comment rewritten, one `finiteArrayQ`, `Normal[Diagonal[mat]]` computed once.
- C2, tests that did not read as physics: `rho0`, `rhoT`, and expected values derived from their inputs.
- D1, `MatrixExp[qo, qs]` non-physical on a mixed state: fixed through the generator M ⊗ 1 + 1 ⊗ M* and tested.
- D2, invariants never asserted: unitarity, the group law and [H, Z_k] = 0 asserted, and unitarity at θ = 10⁶.
- E1, only base cases: the concurrence of e^(-igZZ) on `++`, the spin-1 J_x in its eigenbasis, and collective dephasing.
- F1 and F2, numeric where symbolic: symbolic dephasing in γ, t and ρ; a symbolic 3-qubit Ising with the group law; the symbolic closed form with γ declared.
- G1, no brief before the code: a brief is now in the test header, written after the prototype and before the tests. It cannot be moved before the code, and the report says so.
- G2, no asymptotic regime: t -> ∞ gives the pinching map, γ -> 0 the unitary, and θ reaches 10⁶.
- G3, a reference that was the implementation itself: the 40-digit reference with the (|θh_k|/2 + 1)ε allowance; the accuracy claim rests on it and on unitarity.
- G4, missing edge cases: overflow, the complex subnormal part, non-finite entries, the substitution that decides the route, and `MatrixExp[qo, qs]`.

Its out-of-scope note, the QAOA mixer, is section 8.4. Main then moved and the machine rebooted; the prototype was rebuilt on `a1799b61` from `apply_proto.py`, and every number in this report was measured there.

**`/wl-quality` round 2 on the `a1799b61` version: 11 open issues.** What changed for each:

- A1, the underflow mask rested on a claim the kernel refutes, since `Exp` on a packed array is silent at underflow and overflow, and it zeroed subnormal Boltzmann factors, so e^(-H) for energies 720 and 721 had trace 0: the mask is gone, `DiagonalExp-Gibbs-shift-invariance` holds the Gibbs weights at those energies to the weights at 0 and 1, and the overflow test compares at 40 digits.
- A2, a hand-written commutator: `Commutator[x, y, {Dot, 4}]`, documented for square matrices.
- B1, section 5.4 tested the container type: the braces now depend on the rank of the tensor (`stateTensorArray`), which is the reason they exist.
- C1, a second layer of diagonal helpers: `diagonalExponential` and `maskedExp` are gone, and one `With` with a condition reads the diagonal once. The critic also asked that `f[qo]`'s diagonal branch share one detection rule and one sparse builder with the exponential; that is decision 1, with D2.
- C2, the d² by d² generator M ⊗ 1 + 1 ⊗ M* for a mixed state, 2d³ stored entries for a dense M: replaced by u ρ u† with u = e^M.
- D1, `MatrixExp[qo, qs]` returned an invalid state for an operator on part of the register, on main too: such an operator now goes through `Exp[qo][qs]`, tested by `DiagonalExp-exponential-on-part-of-the-register`.
- D2, Euler's formula fails by 10⁻³ on the resonant matrix, identically on main: stated in decision 1 with the single rule the critic proposed. It is not prototyped, because `f[qo]`'s Schur route would still set the ±10⁻³ eigenvalues to zero.
- E1, the mixed-state test could not tell e^(M†) from e^(Mᵀ) or e^(M*): it now covers Y and the non-Hermitian no-jump generator, whose trace a + (1 - a) e^(-g) it asserts.
- E2, no Gibbs physics in the real-exponent regime: the chain's partition function for n = 2..5, its two ground states at βJ = 250 through the overflow, and the shift invariance with subnormal factors.
- G1, the physics rows pass on main too: the brief now says so, and the rows meant to cover the entry-by-entry route assert that the generator is diagonal with no tolerance.
- G2, the general symbolic row fixed n = 3: the closed form for n = 2..5, and the factorization into the QAOA cost-layer gates R_ZZ(2θJ_ij) and R_Z(2θb_k) for n = 2..4.


**`/wl-quality` round 3, the last the loop allows, on the version after round 2: 13 open issues, so the loop did not converge.** Everything below was done after it and has not been reviewed by a fresh critic:

- A1, [H, Z_k] was tested on matrices with a hard-coded dimension: QF's own `Commutator` on the operators.
- B1, tests took `Normal` of 4096 by 4096 matrices: they compare `SparseArray`s (`IdentityMatrix[4096, SparseArray]`, `"NonzeroPositions"`).
- B2, `MapThread` over arithmetic that is already listable: `Clip`.
- C1, `f[qo]`'s diagonal branch keeps its own test of "diagonal" and its own sparse builder, whose exact zero background leaves `Normal` unpacked: not changed, since it changes `f[qo]`'s output; section 8.6 measures it, and it belongs with decision 1.
- C2, the test Liouvillian was spelled out four times: `deLiouvillian`.
- D1, complete positivity of the dephasing map was never asserted: `DiagonalExp-dephasing-complete-positivity`.
- D2, `MatrixExp[qo, qs]` extended the state for an operator on a qudit the state lacks, which risk 2 did not list, and returned an invalid state on a dimension mismatch: the whole-register route now requires the state's qudits and dimensions exactly, everything else goes through `Exp[qo][qs]`, a `ConfirmBy` turns a non-array into a `Failure`, and `DiagonalExp-exponential-beyond-the-register` and risk 2 cover both cases. This is the one kernel change after round 3.
- D3, the higher-precision branch of `diagonalExp` and the whole-register action on a pure state were reached by no test: `DiagonalExp-thirty-digit-diagonal` and `DiagonalExp-exponential-on-a-pure-register`.
- E1, the QAOA use case was tested as the cost layer alone, whose expectation values do not depend on γ: `DiagonalExp-QAOA-depth-one-ring` with symbolic γ and β, and `DiagonalExp-QAOA-scan-with-gamma-declared`.
- E2, the `MatrixExp[qo, qs]` fixes were tested on one qubit only: the entangling -iθZZ on a generic two-qubit ρ, and one half of a Bell pair dephased, from the mixed and from the pure state.
- F1, the β → ∞ limit was a number: the symbolic `Limit`, with βJ = 250 kept as the overflow case, and the correlations tanh(βJ)^(r - 1) of the open chain.
- G1, the brief was written after the code, so every symbolic row passed on main: the rows the critic named as missing (D1, D2, E1, E2) are added, and the brief still dates from after the code.
- G2, the numerical reference was not independent of the stored entries: `DiagonalExp-dephasing-against-QuantumEvolve` integrates the Lindblad equation with `NDSolve`, and the 40-digit test is now called a rounding test of the stored entries.

## 12. Files in `exp-diagonal-shortcut-prototype/`

- `proto-Utilities.diff`, `proto-QuantumOperator.diff`: sections 5.1-5.3. `protos-StateMatrix.diff`: section 5.4.
- `apply_proto.py`: applies the prototype to a paclet exported from `a1799b61` (`git archive a1799b61 QuantumFramework Tests | tar -x -C <dir>`, then `python3 apply_proto.py <dir>`, adding `sparse-state-matrix` for section 5.4), asserting that each edited line matches exactly once.
- `QuantumOperatorDiagonalExp.wlt`: the 38 tests (section 9); `file-base.out` and `file-proto.out` are its runs on main and on the prototype.
- `bench.wls` and `bench.out`: sections 3 (detection costs) and 4.1-4.3, 4.5; `ab.wls` and `ab.out`: the alternate runs of section 4.5's transverse-field control.
- `accuracy.wls`, `accuracy-base.out`, `accuracy-proto.out`: sections 4.4, 8.4 and 8.5.
- `correct.wls`, `correct-base.out`, `correct-proto.out`: section 6.
- `followon.wls` and `followon-base.out`: sections 8.1-8.3.
- `tolerance.wls` and `tolerance.out`, `resonant.wls` and `resonant.out`, `liouvillian-matrix.wls` and `liouvillian-matrix.out`: section 3. `reach.wls` and `reach.out`: section 5.4. `structured.wls` and `structured.out`: risk 7. `nested.wls` and `nested.out`: section 5.3. `edge.wls` and `edge.out`: the exact forms in risk 1, `MatrixExp`'s subnormal result in section 5.1, the e^800 row of section 9, and the guard of section 5.2.
- `suite-proto.out`, `suite-protos.out`: the full suites (section 9); `file-protos.out`: the new tests with section 5.4; `file-bdbe7165.out`: the new tests on main at `bdbe7165` with every edit but the guard. `run-all.sh`: the chain that produced every output here, one kernel at a time.
- `earlier/`: the `2ff7d0ce` and `81df7633` versions of the prototype, with their scripts and outputs.
