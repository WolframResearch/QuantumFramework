# Independent WL quality review, round 1 — applied permanent selector

Scope: `/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/QuantumFramework/Kernel/SecondQuantization/Interferometer.m` and `/Users/mohammadb/.codex/worktrees/56d8/QuantumFramework/Tests/QuantumInterferometer.wlt`. No production files edited. Complete target files and complete quality critic mandate read. Review concerns the applied selector and its regression tests; it does not recertify unrelated Gaussian/coherent algorithms.

Anchor: repository HEAD `b92953cd8b5913fa76db7a1668533fa038a4d0d5`; current manifest version 2.1.3. Kernel SHA256 `cdd11b786204bbfe98ff963e05d34a6998005cbe0c37ae8e52903edd08ccbe6f`; test SHA256 `5c9eb9145dbbcdf513d2826229b5cd588e6ce1715c5a47d01ba7915bd1d4c22a`. External QF audit anchor is older (`f77a6640`); live source, not stale audit line numbers, is the contract here.

Brief: passive bosonic scattering with W[j,k] = <1_j|U|1_k>, integer nonnegative occupations S,T of m modes, and normalized amplitude Perm(W[T,S])/sqrt(product S! product T!). Change only which existing numerical permanent route is chosen. Preserve exact and symbolic expressions, vacuum, number conservation, HOM interference, orientation, factorial normalization, and numerical accuracy. Calibration weights are empirical and architecture/version dependent.

Regimes recovered from the brief and executable artifacts:

| Regime | Evidence |
|---|---|
| General symbolic | Unrestricted symbolic 2x2 W with three photons gives b c^2 + 2 a c d exactly; existing symbolic-photon tests remain unchanged (tests 328–387). |
| Solvable special | Exact HOM output {0,1/2,1/2}; exact Fourier amplitude; decoupled phase networks. |
| Large-size / decoupling limit | 12-photon output collision, 16-photon strong-collision route, 24-photon phase bunching; 128-mode cost boundary. No physical approximation or asymptotic formula is introduced by this dispatch change. |
| Numerical reference | 1741 exact four-mode occupation pairs at machine/35-digit precision; whole output-sector norm; sparse/mixed arithmetic; independently retained 85-digit references and rationalized relative errors. |
| Failure / edge | Vacuum, unequal total photons, exact/symbolic guard, no-collision K=2^n boundary; pre-existing wrong lengths/nonunitary/negative patterns remain covered. |

A — PASS. Documented/current symbols and signatures.
- New predicate uses core WL symbols only. Source 629–637 defines exactly `(w,s,t)` and returns a Boolean on the validated public domain; source 535–540 and 606–610 establish that integer patterns have equal mode length and nonnegative entries before reaching this private helper.
- Public `QuantumInterferometer["Fourier"[m]]` constructor is defined at source 464–465, and `qi["Amplitude",s->t]` at 535–540. The test aliases at 611–612 bind actual PackagePrivate symbols, consistent with existing hafnian alias 533. No invented public option or return shape is used.
- Live WL 15.0 ARM probe: exact symbolic W has Precision Infinity; mixed 50-digit/machine W has MachinePrecision. No messages. Separate actual-package suite records all 60 tests successful and post-run actual-dispatch checks {True,True} in `/tmp/qi-applied-verifier-suite.log`.
- Defects: none. Minimal fix: none.

B — PASS. Idiomatic WL.
- Source 629–637 uses nested lexical `With` for constants, `And` short-circuit guards, listable `s+1`, `Times@@`, `Min`, and `Count`. No imperative loop, mutable cache, dense conversion, or redundant numeric mapping was added.
- `p =!= Infinity` precedes every cost computation. `k < 2^n` precedes row counting and score calculation. Live Trace with exact symbolic W and occupations {1000000,0} found no `Total`, `Min`, `Count`, or `2^n` evaluation: `{}`.
- New tests 708–724 use `Table` to enumerate a genuine Cartesian grid; `MapIndexed` constructs row/column lists directly. This is appropriate structured test data, not list mutation. Existing compiled `Do` loops are unchanged and outside the selector change.
- Defects: none. Minimal fix: none.

C — PASS. Terse physics reasoning.
- Source 619–622 retains the normalized bosonic permanent/grouped Ryser argument. New helper separates the calibrated routing decision from the existing amplitude formula. Source 625–628 explains whole-matrix promotion, occupied-row dependence, weights and WL 15.0 ARM calibration without asserting universal optimality.
- Source-to-HEAD byte checks: `repeatedPermanent` is unchanged; `permanentAmplitude` differs only by replacement of its predicate. Precision promotion, transpose orientation and factorial normalization are byte-identical.
- In a<=b ties, source 632 counts occupied t rows, matching the original s-oriented branch at 644. In a>b it counts s, matching the transposed branch. 17,361 live exact/machine/35-digit small-occupation selector cases agreed with the literal approved plan; separate verifier recorded 388/388 benchmark policy matches.
- Defects: none. Minimal fix: none.

D — PASS. Applicable invariants.
- New production tests 627–645 assert vacuum, total-photon conservation and exact HOM. Tests 694–701 check independent Sqrt[n!] product normalization for output collision. Tests 730–743 check whole-sector normalization and sparse diagonal unitary phases.
- Independent source-fragment live results: vacuum/mismatch/HOM `{1,0,0,1/2}`; unrestricted symbolic three-photon formula True; four-photon machine sector norm error `1.3322676295501878e-15`; 35-digit sector norm residual `0``34.82362025148232`; both collision orientations residual `0``37.473736691237384`; 35-digit sparse repeated-route phase residual `{True,0``33.59088569519845}`; original W remains structurally unchanged; 24-photon phase error `7.108895957933346e-16`.
- High-precision residuals above are significance zeros, not proofs of exact equality. Hermiticity/positivity/trace tests are inapplicable to this passive transition-amplitude selector; no Hamiltonian/density operator is added or changed.
- Defects: none. Minimal fix: none.

E — PASS. Nontrivial structure.
- Artifact does not stop at vacuum, one photon or HOM. Tests exercise collisions on input/output, large photon counts 12/16/24, large mode count 128, arbitrary precision, sparse/mixed matrices, and all 1741 four-mode input/output pairs for n=1..4 at two precisions. Exact independent expectations include the permanent formula with explicit occupation expansion, HOM, diagonal phases, and Sqrt[n!] row-product identity.
- Static checks confirm the entire original 47-test prefix is byte-for-byte unchanged, 60 unique TestIDs exist, and exactly 13 IDs were appended. New tests contain no `Get[`, audit-directory reference, `InstallProposal`, `ProposedAmplitude`, or `SourceHelpers`. Test discovery remains the existing Tests/*.wlt runner. No new wall-time assertion was added.
- Defects: none. Minimal fix: none.

F — PASS. Symbolic-first derivation.
- Exact/symbolic Precision Infinity is rejected by selector before costs and continues to built-in exact Permanent at source 646. Numeric N/SetPrecision occurs only within numerical routing or test comparisons. Empirical weight 1/10 is exact, avoiding accidental score-rounding boundaries.
- Tests 614–645 retain exact and symbolic W and exact HOM; 707–717 derive the four-mode reference exactly before machine/35-digit evaluation. Independent symbolic probe establishes the generic three-photon polynomial b c^2+2 a c d without substituting a special angle or numeric matrix.
- Defects: none. Minimal fix: none.

G — PASS. Regime coverage and limits.
- The regime table above has executable coverage for every row. Current probe adds arbitrary-precision whole-sector normalization and sparse repeated-route checks beyond the integrated suite. Separate actual-kernel 60/60 run confirms old symbolic, Fock-circuit, sampling, squeezing and edge regressions remain successful.
- Limitation to report, not an implementation defect under the approved rule: calibration does not guarantee the fastest route for every representation. Verifier found the 128-mode, 12-photon balanced case about 20.4 ms new vs 12.0 ms old on an unpacked machine Fourier matrix; explicit Complex packing changed this to about 1.33 ms new vs 12.39 ms old (`/tmp/qi-applied-verifier-storage.json`, `...-storage-complex.json`). The source honestly calls these calibrated work scores and names WL 15.0 ARM. Do not claim universal speedup or representation-independent optimum.
- Environment limitation: `/tmp/qi-applied-verifier-suite.log` and `/tmp/qi-applied-tests.log` both report `TOTAL 60 60`, `ALERT LoadMessages 6`, `STATUS RED` from existing offline ImageMetadataTools/Memoize load messages. This is a successful 60-test result in a RED standard runner environment, not GREEN. Isolated core helper probe had no messages and did not install/update dependencies.
- Defects: none within approved scope. Minimal fix: none. Preserve the two limitations in the user-facing completion report.

Reproduction: `/Applications/Wolfram 15.0.app/Contents/MacOS/WolframKernel -noinit -script /tmp/qi-applied-quality-probe.wls` (quote the executable path containing a space). Source fragment `/tmp/qi-applied-quality-source.wl` is verbatim source 615–647; live output `/tmp/qi-applied-quality-live.log` and `/tmp/qi-applied-quality-live.wl`. The probe imports no audit prototype and applies no integration hook. The report does not count prior prototype logs as evidence that production tests ran: actual-package runner logs and source-bound probe are separate.

OPEN ISSUES: 0
