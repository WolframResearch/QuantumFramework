# Independent WL quality review, round 2/3

All sixteen named public targets were read in full. Prior reviewer reports and `reviews/` were not consulted. Applicable `AGENTS.md` locations were searched; no nonempty applicable instructions were found. Targets were not modified. All regeneration was performed in an isolated temporary copy, using the supplied Python environment, installed Wolfram 15.0 kernel, hash-pinned competition dependencies and public factors. The binding below identifies the original reviewed bytes, not regenerated report timestamps.

The operation under review is the four-address, three-word signed-angle XOR oracle with measurement-assisted scratch cleanup. Full many-body state propagation is outside its claim. This review applies the independent critic mandate and the native/style rubrics; no missing kernel or dependency limited the checks.

## A. Documented signatures, options and returns — PASS

`LoaderIdentities.wl:1-15` defines only two project exports, with explicit integer/shape guards and a `Failure` fallback; every external WL symbol is core. Installed Wolfram 15.0 `PolynomialMod` documentation supports reduction by an integer and a list of polynomial/integer moduli, including Boolean idempotence. Installed `BinaryReadList` documentation supports a filename, element type and byte-order option, returning the read list and closing a filename-opened stream. These match `LoaderIdentities.wl:20,27-30` and `run.wls:6`.

Independent fresh-kernel probes returned `BitGet[1,0]=1`, `BitGet[1,1]=0`, `Listable` for both `BitGet` and `BitXor`, and `ByteOrdering` among `BinaryReadList` options. `Options[PolynomialMod]` returned `{CoefficientDomain -> Rationals, Modulus -> 0, Trig -> False}`. The supplied documented binary read executes over all 298200 words. No recalled QuantumFramework API is involved. Defects: none. Smallest fix: none.

## B. Idiomatic functional WL — PASS

`LoaderIdentities.wl:9-15` uses guarded pattern dispatch and `With` for the base and differences. Lines 20-21 use `MapThread` for aligned symbolic substitutions; lines 31-35 use `Outer` and listable bit operations. `run.wls:9-15` uses structural partition/transpose operations, `FoldList` for forward accumulation and `Fold` for reversal. A source search found no `Do`, `For`, `While`, `AppendTo`, index-driven `Table`, `Quiet` or `Print` in either WL target. The `Module` at line 17 is justified by unique formal symbols, not merely constant storage. The physical word arrays are dense integer tables; sparse storage is not their natural representation. Defects: none. Smallest fix: none.

## C. Terse physics reasoning — PASS

`DOCUMENTATION.md:7-38` presents interpolation, signed-bit compression, register allocation, sign gauge and measurement phase cancellation in that order. `LoaderIdentities.wl:19-35` follows the same argument using named objects (`interpolation`, `signed`, `vectors`, `xor`), and `run.wls:9-15` directly expresses chaining and reversal. All definitions participate in a reported identity or data check. The prose distinguishes the Q128 ANGLE allocation from the 661-qubit full block encoding and bounds what the measured evidence proves. Defects: none. Smallest fix: none.

## D. Applicable physical invariants — PASS

Fresh execution of `run.wls` passed all 17 symbolic/edge identities and six full-data checks. The latter reconstruct all 3976 rows/298200 words, check signed-bit redundancy, batch accumulation, reverse clearing and all 25 Q128 profiles. Independent symbolic derivations returned `{a,b,c,d}` from the four Boolean addresses for arbitrary coefficients, `0` for Boolean XOR involution, `1` for `(-1)^(m f) (-1)^(m f)` with integer assumptions, and `0` for `(s x)(s y)-x y` under `s^2==1`. Independent sparse permutation matrices for all 16 one-bit four-address truth tables each satisfied `U.Transpose[U] == IdentityMatrix[8,SparseArray]` exactly. The logical XOR permutation therefore preserves normalization; the target's exact phase-cancellation and fanout identities address coherent phases.

The clean generator reproduced maximum sign-adjusted vector drift `7.771561172376096e-16` and unit-norm error `6.661338147750939e-16`, with identical fixed-array hashes. Official metric gates passed. The fresh unmodified official verifier passed 254464 shots at 64 shots for each legal address, with all 25 circuits traversed forward and backward. Its source `scoreinfra/circuits.py:402-416,438-559` explicitly rejects value disagreement, changed selection, residual workspace and uncorrected phase after each circuit. Fresh dirty-target probes additionally checked 795200 released and 4704 synthetic circuit applications with zero value, phase, workspace or selection errors.

Hamiltonian Hermiticity, conserved commutators, density-matrix positivity/trace and many-body eigen-relations are inapplicable to this declared reversible word-loading operation; no such propagated operator or density matrix is produced. Measurement tests remain sampled, as the artifacts explicitly state; no exhaustive branch certification was inferred. Defects: none. Smallest fix: none.

## E. Nontrivial examples — PASS

The possible hiding cases are zero/constant lookup tables, a single compact bit or a clean target. `LoaderIdentities.wl:38` instead checks an exact nonconstant four-row, three-word table containing both signed endpoints. `probe_loader.py:22-42,54-65` starts arbitrary ten-bit dirty targets, applies each oracle twice and checks phases and scratch; zero, constant and nonconstant tables run at 16, 20, 32 and 128 rows. The released case adds 3976 rows and 25 chained batches. The non-power-of-two 20-row case exercises a real address boundary. These tests all reran successfully. Defects: none. Smallest fix: none.

## F. Symbolic-first reasoning — PASS

`LoaderIdentities.wl:17-35` retains arbitrary four-row coefficients, arbitrary integer measurement/sign parameters and a real Givens angle; reductions are exact. All 512 signed words are enumerated exactly, before released numeric data checks. `generator.py:37-60` necessarily performs a numerical reparameterization of fixed released angle words; its grid rounding and inherited branch threshold are constrained by sign-adjusted drift and official accuracy gates, rather than substituting numbers for the GF(2) or phase derivations. No premature numerical symbolic substitute was found. Defects: none. Smallest fix: none.

## G. Declared regimes and limits — PASS

| Regime | Fresh executed evidence |
|---|---|
| General symbolic | Arbitrary four-row coefficients reconstructed exactly; integer HMR phase, symbolic sign projector, real-angle norm and dirty-target fanout all passed. |
| Exactly solvable | All 512 signed words compact/expand exactly; nonconstant exact three-word table and XOR involution passed. |
| Limiting | Zero and constant tables passed in WL and emitted-circuit probes at four address sizes. |
| Numerical reference | Every released word and official metric gate passed; full default 64-shot verifier reproduced cost `2051023203.3721037`, `65244` Toffolis, `661` full qubits and `128` ANGLE qubits. |
| Failure/edge | Fourteen malformed shape/type/range/batch cases rejected in emitted-code probes, plus WL invalid signed/float/empty/negative cases; endpoint, dirty-target, address, scratch, phase and reversal checks passed. |

The regimes table is explicit in `DOCUMENTATION.md:56-65`; it has no missing applicable row. The zero/constant table is the relevant decoupling limit for an exact finite XOR oracle; an unrelated continuum or large-system limit is not required.

The sibling verified costs independently reproduce inverse-square-root, DQI and DFTHC ratios `3.9385957568764147`, `1.2581323878369532`, and `1.6360612511835888`; their geometric mean is `2.0088882102784944`, agreeing with `RESULTS.md` and `verified.json`. Incumbent improvement `1.6154555940004633` and reduction `38.09795801792165%` also reproduce. These are correctly labeled local estimates, not a new leaderboard score. Defects: none. Smallest fix: none.

## Reproduction evidence

The following commands completed with exit code 0 in the clean copy:

```text
python generator.py --data $DATA --public-factors $PUBLIC_FACTORS
python verify.py --data $DATA
python probe_loader.py
WolframKernel -noinit -noprompt -script run.wls
WolframKernel -noinit -noprompt -script ../independent.wls
```

`comparison.json` records that regenerated build content equals the reviewed build except elapsed provenance; probes equal the reviewed report exactly; WL results agree; and the official verification agrees except elapsed provenance and the build hash that transitively includes it. All original embedded source hashes agree, and all sixteen original target hashes were rechecked unchanged at close. No speed claim was evaluated from elapsed fields.

`independent.wls` and `independent-kernel.json` retain the independent exact commands and results. Regenerated public reports and generated binaries are in the isolated `dfthc_c30/` copy.

## SHA256 binding of original targets

```text
21424e0c63d05950f248f66a063cf4354f81b5e164cb4e719269c4011ae91201  README.md
134ff3f328070e7de0d9337a32bb2e6fc077a8472d98e1618b2b90f1f974ddf6  DOCUMENTATION.md
23197cd93b7c69b3866af56d9e60e3e9880ca34466e489b8edc75c274b6c1b62  RESULTS.md
db94bee88216f449856086cc25b15d53f95ace70fdf2326d8c405b92843ffef5  PLAN.md
2f4134e8910d370b0a2f14189857c13cf5bae4c60ea9d34551d40a5107731581  NOTICE.md
50e6751797c50dedd75ef1b8a0d9e42f5f8472e9fbce91f34718e9f97b0c780a  LICENSE
1759ff52b36e4d67075faecffde2056e5c9b602526f2a01c6c718411d4875e6b  loader.py
d1c34e6a9397ec960a91ecdaf1b1e67953efd87ead0f9a009453b21daa3ee643  generator.py
6d09e92c31b74724efeb384c8364b82071f2dbd232f011880dc3a1ac984c959a  verify.py
ec13b40da6c78c64798db4623c5012e879d4cfe8b9ea2de0118c39f36ad3806a  probe_loader.py
183c80d935f9ce05a255c61fee23771ba75594fbc3890a9408c3ea8d39584e9a  LoaderIdentities.wl
d87dc6c1eb27ad11d01d059d345c1efff464e7dce909d74b3786188eb4c0774a  run.wls
3d911ccb95018b64a9389ced0e1400d8f6f3cd3a61caf67c6fb9107456057a54  build.json
8412005457b6ee96e17460309d3906b15d3da7843022f18a33a0ae6c787a2ef5  verified.json
a061aa0e837c55c6946e6d0566b788080765b823620e679110c7c4214365bf73  probes.json
ee8f9a6e8245456c00b0bb141c2df670d9f9e85c8c176c0da685bf6261c0bbe0  wl-checks.json
```

OPEN ISSUES: 0
