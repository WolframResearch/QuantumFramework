# Independent wl-quality review, round 1 of 3

Scope: the sixteen public files directly in `OngoingProjects/VerifiableCircuitCompilation/dfthc_c30`, bound below. The review excludes `reviews/` and the not-yet-written final evidence manifest. I read the complete critic mandate first, the target files, the applicable empty AGENTS file, both WL idiom rubrics, the supplied installed-reference conversions, current live WL definitions, and the relevant unmodified official verifier implementation. No prior reviewer report was consulted. I did not edit the targets; the recorded end-of-review SHA comparison confirms that all sixteen were unchanged.

Verdict: one small idiom finding in `run.wls:10`; the physics, symbolic reasoning, required regimes, and reported numerical results pass. This is not a correctness failure.

## Evidence and reproduction

All regenerating commands ran in an isolated temporary copy with the required sibling JSON dependencies, never against the target directory. Existing public reports remained available as the reviewed snapshot. Runtime: Python 3.13.1, NumPy 2.5.4, kickmix 0.1; Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. `PYTHON`, `KERNEL`, `DATA` and `PUBLIC_FACTORS` below denote the supplied local runtimes and hash-pinned public input paths; these are not public path requirements.

Commands, from the isolated `dfthc_c30` directory:

```sh
"$PYTHON" generator.py --data "$DATA" --public-factors "$PUBLIC_FACTORS" --output local/rebuilt --report ../rebuilt.json
"$PYTHON" verify.py --data "$DATA" --submission local/rebuilt --report ../verified-rebuilt.json
"$PYTHON" probe_loader.py --submission local/rebuilt
"$KERNEL" -noinit -noprompt -script run.wls
"$KERNEL" -noinit -noprompt -script ../quality-probe.wls
"$KERNEL" -noinit -noprompt -script ../listability-probe.wls
"$PYTHON" ../numeric-quality-probe.py
```

Every command exited 0. No WL probe used `Quiet`, message suppression, or `Print`; the additional scripts end in the expression being checked and export their evidence. The fresh build is exactly equal to the reviewed `build.json` after removing only `elapsed_seconds`, including all 25 circuit hashes, the factorization hash, all fixed-array hashes, word hash, and all profiles. Elapsed provenance fields are not interpreted as performance evidence.

Local evidence beside this report:

- `rebuilt.json`, `verified-rebuilt.json`, and `reproduction-summary.json`;
- `dfthc_c30/probes.json` and `dfthc_c30/wl-checks.json`, freshly regenerated in the copy;
- `quality-probe.wls` / `quality-probe.json`;
- `listability-probe.wls` / `listability-probe.json`;
- `numeric-quality-probe.py` / `numeric-quality-probe.json`;
- `target-sha256.json`, the reviewed target binding.

The official rerun returns PASS, 254464 shots, 64 shots for each of 3976 addresses, 25 circuits in both directions, Q128 for ANGLE, 65244 full-block Toffolis, 661 full-block qubits, and cost 2051023203.3721037. This agrees with the unmodified official `circuits.py:438–552` implementation, which checks value consistency, phase, selection, and clean workspace after every forward and reverse circuit. Its default shot count is 64 at `circuits.py:399`.

## A. Documented symbols — PASS

Evidence: `LoaderIdentities.wl:1–45`, `run.wls:1–23`, `DOCUMENTATION.md:39–54`, the installed `PolynomialMod` and `BinaryReadList` reference conversions, and the live `quality-probe.json` usage/options/attribute results.

No third-party or QuantumFramework WL symbol occurs. The two locally defined exported symbols have usage messages and complete implementation authority. `FourAddress` returns the reconstructed four-row integer matrix on its stated domain and `Failure` otherwise; the additional live probe confirms endpoints and rejects ragged, empty, wrong-row-count, floating, and Boolean tables.

Live authority confirms `PolynomialMod[poly,{m1,m2,...}]` and the relevant options `CoefficientDomain -> Rationals`, `Modulus -> 0`, `Trig -> False`; the exact Boolean-reduction probe passes. Installed `BinaryReadList` documentation confirms the filename-taking form and automatic stream closure; live options include `ByteOrdering`, and the full-data run reads the little-endian unsigned 16-bit cache with dimensions `{3976,75}` and its pinned SHA256. Live `BitGet` usage is zero-indexed, and `BitGet` and `BitXor` have `Listable`; actual nested-array execution passes. No unsupported signature, option, or return-shape claim was found.

Smallest meaningful fix: none.

## B. Idiomatic functional WL — FINDING

B1 (minor idiom defect): `run.wls:10` explicitly threads a function already declared `Listable`:

```wolfram
deltas = Prepend[MapThread[BitXor,{Rest[batches],Most[batches]}],First[batches]];
```

The native-style rubric specifically calls out mapping functions already `Listable`. This extra traversal is unnecessary and obscures that the operation is a direct array XOR. The live `Attributes[BitXor]` result is `{Flat,Listable,OneIdentity,Orderless,Protected}`. The real-kernel full-data comparison in `listability-probe.wls` evaluates

```wolfram
MapThread[BitXor,{Rest[batches],Most[batches]}] ===
  BitXor[Rest[batches],Most[batches]]
```

to `True`, with both dimensions exactly `{24,3976,3}`.

Smallest meaningful fix:

```wolfram
deltas = Prepend[BitXor[Rest[batches],Most[batches]],First[batches]];
```

There is no semantic difference on the validated rectangular batch arrays. I made no target edit.

Other B checks pass: no `Do`, `For`, `While`, `AppendTo`, indexed table construction, or avoidable finite-case branching appears in either WL target. `FourAddress` uses `With`; the symbolic `Module` provides fresh symbolic indeterminates as well as intermediate expressions and is justified. `FoldList` and `Fold` express accumulation and reversal directly (`run.wls:11,15`). Dense integer lookup arrays are the natural representation because every entry is required; no sparse-matrix opportunity was missed. The `MapThread` in `LoaderIdentities.wl:20` applies a custom two-argument reduction, not a preexisting listable function, and is appropriate.

## C. Terse physics reasoning — PASS

Evidence: `DOCUMENTATION.md:3–37`, `LoaderIdentities.wl:6–42`, `run.wls:9–17`, `generator.py:38–60`, and `loader.py:49–132,188–225`.

The argument proceeds from four-address GF(2) interpolation to redundant-sign-bit compression, row-sign projector preservation, measurement-phase cancellation, then full-data chaining and reversal. Intermediate WL names identify these objects (`interpolation`, `signed`, `vectors`, `batches`, `deltas`, `forward`) rather than arbitrary temporary slots. All symbolic locals are used. The Python eraser names the low products, traversal carries, coefficients, and measured parity. The small amount of setup and export scaffolding is proportionate to a reproducible verifier. No dead definition or repeated physics structure requiring a separate quality fix was found.

Smallest meaningful fix: none beyond B1.

## D. Physics invariants — PASS

Evidence: `LoaderIdentities.wl:25–35`, `generator.py:53–59,79–85`, `loader.py:69–70,208–225`, `probe_loader.py:14–43`, fresh full official verification, and both independent probes.

The exact WL checks assert arbitrary-coefficient lookup reconstruction, integer HMR-phase cancellation, ANF phase linearity, dirty-target fanout cancellation, Boolean XOR involution, row-sign projector equality, unit Givens norm, all 512 sign-extension words, and compact/expand identity. All 17 target symbolic/edge checks and all six full-data checks return `True`. The additional real-kernel check

```wolfram
FullSimplify[Transpose[{{Cos[theta],-Sin[theta]},
  {Sin[theta],Cos[theta]}}].{{Cos[theta],-Sin[theta]},
  {Sin[theta],Cos[theta]}} == IdentityMatrix[2], Element[theta,Reals]]
```

returns `True`. `PolynomialMod[(xorPolynomial /. y -> xorPolynomial)-y,{x^2-x,y^2-y,2}] === 0` returns `True` with `xorPolynomial = x+y-2 x y`, exactly as executed in the supplied probe. Together with exact phase repair, the basis XOR involution supplies the applicable effective-oracle unitarity invariant.

The independent numeric probe decodes the original public words and freshly generated signed words using the pinned official Givens routine. Across every one of the 3976 vectors, it finds 2011 negative and 1965 positive row signs, maximum component drift after signs `7.771561172376096e-16`, unit-norm error `6.661338147750939e-16`, and directly computed projector drift `8.881784197001252e-16`. Exact array comparisons confirm unchanged weights, offsets, shifts, and tree. All 298200 angle words satisfy the signed-bit constraint. The fresh official metric gate passes with error_1 `0.013911731985012415` Ha and error_2 `0.01998067964067221` Ha.

The supplementary emitted-circuit run passes 795200 released and 4704 synthetic applications, including arbitrary dirty target words, involution, phase, address preservation, and zero scratch. The full official run additionally covers chaining and reversal at every legal address with 64 sampled measurement outcomes per address.

Hamiltonian propagation, density-matrix positivity/trace, and many-body eigen-relations are inapplicable: the declared artifact implements an address/word oracle and a row-sign gauge, not a many-body propagator, density operator, or eigensystem (`DOCUMENTATION.md:3–5`). No invariant failure or unasserted load-bearing physical assumption remains in that scope. Sampling is expressly not claimed to enumerate all measurement branches.

Smallest meaningful fix: none.

## E. Nontrivial examples — PASS

The base cases that could conceal errors are a zero/constant table, clean word targets, a single address, or no sign changes. The artifact goes beyond all of them.

Evidence: `LoaderIdentities.wl:38` checks an exact nonconstant three-word table containing both signed endpoints; `probe_loader.py:25,32–41,56–63` tests dirty targets, two oracle applications, and nonconstant emitted circuits at 16, 20, 32, and 128 address rows. The non-power-of-two case exercises padding/completion rather than only rectangular powers of two. The real FeMoCo case uses 3976 addresses, 75 angles, 25 circuits, and both row-sign sectors, with complete official forward/reverse coverage and supplementary dirty targets. All reruns pass.

Smallest meaningful fix: none.

## F. Symbolic-first derivations — PASS

Evidence: `DOCUMENTATION.md:7–11,21–36`, `LoaderIdentities.wl:17–35`, `generator.py:38–60,76–85`.

Four arbitrary coefficients remain symbolic through GF(2) reduction. The HMR phase identity retains arbitrary integer measurement and lookup values; Givens normalization retains a real angle; the row-sign projector identity retains symbolic vector components; fanout and XOR involution use exact Boolean polynomials. Exact signed-word identities precede physical numeric checks. The floating inverse-tree gauge is necessary to decode the already supplied numeric factorization and enforces the stated vector-drift threshold and official metric gates; it does not replace the symbolic projector argument or claim a general symbolic factorization.

The cost and geometric-mean figures are inherently numerical resource summaries, not prematurely specialized physics parameters. Independent arithmetic using the sibling verified JSON files reproduces cost ratios `3.9385957568764147`, `1.2581323878369532`, `1.6360612511835888`, and combined local estimate `2.0088882102784944`. It also reproduces `47.55854233515931 * 65244 * 661 = 2051023203.3721037`, prior-incumbent improvement `1.6154555940004633`, cost reduction `38.09795801792165%`, and two-body headroom `0.0193203593277895` mHa. No new leaderboard result or speed claim is inferred.

Smallest meaningful fix: none.

## G. Regimes and limits — PASS

The required table exists at `DOCUMENTATION.md:56–64` and is reiterated in `PLAN.md:14–20`. Each row has executable evidence:

| Required regime | Fresh evidence | Result |
|---|---|---|
| General symbolic | `SymbolicChecks[]`: arbitrary four-address coefficients, integer HMR phase, symbolic projector and Givens norm; independent real-angle rotation-unitarity probe | PASS |
| Exactly solvable | All 512 signed words, compact/expand identity, endpoint and nonconstant three-word tables, exact XOR involution and fanout | PASS |
| Limiting | Zero and constant tables in WL and emitted circuits at 16, 20, 32, 128 rows; differences vanish while arbitrary dirty targets are handled | PASS |
| Numerical reference | All 298200 physical words; independent original-versus-gauged vector/projector comparisons; full unmodified official 3976-address/64-shot/25-circuit forward/reverse run; exact profile/hash reproduction | PASS |
| Failure and edge | 14 malformed input/batch cases rejected by supplementary script; additional live WL empty/ragged/wrong-shape/type checks; signed endpoints; all legal addresses, Q128, dirty-target, phase, selection and scratch checks | PASS |

The zero/constant table is the relevant decoupling limit of this discrete fixed-width oracle. No continuous classical or thermodynamic limit is needed for the stated operation. Measurement sampling remains a declared limit, not an unverified exhaustive-branch claim.

Smallest meaningful fix: none.

## Target SHA256 bindings

| Target | SHA256 |
|---|---|
| `loader.py` | `1759ff52b36e4d67075faecffde2056e5c9b602526f2a01c6c718411d4875e6b` |
| `generator.py` | `d1c34e6a9397ec960a91ecdaf1b1e67953efd87ead0f9a009453b21daa3ee643` |
| `verify.py` | `6d09e92c31b74724efeb384c8364b82071f2dbd232f011880dc3a1ac984c959a` |
| `probe_loader.py` | `ec13b40da6c78c64798db4623c5012e879d4cfe8b9ea2de0118c39f36ad3806a` |
| `LoaderIdentities.wl` | `183c80d935f9ce05a255c61fee23771ba75594fbc3890a9408c3ea8d39584e9a` |
| `run.wls` | `9d2640ed8b541ed75b57fc4abd0ed419bca8c52a34b9043b4b38b3bf9f48387e` |
| `build.json` | `3d911ccb95018b64a9389ced0e1400d8f6f3cd3a61caf67c6fb9107456057a54` |
| `verified.json` | `8412005457b6ee96e17460309d3906b15d3da7843022f18a33a0ae6c787a2ef5` |
| `probes.json` | `a061aa0e837c55c6946e6d0566b788080765b823620e679110c7c4214365bf73` |
| `wl-checks.json` | `9c7b3dc1adce26cc040fa52a987af13d8fc92135ef423b343168e9af155badc2` |
| `README.md` | `21424e0c63d05950f248f66a063cf4354f81b5e164cb4e719269c4011ae91201` |
| `DOCUMENTATION.md` | `134ff3f328070e7de0d9337a32bb2e6fc077a8472d98e1618b2b90f1f974ddf6` |
| `RESULTS.md` | `23197cd93b7c69b3866af56d9e60e3e9880ca34466e489b8edc75c274b6c1b62` |
| `PLAN.md` | `db94bee88216f449856086cc25b15d53f95ace70fdf2326d8c405b92843ffef5` |
| `NOTICE.md` | `2f4134e8910d370b0a2f14189857c13cf5bae4c60ea9d34551d40a5107731581` |
| `LICENSE` | `50e6751797c50dedd75ef1b8a0d9e42f5f8472e9fbce91f34718e9f97b0c780a` |

No target changes were made. The sole unresolved issue is B1.

OPEN ISSUES: 1
