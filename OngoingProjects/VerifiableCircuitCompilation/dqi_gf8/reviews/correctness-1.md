# Independent DQI GF(8) correctness verification — ROUND 1/3

## VERDICT

PASS. Rebuilt the full candidate independently, reproduced its exact binary hash and resource counts, and reran the unmodified default full official verifier with one worker. All 102,400 candidate shots passed with zero value errors and zero phase errors. No target or sibling feasibility source/evidence bytes changed during verification.

The retained baseline report was inspected and hash checked; its full verifier was **not rerun** for this integration review. The reported relative improvement is compared with that retained historical record. Kernel/generated-code timings are observations, not portable speed claims. No notebook was supplied or built.

## REPRODUCED

### Scope and evidence isolation

Read in full `__init__.py`, `gf8_backend.py`, `generator.py`, `verify.py`, `test_integration.py`, `run-certificates.wls`, `certificates.json`, `build.json`, `verified.json`, and `README.md`, plus directly required `other_tasks/FieldArithmetic.wl`, `symbolic.json`, `field_probes.py`, `dqi-baseline.json`, and relevant profiler records. Read the complete external DQI decoder and official verifier, native arithmetic generator implementations for multiplication/inverse/division, their Python binding implementations and relevant binding declarations. No applicable AGENTS.md was present in the repository ancestry/project target.

All generator and verifier output went into `$MIRROR`. SHA-256 manifest comparison confirms every originally snapshotted `.py`, `.wl`, `.wls`, `.json`, and `.md` file in target/sibling feasibility folders is unchanged. The original ignored circuit was separately streamed through SHA-256 and matches the regenerated binary.

### Commands and runtimes

`$KERNEL`, `$PYTHON`, `$DATA` and `$MIRROR` identify the exercised direct kernel, pinned runtime, external competition bundle and temporary mirror. Personal absolute paths are omitted.

Working directory: `$MIRROR`.

```sh
"$KERNEL" -noinit -noprompt -script dqi_gf8/run-certificates.wls
QF_INVSQRT_DATA="$DATA" "$PYTHON" -m unittest dqi_gf8.test_integration -v
"$PYTHON" dqi_gf8/generator.py --data "$DATA"
"$PYTHON" dqi_gf8/verify.py --data "$DATA" --workers 1
"$KERNEL" -noinit -noprompt -script independent-probes.wls
"$PYTHON" independent-probes.py
"$PYTHON" binding-probes.py
```

Established commands completed successfully. Logs: `wl-certificates.log`, `python-tests.log`, `full-build.log`, `full-verify.log`, `independent-wl.log`, `independent-python.log`, `binding-probes.log`, and `boundary-decoders.log`. Diagnostics were captured without Quiet or message suppression. Supplemental WL code uses Echo.

Kernel: 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. Python: 3.13.1, Clang 16.0.0. NumPy: 2.5.4. Native kickmix: 0.1; source HEAD `638777fac4379bc62e2a1e4a52cf1a1753895d48`, source checkout clean.

### Exact claims, signatures and boundaries

- `FieldArithmetic.wl:12–19,21–56`: regenerated GF(8) six-product polynomial identity and seven-subproduct native GF(4096) tower plan; complete 12×12×12 binary multiplication tensor. Independent polynomial-only Python arithmetic checked all 64 GF(8) input pairs and all 1,728 native binary coefficients without using native GF multiplication for expected values.
- `run-certificates.wls:7–29`: all 11 exact checks PASS. Regenerated `certificates.json` is **byte-for-byte identical** to the target (SHA-256 `124e09a628d9e187b7eb9782393543e58097cd1c4502607cb3874d180b837c14`). The arbitrary positive-integer chain identity, scaffold exponents, output exponent, all nonzero inverses, zero extension, exact forward/inverse maps and Frobenius period were reproduced.
- Live field metadata confirms characteristic 2, degree 12, size 4,096, polynomial representation and native modulus x^12+x^6+x^4+x+1. Live `Inverse` options include Modulus; modular matrix inversion was exercised. Coefficients are low power first; zero has an empty coefficient list, correctly padded by the implementation. Installed authoritative FiniteField/FiniteFieldElement documentation agrees.
- Map shapes: seven 3×12 evaluation maps and seven 12×3 output maps. Every configured Frobenius map was simulated on all 4,096 field words by the integration tests.
- `gf8_backend.py:140–182`: direct simulator readout independently confirmed output power 4094 and ordered scaffold powers 3,15,1020,2046 on every one of the 4,096 inputs, including zero. Inverse returns `(target[12], scaffold[48])`; division returns `(target[12], scaffold[60])`, beginning with the inverse followed by the 48 chain bits. The native C++ ordered addition-chain contract matches these words.
- `test_integration.py:17–133`: all six tests PASS (0.972 s observed). Coverage is 8,192 dirty-target multiplication shots for each unconditional/quantum-control mode; all 4,096 inverse inputs with independent randomly drawn masks; every divisor in both disabled/enabled classical-condition branches; every configured Frobenius map; false control and malformed inputs; and the reduced 2,048-shot official decoder. Input/control/mask preservation, expected phase and zero scratch assertions passed.
- Supplemental inverse-phase test crossed all 4,096 inputs with each of the 12 independent mask basis bits (49,152 phase cases). Expected inverses came from independent polynomial arithmetic. Source/mask values were preserved and scratch was zero. Native inverse and division cleanup added zero static Toffolis; complete compute/cleanup totals remain 210 and 252 respectively. Arbitrary masks combine these basis bits by parity; the supplied test does not enumerate all 4,096² input-mask pairs or every measurement branch.
- Wrong-width/field, invalid target, singular/wrong-shaped/nonbinary map, invalid dimensions, and empty inverse were rejected. WL TowerPlan rejects empty/overfull signatures, 0,2,3.0,"3", a symbolic argument and an empty list. Constructor binding probes reject changed certificate source, field source and symbolic data, and false certificate flags.
- Four additional minimum/odd-dimension decoder cases `(m,n)=(3,2),(4,3),(5,2),(7,5)` each passed the unmodified official verifier for 2,048 shots with zero value/phase errors.
- `generator.py:20–41`: external reference SHA check passes; default n_iters is the imported unchanged function (222 at n=70). Only the three scoped arithmetic computation entry points are adapted; native masks, phase routine, cleanup and decoding algorithm remain delegated/imported. This is a narrow adapter requiring disjoint native-width registers, as documented, rather than a complete native builder API replacement.

Three cold-cache `TowerPlan[3]` evaluations, each preceded by `ClearSystemCache[]`, observed 0.016651, 0.015712 and 0.016707 seconds (kernel startup excluded). These are observations; exact resource claims below derive from the independent full circuit build and official execution.

### Full decoder and official cost

The independent build matches every stored build field except observed elapsed time (37.446372999984305 s versus 39.97367529204348 s). Exact circuit bytes: 1,178,484,163. SHA-256:

`bc54c73b8725acbf6d246df5205f6a2dff2ea5783f335154833e57d9d93c8668`

| Claim | Fresh result |
|---|---:|
| m / n / degree / r | 4095 / 70 / 12 / 2016 |
| Qubits | 1903 |
| Static Toffolis | 4503030 |
| Operations | 117029650 |
| Classical bits | 49154 |
| GF(8)-backed multiplications | 58,211 |
| Controlled products | 32,408 |
| Division scaffolds | 888 |
| Inverse scaffolds | 4,095 |
| Official verified shots | 102400 |
| Value / phase errors | 0 / 0 |
| Mean executed Toffolis | 4502782.126171875 |
| Repetition factor | 1.1912604648115832 |
| Effective error bound | 0.0008660587301657926 |
| Reaction depth | 174852 |
| Official cost | 19425188366155.996 |
| Historical baseline cost | 24439458623294.45 |
| Improvement ratio | 1.2581323878369532 |
| Cost reduction, percent | 20.517108559676135 |

All original terminal report result fields, including per-operation totals, exactly reproduce; only the observed duration and the regenerated build.json binding differ. The build.json binding differs because its recorded duration differs. Full run duration observed: 460.92019191698637 seconds. Default 102,400 shots, batch size 256, log2(delta)=-128 and one worker were used, without a seed override. Fixed salt is `qf-other-tasks-feasibility-v1`; circuit hashes determine separate official samples for candidate and baseline.

Independently checked arithmetic decomposition `32408 + 888*6 + 4095*5 = 58211`, static model `5725461 - 58211*(63-42) = 4503030`, official formula `R*Q^2*T`, retained baseline ratio and percentage reduction. The retained baseline is Q=1,893, mean Toffolis=5,725,107.220107422 and cost=24,439,458,623,294.45; those historical numbers were not freshly simulated.

All original report/source/circuit bindings and post-run external official input checks pass. No external reference/verifier code changed. No public artifact contains personal absolute paths or copied private/external source; personal absolute paths have been redacted from this public review record. No commit, push, branch change or submission was performed.

## FINDINGS

No unresolved correctness findings in the stated scope. No missing runtime dependency, unexecuted claimed candidate result, kernel limitation or timeout reduced the substantive verification coverage. Official sampling and sampled measurement branches are the explicitly documented coverage limits; no global optimizer, asymptotic speedup or new Kaggle submission score is established or claimed.

## RECOMMENDED FIXES

None required for correctness. Preserve the distinction between the freshly rerun candidate and retained historical baseline when summarizing the result, and preserve the documented sample coverage rather than describing all input-mask pairs or all measurement branches as exhaustively simulated.

OPEN ISSUES: 0
