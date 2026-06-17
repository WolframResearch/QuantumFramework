# Version anchor

The `file:line` citations throughout this reference set are **point-in-time**: they were read
and verified against one specific snapshot of the paclet. Use them as a fast index, not as
ground truth. When an exact line matters, re-check it against the kernel you actually have.

| Field | Value |
|---|---|
| Paclet | `Wolfram/QuantumFramework` |
| Version | `2.0.0` |
| Verified against commit | `a601a187` (`WolframResearch/QuantumFramework`, 2026-06-17) |
| Catalog counts (live) | `$QuantumOperatorNames` 108 · `$QuantumCircuitOperatorNames` 44 · `$PauliStabilizerNames` 9 |

## How to re-verify a cite on your install

The catalog and most behavior notes are stable across patch releases; only the line numbers
drift. To confirm a specific `file:line` (e.g. `QuantumState/Properties.m:461`):

1. **Find the kernel.** The implementation lives under `<paclet>/Kernel/`. Locate your paclet
   with `PacletObject["Wolfram/QuantumFramework"]["Location"]`, or load a local checkout:

   ```wolfram
   PacletDirectoryLoad["/path/to/QuantumFramework/QuantumFramework"];
   Needs["Wolfram`QuantumFramework`"]
   ```

2. **Read the source, not `Usage.m`.** `Kernel/Usage.m` is auto-generated from the doc pages
   and is never evidence of behavior. Read the actual `Kernel/<subdir>/<file>.m`.

3. **Or test it.** A two-line `wolframscript -file` against the live kernel settles any
   signature or return-shape question faster than trusting a cite. Never `Quiet`; never `Print`
   in committed code.

If your version is newer and a symbol moved, the *behavior* note in this set is still the best
starting point; just treat the line number as approximate and grep for the symbol.
