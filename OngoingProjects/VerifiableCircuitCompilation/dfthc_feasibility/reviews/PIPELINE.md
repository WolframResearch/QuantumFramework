# Bounded review history

Deliverable class: code; no notebook build required. The wl-ship graph used
three correctness visits and two quality-to-correctness back edges, within
its three-visit limit. Every numbered review used a fresh agent.

| Node | Round(s) | Exit and repair |
|---|---|---|
| Correctness visit 1 | 1/3 | One issue: exact numeric radicals/constants were rejected. Fixed the guard. |
| Correctness visit 1 | 2/3 | One issue: impossible scratch allocation leaked Missing into arithmetic. Fixed structured failure and propagation. |
| Correctness visit 1 | 3/3 | Zero issues. Exact numerical ordering was also made explicit before the final reviewed snapshot. |
| Quality | 1/3 | One issue: mapped Listable bit operations. Replaced with direct array threading; returned to correctness. |
| Correctness visit 2 | 1/3 | One issue: zero-imaginary complex representations lost spectral modes. Normalize accepted real values with Re; added regressions. |
| Correctness visit 2 | 2/3 | Zero issues; full reports reproduced. |
| Quality | 2/3 | One issue: no owned Q128/over-cap regression and no assertion of within_cap in legal rows. Added explicit assertions, corrected the regimes table, and returned to correctness. |
| Correctness visit 3 | 1/3 | Zero issues on the final source, including 32 built-in checks and 482 cost rows. |
| Quality | 3/3 | Zero issues on all axes A–G on the same final source. |

The final correctness report is `correctness-v3-1.md`; the final quality
report is `quality-3.md`. Earlier reports deliberately retain their issues
and historical source hashes. The evidence manifest binds the final files.
No source edit changed the reference spectrum, public metric acceptance, or
released resource projection. No full circuit verification was performed.
