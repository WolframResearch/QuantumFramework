# Evidence status for the QF formal-semantics document

This is the claim-by-claim evidence record for `QF-formal-semantics.md`, kept here beside the check scripts rather than inside the shared document. It is a point-in-time audit against kernel `2d3c29e8` (2026-09-22), not a certification: every claim in the document must be re-confirmed against the live kernel source, line by line, before it is relied on, since the kernel moves and a status here reflects only the anchor commit.

Each main-text claim is one of: a definition; a source-derived contract, read off the kernel; a proposition, proved in the text; a tested instance, one of the check scripts beside this file passing on named inputs; or untested. A passing check on named inputs supports the claim it tests and does not establish a universal statement.

| Claim | Status | Evidence |
|---|---|---|
| Return heads (Appendix A) | tested instance | `return-shapes-check.wls`, `section3-audit-check.wls` |
| Composition: sorting, boundary, dimension check, silent picture fallback, merging the operands' parameters | source-derived contract; tested instance | `section3-audit-check.wls` |
| Application: padding to the full input order, untouched subsystems, the qutrit failure | tested instance | `section3-audit-check.wls` |
| $\operatorname{rep}(A\circ B) = \operatorname{rep}(A)\operatorname{rep}(B)$ across bases | source-derived contract; tested instance | the kernel suite's cross-basis composition tests |
| Adjoint and trace are coordinate operations, physical only on orthonormal, respectively equal, bases | proposition (Gram factors); tested instance | `section3-audit-check.wls` |
| Maps on density matrices as operator-space matrices (`"Left"`, `"Right"`, `"Liouvillian"`), channel action | source-derived contract; tested instance | `section4-audit-check.wls` |
| Post-measurement state reduces to the target; branch states on the full register | tested instance | `section3-audit-check.wls` |
| Measurement list constructor forms $\sqrt{E_m}$ from effects; channel list constructor keeps Kraus operators | source-derived contract; tested instance | `section4-audit-check.wls` |
| Phase-space transform preserves composition and application; inverse and rate matrix | proposition (similarity of operator-space matrices); tested instance | `phase-space-functor-check.wls`, `section4-audit-check.wls` |
| Even-$d$ Wigner elements, Gram $4d\cdot I$, unequal traces; `QuantumWignerTransform` and `QuantumPhaseSpaceTransform` differ on a channel | tested instance ($d = 2, 3, 4$) | `section4-audit-check.wls` |
| Tableau conversions keep the global phase; state to tableau to state is the identity with phase, tableau to state to tableau returns the same group | source-derived contract; tested instance | `tableau-representation-check.wls` |
| Method agreement | tested instance (Schrödinger against tensor network on one circuit; tableau against dense on two) | the kernel suite |
| Distances and entanglement: default measure `"Fidelity"` is the distance $1 - \operatorname{Re}\operatorname{Tr}\sqrt{\rho_a\rho_b}$; default monotone concurrence, `QuantumEntangledQ` realignment | source-derived contract; tested instance | `derived-functionals-check.wls` |
| The eight live defects (two earlier ones now fixed) | source-derived fact, reproduced at `2d3c29e8` | `defect-repro.wls`, `section3-audit-check.wls`, `section4-audit-check.wls` |
