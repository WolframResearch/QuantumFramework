# Opinion on the Codex audit of Modeling-Slides.md

> **Status, rechecked 2026-11-10 against the current Modeling-Slides.md.** Already applied in the deck: the cat in a warm bath (evolves ρ_cat − ρ_mix), CKW ("three-tangle vanishes as W's does"), the Helstrom case with two negative eigenvalues, the QASM commas, the non-orthonormal basis wording, and the natural-gradient pole wording. The Young-diagram degeneracy, capacity sentence, repetition-code identity and bound-entanglement wording were not rechecked. Read the corrections below as history, not as open items.

The audit is careful and mostly right. Its physics corrections are real in most cases, it overreaches in a few, and its runtime picture is wrong: the teleportation "stall" does not reproduce, and the 81 "unverified" examples run. The editing direction (shorter prose, the result-bearing formula in the text, fewer panels, complete keys, shared colour scales) suits a catalogue deck. Its requests for added evidence (error bars, residual plots, acceptance curves, extra fidelity panels) pull against "minimal code, one figure" and should be taken selectively.

Both audit files (Report.md, 241 lines; Slide-by-slide.md, 192 lines) were read in full. The deck checked is the current file (SHA-256 37bbd8a6..., the same one Codex read). Kernel checks used `wolframscript` against the local QF paclet; the scripts are in this session's scratchpad.

## Corrections that are right and matter

**Cat in a warm bath (slide 15).** This is the most important one. The cell normalises the trace norm of the cross term C, but C has trace 2e^{-2α²} and the evolution preserves trace, so the plotted ratio can never fall below e^{-2α²}. The kernel confirms it: at α=1, γ=0.01, k_BT=20 the curve ends at 0.140 against a floor of 0.135, while the dashed point-lobe prediction is 4×10⁻⁵; the trace of the evolved C stays at e^{-2} to twelve digits. At α=2 the floor is 3.4×10⁻⁴, inside the plotted range, so it shapes that panel too. Codex's fix (evolve the traceless difference ρ_cat − ρ_mix) is the right one. With that fix the closing sentence turns out to be true after all: the traceless difference also decays far more slowly than the point-lobe prediction once most of the coherence is gone (at α=2, γ=0.01, k_BT=20 it ends at 8×10⁻⁴ against a predicted 3×10⁻¹⁸), so the error was only in the plotted quantity, not in the sentence.

**CKW monogamy (slide 106).** "A superposition that, like W, has no three-way entanglement" is wrong: W has zero three-tangle and is genuinely tripartite entangled. Say the three-tangle vanishes, as it does for W. Codex is right; the fix is one sentence.

**Plain against natural gradient (slide 66).** "At the pole θ=π/2, where both coordinate derivatives vanish" is wrong: ∂_θE = −2cos2φ there. The kernel shows where plain descent actually stops from (π/3, π/12): (θ,φ) → (π/2, π/4) at both learning rates, the one point on the pole where cos2φ=0, a saddle with E=0. The sentence should name that point.

**Exchange energies from Young diagrams (slide 133).** The question says "for every n and d" and the text gives degeneracy f^λ s_λ(d). At n=6, d=3 the diagrams (4,1,1) and (3,3) both have content sum 3; the kernel's eigenvalue tally gives degeneracy 150 = 100 + 50 at E=3J. The plotted panels (n ≤ 5) are correct; the general sentence needs the sum over diagrams with equal content.

**Repetition-code threshold (slide 65).** The deck says `Reduce` gives the threshold "for every pair of neighbouring codes" but checks three pairs. Codex's closed form P_L(2m+3) − P_L(2m+1) = C(2m+1,m)[p(1−p)]^{m+1}(2p−1) is correct (verified symbolically for m = 0..6) and proves it for all odd n in one line. Worth adding.

**Smaller true gaps.** Helstrom (slide 120): "whenever Λ has a negative eigenvalue the best rank-one projector touches the bound" fails when both eigenvalues are negative (then E₂ = I); the plotted priors p ≥ 1/2 never reach that case, so it is a wording fix. Two-qubit bound entanglement (slide 118): "nothing is left to be bound" also needs the fact that every two-qubit NPT state is distillable. Nonorthogonal basis (slide 130): with non-orthonormal F the map FCF⁻¹ is not Hermitian, so "observable" is loose there. OpenQASM (slide 73): `swap q[0] q[3];` and the controlled lines lack the commas the standard requires; QF's importer tolerates it, a reader copying it elsewhere will not.

## Right in substance, but the fix overreaches

**Capacity (slide 109).** The sentence "the maximally mixed input maximizes I_c for both channels" is false where I_c(I/2) < 0 (at p=1 it is −1, a pure input gives 0), and "commutes with every Pauli" does not justify it. But the plotted curve max(0, I_c(I/2)) is the true one-use rate: maximizing over the input spectrum in the kernel, I/2 is the maximizer at every p where the maximum is positive (p ≤ 0.25 sampled, up to the zero near 0.252), and above that the maximum is 0. So Codex's proposal to relabel the curve as only "an achievable lower bound at I/2" understates what is drawn. Fix the sentence (I/2 maximizes I_c wherever it is positive; past the zero Q⁽¹⁾ = 0), keep the curve.

**Tomography reference (slide 34).** The dashed pure-state reference is derived for projecting the linear estimate onto the Bloch ball, while the points use QF's likelihood estimator. Codex is right that this is a comparison, not the estimator's own formula, but the deck already describes the estimator honestly and says the points "scatter about" the reference. Renaming the dashed curve "projected linear-inversion reference" is enough.

**DRAG zero (slide 138).** True in principle: a 40 ns, 1 ns-sampled window has no exact zero at ω=1/β. In practice the sampled power at ω=0.2 is 8×10⁻⁶ against a peak of 157 and a plot floor of 10⁻³, two orders below anything visible. The slide's "puts a zero at ω=1/β" describes the closed-form line, which does have the zero. At most say "the closed form has a zero".

**Scope labels (Shor, Bianchi, curved gamma matrices, data processing, phase certificate, golden-rule window, thermal Jaynes-Cummings).** Each is a fair tightening of a quantifier, none changes a figure. The thermal Jaynes-Cummings slide already asks for long-time averages, so that one is nearly moot. Do these in the prose cut, not as separate fixes.

## Convention, not error

Fidelity: QF's `QuantumSimilarity[..., "Fidelity"]` returns root fidelity; slide 114 already squares it, and elsewhere the fix is naming which one is plotted, not a correction. Zeno threshold `\Omega/2\sqrt2`: the value is right (Ω/(2√2)), the LaTeX just reads ambiguously. Schrödinger-Newton E_total = E/3: correct (virial: 2T + W = 0 gives eigenvalue 3W/2 and total energy W/2), and only matters if the slide calls the eigenvalue an energy of the body.

## Runtime claims

The teleportation cell does not stall. Run alone in a fresh kernel against the local QF it finishes in 6.6 s with no messages, and gives Bob's Bloch vector as the input, (0,0,cos β), (sin β cos ϕ,0,0) and 0 for the four correction choices, exactly as the slide says. Codex's harness started the kernel with `-noinit`, loaded TensorNetworks and Arrays from the paclet repository instead of the deck's own loading cell, and ran everything inside a custom context with a 180 s limit that did not fire; the hang belongs to that setup, not to the slide. A sequential run of the current deck in one kernel (setup plus all 139 examples, in order, as Codex ran them, with a 300 s limit per cell) finished in about 8.5 minutes: every cell returned, none timed out, and none raised a message, with teleportation at 2.9 s. So the 81 examples Codex marks "not run" all evaluate cleanly. That means they run, not that every claim in them is right; the physics points above for slides 65, 66, 106, 109, 118, 120, 130 and 133 sit in that unrun range and stand on their own.

## What to do with it

Fix the cat observable first (it changes what the figure means), then the five one-sentence physics errors (CKW, natural gradient pole, Young-diagram degeneracy, Helstrom case, capacity sentence), and add the one-line repetition-code identity. Then apply the prose cut, keeping one sentence per slide that says which QF or WL construct does the work, since showing that is what the catalogue is for; Codex's template moves that sentence to notes, which would remove the deck's reason to exist. Take the panel reductions and shared colour scales; skip most of the added-evidence requests.
