# Released-design audit and DFTHC feasibility

Scope: audit and necessary-condition/resource screening only. No factor training,
new submission circuits, commit, or push. Existing verified results stay fixed.

Graph: pin source versions and competition hashes -> check pair packing against
the official full-tensor residual -> one shared eigenspectrum -> twelve shape
screens and precision/batch cost sweep -> fresh wl-verify -> fresh wl-quality
-> recommendation. A packing failure returns to packing; an ambiguous numerical
threshold gets one refinement. Each independent review has a three-round cap.
Executable changes from quality return to correctness.

Use the official fixed shifts and float32-rounded FCIDUMP. Pack symmetric pairs
with sqrt(2) on off-diagonal coordinates. The approximation is PSD and has rank
at most R min(B,C); its unconstrained PSD tail is only a necessary error bound.
Passing this bound does not establish DFTHC representability, quantized accuracy,
the one-body residual, normalization, or gap acceptance.

Costs include official state preparation, both Givens passes, ordinary unary
angle lookup, shared scratch, and the 128-qubit ANGLE cap. New-shape costs are
conditional projections, never verified circuit scores. Hold lambda_eff at the
incumbent value solely to compare resources. Audit any directly compatible
public word table through the official metric evaluator before relying on it.

Keep raw integral matrices, spectra, downloaded factors and source in ignored
local storage or temporary audit storage. Public evidence contains code,
aggregate results, source hashes and attribution, with no personal paths.
