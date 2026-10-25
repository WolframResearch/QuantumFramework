# Standalone fixed-factor C30 loader

Build the attributed four-address implementation with no downloaded source
module dependencies. Preserve public W, H2_I, the official fixed shifts,
and the supplied tree; allow only the checked row-sign Givens gauge. No fitting.

Graph: pinned inputs -> sign-gauge and table invariants -> standalone build
-> full official verifier and independent dirty-target probes -> WL symbolic
and full-table checks -> fresh wl-verify -> fresh wl-quality -> report.
Executable quality changes return to correctness. Correctness has at most
three visits and at most three rounds per visit; quality has at most three
rounds. This is code, so there is no notebook build.

Required regimes: arbitrary four-case symbolic GF(2) interpolation; exact
signed-word and phase-repair identities; zero/constant-table limiting cases;
all released words and official numerical metrics; invalid inputs, dirty
targets, address boundaries, phase cleanup and Q128. The official verifier
uses all 3976 legal addresses, 64 measurement shots per address, and all
25 circuits in both directions. Its sampling does not certify every possible
measurement branch.

Push authorization covers the preceding reviewed audit. New loader work is
kept local for this build and verification stage.
