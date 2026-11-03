"""Assemble drafts/exact/P-section.md from the input fragment (titles, Topic, Appears in,
formula, Conditions, non-formula records), entries.py (text), cells/NN.wl (code, verbatim)
and out_cells.txt (kernel output of exactly those cells)."""
import re, pathlib, sys
here = pathlib.Path(__file__).resolve().parent
drafts = here.parent.parent
src = (drafts / "P-section.md").read_text()
sys.path.insert(0, str(here))
from entries import E

body, nonformula = src.split("## Non-formula records", 1)
blocks = re.split(r"(?m)^### ", body)[1:]
assert len(blocks) == 95, len(blocks)

outs = {}
for line in (here / "out_cells.txt").read_text().splitlines():
    m = re.match(r"^(\d\d) \[[^\]]*\] :: (.*)$", line)
    if m:
        outs[int(m.group(1))] = m.group(2)

SHORT = {
    11: "{E^(al*SuperDagger[\\[FormalA]1])**E^(-(\\[FormalA]1*Conjugate[al]))/E^((al*Conjugate[al])/2), E^((al*Conjugate[al])/2)*E^(-(\\[FormalA]1*Conjugate[al]))**E^(al*SuperDagger[\\[FormalA]1]), <same two forms again>}",
    35: "{E^((E^(I*ph)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\\[FormalA]1], 2]*Tanh[rr])/2)**Sech[rr]^(1/2 + SuperDagger[\\[FormalA]1]**\\[FormalA]1)**E^(-1/2*(\\[FormalA]1**\\[FormalA]1*Tanh[rr])/E^(I*ph)), E^(-1/2*(\\[FormalA]1**\\[FormalA]1*Tanh[rr])/E^(I*ph))**Sech[rr]^(1/2 - \\[FormalA]1**SuperDagger[\\[FormalA]1])**E^((E^(I*ph)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\\[FormalA]1], 2]*Tanh[rr])/2)}",
    40: "{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\\[FormalA]1]**SuperDagger[\\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\\[FormalA]1]**\\[FormalA]1)**E^((\\[FormalA]1**\\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), <the same three factors with xi -> rr E^(I th)>, Sqrt[Sech[Abs[xi]]]}",
    46: "{E^(-1/2*(z*SuperDagger[\\[FormalA]1]**SuperDagger[\\[FormalA]1]))**(\\[FormalA]1 - z*SuperDagger[\\[FormalA]1]), E^(<gen2>)**(\\[FormalA]1 - z11*SuperDagger[\\[FormalA]1] - z12*SuperDagger[\\[FormalA]2]), E^(<gen2>)**(\\[FormalA]2 - z12*SuperDagger[\\[FormalA]1] - z22*SuperDagger[\\[FormalA]2])}",
    62: "{1/4 + SuperDagger[\\[FormalA]1]**\\[FormalA]1/2, (-1/4*I)*GeneralizedPower[NonCommutativeMultiply, \\[FormalA]1, 2] + (I/4)*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\\[FormalA]1], 2], -1/4*GeneralizedPower[NonCommutativeMultiply, \\[FormalA]1, 2] - GeneralizedPower[NonCommutativeMultiply, SuperDagger[\\[FormalA]1], 2]/4}",
}

def field(block, name):
    m = re.search(r"(?m)^" + name + r": (.*)$", block)
    return m.group(1) if m else None

out = []
out.append("# Perelomov section (P-001 to P-049, P2-001 to P2-058), exact computation only\n")
out.append("Every check below is exact: QF's symbolic second-quantization algebra, exact symbolic Wolfram computation, exact spin, 2x2 group, coefficient or Jordan-Wigner fermion matrices, or a hand derivation. No truncated Fock matrix is used anywhere. Each `wl` cell runs as written after the loading lines\n")
out.append("```wl\nNeeds[\"Wolfram`QuantumFramework`\"];\nNeeds[\"Wolfram`QuantumFramework`SecondQuantization`\"];\n```\n")
out.append("and its Output is copied from that run (`drafts/exact/P-work/driver.wl` runs every cell in `cells/` in a cleared context; outputs in `out_cells.txt`). The exact checks without a QF route are in `drafts/exact/P-work/wlchecks.wl` (outputs in `out_wl.txt`).\n")
out.append("""## Errors found

No new error in the atlas transcription. Every source correction the atlas records holds under exact computation:
- P-029: the printed middle exponents of the SU(1,1) disentangling are reversed; the corrected form holds exactly in the faithful 2x2 matrices and in QF's single-mode realization, and the printed form fails in both.
- P-019: the factor 2 in the $P\\to W$ and $W\\to Q$ smoothing kernels is required (exact Gaussian integrals).
- P2-013: the corrected phase equation solves the Schrodinger equation exactly; the printed sign leaves a residual proportional to $k\\,\\mathrm{Im}(\\bar\\beta\\zeta)$.
- P2-017: the printed normalization is wrong already for real parameters; the corrected one gives norm 1 exactly for complex $\\zeta,\\alpha$ (exact Bargmann integral).
- P2-023: the trace-log expansion needs a plus sign (exact series).
- P2-054: the extra $e^{|\\gamma|^2}$ in the printed generating function is wrong (QF exact matrix elements).
- P2-038: the corrected sign in the Gaussian kernel holds (normalized Gaussian, exact integral).

Convention gaps the atlas leaves implicit (formulas are right only under a specific index reading):
- P-045: with $\\hat X_{ij}=\\hat c_i\\hat c_j$, the displayed $A$ equals $\\langle \\hat c_j\\hat c_i\\rangle=-\\langle\\hat c_i\\hat c_j\\rangle$, and $B_{ij}=\\langle\\hat c_j^\\dagger\\hat c_i\\rangle-\\delta_{ij}/2$. Read naively as $A_{ij}=\\langle\\hat X_{ij}\\rangle$, $A$ has the wrong sign.
- P-039: $A_{ij}=\\langle\\hat a_i\\hat a_j\\rangle$ holds, but $B_{ij}=\\langle\\hat a_j^\\dagger\\hat a_i\\rangle+\\delta_{ij}/2$, the transpose of $\\langle\\hat X^i{}_j\\rangle$ under the atlas's definition of $\\hat X^i{}_j$.
- P2-012: $\\hat K_2=\\tfrac14(\\hat p\\hat q+\\hat q\\hat p)$ equals $-(\\hat K_+-\\hat K_-)/(2i)$ for $\\hat K_\\pm=\\tfrac12\\hat a^{\\dagger2},\\tfrac12\\hat a^2$. This is a consistent su(1,1) basis only together with the source's $\\hat K_1$, which flips sign too (the relabeling $\\hat a\\to i\\hat a$).

QF defect seen while writing the cells: `BosonicNormalOrder` given a List of expressions returns pure functions for some elements instead of mapping over the list; for example `BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}]` returns `{{0 & }, SuperDagger[\\[FormalA]1]**\\[FormalA]1}`, while `BosonicNormalOrder[a ** ad, {a, ad}]` correctly returns `1 + SuperDagger[\\[FormalA]1]**\\[FormalA]1`. The cells below map it over lists explicitly.

---
""")

for i, blk in enumerate(blocks, start=1):
    e = E[i]
    title = blk.splitlines()[0].strip()
    topic, appears = field(blk, "Topic"), field(blk, "Appears in")
    formula = re.search(r"(?m)^\$\$.*\$\$$", blk).group(0)
    cond = field(blk, "Conditions")
    s = [f"### {title}", f"Topic: {topic}", f"Appears in: {appears}", "", e["desc"], "", formula, ""]
    if cond:
        s.append(f"Conditions: {cond}")
    s.append(f"Audit: {e['audit']}")
    if e["cell"] is True:
        code = (here / "cells" / f"{i:02d}.wl").read_text().rstrip("\n")
        o = SHORT.get(i, outs[i])
        s.append("QF code (exact):")
        s.append("```wl")
        s.append(code)
        s.append("```")
        s.append(f"Output: `{o}`. {e['match']}")
    else:
        s.append(f"QF code (exact): {e['cell']}")
    out.append("\n".join(s) + "\n")

out.append("## Non-formula records" + nonformula)
text = "\n".join(out)
assert "—" not in text and "–" not in text
(drafts / "exact" / "P-section.md").write_text(text)
print("entries:", len(blocks), "with QF cell:", sum(1 for k in E if E[k]["cell"] is True))
