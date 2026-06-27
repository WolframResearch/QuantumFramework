#!/usr/bin/env python3
"""Extract the Colab-specific override pieces from an existing introduction-to-QIS-colab.ipynb
so build_colab.py can regenerate the notebook from introduction-to-QIS-python.md.

Overrides captured (the parts that are NOT a pure function of python.md):
  - badge_cell  : the markdown cell holding the "Open in Colab" badge + how-to-run intro.
  - title_header: the title/subtitle/author/affiliation/hero-image markdown that replaces
                  python.md's md2nb front matter (everything in the title cell BEFORE the
                  "### Setting the Stage" prose, which itself comes from python.md).
  - setup_cell  : the free-Wolfram-Engine setup code cell (install/activate/kernel/helpers).
  - animations  : {fig-NNN: wlgif-source} for the 13 cells where a wlshow(Manipulate[...])
                  in python.md becomes an animated wlgif(Table[...]) in the notebook.

The same PlotLabel sigma/omega Subscript -> SubscriptBox typesetting fix that is applied to
python.md is applied here to the animation sources, so the overrides stay consistent.
"""
import json, re, sys

SRC_NB = sys.argv[1] if len(sys.argv) > 1 else "introduction-to-QIS-colab.ipynb"
OUT    = sys.argv[2] if len(sys.argv) > 2 else "colab_overrides.json"

SUBSCRIPT = re.compile(r'Subscript\[\\\[(Sigma|Omega)\],\s*([xyz])\]')
def plotlabel_fix(s):
    return SUBSCRIPT.sub(lambda m: r'\!\(\*SubscriptBox[\(\[%s]\), \(%s\)]\)' % (m.group(1), m.group(2)), s)

nb = json.load(open(SRC_NB))
def src(c): return "".join(c["source"])
cells = nb["cells"]

badge_cell  = src(cells[0])
setup_cell  = src(cells[2])

# title cell = header + python.md intro prose; keep only the header part.
title_full = src(cells[1])
marker = "### Setting the Stage"
title_header = title_full[:title_full.index(marker)].rstrip() + "\n" if marker in title_full else title_full

animations = {}
for c in cells:
    if c["cell_type"] != "code":
        continue
    m = re.search(r'wlgif\(r%s.*?%s,\s*"(fig-\d+)"\)' % ("'''", "'''"), src(c), re.S)
    if m:
        animations[m.group(1)] = plotlabel_fix(src(c))

json.dump({"badge_cell": badge_cell, "title_header": title_header,
           "setup_cell": setup_cell, "animations": animations},
          open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"wrote {OUT}: badge({len(badge_cell)}c) title_header({len(title_header)}c) "
      f"setup({len(setup_cell)}c) animations({len(animations)})")
