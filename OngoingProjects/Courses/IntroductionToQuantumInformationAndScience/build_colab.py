#!/usr/bin/env python3
"""Build introduction-to-QIS-colab.ipynb from introduction-to-QIS-python.md (self-contained).

python.md is the source of truth: prose, Wolfram-Language code, the text results (embedded as
"*Output:*" fenced blocks) and the figures (embedded as ![fig-NNN](qis-figs/fig-NNN.png) refs,
with the actual images in qis-figs/). The notebook adds Colab-specific pieces that are not a
pure function of python.md and live in colab_overrides.json (see extract_overrides.py):
  - the Open-in-Colab badge cell,
  - a cleaned title header (replacing python.md's md2nb front matter),
  - the free-Wolfram-Engine setup cell,
  - 13 wlshow(Manipulate[...]) cells rendered as animated wlgif(Table[...]) (fig-NNN.gif).

Output cells: text results come straight from the "*Output:*" blocks (-> stdout stream); figures
are read from qis-figs/ and embedded as base64 (image/png, or image/gif for the 13 animations).
No Wolfram kernel is needed to build the notebook.
"""
import argparse, base64, json, os, re

HERE   = os.path.dirname(os.path.abspath(__file__))
MD     = os.path.join(HERE, "introduction-to-QIS-python.md")
OVR    = os.path.join(HERE, "colab_overrides.json")
FIGDIR = os.path.join(HERE, "qis-figs")

# a ```python``` code block, optionally followed by its "*Output:*" text block or an image ref
UNIT = re.compile(
    r'```python[ \t]*\n(?P<code>.*?)\n```'
    r'(?:\s*\*Output:\*\s*\n+```[a-z]*\n(?P<txt>.*?)\n```'
    r'|\s*!\[[^\]]*\]\(qis-figs/fig-\d+\.(?:png|gif)\)(?:\n\*fig-\d+[^*\n]*\*)?)?',
    re.S)

def md_cell(text):   return {"cell_type":"markdown","metadata":{},"source":_lines(text)}
def code_cell(text, outs):
    return {"cell_type":"code","metadata":{},"execution_count":None,"outputs":outs,"source":_lines(text)}
def _lines(t): return (t.rstrip("\n")).splitlines(keepends=True) or [t.rstrip("\n")]

def stream_out(text):
    return {"output_type":"stream","name":"stdout","text":(text+"\n").splitlines(keepends=True)}
def image_out(path, mime):
    b64 = base64.b64encode(open(path,"rb").read()).decode()
    return {"output_type":"display_data","metadata":{},"data":{mime:b64}}

def build():
    ovr   = json.load(open(OVR))
    anims = ovr["animations"]
    md    = open(MD, encoding="utf-8").read()

    intro = md[md.index("### Setting the Stage"): md.index("Let’s start!")+len("Let’s start!")].rstrip()+"\n"
    body  = md[md.index("## Part I"):]

    cells = [md_cell(ovr["badge_cell"]),
             md_cell(ovr["title_header"].rstrip()+"\n\n"+intro),
             code_cell(ovr["setup_cell"], [])]

    missing = []
    pos = 0
    for m in UNIT.finditer(body):
        prose = body[pos:m.start()].strip("\n")
        if prose.strip():
            cells.append(md_cell(prose))
        pos = m.end()

        code = m.group("code")
        figm = re.search(r'"(fig-\d+)"', code)
        fig  = figm.group(1) if figm else None
        is_anim = bool(fig) and fig in anims
        if is_anim:
            code = anims[fig]

        outs = []
        if m.group("txt") is not None:
            outs.append(stream_out(m.group("txt")))
        elif fig:
            ext  = "gif" if is_anim else "png"
            mime = "image/gif" if is_anim else "image/png"
            p = os.path.join(FIGDIR, f"{fig}.{ext}")
            if os.path.exists(p):
                outs.append(image_out(p, mime))
            else:
                missing.append(f"{fig}.{ext}")
        cells.append(code_cell(code, outs))

    tail = body[pos:].strip("\n")
    if tail.strip():
        cells.append(md_cell(tail))

    if missing:
        print("[warn] missing figure files:", missing)
    nstream = sum(1 for c in cells if c["cell_type"]=="code" and any(o["output_type"]=="stream" for o in c["outputs"]))
    nimg    = sum(1 for c in cells if c["cell_type"]=="code" and any(o["output_type"]=="display_data" for o in c["outputs"]))
    print(f"[build] {len(cells)} cells | {nstream} text outputs | {nimg} image outputs")
    return {"cells": cells,
            "metadata": {"language_info": {"name": "python"},
                         "kernelspec": {"name": "python3", "display_name": "Python 3"}},
            "nbformat": 4, "nbformat_minor": 5}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o","--out", default=os.path.join(HERE,"introduction-to-QIS-colab.ipynb"))
    a = ap.parse_args()
    json.dump(build(), open(a.out,"w"), indent=1, ensure_ascii=False)
    print("wrote", a.out)
