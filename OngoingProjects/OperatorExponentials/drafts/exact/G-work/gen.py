import re, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entries import E
base = __import__("os").path.normpath(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", ".."))
src = open(f"{base}/G-section.md").read()
work = f"{base}/exact/G-work"
# outputs
outs = {}
cur = None
for line in open(f"{work}/out-cells.txt").read().split("\n"):
    m = re.match(r"^#### (\d{3}) .*msgs=(.*)$", line)
    if m:
        cur = int(m.group(1)); assert m.group(2) == "{}", (cur, m.group(2)); outs[cur] = []
    elif cur is not None and line != "":
        outs[cur].append(line)
outs = {k: " ".join(v) for k, v in outs.items()}
body, rest = src.split("## Non-formula records", 1)
nonformula = "## Non-formula records" + rest.split("## QF notes", 1)[0]
chunks = re.split(r"(?m)^### ", body)[1:]
assert len(chunks) == 103, len(chunks)
def code(s):
    return f"``` {s} ```" if "`" in s else f"`{s}`"
parts = []
for i, ch in enumerate(chunks, 1):
    lines = ch.rstrip("\n").split("\n")
    title, topic, appears = lines[0], lines[1], lines[2]
    assert topic.startswith("Topic:") and appears.startswith("Appears in:")
    formula = [l for l in lines if l.startswith("$$")]
    assert len(formula) == 1, i
    cond = [l for l in lines if l.startswith("Conditions:")]
    d = E[i]
    out = [f"### {title}", topic, appears, "", d["desc"], "", formula[0], ""]
    out += cond
    out.append(f"Audit: {d['audit']}. Checked by: {d['checked']}.")
    cell = f"{work}/cells/{i:03d}.wl"
    if os.path.exists(cell):
        assert d["none"] is None, i
        out.append("QF code (exact):")
        out.append("```wl")
        out.append(open(cell).read().rstrip("\n"))
        out.append("```")
        out.append(f"Output: {code(outs[i])}. {d['match']}")
    else:
        assert d["none"] is not None, i
        out.append(f"QF code (exact): none. {d['none']}")
    parts.append("\n".join(out))
head = """# Section G: Gerry and Knight, Introductory Quantum Optics (atlas lines 4361 to 5648)

## Errors found

No error in the atlas formulas of this section; every entry was rechecked with exact algebra only. The exact checks confirm the atlas's corrections of three printed equations: the antinormal characteristic function of Eq. (3.128c) (a creation-first product is the normal-ordered function, which differs from the antinormal one by e^{|lambda|^2}; G-028), the odd-number phase of Eq. (7.142) (it is -i, not +i; G-077), and the missing number shift in a e^{-2iKt n} (off by e^{2iKt}; G-079). One remark stays unchecked: the interferometer fringe sign in G-050 depends on the book's Eq. (11.2) layout, which the atlas does not reproduce.

Every `wl` cell below runs as written after these loading lines, with `{a, ad} = FieldVariables[{1}]` (or `FieldVariables[{1, 2}]`, `FieldVariables[{1, 2, 3}]` for several modes) at its top:

```wl
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

The cells are the files `drafts/exact/G-work/cells/NNN.wl` (NNN = entry position), run by `drafts/exact/G-work/driver.wl` through `runner.wls`; their raw output is `drafts/exact/G-work/out-cells.txt`, and every cell ran without a message. Exact WL checks for entries without a QF cell are in `drafts/exact/G-work/wl/NNN.wl`, output `out-wl.txt`. No truncated Fock matrix is used anywhere. Symbols such as `al`, `be`, `xi`, `r`, `th` are left symbolic; `r > 0` and real `th` are passed as assumptions where a cell simplifies."""
tail = """## QF notes from this section

- No wrong value from an exact QF call. `BosonicNormalOrder`, `BosonicExpOrder` (one- and two-mode, normal and antinormal), `BosonicExpSimplify` (one-, two- and three-mode conjugations, products of displacements, conjugation of an exponential), `BosonicExpBraid`, `BosonicVEV`, `BosonicMatrixElement` (symbolic and integer Fock indices, one and several modes, coherent basis), and the symbolic `DisplacementOperator`/`SqueezeOperator` all agree with the atlas.
- Missing features met here, all returning the input unevaluated rather than a wrong value: `BosonicMatrixElement` with an exponential of `ad ** ad` in the coherent basis; `BosonicMatrixElement[{m, m}, DisplacementOperator[al]]` and `BosonicMatrixElement[{2 m, 0}, SqueezeOperator[xi]]` with symbolic m (integer m works); exponentials of `(ad ** a) ** (ad ** a)` (Kerr); `BosonicVEV` of a product of exponentials (route through `BosonicExpOrder` and `BosonicMatrixElement[{0, 0}, ...]`); `BosonicExpSimplify` of the Kerr conjugation, whose commutator series does not close.
- `BosonicExpSimplify[SqueezeOperator[xi, Infinity] ** SqueezeOperator[-xi, Infinity]]` returns `E^(X/2 + (-X)/2)` with the two halves not combined; expanding the exponent gives 1. Not a wrong value.
- `BosonicNormalOrder` needs the variable list as its second argument; with one argument it returns unevaluated.
"""
doc = head + "\n\n" + "\n\n".join(parts) + "\n\n" + nonformula.rstrip("\n") + "\n\n" + tail
assert "—" not in doc and "–" not in doc
open(f"{base}/exact/G-section.md", "w").write(doc)
print(len(parts), sum(1 for i in range(1, 104) if os.path.exists(f"{work}/cells/{i:03d}.wl")))
