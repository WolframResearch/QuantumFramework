#!/usr/bin/env python3
"""backbone_lint.py: the deterministic checker for backbone pages (key v3, nine rules).

Derives the dependency graph from a page's node headers and checks the mechanical
invariants; the authored metadata is diffed against the derived structure instead of
being trusted. Tier 1 findings (ERROR) are decidable graph facts; tier 2 (WARNING)
are mechanical-but-noisy; NOTICE covers declared deviations. --render emits an SVG
generated from the derived graph, so drawing and page cannot drift.

Usage:
  python3 backbone_lint.py PAGE.md [--json] [--render OUT.svg]
"""

import json
import re
import sys
from collections import OrderedDict

LABEL_RX = r"(?:[RCSEN]\d+|S0|[HEN])[\u2032\u2033]?"
HEADER_RX = re.compile(r"^\*\*(" + LABEL_RX + r"|[A-Z][^*(]{0,12})\s*\(([^)]*)\)\s*\.?\*\*", re.M)
MATH_RX = re.compile(r"\$\$([^$]+)\$\$|\$([^$]+)\$")


def math_spans(text):
    return [m.group(1) or m.group(2) for m in MATH_RX.finditer(text)]
ITALIC_RX = re.compile(r"(?<!\*)\*([^*\n][^*]*?)\*(?!\*)", re.S)

PROV_WORDS = ["derived in place", "cited", "conjectured"]
VERIF_WORDS = ["exact decision", "analytic cross-check", "numerical check",
               "verification: none", "verified", "verification:", "checked"]


def norm_label(raw):
    s = raw.replace("$", "").strip()
    s = s.replace("'''", "\u2033").replace("''", "\u2033").replace("'", "\u2032")
    return s


def label_tokens(text):
    toks = []
    for m in re.finditer(r"\b(" + LABEL_RX + r")\b", text.replace("$", "")):
        toks.append(norm_label(m.group(1)))
    return toks


class Node:
    def __init__(self, label, kindtext, body):
        self.label = label
        self.kindtext = kindtext.strip()
        self.body = body
        self.parents = []          # list of (label, weakened_flag)
        self.tail = None
        self.kind = self.classify()

    def classify(self):
        kt = self.kindtext.lower()
        if "\u2032" in self.label or "\u2033" in self.label:
            return "S"  # a primed node is a weakened statement, whatever its base
        if "non-assumption" in kt:
            return "N"
        if "assumption" in kt:
            return "R"
        if "choice" in kt:
            return "C"
        if "contract" in kt:
            return "E"
        if "cited" in kt and self.label.startswith("S0"):
            return "S0"
        if self.label.startswith("H"):
            return "H"
        if self.label.startswith("S"):
            return "S"
        if self.label.startswith("N"):
            return "N"
        if self.label.startswith("R"):
            return "R"
        if self.label.startswith("C"):
            return "C"
        if self.label.startswith("E"):
            return "E"
        return "?"


def parse_page(text):
    nodes = OrderedDict()
    preamble = text.split("\n## ", 1)[0]
    matches = list(HEADER_RX.finditer(text))
    problems = []
    for i, m in enumerate(matches):
        raw_label, kindtext = m.group(1), m.group(2)
        label = norm_label(raw_label)
        if not re.fullmatch(LABEL_RX, label):
            continue  # bold text that merely resembles a header
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        stop = text.find("\n## ", m.end())
        if stop != -1 and stop < end:
            end = stop
        body = text[m.end():end].strip()
        node = Node(label, kindtext, body)
        # parents: label tokens in the kind text, excluding self
        for tok in label_tokens(kindtext):
            if tok == label:
                continue
            weakened = bool(re.search(re.escape(tok.replace("\u2032", "'")) + r"\S*\s+weakened",
                                      kindtext.replace("\u2032", "'")))
            node.parents.append((tok, weakened))
        # a primed node implicitly derives from its base node
        base = label.rstrip("\u2032\u2033")
        if base != label and base not in [p for p, _ in node.parents]:
            node.parents.append((base, False))
        # tail: final italic span that closes the body
        spans = list(ITALIC_RX.finditer(body))
        if spans:
            last = spans[-1]
            if len(body) - last.end() <= 2:
                node.tail = last.group(1).strip()
        if label in nodes:
            problems.append(("ERROR", label, "duplicate node label"))
        nodes[label] = node
    return preamble, nodes, problems


def ancestors(nodes, label, seen=None):
    seen = set() if seen is None else seen
    for p, _w in nodes.get(label, Node(label, "", "")).parents:
        if p not in seen:
            seen.add(p)
            ancestors(nodes, p, seen)
    return seen


def declared_deviations(preamble):
    return set(int(n) for n in re.findall(r"[Dd]eclared deviation[^.]*rule\s*(\d+)", preamble))


def lint(path):
    text = open(path, encoding="utf-8").read()
    preamble, nodes, findings = parse_page(text)
    deviations = declared_deviations(preamble)
    if not nodes:
        findings.append(("ERROR", "-", "no parseable node headers found (is this a backbone page?)"))
        return findings, nodes

    # referential integrity and acyclicity
    for n in nodes.values():
        for p, _w in n.parents:
            if p not in nodes:
                findings.append(("ERROR", n.label, f"parent {p} does not exist on the page"))
    state = {}

    def cyclic(lbl, stack):
        state[lbl] = 1
        for p, _w in nodes[lbl].parents:
            if p in nodes:
                if state.get(p) == 1 or (state.get(p) is None and cyclic(p, stack)):
                    return True
        state[lbl] = 2
        return False
    for lbl in nodes:
        if state.get(lbl) is None and cyclic(lbl, []):
            findings.append(("ERROR", lbl, "cycle in the derives-from relation"))
            break

    # label-prefix versus kind consistency (primed variants are weakened statements, exempt)
    for n in nodes.values():
        if "\u2032" in n.label or "\u2033" in n.label:
            continue
        prefix = re.match(r"[A-Z]+", n.label).group(0)[0]
        expect = {"R": "R", "C": "C", "H": "H", "N": "N", "E": "E"}.get(prefix)
        if n.label == "S0" and n.kind not in ("S0", "E"):
            findings.append(("ERROR", n.label, "S0 must be marked cited (or E with a scope statement)"))
        elif expect and n.kind != expect and not n.label.startswith("S"):
            findings.append(("ERROR", n.label, f"label says {expect}-kind but header reads '{n.kindtext}'"))

    # rule 2: every root consumed by a non-weakening statement
    for n in nodes.values():
        if n.kind != "R":
            continue
        ok = any(any(p == n.label and not w for p, w in m.parents)
                 for m in nodes.values() if m.label != n.label)
        if ok:
            continue
        sev = "NOTICE" if 2 in deviations else "ERROR"
        findings.append((sev, n.label, "rule 2: no statement consumes this root except its own weakening"))

    # tails: presence, provenance, verification, naked standard, cited sources in parents
    for n in nodes.values():
        if n.kind in ("S", "S0", "H", "E"):
            if not n.tail:
                findings.append(("ERROR", n.label, "statement has no italic tail"))
                continue
            low = n.tail.lower()
            if not any(w in low for w in PROV_WORDS):
                findings.append(("ERROR", n.label, "tail names no provenance (derived in place / cited / conjectured)"))
            if not any(w in low for w in VERIF_WORDS):
                findings.append(("ERROR", n.label, "tail names no verification (exact decision / analytic cross-check / numerical check / none)"))
            if re.search(r"\bstandard\b", low) and "cited" not in low:
                findings.append(("ERROR", n.label, "naked 'standard' in tail: pin a source"))
            prov_half = re.split(r"[Vv]erif", n.tail)[0]
            cited_here = set(label_tokens(prov_half)) - {n.label}
            allowed = set(p for p, _w in n.parents)
            for tok in cited_here:
                if tok not in allowed:
                    findings.append(("ERROR", n.label,
                                     f"tail derives from {tok}, which is not among the header's parents"))

    # over-listed parents: the tail's provenance half names sources but omits this parent,
    # and the body never mentions it either (the B2/B3 signature)
    for n in nodes.values():
        if not n.tail:
            continue
        prov_half = re.split(r"[Vv]erif", n.tail)[0]
        prov_named = set(label_tokens(prov_half))
        if not prov_named:
            continue
        mentioned = set(label_tokens(n.body))
        for p, w in n.parents:
            if w or p in prov_named or p in mentioned:
                continue
            findings.append(("WARNING", n.label,
                             f"parent {p} is listed but neither the derivation tail nor the body mentions it; padded from-list?"))

    # audit line at H: recompute closure, compare
    heads = [n for n in nodes.values() if n.kind == "H"]
    if not heads:
        findings.append(("ERROR", "-", "no headline node found"))
    for h in heads:
        m = re.search(r"Audit:(.+?)(?:\*|$)", h.body, re.S)
        if not m:
            findings.append(("ERROR", h.label, "rule 8: headline has no audit line"))
            continue
        audit_toks = set(label_tokens(m.group(1)))
        closure = ancestors(nodes, h.label)
        must_name = {c for c in closure if nodes.get(c) and nodes[c].kind in ("R", "C", "S0", "E")}
        for c in sorted(must_name - audit_toks):
            findings.append(("ERROR", h.label, f"audit line omits {c}, which is in the closure"))
        for c in sorted(audit_toks - closure - {h.label}):
            if c in nodes:
                findings.append(("ERROR", h.label, f"audit line names {c}, which is not in the closure"))
        conj = [c for c in closure | {h.label}
                if nodes.get(c) and nodes[c].tail and "conjectured" in nodes[c].tail.lower()]
        says_none = re.search(r"no conjecture", m.group(1), re.I)
        if conj and says_none:
            findings.append(("ERROR", h.label, f"audit claims no conjectures but {', '.join(conj)} carries one"))
        if conj and not says_none and not re.search(r"conjectur", m.group(1), re.I):
            findings.append(("ERROR", h.label, f"conjecture in closure ({', '.join(conj)}) not surfaced in audit"))

    # symbol closure (tier 2): defined with '=' in one node, used elsewhere without an edge path
    FORMAT_CMDS = {"\\frac", "\\tfrac", "\\sum", "\\exp", "\\max", "\\min", "\\sqrt",
                   "\\hat", "\\bar", "\\dot", "\\vec", "\\left", "\\right", "\\mathcal",
                   "\\mathrm", "\\langle", "\\rangle", "\\otimes", "\\quad", "\\qquad",
                   "\\text", "\\dagger", "\\infty", "\\prime", "\\log", "\\ln"}
    LETTER_NOISE = {"e", "t", "i", "m", "s", "T"}

    def symbols_defined(n):
        out = set()
        for span in math_spans(n.body):
            for sm in re.finditer(r"(?:^|[\s(,;{])(\\[A-Za-z]+|[A-Za-z\u2113])(?:_\{?[A-Za-z0-9]+\}?)?\s*=\s*([^,;$]*)", span):
                sym, rhs = sm.group(1), sm.group(2).strip()
                if sym in FORMAT_CMDS or sym in LETTER_NOISE:
                    continue
                first = rhs.split()[0] if rhs.split() else ""
                if re.fullmatch(r"[-+]?\d+(?:\.\d+)?\s*\.?", rhs) or first == "0":
                    continue  # a numeric pinning (eta = 0.7) or set-to-zero, not a definition
                out.add(sym)
        return out

    def symbols_used(n):
        out = set()
        for span in math_spans(n.body):
            for sm in re.finditer(r"\\[A-Za-z]+|[A-Za-z\u2113]", span):
                tok = sm.group(0)
                if tok in FORMAT_CMDS or tok in LETTER_NOISE:
                    continue
                out.add(tok)
        return out

    definer = {}
    for n in nodes.values():
        for s in symbols_defined(n):
            definer.setdefault(s, n.label)
    usage_count = {s: sum(1 for n in nodes.values() if s in symbols_used(n)) for s in definer}
    for n in nodes.values():
        if n.kind not in ("S", "S0", "H", "E"):
            continue  # roots, choices, and non-assumptions share the notation namespace
        anc = ancestors(nodes, n.label)
        for s in symbols_used(n):
            d = definer.get(s)
            if not d or d == n.label:
                continue
            if usage_count.get(s, 0) > max(3, len(nodes) // 2):
                continue  # background notation
            if d not in anc and n.label not in ancestors(nodes, d):
                findings.append(("WARNING", n.label,
                                 f"uses {s} defined in {d}, but no edge path connects them (missing from-list entry?)"))
    return findings, nodes


# ---------------- SVG generation from the derived graph ----------------

KIND_STYLE = {
    "R": ("#E7EEF9", "#5B7DB1", 1.5, None), "C": ("#FDF2DA", "#C09A4A", 1.5, None),
    "S0": ("#E9F4EA", "#5E9C6B", 1.5, "double"), "E": ("#E9F4EA", "#5E9C6B", 1.5, "double"),
    "S": ("#F4F4F4", "#909090", 1.5, None), "H": ("#FFFFFF", "#333333", 2.5, None),
    "N": ("#FAFAFA", "#999999", 1.5, "dotted"), "?": ("#FFFFFF", "#CC0000", 1.5, None),
}


def gloss_snippet(n, limit=80):
    plain = n.body.split("$")[0]
    plain = re.sub(r"\*+", "", plain).strip()
    if len(plain) <= limit:
        return plain.rstrip(",;: ")
    cut = plain[:limit]
    # prefer a clause boundary, else a word boundary, then mark the truncation
    m = max(cut.rfind(", "), cut.rfind("; "), cut.rfind(": "))
    if m > limit // 2:
        return cut[:m].rstrip() + " \u2026"
    return cut[:cut.rfind(" ")].rstrip(",;: ") + " \u2026"


def render_svg(nodes, title, out_path):
    layer = {}

    def depth(lbl, seen=()):
        if lbl in layer:
            return layer[lbl]
        ps = [p for p, _w in nodes[lbl].parents if p in nodes and p not in seen]
        d = 0 if not ps else 1 + max(depth(p, seen + (lbl,)) for p in ps)
        layer[lbl] = d
        return d
    for lbl in nodes:
        depth(lbl)
    maxd = max(layer.values()) if layer else 0
    rows = {}
    for lbl, d in layer.items():
        rows.setdefault(d, []).append(lbl)
    W = max(980, 200 * max(len(r) for r in rows.values()) + 60)
    H_px = 170 + 120 * (maxd + 1) + 60
    bw, bh = 180, 58
    pos = {}
    for d in range(maxd + 1):
        row = rows.get(d, [])
        for i, lbl in enumerate(sorted(row)):
            x = (i + 1) * W / (len(row) + 1)
            pos[lbl] = (x, 90 + 120 * d)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H_px}" viewBox="0 0 {W} {H_px}" font-family="Helvetica, Arial, sans-serif">',
             f'<rect width="{W}" height="{H_px}" fill="#ffffff"/>',
             f'<text x="{W/2}" y="30" text-anchor="middle" font-size="15" font-weight="bold" fill="#222">{title}</text>',
             f'<text x="{W/2}" y="48" text-anchor="middle" font-size="10.5" fill="#777">generated by backbone_lint.py from the page\u2019s from-lists</text>',
             '<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#666"/></marker>'
             '<marker id="arrd" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#A6699C"/></marker></defs>']
    for n in nodes.values():
        x2, y2 = pos[n.label]
        for p, w in n.parents:
            if p not in pos:
                continue
            x1, y1 = pos[p]
            dash = ' stroke-dasharray="6 4"' if w else ""
            color, mark = ("#A6699C", "arrd") if w else ("#666666", "arr")
            parts.append(f'<line x1="{x1:.0f}" y1="{y1 + bh/2:.0f}" x2="{x2:.0f}" y2="{y2 - bh/2:.0f}" '
                         f'stroke="{color}" stroke-width="1.6"{dash} marker-end="url(#{mark})"/>')
    for n in nodes.values():
        x, y = pos[n.label]
        fill, stroke, sw, extra = KIND_STYLE.get(n.kind, KIND_STYLE["?"])
        dash = ' stroke-dasharray="4 3"' if extra == "dotted" else ""
        parts.append(f'<rect x="{x - bw/2:.0f}" y="{y - bh/2:.0f}" width="{bw}" height="{bh}" rx="8" '
                     f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{dash}/>')
        if extra == "double":
            parts.append(f'<rect x="{x - bw/2 + 4:.0f}" y="{y - bh/2 + 4:.0f}" width="{bw - 8}" height="{bh - 8}" rx="6" fill="none" stroke="{stroke}" stroke-width="0.8"/>')
        kindname = {"R": "assumption", "C": "choice", "S0": "equation, cited", "E": "contract",
                    "S": "statement", "H": "headline", "N": "non-assumption"}.get(n.kind, "?")
        parts.append(f'<text x="{x:.0f}" y="{y - 12:.0f}" text-anchor="middle" font-size="12" font-weight="bold" fill="#222">{n.label} \u00b7 {kindname}</text>')
        snip = gloss_snippet(n)
        if snip:
            half = len(snip) // 2
            cut = snip.rfind(" ", 0, half + 12)
            l1, l2 = (snip[:cut], snip[cut + 1:]) if cut > 0 else (snip, "")
            parts.append(f'<text x="{x:.0f}" y="{y + 3:.0f}" text-anchor="middle" font-size="9.5" fill="#333">{l1}</text>')
            if l2:
                parts.append(f'<text x="{x:.0f}" y="{y + 16:.0f}" text-anchor="middle" font-size="9.5" fill="#333">{l2}</text>')
        if n.kind == "H":
            anc = ancestors(nodes, n.label)
            named = sorted(a for a in anc if a in nodes and nodes[a].kind in ("R", "C", "S0", "E"))
            conj = any(nodes[a].tail and "conjectured" in nodes[a].tail.lower()
                       for a in anc | {n.label} if a in nodes)
            caption = "audit: rests on " + " ".join(named) + (" \u00b7 conjectures in closure" if conj else " \u00b7 no conjectures")
            parts.append(f'<text x="{x:.0f}" y="{y + bh/2 + 16:.0f}" text-anchor="middle" font-size="9" font-style="italic" fill="#888">{caption}</text>')
    parts.append("</svg>")
    open(out_path, "w", encoding="utf-8").write("\n".join(parts))


def main():
    args = sys.argv[1:]
    as_json = "--json" in args
    render_to = None
    if "--render" in args:
        render_to = args[args.index("--render") + 1]
    paths = [a for a in args if a.endswith(".md")]
    if not paths:
        print(__doc__)
        sys.exit(2)
    worst = 0
    for path in paths:
        findings, nodes = lint(path)
        if render_to and nodes:
            title = "Dependency graph \u00b7 " + (re.search(r"^# (.+)$", open(path, encoding="utf-8").read(), re.M) or [None, path]).group(1)
            render_svg(nodes, title, render_to)
        if as_json:
            print(json.dumps({"page": path, "findings": [
                {"severity": s, "node": n, "message": msg} for s, n, msg in findings]}, ensure_ascii=False, indent=1))
        else:
            print(f"== {path}: {len(nodes)} nodes parsed")
            for s, n, msg in findings:
                print(f"  {s:7s} [{n}] {msg}")
            if not findings:
                print("  clean")
        if any(s == "ERROR" for s, _n, _m in findings):
            worst = 1
    sys.exit(worst)


if __name__ == "__main__":
    main()
