import json, re
from collections import OrderedDict
E = json.load(open("exact-merge/entries.json"))
byk = {e["key"]: e for e in E}
groups = [["P002","C002","G003"],["X001","P011","C005","G016"],["X002","P009","G019"],["X003","P003","C003","G018"],
 ["P013","C006","G014"],["C007","C009","G020"],["P014","C010","G024"],["P006","G005","B002"],
 ["X004","P040","G054","C045"],["X005","G052"],["X006","G056"],["P039","G053"],
 ["X007","P042","G069"],["X008","G071"],["P041","G068"],
 ["X009","P030","C038"],["X010","P035"],["P029","C035"],["P034","C042"],["P032","P033","C039"],["P036","P037","C044"],
 ["X016","B071","G094"],["X017","B005"],["X018","B070"],
 ["P026","G029"],["P025","P017","C016"],["P024","G034"],["P027","G030"],
 ["X014","P046","B031"],["X015","B041"],["X011","P074","C082"],["X019","P043"],["P050","B039"],["P061","G023"]]
secondary = {k: g[0] for g in groups for k in g[1:]}
gmap = {g[0]: g for g in groups}
topics = ["BCH and Hadamard","Displacement","Coherent states (general)","Single-mode squeezing","Two-mode squeezing",
 "SU(2) and passive mode mixing","Bosonic canonical transformations","Fermionic canonical transformations",
 "Gaussian vacua and overlaps","Quadratic propagators and Riccati","Forced and parametric oscillators",
 "Normal symbols and normal-order exponentials","Weyl and antinormal symbols, characteristic functions",
 "Grassmann calculus and fermionic coherent states","Gaussian integrals and Pfaffians","Thermal states and traces",
 "Wick, Dyson and auxiliary fields","Nonlinear exponentials (number shift, Kerr, cats)",
 "Fock matrix elements and generating functions","Infinite modes and domains",
 "Lie-group decompositions and generalized coherent states"]
def topic_of(e):
    t = e["topic"]
    return "Bosonic canonical transformations" if t.startswith("Bosonic canonical") else t
def body(e):
    lines = e["text"].split("\n")[1:]
    return "\n".join(l for l in lines if not l.startswith("Topic:") and not l.startswith("Appears in:")).strip()
def hasqf(e):
    return bool(re.search(r"^QF code \(exact\):\s*$", e["text"], re.M))
def slug(s): return re.sub(r"[^a-z0-9]+","-",s.lower()).strip("-")
out = {t: [] for t in topics}; counts = {"entries":0,"qf":0}
idmap = {}
for e in E:
    if e["key"] in secondary: continue
    members = [byk[k] for k in gmap.get(e["key"], [e["key"]])]
    title = e["title"]; anchor = slug(title)
    apps = "; ".join(m["appears"] for m in members)
    txt = [f'<a id="{anchor}"></a>', f"### {title}", f"Appears in: {apps}", "", body(e)]
    for m in members[1:]:
        txt += ["", f"**Also stated in {m['appears'].split(' ')[0].rstrip(',;')} ({m['title']}):**", "", body(m)]
    out[topic_of(e)].append("\n".join(txt))
    counts["entries"] += 1; counts["qf"] += any(hasqf(m) for m in members)
    for i in re.findall(r"\b([XBPCG]2?-\d{3})\b", apps): idmap.setdefault(i, []).append((title, anchor))
nf = []
for s in "XBPCG":
    t = open(f"exact/{s}-section.md").read()
    sec = re.search(r"## Non-formula records\n(.*?)(\n## |\Z)", t, re.S).group(1).strip()
    for l in sec.split("\n"):
        if l.startswith("- "):
            nf.append(l)
            for i in re.findall(r"\b([XBPCG]2?-\d{3})\b", l.split(":")[0]): idmap.setdefault(i, []).append(("Non-formula records", "non-formula-records"))
json.dump({"counts":counts,"idmap":idmap,"nf":nf,"sections":{t:out[t] for t in topics}}, open("exact-merge/merged.json","w"))
missing = [t for t in out if t not in topics]
print(counts, len(nf), sum(len(v) for v in out.values()))
print({t: len(out[t]) for t in topics})
