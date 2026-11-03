import json, re
d = json.load(open("exact-merge/merged.json"))
secs = d["sections"]; c = d["counts"]
def slug(s): return re.sub(r"[^a-z0-9]+","-",s.lower()).strip("-")
unchecked = []
for t, ents in secs.items():
    for e in ents:
        if re.search(r"Audit:\s*Unchecked", e):
            unchecked.append(re.search(r"^### (.*)$", e, re.M).group(1))
head = open("exact-merge/head.md").read()
head = head.replace("{N}", str(c["entries"])).replace("{QF}", str(c["qf"])).replace("{NF}", str(len(d["nf"])))
parts = [head]
parts.append("## Contents\n\n" + "\n".join(f"- [{t}](#{slug(t)}) ({len(v)})" for t, v in secs.items()) + "\n- [Non-formula records](#non-formula-records)\n- [Index of atlas IDs](#index-of-atlas-ids)\n")
parts.append("## Entries marked Unchecked\n\n" + "\n".join(f"- {u}" for u in unchecked) + "\n")
for t, v in secs.items():
    parts.append(f'<a id="{slug(t)}"></a>\n## {t}\n\n' + "\n\n".join(v) + "\n")
parts.append('<a id="non-formula-records"></a>\n## Non-formula records\n\nThese atlas records state domains, theorems in prose, or coverage notes and carry no formula to compute.\n\n' + "\n".join(d["nf"]) + "\n")
def key(i):
    m = re.match(r"([XBPCG])(2?)-(\d+)", i); return ("XBPCG".index(m.group(1)), m.group(2), int(m.group(3)))
idx = []
for i in sorted(d["idmap"], key=key):
    seen = []
    for title, a in d["idmap"][i]:
        if (title, a) not in seen: seen.append((title, a))
    idx.append(f"- {i}: " + "; ".join(f"[{t}](#{a})" for t, a in seen))
parts.append('<a id="index-of-atlas-ids"></a>\n## Index of atlas IDs\n\n' + "\n".join(idx) + "\n")
doc = "\n".join(parts)
doc = doc.replace(" — ", ": ").replace("—", ", ").replace(" – ", " to ")
open("../Operator-Exponentials-QF-Audit.md", "w").write(doc)
print(len(doc), len(unchecked), doc.count("—"), doc.count("–"))
