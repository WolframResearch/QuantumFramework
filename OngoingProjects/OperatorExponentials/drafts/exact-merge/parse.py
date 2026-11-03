import re, json, sys
secs = ["X","B","P","C","G"]
entries = []
for s in secs:
    txt = open(f"exact/{s}-section.md").read()
    body = txt.split("\n## Non-formula records")[0]
    parts = re.split(r"\n(?=### )", body)
    for i,p in enumerate(parts):
        if not p.startswith("### "): continue
        title = p.split("\n",1)[0][4:].strip()
        topic = re.search(r"^Topic:\s*(.*)$", p, re.M)
        app = re.search(r"^Appears in:\s*(.*)$", p, re.M)
        entries.append(dict(key=f"{s}{i:03d}", sec=s, title=title, topic=topic.group(1).strip() if topic else "?", appears=app.group(1).strip() if app else "?", text=p.strip()))
json.dump(entries, open("exact-merge/entries.json","w"))
from collections import defaultdict
g = defaultdict(list)
for e in entries: g[e["topic"]].append(e)
for t in sorted(g):
    print(f"\n## {t} ({len(g[t])})")
    for e in g[t]:
        print(f"{e['key']} | {e['title']} | {e['appears'][:70]}")
