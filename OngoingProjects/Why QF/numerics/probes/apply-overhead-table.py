"""Summarizes apply-overhead.tsv: the median of the repeats for each row and environment, the ratio,
and whether every run gave the same result. Usage: python3 -I apply-overhead-table.py [apply-overhead.tsv]"""
import sys, statistics, collections, pathlib
path = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).with_name("apply-overhead.tsv"))
runs = collections.defaultdict(list)
for line in path.read_text().splitlines():
    if not line or line.startswith("#"):
        continue
    rep, env, row, *fields = line.split("\t")
    runs[(row, env)].append(fields)
def num(x):
    try: return float(x)
    except ValueError: return None
two = {"apply-oracle10-zero", "apply-oracle10-uniform", "apply-qft14-zero", "apply-qft14-random"}
label = {
    "first-state": "A kernel's first QuantumState[\"0000000000\"] (ms)",
    "oracle10-new-states": "Phase oracle, 10 variables, applied to a new random state (s, median of five)",
    "apply-oracle10-zero": "Phase oracle, 10 variables, applied to the zero state (s, first / second)",
    "apply-oracle10-uniform": "Phase oracle, 10 variables, applied to the uniform state (s, first / second)",
    "apply-qft14-zero": "QFT, 14 qubits, applied to the zero state (s, first / second)",
    "apply-qft14-random": "QFT, 14 qubits, applied to a random state (s, first / second)",
    "quench-12": "Trotter quench, 4 steps, 12 qubits, in a fresh kernel (s)",
    "quench-16": "Trotter quench, 4 steps, 16 qubits, in a fresh kernel (s)",
    "quench-20": "Trotter quench, 4 steps, 20 qubits, in a fresh kernel (s)",
    "grover8-greedy": "Grover search of QF-Showcase, 8 qubits, \"Path\" -> \"Greedy\" (s)",
    "grover8-default": "Grover search of QF-Showcase, 8 qubits, default path (s)",
}
def fmt(m): return f"{m:.3g}" if m < 10 else f"{m:.1f}"
print("| Workload | Before | After | Before / after | Same result |")
print("|---|---|---|---|---|")
for row in label:
    cells, ratios, checks = [], [], set()
    for i in ([0, 1] if row in two else [0]):
        meds = []
        for env in ("before", "after"):
            vals = [num(f[i]) for f in runs.get((row, env), []) if len(f) > i]
            meds.append(statistics.median(vals) if vals and None not in vals else None)
        cells.append(meds)
        ratios.append(f"{meds[0] / meds[1]:.1f}x" if None not in meds else "")
    for env in ("before", "after"):
        for f in runs.get((row, env), []):
            if row in two and len(f) > 2: checks.add(f[2])
            elif row.startswith("quench") and len(f) > 1: checks.add(f[1])
            elif row == "oracle10-new-states" and len(f) > 1: checks.add(f[1])
            elif row.startswith("grover") and len(f) > 1: checks.add(f[1])
    before = " / ".join(fmt(c[0]) if c[0] is not None else "failed" for c in cells)
    after = " / ".join(fmt(c[1]) if c[1] is not None else "failed" for c in cells)
    values = [num(c) for c in checks]
    if row == "first-state":
        same = ""
    elif checks and None not in values:
        same = "yes" if max(values) - min(values) < 1e-12 else "NO: " + ", ".join(sorted(checks))
    else:
        same = "yes" if len(checks) == 1 else "NO: " + ", ".join(sorted(checks))
    print(f"| {label[row]} | {before} | {after} | {' / '.join(ratios)} | {same} |")
