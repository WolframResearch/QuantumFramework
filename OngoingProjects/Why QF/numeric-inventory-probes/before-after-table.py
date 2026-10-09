"""Summarizes before-after.tsv into a markdown table: the median of the repeats for each row and
environment, with a check that every environment gave the same result.
Usage: python3 -I before-after-table.py [before-after.tsv]"""
import sys, statistics, collections, pathlib

path = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).with_name("before-after.tsv"))
runs = collections.defaultdict(list)    # (row, env) -> list of field lists
checks = collections.defaultdict(set)   # row -> set of (env, check)
loads = []
for line in path.read_text().splitlines():
    if not line or line.startswith("#"):
        continue
    rep, env, versions, load, row, *fields = line.split("\t")
    runs[(row, env)].append(fields)
    loads.append(float(load))

def number(text):
    try:
        return float(text)
    except ValueError:
        return None

limits = {"fourier-basis-16": 60}  # seconds before a row times out; 240 for the others

def median(values, row=None):
    nums = [number(v) for v in values]
    if any(n is None for n in nums):
        failures = sorted({v for v, n in zip(values, nums) if n is None})
        return {"memory": "out of 4 GB", "timeout": f"over {limits.get(row, 240)} s"}.get(failures[0], failures[0])
    m = statistics.median(nums)
    return f"{m:.3g}" if m < 10 else f"{m:.1f}"

def ratio(row):
    """before over after, from the medians; for a pair, of each of the two"""
    def med(env, i):
        fields = runs.get((row, env), [])
        nums = [number(f[i]) for f in fields]
        return None if not nums or any(n is None for n in nums) else statistics.median(nums)
    parts = []
    for i in ([0, 1] if row in two_numbers else [0]):
        b, a = med("before", i), med("after", i)
        if a is None:
            parts.append("")
        elif b is None:
            parts.append("before failed")
        else:
            parts.append(f"{b / a:.1f}x")
    return " / ".join(parts)

two_numbers = {"apply-qft14-zero", "apply-qft14-random", "apply-oracle10-zero", "apply-oracle10-uniform"}
label = {
    "build-qft14": "QFT, 14 qubits: build (s)",
    "build-oracle10": "Phase oracle, 10 variables: build (s)",
    "controlled-gate": "Controlled phase gate of the QFT: build (ms per gate)",
    "apply-qft14-zero": "QFT, 14 qubits, applied to the zero state with qc[] (s, first / second)",
    "apply-qft14-random": "QFT, 14 qubits, applied to a random state (s, first / second)",
    "apply-oracle10-zero": "Phase oracle, 10 variables, applied to the zero state (s, first / second)",
    "apply-oracle10-uniform": "Phase oracle, 10 variables, applied to the uniform state (s, first / second)",
    "quench-12": "Trotter quench, 4 steps, 12 qubits, default method (s)",
    "quench-16": "Trotter quench, 4 steps, 16 qubits, default method (s)",
    "quench-20": "Trotter quench, 4 steps, 20 qubits, default method (s)",
    "fourier-basis-16": "Change of basis into QuantumBasis[\"Fourier\"[16]] (s)",
}
order = list(label)
envs = ["before", "middle", "after"]
print("| Workload | Before | QF before, new TN and Arrays | After | Before / after | Same result |")
print("|---|---|---|---|---|---|")
for row in order:
    cells = []
    results = set()
    for env in envs:
        fields = runs.get((row, env), [])
        if not fields:
            cells.append("not run")
            continue
        if row in two_numbers:
            cells.append(median([f[0] for f in fields], row) + " / " + median([f[1] for f in fields], row))
            results |= {f[2] for f in fields if len(f) > 2}
        else:
            cells.append(median([f[0] for f in fields], row))
            if row.startswith("quench") or row == "fourier-basis-16":
                results |= {f[1] for f in fields if len(f) > 1 and f[1] not in ("none",)}
    results.discard("None")
    values = [number(r) for r in results]
    if results and all(v is not None for v in values):
        same = "yes" if max(values) - min(values) < 1e-12 else "NO: " + ", ".join(sorted(results))
    else:
        same = "yes" if len(results) <= 1 else "NO: " + ", ".join(sorted(results))
    if row.startswith("build") or row == "controlled-gate":
        same = ""
    print(f"| {label[row]} | " + " | ".join(cells) + f" | {ratio(row)} | {same} |")
if loads:
    print(f"\nload average during the runs: {min(loads):.1f} to {max(loads):.1f}, median {statistics.median(loads):.1f}")
