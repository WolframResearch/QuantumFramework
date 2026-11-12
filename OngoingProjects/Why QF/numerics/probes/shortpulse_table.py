"""Markdown tables for the short-pulse head-to-head, from shortpulse_qutip.tsv, shortpulse_qf.tsv
and shortpulse_reference.tsv (all next to this file). Prints to stdout.

Run: python3 -I shortpulse_table.py
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        lines = fh.read().splitlines()
    versions = [l[2:] for l in lines if l.startswith("# ")]
    rows = list(csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t"))
    return versions, rows


vq, qutip = load("shortpulse_qutip.tsv")
vf, qf = load("shortpulse_qf.tsv")
vr, ref = load("shortpulse_reference.tsv")

# exact (P1, P2) of the transmon from the WorkingPrecision-40 NDSolve reference
exact4 = {r["variant"]: [float(x) for x in r["result"].split(", ")][1:] for r in ref}


def fnum(x, digits=6):
    return "" if x in ("", None) else f"{float(x):.{digits}g}"


def ferr(x):
    return "" if x in ("", None) else f"{float(x):.1e}"


def fsec(x):
    if x in ("", None):
        return ""
    s = float(x)
    return f"{s * 1000:.2g} ms" if s < 1 else f"{s:.3g} s"


def scored(r):
    """(result, error) of a row; a transmon row is scored by the larger error of P1 and P2."""
    if r is None:
        return "", ""
    if r["case"] == "4 transmon":
        if not r["result"]:
            return "failed", ""
        p = [float(x) for x in r["result"].split(", ")]
        e = exact4[r["variant"].split(",")[0]]
        return f"P1 {p[0]:.6g}, P2 {p[1]:.3g}", ferr(max(abs(p[0] - e[0]), abs(p[1] - e[1])))
    if not r["result"]:
        return "failed", ""
    return fnum(r["result"]), ferr(r["abs_error"])


def pick(rows, case, variant, setting):
    return next((r for r in rows if r["case"] == case and r["variant"] == variant and r["setting"].startswith(setting)), None)


def advised(case, variant):
    """The QuTiP row with the docstring's max_step; when max_step alone exceeded nsteps, the rerun with nsteps raised, marked +."""
    ms = pick(qutip, case, variant, "max_step 0.49") or next(
        (r for r in qutip if r["case"] == case and r["variant"] == variant and r["setting"].startswith("max_step") and "nsteps" not in r["setting"] and "docs" not in r["setting"]), None)
    if ms is not None and not ms["result"]:
        raised = next((r for r in qutip if r["case"] == case and r["variant"] == variant and "nsteps" in r["setting"]), None)
        if raised is not None:
            res, err = scored(raised)
            return res + " +", err, fsec(raised["seconds"])
    res, err = scored(ms)
    return res, err, fsec(ms["seconds"]) if ms else ""


def exact_of(case, variant):
    if case == "4 transmon":
        e = exact4[variant.split(",")[0]]
        return f"P1 {e[0]:.6g}, P2 {e[1]:.3g}"
    r = pick(qutip, case, variant, "default")
    return fnum(r["exact"], 6) if r else ""


main = [
    ("1 docs", "pulse 0.7<t<0.75", "QuTiP docs example, <n>(1)"),
    ("2a #1265", "wait 0 (c_x)", "#1265, pi/2 + pi/2, no wait, Pe"),
    ("2a #1265", "wait 40 (c_x_bug)", "#1265, pi/2, wait 40, pi/2, Pe"),
    ("2b #1945", "tlist linspace(0,700,268)", "#1945, 268 output times, max error of <sz>"),
    ("2c #2552", "index93", "#2552 index93, sampled pulses, P0"),
    ("2d #2814 train", "5 Gaussian pi/5, w 0.04", "#2814-style train of 5 pulses, P1"),
] + [("3 pi pulse", f"{s}, w {w}", f"pi pulse, {s}, w = {w}") for w in ("1.0", "0.1", "0.01", "0.001")
     for s in ("square Piecewise", "square UnitStep", "Gaussian sigma=w/4", "sech tau=w/4")] + [
    ("4 transmon", f"{k}, gap {g}", f"transmon, {k} Gaussian, gap {g} ns") for k in ("windowed", "untruncated") for g in (0, 1000)] + [
    ("5 square from t0, then Gaussian", f, f"square pi/2 from t = 0 as {f}, Gaussian pi/2 at 50") for f in
    ("Boole[t>=0] Boole[t<=1]", "UnitStep[t] UnitStep[1-t]", "Boole[0<=t<=1]")]

print("Versions:\n")
for source, vs in (("QuTiP probe", vq), ("QF probe", vf), ("Transmon reference", vr)):
    for v in vs:
        print(f"- {source}: {v}")
print()
print("| Case | Exact | QF default | error | time | QuTiP default | error | time | QuTiP, max_step 0.49 w | error | time |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for case, variant, label in main:
    f_res, f_err = scored(pick(qf, case, variant, "default"))
    f_t = fsec((pick(qf, case, variant, "default") or {}).get("seconds"))
    q = pick(qutip, case, variant, "default")
    q_res, q_err = scored(q)
    a_res, a_err, a_t = advised(case, variant)
    print(f"| {label} | {exact_of(case, variant)} | {f_res} | {f_err} | {f_t} | {q_res} | {q_err} | {fsec(q['seconds']) if q else ''} | {a_res} | {a_err} | {a_t} |")

print("\nQF with and without its restart events (NDSolveValue at its defaults on the same equations):\n")
print("| Case | Restart times | QF default | error | QF, restart events removed | error |")
print("|---|---|---|---|---|---|")
for case, variant, label in main + [("2a #1265", r["variant"], "#1265, " + r["variant"]) for r in qf if r["case"] == "2a #1265" and r["setting"] == "default" and r["variant"] not in ("wait 0 (c_x)", "wait 40 (c_x_bug)")] + [
        ("4 transmon", f"{k}, gap {g}", f"transmon, {k} Gaussian, gap {g} ns") for k in ("windowed", "untruncated") for g in (100, 10000)]:
    d = pick(qf, case, variant, "default")
    c = pick(qf, case, variant, "control")
    if d is None:
        continue
    restarts = next((n.split()[-1] for n in d["note"].split("; ") if n.startswith("restart times")), "")
    print(f"| {label} | {restarts} | {' | '.join(scored(d))} | {' | '.join(scored(c))} |")

print("\n#1265, Pe at t = 99.5 against the wait between the two pi/2 pulses (exact 1):\n")
print("| Wait | QF default | error | QuTiP default | error | QuTiP, max_step 9.8 | error |")
print("|---|---|---|---|---|---|---|")
for r in [r for r in qutip if r["case"] == "2a #1265" and r["setting"] == "default"]:
    f_res, f_err = scored(pick(qf, "2a #1265", r["variant"], "default"))
    a_res, a_err, _ = advised("2a #1265", r["variant"])
    q_res, q_err = scored(r)
    print(f"| {r['variant'].replace('wait ', '')} | {f_res} | {f_err} | {q_res} | {q_err} | {a_res} | {a_err} |")

print("\nOther rows:\n")
print("| Case | Package, setting | Result | Error | Time | Note |")
print("|---|---|---|---|---|---|")
shown = {(c, v) for c, v, _ in main} | {("2a #1265", r["variant"]) for r in qutip if r["case"] == "2a #1265"}
for r in qutip + qf:
    if r["setting"].startswith("control"):
        continue
    extra = ((r["case"], r["variant"]) not in shown and not r["case"].startswith("4")) or "tlist linspace(0,100,201)" in r["setting"] \
        or "docs example" in r["setting"] or r["setting"].startswith("diagnostic") or "StateVector" in r["setting"]
    if r["case"].startswith("4") and r["variant"] in ("windowed", "untruncated"):
        extra = True
    if extra:
        res, err = scored(r) if r["result"] or r["case"] != "4 transmon" else ("failed", "")
        if r["package"].startswith("SciPy"):
            p = [float(x) for x in r["result"].split(", ")]
            e = exact4[r["variant"]]
            res, err = f"P1 {p[0]:.12g}, P2 {p[1]:.6g}", ferr(max(abs(p[0] - e[0]), abs(p[1] - e[1])))
        print(f"| {r['case']}, {r['variant']} | {r['package']}, {r['setting']} | {res} | {err} | {fsec(r['seconds'])} | {r['note'][:110]} |")
