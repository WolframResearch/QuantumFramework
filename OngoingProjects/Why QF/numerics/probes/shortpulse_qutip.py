"""Short-pulse head-to-head, QuTiP side.

Each case runs twice: with every solver option at its default, and with the
max_step that the sesolve/mesolve docstring advises ("less than half the width
of the thinnest pulse"), taken as 0.49 w. When that max_step alone exceeds the
default nsteps (2500 internal steps per tlist interval), the run is repeated
with nsteps raised to fit. Every result is scored against its exact value.

Writes shortpulse_qutip.tsv next to this file. Run with the venv that holds
qutip 5.3.1:  <venv>/bin/python -I shortpulse_qutip.py
"""

import math
import os
import platform
import statistics
import sys
import time
import warnings

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.special import erf

import qutip
from qutip import QobjEvo, basis, destroy, mesolve, num, sesolve, sigmax, sigmay, sigmaz

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "shortpulse_qutip.tsv")
DATA = os.path.join(HERE, "shortpulse_data2552")
REPS = 3
ADVISED = 0.49  # max_step as a fraction of the pulse width

rows = []


def run(f):
    """Median wall time of REPS calls to f, the last value, and any warning or error text."""
    times, value, notes = [], None, []
    for _ in range(REPS):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            t0 = time.perf_counter()
            try:
                value = f()
            except Exception as e:
                times.append(time.perf_counter() - t0)
                notes.append(f"{type(e).__name__}: {str(e).splitlines()[0][:100]}")
                value = None
                break
            times.append(time.perf_counter() - t0)
        notes += [f"{w.category.__name__}: {str(w.message).splitlines()[0][:70]}" for w in caught]
    return value, statistics.median(times), "; ".join(dict.fromkeys(notes))


def record(case, variant, setting, observable, value, exact, seconds, note=""):
    err = "" if value is None or exact is None else f"{abs(value - exact):.3e}"
    rows.append([case, variant, "QuTiP 5.3.1", setting, observable,
                 "" if value is None else f"{value:.12g}",
                 "" if exact is None else f"{exact:.12g}", err, f"{seconds:.4f}", note])
    print("\t".join(rows[-1]), flush=True)


def settings(width, span):
    """Default options, the advised max_step, and the advised max_step with nsteps raised to fit span."""
    ms = ADVISED * width
    return [("default", None),
            (f"max_step {ms:.4g}", {"max_step": ms}),
            (f"max_step {ms:.4g}, nsteps {math.ceil(2 * span / ms)}", {"max_step": ms, "nsteps": math.ceil(2 * span / ms)})]


def sweep(case, variant, observable, solve, exact, width, span, note=""):
    """Run solve(options) at each setting; the nsteps row runs only when the max_step row failed."""
    failed = False
    for name, opts in settings(width, span):
        if "nsteps" in name and not failed:
            continue
        value, sec, msg = run(lambda: solve(opts))
        failed = value is None
        record(case, variant, name, observable, value, exact, sec, "; ".join(filter(None, [note, msg])))


# 1. QuTiP docs, "Working with pulses": pulse on 0.7 < t < 0.75, sigmaz + pulse sigmax, psi0 = |1>.
def docs_solve(opts):
    def pulse(t):
        return 10 * np.pi * (0.7 < t < 0.75)

    tlist = np.linspace(0, 1, 201)
    H = [sigmaz(), [sigmax(), pulse]]
    psi0 = basis(2, 1)
    return sesolve(H, psi0, tlist, e_ops=num(2), options=opts).expect[0][-1]


a = 10 * np.pi
W = math.sqrt(1 + a**2)
docs_exact = math.cos(0.05 * W) ** 2 + math.sin(0.05 * W) ** 2 / W**2
sweep("1 docs", "pulse 0.7<t<0.75", "<n>(1)", docs_solve, docs_exact, 0.05, 1)
value, sec, msg = run(lambda: docs_solve({"max_step": 0.01}))
record("1 docs", "pulse 0.7<t<0.75", "max_step 0.01 (docs example)", "<n>(1)", value, docs_exact, sec, msg)


# 2a. Issue #1265: two cosine pi/2 pulses (width 20) with a wait between them, mesolve with
# T1 = T2 = inf (zero collapse operators), Pe = (1 - <sz>)/2 at the last tlist point 99.5. Exact Pe = 1.
# Edits for QuTiP 5: e_ops passed by keyword, progress_bar dropped, lambda (t, arg) renamed (t, args).
piLen = 20
piAmp = np.pi * 1.0 / 20
T1, T2 = float("inf"), float("inf")
Hq = qutip.Qobj(np.diag([0, 0]), dims=[[2], [2]])
Hlbld = [np.sqrt(1.0 / T1) * qutip.destroy(2), np.sqrt(2.0 / T2) * qutip.num(2)]


def cosine_pulse(t, t0=piLen / 2.0, amp=piAmp, w=piLen, phase=0, df=0):
    x = t - t0
    return amp * 0.5 * (1 + np.cos(2 * np.pi * x / w)) * ((x + w / 2.0) > 0) * ((-x + w / 2.0) > 0) * \
        np.exp(1j * phase - 2j * np.pi * df * x)


def seq1265(wait):
    """Second pulse starts `wait` after the first ends: wait 0 is c_x, wait 40 is c_x_bug."""
    return lambda t, args: cosine_pulse(t, t0=piLen / 2.0, amp=0.5 * piAmp) + \
        cosine_pulse(t, t0=wait + 1.0 * piLen + piLen / 2.0, amp=0.5 * piAmp)


def solve1265(wait, opts):
    H = [Hq, [sigmax(), seq1265(wait)]]
    r = mesolve(H, basis(2, 0), np.arange(0, 100, 0.5), Hlbld,
                e_ops=[sigmax(), sigmay(), sigmaz()], options=opts)
    return (1 - r.expect[2][-1].real) / 2


for wait in [0, 10, 20, 23, 23.0017, 23.0018, 24, 30, 40, 50]:
    label = {0: " (c_x)", 40: " (c_x_bug)", 23.0018: " (BoxiLi threshold)"}.get(wait, "")
    sweep("2a #1265", f"wait {wait}{label}", "Pe(99.5)", lambda o, w=wait: solve1265(w, o), 1.0, piLen, 100)


# 2b. Issue #1945: square pulse on [100, 600] of 2 pi omega_R sigmax, mesolve, <sz> on linspace(0, 700, N).
# Scored by the largest error over the output grid against cos(2 Omega (t - 100)) inside the pulse.
def func(t, *args):
    omega_R = 2e-3
    if t < 100:
        return 0
    elif t > 600:
        return 0
    else:
        return 2 * np.pi * omega_R


def err1945(n, opts):
    times = np.linspace(0, 700, n)
    ez = mesolve([[sigmax(), func]], basis(2, 0), times, e_ops=[sigmaz()], options=opts).expect[0]
    exact = np.cos(2 * (2 * np.pi * 2e-3) * np.clip(times - 100, 0, 500))
    return float(np.max(np.abs(ez - exact)))


for n in [267, 268, 800]:
    sweep("2b #1945", f"tlist linspace(0,700,{n})", "max_t |<sz> - exact|", lambda o, n=n: err1945(n, o), 0.0, 500, 700)


# 2c. Issue #2552: two pi/2 pulses given as sampled arrays (cubic-spline coefficients), mesolve,
# |<0|psi_final>|^2. The drive is along sigmax up to an imaginary part below 1e-14, so the exact
# value is cos^2(Theta) with Theta the integral of the spline, computed exactly by Simpson's rule
# on each knot interval (exact for a cubic).
sx, sy = sigmax(), sigmay()
psi0_2552 = qutip.fock(2, 0)
Omega0 = 10e6 * 2 * np.pi
for name in ["index92", "index93", "index94"]:
    coef = np.loadtxt(os.path.join(DATA, f"{name}_coef.txt"), dtype=np.complex128)
    tlist = np.loadtxt(os.path.join(DATA, f"{name}_tlist.txt"), dtype=np.float64)

    def solve2552(opts, coef=coef, tlist=tlist):
        H = QobjEvo([[sx, (coef * Omega0).real / 2], [sy, (coef * Omega0).imag / 2]], tlist=tlist)
        result = mesolve(H, psi0_2552, tlist=tlist, options=opts)
        return float(np.abs(psi0_2552.dag() @ result.final_state) ** 2)

    spline = qutip.coefficient((coef * Omega0).real / 2, tlist=tlist)
    mids = (tlist[:-1] + tlist[1:]) / 2
    f0 = np.array([spline(t).real for t in tlist])
    fm = np.array([spline(t).real for t in mids])
    theta = float(np.sum(np.diff(tlist) / 6 * (f0[:-1] + 4 * fm + f0[1:])))
    sweep("2c #2552", name, "|<0|psi(T)>|^2", solve2552, math.cos(theta) ** 2, 9.95e-8, tlist[-1],
          f"Theta - pi/2 = {theta - np.pi / 2:.3e}")


# 2d. Issue #2814 has no code; its sketch is a pulse train with long waits. Five Gaussian pi/5 pulses,
# sigma = 0.01 (w = 4 sigma = 0.04), centred at 10, 30, 50, 70, 90 on [0, 100]: with k captured,
# P1 = sin^2(k pi/10), so P1 = 1 needs all five and P1 identifies k.
def train(t):
    s = 0.01
    return sum((np.pi / 5) / (np.sqrt(2 * np.pi) * s) * np.exp(-(t - c) ** 2 / (2 * s**2)) for c in (10, 30, 50, 70, 90))


def solve_train(opts):
    return sesolve([[sigmax() / 2, train]], basis(2, 0), [0, 100], e_ops=[num(2)], options=opts).expect[0][-1].real


def pulses_captured(p1):
    return f"pulses captured {10 / np.pi * np.arcsin(np.sqrt(min(max(p1, 0), 1))):.2f} of 5" if p1 is not None else ""


for name, opts in settings(0.04, 100):
    value, sec, msg = run(lambda: solve_train(opts))
    record("2d #2814 train", "5 Gaussian pi/5, w 0.04", name, "P1(100)", value, 1.0, sec,
           "; ".join(filter(None, [pulses_captured(value), msg])))


# 3. A pi pulse at t = 50 on [0, 100] under H = Omega(t) X/2, four shapes of width w. Exact P1(100) = 1.
def shapes(w):
    s = w / 4
    return {
        "square Piecewise": lambda t: np.pi / w if abs(t - 50) < w / 2 else 0.0,
        "square UnitStep": lambda t: np.pi / w * np.heaviside(t - (50 - w / 2), 1.0) * np.heaviside((50 + w / 2) - t, 1.0),
        "Gaussian sigma=w/4": lambda t: np.pi / (np.sqrt(2 * np.pi) * s) * np.exp(-(t - 50) ** 2 / (2 * s**2)),
        # sech(x) = 2 e^-|x| / (1 + e^-2|x|), the overflow-free form of 1/cosh(x)
        "sech tau=w/4": lambda t: (1 / s) * 2 * np.exp(-abs(t - 50) / s) / (1 + np.exp(-2 * abs(t - 50) / s)),
    }


def solve_pi(om, opts, tlist=(0, 100)):
    return sesolve([[sigmax() / 2, om]], basis(2, 0), list(tlist), e_ops=[num(2)], options=opts).expect[0][-1].real


for w in [1.0, 0.1, 0.01, 0.001]:
    for shape, om in shapes(w).items():
        sweep("3 pi pulse", f"{shape}, w {w}", "P1(100)", lambda o, om=om: solve_pi(om, o), 1.0, w, 100)
        value, sec, msg = run(lambda om=om: solve_pi(om, None, np.linspace(0, 100, 201)))
        record("3 pi pulse", f"{shape}, w {w}", "default, tlist linspace(0,100,201)", "P1(100)", value, 1.0, sec, msg)


# 4. Three-level transmon, alpha/2pi = -200 MHz (time in ns), rotating frame,
# H = alpha/2 a+a+aa + Omega(t)/2 (a + a+), Gaussian pi pulse after an idle gap, no DRAG.
# "windowed": sigma = T/4, T = 16, truncated to [gap, gap + T], area pi over the window, t_end = gap + T.
# "untruncated": sigma = 4, centre gap + 48, t_end = gap + 96.
# The exact value comes from shortpulse_reference.wls; |0> is stationary during the idle gap, so it
# does not depend on the gap. SciPy DOP853 at rtol 1e-13 over the pulse alone is a second reference.
alpha = -2 * np.pi * 0.2
T, sig = 16.0, 4.0
A_win = np.pi / (sig * np.sqrt(2 * np.pi) * erf(T / (2 * np.sqrt(2) * sig)))
b = destroy(3)
H0 = alpha / 2 * b.dag() * b.dag() * b * b
Hd = (b + b.dag()) / 2
proj = [basis(3, k) * basis(3, k).dag() for k in range(3)]


def transmon_pulse(kind, gap):
    if kind == "windowed":
        return (lambda t: A_win * np.exp(-(t - (gap + T / 2)) ** 2 / (2 * sig**2)) if gap <= t <= gap + T else 0.0), gap + T
    return (lambda t: np.pi / (np.sqrt(2 * np.pi) * sig) * np.exp(-(t - (gap + 48)) ** 2 / (2 * sig**2))), gap + 96


def solve_transmon(kind, gap, opts):
    om, t_end = transmon_pulse(kind, gap)
    r = sesolve([H0, [Hd, om]], basis(3, 0), [0, t_end], e_ops=proj, options=opts)
    return [float(r.expect[k][-1].real) for k in range(3)]


for kind in ["windowed", "untruncated"]:
    om, t_end = transmon_pulse(kind, 0)
    H0d, Hdd = H0.full(), Hd.full()
    sol = solve_ivp(lambda t, y: -1j * (H0d + om(t) * Hdd) @ y, (0, t_end), np.array([1, 0, 0], complex),
                    method="DOP853", rtol=1e-13, atol=1e-15, max_step=0.01)
    pref = np.abs(sol.y[:, -1]) ** 2
    rows.append(["4 transmon", f"{kind}", "SciPy DOP853 reference", "rtol 1e-13, max_step 0.01, pulse only",
                 "P1, P2", f"{pref[1]:.15g}, {pref[2]:.15g}", "", "", "", ""])
    print("\t".join(rows[-1]), flush=True)
    for gap in [0, 100, 1000, 10000]:
        failed = False
        for name, opts in settings(T if kind == "windowed" else 4 * sig, gap + t_end):
            if "nsteps" in name and not failed:
                continue
            value, sec, msg = run(lambda o=opts: solve_transmon(kind, gap, o))
            failed = value is None
            rows.append(["4 transmon", f"{kind}, gap {gap}", "QuTiP 5.3.1", name, "P1, P2",
                         "" if value is None else f"{value[1]:.12g}, {value[2]:.12g}", "", "", f"{sec:.4f}", msg])
            print("\t".join(rows[-1]), flush=True)


# 5. A pi/2 square pulse on [0, 1], then a Gaussian pi/2 pulse (sigma = 0.025, w = 0.1) at t = 50, under
# H = Omega(t) X/2: exact P1(100) = 1, and a skipped Gaussian leaves 1/2. The square is written three ways,
# which QF reads differently (shortpulse_qf.wls); here they are the same function.
def edge_gauss(t):
    return (np.pi / 2) / (np.sqrt(2 * np.pi) * 0.025) * np.exp(-(t - 50) ** 2 / (2 * 0.025**2))


edge_forms = {
    "Boole[t>=0] Boole[t<=1]": lambda t: np.pi / 2 * (t >= 0) * (t <= 1) + edge_gauss(t),
    "UnitStep[t] UnitStep[1-t]": lambda t: np.pi / 2 * np.heaviside(t, 1.0) * np.heaviside(1 - t, 1.0) + edge_gauss(t),
    "Boole[0<=t<=1]": lambda t: np.pi / 2 * (0 <= t <= 1) + edge_gauss(t),
}
for form, om in edge_forms.items():
    sweep("5 square from t0, then Gaussian", form, "P1(100)", lambda o, om=om: solve_pi(om, o), 1.0, 0.1, 100)


versions = (f"# QuTiP {qutip.__version__}, SciPy {scipy.__version__}, NumPy {np.__version__}, "
            f"Python {sys.version.split()[0]}, {platform.platform()}, {time.strftime('%Y-%m-%d %H:%M')}")
header = ["case", "variant", "package", "setting", "observable", "result", "exact", "abs_error", "seconds", "note"]
with open(OUT, "w") as fh:
    fh.write(versions + "\n" + "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
print(versions)
