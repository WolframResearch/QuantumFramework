import numpy as np
import qutip as qt

import qutip_qoc as qoc

# objective
initial = qt.qeye(2) # identity
target  = qt.gates.hadamard_transform()

# energy splitting, tunneling, amplitude damping
omega, delta, gamma = 0.1, 1.0, 0.1
sx, sy, sz = qt.sigmax(), qt.sigmay(), qt.sigmaz()

Hc = [sx, sy, sz] # control operator
Hc = [qt.liouvillian(H) for H in Hc]

Hd = 1 / 2 * (omega * sz + delta * sx) # drift term
Hd = qt.liouvillian(H=Hd, c_ops=[np.sqrt(gamma) * qt.sigmam()])

# combined operator list
H = [Hd, Hc[0], Hc[1], Hc[2]]

# pulse time interval
times = np.linspace(0, np.pi/2, 100)

# run the optimization
res_grape = qoc.optimize_pulses(
    objectives=qoc.Objective(initial, H, target),
    control_parameters={
        "ctrl_x": {
            "guess" : np.sin(times),
            "bounds": [-1, 1]
        },
        "ctrl_y": {
            "guess" : np.cos(times),
            "bounds": [-1, 1]
        },
        "ctrl_z": {
            "guess" : np.tanh(times),
            "bounds": [-1, 1]
        }
    },
    tlist=times,
    algorithm_kwargs={
        "alg": "GRAPE",
        "fid_err_targ": 0.01
    }
)

# ------------------------------------------- CRAB EXAMPLE -------------------------------------------

# c0 * sin(c2*t) + c1 * cos(c2*t) + ...
n_params = 3 # adjust in steps of 3

# run the optimization
res_crab = qoc.optimize_pulses(
    objectives=qoc.Objective(initial, H, target),
    control_parameters={
        "ctrl_x": {
            "guess" : [1 for _ in range(n_params)],
            "bounds": [(-1, 1)] * n_params
        },
        "ctrl_y": {
            "guess" : [1 for _ in range(n_params)],
            "bounds": [(-1, 1)] * n_params
        },
        "ctrl_z": {
            "guess" : [1 for _ in range(n_params)],
            "bounds": [(-1, 1)] * n_params
        }
    },
    tlist=times,
    algorithm_kwargs={
        "alg": "CRAB",
        "fid_err_targ": 0.01,
        "fix_frequency": False
    }
)

# ------------------------------------------- GOAT EXAMPLE -------------------------------------------

def sin(t, c):
    return c[0] * np.sin(c[1] * t)

def grad_sin(t, c, idx):
    if idx == 0: # w.r.t. c0
        return np.sin(c[1] * t)
    if idx == 1: # w.r.t. c1
        return c[0] * np.cos(c[1] * t) * t
    if idx == 2: # w.r.t. time
        return c[0] * np.cos(c[1] * t) * c[1]

# similar to qutip.QobjEvo
H = [Hd] + [[hc, sin, {"grad": grad_sin}] for hc in Hc]

ctrl_parameters = {
    id: {
        "guess":  [1, 0], # c0 and c1
        "bounds": [(-1, 1), (0, 2*np.pi)]
    } for id in ['x', 'y', 'z']
}

# magic kwrd to treat time as optimization variable
ctrl_parameters["__time__"] = {
    "guess": times[len(times) // 2],
    "bounds": [times[0], times[-1]],
}

# run the optimization
res_goat = qoc.optimize_pulses(
    objectives=qoc.Objective(initial, H, target),
    control_parameters=ctrl_parameters,
    tlist=times,
    algorithm_kwargs={
        "alg": "GOAT",
        "fid_err_targ": 0.01,
    }
)

# ------------------------------------------- JOPT EXAMPLE -------------------------------------------

import jax

@jax.jit
def sin_y (t, d, **kwargs):
    return d[0] * jax.numpy.sin(d[1] * t)

@jax.jit
def sin_z (t, e, **kwargs):
    return e[0] * jax.numpy.sin(e[1] * t)

import jax
    
@jax.jit
def sin_x (t, c, **kwargs):
    return c[0] * jax.numpy.sin(c[1] * t)

# same for sin_y and sin_z ...

H = [Hd] + [[Hc[0], sin_x], [Hc[1], sin_y], [Hc[2], sin_z]]

res_jopt = qoc.optimize_pulses(
    objectives=qoc.Objective(initial, H, target),
    control_parameters=ctrl_parameters, 
    tlist=times,
    algorithm_kwargs={
        "alg": "JOPT",
        "fid_err_targ": 0.01,
    }
)

# ------------------------------------------- PLOTS -------------------------------------------

import matplotlib.pyplot as plt
from mpl_setup import *

fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.54))

goat_range = times < res_goat.optimized_params[-1]
jopt_range = times < res_jopt.optimized_params[-1]

for i in range(3):
    ax[i].plot(times, res_grape.optimized_controls[i], ':', label='GRAPE')
    ax[i].plot(times, res_crab.optimized_controls[i], '-.', label='CRAB')
    ax[i].plot(times[goat_range], np.array(res_goat.optimized_controls[i])[goat_range], '-', label='GOAT')
    ax[i].plot(times[jopt_range], np.array(res_jopt.optimized_controls[i])[jopt_range], '--', label='JOPT')

    ax[i].set_xlabel(r"Time $t$")

ax[0].legend(loc=0)
ax[0].set_ylabel(r"Pulse amplitude $c_x(t)$", labelpad=-5)
ax[1].set_ylabel(r"Pulse amplitude $c_y(t)$", labelpad=-5)
ax[2].set_ylabel(r"Pulse amplitude $c_z(t)$", labelpad=-5)
ax[2].set_ylim(-0.2, 1.1) # ensure equal spacing between subplots

plt.savefig("qoc_pulse_plot.pdf")
plt.show()



# ------------------------------------------- FIDELITIES -------------------------------------------

print('GRAPE: ', res_grape.fid_err)
print(res_grape.total_seconds, ' seconds')
print()
print('CRAB : ', res_crab.fid_err)
print(res_crab.total_seconds, ' seconds')
print()
print('GOAT : ', res_goat.fid_err)
print(res_goat.total_seconds, ' seconds')
print()
print('JOPT : ', res_jopt.fid_err)
print(res_jopt.total_seconds, ' seconds')