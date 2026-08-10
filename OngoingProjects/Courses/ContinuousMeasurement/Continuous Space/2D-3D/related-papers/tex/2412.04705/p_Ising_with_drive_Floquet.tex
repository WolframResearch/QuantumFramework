import matplotlib.pyplot as plt
import mpl_setup
import numpy as np

import qutip as qt
import time

PI = np.pi


## ====================================================
#                    Input Parameters
## ====================================================

# ----------------------------------------
# Defining the Hamiltonian constants
# ----------------------------------------
# Number of spins
N = 4
# Energy-splitting
g0 = 1
# Coupling strength
J0 = 1.4
# Drive strength
A = 1.  # Drive norm.

# Drive frequency and period
omega = 1.0 * 2 * PI  # drive angular frequency
T = 2*PI / omega  # Drive period

# ---------------------------------
# Numerical sampling
# ---------------------------------
N_T = 10
dt = 0.01
tlist = np.arange(0, N_T * T, dt)
tlist_0 = np.arange(0, T, dt)  # One period tlist


# ======================================================
# Creating Hamiltonian
# ======================================================

def setup_Ising_drive(N, g0, J0, A, omega, data_type='CSR'):
    """
    # N    : number of spins
    # g0   : splitting,
    # J0   : couplings
    # A    : drive amplitude
    # omega: drive frequency
    """
    with qt.CoreOptions(default_dtype=data_type):

        sx_list, sy_list, sz_list = [], [], []
        for i in range(N):
            op_list = [qt.qeye(2)] * N 
            op_list[i] = qt.sigmax().to(data_type)
            sx_list.append(qt.tensor(op_list))
            op_list[i] = qt.sigmay().to(data_type)
            sy_list.append(qt.tensor(op_list))
            op_list[i] = qt.sigmaz().to(data_type)
            sz_list.append(qt.tensor(op_list))

        # Hamiltonian - Energy splitting terms
        H_0 = 0.
        for i in range(N):
            H_0 += g0 * sz_list[i]

        # Interaction terms
        H_1 = qt.qzero_like(H_0)
        for n in range(N - 1):
            H_1 += -J0 * sx_list[n] * sx_list[n + 1]

        # Driving terms
        if A > 0:
            H_d = 0.
            for i in range(N):
                H_d += A * sx_list[i]
            args = {'w': omega}
            H = [H_0, H_1, [H_d, lambda t, w: np.sin(w * t)]]
        else:
            args = {}
            H = [H_0, H_1]

        # Defining initial conditions
        state_list = [qt.basis(2, 1)] * (N-1)
        state_list.append(qt.basis(2, 0))
        psi0 = qt.tensor(state_list)

        # Defining expectation operator
        e_ops = sz_list
        return H, psi0, e_ops, args

options = {"progress_bar": False,
           "store_floquet_states": True}

H, psi0, e_ops, args = setup_Ising_drive(N, g0, J0, A, omega)


# ======================================================
# Initialization of arrays
# ======================================================
computational_time_floquet = np.ones(tlist.shape) * np.nan
computational_time_sesolve = np.ones(tlist.shape) * np.nan

expect_floquet = np.zeros((N, len(tlist)))
expect_sesolve = np.zeros((N, len(tlist)))


# ======================================================
# Starting simulation
# ======================================================
for n, t in enumerate(tlist):
    # --------------------------------
    # Floquet basis
    # --------------------------------
    tic_f = time.perf_counter()
    if t < T:
        # find the floquet modes for the time-dependent hamiltonian
        floquetbasis = qt.FloquetBasis(H, T, args)
    else:
        floquetbasis = qt.FloquetBasis(H, T, args, precompute=tlist_0)

    # Decomposing inital state into Floquet modes
    f_coeff = floquetbasis.to_floquet_basis(psi0)
    # Obtain evolved state in the original basis
    psi_t = floquetbasis.from_floquet_basis(f_coeff, t)
    p_ex = qt.expect(e_ops, psi_t)
    toc_f = time.perf_counter()

    # Saving data
    computational_time_floquet[n] = toc_f - tic_f

    # --------------------------------
    # sesolve
    # --------------------------------
    tic_f = time.perf_counter()
    output = qt.sesolve(H, psi0, tlist[:n + 1], e_ops=e_ops, args=args)
    p_ex_r = output.expect
    toc_f = time.perf_counter()
    # Saving data
    computational_time_sesolve[n] = toc_f - tic_f

    for i in range(N):
        expect_floquet[i, n] = p_ex[i]
        expect_sesolve[i, n] = p_ex_r[i][-1]


# ======================================================
# Plotting results
# ======================================================
sty = ["-", "--"]
fig, axs = plt.subplots(1, 2, figsize=(13.6, 4.5))
ax = axs[0]
ax.plot(tlist / T, computational_time_floquet, sty[0])
ax.plot(tlist / T, computational_time_sesolve, sty[1])

axs[0].set_yscale("log")
axs[0].set_yticks(np.logspace(-3, -1, 3))
lg = []
for i in range(N):
    axs[1].plot(tlist / T, np.real(expect_floquet[i, :]), sty[0], lw='3.5')
    lg.append(f"n={i+1}")
for i in range(N):
    axs[1].plot(tlist / T, np.real(expect_sesolve[i, :]), sty[1], lw='2.5', color="GREEN")

axs[0].set_xlabel(r'$t \, / \, T$')
axs[1].set_xlabel(r'$t \, / \, T$')
ax.set_ylabel(r'$Computational \; Time,\; [s]$')
axs[1].set_ylabel(r'$\langle \sigma_z^{{(n)}} \rangle$')
ax.legend((r"Floquet", r"sesolve"), frameon=False)
axs[1].legend(lg, frameon=False, ncol=2)
xticks = np.rint(np.linspace(0, N_T, N_T + 1, endpoint=True))
axs[0].set_xticks(xticks)
axs[0].set_xlim([0, N_T])
axs[1].set_xticks(xticks)
axs[1].set_xlim([0, N_T])
txt = f'$N={N}$, $g_0={g0}$, $J_0={J0}$, $A={A}$'
axs[0].text(.0, 1.05, txt, transform=axs[0].transAxes)
axs[0].text(-.15, 1.05, r'(a)', transform=axs[0].transAxes)
axs[1].text(-.15, 1.05, r'(b)', transform=axs[1].transAxes)

fig.show()
fig.savefig('./floquet_speed_test_ising.pdf')
