import numpy as np
import scipy.sparse as sp
import qutip as qt
from qutip.core.energy_restricted import EnrSpace

from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_setup import * 

import time as time
from numpy.random import RandomState
prng = RandomState(int(time.time()))


# ---------- Utility Functions ----------

def _enr_project_state(dims, excitations, state):
    """
    Generate a projector onto a Fock state representation in a excitation-number restricted
    state space.  
    """
    nstates, state2idx, idx2state = qt.enr_state_dictionaries(dims, excitations)
    a_op = sp.lil_matrix((nstates, nstates), dtype=np.complex128)

    try:
        n = state2idx[tuple(state)]
    except KeyError: 
        raise ValueError("The state tuple %s is not in the restricted "
                         "state space" % str(tuple(state)))
    a_op[n, n] = 1

    enr_dims = [EnrSpace(dims, excitations)] * 2
    return qt.Qobj(a_op, dims=enr_dims)


def _enr_project_sub_system(dims, excitations, mode, s0):
    """
    Projector onto a particular state in position "mode" with excitation s0
    E.g., mode=0, s0=1, excitations=2, returns |1,0><1,0| + |1,1><1,1|
    """
    nstates, state2idx, idx2state = qt.enr_state_dictionaries(dims, excitations)
    a_op = sp.lil_matrix((nstates, nstates), dtype=np.complex128)

    for n1, state1 in enumerate(list(idx2state.values())):
        for idx, s in enumerate(state1):
            # if s > 0, the annihilation operator of mode idx has a non-zero
            # entry with one less excitation in mode idx in the final state
            if idx == mode and s0==s:
                a_op[n1, n1] = 1

    enr_dims = [EnrSpace(dims, excitations)] * 2
    return qt.Qobj(a_op, dims=enr_dims)


def _enr_projector_subdim(dims, excitations, subrange, n_exc):
    """
    Generate the projector onto excitation subspace
    """
    nstates, state2idx, idx2state = qt.enr_state_dictionaries(dims, excitations)
    a_op = sp.lil_matrix((nstates, nstates), dtype=np.complex128)

    for n1, state1 in enumerate(list(idx2state.values())):
        if sum([state1[i] for i in subrange]) == n_exc:
            a_op[n1, n1] = 1

    enr_dims = [EnrSpace(dims, excitations)] * 2
    return qt.Qobj(a_op, dims=enr_dims)


# ---------- Solve problem ----------

def Waveguide_Solve(d, psi0, deltat, T, phi, N, wq=0, g=1, Ntraj=500):
    """
    A basic implementation of the space-discretized waveguide model from 
    Regidor et al, PHYSICAL REVIEW RESEARCH 3, 023030 (2021).
    Note this only implements a single qubit coupled to a waveguide
    truncated by a mirror at one end, and only in the single-excitation limit.

    Will evolve the state ("psi0") until time ("T") in steps ("deltat"),
    and average over ("Ntraj") iterations

    Parameters
    ----------

    d : list of [:obj:`.Qobj`] of length ("N+1"). d[0:-1] must be the
        annhilation operators defining the ("N") discrete waveguide modes
        and d[-1] must be the annilation operator for the system (qubit)

    psi0 : :obj:`.Qobj`
        initial state vector (ket). 

    deltat : float, simultaneously defines the time step used in the discrete
        propagation below and the length of the waveguide via L = deltat * N

    T : float, total integration time.

    phi : float, phase acquired at the mirror truncating the waveguide.

    N : int, number of modes.

    wq : float, qubit energy

    g : float, coupling strength

    Ntraj: int, number of trajectories used for averaging.

    Returns
    -------

    output_states: list of [:obj:`.Qobj`] representing the qubit+waveguide 
        states at all times in [0, T] in time-steps of deltat
    """

    # list of annihilation operators for cavity modes
    a = d[0:-1]

    # atomic annihilation operator
    sm = d[-1]

    #example of states:
    #a[0].dag()|0> = [1,0,0,0,0]
    #a[-1].dag()|0> = [0,0,0,1,0]
    #sm.dag()|0> = [0,0,0,0,1]

    H0 = wq * sm.dag() * sm

    #Atom-Waveguide couplings
    #left moving box
    Hint_ac = (
        (np.sqrt(g) / np.sqrt(2 * deltat)) *
        (a[0].dag() * sm + sm.dag() * a[0])
    )

    #right moving box.
    Hint_ac += (
        (np.sqrt(g) / np.sqrt(2 * deltat)) * 
        (np.exp(-1.0j * phi) * a[-1].dag() * sm +  
         np.exp(1.0j * phi) * sm.dag() * a[-1])
    )

    H = H0 + Hint_ac
   
    #Utility functions:
    #projector onto single excitation in the last right moving box
    coll_a = _enr_project_sub_system(dims, excitations, N-1, 1)  

    #projector onto the empty state in the right-moving box
    coll_0 = _enr_project_sub_system(dims, excitations, N-1, 0) 

    #propagator for interaction evolution
    U = (-1.0j * deltat * H).expm()

    #projector onto 1 photon subspace
    proj_1ss = _enr_projector_subdim(dims, excitations, np.arange(0,N), 1)

    #construct propagator for moving state between boxes
    U_shift1 = 0
    #Translates all excitations one space to the left every time step.
    #If more than one total excitations, this is insufficient.
    for n in range(N-1):
        U_shift1 += a[n+1].dag() * a[n]
    U_shift = U_shift1 * proj_1ss
    #this adds on projectors onto the 0 photon subspace, to make sure they are not projected out by U_shift:
    U_shift += (
        _enr_project_state(dims,excitations, [0 if m == (N) else 0 for m in range(N+1)]) +
        _enr_project_state(dims,excitations, [1 if m == (N) else 0 for m in range(N+1)])
    )

    # Do the simulation:

    Nsteps = int(np.round(T / deltat))
    output_states = [0 for _ in range(Nsteps+1)]
    output_states[0] = psi0 * psi0.dag()

    for traj in tqdm(range(Ntraj)):

        psi1 = psi0
        nt = 0
        for _ in range(Nsteps):
            nt += 1
            psi1 = U * psi1

            Pout = qt.expect(coll_a, psi1)
            numpa = prng.rand()

            if numpa <= Pout:
                psi1 = coll_a * psi1
                psi1 = a[-1] * psi1
            else:
                psi1 = (coll_0) * psi1

            psi1 = U_shift * psi1
            psi1 = psi1.unit()

            output_states[nt] += psi1 * psi1.dag() / Ntraj

    return output_states


N = 21      # number of waveguide modes
M = 2       # normal, assumed, mode fock-state truncation
wq = 0.     # qubit frequency
g = 1       # coupling constant

#Physical dimensions of the qubit+waveguide modes
dims = [M] * N + [2]
#Restricted number of excitations in the qubit+waveguide modes
excitations = 1           # total number of excitations
initial_excitiations = 1  # initial number of excitations

#Construct ENR operators and states
d = qt.enr_destroy(dims, excitations)
psi0 = qt.enr_fock(dims, excitations, [initial_excitiations if m == N else 0 for m in range(N+1)])

deltat = 1 / (N-1)  #set up so that it arrives back to interact with system at t=(N-1)deltat
T=5  #total integration time

phi = np.pi #phase change at mirror
output_states_phipi = Waveguide_Solve(d, psi0, deltat, T, phi, N,
                                      wq, g,Ntraj = 4000)

phi = 0 #phase change at mirror
output_states_phi0 = Waveguide_Solve(d, psi0, deltat, T, phi, N,
                                     wq, g,Ntraj = 4000)


# ---------- Generate plots ----------

Nsteps = int(np.round(T / deltat))
times = np.linspace(0, T, Nsteps + 1)

plt.figure()
#plt.plot(times, [1/np.e for t in times])
plt.plot(times, qt.expect(d[-1].dag()*d[-1], output_states_phipi), label= r"With mirror, $\phi=\pi$")
plt.plot(times, qt.expect(d[-1].dag()*d[-1], output_states_phi0),'--', label= r"With mirror, $\phi=0$")
plt.plot(times, [(np.cos(np.sqrt(g*deltat/2))**(2*n) * np.exp(-(g/2) * n * deltat))  for n in range(Nsteps+1)], ':',label="Without mirror")

plt.legend(loc=0)
plt.xlabel(r"$t$ $(1/\gamma$)", fontsize=18)
plt.ylabel("TLS Population", fontsize=18)

plt.savefig("enr_waveguide_mirror.pdf")
plt.show()

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
e_ops = [x.dag() * x for x in d[0:-1]]
mode_occup_listphi0 = []
for tk,t in enumerate(times):
    mode_occup = []
    for kk in range(N):
        mode_occup.append(qt.expect(e_ops[kk],output_states_phi0[len(times)-tk-1]))
    mode_occup_listphi0.append(mode_occup)

mode_occup_listphipi = []
for tk,t in enumerate(times):
    mode_occup = []
    for kk in range(N):
        mode_occup.append(qt.expect(e_ops[kk],output_states_phipi[len(times)-tk-1]))
    mode_occup_listphipi.append(mode_occup)

X, Y = np.meshgrid(range(len(mode_occup_listphi0[0])),times[::-1])


fig = plt.figure(figsize=plt.figaspect(1.5))
fig.subplots_adjust(wspace=0, hspace=-0.15)

ax = fig.add_subplot(211, projection='3d')
ax.zaxis.set_rotate_label(False) 
ax.zaxis.set_label_position('top')
ax.plot_surface(X, -Y,np.array(mode_occup_listphi0), rstride=1, cstride=1, cmap=cm.viridis, alpha=0.9)
ax.set_yticks([-2,-4],["2","4"])
ax.tick_params(axis="y", pad=6)
ax.tick_params(axis="z", pad=8)
ax.tick_params(axis="x", pad=4)
ax.set_ylabel("t $(1/\gamma)$")
ax.set_xlabel("n")
#ax.set_zlabel("$B_n^{\dagger}B_n$")

ax.set_box_aspect((4,4,2))
ax.yaxis.labelpad=10
ax.zaxis.labelpad=5
ax.xaxis.labelpad=10
#plt.savefig("waveguide_phi0.pdf")


ax = fig.add_subplot(212, projection='3d')
ax.zaxis.set_rotate_label(False) 
ax.plot_surface(X, -Y,np.array(mode_occup_listphipi), rstride=1, cstride=1, cmap=cm.viridis, alpha=0.9)
ax.set_yticks([-2,-4],["2","4"])
ax.tick_params(axis="y", pad=6)
ax.tick_params(axis="z", pad=8)
ax.tick_params(axis="x", pad=4)
ax.set_ylabel("t $(1/\gamma)$")
ax.set_xlabel("n")

ax.set_box_aspect((4,4,2))
ax.yaxis.labelpad=10
ax.zaxis.labelpad=5
ax.xaxis.labelpad=10
fig.savefig('waveguide_occup.pdf', transparent=True)