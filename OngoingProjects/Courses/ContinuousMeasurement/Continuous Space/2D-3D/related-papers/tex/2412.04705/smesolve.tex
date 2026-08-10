import matplotlib.pyplot as plt
from mpl_setup import * 

import qutip as qt
import numpy as np


# parameters
N = 20                 # Hilbert space dimension
delta = 5 * 2 * np.pi  # cavity detuning
kappa = 1              # cavity decay rate
A = 4                  # intensity of initial state
num_traj = 500

# operators
a = qt.destroy(N)
x = a + a.dag()
H = delta * a.dag() * a

rho_0 = qt.coherent(N, np.sqrt(A))
times = np.arange(0, 1, 0.0025)

stoc_solution = qt.smesolve(
    H, rho_0, times, c_ops=[], sc_ops=[np.sqrt(kappa) * a], e_ops=[x],
    ntraj=num_traj, options={"dt": 0.00125, "store_measurement": True}
)

me_solution = qt.mesolve(H, rho_0, times, c_ops=[np.sqrt(kappa) * a], e_ops=[x])

plt.figure()
plt.plot(times[1:], np.array(stoc_solution.measurement).mean(axis=0)[0, :].real,
         lw=2, label=r'$J_x$')
plt.plot(times, stoc_solution.expect[0],
         label=r'$\langle x \rangle$')
plt.plot(times, me_solution.expect[0], '--', color=GRAY,
         label=r'$\langle x \rangle$ mesolve')

plt.xlabel('Time', fontsize=18)
plt.legend()
plt.savefig("smesolve.pdf")
plt.show()