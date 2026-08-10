import qutip as qt
import numpy as np

ntraj = 100
epsilon1 = 1.0
epsilon2 = 1.0
g = 0.1

sx1 = qt.sigmax() & qt.qeye(2)
sx2 = qt.qeye(2) & qt.sigmax()
sz1 = qt.sigmaz() & qt.qeye(2)
sz2 = qt.qeye(2) & qt.sigmaz()

H = (epsilon1 / 2) * sz1 + (epsilon2 / 2) * sz2 + g * sx1 * sx2

sm1 = qt.sigmam() & qt.qeye(2)
sm2 = qt.qeye(2) & qt.sigmam()
gam = 0.1  # dissipation rate
c_ops = [np.sqrt(gam) * sm1, np.sqrt(gam) * sm2]
psi0 = qt.basis(2, 0) & qt.basis(2, 1)
tlist = np.linspace(0, 40, 100)

result_mc = qt.mcsolve(H, psi0, tlist, c_ops, e_ops=[sz1, sz2], ntraj=ntraj,
                       options={"map": "parallel", "keep_runs_results": True})
result_me = qt.mesolve(H, psi0, tlist, c_ops, e_ops=[sz1, sz2])

import matplotlib.pyplot as plt
from mpl_setup import * 
plt.figure()
plt.plot(tlist, result_me.expect[0],'-', label=r'mesolve')
plt.plot(tlist, result_mc.average_expect[0],'--',alpha=1,label=r'mcsolve average')
plt.fill_between(
    tlist,
    (result_mc.average_expect[0] - result_mc.std_expect[0] / np.sqrt(ntraj)) ,
    (result_mc.average_expect[0] + result_mc.std_expect[0] / np.sqrt(ntraj)),
    alpha=0.5, color="ORANGE"
)
plt.plot(tlist, result_mc.runs_expect[0][4],'-.',label=r'mcsolve 1 run')
plt.plot(tlist, result_mc.runs_expect[0][7],':',label=r'mcsolve 1 run')

plt.xlabel('Time', fontsize=18)
plt.ylabel(r'$\langle \sigma_z^{(1)} \rangle$', fontsize=18)
plt.legend()
plt.savefig("mcsolve.pdf")
plt.show()
