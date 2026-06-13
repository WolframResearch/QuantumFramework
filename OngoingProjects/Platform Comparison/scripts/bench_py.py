import time, numpy as np
def timeit(f, reps=3):
    best=1e9
    for _ in range(reps):
        t=time.perf_counter(); f(); best=min(best, time.perf_counter()-t)
    return best

# ---- Benchmark A: statevector sim, fixed circuit: (H all; CNOT chain; RZ(0.3) all) x3 ----
def bench_statevector(n, reps_layers=3):
    res={}
    # Qiskit Aer
    try:
        from qiskit import QuantumCircuit
        from qiskit_aer import AerSimulator
        def run():
            qc=QuantumCircuit(n)
            for _ in range(reps_layers):
                for q in range(n): qc.h(q)
                for q in range(n-1): qc.cx(q,q+1)
                for q in range(n): qc.rz(0.3,q)
            qc.save_statevector()
            sim=AerSimulator(method='statevector')
            sim.run(qc).result().data()['statevector']
        res['qiskit_aer']=timeit(run)
    except Exception as e: res['qiskit_aer']=f'ERR {e}'
    # Cirq
    try:
        import cirq
        qs=cirq.LineQubit.range(n)
        def run():
            c=cirq.Circuit()
            for _ in range(reps_layers):
                c.append(cirq.H(q) for q in qs)
                c.append(cirq.CNOT(qs[i],qs[i+1]) for i in range(n-1))
                c.append(cirq.rz(0.3)(q) for q in qs)
            cirq.Simulator().simulate(c)
        res['cirq']=timeit(run)
    except Exception as e: res['cirq']=f'ERR {e}'
    return res

# ---- Benchmark C(stim part): stabilizer sim, m random Clifford gates on n qubits ----
def bench_stim(n, m, seed=1):
    import stim
    rng=np.random.default_rng(seed)
    ops=[]
    for _ in range(m):
        if rng.random()<0.5:
            ops.append(('H', int(rng.integers(n))))
        else:
            a,b=rng.integers(n),rng.integers(n)
            while b==a: b=rng.integers(n)
            ops.append(('CNOT', int(a), int(b)))
    def run():
        s=stim.TableauSimulator()
        for op in ops:
            if op[0]=='H': s.h(op[1])
            else: s.cnot(op[1],op[2])
    return timeit(run), ops

# ---- Benchmark B: amplitude-damped qubit, QuTiP mesolve; return final excited-state pop + time ----
def bench_qutip_lindblad():
    import qutip as qt
    # H = (omega/2) sigma_x drive; collapse sqrt(gamma) sigma_minus; start excited
    omega=2*np.pi*1.0; gamma=0.3
    H=0.5*omega*qt.sigmax()
    c_ops=[np.sqrt(gamma)*qt.sigmam()]
    psi0=qt.basis(2,1)  # excited
    tlist=np.linspace(0,2,201)
    res={}
    def run():
        r=qt.mesolve(H,psi0,tlist,c_ops=c_ops,e_ops=[qt.num(2)])
        return r.expect[0][-1]
    res['time']=timeit(run)
    res['final_excited_pop']=float(run())
    return res

if __name__=='__main__':
    import json
    out={}
    for n in (12,16):
        out[f'statevector_n{n}']=bench_statevector(n)
    for (n,m) in ((100,2000),(500,10000)):
        t,_=bench_stim(n,m); out[f'stim_n{n}_m{m}']=t
    out['qutip_lindblad']=bench_qutip_lindblad()
    print(json.dumps(out, indent=2, default=str))
