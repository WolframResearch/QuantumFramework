import json, time, stim
def timeit(f, reps=5):
    best=1e9
    for _ in range(reps):
        t=time.perf_counter(); f(); best=min(best,time.perf_counter()-t)
    return best
for (n,) in [(100,),(500,)]:
    ops=json.load(open(f"/tmp/stab_ops_{n}.json"))
    def run():
        s=stim.TableauSimulator()
        for op in ops:
            if op[0]=='H': s.h(op[1])
            elif op[0]=='S': s.s(op[1])
            else: s.cnot(op[1],op[2])
    print(f"stim n={n} m={len(ops)}: {timeit(run)*1000:.3f} ms")
