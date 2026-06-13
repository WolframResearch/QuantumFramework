import json
def gen(n, m):
    # deterministic, identical across WL/Julia/Python: mix of H, S, CNOT
    ops=[]
    for i in range(m):
        r=i%3
        if r==0: ops.append(["H", i%n])
        elif r==1: ops.append(["S", (i*7+1)%n])
        else:
            a=i%n; b=(i+ (i%5) +1)%n
            if b==a: b=(a+1)%n
            ops.append(["CNOT", a, b])
    return ops
for (n,m) in ((100,2000),(500,10000)):
    json.dump(gen(n,m), open(f"/tmp/stab_ops_{n}.json","w"))
    print(f"wrote /tmp/stab_ops_{n}.json")
