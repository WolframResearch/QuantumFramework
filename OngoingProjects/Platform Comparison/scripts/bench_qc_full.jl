using QuantumClifford, BenchmarkTools, JSON
function runsim(n,ops)
    s=MixedDestabilizer(one(Stabilizer,n))
    for op in ops
        g=op[1]
        if g=="H"; apply!(s, sHadamard(op[2]+1))
        elseif g=="S"; apply!(s, sPhase(op[2]+1))
        else; apply!(s, sCNOT(op[2]+1, op[3]+1)); end
    end
    s
end
for (n,m) in ((100,2000),(500,10000),(1000,20000))
    ops=JSON.parsefile("/tmp/stab_ops_$(n).json")
    runsim(n,ops)
    t=@belapsed runsim($n,$ops)
    println("RESULT|qcjl_n$(n)|", round(t*1000,digits=3))
end
