using QuantumClifford, BenchmarkTools
function genops(n,m)
    ops=Vector{Tuple{Symbol,Int,Int}}()
    for i in 0:m-1
        r=i%3
        if r==0; push!(ops,(:H, i%n, 0))
        elseif r==1; push!(ops,(:S, (i*7+1)%n, 0))
        else
            a=i%n; b=(i+(i%5)+1)%n
            if b==a; b=(a+1)%n; end
            push!(ops,(:CNOT,a,b))
        end
    end
    ops
end
function runsim(n,ops)
    s=MixedDestabilizer(one(Stabilizer,n))
    for op in ops
        if op[1]==:H; apply!(s, sHadamard(op[2]+1))
        elseif op[1]==:S; apply!(s, sPhase(op[2]+1))
        else; apply!(s, sCNOT(op[2]+1, op[3]+1)); end
    end
    s
end
for (n,m) in ((100,2000),(500,10000))
    ops=genops(n,m)
    runsim(n,ops)
    t=@belapsed runsim($n,$ops)
    println("RESULT|qc_jl_n$(n)_m$(m)_ms|", t*1000)
end
