# Phase-resolved QuantumClifford.jl benchmark. Companion to bench_phases_qf.wls /
# bench_phases_stim.py; shared schema  RESULT|qcjl|<phase>|n=<n>|k=<k>|<ms>.
#
# Phases:
#   ingest       op tuples -> typed gate vector (sHadamard/sPhase/sCNOT objects).
#                Native Julia, no foreign-function boundary.
#   simulate     apply! fold of the pre-built gate vector over MixedDestabilizer.
#   materialize  stab_to_gf2 on the full tableau view + phases: dense Bool matrix.
#   measure1     projective Z measurement on qubit 1 of the evolved state.
#
# The JSON parse happens outside all timed regions, as on every engine.
using QuantumClifford, BenchmarkTools, JSON

function build(ops)
    gates = Vector{QuantumClifford.AbstractOperation}()
    sizehint!(gates, length(ops))
    for op in ops
        g = op[1]
        if g == "H"; push!(gates, sHadamard(op[2] + 1))
        elseif g == "S"; push!(gates, sPhase(op[2] + 1))
        else; push!(gates, sCNOT(op[2] + 1, op[3] + 1)); end
    end
    gates
end

function runsim(n, gates)
    s = MixedDestabilizer(one(Stabilizer, n))
    for g in gates
        apply!(s, g)
    end
    s
end

pr(phase, n, k, ms) = println("RESULT|qcjl|", phase, "|n=", n, "|k=", k, "|", round(ms, digits=4))

for (n, m) in ((100, 2000), (500, 10000), (1000, 20000))
    ops = JSON.parsefile(joinpath(@__DIR__, "stab_ops_$(n).json"))
    k = length(ops)
    @assert k == m
    gates = build(ops)
    s = runsim(n, gates)  # warm compile + evolved state for materialize/measure

    pr("ingest", n, k, 1000 * @belapsed build($ops))
    pr("simulate", n, k, 1000 * @belapsed runsim($n, $gates))
    pr("materialize", n, k, 1000 * @belapsed (stab_to_gf2(stabilizerview($s)), phases(stabilizerview($s))))
    pr("measure1", n, k, 1000 * @belapsed projectZ!(copy($s), 1))
end
println("DONE|qcjl")
