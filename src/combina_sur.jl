export CombinaSUR
export status, setup!, solve!

mutable struct CombinaSUR <: AbstractCombinaSolver
    solver_status::Int
    binapprox::BinApprox
end

function CombinaSUR(binapprox::BinApprox)
    solver_status = 0
    return CombinaSUR(solver_status, binapprox)
end

function run!(s::CombinaSUR)
    @info("Running Sum-Up-Rounding...")
    start_time = time()

    nt = s.binapprox.n_t
    nc = s.binapprox.n_c
    b_bin = s.binapprox.b_bin
    eta = s.binapprox.eta
    eta_i = zeros(nc)

    for i = 1:nt
        b_active = 1
        for j = 1:nc
            eta_i[j] += s.binapprox.b_rel[i, j] * s.binapprox.dt[i]
            if eta_i[j] > eta_i[b_active]
                b_active = j
            end
        end
        b_bin[i, b_active] = 1
        eta_i[b_active] -= s.binapprox.dt[i]
        eta = max(eta, maximum(abs.(eta_i)))
    end

    #s.binapprox.b_bin = b_bin
    #s.binapprox.eta = eta
    runtime = time() - start_time
    @info("done")
    @info("    Best solution: $(eta)")
    @info("    Total runtime: $(runtime) s")
    return s
end