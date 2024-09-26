export CombinaSUR
export status, solve!

# TODO preprocessing
# TODO move solver_statuses out

mutable struct CombinaSUR
    sur_status::Int
    binapprox::BinApprox
end

# constructor
function CombinaSUR(binapprox::BinApprox)
    sur_status = 0
    return CombinaSUR(sur_status, binapprox)
end

function setup_sur!(s::CombinaSUR)
    # check and allocate solution b_bin
    if isnothing(s.binapprox.b_bin)
        nt = s.binapprox.n_t
        nc = s.binapprox.n_c
        s.binapprox.b_bin = zeros(nt, nc)
    else
        if size(s.binapprox.b_bin, 1) != nt
            error("b_bin must be either nothing or an array with a suitable number of rows.")
        end
        if size(s.binapprox.b_bin, 2) != nc
            error("b_bin must be either nothing or an array with a suitable number of columns.")
        end
    end
    # check and allocate solution eta
    if isnothing(s.binapprox.eta)
        s.binapprox.eta = -1
    else
        if length(s.binapprox.eta) != 1
            error("eta must be either nothing or a scalar.")
        end
    end
    # set SUR status
    set_sur_status!(s, 1)
    return s
end

function solve!(s::CombinaSUR)
    """
        Solve the combinatorial integral approximation problem.
    """
    if s.sur_status < 1
        setup_sur!(s)
    end
    run_sur!(s)
    set_solution!(s)
    set_sur_status!(s, 2)
    return s
end

function set_sur_status!(s::CombinaSUR, status::Int)
    s.sur_status = status
    return s
end

function set_solution!(s::CombinaSUR)
    return s
end

function run_sur!(s::CombinaSUR)
    @info("Running Sum-up-rounding...")
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

function status(s::CombinaSUR)
    solver_statuses = Dict{Int, String}(
        -1 => "Exception",
        0 => "Built",
        1 => "Initialized",
        2 => "Solution found",
    )
    if s.sur_status in keys(solver_statuses)
        return solver_statuses[s.sur_status]
    else
        error("Unknown solver status, this should not happen.")
    end
end