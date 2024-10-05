export AbstractCombinaSolver
export status, setup!, solve!

# TODO preprocessing
# TODO move solver_statuses out

abstract type AbstractCombinaSolver end
# to be always instantiated with properties
#   solver_status::Int
#   binapprox::BinApprox

function pre_setup!(s::S) where {S <: AbstractCombinaSolver}
    # check and allocate solution b_bin
    nt = s.binapprox.n_t
    nc = s.binapprox.n_c
    if isnothing(s.binapprox.b_bin)
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
    return s
end

function core_setup!(s::S) where {S <: AbstractCombinaSolver}
    return s
end

function post_setup!(s::S) where {S <: AbstractCombinaSolver}
    set_solver_status!(s, 1)
    return s
end

function setup!(s::S) where {S <: AbstractCombinaSolver}
    pre_setup!(s)
    core_setup!(s)
    post_setup!(s)
    return s
end

function solve!(s::S) where {S <: AbstractCombinaSolver}
    """
        Solve the combinatorial integral approximation problem.
    """
    if s.solver_status < 1
        setup!(s)
    end
    run!(s)
    set_solution!(s)
    set_solver_status!(s, 2)
    return s
end

function set_solver_status!(s::S, status::Int) where {S <: AbstractCombinaSolver}
    s.solver_status = status
    return s
end

function set_solution!(s::S) where {S <: AbstractCombinaSolver}
    return s
end

function run!(s::S) where {S <: AbstractCombinaSolver}
    error("Method not implemented...")
end

function status(s::S) where {S <: AbstractCombinaSolver}
    solver_statuses = Dict{Int, String}(
        -1 => "Exception",
        0 => "Built",
        1 => "Initialized",
        2 => "Solution found",
    )
    if s.solver_status in keys(solver_statuses)
        return solver_statuses[s.solver_status]
    else
        error("Unknown solver status, this should not happen.")
    end
end