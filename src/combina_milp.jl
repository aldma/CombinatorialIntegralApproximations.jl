export CombinaMILP
export status, setup!, solve!

# TODO preprocessing

using JuMP
import HiGHS

mutable struct CombinaMILP <: AbstractCombinaSolver
    solver_status::Int
    binapprox::BinApprox
    model::GenericModel
end

# constructor
function CombinaMILP(binapprox::BinApprox)
    solver_status = 0
    model = Model(HiGHS.Optimizer)
    return CombinaMILP(solver_status, binapprox, model)
end

function core_setup!(s::CombinaMILP)
    build_jump_model!(s)
    return s
end

function build_jump_model!(s::CombinaMILP)
    n_t = s.binapprox.n_t
    n_c = s.binapprox.n_c
    b_rel = s.binapprox.b_rel
    dt = s.binapprox.dt
    n_max_switches = nothing # TODO add in BinApprox

    # TODO expose options
    set_attribute(s.model, "presolve", "on")
    set_attribute(s.model, "time_limit", 3600.0)
    set_attribute(s.model, "log_to_console", true)

    # setup model variables
    @variable(s.model, eta_sym)
    @variable(s.model, b_bin_sym[1:n_t, 1:n_c], Bin)
    if !isnothing(n_max_switches)
        @variable(s.model, s_sym[1:(n_t - 1), 1:n_c])
    end

    # setup objective
    @objective(s.model, Min, eta_sym)

    # setup constraints
    #   sos1 constraints
    @constraint(s.model, sos1eq[k = 1:n_t], sum(b_bin_sym[k, i] for i = 1:n_c) == 1.0)
    #   approximation inequalities
    @constraint(
        s.model,
        cia_plus[i = 1:n_c, k = 1:n_t],
        eta_sym >= sum((b_rel[j, i] - b_bin_sym[j, i]) * dt[j] for j = 1:k)
    )
    @constraint(
        s.model,
        cia_mnus[i = 1:n_c, k = 1:n_t],
        eta_sym >= -sum((b_rel[j, i] - b_bin_sym[j, i]) * dt[j] for j = 1:k)
    )
    # maximum number of switches
    if !isnothing(n_max_switches)
        @constraint(
            s.model,
            s_def1[i = 1:n_c, k = 1:(n_t - 1)],
            s_sym[k, i] >= b_bin_sym[k + 1, i] - b_bin_sym[k, i]
        )
        @constraint(
            s.model,
            s_def2[i = 1:n_c, k = 1:(n_t - 1)],
            s_sym[k, i] >= b_bin_sym[k, i] - b_bin_sym[k + 1, i]
        )
        @constraint(
            s.model,
            s_def3[i = 1:n_c, k = 1:(n_t - 1)],
            s_sym[k, i] <= b_bin_sym[k + 1, i] + b_bin_sym[k, i]
        )
        @constraint(
            s.model,
            s_def4[i = 1:n_c, k = 1:(n_t - 1)],
            s_sym[k, i] <= 2 - b_bin_sym[k + 1, i] - b_bin_sym[k, i]
        )
        # direct way
        @constraint(
            s.model,
            maxsw[i = 1:n_c],
            sum(s_sym[k, i] for k = 1:(n_t - 1)) <= n_max_switches[i]
        )
        # TODO pycombina way?
    end
    return s
end

function run!(s::CombinaMILP)
    @info("Running CombinaMILP...")

    start_time = time()
    optimize!(s.model)
    runtime = time() - start_time

    if !has_values(s.model)
        @error("Problem not really solved..")
    end
    b_bin_sym = s.model.obj_dict[:b_bin_sym]
    eta_sym = s.model.obj_dict[:eta_sym]
    b_bin = value.(b_bin_sym)
    eta = value.(eta_sym)
    # honor bounds
    b_bin .= max.(0, min.(b_bin, 1))
    eta = max(0, eta)

    @info("done")
    @info("    Best solution: $(eta)")
    @info("    Total runtime: $(runtime) s")

    s.binapprox.b_bin .= b_bin
    s.binapprox.eta = eta
    return s
end
