export BinApprox

# TODO use keyword arguments and fill in
# TODO allow single control input (and fill in)

# TODO b_valid
# TODO b_adjacencies
# TODO dwell time tolerance
# TODO n max switches
# TODO min_up_times
# TODO min_down_times
# TODO max_up_times
# TODO total_max_up_times
# TODO b_bin_pre, set_b_bin_pre
# TODO set_valid_control_transitions
# TODO set_valid_controls_for_interval

mutable struct BinApprox
    t::AbstractVector
    b_rel::AbstractArray
    dt::AbstractVector
    n_t::Int
    n_c::Int
    eta::Union{Real, Nothing}
    b_bin::Union{AbstractArray, Nothing}
    cia_norm::Symbol
    reduce_problem_size_before_solve::Bool
    binary_threshold::Real
    clamped::Int
end

# constructor
function BinApprox(
    t::S,
    b_rel::V;
    binary_threshold::T = 1e-3,
    reduce_problem_size_before_solve::Bool = false,
    cia_norm::Symbol = :max_norm,
) where {T <: Real, S <: AbstractVector, V <: AbstractArray}
    if binary_threshold < 0
        error("Input binary_threshold must be nonnegative.")
    elseif binary_threshold >= 1 / 2
        error("Input binary_threshold must be less than 1/2.")
    end
    # set time points t
    if length(size(t)) > 1
        error("Input t must be a vector.")
    end
    n_t = length(t) - 1
    dt = similar(t, n_t)
    dt .= t[2:(n_t + 1)] .- t[1:n_t]
    if !all(dt .> 0)
        error("Values in t must be strictly increasing.")
    end
    # set relaxed binaries b_rel
    # b_rel [n_t x n_c]
    if size(b_rel, 1) != n_t
        b_rel = transpose(b_rel)
    end
    if size(b_rel, 1) != n_t
        error("One dimension of b_rel must be |t|-1.")
    end
    if !all((b_rel .>= 0) .& (b_rel .<= 1))
        error("All elements of the relaxed binary input must be 0 <= b <= 1.")
    end
    clamped_down = b_rel .< binary_threshold
    clamped_up = b_rel .> 1.0 - binary_threshold
    b_rel[clamped_down] .= 0
    b_rel[clamped_up] .= 1
    clamped = sum(clamped_down) + sum(clamped_up)

    # number of control intervals
    n_c = size(b_rel, 2)

    # check_sos1_constraint_fulfilled
    sums = sum(b_rel, dims = 2)
    tol = clamped * binary_threshold + (n_c + 1) * eps(eltype(sums))
    if any(sums .> 1.0 + tol) || any(sums .< 1.0 - tol)
        error("The sum of relaxed binary controls per time point must be one.")
    end

    BA = BinApprox(
        t,
        b_rel,
        dt,
        n_t,
        n_c,
        nothing,
        nothing,
        cia_norm,
        reduce_problem_size_before_solve,
        binary_threshold,
        clamped,
    )
    set_cia_norm!(BA, cia_norm)
    return BA
end

function set_cia_norm!(BA::BinApprox, cia_norm::Symbol)
    #=
        Specify the norm applied to || integral_t_0^t_f (b_rel-b_bin)*dt ||

        Allowed choices:

        - "max_norm": maximum over all time points and controls
        - "column_sum_norm": maximum over all time points of the sum (over all controls) of absolute accumulated values (b_rel-b_bin)*dt
        - "row_sum_norm": maximum over all controls of the sum (over all time points) of absolute accumulated values (b_rel-b_bin)*dt

        Usage::

            >>> from pycombina import BinApprox

            >>> t = [0, ..., 9.0, 9.5, 10.0]
            >>> b_rel = [[0.0      , ..., 0.558401, 0.558401, 0.558401],
            ...          [0.0      , ..., 0.0     , 0.0     , 0.0     ],
            ...          [1.0      , ..., 0.441599, 0.441599, 0.441599]])

            >>> binapprox = BinApprox(t, b_rel)

            >>> # Induce max_norm:
            >>> binapprox.set_cia_norm("max_norm")

        :param cia_norm: norm choice.
    =#
    available_norms = [:max_norm, :column_sum_norm, :row_sum_norm]
    if !(cia_norm in available_norms)
        error("cia_norm must be set one of the available norms: $(available_norms).")
    end
    BA.cia_norm = cia_norm
end

#=function determine_number_of_control_intervals!(P::BinApprox)
    P.n_t = length(P.t) - 1
end

function BinApprox.determine_number_of_control_intervals!(P::BinApprox)
    P.n_c
end

function BinApprox.compute_time_grid_from_time_points(P::BinApprox)
    P.n_c
end

function BinApprox._compute_dwell_time_tolerance(P::BinApprox)
    P.n_c
end

function BinApprox.set_b_bin(P::BinApprox)
    P.n_c
end

function BinApprox.set_eta(P::BinApprox)
    P.n_c
end
=#