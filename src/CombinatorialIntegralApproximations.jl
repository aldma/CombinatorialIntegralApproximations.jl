module CombinatorialIntegralApproximations

using Logging

include("binary_approximation.jl")

include("combina_abstract_solver.jl")
include("combina_sur.jl")
include("combina_milp.jl")

end
