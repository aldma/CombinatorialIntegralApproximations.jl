# CombinatorialIntegralApproximations.jl
__Solving binary approximation problems in Julia__

[![Build Status](https://github.com/aldma/CombinatorialIntegralApproximations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aldma/CombinatorialIntegralApproximations.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package provides modelling tools and solution algorithms for combinatorial approximation problems arising, for instance, in mixed-integer optimal control. Implemented algorithms include:

- Sum-Up and Rounding (SUR)
- Mixed-integer linear programming (MILP)

Check out the [CIA paper](https://doi.org/10.1007/s00186-011-0355-4) for a theoretical background and the [pycombina paper](https://doi.org/10.1016/j.ifacol.2020.12.1799) for an algorithmic overview.

Some algorithms rely on:
 - [JuMP](https://jump.dev/JuMP.jl/stable/) for modelling or reformulating problems,
 - [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl) for solving those problems numerically.

## Install

Use `]` to enter `pkg>` mode of Julia, then
```julia
pkg> add CombinatorialIntegralApproximations
```

## Documentation

Coming soon

## Credits

`CombinatorialIntegralApproximations.jl` is a pure Julia implementation of the software package [pycombina](https://github.com/adbuerger/pycombina), developed in python (with the `CombinaBnB` solver written in C++).

# Bug reports and discussions

Contributions are welcome in the form of issues notification or pull requests.
We recommend looking at already implemented algorithms to get inspiration on how to structure new ones.
If you think you found a bug, please open an [issue](https://github.com/aldma/CombinatorialIntegralApproximations.jl/issues).
Focused suggestions and requests can also be opened as issues.
Before opening a pull request, it is recommended to start an issue or a discussion on the topic.
