using CombinatorialIntegralApproximations
using Test

@testset "Input test" begin
    @testset "Test single control ok" begin
        t = [0, 1, 3, 4, 7]
        nt = length(t) - 1
        nc = 3
        b_rel = zeros(nt, nc)
        b_rel[1, :] .= [0.1, 0.3, 0.6]
        b_rel[2, :] .= [0.2, 0.5, 0.3]
        b_rel[3, :] .= [0.4, 0.2, 0.4]
        b_rel[4, :] .= [0.1, 0.7, 0.2]
        BA = BinApprox(t, b_rel)
        @test BA.n_t == nt
        @test BA.n_c == nc
        @test BA.t == t
        @test BA.b_rel == b_rel
        @test BA.dt == diff(t)
        @test isnothing(BA.eta)
        @test isnothing(BA.b_bin)
        @test BA.cia_norm == :max_norm
    end

    @testset "Test single control sos1 violated" begin
        t = [0, 1, 2, 3]
        b_rel = [0.1, 0.3, 0.2]
        msg = str -> occursin("sum of relaxed binary controls", str)
        @test_throws msg BinApprox(t, b_rel)
    end

    @testset "Test single control invalid dimensions" begin
        t = [0, 1, 2, 3]
        b_rel = [0.1, 0.3, 0.2, 0.5]
        msg = str -> occursin("dimension of b_rel must", str)
        @test_throws msg BinApprox(t, b_rel)
    end

    @testset "Test single control t not increasing" begin
        t = [0, 2, 1, 3]
        b_rel = [0.1, 0.3, 0.2]
        msg = str -> occursin("must be strictly increasing", str)
        @test_throws msg BinApprox(t, b_rel)
    end

    @testset "Test single control b_rel not relaxed binary solution" begin
        t = [0, 1, 2, 3]
        b_rel = [0.1, 0.3, 1.2]
        msg = str -> occursin("relaxed binary input must be", str)
        @test_throws msg BinApprox(t, b_rel)
    end

    @testset "Test manual extend sos1 fulfilled" begin
        for run_idx = 1:50
            nt = Int(round(10 + rand() * (1000 - 10)))
            nc = Int(round(1 + rand() * (100 - 1)))
            t = sort(rand(nt + 1) .* nt)
            b_rel = zeros(nt, nc)
            for j = 1:nt
                for i = 1:(nc - 1)
                    b_rel[j, i] = (1 .- sum(b_rel[j, 1:(i - 1)])) .* rand()
                end
                b_rel[j, nc] = 1 .- sum(b_rel[j, :])
            end
            @assert maximum(abs.(sum(b_rel, dims = 2) .- 1)) == 0
            BA = BinApprox(t, b_rel)
            @test BA.n_t == nt
            @test BA.n_c == nc
        end
    end

    @testset "Test manual scale sos1 fulfilled" begin
        for run_idx = 1:50
            nt = Int(round(10 + rand() * (1000 - 10)))
            nc = Int(round(2 + rand() * (100 - 2)))
            t = sort(rand(nt + 1) .* nt)
            b_rel = rand(nt, nc) ./ nc
            for j = 1:nt
                b_rel[j, :] .= b_rel[j, :] ./ sum(b_rel[j, :])
            end
            T = eltype(b_rel)
            @assert maximum(abs.(sum(b_rel, dims = 2) .- 1)) < 10 * eps(T)
            BA = BinApprox(t, b_rel)
            @test BA.n_t == nt
            @test BA.n_c == nc
        end
    end

    @testset "Test multiple controls sos1 fulfilled manual" begin
        t = [0, 1, 2, 3]
        nt = length(t) - 1
        nc = 2
        b_rel = zeros(nt, nc)
        b_rel[:, 1] .= [0.1, 0.3, 0.3]
        b_rel[:, 2] .= [0.9, 0.7, 0.7]
        @assert maximum(abs.(sum(b_rel, dims = 2) .- 1)) == 0
        BA = BinApprox(t, b_rel)
        @test BA.n_t == nt
        @test BA.n_c == nc
        t = [0, 2, 1, 3]
        @assert nt == length(t) - 1
        msg = str -> occursin("must be strictly increasing", str)
        @test_throws msg BinApprox(t, b_rel)
    end

    @testset "Test multiple controls sos1 violated" begin
        t = [0, 1, 2, 3]
        nt = length(t) - 1
        nc = 2
        b_rel = zeros(nt, nc)
        b_rel[:, 1] .= [0.1, 0.3, 0.3]
        b_rel[:, 2] .= [0.3, 0.4, 0.5]
        msg = str -> occursin("sum of relaxed binary controls", str)
        @test_throws msg BinApprox(t, b_rel)
    end
end

include("utils.jl")

@testset "Rounding test" begin
    b_rel_single = get_single_input()
    nt = length(b_rel_single)
    nc = 2
    b_rel = zeros(nt, nc)
    b_rel[:, 1] .= b_rel_single
    b_rel[:, 2] .= 1 .- b_rel_single
    t = collect(range(start = 0, stop = 240, length = nt + 1))
    binapprox = BinApprox(t, b_rel)
    @testset "SUR" begin
        s = CombinaSUR(binapprox)
        @test s.sur_status == 0
        @test occursin("Built", status(s))
        @test isnothing(binapprox.b_bin)
        @test isnothing(binapprox.eta)
        CombinatorialIntegralApproximations.setup_sur!(s)
        @test s.sur_status == 1
        @test occursin("Initialized", status(s))
        @test size(binapprox.b_bin, 1) == nt
        @test size(binapprox.b_bin, 2) == nc
        @test length(binapprox.eta) == 1
        @test size(s.binapprox.b_bin, 1) == nt
        @test size(s.binapprox.b_bin, 2) == nc
        @test length(s.binapprox.eta) == 1
        solve!(s)
        @test s.sur_status == 2
        @test occursin("Solution found", status(s))
        @test size(binapprox.b_bin, 1) == nt
        @test size(binapprox.b_bin, 2) == nc
        @test length(binapprox.eta) == 1
        @test size(s.binapprox.b_bin, 1) == nt
        @test size(s.binapprox.b_bin, 2) == nc
        @test length(s.binapprox.eta) == 1
        # TODO check values too
    end
end