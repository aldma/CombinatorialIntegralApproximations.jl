using CombinatorialIntegralApproximations
using Plots

include("../test/utils.jl")

b_rel_single = get_single_input()

nt = length(b_rel_single)
nc = 2
b_rel = zeros(nt, nc)
b_rel[:, 1] .= b_rel_single
b_rel[:, 2] .= 1 .- b_rel_single

t = collect(range(start = 0, stop = 240, length = nt + 1))

binapprox = BinApprox(t, b_rel)
s = CombinaSUR(binapprox)
status(s)
CombinatorialIntegralApproximations.setup!(s)
status(s)
solve!(s)
status(s)

b_bin = s.binapprox.b_bin

# TODO plotting function in utils

t_plt = t # [nt+1]
a_plt = [b_rel; b_rel[nt, :]']
b_plt = [b_bin; b_bin[nt, :]']
plot(t_plt, [a_plt[:, j] for j = 1:nc], seriestype = :scatter)
plot!(t_plt, [b_plt[:, j] for j = 1:nc], seriestype = :path, linetype = :steppost)
savefig("rounding_example_1.pdf")

s2 = CombinaMILP(binapprox)
solve!(s2)
b_bin2 = s2.binapprox.b_bin

b_tmp = b_bin - b_bin2
b_plt .= [b_tmp; b_tmp[nt, :]']
plot(t_plt, [a_plt[:, j] for j = 1:nc], seriestype = :scatter)
plot!(t_plt, [b_plt[:, j] for j = 1:nc], seriestype = :path, linetype = :steppost)
savefig("rounding_example_diff.pdf")