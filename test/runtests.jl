using Combina
using Test

@testset "Combina.jl" begin
    # Write your tests here.
    @test Combina.greet_your_package_name() == "Hello Combina!"
    @test Combina.greet_your_package_name() != "Hello world!"
end
