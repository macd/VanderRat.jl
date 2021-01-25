using VanderRat
using Test
using LinearAlgebra

include("tests.jl")

@testset "VanderRat.jl" begin
    @test nvartest(k=1000, n=2) < 1e-13
    @test nvartest(k=3000, n=3) < 1e-5
    @test doreduce() < 1e-4
end
