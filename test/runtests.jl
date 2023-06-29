using VanderRat
using Test
using LinearAlgebra

include("tests.jl")

@testset "VanderRat.jl" begin
    @test nvartest(k=1000,  n=2) < 1e-12      # default N=20, M=20
    @test nvartest(k=10000,  n=3) < 1e-10
    @test nvartest(k=10000, n=4, M=5, N=5) < 1e-3
    @test doreduce() < 1e-4
end
