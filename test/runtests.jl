using VanderRat
using Test
using LinearAlgebra

include("tests.jl")

@testset "VanderRat.jl" begin
    @test bivartest()
end
