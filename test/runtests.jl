using VanderRat
using Test
using LinearAlgebra

include("tests.jl")

fun_a5_2(x) = log(2.25 - x'*x)
hyptan(x) = tanh(5*(x[1] − x[2]))
sexp(x) = exp(-x'*x / 1000)
nvarf(x) = exp(prod(x)) / prod([(y^2 - 1.44) for y in x])

# as the dimensionality goes up, so does the need for data (or to at least
# lower the order of M,N the approximant
@testset "VanderRat.jl" begin

    @test nvartest(nvarf, k=1000,   n=2) < 1e-11    # default N=20, M=20
    @test nvartest(nvarf, k=10000,  n=3) < 1e-10
    @test nvartest(nvarf, k=10000,  n=4, M=10, N=10) < 5e-2
    
    @test nvartest(hyptan, k=1000,  n=2) < 1e-6
    @test nvartest(fun_a5_2, k=1000, n=2) < 1e-6 # default N=20, M=20
    @test nvartest(sexp, k=10000,  n=2, M=30, N=20) < 1e-6 # default N=20, M=20
    @test doreduce() < 1e-4

end
