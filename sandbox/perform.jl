using LinearAlgebra
using VanderRat
using BenchmarkTools

nvarf(x) = exp(prod(x)) / prod([(y^2 - 1.44) for y in x])

k=1000
n=2
M=20
N=20

x = 2*rand(k, n) .- 1.0
ff = nvarf.(eachrow(x))

@btime vf = VanderRat.vrfit($x, $ff, $M, $N);

#=
julia @v1.7> include("perform.jl")
  94.061 ms (53411 allocations: 436.18 MiB)

Doesn't seem too bad.  Should profile at some point to see
if there is anything low hanging to fix
=#
