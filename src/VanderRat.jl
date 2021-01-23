module VanderRat
using LinearAlgebra

export polyfitA, polyvalA, vandermat, vrfit, veval, vrval, degreduce

struct VanderFit{T}
    VL::Matrix{T}   # we only save this for resizing M and N
    RL::Matrix{T}   # recursion matrix
    M::Int          # order of numerator
    N::Int          # order of denominator
    a::Vector{T}    # numerator coefs for size M, N
    b::Vector{T}    # denominator coefs for size M, N
end

include("utils.jl")
include("vanderational.jl")
include("degreduce.jl")
include("poly.jl")

end
