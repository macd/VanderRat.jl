module VanderRat
using LinearAlgebra

export VanderFit, vandermat, vrfit, veval, vrval, calcreduce, degreduce!, polyfitA, polyvalA

mutable struct VanderFit{T}
    VL::Matrix{T}   # we only save this for resizing M and N
    RL::Matrix{T}   # recursion matrix
    M::Int          # (current) order of numerator
    N::Int          # (current) order of denominator
    a::Vector{T}    # numerator coefs for size M, N
    b::Vector{T}    # denominator coefs for size M, N
    Mmax::Int       # M of first created VanderFit 
    Nmax::Int       # N of first created VanderFit 
    n::Int          # vector dimension. Needed for resizing
    f::Vector{T}    # function values.  Needed for resizing
end

function VanderFit(vf::VanderFit)
    VanderFit(vf.VL, vf.RL, vf.M, vf.N, copy(vf.a), copy(vf.b), vf.Mmax, vf.Nmax, vf.n, copy(vf.f))
end

include("utils.jl")
include("vanderational.jl")
include("degreduce.jl")
include("poly.jl")

end
