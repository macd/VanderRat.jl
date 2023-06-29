module VanderRat
using LinearAlgebra

export VanderFit, vandermat, vrfit, veval, vrval, calcreduce, chebpts,
       degreduce!, polyfitA, polyvalA, polyfitD, polyvalD


# So VL is essentially the polynomial basis that was constructed for the K 
# points in n-space (from the coordinates of the sample points). It is used to 
# calculate the a and b vectors, so it is actually not needed for the estimation. 
# However, it is needed when we reduce the size of the polynomials so we can
# recalculate a and b.
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

# Copy constructor
function VanderFit(vf::VanderFit)
    VanderFit(vf.VL, vf.RL, vf.M, vf.N, copy(vf.a), copy(vf.b),
              vf.Mmax, vf.Nmax, vf.n, copy(vf.f))
end

include("utils.jl")
include("vanderational.jl")
include("degreduce.jl")
include("poly.jl")

end
