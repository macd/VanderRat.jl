
# copy-pasta... don't want big dependencies
function chebpts(n, ::Type{T}=Float64) where {T <: AbstractFloat}
    m = n - 1
    x = sinpi.([-m:2:m;] ./ 2m)
    return x
end;

function chebpts(n::Int, l::T, u::T) where {T <: AbstractFloat}
    x = chebpts(n) * ((u - l) / T(2))
    s = l - x[1]
    return x .+ s
end

