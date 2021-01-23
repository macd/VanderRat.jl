
# copy-pasta... don't want big dependencies
function chebpts(n)
    m = n - 1
    x = sinpi.([-m:2:m;] ./ 2m)
    return x
end;

function chebpts(n::Int, l::Float64, u::Float64)
    x = chebpts(n) * ((u - l) / 2.0)
    s = l - x[1]
    return x .+ s
end

