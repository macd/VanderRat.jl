
# function testvanderrat(d=300, n=3, M=3, N=3)
#     x = 2.0 * rand(d, n) .- 1.0
#     f = exp.(prod(x[i, :]) for i in 1:d)
#     vf = vrfit(x, f, M, N)

#     xx = rand(n)
#     println(exp(prod(xx)), "   ", vrval(xx, vf))
# end


# # two dimensionl with chebyshev tensor product grid
# function tvdr(p=5, M=2, N=2)
#     pts = chebpts(p)
#     d = p ^ 2
#     x = permutedims(reduce(hcat, [[j, i] for i in pts for j in pts]))
#     f = exp.(prod(x[i, :]) for i in 1:d)
#     vf = vrfit(x, f, M, N)

#     #xx = rand(2)
#     xx = [1.0, 1.0]
#     println(exp(prod(xx)), "   ", vrval(xx, vf))
# end

# One of the functions from the paper
nvarf(x) = exp(prod(x)) / prod([(y^2 - 1.44) for y in x])

function nvartest(;k=1000, n=2, M=20, N=20)
    x = rand(k, n)
    f = nvarf.(eachrow(x))
    vf = VanderRat.vrfit(x, f, M, N)
    xx = rand(n)
    return abs(nvarf(xx) - vrval(xx, vf))
end

function doreduce(;k=1000, n=2, M=20, N=20, eta=1e-12)
    x = rand(k, n)
    f = nvarf.(eachrow(x))
    vf = VanderRat.vrfit(x, f, M, N)
    return degreduce(vf.VL, f, M, N, n, eta=eta)
end
