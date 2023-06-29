
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

# Random sample points and random test points. What's the worry?
function nvartest(;k=1000, n=2, M=20, N=20)
    x = rand(k, n)
    f = nvarf.(eachrow(x))
    vf = VanderRat.vrfit(x, f, M, N)
    # test at 5 random locations
    xx = rand(5, n)
    return maximum(abs.(nvarf.(eachrow(xx)) - vrval.(eachrow(xx), (vf,))))
end

# test that degreduce! actually does something
function doreduce()
    k=10000; n=2; M=20; N=20
    x = rand(k, n)
    f = nvarf.(eachrow(x))
    vf = VanderRat.vrfit(x, f, M, N)
    vff = VanderFit(vf)
    M, N = calcreduce(vff)
    println("M: ", M, "  N: ", N)

    xx = rand(2)

    # println("\nOriginal 20, 20")
    # println(abs(nvarf(xx) - vrval(xx, vff)))
    
    # println("\nreducing to 12, 9")
    # degreduce!(vff, 12, 9)
    # println(abs(nvarf(xx) - vrval(xx, vff)))

    # println("\nreducing to 8, 6")
    degreduce!(vff, 8, 6)
    #println(abs(nvarf(xx) - vrval(xx, vff)))
    
    return abs(nvarf(xx) - vrval(xx, vff))
end
