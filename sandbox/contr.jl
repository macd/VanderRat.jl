using LinearAlgebra
using VanderRat
using PyPlot

# One of the functions from the paper
nvarf(x) = exp(prod(x)) / prod([(y^2 - 1.44) for y in x])

function meshgrid(x, y)
    xx = permutedims(repeat(x, 1, length(y)))
    yy = repeat(y, 1, length(x))
    return xx, yy
end

# square meshgrid. Takes x [low to high]
function meshgrid(x)
    return meshgrid(x, reverse(x))
end

function finst(;k=1000, n=2, M=20, N=20)
    x = 2*rand(k, n) .- 1.0
    ff = nvarf.(eachrow(x))
    vf = VanderRat.vrfit(x, ff, M, N)
    show_contrs(vf)
    readline()
    vff = VanderFit(vf)
    degreduce!(vff, 12, 9)
    show_contrs(vff)

    # the values at a random test point
    xx = rand(n)
    println(nvarf(xx), "  ", vrval(xx, vf), "   ", vrval(xx, vff))
end

function show_contrs(vf)
    f(x) = veval(x, vf.b, vf.RL, vf.N)
    g(x) = veval(x, vf.a, vf.RL, vf.M)
    r = [-1.6:0.1:1.6;]
    xx, yy = meshgrid(r)

    zz = [vec(xx) vec(yy)]
    
    ZF = f.(eachrow(zz))
    ZF = reshape(ZF, 33, 33)
    
    ZG = g.(eachrow(zz))
    ZG = reshape(ZG, 33, 33)
    
    contour(xx, yy, ZF, [0.0])
    contour(xx, yy, ZG, [0.0])
end
