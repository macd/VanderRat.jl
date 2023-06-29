using LinearAlgebra
using Random
using VanderRat
using PyPlot

# One of the functions from the paper
nvarf(x) = exp(prod(x)) / prod([(y^2 - 1.44) for y in x])

# assume low to high in both x and y
function meshgrid(x, y)
    xx = repeat(x', length(y), 1)
    yy = repeat(reverse(y), 1, length(x))
    return xx, yy
end

# square meshgrid. Takes x [low to high]
function meshgrid(x)
    return meshgrid(x, x)
end

function finst(;k=1000, n=2, M=20, N=20, seed=nothing)
    if seed === nothing
        seed = rand(0:99999)
    end
    Random.seed!(seed)
    x = 2*rand(k, n) .- 1.0
    ff = nvarf.(eachrow(x))
    vf = VanderRat.vrfit(x, ff, M, N)
    show_contrs(vf, seed)
    
    # readline()
    # vff = VanderFit(vf)
    # degreduce!(vff, 12, 9)
    # show_contrs(vff)

    # # the values at a random test point
    # xx = rand(n)
    # println(nvarf(xx), "  ", vrval(xx, vf), "   ", vrval(xx, vff))
end



function show_contrs(vf, seed)
    f(x) = veval(x, vf.b, vf.RL, vf.N)
    g(x) = veval(x, vf.a, vf.RL, vf.M)
    r = [-1.4:0.005:1.4;]
    #r = [-1.8:0.005:1.8;] .+ 2.0
    #r = [-4:0.005:4;]
    n = length(r)
    xx, yy = meshgrid(r)

    zz = [vec(xx) vec(yy)]
    
    ZF = f.(vec(eachrow(zz)))
    ZF = reshape(ZF, n, n)
    
    ZG = g.(eachrow(zz))
    ZG = reshape(ZG, n, n)
    
    contour(xx, yy, ZF, [0.0], colors="r")
    readline()
    contour(xx, yy, ZG, [0.0], colors="g", alpha=0.5)
    title("Numerator and denominator zeros with seed $seed")
    xlabel("Red: denominator     Green: numerator     Brown: overlapped")
end


function fcheb(;k=32, n=2, M=20, N=20)
    x = VanderRat.chebpts(k)
    xx, yy = meshgrid(x)
    zz = [vec(xx) vec(yy)]
    ff = nvarf.(eachrow(zz))
    vf = VanderRat.vrfit(zz, ff, M, N)
    println(typeof(vf.b), "  ", typeof(vf.a))
    show_contrs(vf, k)
    
    readline()
    vff = VanderFit(vf)
    degreduce!(vff, 12, 9)
    show_contrs(vff, k)

    # the values at a random test point
    xx = rand(n)
    println(nvarf(xx), "  ", vrval(xx, vf), "   ", vrval(xx, vff))
end
