# VanderRat

This package is a Julia implementation of some of the algorithms presented 
paper [Practical algorithms for multivariate rational approximation](https://doi.org/10.1016/j.cpc.2020.107663) (It does not implement the Algorithm 4.1)

## The algorithm

## Example

```julia
using LinearAlgebra
using VanderRat
using BenchmarkTools

# Example A.5.1 in [1]
nvarf(x) = exp(prod(x)) / prod([(y^2 - 1.44) for y in x])

k=1000
n=2
M=20
N=20

x = 2*rand(k, n) .- 1.0
ff = nvarf.(eachrow(x))

@btime vf = VanderRat.vrfit($x, $ff, $M, $N);

#julia 1.14> @btime vf = VanderRat.vrfit($x, $ff, $M, $N);
#  60.621 ms (3762 allocations: 30.12 MiB)

# Random sample points and random test points. What's the worry?
# when n=3, error is a min at M=N=11 (only tested M=N) Error gets
# larger for higher orders eventually triggering a rank deficient error
# because number of sample points is less than the number of basis functions
# We now give a warning when that happens.
function nvartest(func; llim=-1.0, ulim=1.0, k=1000, n=2, M=20, N=20)
    x = llim .+ (ulim - llim)*rand(k, n)
    f = func.(eachrow(x))
    vf = VanderRat.vrfit(x, f, M, N)
    # test at 5 random locations
    xx = llim .+ (ulim - llim)*rand(5, n)
    return maximum(abs.(func.(eachrow(xx)) - vrval.(eachrow(xx), (vf,))))
end

#julia 1.14> nvartest(nvarf)
#5.995204332975845e-15
```

## Notes

This package also contains an implementation of [Vandermonde with
Arnoldi](https://epubs.siam.org/doi/pdf/10.1137/19M130100X)

## References

[1] [Practical algorithms for multivariate rational approximation](https://doi.org/10.1016/j.cpc.2020.107663)

[2]


