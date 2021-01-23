
# x is the matrix of K n-space points of dimension (K, n)
# L is the max degree monomial
function vandermat(x::Matrix{T}, L::Int) where {T}
    K, n = size(x)
    U = binomial(n + L, L)
    
    i = 2
    ii = ones(Int, n+1)
        
    Q = Matrix{T}(undef, K, U)
    r = zeros(T, U, U)
    r[1, 1] = one(T)
    
    Q[:, 1] .= one(T) / sqrt(K)
    
    @inbounds @views for d in 1:L
        for j in 1:n
            i_star = i
            for k in ii[j]:ii[n+1]
                Q[:, i] = x[:, j] .* Q[:, k]
                for l in 1:i-1
                    r[l, i] = Q[:, l]' * Q[:, i]
                    Q[:, i] = Q[:, i] - r[l, i] * Q[:, l]
                end
                r[i, i] = norm(Q[:, i])
                if isapprox(r[i, i], zero(T), atol=1e-12)
                    println("rank deficient")
                    # and leave column i un-normalized
                else
                    Q[:, i] ./= r[i, i]
                end
                i += 1
            end
            ii[j] = i_star
        end
        ii[n+1] = i - 1
    end
    return Q, r
end


# x is the (K, n) matrix of K points of n-dimensional space.
# f is the K vector of values of the funcion at the x locations
# M is the degree of the numerator and N is the degree of the denominator
function vrfit(x, f, M, N)
    _, n = size(x)
    Lmin, Lmax = extrema([M, N])
    VL, RL = vandermat(x, Lmax)
    U = binomial(n + Lmin, Lmin)

    VM, RM = M >= N ? (VL, RL) : (view(VL, :, 1:U), view(RL, :, 1:U))
    VN, RN = N >= M ? (VL, RL) : (view(VL, :, 1:U), view(RL, :, 1:U))
    
    fVN = Diagonal(f) * VN
    Z = pinv(VM) * fVN
    W = VM * Z - fVN
    G = svd(W)
    b = G.Vt[end, :]
    a = Z * b
    return VanderFit(VL, RL, M, N, a, b)
end


function veval(x::Vector{T}, c::Vector{T}, r::Matrix{T}, L::Int) where {T}
    n = length(x)
    U = binomial(n + L, L)

    i = 2
    ii = ones(Int, n+1)

    y = Vector{T}(undef, U)
    y[1] = one(T)
    
    @inbounds for d in 1:L
        for j in 1:n
            i_star = i
            for k in ii[j]:ii[n+1]
                y[i] = x[j] * y[k]
                for l in 1:i-1
                    y[i] = y[i] - r[l, i] * y[l]
                end
                y[i] = y[i] / r[i, i]
                i += 1
            end
            ii[j] = i_star
        end
        ii[n+1] = i - 1
    end
    return sum(c .* y)
end


function vrval(x::Vector{T}, vf::VanderFit{T}) where {T}
    n = veval(x, vf.a, vf.RL, vf.M)
    d = veval(x, vf.b, vf.RL, vf.N)
    return n / d
end
