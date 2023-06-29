
# This is algorithm 3.1 from the paper.
# x is the matrix of K by n (K:number of samples, n: dimesionality of space)
# set of points.
# L is the max degree monomial for the approximation
function vandermat(x::AbstractMatrix{T}, L::Int) where {T}
    K, n = size(x)
    U = binomial(n + L, L)
    
    i = 2
    I = ones(Int, n+1)
        
    Q = Matrix{T}(undef, K, U)
    r = zeros(T, U, U)
    r[1, 1] = one(T)
    
    Q[:, 1] .= one(T) / sqrt(K)
    q = zeros(T, K)
    s = zeros(T, K+1)
    @inbounds @views for d in 1:L
        for j in 1:n
            i_star = i
            for k in I[j]:I[n+1]
                q .= x[:, j] .* Q[:, k]

                # My original version
                # for l in 1:i-1
                #     r[l, i] = Q[:, l]' * q
                #     q .= q - r[l, i] * Q[:, l]
                # end

                # 2 rounds of Orthogonalization
                for t in 1:2
                    s[1:i-1] .= Q[:, 1:i-1]' * q
                    q .= q - Q[:, 1:i-1] * s[1:i-1]
                    r[1:i-1, i] .+= s[1:i-1]
                end
                
                r[i, i] = norm(q)
                if isapprox(r[i, i], zero(T), atol=1e-12)
                    error("vandermat: x is near rank deficient")
                end
                Q[:, i] .= q ./ r[i, i]
                i += 1
            end
            I[j] = i_star
        end
        I[n+1] = i - 1
    end
    return Q, r
end

# So this top level algo does not have pseudo code in the paper but is
# described in "Section 3.1 Basic Algorithm"
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
    
    # from text between Eqn. 3 and Eqn 4
    Z = pinv(VM) * fVN
    
    # text just after Eqn 4
    W = VM * Z - fVN
    G = svd(W)
    b = G.Vt[end, :]
    a = Z * b
    return VanderFit(VL, RL, M, N, a, b, M, N, n, f)
end

# Algorithm 3.2 for evaluating the polynomials. It is sort of a generalized
# Clenshaw algorithm.
function veval(x::AbstractVector{T}, c::W, r::AbstractMatrix{T}, L::Int) where {T, W <: AbstractVector{T}}
    n = length(x)
    U = binomial(n + L, L)

    i = 2
    I = ones(Int, n+1)

    y = Vector{T}(undef, U)
    y[1] = one(T)
    
    @inbounds for d in 1:L
        for j in 1:n
            i_star = i
            for k in I[j]:I[n+1]
                y[i] = x[j] * y[k]
                for l in 1:i-1
                    y[i] = y[i] - r[l, i] * y[l]
                end
                y[i] = y[i] / r[i, i]
                i += 1
            end
            I[j] = i_star
        end
        I[n+1] = i - 1
    end
    return sum(c .* y)
end


# x is a single point in n space
function vrval(x::AbstractVector{T}, vf::VanderFit{T}) where {T}
    n = veval(x, vf.a, vf.RL, vf.M)
    d = veval(x, vf.b, vf.RL, vf.N)
    return n / d
end
