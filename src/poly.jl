# These are the single dimensional polynomial fitting routines
# from "Vandermonde with Arnoldi"

function polyfitA(x::S, f::S, n::Int) where {T, S <: AbstractVector{T}}
    m = length(x)
    H = zeros(T, n+1, n)
    Q = Array{T, 2}(undef, m, n + 1)
    Q[:, 1] .= one(T)
    q = zeros(T, m)
    @inbounds @views for k = 1:n
        q .= x .* Q[:, k]
        for j = 1:k
            H[j, k] = Q[:, j]' * q / m
            q .-= H[j, k] * Q[:, j]
        end
        H[k+1, k]  = norm(q) / sqrt(m)
        Q[:, k+1] .= q / H[k+1, k]
    end
    d = Q \ f
    return d, H
end

function polyvalA(d::S, H::U, s::S) where {T, S <: AbstractVector{T},
                                              U <: AbstractMatrix{T}}
    m = length(s)
    n = size(H, 2)
    W = Array{T, 2}(undef, m, n + 1)
    W[:, 1] .= one(T)
    w = zeros(T, m)
    
    @inbounds @views for k = 1:n
        w .= s .* W[:, k]
        for j = 1:k
            w .-= H[j, k] * W[:, j]
        end
        W[:, k+1] = w / H[k+1,k]
    end
    return W * d
end

# Algorithm 2 from "Convergence and Near-optimal Sampling..." by
# Zhu and Nakatsukasa (2023). This is a refactor of the algorithm from
# "Vandermonde with Arnoldi"
# Biggest change is from mgs to 2 times cgs, but to accomplish this
# we need to rewrite the inner loop as a matrix operation. Then, around
# that we add an iteration. For this to work, we need a temp variable,
# here called s.
# Also, changing all the variable names is very annoying.
# NB: The N in the algorithm must be n+1. That should be stated.
function polyfitD(x::S, f::S, n::Int; iter=2) where {T, S <: AbstractVector{T}}
    m = length(x)
    H = zeros(T, n+1, n)
    Q = Array{T, 2}(undef, m, n + 1)
    Q[:, 1] .= one(T)
    v = zeros(T, m)
    s = zeros(T, m)
    @inbounds @views for l = 1:n
        v .= x .* Q[:, l]
        for t in 1:iter
            s[1:l] .= Q[:, 1:l]' * v / m
            v .-= Q[:, 1:l] * s[1:l]
            H[1:l, l] .+= s[1:l]
        end
        H[l+1, l]  = norm(v) / sqrt(m)
        Q[:, l+1] .= v / H[l+1, l]
    end
    d = Q \ f
    return d, H
end


function polyvalD(d::S, H::V, y::S) where {T, S <: AbstractVector{T},
                                              V <: AbstractMatrix{T}}
    m = length(y)
    n = size(H, 2)
    U = Array{T, 2}(undef, m, n + 1)
    U[:, 1] .= one(T)
    v = zeros(T, m)
    @inbounds @views for l = 1:n
        v .= y .* U[:, l]
        v .-= U[:, 1:l] * H[1:l, l]
        U[:, l+1] .= v / H[l+1, l]
    end
    return U * d
end
