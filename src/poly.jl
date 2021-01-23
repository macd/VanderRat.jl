# These are the single dimensional polynomial fitting routines
# from "Vandermonde with Arnoldi" 
function polyfitA(x::S, f::S, n) where {T, S <: AbstractVector{T}}
    m = length(x)
    H = zeros(T, n+1, n)
    Q = Array{T, 2}(undef, m, n + 1)
    Q[:, 1] .= one(T)
    @inbounds @views for k = 1:n
        q = x .* Q[:, k]
        for j = 1:k
            H[j, k] = Q[:, j]' * q/m
            q = q - H[j, k] * Q[:, j]
        end
        H[k+1, k] = norm(q) / sqrt(m)
        Q[:, k+1] = q / H[k+1, k]
    end
    d = Q \ f
    return (d, H, Q)
end

function polyvalA(d::S, H::U, s::S) where {T, S <: AbstractVector{T}, U <: AbstractMatrix{T}}
    m = length(s)
    n = size(H, 2)
    W = Array{T, 2}(undef, m, n + 1)
    W[:, 1] .= one(T)
    
    @inbounds @views for k = 1:n
        w = s .* W[:, k]
        for j = 1:k
            w = w - H[j, k] * W[:, j]
        end
        W[:, k+1] = w / H[k+1,k]
    end
    return W * d
end
