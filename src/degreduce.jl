# recall that L = max(M, N)
#function calcreduce(VL::AbstractMatrix{T}, f::AbstractVector{T}, M, N, n; eta=1e-12) where {T}
function calcreduce(vf::VanderFit; eta=1e-12) where {T}
    UM = binomial(vf.Mmax + vf.n, vf.Mmax)
    UN = binomial(vf.Nmax + vf.n, vf.Nmax)
    # Reduce denominator degree
    FVN = Diagonal(f) * vf.VL[:, 1:UN]
    MT = vf.Mmax
    while true
        UMT1 = binomial(MT-1+n, MT-1)
        @views u, s, vt = svd(VL[:, 1:UMT1]' * FVN)
        smin, smax = extrema(s)
        println("d: ", smin, "  ", smax, "   ", MT)
        if smin < eta * smax
            MT -= 1
        else
            break
        end
    end
    
    # Reduce numerator degree using the new VM
    UMT = binomial(MT+n, MT)
    FIVM = Diagonal(f .^ -1) * VL[:, 1:UMT]
    NT = vf.Nmax
    while true
        UNT1 = binomial(NT-1+n, NT-1)
        @views u, s, vt = svd(VL[:, UNT1]' * FIVM)
        smin, smax = extrema(s)
        println("n: ", smin, "  ", smax, "   ", NT)
        if smin < eta * smax
            NT -= 1
        else
            break
        end
    end

    return MT, NT
end

function degreduce!(vf, M, N)
    if M > vf.M || N > vf.N
        println("Error: Trying to resize to larger degree")
        return nothing
    end
    UM = binomial(vf.n + M, M)
    UN = binomial(vf.n + N, N)

    VM, RM = (view(vf.VL, :, 1:UM), view(vf.RL, :, 1:UM))
    VN, RN = (view(vf.VL, :, 1:UN), view(vf.RL, :, 1:UN))
    
    fVN = Diagonal(vf.f) * VN
    Z = pinv(VM) * fVN
    W = VM * Z - fVN
    G = svd(W)
    vf.b = G.Vt[end, :]
    vf.a = Z * vf.b
    vf.M = M
    vf.N = N
    nothing
end
