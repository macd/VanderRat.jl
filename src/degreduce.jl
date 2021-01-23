

# recall that L = max(M, N)
function degreduce(VL, f, M, N, eta=1e-12)
    
    # Reduce denominator degree
    FVN = md(f, VL[:, 1:N])
    MT = M
    while true
        @views u, s, vt = svd(VL[:, 1:MT-1]' * FVN)
        smin, smax = extrema(s)
        println("d: ", smin, "  ", smax, "   ", MT)
        if smin < eta * smax
            MT -= 1
        else
            break
        end
    end
    
    # Reduce numerator degree
    FIVM = md(f .^ -1, VL[:, 1:M])
    NT = N
    while true
        @views u, s, vt = svd(VL[:, 1:NT-1]' * FIVM)
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
