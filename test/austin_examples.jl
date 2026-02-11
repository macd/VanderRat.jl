"""
    examples(CAS::Integer)

Julia translation of MATLAB mLF `examples.m`.

Returns:
- `H`: callable;
      supports `H(X::AbstractMatrix)`
      (N×n, rows are points) -> Vector{Float64}
      and `H(x::AbstractVector)` (length n) -> Float64

- `info`: Dict{Symbol,Any} with metadata and grids (p_c, p_r, ip, bounds, etc.)
"""
function austin_examples(CAS::Integer)
    Random.seed!(137)
    rnd  = 0.0
    Nip  = 20
    xbnd = [-1.0, 1.0]
    mpts = 1
    tag  = String[]

    # --- per-case definitions ---

    # bounds can be either a single 2-vector or a vector-of-2-vectors
    xbnd_cell = nothing  # will set to Vector{Vector{Float64}} when needed

    # Helper to make H accept both matrix and vector
    Hmat = nothing  # will be set to a function X::AbstractMatrix -> Vector{Float64}

    if CAS == 1
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.(exp.(x1 .* x2) ./ ((x1.^2 .- 1.44) .* (x2.^2 .- 1.44)))
        end

    elseif CAS == 2
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.(log2.(2.25 .- x1.^2 .- x2.^2))
        end

    elseif CAS == 3
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.(tanh.(5 .* (x1 .- x2)))
        end

    elseif CAS == 4
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.(exp.(-(x1.^2 .+ x2.^2) ./ 1000))
        end

    elseif CAS == 5
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.(abs.(x1 .- x2).^3)
        end

    elseif CAS == 6
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.((x1 .+ x2.^3) ./ (x1 .* x2.^2 .+ 2))
        end
        xbnd = [1e-10, 1.0]

    elseif CAS == 7
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            num = x1.^2 .+ x2.^2 .+ x1 .- x2 .- 1
            den = (x1 .- 1.1) .* (x2 .- 1.1)
            Float64.(num ./ den)
        end

    elseif CAS == 8
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            num = x1.^4 .+ x2.^4 .+ (x1.^2 .* x2.^2) .+ (x1 .* x2)
            den = (x1 .- 1.1) .* (x2 .- 1.1)
            Float64.(num ./ den)
        end

    elseif CAS == 9
        n = 4
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]; x3 = @view X[:,3]; x4 = @view X[:,4]
            num = x1.^2 .+ x2.^2 .+ x1 .- x2 .+ 1
            den = (x3 .- 1.5) .* (x4 .- 1.5)
            Float64.(num ./ den)
        end

    elseif CAS == 10
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            num = x1.^2 .+ x2.^2 .+ x1 .- x2 .- 1
            den = x1.^3 .+ x2.^3 .+ 4
            Float64.(num ./ den)
        end

    elseif CAS == 11
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            Float64.((x1.^3 .+ x2.^3) ./ (x1.^2 .+ x2.^2 .+ 3))
        end

    elseif CAS == 12
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            num = x1.^4 .+ x2.^4 .+ x1.^2 .* x2.^2 .+ x1 .* x2
            den = x1.^2 .* x2.^2 .- 2 .* x1.^2 .- 2 .* x2.^2 .+ 4
            Float64.(num ./ den)
        end

    elseif CAS == 13
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            num = x1.^3 .+ x2.^3
            den = x1.^2 .* x2.^2 .- 2 .* x1.^2 .- 2 .* x2.^2 .+ 4
            Float64.(num ./ den)
        end

    elseif CAS == 14
        n = 2
        Hmat = X -> begin
            x1 = @view X[:,1]; x2 = @view X[:,2]
            num = x1.^4 .+ x2.^4 .+ x1.^2 .* x2.^2 .+ x1 .* x2
            den = x1.^3 .+ x2.^3 .+ 4
            Float64.(num ./ den)
        end

    elseif CAS == 15
        # This needs work on correctly defining the bounds
        n = 3
        Hmat = X -> begin
            E = @view X[:,1]; G = @view X[:,2]; M = @view X[:,3]
            num = 2 .* sqrt(2.0) .* M .* G .* sqrt.(M.^2 .* (M.^2 .+ G.^2))
            den = (pi .* sqrt.(M.^2 .+ sqrt.(M.^2 .* (M.^2 .+ G.^2)))) .* ((E.^2 .- M.^2).^2 .+ M.^2 .* G.^2)
            Float64.(num ./ den)
        end


    elseif CAS == 16
        n = 4
        Hmat = X -> begin
            x1 = @view X[:,1]; x2=@view X[:,2]; x3=@view X[:,3]; x4=@view X[:,4]
            num = atan.(x1) .+ atan.(x2) .+ atan.(x3) .+ atan.(x4)
            den = x1.^2 .* x2.^2 .- x1.^2 .- x2.^2 .+ 1
            Float64.(num ./ den)
        end
        xbnd = [-0.95, 0.95]

    elseif CAS == 17
        n = 4
        Hmat = X -> begin
            x1=@view X[:,1]; x2=@view X[:,2]; x3=@view X[:,3]; x4=@view X[:,4]
            num = exp.(x1 .* x2 .* x3 .* x4)
            den = x1.^2 .+ x2.^2 .- x3 .* x4 .+ 3
            Float64.(num ./ den)
        end

    elseif CAS == 18
        n = 4
        Hmat = X -> begin
            x1=@view X[:,1]; x2=@view X[:,2]; x3=@view X[:,3]; x4=@view X[:,4]
            Float64.(10 .* (sinc.(x1).*sinc.(x2).*sinc.(x3).*sinc.(x4)))
        end
        xbnd = [1e-6, 4*pi]

    elseif CAS == 19
        n = 2
        Hmat = X -> begin
            x1=@view X[:,1]; x2=@view X[:,2]
            Float64.(10 .* (sin.(x1).*sin.(x2)) ./ (x1.*x2))
        end
        xbnd = [1e-6, 4*pi]

    elseif CAS == 20
        n = 2
        Hmat = X -> begin
            x1=@view X[:,1]; x2=@view X[:,2]
            Float64.(x1.^2 .+ x2.^2 .+ x1.*x2 .- x2 .+ 1)
        end

    else
        error("examples: case $CAS not implemented in this translation block")
    end

    xlim = Vector{Vector{Float64}}(undef, n)
    if xbnd_cell === nothing
        for ii in 1:n
            xlim[ii] = [float(xbnd[1]), float(xbnd[2])]
        end
    else
        @assert length(xbnd_cell) == n "xbnd cell length must match n"
        for ii in 1:n
            xlim[ii] = [float(xbnd_cell[ii][1]), float(xbnd_cell[ii][2])]
        end
    end


    # --- wrap H to accept vector points too ---
    H = function (x)
        if x isa AbstractVector
            @assert length(x) == n "H(x): length(x) must be $n"
            X = reshape(Float64.(collect(x)), 1, n)
            return Hmat(X)[1]
        elseif x isa AbstractMatrix
            @assert size(x,2) == n "H(X): size(X,2) must be $n"
            return Hmat(x)
        else
            error("H expects a vector (length n) or a matrix (N×n).")
        end
    end

    info = Dict{Symbol,Any}()
    info[:n]      = n
    info[:mpts]   = mpts
    info[:bound]  = xlim

    return H, info
end
