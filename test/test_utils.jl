function run_case(cas; K=1000, M=20, N=20, EPS=1e-6, ntest=500)
    H, info = austin_examples(cas)
    n = info[:n]
    bounds = info[:bound]

    # training samples
    X = hcat([rand(K) .* (b[2] - b[1]) .+ b[1] for b in bounds]...)
    f = H(X)

    vf = vrfit(X, f, M, N)

    # test samples
    Xt = hcat([rand(ntest) .* (b[2] - b[1]) .+ b[1] for b in bounds]...)
    ft = H(Xt)
    rt = vrval.(eachrow(Xt), Ref(vf))

    return maximum(abs.(ft - rt))
end

function run_case_allow_rank_deficiency(cas; kwargs...)
    try
        err = run_case(cas; kwargs...)
        return (:ok, err)
    catch e
        if occursin("near rank deficient", sprint(showerror, e))
            return (:rank_deficient, NaN)
        else
            rethrow(e)
        end
    end
end
