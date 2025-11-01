module GeoMatch

using GeoInterface 
using DataFrames, StatsBase, StatsModels, GLM
using Distances
using NearestNeighbors
using LinearAlgebra
using Random
using Statistics

export MatchConfig, MatchResult, match_spatial, balance_table, love_plot!, plot_matches!, html_report

############################
## core data structures
############################

struct MatchConfig
    method::Symbol                # :geoNN | :kernel | :caliper
    k::Int
    radius_km::Union{Nothing,Float64}
    kernel::Symbol                # :gaussian | :triangular | :epanechnikov
    bandwidth::Union{Nothing,Float64}
    hybrid::Bool
    ps_formula::Union{Nothing,FormulaTerm}
    region_col::Union{Nothing,Symbol}
    distance_metric::Symbol       # :haversine | :euclidean
    crs_target::Union{Nothing,String}
    replace::Bool
    seed::Union{Nothing,Int}
end

struct MatchResult
    pairs::DataFrame              # columns: :t_id, :c_id, :distance, :weight
    dropped_treated::Vector{Int}
    dropped_control::Vector{Int}
    config::MatchConfig
    summary::Dict{String,Any}
end

############################
## utilities: geometry & distance
############################

## note: users should ensure geometries are points or have centroids computed upstream

function _coords_lonlat(gdf::AbstractDataFrame)
    ## extract lon/lat arrays from geodataframe geometries
    ## assume geometries are points in lon-lat (wgs84); upstream reprojection handled separately
    lon = Float64[]
    lat = Float64[]
    for geom in gdf.geometry
        p = GeoInterface.coordinates(geom)
        push!(lon, p[1]); push!(lat, p[2])
    end
    return hcat(lon, lat)
end

function _haversine_matrix(x::Matrix{Float64}, y::Matrix{Float64})
    ## compute pairwise haversine distance (km) between rows of x and y
    n, m = size(x,1), size(y,1)
    D = Matrix{Float64}(undef, n, m)
    for i in 1:n
        for j in 1:m
            D[i,j] = haversine((x[i,2], x[i,1]), (y[j,2], y[j,1])) / 1000.0
        end
    end
    return D
end

function _distance_matrix(treated::AbstractDataFrame, control::AbstractDataFrame; metric::Symbol=:haversine)
    ## dispatch to distance backend; mvp supports haversine only
    if metric == :haversine
        xt = _coords_lonlat(treated)
        xc = _coords_lonlat(control)
        return _haversine_matrix(xt, xc)
    else
        error("distance metric not implemented in mvp")
    end
end

############################
## candidate filtering
############################

function _apply_constraints!(D::Matrix{Float64}; radius_km::Union{Nothing,Float64}=nothing,
                             region_t::Union{Nothing,Vector}=nothing,
                             region_c::Union{Nothing,Vector}=nothing)
    ## set distances to +inf if violating radius or region constraints
    n, m = size(D)
    if radius_km !== nothing
        for i in 1:n, j in 1:m
            if D[i,j] > radius_km
                D[i,j] = Inf
            end
        end
    end
    if region_t !== nothing && region_c !== nothing
        for i in 1:n, j in 1:m
            if region_t[i] != region_c[j]
                D[i,j] = Inf
            end
        end
    end
    return D
end

############################
## matching strategies
############################

function _knn_match(D::Matrix{Float64}; k::Int=1, replace::Bool=false)
    ## perform constrained k-nn matching on distance matrix with +inf as invalid
    n, m = size(D)
    t_ids = Int[]
    c_ids = Int[]
    dists = Float64[]
    used = falses(m)
    for i in 1:n
        row = D[i, :]
        order = sortperm(eachindex(row); by=j->row[j])
        picked = 0
        for j in order
            if !isfinite(row[j]); continue; end
            if !replace && used[j]; continue; end
            push!(t_ids, i); push!(c_ids, j); push!(dists, row[j])
            used[j] = true
            picked += 1
            if picked == k; break; end
        end
    end
    df = DataFrame(t_id=t_ids, c_id=c_ids, distance=dists)
    df.weight = ones(size(df,1))
    return df
end

function _kernel_weights(dist::Vector{Float64}; kernel::Symbol=:gaussian, bandwidth::Float64)
    ## compute nonnegative kernel weights given distances and bandwidth
    w = similar(dist)
    if kernel == :gaussian
        @inbounds for i in eachindex(dist)
            z = dist[i] / bandwidth
            w[i] = exp(-0.5*z*z)
        end
    elseif kernel == :triangular
        @inbounds for i in eachindex(dist)
            z = dist[i] / bandwidth
            w[i] = max(0.0, 1.0 - abs(z))
        end
    else
        error("kernel not implemented")
    end
    return w ./ sum(w)
end

function _kernel_match(D::Matrix{Float64}; bandwidth::Float64, k::Int=5)
    ## for each treated, take k nearest valid controls then assign kernel weights
    n, m = size(D)
    t_ids = Int[]; c_ids = Int[]; dists = Float64[]; weights = Float64[]
    for i in 1:n
        row = D[i, :]
        idx = sortperm(eachindex(row); by=j->row[j])
        cand = [j for j in idx if isfinite(row[j])][1:min(k, count(isfinite, row))]
        if isempty(cand); continue; end
        dist_i = row[cand]
        w = _kernel_weights(dist_i; bandwidth=bandwidth)
        for (cj, wj, dj) in zip(cand, w, dist_i)
            push!(t_ids, i); push!(c_ids, cj); push!(dists, dj); push!(weights, wj)
        end
    end
    return DataFrame(t_id=t_ids, c_id=c_ids, distance=dists, weight=weights)
end

############################
## hybrid with propensity score
############################

function _estimate_ps(df::DataFrame, formula::FormulaTerm)
    ## estimate propensity scores via glm (logit)
    m = glm(formula, df, Binomial(), LogitLink())
    return predict(m)
end

function _hybrid_refine(pairs::DataFrame, ps_t::Vector{Float64}, ps_c::Vector{Float64})
    ## reorder or filter pairs by |ps_t - ps_c| within treated groups
    pairs.ps_gap = similar(pairs.distance)
    @inbounds for r in 1:size(pairs,1)
        pairs.ps_gap[r] = abs(ps_t[pairs.t_id[r]] - ps_c[pairs.c_id[r]])
    end
    sort!(pairs, [:t_id, :ps_gap])
    return pairs
end

############################
## public api: match_spatial
############################

function match_spatial(treated::AbstractDataFrame, control::AbstractDataFrame;
                       covars::Vector{Symbol},
                       method::AbstractString="geoNN",
                       k::Int=1,
                       radius_km::Union{Nothing,Float64}=nothing,
                       kernel::AbstractString="gaussian",
                       bandwidth::Union{Nothing,Float64}=nothing,
                       caliper::Union{Nothing,Float64}=nothing,   # reserved for ps caliper
                       hybrid::Bool=false,
                       ps_formula=nothing,
                       region_col::Union{Nothing,Symbol}=nothing,
                       distance_metric::Symbol=:haversine,
                       crs_target::Union{Nothing,String}=nothing,
                       replace::Bool=false,
                       seed::Union{Nothing,Int}=nothing)

    ## mvp: haversine only; kernel bandwidth required for kernel method
    rng = seed === nothing ? Random.GLOBAL_RNG : MersenneTwister(seed)

    D = _distance_matrix(treated, control; metric=distance_metric)

    region_t = region_col === nothing ? nothing : Vector{Any}(treated[:, region_col])
    region_c = region_col === nothing ? nothing : Vector{Any}(control[:, region_col])
    _apply_constraints!(D; radius_km=radius_km, region_t=region_t, region_c=region_c)

    pairs = if method == "geoNN"
        _knn_match(D; k=k, replace=replace)
    elseif method == "kernel"
        bandwidth === nothing && error("bandwidth is required for kernel matching")
        _kernel_match(D; bandwidth=bandwidth, k=k)
    else
        error("method not implemented in mvp")
    end

    if hybrid
        ps_formula === nothing && error("ps_formula is required for hybrid")
        ## build combined df to fit ps for control; assume treated has indicator = 1, control = 0
        df_t = DataFrame(treated[:, covars]); df_t.:treat = 1
        df_c = DataFrame(control[:, covars]); df_c.:treat = 0
        df = vcat(df_t, df_c; cols=:union)
        ps = _estimate_ps(df, ps_formula)
        ps_t = ps[1:size(df_t,1)]
        ps_c = ps[(size(df_t,1)+1):end]
        pairs = _hybrid_refine(pairs, ps_t, ps_c)
    end

    dropped_t = setdiff(collect(1:size(treated,1)), unique(pairs.t_id))
    dropped_c = setdiff(collect(1:size(control,1)), unique(pairs.c_id))

    summary = Dict(
        "method" => method,
        "n_treated" => size(treated,1),
        "n_control" => size(control,1),
        "n_pairs" => size(pairs,1),
        "radius_km" => radius_km,
        "hybrid" => hybrid
    )

    return MatchResult(pairs, dropped_t, dropped_c,
                       MatchConfig(Symbol(method), k, radius_km,
                                   Symbol(kernel), bandwidth, hybrid, ps_formula,
                                   region_col, distance_metric, crs_target, replace, seed),
                       summary)
end

############################
## diagnostics: balance
############################

function _smd(x_t::AbstractVector, x_c::AbstractVector, w::AbstractVector)
    ## compute standardized mean difference with optional weights for controls
    mt = mean(x_t)
    mc = sum(w .* x_c) / sum(w)
    st = std(x_t); sc = sqrt(sum(w .* (x_c .- mc).^2) / sum(w))
    return (mt - mc) / sqrt((st^2 + sc^2) / 2)
end

function balance_table(res::MatchResult, treated::AbstractDataFrame, control::AbstractDataFrame, covars::Vector{Symbol})
    ## return dataframe with smd before/after matching
    ## build weights: for each pair, control-side weights sum within treated
    w = zeros(size(control,1))
    for r in 1:size(res.pairs,1)
        w[res.pairs.c_id[r]] += res.pairs.weight[r]
    end
    w .+= 1e-12

    rows = String[]; smd_before = Float64[]; smd_after = Float64[]
    for v in covars
        x_t = Vector{Float64}(treated[:, v])
        x_c = Vector{Float64}(control[:, v])
        push!(rows, String(v))
        push!(smd_before, (mean(x_t) - mean(x_c)) / sqrt((var(x_t) + var(x_c)) / 2))
        push!(smd_after, _smd(x_t, x_c, w))
    end
    return DataFrame(var=rows, smd_before=smd_before, smd_after=smd_after)
end

############################
## visualization stubs
############################

function love_plot!(ax, tbl::DataFrame; threshold=0.1)
    ## plot smd_before and smd_after on a shared axis (user passes axis from makie)
    ## mvp: left to user; here we return the processed arrays for plotting
    return (tbl.var, tbl.smd_before, tbl.smd_after, threshold)
end

function plot_matches!(ax, treated::AbstractDataFrame, control::AbstractDataFrame, res::MatchResult)
    ## mvp: return coordinates for lines; user draws with geomakie
    ## returns vectors: (lon_t, lat_t, lon_c, lat_c, weight)
    xt = _coords_lonlat(treated); xc = _coords_lonlat(control)
    lon_t = Float64[]; lat_t = Float64[]; lon_c = Float64[]; lat_c = Float64[]; w = Float64[]
    for r in 1:size(res.pairs,1)
        i = res.pairs.t_id[r]; j = res.pairs.c_id[r]
        push!(lon_t, xt[i,1]); push!(lat_t, xt[i,2])
        push!(lon_c, xc[j,1]); push!(lat_c, xc[j,2])
        push!(w, res.pairs.weight[r])
    end
    return (lon_t, lat_t, lon_c, lat_c, w)
end

############################
## html report (minimal)
############################

function html_report(path::AbstractString, res::MatchResult, balance::DataFrame)
    ## write a very simple html file summarizing matching
    open(path, "w") do io
        write(io, """
        <html><head><meta charset="utf-8"><title>GeoMatch Report</title></head>
        <body>
          <h1>GeoMatch Report</h1>
          <h2>Summary</h2>
          <pre>$(res.summary)</pre>
          <h2>Balance (SMD)</h2>
          <table border="1" cellpadding="6" cellspacing="0">
            <tr><th>var</th><th>smd_before</th><th>smd_after</th></tr>
            $(join([ "<tr><td>$(balance.var[i])</td><td>$(round(balance.smd_before[i], digits=3))</td><td>$(round(balance.smd_after[i], digits=3))</td></tr>" for i in 1:size(balance,1) ], ""))
          </table>
        </body></html>
        """)
    end
    return path
end

end # module