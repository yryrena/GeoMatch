## kernel matching example
import Pkg; Pkg.activate(".")

using DataFrames
using GeoInterface, GeometryBasics
using GeoMatch
using Random

## set seed for reproducibility
Random.seed!(123)

## build treated table (small cluster near some center)
t = DataFrame(
    id = 1:8,
    x  = 116.3 .+ 0.03 .* rand(8),
    y  = 39.85 .+ 0.03 .* rand(8),
    income = 50 .+ 5 .* rand(8),
    pop    = 20 .+ 2 .* rand(8),
)
t.geometry = [GeometryBasics.Point(t.x[i], t.y[i]) for i in 1:nrow(t)]

## build control table (wider area)
c = DataFrame(
    id = 101:260,
    x  = 116.1 .+ 0.6 .* rand(160),
    y  = 39.7  .+ 0.6 .* rand(160),
    income = 48 .+ 6 .* rand(160),
    pop    = 22 .+ 3 .* rand(160),
)
c.geometry = [GeometryBasics.Point(c.x[i], c.y[i]) for i in 1:nrow(c)]

## run kernel matching
res = GeoMatch.match_spatial(
    t, c;
    covars=[:income, :pop],
    method="kernel",
    k=5,                      ## take 5 nearest controls as kernel candidates
    radius_km=80.0,
    bandwidth=12.0,           ## kernel bandwidth in km
    distance_metric=:haversine
)

## balance diagnostics
bt = GeoMatch.balance_table(res, t, c, [:income, :pop])

## write a minimal html report
GeoMatch.html_report("geomatch_kernel.html", res, bt)
println("report written to geomatch_kernel.html")
 
println(first(res.pairs, 6))