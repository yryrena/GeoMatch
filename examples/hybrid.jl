## hybrid matching example (geo + propensity score refinement) 
import Pkg; Pkg.activate(".")

using DataFrames
using GeoInterface, GeometryBasics
using GeoMatch
using StatsModels
using Random

## set seed for reproducibility
Random.seed!(123)

## treated table
t = DataFrame(
    id = 1:10,
    x  = 116.3 .+ 0.03 .* rand(10),
    y  = 39.85 .+ 0.03 .* rand(10),
    income = 50 .+ 5 .* rand(10),
    pop    = 20 .+ 2 .* rand(10),
    urban  = rand(Bool, 10),
)
t.geometry = [GeometryBasics.Point(t.x[i], t.y[i]) for i in 1:nrow(t)]

## control table
c = DataFrame(
    id = 101:300,
    x  = 116.1 .+ 0.6 .* rand(200),
    y  = 39.7  .+ 0.6 .* rand(200),
    income = 48 .+ 6 .* rand(200),
    pop    = 22 .+ 3 .* rand(200),
    urban  = rand(Bool, 200),
)
c.geometry = [GeometryBasics.Point(c.x[i], c.y[i]) for i in 1:nrow(c)]

## run geo k-nn matching with hybrid propensity score refinement
res = GeoMatch.match_spatial(
    t, c;
    covars=[:income, :pop, :urban],
    method="geoNN",
    k=2,
    radius_km=80.0,
    distance_metric=:haversine,
    hybrid=true,
    ps_formula=@formula(treat ~ income + pop + urban)
)

## balance diagnostics
bt = GeoMatch.balance_table(res, t, c, [:income, :pop])

## write a minimal html report
GeoMatch.html_report("geomatch_hybrid.html", res, bt)
println("report written to geomatch_hybrid.html")
 
println(first(res.pairs, 6))
