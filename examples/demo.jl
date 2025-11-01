## minimal demo
## comments use lowercase and start with ##
import Pkg; Pkg.activate(".")
using DataFrames, GeoInterface, GeometryBasics
using GeoMatch

## build treated table with wgs84 lon-lat
t = DataFrame(id=1:3,
              x=[116.4,116.5,116.6],
              y=[39.9,39.92,39.95],
              income=50 .+ 10rand(3),
              pop=20 .+ 5rand(3))
t.geometry = [GeometryBasics.Point(t.x[i], t.y[i]) for i in 1:nrow(t)]

## build control table
c = DataFrame(id=101:130,
              x=116.1 .+ rand(30),
              y=39.7 .+ 0.6rand(30),
              income=48 .+ 10rand(30),
              pop=22 .+ 5rand(30))
c.geometry = [GeometryBasics.Point(c.x[i], c.y[i]) for i in 1:nrow(c)]

## run matching
res = GeoMatch.match_spatial(t, c;
    covars=[:income, :pop],
    method="geoNN",
    k=2,
    radius_km=50.0,
    distance_metric=:haversine,
    replace=false)

## show first pairs
println(first(res.pairs, 6))

## balance and report
bt = GeoMatch.balance_table(res, t, c, [:income, :pop])
GeoMatch.html_report("geomatch_report.html", res, bt)
println("report written to geomatch_report.html")
