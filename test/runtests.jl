## minimal tests for mvp
## comments use lowercase and start with ##
using Test
using GeoMatch
using DataFrames, GeoInterface, GeometryBasics

@testset "geomatch mvp" begin
    t = DataFrame(id = 1:2, x = [0.0, 0.1], y = [0.0, 0.1], a = [1.0, 2.0])
    t.geometry = [GeometryBasics.Point(t.x[i], t.y[i]) for i = 1:nrow(t)]

    c = DataFrame(id = 101:104, x = rand(4), y = rand(4), a = rand(4))
    c.geometry = [GeometryBasics.Point(c.x[i], c.y[i]) for i = 1:nrow(c)]

    res = GeoMatch.match_spatial(
        t,
        c;
        covars = [:a],
        method = "geoNN",
        k = 1,
        radius_km = 5000.0,
    )
    @test nrow(res.pairs) == 2

    bt = GeoMatch.balance_table(res, t, c, [:a])
    @test all(ismissing.(bt.smd_before) .== false)
end
