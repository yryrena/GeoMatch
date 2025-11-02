# GeoMatch.jl

[![ci](https://github.com/yryrena/GeoMatch/actions/workflows/ci.yml/badge.svg)](https://github.com/yryrena/GeoMatch/actions/workflows/ci.yml)

[![ascii guard](https://github.com/yryrena/GeoMatch/actions/workflows/ascii-guard.yml/badge.svg)](https://github.com/yryrena/GeoMatch/actions/workflows/ascii-guard.yml)
[![julia tests](https://github.com/yryrena/GeoMatch/actions/workflows/julia-tests.yml/badge.svg)](https://github.com/yryrena/GeoMatch/actions/workflows/julia-tests.yml)

Minimal, working proof-of-concept for spatial matching in Julia.

Given a **treated** table and a **control** table with point geometries (WGS84 lon-lat), GeoMatch builds matched pairs or weighted sets based on great-circle distance, with optional covariate balancing diagnostics.

MVP scope: haversine distance, K-nearest neighbors or kernel matching, basic balance table, and a tiny HTML report.

------

## Features

- Haversine distance on lon-lat coordinates  
- K-nearest neighbor matching with or without replacement  
- Kernel matching with Gaussian or triangular kernels  
- Optional radius caliper and same-region constraint  
- Simple balance diagnostics (standardized mean differences)  
- Minimal HTML report generator  

------

## Status

- Tested on Julia `1.11`
- Lightweight dependencies: `DataFrames`, `Distances`, `NearestNeighbors`, `StatsBase`, `StatsModels`, `GLM`, `GeoInterface`
- Geometries are passed via a `geometry` column containing `GeometryBasics.Point` values

------

## Installation

From a local checkout:

```julia
julia> ]
(@v1.11) pkg> dev /path/to/GeoMatch
(@v1.11) pkg> test GeoMatch
```

From GitHub: 

```julia
julia> ]
(@v1.11) pkg> add https://github.com/yryrena/GeoMatch.git
(@v1.11) pkg> test GeoMatch
```

> Recommendation: do not commit `Manifest.toml` for a package repo. 

---

## Quick Start

```julia
using DataFrames
using GeometryBasics          ## for Point(x, y)
using GeoInterface            ## used internally; not strictly required in user code
using GeoMatch

## treated table
t = DataFrame(
    id = 1:3,
    x  = [116.4, 116.5, 116.6],      ## longitude
    y  = [39.90, 39.92, 39.95],      ## latitude
    income = 50 .+ 10*rand(3),
    pop    = 20 .+  5*rand(3)
)
t.geometry = [GeometryBasics.Point(t.x[i], t.y[i]) for i in 1:nrow(t)]

## control table
c = DataFrame(
    id = 101:130,
    x  = 116.1 .+ rand(30),
    y  = 39.7  .+ 0.6*rand(30),
    income = 48 .+ 10*rand(30),
    pop    = 22 .+  5*rand(30)
)
c.geometry = [GeometryBasics.Point(c.x[i], c.y[i]) for i in 1:nrow(c)]

## run matching
res = GeoMatch.match_spatial(
    t, c;
    covars = [:income, :pop],
    method = "geoNN",              ## "geoNN" or "kernel"
    k = 2,
    radius_km = 50.0,
    distance_metric = :haversine,
    replace = false
)

## show first pairs
first(res.pairs, 6) |> println

## balance diagnostics
bt = GeoMatch.balance_table(res, t, c, [:income, :pop])

## minimal HTML report
GeoMatch.html_report("geomatch_report.html", res, bt)
println("report written to geomatch_report.html")
```

Example output for `res.pairs`:

```text
6x4 DataFrame
 Row | t_id  c_id  distance  weight
     | Int64 Int64 Float64   Float64
-----+-------------------------------
   1 |     1    16   6.88871     1.0
   2 |     1    19  11.76730     1.0
   3 |     2    10   8.28476     1.0
   4 |     2    22   8.51346     1.0
   5 |     3    26   5.95047     1.0
   6 |     3     7   8.30584     1.0
```
```bash
julia --project=. examples/kernel.jl
julia --project=. examples/hybrid.jl
```

- `kernel.jl` shows spatial kernel matching with a Gaussian kernel (`method="kernel"`, `k=5`, `bandwidth=12.0` km).
- `hybrid.jl` shows geo K-NN with propensity-score refinement (`hybrid=true`, `ps_formula=@formula(treat ~ income + pop + urban)`).

Reports are saved as `geomatch_kernel.html` and `geomatch_hybrid.html`.

------

## API

#### `match_spatial(treated, control; kwargs...) -> MatchResult`

**Inputs**

- `treated::AbstractDataFrame`, `control::AbstractDataFrame`
  Must each contain a `geometry` column with `GeometryBasics.Point` values whose `.x` and `.y` represent longitude and latitude (WGS84).
- `covars::Vector{Symbol}`
   Covariates used for balance reporting and, when `hybrid = true`, in the propensity score model.
- `method::String = "geoNN"`
  - `"geoNN"`: distance-based KNN matching
  - `"kernel"`: distance-weighted matching using kernels
- `k::Int = 1`
   Number of matches per treated unit (for `"kernel"`, k nearest controls are used to compute weights).
- `radius_km::Union{Nothing, Float64} = nothing`
   Maximum allowed distance. Pairs beyond this are discarded.
- `kernel::String = "gaussian"`
   `"gaussian"` or `"triangular"` (only used with `method = "kernel"`).
- `bandwidth::Union{Nothing, Float64} = nothing`
   Required if `method = "kernel"`.
- `hybrid::Bool = false`
   If true, estimate a logit propensity score and sort within treated units by absolute PS gap.
- `ps_formula = nothing`
   A `StatsModels` formula, e.g. `@formula(treat ~ income + pop)`. Required if `hybrid = true`.
- `region_col::Union{Nothing, Symbol} = nothing`
   If provided, only match within equal region labels across treated and control.
- `distance_metric::Symbol = :haversine`
   Currently only `:haversine` is implemented.
- `crs_target::Union{Nothing, String} = nothing`
   Reserved for future CRS handling. Not used in the MVP.
- `replace::Bool = false`
   If false, each control unit can be used at most once in KNN.
- `seed::Union{Nothing, Int} = nothing`
   Optional RNG seed for reproducibility where randomization is introduced (not used by deterministic KNN).

**Output: `MatchResult`**

- `pairs::DataFrame` with columns:
  - `:t_id` (index into treated rows, 1-based)
  - `:c_id` (index into control rows, 1-based)
  - `:distance` (km)
  - `:weight` (1.0 for KNN; kernel weights for `"kernel"`)
- `dropped_treated::Vector{Int}` indexes of treated units without valid matches
- `dropped_control::Vector{Int}` indexes of control units never used
- `config::MatchConfig` echo of settings
- `summary::Dict{String, Any}` basic counts and settings

------

#### `balance_table(res, treated, control, covars) -> DataFrame`

Returns a table with standardized mean differences before and after matching:

- Columns: `:var`, `:smd_before`, `:smd_after`
- For `"kernel"`, control side is weighted by the kernel weights; for `"geoNN"`, weights are 1 per selected pair (summing over reuse if `replace = true`).

------

#### `html_report(path, res, balance) -> String`

Writes a minimal static HTML file with a summary and the balance table. Returns the written `path`.

------

## Tests

A minimal test suite is included:

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test builds tiny treated/control tables, runs `geoNN` matching, and checks that:

- at least two pairs are formed,
- the balance table contains nonmissing SMD values.

------

## Design Notes

- Geometry creation: prefer `GeometryBasics.Point(x, y)` for lon-lat points. Prior `GeoInterface.point((x, y))` patterns may not work with recent `GeoInterface` versions.
- Distance: haversine is computed in kilometers on lon-lat. Projections and alternative metrics are not implemented in the MVP.
- Determinism: KNN matching is deterministic given a fixed dataset and constraints. Kernel matching uses deterministic weights for the k nearest valid controls.
- Performance: the MVP computes a dense distance matrix. For large `n` and `m`, consider tiling or approximate neighbors in future versions.

---

## Roadmap

- Optional CRS transforms and planar distances
- Exact calipers on covariates and propensity score
- Ratio matching and optimal pair assignment
- Richer diagnostics and plots
- Benchmarks and large-n strategies
