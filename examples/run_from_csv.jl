## examples/run_from_csv.jl
## read csvs -> build geometry column -> call match_spatial -> write pairs.csv and balance.html (+ optional love plot)

using CSV
using DataFrames
using GeometryBasics
using GeoMatch

const _HAVE_PLOTS = try
    @eval using Plots
    true
catch
    false
end

const _HAVE_PT = try
    @eval using PrettyTables
    true
catch
    false
end

## helpers
function normalize_name(x)::Symbol
    s = String(x)
    s = replace(s, '\uFEFF' => "", '\u00A0' => " ")
    s = strip(lowercase(s))
    s = replace(s, r"\p{Zs}+" => " ")
    s = replace(s, r"[ \t]+" => "_")
    return Symbol(s)
end

function force_symbol_names!(df::DataFrame)
    # ensure all column names are Symbols
    if any(n -> n isa String, names(df))
        rename!(df, Dict(n => Symbol(n) for n in names(df)))
    end
end

function first_existing(paths::Vector{String})
    for p in paths
        if isfile(p)
            return p
        end
    end
    return paths[1]
end

function _write_balance_html(df::DataFrame, outfile::AbstractString)
    if _HAVE_PT
        io = open(outfile, "w")
        try
            PrettyTables.pretty_table(
                io, df;
                backend = :html,
                stand_alone = true,       
                wrap_table_in_div = true   
            )
        finally
            close(io)
        end
    else
        CSV.write(replace(outfile, ".html" => ".csv"), df)
    end
end


## env overrides (empty means auto resolve)
T_CSV   = get(ENV, "T_CSV",   "")
C_CSV   = get(ENV, "C_CSV",   "")
GEO_LON = get(ENV, "GEO_LON", "lon")
GEO_LAT = get(ENV, "GEO_LAT", "lat")
COVARS  = split(get(ENV, "COVARS", ""), ',')
METHOD  = get(ENV, "METHOD", "geoNN")   # "geoNN" or "kernel"
K       = parse(Int, get(ENV, "K", "1"))
RADIUS  = let s = get(ENV, "RADIUS_KM", ""); isempty(s) ? nothing : parse(Float64, s) end
BANDW   = let s = get(ENV, "BANDWIDTH", ""); isempty(s) ? nothing : parse(Float64, s) end
REPLACE = get(ENV, "REPLACE", "false") in ("1","true","yes")

if isempty(T_CSV)
    T_CSV = first_existing(["data/treated.csv", "treated.csv"])
end
if isempty(C_CSV)
    C_CSV = first_existing(["data/control.csv", "control.csv"])
end

println("reading CSVs...")
t = CSV.read(T_CSV, DataFrame)
c = CSV.read(C_CSV, DataFrame)

# normalize dataframe column names (content-wise) then force Symbol names
rename!(t, Dict(n => normalize_name(n) for n in names(t)))
rename!(c, Dict(n => normalize_name(n) for n in names(c)))
force_symbol_names!(t)
force_symbol_names!(c)

# normalize user-provided names in the same way
lon_col = normalize_name(GEO_LON)
lat_col = normalize_name(GEO_LAT)
cov_syms = isempty(COVARS[1]) ? Symbol[] : [normalize_name(v) for v in COVARS]

# debug
println("treated propertynames = ", propertynames(t))
println("control propertynames = ", propertynames(c))
println("expect lon=", lon_col, " lat=", lat_col, " covars=", cov_syms)

# validate required columns exist 
@assert lon_col in propertynames(t) "## treated is missing column $(GEO_LON)"
@assert lat_col in propertynames(t) "## treated is missing column $(GEO_LAT)"
@assert lon_col in propertynames(c) "## control is missing column $(GEO_LON)"
@assert lat_col in propertynames(c) "## control is missing column $(GEO_LAT)"

# build geometry columns
t.geometry = Point2f.(Float32.(t[!, lon_col]), Float32.(t[!, lat_col]))
c.geometry = Point2f.(Float32.(c[!, lon_col]), Float32.(c[!, lat_col]))

if METHOD == "kernel"
    @assert BANDW !== nothing "## bandwidth is required for METHOD=kernel"
end

println("running match_spatial...")
res = GeoMatch.match_spatial(
    t, c;
    covars = cov_syms,
    method = METHOD,
    k = K,
    radius_km = RADIUS,
    kernel = "gaussian",
    bandwidth = BANDW,
    replace = REPLACE,
)

pairs_out = "pairs.csv"
CSV.write(pairs_out, res.pairs)
println("wrote: $(pairs_out)")

if !isempty(cov_syms)
    bt = GeoMatch.balance_table(res, t, c, cov_syms)
    _write_balance_html(bt, "balance.html")
    println("wrote: balance.html")
    if _HAVE_PLOTS
        try
            GeoMatch.love_plot(bt; outfile = "balance_loveplot.png")
            println("wrote: balance_loveplot.png")
        catch e
            @warn "love_plot failed, skipping" exception = e
        end
    else
        println("plots.jl not installed; skip love plot (pkg> add Plots to enable)")
    end
else
    println("no covariates provided; skip balance table and love plot")
end

println("done.")
