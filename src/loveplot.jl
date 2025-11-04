## src/loveplot.jl
## love plot for standardized mean differences (before vs after)

module LovePlot

using DataFrames
using CategoricalArrays    

function love_plot(bt::DataFrame; outfile::AbstractString = "balance_loveplot.png")
    try
        @eval using Plots
    catch
        error("Plots.jl is required for love_plot; run: pkg> add Plots")
    end

    # be tolerant to name types (Symbol/String) and case
    required = [:var, :smd_before, :smd_after]
    colsyms = Symbol.(names(bt))
    missing_cols = [c for c in required if !(c in colsyms)]
    !isempty(missing_cols) && error("balance table missing columns: $(missing_cols)")

    # build long-format data
    df = DataFrame(
        var  = repeat(bt.var, inner = 2),
        smd  = vcat(abs.(bt.smd_before), abs.(bt.smd_after)),
        when = vcat(fill("before", nrow(bt)), fill("after", nrow(bt))),
    )

    # order variables by max |SMD| to show largest first
    order  = sortperm(max.(abs.(bt.smd_before), abs.(bt.smd_after)); rev = true)
    levels = bt.var[order]

    # plot
    plt = plot(
        legend = :bottomright,
        xlabel = "abs standardized mean difference",
        ylabel = "",
        size   = (800, 400 + 18 * length(levels)),
    )

    for w in ("before", "after")
        sub = df[df.when .== w, :]
        sub.var = categorical(sub.var, levels = levels)
        scatter!(sub.smd, string.(sub.var); label = w, markersize = 6)
    end

    vline!([0.1, 0.2]; linestyle = :dash, label = ["0.1 threshold" "0.2 threshold"])
    png(plt, outfile)
    return outfile
end

end   ## module
